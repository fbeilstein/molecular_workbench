from typing import Dict, List, Tuple
import numpy as np
from pyscf import gto
from qm.analyzers.planar_pi import extract_step0_planar_pi
from qm.analyzers.valence import ValenceLocalizer
from qm.classification import count_core_orbitals
from qm.renderer import Renderer
from qm.models import Orbital

HARTREE_TO_EV = 27.2114


def _frontier_label_and_key(mo_idx: int, homo_idx: int, lumo_idx: int) -> Tuple[str, str]:
    if mo_idx == homo_idx:
        return "HOMO", "homo"
    if mo_idx < homo_idx:
        d = homo_idx - mo_idx
        return f"HOMO-{d}", f"homo_minus_{d}"
    if mo_idx == lumo_idx:
        return "LUMO", "lumo"
    d = mo_idx - lumo_idx
    return f"LUMO+{d}", f"lumo_plus_{d}"


def _group_degenerate_clusters(indices: List[int], mo_energy: np.ndarray, tol_ha: float = 1e-4) -> List[List[int]]:
    if not indices:
        return []
    clusters: List[List[int]] = [[indices[0]]]
    for idx in indices[1:]:
        if abs(float(mo_energy[idx]) - float(mo_energy[clusters[-1][-1]])) < tol_ha:
            clusters[-1].append(idx)
        else:
            clusters.append([idx])
    return clusters


def _align_degenerate_canonical_subspaces(
    can_coeff: np.ndarray,
    mo_energy: np.ndarray,
    ovlp: np.ndarray,
    chem_orbitals: List[Orbital],
    n_core: int,
    nocc: int,
    threshold: float = 0.90,
) -> np.ndarray:
    """
    Rotate strictly degenerate canonical MO subspaces (|Δε| < 1e-4 Ha, dim >= 2)
    so their orthonormal basis aligns one-to-one with matching chemical orbitals
    (e.g. symmetric benzene e1g/e2u or linear ethyne pi_u/pi_g).
    Every rotated basis vector remains an exact canonical Fock eigenfunction.
    """
    C_aligned = can_coeff.copy()
    n_mo = C_aligned.shape[1]

    occ_clusters = _group_degenerate_clusters(list(range(n_core, nocc)), mo_energy)
    virt_clusters = _group_degenerate_clusters(list(range(nocc, n_mo)), mo_energy)

    for idxs in occ_clusters + virt_clusters:
        g = len(idxs)
        if g < 2:
            continue
        is_virt = idxs[0] >= nocc
        U0 = C_aligned[:, idxs]  # nao x g, S-orthonormal

        matching_vecs = []
        for orb in chem_orbitals:
            orb_is_anti = ("star" in orb.type) or ("*" in orb.type)
            if orb_is_anti != is_virt:
                continue
            proj = U0.T @ ovlp @ orb.mo_coeff
            w_sub = float(np.dot(proj, proj))
            if w_sub >= threshold:
                matching_vecs.append((orb.energy_ev if orb.energy_ev is not None else 0.0, proj))

        if not matching_vecs:
            continue

        matching_vecs.sort(key=lambda item: item[0])
        selected_projs = []
        for _, proj in matching_vecs:
            cand = np.column_stack(selected_projs + [proj]) if selected_projs else proj.reshape(-1, 1)
            gram = cand.T @ cand
            if float(np.min(np.linalg.eigvalsh(gram))) > 1e-3:
                selected_projs.append(proj)
            if len(selected_projs) == g:
                break

        if not selected_projs:
            continue

        A = np.column_stack(selected_projs)  # g x p
        ev, V = np.linalg.eigh(A.T @ A)
        Q_match = A @ (V @ np.diag(1.0 / np.sqrt(np.maximum(ev, 1e-8))) @ V.T)  # g x p

        p = Q_match.shape[1]
        if p < g:
            P_rem = np.eye(g) - Q_match @ Q_match.T
            u_svd, _, _ = np.linalg.svd(P_rem)
            Q_perp = u_svd[:, : (g - p)]
            Q = np.column_stack([Q_match, Q_perp]) if is_virt else np.column_stack([Q_perp, Q_match])
        else:
            Q = Q_match

        C_aligned[:, idxs] = U0 @ Q

    return C_aligned


class QuantumPipeline:
    def __init__(
        self,
        mol: gto.M,
        mf,
        smiles: str,
        name: str,
        out_dir: str,
        grid_points: int = 50,
        frontier_depth: int = 0,
        snap_threshold: float = 0.90,
    ):
        self.mol = mol
        self.mf = mf
        self.smiles = smiles
        self.name = name
        self.out_dir = out_dir
        self.grid_points = grid_points
        self.frontier_depth = max(0, int(frontier_depth))
        self.snap_threshold = float(snap_threshold)
        self.orbitals: List[Orbital] = []

    def run(self):
        print(f"--- Running Unified Quantum Pipeline for {self.name} (frontier_depth={self.frontier_depth}) ---")

        all_atom_ids = [f"{self.mol.atom_symbol(i)}{i}" for i in range(self.mol.natm)]
        nocc = self.mol.nelectron // 2 if self.mol.spin == 0 else self.mol.nelec[0]

        raw_coeff = (
            self.mf.mo_coeff[0]
            if (self.mol.spin > 0 and isinstance(self.mf.mo_coeff, (tuple, list, np.ndarray)) and np.ndim(self.mf.mo_coeff) == 3)
            else self.mf.mo_coeff
        )
        mo_energy = (
            self.mf.mo_energy[0]
            if (self.mol.spin > 0 and isinstance(self.mf.mo_energy, (tuple, list, np.ndarray)) and np.ndim(self.mf.mo_energy) == 2)
            else self.mf.mo_energy
        )
        raw_fock = self.mf.get_fock()
        fock = (
            raw_fock[0]
            if (self.mol.spin > 0 and isinstance(raw_fock, (tuple, list, np.ndarray)) and np.ndim(raw_fock) == 3)
            else raw_fock
        )
        ovlp = self.mol.intor_symmetric("int1e_ovlp")
        n_mo = raw_coeff.shape[1]
        n_core = count_core_orbitals(self.mol, mo_energy, nocc)

        # 1. Step-0 Planar Conjugated π Subspace Separation (unified engine shared with server_irc.py)
        C_non_pi_occ, C_non_pi_virt, pi_meta, pi_vecs = extract_step0_planar_pi(
            self.mol, self.mf, atom_labels=all_atom_ids, smiles=self.smiles
        )

        chem_orbitals: List[Orbital] = []
        for key, meta in pi_meta.items():
            vec = pi_vecs[key].copy()
            nrm = np.sqrt(float(vec.T @ ovlp @ vec))
            if nrm > 1e-8:
                vec /= nrm
            orb_type = "delocalized_pi" if meta.get("type") == "delocalized_pi" else "delocalized_pistar"
            chem_orbitals.append(
                Orbital(
                    type=orb_type,
                    atoms=meta.get("atoms", []),
                    mo_coeff=vec,
                    energy_ev=meta.get("energy_ev"),
                    canonical_label=meta.get("label"),
                )
            )

        # 2. Valence Localizer (IBO on non-pi subspaces) with diagonal Fock energies on ALL orbitals
        valence_loc = ValenceLocalizer(self.mol, raw_coeff, fock=fock)
        for orb in valence_loc.localize(C_non_pi_occ, C_non_pi_virt):
            vec = orb.mo_coeff.copy()
            nrm = np.sqrt(float(vec.T @ ovlp @ vec))
            if nrm > 1e-8:
                vec /= nrm
            orb.mo_coeff = vec
            chem_orbitals.append(orb)

        # 3. Align strictly degenerate canonical subspaces (|Δε| < 1e-4 Ha) before Frontier extraction & overlap check
        can_coeff = _align_degenerate_canonical_subspaces(
            raw_coeff, mo_energy, ovlp, chem_orbitals, n_core, nocc, threshold=self.snap_threshold
        )

        # 4. Frontier Orbitals (HOMO-d .. HOMO, LUMO .. LUMO+d; default d=0 -> single HOMO & LUMO)
        homo_idx = nocc - 1
        lumo_idx = nocc
        frontier_orbitals: List[Orbital] = []
        frontier_by_mo: Dict[int, Orbital] = {}

        occ_start = max(n_core, homo_idx - self.frontier_depth)
        for m in range(occ_start, homo_idx + 1):
            lbl, fkey = _frontier_label_and_key(m, homo_idx, lumo_idx)
            orb_type = "homo" if m == homo_idx else "frontier"
            f_orb = Orbital(
                type=orb_type,
                atoms=[],
                mo_coeff=can_coeff[:, m].copy(),
                energy_ev=round(float(mo_energy[m] * HARTREE_TO_EV), 2),
                canonical_label=lbl,
                frontier_key=fkey,
            )
            frontier_orbitals.append(f_orb)
            frontier_by_mo[m] = f_orb

        virt_end = min(n_mo, lumo_idx + self.frontier_depth + 1)
        for m in range(lumo_idx, virt_end):
            lbl, fkey = _frontier_label_and_key(m, homo_idx, lumo_idx)
            orb_type = "lumo" if m == lumo_idx else "frontier"
            f_orb = Orbital(
                type=orb_type,
                atoms=[],
                mo_coeff=can_coeff[:, m].copy(),
                energy_ev=round(float(mo_energy[m] * HARTREE_TO_EV), 2),
                canonical_label=lbl,
                frontier_key=fkey,
            )
            frontier_orbitals.append(f_orb)
            frontier_by_mo[m] = f_orb

        # 5. Post-Localization Overlap-Gated Snap ("Snap to Grid" when w >= 0.90)
        # Check overlap against canonical frontier manifold (up to max(8, frontier_depth) around HOMO/LUMO).
        max_tag_depth = max(8, self.frontier_depth)
        occ_cand = list(range(max(n_core, homo_idx - max_tag_depth), nocc))
        virt_cand = list(range(lumo_idx, min(n_mo, lumo_idx + max_tag_depth + 1)))

        for orb in chem_orbitals:
            if orb.type == "core":
                continue
            is_anti = ("star" in orb.type) or ("*" in orb.type)
            cand_indices = virt_cand if is_anti else occ_cand
            if not cand_indices:
                continue

            C_cand = can_coeff[:, cand_indices]
            overlaps = C_cand.T @ ovlp @ orb.mo_coeff
            weights = overlaps ** 2
            best_rel = int(np.argmax(weights))
            best_w = float(weights[best_rel])

            if best_w >= self.snap_threshold:
                best_m = cand_indices[best_rel]
                sign = 1.0 if float(overlaps[best_rel]) >= 0.0 else -1.0
                snapped_vec = (sign * can_coeff[:, best_m]).copy()
                orb.mo_coeff = snapped_vec
                orb.energy_ev = round(float(mo_energy[best_m] * HARTREE_TO_EV), 2)
                tag, _ = _frontier_label_and_key(best_m, homo_idx, lumo_idx)
                orb.frontier_tag = tag
                if best_m in frontier_by_mo:
                    frontier_by_mo[best_m].mo_coeff = snapped_vec.copy()

        self.orbitals = frontier_orbitals + chem_orbitals

        # 6. Render
        renderer = Renderer(self.mol, self.name, self.out_dir, self.grid_points)
        manifest = renderer.render(self.orbitals)

        return manifest
