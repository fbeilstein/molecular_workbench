from typing import List, Optional, Tuple
import numpy as np
from pyscf import gto, lo
from qm.constants import HARTREE_TO_EV
from qm.classification import (
    _atom_populations,
    _canonical_atom_sort_key,
    _classify_orbital,
    _classify_virtual_orbital,
)
from qm.models import Orbital


class ValenceLocalizer:
    """Unified Occupied + Valence Virtual IBO Localizer shared by server.py and server_irc.py."""

    def __init__(
        self,
        mol: gto.M,
        mo_coeff: np.ndarray,
        fock: Optional[np.ndarray] = None,
        atom_labels: Optional[List[str]] = None,
        occ_all_coeff: Optional[np.ndarray] = None,
    ):
        self.mol = mol
        self.mo_coeff = mo_coeff
        self.fock = fock
        self.ovlp = mol.intor_symmetric("int1e_ovlp")
        self.atom_labels = [mol.atom_symbol(i) for i in range(mol.natm)]
        self.atom_ids = (
            list(atom_labels)
            if (atom_labels is not None and len(atom_labels) == mol.natm)
            else [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
        )
        if occ_all_coeff is not None:
            self.occ_all_coeff = occ_all_coeff
        else:
            nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
            self.occ_all_coeff = mo_coeff[:, :nocc] if mo_coeff is not None else None

    def _eval_energy_ev(self, vec: np.ndarray) -> Optional[float]:
        if self.fock is None:
            return None
        denom = float(vec.T @ self.ovlp @ vec)
        if denom < 1e-8:
            return None
        return round(float((vec.T @ self.fock @ vec) / denom) * HARTREE_TO_EV, 2)

    def _iao_atom_populations(
        self,
        C_vecs: np.ndarray,
        Q_iao: Optional[np.ndarray],
        ref_ao_atom: Optional[List[int]],
    ) -> np.ndarray:
        """Compute non-negative fractional atom populations (each row sums to 1.0)."""
        n_orb = C_vecs.shape[1]
        pops = np.zeros((n_orb, self.mol.natm))
        if n_orb == 0:
            return pops
        if Q_iao is not None and ref_ao_atom is not None:
            C_in_iao = Q_iao.T @ self.ovlp @ C_vecs
            sq = C_in_iao ** 2
            for mu, a_idx in enumerate(ref_ao_atom):
                pops[:, a_idx] += sq[mu, :]
        else:
            ev_s, U_s = np.linalg.eigh(self.ovlp)
            S_half = U_s @ np.diag(np.sqrt(np.maximum(ev_s, 1e-12))) @ U_s.T
            C_low = S_half @ C_vecs
            sq = C_low ** 2
            ao_atom = [int(l.split()[0]) for l in self.mol.ao_labels()]
            for mu, a_idx in enumerate(ao_atom):
                pops[:, a_idx] += sq[mu, :]
        row_sums = pops.sum(axis=1, keepdims=True)
        pops = np.where(row_sums > 1e-10, pops / np.maximum(row_sums, 1e-10), pops)
        return pops

    def _purify_multiple_bonds(
        self,
        C_loc: np.ndarray,
        frac_pops: np.ndarray,
        Q_iao: Optional[np.ndarray],
        ref_ao_atom: Optional[List[int]],
        is_virt: bool = False,
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Multiple-Bond σ/π (and σ*/π*) Subspace Fock Purification.

        Eliminates bent "banana-bond" mixing on double/triple bonds by diagonalizing
        the m x m Fock block strictly within each heavy-atom pair's localized bond subspace.
        """
        if self.fock is None or C_loc.shape[1] < 2:
            return C_loc, frac_pops

        C_out = C_loc.copy()
        pops_out = frac_pops.copy()
        n_orb = C_out.shape[1]

        pair_groups = {}
        for k in range(n_orb):
            top2 = np.argsort(-pops_out[k])[:2]
            if len(top2) < 2:
                continue
            a1, a2 = int(min(top2)), int(max(top2))
            if (
                self.mol.atom_symbol(a1) != "H"
                and self.mol.atom_symbol(a2) != "H"
                and pops_out[k, a1] >= 0.15
                and pops_out[k, a2] >= 0.15
                and (pops_out[k, a1] + pops_out[k, a2]) >= 0.72
            ):
                pair_groups.setdefault((a1, a2), []).append(k)

        coords_b = self.mol.atom_coords()
        for (a1, a2), cols in pair_groups.items():
            if len(cols) < 2:
                continue
            C_pair = C_out[:, cols]
            # 1. Exact S_ao Löwdin orthonormalization within the pair subspace
            S_sub = C_pair.T @ self.ovlp @ C_pair
            ev_s, U_s = np.linalg.eigh(S_sub)
            C_pair = C_pair @ (U_s @ np.diag(1.0 / np.sqrt(np.maximum(ev_s, 1e-10))) @ U_s.T)

            # 2. Diagonalize m x m Fock block:
            #    Occupied: lowest root = pure σ, upper roots = pure π
            #    Virtual:  lower roots = pure π*, highest root = pure σ*
            F_sub = C_pair.T @ self.fock @ C_pair
            ev_f, U_f = np.linalg.eigh(F_sub)
            C_purified = C_pair @ U_f

            # 3. Deterministic geometric gauge for 2D degenerate triple-bond π subspace (m == 3)
            #    Pick the nearest non-collinear heavy neighbor (>7° off bond axis) so a nearly
            #    linear bond (e.g. acylium ∠C5-C6-O8 ≈ 178°) never divides by sin(2°) noise.
            if C_purified.shape[1] == 3:
                pi_sl = slice(0, 2) if is_virt else slice(1, 3)
                dE_pi = abs(float(ev_f[pi_sl.stop - 1]) - float(ev_f[pi_sl.start])) * HARTREE_TO_EV
                if dE_pi < 3.0:
                    bond_v = coords_b[a2] - coords_b[a1]
                    bond_hat = bond_v / max(float(np.linalg.norm(bond_v)), 1e-8)
                    best_nb, best_d = None, 999.0
                    for nb in range(self.mol.natm):
                        if nb not in (a1, a2) and self.mol.atom_symbol(nb) != "H":
                            rv = coords_b[nb] - coords_b[a1]
                            rv_n = float(np.linalg.norm(rv))
                            if rv_n < 1e-6:
                                continue
                            sin_ang = float(np.linalg.norm(np.cross(bond_hat, rv / rv_n)))
                            if sin_ang < 0.12:
                                continue
                            d_nb = min(rv_n, float(np.linalg.norm(coords_b[nb] - coords_b[a2])))
                            if d_nb < best_d:
                                best_d, best_nb = d_nb, nb
                    ref_v = (
                        coords_b[best_nb] - coords_b[a1]
                        if best_nb is not None
                        else np.array([1.0, 0.0, 0.0])
                    )
                    plane_norm = np.cross(bond_hat, ref_v)
                    if float(np.linalg.norm(plane_norm)) < 1e-4:
                        plane_norm = np.cross(bond_hat, np.array([0.0, 1.0, 0.0]))
                    plane_norm /= np.linalg.norm(plane_norm)
                    probe_pt = (0.5 * (coords_b[a1] + coords_b[a2]) + 0.5 * plane_norm).reshape(1, 3)
                    if is_virt:
                        probe_pt = (
                            0.25 * coords_b[a1] + 0.75 * coords_b[a2] + 0.5 * plane_norm
                        ).reshape(1, 3)
                    vals_n = (self.mol.eval_gto("GTOval", probe_pt) @ C_purified[:, pi_sl])[0]
                    theta = float(np.arctan2(vals_n[0], vals_n[1]))
                    R2 = np.array([[np.cos(theta), np.sin(theta)], [-np.sin(theta), np.cos(theta)]])
                    C_purified[:, pi_sl] = C_purified[:, pi_sl] @ R2

            C_out[:, cols] = C_purified
            pops_out[cols, :] = self._iao_atom_populations(C_purified, Q_iao, ref_ao_atom)

        return C_out, pops_out

    def _extract_vvo_subspace(
        self,
        virt_coeff: np.ndarray,
        Q_iao: Optional[np.ndarray],
    ) -> np.ndarray:
        """Extract exact minimal-valence virtual subspace (VVOs) from virt_coeff via IAO SVD.

        In minimal bases (sto-3g), nao == n_iao so the full non-pi virtual space is already
        exact valence. In extended bases (e.g. 6-31g*), SVD projection onto Q_iao separates
        true chemical antibonds (s ~ 1.0) from diffuse/polarization Rydberg virtuals (s ~ 0.0).
        """
        if virt_coeff is None or virt_coeff.shape[1] == 0:
            return np.zeros((self.mol.nao, 0))
        if Q_iao is None or self.mol.nao <= Q_iao.shape[1]:
            # Minimal basis: ensure exact S-orthonormality
            S_v = virt_coeff.T @ self.ovlp @ virt_coeff
            ev, U = np.linalg.eigh(S_v)
            keep = ev > 1e-8
            return virt_coeff @ (U[:, keep] / np.sqrt(ev[keep]))

        M = Q_iao.T @ self.ovlp @ virt_coeff
        _, s_vals, Vt = np.linalg.svd(M, full_matrices=False)
        n_vvo = int(np.sum(s_vals > 0.50))
        if self.occ_all_coeff is not None:
            max_vvo = max(0, Q_iao.shape[1] - self.occ_all_coeff.shape[1])
            n_vvo = min(n_vvo, max_vvo)
        if n_vvo <= 0:
            return np.zeros((self.mol.nao, 0))

        C_vvo = virt_coeff @ Vt[:n_vvo].T
        S_v = C_vvo.T @ self.ovlp @ C_vvo
        ev, U = np.linalg.eigh(S_v)
        keep = ev > 1e-8
        return C_vvo @ (U[:, keep] / np.sqrt(ev[keep]))

    def localize(self, cols_to_localize, virt_cols) -> List[Orbital]:
        """
        Localize non-conjugated occupied valence MOs AND valence virtual MOs (VVOs)
        using Knizia's IAO/IBO localization + multiple-bond Fock purification +
        intra-atomic LP Fock purification.
        """
        orbitals: List[Orbital] = []
        occ_coeff = None
        if cols_to_localize is not None:
            if isinstance(cols_to_localize, np.ndarray) and cols_to_localize.ndim == 2:
                if cols_to_localize.shape[1] > 0:
                    occ_coeff = cols_to_localize
            elif isinstance(cols_to_localize, (list, tuple)) and len(cols_to_localize) > 0:
                occ_coeff = self.mo_coeff[:, cols_to_localize]

        virt_coeff = None
        if virt_cols is not None:
            if isinstance(virt_cols, np.ndarray) and virt_cols.ndim == 2:
                if virt_cols.shape[1] > 0:
                    virt_coeff = virt_cols
            elif isinstance(virt_cols, (list, tuple)) and len(virt_cols) > 0:
                virt_coeff = self.mo_coeff[:, virt_cols]

        if occ_coeff is None and virt_coeff is None:
            return orbitals

        try:
            # Build shared full-molecule Knizia IAOs once for both occupied and virtual subspaces
            iao_source = (
                self.occ_all_coeff
                if (self.occ_all_coeff is not None and self.occ_all_coeff.shape[1] > 0)
                else occ_coeff
            )
            iao_raw = None
            Q_iao = None
            ref_ao_atom = None
            if iao_source is not None and iao_source.shape[1] > 0:
                iao_raw = lo.iao.iao(self.mol, iao_source)
                Q_iao = lo.vec_lowdin(iao_raw, self.ovlp)
                ref_mol = lo.iao.reference_mol(self.mol)
                ref_ao_atom = [int(l.split()[0]) for l in ref_mol.ao_labels()]

            # ------------------------------------------------------------
            # 1. Occupied Valence IBOs + Multiple-Bond & LP Fock Purification
            # ------------------------------------------------------------
            if occ_coeff is not None and occ_coeff.shape[1] > 0:
                loc_occ = lo.ibo.ibo(self.mol, occ_coeff, iaos=iao_raw, verbose=0)
                occ_frac_pops = self._iao_atom_populations(loc_occ, Q_iao, ref_ao_atom)
                loc_occ, occ_frac_pops = self._purify_multiple_bonds(
                    loc_occ, occ_frac_pops, Q_iao, ref_ao_atom, is_virt=False
                )

                seen_pi = set()
                for i in range(loc_occ.shape[1]):
                    mo = loc_occ[:, i]
                    pop = occ_frac_pops[i] * 2.0
                    info = _classify_orbital(
                        self.mol, mo, pop, self.atom_labels, self.atom_ids, self.ovlp
                    )
                    orb_type = info.get("type", "unknown")
                    atoms = info.get("atoms", [info.get("atom", "?")])
                    if len(atoms) >= 2:
                        atoms = sorted(atoms, key=_canonical_atom_sort_key)

                    # Limit hypervalent M=O bonds to at most 1 pi bond
                    if orb_type == "pi" and len(atoms) == 2:
                        pair = tuple(sorted(atoms))
                        s1 = "".join(filter(str.isalpha, pair[0]))
                        s2 = "".join(filter(str.isalpha, pair[1]))
                        hyper = {"S", "P", "Cl", "Br", "I", "N"}
                        if (s1 in hyper and s2 == "O") or (s2 in hyper and s1 == "O"):
                            if pair in seen_pi:
                                orb_type = "lone_pair"
                                atoms = [pair[0] if s1 == "O" else pair[1]]
                            else:
                                seen_pi.add(pair)

                    orbitals.append(
                        Orbital(
                            type=orb_type,
                            atoms=atoms,
                            mo_coeff=mo.copy(),
                            energy_ev=self._eval_energy_ev(mo),
                        )
                    )

                # Intra-atomic Lone-Pair Subspace Fock Purification ("No Rabbit Ears")
                if self.fock is not None:
                    lp_groups = {}
                    for idx_orb, orb in enumerate(orbitals):
                        if orb.type == "lone_pair" and orb.atoms:
                            lp_groups.setdefault(orb.atoms[0], []).append(idx_orb)

                    for _, orb_indices in lp_groups.items():
                        if len(orb_indices) < 2:
                            continue
                        C_sub = np.column_stack([orbitals[k].mo_coeff for k in orb_indices])
                        S_sub = C_sub.T @ self.ovlp @ C_sub
                        ev_s, U_s = np.linalg.eigh(S_sub)
                        C_orth = C_sub @ (
                            U_s @ np.diag(1.0 / np.sqrt(np.maximum(ev_s, 1e-10))) @ U_s.T
                        )
                        F_sub = C_orth.T @ self.fock @ C_orth
                        ev_f, U_f = np.linalg.eigh(F_sub)
                        C_pur = C_orth @ U_f
                        for j, k in enumerate(orb_indices):
                            orbitals[k].mo_coeff = C_pur[:, j].copy()
                            orbitals[k].energy_ev = round(float(ev_f[j] * HARTREE_TO_EV), 2)

            # ------------------------------------------------------------
            # 2. Orthonormal Valence Virtual Orbitals (VVOs) via IAO + IBO
            #    Note: Keep raw p=4 IBO virtuals diabatic (do NOT diagonalize
            #    2x2 Fock blocks across π* pairs, which forces adiabatic
            #    avoided-crossing 90° swaps when a reactive π* climbs past a
            #    spectator π* into σ*). Raw p=4 IBOs already keep σ* pure.
            # ------------------------------------------------------------
            if virt_coeff is not None and virt_coeff.shape[1] > 0:
                C_vvo = self._extract_vvo_subspace(virt_coeff, Q_iao)
                if C_vvo.shape[1] > 0:
                    loc_vvo = lo.ibo.ibo(self.mol, C_vvo, iaos=iao_raw, verbose=0)
                    vvo_frac_pops = self._iao_atom_populations(loc_vvo, Q_iao, ref_ao_atom)

                    for i in range(loc_vvo.shape[1]):
                        v_mo = loc_vvo[:, i]
                        v_info = _classify_virtual_orbital(
                            self.mol, v_mo, vvo_frac_pops[i], self.atom_labels, self.atom_ids, self.ovlp
                        )
                        v_type = v_info.get("type", "sigma_star")
                        v_atoms = v_info.get("atoms", [v_info.get("atom", "?")])
                        can_lbl = f"p*({v_atoms[0]})" if (v_type == "p_star" and v_atoms) else None
                        orbitals.append(
                            Orbital(
                                type=v_type,
                                atoms=v_atoms,
                                mo_coeff=v_mo.copy(),
                                energy_ev=self._eval_energy_ev(v_mo),
                                canonical_label=can_lbl,
                                sym_score=1.0,
                            )
                        )

        except Exception as e:
            print(f"Valence localization failed: {e}")

        return orbitals
