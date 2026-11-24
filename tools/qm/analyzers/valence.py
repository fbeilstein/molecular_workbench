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

    def _decouple_virtual_sigma_star(
        self,
        C_loc: np.ndarray,
        frac_pops: np.ndarray,
        Q_iao: Optional[np.ndarray],
        ref_ao_atom: Optional[List[int]],
    ) -> Tuple[np.ndarray, np.ndarray]:
        """Decouple pure σ* from π* on 2-atom virtual multiple-bond groups via minimal rotation.

        On homonuclear multiple bonds (e.g. N≡N, C≡C), σ* and π* share identical 50/50 IAO
        populations, leaving the Pipek-Mezey IBO objective invariant to σ*/π* mixing.
        Projecting out ONLY the highest-energy σ* Fock eigenvector via Löwdin minimal-rotation
        projection yields pure σ* and pure π* while preserving diabatic π* orientation
        (no 2x2 Fock diagonalization inside the π* plane).
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

        for _, cols in pair_groups.items():
            m = len(cols)
            if m < 2 or m > 3:
                continue
            C_pair = C_out[:, cols]
            S_sub = C_pair.T @ self.ovlp @ C_pair
            ev_s, U_s = np.linalg.eigh(S_sub)
            C_pair = C_pair @ (U_s @ np.diag(1.0 / np.sqrt(np.maximum(ev_s, 1e-10))) @ U_s.T)

            F_sub = C_pair.T @ self.fock @ C_pair
            ev_f, U_f = np.linalg.eigh(F_sub)
            gap_ev = (float(ev_f[-1]) - float(ev_f[-2])) * HARTREE_TO_EV
            if gap_ev < 3.0:
                continue

            u_sig = U_f[:, -1].copy()
            j_sig = int(np.argmax(np.abs(u_sig)))
            if u_sig[j_sig] < 0.0:
                u_sig = -u_sig

            # Skip if σ* is already decoupled to machine precision
            if float(abs(u_sig[j_sig])) > 1.0 - 1e-6:
                continue

            keep = [j for j in range(m) if j != j_sig]
            I_keep = np.eye(m)[:, keep]
            V = I_keep - np.outer(u_sig, u_sig[keep])
            G = V.T @ V
            ev_g, U_g = np.linalg.eigh(G)
            V_orth = V @ (U_g @ np.diag(1.0 / np.sqrt(np.maximum(ev_g, 1e-10))) @ U_g.T)

            U_rot = np.zeros((m, m))
            U_rot[:, j_sig] = u_sig
            for idx_k, orig_j in enumerate(keep):
                U_rot[:, orig_j] = V_orth[:, idx_k]

            C_decoupled = C_pair @ U_rot
            C_out[:, cols] = C_decoupled
            pops_out[cols, :] = self._iao_atom_populations(C_decoupled, Q_iao, ref_ao_atom)

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

                # Cumulated 3c-4e Coplanar π / Lone-Pair Decoupling (e.g. O=C=O):
                # A central atom C only has one 2p orbital per π plane, so two coplanar π-classified
                # IBOs on (C-A) and (C-B) represent 1 covalent π bond + 1 terminal lone pair LP,
                # never two coplanar π bonds (which would create triple bonds on both A and B).
                self._decouple_coplanar_3c4e_pi(orbitals, Q_iao, ref_ao_atom)

            # ------------------------------------------------------------
            # 2. Orthonormal Valence Virtual Orbitals (VVOs) via IAO + IBO
            #    Note: Keep π* virtuals diabatic (do NOT diagonalize 2x2 Fock
            #    blocks inside a degenerate π* pair, which causes adiabatic
            #    avoided-crossing swaps during reactions). Only decouple the
            #    highest-energy σ* Fock eigenvector from π* on 2-atom multiple
            #    bonds via minimal-rotation orthogonal projection (fixes 50/50
            #    homonuclear σ*/π* mixing such as N≡N without rotating π*).
            # ------------------------------------------------------------
            if virt_coeff is not None and virt_coeff.shape[1] > 0:
                C_vvo = self._extract_vvo_subspace(virt_coeff, Q_iao)
                if C_vvo.shape[1] > 0:
                    loc_vvo = lo.ibo.ibo(self.mol, C_vvo, iaos=iao_raw, verbose=0)
                    vvo_frac_pops = self._iao_atom_populations(loc_vvo, Q_iao, ref_ao_atom)
                    loc_vvo, vvo_frac_pops = self._decouple_virtual_sigma_star(
                        loc_vvo, vvo_frac_pops, Q_iao, ref_ao_atom
                    )

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
                    self._align_coplanar_pistar_atoms(orbitals)

        except Exception as e:
            print(f"Valence localization failed: {e}")

        return orbitals

    def _decouple_coplanar_3c4e_pi(
        self,
        orbitals: List[Orbital],
        Q_iao: Optional[np.ndarray],
        ref_ao_atom: Optional[List[int]],
    ) -> None:
        """Decouple coplanar 3c-4e π IBO pairs sharing a central atom into 1 π bond + 1 lone pair."""
        ao_lbls = self.mol.ao_labels()
        pi_entries = []
        for idx_o, orb in enumerate(orbitals):
            if orb.type == "pi" and len(orb.atoms) == 2:
                try:
                    i0 = self.atom_ids.index(orb.atoms[0])
                    i1 = self.atom_ids.index(orb.atoms[1])
                    pi_entries.append((idx_o, i0, i1))
                except ValueError:
                    continue
        if len(pi_entries) < 2:
            return

        from qm.analyzers.planar_pi.geometry import detect_planar_pi_fragments
        frags = detect_planar_pi_fragments(self.mol)
        fallback_normal = frags[0][1] if frags else np.array([0.0, 1.0, 0.0])
        coords_bohr = self.mol.atom_coords()

        def _local_out_normal(c_i: int, a_i: int, b_i: int) -> np.ndarray:
            vA = coords_bohr[a_i] - coords_bohr[c_i]
            vB = coords_bohr[b_i] - coords_bohr[c_i]
            nA, nB = float(np.linalg.norm(vA)), float(np.linalg.norm(vB))
            if nA > 1e-6 and nB > 1e-6:
                cross = np.cross(vA / nA, vB / nB)
                cn = float(np.linalg.norm(cross))
                if cn > 0.05:
                    return cross / cn
            v_ax = coords_bohr[b_i] - coords_bohr[a_i]
            ax_n = float(np.linalg.norm(v_ax))
            if ax_n > 1e-6:
                v_ax = v_ax / ax_n
                heavy_others = [
                    i for i in range(self.mol.natm)
                    if i not in (c_i, a_i, b_i) and self.mol.atom_charge(i) > 1
                ]
                others = heavy_others or [i for i in range(self.mol.natm) if i not in (c_i, a_i, b_i)]
                if others:
                    nearest = min(others, key=lambda i: float(np.linalg.norm(coords_bohr[i] - coords_bohr[c_i])))
                    cross = np.cross(v_ax, coords_bohr[nearest] - coords_bohr[c_i])
                    cn = float(np.linalg.norm(cross))
                    if cn > 1e-4:
                        return cross / cn
            return fallback_normal

        by_center = {}
        used_idx = set()
        for a_i in range(len(pi_entries)):
            idx1, a1, b1 = pi_entries[a_i]
            if idx1 in used_idx:
                continue
            for b_i in range(a_i + 1, len(pi_entries)):
                idx2, a2, b2 = pi_entries[b_i]
                if idx2 in used_idx:
                    continue
                shared = {a1, b1} & {a2, b2}
                if len(shared) != 1:
                    continue
                c_idx = list(shared)[0]
                x1 = b1 if a1 == c_idx else a1
                x2 = b2 if a2 == c_idx else a2
                if x1 == x2:
                    continue
                c_p = [i for i, l in enumerate(ao_lbls) if int(l.split()[0]) == c_idx and "2p" in l]
                if len(c_p) != 3:
                    continue
                u1 = orbitals[idx1].mo_coeff[c_p]
                u2 = orbitals[idx2].mo_coeff[c_p]
                n1, n2 = float(np.linalg.norm(u1)), float(np.linalg.norm(u2))
                if n1 > 1e-6 and n2 > 1e-6 and abs(float(np.dot(u1 / n1, u2 / n2))) > 0.70:
                    used_idx.add(idx1)
                    used_idx.add(idx2)
                    ref_normal = _local_out_normal(c_idx, min(x1, x2), max(x1, x2))
                    out_align = abs(float(np.dot(u1 / n1, ref_normal)))
                    pair = (
                        ((idx1, x1), (idx2, x2), out_align)
                        if x1 < x2
                        else ((idx2, x2), (idx1, x1), out_align)
                    )
                    by_center.setdefault((c_idx, min(x1, x2), max(x1, x2)), []).append(pair)
                    break

        for (c_idx, xA, xB), pairs in by_center.items():
            assignments = []
            if len(pairs) == 1:
                (iA, _), (iB, _), _ = pairs[0]
                p_mat = self._iao_atom_populations(
                    np.column_stack([orbitals[iA].mo_coeff, orbitals[iB].mo_coeff]), Q_iao, ref_ao_atom
                )
                # Prefer xA (lower index) when symmetric within 1% so in-plane choice stays frame-invariant
                if p_mat[0, c_idx] >= p_mat[1, c_idx] - 0.01:
                    assignments.append((iA, xA, iB, xB))
                else:
                    assignments.append((iB, xB, iA, xA))
            elif len(pairs) == 2:
                pairs_sorted = sorted(pairs, key=lambda item: -item[2])
                (iA_out, _), (iB_out, _), _ = pairs_sorted[0]
                (iA_in, _), (iB_in, _), _ = pairs_sorted[1]
                p_out = self._iao_atom_populations(
                    np.column_stack([orbitals[iA_out].mo_coeff, orbitals[iB_out].mo_coeff]), Q_iao, ref_ao_atom
                )
                # Prefer xB (higher index) for out-of-plane and xA for in-plane when symmetric within 1%
                if p_out[0, c_idx] > p_out[1, c_idx] + 0.01:
                    assignments.append((iA_out, xA, iB_out, xB))
                    assignments.append((iB_in, xB, iA_in, xA))
                else:
                    assignments.append((iB_out, xB, iA_out, xA))
                    assignments.append((iA_in, xA, iB_in, xB))

            th = np.radians(10.0)
            for i_bond, x_bond, i_lp, x_lp in assignments:
                vb0 = orbitals[i_bond].mo_coeff.copy()
                vl0 = orbitals[i_lp].mo_coeff.copy()
                cands = np.column_stack([
                    np.cos(th) * vb0 + np.sin(th) * vl0,
                    np.cos(-th) * vb0 + np.sin(-th) * vl0,
                ])
                c_pops = self._iao_atom_populations(cands, Q_iao, ref_ao_atom)[:, c_idx]
                best_th = th if c_pops[0] >= c_pops[1] else -th
                vb = np.cos(best_th) * vb0 + np.sin(best_th) * vl0
                vlp = -np.sin(best_th) * vb0 + np.cos(best_th) * vl0
                orbitals[i_bond].mo_coeff = vb
                orbitals[i_bond].energy_ev = self._eval_energy_ev(vb)
                orbitals[i_bond].atoms = sorted([self.atom_ids[c_idx], self.atom_ids[x_bond]], key=_canonical_atom_sort_key)
                orbitals[i_lp].mo_coeff = vlp
                orbitals[i_lp].energy_ev = self._eval_energy_ev(vlp)
                orbitals[i_lp].type = "lone_pair"
                orbitals[i_lp].atoms = [self.atom_ids[x_lp]]

    def _align_coplanar_pistar_atoms(self, orbitals: List[Orbital]) -> None:
        """Match perpendicular π* antibonds on cumulated A=C=B centers to their coplanar occupied π bond atoms."""
        ao_lbls = self.mol.ao_labels()
        occ_pi = [o for o in orbitals if o.type == "pi" and len(o.atoms) == 2]
        virt_pi = [o for o in orbitals if o.type == "pi_star" and len(o.atoms) == 2]
        if len(occ_pi) < 2 or len(virt_pi) < 2:
            return
        for vp in virt_pi:
            for c_id in vp.atoms:
                try:
                    c_idx = self.atom_ids.index(c_id)
                except ValueError:
                    continue
                c_p = [i for i, l in enumerate(ao_lbls) if int(l.split()[0]) == c_idx and "2p" in l]
                if len(c_p) != 3:
                    continue
                uv = vp.mo_coeff[c_p]
                nv = float(np.linalg.norm(uv))
                if nv < 1e-6:
                    continue
                cand_occ = [o for o in occ_pi if c_id in o.atoms]
                if len(cand_occ) == 2 and cand_occ[0].atoms != cand_occ[1].atoms:
                    best_o, best_dot = None, -1.0
                    for o in cand_occ:
                        uo = o.mo_coeff[c_p]
                        no = float(np.linalg.norm(uo))
                        if no > 1e-6:
                            d = abs(float(np.dot(uv / nv, uo / no)))
                            if d > best_dot:
                                best_dot, best_o = d, o
                    if best_o is not None and best_dot > 0.70:
                        vp.atoms = list(best_o.atoms)
                        break
