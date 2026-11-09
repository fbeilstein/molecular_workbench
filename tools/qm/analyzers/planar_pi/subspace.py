"""
Planar Pi Subspace Separation Module — Step-0 Occupied & Virtual Conjugated π Extraction.
"""

from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np
from pyscf import lo

from qm.classification import _atom_populations, _classify_orbital, count_core_orbitals
from qm.constants import (
    HARTREE_TO_EV,
    PLANAR_PI_MAX_PLANE_DIST_BOHR,
    PLANAR_PI_MIN_SVD_WEIGHT,
)
from .geometry import detect_planar_pi_fragments


def _snap_degenerate_subspace_pairs(
    mol: Any, C_sub: np.ndarray, energies_ha: np.ndarray, ref_atom: int, ovlp: np.ndarray
) -> np.ndarray:
    """Break arbitrary 2D degeneracy (e.g. symmetric benzene e1g/e2u) using target atom population."""
    if C_sub.shape[1] < 2:
        return C_sub
    ao_labels = mol.ao_labels(fmt=False)
    P_A = np.zeros_like(ovlp)
    for mu, (atom_idx, *_) in enumerate(ao_labels):
        if atom_idx == ref_atom:
            P_A[mu, :] += 0.5 * ovlp[mu, :]
            P_A[:, mu] += 0.5 * ovlp[:, mu]

    C_out = C_sub.copy()
    i = 0
    while i < C_out.shape[1] - 1:
        if abs(float(energies_ha[i]) - float(energies_ha[i + 1])) < 1e-3:
            pair = C_out[:, [i, i + 1]]
            _, evecs = np.linalg.eigh(pair.T @ P_A @ pair)
            C_out[:, [i, i + 1]] = pair @ evecs
            i += 2
        else:
            i += 1
    return C_out


def _extract_localized_pistar_projections(
    mol: Any,
    C_non_pi_occ: np.ndarray,
    C_virt: np.ndarray,
    ovlp: np.ndarray,
    ao_labels: List[Any],
    all_atom_ids: List[str],
) -> List[np.ndarray]:
    """Find external localized 2-center pi bonds from non-pi occupied space and project antibonds."""
    loc_pistar_vecs: List[np.ndarray] = []
    if C_non_pi_occ.shape[1] == 0:
        return loc_pistar_vecs

    try:
        iaos = lo.iao.iao(mol, C_non_pi_occ)
        loc_occ = lo.ibo.ibo(mol, C_non_pi_occ, iaos=iaos)
        atom_labels_list = [mol.atom_symbol(i) for i in range(mol.natm)]
        for i in range(loc_occ.shape[1]):
            mo = loc_occ[:, i]
            pop = _atom_populations(mol, mo, ovlp)
            info = _classify_orbital(mol, mo, pop, atom_labels_list, all_atom_ids, ovlp)
            if info.get("type") == "pi" and len(info.get("atoms", [])) == 2:
                atoms = info["atoms"]
                idx1 = all_atom_ids.index(atoms[0])
                idx2 = all_atom_ids.index(atoms[1])
                c_A = np.zeros_like(mo)
                c_B = np.zeros_like(mo)
                for mu, (atom_idx, *_) in enumerate(ao_labels):
                    if atom_idx == idx1:
                        c_A[mu] = mo[mu]
                    elif atom_idx == idx2:
                        c_B[mu] = mo[mu]
                c_anti = c_A - c_B
                c_v = C_virt @ (C_virt.T @ ovlp @ c_anti)
                nrm = np.sqrt(float(c_v.T @ ovlp @ c_v))
                if nrm > 1e-4:
                    c_v /= nrm
                    loc_pistar_vecs.append(c_v)
    except Exception:
        pass
    return loc_pistar_vecs


def _extract_iao_valence_p_vecs(
    mol: Any,
    fg: List[int],
    normal: np.ndarray,
    ovlp: np.ndarray,
    ao_labels: List[Any],
    iao_vecs: Optional[np.ndarray],
    iao_ao_labels: Optional[List[Any]],
) -> List[np.ndarray]:
    """Build minimal valence p_perp vectors across fragment atoms via IAO or valence AOs."""
    p_vecs: List[np.ndarray] = []
    for ai in fg:
        z = mol.atom_charge(ai)
        val_n = "2p" if z <= 10 else ("3p" if z <= 18 else "4p")
        v = None
        if iao_vecs is not None and iao_ao_labels is not None:
            found_iao = {}
            for mu_iao, (a, _el, nl, m) in enumerate(iao_ao_labels):
                if a == ai and val_n in nl and m in ("x", "y", "z"):
                    found_iao[m] = mu_iao
            if len(found_iao) == 3:
                v_cand = (
                    iao_vecs[:, found_iao["x"]] * normal[0]
                    + iao_vecs[:, found_iao["y"]] * normal[1]
                    + iao_vecs[:, found_iao["z"]] * normal[2]
                )
                nrm = np.sqrt(float(v_cand.T @ ovlp @ v_cand))
                if nrm > 1e-8:
                    v = v_cand / nrm

        if v is None:
            v_ao = np.zeros(mol.nao)
            found = {}
            for mu, (a, _el, nl, m) in enumerate(ao_labels):
                if a == ai and val_n in nl and m in ("x", "y", "z"):
                    found[m] = mu
            if len(found) == 3:
                v_ao[found["x"]] = normal[0]
                v_ao[found["y"]] = normal[1]
                v_ao[found["z"]] = normal[2]
                nrm = np.sqrt(float(v_ao.T @ ovlp @ v_ao))
                if nrm > 1e-8:
                    v = v_ao / nrm

        if v is not None:
            p_vecs.append(v)
    return p_vecs


def extract_step0_planar_pi(
    mol: Any,
    mf: Any,
    atom_labels: Optional[List[str]] = None,
    smiles: Optional[str] = None,
    max_plane_dist_bohr: float = PLANAR_PI_MAX_PLANE_DIST_BOHR,
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]:
    """
    Step-0 Planar Conjugated π Subspace Separation.
    Un-mixes pure out-of-plane p_perp conjugated π and π* subspaces from hyperconjugation tails via SVD,
    canonicalizes them with the Fock operator, and excludes them from the IBO localization subspaces.
    Returns (C_non_pi_occ, C_non_pi_virt, canonical_pi_meta, canonical_pi_vecs).
    """
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    mo_coeff = mf.mo_coeff[0] if mol.spin > 0 else mf.mo_coeff
    mo_energy = mf.mo_energy[0] if mol.spin > 0 else mf.mo_energy

    n_core = count_core_orbitals(mol, mo_energy, nocc)
    C_non_pi_occ = mo_coeff[:, n_core:nocc].copy()
    C_virt = mo_coeff[:, nocc:].copy()
    C_non_pi_virt = C_virt.copy()

    all_atom_ids = (
        list(atom_labels)
        if (atom_labels and len(atom_labels) == mol.natm)
        else [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
    )

    ovlp = mol.intor_symmetric("int1e_ovlp")
    raw_fock = mf.get_fock()
    fock = (
        raw_fock[0]
        if (mol.spin > 0 and isinstance(raw_fock, (tuple, list, np.ndarray)) and np.ndim(raw_fock) == 3)
        else raw_fock
    )
    ao_labels = mol.ao_labels(fmt=False)
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]
    coords_b = mol.atom_coords()

    canonical_pi_meta: Dict[str, Dict[str, Any]] = {}
    canonical_pi_vecs: Dict[str, np.ndarray] = {}

    planar_frags = detect_planar_pi_fragments(mol, smiles=smiles, max_plane_dist_bohr=max_plane_dist_bohr)
    frag_occ_data = []

    for fg, normal, cyc_nodes in planar_frags:
        if C_non_pi_occ.shape[1] == 0:
            break
        fg_set = set(fg)
        p_shells: Dict[Tuple[int, str], Dict[str, int]] = {}
        for mu, (ai, _el, nl, m) in enumerate(ao_labels):
            if ai in fg_set and "p" in nl and m in ("x", "y", "z"):
                p_shells.setdefault((ai, nl), {})[m] = mu

        p_vecs = []
        for (_ai, _nl), m_map in sorted(p_shells.items()):
            if len(m_map) == 3:
                v = np.zeros(mol.nao)
                v[m_map["x"]] = normal[0]
                v[m_map["y"]] = normal[1]
                v[m_map["z"]] = normal[2]
                nrm = np.sqrt(float(v.T @ ovlp @ v))
                if nrm > 1e-8:
                    p_vecs.append(v / nrm)
        if not p_vecs:
            continue

        B_p = np.column_stack(p_vecs)
        ev_b, U_b = np.linalg.eigh(B_p.T @ ovlp @ B_p)
        Q_pi = B_p @ (U_b @ np.diag(1.0 / np.sqrt(np.maximum(ev_b, 1e-8))) @ U_b.T)

        # 1. Occupied conjugated pi subspace separation
        M_occ = Q_pi.T @ ovlp @ C_non_pi_occ
        _, s_occ, Vt_occ = np.linalg.svd(M_occ, full_matrices=True)
        n_pi = int(np.sum((s_occ ** 2) >= PLANAR_PI_MIN_SVD_WEIGHT))
        if n_pi == 0:
            continue

        C_occ_rot = C_non_pi_occ @ Vt_occ.T
        C_pi_occ = C_occ_rot[:, :n_pi]
        C_non_pi_occ = C_occ_rot[:, n_pi:]

        e_po, U_po = np.linalg.eigh(C_pi_occ.T @ fock @ C_pi_occ)
        C_pi_occ = C_pi_occ @ U_po
        C_pi_occ = _snap_degenerate_subspace_pairs(mol, C_pi_occ, e_po, fg[0], ovlp)

        # Active conjugated atoms check
        pi_atom_pop = np.zeros(mol.natm)
        occ_col_pops = []
        for k in range(n_pi):
            vp = C_pi_occ[:, k]
            ps = np.outer(vp, vp) @ ovlp
            ap = np.zeros(mol.natm)
            for mu in range(mol.nao):
                ap[ao_atom[mu]] += ps[mu, mu] * 2.0
            occ_col_pops.append(ap)
            pi_atom_pop += ap

        active_fg = [i for i in fg if pi_atom_pop[i] >= 0.10]
        if len(active_fg) < 3:
            active_fg = list(fg)
        atom_ids = [all_atom_ids[i] for i in active_fg]
        frag_lbl = "-".join(atom_ids)

        cyc_active = set(cyc_nodes) & set(active_fg)
        exo_ids = [all_atom_ids[i] for i in active_fg if i not in cyc_active]
        if len(cyc_active) == 6:
            frag_short = f"Ar-{'-'.join(exo_ids)}" if exo_ids else "Ar"
        elif len(atom_ids) <= 4:
            frag_short = frag_lbl
        else:
            frag_short = f"{atom_ids[0]}..{atom_ids[-1]}"

        for i in range(n_pi):
            vec = C_pi_occ[:, i].copy()
            norm = np.sqrt(float(vec.T @ ovlp @ vec))
            if norm > 1e-6:
                vec /= norm
            dom_ai = max(active_fg, key=lambda ai: float(occ_col_pops[i][ai]))
            pt_up = (coords_b[dom_ai] + 0.5 * normal).reshape(1, 3)
            if float((mol.eval_gto("GTOval", pt_up) @ vec)[0]) < 0.0:
                vec = -vec

            key = f"DELOC_PI_{frag_lbl}_{i+1:02d}"
            node_str = f"{i} node{'s' if i != 1 else ''}" if i > 0 else "all-bonding, 0 nodes"
            canonical_pi_vecs[key] = vec
            canonical_pi_meta[key] = {
                "name": key,
                "label": f"π_deloc({frag_short}) [π{i+1}]",
                "description": f"Delocalized conjugated π{i+1} ({node_str})",
                "category": "conjugated_pi",
                "type": "delocalized_pi",
                "chem_key": f"deloc_pi_{frag_lbl}_{i+1}",
                "mo_index": i + 1,
                "energy_ev": round(float(e_po[i] * HARTREE_TO_EV), 2),
                "atoms": atom_ids,
            }

        frag_occ_data.append((fg, normal, cyc_nodes, active_fg, atom_ids, frag_lbl, frag_short, n_pi))

    # 2. Localized pi* protection
    loc_pistar_vecs = _extract_localized_pistar_projections(mol, C_non_pi_occ, C_virt, ovlp, ao_labels, all_atom_ids)
    if loc_pistar_vecs:
        V_loc = np.column_stack(loc_pistar_vecs)
        P_perp_loc = np.eye(mol.nao) - V_loc @ np.linalg.inv(V_loc.T @ ovlp @ V_loc) @ V_loc.T @ ovlp
        C_virt_deloc_raw = P_perp_loc @ C_virt
        metric = C_virt_deloc_raw.T @ ovlp @ C_virt_deloc_raw
        ev, U_m = np.linalg.eigh(metric)
        keep = ev > 1e-6
        C_virt_deloc = C_virt_deloc_raw @ (U_m[:, keep] / np.sqrt(ev[keep]))
    else:
        C_virt_deloc = C_virt.copy()

    # 3. Delocalized pi* subspace extraction for each planar fragment
    iao_vecs, iao_ao_labels = None, None
    if nocc > 0:
        try:
            occ_for_iao = mo_coeff[:, :nocc]
            iao_raw = lo.iao.iao(mol, occ_for_iao)
            iao_vecs = lo.vec_lowdin(iao_raw, ovlp)
            ref_mol = lo.iao.reference_mol(mol)
            iao_ao_labels = ref_mol.ao_labels(fmt=False)
        except Exception:
            iao_vecs, iao_ao_labels = None, None

    all_deloc_pi_virt = []
    for fg, normal, cyc_nodes, active_fg, atom_ids, frag_lbl, frag_short, n_pi in frag_occ_data:
        n_pistar = len(fg) - n_pi
        if n_pistar <= 0 or C_virt_deloc.shape[1] == 0:
            continue

        p_vecs = _extract_iao_valence_p_vecs(mol, fg, normal, ovlp, ao_labels, iao_vecs, iao_ao_labels)
        if not p_vecs:
            continue

        B_p = np.column_stack(p_vecs)
        ev_b, U_b = np.linalg.eigh(B_p.T @ ovlp @ B_p)
        Q_pi_val = B_p @ (U_b @ np.diag(1.0 / np.sqrt(np.maximum(ev_b, 1e-8))) @ U_b.T)

        M_virt = Q_pi_val.T @ ovlp @ C_virt_deloc
        _, _, Vt_virt = np.linalg.svd(M_virt, full_matrices=False)
        C_pi_virt = C_virt_deloc @ Vt_virt[:n_pistar].T
        e_pv, U_pv = np.linalg.eigh(C_pi_virt.T @ fock @ C_pi_virt)
        C_pi_virt = C_pi_virt @ U_pv
        C_pi_virt = _snap_degenerate_subspace_pairs(mol, C_pi_virt, e_pv, fg[0], ovlp)

        for i in range(C_pi_virt.shape[1]):
            vec = C_pi_virt[:, i].copy()
            norm = np.sqrt(float(vec.T @ ovlp @ vec))
            if norm > 1e-6:
                vec /= norm
            dom_ai = active_fg[0]
            pt_up = (coords_b[dom_ai] + 0.5 * normal).reshape(1, 3)
            if float((mol.eval_gto("GTOval", pt_up) @ vec)[0]) < 0.0:
                vec = -vec

            key = f"DELOC_PISTAR_{frag_lbl}_{i+1:02d}"
            canonical_pi_vecs[key] = vec
            canonical_pi_meta[key] = {
                "name": key,
                "label": f"π*_deloc({frag_short}) [π*{i+1}]",
                "description": f"Delocalized antibonding π*{i+1}",
                "category": "conjugated_pi_star",
                "type": "delocalized_pistar",
                "chem_key": f"deloc_pistar_{frag_lbl}_{i+1}",
                "mo_index": i + 1,
                "energy_ev": round(float(e_pv[i] * HARTREE_TO_EV), 2),
                "atoms": atom_ids,
            }
            all_deloc_pi_virt.append(vec)

    # 4. Project ONLY delocalized pi* out of C_virt to get C_non_pi_virt for ValenceLocalizer
    if all_deloc_pi_virt:
        V_deloc = np.column_stack(all_deloc_pi_virt)
        P_perp_deloc = np.eye(mol.nao) - V_deloc @ np.linalg.inv(V_deloc.T @ ovlp @ V_deloc) @ V_deloc.T @ ovlp
        C_non_pi_virt_raw = P_perp_deloc @ C_virt
        metric2 = C_non_pi_virt_raw.T @ ovlp @ C_non_pi_virt_raw
        ev2, U_m2 = np.linalg.eigh(metric2)
        keep2 = ev2 > 1e-6
        C_non_pi_virt = C_non_pi_virt_raw @ (U_m2[:, keep2] / np.sqrt(ev2[keep2]))
    else:
        C_non_pi_virt = C_virt.copy()

    return C_non_pi_occ, C_non_pi_virt, canonical_pi_meta, canonical_pi_vecs

