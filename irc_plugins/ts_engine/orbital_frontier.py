import os
import sys
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from pyscf import gto

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

BOHR_TO_ANG = 0.52917721092
ANG_TO_BOHR = 1.8897259886

_COV_RADII = {
    "H": 0.32, "C": 0.77, "N": 0.75, "O": 0.73, "F": 0.71,
    "S": 0.99, "P": 1.06, "Cl": 0.99, "Br": 1.14, "I": 1.33,
}


def detect_planar_pi_fragments(
    mol: Any, max_plane_dist_bohr: float = 0.50
) -> List[Tuple[List[int], np.ndarray, List[int]]]:
    """
    Detect connected planar conjugated heavy-atom fragments directly from 3D coordinates.
    Matches tools/qm/analyzers/topology.py Step-0 planarity criterion (max_dist <= 0.50 Bohr),
    iteratively pruning non-cyclic out-of-plane substituents while preserving planar rings.
    Returns list of (sorted_atom_indices, unit_normal, sorted_cycle_atom_indices).
    """
    coords_b = mol.atom_coords()
    coords_a = coords_b * BOHR_TO_ANG
    n = mol.natm
    syms = [mol.atom_symbol(i) for i in range(n)]

    full_adj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            rc = (_COV_RADII.get(syms[i], 0.80) + _COV_RADII.get(syms[j], 0.80)) * 1.25
            if float(np.linalg.norm(coords_a[i] - coords_a[j])) <= rc:
                full_adj[i].append(j)
                full_adj[j].append(i)

    core_unsat = {
        i for i, s in enumerate(syms)
        if s in ("C", "N") and 1 <= len(full_adj[i]) <= 3
    }
    cand = set(core_unsat)
    for i, s in enumerate(syms):
        if s in ("O", "S") and len(full_adj[i]) in (1, 2):
            if any(
                j in core_unsat and float(np.linalg.norm(coords_a[i] - coords_a[j])) < 1.42
                for j in full_adj[i]
            ):
                cand.add(i)

    sub_adj = {i: [j for j in full_adj[i] if j in cand] for i in cand}

    def _cycle_nodes(nodes: Any) -> set:
        rem = set(nodes)
        changed = True
        while changed:
            changed = False
            for u in list(rem):
                if sum(1 for v in sub_adj[u] if v in rem) < 2:
                    rem.remove(u)
                    changed = True
        return rem

    def _clean_stubs(nodes: Any) -> List[int]:
        rem = set(nodes)
        changed = True
        while changed:
            changed = False
            cyc = _cycle_nodes(rem)
            for u in list(rem):
                nbs = [v for v in sub_adj[u] if v in rem]
                if len(nbs) == 0:
                    rem.remove(u)
                    changed = True
                elif len(nbs) == 1 and u not in cyc:
                    if syms[u] == "C" and float(np.linalg.norm(coords_a[u] - coords_a[nbs[0]])) > 1.39:
                        rem.remove(u)
                        changed = True
        return sorted(rem)

    visited: set = set()
    comps: List[List[int]] = []
    for u in sorted(cand):
        if u not in visited:
            comp, q = [], [u]
            visited.add(u)
            while q:
                cur = q.pop(0)
                comp.append(cur)
                for nb in sub_adj[cur]:
                    if nb not in visited:
                        visited.add(nb)
                        q.append(nb)
            comp = _clean_stubs(comp)
            if len(comp) >= 3:
                comps.append(comp)

    planar_frags: List[Tuple[List[int], np.ndarray, List[int]]] = []
    for comp in comps:
        cur = list(comp)
        while len(cur) >= 3:
            cb = coords_b[cur]
            centered = cb - cb.mean(axis=0)
            _, evecs = np.linalg.eigh(centered.T @ centered)
            normal = evecs[:, 0]
            if normal[2] < 0 or (abs(normal[2]) < 1e-6 and normal[1] < 0):
                normal = -normal
            dists = np.abs(centered @ normal)
            if float(np.max(dists)) <= max_plane_dist_bohr:
                planar_frags.append((sorted(cur), normal, sorted(_cycle_nodes(cur))))
                break
            cyc = _cycle_nodes(cur)
            exo = [k for k, u in enumerate(cur) if u not in cyc]
            if not exo:
                break
            worst_k = max(exo, key=lambda k: float(dists[k]))
            cur.pop(worst_k)
            cur = _clean_stubs(cur)

    return planar_frags


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


def extract_step0_planar_pi(
    mol: Any, mf: Any, atom_labels: Optional[List[str]] = None
) -> Tuple[np.ndarray, np.ndarray, Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]:
    """
    Step-0 Planar Conjugated π Subspace Separation (matching server.py / tools/qm/pipeline.py).
    Un-mixes pure out-of-plane p_perp conjugated π and π* subspaces from hyperconjugation tails via SVD,
    canonicalizes them with the Fock operator, and excludes them from the IBO localization subspaces.
    Returns (C_non_pi_occ, C_non_pi_virt, canonical_pi_meta, canonical_pi_vecs).
    """
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    mo_coeff = mf.mo_coeff[0] if mol.spin > 0 else mf.mo_coeff
    mo_energy = mf.mo_energy[0] if mol.spin > 0 else mf.mo_energy

    core_mask = (mo_energy[:nocc] * 27.2114) < -50.0
    C_non_pi_occ = mo_coeff[:, :nocc][:, ~core_mask].copy()
    C_non_pi_virt = mo_coeff[:, nocc:].copy()

    all_atom_ids = (
        list(atom_labels)
        if (atom_labels and len(atom_labels) == mol.natm)
        else [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
    )

    ovlp = mol.intor_symmetric("int1e_ovlp")
    fock = mf.get_fock()
    ao_labels = mol.ao_labels(fmt=False)
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]
    coords_b = mol.atom_coords()

    canonical_pi_meta: Dict[str, Dict[str, Any]] = {}
    canonical_pi_vecs: Dict[str, np.ndarray] = {}

    planar_frags = detect_planar_pi_fragments(mol, max_plane_dist_bohr=0.50)

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
        n_pi = int(np.sum((s_occ ** 2) >= 0.55))
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
                "energy_ev": round(float(e_po[i] * 27.2114), 2),
                "atoms": atom_ids,
            }

        # 2. Virtual conjugated pi* subspace separation
        if C_non_pi_virt.shape[1] > 0:
            M_virt = Q_pi.T @ ovlp @ C_non_pi_virt
            _, s_virt, Vt_virt = np.linalg.svd(M_virt, full_matrices=True)
            n_pv = min(n_pi, int(np.sum((s_virt ** 2) >= 0.35)))
            if n_pv > 0:
                C_virt_rot = C_non_pi_virt @ Vt_virt.T
                C_pi_virt = C_virt_rot[:, :n_pv]
                C_non_pi_virt = C_virt_rot[:, n_pv:]

                e_pv, U_pv = np.linalg.eigh(C_pi_virt.T @ fock @ C_pi_virt)
                C_pi_virt = C_pi_virt @ U_pv
                C_pi_virt = _snap_degenerate_subspace_pairs(mol, C_pi_virt, e_pv, fg[0], ovlp)

                for i in range(n_pv):
                    vec = C_pi_virt[:, i].copy()
                    norm = np.sqrt(float(vec.T @ ovlp @ vec))
                    if norm > 1e-6:
                        vec /= norm
                    ps = np.outer(vec, vec) @ ovlp
                    ap = np.zeros(mol.natm)
                    for mu in range(mol.nao):
                        ap[ao_atom[mu]] += ps[mu, mu]
                    dom_ai = max(active_fg, key=lambda ai: float(ap[ai]))
                    pt_up = (coords_b[dom_ai] + 0.5 * normal).reshape(1, 3)
                    if float((mol.eval_gto("GTOval", pt_up) @ vec)[0]) < 0.0:
                        vec = -vec

                    key = f"DELOC_PISTAR_{frag_lbl}_{i+1:02d}"
                    canonical_pi_vecs[key] = vec
                    canonical_pi_meta[key] = {
                        "name": key,
                        "label": f"π*_deloc({frag_short}) [π*{i+1}]",
                        "description": f"Delocalized conjugated π*{i+1}",
                        "category": "antibonding",
                        "type": "delocalized_pistar",
                        "chem_key": f"deloc_pistar_{frag_lbl}_{i+1}",
                        "mo_index": 100 + i + 1,
                        "energy_ev": round(float(e_pv[i] * 27.2114), 2),
                        "atoms": atom_ids,
                    }

    return C_non_pi_occ, C_non_pi_virt, canonical_pi_meta, canonical_pi_vecs


def extract_frontier_orbitals(
    mol: Any,
    mf: Any,
    prev_mol: Optional[Any] = None,
    prev_mo_coeffs: Optional[Dict[str, np.ndarray]] = None,
    ibos_m: Optional[Dict[str, Any]] = None,
    star_m: Optional[Dict[str, Any]] = None,
    ibos_v: Optional[Dict[str, np.ndarray]] = None,
    atom_labels: Optional[List[str]] = None,
    step0_pi: Optional[Tuple[Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]] = None,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]:
    """Extract canonical frontier MOs and Step-0 planar conjugated pi systems."""
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    mo_coeff = mf.mo_coeff[0] if mol.spin > 0 else mf.mo_coeff
    mo_energy = mf.mo_energy[0] if mol.spin > 0 else mf.mo_energy

    all_atom_ids = (
        list(atom_labels)
        if (atom_labels and len(atom_labels) == mol.natm)
        else [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
    )

    ovlp = mol.intor_symmetric("int1e_ovlp")
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]
    s_cross = gto.intor_cross("int1e_ovlp", prev_mol, mol) if (prev_mol and prev_mo_coeffs) else None

    target_keys = [
        ("HOMO-3", nocc - 4, "HOMO-3"),
        ("HOMO-2", nocc - 3, "HOMO-2"),
        ("HOMO-1", nocc - 2, "HOMO-1"),
        ("HOMO", nocc - 1, "HOMO (Frontier Donor)"),
        ("LUMO", nocc, "LUMO (Frontier Acceptor)"),
        ("LUMO+1", nocc + 1, "LUMO+1"),
        ("LUMO+2", nocc + 2, "LUMO+2"),
    ]

    orbitals_meta: Dict[str, Dict[str, Any]] = {}
    frontier_vectors: Dict[str, np.ndarray] = {}

    for name, mo_idx, label in target_keys:
        if 0 <= mo_idx < mo_coeff.shape[1]:
            vec = mo_coeff[:, mo_idx].copy()
            if s_cross is not None and prev_mo_coeffs and name in prev_mo_coeffs:
                if float(prev_mo_coeffs[name] @ s_cross @ vec) < 0.0:
                    vec = -vec

            ps = np.outer(vec, vec) @ ovlp
            ap = np.zeros(mol.natm)
            for mu in range(mol.nao):
                ap[ao_atom[mu]] += ps[mu, mu]
            dom_atoms = [all_atom_ids[i] for i in np.argsort(-ap) if ap[i] >= 0.12][:4]

            frontier_vectors[name] = vec
            e_ev = float(mo_energy[mo_idx] * 27.2114)
            orbitals_meta[name] = {
                "name": name,
                "label": label,
                "category": "frontier",
                "type": "frontier",
                "mo_index": int(mo_idx),
                "energy_ev": round(e_ev, 3),
                "is_occupied": bool(mo_idx < nocc),
                "atoms": dom_atoms,
            }

    if step0_pi is not None:
        canonical_pi_meta, canonical_pi_vecs = step0_pi
    else:
        _, _, canonical_pi_meta, canonical_pi_vecs = extract_step0_planar_pi(mol, mf, atom_labels=atom_labels)

    if s_cross is not None and prev_mo_coeffs:
        for k, v in list(canonical_pi_vecs.items()):
            if k in prev_mo_coeffs and float(prev_mo_coeffs[k] @ s_cross @ v) < 0.0:
                canonical_pi_vecs[k] = -v

    all_vecs = {**frontier_vectors, **canonical_pi_vecs}
    return orbitals_meta, canonical_pi_meta, all_vecs
