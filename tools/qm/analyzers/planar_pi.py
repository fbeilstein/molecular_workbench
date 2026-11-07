"""
Planar Pi Engine — Unified Step-0 Planar Conjugated π Subspace Separation

Shared across single-molecule workbench (server.py) and reaction IRC trajectories (server_irc.py).
Uses Planar BFS Flood-Fill to detect maximal contiguous planar conjugated fragments
without arbitrary ring-versus-chain assumptions, separates out-of-plane p_perp subspaces via SVD,
canonicalizes with the Fock operator, breaks 2D degeneracy, and leaves non-pi subspaces for IBO.
"""

import os
import sys
from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np
from pyscf import gto, lo

try:
    from qm.classification import _atom_populations, _classify_orbital
except ImportError:
    _tools_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
    if _tools_dir not in sys.path:
        sys.path.insert(0, _tools_dir)
    from qm.classification import _atom_populations, _classify_orbital

BOHR_TO_ANG = 0.52917721092
ANG_TO_BOHR = 1.8897259886

_COV_RADII = {
    "H": 0.32, "C": 0.77, "N": 0.75, "O": 0.73, "F": 0.71,
    "S": 0.99, "P": 1.06, "Cl": 0.99, "Br": 1.14, "I": 1.33,
}


def get_local_normal(coords_bohr: np.ndarray, atom_idx: int, neighbors: List[int]) -> Optional[np.ndarray]:
    """Compute unit normal vector of local coordination plane at atom_idx."""
    if len(neighbors) >= 2:
        best_norm, best_l = None, 0.0
        for i in range(len(neighbors)):
            for j in range(i + 1, len(neighbors)):
                v1 = coords_bohr[neighbors[i]] - coords_bohr[atom_idx]
                v2 = coords_bohr[neighbors[j]] - coords_bohr[atom_idx]
                norm = np.cross(v1, v2)
                l = float(np.linalg.norm(norm))
                if l > best_l:
                    best_l = l
                    best_norm = norm / l
        return best_norm
    return None


def detect_planar_pi_fragments(
    mol: Any,
    smiles: Optional[str] = None,
    max_plane_dist_bohr: float = 0.50,
    min_cos_align: float = 0.88,
) -> List[Tuple[List[int], np.ndarray, List[int]]]:
    """
    Planar BFS Flood-Fill: Detects maximal planar conjugated fragments directly from 3D geometry.
    Expands contiguous conjugated atoms via BFS as long as global fragment planarity 
    (max_dist <= max_plane_dist_bohr) and local normal alignment (cos >= min_cos_align) hold.
    Returns list of (sorted_atom_indices, unit_normal, sorted_cycle_atom_indices).
    """
    coords_b = mol.atom_coords()
    coords_a = coords_b * BOHR_TO_ANG
    n = mol.natm
    syms = [mol.atom_symbol(i) for i in range(n)]

    # 1. 3D covalent radii connectivity directly from atomic coordinates
    full_adj: Dict[int, List[int]] = {i: [] for i in range(n)}
    for i in range(n):
        for j in range(i + 1, n):
            rc = (_COV_RADII.get(syms[i], 0.80) + _COV_RADII.get(syms[j], 0.80)) * 1.25
            if float(np.linalg.norm(coords_a[i] - coords_a[j])) <= rc:
                full_adj[i].append(j)
                full_adj[j].append(i)

    # Candidate unsaturated backbone atoms (C, N with coordination <= 3, plus conjugated O, S)
    core_unsat = {
        i for i, s in enumerate(syms)
        if s in ("C", "N") and 1 <= len(full_adj[i]) <= 3
    }
    cand = set(core_unsat)
    for i, s in enumerate(syms):
        if s in ("O", "S") and len(full_adj[i]) in (1, 2):
            if any(
                j in core_unsat and float(np.linalg.norm(coords_a[i] - coords_a[j])) < 1.45
                for j in full_adj[i]
            ):
                cand.add(i)

    adj: Dict[int, Set[int]] = {u: set(v for v in full_adj[u] if v in cand) for u in cand}

    if not adj:
        return []

    # Precompute local normal vector for each atom
    local_normals: Dict[int, np.ndarray] = {}
    for i in adj:
        n_vec = get_local_normal(coords_b, i, list(adj[i]))
        if n_vec is not None:
            local_normals[i] = n_vec

    # Form 3-atom seeds
    seeds: List[Set[int]] = []
    for u in adj:
        nbrs = sorted(list(adj[u]))
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                v, w = nbrs[i], nbrs[j]
                triad = [v, u, w]
                pts = coords_b[triad]
                centered = pts - pts.mean(axis=0)
                _, evecs = np.linalg.eigh(centered.T @ centered)
                fnorm = evecs[:, 0]
                aligned = True
                for a in triad:
                    if a in local_normals:
                        if abs(float(np.dot(local_normals[a], fnorm))) < min_cos_align:
                            aligned = False
                            break
                if aligned:
                    seeds.append(set(triad))

    for u in adj:
        for v in adj[u]:
            for w in adj.get(v, []):
                if w != u and w not in adj[u]:
                    triad = [u, v, w]
                    pts = coords_b[triad]
                    centered = pts - pts.mean(axis=0)
                    _, evecs = np.linalg.eigh(centered.T @ centered)
                    fnorm = evecs[:, 0]
                    aligned = True
                    for a in triad:
                        if a in local_normals:
                            if abs(float(np.dot(local_normals[a], fnorm))) < min_cos_align:
                                aligned = False
                                break
                    if aligned:
                        seeds.append(set(triad))

    # Planar BFS Flood-Fill from each seed
    candidates: List[Tuple[Set[int], np.ndarray]] = []
    for seed in seeds:
        curr = set(seed)
        frontier: Set[int] = set()
        for a in curr:
            for nbr in adj.get(a, []):
                if nbr not in curr:
                    frontier.add(nbr)

        while frontier:
            best_nbr = None
            best_md = 999.0
            for nbr in list(frontier):
                cand_set = curr | {nbr}
                pts = coords_b[list(cand_set)]
                centered = pts - pts.mean(axis=0)
                _, evecs = np.linalg.eigh(centered.T @ centered)
                fnorm = evecs[:, 0]
                dists = np.abs(centered @ fnorm)
                md = float(np.max(dists))
                if md <= max_plane_dist_bohr:
                    if nbr in local_normals and abs(float(np.dot(local_normals[nbr], fnorm))) < min_cos_align:
                        continue
                    if md < best_md:
                        best_md = md
                        best_nbr = nbr

            if best_nbr is not None:
                curr.add(best_nbr)
                frontier.remove(best_nbr)
                for nbr in adj.get(best_nbr, []):
                    if nbr not in curr:
                        frontier.add(nbr)
            else:
                break

        if len(curr) >= 3:
            pts = coords_b[list(curr)]
            centered = pts - pts.mean(axis=0)
            _, evecs = np.linalg.eigh(centered.T @ centered)
            normal = evecs[:, 0]
            if normal[2] < 0 or (abs(normal[2]) < 1e-6 and normal[1] < 0):
                normal = -normal
            candidates.append((curr, normal))

    # Deduplicate & select maximal non-overlapping planar fragments
    candidates.sort(key=lambda item: len(item[0]), reverse=True)
    final_frags: List[Tuple[List[int], np.ndarray, List[int]]] = []
    for s, norm in candidates:
        if any(len(s & set(f)) >= 0.70 * len(s) for f, _, _ in final_frags):
            continue
        final_frags.append((sorted(list(s)), norm, sorted(list(s))))

    return final_frags


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
    mol: Any,
    mf: Any,
    atom_labels: Optional[List[str]] = None,
    smiles: Optional[str] = None,
    max_plane_dist_bohr: float = 0.50,
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

    core_mask = (mo_energy[:nocc] * 27.2114) < -50.0
    C_non_pi_occ = mo_coeff[:, :nocc][:, ~core_mask].copy()
    C_virt = mo_coeff[:, nocc:].copy()
    C_non_pi_virt = C_virt.copy()

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

        frag_occ_data.append((fg, normal, cyc_nodes, active_fg, atom_ids, frag_lbl, frag_short, n_pi))

    # 2. Localized pi* protection:
    # Identify any external 2-center pi bonds from the remaining non-delocalized occupied space C_non_pi_occ.
    # Construct their projected antibonds onto C_virt, and project them out to create a pure delocalized virtual search space C_virt_deloc.
    loc_pistar_vecs = []
    if C_non_pi_occ.shape[1] > 0:
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
    all_deloc_pi_virt = []
    for fg, normal, cyc_nodes, active_fg, atom_ids, frag_lbl, frag_short, n_pi in frag_occ_data:
        n_pistar = len(fg) - n_pi
        if n_pistar <= 0 or C_virt_deloc.shape[1] == 0:
            continue

        # Build compact valence p-shell projector for the fragment
        p_vecs = []
        for ai in fg:
            z = mol.atom_charge(ai)
            val_n = "2p" if z <= 10 else ("3p" if z <= 18 else "4p")
            v = np.zeros(mol.nao)
            found = {}
            for mu, (a, el, nl, m) in enumerate(ao_labels):
                if a == ai and val_n in nl and m in ("x", "y", "z"):
                    found[m] = mu
            if len(found) == 3:
                v[found["x"]] = normal[0]
                v[found["y"]] = normal[1]
                v[found["z"]] = normal[2]
                nrm = np.sqrt(float(v.T @ ovlp @ v))
                if nrm > 1e-8:
                    p_vecs.append(v / nrm)

        if p_vecs:
            B_p = np.column_stack(p_vecs)
            ev_b, U_b = np.linalg.eigh(B_p.T @ ovlp @ B_p)
            Q_pi_val = B_p @ (U_b @ np.diag(1.0 / np.sqrt(np.maximum(ev_b, 1e-8))) @ U_b.T)
        else:
            continue

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
                "energy_ev": round(float(e_pv[i] * 27.2114), 2),
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

