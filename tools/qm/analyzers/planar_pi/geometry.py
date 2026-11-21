"""
Planar Pi Geometry Module — Planar BFS Flood-Fill & Geometric Coplanarity.
"""

from typing import Any, Dict, List, Optional, Set, Tuple
import numpy as np

from qm.constants import (
    ANG_TO_BOHR,
    BOHR_TO_ANG,
    COVALENT_RADII_ANG,
    CONJUGATED_HETERO_BOND_MAX_ANG,
    PLANAR_PI_MAX_PLANE_DIST_BOHR,
    PLANAR_PI_MIN_COS_ALIGN,
)

_COV_RADII = COVALENT_RADII_ANG


def get_local_normal(coords_bohr: np.ndarray, atom_idx: int, neighbors: List[int]) -> Optional[np.ndarray]:
    """Compute unit normal vector of local coordination plane at atom_idx."""
    if len(neighbors) < 2:
        return None
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


def _find_planar_seeds(
    adj: Dict[int, Set[int]],
    coords_b: np.ndarray,
    local_normals: Dict[int, np.ndarray],
    min_cos_align: float,
) -> List[Set[int]]:
    """Generate 3-atom planar seed triads whose normals align with local coordination planes."""
    seeds: List[Set[int]] = []

    def _check_and_add_triad(triad: List[int]) -> None:
        pts = coords_b[triad]
        centered = pts - pts.mean(axis=0)
        _, evecs = np.linalg.eigh(centered.T @ centered)
        fnorm = evecs[:, 0]
        for a in triad:
            if a in local_normals and abs(float(np.dot(local_normals[a], fnorm))) < min_cos_align:
                return
        seeds.append(set(triad))

    for u in adj:
        nbrs = sorted(list(adj[u]))
        for i in range(len(nbrs)):
            for j in range(i + 1, len(nbrs)):
                _check_and_add_triad([nbrs[i], u, nbrs[j]])

    for u in adj:
        for v in adj[u]:
            for w in adj.get(v, []):
                if w != u and w not in adj[u]:
                    _check_and_add_triad([u, v, w])

    return seeds


def _expand_seed_bfs(
    seed: Set[int],
    adj: Dict[int, Set[int]],
    coords_b: np.ndarray,
    local_normals: Dict[int, np.ndarray],
    max_plane_dist_bohr: float,
    min_cos_align: float,
) -> Optional[Tuple[Set[int], np.ndarray]]:
    """Flood-fill candidate atoms via BFS while maintaining plane distance and normal alignment."""
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

    if len(curr) < 3:
        return None

    pts = coords_b[list(curr)]
    centered = pts - pts.mean(axis=0)
    _, evecs = np.linalg.eigh(centered.T @ centered)
    normal = evecs[:, 0]
    if normal[2] < 0 or (abs(normal[2]) < 1e-6 and normal[1] < 0):
        normal = -normal
    return curr, normal


def detect_planar_pi_fragments(
    mol: Any,
    smiles: Optional[str] = None,
    max_plane_dist_bohr: float = PLANAR_PI_MAX_PLANE_DIST_BOHR,
    min_cos_align: float = PLANAR_PI_MIN_COS_ALIGN,
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
    # Exclude carboxylic acid / ester / amide carbons (carbons with >= 2 heteroatom neighbors like -C(=O)OR, -C(=O)Cl)
    # to keep their carbonyl double bonds cleanly localized and prevent false flood-fill into adjacent rings.
    hetero_neighbors = {
        i: sum(1 for nbr in full_adj[i] if syms[nbr] in ("O", "S", "N", "Cl", "F", "Br"))
        for i in range(n)
    }
    core_unsat = {
        i for i, s in enumerate(syms)
        if s in ("C", "N") and 1 <= len(full_adj[i]) <= 3 and hetero_neighbors[i] < 2
    }
    cand = set(core_unsat)
    for i, s in enumerate(syms):
        if s in ("O", "S") and len(full_adj[i]) in (1, 2):
            if any(
                j in core_unsat and float(np.linalg.norm(coords_a[i] - coords_a[j])) < CONJUGATED_HETERO_BOND_MAX_ANG
                for j in full_adj[i]
            ):
                cand.add(i)

    adj: Dict[int, Set[int]] = {u: set(v for v in full_adj[u] if v in cand) for u in cand}
    final_frags: List[Tuple[List[int], np.ndarray, List[int]]] = []

    if adj:
        # Precompute local normal vector for each atom
        local_normals: Dict[int, np.ndarray] = {}
        for i in adj:
            n_vec = get_local_normal(coords_b, i, list(adj[i]))
            if n_vec is not None:
                local_normals[i] = n_vec

        # Form 3-atom seeds & BFS flood-fill
        seeds = _find_planar_seeds(adj, coords_b, local_normals, min_cos_align)
        candidates: List[Tuple[Set[int], np.ndarray]] = []
        for seed in seeds:
            res = _expand_seed_bfs(seed, adj, coords_b, local_normals, max_plane_dist_bohr, min_cos_align)
            if res is not None:
                candidates.append(res)

        # Deduplicate & select maximal non-overlapping planar fragments
        candidates.sort(key=lambda item: len(item[0]), reverse=True)
        for s, norm in candidates:
            if len(s) <= 3 and sum(1 for ai in s if syms[ai] in ("O", "S", "N")) >= 2:
                continue
            if any(len(s & set(f)) >= 0.70 * len(s) for f, _, _ in final_frags):
                continue
            final_frags.append((sorted(list(s)), norm, sorted(list(s))))

    # 2. Dedicated planar 3-center (or 4-center) O=C-X (X = O, N) carbonyl conjugation detector
    #    Captures esters (-C(=O)-OR), carboxylic acids (-C(=O)-OH), carboxylates (-COO-),
    #    amides (-C(=O)-NR2), ureas, and carbonates as standalone conjugated pi systems without
    #    allowing them to falsely flood-fill into sterically twisted adjacent aromatic rings.
    used_atoms: Set[int] = set()
    for f_atoms, _, _ in final_frags:
        used_atoms.update(f_atoms)

    def _angle_sum_deg(center_idx: int, nbrs_list: List[int]) -> float:
        total = 0.0
        m = len(nbrs_list)
        for ii in range(m):
            for jj in range(ii + 1, m):
                v1 = coords_a[nbrs_list[ii]] - coords_a[center_idx]
                v2 = coords_a[nbrs_list[jj]] - coords_a[center_idx]
                n1, n2 = float(np.linalg.norm(v1)), float(np.linalg.norm(v2))
                if n1 > 1e-6 and n2 > 1e-6:
                    total += float(np.degrees(np.arccos(np.clip(np.dot(v1, v2) / (n1 * n2), -1.0, 1.0))))
        return total

    for c_idx in range(n):
        if syms[c_idx] != "C" or c_idx in used_atoms:
            continue
        c_nbrs = full_adj[c_idx]
        if len(c_nbrs) != 3:
            continue
        # Must be planar trigonal sp2 carbon (not a pyramidal tetrahedral intermediate)
        if _angle_sum_deg(c_idx, c_nbrs) < 350.0:
            continue
        c_norm = get_local_normal(coords_b, c_idx, c_nbrs)
        if c_norm is None:
            continue

        term_o = [
            nb for nb in c_nbrs
            if syms[nb] == "O"
            and len(full_adj[nb]) == 1
            and float(np.linalg.norm(coords_a[c_idx] - coords_a[nb])) < 1.32
        ]
        if not term_o:
            continue

        conj_x = []
        for nb in c_nbrs:
            if nb in term_o:
                continue
            d_cx = float(np.linalg.norm(coords_a[c_idx] - coords_a[nb]))
            if d_cx >= 1.42:
                continue
            if syms[nb] == "O":
                # 1-coord carboxylate O- or 2-coord neutral ester/acid -OR/-OH (excludes 3-coord protonated oxonium!)
                if len(full_adj[nb]) not in (1, 2):
                    continue
            elif syms[nb] == "N":
                if len(full_adj[nb]) not in (2, 3):
                    continue
                if len(full_adj[nb]) == 3 and _angle_sum_deg(nb, full_adj[nb]) < 345.0:
                    continue
            else:
                continue

            # Check coplanarity of X's local plane with the carbonyl plane
            if len(full_adj[nb]) >= 2:
                x_norm = get_local_normal(coords_b, nb, full_adj[nb])
                if x_norm is None or abs(float(np.dot(c_norm, x_norm))) < min_cos_align:
                    continue
            conj_x.append(nb)

        # Require at least 2 conjugated O/N atoms on C (e.g. 1 carbonyl =O + >=1 ester/acid/amide O/N, or 2 carboxylate O's)
        conj_atoms = sorted(set([c_idx] + term_o + conj_x))
        if len(conj_atoms) >= 3 and not (set(conj_atoms) & used_atoms):
            norm_out = c_norm.copy()
            if norm_out[2] < 0 or (abs(norm_out[2]) < 1e-6 and norm_out[1] < 0):
                norm_out = -norm_out
            # Note: pass empty cycle list [] so frag_short uses atom labels (e.g. C6-O8-O24)
            final_frags.append((conj_atoms, norm_out, []))
            used_atoms.update(conj_atoms)

    return final_frags

