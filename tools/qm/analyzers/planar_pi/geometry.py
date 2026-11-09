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
    core_unsat = {
        i for i, s in enumerate(syms)
        if s in ("C", "N") and 1 <= len(full_adj[i]) <= 3
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
    if not adj:
        return []

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
    final_frags: List[Tuple[List[int], np.ndarray, List[int]]] = []
    for s, norm in candidates:
        if any(len(s & set(f)) >= 0.70 * len(s) for f, _, _ in final_frags):
            continue
        final_frags.append((sorted(list(s)), norm, sorted(list(s))))

    return final_frags

