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
    """Compute unit normal vector of local coordination plane at atom_idx (returns None for sp linear centers)."""
    if len(neighbors) < 2:
        return None
    best_norm, best_sin = None, 0.0
    for i in range(len(neighbors)):
        for j in range(i + 1, len(neighbors)):
            v1 = coords_bohr[neighbors[i]] - coords_bohr[atom_idx]
            v2 = coords_bohr[neighbors[j]] - coords_bohr[atom_idx]
            n1 = float(np.linalg.norm(v1))
            n2 = float(np.linalg.norm(v2))
            if n1 < 1e-6 or n2 < 1e-6:
                continue
            norm = np.cross(v1, v2)
            l = float(np.linalg.norm(norm))
            sin_ang = l / (n1 * n2)
            if sin_ang > best_sin and l > 1e-6:
                best_sin = sin_ang
                best_norm = norm / l
    # Linear sp centers (∠ ~ 180°, sin(θ) < 0.15) do not define a unique local plane normal
    if best_sin < 0.15:
        return None
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
        evals, evecs = np.linalg.eigh(centered.T @ centered)
        # Collinear triads (e.g. O=C=O, C-C≡N, C-C≡C) span only 1D and do not define a 2D plane seed
        if float(evals[1]) < 0.15:
            return
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
    max_plane_dist_bohr: float = 0.45,
    min_cos_align: float = 0.92,
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

    # Candidate unsaturated planar backbone atoms:
    # sp2/sp C and N (excluding pyramidal sp3 N like amines and tetrahedral TS carbons),
    # plus 1- or 2-coordinate O/S within conjugated bond distance (< 1.42 Å; excludes 3-coord oxonium).
    max_hetero_bond = min(CONJUGATED_HETERO_BOND_MAX_ANG, 1.42)
    core_unsat: Set[int] = set()
    for i, s in enumerate(syms):
        if s in ("C", "N") and 1 <= len(full_adj[i]) <= 3:
            if len(full_adj[i]) == 3 and _angle_sum_deg(i, full_adj[i]) < 345.0:
                continue
            core_unsat.add(i)

    cand = set(core_unsat)
    for i, s in enumerate(syms):
        if s in ("O", "S") and len(full_adj[i]) in (1, 2):
            if any(
                j in core_unsat and float(np.linalg.norm(coords_a[i] - coords_a[j])) < max_hetero_bond
                for j in full_adj[i]
            ):
                cand.add(i)

    adj: Dict[int, Set[int]] = {u: set(v for v in full_adj[u] if v in cand) for u in cand}
    if not adj:
        return []

    # Precompute local normal vector using full_adj for 2+/3-coord centers so sp2 carbonyls
    # use their true trigonal plane normal (even if one substituent like protonated O+ is non-conjugated).
    local_normals: Dict[int, np.ndarray] = {}
    for i in adj:
        nbs_for_norm = full_adj[i] if len(full_adj[i]) >= 2 else list(adj[i])
        n_vec = get_local_normal(coords_b, i, nbs_for_norm)
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
        sorted_s = sorted(list(s))
        final_frags.append((sorted_s, norm, sorted_s))

    return final_frags

