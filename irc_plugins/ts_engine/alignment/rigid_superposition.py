"""
Rigid Superposition Module
Implements Kabsch SE(3) superposition minimizing RMSD without reflection inversion.
"""

from typing import Tuple
import numpy as np


def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute optimal rotation matrix R in SO(3) and translation vector t
    minimizing RMSD(P, Q @ R.T + t).
    Guarantees det(R) = +1 (strictly forbids enantiomeric reflection).
    """
    cP = np.mean(P, axis=0)
    cQ = np.mean(Q, axis=0)
    H = (Q - cQ).T @ (P - cP)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    t = cP - cQ @ R.T
    return R, t


def compute_rmsd(P: np.ndarray, Q: np.ndarray) -> float:
    """Compute Root Mean Square Deviation between two sets of 3D coordinates."""
    if len(P) == 0 or len(Q) == 0 or len(P) != len(Q):
        return 0.0
    return float(np.sqrt(np.mean(np.sum((P - Q) ** 2, axis=1))))


def align_coordinates(coords: np.ndarray, R: np.ndarray, t: np.ndarray) -> np.ndarray:
    """Apply SE(3) transformation: coords_aligned = coords @ R.T + t."""
    if len(coords) == 0:
        return coords.copy()
    return coords @ R.T + t
