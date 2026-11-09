"""
Planar Pi Engine — Unified Step-0 Planar Conjugated π Subspace Separation.

Shared across single-molecule workbench (server.py) and reaction IRC trajectories (server_irc.py).
Uses Planar BFS Flood-Fill to detect maximal contiguous planar conjugated fragments
without arbitrary ring-versus-chain assumptions, separates out-of-plane p_perp subspaces via SVD,
canonicalizes with the Fock operator, breaks 2D degeneracy, and leaves non-pi subspaces for IBO.
"""

from qm.constants import BOHR_TO_ANG, ANG_TO_BOHR
from .geometry import get_local_normal, detect_planar_pi_fragments
from .subspace import extract_step0_planar_pi, _snap_degenerate_subspace_pairs
from .frontier import extract_frontier_orbitals

__all__ = [
    "get_local_normal",
    "detect_planar_pi_fragments",
    "_snap_degenerate_subspace_pairs",
    "extract_step0_planar_pi",
    "extract_frontier_orbitals",
    "BOHR_TO_ANG",
    "ANG_TO_BOHR",
]
