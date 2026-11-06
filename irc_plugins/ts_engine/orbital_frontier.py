"""
Orbital Frontier & Planar Pi Subspace Extraction (Re-export Wrapper)

Unified with tools/qm/analyzers/planar_pi.py so that both server.py and server_irc.py
share the exact same planar BFS flood-fill and Step-0 subspace separation engine.
"""

import os
import sys

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

from qm.analyzers.planar_pi import (
    detect_planar_pi_fragments,
    extract_step0_planar_pi,
    extract_frontier_orbitals,
    _snap_degenerate_subspace_pairs,
    BOHR_TO_ANG,
    ANG_TO_BOHR,
)

__all__ = [
    "detect_planar_pi_fragments",
    "extract_step0_planar_pi",
    "extract_frontier_orbitals",
    "_snap_degenerate_subspace_pairs",
    "BOHR_TO_ANG",
    "ANG_TO_BOHR",
]
