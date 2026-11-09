"""
Common definitions, paths, and constants for Orbital Tracker.
"""

import os

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, "trajectories")
NPZ_PATH = os.path.join(TRAJECTORIES_DIR, "orbital_tracks.npz")
REGISTRY_PATH = os.path.join(TRAJECTORIES_DIR, "orbital_registry.json")

FRONTIER_KEYS = ["HOMO-3", "HOMO-2", "HOMO-1", "HOMO", "LUMO", "LUMO+1", "LUMO+2"]
FRONTIER_TRACK_MAP = {
    "HOMO-3": "trk_homo_3",
    "HOMO-2": "trk_homo_2",
    "HOMO-1": "trk_homo_1",
    "HOMO": "trk_homo",
    "LUMO": "trk_lumo",
    "LUMO+1": "trk_lumo_1",
    "LUMO+2": "trk_lumo_2",
}

SUBSCRIPTS = ["₁", "₂", "₃", "₄", "₅"]

