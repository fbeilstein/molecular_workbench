"""
Shared physical constants, conversion factors, and default parameters for IRC reaction workflows.
"""

# Physical and unit conversions
BOHR_TO_ANG: float = 0.529177210903
ANG_TO_BOHR: float = 1.0 / BOHR_TO_ANG
EV_TO_KCAL: float = 23.060548
HARTREE_TO_EV: float = 27.211386245988
HARTREE_TO_KCAL: float = 627.509474

# Docking and reaction coordinate geometry adjustments (Angstrom)
DOCK_EXTRA: float = 1.0       # Pre-reaction complex: forming bond = d0 + 1.0 A
BREAK_EXTRA: float = 1.4      # End of drive: breaking bond = d0 + 1.4 A

# Pipeline & Server defaults
DEFAULT_PORT_IRC: int = 8082
DEFAULT_SOLVENT: str = "water"
DEFAULT_IRC_GRID_POINTS: int = 32
DEFAULT_IRC_MARGIN_ANG: float = 2.5
DEFAULT_IRC_MAX_STEPS: int = 35
DEFAULT_IRC_ISOVALUE: float = 0.045

# Centralized Directory Architecture
import os
ROOT_DIR: str = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
CACHE_DIR_NAME: str = "irc_cache"
CACHE_DIR: str = os.path.join(ROOT_DIR, CACHE_DIR_NAME)
OUT_DIR: str = CACHE_DIR  # Alias for backward compatibility
STEPS_DIR: str = os.path.join(CACHE_DIR, "steps")
CASCADE_DIR: str = os.path.join(CACHE_DIR, "cascade")
CUBES_DIR: str = os.path.join(CACHE_DIR, "cubes")
BUNDLES_DIR: str = os.path.join(CACHE_DIR, "bundles")
PROJECT_FILE: str = os.path.join(CACHE_DIR, "project.json")
LEGACY_STATE_FILE: str = os.path.join(CACHE_DIR, "reaction_state.json")
LEGACY_TRAJ_FILE: str = os.path.join(CACHE_DIR, "irc_trajectory.json")
MASTER_PATHWAY_FILE: str = os.path.join(CACHE_DIR, "master_pathway.json")
PROLONGED_PATH: str = os.path.join(CASCADE_DIR, "prolonged_pathway.json")
NPZ_PATH: str = os.path.join(CASCADE_DIR, "orbital_tracks.npz")
REGISTRY_PATH: str = os.path.join(CASCADE_DIR, "orbital_registry.json")

