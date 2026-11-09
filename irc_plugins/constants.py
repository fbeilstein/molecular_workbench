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

