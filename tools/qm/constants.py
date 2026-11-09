"""
Centralized physical constants, conversion factors, chemical properties,
and default numerical thresholds for quantum-mechanical orbital pipeline.
"""

from typing import Dict, Tuple

# ==============================================================================
# Physical Units & Conversion Factors
# ==============================================================================

BOHR_TO_ANG = 0.52917721092
ANG_TO_BOHR = 1.8897259886
ANG2BOHR = ANG_TO_BOHR

HARTREE_TO_EV = 27.2114
HARTREE_TO_KCAL = 627.5095
EV_TO_KCAL = 23.060548

# ==============================================================================
# Chemical Elements & Radii (Angstroms)
# ==============================================================================

# Standard single-bond covalent radii (Cordero et al. / PySCF compatible)
COVALENT_RADII_ANG: Dict[str, float] = {
    "H": 0.31,
    "He": 0.28,
    "Li": 1.28,
    "Be": 0.96,
    "B": 0.84,
    "C": 0.76,
    "N": 0.71,
    "O": 0.66,
    "F": 0.57,
    "Ne": 0.58,
    "Na": 1.66,
    "Mg": 1.41,
    "Al": 1.21,
    "Si": 1.11,
    "P": 1.07,
    "S": 1.05,
    "Cl": 1.02,
    "Ar": 1.06,
    "K": 2.03,
    "Ca": 1.76,
    "Sc": 1.70,
    "Ti": 1.60,
    "V": 1.53,
    "Cr": 1.39,
    "Mn": 1.39,
    "Fe": 1.32,
    "Co": 1.26,
    "Ni": 1.24,
    "Cu": 1.32,
    "Zn": 1.22,
    "Ga": 1.22,
    "Ge": 1.20,
    "As": 1.19,
    "Se": 1.20,
    "Br": 1.20,
    "Kr": 1.16,
    "Rb": 2.20,
    "Sr": 1.95,
    "I": 1.39,
    "Xe": 1.40,
}

# Maximum normal valences for bond connectivity heuristics
MAX_VALENCE_DICT: Dict[str, int] = {
    "H": 1,
    "F": 1,
    "Cl": 1,
    "Br": 1,
    "I": 1,
    "O": 3,
    "N": 4,
    "C": 4,
    "S": 6,
    "P": 5,
    "Si": 4,
}

# Number of core electrons per element (for closed-shell core orbital identification)
CORE_ELECTRONS: Dict[str, int] = {
    "H": 0, "He": 0,
    "Li": 2, "Be": 2, "B": 2, "C": 2, "N": 2, "O": 2, "F": 2, "Ne": 2,
    "Na": 10, "Mg": 10, "Al": 10, "Si": 10, "P": 10, "S": 10, "Cl": 10, "Ar": 10,
    "K": 18, "Ca": 18,
    "Sc": 18, "Ti": 18, "V": 18, "Cr": 18, "Mn": 18, "Fe": 18, "Co": 18, "Ni": 18, "Cu": 18, "Zn": 18,
    "Ga": 28, "Ge": 28, "As": 28, "Se": 28, "Br": 28, "Kr": 28,
    "Rb": 36, "Sr": 36, "I": 46, "Xe": 46,
}

# ==============================================================================
# Quantum Orbital Classification Thresholds
# ==============================================================================

# Step 0: Planar Conjugated Pi Subspace
PLANAR_PI_MAX_PLANE_DIST_BOHR: float = 0.50     # ~0.26 Angstrom maximum out-of-plane deviation
PLANAR_PI_MIN_COS_ALIGN: float = 0.88          # cos(28 deg) minimum normal alignment
PLANAR_PI_MIN_SVD_WEIGHT: float = 0.55         # Minimum singular value weight for occupied pi inclusion
CONJUGATED_HETERO_BOND_MAX_ANG: float = 1.45   # Maximum bond length for conjugated C/N-heteroatom

# Step 1: Population Analysis
LP_DOMINANT_POPULATION_THRESHOLD: float = 0.88  # Isolated lone pair >= 88% on one heavy atom
LP_HYPERCONJUGATION_PRIMARY_MIN: float = 0.78   # Polarized lone pair >= 78% on primary atom...
LP_HYPERCONJUGATION_SECONDARY_MAX: float = 0.14 # ...and < 14% on any secondary atom
TWO_CENTER_BOND_MIN_POPULATION: float = 0.10    # Minimum contribution per atom for 2-center bond

# Step 2: Nodal Plane Discrete Contour Sampling
NODAL_SCAN_RADII_BOHR: Tuple[float, float] = (0.45, 0.85)  # Concentric circles (~0.24 A and ~0.45 A)
NODAL_SCAN_SAMPLE_POINTS: int = 12                         # Number of angular points per circle
NODAL_SIGNIFICANT_AMP_RATIO: float = 0.10                  # Points below 10% of max amplitude dropped
NODAL_MIN_SIGNIFICANT_POINTS: int = 4                      # Minimum valid points required for pi test
NODAL_AMPLITUDE_FLOOR: float = 1e-5                        # Disregard near-zero wavefunctions
SP3_MAX_NEIGHBOR_DIST_ANG: float = 1.95                    # Distance threshold for coordination count

# Post-Localization Snap ("Snap to Grid")
FRONTIER_SNAP_THRESHOLD: float = 0.90                      # Overlap threshold (w >= 0.90) for canonical snap
FRONTIER_DEFAULT_DEPTH: int = 0                            # Default depth around HOMO/LUMO (HOMO..LUMO)

# Default SCF and Geometry Optimization Parameters
DEFAULT_BASIS: str = "6-31g*"
DEFAULT_XC_METHOD: str = "b3lyp"
DEFAULT_GRID_POINTS: int = 50                              # 50x50x50 volumetric grid
DEFAULT_ORBITAL_ISOVALUE: float = 0.045                    # Standard isovalue for Marching Cubes
DEFAULT_PORT_BENCHMARK: int = 8080
DEFAULT_PORT_IRC: int = 8082

