"""
Alignment Engine Package
Modular chemical graph isomorphism mapping, rigid superposition, and geodesic bridge interpolation.
"""

from .rigid_superposition import kabsch_fit, compute_rmsd, align_coordinates
from .graph_mapper import find_step_chemical_overlaps, extract_fragment_molblock
from .bridge_interpolator import generate_bridge_frames, build_multistep_bridge
from .window_optimizer import find_optimal_overlap_window

__all__ = [
    "kabsch_fit",
    "compute_rmsd",
    "align_coordinates",
    "find_step_chemical_overlaps",
    "extract_fragment_molblock",
    "generate_bridge_frames",
    "build_multistep_bridge",
    "find_optimal_overlap_window",
]
