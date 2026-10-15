import os
import sys
from typing import Any, List, Tuple
from pyscf import gto, scf

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

# Re-export grid, isosurface, and cube generation routines
from .orbital_grid import (
    ANG_TO_BOHR,
    BOHR_TO_ANG,
    build_global_grid,
    evaluate_mo_grid,
    extract_isosurface_mesh,
    generate_cube_text
)

# Re-export IBO and density flux calculation routines
from .orbital_ibo import (
    compute_density_flux,
    compute_frame_ibos,
    classify_ibo_quick as _classify_ibo_quick
)

# Re-export frontier and delocalized pi systems extraction
from .orbital_frontier import extract_frontier_orbitals


def compute_frame_scf(
    symbols: List[str],
    positions: List[List[float]],
    charge: int = -1,
    spin: int = 0,
    basis: str = "sto-3g"
) -> Tuple[Any, Any, List[float]]:
    """Run RHF/UHF on a frame geometry and return (mol, mf, mulliken_charges)."""
    atom_str = "; ".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(symbols, positions))
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol) if spin > 0 else scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    try:
        mull = mf.mulliken_pop(verbose=0)[1].round(3).tolist()
    except Exception:
        mull = [0.0] * len(symbols)
    return mol, mf, mull


__all__ = [
    "BOHR_TO_ANG",
    "ANG_TO_BOHR",
    "build_global_grid",
    "compute_frame_scf",
    "compute_density_flux",
    "compute_frame_ibos",
    "extract_frontier_orbitals",
    "evaluate_mo_grid",
    "extract_isosurface_mesh",
    "generate_cube_text",
    "_classify_ibo_quick"
]
