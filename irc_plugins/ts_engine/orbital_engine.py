import os
import sys
from typing import Any, List, Optional, Tuple
from pyscf import gto, scf
from ..orbital_selection import SOLVENT_DIELECTRICS

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
from .orbital_frontier import (
    detect_planar_pi_fragments,
    extract_step0_planar_pi,
    extract_frontier_orbitals,
)


def compute_frame_scf(
    symbols: List[str],
    positions: List[List[float]],
    charge: int = -1,
    spin: int = 0,
    basis: str = "sto-3g",
    dm0: Any = None,
    solvent: Optional[str] = None,
) -> Tuple[Any, Any, List[float]]:
    """Run RHF/UHF on a frame geometry (with optional ddCOSMO dielectric solvation) and return (mol, mf, mulliken_charges)."""
    atom_str = "; ".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(symbols, positions))
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol) if spin > 0 else scf.RHF(mol)
    solv_key = str(solvent).strip().lower() if solvent else ""
    eps = SOLVENT_DIELECTRICS.get(solv_key) if solv_key not in ("", "gas", "none", "vacuum", "gas_phase") else None
    if eps is not None:
        mf = mf.ddCOSMO()
        mf.with_solvent.eps = float(eps)
        mf.with_solvent.grids.level = 0
        mf.with_solvent.lebedev_order = 11
        mf.with_solvent.lmax = 4
        mf.conv_tol = 1e-6
    mf.verbose = 0
    if dm0 is not None:
        mf.kernel(dm0=dm0)
    else:
        mf.kernel()
    if not mf.converged:
        mf = mf.newton()
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
    "detect_planar_pi_fragments",
    "extract_step0_planar_pi",
    "extract_frontier_orbitals",
    "evaluate_mo_grid",
    "extract_isosurface_mesh",
    "generate_cube_text",
    "_classify_ibo_quick"
]
