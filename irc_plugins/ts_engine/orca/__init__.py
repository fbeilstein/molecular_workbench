"""
ORCA Engine Package
Modular execution, trajectory parsing, and calculation workflows for ORCA DFT.
"""

from .executor import ORCA_BIN, run_orca_calc
from .parser import parse_orca_trj_file, build_trajectory_result
from .calculations import (
    run_orca_dft_ts,
    eval_orca_dft_profile,
    run_orca_native_irc,
    run_orca_neb_ci,
    run_orca_thermo_profile,
)

__all__ = [
    "ORCA_BIN",
    "run_orca_calc",
    "parse_orca_trj_file",
    "build_trajectory_result",
    "run_orca_dft_ts",
    "eval_orca_dft_profile",
    "run_orca_native_irc",
    "run_orca_neb_ci",
    "run_orca_thermo_profile",
]
