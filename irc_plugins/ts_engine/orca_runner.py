"""
ORCA Runner Facade
Maintains backward compatibility by re-exporting the modular ORCA engine.
Implementation has been decomposed into irc_plugins.ts_engine.orca subpackage.
"""

from .orca import (
    ORCA_BIN,
    run_orca_calc,
    parse_orca_trj_file,
    build_trajectory_result,
    run_orca_dft_ts,
    eval_orca_dft_profile,
    run_orca_native_irc,
    run_orca_neb_ci,
    run_orca_thermo_profile,
)

# Private aliases for backward compatibility
_run_orca_calc = run_orca_calc
_parse_orca_trj_file = parse_orca_trj_file
_build_trajectory_result = build_trajectory_result

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
    "_run_orca_calc",
    "_parse_orca_trj_file",
    "_build_trajectory_result",
]
