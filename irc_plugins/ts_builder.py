"""
Facade re-exporting transition-state engine components for backward compatibility.
Implementation logic has been modularized into `irc_plugins.ts_engine`:
  - `ts_engine.common`: Periodic table, vector helpers, Kabsch alignment, molblock conversions.
  - `ts_engine.reaction_model`: Chemical graph tracking, bond order classifications, attack vectors.
  - `ts_engine.docking`: ETKDG conformer generation, multi-start stochastic grid docking.
  - `ts_engine.refinement`: In-place fragment geometry refinement (MMFF94 / GFN2-xTB).
  - `ts_engine.calculators`: GFN2-xTB factory, harmonic restraint calculator with spectator locks.
  - `ts_engine.strategies`: Multi-strategy TS optimizers (Drive->Sella, CI-NEB, Dimer, Seeded Sella).
  - `ts_engine.vibrations`: Normal-mode vibrational frequency projection & analysis.
"""

from .ts_engine import (
    ReactionModel,
    embed_fragments,
    dock_fragments,
    auto_dock_molblocks,
    refine_fragment_shape,
    molblocks_to_positions,
    positions_to_molblocks,
    set_positions,
    align_to_reference,
    make_xtb,
    _make_restrained_calc,
    check_geometry_sanity,
    build_reaction_eigenvector,
    run_rits_ts,
    run_lego_ts,
    run_guided_manual_ts,
    run_genetic_swarm_ts,
    run_orca_dft_ts,
    run_dimer_ts,
    run_drive_sella_ts,
    run_seeded_sella_ts,
    run_neb_ts,
    analyse_mode,
    rcov,
    rvdw,
    EV_TO_KCAL,
    DOCK_EXTRA,
    BREAK_EXTRA,
)

