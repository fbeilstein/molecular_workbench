import numpy as np
from ase import Atoms
from ase.optimize import LBFGS
from sella import Sella
from .calculators import make_xtb, _make_restrained_calc

def run_guided_manual_ts(atoms, model, log):
    """Two-Stage Transition-Zone Valley Relaxation for crude manual positioning."""
    log("=== Strategy: Quantum Valley Relaxation (Crude Manual Placement) ===")
    log("  Stage 1: Relaxing into reaction valley under quantum forces (angles free)...")

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)
    targets = []
    for i, j in model.formed:
        targets.append((i, j, model.d0(i, j) + 0.65))
    for i, j in model.broken:
        targets.append((i, j, model.d0(i, j) + 0.35))

    pos_init = atoms.positions.copy()
    spec = []
    for i, j in model.spectator_bonds:
        d_eq = float(np.linalg.norm(pos_init[i] - pos_init[j]))
        spec.append((i, j, d_eq))

    calc = _make_restrained_calc(base, targets, k=35.0, spectator_restraints=spec, model=model)
    atoms.calc = calc
    opt_stage1 = LBFGS(atoms, logfile=None, maxstep=0.15)
    opt_stage1.run(fmax=0.05, steps=70)
    log(f"  Stage 1 finished in {opt_stage1.nsteps} steps. Fragments aligned into reaction valley.")

    log("  Stage 2: Polishing saddle point with Sella DLC (guarded)...")
    from .calculators import guard_saddle_polish
    from .vibrations import analyse_mode
    atoms.calc = make_xtb(model.charge, model.multiplicity)
    ts_rep = None
    try:
        ts_rep = analyse_mode(atoms, model, log)
    except Exception:
        pass

    if ts_rep and ts_rep.get("n_imag") == 1 and ts_rep.get("mode_consistent", True):
        log(f"  ✓ Saddle point confirmed by vibrational analysis ({ts_rep['imag_freqs'][0]:.0f}i cm-1).")
        return [{"step": 1, "energy_ev": float(atoms.get_potential_energy()), "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))}]

    if ts_rep and ts_rep.get("n_imag") == 0:
        log("  ℹ️ Local minimum detected in scan: steering along reaction coordinate into saddle region...")

    return guard_saddle_polish(atoms, model, log, max_steps=20, max_drift=0.15)

