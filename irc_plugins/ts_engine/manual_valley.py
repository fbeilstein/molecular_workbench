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

    calc = _make_restrained_calc(base, targets, k=35.0, spectator_restraints=spec)
    atoms.calc = calc
    opt_stage1 = LBFGS(atoms, logfile=None, maxstep=0.15)
    opt_stage1.run(fmax=0.05, steps=70)
    log(f"  Stage 1 finished in {opt_stage1.nsteps} steps. Fragments aligned into reaction valley.")

    log("  Stage 2: Releasing restraints and polishing saddle point with Sella DLC...")
    atoms.calc = make_xtb(model.charge, model.multiplicity)
    opt_sella = Sella(atoms, internal=True, order=1, logfile=None)
    profile = []

    def on_step():
        try:
            profile.append({
                "step": opt_sella.nsteps,
                "energy_ev": float(atoms.get_potential_energy()),
                "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
            })
        except Exception:
            pass

    opt_sella.attach(on_step)
    converged = opt_sella.run(fmax=0.03, steps=60)
    log(f"  Sella Polish {'Converged' if converged else 'Finished'} in {opt_sella.nsteps} steps.")
    return profile
