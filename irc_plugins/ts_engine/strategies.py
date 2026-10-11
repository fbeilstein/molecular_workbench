import numpy as np
from .common import EV_TO_KCAL, DOCK_EXTRA, BREAK_EXTRA
from .calculators import make_xtb, _make_restrained_calc, build_reaction_eigenvector


def _polish_with_sella(atoms, model, log, max_steps=120):
    """Polish a candidate saddle point using Sella in internal coordinates."""
    from sella import Sella
    atoms.calc = make_xtb(model.charge, model.multiplicity)
    try:
        opt = Sella(atoms, internal=True, order=1, logfile=None)
        opt.run(fmax=0.03, steps=max_steps)
    except Exception as e:
        log(f"  Sella polish note: {e}")


# ---------------------------------------------------------------------------
# Strategy 1: Constrained Coordinate Drive across barrier peak -> Sella
# ---------------------------------------------------------------------------
def _build_drive_targets(model, pos, pairs):
    d_cur = np.array([np.linalg.norm(pos[i] - pos[j]) for i, j in pairs])
    d_start, d_end = [], []
    for (i, j), dc in zip(pairs, d_cur):
        d0 = model.d0(i, j)
        if (i, j) in model.formed:
            d_start.append(min(max(dc, d0 + DOCK_EXTRA), d0 + 2.5))
            d_end.append(d0)
        else:
            d_start.append(min(dc, d0 + 0.1))
            d_end.append(d0 + BREAK_EXTRA)
    return np.array(d_start), np.array(d_end)


def run_drive_sella_ts(atoms, model, log, n_points=11, k=25.0):
    """Strategy: Constrained Coordinate Drive across barrier peak -> Sella."""
    from ase.optimize import LBFGS
    log("\n=== Strategy: Constrained Coordinate Drive -> Sella ===")

    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    d_start, d_end = _build_drive_targets(model, pos, pairs)
    spec_restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)

    frames, energies, lams = [], [], []
    i_max = 0

    for n, lam in enumerate(np.linspace(0.0, 1.0, n_points)):
        targets = d_start + lam * (d_end - d_start)
        calc = _make_restrained_calc(base, [(i, j, t) for (i, j), t in zip(pairs, targets)], k,
                                     spectator_restraints=spec_restraints)
        atoms.calc = calc
        opt = LBFGS(atoms, logfile=None, maxstep=0.15)
        opt.run(fmax=0.08, steps=100)

        e = calc.e_base
        frames.append(atoms.get_positions().copy())
        energies.append(e)
        lams.append(float(lam))
        if e > energies[i_max]:
            i_max = n

        rel = (e - energies[0]) * EV_TO_KCAL
        log(f"  lambda={lam:4.2f}  dE={rel:+7.2f} kcal/mol")

        passed_barrier = (i_max >= 2 and (energies[i_max] - energies[0]) * EV_TO_KCAL > 3.0
                          and (energies[i_max] - e) * EV_TO_KCAL > 6.0)
        if passed_barrier:
            log("  (passed barrier maximum - stopping scan early)")
            break

    atoms.set_positions(frames[i_max])
    log(f"  Barrier peak selected at lambda={lams[i_max]:.2f}. Refining with Sella...")
    _polish_with_sella(atoms, model, log, max_steps=120)

    profile = [[l, (e - energies[0]) * EV_TO_KCAL] for l, e in zip(lams, energies)]
    return profile


# ---------------------------------------------------------------------------
# Strategy 2: Climbing Image NEB (CI-NEB with IDPP)
# ---------------------------------------------------------------------------
def _build_relaxed_product_endpoint(atoms, model):
    """Build product endpoint P via restrained GFN2-xTB relaxation."""
    from ase.optimize import LBFGS
    p_atoms = atoms.copy()
    p_targets = [(i, j, model.d0(i, j)) for i, j in model.formed]
    p_targets += [(i, j, model.d0(i, j) + 1.25) for i, j in model.broken]
    spec_restraints = [(i, j, float(np.linalg.norm(atoms.positions[i] - atoms.positions[j])))
                       for i, j in model.spectator_bonds]

    p_base = p_atoms.copy()
    p_base.calc = make_xtb(model.charge, model.multiplicity)
    p_atoms.calc = _make_restrained_calc(p_base, p_targets, k=30.0, spectator_restraints=spec_restraints)
    p_opt = LBFGS(p_atoms, logfile=None, maxstep=0.2)
    p_opt.run(fmax=0.08, steps=50)
    return p_atoms.get_positions()


def run_neb_ts(atoms, model, log, n_images=7):
    """Strategy: Climbing Image Nudged Elastic Band (CI-NEB with IDPP)."""
    from ase.optimize import MDMin
    from ase.mep import NEB
    from ase.mep.neb import idpp_interpolate
    log("\n=== Strategy: Climbing Image NEB (CI-NEB) ===")

    pos_p = _build_relaxed_product_endpoint(atoms, model)

    images = [atoms.copy() for _ in range(n_images)]
    images[0].set_positions(atoms.get_positions())
    images[-1].set_positions(pos_p)
    idpp_interpolate(images, mic=False)

    for img in images:
        img.calc = make_xtb(model.charge, model.multiplicity)

    neb = NEB(images, climb=True, method='improvedtangent')
    opt = MDMin(neb, logfile=None)
    opt.run(fmax=0.10, steps=80)

    energies = [img.get_potential_energy() for img in images]
    mid_idx = int(np.argmax(energies[1:-1])) + 1
    atoms.set_positions(images[mid_idx].get_positions())
    log(f"  CI-NEB finished. Image {mid_idx} selected as barrier top (E = {energies[mid_idx]:.4f} eV). Refining with Sella...")

    _polish_with_sella(atoms, model, log, max_steps=100)
    profile = [[float(idx), float((e - energies[0]) * EV_TO_KCAL)] for idx, e in enumerate(energies)]
    return profile


# ---------------------------------------------------------------------------
# Strategy 3: Dimer Method (ase.mep.dimer)
# ---------------------------------------------------------------------------
def _precondition_dimer_guess(atoms, model, log):
    from ase.optimize import LBFGS
    pos = atoms.get_positions()
    far_forming = [pair for pair in model.formed
                   if np.linalg.norm(pos[pair[0]] - pos[pair[1]]) > model.d0(pair[0], pair[1]) + 1.0]
    if not far_forming:
        return

    log("  Initial manual placement has large forming distance; pre-contracting into the TS zone...")
    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)
    targets = [(i, j, model.d0(i, j) + 0.65) for i, j in model.formed]
    spec = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]
    atoms.calc = _make_restrained_calc(base, targets, 35.0, spectator_restraints=spec)
    pre_opt = LBFGS(atoms, logfile=None, maxstep=0.2)
    pre_opt.run(fmax=0.1, steps=25)
    atoms.set_positions(atoms.get_positions())
    log("  Pre-contraction complete. Initializing Dimer saddle search...")


def run_dimer_ts(atoms, model, log, steps=100, fmax=0.04):
    """Strategy: Dimer Method (ase.mep.dimer)."""
    from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate
    log("\n=== Strategy: Dimer Method (ase.mep.dimer) ===")
    log("Dynamically searching and following lowest-curvature reaction mode using force replicas...")

    _precondition_dimer_guess(atoms, model, log)

    atoms.calc = make_xtb(model.charge, model.multiplicity)
    d_control = DimerControl(
        initial_eigenmode_method='displacement',
        displacement_method='vector',
        dimer_separation=0.005,
        maximum_translation=0.08,
        trial_trans_step=0.005,
        f_rot_min=0.05,
        f_rot_max=0.5,
        max_num_rot=3
    )
    d_atoms = MinModeAtoms(atoms, d_control)

    v0 = build_reaction_eigenvector(atoms, model)
    if v0 is not None:
        try:
            d_atoms.initialize_eigenmodes(eigenmodes=[v0])
            log("  Seeded dimer eigenmode along forming/breaking bonds.")
        except Exception:
            pass

    opt = MinModeTranslate(d_atoms, logfile=None)

    def _log_dimer():
        try:
            f = float(np.max(np.linalg.norm(d_atoms.get_forces(), axis=1)))
            e = atoms.get_potential_energy()
            if opt.nsteps % 5 == 0 or f < fmax:
                log(f"  Dimer step {opt.nsteps:3d}: E = {e:.4f} eV, max|F| = {f:.4f} eV/A")
        except Exception:
            pass

    opt.attach(_log_dimer)
    try:
        opt.run(fmax=fmax, steps=steps)
        log("  Dimer stage finished.")
    except Exception as e:
        log(f"  Dimer search note: {e}")

    log("  Polishing Dimer saddle point with Sella...")
    _polish_with_sella(atoms, model, log, max_steps=30)
    return True


# ---------------------------------------------------------------------------
# Strategy 4: Pre-Relaxation -> Sella
# ---------------------------------------------------------------------------
def run_seeded_sella_ts(atoms, model, log):
    """Strategy: Constrained Pre-Relaxation -> Sella."""
    from ase.optimize import LBFGS
    log("\n=== Strategy: Pre-Relaxation -> Sella ===")
    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in pairs]
    spec_restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity)
    atoms.calc = _make_restrained_calc(base, restraints, 25.0, spectator_restraints=spec_restraints)
    opt = LBFGS(atoms, logfile=None, maxstep=0.15)
    opt.run(fmax=0.06, steps=120)
    log(f"  Pre-relaxation complete. Relaxed energy = {atoms.calc.e_base:.4f} eV")

    _polish_with_sella(atoms, model, log, max_steps=150)
    return []
