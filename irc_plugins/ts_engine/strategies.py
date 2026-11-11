import numpy as np
from .common import EV_TO_KCAL, DOCK_EXTRA, BREAK_EXTRA
from .calculators import make_xtb, _make_restrained_calc, build_reaction_eigenvector


def _polish_with_sella(atoms, model, log, max_steps=60):
    """Polish a candidate saddle point using Sella with TS trust zone guarding."""
    from .calculators import guard_saddle_polish
    from .vibrations import analyse_mode
    solvent = getattr(model, "solvent", "water")
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
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

    return guard_saddle_polish(atoms, model, log, max_steps=max_steps, max_drift=0.35)




# ---------------------------------------------------------------------------
# Strategy 1: Constrained Coordinate Drive across barrier peak -> Sella
# ---------------------------------------------------------------------------
def _build_drive_targets(model, pos, pairs, scan_start=None, scan_end=None):
    d_cur = np.array([np.linalg.norm(pos[i] - pos[j]) for i, j in pairs])
    d_start, d_end = [], []
    for (i, j), dc in zip(pairs, d_cur):
        d0 = float(model.d0(i, j))
        if (i, j) in model.formed:
            d_s = round(min(max(dc, d0 + 1.0), d0 + 2.5), 3)
            d_e = round(d0, 3)
        else:
            d_s = round(dc if dc < 1.35 * d0 else d0, 3)
            d_e = round(max(d_s + 0.42 * d0, 1.55 * d0), 3)
        d_start.append(d_s)
        d_end.append(d_e)

    if scan_start is not None:
        try:
            val = float(scan_start)
            if val > 0.2:
                d_start = np.array([val] * len(pairs))
        except (ValueError, TypeError):
            pass

    if scan_end is not None:
        try:
            val = float(scan_end)
            if val > 0.2:
                d_end = np.array([val] * len(pairs))
        except (ValueError, TypeError):
            pass

    return np.array(d_start), np.array(d_end)


def run_drive_sella_ts(atoms, model, log, n_points=None, k=60.0, scan_start=None, scan_end=None):
    """Strategy: Constrained Coordinate Drive across barrier peak -> Sella."""
    from ase.optimize import LBFGS
    solvent = getattr(model, "solvent", "water")
    log(f"\n=== Strategy: Constrained Coordinate Drive in {solvent.upper()} -> Sella ===")

    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    d_start, d_end = _build_drive_targets(model, pos, pairs, scan_start=scan_start, scan_end=scan_end)
    if not n_points or int(n_points) <= 0:
        if len(d_start) > 0:
            n_points = int(max(15, min(40, round(abs(d_end[0] - d_start[0]) / 0.04) + 1)))
        else:
            n_points = 21
    else:
        n_points = int(n_points)
    spec_restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

    frames, energies, lams, dist_list = [], [], [], []
    i_max = 0

    # Ensure current starting positions match d_start on reacting pairs
    init_pos = atoms.get_positions().copy()
    for idx, (i, j) in enumerate(pairs):
        vec = init_pos[j] - init_pos[i]
        norm = np.linalg.norm(vec)
        if norm > 1e-6:
            init_pos[j] = init_pos[i] + (vec / norm) * d_start[idx]
    atoms.set_positions(init_pos)

    log(f"  Scanning reaction coordinate across {n_points} points (d: {d_start[0]:.2f} -> {d_end[0]:.2f} A)...")

    from ase.constraints import FixBondLength

    for n, lam in enumerate(np.linspace(0.0, 1.0, n_points)):
        targets = d_start + lam * (d_end - d_start)

        # Clear prior constraints to allow setting new distance
        atoms.set_constraint(None)

        # Enforce exact target distances along reacting coordinate(s)
        cur_pos = atoms.get_positions().copy()
        for idx, (i, j) in enumerate(pairs):
            vec = cur_pos[j] - cur_pos[i]
            norm = np.linalg.norm(vec)
            if norm > 1e-6:
                cur_pos[j] = cur_pos[i] + (vec / norm) * targets[idx]
        atoms.set_positions(cur_pos)

        # Fix reacting bond lengths rigidly while relaxing all orthogonal DOFs
        calc = _make_restrained_calc(
            base,
            [],
            0.0,
            spectator_restraints=spec_restraints,
        )
        atoms.calc = calc
        atoms.set_constraint([FixBondLength(i, j) for i, j in pairs])
        opt = LBFGS(atoms, logfile=None, maxstep=0.08)
        opt.run(fmax=0.04, steps=45)

        e = calc.e_base
        cur_p = atoms.get_positions().copy()
        frames.append(cur_p)
        energies.append(e)
        lams.append(float(lam))
        d_val = float(np.linalg.norm(cur_p[pairs[0][0]] - cur_p[pairs[0][1]])) if pairs else 0.0
        dist_list.append(d_val)

        if e > energies[i_max]:
            i_max = n

        rel = (e - energies[0]) * EV_TO_KCAL
        log(f"  [{n+1:2d}/{n_points}] d = {d_val:.3f} A  dE = {rel:+6.2f} kcal/mol")

        # Early break if we have passed the peak by > 2.5 kcal/mol after at least 3 points
        if i_max >= 2 and (energies[i_max] - e) * EV_TO_KCAL > 2.5:
            log("  ✓ Passed barrier maximum - stopping scan early.")
            break

    # Identify saddle peak and reactant minimum well
    i_min = int(np.argmin(energies[:i_max + 1])) if i_max > 0 else 0
    barrier_kcal = (energies[i_max] - energies[i_min]) * EV_TO_KCAL
    d_peak = dist_list[i_max] if dist_list else d_start[0]
    d_well = dist_list[i_min] if dist_list else d_start[0]

    atoms.set_positions(frames[i_max])
    atoms.set_constraint(None)
    log(f"  ✓ Saddle point located at d = {d_peak:.3f} A (Well at {d_well:.3f} A, ΔE‡ = {barrier_kcal:.2f} kcal/mol).")
    log("  Polishing saddle point with guarded Sella (max drift 0.05 A)...")

    from .calculators import guard_saddle_polish
    guard_saddle_polish(atoms, model, log, max_steps=40, max_drift=0.05)

    profile = [[l, (e - energies[i_min]) * EV_TO_KCAL] for l, e in zip(lams, energies)]
    return profile


# ---------------------------------------------------------------------------
# Strategy 2: Climbing Image NEB (CI-NEB with IDPP)
# ---------------------------------------------------------------------------
def _build_relaxed_product_endpoint(atoms, model):
    """Build product endpoint P via true Product conformer Kabsch SO(3) alignment."""
    try:
        from .lego_builder import _embed_aligned_products
        pos_P = _embed_aligned_products(model, atoms.get_positions())
        if pos_P is not None and not np.allclose(pos_P, 0.0):
            return pos_P
    except Exception:
        pass

    from ase.optimize import LBFGS
    solvent = getattr(model, "solvent", "water")
    p_atoms = atoms.copy()
    p_targets = [(i, j, model.d0(i, j)) for i, j in model.formed]
    p_targets += [(i, j, model.d0(i, j) + 1.25) for i, j in model.broken]
    spec_restraints = [(i, j, float(np.linalg.norm(atoms.positions[i] - atoms.positions[j])))
                       for i, j in model.spectator_bonds]

    p_base = p_atoms.copy()
    p_base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    p_atoms.calc = _make_restrained_calc(p_base, p_targets, k=30.0, spectator_restraints=spec_restraints)
    p_opt = LBFGS(p_atoms, logfile=None, maxstep=0.2)
    p_opt.run(fmax=0.08, steps=50)
    return p_atoms.get_positions()


def run_neb_ts(atoms, model, log, n_images=9):
    """Strategy: Climbing Image Nudged Elastic Band (CI-NEB with IDPP)."""
    from ase.optimize import FIRE
    from ase.mep import NEB
    from ase.mep.neb import idpp_interpolate
    solvent = getattr(model, "solvent", "water")
    log(f"\n=== Strategy: Climbing Image NEB (CI-NEB) in {solvent.upper()} ===")

    pos_p = _build_relaxed_product_endpoint(atoms, model)

    images = [atoms.copy() for _ in range(n_images)]
    images[0].set_positions(atoms.get_positions())
    images[-1].set_positions(pos_p)
    idpp_interpolate(images, mic=False)

    for img in images:
        img.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

    neb = NEB(images, climb=True, method='improvedtangent')
    opt = FIRE(neb, logfile=None)
    opt.run(fmax=0.08, steps=120)

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

    solvent = getattr(model, "solvent", "water")
    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    targets = [(i, j, model.d0(i, j) + 0.65) for i, j in model.formed]
    spec = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]
    atoms.calc = _make_restrained_calc(base, targets, 35.0, spectator_restraints=spec, model=model)
    pre_opt = LBFGS(atoms, logfile=None, maxstep=0.2)
    pre_opt.run(fmax=0.1, steps=25)
    atoms.set_positions(atoms.get_positions())
    log("  Pre-contraction complete. Initializing Dimer saddle search...")


def run_dimer_ts(atoms, model, log, steps=100, fmax=0.04):
    """Strategy: Dimer Method (ase.mep.dimer)."""
    from ase.mep.dimer import DimerControl, MinModeAtoms, MinModeTranslate
    solvent = getattr(model, "solvent", "water")
    log(f"\n=== Strategy: Dimer Method (ase.mep.dimer) in {solvent.upper()} ===")
    log("Dynamically searching and following lowest-curvature reaction mode using force replicas...")

    _precondition_dimer_guess(atoms, model, log)

    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
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
    _polish_with_sella(atoms, model, log, max_steps=60)
    return True


# ---------------------------------------------------------------------------
# Strategy 4: Pre-Relaxation -> Sella
# ---------------------------------------------------------------------------
def run_seeded_sella_ts(atoms, model, log):
    """Strategy: Constrained Pre-Relaxation -> Sella."""
    from ase.optimize import LBFGS
    solvent = getattr(model, "solvent", "water")
    log(f"\n=== Strategy: Pre-Relaxation in {solvent.upper()} -> Sella ===")
    pairs = model.formed + model.broken
    pos = atoms.get_positions()
    restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in pairs]
    spec_restraints = [(i, j, float(np.linalg.norm(pos[i] - pos[j]))) for i, j in model.spectator_bonds]

    base = atoms.copy()
    base.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    atoms.calc = _make_restrained_calc(base, restraints, 25.0, spectator_restraints=spec_restraints)
    opt = LBFGS(atoms, logfile=None, maxstep=0.15)
    opt.run(fmax=0.06, steps=120)
    log(f"  Pre-relaxation complete. Relaxed energy = {atoms.calc.e_base:.4f} eV")

    _polish_with_sella(atoms, model, log, max_steps=150)
    return []
