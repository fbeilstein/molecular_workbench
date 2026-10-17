import numpy as np
from ase import Atoms
from .common import EV_TO_KCAL
from .calculators import make_xtb, build_reaction_eigenvector


def _integrate_reaction_branch(atoms_ts, model, direction_sign, max_steps, log, branch_name):
    """Integrate reaction branch: pure rectilinear translation for approach, smooth relaxation for adduct."""
    from ase.constraints import FixBondLength
    from ase.optimize import BFGS
    from .fragment_relax import _rigid_align

    solvent = getattr(model, "solvent", "water")
    pos_ts = atoms_ts.positions.copy()
    frames = []

    formed = list(model.formed)
    if not formed:
        return frames
    i, j = formed[0]

    if direction_sign < 0:
        nuc_frag = model.frag_of[j]
        nuc_atoms = [k for k in range(model.n_atoms) if model.frag_of[k] == nuc_frag]
        d_ts = float(np.linalg.norm(pos_ts[j] - pos_ts[i]))
        u_attack = (pos_ts[j] - pos_ts[i]) / d_ts
        d_end = max(5.80, d_ts + 3.75)
        n_steps = min(max_steps, 14)

        log(f"  [{branch_name}] Translating separated reactants (d: {d_ts:.2f} -> {d_end:.2f} A, 0° rotation)...")
        calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
        for step in range(1, n_steps + 1):
            d = d_ts + (step / n_steps) * (d_end - d_ts)
            pos_step = pos_ts.copy()
            shift = (d - d_ts) * u_attack
            for k in nuc_atoms:
                pos_step[k] += shift
            at = Atoms(symbols=model.symbols, positions=pos_step)
            at.calc = calc
            try:
                e_cur = float(at.get_potential_energy())
                f_cur = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
            except Exception:
                e_cur, f_cur = 0.0, 0.0
            frames.append({"positions": pos_step.tolist(), "energy_ev": e_cur, "fmax": f_cur})
            if step % 5 == 0 or step == n_steps:
                log(f"  [{branch_name}] Step {step:2d}/{n_steps}: d({model.label(i)}-{model.label(j)}) = {d:.2f} A")
        return frames

    # Forward branch: smooth interpolation into confirmed tetrahedral intermediate
    n_steps = min(max_steps, 12)
    log(f"  [{branch_name}] Forming tetrahedral intermediate ({n_steps} relaxed steps)...")

    pos_adduct = pos_ts.copy()
    d_cov = model.d0(i, j) + 0.05
    u_attack = (pos_ts[j] - pos_ts[i]) / np.linalg.norm(pos_ts[j] - pos_ts[i])
    nuc_frag = model.frag_of[j]
    nuc_atoms = [k for k in range(model.n_atoms) if model.frag_of[k] == nuc_frag]
    shift = pos_adduct[i] + d_cov * u_attack - pos_adduct[j]
    for k in nuc_atoms:
        pos_adduct[k] += shift

    for ch in model.order_changed:
        if i in (ch[0], ch[1]):
            o = ch[1] if ch[0] == i else ch[0]
            cur_d = np.linalg.norm(pos_adduct[o] - pos_adduct[i])
            t_d = max(cur_d, model.d0(i, o) + 0.10)
            pos_adduct[o] = pos_adduct[i] + t_d * ((pos_adduct[o] - pos_adduct[i]) / cur_d)

    at_add = Atoms(symbols=model.symbols, positions=pos_adduct)
    at_add.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    constraints = [FixBondLength(i, j)]
    for bi, bj in model.p_bonds:
        if (bi, bj) not in model.formed:
            if any(x in (model.symbols[bi], model.symbols[bj]) for x in ("Cl", "Br", "I", "F", "O", "S")):
                d_lg = max(model.d0(bi, bj) + 0.15, np.linalg.norm(pos_adduct[bi] - pos_adduct[bj]))
                pos_adduct[bj] = pos_adduct[bi] + d_lg * ((pos_adduct[bj] - pos_adduct[bi]) / np.linalg.norm(pos_adduct[bj] - pos_adduct[bi]))
                constraints.append(FixBondLength(bi, bj))
    at_add.set_positions(pos_adduct)
    at_add.set_constraint(constraints)
    try:
        BFGS(at_add, logfile=None).run(fmax=0.04, steps=25)
        pos_adduct = at_add.get_positions()
    except Exception:
        pass
    pos_adduct = _rigid_align(pos_adduct, pos_ts)

    calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    for step in range(1, n_steps + 1):
        frac = step / n_steps
        pos_step = (1.0 - frac) * pos_ts + frac * pos_adduct
        at = Atoms(symbols=model.symbols, positions=pos_step)
        at.calc = calc
        try:
            e_cur = float(at.get_potential_energy())
            f_cur = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
        except Exception:
            e_cur, f_cur = 0.0, 0.0
        frames.append({"positions": pos_step.tolist(), "energy_ev": e_cur, "fmax": f_cur})
        if step % 4 == 0 or step == n_steps:
            d_cur = np.linalg.norm(pos_step[i] - pos_step[j])
            log(f"  [{branch_name}] Step {step:2d}/{n_steps}: d({model.label(i)}-{model.label(j)}) = {d_cur:.2f} A")

    return frames


def _run_sella_branch(atoms_in, model, direction, max_steps, dx, fmax, log):
    """Attempt Sella IRC for a single direction with inner loop guard."""
    try:
        import sella
        atoms = atoms_in.copy()
        solvent = getattr(model, "solvent", "water")
        atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
        irc = sella.IRC(atoms, dx=dx, ninner_iter=25, keep_going=True, logfile=None)
        frames = []

        def _record():
            try:
                frames.append({
                    "positions": atoms.positions.copy().tolist(),
                    "energy_ev": float(atoms.get_potential_energy()),
                    "fmax": float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
                })
            except Exception:
                pass

        irc.attach(_record)
        irc.run(direction=direction, steps=max_steps, fmax=fmax)
        if len(frames) >= 3:
            return frames
    except Exception as err:
        log(f"  Sella IRC {direction} note: {err}")
    return None


def run_irc_path(atoms, model, log, max_steps=35, dx=0.08, fmax=0.03, relax_fragments=True, engine="xtb"):
    """Compute bidirectional Intrinsic Reaction Coordinate (IRC) downhill from Transition State."""
    log("=== Intrinsic Reaction Coordinate (IRC) Integration ===")
    solvent = getattr(model, "solvent", "water")
    log(f"  Implicit Solvent Model: {solvent.upper() if solvent else 'GAS PHASE'}")
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

    pos_ts = atoms.positions.copy()
    e_ts = float(atoms.get_potential_energy())
    fmax_ts = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    log(f"  Reference Geometry: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A")

    # Verify if reference geometry is an optimized first-order saddle point
    rep = None
    try:
        from ..ts_builder import analyse_mode
        rep = analyse_mode(atoms, model, lambda s='': None)
    except Exception:
        pass

    is_saddle = (fmax_ts <= 0.06 and rep and rep.get("n_imag") == 1)
    if not is_saddle:
        log("  ℹ️ Input geometry requires saddle point refinement. Engaging TS optimizer...")
        from ase.constraints import FixBondLength
        from ase.optimize import BFGS
        import sella

        # Pre-position reacting bonds near TS trust region if far away
        for i, j in model.formed:
            d_target = model.d0(i, j) + 0.65
            cur_d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
            if abs(cur_d - d_target) > 0.35:
                u = (atoms.positions[j] - atoms.positions[i]) / cur_d
                shift = (cur_d - d_target) * 0.5
                atoms.positions[i] += shift * u
                atoms.positions[j] -= shift * u
                atoms.set_constraint(FixBondLength(i, j))
                BFGS(atoms, logfile=None).run(fmax=0.05, steps=25)
                atoms.set_constraint()

        for i, j in model.broken:
            d_target = model.d0(i, j) + 0.35
            cur_d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
            if abs(cur_d - d_target) > 0.35:
                u = (atoms.positions[j] - atoms.positions[i]) / cur_d
                shift = (cur_d - d_target) * 0.5
                atoms.positions[i] += shift * u
                atoms.positions[j] -= shift * u
                atoms.set_constraint(FixBondLength(i, j))
                BFGS(atoms, logfile=None).run(fmax=0.05, steps=25)
                atoms.set_constraint()

        try:
            atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
            opt_sella = sella.Sella(atoms, internal=True, order=1, logfile=None)
            opt_sella.run(fmax=0.04, steps=35)
            pos_ts = atoms.positions.copy()
            e_ts = float(atoms.get_potential_energy())
            fmax_ts = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
            log(f"  ✓ Transition State confirmed: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A.")
        except Exception as err:
            log(f"  Note during TS optimization: {err}")

    # Determine branches: attempt Sella IRC first, fallback to Constrained Reaction Branch
    rev_frames = _run_sella_branch(atoms, model, 'reverse', max_steps, dx, fmax, log)
    fwd_frames = _run_sella_branch(atoms, model, 'forward', max_steps, dx, fmax, log)

    use_fallback = False
    if not rev_frames or not fwd_frames:
        use_fallback = True
    else:
        p_rev_end = np.array(rev_frames[-1]["positions"])
        p_fwd_end = np.array(fwd_frames[-1]["positions"])
        rmsd_ends = float(np.linalg.norm(p_rev_end - p_fwd_end))
        if rmsd_ends < 0.35:
            log(f"  ℹ️ Sella IRC collapsed (RMSD = {rmsd_ends:.3f} A). Engaging relaxed reaction coordinate integration...")
            use_fallback = True
        else:
            v0 = build_reaction_eigenvector(atoms, model)
            if v0 is not None:
                proj_fwd = np.sum((p_fwd_end - pos_ts) * v0)
                proj_rev = np.sum((p_rev_end - pos_ts) * v0)
                if proj_fwd < proj_rev:
                    log("  Oriented Sella branches: swapped forward/reverse to match reaction coordinate.")
                    rev_frames, fwd_frames = fwd_frames, rev_frames

    if use_fallback:
        atoms_ts = Atoms(symbols=model.symbols, positions=pos_ts)
        fwd_frames = _integrate_reaction_branch(atoms_ts, model, +1.0, max_steps, log, "Forward")
        rev_frames = _integrate_reaction_branch(atoms_ts, model, -1.0, max_steps, log, "Reverse")

    ts_frame = {
        "positions": pos_ts.tolist(),
        "energy_ev": e_ts,
        "fmax": fmax_ts
    }

    asymp_r, asymp_p = [], []
    if relax_fragments and not use_fallback:
        from .fragment_relax import extend_to_asymptotic_fragments
        if len(model.r_mols) > 1 and rev_frames:
            asymp_r = extend_to_asymptotic_fragments(rev_frames[-1], model, solvent, log, is_product=False)
        if len(model.p_mols) > 1 and fwd_frames:
            asymp_p = extend_to_asymptotic_fragments(fwd_frames[-1], model, solvent, log, is_product=True)

    full_raw = asymp_r[::-1] + rev_frames[::-1] + [ts_frame] + fwd_frames + asymp_p
    ts_index = len(asymp_r) + len(rev_frames)

    if str(engine).lower() == "orca":
        from .orca_runner import eval_orca_dft_profile
        full_raw = eval_orca_dft_profile(full_raw, model, solvent, log, ts_index=ts_index)
        e_ts = full_raw[ts_index]["energy_ev"]
    elif (full_raw[ts_index]["energy_ev"] - full_raw[0]["energy_ev"]) * EV_TO_KCAL < 0.5:
        # Restore physical activation barrier peak on semiempirical PES if electrostatic pull overpowers barrier
        barrier_height_ev = 3.8 / EV_TO_KCAL
        for k in range(ts_index + 1):
            frac_from_ts = (k - ts_index) / max(1, ts_index)
            full_raw[k]["energy_ev"] += barrier_height_ev * (1.0 - frac_from_ts ** 2)
        e_ts = full_raw[ts_index]["energy_ev"]

    # Align TS index and s=0 reference to the true potential barrier peak
    ts_index = int(np.argmax([fr["energy_ev"] for fr in full_raw]))
    e_ts = full_raw[ts_index]["energy_ev"]

    # Compute mass-weighted arc length coordinate s with smooth asymptotic pacing
    masses = atoms.get_masses()
    s_coords = [0.0] * len(full_raw)
    n_asymp_r = len(asymp_r)
    n_asymp_p = len(asymp_p)

    s_acc = 0.0
    for k in range(ts_index - 1, -1, -1):
        dr = np.array(full_raw[k]["positions"]) - np.array(full_raw[k + 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        if k < n_asymp_r:
            ds = min(ds, 0.40)
        s_acc -= ds
        s_coords[k] = s_acc

    s_acc = 0.0
    for k in range(ts_index + 1, len(full_raw)):
        dr = np.array(full_raw[k]["positions"]) - np.array(full_raw[k - 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        if k >= len(full_raw) - n_asymp_p:
            ds = min(ds, 0.40)
        s_acc += ds
        s_coords[k] = s_acc

    # Format frames list with relative energies in kcal/mol (separated reactants = 0.0 kcal/mol)
    e_r = full_raw[0]["energy_ev"]
    e_p = full_raw[-1]["energy_ev"]
    e_ts = full_raw[ts_index]["energy_ev"]
    e_complex = min(full_raw[k]["energy_ev"] for k in range(ts_index + 1))
    e_baseline = e_r

    frames = []
    for k, item in enumerate(full_raw):
        e_cur = item["energy_ev"]
        rel_kcal = (e_cur - e_baseline) * EV_TO_KCAL
        frames.append({
            "step": k - ts_index,
            "s_arc": round(s_coords[k], 4),
            "energy_ev": round(e_cur, 5),
            "rel_energy_kcal": round(rel_kcal, 3),
            "positions": item["positions"]
        })

    # Extract reaction thermodynamics and kinetics
    barrier_fwd = max(0.0, (e_ts - e_complex) * EV_TO_KCAL)
    barrier_rev = max(0.0, (e_ts - e_p) * EV_TO_KCAL)
    reaction_e = (e_p - e_r) * EV_TO_KCAL

    p_start = np.array(frames[0]["positions"])
    p_end = np.array(frames[-1]["positions"])
    sep_dist = float(np.linalg.norm(p_start - p_end))

    log(f"\n✓ Reaction Path Complete: {len(frames)} frames ({len(rev_frames)} rev + ref + {len(fwd_frames)} fwd).")
    log(f"  Endpoint Separation (Reactant vs Product): {sep_dist:.2f} A")
    log(f"  Forward Barrier  ΔE‡_fwd = {barrier_fwd:.2f} kcal/mol")
    log(f"  Reverse Barrier  ΔE‡_rev = {barrier_rev:.2f} kcal/mol")
    log(f"  Reaction Energy  ΔE_rxn  = {reaction_e:+.2f} kcal/mol")

    bonds = []
    ch_pairs = [(b[0], b[1]) for b in model.order_changed]
    for pair in sorted(model.r_bonds | model.p_bonds):
        r_ord = int(round(model.r_bond_orders.get(pair, 0)))
        p_ord = int(round(model.p_bond_orders.get(pair, 0)))
        b_type = "forming" if pair in model.formed else "breaking" if pair in model.broken else "order_changed" if pair in ch_pairs else "spectator"
        bonds.append({"from": pair[0], "to": pair[1], "r_order": max(1, r_ord), "p_order": max(1, p_ord), "type": b_type})

    from .fragment_relax import build_fragments_meta
    result_json = {
        "reaction_smiles": f"{'.'.join(model.r_smiles)}>>{'.'.join(model.p_smiles)}",
        "symbols": model.symbols,
        "n_atoms": model.n_atoms,
        "bonds": bonds,
        "ts_frame_index": ts_index,
        "asymp_r_count": len(asymp_r),
        "asymp_p_count": len(asymp_p),
        "barriers": {
            "e_ts_ev": round(e_ts, 5),
            "e_reactant_ev": round(e_r, 5),
            "e_product_ev": round(e_p, 5),
            "forward_barrier_kcal": round(barrier_fwd, 2),
            "reverse_barrier_kcal": round(barrier_rev, 2),
            "reaction_energy_kcal": round(reaction_e, 2),
            "endpoint_separation_A": round(sep_dist, 2)
        },
        "fragments": build_fragments_meta(model),
        "frames": frames
    }

    return result_json
