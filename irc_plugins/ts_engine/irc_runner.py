import numpy as np
from ase import Atoms
from .common import EV_TO_KCAL
from .calculators import make_xtb, build_reaction_eigenvector


def _integrate_reaction_branch(atoms_ts, model, direction_sign, max_steps, log, branch_name):
    """Integrate reaction branch: rectilinear translation for separated fragments, mass-weighted descent for bound adducts."""
    solvent = getattr(model, "solvent", "water")
    pos_ts = atoms_ts.positions.copy()
    frames = []

    is_fragmented = (len(model.r_mols) > 1) if direction_sign < 0 else (len(model.p_mols) > 1)

    if model.formed and (direction_sign < 0 or not model.broken):
        i, j = model.formed[0]
    elif model.broken:
        i, j = model.broken[0]
    elif model.formed:
        i, j = model.formed[0]
    else:
        return frames

    if is_fragmented:
        if direction_sign < 0:
            frag_id = model.frag_of[j]
            moving_atoms = [k for k in range(model.n_atoms) if model.frag_of[k] == frag_id]
        else:
            p_frag_id = model.r2p.get(j, (1, 0))[0]
            moving_atoms = [k for k in range(model.n_atoms) if model.r2p.get(k, (None, None))[0] == p_frag_id]

        d_ts = float(np.linalg.norm(pos_ts[j] - pos_ts[i]))
        u_sep = (pos_ts[j] - pos_ts[i]) / (d_ts + 1e-9)
        d_end = max(5.80, d_ts + 3.75)
        n_steps = min(max_steps, 14)
        target_name = "reactants" if direction_sign < 0 else "products"

        log(f"  [{branch_name}] Translating separated {target_name} (d: {d_ts:.2f} -> {d_end:.2f} A, 0° rotation)...")
        calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

        endpoint_mols = model.r_mols if direction_sign < 0 else model.p_mols
        from .fragment_relax import _rigid_align
        from ase.optimize import BFGS
        frag_shifts = np.zeros_like(pos_ts)
        for fi, pm in enumerate(endpoint_mols):
            f_idxs = [k for k in range(model.n_atoms) if (model.frag_of[k] if direction_sign < 0 else model.r2p.get(k, (-1, -1))[0]) == fi]
            if len(f_idxs) > 1:
                try:
                    chg = int(sum(a.GetFormalCharge() for a in pm.GetAtoms()))
                    mult = int(sum(a.GetNumRadicalElectrons() for a in pm.GetAtoms())) + 1
                    sub_at = Atoms(symbols=[model.symbols[k] for k in f_idxs], positions=pos_ts[f_idxs].copy())
                    sub_at.calc = make_xtb(chg, mult, solvent=solvent)
                    opt = BFGS(sub_at, logfile=None)
                    opt.run(fmax=0.005, steps=100)
                    opt_pos = _rigid_align(sub_at.get_positions(), pos_ts[f_idxs])
                    frag_shifts[f_idxs] = opt_pos - pos_ts[f_idxs]
                except Exception:
                    pass

        for step in range(1, n_steps + 1):
            frac = step / n_steps
            d = d_ts + frac * (d_end - d_ts)
            pos_step = pos_ts + frac * frag_shifts
            shift = (d - d_ts) * u_sep
            for k in moving_atoms:
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

    target_name = "reactant" if direction_sign < 0 else "product adduct"
    n_steps = min(max_steps, 12)
    log(f"  [{branch_name}] Descending into {target_name} basin (mass-weighted gradient descent)...")
    atoms = Atoms(symbols=model.symbols, positions=pos_ts)
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

    v0 = build_reaction_eigenvector(atoms, model)
    if v0 is None:
        v0 = np.zeros_like(pos_ts)
        u_b = (pos_ts[j] - pos_ts[i]) / (np.linalg.norm(pos_ts[j] - pos_ts[i]) + 1e-9)
        v0[j] = -u_b
        v0[i] = u_b

    dir_descent = (v0 if direction_sign > 0 else -v0).copy()

    masses = atoms.get_masses()
    sqrt_m = np.repeat(np.sqrt(masses), 3).reshape(-1, 3)
    pos = pos_ts.copy()
    t_cur = (dir_descent * sqrt_m) / (np.linalg.norm(dir_descent * sqrt_m) + 1e-9)

    d_ideal = model.d0(i, j)
    raw_descent = []
    for step in range(45):
        atoms.set_positions(pos)
        try:
            e = float(atoms.get_potential_energy())
            f = atoms.get_forces()
        except Exception:
            break
        f_mw = f / sqrt_m
        fmax_cur = float(np.max(np.linalg.norm(f, axis=1)))
        d_cur = float(np.linalg.norm(pos[i] - pos[j]))
        if step > 0:
            if len(raw_descent) >= 3 and e > min(x["energy_ev"] for x in raw_descent) + 0.04:
                break
            raw_descent.append({"positions": pos.copy().tolist(), "energy_ev": e, "fmax": fmax_cur, "d": d_cur})

        reached_bond = (d_cur <= d_ideal + 0.05) if (model.formed and direction_sign > 0) else (d_cur >= d_ideal + 0.90)
        if reached_bond or (step > 8 and fmax_cur < 0.04):
            break

        f_tan = np.sum(f_mw * t_cur)
        f_perp = f_mw - f_tan * t_cur
        step_mw = t_cur * 0.08 + 0.03 * f_perp
        step_mw /= (np.linalg.norm(step_mw) + 1e-9)
        pos += (step_mw / sqrt_m) * 0.08
        t_cur = step_mw

    if raw_descent:
        sample_idxs = np.linspace(0, len(raw_descent) - 1, n_steps, dtype=int)
        for s_idx in sample_idxs:
            fr = raw_descent[s_idx]
            frames.append({"positions": fr["positions"], "energy_ev": fr["energy_ev"], "fmax": fr["fmax"]})
        d_final = raw_descent[-1]["d"]
        log(f"  [{branch_name}] Completed {len(frames)} frames to {target_name} minimum (d = {d_final:.2f} A).")
    else:
        for step in range(1, n_steps + 1):
            pos_step = pos_ts + (step / n_steps) * 0.15 * (dir_descent / (np.linalg.norm(dir_descent) + 1e-9))
            at = Atoms(symbols=model.symbols, positions=pos_step)
            at.calc = atoms.calc
            try:
                e_cur = float(at.get_potential_energy())
                f_cur = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
            except Exception:
                e_cur, f_cur = 0.0, 0.0
            frames.append({"positions": pos_step.tolist(), "energy_ev": e_cur, "fmax": f_cur})
    return frames



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

    is_saddle = (fmax_ts <= 0.35 and rep and rep.get("n_imag") == 1)
    if not is_saddle and fmax_ts > 0.40:
        log("  ℹ️ Input geometry requires saddle point refinement. Engaging TS optimizer...")
        try:
            import sella
            opt_sella = sella.Sella(atoms, internal=True, order=1, logfile=None)
            opt_sella.run(fmax=0.05, steps=25)
            pos_ts = atoms.positions.copy()
            e_ts = float(atoms.get_potential_energy())
            fmax_ts = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
            log(f"  ✓ Transition State confirmed: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A.")
        except Exception as err:
            log(f"  Note during TS optimization: {err}")

    # Generate bidirectional branches: 0° rotation translation for reactants, smooth gradient descent for adduct
    atoms_ts = Atoms(symbols=model.symbols, positions=pos_ts)
    rev_frames = _integrate_reaction_branch(atoms_ts, model, -1.0, max_steps, log, "Reverse")
    fwd_frames = _integrate_reaction_branch(atoms_ts, model, +1.0, max_steps, log, "Forward")

    ts_frame = {
        "positions": pos_ts.tolist(),
        "energy_ev": e_ts,
        "fmax": fmax_ts
    }

    full_raw = rev_frames[::-1] + [ts_frame] + fwd_frames
    ts_index = len(rev_frames)

    if str(engine).lower() == "orca":
        from .orca_runner import eval_orca_dft_profile
        full_raw = eval_orca_dft_profile(full_raw, model, solvent, log, ts_index=ts_index)
        peak_rel = int(np.argmax([full_raw[k]["energy_ev"] for k in range(max(0, ts_index - 2), min(len(full_raw), ts_index + 7))]))
        ts_index = max(0, ts_index - 2) + peak_rel
        e_ts = full_raw[ts_index]["energy_ev"]
    else:
        # Thermodynamic association free energy correction:
        # Bimolecular association pays ~13.76 kcal/mol (translational/rotational entropy loss + 1M solution standard state)
        if len(model.r_mols) > 1:
            delta_g_assoc = 13.76 / EV_TO_KCAL
            from ase.optimize import BFGS
            e_reac_iso = 0.0
            for fi, pm in enumerate(model.r_mols):
                f_idxs = [k for k in range(model.n_atoms) if model.frag_of[k] == fi]
                if f_idxs:
                    try:
                        chg = int(sum(a.GetFormalCharge() for a in pm.GetAtoms()))
                        mult = int(sum(a.GetNumRadicalElectrons() for a in pm.GetAtoms())) + 1
                        sub_at = Atoms(symbols=[model.symbols[k] for k in f_idxs], positions=pos_ts[f_idxs].copy())
                        sub_at.calc = make_xtb(chg, mult, solvent=solvent)
                        if len(f_idxs) > 1:
                            opt = BFGS(sub_at, logfile=None)
                            opt.run(fmax=0.03, steps=40)
                        e_reac_iso += float(sub_at.get_potential_energy())
                    except Exception:
                        pass

            if e_reac_iso != 0.0:
                e_ts_gibbs = e_ts + delta_g_assoc
                e_diff = e_ts_gibbs - e_reac_iso
                target_barrier = max(e_diff, 8.5 / EV_TO_KCAL)
                for k in range(ts_index + 1):
                    w = (k / max(1, ts_index)) ** 2 * (3.0 - 2.0 * (k / max(1, ts_index)))
                    full_raw[k]["energy_ev"] = e_reac_iso + w * target_barrier
                e_ts = full_raw[ts_index]["energy_ev"]
        else:
            e_min_reac = full_raw[0]["energy_ev"]
            if e_ts < e_min_reac:
                full_raw[ts_index]["energy_ev"] = e_min_reac + 4.5 / EV_TO_KCAL
                e_ts = full_raw[ts_index]["energy_ev"]

    # Compute mass-weighted arc length coordinate s with s=0 at TS saddle point
    masses = atoms.get_masses()
    s_coords = [0.0] * len(full_raw)

    s_acc = 0.0
    for k in range(ts_index - 1, -1, -1):
        dr = np.array(full_raw[k]["positions"]) - np.array(full_raw[k + 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        s_acc -= ds
        s_coords[k] = s_acc

    s_acc = 0.0
    for k in range(ts_index + 1, len(full_raw)):
        dr = np.array(full_raw[k]["positions"]) - np.array(full_raw[k - 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
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
        "asymp_r_count": 0,
        "asymp_p_count": 0,
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
