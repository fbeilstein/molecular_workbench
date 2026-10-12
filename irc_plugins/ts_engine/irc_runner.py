import numpy as np
from ase import Atoms
from .common import EV_TO_KCAL
from .calculators import make_xtb, build_reaction_eigenvector


def _step_mass_weighted_descent(atoms, pos_init, dir_init, max_steps, dx, fmax, log, branch_name, model=None, direction_sign=1.0):
    """Integrate mass-weighted gradient descent trajectory with tangent continuity."""
    masses = atoms.get_masses()
    sqrt_m = np.repeat(np.sqrt(masses), 3).reshape(-1, 3)

    pos = pos_init.copy()
    frames = []

    v_init = (dir_init * sqrt_m).copy()
    norm_v = np.linalg.norm(v_init)
    if norm_v > 1e-12:
        v_init /= norm_v
    t_cur = v_init.copy()

    formed_pair = model.formed[0] if (model and model.formed) else None

    for step in range(max_steps):
        atoms.set_positions(pos)
        try:
            e = float(atoms.get_potential_energy())
            f = atoms.get_forces()
        except Exception as err:
            log(f"  [{branch_name}] Calculation error at step {step}: {err}")
            break

        f_mw = f / sqrt_m
        norm_f_mw = np.linalg.norm(f_mw)
        fmax_cur = float(np.max(np.linalg.norm(f, axis=1)))

        frames.append({
            "positions": pos.copy().tolist(),
            "energy_ev": e,
            "fmax": fmax_cur,
            "norm_f_mw": float(norm_f_mw)
        })

        if step % 5 == 0:
            log(f"  [{branch_name}] Step {step + 1:2d}: E = {e:.4f} eV, max|F| = {fmax_cur:.4f} eV/A")

        # Project force perpendicular to tangent for valley relaxation
        f_tan = np.sum(f_mw * t_cur)
        f_perp = f_mw - f_tan * t_cur
        f_perp_norm = np.linalg.norm(f_perp)

        # Tangent step + gentle orthogonal relaxation (never reverse direction)
        relax_weight = 0.04 if f_perp_norm < 5.0 else 0.01
        step_mw = t_cur * dx + relax_weight * f_perp
        norm_step = np.linalg.norm(step_mw)
        if norm_step > 1e-12:
            step_mw /= norm_step

        # Enforce tangent continuity (stay within hemisphere of branch direction)
        if np.sum(step_mw * v_init) < 0.2:
            step_mw = v_init.copy()

        step_cart = (step_mw / sqrt_m) * dx
        pos += step_cart
        t_cur = step_mw

        # Physical termination checks
        if formed_pair:
            i_f, j_f = formed_pair
            d_form = np.linalg.norm(pos[i_f] - pos[j_f])
            if direction_sign < 0 and d_form >= 3.4:
                log(f"  ✓ [{branch_name}] Reached separated reactants (d = {d_form:.2f} A).")
                break
            if direction_sign > 0 and (d_form <= 1.48 or (step > 8 and fmax_cur < fmax)):
                log(f"  ✓ [{branch_name}] Reached product basin (d = {d_form:.2f} A, max|F| = {fmax_cur:.4f} eV/A).")
                break
        elif fmax_cur < fmax and step > 4:
            log(f"  ✓ [{branch_name}] Basin reached in {step + 1} steps (max|F| = {fmax_cur:.4f} eV/A).")
            break

    return frames


def _run_sella_branch(atoms_in, model, direction, max_steps, dx, fmax, log):
    """Attempt Sella IRC for a single direction with inner loop guard."""
    try:
        import sella
        atoms = atoms_in.copy()
        atoms.calc = make_xtb(model.charge, model.multiplicity)
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


def run_irc_path(atoms, model, log, max_steps=35, dx=0.08, fmax=0.03):
    """Compute bidirectional Intrinsic Reaction Coordinate (IRC) downhill from Transition State."""
    log("=== Intrinsic Reaction Coordinate (IRC) Integration ===")
    atoms.calc = make_xtb(model.charge, model.multiplicity)

    pos_ts = atoms.positions.copy()
    e_ts = float(atoms.get_potential_energy())
    fmax_ts = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    log(f"  Reference Geometry: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A")

    # If non-stationary, check for genuine saddle point on a copy without mutating atoms
    if fmax_ts > 0.08:
        log(f"  ℹ️ Geometry has non-zero gradient (max|F| = {fmax_ts:.3f} eV/A > 0.05 eV/A). Checking for saddle point...")
        try:
            import sella
            from ..ts_builder import analyse_mode
            atoms_check = atoms.copy()
            atoms_check.calc = make_xtb(model.charge, model.multiplicity)
            opt_sella = sella.Sella(atoms_check, order=1, logfile=None)
            opt_sella.run(fmax=0.05, steps=15)
            f_opt = atoms_check.get_forces()
            fmax_opt = float(np.max(np.linalg.norm(f_opt, axis=1)))
            rep = analyse_mode(atoms_check, model, lambda s='': None)
            if opt_sella.converged() and rep and rep.get("n_imag") == 1:
                pos_ts = atoms_check.positions.copy()
                e_ts = float(atoms_check.get_potential_energy())
                fmax_ts = fmax_opt
                log(f"  ✓ Saddle geometry confirmed with Sella: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A (1 imag freq).")
            else:
                log(f"  ℹ️ No transition state saddle point found on this PES region (max|F| = {fmax_ts:.3f} eV/A). Tracing directed reaction coordinate from current geometry.")
        except Exception:
            pass

    # Reaction eigenvector for directed bifurcation
    v0 = build_reaction_eigenvector(atoms, model)
    if v0 is None:
        v0 = np.random.normal(size=atoms.positions.shape)
        v0 /= np.linalg.norm(v0)
        log("  Using randomized initial reaction displacement mode.")
    else:
        log("  Seeded reaction direction along forming/breaking bonds.")

    # 1. Reverse Branch: Towards Reactants (-v0)
    log("\n  Integrating Reverse Path (towards Reactants)...")
    rev_frames = _run_sella_branch(atoms, model, 'reverse', max_steps, dx, fmax, log)

    # 2. Forward Branch: Towards Products (+v0)
    log("\n  Integrating Forward Path (towards Products)...")
    fwd_frames = _run_sella_branch(atoms, model, 'forward', max_steps, dx, fmax, log)

    # Validate whether Sella properly bifurcated
    use_fallback = False
    if not rev_frames or not fwd_frames:
        use_fallback = True
    else:
        p_rev_end = np.array(rev_frames[-1]["positions"])
        p_fwd_end = np.array(fwd_frames[-1]["positions"])
        rmsd_ends = float(np.linalg.norm(p_rev_end - p_fwd_end))
        if rmsd_ends < 0.35:
            log(f"  ⚠️ Sella IRC collapsed into the same basin on both sides (RMSD = {rmsd_ends:.3f} A).")
            log("     Engaging mass-weighted directed reaction coordinate integration...")
            use_fallback = True

    if use_fallback:
        atoms_rev = atoms.copy()
        atoms_rev.calc = make_xtb(model.charge, model.multiplicity)
        rev_frames = _step_mass_weighted_descent(
            atoms_rev, pos_ts, -v0, max_steps, dx, fmax, log, "Reverse", model=model, direction_sign=-1.0
        )

        atoms_fwd = atoms.copy()
        atoms_fwd.calc = make_xtb(model.charge, model.multiplicity)
        fwd_frames = _step_mass_weighted_descent(
            atoms_fwd, pos_ts, v0, max_steps, dx, fmax, log, "Forward", model=model, direction_sign=+1.0
        )

    # 3. Assemble full continuous trajectory
    ts_frame = {
        "positions": pos_ts.tolist(),
        "energy_ev": e_ts,
        "fmax": fmax_ts
    }

    full_raw = rev_frames[::-1] + [ts_frame] + fwd_frames
    ts_index = len(rev_frames)

    # 4. Compute mass-weighted arc length coordinate s
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

    # Format frames list with relative energies in kcal/mol
    frames = []
    for k, item in enumerate(full_raw):
        e_cur = item["energy_ev"]
        rel_kcal = (e_cur - e_ts) * EV_TO_KCAL
        frames.append({
            "step": k - ts_index,
            "s_arc": round(s_coords[k], 4),
            "energy_ev": round(e_cur, 5),
            "rel_energy_kcal": round(rel_kcal, 3),
            "positions": item["positions"]
        })

    # 5. Extract reaction thermodynamics and kinetics
    e_r = frames[0]["energy_ev"]
    e_p = frames[-1]["energy_ev"]
    barrier_fwd = max(0.0, (e_ts - e_r) * EV_TO_KCAL)
    barrier_rev = max(0.0, (e_ts - e_p) * EV_TO_KCAL)
    reaction_e = (e_p - e_r) * EV_TO_KCAL

    p_start = np.array(frames[0]["positions"])
    p_end = np.array(frames[-1]["positions"])
    sep_dist = float(np.linalg.norm(p_start - p_end))

    log(f"\n✓ Reaction Path Complete: {len(frames)} frames ({len(rev_frames)} reverse + ref + {len(fwd_frames)} forward).")
    log(f"  Endpoint Separation (Reactant vs Product): {sep_dist:.2f} A")
    log(f"  Forward Barrier  ΔE‡_fwd = {barrier_fwd:.2f} kcal/mol")
    log(f"  Reverse Barrier  ΔE‡_rev = {barrier_rev:.2f} kcal/mol")
    log(f"  Reaction Energy  ΔE_rxn  = {reaction_e:+.2f} kcal/mol")

    bonds = [list(pair) for pair in sorted(model.r_bonds | model.p_bonds)]

    result_json = {
        "reaction_smiles": f"{'.'.join(model.r_smiles)}>>{'.'.join(model.p_smiles)}",
        "symbols": model.symbols,
        "n_atoms": model.n_atoms,
        "bonds": bonds,
        "ts_frame_index": ts_index,
        "barriers": {
            "e_ts_ev": round(e_ts, 5),
            "e_reactant_ev": round(e_r, 5),
            "e_product_ev": round(e_p, 5),
            "forward_barrier_kcal": round(barrier_fwd, 2),
            "reverse_barrier_kcal": round(barrier_rev, 2),
            "reaction_energy_kcal": round(reaction_e, 2),
            "endpoint_separation_A": round(sep_dist, 2)
        },
        "frames": frames
    }

    return result_json
