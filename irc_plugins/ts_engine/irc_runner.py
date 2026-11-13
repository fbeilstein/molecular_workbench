import numpy as np
from ase import Atoms
from ase.constraints import FixBondLength
from ase.optimize import BFGS

from .common import EV_TO_KCAL
from .calculators import make_xtb, build_reaction_eigenvector


def _endpoint_frag_map(model, direction_sign: float):
    """Return per-atom fragment index mapping on the reactant (<0) or product (>0) side."""
    if direction_sign < 0:
        return {k: int(model.frag_of[k]) for k in range(model.n_atoms)}
    return {k: int(model.r2p.get(k, (-1, -1))[0]) for k in range(model.n_atoms)}


def _select_interfragment_pair(model, direction_sign: float, pos: np.ndarray):
    """Select the primary forming/breaking atom pair separating distinct endpoint fragments."""
    frag_map = _endpoint_frag_map(model, direction_sign)
    cands = (list(model.formed) + list(model.broken)) if direction_sign < 0 else (list(model.broken) + list(model.formed))
    for i, j in cands:
        if frag_map.get(i, -1) != frag_map.get(j, -1) and frag_map.get(i, -1) >= 0 and frag_map.get(j, -1) >= 0:
            return i, j

    # Fallback: closest pair across fragment 0 and fragment 1
    f0 = [k for k, fid in frag_map.items() if fid == 0]
    f1 = [k for k, fid in frag_map.items() if fid == 1]
    if f0 and f1:
        sub_d = np.linalg.norm(pos[f0][:, None, :] - pos[f1][None, :, :], axis=-1)
        a_idx, b_idx = np.unravel_index(int(np.argmin(sub_d)), sub_d.shape)
        return f0[int(a_idx)], f1[int(b_idx)]
    return (cands[0][0], cands[0][1]) if cands else (None, None)


def _select_bond_pair(model, direction_sign: float):
    """Select the primary forming/breaking bond pair for the given reaction branch."""
    if model.formed and (direction_sign < 0 or not model.broken):
        return model.formed[0]
    if model.broken:
        return model.broken[0]
    if model.formed:
        return model.formed[0]
    return None, None


def _topology_preserved(pos_old: np.ndarray, pos_new: np.ndarray, model) -> bool:
    """Guard TS optimization against collapsing into separated products/reactants."""
    pairs = list(model.formed) + list(model.broken)
    for i, j in pairs:
        d_old = float(np.linalg.norm(pos_old[i] - pos_old[j]))
        d_new = float(np.linalg.norm(pos_new[i] - pos_new[j]))
        if d_new > max(3.6, d_old * 1.55) or d_new < max(0.85, d_old * 0.65):
            return False
    return True


def _branch_progress(pos_ts: np.ndarray, pos_pt: np.ndarray, model) -> float:
    """Signed geometric progress relative to pos_ts (>0 toward products, <0 toward reactants)."""
    prog = 0.0
    for i, j in getattr(model, "formed", ()) or ():
        prog += float(np.linalg.norm(pos_ts[i] - pos_ts[j]) - np.linalg.norm(pos_pt[i] - pos_pt[j]))
    for i, j in getattr(model, "broken", ()) or ():
        prog += float(np.linalg.norm(pos_pt[i] - pos_pt[j]) - np.linalg.norm(pos_ts[i] - pos_ts[j]))
    return prog


def _fallback_rxn_direction(pos: np.ndarray, model) -> np.ndarray:
    """Fallback reaction tangent built from formed/broken bond vectors."""
    u = build_reaction_eigenvector(Atoms(symbols=model.symbols, positions=pos), model)
    if u is not None and np.linalg.norm(u) > 1e-8:
        return u
    pair = _select_bond_pair(model, +1.0)
    v = np.zeros_like(pos)
    if pair[0] is not None:
        i, j = pair
        sign = -1.0 if pair in model.formed else 1.0
        dv = pos[j] - pos[i]
        u_b = dv / max(float(np.linalg.norm(dv)), 1e-9)
        v[j] = sign * u_b
        v[i] = -sign * u_b
    return v


def _prepare_ts_and_mode(atoms: Atoms, model, solvent: str, dx: float, log):
    """Optionally polish TS saddle without topology drift and extract mass-weighted IRC mode."""
    pos_ts = atoms.positions.copy()
    e_ts = float(atoms.get_potential_energy())
    fmax_ts = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    log(f"  Reference Geometry: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A")

    # Gentle Sella refinement when TS still has significant residual gradient
    if fmax_ts > 0.015:
        steps_lim = 30 if fmax_ts <= 0.25 else 50
        target_f = 0.015 if fmax_ts <= 0.25 else 0.025
        log("  Refining Transition State saddle geometry (Sella order=1)...")
        try:
            import sella
            at_opt = Atoms(symbols=model.symbols, positions=pos_ts.copy())
            at_opt.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
            opt_sella = sella.Sella(at_opt, internal=True, order=1, logfile=None)
            opt_sella.run(fmax=target_f, steps=steps_lim)
            pos_opt = at_opt.get_positions().copy()
            fmax_opt = float(np.max(np.linalg.norm(at_opt.get_forces(), axis=1)))
            converged_or_local = fmax_opt <= 0.035 or (fmax_ts <= 0.20 and fmax_opt < fmax_ts)
            if converged_or_local and _topology_preserved(pos_ts, pos_opt, model):
                pos_ts = pos_opt
                e_ts = float(at_opt.get_potential_energy())
                fmax_ts = fmax_opt
                atoms.set_positions(pos_ts)
                log(f"  ✓ Transition State refined: E = {e_ts:.4f} eV, max|F| = {fmax_ts:.4f} eV/A")
            else:
                log("  Note: Kept input TS geometry (Sella refinement drifted or remained unconverged).")
        except Exception as err:
            log(f"  Note during TS refinement: {err}")

    # Compute exact Hessian eigenmodes once at TS and orient along reactants -> products
    sqrt_m = np.sqrt(atoms.get_masses())[:, None]
    u_rxn = _fallback_rxn_direction(pos_ts, model)
    u_rxn_mw = (u_rxn * sqrt_m).ravel()
    u_rxn_mw /= max(float(np.linalg.norm(u_rxn_mw)), 1e-12)

    v_step = None
    irc_dx = max(float(dx), 0.10)
    try:
        from .vibrations import _compute_vibrations
        at_vib = Atoms(symbols=model.symbols, positions=pos_ts.copy())
        at_vib.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
        freqs, modes = _compute_vibrations(at_vib)
        imag = [(k, float(f.imag)) for k, f in enumerate(freqs) if getattr(f, "imag", 0.0) > 20.0]
        if imag:
            log("  Imaginary frequencies: " + ", ".join(f"{fi:.0f}i" for _, fi in imag) + " cm-1")
            best_k, best_score = None, -1.0
            for k, fi in imag:
                vm = (np.real(modes[k]).reshape(-1, 3) * sqrt_m).ravel()
                ov = abs(float(np.dot(vm, u_rxn_mw))) / max(float(np.linalg.norm(vm)), 1e-12)
                if ov >= 0.20:
                    score = fi * ov
                    if score > best_score:
                        best_score = score
                        best_k = k
            if best_k is not None:
                v_mw = np.real(modes[best_k]).reshape(-1, 3) * sqrt_m
                v_mw /= max(float(np.linalg.norm(v_mw)), 1e-12)
                if float(np.dot(v_mw.ravel(), u_rxn_mw)) < 0.0:
                    v_mw = -v_mw
                v_step = (irc_dx * v_mw / sqrt_m).ravel()
            else:
                log("  Warning: Imaginary modes have <0.20 overlap with reactive bonds; using bond-vector tangent.")
        else:
            log("  Warning: No imaginary mode > 20i cm-1 found at TS; using bond-vector tangent.")
    except Exception as err:
        log(f"  Note during vibrational mode extraction: {err}")

    if v_step is None:
        v_step = (irc_dx * u_rxn_mw.reshape(-1, 3) / sqrt_m).ravel()

    return pos_ts, e_ts, fmax_ts, v_step, irc_dx


def _integrate_sella_branch(
    pos_ts: np.ndarray,
    e_ts: float,
    v_step: np.ndarray,
    model,
    direction_sign: float,
    d_target: float,
    max_steps: int,
    dx: float,
    fmax: float,
    log,
    branch_name: str,
):
    """Integrate downhill from TS using mode-seeded sella.IRC plus optional BFGS basin polish."""
    import sella

    solvent = getattr(model, "solvent", "water")
    direction = "reverse" if direction_sign < 0 else "forward"
    is_fragmented = (len(model.r_mols) > 1) if direction_sign < 0 else (len(model.p_mols) > 1)
    pair_i, pair_j = _select_interfragment_pair(model, direction_sign, pos_ts) if is_fragmented else _select_bond_pair(model, direction_sign)

    polish_fmax = min(float(fmax), 0.02)
    n_irc_steps = max(int(max_steps), 40)
    log(f"  [{branch_name}] Integrating downhill IRC (sella.IRC, dx={dx:.2f}, tol={polish_fmax:.3f} eV/A)...")

    at = Atoms(symbols=model.symbols, positions=pos_ts.copy())
    at.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

    pts = []
    try:
        irc = sella.IRC(at, logfile=None, dx=dx, ninner_iter=15, keep_going=True)
        irc.pes.kick(0, True, **irc.peskwargs)
        irc.H0 = irc.pes.get_H().asarray().copy()
        irc.pescurr = irc.pes.curr.copy()
        irc.peslast = irc.pes.last.copy()
        irc.v0ts = v_step.copy()

        prev_e = e_ts
        d_max_sep = 0.0
        flat_count = 0
        for step_i, _ in enumerate(irc.irun(fmax=polish_fmax, steps=n_irc_steps, direction=direction)):
            if step_i == 0:
                continue  # skip duplicate TS point
            p = at.get_positions().copy()
            if direction_sign * _branch_progress(pos_ts, p, model) < -0.02:
                log(f"  [{branch_name}] Stopped sella.IRC at step {step_i} (wrong-way drift toward opposite basin).")
                break
            e = float(at.get_potential_energy())
            fm = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
            de_kcal = (e - prev_e) * EV_TO_KCAL
            if de_kcal > 0.05:
                break
            pts.append({"positions": p, "energy_ev": e, "fmax": fm})
            prev_e = e

            if is_fragmented and pair_i is not None:
                d_cur = float(np.linalg.norm(p[pair_i] - p[pair_j]))
                if d_cur > d_max_sep:
                    d_max_sep = d_cur

                # If fragment cleaved (>2.8 A) and starts drifting back, terminate branch
                if d_max_sep > 2.80 and d_cur < (d_max_sep - 0.10):
                    log(f"  [{branch_name}] Stopped sella.IRC: fragment reached maximum separation ({d_max_sep:.2f} A) and began curling back.")
                    break

                # If energy flattens in asymptotic dissociation region, terminate branch
                if d_cur > 2.80 and abs(de_kcal) < 0.02:
                    flat_count += 1
                    if flat_count >= 3:
                        log(f"  [{branch_name}] Stopped sella.IRC: energy plateaued in asymptotic dissociation region (d = {d_cur:.2f} A).")
                        break
                else:
                    flat_count = 0

                if d_cur >= d_target - 0.1:
                    break
    except Exception as err:
        log(f"  [{branch_name}] sella.IRC note: {err}")

    if pts and max(direction_sign * _branch_progress(pos_ts, pt["positions"], model) for pt in pts) <= 0.01:
        pts = []

    # Fallback kick for bound unimolecular branches when sella.IRC yields nothing;
    # barrierless multi-fragment exits start from pos_ts via constrained separation instead.
    if not pts and not is_fragmented:
        sign = -1.0 if direction_sign < 0 else 1.0
        p0 = pos_ts + sign * v_step.reshape(-1, 3)
        at_fb = Atoms(symbols=model.symbols, positions=p0)
        at_fb.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
        try:
            pts.append({
                "positions": p0,
                "energy_ev": float(at_fb.get_potential_energy()),
                "fmax": float(np.max(np.linalg.norm(at_fb.get_forces(), axis=1))),
            })
        except Exception:
            pass

    # Polish bound unimolecular basin to converged gradient. Disconnected fragments must
    # NOT be polished unconstrained here (they are relaxed independently via fragment_relax).
    if pts and not is_fragmented:
        end_pos = pts[-1]["positions"].copy()
        end_e = pts[-1]["energy_ev"]
        end_fmax = pts[-1]["fmax"]
        allow_polish = (end_fmax > polish_fmax)

        if allow_polish and end_fmax > polish_fmax:
            at_end = Atoms(symbols=model.symbols, positions=end_pos)
            at_end.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
            try:
                opt = BFGS(at_end, logfile=None)
                last_saved = end_pos.copy()
                for _ in opt.irun(fmax=polish_fmax, steps=55):
                    p_cur = at_end.get_positions().copy()
                    if direction_sign * _branch_progress(pos_ts, p_cur, model) < -0.02:
                        log(f"  [{branch_name}] Stopped BFGS polish (wrong-way drift toward opposite basin).")
                        break
                    e_cur = float(at_end.get_potential_energy())
                    fm_cur = float(np.max(np.linalg.norm(at_end.get_forces(), axis=1)))
                    if e_cur <= end_e and float(np.sqrt(np.mean((p_cur - last_saved) ** 2))) >= 0.03:
                        pts.append({"positions": p_cur, "energy_ev": e_cur, "fmax": fm_cur})
                        last_saved = p_cur
                        end_e = e_cur
                p_fin = at_end.get_positions().copy()
                e_fin = float(at_end.get_potential_energy())
                fm_fin = float(np.max(np.linalg.norm(at_end.get_forces(), axis=1)))
                if (
                    direction_sign * _branch_progress(pos_ts, p_fin, model) >= -0.02
                    and e_fin <= end_e + 1e-5
                    and float(np.sqrt(np.mean((p_fin - pts[-1]["positions"]) ** 2))) > 0.005
                ):
                    pts.append({"positions": p_fin, "energy_ev": e_fin, "fmax": fm_fin})
            except Exception as err:
                log(f"  [{branch_name}] BFGS endpoint polish note: {err}")

    if pts:
        d_fin = float(np.linalg.norm(pts[-1]["positions"][pair_i] - pts[-1]["positions"][pair_j])) if pair_i is not None else 0.0
        log(f"  [{branch_name}] {len(pts)} IRC frames -> basin d = {d_fin:.2f} A (max|F| = {pts[-1]['fmax']:.4f} eV/A).")
    return pts


def _separate_fragments_constrained(
    enc_pos: np.ndarray,
    model,
    direction_sign: float,
    d_target: float,
    relax_fragments: bool,
    log,
    branch_name: str,
):
    """Extend multi-fragment endpoint from IRC encounter basin out to asymptotic separation."""
    frag_map = _endpoint_frag_map(model, direction_sign)
    pair_i, pair_j = _select_interfragment_pair(model, direction_sign, enc_pos)
    if pair_i is None or pair_j is None or frag_map.get(pair_i) == frag_map.get(pair_j):
        return []

    d_cur = float(np.linalg.norm(enc_pos[pair_j] - enc_pos[pair_i]))
    if d_cur >= d_target - 0.15:
        return []

    moving = [k for k, fid in frag_map.items() if fid == frag_map[pair_j]]
    solvent = getattr(model, "solvent", "water")
    calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)
    n_sep = int(np.clip(round((d_target - d_cur) / 0.45), 6, 10))
    mode_tag = "constrained BFGS" if relax_fragments else "rigid"
    log(f"  [{branch_name}] Separating fragments ({mode_tag}, d: {d_cur:.2f} -> {d_target:.2f} A, {n_sep} frames)...")

    sep_pts = []
    prev_pos = enc_pos.copy()
    for step_idx in range(1, n_sep + 1):
        d_k = d_cur + (step_idx / n_sep) * (d_target - d_cur)
        p_init = prev_pos.copy()
        vec = p_init[pair_j] - p_init[pair_i]
        d_now = float(np.linalg.norm(vec))
        shift = (d_k - d_now) * (vec / max(d_now, 1e-8))
        for k in moving:
            p_init[k] += shift

        at = Atoms(symbols=model.symbols, positions=p_init)
        at.calc = calc
        p_final = p_init.copy()
        if relax_fragments:
            try:
                at.set_constraint(FixBondLength(pair_i, pair_j))
                opt = BFGS(at, logfile=None)
                opt.run(fmax=0.04, steps=18)
                p_cand = at.get_positions().copy()
                # Check against nuclear clashes before accepting relaxed coordinates
                d_mat = np.linalg.norm(p_cand[:, None, :] - p_cand[None, :, :], axis=-1) + np.eye(model.n_atoms) * 99.0
                if float(np.min(d_mat)) > 0.70:
                    p_final = p_cand
            except Exception:
                pass
            finally:
                at.set_constraint(None)
                at.set_positions(p_final)

        try:
            e_val = float(at.get_potential_energy())
            f_val = float(np.max(np.linalg.norm(at.get_forces(), axis=1)))
        except Exception:
            e_val = sep_pts[-1]["energy_ev"] if sep_pts else 0.0
            f_val = 0.0

        sep_pts.append({"positions": p_final, "energy_ev": e_val, "fmax": f_val})
        prev_pos = p_final

    return sep_pts


def _format_raw_frames(pts):
    return [{"positions": np.asarray(p["positions"]).tolist(), "energy_ev": float(p["energy_ev"]), "fmax": float(p.get("fmax", 0.0))} for p in pts]


def run_irc_path(
    atoms,
    model,
    log,
    max_steps=35,
    dx=0.08,
    fmax=0.03,
    relax_fragments=True,
    engine="xtb",
    **kwargs,
):
    """Compute bidirectional Intrinsic Reaction Coordinate (IRC) downhill from Transition State."""
    log("=== Intrinsic Reaction Coordinate (IRC) Integration ===")
    solvent = getattr(model, "solvent", "water")
    log(f"  Implicit Solvent Model: {solvent.upper() if solvent else 'GAS PHASE'}")
    atoms = Atoms(symbols=model.symbols, positions=atoms.positions.copy())
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=solvent)

    pos_ts, e_ts, fmax_ts, v_step, irc_dx = _prepare_ts_and_mode(atoms, model, solvent, dx, log)

    # Determine asymptotic separation targets when reactants or products have >= 2 fragments
    r_pair = _select_interfragment_pair(model, -1.0, pos_ts) if len(model.r_mols) > 1 else (None, None)
    p_pair = _select_interfragment_pair(model, +1.0, pos_ts) if len(model.p_mols) > 1 else (None, None)
    d_target_r = max(6.20, float(np.linalg.norm(pos_ts[r_pair[0]] - pos_ts[r_pair[1]])) + 3.20) if r_pair[0] is not None else 6.20
    d_target_p = max(6.20, float(np.linalg.norm(pos_ts[p_pair[0]] - pos_ts[p_pair[1]])) + 3.20) if p_pair[0] is not None else 6.20
    d_target_r = min(d_target_r, 7.50)
    d_target_p = min(d_target_p, 7.50)

    rev_irc = _integrate_sella_branch(pos_ts, e_ts, v_step, model, -1.0, d_target_r, max_steps, irc_dx, fmax, log, "Reverse")
    fwd_irc = _integrate_sella_branch(pos_ts, e_ts, v_step, model, +1.0, d_target_p, max_steps, irc_dx, fmax, log, "Forward")

    rev_sep = []
    if len(model.r_mols) > 1:
        enc_r = rev_irc[-1]["positions"] if rev_irc else pos_ts
        rev_sep = _separate_fragments_constrained(enc_r, model, -1.0, d_target_r, relax_fragments, log, "Reverse")

    fwd_sep = []
    if len(model.p_mols) > 1:
        enc_p = fwd_irc[-1]["positions"] if fwd_irc else pos_ts
        fwd_sep = _separate_fragments_constrained(enc_p, model, +1.0, d_target_p, relax_fragments, log, "Forward")

    ts_frame = {"positions": pos_ts.tolist(), "energy_ev": e_ts, "fmax": fmax_ts}
    rev_sep_raw = _format_raw_frames(rev_sep[::-1])
    rev_irc_raw = _format_raw_frames(rev_irc[::-1])
    fwd_irc_raw = _format_raw_frames(fwd_irc)
    fwd_sep_raw = _format_raw_frames(fwd_sep)

    # Path order: separated reactants -> encounter R -> TS -> encounter P -> separated products
    full_raw = rev_sep_raw + rev_irc_raw + [ts_frame] + fwd_irc_raw + fwd_sep_raw
    ts_index = len(rev_sep_raw) + len(rev_irc_raw)
    asymp_r_count = len(rev_sep_raw)
    asymp_p_count = len(fwd_sep_raw)

    if str(engine).lower() == "orca":
        from .orca_runner import eval_orca_dft_profile
        full_raw = eval_orca_dft_profile(full_raw, model, solvent, log, ts_index=ts_index)

    # Preserve TS saddle anchor without overwriting with raw np.argmax (essential for submerged barriers)
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

    e_r = full_raw[0]["energy_ev"]
    e_p = full_raw[-1]["energy_ev"]
    e_complex = min(full_raw[k]["energy_ev"] for k in range(ts_index + 1))

    frames = []
    for k, item in enumerate(full_raw):
        e_cur = item["energy_ev"]
        rel_kcal = (e_cur - e_r) * EV_TO_KCAL
        frames.append({
            "step": k - ts_index,
            "s_arc": round(s_coords[k], 4),
            "energy_ev": round(e_cur, 5),
            "rel_energy_kcal": round(rel_kcal, 3),
            "positions": item["positions"],
        })

    barrier_fwd = max(0.0, (e_ts - e_complex) * EV_TO_KCAL)
    barrier_rev = max(0.0, (e_ts - e_p) * EV_TO_KCAL)
    reaction_e = (e_p - e_r) * EV_TO_KCAL

    p_start = np.array(frames[0]["positions"])
    p_end = np.array(frames[-1]["positions"])
    sep_dist = float(np.linalg.norm(p_start - p_end))

    log(f"\n✓ Reaction Path Complete: {len(frames)} frames ({asymp_r_count} sep_R + {len(rev_irc_raw)} irc_R + TS + {len(fwd_irc_raw)} irc_P + {asymp_p_count} sep_P).")
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
    return {
        "reaction_smiles": f"{'.'.join(model.r_smiles)}>>{'.'.join(model.p_smiles)}",
        "symbols": model.symbols,
        "n_atoms": model.n_atoms,
        "bonds": bonds,
        "ts_frame_index": ts_index,
        "asymp_r_count": asymp_r_count,
        "asymp_p_count": asymp_p_count,
        "barriers": {
            "e_ts_ev": round(e_ts, 5),
            "e_reactant_ev": round(e_r, 5),
            "e_product_ev": round(e_p, 5),
            "forward_barrier_kcal": round(barrier_fwd, 2),
            "reverse_barrier_kcal": round(barrier_rev, 2),
            "reaction_energy_kcal": round(reaction_e, 2),
            "endpoint_separation_A": round(sep_dist, 2),
        },
        "fragments": build_fragments_meta(model),
        "frames": frames,
    }
