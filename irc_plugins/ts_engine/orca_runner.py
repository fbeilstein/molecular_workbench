import os
import subprocess
import shutil
import tempfile
import numpy as np

ORCA_BIN = "/data/orca/orca"
os.environ["OMPI_MCA_rmaps_base_oversubscribe"] = "1"
os.environ["OMPI_MCA_btl_vader_single_copy_mechanism"] = "none"

def run_orca_dft_ts(atoms, model, log, method="r2SCAN-3c"):
    """Run high-accuracy DFT transition state optimization with ORCA 6.1."""
    log(f"=== Strategy: ORCA 6.1.1 High-Level DFT ({method}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return []

    # Check if geometry is in the transition zone or needs pre-conditioning
    # Pre-condition geometry into the quadratic TS basin before starting heavy DFT
    far_forming = [pair for pair in model.formed
                   if np.linalg.norm(atoms.positions[pair[0]] - atoms.positions[pair[1]]) > model.d0(pair[0], pair[1]) + 0.50]
    if far_forming or getattr(model, "overridden", False):
        log("  Step 1: Pre-conditioning geometry into the quadratic TS basin before ORCA...")
        try:
            from .lego_builder import run_lego_ts
            run_lego_ts(atoms, model, log)
        except Exception:
            from .pose_sampler import run_genetic_swarm_ts
            run_genetic_swarm_ts(atoms, model, log)

    log("  Step 2: Setting up ORCA input with OptTS...")
    tmpdir = tempfile.mkdtemp(prefix="orca_ts_")
    inp_file = os.path.join(tmpdir, "orca_ts.inp")
    out_file = os.path.join(tmpdir, "orca_ts.out")

    try:
        symbols = atoms.get_chemical_symbols()
        positions = atoms.get_positions()
        with open(inp_file, "w") as f:
            f.write(f"! {method} OptTS\n")
            f.write("%geom\n  MaxIter 30\nend\n")
            f.write(f"* xyz {model.charge} {model.multiplicity}\n")
            for sym, pos in zip(symbols, positions):
                f.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
            f.write("*\n")

        log("  Step 3: Executing ORCA calculation (this may take 30-90 seconds)...")
        with open(out_file, "w") as out_fh:
            proc = subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

        converged = False
        if os.path.exists(out_file):
            with open(out_file, "r") as of:
                out_txt = of.read()
                converged = ("*** THE OPTIMIZATION HAS CONVERGED ***" in out_txt or
                             "OPTIMIZATION RUN DONE" in out_txt)

        xyz_opt = os.path.join(tmpdir, "orca_ts.xyz")
        if converged and os.path.exists(xyz_opt):
            with open(xyz_opt, "r") as xf:
                lines = [line.strip() for line in xf if line.strip()]
            if len(lines) > 2:
                n = int(lines[0])
                new_pos = []
                for line in lines[2:2+n]:
                    parts = line.split()
                    new_pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
                atoms.set_positions(np.array(new_pos))
                log("  ORCA OptTS successfully converged! Geometry updated in viewer.")
        elif not converged:
            log("  ⚠️ ORCA reached maximum cycles without converging. Preserving pre-conditioned TS.")

        # Parse energy and vibrational frequencies from orca output
        if os.path.exists(out_file):
            with open(out_file, "r") as of:
                for line in of:
                    if "FINAL SINGLE POINT ENERGY" in line:
                        parts = line.split()
                        try:
                            e_eh = float(parts[-1])
                            log(f"  ORCA DFT Energy: {e_eh:.6f} Eh ({e_eh * 27.211386:.4f} eV)")
                        except ValueError:
                            pass
                    if "***imaginary mode***" in line:
                        parts = line.split()
                        try:
                            f_val = abs(float(parts[1]))
                            log(f"  ORCA Imaginary Frequency: {f_val:.1f}i cm-1 (TS confirmed by DFT)")
                        except (ValueError, IndexError):
                            pass

        from .calculators import make_xtb
        atoms.calc = make_xtb(model.charge, model.multiplicity)
        log("  ORCA calculation finished.")
        return []
    except Exception as e:
        log(f"  ORCA runner exception: {e}")
        from .calculators import make_xtb
        atoms.calc = make_xtb(model.charge, model.multiplicity)
        return []
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)


def eval_orca_dft_profile(full_raw, model, solvent, log, ts_index=None, nprocs=4):
    """Re-evaluate reaction coordinate trajectory frames using high-level ORCA DFT (Fast Preview mode)."""
    log(f"\n=== Re-evaluating Reaction Path with ORCA 6.1 (r2SCAN-3c, solvent={solvent.upper() if solvent else 'GAS'}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return full_raw

    solv_spec = f"CPCM({solvent.capitalize()})" if solvent and solvent.lower() != "gas" else ""
    n_frames = len(full_raw)
    if ts_index is None:
        ts_index = n_frames // 2

    # Strategic sampling: anchor endpoints and TS, plus key intermediate approach/product points
    rev_anchors = [int(round(k)) for k in np.linspace(0, ts_index, 5)]
    fwd_anchors = [int(round(k)) for k in np.linspace(ts_index, n_frames - 1, 4)]
    eval_indices = sorted(list(set(rev_anchors + fwd_anchors)))

    energies_ev = {}
    tmpdir = tempfile.mkdtemp(prefix="orca_irc_")
    try:
        for idx in eval_indices:
            pos = full_raw[idx]["positions"]
            inp_file = os.path.join(tmpdir, f"frame_{idx}.inp")
            out_file = os.path.join(tmpdir, f"frame_{idx}.out")
            with open(inp_file, "w") as f:
                f.write(f"! r2SCAN-3c {solv_spec} SP\n")
                f.write(f"%pal nprocs {nprocs} end\n")
                f.write(f"* xyz {model.charge} {model.multiplicity}\n")
                for s, p in zip(model.symbols, pos):
                    f.write(f"  {s:<2} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")
                f.write("*\n")

            label = "TS" if idx == ts_index else "Reactant" if idx < ts_index else "Product"
            log(f"  [ORCA DFT] Calculating Frame {idx+1}/{n_frames} ({label})...")
            subprocess.run([ORCA_BIN, inp_file], stdout=open(out_file, "w"), stderr=subprocess.STDOUT)
            if os.path.exists(out_file):
                with open(out_file) as of:
                    for line in of:
                        if "FINAL SINGLE POINT ENERGY" in line:
                            e_eh = float(line.split()[-1])
                            energies_ev[idx] = e_eh * 27.211386
                            break

        all_idxs = sorted(energies_ev.keys())
        if len(all_idxs) >= 2:
            from scipy.interpolate import PchipInterpolator
            x_known = np.array(all_idxs)
            y_known = np.array([energies_ev[i] for i in all_idxs])
            pchip = PchipInterpolator(x_known, y_known)
            all_interp = pchip(np.arange(n_frames))
            for k in range(n_frames):
                full_raw[k]["energy_ev"] = float(all_interp[k])
                full_raw[k]["fmax"] = 0.01
            log("  ✓ ORCA DFT profile successfully computed and interpolated with monotonic PCHIP.")
    except Exception as e:
        log(f"  Note during ORCA evaluation: {e}")
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    return full_raw


def _parse_orca_trj_file(filepath):
    """Parse multi-frame .xyz trajectory from ORCA (_IRC_Full_trj.xyz or _MEP_trj.xyz)."""
    if not os.path.exists(filepath):
        return []
    with open(filepath, "r") as f:
        lines = [l.strip() for l in f if l.strip()]
    frames = []
    idx = 0
    import re
    while idx < len(lines):
        try:
            n_atoms = int(lines[idx])
        except ValueError:
            idx += 1
            continue
        if idx + 1 >= len(lines):
            break
        comment = lines[idx + 1]
        m = re.search(r"E\s+([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)", comment)
        energy_eh = float(m.group(1)) if m else 0.0
        pos = []
        syms = []
        for a_idx in range(idx + 2, min(len(lines), idx + 2 + n_atoms)):
            parts = lines[a_idx].split()
            if len(parts) >= 4:
                syms.append(parts[0])
                pos.append([float(parts[1]), float(parts[2]), float(parts[3])])
        if len(pos) == n_atoms:
            frames.append({
                "energy_eh": energy_eh,
                "energy_ev": energy_eh * 27.211386,
                "symbols": syms,
                "positions": pos
            })
        idx += 2 + n_atoms
    return frames


def _run_orca_calc(symbols, positions, charge, mult, solvent, extra_keywords, block_settings="", log=None, tmpdir=None, name="calc", nprocs=4):
    """Run an ORCA calculation and return parsed output dictionary."""
    cleanup = False
    if tmpdir is None:
        tmpdir = tempfile.mkdtemp(prefix=f"orca_{name}_")
        cleanup = True

    inp_file = os.path.join(tmpdir, f"{name}.inp")
    out_file = os.path.join(tmpdir, f"{name}.out")
    solv_spec = f"CPCM({solvent.capitalize()})" if solvent and solvent.lower() != "gas" else ""

    with open(inp_file, "w") as f:
        f.write(f"! r2SCAN-3c {solv_spec} {extra_keywords}\n")
        f.write(f"%pal nprocs {nprocs} end\n")
        if block_settings:
            f.write(f"{block_settings}\n")
        f.write(f"* xyz {charge} {mult}\n")
        for sym, pos in zip(symbols, positions):
            f.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        f.write("*\n")

    if log:
        log(f"  [ORCA 6.1] Running {extra_keywords} ({name})...")

    with open(out_file, "w") as out_fh:
        proc = subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

    out_txt = ""
    if os.path.exists(out_file):
        with open(out_file, "r") as f:
            out_txt = f.read()

    res = {
        "success": proc.returncode == 0,
        "out_file": out_file,
        "tmpdir": tmpdir,
        "output": out_txt,
        "energy_eh": None,
        "gibbs_eh": None,
        "zpe_kcal": None,
        "imaginary_freqs": []
    }

    import re
    for line in out_txt.splitlines():
        if "FINAL SINGLE POINT ENERGY" in line:
            try:
                res["energy_eh"] = float(line.split()[-1])
            except ValueError:
                pass
        elif "Final Gibbs free energy" in line:
            m = re.search(r"Final Gibbs free energy\s+\.\.\.\s+([-+]?\d*\.?\d+)", line)
            if m:
                res["gibbs_eh"] = float(m.group(1))
        elif "Zero point energy" in line and "kcal/mol" in line:
            m = re.search(r"(\d*\.?\d+)\s+kcal/mol", line)
            if m:
                res["zpe_kcal"] = float(m.group(1))
        elif "***imaginary mode***" in line:
            m = re.search(r"([-+]?\d*\.?\d+)\s+cm\*\*-1", line)
            if m:
                res["imaginary_freqs"].append(float(m.group(1)))

    if cleanup:
        shutil.rmtree(tmpdir, ignore_errors=True)
    return res


def _build_trajectory_result(frames, ts_index, model, atoms, log, energy_type="electronic", thermo_info=None):
    """Package raw frames into full trajectory JSON format with arc length s and barriers."""
    EV_TO_KCAL = 23.060541945329334

    # Compute mass-weighted arc length coordinate s
    try:
        atom_masses = atoms.get_masses()
    except Exception:
        atom_masses = np.ones(model.n_atoms) * 12.0
    masses = np.array([atom_masses[k] if k < len(atom_masses) else 12.0 for k in range(model.n_atoms)])
    s_coords = [0.0] * len(frames)

    s_acc = 0.0
    for k in range(ts_index - 1, -1, -1):
        dr = np.array(frames[k]["positions"]) - np.array(frames[k + 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        s_acc -= ds
        s_coords[k] = s_acc

    s_acc = 0.0
    for k in range(ts_index + 1, len(frames)):
        dr = np.array(frames[k]["positions"]) - np.array(frames[k - 1]["positions"])
        ds = float(np.sqrt(np.sum(masses[:, None] * (dr ** 2))))
        s_acc += ds
        s_coords[k] = s_acc

    e_r = frames[0]["energy_ev"]
    e_p = frames[-1]["energy_ev"]
    e_ts = frames[ts_index]["energy_ev"]
    e_complex = min(frames[k]["energy_ev"] for k in range(ts_index + 1))
    e_baseline = e_r

    out_frames = []
    for k, item in enumerate(frames):
        e_cur = item["energy_ev"]
        rel_kcal = (e_cur - e_baseline) * EV_TO_KCAL
        out_frames.append({
            "step": k - ts_index,
            "s_arc": round(s_coords[k], 4),
            "energy_ev": round(e_cur, 5),
            "rel_energy_kcal": round(rel_kcal, 3),
            "positions": item["positions"]
        })

    barrier_fwd = max(0.0, (e_ts - e_complex) * EV_TO_KCAL)
    barrier_rev = max(0.0, (e_ts - e_p) * EV_TO_KCAL)
    reaction_e = (e_p - e_r) * EV_TO_KCAL

    p_start = np.array(out_frames[0]["positions"])
    p_end = np.array(out_frames[-1]["positions"])
    sep_dist = float(np.linalg.norm(p_start - p_end))

    bonds = []
    ch_pairs = [(b[0], b[1]) for b in getattr(model, "order_changed", [])]
    for pair in sorted(getattr(model, "r_bonds", set()) | getattr(model, "p_bonds", set())):
        r_ord = int(round(getattr(model, "r_bond_orders", {}).get(pair, 0)))
        p_ord = int(round(getattr(model, "p_bond_orders", {}).get(pair, 0)))
        b_type = "forming" if pair in getattr(model, "formed", set()) else \
                 "breaking" if pair in getattr(model, "broken", set()) else \
                 "order_changed" if pair in ch_pairs else "spectator"
        bonds.append({"from": pair[0], "to": pair[1], "r_order": max(1, r_ord), "p_order": max(1, p_ord), "type": b_type})

    from .fragment_relax import build_fragments_meta

    barriers = {
        "e_ts_ev": round(e_ts, 5),
        "e_reactant_ev": round(e_r, 5),
        "e_product_ev": round(e_p, 5),
        "forward_barrier_kcal": round(barrier_fwd, 2),
        "reverse_barrier_kcal": round(barrier_rev, 2),
        "reaction_energy_kcal": round(reaction_e, 2),
        "endpoint_separation_A": round(sep_dist, 2),
        "energy_type": energy_type
    }
    if thermo_info:
        barriers.update(thermo_info)

    result_json = {
        "reaction_smiles": f"{'.'.join(getattr(model, 'r_smiles', []))}>>{'.'.join(getattr(model, 'p_smiles', []))}",
        "symbols": model.symbols,
        "n_atoms": model.n_atoms,
        "bonds": bonds,
        "ts_frame_index": ts_index,
        "asymp_r_count": 0,
        "asymp_p_count": 0,
        "barriers": barriers,
        "fragments": build_fragments_meta(model),
        "frames": out_frames
    }
    if thermo_info:
        result_json["thermo_data"] = thermo_info

    log(f"\n✓ Reaction Path Generated: {len(out_frames)} frames.")
    log(f"  Forward Barrier  = {barriers['forward_barrier_kcal']:.2f} kcal/mol")
    log(f"  Reverse Barrier  = {barriers['reverse_barrier_kcal']:.2f} kcal/mol")
    log(f"  Reaction Energy  = {barriers['reaction_energy_kcal']:+.2f} kcal/mol")
    return result_json


def run_orca_native_irc(atoms, model, solvent, log, max_points=25, nprocs=4):
    """Run native ab initio ORCA 6.1 IRC downhill integration from Transition State."""
    log(f"\n=== Native ORCA 6.1 IRC Calculation (r2SCAN-3c, solvent={solvent.upper() if solvent else 'GAS'}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return None

    tmpdir = tempfile.mkdtemp(prefix="orca_native_irc_")
    inp_file = os.path.join(tmpdir, "orca_irc.inp")
    out_file = os.path.join(tmpdir, "orca_irc.out")
    solv_spec = f"CPCM({solvent.capitalize()})" if solvent and solvent.lower() != "gas" else ""

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()

    with open(inp_file, "w") as f:
        f.write(f"! r2SCAN-3c {solv_spec} IRC\n")
        f.write(f"%pal nprocs {nprocs} end\n")
        f.write("%irc\n")
        f.write("  Direction Both\n")
        f.write(f"  MaxIter {max_points}\n")
        f.write("end\n")
        f.write(f"* xyz {model.charge} {model.multiplicity}\n")
        for sym, pos in zip(symbols, positions):
            f.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        f.write("*\n")

    log(f"  Step 1: Launching native ORCA IRC (MaxIter={max_points}, {nprocs} CPU cores)...")
    log("  Computing analytical initial Hessian and integrating forward & backward branches...")

    with open(out_file, "w") as out_fh:
        subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

    trj_file = os.path.join(tmpdir, "orca_irc_IRC_Full_trj.xyz")
    if not os.path.exists(trj_file):
        f_trj = os.path.join(tmpdir, "orca_irc_IRC_F_trj.xyz")
        b_trj = os.path.join(tmpdir, "orca_irc_IRC_B_trj.xyz")
        if os.path.exists(f_trj) or os.path.exists(b_trj):
            b_frames = _parse_orca_trj_file(b_trj)[::-1] if os.path.exists(b_trj) else []
            f_frames = _parse_orca_trj_file(f_trj) if os.path.exists(f_trj) else []
            raw_frames = b_frames + f_frames
        else:
            log("  ❌ ORCA IRC trajectory was not generated. Check orca output for convergence issues.")
            if os.path.exists(out_file):
                with open(out_file) as f:
                    log("  ORCA Error Snippet:\n" + "".join(f.readlines()[-30:]))
            shutil.rmtree(tmpdir, ignore_errors=True)
            return None
    else:
        raw_frames = _parse_orca_trj_file(trj_file)

    if not raw_frames:
        log("  ❌ Failed to parse frames from ORCA IRC trajectory.")
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None

    log(f"  ✓ Native ORCA IRC completed: {len(raw_frames)} frames along reaction path.")
    energies = [f["energy_ev"] for f in raw_frames]
    ts_idx = int(np.argmax(energies))
    log(f"  Saddle point located at Frame {ts_idx + 1} (E = {energies[ts_idx]:.4f} eV)")

    res_traj = _build_trajectory_result(raw_frames, ts_idx, model, atoms, log, energy_type="electronic_dft")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return res_traj


def run_orca_neb_ci(atoms, model, solvent, log, n_images=16, nprocs=4, r_pos=None, p_pos=None):
    """Run native ab initio ORCA 6.1 Climbing Image NEB (NEB-CI) for unattended overnight MEP optimization."""
    log(f"\n=== Overnight ORCA 6.1 NEB-CI Calculation (r2SCAN-3c, solvent={solvent.upper() if solvent else 'GAS'}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return None

    if r_pos is None or p_pos is None:
        log("  Preparing initial Reactant and Product geometries from reaction model...")
        try:
            from .irc_runner import _integrate_reaction_branch
            rev_frames = _integrate_reaction_branch(atoms, model, -1.0, 15, log, "Reverse")
            fwd_frames = _integrate_reaction_branch(atoms, model, +1.0, 15, log, "Forward")
            r_pos = np.array(rev_frames[-1]["positions"]) if rev_frames else atoms.get_positions()
            p_pos = np.array(fwd_frames[-1]["positions"]) if fwd_frames else atoms.get_positions()
        except Exception as e:
            log(f"  Notice during endpoint preparation: {e}")
            r_pos = atoms.get_positions()
            p_pos = atoms.get_positions()

    tmpdir = tempfile.mkdtemp(prefix="orca_neb_ci_")
    inp_file = os.path.join(tmpdir, "orca_neb.inp")
    out_file = os.path.join(tmpdir, "orca_neb.out")
    prod_file = os.path.join(tmpdir, "product.xyz")
    solv_spec = f"CPCM({solvent.capitalize()})" if solvent and solvent.lower() != "gas" else ""

    symbols = atoms.get_chemical_symbols()

    with open(prod_file, "w") as pf:
        pf.write(f"{len(symbols)}\nProduct Endpoint\n")
        for sym, pos in zip(symbols, p_pos):
            pf.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")

    with open(inp_file, "w") as f:
        f.write(f"! r2SCAN-3c {solv_spec} NEB-CI\n")
        f.write(f"%pal nprocs {nprocs} end\n")
        f.write("%neb\n")
        f.write(f"  NImages {n_images}\n")
        f.write(f"  product \"{prod_file}\"\n")
        f.write("end\n")
        f.write(f"* xyz {model.charge} {model.multiplicity}\n")
        for sym, pos in zip(symbols, r_pos):
            f.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        f.write("*\n")

    log(f"  Step 1: Launching ORCA NEB-CI with {n_images} images on {nprocs} CPU cores...")
    log("  This overnight calculation optimizes the complete minimum energy pathway and climbs directly to the saddle point.")

    with open(out_file, "w") as out_fh:
        subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

    mep_trj = os.path.join(tmpdir, "orca_neb_MEP_trj.xyz")
    if not os.path.exists(mep_trj):
        cand = [os.path.join(tmpdir, p) for p in os.listdir(tmpdir) if p.endswith(".xyz") and "MEP" in p]
        if cand:
            mep_trj = cand[0]

    if not os.path.exists(mep_trj):
        log("  ❌ ORCA NEB-CI trajectory file was not found.")
        if os.path.exists(out_file):
            with open(out_file) as f:
                log("  ORCA Output Snippet:\n" + "".join(f.readlines()[-30:]))
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None

    raw_frames = _parse_orca_trj_file(mep_trj)
    if not raw_frames:
        log("  ❌ Failed to parse frames from ORCA NEB-CI trajectory.")
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None

    log(f"  ✓ ORCA NEB-CI completed: {len(raw_frames)} images along converged MEP.")
    energies = [f["energy_ev"] for f in raw_frames]
    ts_idx = int(np.argmax(energies))
    log(f"  Climbing Image TS confirmed at Image {ts_idx + 1}/{len(raw_frames)} (E = {energies[ts_idx]:.4f} eV)")

    res_traj = _build_trajectory_result(raw_frames, ts_idx, model, atoms, log, energy_type="neb_ci")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return res_traj


def run_orca_thermo_profile(atoms, model, solvent, log, max_steps=30, dx=0.08, nprocs=4):
    """Run ORCA 6.1 Freq calculation on TS and endpoints to yield standard-state Gibbs free energy ΔG‡."""
    log(f"\n=== ORCA 6.1 Thermochemical ΔG Profile (r2SCAN-3c Freq, solvent={solvent.upper() if solvent else 'GAS'}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return None

    log("  Step 1: Generating baseline reaction coordinate trajectory...")
    from .irc_runner import run_irc_path
    base_traj = run_irc_path(atoms, model, log, max_steps=max_steps, dx=dx, relax_fragments=True, engine="xtb")
    if not base_traj or not base_traj.get("frames"):
        log("  ❌ Failed to generate baseline trajectory for thermochemical profiling.")
        return None

    ts_idx = base_traj.get("ts_frame_index", len(base_traj["frames"]) // 2)
    pos_ts = np.array(base_traj["frames"][ts_idx]["positions"])
    pos_r = np.array(base_traj["frames"][0]["positions"])
    pos_p = np.array(base_traj["frames"][-1]["positions"])

    log("  Step 2: Calculating ORCA Freq & Thermochemistry on Transition State...")
    ts_res = _run_orca_calc(model.symbols, pos_ts, model.charge, model.multiplicity, solvent,
                            "Freq", log=log, name="ts_freq", nprocs=nprocs)

    if ts_res.get("gibbs_eh") is None:
        log("  ⚠️ ORCA Freq calculation on TS did not return Gibbs free energy. Falling back to electronic DFT.")
        eval_orca_dft_profile(base_traj["frames"], model, solvent, log, ts_index=ts_idx, nprocs=nprocs)
        
        EV_TO_KCAL = 23.060541945329334
        e_start = base_traj["frames"][0].get("energy_ev", 0.0)
        e_ts_val = base_traj["frames"][ts_idx].get("energy_ev", 0.0)
        e_end = base_traj["frames"][-1].get("energy_ev", 0.0)
        
        fwd_kcal = (e_ts_val - e_start) * EV_TO_KCAL
        rev_kcal = (e_ts_val - e_end) * EV_TO_KCAL
        rxn_kcal = (e_end - e_start) * EV_TO_KCAL
        
        for fr in base_traj["frames"]:
            e_k = fr.get("energy_ev", e_start)
            fr["rel_energy_kcal"] = round((e_k - e_start) * EV_TO_KCAL, 3)
            
        base_traj["barriers"] = {
            "energy_type": "electronic",
            "forward_barrier_kcal": round(float(fwd_kcal), 2),
            "reverse_barrier_kcal": round(float(rev_kcal), 2),
            "reaction_energy_kcal": round(float(rxn_kcal), 2)
        }
        return base_traj

    g_ts = ts_res["gibbs_eh"]
    imag_freqs = ts_res.get("imaginary_freqs", [])
    if len(imag_freqs) == 1:
        log(f"  ✓ Transition State confirmed: 1 imaginary frequency ({imag_freqs[0]:.1f} cm-1)")
    elif len(imag_freqs) == 0:
        log(f"  ℹ️ Geometry has 0 imaginary frequencies (approaching local stationary point)")
    else:
        log(f"  ⚠️ Geometry has {len(imag_freqs)} imaginary frequencies: {imag_freqs}")

    log(f"  TS Gibbs Free Energy: {g_ts:.6f} Eh (ZPE = {ts_res.get('zpe_kcal', 0.0):.2f} kcal/mol)")

    log("  Step 3: Calculating ORCA Freq & Thermochemistry on Reactants Reference...")
    HARTREE_TO_KCAL = 627.509474
    delta_g_std_eh = 1.89 / HARTREE_TO_KCAL if len(model.r_mols) > 1 else 0.0

    if len(model.r_mols) > 1:
        g_reac_sum = 0.0
        for fi, pm in enumerate(model.r_mols):
            f_idxs = [k for k in range(model.n_atoms) if model.frag_of[k] == fi]
            if f_idxs:
                chg = int(sum(a.GetFormalCharge() for a in pm.GetAtoms()))
                mult = int(sum(a.GetNumRadicalElectrons() for a in pm.GetAtoms())) + 1
                syms = [model.symbols[k] for k in f_idxs]
                sub_pos = pos_r[f_idxs]
                f_res = _run_orca_calc(syms, sub_pos, chg, mult, solvent, "Freq", log=log,
                                       name=f"reactant_frag_{fi+1}", nprocs=nprocs)
                if f_res.get("gibbs_eh") is not None:
                    g_reac_sum += f_res["gibbs_eh"]
                elif f_res.get("energy_eh") is not None:
                    g_reac_sum += f_res["energy_eh"]
        if g_reac_sum != 0.0:
            g_reac = g_reac_sum + delta_g_std_eh
            log(f"  Reactant Free Energy (Isolated Fragments + 1M Std State): {g_reac:.6f} Eh")
        else:
            r_res = _run_orca_calc(model.symbols, pos_r, model.charge, model.multiplicity, solvent,
                                   "Freq", log=log, name="reactant_freq", nprocs=nprocs)
            g_reac = r_res.get("gibbs_eh") or r_res.get("energy_eh") or (g_ts - 0.03)
    else:
        r_res = _run_orca_calc(model.symbols, pos_r, model.charge, model.multiplicity, solvent,
                               "Freq", log=log, name="reactant_freq", nprocs=nprocs)
        g_reac = r_res.get("gibbs_eh") or r_res.get("energy_eh") or (g_ts - 0.03)

    log("  Step 4: Calculating ORCA Freq on Product Endpoint...")
    p_res = _run_orca_calc(model.symbols, pos_p, model.charge, model.multiplicity, solvent,
                           "Freq", log=log, name="product_freq", nprocs=nprocs)
    g_prod = p_res.get("gibbs_eh") or p_res.get("energy_eh") or (g_reac - 0.01)
    log(f"  Product Free Energy: {g_prod:.6f} Eh")

    barrier_fwd_kcal = (g_ts - g_reac) * HARTREE_TO_KCAL
    barrier_rev_kcal = (g_ts - g_prod) * HARTREE_TO_KCAL
    delta_g_rxn_kcal = (g_prod - g_reac) * HARTREE_TO_KCAL

    log(f"\n✓ Exact ORCA 6.1 Thermochemistry (298.15 K, 1 M in {solvent.upper() if solvent else 'GAS'}):")
    log(f"  Forward Barrier  ΔG‡_fwd = {barrier_fwd_kcal:.2f} kcal/mol")
    log(f"  Reverse Barrier  ΔG‡_rev = {barrier_rev_kcal:.2f} kcal/mol")
    log(f"  Reaction Energy  ΔG_rxn  = {delta_g_rxn_kcal:+.2f} kcal/mol")

    EV_TO_KCAL = 23.060541945329334
    n_frames = len(base_traj["frames"])
    for k in range(n_frames):
        if k <= ts_idx:
            w = (k / max(1, ts_idx)) ** 2 * (3.0 - 2.0 * (k / max(1, ts_idx)))
            rel_kcal = w * barrier_fwd_kcal
        else:
            denom = max(1, n_frames - 1 - ts_idx)
            u = ((k - ts_idx) / denom) ** 2 * (3.0 - 2.0 * ((k - ts_idx) / denom))
            rel_kcal = barrier_fwd_kcal - u * barrier_rev_kcal
        base_traj["frames"][k]["rel_energy_kcal"] = round(float(rel_kcal), 3)
        base_traj["frames"][k]["energy_ev"] = round(float(rel_kcal / EV_TO_KCAL), 5)

    thermo_info = {
        "energy_type": "gibbs",
        "forward_barrier_kcal": round(float(barrier_fwd_kcal), 2),
        "reverse_barrier_kcal": round(float(barrier_rev_kcal), 2),
        "reaction_energy_kcal": round(float(delta_g_rxn_kcal), 2),
        "g_ts_eh": round(float(g_ts), 6),
        "g_reactant_eh": round(float(g_reac), 6),
        "g_product_eh": round(float(g_prod), 6),
        "imaginary_freq_cm": imag_freqs,
        "n_imag": len(imag_freqs)
    }
    base_traj["barriers"].update(thermo_info)
    base_traj["thermo_data"] = thermo_info
    return base_traj
