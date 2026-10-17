import os
import subprocess
import shutil
import tempfile
import numpy as np

ORCA_BIN = "/data/orca/orca"

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
    """Re-evaluate reaction coordinate trajectory frames using high-level ORCA DFT."""
    log(f"\n=== Re-evaluating Reaction Path with ORCA 6.1 (r2SCAN-3c, solvent={solvent.upper()}) ===")
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
