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
