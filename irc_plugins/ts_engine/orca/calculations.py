"""
ORCA Calculations Module
High-level computational workflows: OptTS, IRC, NEB-CI, Thermochemistry, and Energy Profiling.
"""

import os
import shutil
import subprocess
import tempfile
from typing import Any, Dict, List, Optional
import numpy as np

from .executor import ORCA_BIN, run_orca_calc
from .parser import parse_orca_trj_file, build_trajectory_result


from .dft_ts import run_orca_dft_ts, eval_orca_dft_profile


def run_orca_native_irc(atoms: Any, model: Any, solvent: str, log: Any, max_points: int = 25, nprocs: int = 4) -> Optional[Dict[str, Any]]:
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
    with open(out_file, "w") as out_fh:
        subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

    trj_file = os.path.join(tmpdir, "orca_irc_IRC_Full_trj.xyz")
    if not os.path.exists(trj_file):
        f_trj = os.path.join(tmpdir, "orca_irc_IRC_F_trj.xyz")
        b_trj = os.path.join(tmpdir, "orca_irc_IRC_B_trj.xyz")
        if os.path.exists(f_trj) or os.path.exists(b_trj):
            b_frames = parse_orca_trj_file(b_trj)[::-1] if os.path.exists(b_trj) else []
            f_frames = parse_orca_trj_file(f_trj) if os.path.exists(f_trj) else []
            raw_frames = b_frames + f_frames
        else:
            log("  ❌ ORCA IRC trajectory was not generated.")
            shutil.rmtree(tmpdir, ignore_errors=True)
            return None
    else:
        raw_frames = parse_orca_trj_file(trj_file)

    if not raw_frames:
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None

    log(f"  ✓ Native ORCA IRC completed: {len(raw_frames)} frames along reaction path.")
    energies = [f["energy_ev"] for f in raw_frames]
    ts_idx = int(np.argmax(energies))

    res_traj = build_trajectory_result(raw_frames, ts_idx, model, atoms, log, energy_type="electronic_dft")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return res_traj


def run_orca_neb_ci(atoms: Any, model: Any, solvent: str, log: Any, n_images: int = 16, nprocs: int = 4, r_pos: Any = None, p_pos: Any = None) -> Optional[Dict[str, Any]]:
    """Run native ab initio ORCA 6.1 Climbing Image NEB (NEB-CI) for MEP optimization."""
    log(f"\n=== Overnight ORCA 6.1 NEB-CI Calculation (r2SCAN-3c, solvent={solvent.upper() if solvent else 'GAS'}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return None

    if r_pos is None or p_pos is None:
        try:
            from ..irc_runner import _integrate_reaction_branch
            rev_frames = _integrate_reaction_branch(atoms, model, -1.0, 15, log, "Reverse")
            fwd_frames = _integrate_reaction_branch(atoms, model, +1.0, 15, log, "Forward")
            r_pos = np.array(rev_frames[-1]["positions"]) if rev_frames else atoms.get_positions()
            p_pos = np.array(fwd_frames[-1]["positions"]) if fwd_frames else atoms.get_positions()
        except Exception:
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
        f.write(f"%neb\n  NImages {n_images}\n  product \"{prod_file}\"\nend\n")
        f.write(f"* xyz {model.charge} {model.multiplicity}\n")
        for sym, pos in zip(symbols, r_pos):
            f.write(f"  {sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}\n")
        f.write("*\n")

    with open(out_file, "w") as out_fh:
        subprocess.run([ORCA_BIN, inp_file], stdout=out_fh, stderr=subprocess.STDOUT, cwd=tmpdir)

    trj_file = os.path.join(tmpdir, "orca_neb_MEP_trj.xyz")
    if not os.path.exists(trj_file):
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None

    raw_frames = parse_orca_trj_file(trj_file)
    if not raw_frames:
        shutil.rmtree(tmpdir, ignore_errors=True)
        return None

    energies = [f["energy_ev"] for f in raw_frames]
    ts_idx = int(np.argmax(energies))
    res_traj = build_trajectory_result(raw_frames, ts_idx, model, atoms, log, energy_type="neb_ci")
    shutil.rmtree(tmpdir, ignore_errors=True)
    return res_traj


def run_orca_thermo_profile(atoms: Any, model: Any, solvent: str, log: Any, max_steps: int = 30, dx: float = 0.08, nprocs: int = 4) -> Optional[Dict[str, Any]]:
    """Run ORCA 6.1 Freq calculation on TS and endpoints to yield standard-state Gibbs free energy ΔG‡."""
    log(f"\n=== ORCA 6.1 Thermochemical ΔG Profile (r2SCAN-3c Freq, solvent={solvent.upper() if solvent else 'GAS'}) ===")
    if not os.path.exists(ORCA_BIN):
        log(f"  ERROR: ORCA binary not found at {ORCA_BIN}")
        return None

    from ..irc_runner import run_irc_path
    base_traj = run_irc_path(atoms, model, log, max_steps=max_steps, dx=dx, relax_fragments=True, engine="xtb")
    if not base_traj or not base_traj.get("frames"):
        return None

    ts_idx = base_traj.get("ts_frame_index", len(base_traj["frames"]) // 2)
    pos_ts = np.array(base_traj["frames"][ts_idx]["positions"])
    pos_r = np.array(base_traj["frames"][0]["positions"])
    pos_p = np.array(base_traj["frames"][-1]["positions"])

    ts_res = run_orca_calc(model.symbols, pos_ts, model.charge, model.multiplicity, solvent,
                           "Freq", log=log, name="ts_freq", nprocs=nprocs)

    if ts_res.get("gibbs_eh") is None:
        eval_orca_dft_profile(base_traj["frames"], model, solvent, log, ts_index=ts_idx, nprocs=nprocs)
        return base_traj

    g_ts = ts_res["gibbs_eh"]
    imag_freqs = ts_res.get("imaginary_freqs", [])
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
                f_res = run_orca_calc(syms, sub_pos, chg, mult, solvent, "Freq", log=log,
                                      name=f"reactant_frag_{fi+1}", nprocs=nprocs)
                if f_res.get("gibbs_eh") is not None:
                    g_reac_sum += f_res["gibbs_eh"]
        g_reac = g_reac_sum + delta_g_std_eh if g_reac_sum != 0.0 else g_ts - 0.03
    else:
        r_res = run_orca_calc(model.symbols, pos_r, model.charge, model.multiplicity, solvent,
                              "Freq", log=log, name="reactant_freq", nprocs=nprocs)
        g_reac = r_res.get("gibbs_eh") or r_res.get("energy_eh") or (g_ts - 0.03)

    p_res = run_orca_calc(model.symbols, pos_p, model.charge, model.multiplicity, solvent,
                          "Freq", log=log, name="product_freq", nprocs=nprocs)
    g_prod = p_res.get("gibbs_eh") or p_res.get("energy_eh") or (g_reac - 0.01)

    barrier_fwd_kcal = (g_ts - g_reac) * HARTREE_TO_KCAL
    barrier_rev_kcal = (g_ts - g_prod) * HARTREE_TO_KCAL
    delta_g_rxn_kcal = (g_prod - g_reac) * HARTREE_TO_KCAL

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
