#!/usr/bin/env python3
"""
rxn_path_xtb — xTB-based reaction path finder (Tier 1)

Uses xTB's built-in RMSD-push/pull path finder to generate
reaction trajectories from matched reactant/product complex XYZ files.

Pipeline:
  1. Pre-optimize both endpoints at GFN2-xTB level (critical!)
  2. Run xTB --path to find reaction path + TS guess
  3. Parse trajectory + energies into unified output format

Usage:
    python rxn_path_xtb.py reactant.xyz product.xyz --charge -1 -o output/
"""

import argparse
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile

# ── Configuration ────────────────────────────────────────────────────────────

XTB_BIN = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'xtb', 'bin', 'xtb')
XTB_PARAM = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         'xtb', 'share', 'xtb')

# Default path finder parameters (from xTB docs, good starting point)
DEFAULT_PATH_PARAMS = {
    'nrun': 1,
    'npoint': 25,
    'anopt': 10,
    'kpush': 0.003,
    'kpull': -0.015,
    'ppull': 0.05,
    'alp': 1.2,
}


# ── Helpers ──────────────────────────────────────────────────────────────────

def _read_xyz(path):
    """Read XYZ file → list of (symbol, x, y, z)."""
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def _write_xyz(atoms, path, comment=''):
    """Write list of (sym, x, y, z) to XYZ file."""
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for sym, x, y, z in atoms:
            f.write(f"{sym:2s} {x:14.8f} {y:14.8f} {z:14.8f}\n")


def _run_xtb(args, cwd, timeout=600):
    """
    Run xTB with given arguments. Returns (returncode, stdout).
    Sets XTBPATH for parameter files.
    """
    env = os.environ.copy()
    env['XTBPATH'] = XTB_PARAM
    env['OMP_NUM_THREADS'] = '4'
    env['MKL_NUM_THREADS'] = '4'

    cmd = [XTB_BIN] + args
    try:
        proc = subprocess.run(
            cmd, cwd=cwd, capture_output=True, text=True,
            timeout=timeout, env=env
        )
        return proc.returncode, proc.stdout + proc.stderr
    except subprocess.TimeoutExpired:
        return -1, "TIMEOUT: xTB exceeded time limit"


# ── Step 1: Pre-optimize endpoints ──────────────────────────────────────────

# Covalent radii (Å) for fragment detection
_COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84,
    'C': 0.76, 'N': 0.71, 'O': 0.66, 'F': 0.57, 'Ne': 0.58,
    'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07,
    'S': 1.05, 'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76,
    'Br': 1.20, 'I': 1.39,
}


def _detect_fragments(atoms, scale=1.3):
    """
    Detect molecular fragments in a complex using bonded connectivity.
    Returns list of sets of atom indices, one set per fragment.
    """
    import math
    n = len(atoms)
    # Build adjacency
    adj = [set() for _ in range(n)]
    for i in range(n):
        for j in range(i + 1, n):
            si, xi, yi, zi = atoms[i]
            sj, xj, yj, zj = atoms[j]
            ri = _COVALENT_RADII.get(si, 1.0)
            rj = _COVALENT_RADII.get(sj, 1.0)
            dist = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
            if dist < scale * (ri + rj):
                adj[i].add(j)
                adj[j].add(i)

    # BFS to find connected components
    visited = [False] * n
    fragments = []
    for start in range(n):
        if visited[start]:
            continue
        frag = set()
        queue = [start]
        while queue:
            node = queue.pop(0)
            if visited[node]:
                continue
            visited[node] = True
            frag.add(node)
            for nb in adj[node]:
                if not visited[nb]:
                    queue.append(nb)
        fragments.append(frag)
    return fragments


def _find_closest_interfragment_pairs(atoms, fragments, n_pairs=3):
    """
    Find the N closest atom pairs between different fragments.
    Returns list of (atom_i, atom_j, distance) tuples with 1-based indices.
    """
    import math
    pairs = []
    for fi, frag_i in enumerate(fragments):
        for fj, frag_j in enumerate(fragments):
            if fj <= fi:
                continue
            for ai in frag_i:
                for aj in frag_j:
                    si, xi, yi, zi = atoms[ai]
                    sj, xj, yj, zj = atoms[aj]
                    dist = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
                    pairs.append((ai + 1, aj + 1, dist))  # 1-based for xTB
    pairs.sort(key=lambda x: x[2])
    return pairs[:n_pairs]


def optimize_endpoint(xyz_path, charge=0, solvent='water', cwd=None,
                      keep_fragments_separated=False):
    """
    Optimize a molecular geometry with GFN2-xTB.

    If keep_fragments_separated=True, detects molecular fragments and adds
    interfragment distance constraints to prevent recombination of ion pairs
    (e.g. E1 products: carbocation + leaving group). Internal geometry of
    each fragment relaxes freely, so carbocations naturally adopt sp2 planar.

    Returns path to optimized XYZ file, or original if optimization fails.
    """
    if cwd is None:
        cwd = os.path.dirname(xyz_path)

    basename = os.path.splitext(os.path.basename(xyz_path))[0]
    work_dir = os.path.join(cwd, f'_opt_{basename}')
    os.makedirs(work_dir, exist_ok=True)
    input_xyz = os.path.join(work_dir, 'input.xyz')
    shutil.copy2(xyz_path, input_xyz)

    args = ['input.xyz', '--opt', 'tight', '--gfn', '2',
            '--chrg', str(charge)]
    if solvent:
        args += ['--alpb', solvent]

    # For dissociation products (ion pairs), constrain interfragment distances
    # so oppositely-charged fragments don't recombine during optimization
    if keep_fragments_separated:
        atoms = _read_xyz(input_xyz)
        fragments = _detect_fragments(atoms, scale=1.15)
        if len(fragments) > 1:
            pairs = _find_closest_interfragment_pairs(atoms, fragments,
                                                      n_pairs=len(fragments))
            if pairs:
                inp_path = os.path.join(work_dir, 'constrain.inp')
                with open(inp_path, 'w') as f:
                    f.write("$constrain\n")
                    f.write("  force constant=0.5\n")
                    for ai, aj, dist in pairs:
                        # Pin at current distance (or slightly beyond) to
                        # prevent collapse while allowing internal relaxation
                        pin_dist = max(dist, 3.5)
                        f.write(f"  distance: {ai}, {aj}, {pin_dist:.4f}\n")
                        print(f"  Constraining atoms {ai}–{aj} at {pin_dist:.2f} Å "
                              f"(current {dist:.2f} Å)")
                    f.write("$end\n")
                args += ['--input', 'constrain.inp']
                print(f"  {len(fragments)} fragments detected — "
                      f"constraining {len(pairs)} interfragment distances")

    print(f"  Optimizing {basename} at GFN2-xTB...")
    rc, output = _run_xtb(args, work_dir, timeout=300)

    opt_xyz = os.path.join(work_dir, 'xtbopt.xyz')
    if rc == 0 and os.path.exists(opt_xyz):
        print(f"  ✓ Optimized {basename}")
        return opt_xyz

    # Try looser convergence
    args_loose = list(args)
    try:
        idx = args_loose.index('tight')
        args_loose[idx] = 'normal'
    except ValueError:
        pass
    rc2, _ = _run_xtb(args_loose, work_dir, timeout=300)
    if rc2 == 0 and os.path.exists(opt_xyz):
        print(f"  ✓ Optimized {basename} (normal convergence)")
        return opt_xyz

    print(f"  ⚠ Optimization failed for {basename}, using original")
    return xyz_path


# ── Step 2: Run xTB path finder ─────────────────────────────────────────────

def run_xtb_path(reactant_xyz, product_xyz, charge=0, solvent='water',
                 work_dir=None, path_params=None):
    """
    Run the xTB RMSD-push/pull path finder.
    Automatically retries with escalating push/pull parameters if first attempt
    finds no barrier.

    Args:
        reactant_xyz: path to optimized reactant complex XYZ
        product_xyz: path to optimized product complex XYZ
        charge: total system charge
        solvent: solvent name for ALPB (None for gas phase)
        work_dir: directory for calculation files
        path_params: dict overriding DEFAULT_PATH_PARAMS

    Returns:
        dict with keys: success, trajectory_xyz, ts_xyz, energies, log
    """
    import math

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix='xtb_path_')
    os.makedirs(work_dir, exist_ok=True)

    # Calculate RMSD between endpoints for adaptive parameter tuning
    r_atoms = _read_xyz(reactant_xyz)
    p_atoms = _read_xyz(product_xyz)
    rmsd = 0.0
    for (_, rx, ry, rz), (_, px, py, pz) in zip(r_atoms, p_atoms):
        rmsd += (rx - px)**2 + (ry - py)**2 + (rz - pz)**2
    rmsd = math.sqrt(rmsd / len(r_atoms))

    # Define parameter sets to try (escalating aggressiveness)
    param_sets = []
    if path_params:
        # User-specified params first
        p = DEFAULT_PATH_PARAMS.copy()
        p.update(path_params)
        param_sets.append(('user', p))

    # Auto-tuned based on RMSD
    if rmsd < 1.0:
        # Small structural change: gentle push/pull
        param_sets.append(('gentle', {
            'nrun': 1, 'npoint': 25, 'anopt': 10,
            'kpush': 0.003, 'kpull': -0.015, 'ppull': 0.05, 'alp': 1.2,
        }))
    else:
        # Larger structural change: stronger push/pull
        param_sets.append(('moderate', {
            'nrun': 1, 'npoint': 30, 'anopt': 10,
            'kpush': 0.005, 'kpull': -0.025, 'ppull': 0.05, 'alp': 1.0,
        }))

    # Aggressive fallback
    param_sets.append(('aggressive', {
        'nrun': 1, 'npoint': 40, 'anopt': 15,
        'kpush': 0.010, 'kpull': -0.050, 'ppull': 0.10, 'alp': 0.8,
    }))

    # Copy structures into working directory
    start_xyz = os.path.join(work_dir, 'start.xyz')
    end_xyz = os.path.join(work_dir, 'end.xyz')
    shutil.copy2(reactant_xyz, start_xyz)
    shutil.copy2(product_xyz, end_xyz)

    print(f"  R/P RMSD = {rmsd:.3f} Å")

    best_result = None

    for attempt_name, params in param_sets:
        # Write path input file
        inp_path = os.path.join(work_dir, 'path.inp')
        with open(inp_path, 'w') as f:
            f.write("$path\n")
            for k, v in params.items():
                f.write(f"  {k}={v}\n")
            f.write("$end\n")

        # Run xTB path
        args = ['start.xyz', '--path', 'end.xyz', '--input', 'path.inp',
                '--chrg', str(charge), '--gfn', '2']
        if solvent:
            args += ['--alpb', solvent]

        print(f"  Running xTB path ({attempt_name}, kpush={params['kpush']})...")
        rc, output = _run_xtb(args, work_dir, timeout=600)

        # Save log
        log_path = os.path.join(work_dir, 'xtb_path.log')
        with open(log_path, 'w') as f:
            f.write(output)

        # Parse results
        result = {
            'success': False,
            'trajectory_xyz': None,
            'ts_xyz': None,
            'energies': None,
            'barrier_forward_kcal': None,
            'barrier_reverse_kcal': None,
            'reaction_energy_kcal': None,
            'method': 'xtb_path',
            'log': log_path,
            'work_dir': work_dir,
        }

        if rc != 0:
            print(f"  ✗ xTB path failed ({attempt_name}, exit code {rc})")
            best_result = best_result or result
            continue

        # Check for trajectory files
        # The final optimized path is always written to xtbpath.xyz
        final_path = os.path.join(work_dir, 'xtbpath.xyz')
        if os.path.exists(final_path) and os.path.getsize(final_path) > 0:
            result['trajectory_xyz'] = final_path
        else:
            # Fallback if xtbpath.xyz is missing
            m = re.search(r'path\s+(\d+)\s+taken', output)
            if m:
                chosen_idx = int(m.group(1)) - 1
                chosen_path = os.path.join(work_dir, f'xtbpath_{chosen_idx}.xyz')
                if os.path.exists(chosen_path):
                    result['trajectory_xyz'] = chosen_path

        # TS geometry
        ts_path = os.path.join(work_dir, 'xtbpath_ts.xyz')
        if os.path.exists(ts_path):
            result['ts_xyz'] = ts_path

        # Parse energy data from output
        energies = _parse_path_energies(output)
        if energies:
            result['energies'] = energies

        # Parse barrier info
        fw_match = re.search(r'forward\s+barrier\s*\(kcal\)\s*:\s*([\d.]+)', output)
        bw_match = re.search(r'backward\s+barrier\s*\(kcal\)\s*:\s*([\d.]+)', output)
        re_match = re.search(r'reaction\s+energy\s*\(kcal\)\s*:\s*([-\d.]+)', output)

        if fw_match:
            result['barrier_forward_kcal'] = float(fw_match.group(1))
        if bw_match:
            result['barrier_reverse_kcal'] = float(bw_match.group(1))
        if re_match:
            result['reaction_energy_kcal'] = float(re_match.group(1))

        # Check if we got a meaningful barrier
        fw = result.get('barrier_forward_kcal') or 0
        if result['trajectory_xyz'] and fw > 0.5:
            result['success'] = True
            print(f"  ✓ Path found: barrier = {result['barrier_forward_kcal']:.1f} kcal/mol")
            return result  # Success — no need to retry

        # Zero or low barrier — save as best so far but try next params
        if result['trajectory_xyz']:
            print(f"  ⚠ Barrier ≈ 0, trying stronger parameters...")
            best_result = result
        else:
            best_result = best_result or result

    # Exhausted all parameter sets
    if best_result and best_result.get('trajectory_xyz'):
        best_result['success'] = True
        fw = best_result.get('barrier_forward_kcal')
        if fw is not None and fw > 0.1:
            print(f"  ✓ Path found: barrier = {fw:.1f} kcal/mol")
        else:
            print(f"  ⚠ Path found but barrier ≈ 0 (may be barrierless at GFN2-xTB level)")
        return best_result

    return best_result or result


def _parse_path_energies(output):
    """
    Parse the path data table from xTB output.
    Returns list of dicts: [{point, drms, energy_kcal, pmode_ovlp, pmode_grad}, ...]
    """
    energies = []

    # Look for the "path data" section
    # Format: point  drms  energy  pmode ovlp  pmode grad
    in_path_data = False
    for line in output.split('\n'):
        if 'path data' in line.lower():
            in_path_data = True
            continue
        if in_path_data and line.strip().startswith('point'):
            continue  # header line
        if in_path_data and line.strip():
            parts = line.split()
            if len(parts) >= 3:
                try:
                    point = int(parts[0])
                    drms = float(parts[1])
                    energy = float(parts[2])
                    energies.append({
                        'point': point,
                        'drms': drms,
                        'energy_kcal': energy,
                    })
                except (ValueError, IndexError):
                    if energies:  # we were reading data, now we hit something else
                        break
            else:
                if energies:
                    break

    return energies if energies else None


# ── Step 3: Normalize trajectory output ──────────────────────────────────────

def normalize_trajectory(result, out_dir, name):
    """
    Convert xTB path output to unified format:
    - {name}_MEP_trj.xyz  (multi-frame trajectory)
    - {name}_energy_profile.json
    """
    os.makedirs(out_dir, exist_ok=True)

    trj_out = os.path.join(out_dir, f'{name}_MEP_trj.xyz')
    profile_out = os.path.join(out_dir, f'{name}_energy_profile.json')

    # Copy trajectory
    if result.get('trajectory_xyz') and os.path.exists(result['trajectory_xyz']):
        shutil.copy2(result['trajectory_xyz'], trj_out)
        print(f"  ✓ Trajectory: {trj_out}")
    else:
        print(f"  ✗ No trajectory to normalize")
        return None, None

    # Copy TS geometry
    if result.get('ts_xyz') and os.path.exists(result['ts_xyz']):
        ts_out = os.path.join(out_dir, f'{name}_TS.xyz')
        shutil.copy2(result['ts_xyz'], ts_out)

    # Build energy profile JSON
    profile = {
        'method': result.get('method', 'xtb_path'),
        'barrier_forward_kcal': result.get('barrier_forward_kcal'),
        'barrier_reverse_kcal': result.get('barrier_reverse_kcal'),
        'reaction_energy_kcal': result.get('reaction_energy_kcal'),
        'frames': [],
    }

    # If we have parsed energies from the path data table, use those
    # (they're the clean, properly interpolated path points)
    if result.get('energies'):
        for j, e in enumerate(result['energies']):
            profile['frames'].append({
                'index': j,
                'energy_kcal': e.get('energy_kcal'),
                'drms': e.get('drms'),
            })
    else:
        # Fallback: parse energies from the trajectory file comments
        if os.path.exists(trj_out):
            with open(trj_out) as f:
                content = f.read()
            lines = content.strip().split('\n')
            frame_idx = 0
            i = 0
            while i < len(lines):
                try:
                    n_atoms = int(lines[i].strip())
                    comment = lines[i + 1].strip() if i + 1 < len(lines) else ''
                    energy = None
                    e_match = re.search(r'energy:\s*([-\d.]+)', comment)
                    if e_match:
                        energy = float(e_match.group(1))
                    profile['frames'].append({
                        'index': frame_idx,
                        'energy_kcal': energy,
                        'comment': comment,
                    })
                    frame_idx += 1
                    i += n_atoms + 2
                except (ValueError, IndexError):
                    break

    with open(profile_out, 'w') as f:
        json.dump(profile, f, indent=2)
    print(f"  ✓ Energy profile: {profile_out}")

    return trj_out, profile_out


# ── Step 2b: Relaxed scan fallback ──────────────────────────────────────────

def run_xtb_scan(reactant_xyz, driving_coords, charge=0, solvent='water',
                 work_dir=None, n_points=20):
    """
    Run a relaxed coordinate scan along the forming bond.
    Fallback for when --path can't connect R and P.

    Scans the forming bond from its current distance down to ~1.5Å (bonded),
    generating a trajectory of constrained optimizations.

    Args:
        reactant_xyz: path to (optimized) reactant complex XYZ
        driving_coords: dict with 'complex_formed_1based' and 'complex_broken_1based'
        charge: total system charge
        solvent: solvent name for ALPB
        work_dir: directory for calculation files
        n_points: number of scan points

    Returns:
        dict similar to run_xtb_path result
    """
    import math

    if work_dir is None:
        work_dir = tempfile.mkdtemp(prefix='xtb_scan_')
    os.makedirs(work_dir, exist_ok=True)

    atoms = _read_xyz(reactant_xyz)
    formed = driving_coords.get('complex_formed_1based', [])
    broken = driving_coords.get('complex_broken_1based', [])
    pt = driving_coords.get('complex_pt_1based', [])

    if not formed and not broken and not pt:
        return {'success': False, 'method': 'xtb_scan',
                'error': 'No driving coordinates for scan'}

    scans = []
    scan_labels = []
    
    if pt:
        # Proton transfer: scan the donor-acceptor distance AND the breaking proton bond simultaneously
        scan_pair = pt[0]  # 1-based (donor, acceptor)
        donor, acceptor = scan_pair[0], scan_pair[1]
        
        # Find the transferring proton
        proton = None
        for fbond in formed:
            for bbond in broken:
                shared = set(fbond).intersection(set(bbond))
                if shared:
                    proton = list(shared)[0]
                    break
            if proton: break
            
        if proton:
            # 1. Scan heavy atoms
            i, j = donor - 1, acceptor - 1
            si, xi, yi, zi = atoms[i]
            sj, xj, yj, zj = atoms[j]
            d_current = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
            d_start_heavy = min(max(d_current, 3.0), 3.5)
            d_end_heavy = 2.4  # H-bond distance
            scans.append({'pair': (donor, acceptor), 'start': d_start_heavy, 'end': d_end_heavy})
            scan_labels.append(f"PT heavy {si}({donor})-{sj}({acceptor})")
            
            # 2. Scan proton breaking
            i, j = donor - 1, proton - 1
            si, xi, yi, zi = atoms[i]
            sj, xj, yj, zj = atoms[j]
            d_start_pt = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
            d_end_pt = 1.45  # Perfect distance to deposit proton onto acceptor when heavy atoms are at 2.40A
            scans.append({'pair': (donor, proton), 'start': d_start_pt, 'end': d_end_pt})
            scan_labels.append(f"PT proton {si}({donor})-{sj}({proton})")
        else:
            # Fallback if proton not found
            i, j = donor - 1, acceptor - 1
            si, xi, yi, zi = atoms[i]
            sj, xj, yj, zj = atoms[j]
            d_current = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
            scans.append({'pair': (donor, acceptor), 'start': min(max(d_current, 3.0), 3.5), 'end': 2.4})
            scan_labels.append(f"PT {si}({donor})-{sj}({acceptor})")
            
    elif formed:
        scan_pair = formed[0]  # 1-based
        i, j = scan_pair[0] - 1, scan_pair[1] - 1
        si, xi, yi, zi = atoms[i]
        sj, xj, yj, zj = atoms[j]
        d_start = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
        ri = _COVALENT_RADII.get(si, 0.8)
        rj = _COVALENT_RADII.get(sj, 0.8)
        d_end = (ri + rj) * 1.1
        scans.append({'pair': scan_pair, 'start': d_start, 'end': d_end})
        scan_labels.append(f"forming {si}({scan_pair[0]})-{sj}({scan_pair[1]})")
        
    elif broken:
        scan_pair = broken[0]
        i, j = scan_pair[0] - 1, scan_pair[1] - 1
        si, xi, yi, zi = atoms[i]
        sj, xj, yj, zj = atoms[j]
        d_start = math.sqrt((xi - xj)**2 + (yi - yj)**2 + (zi - zj)**2)
        ri = _COVALENT_RADII.get(si, 0.8)
        rj = _COVALENT_RADII.get(sj, 0.8)
        d_end = max(8.0, (ri + rj) * 2.5)
        scans.append({'pair': scan_pair, 'start': d_start, 'end': d_end})
        scan_labels.append(f"breaking {si}({scan_pair[0]})-{sj}({scan_pair[1]})")

    print(f"  Scans:")
    for lbl, s in zip(scan_labels, scans):
        print(f"    - {lbl}, {s['start']:.2f} → {s['end']:.2f} Å")

    # Copy input
    input_xyz = os.path.join(work_dir, 'input.xyz')
    shutil.copy2(reactant_xyz, input_xyz)

    # Write scan input
    inp_path = os.path.join(work_dir, 'scan.inp')
    with open(inp_path, 'w') as f:
        f.write("$constrain\n")
        f.write("  force constant=1.0\n")
        for s in scans:
            f.write(f"  distance: {s['pair'][0]}, {s['pair'][1]}, auto\n")
            
        # Add angle constraint to ensure linear proton transfer
        if pt and proton:
            f.write(f"  angle: {donor}, {proton}, {acceptor}, auto\n")
            
        # Prevent ALL spontaneous side-reactions (like C-O SN1 heterolysis) 
        # by freezing the entire spectator backbone in its reactant state.
        reactive_atoms = set()
        for fbond in formed: reactive_atoms.update(fbond)
        for bbond in broken: reactive_atoms.update(bbond)
        
        for i, (sym_i, xi, yi, zi) in enumerate(atoms):
            for j in range(i + 1, len(atoms)):
                sym_j, xj, yj, zj = atoms[j]
                dist = math.sqrt((xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2)
                
                ri = _COVALENT_RADII.get(sym_i, 0.8)
                rj = _COVALENT_RADII.get(sym_j, 0.8)
                
                # If they are covalently bonded
                if dist < (ri + rj) * 1.3:
                    idx_i, idx_j = i + 1, j + 1
                    
                    # Allow the atoms explicitly participating in the mechanism to move freely
                    if (idx_i in reactive_atoms) and (idx_j in reactive_atoms):
                        continue
                        
                    # Never freeze the exact scan coordinates
                    is_scan_coord = False
                    for s in scans:
                        sp = s['pair']
                        if (idx_i == sp[0] and idx_j == sp[1]) or (idx_i == sp[1] and idx_j == sp[0]):
                            is_scan_coord = True
                    if is_scan_coord:
                        continue
                        
                    f.write(f"  distance: {idx_i}, {idx_j}, auto\n")
        f.write("$end\n")
        f.write("$scan\n")
        if len(scans) > 1:
            f.write("  mode=simultaneous\n")
        for idx, s in enumerate(scans):
            f.write(f"  {idx+1}: {s['start']:.4f}, {s['end']:.4f}, {n_points}\n")
        f.write("$end\n")
        f.write("$opt\n")
        f.write("  optlevel=normal\n")
        f.write("$end\n")

    args = ['input.xyz', '--opt', '--gfn', '2',
            '--chrg', str(charge), '--input', 'scan.inp']
            
    # CRITICAL: Do NOT use ALPB solvent during the scan for reactions with charged 
    # intermediates (like tBuOH2+). Implicit continuum solvents over-stabilize bare ions 
    # (like H+) causing spontaneous, artifactual dissociation before the actual TS is reached.
    # We rely on gas-phase scanning to produce the visually correct, drawn mechanism.
    # if solvent:
    #     args += ['--alpb', solvent]

    print(f"  Running xTB relaxed scan...")
    rc, output = _run_xtb(args, work_dir, timeout=600)

    # Save log
    log_path = os.path.join(work_dir, 'xtb_scan.log')
    with open(log_path, 'w') as f:
        f.write(output)

    result = {
        'success': False,
        'trajectory_xyz': None,
        'ts_xyz': None,
        'energies': None,
        'barrier_forward_kcal': None,
        'barrier_reverse_kcal': None,
        'reaction_energy_kcal': None,
        'method': 'xtb_scan',
        'log': log_path,
        'work_dir': work_dir,
    }

    # Check for scan trajectory
    scan_traj = os.path.join(work_dir, 'xtbscan.log')
    if os.path.exists(scan_traj) and os.path.getsize(scan_traj) > 0:
        result['trajectory_xyz'] = scan_traj

        # Parse energies from the trajectory
        with open(scan_traj) as f:
            content = f.read()
        lines = content.strip().split('\n')
        energies = []
        idx = 0
        i = 0
        while i < len(lines):
            try:
                n_atoms = int(lines[i].strip())
                comment = lines[i + 1].strip() if i + 1 < len(lines) else ''
                # xTB scan log has "energy:" in comment
                e_match = re.search(r'energy:\s*([-\d.]+)', comment)
                if e_match:
                    energies.append(float(e_match.group(1)))
                i += n_atoms + 2
                idx += 1
            except (ValueError, IndexError):
                break

        if energies:
            # Convert from Hartree to kcal/mol relative to first point
            e_min = energies[0]
            energies_kcal = [(e - e_min) * 627.509 for e in energies]
            barrier = max(energies_kcal)
            result['energies'] = [{'point': k, 'energy_kcal': e}
                                  for k, e in enumerate(energies_kcal)]
            result['barrier_forward_kcal'] = barrier
            result['reaction_energy_kcal'] = energies_kcal[-1] - energies_kcal[0]

        result['success'] = True
        print(f"  ✓ Scan complete: {idx} frames")
        if result.get('barrier_forward_kcal'):
            print(f"    Max energy: {result['barrier_forward_kcal']:.1f} kcal/mol")
    else:
        print(f"  ✗ Scan failed (exit code {rc})")

    return result


# ── Full pipeline ────────────────────────────────────────────────────────────

def run_reaction_path(reactant_xyz, product_xyz, charge=0, solvent='water',
                      out_dir='.', name='reaction', path_params=None,
                      driving_coords=None):
    """
    Full Tier-1 reaction path pipeline:
    1. Pre-optimize endpoints
    2. Run xTB --path
    3. If path fails and driving_coords available, run relaxed scan
    4. Normalize output

    Args:
        driving_coords: dict with 'complex_formed_1based' and 'complex_broken_1based'
                        from reaction metadata JSON

    Returns dict with trajectory path, energy profile, etc.
    """
    work_dir = os.path.join(out_dir, 'xtb_path')
    os.makedirs(work_dir, exist_ok=True)

    print(f"\n{'═'*60}")
    print(f" xTB Reaction Path Finder")
    print(f" {name}: {os.path.basename(reactant_xyz)} → {os.path.basename(product_xyz)}")
    if driving_coords:
        formed = driving_coords.get('complex_formed_1based', [])
        broken = driving_coords.get('complex_broken_1based', [])
        pt = driving_coords.get('complex_pt_1based', [])
        if formed or broken or pt:
            parts = []
            if formed:
                parts.append(f"forming: {formed}")
            if broken:
                parts.append(f"breaking: {broken}")
            if pt:
                parts.append(f"proton transfer: {pt}")
            print(f" Driving coords: {', '.join(parts)}")
    print(f"{'═'*60}")

    # Step 1: Optimize endpoints
    # Detect dissociation: if product has more fragments than reactant,
    # OR if driving_coords indicate bonds are breaking without forming,
    # constrain product optimization to prevent ion-pair recombination
    r_atoms_raw = _read_xyz(reactant_xyz)
    p_atoms_raw = _read_xyz(product_xyz)
    r_frags = _detect_fragments(r_atoms_raw)
    # Use a tighter scale for product fragment detection — leaving groups
    # may be placed just barely beyond the default covalent threshold
    p_frags = _detect_fragments(p_atoms_raw, scale=1.15)
    is_dissociation = len(p_frags) > len(r_frags) or len(p_frags) > 1
    # Also detect from driving coordinates: if we have broken bonds
    # but no formed bonds, it's a pure dissociation (E1 step 1, etc.)
    if not is_dissociation and driving_coords:
        has_broken = bool(driving_coords.get('complex_broken_1based'))
        has_formed = bool(driving_coords.get('complex_formed_1based'))
        if has_broken and not has_formed:
            is_dissociation = True
    if is_dissociation:
        print(f"\n  ⚡ Dissociation detected: {len(r_frags)} fragment(s) → "
              f"{len(p_frags)} fragment(s)")
        print(f"    Product optimization SKIPPED — fragments already optimized individually")
        print(f"    (unconstrained optimization would recombine the ion pair)")

    print(f"\n Step 1: Optimize endpoints")
    print(f"{'─'*40}")
    print(f"  ⚡ Using pre-optimized geometries from rxn_prep (fragments were optimized individually)")
    print(f"  ⚡ Skipping unconstrained complex optimization to preserve reactive attack angles")
    opt_r = reactant_xyz
    opt_p = product_xyz

    # Verify atom counts still match after optimization
    r_atoms = _read_xyz(opt_r)
    p_atoms = _read_xyz(opt_p)
    if len(r_atoms) != len(p_atoms):
        print(f"  ✗ Atom count mismatch after optimization: R={len(r_atoms)}, P={len(p_atoms)}")
        return {'success': False, 'error': 'Atom count mismatch'}

    if not all(r[0] == p[0] for r, p in zip(r_atoms, p_atoms)):
        print(f"  ✗ Atom type mismatch after optimization")
        return {'success': False, 'error': 'Atom type mismatch'}

    # Step 2: Run path finder or relaxed scan
    has_formed = bool(driving_coords and driving_coords.get('complex_formed_1based'))
    has_broken = bool(driving_coords and driving_coords.get('complex_broken_1based'))
    has_pt = bool(driving_coords and driving_coords.get('complex_pt_1based'))
    should_scan = has_formed or has_broken or has_pt

    result = {'success': False}
    
    if should_scan and driving_coords:
        print(f"\n Step 2: Relaxed scan")
        print(f"{'─'*40}")
        scan_dir = os.path.join(out_dir, 'xtb_scan')
        result = run_xtb_scan(opt_r, driving_coords, charge=charge,
                              solvent=solvent, work_dir=scan_dir)
        
        if not result.get('success'):
            print(f"  ⚠ Scan failed or returned zero barrier, falling back to xTB path finder...")

    if not result.get('success'):
        print(f"\n Step 2: xTB path finder")
        print(f"{'─'*40}")
        result = run_xtb_path(opt_r, opt_p, charge=charge, solvent=solvent,
                              work_dir=work_dir, path_params=path_params)

        # Step 2b: If path failed and we have driving coords, try relaxed scan
        # Use ORIGINAL (unoptimized) reactant so fragment distances are preserved
        if not result.get('success') or ((result.get('barrier_forward_kcal') or 0) < 0.5):
            if driving_coords and (driving_coords.get('complex_formed_1based') or
                                   driving_coords.get('complex_broken_1based') or
                                   driving_coords.get('complex_pt_1based')):
                print(f"\n Step 2b: Relaxed scan fallback")
                print(f"{'─'*40}")
                scan_dir = os.path.join(out_dir, 'xtb_scan')
                scan_result = run_xtb_scan(reactant_xyz, driving_coords, charge=charge,
                                           solvent=solvent, work_dir=scan_dir)
                if scan_result.get('success'):
                    result = scan_result

    # Step 2c: If EVERYTHING failed (not just barrierless), try DE-GSM
    # A successful scan with barrier=0 is a valid barrierless result, NOT a failure
    if not result.get('success'):
        gsm_bin = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'gsm', 'gsm.orca')
        if os.path.exists(gsm_bin):
            print(f"\n Step 2c: DE-GSM fallback")
            print(f"{'─'*40}")
            try:
                from rxn_path_gsm import run_gsm
                gsm_result = run_gsm(
                    opt_r, opt_p,
                    charge=charge,
                    name=name,
                    out_dir=out_dir,
                    nnodes=15,
                    driving_coords=driving_coords,
                )
                if gsm_result.get('success'):
                    result = gsm_result
            except Exception as e:
                print(f"  ✗ GSM error: {e}")

    # Step 3: Normalize output
    if result.get('trajectory_xyz'):
        print(f"\n Step 3: Normalize output")
        print(f"{'─'*40}")
        trj, profile = normalize_trajectory(result, out_dir, name)
        result['normalized_trajectory'] = trj
        result['energy_profile_json'] = profile

    # Summary
    print(f"\n{'═'*60}")
    if result['success']:
        print(f" ✓ Reaction path found! (method: {result.get('method', 'unknown')})")
        if result.get('barrier_forward_kcal') is not None:
            print(f"   Forward barrier:  {result['barrier_forward_kcal']:6.1f} kcal/mol")
        if result.get('barrier_reverse_kcal') is not None:
            print(f"   Reverse barrier:  {result['barrier_reverse_kcal']:6.1f} kcal/mol")
        if result.get('reaction_energy_kcal') is not None:
            print(f"   Reaction energy:  {result['reaction_energy_kcal']:6.1f} kcal/mol")
    else:
        print(f" ✗ Path finder failed")
    print(f"{'═'*60}\n")

    return result


# ── CLI ──────────────────────────────────────────────────────────────────────

def _load_driving_coords(meta_path):
    """Load driving coordinates from reaction metadata JSON."""
    try:
        with open(meta_path) as f:
            meta = json.load(f)
        dc = meta.get('driving_coordinates', {})
        if dc.get('complex_formed_1based') or dc.get('complex_broken_1based') or dc.get('complex_pt_1based'):
            return dc
    except Exception:
        pass
    return None


def main():
    parser = argparse.ArgumentParser(
        description='xTB reaction path finder (Tier 1)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s reactant.xyz product.xyz --charge -1 -o output/
  %(prog)s r.xyz p.xyz --charge 0 --solvent water --name sn2
  %(prog)s r.xyz p.xyz --charge -1 --meta reaction.json
        """
    )
    parser.add_argument('reactant_xyz', help='Reactant complex .xyz')
    parser.add_argument('product_xyz', help='Product complex .xyz')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory')
    parser.add_argument('--name', '-n', default='reaction', help='Base name for output files')
    parser.add_argument('--charge', type=int, default=0, help='Total system charge')
    parser.add_argument('--solvent', default='water', help='Solvent for ALPB (or "none")')
    parser.add_argument('--npoint', type=int, default=25, help='Number of path points')
    parser.add_argument('--nrun', type=int, default=1, help='Number of path trial runs')
    parser.add_argument('--meta', default=None,
                        help='Reaction metadata JSON (for driving coordinates)')

    args = parser.parse_args()

    for f in [args.reactant_xyz, args.product_xyz]:
        if not os.path.exists(f):
            print(f"ERROR: {f} not found", file=sys.stderr)
            sys.exit(1)

    solvent = None if args.solvent.lower() == 'none' else args.solvent
    path_params = {}
    if args.npoint != 25:
        path_params['npoint'] = args.npoint
    if args.nrun != 1:
        path_params['nrun'] = args.nrun

    # Load driving coordinates from metadata if available
    driving_coords = None
    if args.meta and os.path.exists(args.meta):
        driving_coords = _load_driving_coords(args.meta)
        if driving_coords:
            print(f"Loaded driving coordinates from {args.meta}")
    else:
        # Auto-detect metadata in output directory
        auto_meta = os.path.join(args.output_dir, f"{args.name}_reaction.json")
        if os.path.exists(auto_meta):
            driving_coords = _load_driving_coords(auto_meta)
            if driving_coords:
                print(f"Auto-loaded driving coordinates from {auto_meta}")

    result = run_reaction_path(
        args.reactant_xyz, args.product_xyz,
        charge=args.charge, solvent=solvent,
        out_dir=args.output_dir, name=args.name,
        path_params=path_params or None,
        driving_coords=driving_coords,
    )

    sys.exit(0 if result['success'] else 1)


if __name__ == '__main__':
    main()
