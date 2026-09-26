#!/usr/bin/env python3
"""
rxn_orbitals — Intrinsic Bond Orbital (IBO) movie along the MEP

Computes IBOs of the full complex at each frame. IBOs are localized
to individual bonds and lone pairs, and smoothly evolve during the
reaction, showing electron flow during bond formation/breaking.

Usage:
    python rxn_orbitals.py output/sn2/ --name sn2 --charge -1
"""
import argparse, json, os, sys, time
import numpy as np

os.environ['HDF5_USE_FILE_LOCKING'] = 'FALSE'

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chem_constants import (
    ANG2BOHR, GRID_POINTS, GRID_MARGIN, BASIS, METHOD,
    CORE_ENERGY_THRESH, ATOMIC_NUMBERS, CORE_ELECTRONS,
)
from cube_writer import make_grid, write_cube_file, eval_mo_on_grid, make_global_grid

def parse_trajectory(path):
    frames = []
    with open(path) as f:
        lines = f.readlines()
    i = 0
    while i < len(lines):
        try:
            n = int(lines[i].strip())
            comment = lines[i+1].strip() if i+1 < len(lines) else ''
            atoms = []
            for j in range(n):
                p = lines[i+2+j].split()
                atoms.append((p[0], float(p[1]), float(p[2]), float(p[3])))
            frames.append({'atoms': atoms, 'comment': comment})
            i += n + 2
        except (ValueError, IndexError):
            break
    return frames

def classify_ibo(mol, ibo_vec):
    """Classify one IBO → (type, label, pop_vector)."""
    ovlp = mol.intor_symmetric('int1e_ovlp')
    ao_labels = mol.ao_labels()
    ao_atom = [int(l.split()[0]) for l in ao_labels]
    dm = np.outer(ibo_vec, ibo_vec)
    ps = dm @ ovlp
    pop = np.zeros(mol.natm)
    for mu in range(mol.nao):
        pop[ao_atom[mu]] += ps[mu, mu]
    pop *= 2  # 2 electrons

    atom_syms = [mol.atom_symbol(i) for i in range(mol.natm)]
    # Find atoms with significant population
    ranked = sorted(range(mol.natm), key=lambda i: -pop[i])
    top = [(i, pop[i]) for i in ranked if pop[i] > 0.15]

    if len(top) == 1 or (len(top) >= 1 and top[0][1] > 1.85 and (len(top) < 2 or top[1][1] < 0.15)):
        ai = top[0][0]
        label = f'LP({atom_syms[ai]}{ai+1})'
        typ = 'lp'
    elif len(top) >= 2:
        a1, a2 = top[0][0], top[1][0]
        # Sort indices to prevent label flipping (e.g., C1-C2 to C2-C1) from tiny population shifts
        if a1 > a2:
            a1, a2 = a2, a1
        s1, s2 = atom_syms[a1], atom_syms[a2]
        label = f'{s1}{a1+1}-{s2}{a2+1}'
        typ = 'bond'
    else:
        label = 'core'
        typ = 'core'
    return typ, label, pop

def compute_ibo_energies(mol, mf, ibo_coeff):
    """Compute orbital energy of each IBO from diagonal of Fock matrix."""
    fock = mf.get_fock()
    # UKS returns (fock_alpha, fock_beta) — use alpha since IBOs are from alpha MOs
    if isinstance(fock, (tuple, list)) or (isinstance(fock, np.ndarray) and fock.ndim == 3):
        fock = fock[0]
    ibo_fock = ibo_coeff.T @ fock @ ibo_coeff
    return np.diag(ibo_fock) * 27.2114  # Hartree → eV

def is_core_ibo(energy_ev):
    """Core IBOs have very negative orbital energies (< -40 eV)."""
    return energy_ev < CORE_ENERGY_THRESH

def match_ibos(pops_prev, pops_curr):
    """Match current IBOs to previous frame by population fingerprint.
    Returns permutation: result[i] = index in curr that matches prev[i]."""
    from scipy.optimize import linear_sum_assignment
    n = len(pops_prev)
    m = len(pops_curr)
    cost = np.zeros((n, m))
    for i in range(n):
        for j in range(m):
            cost[i,j] = np.sum((pops_prev[i] - pops_curr[j])**2)
    row_ind, col_ind = linear_sum_assignment(cost)
    perm = list(col_ind)
    return perm

def fix_phase(ibo_vec, prev_vec, mol):
    """Ensure consistent phase by checking sign of overlap with prev."""
    ovlp = mol.intor_symmetric('int1e_ovlp')
    s = prev_vec @ ovlp @ ibo_vec
    if s < 0:
        return -ibo_vec
    return ibo_vec

def compute_ibo_movie(job_dir, name, charge=0, spin=0):
    from pyscf import gto, dft, lo

    trj_path = os.path.join(job_dir, f'{name}_MEP_trj.xyz')
    if not os.path.exists(trj_path):
        print(f'  ✗ Trajectory not found: {trj_path}'); return None

    frames = parse_trajectory(trj_path)
    n_frames = len(frames)
    n_atoms = len(frames[0]['atoms'])
    print(f'  {n_frames} frames, {n_atoms} atoms')

    orb_dir = os.path.join(job_dir, 'mep_orbitals')
    os.makedirs(orb_dir, exist_ok=True)

    # Phase 1: compute IBOs at all frames, classify, track
    all_ibos = []    # [frame_idx] → list of {vec, type, label, pop, is_core}
    prev_pops = None
    prev_vecs = None
    all_lumos = []
    prev_lumo = None
    chkfile = os.path.join(orb_dir, '.chk_ibo.chk')

    # Auto-detect odd electron count and adjust spin (check once on first frame)
    frame_spin = spin
    first_atoms = frames[0]['atoms']
    total_nuclear_charge = sum(ATOMIC_NUMBERS.get(s, 0) for s, x, y, z in first_atoms)
    n_electrons = total_nuclear_charge - charge
    if n_electrons % 2 != 0 and frame_spin == 0:
        frame_spin = 1  # doublet for odd-electron systems
        print(f'  ⚠ Odd electron count ({n_electrons}), using spin=1 (doublet)')

    for fi, frame in enumerate(frames):
        t0 = time.time()
        atoms = frame['atoms']
        atom_str = '; '.join(f'{s} {x} {y} {z}' for s,x,y,z in atoms)

        mol = gto.M(atom=atom_str, basis=BASIS, charge=charge, spin=frame_spin, verbose=0)

        if frame_spin > 0:
            mf = dft.UKS(mol)
        else:
            mf = dft.RKS(mol)
        mf.xc = METHOD
        mf.verbose = 0
        mf.chkfile = chkfile
        if fi > 0 and os.path.exists(chkfile):
            mf.init_guess = 'chkfile'
        mf.kernel()

        if frame_spin > 0:
            # Open-shell: use alpha MOs
            nocc = mol.nelec[0]
            occ = mf.mo_coeff[0][:, :nocc]
        else:
            nocc = mol.nelectron // 2
            occ = mf.mo_coeff[:, :nocc]
        ibo_coeff = lo.ibo.ibo(mol, occ, verbose=0)

        # Compute IBO energies and classify
        ibo_energies = compute_ibo_energies(mol, mf, ibo_coeff)
        frame_ibos = []
        pops = []
        for k in range(ibo_coeff.shape[1]):
            vec = ibo_coeff[:, k]
            typ, label, pop = classify_ibo(mol, vec)
            core = is_core_ibo(ibo_energies[k])
            if core:
                typ = 'core'
            frame_ibos.append({'vec': vec, 'type': typ, 'label': label,
                               'pop': pop, 'is_core': core, 'energy_ev': float(ibo_energies[k])})
            pops.append(pop)

        pops = np.array(pops)

        # Match to previous frame
        if prev_pops is not None:
            perm = match_ibos(prev_pops, pops)
            frame_ibos = [frame_ibos[perm[i]] for i in range(len(perm))]
            pops = pops[perm]
            # Fix phase
            for k in range(len(frame_ibos)):
                if prev_vecs is not None and k < len(prev_vecs):
                    frame_ibos[k]['vec'] = fix_phase(frame_ibos[k]['vec'],
                                                      prev_vecs[k], mol)

        prev_pops = pops
        prev_vecs = [ib['vec'] for ib in frame_ibos]
        all_ibos.append(frame_ibos)
        
        # Track LUMO (lowest unoccupied alpha MO)
        lumo_vec = mf.mo_coeff[0][:, nocc] if frame_spin > 0 else mf.mo_coeff[:, nocc]
        
        # Phase match LUMO to prevent flipping
        if prev_lumo is not None:
            lumo_vec = fix_phase(lumo_vec, prev_lumo, mol)
        prev_lumo = lumo_vec
        all_lumos.append(lumo_vec)

        dt = time.time() - t0

        n_val = sum(1 for ib in frame_ibos if not ib['is_core'])
        print(f'  ✓ Frame {fi:3d}/{n_frames-1}  '
              f'{n_val} valence IBOs + 1 LUMO ({dt:.1f}s)')

    # Phase 2: determine significance using explicit driving coordinates
    first = all_ibos[0]
    last = all_ibos[-1]
    n_ibos = len(first)
    
    import re
    formed_pairs = []
    broken_pairs = []
    rxn_json_path = os.path.join(job_dir, f'{name}_reaction.json')
    if os.path.exists(rxn_json_path):
        try:
            with open(rxn_json_path) as f:
                rxn_data = json.load(f)
                dc = rxn_data.get('driving_coordinates', {})
                for pair in dc.get('complex_formed_1based', []):
                    formed_pairs.append(set(pair))
                for pair in dc.get('complex_broken_1based', []):
                    broken_pairs.append(set(pair))
        except Exception as e:
            print(f"  ⚠ Failed to load explicit reaction driving coords: {e}")

    ibo_meta = []
    for k in range(n_ibos):
        if first[k]['is_core']:
            continue
        label_start = first[k]['label']
        label_end = last[k]['label']
        
        # Extract 1-based atom indices from the labels
        start_atoms = set(int(x) for x in re.findall(r'\d+', label_start))
        end_atoms = set(int(x) for x in re.findall(r'\d+', label_end))
        
        significant = False
        if formed_pairs or broken_pairs:
            # Explicit matching: the orbital IS the bond that breaks or forms
            if any(start_atoms == p for p in broken_pairs):
                significant = True
            if any(end_atoms == p for p in formed_pairs):
                significant = True
        else:
            # Fallback if no explicit mechanism was provided
            pop_diff = np.max(np.abs(first[k]['pop'] - last[k]['pop']))
            significant = bool(label_start != label_end) or bool(pop_diff > 0.25)
        
        # Label formatting
        cat_sig = bool(label_start != label_end)
        if cat_sig:
            label = f'{label_start} → {label_end}'
        elif significant:
            label = f'{label_start} (reactive)'
        else:
            label = label_start
            
        ibo_meta.append({
            'ibo_index': k,
            'label': label,
            'label_start': label_start,
            'label_end': label_end,
            'significant': significant,
            'type_start': first[k]['type'],
            'type_end': last[k]['type'],
        })

    sig_count = sum(1 for m in ibo_meta if m['significant'])
    print(f'\n  {len(ibo_meta)} valence IBOs, {sig_count} significant (character change)')
    for m in ibo_meta:
        tag = ' ★' if m['significant'] else ''
        print(f"    [{m['ibo_index']:2d}] {m['label']}{tag}")
        
    # Always include LUMO for reference, but NEVER mark it as significant 
    # so it doesn't clutter the user's localized reaction view.
    ibo_meta.append({
        'ibo_index': 'lumo',
        'label': 'LUMO (p-orbital)',
        'label_start': 'LUMO',
        'label_end': 'LUMO',
        'significant': False,
        'type_start': 'lumo',
        'type_end': 'lumo'
    })
    print(f'\n  Writing cube files...')
    manifest_frames = []
    t_cubes = time.time()
    
    # Compute a universal grid bounding box across all frames
    global_origin, global_step, global_npts = make_global_grid([f['atoms'] for f in frames])

    for fi, frame in enumerate(frames):
        atoms = frame['atoms']
        origin, step, npts = global_origin, global_step, global_npts
        atom_str = '; '.join(f'{s} {x} {y} {z}' for s,x,y,z in atoms)
        mol = gto.M(atom=atom_str, basis=BASIS, charge=charge, spin=frame_spin, verbose=0)

        frame_entry = {'frame_index': fi}
        for meta in ibo_meta:
            k = meta['ibo_index']
            key = f'ibo_{k}'
            if k == 'lumo':
                vec = all_lumos[fi]
            else:
                vec = all_ibos[fi][k]['vec']
            
            vals = eval_mo_on_grid(mol, vec, origin, step, npts)
            cube_name = f'{name}_f{fi:03d}_{key}.cube'
            write_cube_file(os.path.join(orb_dir, cube_name),
                       f'IBO {k}: {meta["label"]}',
                       atoms, origin, step, npts, vals)
            frame_entry[key] = {'file': cube_name}

        manifest_frames.append(frame_entry)
        if fi % 4 == 0 or fi == n_frames-1:
            print(f'    Frame {fi}/{n_frames-1}')

    dt_cubes = time.time() - t_cubes

    # Assign colors
    COLORS = ['#ff5c5c','#55aaff','#33cc66','#ff9933','#c77dff','#ffdd57','#ff6b9d','#00d4aa']
    for ci, m in enumerate(ibo_meta):
        m['color'] = COLORS[ci % len(COLORS)]

    manifest = {
        'name': name, 'method': METHOD, 'basis': BASIS, 'charge': charge,
        'grid_points': GRID_POINTS, 'n_frames': n_frames, 'n_atoms': n_atoms,
        'mode': 'ibo',
        'ibos': ibo_meta,
        'orbital_keys': [f'ibo_{m["ibo_index"]}' for m in ibo_meta],
        'frames': manifest_frames,
    }

    manifest_path = os.path.join(orb_dir, f'{name}_mep_orbitals.json')
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    top = os.path.join(job_dir, f'{name}_mep_orbitals.json')
    with open(top, 'w') as f:
        json.dump(manifest, f, indent=2)

    if os.path.exists(chkfile):
        os.remove(chkfile)

    n_cubes = len(ibo_meta) * n_frames
    print(f'\n  ✓ {n_cubes} cube files ({len(ibo_meta)} IBOs × {n_frames} frames)')
    print(f'  ✓ Cubes written in {dt_cubes:.0f}s')
    print(f'  ✓ Manifest: {manifest_path}')
    return manifest_path

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('job_dir')
    ap.add_argument('--name', required=True)
    ap.add_argument('--charge', type=int, default=0)
    ap.add_argument('--spin', type=int, default=0)
    args = ap.parse_args()
    print(f"\n{'═'*60}")
    print(f' IBO Movie — Intrinsic Bond Orbitals along MEP')
    print(f' {args.name}: {args.job_dir}')
    print(f"{'═'*60}\n")
    r = compute_ibo_movie(args.job_dir, args.name, args.charge, args.spin)
    if r:
        print(f"\n{'═'*60}\n ✓ IBO movie complete\n{'═'*60}\n")
    else:
        print('\n ✗ Failed'); sys.exit(1)

if __name__ == '__main__':
    main()
