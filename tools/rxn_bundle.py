#!/usr/bin/env python3
"""
rxn_bundle — Assemble reaction output into a clean bundle structure

Organizes flat output files into a structured bundle directory with
a master manifest (bundle.json) for consumption by slide widgets.

Usage:
    python rxn_bundle.py output/sn2/ --name sn2
"""

import argparse
import glob
import json
import os
import shutil
import sys


def assemble_bundle(job_dir, name):
    """Organize output files into a clean bundle structure.

    Creates:
        {job_dir}/molecules/       — individual fragment XYZ, SVG, orbitals
        {job_dir}/mep_orbitals/    — keyframe orbital cubes (already exists if computed)
        {job_dir}/{name}_bundle.json — master manifest

    Returns path to bundle.json
    """
    print(f"  Assembling bundle for '{name}' in {job_dir}")

    # ── Create subdirectories ────────────────────────────────────────
    mol_dir = os.path.join(job_dir, 'molecules')
    mep_dir = os.path.join(job_dir, 'mep_orbitals')
    os.makedirs(mol_dir, exist_ok=True)
    os.makedirs(mep_dir, exist_ok=True)

    # ── Read reaction metadata ───────────────────────────────────────
    meta_path = os.path.join(job_dir, f"{name}_reaction.json")
    meta = {}
    if os.path.exists(meta_path):
        with open(meta_path) as f:
            meta = json.load(f)

    reactants_smi = meta.get('reactants', [])
    products_smi = meta.get('products', [])
    total_charge = meta.get('total_charge', 0)
    driving_coords = meta.get('driving_coordinates', {})

    # ── Collect molecule fragments ───────────────────────────────────
    molecules = []

    # Find all individual fragment XYZ files (not complex, not trajectory)
    for f in sorted(os.listdir(job_dir)):
        if not f.endswith('.xyz'):
            continue
        if 'complex' in f or 'MEP' in f or 'TS' in f:
            continue
        if not f.startswith(name + '_'):
            continue

        frag_name = f.replace('.xyz', '')
        # Determine role: r1, r2, ... or p1, p2, ...
        suffix = frag_name.replace(name + '_', '')
        if not (suffix.startswith('r') or suffix.startswith('p')):
            continue

        role = 'reactant' if suffix.startswith('r') else 'product'
        idx = int(suffix[1:]) - 1 if suffix[1:].isdigit() else 0
        smiles = ''
        if role == 'reactant' and idx < len(reactants_smi):
            smiles = reactants_smi[idx]
        elif role == 'product' and idx < len(products_smi):
            smiles = products_smi[idx]

        mol_info = {
            'key': suffix,
            'role': role,
            'smiles': smiles,
            'xyz': f'molecules/{f}',
        }

        # Copy XYZ to molecules/
        src = os.path.join(job_dir, f)
        dst = os.path.join(mol_dir, f)
        if not os.path.exists(dst):
            shutil.copy2(src, dst)

        # Copy SVG if exists
        svg_name = f"{frag_name}.svg"
        svg_src = os.path.join(job_dir, svg_name)
        if os.path.exists(svg_src):
            svg_dst = os.path.join(mol_dir, svg_name)
            if not os.path.exists(svg_dst):
                shutil.copy2(svg_src, svg_dst)
            mol_info['svg'] = f'molecules/{svg_name}'

        # Copy orbital manifest and cubes if exist
        orb_name = f"{frag_name}_orbitals.json"
        orb_src = os.path.join(job_dir, orb_name)
        if os.path.exists(orb_src):
            orb_dst = os.path.join(mol_dir, orb_name)
            if not os.path.exists(orb_dst):
                shutil.copy2(orb_src, orb_dst)
            mol_info['orbitals'] = f'molecules/{orb_name}'

            # Copy associated cube files
            with open(orb_src) as of:
                orb_data = json.load(of)
            _copy_cube_files(orb_data, job_dir, mol_dir)

        molecules.append(mol_info)

    print(f"  {len(molecules)} molecule fragments")

    # ── Trajectory info ──────────────────────────────────────────────
    trj_path = os.path.join(job_dir, f"{name}_MEP_trj.xyz")
    profile_path = os.path.join(job_dir, f"{name}_energy_profile.json")
    ts_path = os.path.join(job_dir, f"{name}_TS.xyz")

    trajectory = None
    if os.path.exists(trj_path):
        # Count frames
        with open(trj_path) as f:
            first_line = f.readline().strip()
            n_atoms = int(first_line) if first_line.isdigit() else 0

        n_frames = 0
        if n_atoms > 0:
            total_lines = sum(1 for _ in open(trj_path))
            n_frames = total_lines // (n_atoms + 2)

        trajectory = {
            'xyz': f"{name}_MEP_trj.xyz",
            'n_frames': n_frames,
            'n_atoms': n_atoms,
        }

        if os.path.exists(profile_path):
            trajectory['energy_profile'] = f"{name}_energy_profile.json"
            with open(profile_path) as f:
                ep = json.load(f)
            trajectory['method'] = ep.get('method', 'unknown')

        if os.path.exists(ts_path):
            trajectory['ts_xyz'] = f"{name}_TS.xyz"

    # ── Energy data ──────────────────────────────────────────────────
    barrier_fwd = None
    barrier_rev = None
    rxn_energy = None
    method = None
    if os.path.exists(profile_path):
        with open(profile_path) as f:
            ep = json.load(f)
        barrier_fwd = ep.get('barrier_forward_kcal')
        barrier_rev = ep.get('barrier_reverse_kcal')
        rxn_energy = ep.get('reaction_energy_kcal')
        method = ep.get('method')

    # ── MEP orbitals ─────────────────────────────────────────────────
    mep_orb_path = os.path.join(job_dir, f"{name}_mep_orbitals.json")
    mep_orbitals = None
    if os.path.exists(mep_orb_path):
        mep_orbitals = f"{name}_mep_orbitals.json"

    # ── Build bundle manifest ────────────────────────────────────────
    # Build title from SMILES
    r_str = ' + '.join(reactants_smi) if reactants_smi else '?'
    p_str = ' + '.join(products_smi) if products_smi else '?'
    title = f"{r_str} → {p_str}"

    bundle = {
        'name': name,
        'title': title,
        'reaction_smiles': f"{'.' .join(reactants_smi)}>>{'.' .join(products_smi)}",
        'charge': total_charge,
        'method': method,
        'barrier_forward_kcal': barrier_fwd,
        'barrier_reverse_kcal': barrier_rev,
        'reaction_energy_kcal': rxn_energy,
        'trajectory': trajectory,
        'molecules': molecules,
        'mep_orbitals': mep_orbitals,
        'driving_coordinates': driving_coords if driving_coords else None,
        'complexes': {
            'reactant': f"{name}_reactant_complex.xyz",
            'product': f"{name}_product_complex.xyz",
        },
        'metadata': f"{name}_reaction.json",
    }

    bundle_path = os.path.join(job_dir, f"{name}_bundle.json")
    with open(bundle_path, 'w') as f:
        json.dump(bundle, f, indent=2)

    print(f"  ✓ Bundle manifest: {bundle_path}")

    # Summary
    print(f"\n  Bundle contents:")
    print(f"    Molecules:    {len(molecules)}")
    if trajectory:
        print(f"    Trajectory:   {trajectory['n_frames']} frames × {trajectory['n_atoms']} atoms")
    if barrier_fwd is not None:
        print(f"    Barrier:      {barrier_fwd:.1f} kcal/mol")
    if mep_orbitals:
        with open(mep_orb_path) as f:
            mep_data = json.load(f)
        if mep_data.get('mode') == 'per_frame':
            print(f"    MEP orbitals: {len(mep_data.get('frames', []))} frames (per-frame movie)")
        else:
            print(f"    MEP orbitals: {len(mep_data.get('keyframes', []))} keyframes")

    return bundle_path


def _copy_cube_files(orb_data, src_dir, dst_dir):
    """Copy cube files referenced in an orbital manifest."""
    for section in ['canonical', 'localized']:
        data = orb_data.get(section, {})
        if isinstance(data, dict):
            for key, val in data.items():
                if isinstance(val, dict) and 'file' in val:
                    _do_copy(val['file'], src_dir, dst_dir)
                elif isinstance(val, list):
                    for item in val:
                        if isinstance(item, dict) and 'file' in item:
                            _do_copy(item['file'], src_dir, dst_dir)


def _do_copy(filename, src_dir, dst_dir):
    """Copy a single file if it exists and not already copied."""
    src = os.path.join(src_dir, filename)
    dst = os.path.join(dst_dir, filename)
    if os.path.exists(src) and not os.path.exists(dst):
        shutil.copy2(src, dst)


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description='Assemble reaction bundle')
    parser.add_argument('job_dir', help='Job output directory')
    parser.add_argument('--name', required=True, help='Reaction name')
    args = parser.parse_args()

    print(f"\n{'═' * 60}")
    print(f" Reaction Bundle Assembler")
    print(f"{'═' * 60}\n")

    result = assemble_bundle(args.job_dir, args.name)

    print(f"\n{'═' * 60}")
    print(f" ✓ Bundle assembled: {result}")
    print(f"{'═' * 60}\n")


if __name__ == '__main__':
    main()
