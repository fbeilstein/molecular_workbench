#!/usr/bin/env python3
"""
rxn_export — Package a completed reaction job into a data-only .rxnbundle.zip

The archive contains ONLY data: JSON manifests, cube files, XYZ text.
All compression is handled by ZIP_DEFLATED — no double-wrapping with gzip.

Usage:
    python rxn_export.py output/sn2/
    python rxn_export.py output/sn2/ --clean   # delete intermediates after export
"""
import argparse, json, os, shutil, sys, zipfile, glob

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cube_converter import CubeConverter


def load_json(path):
    with open(path) as f:
        return json.load(f)


def load_text(path):
    with open(path) as f:
        return f.read()


def parse_xyz_frames(text):
    """Parse multi-frame XYZ into list of dicts."""
    frames, lines, i = [], text.strip().split('\n'), 0
    while i < len(lines):
        try:
            n = int(lines[i].strip())
            comment = lines[i+1].strip() if i+1 < len(lines) else ''
            atoms = []
            for j in range(n):
                p = lines[i+2+j].split()
                atoms.append({'sym':p[0],'x':float(p[1]),'y':float(p[2]),'z':float(p[3])})
            frames.append({'atoms': atoms, 'comment': comment})
            i += n + 2
        except (ValueError, IndexError):
            break
    return frames






def _get_charge_from_smiles(smiles):
    mol_charge = 0
    if '-]' in smiles: mol_charge = -smiles.count('-]')
    elif '+]' in smiles: mol_charge = smiles.count('+]')
    try:
        import rdkit.Chem as Chem
        m = Chem.MolFromSmiles(smiles)
        if m:
            mol_charge = sum(atom.GetFormalCharge() for atom in m.GetAtoms())
    except:
        pass
    return mol_charge

def _process_molecule_orbitals(mol_info, job_dir, cube_dir, zf):
    key, role = mol_info['key'], mol_info['role']
    smiles = mol_info.get('smiles', '')
    mol_data = {'key': key, 'role': role, 'smiles': smiles, 'orbitals': []}

    xyz_path = __import__('os').path.join(job_dir, mol_info.get('xyz', ''))
    if __import__('os').path.exists(xyz_path):
        mol_data['xyz'] = load_text(xyz_path)
    
    if not mol_data.get('xyz'):
        return mol_data
        
    mol_charge = _get_charge_from_smiles(smiles)
    
    try:
        from mol_orbitals import compute_all_orbitals
        orb_manifest = compute_all_orbitals(
            xyz_file=xyz_path, name=key, out_dir=cube_dir, 
            charge=mol_charge, method='b3lyp', smiles=smiles
        )
        
        for cf in __import__('os').listdir(cube_dir):
            if cf.endswith('.cube'):
                CubeConverter.process_cube(__import__('os').path.join(cube_dir, cf), adaptive=True)
                
        from orbital_packager import package_orbitals
        orbital_groups, esp_surface = package_orbitals(orb_manifest, cube_dir, zf)
        
        mol_data['orbitals'] = orbital_groups
        if esp_surface:
            mol_data['esp_surface'] = esp_surface
            
        print(f"    ✓ {key} ({smiles}): {sum(len(g['items']) for g in mol_data['orbitals'])} orbitals")
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"    ✗ {key} ({smiles}): computation failed ({e})")
        
    return mol_data

def export_bundle(job_dir, output_path=None, significant_only=False, clean=False):
    """Package job_dir into a data-only .rxnbundle.zip.

    Cube files are stored raw — ZIP_DEFLATED handles all compression
    in one pass for better ratios than double gzip+zip.
    """
    bundle_files = [f for f in os.listdir(job_dir) if f.endswith('_bundle.json')]
    if not bundle_files:
        print(f'  ✗ No _bundle.json in {job_dir}')
        return None

    bundle = load_json(os.path.join(job_dir, bundle_files[0]))
    name = bundle['name']
    if output_path is None:
        output_path = os.path.join(job_dir, f'{name}.rxnbundle.zip')

    print(f'  Exporting: {name}\n  Source: {job_dir}\n  Output: {output_path}')

    manifest = {
        'name': name,
        'title': bundle.get('title', name),
        'smiles': bundle.get('smiles', ''),
        'reaction_smiles': bundle.get('smiles', ''),
        'charge': bundle.get('charge', 0),
        'engine': bundle.get('engine', 'xtb'),
        'method': bundle.get('method', 'unknown'),
        'barrier_forward_kcal': bundle.get('barrier_forward_kcal'),
        'barrier_reverse_kcal': bundle.get('barrier_reverse_kcal'),
        'reaction_energy_kcal': bundle.get('reaction_energy_kcal'),
        'molecules': [],
        'has_trajectory': False,
        'has_orbitals': False,
    }

    with zipfile.ZipFile(output_path, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        # --- Reaction SVG (single file for the whole reaction) ---
        rxn_svg_path = os.path.join(job_dir, f'{name}_reaction.svg')
        if os.path.exists(rxn_svg_path):
            zf.write(rxn_svg_path, 'reaction.svg')
            manifest['has_reaction_svg'] = True
            
        # --- Original Ketcher Layout ---
        ket_path = os.path.join(job_dir, f'{name}.ket')
        if os.path.exists(ket_path):
            zf.write(ket_path, 'layout.ket')
            manifest['ket_file'] = 'layout.ket'

        # --- Molecules + Orbitals (HOMO/LUMO + Localized) ---
        print('  Computing molecule orbitals...')
        from mol_orbitals import compute_all_orbitals
        
        cube_dir = os.path.join(job_dir, 'orbitals_raw')
        os.makedirs(cube_dir, exist_ok=True)
        for mol_info in bundle.get('molecules', []):
            mol_data = _process_molecule_orbitals(mol_info, job_dir, cube_dir, zf)
            manifest['molecules'].append(mol_data)
            zf.writestr(f'molecules/{mol_data["key"]}.json', json.dumps(mol_data, indent=2))

        zf.writestr('manifest.json', json.dumps(manifest, indent=2))

    sz = os.path.getsize(output_path) / (1024*1024)
    print(f'\n  ✓ {output_path} ({sz:.1f} MB)')
    print(f'  ✓ {len(manifest["molecules"])} molecules')
    if manifest['has_trajectory']:
        print(f'  ✓ {manifest.get("n_frames","?")} frames')
    if manifest['has_orbitals']:
        print(f'  ✓ {manifest.get("n_orbital_cubes","?")} cubes')
    return output_path


def clean_intermediates(job_dir):
    """Remove intermediate files, keeping only: .rxnbundle.zip, .svg, .cdxml."""
    keep_ext = {'.zip', '.svg', '.cdxml', '.log'}
    keep_special = {'compute_', 'orbitals_raw'}  # keep compute scripts and raw cubes

    removed = 0
    for f in os.listdir(job_dir):
        fp = os.path.join(job_dir, f)
        _, ext = os.path.splitext(f)

        if ext in keep_ext:
            continue
        if any(f.startswith(p) for p in keep_special):
            continue

        if os.path.isdir(fp):
            shutil.rmtree(fp)
            print(f'    rm -rf {f}/')
            removed += 1
        else:
            os.remove(fp)
            removed += 1

    print(f'  ✓ Cleaned {removed} intermediate files/dirs')


def main():
    ap = argparse.ArgumentParser(description='Package reaction data into .rxnbundle.zip')
    ap.add_argument('job_dir')
    ap.add_argument('-o', '--output')
    ap.add_argument('--significant-only', action='store_true')
    ap.add_argument('--clean', action='store_true',
                    help='Delete intermediate files after export (keep .zip, .svg, .cdxml)')
    args = ap.parse_args()
    print(f"\n{'═'*50}\n Reaction Bundle Export\n{'═'*50}\n")
    r = export_bundle(args.job_dir, args.output, args.significant_only, args.clean)
    if r:
        print(f"\n{'═'*50}\n ✓ Bundle: {r}\n{'═'*50}\n")
        if args.clean:
            print('  Cleaning intermediates...')
            clean_intermediates(args.job_dir)
            # Also clean up unneeded cube files in mep_orbitals if they exist
            orb_dir = os.path.join(args.job_dir, 'mep_orbitals')
            if os.path.exists(orb_dir):
                for c in glob.glob(os.path.join(orb_dir, '*.cube')):
                    os.remove(c)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
