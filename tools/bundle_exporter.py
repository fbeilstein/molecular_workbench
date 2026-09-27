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
        'charge': bundle.get('charge', 0),
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

        # --- Molecules + Orbitals (HOMO/LUMO + Localized) ---
        print('  Computing molecule orbitals...')
        from mol_orbitals import compute_all_orbitals
        import tempfile
        for mol_info in bundle.get('molecules', []):
            key, role = mol_info['key'], mol_info['role']
            smiles = mol_info.get('smiles', '')
            mol_data = {'key': key, 'role': role, 'smiles': smiles, 'orbitals': []}

            xyz_path = os.path.join(job_dir, mol_info.get('xyz', ''))
            if os.path.exists(xyz_path):
                mol_data['xyz'] = load_text(xyz_path)

            if mol_data.get('xyz'):
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
                
                with tempfile.TemporaryDirectory() as tmp_dir:
                    try:
                        orb_manifest = compute_all_orbitals(
                            xyz_file=xyz_path, name=key, out_dir=tmp_dir, 
                            charge=mol_charge, method='b3lyp', smiles=smiles
                        )
                        
                        # Process cubes into .json surfaces
                        for cube_file in glob.glob(os.path.join(tmp_dir, '*.cube')):
                            CubeConverter.process_cube(cube_file, adaptive=True)
                            os.remove(cube_file)
                        
                        # Add Frontier Orbitals (HOMO + LUMO only)
                        frontier_group = {'name': 'Frontier', 'items': []}
                        for lbl in ['homo', 'lumo']:
                            if lbl in orb_manifest['canonical']:
                                info = orb_manifest['canonical'][lbl]
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                frontier_group['items'].append({
                                    'label': lbl.upper(),
                                    'file': f'molecules/{base_file}',
                                    'energy_ev': info['energy_ev'],
                                    'color': '#ff5c5c' if lbl == 'homo' else '#55aaff'
                                })
                        if frontier_group['items']:
                            mol_data['orbitals'].append(frontier_group)

                        # Add canonical π System (auto-detected from ring planes)
                        has_pi_system = False
                        pi_system = orb_manifest['canonical'].get('pi_system', [])
                        if pi_system:
                            has_pi_system = True
                            pi_group = {'name': 'π System (Canonical)', 'items': []}
                            PI_COLORS = ['#e040fb', '#ab47bc', '#7b1fa2', '#4a148c', '#ea80fc']
                            for ci, info in enumerate(pi_system):
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                pi_group['items'].append({
                                    'label': f"π ({info['homo_label']}, {info['energy_ev']} eV)",
                                    'file': f'molecules/{base_file}',
                                    'energy_ev': info['energy_ev'],
                                    'color': PI_COLORS[ci % len(PI_COLORS)]
                                })
                            mol_data['orbitals'].append(pi_group)

                        # Add Localized Orbitals
                        loc = orb_manifest['localized']
                        
                        def _sort_bonds(bonds):
                            return sorted(bonds, key=lambda x: (min(x['atoms']), max(x['atoms'])))
                            
                        if 'sigma' in loc: loc['sigma'] = _sort_bonds(loc['sigma'])
                        if 'sigma_star' in loc: loc['sigma_star'] = _sort_bonds(loc['sigma_star'])
                        if 'pi' in loc: loc['pi'] = _sort_bonds(loc['pi'])
                        if 'pi_star' in loc: loc['pi_star'] = _sort_bonds(loc['pi_star'])
                        
                        if loc.get('sigma'):
                            sigma_group = {'name': 'σ Bonds', 'items': []}
                            for info in loc['sigma']:
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                if len(info['atoms']) > 2:
                                    lbl = "deloc-σ(" + ",".join(info['atoms']) + ")"
                                else:
                                    lbl = f"σ({info['atoms'][0]}–{info['atoms'][1]})"
                                    
                                sigma_group['items'].append({
                                    'label': lbl,
                                    'file': f'molecules/{base_file}',
                                    'color': '#44cc77'
                                })
                            mol_data['orbitals'].append(sigma_group)
                            
                        # Show localized π bonds (now with multi-center support)
                        if loc.get('pi'):
                            pi_group = {'name': 'π Bonds', 'items': []}
                            for info in loc['pi']:
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                
                                if 'canonical_label' in info:
                                    lbl = info['canonical_label']
                                elif len(info['atoms']) > 2:
                                    lbl = "deloc-π(" + ",".join(info['atoms']) + ")"
                                else:
                                    lbl = f"π({info['atoms'][0]}={info['atoms'][1]})"
                                    
                                pi_group['items'].append({
                                    'label': lbl,
                                    'file': f'molecules/{base_file}',
                                    'color': '#cc44bb'
                                })
                            mol_data['orbitals'].append(pi_group)
                        if loc.get('sigma_star'):
                            sigma_star_group = {'name': 'σ* Anti-Bonds', 'items': []}
                            for info in loc['sigma_star']:
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                if len(info['atoms']) > 2:
                                    lbl = "σ*(" + ",".join(info['atoms']) + ")"
                                else:
                                    lbl = f"σ*({info['atoms'][0]}–{info['atoms'][1]})"
                                sigma_star_group['items'].append({
                                    'label': lbl,
                                    'file': f'molecules/{base_file}',
                                    'color': '#ffaa00'
                                })
                            mol_data['orbitals'].append(sigma_star_group)

                        if loc.get('pi_star'):
                            pi_star_group = {'name': 'π* Anti-Bonds', 'items': []}
                            for info in loc['pi_star']:
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                if 'canonical_label' in info:
                                    lbl = info['canonical_label']
                                elif len(info['atoms']) > 2:
                                    lbl = "π*(" + ",".join(info['atoms']) + ")"
                                else:
                                    lbl = f"π*({info['atoms'][0]}–{info['atoms'][1]})"
                                pi_star_group['items'].append({
                                    'label': lbl,
                                    'file': f'molecules/{base_file}',
                                    'color': '#ff00aa'
                                })
                            mol_data['orbitals'].append(pi_star_group)
                        if loc.get('lone_pairs'):
                            lp_group = {'name': 'Lone Pairs', 'items': []}
                            for info in loc['lone_pairs']:
                                base_file = info['file'].replace('.cube', '.json')
                                for sign in ['pos', 'neg']:
                                    sign_file = info['file'].replace('.cube', f'_{sign}.json')
                                    if os.path.exists(os.path.join(tmp_dir, sign_file)):
                                        zf.write(os.path.join(tmp_dir, sign_file), f'molecules/{sign_file}')
                                lp_group['items'].append({
                                    'label': f"LP({info['atom']} #{info['index']})",
                                    'file': f'molecules/{base_file}',
                                    'color': '#eeee22'
                                })
                            mol_data['orbitals'].append(lp_group)

                        # Add ESP surface if available
                        esp_info = orb_manifest.get('esp_surface')
                        if esp_info:
                            esp_file = esp_info['file']
                            esp_path = os.path.join(tmp_dir, esp_file)
                            if os.path.exists(esp_path):
                                zf.write(esp_path, f'molecules/{esp_file}')
                                mol_data['esp_surface'] = {
                                    'file': f'molecules/{esp_file}',
                                    'esp_min': esp_info['esp_min'],
                                    'esp_max': esp_info['esp_max'],
                                }
                            
                        print(f"    ✓ {key} ({smiles}): {sum(len(g['items']) for g in mol_data['orbitals'])} orbitals")
                    except Exception as e:
                        print(f"    ✗ {key} ({smiles}): computation failed ({e})")

            manifest['molecules'].append(mol_data)
            zf.writestr(f'molecules/{key}.json', json.dumps(mol_data, indent=2))

        # --- Trajectory ---
        trj_info = bundle.get('trajectory')
        if trj_info:
            trj_path = os.path.join(job_dir, trj_info.get('xyz', ''))
            if os.path.exists(trj_path):
                frames = parse_xyz_frames(load_text(trj_path))
                zf.writestr('trajectory/frames.json', json.dumps(frames, indent=1))
                manifest['has_trajectory'] = True
                manifest['n_frames'] = len(frames)
                manifest['n_atoms'] = len(frames[0]['atoms']) if frames else 0
            ep_path = os.path.join(job_dir, trj_info.get('energy_profile', ''))
            if os.path.exists(ep_path):
                zf.write(ep_path, 'trajectory/energy.json')

        # --- IBO Orbitals (raw cube, no gzip) ---
        orb_path = os.path.join(job_dir, f'{name}_mep_orbitals.json')
        if os.path.exists(orb_path):
            orb_manifest = load_json(orb_path)
            manifest['has_orbitals'] = True
            ibos = orb_manifest.get('ibos', [])
            keys = orb_manifest.get('orbital_keys', [])
            if significant_only:
                filtered_keys, filtered_ibos = [], []
                for k, m in zip(keys, ibos):
                    if m.get('significant'):
                        filtered_keys.append(k)
                        filtered_ibos.append(m)
                keys, ibos = filtered_keys, filtered_ibos
            orb_export = {
                'mode': orb_manifest.get('mode'),
                'method': orb_manifest.get('method'),
                'basis': orb_manifest.get('basis'),
                'n_frames': orb_manifest.get('n_frames'),
                'ibos': ibos, 'orbital_keys': keys, 'frames': [],
            }
            orb_dir = os.path.join(job_dir, 'mep_orbitals')
            n_cubes = 0
            for frame in orb_manifest.get('frames', []):
                fe = {'frame_index': frame['frame_index']}
                for key in keys:
                    info = frame.get(key)
                    if not info:
                        continue
                    cp = os.path.join(orb_dir, info['file'])
                    if os.path.exists(cp):
                        # Convert to surface — adaptive threshold for consistent volume across frames
                        CubeConverter.process_cube(cp, adaptive=True)
                        
                        base_file = info['file'].replace('.cube', '.json')
                        for sign in ['pos', 'neg']:
                            sign_file = info['file'].replace('.cube', f'_{sign}.json')
                            cp_sign = os.path.join(orb_dir, sign_file)
                            if os.path.exists(cp_sign):
                                zf.write(cp_sign, f"orbitals/{sign_file}")
                                if clean: os.remove(cp_sign)  # Cleanup the JSON file
                        
                        fe[key] = {'file': f"orbitals/{base_file}"}
                        n_cubes += 1
                        if clean: os.remove(cp) # cleanup the original cube
                orb_export['frames'].append(fe)
            zf.writestr('trajectory/orbitals.json', json.dumps(orb_export, indent=1))
            manifest['n_orbital_cubes'] = n_cubes

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
    keep_special = {'compute_'}  # keep compute scripts

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
