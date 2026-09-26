import os
import json
import zipfile
import tempfile
import argparse
import shutil
import re

def load_json_from_zip(zf, filename):
    try:
        with zf.open(filename) as f:
            return json.load(f)
    except KeyError:
        return None

def write_json_to_zip(zf, filename, data):
    zf.writestr(filename, json.dumps(data, indent=2))

def parse_frame_idx(filename):
    # e.g., orbitals/sn1_2_f005_ibo_1_pos.json
    m = re.search(r'_f(\d+)_', filename)
    if m:
        return int(m.group(1))
    return None

def rename_frame_idx(filename, new_idx):
    return re.sub(r'_f\d+_', f'_f{new_idx:03d}_', filename)

def trim_bundle(in_zip, out_zip, start_idx, end_idx):
    print(f"Trimming {in_zip} from frame {start_idx} to {end_idx}...")
    with zipfile.ZipFile(in_zip, 'r') as z_in, zipfile.ZipFile(out_zip, 'w', compression=zipfile.ZIP_DEFLATED) as z_out:
        manifest = load_json_from_zip(z_in, 'manifest.json')
        if not manifest:
            raise ValueError("Invalid bundle: missing manifest.json")
            
        n_frames = manifest.get('n_frames', 0)
        if end_idx >= n_frames:
            end_idx = n_frames - 1
            
        keep_count = end_idx - start_idx + 1
        print(f"  Keeping {keep_count} frames")
        
        manifest['n_frames'] = keep_count
        
        frames = load_json_from_zip(z_in, 'trajectory/frames.json')
        if frames:
            frames = frames[start_idx:end_idx+1]
            write_json_to_zip(z_out, 'trajectory/frames.json', frames)
            
        orbitals = load_json_from_zip(z_in, 'trajectory/orbitals.json')
        if orbitals and 'frames' in orbitals:
            orb_frames = orbitals['frames'][start_idx:end_idx+1]
            for i, f in enumerate(orb_frames):
                f['frame_index'] = i
                # Also update internal file references
                for k, v in f.items():
                    if isinstance(v, dict) and 'file' in v:
                        v['file'] = rename_frame_idx(v['file'], i)
            orbitals['frames'] = orb_frames
            orbitals['n_frames'] = keep_count
            write_json_to_zip(z_out, 'trajectory/orbitals.json', orbitals)
            
        for item in z_in.infolist():
            if item.filename in ['trajectory/frames.json', 'trajectory/orbitals.json', 'manifest.json']:
                continue
                
            if item.filename.startswith('orbitals/'):
                f_idx = parse_frame_idx(item.filename)
                if f_idx is not None:
                    if start_idx <= f_idx <= end_idx:
                        new_name = rename_frame_idx(item.filename, f_idx - start_idx)
                        z_out.writestr(new_name, z_in.read(item.filename))
                    continue
                    
            z_out.writestr(item, z_in.read(item.filename))
            
        write_json_to_zip(z_out, 'manifest.json', manifest)
    print(f"Saved trimmed bundle to {out_zip}")

def merge_bundles(in_zips, out_zip):
    print(f"Merging {len(in_zips)} bundles...")
    
    merged_manifest = None
    merged_frames = []
    merged_orb_frames = []
    
    orb_keys = []
    ibos_meta = []
    
    frame_offset = 0
    
    with zipfile.ZipFile(out_zip, 'w', compression=zipfile.ZIP_DEFLATED) as z_out:
        for zip_idx, in_zip in enumerate(in_zips):
            print(f"  Processing {in_zip} (offset={frame_offset})")
            with zipfile.ZipFile(in_zip, 'r') as z_in:
                manifest = load_json_from_zip(z_in, 'manifest.json')
                n_frames = manifest.get('n_frames', 0)
                
                if zip_idx == 0:
                    merged_manifest = manifest.copy()
                    # Copy non-trajectory files from the first bundle
                    for item in z_in.infolist():
                        if not item.filename.startswith('orbitals/') and not item.filename.startswith('trajectory/') and item.filename != 'manifest.json':
                            z_out.writestr(item, z_in.read(item.filename))
                
                frames = load_json_from_zip(z_in, 'trajectory/frames.json')
                if frames:
                    merged_frames.extend(frames)
                    
                orbitals = load_json_from_zip(z_in, 'trajectory/orbitals.json')
                if orbitals:
                    for k, m in zip(orbitals.get('orbital_keys', []), orbitals.get('ibos', [])):
                        if k not in orb_keys:
                            orb_keys.append(k)
                            ibos_meta.append(m)
                            
                    if 'frames' in orbitals:
                        for f in orbitals['frames']:
                            f['frame_index'] += frame_offset
                            for k, v in f.items():
                                if isinstance(v, dict) and 'file' in v:
                                    v['file'] = rename_frame_idx(v['file'], f['frame_index'])
                        merged_orb_frames.extend(orbitals['frames'])
                
                for item in z_in.infolist():
                    if item.filename.startswith('orbitals/'):
                        f_idx = parse_frame_idx(item.filename)
                        if f_idx is not None:
                            new_name = rename_frame_idx(item.filename, f_idx + frame_offset)
                            z_out.writestr(new_name, z_in.read(item.filename))
                            
            frame_offset += n_frames
            
        merged_manifest['n_frames'] = frame_offset
        merged_manifest['name'] = merged_manifest['name'] + "_merged"
        write_json_to_zip(z_out, 'manifest.json', merged_manifest)
        
        if merged_frames:
            write_json_to_zip(z_out, 'trajectory/frames.json', merged_frames)
            
        if merged_orb_frames:
            base_orb = load_json_from_zip(zipfile.ZipFile(in_zips[0], 'r'), 'trajectory/orbitals.json')
            if not base_orb:
                base_orb = {}
            base_orb['n_frames'] = frame_offset
            base_orb['frames'] = merged_orb_frames
            base_orb['orbital_keys'] = orb_keys
            base_orb['ibos'] = ibos_meta
            write_json_to_zip(z_out, 'trajectory/orbitals.json', base_orb)

    print(f"Saved merged bundle ({frame_offset} total frames) to {out_zip}")

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Reaction Trajectory Editor')
    subparsers = parser.add_subparsers(dest='command', required=True)
    
    p_trim = subparsers.add_parser('trim')
    p_trim.add_argument('input_zip')
    p_trim.add_argument('-o', '--output', required=True)
    p_trim.add_argument('--start', type=int, required=True)
    p_trim.add_argument('--end', type=int, required=True)
    
    p_merge = subparsers.add_parser('merge')
    p_merge.add_argument('input_zips', nargs='+')
    p_merge.add_argument('-o', '--output', required=True)
    
    args = parser.parse_args()
    
    if args.command == 'trim':
        trim_bundle(args.input_zip, args.output, args.start, args.end)
    elif args.command == 'merge':
        merge_bundles(args.input_zips, args.output)
