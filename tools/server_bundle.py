import json
import os
import sys
import zipfile

TOOLS_DIR = os.path.dirname(os.path.abspath(__file__))
OUTPUT_DIR = os.path.join(os.path.dirname(TOOLS_DIR), 'output')


def load_cube_converter():
    if TOOLS_DIR not in sys.path:
        sys.path.insert(0, TOOLS_DIR)
    try:
        from cube_converter import CubeConverter
        return CubeConverter
    except ImportError:
        return None


def is_bundle_metadata(filename):
    if filename == 'manifest.json' or filename.endswith('.xyz'):
        return True
    return filename.endswith('.json') and not filename.endswith('_pos.json') and not filename.endswith('_neg.json')


def purge_cube_dir(dir_path):
    if not os.path.exists(dir_path):
        return 0
    count = 0
    for f in os.listdir(dir_path):
        if not f.endswith('.cube'):
            continue
        try:
            os.remove(os.path.join(dir_path, f))
            count += 1
        except OSError:
            pass
    return count


def try_remesh(item, name, isovalue, converter, remeshed_meshes):
    if not (isovalue and converter):
        return None
    if not (item.filename.endswith('_pos.json') or item.filename.endswith('_neg.json')):
        return None
    if item.filename.endswith('esp_pos.json') or item.filename.endswith('esp_neg.json'):
        return None

    base_name = os.path.basename(item.filename).replace('_pos.json', '').replace('_neg.json', '')
    raw_cube = os.path.join(OUTPUT_DIR, name, 'orbitals_raw', f'{base_name}.cube')
    if not os.path.exists(raw_cube):
        raw_cube = os.path.join(OUTPUT_DIR, name, 'mep_orbitals', f'{base_name}.cube')
    if not os.path.exists(raw_cube):
        return None

    sign = 1 if item.filename.endswith('_pos.json') else -1
    cache_key = (raw_cube, sign * float(isovalue))

    if cache_key not in remeshed_meshes:
        try:
            vol, origin, step = converter.read_cube(raw_cube)
            remeshed_meshes[cache_key] = converter.extract_mesh(vol, origin, step, cache_key[1])
        except Exception as e:
            remeshed_meshes[cache_key] = None
            print(f"Error remeshing {raw_cube}: {e}")

    return remeshed_meshes.get(cache_key)


def package_custom_bundle(zip_path, custom_zip_path, keep_files, isovalue, name):
    converter = load_cube_converter() if isovalue else None
    remeshed_meshes = {}
    seen_files = set()

    with zipfile.ZipFile(zip_path, 'r') as zin, zipfile.ZipFile(custom_zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename in seen_files:
                continue
            seen_files.add(item.filename)

            keep_item = is_bundle_metadata(item.filename) or item.filename in keep_files
            if not keep_item:
                continue

            mesh = try_remesh(item, name, isovalue, converter, remeshed_meshes)
            if mesh:
                zout.writestr(item, json.dumps(mesh))
            else:
                zout.writestr(item, zin.read(item.filename))


def read_bundle_entry(zip_path, inner_path):
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            return z.read(inner_path), None
    except KeyError:
        if inner_path.endswith('_pos.json') or inner_path.endswith('_neg.json'):
            return b'{"vertices": [], "faces": []}', None
        return None, f'File not found in bundle: {inner_path}'
