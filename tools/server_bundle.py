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
    except (ImportError, ModuleNotFoundError):
        import glob
        v_sites = glob.glob(os.path.join(os.path.dirname(TOOLS_DIR), '.venv', 'lib', 'python*', 'site-packages'))
        if v_sites and v_sites[0] not in sys.path:
            sys.path.insert(0, v_sites[0])
        try:
            from cube_converter import CubeConverter
            return CubeConverter
        except Exception:
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


def locate_raw_cube(name, base_name):
    import urllib.parse
    for n in [name, urllib.parse.unquote(str(name))]:
        for sub in ['orbitals_raw', 'mep_orbitals', '']:
            d = os.path.join(OUTPUT_DIR, n, sub) if sub else os.path.join(OUTPUT_DIR, n)
            if not os.path.isdir(d):
                continue
            direct = os.path.join(d, f'{base_name}.cube')
            if os.path.exists(direct):
                return direct
            if '_' in base_name:
                prefix, s_idx = base_name.rsplit('_', 1)
                cands = sorted(
                    [c for c in os.listdir(d) if c.startswith(prefix + '_') and c.endswith('.cube')],
                    key=lambda x: int(x.rsplit('_', 1)[1].split('.')[0]) if x.rsplit('_', 1)[1].split('.')[0].isdigit() else 0
                )
                if len(cands) == 1:
                    return os.path.join(d, cands[0])
                elif len(cands) > 1 and s_idx.isdigit():
                    num = int(s_idx)
                    for c in cands:
                        c_num = int(c.rsplit('_', 1)[1].split('.')[0])
                        if c_num in (num, num + 1, num + 2, num - 1):
                            return os.path.join(d, c)
    return None


def try_remesh(item, name, isovalue, converter, remeshed_meshes):
    if not (isovalue and converter):
        return None
    if not (item.filename.endswith('_pos.json') or item.filename.endswith('_neg.json')):
        return None
    if item.filename.endswith('esp_pos.json') or item.filename.endswith('esp_neg.json'):
        return None

    base_name = os.path.basename(item.filename).replace('_pos.json', '').replace('_neg.json', '')
    raw_cube = locate_raw_cube(name, base_name)

    if not raw_cube:
        print(f"[try_remesh] Raw cube not found for {base_name} in {name}")
        return None

    sign = 1 if item.filename.endswith('_pos.json') else -1
    cache_key = (raw_cube, sign * float(isovalue))

    if cache_key not in remeshed_meshes:
        try:
            vol, origin, step = converter.read_cube(raw_cube)
            mesh = converter.extract_mesh(vol, origin, step, cache_key[1])
            remeshed_meshes[cache_key] = mesh if mesh is not None else {'vertices': [], 'faces': [], 'normals': []}
        except Exception as e:
            remeshed_meshes[cache_key] = {'vertices': [], 'faces': [], 'normals': []}
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

            fn = item.filename
            fn_base = fn.replace('molecules/', '')
            matched_keep = False
            if keep_files is not None:
                matched_keep = (fn in keep_files) or (fn_base in keep_files) or (f'molecules/{fn}' in keep_files)

            keep_item = (keep_files is None) or is_bundle_metadata(item.filename) or matched_keep
            if not keep_item:
                continue

            mesh = try_remesh(item, name, isovalue, converter, remeshed_meshes)
            if mesh is not None:
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
