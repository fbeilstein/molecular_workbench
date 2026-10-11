#!/usr/bin/env python3
"""
Molecular Workbench — Server

Serves the UI and dispatches API calls.

Workflow:
    1. Draw reactants and products
    2. Compute Reaction -> automatically maps atoms, generates 3D geometries, finds path, generates orbitals
    3. Verify Bundle → preview results in UI

Launch:
    python server.py
    python server.py --port 8080
"""

import http.server
import json
import os
import signal
import subprocess
import sys
import threading
import urllib.parse
import zipfile

# ── Path setup ───────────────────────────────────────────────────────────────

PACKAGE_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_DIR = PACKAGE_DIR
TOOLS_DIR = os.path.join(PACKAGE_DIR, 'tools')
OUTPUT_DIR = os.path.join(PACKAGE_DIR, 'output')
KETCHER_DIR = os.path.join(PROJECT_DIR, '.ketcher')
VENV_PYTHON = os.path.join(PROJECT_DIR, '.venv', 'bin', 'python')

if not os.path.exists(VENV_PYTHON):
    VENV_PYTHON = sys.executable

os.makedirs(OUTPUT_DIR, exist_ok=True)


# ── Ketcher download ────────────────────────────────────────────────────────

def _ensure_ketcher():
    """Download Ketcher standalone if not present."""
    if os.path.exists(os.path.join(KETCHER_DIR, 'index.html')):
        return
    print("Downloading Ketcher editor...")
    os.makedirs(KETCHER_DIR, exist_ok=True)
    url = "https://github.com/nicedouble/KetcherFinder/releases/download/v1.0/ketcher-standalone-2.27.0.tgz"
    import tarfile, urllib.request, io, shutil
    try:
        resp = urllib.request.urlopen(url, timeout=60)
        data = resp.read()
        with tarfile.open(fileobj=io.BytesIO(data), mode='r:gz') as tar:
            tar.extractall(KETCHER_DIR)
        pkg = os.path.join(KETCHER_DIR, 'package')
        if os.path.isdir(pkg):
            for item in os.listdir(pkg):
                shutil.move(os.path.join(pkg, item), KETCHER_DIR)
            os.rmdir(pkg)
        print("Ketcher downloaded.")
    except Exception as e:
        print(f"Warning: Could not download Ketcher: {e}")


# ── MIME types ───────────────────────────────────────────────────────────────

MIME = {
    '.html': 'text/html', '.css': 'text/css', '.js': 'application/javascript',
    '.json': 'application/json', '.svg': 'image/svg+xml', '.png': 'image/png',
    '.wasm': 'application/wasm', '.cube': 'text/plain', '.xyz': 'text/plain',
    '.sdf': 'chemical/x-mdl-sdfile', '.mol': 'chemical/x-mdl-molfile',
    '.cdxml': 'application/xml', '.sh': 'text/plain',
}


def _mime(path):
    ext = os.path.splitext(path)[1].lower()
    return MIME.get(ext, 'application/octet-stream')


# ── Utility Helpers ──────────────────────────────────────────────────────────

def _try_remesh(item, name, isovalue, converter, remeshed_meshes):
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


def _kill_job(job):
    if not hasattr(job, 'pid'):
        return
    try:
        os.killpg(os.getpgid(job.pid), signal.SIGTERM)
    except ProcessLookupError:
        pass


def _is_bundle_metadata(filename):
    if filename == 'manifest.json' or filename.endswith('.xyz'):
        return True
    return filename.endswith('.json') and not filename.endswith('_pos.json') and not filename.endswith('_neg.json')


def _purge_cube_dir(dir_path):
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


def _load_cube_converter():
    sys.path.insert(0, TOOLS_DIR)
    try:
        from cube_converter import CubeConverter
        return CubeConverter
    except ImportError:
        return None


def _execute_script_job(log_path, script_path, active_jobs, name):
    with open(log_path, 'w') as lf:
        p = subprocess.Popen(["bash", script_path], stdout=lf, stderr=subprocess.STDOUT, start_new_session=True)
        active_jobs[name] = p
        p.wait()
        active_jobs.pop(name, None)


def _package_custom_bundle(zip_path, custom_zip_path, keep_files, isovalue, name):
    converter = _load_cube_converter() if isovalue else None
    remeshed_meshes = {}
    seen_files = set()

    with zipfile.ZipFile(zip_path, 'r') as zin, zipfile.ZipFile(custom_zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zout:
        for item in zin.infolist():
            if item.filename in seen_files:
                continue
            seen_files.add(item.filename)

            keep_item = _is_bundle_metadata(item.filename) or item.filename in keep_files
            if not keep_item:
                continue

            mesh = _try_remesh(item, name, isovalue, converter, remeshed_meshes)
            if mesh:
                zout.writestr(item, json.dumps(mesh))
            else:
                zout.writestr(item, zin.read(item.filename))


def _read_bundle_entry(zip_path, inner_path):
    try:
        with zipfile.ZipFile(zip_path, 'r') as z:
            return z.read(inner_path), None
    except KeyError:
        if inner_path.endswith('_pos.json') or inner_path.endswith('_neg.json'):
            return b'{"vertices": [], "faces": []}', None
        return None, f'File not found in bundle: {inner_path}'


# ── Request handler ──────────────────────────────────────────────────────────

class WorkbenchHandler(http.server.BaseHTTPRequestHandler):

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        qs = urllib.parse.parse_qs(parsed.query)

        if path in ('/', '/index.html'):
            return self._serve_file(os.path.join(PACKAGE_DIR, 'templates', 'index.html'))
        if path.startswith('/static/'):
            return self._serve_file(os.path.join(PACKAGE_DIR, path.lstrip('/')))
        if path.startswith('/ketcher/'):
            return self._serve_file(os.path.join(KETCHER_DIR, path[len('/ketcher/'):]))
        if path.startswith('/api/bundle/'):
            return self._serve_bundle_file(path[len('/api/bundle/'):])
        if path.startswith('/output/'):
            return self._serve_file(os.path.join(OUTPUT_DIR, path[len('/output/'):]))
        if path == '/api/dirs':
            dirs = sorted(d for d in os.listdir(OUTPUT_DIR) if os.path.isdir(os.path.join(OUTPUT_DIR, d)))
            return self._json_response(dirs)
        if path == '/api/job/status':
            return self._api_job_status(qs)

        self.send_error(404)

    def do_POST(self):
        path = urllib.parse.urlparse(self.path).path
        body = self._read_json()

        routes = {
            '/api/rxn/compute-script': self._api_compute_script,
            '/api/job/stop': self._api_job_stop,
            '/api/bundle/package': self._api_bundle_package,
            '/api/orbital/remesh': self._api_orbital_remesh,
            '/api/bundle/purge-cubes': self._api_purge_cubes,
        }
        handler = routes.get(path)
        if handler:
            return handler(body)
        self.send_error(404)

    # ── API: Compute Script ───────────────────────────────────────────

    def _api_compute_script(self, body):
        sys.path.insert(0, TOOLS_DIR)
        from script_builder import build_compute_script

        script_path, name, job_dir = build_compute_script(body, OUTPUT_DIR, VENV_PYTHON, TOOLS_DIR)
        log_path = os.path.join(job_dir, 'compute.log')

        if not hasattr(self.server, 'active_jobs'):
            self.server.active_jobs = {}
        self.server.active_jobs[name] = "starting"

        t = threading.Thread(
            target=_execute_script_job,
            args=(log_path, script_path, self.server.active_jobs, name),
            daemon=True
        )
        t.start()

        self._json_response({'status': 'started', 'name': name})

    def _api_job_status(self, qs):
        name = qs.get('name', [''])[0]
        jobs = getattr(self.server, 'active_jobs', {})
        is_running = name in jobs
        job_dir = os.path.join(OUTPUT_DIR, name)
        log_path = os.path.join(job_dir, 'compute.log')
        bundle_path = os.path.join(job_dir, f"{name}.rxnbundle.zip")

        log_content = ""
        if os.path.exists(log_path):
            with open(log_path, 'r', errors='ignore') as f:
                log_content = f.read()

        self._json_response({
            'running': is_running,
            'log': log_content,
            'ready': os.path.exists(bundle_path)
        })

    def _api_job_stop(self, body):
        name = body.get('name', '')
        jobs = getattr(self.server, 'active_jobs', {})
        if name not in jobs:
            return self._json_response({'status': 'not_running'})

        _kill_job(jobs[name])
        jobs.pop(name, None)
        self._json_response({'status': 'stopped'})

    def _serve_bundle_file(self, path_in_bundle):
        parts = path_in_bundle.split('/', 1)
        if len(parts) < 2:
            return self.send_error(400)
        job_name, inner_path = parts
        zip_path = os.path.join(OUTPUT_DIR, job_name, f'{job_name}.rxnbundle.zip')
        if not os.path.exists(zip_path):
            return self.send_error(404, f'No bundle: {zip_path}')

        data, err = _read_bundle_entry(zip_path, inner_path)
        if err:
            return self.send_error(404, err)

        ext = os.path.splitext(inner_path)[1]
        ct = {'.json': 'application/json', '.gz': 'application/gzip',
              '.svg': 'image/svg+xml', '.cube': 'text/plain'}.get(ext, 'application/octet-stream')
        self.send_response(200)
        self.send_header('Content-Type', ct)
        self.send_header('Content-Length', str(len(data)))
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        self.end_headers()
        self.wfile.write(data)

    def _api_bundle_package(self, body):
        name = body.get('name')
        if not name:
            return self.send_error(400)

        zip_path = os.path.join(OUTPUT_DIR, name, f'{name}.rxnbundle.zip')
        if not os.path.exists(zip_path):
            return self.send_error(404, 'Bundle not found')

        custom_zip_path = os.path.join(OUTPUT_DIR, name, f'{name}_custom.rxnbundle.zip')
        _package_custom_bundle(zip_path, custom_zip_path, set(body.get('files', [])), body.get('isovalue'), name)
        self._json_response({'download_url': f'/output/{name}/{name}_custom.rxnbundle.zip'})

    def _api_orbital_remesh(self, body):
        name = body.get('name')
        cube_file_path = body.get('cube_file')
        isovalue = body.get('isovalue', 0.045)

        if not name or not cube_file_path:
            return self.send_error(400)

        base_name = os.path.basename(cube_file_path).replace('.json', '').replace('_pos', '').replace('_neg', '')
        raw_cube = os.path.join(OUTPUT_DIR, name, 'orbitals_raw', f'{base_name}.cube')
        if not os.path.exists(raw_cube):
            raw_cube = os.path.join(OUTPUT_DIR, name, 'mep_orbitals', f'{base_name}.cube')
            if not os.path.exists(raw_cube):
                return self._json_response({'error': f'Cube file not found: {base_name}.cube'})

        converter = _load_cube_converter()
        if not converter:
            return self._json_response({'error': 'CubeConverter tool unavailable'})

        try:
            vol, origin, step = converter.read_cube(raw_cube)
            pos = converter.extract_mesh(vol, origin, step, isovalue)
            neg = converter.extract_mesh(vol, origin, step, -isovalue)
            self._json_response({'pos': pos, 'neg': neg})
        except Exception as e:
            self._json_response({'error': str(e)})

    def _api_purge_cubes(self, body):
        name = body.get('name')
        if not name:
            return self.send_error(400)

        count = _purge_cube_dir(os.path.join(OUTPUT_DIR, name, 'orbitals_raw'))
        count += _purge_cube_dir(os.path.join(OUTPUT_DIR, name, 'mep_orbitals'))
        self._json_response({'status': 'ok', 'purged': count})

    # ── Utility methods ───────────────────────────────────────────────

    def _serve_file(self, filepath):
        if not os.path.isfile(filepath):
            return self.send_error(404)
        self.send_response(200)
        self.send_header('Content-Type', _mime(filepath))
        self.send_header('Content-Length', os.path.getsize(filepath))
        self.send_header('Access-Control-Allow-Origin', '*')
        self.end_headers()
        with open(filepath, 'rb') as f:
            self.wfile.write(f.read())

    def _read_json(self):
        length = int(self.headers.get('Content-Length', 0))
        if length == 0:
            return {}
        return json.loads(self.rfile.read(length))

    def _json_response(self, obj, status=200):
        data = json.dumps(obj).encode()
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.send_header('Content-Length', len(data))
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
        self.end_headers()
        self.wfile.write(data)

    def log_message(self, fmt, *args):
        pass  # Quiet


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Molecular Workbench')
    parser.add_argument('--port', type=int, default=8080)
    parser.add_argument('--host', default='0.0.0.0')
    args = parser.parse_args()

    _ensure_ketcher()

    class ReusableServer(http.server.HTTPServer):
        allow_reuse_address = True

    try:
        server = ReusableServer((args.host, args.port), WorkbenchHandler)
    except Exception as e:
        print('Failed to bind server:', e)
        raise

    print(f"  🧪 Molecular Workbench")
    print(f"  → http://localhost:{args.port}")
    print(f"  Tools:  {TOOLS_DIR}")
    print(f"  Output: {OUTPUT_DIR}")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nShutting down.")
        server.shutdown()


if __name__ == '__main__':
    main()
