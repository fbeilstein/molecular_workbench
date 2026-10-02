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
    python server.py --port 8081
"""

import http.server
import json
import os
import sys
import urllib.parse
import time

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


# ── Request handler ──────────────────────────────────────────────────────────

def _try_remesh(item, name, isovalue, converter, remeshed_meshes):
    if not (isovalue and converter): return None
    if not (item.filename.endswith('_pos.json') or item.filename.endswith('_neg.json')): return None
    if item.filename.endswith('esp_pos.json') or item.filename.endswith('esp_neg.json'): return None
    
    base_name = __import__('os').path.basename(item.filename).replace('_pos.json', '').replace('_neg.json', '')
    raw_cube = __import__('os').path.join(OUTPUT_DIR, name, 'orbitals_raw', f'{base_name}.cube')
    if not __import__('os').path.exists(raw_cube):
        raw_cube = __import__('os').path.join(OUTPUT_DIR, name, 'mep_orbitals', f'{base_name}.cube')
    if not __import__('os').path.exists(raw_cube):
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

class WorkbenchHandler(http.server.BaseHTTPRequestHandler):

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        qs = urllib.parse.parse_qs(parsed.query)

        if path == '/' or path == '/index.html':
            self._serve_file(os.path.join(PACKAGE_DIR, 'templates', 'index.html'))

        elif path.startswith('/static/'):
            self._serve_file(os.path.join(PACKAGE_DIR, path.lstrip('/')))

        elif path.startswith('/ketcher/'):
            rel = path[len('/ketcher/'):]
            self._serve_file(os.path.join(KETCHER_DIR, rel))

        elif path.startswith('/api/bundle/'):
            self._serve_bundle_file(path[len('/api/bundle/'):])

        elif path.startswith('/output/'):
            rel = path[len('/output/'):]
            self._serve_file(os.path.join(OUTPUT_DIR, rel))

        elif path == '/api/dirs':
            dirs = sorted(d for d in os.listdir(OUTPUT_DIR)
                          if os.path.isdir(os.path.join(OUTPUT_DIR, d)))
            self._json_response(dirs)

        elif path == '/api/job/status':
            self._api_job_status(qs)

        else:
            self.send_error(404)

    def do_POST(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        body = self._read_json()

        if path == '/api/rxn/compute-script':
            self._api_compute_script(body)
        elif path == '/api/job/stop':
            self._api_job_stop(body)
        elif path == '/api/bundle/package':
            self._api_bundle_package(body)
        elif path == '/api/orbital/remesh':
            self._api_orbital_remesh(body)
        elif path == '/api/bundle/purge-cubes':
            self._api_purge_cubes(body)
        elif path == '/api/bundle/purge-cubes':
            self._api_purge_cubes(body)
            self._api_orbital_remesh(body)
        else:
            self.send_error(404)

    # ── API: Compute Script ───────────────────────────────────────────

    def _api_compute_script(self, body):
        import sys
        sys.path.insert(0, TOOLS_DIR)
        from script_builder import build_compute_script
        
        script_path, name, job_dir = build_compute_script(body, OUTPUT_DIR, VENV_PYTHON, TOOLS_DIR)

        log_path = os.path.join(job_dir, 'compute.log')
        import subprocess, threading
        if not hasattr(self.server, 'active_jobs'):
            self.server.active_jobs = {}
            
        self.server.active_jobs[name] = "starting"
        
        def run_job():
            with open(log_path, 'w') as lf:
                p = subprocess.Popen(["bash", script_path], stdout=lf, stderr=subprocess.STDOUT)
                self.server.active_jobs[name] = p
                p.wait()
                self.server.active_jobs.pop(name, None)
                
        threading.Thread(target=run_job).start()

        self._json_response({
            'status': 'started',
            'name': name
        })

    def _api_job_status(self, query):
        name = query.get('name', [''])[0]
        job_dir = os.path.join(OUTPUT_DIR, name)
        log_path = os.path.join(job_dir, 'compute.log')
        
        is_running = name in getattr(self.server, 'active_jobs', {})
        log_content = ""
        if os.path.exists(log_path):
            with open(log_path, 'r') as f:
                lines = f.readlines()
                log_content = ''.join(lines[-50:])
                
        bundle_ready = os.path.exists(os.path.join(job_dir, f'{name}.rxnbundle.zip'))
        
        self._json_response({
            'running': is_running,
            'log': log_content,
            'ready': bundle_ready
        })

    def _api_job_stop(self, body):
        name = body.get('name', '')
        jobs = getattr(self.server, 'active_jobs', {})
        if name in jobs:
            jobs[name].terminate()
            jobs.pop(name, None)
            self._json_response({'status': 'stopped'})
        else:
            self._json_response({'status': 'not_running'})

    def _serve_bundle_file(self, path_in_bundle):
        """Serve a file from inside a .rxnbundle.zip.
        URL: /api/bundle/<job_name>/<path_inside_zip>
        """
        parts = path_in_bundle.split('/', 1)
        if len(parts) < 2:
            self.send_error(400)
            return
        job_name, inner_path = parts
        zip_path = os.path.join(OUTPUT_DIR, job_name, f'{job_name}.rxnbundle.zip')
        if not os.path.exists(zip_path):
            self.send_error(404, f'No bundle: {zip_path}')
            return
        try:
            import zipfile
            with zipfile.ZipFile(zip_path, 'r') as z:
                data = z.read(inner_path)
            ext = os.path.splitext(inner_path)[1]
            ct = {'.json': 'application/json', '.gz': 'application/gzip',
                  '.svg': 'image/svg+xml', '.cube': 'text/plain'}.get(ext, 'application/octet-stream')
            self.send_response(200)
            self.send_header('Content-Type', ct)
            self.send_header('Content-Length', str(len(data)))
            self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
            self.end_headers()
            self.wfile.write(data)
        except KeyError:
            if inner_path.endswith('_pos.json') or inner_path.endswith('_neg.json'):
                data = b'{"vertices": [], "faces": []}'
                self.send_response(200)
                self.send_header('Content-Type', 'application/json')
                self.send_header('Content-Length', str(len(data)))
                self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
                self.end_headers()
                self.wfile.write(data)
            else:
                self.send_error(404, f'Not in zip: {inner_path}')
        except Exception as e:
            self.send_error(500, str(e))


    def _api_bundle_package(self, body):
        import zipfile, json
        name = body.get('name')
        keep_files = set(body.get('keep_files', []))
        isovalue = body.get('isovalue')
        
        zip_path = os.path.join(OUTPUT_DIR, name, f'{name}.rxnbundle.zip')
        if not os.path.exists(zip_path):
            self._json_response({'error': 'Bundle not found'})
            return
            
        custom_zip_path = os.path.join(OUTPUT_DIR, name, f'{name}_custom.rxnbundle.zip')
        
        # If we need to remesh, import converter
        converter = None
        if isovalue:
            import sys
            sys.path.insert(0, TOOLS_DIR)
            try:
                from cube_converter import CubeConverter
                converter = CubeConverter
            except ImportError:
                pass
                
        # Remesh cache
        remeshed_meshes = {}
        seen_files = set()
        
        with zipfile.ZipFile(zip_path, 'r') as zin, zipfile.ZipFile(custom_zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename in seen_files:
                    continue
                seen_files.add(item.filename)
                
                is_metadata = item.filename == 'manifest.json' or item.filename.endswith('.xyz')
                if item.filename.endswith('.json') and not item.filename.endswith('_pos.json') and not item.filename.endswith('_neg.json'):
                    is_metadata = True
                
                keep_item = is_metadata or item.filename in keep_files
                
                if not keep_item:
                    continue
                
                mesh = _try_remesh(item, name, isovalue, converter, remeshed_meshes)
                if mesh:
                    zout.writestr(item, json.dumps(mesh))
                else:
                    zout.writestr(item, zin.read(item.filename))
        
        self._json_response({'download_url': f'/output/{name}/{name}_custom.rxnbundle.zip'})

    def _api_orbital_remesh(self, body):
        name = body.get('name')
        cube_file_path = body.get('cube_file')  # e.g. "molecules/[18]annulene_homo.json"
        isovalue = body.get('isovalue', 0.045)
        
        if not name or not cube_file_path:
            self.send_error(400)
            return
            
        base_name = os.path.basename(cube_file_path).replace('.json', '').replace('_pos', '').replace('_neg', '')
        
        raw_cube = os.path.join(OUTPUT_DIR, name, 'orbitals_raw', f'{base_name}.cube')
        
        if not os.path.exists(raw_cube):
            raw_cube = os.path.join(OUTPUT_DIR, name, 'mep_orbitals', f'{base_name}.cube')
            if not os.path.exists(raw_cube):
                self._json_response({'error': f'Cube file not found: {base_name}.cube'})
                return
                
        sys.path.insert(0, TOOLS_DIR)
        from cube_converter import CubeConverter
        
        try:
            vol, origin, step = CubeConverter.read_cube(raw_cube)
            pos = CubeConverter.extract_mesh(vol, origin, step, isovalue)
            neg = CubeConverter.extract_mesh(vol, origin, step, -isovalue)
            
            self._json_response({'pos': pos, 'neg': neg})
        except Exception as e:
            self._json_response({'error': str(e)})

    def _api_purge_cubes(self, body):
        import shutil
        name = body.get('name')
        if not name:
            self.send_error(400)
            return
            
        raw_dir = os.path.join(OUTPUT_DIR, name, 'orbitals_raw')
        count = 0
        if os.path.exists(raw_dir):
            for f in os.listdir(raw_dir):
                if f.endswith('.cube'):
                    os.remove(os.path.join(raw_dir, f))
                    count += 1
                    
        mep_dir = os.path.join(OUTPUT_DIR, name, 'mep_orbitals')
        if os.path.exists(mep_dir):
            for f in os.listdir(mep_dir):
                if f.endswith('.cube'):
                    os.remove(os.path.join(mep_dir, f))
                    count += 1
                    
        self._json_response({'status': 'ok', 'purged': count})

    # ── Utility methods ───────────────────────────────────────────────

    def _serve_file(self, filepath):
        if not os.path.isfile(filepath):
            self.send_error(404)
            return
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
    print(f'Attempting to bind to {args.host}:{args.port}')
    import os
    os.system(f'lsof -i :{args.port} || echo "Nothing on {args.port}"')
    try:
        server = ReusableServer((args.host, args.port), WorkbenchHandler)
    except Exception as e:
        print('Failed! ', e)
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

