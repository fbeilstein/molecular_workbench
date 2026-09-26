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


def _detect_charge_from_smiles(smiles):
    """Detect total charge from SMILES string by calculating formal charges."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            mol.UpdatePropertyCache(strict=False)
            return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    except ImportError:
        pass
        
    import re
    total = 0
    # Match patterns like [C+], [NH4+], [Br-], [O-2], [Fe+3], [O+:6]
    for match in re.finditer(r'\[([^\]]*?)([+-])(\d*)(?::\d+)?\]', smiles):
        sign = 1 if match.group(2) == '+' else -1
        magnitude = int(match.group(3)) if match.group(3) else 1
        total += sign * magnitude
    return total

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
        else:
            self.send_error(404)

    # ── API: Compute Script ───────────────────────────────────────────

    def _api_compute_script(self, body):
        rxn_smiles = body.get('smiles', '')
        name = body.get('name', 'molecule')
        charge = body.get('charge', 0)

        smiles_charge = _detect_charge_from_smiles(rxn_smiles)
        if charge == 0 and smiles_charge != 0:
            charge = smiles_charge

        job_dir = os.path.join(OUTPUT_DIR, name)
        os.makedirs(job_dir, exist_ok=True)

        ket_data = body.get('ket', '')
        if ket_data:
            ket_path = os.path.join(job_dir, f"{name}.ket")
            with open(ket_path, 'w') as f:
                f.write(ket_data)

        svg_data = body.get('svg', '')
        if svg_data:
            svg_path = os.path.join(job_dir, f"{name}.svg")
            with open(svg_path, 'w') as f:
                f.write(svg_data)

        script_path = os.path.join(job_dir, f'compute_{name}.sh')

        lines = [
            '#!/bin/bash',
            'set -e',
            'export PYTHONUNBUFFERED=1',
            f'PYTHON="{VENV_PYTHON}"',
            f'TOOLS="{TOOLS_DIR}"',
            f'OUT="{job_dir}"',
            f'echo "═══ Computing: {name} ═══"',
            '',
            '# Step 1: Prepare molecule (3D structure, SVG, MOL)',
            f'$PYTHON $TOOLS/mol_prep.py "{rxn_smiles}" --name {name} -o $OUT --charge {charge}',
            f'obabel $OUT/{name}.mol -O $OUT/{name}.cdxml 2>/dev/null || true',
            '',
            '# Step 2: Create bundle manifest',
            f'cat <<EOF > $OUT/{name}_bundle.json',
            '{',
            f'  "name": "{name}",',
            f'  "smiles": "{rxn_smiles}",',
            f'  "charge": {charge},',
            f'  "method": "b3lyp",',
            '  "molecules": [',
            '    {',
            f'      "key": "{name}",',
            '      "role": "molecule",',
            f'      "smiles": "{rxn_smiles}",',
            f'      "xyz": "{name}.xyz"',
            '    }',
            '  ]',
            '}',
            'EOF',
            '',
            '# Export data-only archive + clean intermediates',
            f'$PYTHON $TOOLS/rxn_export.py $OUT --clean',
            '',
            'echo "DONE"'
        ]

        with open(script_path, 'w') as f:
            f.write('\n'.join(lines) + '\n')
        os.chmod(script_path, 0o755)

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
            self.end_headers()
            self.wfile.write(data)
        except KeyError:
            self.send_error(404, f'Not in zip: {inner_path}')
        except Exception as e:
            self.send_error(500, str(e))


    def _api_bundle_package(self, body):
        import zipfile
        name = body.get('name')
        keep_files = set(body.get('keep_files', []))
        
        zip_path = os.path.join(OUTPUT_DIR, name, f'{name}.rxnbundle.zip')
        if not os.path.exists(zip_path):
            self._json_response({'error': 'Bundle not found'})
            return
            
        custom_zip_path = os.path.join(OUTPUT_DIR, name, f'{name}_custom.rxnbundle.zip')
        with zipfile.ZipFile(zip_path, 'r') as zin, zipfile.ZipFile(custom_zip_path, 'w', compression=zipfile.ZIP_DEFLATED) as zout:
            for item in zin.infolist():
                if item.filename == 'manifest.json' or item.filename.endswith('.xyz') or item.filename in keep_files or item.filename.endswith('.json'): 
                    # Actually if we keep ALL json it defeats the purpose for orbitals.
                    # Wait, manifest is .json. Orbitals are _pos.json, _neg.json.
                    pass
                if item.filename == 'manifest.json' or item.filename.endswith('.xyz') or item.filename in keep_files:
                    zout.writestr(item, zin.read(item.filename))
        
        self._json_response({'download_url': f'/output/{name}/{name}_custom.rxnbundle.zip'})

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

