import http.server
import json
import os
import socketserver
import urllib.parse
from multistep_plugins import get_all_steps, get_step_by_id

PORT = 8083
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
os.makedirs(OUT_DIR, exist_ok=True)
STATE_FILE = os.path.join(OUT_DIR, 'multistep_state.json')
LEGACY_STATE_FILE = os.path.join(ROOT_DIR, 'multistep_state.json')
pipeline_state = {}


def load_saved_state():
    global pipeline_state
    sf = STATE_FILE if os.path.exists(STATE_FILE) else (LEGACY_STATE_FILE if os.path.exists(LEGACY_STATE_FILE) else None)
    if sf and os.path.exists(sf):
        try:
            with open(sf, 'r') as f:
                data = json.load(f)
                pipeline_state.clear()
                pipeline_state.update(data)
                print(f"[Multistep] Loaded persistent state from {sf}: {list(pipeline_state.keys())}")
        except Exception as e:
            print(f"[Multistep] Failed to load state: {e}")


def save_current_state():
    try:
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(STATE_FILE, 'w') as f:
            json.dump(pipeline_state, f, indent=2)
    except Exception as e:
        print(f"[Multistep] Failed to save state: {e}")


class MultistepHandler(http.server.SimpleHTTPRequestHandler):
    def log_message(self, format, *args):
        try:
            super().log_message(format, *args)
        except Exception:
            pass

    def end_headers(self):
        # Prevent caching for live development
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        self.send_header('Access-Control-Allow-Origin', '*')
        self.send_header('Access-Control-Allow-Methods', 'GET, POST, OPTIONS, HEAD')
        self.send_header('Access-Control-Allow-Headers', 'Content-Type')
        super().end_headers()

    def do_OPTIONS(self):
        self.send_response(204)
        self.end_headers()

    def do_HEAD(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path
        if path == '/':
            self.path = '/templates_multistep/index.html'
            return super().do_HEAD()
        if path.startswith(('/static_multistep/', '/static_viewer/', '/trajectories/', '/out/', '/bundles/')) or path.endswith('.html'):
            return super().do_HEAD()
        if path in ('/api/plugins', '/api/state'):
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.end_headers()
            return
        self.send_response(404)
        self.end_headers()

    def _send_json(self, data, status=200):
        try:
            body = json.dumps(data).encode('utf-8')
            self.send_response(status)
            self.send_header('Content-Type', 'application/json')
            self.send_header('Content-Length', str(len(body)))
            self.end_headers()
            self.wfile.write(body)
        except (BrokenPipeError, ConnectionResetError):
            pass

    def _read_json_body(self):
        content_length = int(self.headers.get('Content-Length', 0))
        if not content_length:
            return {}
        post_data = self.rfile.read(content_length)
        return json.loads(post_data.decode()) if post_data else {}

    def _handle_get_plugins(self):
        steps = get_all_steps()
        resp = [
            {
                "step_id": s.step_id,
                "step_name": s.step_name,
                "ui_schema": s.get_ui_schema(pipeline_state)
            }
            for s in steps
        ]
        self._send_json(resp)

    def do_GET(self):
        try:
            parsed = urllib.parse.urlparse(self.path)
            path = parsed.path

            if path == '/':
                self.path = '/templates_multistep/index.html'
                return super().do_GET()

            if path.startswith(('/static_multistep/', '/static_viewer/', '/trajectories/', '/out/', '/bundles/')) or path.endswith('.html'):
                return super().do_GET()

            if path == '/api/plugins':
                return self._handle_get_plugins()

            if path == '/api/state':
                clean_state = {}
                for k, v in pipeline_state.items():
                    if isinstance(v, dict) and "megajson" in v and isinstance(v["megajson"], dict):
                        v_copy = dict(v)
                        mj_copy = dict(v["megajson"])
                        if "highres_cubes" in mj_copy:
                            mj_copy["highres_cubes"] = {
                                ck: (cv if (isinstance(cv, str) and "\n" not in cv) else True)
                                for ck, cv in mj_copy["highres_cubes"].items()
                            }
                        if "molecule_bundles" in mj_copy:
                            mj_copy["molecule_bundles"] = {
                                m_id: {
                                    "filename": data.get("filename", "bundle.zip"),
                                    "orbitals": data.get("orbitals", []),
                                    "molecule_name": data.get("molecule_name", m_id),
                                    "molecule_key": data.get("molecule_key", m_id),
                                }
                                for m_id, data in mj_copy["molecule_bundles"].items()
                                if isinstance(data, dict)
                            }
                        v_copy["megajson"] = mj_copy
                        clean_state[k] = v_copy
                    else:
                        clean_state[k] = v
                return self._send_json(clean_state)

            self.send_response(404)
            self.end_headers()
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._send_json({"error": str(e)}, 500)

    def _handle_step_run(self, step, user_input):
        global pipeline_state
        steps = get_all_steps()
        step_idx = next((i for i, s in enumerate(steps) if s.step_id == step.step_id), -1)

        # Invalidate downstream steps
        if step_idx >= 0:
            for downstream in steps[step_idx + 1:]:
                if downstream.step_id in pipeline_state:
                    del pipeline_state[downstream.step_id]

        result = step.run(pipeline_state, user_input)
        if result.success:
            pipeline_state[step.step_id] = result.data
            save_current_state()

        self._send_json({
            "success": result.success,
            "data": result.data,
            "error": result.error_message,
            "warnings": result.warnings
        })

    def _handle_step_custom(self, step, user_input):
        print(f"CUSTOM_REQ: action={user_input.get('action_type')} f={user_input.get('frame_index')} orb={user_input.get('orbital')} chem={user_input.get('chem_key')} sel={len(user_input.get('selected_tracks', []))}")
        result_data = step.custom_action(user_input, pipeline_state) if hasattr(step, 'custom_action') else {}
        save_current_state()
        self._send_json(result_data)

    def _handle_post_step(self, path):
        parts = path.split('/')
        if len(parts) < 5:
            return self._send_json({"error": "Invalid step endpoint"}, 400)

        step_id, action = parts[3], parts[4]
        step = get_step_by_id(step_id)
        if not step:
            return self._send_json({"error": f"Step '{step_id}' not found"}, 404)

        user_input = self._read_json_body()
        if action == 'run':
            return self._handle_step_run(step, user_input)
        if action in ('custom', 'action'):
            return self._handle_step_custom(step, user_input)

        self._send_json({"error": f"Unknown action: {action}"}, 400)

    def do_POST(self):
        try:
            global pipeline_state
            parsed = urllib.parse.urlparse(self.path)
            path = parsed.path

            if path.startswith('/api/step/'):
                return self._handle_post_step(path)

            if path == '/api/reset':
                pipeline_state.clear()
                save_current_state()
                return self._send_json({"success": True})

            self.send_response(404)
            self.end_headers()
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._send_json({"error": str(e)}, 500)


if __name__ == '__main__':
    load_saved_state()
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    socketserver.ThreadingTCPServer.daemon_threads = True
    with socketserver.ThreadingTCPServer(("", PORT), MultistepHandler) as httpd:
        print(f"Multistep Reaction Pathway Server running at http://localhost:{PORT}")
        httpd.serve_forever()
