import http.server
import json
import os
import socketserver
import urllib.parse
from multistep_plugins import get_all_steps, get_step_by_id

PORT = 8083
STATE_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'multistep_state.json')
pipeline_state = {}


def load_saved_state():
    global pipeline_state
    if os.path.exists(STATE_FILE):
        try:
            with open(STATE_FILE, 'r') as f:
                data = json.load(f)
                pipeline_state.clear()
                pipeline_state.update(data)
                print(f"[Multistep] Loaded persistent state: {list(pipeline_state.keys())}")
        except Exception as e:
            print(f"[Multistep] Failed to load state: {e}")


def save_current_state():
    try:
        with open(STATE_FILE, 'w') as f:
            json.dump(pipeline_state, f, indent=2)
    except Exception as e:
        print(f"[Multistep] Failed to save state: {e}")


class MultistepHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Prevent caching for live development
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()

    def _send_json(self, data, status=200):
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.end_headers()
        self.wfile.write(json.dumps(data).encode())

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
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path == '/':
            self.path = '/templates_multistep/index.html'
            return super().do_GET()

        if path.startswith('/static_multistep/') or path.startswith('/trajectories/'):
            return super().do_GET()

        if path == '/api/plugins':
            return self._handle_get_plugins()

        if path == '/api/state':
            return self._send_json(pipeline_state)

        self.send_response(404)
        self.end_headers()

    def _handle_step_run(self, step, user_input):
        global pipeline_state
        steps = get_all_steps()
        step_idx = steps.index(step)

        # Invalidate downstream steps
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


if __name__ == '__main__':
    load_saved_state()
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    socketserver.ThreadingTCPServer.daemon_threads = True
    with socketserver.ThreadingTCPServer(("", PORT), MultistepHandler) as httpd:
        print(f"Multistep Reaction Pathway Server running at http://localhost:{PORT}")
        httpd.serve_forever()
