import os
import http.server
import socketserver
import json
import urllib.parse
def get_plugins():
    import sys
    for mod_name in list(sys.modules.keys()):
        if mod_name.startswith("irc_plugins"):
            del sys.modules[mod_name]
    import irc_plugins
    return irc_plugins.get_plugins()

PORT = 8082
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
os.makedirs(OUT_DIR, exist_ok=True)
STATE_FILE = os.path.join(OUT_DIR, 'reaction_state.json')
LEGACY_STATE_FILE = os.path.join(ROOT_DIR, 'reaction_state.json')
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
                print(f"Loaded persistent reaction state from {sf}: {list(pipeline_state.keys())}")
        except Exception as e:
            print(f"Failed to load state: {e}")


def save_current_state():
    try:
        os.makedirs(OUT_DIR, exist_ok=True)
        with open(STATE_FILE, 'w') as f:
            json.dump(pipeline_state, f, indent=2)
    except Exception as e:
        print(f"Failed to save state: {e}")


load_saved_state()


class IRCHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Prevent caching for development
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()

    def _send_json(self, data, status=200):
        self.send_response(status)
        self.send_header('Content-Type', 'application/json')
        self.end_headers()
        self.wfile.write(json.dumps(data).encode())

    def _send_text(self, text, status=200):
        self.send_response(status)
        self.send_header('Content-Type', 'text/plain')
        self.end_headers()
        self.wfile.write(text.encode())

    def _read_json_body(self):
        content_length = int(self.headers.get('Content-Length', 0))
        if not content_length:
            return {}
        post_data = self.rfile.read(content_length)
        return json.loads(post_data.decode()) if post_data else {}

    # -----------------------------------------------------------------------
    # GET Handlers
    # -----------------------------------------------------------------------
    def _handle_get_log(self):
        log_path = os.path.join(OUT_DIR, 'optimization.log')
        if not os.path.exists(log_path):
            log_path = os.path.join(ROOT_DIR, 'optimization.log')
        try:
            with open(log_path, 'r') as f:
                content = f.read()
            self._send_text(content)
        except FileNotFoundError:
            self._send_text("")

    def _handle_get_plugins(self):
        plugins = get_plugins()
        resp = [
            {
                "step_id": p.step_id,
                "step_name": p.step_name,
                "ui_schema": p.get_ui_schema(pipeline_state)
            }
            for p in plugins
        ]
        self._send_json(resp)

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        path = parsed.path

        if path == '/':
            self.path = '/templates_irc/index.html'
            return super().do_GET()

        if path.startswith('/static_irc/') or path.startswith('/out/'):
            return super().do_GET()

        if path.startswith('/ketcher/'):
            self.path = '/.ketcher/' + path[len('/ketcher/'):]
            return super().do_GET()

        if path == '/api/log':
            return self._handle_get_log()

        if path == '/api/plugins':
            return self._handle_get_plugins()

        if path == '/api/state':
            return self._send_json(pipeline_state)

        self.send_response(404)
        self.end_headers()

    # -----------------------------------------------------------------------
    # POST Handlers
    # -----------------------------------------------------------------------
    def _handle_step_run(self, plugin, user_input, plugins=None):
        global pipeline_state
        if plugins is None:
            plugins = get_plugins()
        plugin_idx = next((i for i, p in enumerate(plugins) if p.step_id == plugin.step_id), -1)

        # Invalidate downstream steps before running this one
        if plugin_idx >= 0:
            for downstream in plugins[plugin_idx + 1:]:
                if downstream.step_id in pipeline_state:
                    del pipeline_state[downstream.step_id]

        result = plugin.run(pipeline_state, user_input)
        if result.success:
            pipeline_state[plugin.step_id] = result.data
            save_current_state()

        self._send_json({
            "success": result.success,
            "data": result.data,
            "error": result.error_message,
            "warnings": result.warnings
        })

    def _handle_step_custom(self, plugin, user_input):
        result_data = plugin.custom_action(user_input, pipeline_state) if hasattr(plugin, 'custom_action') else {}
        self._send_json(result_data)

    def _handle_post_step(self, path):
        parts = path.split('/')
        if len(parts) < 5:
            return self._send_json({"error": "Invalid step endpoint"}, 400)

        step_id, action = parts[3], parts[4]
        plugins = get_plugins()
        plugin = next((p for p in plugins if p.step_id == step_id), None)
        if not plugin:
            return self._send_json({"error": f"Plugin '{step_id}' not found"}, 404)

        user_input = self._read_json_body()
        if action == 'run':
            return self._handle_step_run(plugin, user_input, plugins)
        if action == 'custom':
            return self._handle_step_custom(plugin, user_input)

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
                return self._send_json({})

            if path == '/api/environment':
                data = self._read_json_body()
                solvent = data.get("solvent", "water")
                pipeline_state["solvent"] = solvent
                save_current_state()
                return self._send_json({"success": True, "solvent": solvent})

            self.send_response(404)
            self.end_headers()
        except Exception as e:
            import traceback
            traceback.print_exc()
            self._send_json({"success": False, "error": str(e)}, 500)


if __name__ == '__main__':
    load_saved_state()
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    socketserver.ThreadingTCPServer.daemon_threads = True
    with socketserver.ThreadingTCPServer(("", PORT), IRCHandler) as httpd:
        print(f"IRC Master Tool server running at http://localhost:{PORT}")
        httpd.serve_forever()
