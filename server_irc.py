import os
import http.server
import socketserver
import json
import urllib.parse
import shutil

import irc_plugins

def get_plugins():
    return irc_plugins.get_plugins()

PORT = 8082
ROOT_DIR = os.path.dirname(os.path.abspath(__file__))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
os.makedirs(OUT_DIR, exist_ok=True)

from irc_plugins import project_manager

pipeline_state = {}


def load_saved_state():
    global pipeline_state
    try:
        data = project_manager.load_active_state()
        pipeline_state.clear()
        pipeline_state.update(data)
        cur_id = project_manager.init_project().get("active_step_id", "step_1")
        print(f"Loaded persistent reaction state for active step [{cur_id}]: {list(pipeline_state.keys())}")
    except Exception as e:
        print(f"Failed to load state: {e}")


def save_current_state():
    try:
        project_manager.save_active_state(pipeline_state)
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

    def translate_path(self, path):
        clean_path = path.split('?', 1)[0].split('#', 1)[0]
        if clean_path.startswith('/ketcher/'):
            subpath = clean_path[len('/ketcher/'):]
            return os.path.join(ROOT_DIR, '.ketcher', *[p for p in subpath.split('/') if p])
        if clean_path == '/':
            return os.path.join(ROOT_DIR, 'templates_irc', 'index.html')
        return super().translate_path(path)

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

        if path == '/' or path.startswith(('/ketcher/', '/static_irc/', '/static_viewer/', '/out/', '/bundles/')) or path.endswith('.html'):
            return super().do_GET()

        if path == '/api/log':
            return self._handle_get_log()

        if path == '/api/plugins':
            return self._handle_get_plugins()

        if path == '/api/state':
            return self._send_json(pipeline_state)

        if path == '/api/project':
            return self._send_json(project_manager.get_project_summary())

        if path == '/api/project/archive/export':
            from irc_plugins import project_archive
            zip_bytes, filename = project_archive.export_project_archive()
            self.send_response(200)
            self.send_header('Content-Type', 'application/zip')
            self.send_header('Content-Disposition', f'attachment; filename="{filename}"')
            self.send_header('Content-Length', str(len(zip_bytes)))
            self.end_headers()
            self.wfile.write(zip_bytes)
            return

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

        saved_docked = pipeline_state.get("step_3_ts", {}).get("docked_reactants")
        if plugin_idx >= 0:
            for downstream in plugins[plugin_idx + 1:]:
                if downstream.step_id in pipeline_state:
                    del pipeline_state[downstream.step_id]

        result = plugin.run(pipeline_state, user_input)
        if result.success:
            pipeline_state[plugin.step_id] = result.data
            if saved_docked and "step_3_ts" not in pipeline_state:
                raw = pipeline_state.get("step_1_input", {}).get("raw_smiles", "")
                r_part = raw.split(">>")[0] if ">>" in raw else ""
                r_smiles = [s.strip() for s in r_part.split(".") if s.strip()]
                try:
                    from rdkit import Chem
                    mols = [Chem.MolFromMolBlock(b, removeHs=False) for b in saved_docked]
                    can_docked = [Chem.CanonSmiles(Chem.MolToSmiles(Chem.RemoveHs(m))) for m in mols if m]
                    can_reac = [Chem.CanonSmiles(s) for s in r_smiles]
                    if sorted(can_docked) == sorted(can_reac):
                        pipeline_state["step_3_ts"] = {"docked_reactants": saved_docked}
                except Exception:
                    pass
            save_current_state()

        self._send_json({
            "success": result.success,
            "data": result.data,
            "error": result.error_message,
            "warnings": result.warnings
        })

    def _handle_step_custom(self, plugin, user_input):
        result_data = plugin.custom_action(user_input, pipeline_state) if hasattr(plugin, 'custom_action') else {}
        if result_data.get("success") and "molblocks" in result_data:
            pipeline_state.setdefault(plugin.step_id, {})["docked_reactants"] = result_data["molblocks"]
            save_current_state()
        if result_data.get("success") and "trajectory" in result_data:
            pipeline_state.setdefault(plugin.step_id, {})["trajectory"] = result_data["trajectory"]
            save_current_state()
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

            if path == '/api/project/step/select':
                data = self._read_json_body()
                target_id = data.get("step_id", "step_1")
                new_state, summary = project_manager.select_step(target_id, pipeline_state)
                pipeline_state.clear()
                pipeline_state.update(new_state)
                return self._send_json({"success": True, "active_step_id": target_id, "project": summary, "state": pipeline_state})

            if path == '/api/project/step/add':
                data = self._read_json_body()
                seed = bool(data.get("seed", True))
                new_id, new_state, summary = project_manager.add_step(seed_from_previous=seed, current_state=pipeline_state)
                pipeline_state.clear()
                pipeline_state.update(new_state)
                return self._send_json({"success": True, "active_step_id": new_id, "project": summary, "state": pipeline_state})

            if path == '/api/project/step/delete':
                data = self._read_json_body()
                target_id = data.get("step_id")
                new_id, new_state, summary = project_manager.delete_step(target_id)
                pipeline_state.clear()
                pipeline_state.update(new_state)
                return self._send_json({"success": True, "active_step_id": new_id, "project": summary, "state": pipeline_state})

            if path == '/api/project/stitch':
                save_current_state()
                success, master, msg = project_manager.stitch_project_trajectories()
                if success:
                    pipeline_state.setdefault("step_5_cascade", {})["megajson"] = master
                    return self._send_json({"success": True, "megajson": master, "message": msg})
                else:
                    return self._send_json({"success": False, "error": msg}, 400)

            if path == '/api/project/archive/import':
                from irc_plugins import project_archive
                content_len = int(self.headers.get('Content-Length', 0))
                raw_body = self.rfile.read(content_len) if content_len > 0 else b""
                if b'PK\x03\x04' in raw_body:
                    pk_start = raw_body.find(b'PK\x03\x04')
                    eocd = raw_body.rfind(b'PK\x05\x06')
                    if eocd != -1 and eocd + 22 <= len(raw_body):
                        comment_len = int.from_bytes(raw_body[eocd + 20 : eocd + 22], "little")
                        raw_body = raw_body[pk_start : eocd + 22 + comment_len]
                    else:
                        raw_body = raw_body[pk_start:]

                success, summary, new_state, msg = project_archive.import_project_archive(raw_body)
                if success:
                    pipeline_state.clear()
                    pipeline_state.update(new_state)
                    return self._send_json({"success": True, "project": summary, "state": pipeline_state, "message": msg})
                else:
                    return self._send_json({"success": False, "error": msg}, 400)

            if path == '/api/reset':
                data = self._read_json_body()
                if data.get("all_steps"):
                    shutil.rmtree(project_manager.STEPS_DIR, ignore_errors=True)
                    if os.path.exists(project_manager.PROJECT_FILE):
                        try:
                            os.remove(project_manager.PROJECT_FILE)
                        except Exception:
                            pass
                    pipeline_state.clear()
                    project_manager.init_project()
                else:
                    pipeline_state.clear()
                    save_current_state()
                return self._send_json({"success": True, "project": project_manager.get_project_summary()})

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
    import sys
    if len(sys.argv) > 1:
        try:
            PORT = int(sys.argv[1])
        except ValueError:
            pass
    load_saved_state()
    socketserver.ThreadingTCPServer.allow_reuse_address = True
    socketserver.ThreadingTCPServer.daemon_threads = True
    with socketserver.ThreadingTCPServer(("", PORT), IRCHandler) as httpd:
        print(f"IRC Master Tool server running at http://localhost:{PORT}")
        httpd.serve_forever()
