import http.server
import socketserver
import json
import os
import urllib.parse
from irc_plugins import get_plugins

PORT = 8082

# Global state for the single-user local app
# Form: { step_id: { data }, ... }
pipeline_state = {}

class IRCHandler(http.server.SimpleHTTPRequestHandler):
    def end_headers(self):
        # Prevent caching for development
        self.send_header('Cache-Control', 'no-cache, no-store, must-revalidate')
        self.send_header('Pragma', 'no-cache')
        self.send_header('Expires', '0')
        super().end_headers()

    def do_GET(self):
        parsed = urllib.parse.urlparse(self.path)
        
        if parsed.path == '/':
            self.path = '/templates_irc/index.html'
            return super().do_GET()
        elif parsed.path.startswith('/static_irc/'):
            return super().do_GET()
            
        elif parsed.path == '/api/plugins':
            plugins = get_plugins()
            resp = [
                {
                    "step_id": p.step_id, 
                    "step_name": p.step_name,
                    "ui_schema": p.get_ui_schema(pipeline_state)
                } 
                for p in plugins
            ]
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(resp).encode())
            return
            
        elif parsed.path == '/api/state':
            self.send_response(200)
            self.send_header('Content-Type', 'application/json')
            self.end_headers()
            self.wfile.write(json.dumps(pipeline_state).encode())
            return

        self.send_response(404)
        self.end_headers()

    def do_POST(self):
        global pipeline_state
        parsed = urllib.parse.urlparse(self.path)
        
        if parsed.path.startswith('/api/step/'):
            parts = parsed.path.split('/')
            if len(parts) >= 5 and parts[4] == 'run':
                step_id = parts[3]
                
                content_length = int(self.headers['Content-Length'])
                post_data = self.rfile.read(content_length)
                user_input = json.loads(post_data.decode()) if post_data else {}
                
                plugins = get_plugins()
                plugin = next((p for p in plugins if p.step_id == step_id), None)
                
                if not plugin:
                    self.send_response(404)
                    self.end_headers()
                    self.wfile.write(json.dumps({"error": "Plugin not found"}).encode())
                    return
                    
                # Invalidate downstream steps before running this one
                plugin_idx = plugins.index(plugin)
                for downstream in plugins[plugin_idx+1:]:
                    if downstream.step_id in pipeline_state:
                        del pipeline_state[downstream.step_id]

                result = plugin.run(pipeline_state, user_input)
                
                if result.success:
                    pipeline_state[step_id] = result.data
                    
                self.send_response(200)
                self.send_header('Content-Type', 'application/json')
                self.end_headers()
                self.wfile.write(json.dumps({
                    "success": result.success,
                    "data": result.data,
                    "error": result.error_message,
                    "warnings": result.warnings
                }).encode())
                return
                
        elif parsed.path == '/api/reset':
            pipeline_state = {}
            self.send_response(200)
            self.end_headers()
            self.wfile.write(b'{}')
            return

        self.send_response(404)
        self.end_headers()

if __name__ == '__main__':
    with socketserver.TCPServer(("", PORT), IRCHandler) as httpd:
        print(f"IRC Master Tool server running at http://localhost:{PORT}")
        httpd.serve_forever()
