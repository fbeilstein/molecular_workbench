import json
from typing import Any, Dict
from .standalone_css import STANDALONE_CSS
from .standalone_js import STANDALONE_JS


def generate_standalone_html(megajson: Dict[str, Any], manifest: Dict[str, Any]) -> str:
    """Create a complete, self-contained standalone HTML viewer running out of the box."""
    # Sanitize JSON payloads for script embedding
    json_data = json.dumps(megajson).replace("</script>", "<\\/script>")
    svg_content = megajson.get("mechanism_svg", "")
    wiring_data = json.dumps(megajson.get("svg_wiring", {})).replace("</script>", "<\\/script>")

    script_body = STANDALONE_JS.replace("__JSON_DATA__", json_data).replace("__WIRING_DATA__", wiring_data)

    total_frames = manifest.get("total_frames", 1)
    max_frame = max(0, total_frames - 1)
    title = manifest.get("title", "Multistep Reaction Pathway Studio")
    n_steps = manifest.get("num_steps", 1)
    res = manifest.get("orbital_resolution", 32)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{title}</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
{STANDALONE_CSS}
    </style>
</head>
<body>
    <div class="header">
        <div>
            <h2>{title}</h2>
            <div class="chips">
                <span class="chip blue">{n_steps} Reaction Steps</span>
                <span class="chip amber">{total_frames} Frames</span>
                <span class="chip emerald">Intrinsic Reaction Coordinates</span>
                <span class="chip">Grid: {res}³</span>
            </div>
        </div>
        <div style="font-size: 0.75rem; color: #64748b; font-family: monospace;">Standalone Bundle Viewer</div>
    </div>

    <!-- Dual Directional Playback Transport Bar -->
    <div class="transport-bar">
        <button id="btn-play-fwd" class="btn btn-primary">▶ Play Fwd</button>
        <button id="btn-play-rev" class="btn btn-amber">◀ Play Rev</button>
        <button id="btn-pause" class="btn">⏸ Pause</button>
        <span id="dir-badge" class="dir-badge dir-fwd">Forward ⏩</span>
        <input type="range" id="scrubber" min="0" max="{max_frame}" value="0" style="flex:1; accent-color: #38bdf8; cursor: pointer;">
        <span id="readout" style="font-size: 0.8rem; font-family: monospace; min-width: 105px; text-align: right;">Frame 0 / {max_frame}</span>
    </div>

    <!-- Main Workspace: 3D Pathway & Synchronized Mechanism SVG -->
    <div class="workspace">
        <div class="viewer-box">
            <div id="viewer3d" style="width: 100%; height: 100%;"></div>
            <div id="v-status" style="position: absolute; bottom: 8px; left: 8px; font-size: 0.72rem; color: #94a3b8; background: rgba(0,0,0,0.6); padding: 2px 6px; border-radius: 4px;">Ready</div>
        </div>
        <div class="mechanism-box">
            <div class="mechanism-header">
                <span>Reaction Mechanism Diagram</span>
                <span style="font-size: 0.72rem; color: #38bdf8;">Click arrows to play segment • Right-click molecules for orbitals</span>
            </div>
            <div id="svg-view" class="svg-container">
                {svg_content}
            </div>
        </div>
    </div>

    <!-- PES Free Energy Chart -->
    <div class="chart-box">
        <canvas id="pes-canvas"></canvas>
    </div>

    <script>
{script_body}
    </script>
</body>
</html>"""
