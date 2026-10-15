import json
from typing import Any, Dict


def generate_standalone_html(megajson: Dict[str, Any], manifest: Dict[str, Any]) -> str:
    """Create a complete, self-contained standalone HTML viewer running out of the box."""
    # Sanitize JSON payloads for script embedding
    json_data = json.dumps(megajson).replace("</script>", "<\\/script>")
    svg_content = megajson.get("mechanism_svg", "")
    wiring_data = json.dumps(megajson.get("svg_wiring", {})).replace("</script>", "<\\/script>")

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="utf-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{manifest.get("title", "Multistep Reaction Pathway Studio")}</title>
    <link rel="preconnect" href="https://fonts.googleapis.com">
    <link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
    <link href="https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&family=JetBrains+Mono:wght@400;500;600&display=swap" rel="stylesheet">
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        :root {{
            --bg-body: #0b0f17;
            --bg-surface: #111827;
            --border-color: rgba(255, 255, 255, 0.08);
            --primary: #38bdf8;
            --text-main: #f8fafc;
            --text-muted: #94a3b8;
        }}
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: 'Inter', -apple-system, sans-serif;
            background: var(--bg-body);
            color: var(--text-main);
            min-height: 100vh;
            display: flex;
            flex-direction: column;
            padding: 1rem 1.5rem;
            gap: 1rem;
        }}
        .header {{
            display: flex;
            justify-content: space-between;
            align-items: center;
            background: rgba(17, 24, 39, 0.85);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.75rem 1.25rem;
        }}
        .chips {{ display: flex; gap: 0.5rem; margin-top: 0.25rem; }}
        .chip {{
            background: rgba(255,255,255,0.06);
            border: 1px solid var(--border-color);
            padding: 2px 8px;
            border-radius: 9999px;
            font-size: 0.75rem;
        }}
        .chip.blue {{ color: #38bdf8; }}
        .chip.amber {{ color: #f59e0b; }}
        .chip.emerald {{ color: #10b981; }}
        .transport-bar {{
            display: flex;
            align-items: center;
            gap: 0.6rem;
            background: rgba(0,0,0,0.3);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.5rem 0.8rem;
        }}
        .btn {{
            background: rgba(255,255,255,0.08);
            border: 1px solid var(--border-color);
            color: #fff;
            padding: 6px 12px;
            border-radius: 6px;
            cursor: pointer;
            font-size: 0.82rem;
            font-weight: 600;
            transition: all 0.15s ease;
        }}
        .btn:hover {{ background: rgba(255,255,255,0.16); }}
        .btn-primary {{ background: #2563eb; border-color: #3b82f6; }}
        .btn-primary:hover {{ background: #1d4ed8; }}
        .btn-amber {{ background: #d97706; border-color: #f59e0b; }}
        .dir-badge {{
            font-size: 0.75rem;
            padding: 3px 8px;
            border-radius: 4px;
            font-weight: 700;
        }}
        .dir-fwd {{ background: rgba(56,189,248,0.15); color: #38bdf8; border: 1px solid rgba(56,189,248,0.3); }}
        .dir-rev {{ background: rgba(245,158,11,0.15); color: #f59e0b; border: 1px solid rgba(245,158,11,0.3); }}
        .workspace {{
            display: flex;
            gap: 1rem;
            height: 480px;
        }}
        .viewer-box {{
            flex: 1;
            position: relative;
            background: #090d16;
            border: 1px solid var(--border-color);
            border-radius: 8px;
            overflow: hidden;
        }}
        .mechanism-box {{
            flex: 1;
            position: relative;
            background: rgba(17, 24, 39, 0.75);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            display: flex;
            flex-direction: column;
            overflow: hidden;
        }}
        .mechanism-header {{
            padding: 0.5rem 0.75rem;
            background: rgba(0,0,0,0.25);
            border-bottom: 1px solid var(--border-color);
            font-size: 0.82rem;
            font-weight: 600;
            display: flex;
            justify-content: space-between;
        }}
        .svg-container {{
            flex: 1;
            display: flex;
            align-items: center;
            justify-content: center;
            padding: 0.5rem;
            overflow: auto;
        }}
        .svg-container svg {{
            max-width: 100%;
            max-height: 100%;
            display: block;
        }}
        /* SVG Interactive Classes */
        .svg-clickable-arrow {{ cursor: pointer; transition: all 0.2s ease; }}
        .svg-clickable-arrow:hover {{ filter: drop-shadow(0 0 8px #38bdf8); opacity: 0.9; }}
        .svg-clickable-state {{ cursor: pointer; transition: all 0.2s ease; }}
        .svg-clickable-state:hover {{ filter: drop-shadow(0 0 8px #10b981); }}
        .svg-clickable-molecule {{ cursor: pointer; }}
        .svg-clickable-molecule:hover {{ filter: drop-shadow(0 0 8px #a855f7); }}
        .svg-push-active {{
            opacity: 1.0 !important;
            filter: drop-shadow(0 0 6px #38bdf8) !important;
            transition: all 0.25s ease-in-out;
            animation: pulse-push 1.5s infinite;
        }}
        .svg-push-inactive {{
            opacity: 0.15 !important;
            pointer-events: none;
            transition: opacity 0.25s ease-in-out;
        }}
        @keyframes pulse-push {{
            0%, 100% {{ filter: drop-shadow(0 0 4px #38bdf8); }}
            50% {{ filter: drop-shadow(0 0 10px #0284c7); }}
        }}
        .chart-box {{
            background: rgba(0,0,0,0.3);
            border: 1px solid var(--border-color);
            border-radius: 8px;
            padding: 0.6rem;
            height: 140px;
        }}
        #pes-canvas {{ width: 100%; height: 100%; display: block; }}
    </style>
</head>
<body>
    <div class="header">
        <div>
            <h1 style="font-size: 1.15rem; font-weight: 700;">{manifest.get("title", "Multistep Reaction Pathway")}</h1>
            <div class="chips">
                <span class="chip blue">{len(manifest.get("symbols", []))} atoms</span>
                <span class="chip">{manifest.get("total_frames", 0)} frames</span>
                <span class="chip amber">Barrier: {(manifest.get("overall_barrier_kcal") or 0.0):.1f} kcal/mol</span>
                <span class="chip emerald">Net ΔG: {(manifest.get("overall_reaction_energy_kcal") or 0.0):.1f} kcal/mol</span>
                <span class="chip">Grid: {manifest.get("orbital_resolution", 32)}³</span>
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
        <input type="range" id="scrubber" min="0" max="{manifest.get("total_frames", 1) - 1}" value="0" style="flex:1; accent-color: #38bdf8; cursor: pointer;">
        <span id="readout" style="font-size: 0.8rem; font-family: monospace; min-width: 105px; text-align: right;">Frame 0 / {manifest.get("total_frames", 1) - 1}</span>
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
        const MEGDATA = {json_data};
        const WIRING = {wiring_data};
        const TOTAL_FRAMES = MEGDATA.frames ? MEGDATA.frames.length : 0;
        let curFrame = 0;
        let playDirection = 'forward';
        let isPlaying = false;
        let timer = null;

        // Init 3Dmol viewer
        const viewerEl = document.getElementById('viewer3d');
        const viewer = $3Dmol.createViewer(viewerEl, {{ backgroundColor: '#090d16' }});
        let hasCamInit = false;

        function renderFrame(f) {{
            curFrame = f;
            document.getElementById('scrubber').value = f;
            document.getElementById('readout').textContent = `Frame ${{f}} / ${{TOTAL_FRAMES - 1}}`;
            
            // Render 3D Molecule
            const frame = MEGDATA.frames[f];
            if (frame) {{
                viewer.clear();
                let n = MEGDATA.symbols.length;
                let xyz = n + '\\nFrame ' + f + '\\n';
                for (let i = 0; i < n; i++) {{
                    const s = MEGDATA.symbols[i];
                    const [x, y, z] = frame.coords[i];
                    xyz += s + ' ' + x.toFixed(4) + ' ' + y.toFixed(4) + ' ' + z.toFixed(4) + '\\n';
                }}
                viewer.addModel(xyz, 'xyz');
                viewer.setStyle({{}}, {{
                    stick: {{ radius: 0.14, colorscheme: 'Jmol' }},
                    sphere: {{ scale: 0.28, colorscheme: 'Jmol' }}
                }});
                if (!hasCamInit) {{
                    viewer.zoomTo();
                    hasCamInit = true;
                }}
                viewer.render();
            }}

            // Update Mechanism Push Arrows & Interactive States
            updateSvgPushArrows(f, playDirection);
            drawPes();
        }}

        function updateSvgPushArrows(f, dir) {{
            const svg = document.querySelector('#svg-view svg');
            if (!svg || !WIRING) return;

            Object.entries(WIRING).forEach(([grpId, cfg]) => {{
                const el = svg.getElementById(grpId) || svg.querySelector(`[inkscape\\\\:label="${{grpId}}"]`);
                if (!el) return;

                if (cfg.role === 'push_arrow') {{
                    const sf = cfg.startFrame !== undefined ? cfg.startFrame : 0;
                    const ef = cfg.endFrame !== undefined ? cfg.endFrame : TOTAL_FRAMES - 1;
                    const cDir = cfg.direction || 'both';
                    const dirMatch = (cDir === 'both') || (cDir === dir);
                    const frameMatch = (f >= sf && f <= ef);

                    if (dirMatch && frameMatch) {{
                        el.classList.add('svg-push-active');
                        el.classList.remove('svg-push-inactive');
                    }} else {{
                        el.classList.add('svg-push-inactive');
                        el.classList.remove('svg-push-active');
                    }}
                }}
            }});
        }}

        function wireInteractiveSvg() {{
            const svg = document.querySelector('#svg-view svg');
            if (!svg || !WIRING) return;

            Object.entries(WIRING).forEach(([grpId, cfg]) => {{
                const el = svg.getElementById(grpId) || svg.querySelector(`[inkscape\\\\:label="${{grpId}}"]`);
                if (!el) return;

                if (cfg.role === 'arrow') {{
                    el.classList.add('svg-clickable-arrow');
                    el.title = `Reaction Arrow: Play frames ${{cfg.startFrame || 0}} → ${{cfg.endFrame || TOTAL_FRAMES - 1}}`;
                    el.addEventListener('click', (e) => {{
                        e.stopPropagation();
                        playSegment(cfg.startFrame || 0, cfg.endFrame || TOTAL_FRAMES - 1, cfg.direction || 'forward');
                    }});
                }} else if (cfg.role === 'state') {{
                    el.classList.add('svg-clickable-state');
                    el.title = `State: Jump to frame ${{cfg.targetFrame || 0}}`;
                    el.addEventListener('click', (e) => {{
                        e.stopPropagation();
                        pause();
                        renderFrame(cfg.targetFrame || 0);
                    }});
                }} else if (cfg.role === 'molecule') {{
                    el.classList.add('svg-clickable-molecule');
                    el.title = `Molecule: Right click for orbitals`;
                }}
            }});
        }}

        function setDirection(dir) {{
            playDirection = dir;
            const badge = document.getElementById('dir-badge');
            if (badge) {{
                badge.className = `dir-badge ${{dir === 'forward' ? 'dir-fwd' : 'dir-rev'}}`;
                badge.textContent = dir === 'forward' ? 'Forward ⏩' : 'Reverse ⏪';
            }}
        }}

        function playForward() {{
            pause();
            setDirection('forward');
            if (curFrame >= TOTAL_FRAMES - 1) curFrame = 0;
            isPlaying = true;
            timer = setInterval(() => {{
                if (curFrame >= TOTAL_FRAMES - 1) {{ pause(); }}
                else {{ renderFrame(curFrame + 1); }}
            }}, 100);
        }}

        function playReverse() {{
            pause();
            setDirection('reverse');
            if (curFrame <= 0) curFrame = TOTAL_FRAMES - 1;
            isPlaying = true;
            timer = setInterval(() => {{
                if (curFrame <= 0) {{ pause(); }}
                else {{ renderFrame(curFrame - 1); }}
            }}, 100);
        }}

        function pause() {{
            isPlaying = false;
            if (timer) {{ clearInterval(timer); timer = null; }}
        }}

        function playSegment(sf, ef, dir) {{
            pause();
            setDirection(dir);
            curFrame = (dir === 'forward') ? sf : ef;
            renderFrame(curFrame);
            isPlaying = true;
            timer = setInterval(() => {{
                if (dir === 'forward') {{
                    if (curFrame >= ef) pause();
                    else renderFrame(curFrame + 1);
                }} else {{
                    if (curFrame <= sf) pause();
                    else renderFrame(curFrame - 1);
                }}
            }}, 90);
        }}

        // Draw PES Canvas
        const pesCanvas = document.getElementById('pes-canvas');
        const pesCtx = pesCanvas.getContext('2d');
        function drawPes() {{
            const w = pesCanvas.parentElement.clientWidth;
            const h = pesCanvas.parentElement.clientHeight;
            if (pesCanvas.width !== w) pesCanvas.width = w;
            if (pesCanvas.height !== h) pesCanvas.height = h;

            pesCtx.clearRect(0, 0, w, h);
            const ep = (MEGDATA.energy_profile && MEGDATA.energy_profile.frame_energies_kcal) || [];
            if (!ep.length) return;

            const minE = Math.min(...ep) - 1.0;
            const maxE = Math.max(...ep) + 2.0;
            const padX = 35, padY = 20;
            const pw = w - padX * 2, ph = h - padY * 2;

            // Draw line
            pesCtx.beginPath();
            pesCtx.strokeStyle = '#38bdf8';
            pesCtx.lineWidth = 2.5;
            ep.forEach((e, i) => {{
                const x = padX + (i / (ep.length - 1)) * pw;
                const y = h - padY - ((e - minE) / (maxE - minE)) * ph;
                if (i === 0) pesCtx.moveTo(x, y); else pesCtx.lineTo(x, y);
            }});
            pesCtx.stroke();

            // Cursor
            const cx = padX + (curFrame / (ep.length - 1)) * pw;
            const cy = h - padY - ((ep[curFrame] - minE) / (maxE - minE)) * ph;
            pesCtx.beginPath();
            pesCtx.arc(cx, cy, 6, 0, Math.PI * 2);
            pesCtx.fillStyle = '#ef4444';
            pesCtx.fill();
            pesCtx.strokeStyle = '#fff';
            pesCtx.lineWidth = 2;
            pesCtx.stroke();
        }}

        document.getElementById('btn-play-fwd').onclick = playForward;
        document.getElementById('btn-play-rev').onclick = playReverse;
        document.getElementById('btn-pause').onclick = pause;
        document.getElementById('scrubber').oninput = (e) => {{
            pause();
            renderFrame(+e.target.value);
        }};

        wireInteractiveSvg();
        renderFrame(0);
        window.addEventListener('resize', drawPes);
    </script>
</body>
</html>"""
