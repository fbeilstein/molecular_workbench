/**
 * Popout window renderer for molecule orbitals & ESP.
 */

export function openMoleculePopout(moleculeName, bundleData) {
    const title = `Orbitals - ${moleculeName}`;
    const win = window.open('', '_blank', 'width=950,height=620');
    if (!win) {
        alert('Popup was blocked by your browser. Please allow popups for this site.');
        return;
    }

    const bDataStr = JSON.stringify(bundleData || {});
    win.document.write(`<!DOCTYPE html>
<html>
<head>
    <title>${title}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body { margin: 0; background: #0b0f17; color: #fff; font-family: sans-serif; display: flex; height: 100vh; }
        #sidebar { width: 280px; background: #111827; border-right: 1px solid rgba(255,255,255,0.1); padding: 12px; display: flex; flex-direction: column; gap: 8px; }
        #viewer { flex: 1; position: relative; }
        .item { padding: 6px 10px; background: rgba(255,255,255,0.05); border-radius: 4px; cursor: pointer; font-size: 13px; display: flex; justify-content: space-between; align-items: center; }
        .item.active { background: #0284c7; font-weight: bold; }
        .cat-title { font-size: 11px; text-transform: uppercase; color: #94a3b8; font-weight: bold; margin-top: 6px; }
    </style>
</head>
<body>
    <div id="sidebar">
        <h3 style="margin: 0 0 6px 0; font-size: 14px;">${moleculeName} Orbitals</h3>
        <div id="list" style="flex:1; overflow-y:auto; display:flex; flex-direction:column; gap:4px;"></div>
    </div>
    <div id="viewer"></div>
    <script>
        const bundle = ${bDataStr};
        const viewer = $3Dmol.createViewer(document.getElementById('viewer'), { backgroundColor: '#090d16' });
        const listEl = document.getElementById('list');
        const files = bundle.files || {};

        function renderMol() {
            if (bundle.xyz) {
                viewer.addModel(bundle.xyz, 'xyz');
                viewer.setStyle({}, { stick: { radius: 0.12, colorscheme: 'Jmol' }, sphere: { scale: 0.25, colorscheme: 'Jmol' } });
            }
        }

        function renderOrb(name) {
            viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllSurfaces();
            renderMol();
            const posKey = Object.keys(files).find(k => k === name + '_pos.json' || k.endsWith('/' + name + '_pos.json'));
            const negKey = Object.keys(files).find(k => k === name + '_neg.json' || k.endsWith('/' + name + '_neg.json'));
            if (posKey && negKey) {
                const pos = typeof files[posKey] === 'string' ? JSON.parse(files[posKey]) : files[posKey];
                const neg = typeof files[negKey] === 'string' ? JSON.parse(files[negKey]) : files[negKey];
                viewer.addCustom({ vertexArr: pos.vertices, faceArr: pos.faces, color: $3Dmol.CC.color('#38bdf8'), opacity: 0.80 });
                viewer.addCustom({ vertexArr: neg.vertices, faceArr: neg.faces, color: $3Dmol.CC.color('#f59e0b'), opacity: 0.80 });
            } else if (files[name]) {
                const vol = new $3Dmol.VolumeData(files[name], 'cube');
                viewer.addIsosurface(vol, { isoval: 0.03, color: '#38bdf8', opacity: 0.75 });
                viewer.addIsosurface(vol, { isoval: -0.03, color: '#f59e0b', opacity: 0.75 });
            }
            viewer.zoomTo(); viewer.render();
        }

        function renderEsp() {
            viewer.removeAllModels(); viewer.removeAllShapes(); viewer.removeAllSurfaces();
            renderMol();
            const espKey = Object.keys(files).find(k => k.endsWith('_esp.json'));
            if (espKey) {
                const esp = typeof files[espKey] === 'string' ? JSON.parse(files[espKey]) : files[espKey];
                const absMax = Math.max(Math.abs(esp.esp_min || -30), Math.abs(esp.esp_max || 30), 20);
                const colors = esp.esp_values.map(v => {
                    const t = Math.max(-1, Math.min(1, v / absMax));
                    return t < 0 ? { r: 1.0, g: 1.0 + t, b: 1.0 + t } : { r: 1.0 - t, g: 1.0 - t, b: 1.0 };
                });
                viewer.addCustom({ vertexArr: esp.vertices, faceArr: esp.faces, color: colors, opacity: 0.85 });
            }
            viewer.zoomTo(); viewer.render();
        }

        if (bundle.esp_surface) {
            const espDiv = document.createElement('div');
            espDiv.className = 'item';
            espDiv.innerHTML = '<span>⚡ ESP Surface</span>';
            espDiv.onclick = () => {
                document.querySelectorAll('.item').forEach(d => d.classList.remove('active'));
                espDiv.classList.add('active');
                renderEsp();
            };
            listEl.appendChild(espDiv);
        }

        const tree = bundle.orbitals_tree || [];
        if (tree.length > 0) {
            tree.forEach(cat => {
                const header = document.createElement('div');
                header.className = 'cat-title';
                header.textContent = cat.name;
                listEl.appendChild(header);
                (cat.items || []).forEach(item => {
                    const base = item.file ? item.file.replace(/^molecules\\//, '').replace(/\\.json$/, '') : item.label;
                    const div = document.createElement('div');
                    div.className = 'item';
                    div.innerHTML = '<span>' + item.label + '</span><span style="font-size:11px;opacity:0.7">' + (item.energy_ev ? item.energy_ev + ' eV' : '') + '</span>';
                    div.onclick = () => {
                        document.querySelectorAll('.item').forEach(d => d.classList.remove('active'));
                        div.classList.add('active');
                        renderOrb(base);
                    };
                    listEl.appendChild(div);
                });
            });
        }
        const firstItem = listEl.querySelector('.item');
        if (firstItem) firstItem.click();
        else renderMol();
    </script>
</body>
</html>`);
}
