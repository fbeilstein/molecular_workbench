window.WB = window.WB || {};

window.WB.shiftColor = function(hex) {
    if (!hex || hex.length < 7) return '#cccccc';
    const r = parseInt(hex.slice(1, 3), 16);
    const g = parseInt(hex.slice(3, 5), 16);
    const b = parseInt(hex.slice(5, 7), 16);
    const nr = Math.min(255, 255 - r + 40);
    const ng = Math.min(255, 255 - g + 40);
    const nb = Math.min(255, 255 - b + 40);
    return `#${nr.toString(16).padStart(2, '0')}${ng.toString(16).padStart(2, '0')}${nb.toString(16).padStart(2, '0')}`;
};

window.WB.lighten = function(hex) {
    // Basic lighten for missing function
    return window.WB.shiftColor(hex);
};

window.WB.remeshCurrentOrbital = async function(isovalue) {
    window._customIsovalue = parseFloat(isovalue);
    window._currentOrbitalCache = {}; // clear cache so next redraw refetches/remeshes
    if (window.WB.onChangeHandler) {
        window.WB.onChangeHandler();
    }
};

window.WB.fetchOrbitalMesh = async function(file, sign, base) {
    try {
        const actualFile = file.replace('.json', `_${sign}.json`);
        let mesh = window._currentOrbitalCache[actualFile];
        if (!mesh) {
            if (window._customIsovalue) {
                // Remesh dynamically
                let name = document.getElementById('job-name').value || document.getElementById('run-chooser').value;
                const chooser = document.getElementById('chemical-chooser');
                if (chooser && chooser.value && chooser.value.startsWith('bundle:')) {
                    name = chooser.value.split(':')[1];
                }
                
                const res = await fetch('/api/orbital/remesh', {
                    method: 'POST',
                    headers: {'Content-Type': 'application/json'},
                    body: JSON.stringify({ name: name, cube_file: file, isovalue: window._customIsovalue })
                });
                if (!res.ok) return null;
                const data = await res.json();
                if (!data.error) {
                    window._currentOrbitalCache[file.replace('.json', '_pos.json')] = data.pos;
                    window._currentOrbitalCache[file.replace('.json', '_neg.json')] = data.neg;
                    mesh = window._currentOrbitalCache[actualFile];
                }
            } else {
                // Fetch from bundle
                const res = await fetch(`${base}/${actualFile}`);
                if (!res.ok) return null;
                mesh = await res.json();
                window._currentOrbitalCache[actualFile] = mesh;
            }
        }
        return mesh;
    } catch (e) {
        return null;
    }
};

let currentRenderId = 0;

window.WB.onChangeHandler = async function() {
    const renderId = ++currentRenderId;
    
    const sidebar = document.getElementById('orbital-toggles');
    const shapesToAdd = [];
    const checkedBoxes = sidebar.querySelectorAll('input[data-orb-file]:checked');
    
    // Get the base path for fetching meshes
    const chooser = document.getElementById('chemical-chooser');
    if (!chooser || !chooser.value) return;
    const parts = chooser.value.split(':');
    const bundleName = parts[1];
    const base = `/api/bundle/${bundleName}`;
    
    for (const chk of checkedBoxes) {
        const file = chk.getAttribute('data-orb-file');
        const color = chk.getAttribute('data-orb-color');
        if (!file) continue;
        
        if (file.endsWith('.json')) {
            for (let sign of ['pos', 'neg']) {
                const mesh = await window.WB.fetchOrbitalMesh(file, sign, base);
                if (mesh && mesh.vertices) {
                    const finalColor = (sign === 'neg') ? window.WB.shiftColor(color) : color;
                    const alpha = (sign === 'neg') ? 0.5 : 0.7;
                    shapesToAdd.push({vertexArr: mesh.vertices, faceArr: mesh.faces, color: $3Dmol.CC.color(finalColor), opacity: alpha});
                }
            }
        } else if (file.endsWith('.cube')) {
            let txt = window._currentOrbitalCache[file];
            if (!txt) {
                try {
                    const r = await fetch(`${base}/${file}`);
                    const buf = await r.arrayBuffer();
                    txt = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'))).text();
                    window._currentOrbitalCache[file] = txt;
                } catch(e) { continue; }
            }
            const vol = new $3Dmol.VolumeData(txt, 'cube');
            shapesToAdd.push({vol: vol, isoval: 0.03, color: color, opacity: 0.55, smoothness: 1});
            shapesToAdd.push({vol: vol, isoval: -0.03, color: window.WB.lighten(color), opacity: 0.35, smoothness: 1});
        }
    }

    if (renderId !== currentRenderId) return; // Abort if a newer render started
    
    const viewer = window.WB.viewer;
    if (!viewer) return;
    
    viewer.removeAllModels(); 
    viewer.removeAllShapes();
    viewer.removeAllSurfaces();
    
    if (window._currentMolXyz) {
        viewer.addModel(window._currentMolXyz, 'xyz');
        viewer.setStyle({}, {stick:{radius:0.12}, sphere:{scale:0.25}});
    }
    
    // Check ESP surface
    const espCb = sidebar.querySelector('input[data-esp]:checked');
    if (espCb && window._currentEspCache && window._currentEspCache.vertices) {
        const espCache = window._currentEspCache;
        const espMin = espCache.esp_min;
        const espMax = espCache.esp_max;
        const absMax = Math.max(Math.abs(espMin), Math.abs(espMax), 20);
        const colors = espCache.esp_values.map(v => {
            const t = Math.max(-1, Math.min(1, v / absMax));
            let r, g, b;
            if (t < 0) { r = 1.0; g = 1.0 + t; b = 1.0 + t; } 
            else { r = 1.0 - t; g = 1.0 - t; b = 1.0; }
            return {r: r, g: g, b: b};
        });
        viewer.addCustom({
            vertexArr: espCache.vertices,
            faceArr: espCache.faces,
            normalArr: espCache.normals || undefined,
            color: colors.length ? colors : $3Dmol.CC.color('#aaaaaa'),
            opacity: 0.85,
        });
    }
    
    for (const shape of shapesToAdd) {
        if (shape.vol) {
            viewer.addIsosurface(shape.vol, {isoval: shape.isoval, color: shape.color, opacity: shape.opacity, smoothness: shape.smoothness});
        } else {
            viewer.addCustom(shape);
        }
    }
    viewer.render();
};
