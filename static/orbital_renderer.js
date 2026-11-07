window.WB = window.WB || {};

// ── Unified Kekulé Bond Perception & SDF MolBlock Builder (shared with viewer.html & server_irc.py) ──

const WB_COV_RADII = {
    H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
    P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39
};
const WB_MAX_VALENCE = { H: 1, F: 1, Cl: 1, Br: 1, I: 1, O: 3, N: 4, C: 4, S: 6, P: 5 };

window.WB.parseXyz = function(xyzText) {
    if (!xyzText) return { symbols: [], coords: [] };
    const lines = xyzText.trim().split(/\r?\n/);
    if (lines.length < 3) return { symbols: [], coords: [] };
    const n = parseInt(lines[0].trim(), 10);
    const count = (!Number.isNaN(n) && n > 0) ? Math.min(n, lines.length - 2) : (lines.length - 2);
    const symbols = [];
    const coords = [];
    for (let i = 0; i < count; i++) {
        const parts = lines[i + 2].trim().split(/\s+/);
        if (parts.length >= 4) {
            symbols.push(parts[0]);
            coords.push([parseFloat(parts[1]), parseFloat(parts[2]), parseFloat(parts[3])]);
        }
    }
    return { symbols, coords };
};

window.WB.computeBonds = function(symbols, coords) {
    const bonds = [];
    const n = Math.min(symbols.length, coords.length);
    const valences = new Array(n).fill(0);
    const candidates = [];
    const adj = Array.from({ length: n }, () => []);
    const distMat = {};

    for (let i = 0; i < n; i++) {
        const s1 = symbols[i];
        const r1 = WB_COV_RADII[s1] || 0.75;
        for (let j = i + 1; j < n; j++) {
            const s2 = symbols[j];
            const r2 = WB_COV_RADII[s2] || 0.75;
            if (s1 === 'H' && s2 === 'H' && n > 2) continue;

            const dx = coords[i][0] - coords[j][0];
            const dy = coords[i][1] - coords[j][1];
            const dz = coords[i][2] - coords[j][2];
            const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
            if (d >= 0.4 && d <= 1.28 * (r1 + r2)) {
                adj[i].push(j);
                adj[j].push(i);
                distMat[`${i}_${j}`] = d;
                candidates.push({ i, j, s1, s2, d });
            }
        }
    }

    const getDist = (u, v) => distMat[u < v ? `${u}_${v}` : `${v}_${u}`] || 1.4;
    const ringDoubleBonds = new Set();
    const ringBondsAll = new Set();
    const cycles = [];
    const visited = new Set();
    const path = [];

    const findCycles = (v, start, depth) => {
        if (depth === 6) {
            if (adj[v].includes(start)) cycles.push([...path]);
            return;
        }
        for (const u of adj[v]) {
            if ((symbols[u] === 'C' || symbols[u] === 'N') && !visited.has(u)) {
                visited.add(u);
                path.push(u);
                findCycles(u, start, depth + 1);
                path.pop();
                visited.delete(u);
            }
        }
    };

    for (let i = 0; i < n; i++) {
        if (symbols[i] === 'C' || symbols[i] === 'N') {
            visited.add(i);
            path.push(i);
            findCycles(i, i, 1);
            path.pop();
            visited.delete(i);
        }
    }

    const seenRings = new Set();
    for (const ring of cycles) {
        const sorted = [...ring].sort((a, b) => a - b).join(',');
        if (seenRings.has(sorted)) continue;
        seenRings.add(sorted);

        let isConjugated = true;
        for (let k = 0; k < 6; k++) {
            if (getDist(ring[k], ring[(k + 1) % 6]) > 1.48) {
                isConjugated = false;
                break;
            }
        }
        if (!isConjugated) continue;

        const edges0 = [[ring[0], ring[1]], [ring[2], ring[3]], [ring[4], ring[5]]];
        const edges1 = [[ring[1], ring[2]], [ring[3], ring[4]], [ring[5], ring[0]]];
        const sum0 = edges0.reduce((acc, [u, v]) => acc + getDist(u, v), 0);
        const sum1 = edges1.reduce((acc, [u, v]) => acc + getDist(u, v), 0);
        const chosen = sum0 <= sum1 ? edges0 : edges1;

        for (const [u, v] of chosen) ringDoubleBonds.add(u < v ? `${u}_${v}` : `${v}_${u}`);
        for (let k = 0; k < 6; k++) {
            const u = ring[k], v = ring[(k + 1) % 6];
            ringBondsAll.add(u < v ? `${u}_${v}` : `${v}_${u}`);
        }
    }

    for (const c of candidates) {
        const pairKey = c.i < c.j ? `${c.i}_${c.j}` : `${c.j}_${c.i}`;
        const pair = [c.s1, c.s2].sort().join('-');
        let targetOrder = 1;

        if (ringDoubleBonds.has(pairKey)) targetOrder = 2;
        else if (!ringBondsAll.has(pairKey)) {
            if (pair === 'C-O') targetOrder = c.d < 1.16 ? 3 : (c.d < 1.265 ? 2 : 1);
            else if (pair === 'C-C') targetOrder = c.d < 1.24 ? 3 : (c.d < 1.39 ? 2 : 1);
            else if (pair === 'C-N') targetOrder = c.d < 1.18 ? 3 : (c.d < 1.32 ? 2 : 1);
            else if (pair === 'N-N') targetOrder = c.d < 1.14 ? 3 : (c.d < 1.28 ? 2 : 1);
            else if (pair === 'N-O' || pair === 'O-O') targetOrder = c.d < 1.25 ? 2 : 1;
        }

        bonds.push({ a1: c.i, a2: c.j, order: 1, target: targetOrder, s1: c.s1, s2: c.s2, d: c.d, pairKey });
        valences[c.i]++;
        valences[c.j]++;
    }

    bonds.sort((a, b) => {
        const prioA = a.target === 3 ? 0 : (ringDoubleBonds.has(a.pairKey) ? 1 : (a.target === 2 ? 2 : 3));
        const prioB = b.target === 3 ? 0 : (ringDoubleBonds.has(b.pairKey) ? 1 : (b.target === 2 ? 2 : 3));
        return prioA !== prioB ? prioA - prioB : a.d - b.d;
    });

    for (const b of bonds) {
        if (b.target > 1) {
            const v1Cap = WB_MAX_VALENCE[b.s1] || 4;
            const v2Cap = WB_MAX_VALENCE[b.s2] || 4;
            const maxBoost = Math.min(v1Cap - valences[b.a1], v2Cap - valences[b.a2]);
            const boost = Math.max(0, Math.min(b.target - 1, maxBoost));
            if (boost > 0) {
                b.order += boost;
                valences[b.a1] += boost;
                valences[b.a2] += boost;
            }
        }
    }

    return bonds;
};

window.WB.buildMolBlock = function(symbols, coords) {
    const natoms = Math.min(symbols.length, coords.length);
    const bonds = window.WB.computeBonds(symbols, coords);
    let sdf = `\n  Molecular Workbench\n\n${String(natoms).padStart(3, ' ')}${String(bonds.length).padStart(3, ' ')}  0  0  0  0  0  0  0  0999 V2000\n`;
    for (let i = 0; i < natoms; i++) {
        const x = Number(coords[i][0]).toFixed(4).padStart(10, ' ');
        const y = Number(coords[i][1]).toFixed(4).padStart(10, ' ');
        const z = Number(coords[i][2]).toFixed(4).padStart(10, ' ');
        sdf += `${x}${y}${z} ${symbols[i].padEnd(3, ' ')} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
    }
    for (let i = 0; i < bonds.length; i++) {
        sdf += `${String(bonds[i].a1 + 1).padStart(3, ' ')}${String(bonds[i].a2 + 1).padStart(3, ' ')}${String(bonds[i].order).padStart(3, ' ')}  0  0  0  0\n`;
    }
    sdf += 'M  END\n$$$$\n';
    return sdf;
};

window.WB.renderMoleculeModel = function(viewer, xyzText) {
    if (!viewer || !xyzText) return;
    const { symbols, coords } = window.WB.parseXyz(xyzText);
    if (symbols.length > 0) {
        const sdf = window.WB.buildMolBlock(symbols, coords);
        const m = viewer.addModel(sdf, 'sdf');
        m.setStyle({}, {
            stick: { radius: 0.13, colorscheme: 'Jmol' },
            sphere: { scale: 0.28, colorscheme: 'Jmol' }
        });
    } else {
        viewer.addModel(xyzText, 'xyz');
        viewer.setStyle({}, {
            stick: { radius: 0.13, colorscheme: 'Jmol' },
            sphere: { scale: 0.28, colorscheme: 'Jmol' }
        });
    }
};

window.WB.getOrbitalPhaseColors = function(fileOrKey) {
    const fKey = String(fileOrKey || '').toLowerCase();
    const isAnti = fKey.includes('star') || fKey.includes('*') || fKey.includes('lumo') || fKey.includes('antibond') || fKey.includes('virtual');
    return isAnti
        ? { pos: '#fb923c', neg: '#db2777' }
        : { pos: '#38bdf8', neg: '#1d4ed8' };
};

window.WB.shiftColor = function(hex) {
    return (hex && hex.toLowerCase() === '#fb923c') ? '#db2777' : '#1d4ed8';
};

window.WB.lighten = function(hex) {
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
                if (name) name = decodeURIComponent(name);

                try {
                    const res = await fetch('/api/orbital/remesh', {
                        method: 'POST',
                        headers: {'Content-Type': 'application/json'},
                        body: JSON.stringify({ name: name, cube_file: file, isovalue: window._customIsovalue })
                    });
                    if (res.ok) {
                        const data = await res.json();
                        if (!data.error && data.pos && data.neg) {
                            window._currentOrbitalCache[file.replace('.json', '_pos.json')] = data.pos;
                            window._currentOrbitalCache[file.replace('.json', '_neg.json')] = data.neg;
                            mesh = window._currentOrbitalCache[actualFile];
                        }
                    }
                } catch (err) {
                    // remesh request failed, will fallback to bundle mesh below
                }
            }
            if (!mesh) {
                // Fetch from bundle (fallback or default)
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
        if (!file) continue;
        const pc = window.WB.getOrbitalPhaseColors(file);
        
        if (file.endsWith('.json')) {
            for (const sign of ['pos', 'neg']) {
                const mesh = await window.WB.fetchOrbitalMesh(file, sign, base);
                const verts = mesh && (mesh.vertices || mesh.vertexArr);
                const faces = mesh && (mesh.faces || mesh.faceArr);
                const norms = mesh && (mesh.normals || mesh.normalArr);
                if (verts && faces) {
                    const finalColor = (sign === 'neg') ? pc.neg : pc.pos;
                    shapesToAdd.push({
                        vertexArr: verts,
                        faceArr: faces,
                        normalArr: norms || undefined,
                        color: $3Dmol.CC.color(finalColor),
                        opacity: 0.75
                    });
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
            const iso = window._customIsovalue || 0.035;
            shapesToAdd.push({vol: vol, isoval: iso, color: pc.pos, opacity: 0.75, smoothness: 1});
            shapesToAdd.push({vol: vol, isoval: -iso, color: pc.neg, opacity: 0.75, smoothness: 1});
        }
    }

    if (renderId !== currentRenderId) return; // Abort if a newer render started
    
    const viewer = window.WB.viewer;
    if (!viewer) return;
    
    viewer.removeAllModels(); 
    viewer.removeAllShapes();
    viewer.removeAllSurfaces();
    
    if (window._currentMolXyz) {
        window.WB.renderMoleculeModel(viewer, window._currentMolXyz);
    }
    
    // Check ESP surface
    const espCb = sidebar.querySelector('input[data-esp]:checked');
    const colorbar = document.getElementById('esp-colorbar');
    if (espCb && window._currentEspCache && window._currentEspCache.vertices) {
        const espCache = window._currentEspCache;
        const espMin = espCache.esp_min;
        const espMax = espCache.esp_max;
        const absMax = Math.max(Math.abs(espMin || -30), Math.abs(espMax || 30), 20);
        if (colorbar) {
            colorbar.style.display = 'flex';
            const minEl = document.getElementById('esp-colorbar-min');
            const maxEl = document.getElementById('esp-colorbar-max');
            if (minEl) minEl.textContent = `${(-absMax).toFixed(1)}`;
            if (maxEl) maxEl.textContent = `+${absMax.toFixed(1)}`;
        }
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
            opacity: 0.65,
        });
    } else if (colorbar) {
        colorbar.style.display = 'none';
    }

    // Check Atom Charges
    const chargesCb = sidebar.querySelector('input[data-charges]:checked');
    if (chargesCb && window._currentMolXyz) {
        const lines = window._currentMolXyz.trim().split('\n');
        if (lines.length >= 3) {
            const charges = window._currentMolCharges || [];
            for (let i = 2; i < lines.length; i++) {
                const parts = lines[i].trim().split(/\s+/);
                if (parts.length >= 4) {
                    const idx = i - 2;
                    let q = charges[idx];
                    const sym = parts[0];
                    if (q === undefined) {
                        if (sym === 'H') q = 0.06;
                        else if (sym === 'C') q = -0.06;
                        else if (sym === 'O') q = -0.50;
                        else if (sym === 'N') q = -0.40;
                        else if (sym === 'Cl') q = -0.20;
                        else q = 0.0;
                    }
                    const sign = q > 0 ? '+' : '';
                    viewer.addLabel(`${sym}: ${sign}${Number(q).toFixed(2)}`, {
                        position: { x: parseFloat(parts[1]), y: parseFloat(parts[2]) + 0.35, z: parseFloat(parts[3]) },
                        fontSize: 10,
                        fontColor: '#ffffff',
                        backgroundColor: 'rgba(15, 23, 42, 0.75)',
                        backgroundOpacity: 0.75
                    });
                }
            }
        }
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
