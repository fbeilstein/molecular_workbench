/**
 * 3Dmol viewer and isosurface renderer for Step 5 Bundle Viewer.
 * Preserves camera zoom/orientation across frames; provides zoom & scale controls.
 */
import { getEspColor, getPhaseColors } from '../orbitals/orbitals_types.js';

const COV_RADII = {
    H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
    P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39
};
const MAX_VALENCE = { H: 1, F: 1, Cl: 1, Br: 1, I: 1, O: 3, N: 4, C: 4, S: 6, P: 5 };

function computeBonds(symbols, coords) {
    const bonds = [];
    const n = Math.min(symbols.length, coords.length);
    const valences = new Array(n).fill(0);
    const candidates = [];
    const adj = Array.from({ length: n }, () => []);
    const distMat = {};

    for (let i = 0; i < n; i++) {
        const s1 = symbols[i];
        const r1 = COV_RADII[s1] || 0.75;
        for (let j = i + 1; j < n; j++) {
            const s2 = symbols[j];
            const r2 = COV_RADII[s2] || 0.75;
            if (s1 === 'H' && s2 === 'H' && n > 2) continue;

            const dx = coords[i][0] - coords[j][0];
            const dy = coords[i][1] - coords[j][1];
            const dz = coords[i][2] - coords[j][2];
            const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
            const maxD = 1.28 * (r1 + r2);

            if (d >= 0.4 && d <= maxD) {
                adj[i].push(j);
                adj[j].push(i);
                const key = `${i}_${j}`;
                distMat[key] = d;
                candidates.push({ i, j, s1, s2, d });
            }
        }
    }

    const getDist = (u, v) => {
        const key = u < v ? `${u}_${v}` : `${v}_${u}`;
        return distMat[key] || 1.4;
    };

    // 1. Detect 6-membered aromatic/conjugated rings (primarily C, N)
    const ringDoubleBonds = new Set();
    const ringBondsAll = new Set();
    const ringAtoms = [];
    for (let i = 0; i < n; i++) {
        if (symbols[i] === 'C' || symbols[i] === 'N') {
            ringAtoms.push(i);
        }
    }

    const cycles = [];
    const visited = new Set();
    const path = [];

    const findCycles = (v, start, depth) => {
        if (depth === 6) {
            if (adj[v].includes(start)) {
                cycles.push([...path]);
            }
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

    for (const startAtom of ringAtoms) {
        visited.clear();
        visited.add(startAtom);
        path.length = 0;
        path.push(startAtom);
        findCycles(startAtom, startAtom, 1);
    }

    const uniqueRings = [];
    const seenRings = new Set();
    for (const c of cycles) {
        const sorted = [...c].sort((a, b) => a - b).join(',');
        if (!seenRings.has(sorted)) {
            seenRings.add(sorted);
            uniqueRings.push(c);
        }
    }

    // Kekulize each 6-membered ring:
    // A 6-cycle has two alternating matchings: edges0 and edges1.
    // Pick matching with smaller total bond length sum.
    for (const ring of uniqueRings) {
        let isConjugated = true;
        for (let k = 0; k < 6; k++) {
            const u = ring[k];
            const v = ring[(k + 1) % 6];
            if (getDist(u, v) > 1.48) {
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
        for (const [u, v] of chosen) {
            const key = u < v ? `${u}_${v}` : `${v}_${u}`;
            ringDoubleBonds.add(key);
        }
        for (let k = 0; k < 6; k++) {
            const u = ring[k];
            const v = ring[(k + 1) % 6];
            const key = u < v ? `${u}_${v}` : `${v}_${u}`;
            ringBondsAll.add(key);
        }
    }

    for (const c of candidates) {
        const pairKey = c.i < c.j ? `${c.i}_${c.j}` : `${c.j}_${c.i}`;
        const pair = [c.s1, c.s2].sort().join('-');
        let targetOrder = 1;

        if (ringDoubleBonds.has(pairKey)) {
            targetOrder = 2;
        } else if (!ringBondsAll.has(pairKey)) {
            if (pair === 'C-O') {
                if (c.d < 1.16) targetOrder = 3;
                else if (c.d < 1.265) targetOrder = 2;
            } else if (pair === 'C-C') {
                if (c.d < 1.24) targetOrder = 3;
                else if (c.d < 1.39) targetOrder = 2;
            } else if (pair === 'C-N') {
                if (c.d < 1.18) targetOrder = 3;
                else if (c.d < 1.32) targetOrder = 2;
            } else if (pair === 'N-N') {
                if (c.d < 1.14) targetOrder = 3;
                else if (c.d < 1.28) targetOrder = 2;
            } else if (pair === 'N-O' || pair === 'O-O') {
                if (c.d < 1.25) targetOrder = 2;
            }
        }

        bonds.push({
            a1: c.i, a2: c.j, order: 1, target: targetOrder,
            s1: c.s1, s2: c.s2, d: c.d, pairKey
        });
        valences[c.i]++;
        valences[c.j]++;
    }

    // Priority sort: triple bonds first (target 3), then ring double bonds, then other double bonds, then by length
    bonds.sort((a, b) => {
        const prioA = a.target === 3 ? 0 : (ringDoubleBonds.has(a.pairKey) ? 1 : (a.target === 2 ? 2 : 3));
        const prioB = b.target === 3 ? 0 : (ringDoubleBonds.has(b.pairKey) ? 1 : (b.target === 2 ? 2 : 3));
        if (prioA !== prioB) return prioA - prioB;
        return a.d - b.d;
    });

    // Promote to double/triple bond if atom valences allow
    for (const b of bonds) {
        if (b.target > 1) {
            const v1Cap = MAX_VALENCE[b.s1] || 4;
            const v2Cap = MAX_VALENCE[b.s2] || 4;
            const desiredBoost = b.target - 1;
            const maxBoost = Math.min(v1Cap - valences[b.a1], v2Cap - valences[b.a2]);
            const boost = Math.max(0, Math.min(desiredBoost, maxBoost));
            if (boost > 0) {
                b.order += boost;
                valences[b.a1] += boost;
                valences[b.a2] += boost;
            }
        }
    }

    return bonds;
}

export class Bundle3DViewer {
    constructor(containerEl, options = {}) {
        this.containerEl = containerEl;
        this.statusEl = (options && options.nodeType) ? options : (options && options.statusEl) || null;
        this.viewer = null;
        this.renderSeq = 0;
        this.meshCache = {};
        this.cubeCache = {};
        this.vdoCache = {};
        this.espMeshCache = {};
        this._espInFlight = new Set();
        this.activeFrameIdx = 0;
        this._isPrewarming = false;
        this.hasInitializedCamera = false;
        this.moleculeScale = 1.0;
        this.init();
    }

    clearCache() {
        this.meshCache = {};
        this.cubeCache = {};
        this.vdoCache = {};
        this.espMeshCache = {};
        this._espInFlight.clear();
        this._isPrewarming = false;
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            this.containerEl.innerHTML = '<div class="empty-state">3Dmol.js not loaded.</div>';
            return;
        }
        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: '#090d16',
            defaultcolors: $3Dmol.rasmolElementColors
        });

        this.containerEl.addEventListener('wheel', (e) => {
            e.preventDefault();
            const factor = e.deltaY < 0 ? 1.12 : 0.89;
            this.zoom(factor);
        }, { passive: false });
    }

    zoom(factor) {
        if (!this.viewer) return;
        this.viewer.zoom(factor);
        this.viewer.render();
    }

    zoomIn() {
        this.zoom(1.2);
    }

    zoomOut() {
        this.zoom(0.83);
    }

    resetCamera() {
        if (!this.viewer) return;
        this.viewer.zoomTo();
        this.viewer.render();
    }

    setMoleculeScale(scale) {
        this.moleculeScale = Math.max(0.4, Math.min(3.0, scale));
        if (this.viewer) {
            this.viewer.setStyle({}, {
                stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
                sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
            });
            this.viewer.render();
        }
    }

    setStatus(text, isError = false) {
        if (!this.statusEl || !this.statusEl.style) return;
        this.statusEl.textContent = text;
        this.statusEl.style.color = isError ? '#f87171' : '#94a3b8';
    }

    renderMolecule(symbols, coords, charges = null, showCharges = false, showEspSurface = false, frameIdx = 0) {
        if (!this.viewer || !symbols || !coords) return;
        this.activeFrameIdx = frameIdx;
        this.viewer.removeAllSurfaces();
        this.viewer.removeAllShapes();
        this.viewer.clear();

        const natoms = Math.min(symbols.length, coords.length);
        const n = natoms;
        const bonds = computeBonds(symbols, coords);
        let sdf = `\n  Reaction Pathway Viewer\n\n${String(natoms).padStart(3, ' ')}${String(bonds.length).padStart(3, ' ')}  0  0  0  0  0  0  0  0999 V2000\n`;
        for (let i = 0; i < natoms; i++) {
            const sym = symbols[i];
            const pos = coords[i];
            const x = Number(pos[0]).toFixed(4).padStart(10, ' ');
            const y = Number(pos[1]).toFixed(4).padStart(10, ' ');
            const z = Number(pos[2]).toFixed(4).padStart(10, ' ');
            sdf += `${x}${y}${z} ${sym.padEnd(3, ' ')} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
        }
        for (let i = 0; i < bonds.length; i++) {
            const b = bonds[i];
            sdf += `${String(b.a1 + 1).padStart(3, ' ')}${String(b.a2 + 1).padStart(3, ' ')}${String(b.order).padStart(3, ' ')}  0  0  0  0\n`;
        }
        sdf += 'M  END\n$$$$\n';

        this.viewer.addModel(sdf, 'sdf');
        this.viewer.setStyle({}, {
            stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
            sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
        });

        // 1. ESP Potential Surface (cached WebGL mesh or client VDW extraction)
        if (showEspSurface) {
            const actualCharges = (charges && charges.length >= n) ? charges : symbols.map(s => {
                if (s === 'H') return 0.08;
                if (s === 'C') return -0.05;
                if (s === 'O') return -0.50;
                if (s === 'Cl') return -0.45;
                return 0.0;
            });

            const cachedEsp = this.espMeshCache[frameIdx];
            if (cachedEsp && cachedEsp.vertexArr && cachedEsp.faceArr) {
                this.viewer.addCustom({
                    vertexArr: cachedEsp.vertexArr,
                    faceArr: cachedEsp.faceArr,
                    color: cachedEsp.color,
                    opacity: 0.65
                });
            } else if (!this._espInFlight.has(frameIdx)) {
                this._espInFlight.add(frameIdx);
                const atoms = this.viewer.selectedAtoms({});
                for (let i = 0; i < atoms.length && i < actualCharges.length; i++) {
                    const qVal = actualCharges[i] || 0.0;
                    atoms[i].partialCharge = qVal;
                    atoms[i].properties = atoms[i].properties || {};
                    atoms[i].properties.partialCharge = qVal;
                }
                const p = this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: 0.65,
                    colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
                });
                if (p && typeof p.then === 'function') {
                    p.then(sId => {
                        this._espInFlight.delete(frameIdx);
                        const sObj = this.viewer.surfaces ? this.viewer.surfaces[sId] : null;
                        const grp = sObj?.geo?.geometryGroups?.[0];
                        if (grp && grp.vertexArray && grp.faceArray) {
                            const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
                            const verts = [];
                            for (let vi = 0; vi < vCount; vi++) {
                                verts.push({ x: grp.vertexArray[vi * 3], y: grp.vertexArray[vi * 3 + 1], z: grp.vertexArray[vi * 3 + 2] });
                            }
                            const faces = Array.from(grp.faceArray.slice(0, grp.faces * 3));
                            const cols = [];
                            if (grp.colorArray && grp.colorArray.length >= vCount * 3) {
                                for (let ci = 0; ci < vCount; ci++) {
                                    cols.push({ r: grp.colorArray[ci * 3], g: grp.colorArray[ci * 3 + 1], b: grp.colorArray[ci * 3 + 2] });
                                }
                            }
                            this.espMeshCache[frameIdx] = {
                                vertexArr: verts,
                                faceArr: faces,
                                color: cols.length ? cols : undefined
                            };
                        }
                        if (this.activeFrameIdx === frameIdx) {
                            this.viewer.render();
                        } else {
                            try { this.viewer.removeSurface(sId); } catch (e) {}
                        }
                    }).catch(() => { this._espInFlight.delete(frameIdx); });
                } else {
                    this._espInFlight.delete(frameIdx);
                    if (this.activeFrameIdx === frameIdx) {
                        this.viewer.render();
                    }
                }
            }
        }

        // 2. Atom Numeric Partial Charge Labels
        if (showCharges) {
            const actualCharges = (charges && charges.length >= n) ? charges : symbols.map(s => {
                if (s === 'H') return 0.08;
                if (s === 'C') return -0.05;
                if (s === 'O') return -0.50;
                if (s === 'Cl') return -0.45;
                return 0.0;
            });
            for (let i = 0; i < n; i++) {
                const q = actualCharges[i];
                const [x, y, z] = coords[i];
                this.viewer.addLabel(`${symbols[i]}: ${q > 0 ? '+' : ''}${Number(q).toFixed(2)}`, {
                    position: { x, y: y + 0.35, z },
                    fontSize: 10,
                    fontColor: '#ffffff',
                    backgroundColor: 'rgba(15, 23, 42, 0.75)',
                    backgroundOpacity: 0.75
                });
            }
        }
        this.viewer.render();
    }

    async prewarmEspCache(symbols, frames, trajectoryCharges, onProgress) {
        if (this._isPrewarming || !symbols || !frames || !frames.length) return;
        this._isPrewarming = true;
        try {
            if (!this._bgDiv) {
                this._bgDiv = document.createElement('div');
                this._bgDiv.style.cssText = 'position:fixed;left:-9999px;width:120px;height:120px;opacity:0;pointer-events:none;';
                document.body.appendChild(this._bgDiv);
                this._bgViewer = $3Dmol.createViewer(this._bgDiv);
            }
            const bg = this._bgViewer;
            const total = frames.length, n = symbols.length;

            for (let f = 0; f < total; f++) {
                if (!this._isPrewarming) break;
                if (this.espMeshCache[f]) continue;
                const coords = frames[f]?.coords || frames[f]?.positions;
                const charges = trajectoryCharges?.[f] || frames[f]?.charges;
                if (!coords || !charges || charges.length < n) continue;

                bg.clear();
                let mdl = `${n}\nFrame Model\n`;
                for (let i = 0; i < n; i++) {
                    const [x, y, z] = coords[i];
                    mdl += `${symbols[i]} ${x.toFixed(4)} ${y.toFixed(4)} ${z.toFixed(4)}\n`;
                }
                bg.addModel(mdl, 'xyz');
                const atoms = bg.selectedAtoms({});
                for (let i = 0; i < atoms.length && i < charges.length; i++) {
                    atoms[i].partialCharge = charges[i] || 0.0;
                    atoms[i].properties = { partialCharge: charges[i] || 0.0 };
                }
                const p = bg.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: 0.65,
                    colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
                });
                if (p && typeof p.then === 'function') {
                    const sId = await p;
                    const grp = bg.surfaces?.[sId]?.geo?.geometryGroups?.[0];
                    if (grp && grp.vertexArray && grp.faceArray) {
                        const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
                        const verts = [];
                        for (let vi = 0; vi < vCount; vi++) {
                            verts.push({ x: grp.vertexArray[vi * 3], y: grp.vertexArray[vi * 3 + 1], z: grp.vertexArray[vi * 3 + 2] });
                        }
                        const faces = Array.from(grp.faceArray.slice(0, grp.faces * 3));
                        const cols = [];
                        if (grp.colorArray && grp.colorArray.length >= vCount * 3) {
                            for (let ci = 0; ci < vCount; ci++) {
                                cols.push({ r: grp.colorArray[ci * 3], g: grp.colorArray[ci * 3 + 1], b: grp.colorArray[ci * 3 + 2] });
                            }
                        }
                        this.espMeshCache[f] = { vertexArr: verts, faceArr: faces, color: cols.length ? cols : undefined };
                    }
                    try { bg.removeSurface(sId); } catch (e) {}
                }
                if (onProgress) onProgress(f + 1, total);
            }
        } catch (e) {
            console.warn('[Bundle3DViewer] ESP prewarm error:', e);
        } finally {
            this._isPrewarming = false;
        }
    }

    _drawOrbitalLayer(orb, frameIdx, isovalue) {
        const orbKey = orb.orbital_key || orb.key || orb;
        if (orbKey === 'track_esp_charge') return;
        const pc = getPhaseColors(orb.color || '#06b6d4');
        const mKey = `${frameIdx}_${orbKey}`;
        const isDefaultIso = Math.abs(isovalue - 0.03) < 0.003;
        const isoKey = `${mKey}_${Math.round(isovalue * 1000)}`;

        // 1. If we have volumetric cube data and isovalue is not default 0.03 (or mesh not cached),
        // use 3Dmol's native VolumeData marching cubes for instant arbitrary isovalue rendering:
        const cubeData = this.cubeCache[mKey];
        if (cubeData && (!isDefaultIso || !this.meshCache[mKey])) {
            let vdo = this.vdoCache[mKey];
            if (!vdo) {
                try {
                    vdo = new $3Dmol.VolumeData(cubeData, 'cube');
                    this.vdoCache[mKey] = vdo;
                } catch (e) {}
            }
            if (vdo) {
                try {
                    this.viewer.addIsosurface(vdo, {
                        isoval: Math.abs(isovalue),
                        color: pc.pos,
                        opacity: 0.72,
                        wireframe: false
                    });
                    this.viewer.addIsosurface(vdo, {
                        isoval: -Math.abs(isovalue),
                        color: pc.neg,
                        opacity: 0.72,
                        wireframe: false
                    });
                    return;
                } catch (e) {}
            }
        }

        // 2. Precalculated 3D surface mesh
        const mesh = this.meshCache[isoKey] || this.meshCache[mKey];
        if (mesh) {
            if (mesh.pos && mesh.pos.vertices && mesh.pos.faces) {
                this.viewer.addCustom({
                    vertexArr: mesh.pos.vertices,
                    faceArr: mesh.pos.faces,
                    normalArr: mesh.pos.normals || undefined,
                    color: $3Dmol.CC.color(pc.pos),
                    opacity: 0.72
                });
            }
            if (mesh.neg && mesh.neg.vertices && mesh.neg.faces) {
                this.viewer.addCustom({
                    vertexArr: mesh.neg.vertices,
                    faceArr: mesh.neg.faces,
                    normalArr: mesh.neg.normals || undefined,
                    color: $3Dmol.CC.color(pc.neg),
                    opacity: 0.72
                });
            }
            return;
        }

        // 3. Fallback to cubeData
        if (cubeData) {
            let vdo = this.vdoCache[mKey];
            if (!vdo) {
                try {
                    vdo = new $3Dmol.VolumeData(cubeData, 'cube');
                    this.vdoCache[mKey] = vdo;
                } catch (e) {}
            }
            if (vdo) {
                try {
                    this.viewer.addIsosurface(vdo, {
                        isoval: Math.abs(isovalue),
                        color: pc.pos,
                        opacity: 0.72,
                        wireframe: false
                    });
                    this.viewer.addIsosurface(vdo, {
                        isoval: -Math.abs(isovalue),
                        color: pc.neg,
                        opacity: 0.72,
                        wireframe: false
                    });
                } catch (e) {}
            }
        }
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false, showEspSurface = false) {
        if (!this.viewer || !symbols || !coords) return;
        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx);

        if (!this.hasInitializedCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) {
                this._drawOrbitalLayer(orb, frameIdx, isovalue);
            }
        }
        this.viewer.render();
    }

    async renderScene(options) {
        const {
            symbols, coords, charges, showCharges, showEspSurface,
            activeOrbitals, isovalue, frameIdx,
            hideAllOrbitals, fetchCubeFn, resetCamera
        } = options;

        if (!this.viewer) return;
        const currentSeq = ++this.renderSeq;

        // 1. Fetch missing orbital surfaces/cubes in parallel (skip ESP pseudo-track)
        const missing = [];
        if (!hideAllOrbitals && activeOrbitals && activeOrbitals.length > 0 && fetchCubeFn) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                if (orbKey === 'track_esp_charge') continue;
                const cacheKey = `${frameIdx}_${orbKey}`;
                const isDefaultIso = Math.abs(isovalue - 0.03) < 0.003;
                const isoKey = `${cacheKey}_${Math.round(isovalue * 1000)}`;
                const hasSurface = this.cubeCache[cacheKey] || this.meshCache[isoKey] || this.meshCache[cacheKey];

                if (!hasSurface) {
                    missing.push(fetchCubeFn(frameIdx, orbKey, orb.chem_key || null, isovalue).then(res => {
                        if (res && res.success) {
                            if (res.cube) {
                                this.cubeCache[cacheKey] = res.cube;
                            }
                            if (res.pos || res.neg) {
                                const meshData = { pos: res.pos, neg: res.neg };
                                this.meshCache[cacheKey] = meshData;
                                this.meshCache[isoKey] = meshData;
                            }
                        } else {
                            this.meshCache[isoKey] = null;
                            this.meshCache[cacheKey] = null;
                        }
                    }).catch(() => {
                        this.meshCache[isoKey] = null;
                        this.meshCache[cacheKey] = null;
                    }));
                }
            }
        }

        if (missing.length > 0) {
            this.setStatus(`Frame ${frameIdx} | Loading Orbitals...`);
            await Promise.all(missing);
            if (this.renderSeq !== currentSeq) return;
        }

        // 2. Render atoms and electrostatic surface
        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx);

        // 3. ZoomTo on first render or explicit reset
        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.viewer.render();
            this.hasInitializedCamera = true;
        }

        // 4. Render orbital layers (prefers surfaces)
        if (!hideAllOrbitals && activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) {
                this._drawOrbitalLayer(orb, frameIdx, isovalue);
            }
        }

        if (this.renderSeq === currentSeq) {
            this.viewer.render();
            this.setStatus(`Frame ${frameIdx} | Ready`);
        }
    }
}
