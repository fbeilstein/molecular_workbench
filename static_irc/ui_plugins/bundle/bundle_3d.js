/**
 * 3Dmol viewer and isosurface renderer for Step 5 Bundle Viewer.
 * Preserves camera zoom/orientation across frames; provides zoom & scale controls.
 */
import { getEspColor, getPhaseColors } from '../orbitals/orbitals_types.js';

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
                stick: { radius: 0.14 * this.moleculeScale },
                sphere: { scale: 0.28 * this.moleculeScale }
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
        this.viewer.clear();

        const n = Math.min(symbols.length, coords.length);
        let mdl = `${n}\nFrame Model\n`;
        for (let i = 0; i < n; i++) {
            const sym = symbols[i];
            const [x, y, z] = coords[i];
            mdl += `${sym} ${x.toFixed(4)} ${y.toFixed(4)} ${z.toFixed(4)}\n`;
        }

        this.viewer.addModel(mdl, 'xyz');
        this.viewer.setStyle({}, {
            stick: { radius: 0.14 * this.moleculeScale },
            sphere: { scale: 0.28 * this.moleculeScale }
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
        const mesh = isDefaultIso ? (this.meshCache[isoKey] || this.meshCache[mKey]) : this.meshCache[isoKey];
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
                const hasSurface = this.cubeCache[cacheKey] || (isDefaultIso ? this.meshCache[cacheKey] : this.meshCache[isoKey]);

                if (!hasSurface) {
                    missing.push(fetchCubeFn(frameIdx, orbKey, orb.chem_key || null, isovalue).then(res => {
                        if (res && res.success) {
                            if (res.cube) {
                                this.cubeCache[cacheKey] = res.cube;
                            }
                            if (res.pos || res.neg) {
                                if (isDefaultIso) {
                                    this.meshCache[cacheKey] = { pos: res.pos, neg: res.neg };
                                } else {
                                    this.meshCache[isoKey] = { pos: res.pos, neg: res.neg };
                                }
                            }
                        } else {
                            if (!isDefaultIso) this.meshCache[isoKey] = null;
                            else this.meshCache[cacheKey] = null;
                        }
                    }).catch(() => {
                        if (!isDefaultIso) this.meshCache[isoKey] = null;
                        else this.meshCache[cacheKey] = null;
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
