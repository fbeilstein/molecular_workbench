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
        this.hasInitializedCamera = false;
        this.moleculeScale = 1.0;
        this.init();
    }

    clearCache() {
        this.meshCache = {};
        this.cubeCache = {};
        this.vdoCache = {};
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

    renderMolecule(symbols, coords, charges = null, showCharges = false) {
        if (!this.viewer || !symbols || !coords) return;
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

        if (showCharges && charges && charges.length >= n) {
            const atoms = this.viewer.selectedAtoms({});
            for (let i = 0; i < atoms.length && i < charges.length; i++) {
                atoms[i].partialCharge = charges[i] || 0.0;
                atoms[i].properties = atoms[i].properties || {};
                atoms[i].properties.partialCharge = charges[i] || 0.0;
            }
            this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.65,
                colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
            });
            for (let i = 0; i < n; i++) {
                const q = charges[i];
                const [x, y, z] = coords[i];
                this.viewer.addLabel(`${symbols[i]}: ${q > 0 ? '+' : ''}${q.toFixed(2)}`, {
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

    _drawOrbitalLayer(orb, frameIdx, isovalue) {
        const orbKey = orb.orbital_key || orb.key || orb;
        if (orbKey === 'track_esp_charge') return;
        const pc = getPhaseColors(orb.color || '#06b6d4');
        const mKey = `${frameIdx}_${orbKey}`;

        // 1. Prefer precalculated 3D surface mesh (Direct WebGL vertex/face buffers)
        const mesh = this.meshCache[mKey];
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

        // 2. Fallback to volumetric cube isosurface
        const cubeData = this.cubeCache[mKey];
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

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false) {
        if (!this.viewer || !symbols || !coords) return;
        this.renderMolecule(symbols, coords, charges, showCharges);

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
            symbols, coords, charges, showCharges,
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
                if (!this.meshCache[cacheKey] && !this.cubeCache[cacheKey]) {
                    missing.push(fetchCubeFn(frameIdx, orbKey, orb.chem_key || null, isovalue).then(res => {
                        if (res && res.success) {
                            if (res.pos || res.neg) {
                                this.meshCache[cacheKey] = { pos: res.pos, neg: res.neg };
                            } else if (res.cube) {
                                this.cubeCache[cacheKey] = res.cube;
                            }
                        } else {
                            this.meshCache[cacheKey] = null;
                        }
                    }).catch(() => {
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
        this.renderMolecule(symbols, coords, charges, showCharges);

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
