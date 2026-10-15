/**
 * 3Dmol scene and isosurface renderer for Reaction Orbital Studio.
 * Preserves camera zoom/orientation across frames; provides zoom & scale controls.
 */
import { getEspColor } from './orbitals_types.js';

export class Orbitals3DViewer {
    constructor(containerEl, statusEl) {
        this.containerEl = containerEl;
        this.statusEl = statusEl;
        this.viewer = null;
        this.renderSeq = 0;
        this.cubeCache = {};
        this.vdoCache = {};
        this.hasInitializedCamera = false;
        this.moleculeScale = 1.0;
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            this.containerEl.innerHTML = '<div class="empty-state">3Dmol.js library not loaded.</div>';
            return;
        }
        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: '#0f172a',
            id: 'orbitals-3dmol-canvas'
        });

        // Ensure mouse wheel reliably scales the 3D molecule without page scroll interference
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
        if (!this.statusEl) return;
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
            stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
            sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
        });

        if (showCharges && charges && charges.length >= n) {
            const atoms = this.viewer.selectedAtoms({});
            for (let i = 0; i < atoms.length && i < charges.length; i++) {
                atoms[i].partialCharge = charges[i] || 0.0;
            }
            this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.72,
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

    renderMoleculeFast(symbols, coords) {
        if (!this.viewer || !symbols || !coords) return;
        this.renderMolecule(symbols, coords, null, false);
        if (!this.hasInitializedCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false) {
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
            stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
            sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
        });

        if (showCharges && charges && charges.length >= n) {
            const atoms = this.viewer.selectedAtoms({});
            for (let i = 0; i < atoms.length && i < charges.length; i++) {
                atoms[i].partialCharge = charges[i] || 0.0;
            }
            this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.72,
                colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
            });
        }

        if (activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                const cubeData = this.cubeCache[`${frameIdx}_${orbKey}`];
                if (cubeData) {
                    const vKey = `${frameIdx}_${orbKey}`;
                    let vdo = this.vdoCache[vKey];
                    if (!vdo) {
                        try {
                            vdo = new $3Dmol.VolumeData(cubeData, 'cube');
                            this.vdoCache[vKey] = vdo;
                        } catch (e) {}
                    }
                    if (vdo) {
                        try {
                            this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: orb.color || '#06b6d4', opacity: 0.72 });
                            this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: '#f43f5e', opacity: 0.72 });
                        } catch (e) {}
                    }
                }
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

        this.renderMolecule(symbols, coords, charges, showCharges);

        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (hideAllOrbitals || !activeOrbitals || activeOrbitals.length === 0) {
            this.viewer.removeAllSurfaces();
            this.viewer.render();
            this.setStatus(`Frame ${frameIdx} | Ready`);
            return;
        }

        this.setStatus(`Frame ${frameIdx} | Loading Orbitals...`);

        const missing = [];
        for (const orb of activeOrbitals) {
            const orbKey = orb.orbital_key || orb.key || orb;
            const chemKey = orb.chem_key || null;
            const cacheKey = `${frameIdx}_${orbKey}`;
            if (!this.cubeCache[cacheKey] && fetchCubeFn) {
                missing.push(fetchCubeFn(frameIdx, orbKey, chemKey, isovalue).then(res => {
                    if (res && res.success && res.cube) this.cubeCache[cacheKey] = res.cube;
                }).catch(() => {}));
            }
        }

        if (missing.length > 0) {
            await Promise.all(missing);
            if (this.renderSeq !== currentSeq) return;
        }

        this.viewer.removeAllSurfaces();
        for (const orb of activeOrbitals) {
            const orbKey = orb.orbital_key || orb.key || orb;
            const color = orb.color || '#06b6d4';
            const cubeData = this.cubeCache[`${frameIdx}_${orbKey}`];
            if (cubeData) {
                const vKey = `${frameIdx}_${orbKey}`;
                let vdo = this.vdoCache[vKey];
                if (!vdo) {
                    try {
                        vdo = new $3Dmol.VolumeData(cubeData, 'cube');
                        this.vdoCache[vKey] = vdo;
                    } catch (err) {}
                }
                if (vdo) {
                    try {
                        this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: color, opacity: 0.72 });
                        this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: '#f43f5e', opacity: 0.72 });
                    } catch (err) {}
                }
            }
        }

        if (this.renderSeq === currentSeq) {
            this.viewer.render();
            this.setStatus(`Frame ${frameIdx} | Ready`);
        }
    }

    clearCache() {
        this.cubeCache = {};
        this.vdoCache = {};
    }
}
