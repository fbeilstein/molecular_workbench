/**
 * 3Dmol scene and volumetric isosurface renderer for Step 6 Orbitals.
 * Preserves camera orientation and zoom across frames, provides instant fast scrubbing,
 * and renders distinct phase colors for bonding vs. anti-bonding orbitals.
 */
import { getEspColor, getPhaseColors } from './orbitals_types.js';

export class Orbitals3DViewer {
    constructor(containerEl, statusEl) {
        this.containerEl = containerEl;
        this.statusEl = statusEl;
        this.viewer = null;
        this.cubeCache = {};
        this.vdoCache = {};
        this.hasInitializedCamera = false;
        this.moleculeScale = 1.0;
        this.currentFrameIdx = 0;
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            this.containerEl.innerHTML = '<div style="padding: 1rem; color: #ef4444;">3Dmol.js library not loaded.</div>';
            return;
        }
        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: '#0f172a',
            id: 'orbitals-3dmol-canvas'
        });

        // Mouse wheel scales 3D molecule without page scroll interference
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

    renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx) {
        if (!this.viewer || !activeOrbitals) return;
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
                    } catch (e) {
                        console.error('VolumeData error:', e);
                    }
                }
                if (vdo) {
                    try {
                        const pc = getPhaseColors(orb);
                        this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: pc.pos, opacity: 0.72 });
                        this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: pc.neg, opacity: 0.72 });
                    } catch (e) {
                        console.error('addIsosurface error:', e);
                    }
                }
            }
        }
        this.viewer.render();
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false) {
        this.currentFrameIdx = frameIdx;
        this.renderMolecule(symbols, coords, charges, showCharges);

        if (!this.hasInitializedCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (activeOrbitals && activeOrbitals.length > 0) {
            this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
        }
    }

    async renderScene(options) {
        const {
            symbols, coords, charges, showCharges,
            activeOrbitals, isovalue, frameIdx,
            hideAllOrbitals, fetchCubeFn, resetCamera
        } = options;

        if (!this.viewer) return;
        this.currentFrameIdx = frameIdx;

        // 1. Immediately render molecular coordinates
        this.renderMolecule(symbols, coords, charges, showCharges);

        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (hideAllOrbitals || !activeOrbitals || activeOrbitals.length === 0) {
            this.setStatus(`Frame ${frameIdx} | Ready`);
            return;
        }

        // 2. Render any already-cached orbital surfaces immediately
        this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);

        // 3. Find missing orbital cubes
        const missing = [];
        if (fetchCubeFn) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                const cacheKey = `${frameIdx}_${orbKey}`;
                if (!this.cubeCache[cacheKey]) {
                    missing.push({ orb, orbKey, cacheKey });
                }
            }
        }

        if (missing.length > 0) {
            this.setStatus(`Frame ${frameIdx} | Computing orbital (${missing.map(m => m.orbKey).join(', ')})...`);
            const fetchErrors = [];
            try {
                await Promise.all(missing.map(async (m) => {
                    const res = await fetchCubeFn(frameIdx, m.orbKey, m.orb.chem_key || null, isovalue);
                    if (res && res.success && res.cube) {
                        this.cubeCache[m.cacheKey] = res.cube;
                    } else if (res && res.error) {
                        fetchErrors.push(res.error);
                    }
                }));
            } catch (err) {
                fetchErrors.push(err.message);
                console.error('fetchCube error:', err);
            }

            // 4. Redraw scene once missing cubes are fetched (if user is still on this frame)
            if (this.currentFrameIdx === frameIdx) {
                this.renderMolecule(symbols, coords, charges, showCharges);
                this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
                if (fetchErrors.length > 0) {
                    this.setStatus(`Frame ${frameIdx} | ${fetchErrors.join(' · ')}`, true);
                } else {
                    const names = activeOrbitals.map(o => o.title || o.orbital_key || o.key).join(', ');
                    this.setStatus(`Frame ${frameIdx} | Ready (${names})`);
                }
            }
        } else {
            const names = activeOrbitals.map(o => o.title || o.orbital_key || o.key).join(', ');
            this.setStatus(`Frame ${frameIdx} | Ready (${names})`);
        }
    }

    clearCache() {
        this.cubeCache = {};
        this.vdoCache = {};
    }
}
