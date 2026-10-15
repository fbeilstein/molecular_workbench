/**
 * 3Dmol viewer and isosurface renderer for Step 5 Bundle Viewer.
 * Preserves camera zoom/orientation across frames; provides zoom & scale controls.
 */
import { getEspColor, getPhaseColors } from '../orbitals/orbitals_types.js';

export class Bundle3DViewer {
    constructor(containerEl, statusEl) {
        this.containerEl = containerEl;
        this.statusEl = statusEl;
        this.viewer = null;
        this.renderSeq = 0;
        this.cubeCache = {};
        this.hasInitializedCamera = false;
        this.moleculeScale = 1.0;
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
            id: 'bundle-3dmol-canvas'
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
            for (let i = 0; i < n; i++) {
                const q = charges[i];
                const [x, y, z] = coords[i];
                this.viewer.addSphere({
                    center: { x, y, z },
                    radius: 0.55 * this.moleculeScale,
                    color: getEspColor(q),
                    opacity: 0.55
                });
                this.viewer.addLabel(`${symbols[i]}: ${q > 0 ? '+' : ''}${q.toFixed(2)}`, {
                    position: { x, y: y + 0.35, z },
                    fontSize: 10,
                    fontColor: '#ffffff',
                    backgroundColor: 'rgba(15, 23, 42, 0.75)'
                });
            }
        }
    }

    async renderScene(options) {
        const {
            symbols, coords, charges, showCharges,
            activeOrbitals, isovalue, frameIdx,
            hideAllOrbitals, fetchCubeFn, resetCamera
        } = options;

        if (!this.viewer) return;
        const currentSeq = ++this.renderSeq;

        // 1. Draw atoms
        this.renderMolecule(symbols, coords, charges, showCharges);

        // 2. Only zoomTo on first render or explicit reset
        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        // 3. Check if all orbitals hidden
        if (hideAllOrbitals || !activeOrbitals || activeOrbitals.length === 0) {
            this.viewer.render();
            this.setStatus(`Frame ${frameIdx} | Ready`);
            return;
        }

        this.setStatus(`Frame ${frameIdx} | Loading Orbitals...`);

        // 4. Render active orbital isosurfaces
        for (const orb of activeOrbitals) {
            if (this.renderSeq !== currentSeq) return;

            const orbKey = orb.orbital_key || orb.key || orb;
            const chemKey = orb.chem_key || null;
            const color = orb.color || '#06b6d4';
            const cacheKey = `${frameIdx}_${orbKey}`;

            let cubeData = this.cubeCache[cacheKey];
            if (cubeData === null) continue;

            if (!cubeData && fetchCubeFn) {
                try {
                    const res = await fetchCubeFn(frameIdx, orbKey, chemKey, isovalue);
                    if (this.renderSeq !== currentSeq) return;
                    if (res && res.success && res.cube) {
                        cubeData = res.cube;
                        this.cubeCache[cacheKey] = cubeData;
                    } else {
                        this.cubeCache[cacheKey] = null;
                        continue;
                    }
                } catch (e) {
                    this.cubeCache[cacheKey] = null;
                    continue;
                }
            }

            if (cubeData && this.renderSeq === currentSeq) {
                try {
                    const vdo = new $3Dmol.VolumeData(cubeData, 'cube');
                    const pc = getPhaseColors(color);
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
                } catch (err) {
                    console.warn(`Error drawing isosurface for ${orbKey}:`, err);
                }
            }
        }

        if (this.renderSeq === currentSeq) {
            this.viewer.render();
            this.setStatus(`Frame ${frameIdx} | Ready`);
        }
    }
}
