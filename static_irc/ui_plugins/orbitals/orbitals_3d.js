/**
 * 3D Isosurface and Molecular Geometry Viewer for Step 6 Orbitals.
 * Renders stick/sphere models and 3D Gaussian cube volumetric isosurfaces
 * using 3Dmol.js with distinct colors for bonding vs anti-bonding orbitals.
 */
import { getPhaseColors } from './orbitals_types.js';

export class Orbitals3DViewer {
    constructor(containerEl, statusEl = null) {
        this.containerEl = containerEl;
        this.statusEl = statusEl;
        this.viewer = null;
        this.cubeCache = {};
        this.vdoCache = {};
        this.moleculeScale = 1.0;
        this.currentOrientation = null;
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            console.error('3Dmol library is not loaded');
            this.setStatus('Error: 3Dmol library missing', true);
            return;
        }

        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: 'rgba(15, 23, 42, 1.0)',
            antialias: true
        });
        this.viewer.render();
    }

    saveCamera() {
        if (this.viewer) {
            try {
                this.currentOrientation = this.viewer.getView();
            } catch (e) {}
        }
    }

    restoreCamera() {
        if (this.viewer && this.currentOrientation) {
            try {
                this.viewer.setView(this.currentOrientation);
            } catch (e) {}
        }
    }

    zoom(delta) {
        if (!this.viewer) return;
        this.viewer.zoom(delta);
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

    renderFrameDirect(options) {
        const { symbols, coords, charges, showCharges, activeOrbitals, isovalue, frameIdx } = options;
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
                    fontSize: 10, fontColor: '#ffffff', backgroundColor: 'rgba(15, 23, 42, 0.75)', backgroundOpacity: 0.75
                });
            }
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
                            // DISTINCT COLOR: pass the full orbital object to getPhaseColors
                            const pc = getPhaseColors(orb);
                            this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: pc.pos, opacity: 0.72 });
                            this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: pc.neg, opacity: 0.72 });
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
        this.saveCamera();

        if (hideAllOrbitals || !activeOrbitals || activeOrbitals.length === 0) {
            this.renderMolecule(symbols, coords, charges, showCharges);
            if (resetCamera) this.resetCamera();
            else this.restoreCamera();
            return;
        }

        const missing = [];
        for (const orb of activeOrbitals) {
            const orbKey = orb.orbital_key || orb.key || orb;
            const k = `${frameIdx}_${orbKey}`;
            if (!this.cubeCache[k]) {
                missing.push(orb);
            }
        }

        if (missing.length > 0 && typeof fetchCubeFn === 'function') {
            this.setStatus(`Calculating orbital surface${missing.length > 1 ? 's' : ''}...`);
            for (const orb of missing) {
                const orbKey = orb.orbital_key || orb.key || orb;
                const ck = orb.chem_key || null;
                try {
                    const res = await fetchCubeFn(frameIdx, orbKey, isovalue, ck);
                    if (res && res.cube) {
                        this.cubeCache[`${frameIdx}_${orbKey}`] = res.cube;
                    }
                } catch (err) {
                    console.error('Error fetching cube for', orbKey, err);
                }
            }
            this.setStatus('Ready');
        }

        this.renderFrameDirect({
            symbols, coords, charges, showCharges,
            activeOrbitals, isovalue, frameIdx
        });

        if (resetCamera) {
            this.resetCamera();
        } else {
            this.restoreCamera();
        }
    }
}
