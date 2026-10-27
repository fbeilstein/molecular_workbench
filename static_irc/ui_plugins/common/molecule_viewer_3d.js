/**
 * Unified 3D Molecular & Orbital Viewer for IRC Pipeline.
 * Consolidates all 3Dmol rendering into a single high-performance component.
 * Geometry and bond perception delegated to three_scene_core.js.
 */

import { symbolsCoordsToMolBlock, renderDualSuperposition } from './three_scene_core.js';

export class MoleculeViewer3D {
    constructor(containerEl, options = {}) {
        this.containerEl = containerEl;
        this.options = options;
        this.viewer = null;
        this.statusEl = options.statusEl || null;
        this.cubeCache = options.cubeCache || {};
        this.vdoCache = {};
        this.hasInitializedCamera = false;
        this.moleculeScale = options.scale || 1.0;
        this.currentFrameIdx = 0;
        this.bondTable = options.bonds || null;
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            this.containerEl.innerHTML = '<div style="padding: 1rem; color: #dc2626; font-size: 0.85rem;">3Dmol.js library not loaded.</div>';
            return;
        }

        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: this.options.backgroundColor || '#0f172a',
            id: 'unified-3dmol-canvas'
        });

        this.containerEl.addEventListener('wheel', (e) => {
            e.preventDefault();
            this.zoom(e.deltaY < 0 ? 1.12 : 0.89);
        }, { passive: false });

        requestAnimationFrame(() => {
            if (this.viewer) { this.viewer.resize(); this.viewer.render(); }
        });
        setTimeout(() => {
            if (this.viewer) { this.viewer.resize(); this.viewer.render(); }
        }, 150);

        this._resizeHandler = () => this.resize();
        window.addEventListener('resize', this._resizeHandler);
    }

    zoom(factor) { if (this.viewer) { this.viewer.zoom(factor); this.viewer.render(); } }
    zoomIn() { this.zoom(1.2); }
    zoomOut() { this.zoom(0.83); }
    resetCamera() { if (this.viewer) { this.viewer.zoomTo(); this.viewer.render(); } }

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
        this.statusEl.style.color = isError ? '#ef4444' : '#94a3b8';
    }

    setBondTable(bonds) { this.bondTable = bonds; }
    buildMolBlock(symbols, coords, bonds = null) {
        return symbolsCoordsToMolBlock(symbols, coords, bonds || this.bondTable);
    }

    renderMolecule(symbols, coords, charges = null, showCharges = false, bonds = null) {
        if (!this.viewer || !symbols || !coords) return;
        this.viewer.clear();
        this.viewer.addModel(this.buildMolBlock(symbols, coords, bonds), 'mol');
        this.viewer.setStyle({}, {
            stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
            sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
        });

        if (showCharges && charges && charges.length >= symbols.length) {
            const atoms = this.viewer.selectedAtoms({});
            for (let i = 0; i < atoms.length && i < charges.length; i++) {
                const q = charges[i] || 0.0;
                atoms[i].partialCharge = q;
                atoms[i].properties = Object.assign(atoms[i].properties || {}, { partialCharge: q });
            }
            this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.65,
                colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
            });
            for (let i = 0; i < symbols.length; i++) {
                const [x, y, z] = coords[i];
                this.viewer.addLabel(`${symbols[i]}: ${charges[i] > 0 ? '+' : ''}${charges[i].toFixed(2)}`, {
                    position: { x, y: y + 0.35, z }, fontSize: 10, fontColor: '#ffffff',
                    backgroundColor: 'rgba(15, 23, 42, 0.75)', backgroundOpacity: 0.75
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
                    } catch (e) { console.error('VolumeData error:', e); }
                }
                if (vdo) {
                    try {
                        const isAnti = (orb.orbital_class === 'antibonding') || (orb.type && orb.type.includes('star'));
                        this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: isAnti ? '#fb923c' : '#38bdf8', opacity: 0.72 });
                        this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: isAnti ? '#db2777' : '#1d4ed8', opacity: 0.72 });
                    } catch (e) { console.error('addIsosurface error:', e); }
                }
            }
        }
        this.viewer.render();
    }

    resize() {
        if (this.viewer) { this.viewer.resize(); this.viewer.render(); }
    }

    loadTrajectory(traj) {
        if (!this.viewer || !traj || !traj.frames || !traj.frames.length) return;
        this.currentTraj = traj;
        this.symbols = traj.symbols || [];
        this.bondTable = traj.bonds || null;
        const tsIdx = traj.ts_frame_index || 0;
        this.setFrame(tsIdx);
        this.viewer.zoomTo();
    }

    setFrame(idx) {
        if (!this.viewer || !this.currentTraj || !this.currentTraj.frames) return;
        const frames = this.currentTraj.frames;
        if (!frames.length) return;
        const clampedIdx = Math.max(0, Math.min(frames.length - 1, idx));
        this.currentFrameIdx = clampedIdx;

        const frame = frames[clampedIdx];
        const symbols = frame.symbols || this.symbols;
        const coords = frame.coords || frame.positions || [];
        const bonds = frame.bonds || this.bondTable || null;
        this.renderMolecule(symbols, coords, frame.charges || null, false, bonds);
    }

    renderOverlap(prevSymbols, prevCoords, nextSymbols, nextCoords, rmsd = 0, isAligned = true) {
        renderDualSuperposition(this.viewer, prevSymbols, prevCoords, nextSymbols, nextCoords, this.moleculeScale, isAligned);
        const modeLabel = isAligned ? 'Aligned' : 'Raw Unaligned';
        this.setStatus(`Boundary Overlap (${modeLabel}) | RMSD: ${rmsd.toFixed(4)} Å (Cyan: Step K Product, ${isAligned ? 'Magenta' : 'Yellow'}: Step K+1 Reactant)`);
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false, bonds = null) {
        this.currentFrameIdx = frameIdx;
        this.renderMolecule(symbols, coords, charges, showCharges, bonds);
        if (!this.hasInitializedCamera) { this.viewer.zoomTo(); this.hasInitializedCamera = true; }
        if (activeOrbitals && activeOrbitals.length > 0) this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
    }

    async renderScene(options) {
        const { symbols, coords, charges, showCharges, activeOrbitals, isovalue, frameIdx, hideAllOrbitals, fetchCubeFn, resetCamera, bonds } = options;
        if (!this.viewer) return;
        this.currentFrameIdx = frameIdx;
        this.renderMolecule(symbols, coords, charges, showCharges, bonds);

        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (hideAllOrbitals || !activeOrbitals || activeOrbitals.length === 0) {
            this.setStatus(`Frame ${frameIdx} | Ready`);
            return;
        }

        this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
        const missing = [];
        if (fetchCubeFn) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                if (!this.cubeCache[`${frameIdx}_${orbKey}`]) missing.push({ orb, orbKey, cacheKey: `${frameIdx}_${orbKey}` });
            }
        }

        if (missing.length > 0) {
            this.setStatus(`Frame ${frameIdx} | Computing orbital (${missing.map(m => m.orbKey).join(', ')})...`);
            const fetchErrors = [];
            try {
                await Promise.all(missing.map(async (m) => {
                    const res = await fetchCubeFn(frameIdx, m.orbKey, m.orb.chem_key || null, isovalue);
                    if (res && res.success && res.cube) this.cubeCache[m.cacheKey] = res.cube;
                    else if (res && res.error) fetchErrors.push(res.error);
                }));
            } catch (err) { fetchErrors.push(err.message); }

            if (this.currentFrameIdx === frameIdx) {
                this.renderMolecule(symbols, coords, charges, showCharges, bonds);
                this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
                this.setStatus(fetchErrors.length > 0 ? `Frame ${frameIdx} | ${fetchErrors.join(' · ')}` : `Frame ${frameIdx} | Ready`, fetchErrors.length > 0);
            }
        } else {
            this.setStatus(`Frame ${frameIdx} | Ready`);
        }

    }

    clearCache() {
        this.cubeCache = {};
        this.vdoCache = {};
    }
}
