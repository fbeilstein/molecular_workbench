/**
 * 3Dmol viewer and isosurface renderer for Step 5 & Step 7 Bundle Viewer.
 * Preserves camera zoom/orientation across frames; provides zoom & scale controls.
 * Bond perception and ESP prewarming delegated to three_scene_core.js and bundle_esp_prewarmer.js.
 */

import { getPhaseColors } from '../orbitals/orbitals_types.js';
import { symbolsCoordsToMolBlock } from '../common/three_scene_core.js';
import { prewarmEspCache as runEspPrewarm } from './bundle_esp_prewarmer.js';

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
        this.meshCache = {}; this.cubeCache = {}; this.vdoCache = {}; this.espMeshCache = {};
        this._espInFlight.clear(); this._isPrewarming = false;
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            this.containerEl.innerHTML = '<div class="empty-state">3Dmol.js not loaded.</div>';
            return;
        }
        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: '#090d16', defaultcolors: $3Dmol.rasmolElementColors
        });
        this.containerEl.addEventListener('wheel', (e) => {
            e.preventDefault();
            this.zoom(e.deltaY < 0 ? 1.12 : 0.89);
        }, { passive: false });
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

        const n = Math.min(symbols.length, coords.length);
        const sdf = symbolsCoordsToMolBlock(symbols, coords);
        this.viewer.addModel(sdf, 'sdf');
        this.viewer.setStyle({}, {
            stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
            sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
        });

        if (showEspSurface) {
            const actualCharges = (charges && charges.length >= n) ? charges : symbols.map(s => s === 'H' ? 0.08 : (s === 'O' ? -0.50 : 0.0));
            const cachedEsp = this.espMeshCache[frameIdx];
            if (cachedEsp && cachedEsp.vertexArr && cachedEsp.faceArr) {
                this.viewer.addCustom({ vertexArr: cachedEsp.vertexArr, faceArr: cachedEsp.faceArr, color: cachedEsp.color, opacity: 0.65 });
            } else if (!this._espInFlight.has(frameIdx)) {
                this._espInFlight.add(frameIdx);
                const atoms = this.viewer.selectedAtoms({});
                for (let i = 0; i < atoms.length && i < actualCharges.length; i++) {
                    const qVal = actualCharges[i] || 0.0;
                    atoms[i].partialCharge = qVal;
                    atoms[i].properties = Object.assign(atoms[i].properties || {}, { partialCharge: qVal });
                }
                const p = this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                    opacity: 0.65,
                    colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
                });
                if (p && typeof p.then === 'function') {
                    p.then(sId => {
                        this._espInFlight.delete(frameIdx);
                        const grp = this.viewer.surfaces?.[sId]?.geo?.geometryGroups?.[0];
                        if (grp && grp.vertexArray && grp.faceArray) {
                            const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
                            const verts = [], faces = Array.from(grp.faceArray.slice(0, grp.faces * 3)), cols = [];
                            for (let vi = 0; vi < vCount; vi++) verts.push({ x: grp.vertexArray[vi * 3], y: grp.vertexArray[vi * 3 + 1], z: grp.vertexArray[vi * 3 + 2] });
                            if (grp.colorArray && grp.colorArray.length >= vCount * 3) {
                                for (let ci = 0; ci < vCount; ci++) cols.push({ r: grp.colorArray[ci * 3], g: grp.colorArray[ci * 3 + 1], b: grp.colorArray[ci * 3 + 2] });
                            }
                            this.espMeshCache[frameIdx] = { vertexArr: verts, faceArr: faces, color: cols.length ? cols : undefined };
                        }
                        if (this.activeFrameIdx === frameIdx) this.viewer.render();
                        else try { this.viewer.removeSurface(sId); } catch (e) {}
                    }).catch(() => { this._espInFlight.delete(frameIdx); });
                } else {
                    this._espInFlight.delete(frameIdx);
                    if (this.activeFrameIdx === frameIdx) this.viewer.render();
                }
            }
        }

        if (showCharges) {
            const actualCharges = (charges && charges.length >= n) ? charges : symbols.map(s => s === 'H' ? 0.08 : (s === 'O' ? -0.50 : 0.0));
            for (let i = 0; i < n; i++) {
                const [x, y, z] = coords[i];
                this.viewer.addLabel(`${symbols[i]}: ${actualCharges[i] > 0 ? '+' : ''}${Number(actualCharges[i]).toFixed(2)}`, {
                    position: { x, y: y + 0.35, z }, fontSize: 10, fontColor: '#ffffff',
                    backgroundColor: 'rgba(15, 23, 42, 0.75)', backgroundOpacity: 0.75
                });
            }
        }
        this.viewer.render();
    }

    async prewarmEspCache(symbols, frames, trajectoryCharges, onProgress) {
        if (this._isPrewarming) return;
        this._isPrewarming = true;
        try {
            await runEspPrewarm(symbols, frames, trajectoryCharges, this.espMeshCache, onProgress);
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

        const cubeData = this.cubeCache[mKey];
        if (cubeData && (!isDefaultIso || !this.meshCache[mKey])) {
            let vdo = this.vdoCache[mKey];
            if (!vdo) {
                try { vdo = new $3Dmol.VolumeData(cubeData, 'cube'); this.vdoCache[mKey] = vdo; } catch (e) {}
            }
            if (vdo) {
                try {
                    this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: pc.pos, opacity: 0.72 });
                    this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: pc.neg, opacity: 0.72 });
                    return;
                } catch (e) {}
            }
        }

        const mesh = this.meshCache[isoKey] || this.meshCache[mKey];
        if (mesh) {
            if (mesh.pos?.vertices && mesh.pos?.faces) this.viewer.addCustom({ vertexArr: mesh.pos.vertices, faceArr: mesh.pos.faces, color: $3Dmol.CC.color(pc.pos), opacity: 0.72 });
            if (mesh.neg?.vertices && mesh.neg?.faces) this.viewer.addCustom({ vertexArr: mesh.neg.vertices, faceArr: mesh.neg.faces, color: $3Dmol.CC.color(pc.neg), opacity: 0.72 });
            return;
        }

        if (cubeData) {
            let vdo = this.vdoCache[mKey] || new $3Dmol.VolumeData(cubeData, 'cube');
            this.vdoCache[mKey] = vdo;
            try {
                this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: pc.pos, opacity: 0.72 });
                this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: pc.neg, opacity: 0.72 });
            } catch (e) {}
        }
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false, showEspSurface = false) {
        if (!this.viewer || !symbols || !coords) return;
        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx);
        if (!this.hasInitializedCamera) { this.viewer.zoomTo(); this.hasInitializedCamera = true; }
        if (activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) this._drawOrbitalLayer(orb, frameIdx, isovalue);
        }
        this.viewer.render();
    }

    async renderScene(options) {
        const { symbols, coords, charges, showCharges, showEspSurface, activeOrbitals, isovalue, frameIdx, hideAllOrbitals, fetchCubeFn, resetCamera } = options;
        if (!this.viewer) return;
        const currentSeq = ++this.renderSeq;

        const missing = [];
        if (!hideAllOrbitals && activeOrbitals && activeOrbitals.length > 0 && fetchCubeFn) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                if (orbKey === 'track_esp_charge') continue;
                const cacheKey = `${frameIdx}_${orbKey}`;
                const isoKey = `${cacheKey}_${Math.round(isovalue * 1000)}`;
                if (!this.cubeCache[cacheKey] && !this.meshCache[isoKey] && !this.meshCache[cacheKey]) {
                    missing.push(fetchCubeFn(frameIdx, orbKey, orb.chem_key || null, isovalue).then(res => {
                        if (res?.success) {
                            if (res.cube) this.cubeCache[cacheKey] = res.cube;
                            if (res.pos || res.neg) this.meshCache[cacheKey] = this.meshCache[isoKey] = { pos: res.pos, neg: res.neg };
                        } else { this.meshCache[isoKey] = this.meshCache[cacheKey] = null; }
                    }).catch(() => { this.meshCache[isoKey] = this.meshCache[cacheKey] = null; }));
                }
            }
        }

        if (missing.length > 0) {
            this.setStatus(`Frame ${frameIdx} | Loading Orbitals...`);
            await Promise.all(missing);
            if (this.renderSeq !== currentSeq) return;
        }

        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx);
        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.viewer.render();
            this.hasInitializedCamera = true;
        }

        if (!hideAllOrbitals && activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) this._drawOrbitalLayer(orb, frameIdx, isovalue);
        }

        if (this.renderSeq === currentSeq) {
            this.viewer.render();
            this.setStatus(`Frame ${frameIdx} | Ready`);
        }
    }
}
