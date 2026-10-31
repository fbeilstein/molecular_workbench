/**
 * 3Dmol viewer and isosurface renderer for Step 7 Bundle Viewer.
 * Extends unified MoleculeViewer3D with mesh/cube isosurface rendering and ESP surfaces.
 */

import { MoleculeViewer3D } from '../common/molecule_viewer_3d.js';
import { getPhaseColors } from '../orbitals/orbitals_types.js';
import { prewarmEspCache as runEspPrewarm } from './bundle_esp_prewarmer.js';

export class Bundle3DViewer extends MoleculeViewer3D {
    constructor(containerEl, options = {}) {
        super(containerEl, { backgroundColor: '#090d16', ...options });
        this.meshCache = options.meshCache || {};
        this.espMeshCache = options.espMeshCache || {};
        this._espInFlight = new Set();
        this._isPrewarming = false;
        this.renderSeq = 0;
        this.activeFrameIdx = 0;
        this.init();
    }

    clearCache() {
        super.clearCache();
        this.meshCache = {};
        this.espMeshCache = {};
        this._espInFlight.clear();
        this._isPrewarming = false;
    }

    renderMolecule(symbols, coords, charges = null, showCharges = false, showEspSurface = false, frameIdx = 0, bonds = null) {
        if (!this.viewer || !symbols || !coords) return;
        this.activeFrameIdx = frameIdx;
        this.viewer.removeAllSurfaces();
        this.viewer.removeAllShapes();
        super.renderMolecule(symbols, coords, charges, showCharges, bonds);

        if (showEspSurface) {
            this._renderEspSurface(symbols, coords, charges, frameIdx);
        }
    }

    _renderEspSurface(symbols, coords, charges, frameIdx) {
        const n = Math.min(symbols.length, coords.length);
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

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false, showEspSurface = false, bonds = null) {
        if (!this.viewer || !symbols || !coords) return;
        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx, bonds);
        if (!this.hasInitializedCamera) { this.viewer.zoomTo(); this.hasInitializedCamera = true; }
        if (activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) this._drawOrbitalLayer(orb, frameIdx, isovalue);
        }
        this.viewer.render();
    }

    async renderScene(options) {
        const { symbols, coords, charges, showCharges, showEspSurface, activeOrbitals, isovalue, frameIdx, hideAllOrbitals, fetchCubeFn, resetCamera, bonds } = options;
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

        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx, bonds);
        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
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
