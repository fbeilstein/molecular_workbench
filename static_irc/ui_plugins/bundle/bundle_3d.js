/**
 * 3Dmol viewer and isosurface renderer for Step 7 Bundle Viewer.
 * Extends unified MoleculeViewer3D with mesh/cube isosurface rendering and ESP surfaces.
 */

import { MoleculeViewer3D } from '../common/molecule_viewer_3d.js';
import { getPhaseColors } from '../orbitals/orbitals_types.js';
import { prewarmEspCache as runEspPrewarm } from './bundle_esp_prewarmer.js';

export class Bundle3DViewer extends MoleculeViewer3D {
    constructor(containerEl, options = {}) {
        super(containerEl, { backgroundColor: '#0f141c', ...options });
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

    renderMolecule(symbols, coords, charges = null, showCharges = false, showEspSurface = false, frameIdx = 0, bonds = null, showLabels = null, atomLabels = null) {
        if (!this.viewer || !symbols || !coords) return;
        this.activeFrameIdx = frameIdx;
        this.viewer.removeAllSurfaces();
        this.viewer.removeAllShapes();
        // Pass showCharges=false to super so super does not create a redundant VDW surface.
        // Charge labels and ESP surface are handled independently below.
        super.renderMolecule(symbols, coords, charges, false, bonds, showLabels, atomLabels);

        if (showCharges) {
            this._renderChargeLabels(symbols, coords, charges, atomLabels, showLabels);
        }

        if (showEspSurface) {
            this._renderEspSurface(symbols, coords, charges, frameIdx);
        }
    }

    _renderChargeLabels(symbols, coords, charges, atomLabels, showLabels) {
        if (!this.viewer || !symbols || !coords || !charges) return;
        const labelsVisible = (showLabels !== null && showLabels !== undefined) ? showLabels : this.showLabels;
        for (let i = 0; i < symbols.length && i < charges.length; i++) {
            const [x, y, z] = coords[i];
            const q = charges[i];
            if (q === undefined || q === null) continue;
            const sign = q > 0 ? '+' : '';
            const tag = (atomLabels && atomLabels[i]) ? atomLabels[i] : `${symbols[i]}${i}`;
            const lbl = labelsVisible
                ? `${tag}: ${sign}${Number(q).toFixed(2)}`
                : `${symbols[i]}: ${sign}${Number(q).toFixed(2)}`;
            this.viewer.addLabel(lbl, {
                position: { x, y: y + 0.35, z },
                fontSize: 10,
                fontColor: '#ffffff',
                backgroundColor: 'rgba(15, 23, 42, 0.85)',
                backgroundOpacity: 0.85,
                borderColor: '#c084fc',
                borderThickness: 1,
                inFront: true
            });
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
                    const surf = this.viewer.surfaces?.[sId];
                    const meshObj = surf?.[0] || surf;
                    const grp = meshObj?.geo?.geometryGroups?.[0];
                    if (grp && grp.vertexArray && grp.faceArray) {
                        const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
                        const fCount = grp.faceidx || grp.faceArray.length || 0;
                        const verts = [], faces = Array.from(grp.faceArray.slice(0, fCount)), cols = [];
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
        const pc = getPhaseColors(orb);
        const mKey = `${frameIdx}_${orbKey}`;
        const isDefaultIso = Math.abs(isovalue - 0.035) < 0.003;
        const isoKey = `${mKey}_${Math.round(isovalue * 1000)}`;

        const cubeData = this.cubeCache[mKey];
        if (cubeData && (!isDefaultIso || !this.meshCache[mKey])) {
            let vdo = this.vdoCache[mKey];
            if (!vdo) {
                try { vdo = new $3Dmol.VolumeData(cubeData, 'cube'); this.vdoCache[mKey] = vdo; } catch (e) {}
            }
            if (vdo) {
                try {
                    this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: pc.pos, opacity: 0.75 });
                    this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: pc.neg, opacity: 0.75 });
                    return;
                } catch (e) {}
            }
        }

        const mesh = this.meshCache[isoKey] || (isDefaultIso ? this.meshCache[mKey] : null);
        if (mesh) {
            const pV = mesh.pos?.vertices || mesh.pos?.vertexArr;
            const pF = mesh.pos?.faces || mesh.pos?.faceArr;
            const pN = mesh.pos?.normals || mesh.pos?.normalArr;
            if (pV && pF) this.viewer.addCustom({ vertexArr: pV, faceArr: pF, normalArr: pN || undefined, color: $3Dmol.CC.color(pc.pos), opacity: 0.75 });
            const nV = mesh.neg?.vertices || mesh.neg?.vertexArr;
            const nF = mesh.neg?.faces || mesh.neg?.faceArr;
            const nN = mesh.neg?.normals || mesh.neg?.normalArr;
            if (nV && nF) this.viewer.addCustom({ vertexArr: nV, faceArr: nF, normalArr: nN || undefined, color: $3Dmol.CC.color(pc.neg), opacity: 0.75 });
            return;
        }

        if (cubeData) {
            let vdo = this.vdoCache[mKey] || new $3Dmol.VolumeData(cubeData, 'cube');
            this.vdoCache[mKey] = vdo;
            try {
                this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: pc.pos, opacity: 0.75 });
                this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: pc.neg, opacity: 0.75 });
            } catch (e) {}
        }
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false, showEspSurface = false, bonds = null, showLabels = null, atomLabels = null) {
        if (!this.viewer || !symbols || !coords) return;
        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx, bonds, showLabels, atomLabels);
        if (!this.hasInitializedCamera) { this.viewer.zoomTo(); this.hasInitializedCamera = true; }
        if (activeOrbitals && activeOrbitals.length > 0) {
            for (const orb of activeOrbitals) this._drawOrbitalLayer(orb, frameIdx, isovalue);
        }
        this.viewer.render();
    }

    async renderScene(options) {
        const { symbols, coords, charges, showCharges, showEspSurface, showLabels, activeOrbitals, isovalue, frameIdx, hideAllOrbitals, fetchCubeFn, resetCamera, bonds, atomLabels } = options;
        if (!this.viewer) return;
        const currentSeq = ++this.renderSeq;

        const missing = [];
        if (!hideAllOrbitals && activeOrbitals && activeOrbitals.length > 0 && fetchCubeFn) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                if (orbKey === 'track_esp_charge') continue;
                const cacheKey = `${frameIdx}_${orbKey}`;
                const isDefaultIso = Math.abs(isovalue - 0.035) < 0.003;
                const isoKey = `${cacheKey}_${Math.round(isovalue * 1000)}`;
                const hasData = this.cubeCache[cacheKey] || this.meshCache[isoKey] || (isDefaultIso && this.meshCache[cacheKey]);
                if (!hasData) {
                    missing.push(fetchCubeFn(frameIdx, orbKey, orb.chem_key || null, isovalue).then(res => {
                        if (res?.success) {
                            if (res.cube) this.cubeCache[cacheKey] = res.cube;
                            if (res.pos || res.neg) {
                                this.meshCache[isoKey] = { pos: res.pos, neg: res.neg };
                                if (isDefaultIso) this.meshCache[cacheKey] = this.meshCache[isoKey];
                            }
                        } else {
                            this.meshCache[isoKey] = null;
                            if (isDefaultIso) this.meshCache[cacheKey] = null;
                        }
                    }).catch(() => {
                        this.meshCache[isoKey] = null;
                        if (isDefaultIso) this.meshCache[cacheKey] = null;
                    }));
                }
            }
        }

        if (missing.length > 0) {
            this.setStatus(`Frame ${frameIdx} | Loading Orbitals...`);
            await Promise.all(missing);
            if (this.renderSeq !== currentSeq) return;
        }

        this.renderMolecule(symbols, coords, charges, showCharges, showEspSurface, frameIdx, bonds, showLabels, atomLabels);
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
