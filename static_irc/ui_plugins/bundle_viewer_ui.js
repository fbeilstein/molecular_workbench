/**
 * Step 5 & Step 7: Master Multistep Reaction Studio & Viewer facade.
 * Caching and synchronization delegated to bundle_preload_sync.js.
 */

import { Bundle3DViewer } from './bundle/bundle_3d.js?v=2.1.0';
import { BundlePesChart } from './bundle/bundle_pes.js?v=2.0.6';
import { SvgMechanismManager } from './bundle/bundle_svg.js?v=2.0.4';
import { BundleHighresManager } from './bundle/bundle_highres.js?v=2.0.3';
import { BundleMoleculeModal } from './bundle/bundle_molecule_modal.js?v=2.1.2';
import { exportMaxCompBundleZip } from './bundle/bundle_export.js';
import { BundleTransportController } from './bundle/bundle_transport.js?v=2.1.0';
import { BundleWiringManager } from './bundle/bundle_wiring_ui.js?v=2.0.9';
import { renderOrbitalTracksSidebar } from './bundle/bundle_orbital_panel.js';
import { getBundleViewerTemplate } from './bundle/bundle_viewer_template.js?v=2.1.0';
import {
    preloadAllCachedCubes,
    fetchCube,
    saveSvgToBackend,
    saveWiringToBackend,
} from './bundle/bundle_preload_sync.js';
import { bindBundleEvents } from './bundle/bundle_events.js';

export const BundleViewerUI = {
    state: {
        stepId: 'step_7_bundle', megajson: null, currentFrame: 0,
        isPlaying: false, playDirection: 'forward', playbackSpeed: 1.0, playTimer: null,
        viewer3D: null, pesChart: null, svgManager: null, highresManager: null,
        moleculeModal: null, transport: null, wiringMgr: null, activeOrbitals: new Set(),
        showCharges: false, showEspSurface: false, hideAllOrbitals: false, isovalue: 0.03, orbitalResolution: 32,
        tracks: [], frameDataCache: {}, activeTab: 'interactive'
    },

    async render(container, uiSchema, globalState, stepId) {
        window.bundleViewerUI = this;
        this.state.stepId = stepId;
        const megajson = uiSchema.megajson;
        this.state.megajson = megajson;

        if (!megajson || !megajson.frames) {
            container.innerHTML = `<div class="empty-state" style="padding: 2rem; text-align: center;">No stitched reaction pathway found. Complete earlier steps first.</div>`;
            return;
        }

        const totalFrames = megajson.frames.length;
        const timeline = uiSchema.orbital_timeline || megajson.orbital_timeline || { tracks: [] };
        this.state.tracks = timeline.tracks || [];
        this.state.isovalue = timeline.default_isovalue || 0.03;
        this.state.showCharges = false;
        this.state.showEspSurface = this.state.tracks.some(t => t.orbital_key === 'track_esp_charge');
        this.state.orbitalResolution = uiSchema.orbital_resolution || megajson.orbital_resolution || 32;
        this.state.trajectorySymbols = uiSchema.trajectory_symbols || megajson.trajectory_symbols || megajson.frames.map(f => f.symbols || megajson.symbols);
        this.state.trajectoryBonds = uiSchema.trajectory_bonds || megajson.trajectory_bonds || megajson.frames.map(f => f.bonds || null);
        this.state.trajectoryCharges = uiSchema.trajectory_charges || megajson.trajectory_charges || [];

        const initialActive = new Set();
        this.state.tracks.forEach(trk => {
            const tid = trk.track_id || trk.orbital_key || trk.id;
            if (trk.default_on !== false && tid) initialActive.add(tid);
            if ((trk.default_on || trk.orbital_key === 'track_esp_charge') && trk.orbital_key) initialActive.add(trk.orbital_key);
        });
        this.state.activeOrbitals = initialActive;

        const ep = (megajson.energy_profile?.frame_energies_kcal) || megajson.frames.map(f => f.rel_energy_kcal || 0.0);
        container.innerHTML = getBundleViewerTemplate(megajson, totalFrames, this.state);

        this.initControllers(container, ep);
        this.bindEvents(container);
        this.renderTracksList();

        const initialSvg = uiSchema.mechanism_svg || megajson.mechanism_svg || '';
        if (initialSvg) {
            try { this.state.svgManager.loadSvg(initialSvg); this.renderGroupsTree(); }
            catch (svgErr) { console.warn('[BundleViewer] Svg parse error:', svgErr); }
        }

        try { await this.state.transport.seekFrame(0); } catch (e) {}
        if (megajson.esp_surfaces) Object.assign(this.state.viewer3D.espMeshCache, megajson.esp_surfaces);
        if (this.state.showEspSurface) this.state.viewer3D.prewarmEspCache(this.state.trajectorySymbols || megajson.symbols, megajson.frames, megajson.trajectory_charges);
        preloadAllCachedCubes(this.state, false, () => this.refreshViewer());
    },

    initControllers(container, ep) {
        this.state.transport = new BundleTransportController(this.state, { onFrameSeek: async () => this.refreshViewer() });
        this.state.wiringMgr = new BundleWiringManager(this.state, {
            onSelectGroup: (grpId) => { this.state.svgManager.selectGroup(grpId); this.renderGroupsTree(); },
            onSaveWiring: (grpId, cfg) => { this.state.svgManager.setGroupWiring(grpId, cfg); this.renderGroupsTree(); },
            onOpenMolecule: (grpId, cfg) => this.openMoleculeOrbitals(grpId, cfg)
        });
        this.state.viewer3D = new Bundle3DViewer(container.querySelector('#bundle-viewport-3d'), { isovalue: this.state.isovalue });
        this.state.pesChart = new BundlePesChart(container.querySelector('#bundle-pes-canvas'), {
            trajectory: this.state.megajson, energies: ep, boundaries: this.state.megajson.step_boundaries || [],
            onSeekFrame: (f) => this.state.transport.seekFrame(f)
        });
        this.state.svgManager = new SvgMechanismManager(container.querySelector('#bundle-svg-view'), {
            initialWiring: this.state.megajson.svg_wiring || {},
            onPlaySegment: (sf, ef, dir) => this.state.transport.playSegment(sf, ef, dir),
            onJumpFrame: (tf) => this.state.transport.seekFrame(tf),
            onOpenMolecule: (grpId, cfg) => this.openMoleculeOrbitals(grpId, cfg),
            onGroupSelected: (grpId, cfg) => this.state.wiringMgr.renderWiringPanel(document.querySelector('#group-wiring-form'), grpId, cfg),
            onSaveWiring: (wiring) => saveWiringToBackend(this.state, wiring)
        });
        this.state.highresManager = new BundleHighresManager({
            stepId: this.state.stepId, tracks: this.state.tracks, totalFrames: this.state.megajson.frames.length,
            currentResolution: this.state.orbitalResolution, onCompleted: (newRes) => this.onHighresCompleted(newRes)
        });
        this.state.moleculeModal = new BundleMoleculeModal({
            onBundleUploaded: (molId, bundle, molName) => this.onMoleculeBundleUploaded(molId, bundle, molName)
        });
    },

    renderTracksList() {
        const list = document.querySelector('#bundle-tracks-list');
        renderOrbitalTracksSidebar(list, this.state.tracks, this.state.activeOrbitals, {
            onChanged: () => {
                const espCb = document.querySelector('#chk-bundle-esp-surface');
                if (espCb) espCb.checked = this.state.activeOrbitals.has('track_esp_charge');
                this.state.showEspSurface = this.state.activeOrbitals.has('track_esp_charge');
                this.refreshViewer();
            }
        });
    },

    renderGroupsTree() {
        const treeContainer = document.querySelector('#group-tree-container');
        if (treeContainer) {
            this.state.wiringMgr.renderGroupsTree(
                treeContainer,
                this.state.svgManager.groupsTree,
                this.state.svgManager.selectedGroupId,
                this.state.svgManager.wiring
            );
        }
    },

    bindEvents(container) {
        bindBundleEvents(container, this);
        const isoSlider = container.querySelector('#slider-bundle-isovalue');
        if (isoSlider) {
            let isoTimer = null;
            const onIso = (e, force) => {
                const newIso = parseFloat(e.target.value);
                if (isNaN(newIso) || newIso <= 0) return;
                this.state.isovalue = newIso;
                if (this.state.megajson) this.state.megajson.orbital_isovalue = newIso;
                const lbl = container.querySelector('#label-bundle-isovalue');
                if (lbl) lbl.textContent = newIso.toFixed(3);
                this.refreshViewer();
                clearTimeout(isoTimer);
                if (force) preloadAllCachedCubes(this.state, true, () => this.refreshViewer());
                else isoTimer = setTimeout(() => preloadAllCachedCubes(this.state, true, () => this.refreshViewer()), 250);
            };
            isoSlider.addEventListener('input', (e) => onIso(e, false));
            isoSlider.addEventListener('change', (e) => onIso(e, true));
        }
    },

    refreshViewer() {
        const f = this.state.currentFrame;
        const frame = this.state.megajson.frames[f];
        if (!frame || !frame.coords) return;
        const symbols = frame.symbols || this.state.trajectorySymbols?.[f] || this.state.megajson.symbols;
        const bonds = frame.bonds || this.state.trajectoryBonds?.[f] || null;

        const showEspSurface = this.state.showEspSurface || this.state.activeOrbitals.has('track_esp_charge');
        const showCharges = this.state.showCharges;
        const charges = this.state.frameDataCache[f]?.charges || this.state.trajectoryCharges?.[f] || frame.charges || null;

        if (!charges && (showCharges || showEspSurface)) {
            if (!this._fetchingCharges) this._fetchingCharges = new Set();
            if (!this._fetchingCharges.has(f)) {
                this._fetchingCharges.add(f);
                fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST', headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_charges', frame_index: f })
                }).then(r => r.json()).then(data => {
                    this._fetchingCharges.delete(f);
                    if (data?.charges) {
                        this.state.frameDataCache[f] = Object.assign(this.state.frameDataCache[f] || {}, { charges: data.charges });
                        if (this.state.currentFrame === f) this.refreshViewer();
                    }
                }).catch(() => { this._fetchingCharges.delete(f); });
            }
        }
        if (this.state.svgManager) this.state.svgManager.updatePushArrows(f, this.state.playDirection);

        const activeTracks = this.state.tracks.filter(t => {
            const tid = t.track_id || t.orbital_key || t.id;
            return this.state.activeOrbitals.has(tid) || (t.orbital_key && this.state.activeOrbitals.has(t.orbital_key));
        });
        const moTracks = this.state.hideAllOrbitals ? [] : activeTracks.filter(t => t.orbital_key !== 'track_esp_charge');
        const doRender = () => this.state.viewer3D.renderFastFrame(symbols, frame.coords, moTracks, this.state.isovalue, f, charges, showCharges, showEspSurface, bonds);
        doRender();

        for (const orb of moTracks) {
            const ok = orb.orbital_key || orb.key || orb;
            const ck = orb.chem_key || ok;
            const isoKey = `${f}_${ok}_${Math.round(this.state.isovalue * 1000)}`;
            const hasData = this.state.viewer3D.cubeCache[`${f}_${ok}`] || this.state.viewer3D.meshCache[isoKey] || (Math.abs(this.state.isovalue - 0.03) < 0.003 && this.state.viewer3D.meshCache[`${f}_${ok}`]);

            if (!hasData) {
                fetchCube(this.state, f, ok, ck, this.state.isovalue).then(res => {
                    if (res?.success) {
                        if (res.cube) this.state.viewer3D.cubeCache[`${f}_${ok}`] = res.cube;
                        if (res.pos || res.neg) {
                            this.state.viewer3D.meshCache[`${f}_${ok}`] = this.state.viewer3D.meshCache[isoKey] = { pos: res.pos, neg: res.neg };
                        }
                        if (this.state.currentFrame === f) doRender();
                    }
                }).catch(() => {});
            }
        }
    },

    openMoleculeOrbitals(molId, cfg) {
        const molName = (cfg && (cfg.moleculeName || cfg.molName)) || molId;
        const b = this.state.megajson.molecule_bundles || {};
        const norm = s => String(s || '').toLowerCase().replace(/[\s\-_]+/g, '');
        const target = norm(molName), idTarget = norm(molId);
        const existingBundle = b[molId] || b[molName] || Object.values(b).find(v => {
            const vk = norm(v.molecule_key || v.molecule_id || ''), vm = norm(v.molecule_name || v.name || '');
            return vk === target || vm === target || vk === idTarget || vm === idTarget;
        }) || null;
        this.state.moleculeModal.show(document.body, molId, molName, existingBundle);
    },

    onHighresCompleted(newResolution) {
        this.state.orbitalResolution = newResolution;
        const badge = document.querySelector('#badge-orbital-res');
        if (badge) badge.textContent = `Grid: ${newResolution}³ (High-Res)`;
        if (this.state.viewer3D) this.state.viewer3D.clearCache();
        this.refreshViewer();
    },

    onMoleculeBundleUploaded(molId, bundleData, molName) {
        const b = this.state.megajson.molecule_bundles = this.state.megajson.molecule_bundles || {};
        b[molId] = bundleData;
        if (molName) b[molName] = bundleData;
        if (bundleData.molecule_key) b[bundleData.molecule_key] = bundleData;
    }
};
