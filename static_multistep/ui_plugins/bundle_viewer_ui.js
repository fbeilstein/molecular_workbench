/**
 * Step 5: Master Multistep Reaction Studio & Viewer facade.
 */

import { Bundle3DViewer } from './bundle/bundle_3d.js?v=2.0.5';
import { BundlePesChart } from './bundle/bundle_pes.js?v=2.0.3';
import { SvgMechanismManager } from './bundle/bundle_svg.js?v=2.0.4';
import { BundleHighresManager } from './bundle/bundle_highres.js?v=2.0.3';
import { BundleMoleculeModal } from './bundle/bundle_molecule_modal.js';
import { exportMultiModelXyz, exportMaxCompBundleZip, purgeTemporaryCubes } from './bundle/bundle_export.js';
import { BundleTransportController } from './bundle/bundle_transport.js?v=2.0.4';
import { BundleWiringManager } from './bundle/bundle_wiring_ui.js';
import { renderOrbitalTracksSidebar } from './bundle/bundle_orbital_panel.js';
import { getBundleViewerTemplate } from './bundle/bundle_viewer_template.js';

export const BundleViewerUI = {
    state: {
        stepId: 'step_5_bundle_viewer', megajson: null, currentFrame: 0,
        isPlaying: false, playDirection: 'forward', playbackSpeed: 1.0, playTimer: null,
        viewer3D: null, pesChart: null, svgManager: null, highresManager: null,
        moleculeModal: null, transport: null, wiringMgr: null, activeOrbitals: new Set(),
        showCharges: false, hideAllOrbitals: false, isovalue: 0.03, orbitalResolution: 32,
        tracks: [], frameDataCache: {}, activeTab: 'interactive'
    },

    async render(container, uiSchema, globalState, stepId) {
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
        this.state.showCharges = this.state.tracks.some(t => t.orbital_key === 'track_esp_charge');
        this.state.orbitalResolution = uiSchema.orbital_resolution || megajson.orbital_resolution || 32;

        const initialActive = new Set();
        this.state.tracks.forEach(trk => {
            if (trk.default_on && trk.orbital_key) initialActive.add(trk.orbital_key);
        });
        this.state.activeOrbitals = initialActive;

        const ep = (megajson.energy_profile && megajson.energy_profile.frame_energies_kcal) ||
                   megajson.frames.map(f => f.rel_energy_kcal || 0.0);

        container.innerHTML = getBundleViewerTemplate(megajson, totalFrames, this.state);

        this.initControllers(container, ep);
        this.bindEvents(container);
        this.renderTracksList();

        const initialSvg = megajson.mechanism_svg || '';
        if (initialSvg) {
            try {
                this.state.svgManager.loadSvg(initialSvg);
                this.renderGroupsTree();
            } catch (svgErr) {
                console.warn('[BundleViewer] Could not parse initial mechanism SVG:', svgErr);
            }
        }

        try {
            await this.state.transport.seekFrame(0);
        } catch (seekErr) {
            console.error('[BundleViewer] Failed to render initial frame 0:', seekErr);
        }
        this.preloadAllCachedCubes();
    },

    initControllers(container, ep) {
        this.state.transport = new BundleTransportController(this.state, {
            onFrameSeek: async () => this.refreshViewer()
        });

        this.state.wiringMgr = new BundleWiringManager(this.state, {
            onSelectGroup: (grpId) => this.selectGroup(grpId),
            onSaveWiring: (grpId, cfg) => {
                this.state.svgManager.setGroupWiring(grpId, cfg);
                this.renderGroupsTree();
            },
            onOpenMolecule: (grpId, cfg) => this.openMoleculeOrbitals(grpId, cfg)
        });

        this.state.viewer3D = new Bundle3DViewer(container.querySelector('#bundle-viewport-3d'), {
            isovalue: this.state.isovalue
        });

        this.state.pesChart = new BundlePesChart(container.querySelector('#bundle-pes-canvas'), {
            energies: ep,
            boundaries: this.state.megajson.step_boundaries || [],
            onSeekFrame: (f) => this.state.transport.seekFrame(f)
        });

        this.state.svgManager = new SvgMechanismManager(container.querySelector('#bundle-svg-view'), {
            initialWiring: this.state.megajson.svg_wiring || {},
            onPlaySegment: (sf, ef, dir) => this.state.transport.playSegment(sf, ef, dir),
            onJumpFrame: (tf) => this.state.transport.seekFrame(tf),
            onOpenMolecule: (grpId, cfg) => this.openMoleculeOrbitals(grpId, cfg),
            onGroupSelected: (grpId, cfg) => this.renderWiringPanel(grpId, cfg),
            onSaveWiring: (wiring) => this.saveWiringToBackend(wiring)
        });

        this.state.highresManager = new BundleHighresManager({
            stepId: this.state.stepId,
            tracks: this.state.tracks,
            totalFrames: this.state.megajson.frames.length,
            currentResolution: this.state.orbitalResolution,
            onCompleted: (newRes) => this.onHighresCompleted(newRes)
        });

        this.state.moleculeModal = new BundleMoleculeModal({
            onBundleUploaded: (molId, bundle) => this.onMoleculeBundleUploaded(molId, bundle)
        });
    },

    renderTracksList() {
        const list = document.querySelector('#bundle-tracks-list');
        renderOrbitalTracksSidebar(list, this.state.tracks, this.state.activeOrbitals, {
            onChanged: () => this.refreshViewer()
        });
    },

    renderGroupsTree() {
        const treeContainer = document.querySelector('#group-tree-container');
        this.state.wiringMgr.renderGroupsTree(
            treeContainer,
            this.state.svgManager.groupsTree,
            this.state.svgManager.selectedGroupId,
            this.state.svgManager.wiring
        );
    },

    renderWiringPanel(grpId, cfg) {
        const panel = document.querySelector('#group-wiring-form');
        this.state.wiringMgr.renderWiringPanel(panel, grpId, cfg);
    },

    selectGroup(grpId) {
        this.state.svgManager.selectGroup(grpId);
        this.renderGroupsTree();
    },

    bindEvents(container) {
        const svgInput = container.querySelector('#input-upload-svg');
        svgInput.addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (!file) return;
            const reader = new FileReader();
            reader.onload = async () => {
                try {
                    const svgText = reader.result;
                    const success = this.state.svgManager.loadSvg(svgText);
                    if (success) {
                        if (this.state.megajson) this.state.megajson.mechanism_svg = svgText;
                        this.renderGroupsTree();
                        await this.saveSvgToBackend(svgText);
                        const hint = container.querySelector('#label-wiring-hint');
                        if (hint) {
                            hint.textContent = `✓ Loaded ${file.name}`;
                            hint.style.color = '#10b981';
                        }
                    } else {
                        alert('Could not parse SVG document.');
                    }
                } catch (err) {
                    alert('Error reading SVG file: ' + err.message);
                } finally {
                    svgInput.value = '';
                }
            };
            reader.readAsText(file);
        });

        const tabInteractive = container.querySelector('#tab-mode-interactive');
        const tabWiring = container.querySelector('#tab-mode-wiring');
        const hint = container.querySelector('#label-wiring-hint');
        const wiringStatus = container.querySelector('#label-wiring-status');

        tabInteractive.addEventListener('click', () => {
            tabInteractive.classList.add('active', 'btn-primary');
            tabInteractive.classList.remove('btn-secondary');
            tabWiring.classList.remove('active', 'btn-primary');
            tabWiring.classList.add('btn-secondary');
            hint.textContent = 'Click arrows to play • Right-click molecules';
            if (wiringStatus) {
                wiringStatus.textContent = 'Hierarchy';
                wiringStatus.className = 'chip blue';
            }
            this.state.svgManager.setMode('interactive');
        });

        tabWiring.addEventListener('click', () => {
            tabWiring.classList.add('active', 'btn-primary');
            tabWiring.classList.remove('btn-secondary');
            tabInteractive.classList.remove('active', 'btn-primary');
            tabInteractive.classList.add('btn-secondary');
            hint.textContent = 'Click element in tree/canvas to wire actions';
            if (wiringStatus) {
                wiringStatus.textContent = 'Wiring Mode';
                wiringStatus.className = 'chip amber';
            }
            this.state.svgManager.setMode('edit');
            this.renderGroupsTree();
        });

        const t = this.state.transport;
        container.querySelector('#btn-play-fwd').addEventListener('click', () => t.playForward());
        container.querySelector('#btn-play-rev').addEventListener('click', () => t.playReverse());
        container.querySelector('#btn-playback-pause').addEventListener('click', () => t.pause());

        container.querySelector('#sel-playback-speed').addEventListener('change', (e) => {
            this.state.playbackSpeed = parseFloat(e.target.value);
            if (this.state.isPlaying) {
                t.pause();
                if (this.state.playDirection === 'forward') t.playForward();
                else t.playReverse();
            }
        });

        const scrubber = container.querySelector('#bundle-frame-scrubber');
        if (scrubber) {
            scrubber.addEventListener('input', (e) => {
                t.pause();
                const newF = parseInt(e.target.value, 10);
                if (newF > this.state.currentFrame) t.setDirection('forward');
                else if (newF < this.state.currentFrame) t.setDirection('reverse');
                t.seekFrame(newF);
            });
        }

        container.querySelector('#slider-bundle-isovalue').addEventListener('input', (e) => {
            this.state.isovalue = parseFloat(e.target.value);
            container.querySelector('#label-bundle-isovalue').textContent = this.state.isovalue;
            this.refreshViewer();
        });

        container.querySelector('#btn-recompute-highres').addEventListener('click', () => {
            this.state.highresManager.renderModal(container, this.state.tracks);
        });

        container.querySelector('#btn-export-bundle-zip').addEventListener('click', () => exportMaxCompBundleZip(this.state.stepId));
        container.querySelector('#btn-export-xyz').addEventListener('click', () => exportMultiModelXyz(this.state.megajson));
        container.querySelector('#btn-purge-cubes')?.addEventListener('click', () => purgeTemporaryCubes(this.state.stepId));

        container.querySelector('#btn-bundle-zoom-in')?.addEventListener('click', () => this.state.viewer3D?.zoomIn());
        container.querySelector('#btn-bundle-zoom-out')?.addEventListener('click', () => this.state.viewer3D?.zoomOut());
        container.querySelector('#btn-bundle-zoom-fit')?.addEventListener('click', () => this.state.viewer3D?.resetCamera());

        container.querySelector('#chk-bundle-esp').addEventListener('change', (e) => {
            this.state.showCharges = e.target.checked;
            this.refreshViewer();
        });
        container.querySelector('#chk-bundle-mute').addEventListener('change', (e) => {
            this.state.hideAllOrbitals = e.target.checked;
            this.refreshViewer();
        });
    },

    async preloadAllCachedCubes() {
        if (!this.state.megajson || !this.state.viewer3D) return;
        const total = this.state.megajson.frames.length;
        const targets = this.state.tracks.filter(t => t.orbital_key && t.orbital_key !== 'track_esp_charge');
        if (!targets.length) return;

        const batchSize = 35;
        for (let start = 0; start < total; start += batchSize) {
            const end = Math.min(total, start + batchSize);
            const frames = [];
            for (let i = start; i < end; i++) {
                const need = targets.some(t => !this.state.viewer3D.meshCache[`${i}_${t.orbital_key}`] && !this.state.viewer3D.cubeCache[`${i}_${t.orbital_key}`]);
                if (need) frames.push(i);
            }
            if (frames.length > 0) {
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ action_type: 'get_cached_meshes', frames, tracks: targets, isovalue: this.state.isovalue })
                    });
                    const d = await res.json();
                    if (d && d.success) {
                        if (d.meshes) Object.assign(this.state.viewer3D.meshCache, d.meshes);
                        if (d.cubes) Object.assign(this.state.viewer3D.cubeCache, d.cubes);
                        if (frames.includes(this.state.currentFrame)) this.refreshViewer();
                    }
                } catch (e) {}
            }
        }
    },

    async fetchCube(frameIdx, orbKey, chemKey, isovalue) {
        const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                action_type: 'get_orbital_mesh',
                frame_index: frameIdx,
                orbital: orbKey,
                chem_key: chemKey,
                isovalue: isovalue
            })
        });
        return await res.json();
    },

    refreshViewer() {
        const f = this.state.currentFrame;
        const frame = this.state.megajson.frames[f];
        if (!frame) return;

        const espActive = this.state.activeOrbitals.has('track_esp_charge');
        const showCharges = this.state.showCharges || espActive;
        const charges = (this.state.frameDataCache[f] && this.state.frameDataCache[f].charges) ||
                        (this.state.megajson.trajectory_charges && this.state.megajson.trajectory_charges[f]) || null;
        const activeTracks = this.state.tracks.filter(t => this.state.activeOrbitals.has(t.orbital_key));
        const moTracks = this.state.hideAllOrbitals ? [] : activeTracks.filter(t => t.orbital_key !== 'track_esp_charge');

        this.state.viewer3D.renderFastFrame(
            this.state.megajson.symbols,
            frame.coords,
            moTracks,
            this.state.isovalue,
            f,
            charges,
            showCharges
        );

        for (const orb of moTracks) {
            const ok = orb.orbital_key || orb.key || orb;
            const ck = orb.chem_key || ok;
            const hasData = this.state.viewer3D.meshCache[`${f}_${ok}`] || this.state.viewer3D.cubeCache[`${f}_${ok}`];
            if (!hasData) {
                this.fetchCube(f, ok, ck, this.state.isovalue).then(res => {
                    if (res && res.success) {
                        if (res.pos || res.neg) {
                            this.state.viewer3D.meshCache[`${f}_${ok}`] = { pos: res.pos, neg: res.neg };
                        } else if (res.cube) {
                            this.state.viewer3D.cubeCache[`${f}_${ok}`] = res.cube;
                        }
                        if (this.state.currentFrame === f) {
                            this.state.viewer3D.renderFastFrame(
                                this.state.megajson.symbols,
                                frame.coords,
                                moTracks,
                                this.state.isovalue,
                                f,
                                charges,
                                showCharges
                            );
                        }
                    }
                }).catch(() => {});
            }
        }
    },

    openMoleculeOrbitals(molId, cfg) {
        const molName = (cfg && cfg.moleculeName) || molId;
        const existingBundle = (this.state.megajson.molecule_bundles && this.state.megajson.molecule_bundles[molId]) || null;
        this.state.moleculeModal.show(document.body, molId, molName, existingBundle);
    },

    onHighresCompleted(newResolution) {
        this.state.orbitalResolution = newResolution;
        const badge = document.querySelector('#badge-orbital-res');
        if (badge) badge.textContent = `Grid: ${newResolution}³ (High-Res)`;
        if (this.state.viewer3D) this.state.viewer3D.clearCache();
        this.refreshViewer();
    },

    onMoleculeBundleUploaded(molId, bundleData) {
        if (!this.state.megajson.molecule_bundles) {
            this.state.megajson.molecule_bundles = {};
        }
        this.state.megajson.molecule_bundles[molId] = bundleData;
    },

    async saveSvgToBackend(svgText) {
        try {
            await fetch(`/api/step/${this.state.stepId}/custom`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ action_type: 'save_mechanism_svg', svg_text: svgText }) });
        } catch (err) {}
    },

    async saveWiringToBackend(wiring) {
        try {
            await fetch(`/api/step/${this.state.stepId}/custom`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ action_type: 'save_svg_wiring', wiring }) });
        } catch (err) {}
    }
};
