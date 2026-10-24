/**
 * Step 5: Master Multistep Reaction Studio & Viewer facade.
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
        this.state.trajectoryCharges = uiSchema.trajectory_charges || (megajson && megajson.trajectory_charges) || [];

        const initialActive = new Set();
        this.state.tracks.forEach(trk => {
            if ((trk.default_on || trk.orbital_key === 'track_esp_charge') && trk.orbital_key) initialActive.add(trk.orbital_key);
        });
        this.state.activeOrbitals = initialActive;

        const ep = (megajson.energy_profile && megajson.energy_profile.frame_energies_kcal) ||
                   megajson.frames.map(f => f.rel_energy_kcal || 0.0);

        container.innerHTML = getBundleViewerTemplate(megajson, totalFrames, this.state);

        this.initControllers(container, ep);
        this.bindEvents(container);
        this.renderTracksList();

        const initialSvg = uiSchema.mechanism_svg || megajson.mechanism_svg || '';
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
        if (megajson.esp_surfaces) Object.assign(this.state.viewer3D.espMeshCache, megajson.esp_surfaces);
        if (this.state.showEspSurface) this.state.viewer3D.prewarmEspCache(megajson.symbols, megajson.frames, megajson.trajectory_charges);
        this.preloadAllCachedCubes();
    },

    initControllers(container, ep) {
        this.state.transport = new BundleTransportController(this.state, { onFrameSeek: async () => this.refreshViewer() });
        this.state.wiringMgr = new BundleWiringManager(this.state, {
            onSelectGroup: (grpId) => this.selectGroup(grpId),
            onSaveWiring: (grpId, cfg) => { this.state.svgManager.setGroupWiring(grpId, cfg); this.renderGroupsTree(); },
            onOpenMolecule: (grpId, cfg) => this.openMoleculeOrbitals(grpId, cfg)
        });
        this.state.viewer3D = new Bundle3DViewer(container.querySelector('#bundle-viewport-3d'), { isovalue: this.state.isovalue });
        this.state.pesChart = new BundlePesChart(container.querySelector('#bundle-pes-canvas'), {
            trajectory: this.state.megajson,
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

        const setMode = (isWiring) => {
            tabWiring.classList.toggle('active', isWiring);
            tabWiring.classList.toggle('btn-primary', isWiring);
            tabWiring.classList.toggle('btn-secondary', !isWiring);
            tabInteractive.classList.toggle('active', !isWiring);
            tabInteractive.classList.toggle('btn-primary', !isWiring);
            tabInteractive.classList.toggle('btn-secondary', isWiring);
            hint.textContent = isWiring ? 'Click element in tree/canvas to wire' : 'Click arrows to play • Right-click molecules';
            if (wiringStatus) {
                wiringStatus.textContent = isWiring ? 'Wiring Mode' : 'Hierarchy';
                wiringStatus.className = isWiring ? 'chip amber' : 'chip blue';
            }
            this.state.svgManager.setMode(isWiring ? 'edit' : 'interactive');
            if (isWiring) this.renderGroupsTree();
        };
        tabInteractive.addEventListener('click', () => setMode(false));
        tabWiring.addEventListener('click', () => setMode(true));

        const t = this.state.transport;
        container.querySelector('#btn-step-prev')?.addEventListener('click', () => t.stepBackward());
        container.querySelector('#btn-play-rev')?.addEventListener('click', () => t.playReverse());
        container.querySelector('#btn-playback-pause')?.addEventListener('click', () => t.pause());
        container.querySelector('#btn-play-fwd')?.addEventListener('click', () => t.playForward());
        container.querySelector('#btn-step-next')?.addEventListener('click', () => t.stepForward());

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

        let isoDebounceTimer = null;
        const handleIsovalue = (e, immediateAll = false) => {
            const newIso = parseFloat(e.target.value);
            if (isNaN(newIso) || newIso <= 0) return;
            this.state.isovalue = newIso;
            if (this.state.megajson) this.state.megajson.orbital_isovalue = newIso;
            const lbl = container.querySelector('#label-bundle-isovalue');
            if (lbl) lbl.textContent = newIso.toFixed(3);

            // 1. Immediately re-render active frame
            this.refreshViewer();

            // 2. Schedule recomputing all frames across entire trajectory
            clearTimeout(isoDebounceTimer);
            if (immediateAll) {
                this.preloadAllCachedCubes(true);
            } else {
                isoDebounceTimer = setTimeout(() => {
                    this.preloadAllCachedCubes(true);
                }, 250);
            }
        };

        const isoSlider = container.querySelector('#slider-bundle-isovalue');
        if (isoSlider) {
            isoSlider.addEventListener('input', (e) => handleIsovalue(e, false));
            isoSlider.addEventListener('change', (e) => handleIsovalue(e, true));
        }

        container.querySelector('#btn-recompute-highres').addEventListener('click', () => {
            this.state.highresManager.renderModal(container, this.state.tracks);
        });

        container.querySelector('#btn-export-bundle-zip').addEventListener('click', () => {
            exportMaxCompBundleZip(this.state.stepId, {
                isovalue: this.state.isovalue,
                esp_surfaces: this.state.viewer3D?.espMeshCache || {}
            });
        });

        container.querySelector('#btn-bundle-zoom-in')?.addEventListener('click', () => this.state.viewer3D?.zoomIn());
        container.querySelector('#btn-bundle-zoom-out')?.addEventListener('click', () => this.state.viewer3D?.zoomOut());
        container.querySelector('#btn-bundle-zoom-fit')?.addEventListener('click', () => this.state.viewer3D?.resetCamera());

        container.querySelector('#chk-bundle-charges')?.addEventListener('change', (e) => {
            this.state.showCharges = e.target.checked;
            this.refreshViewer();
        });
        container.querySelector('#chk-bundle-esp-surface')?.addEventListener('change', (e) => {
            this.state.showEspSurface = e.target.checked;
            if (e.target.checked) {
                this.state.activeOrbitals.add('track_esp_charge');
                this.state.viewer3D?.prewarmEspCache(this.state.megajson.symbols, this.state.megajson.frames, this.state.megajson.trajectory_charges);
            } else {
                this.state.activeOrbitals.delete('track_esp_charge');
            }
            this.renderTracksList();
            this.refreshViewer();
        });
        container.querySelector('#chk-bundle-mute').addEventListener('change', (e) => {
            this.state.hideAllOrbitals = e.target.checked;
            this.refreshViewer();
        });
    },

    async preloadAllCachedCubes(force = false) {
        if (!this.state.megajson || !this.state.viewer3D) return;
        const total = this.state.megajson.frames.length;
        const targets = this.state.tracks.filter(t => t.orbital_key && t.orbital_key !== 'track_esp_charge');
        if (!targets.length) return;

        const iso = this.state.isovalue;
        const isoSuffix = `_${Math.round(iso * 1000)}`;
        const seq = ++this._preloadSeq || (this._preloadSeq = 1);

        if (force) {
            this.state.viewer3D?.setStatus(`Remeshing all frames for isovalue ${iso.toFixed(3)}...`);
            fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ action_type: 'set_isovalue', isovalue: iso })
            }).catch(() => {});
        }

        const batchSize = 35;
        for (let start = 0; start < total; start += batchSize) {
            if (this._preloadSeq !== seq) return;
            const end = Math.min(total, start + batchSize);
            const frames = [];
            for (let i = start; i < end; i++) {
                const need = force || targets.some(t => {
                    const ok = t.orbital_key;
                    const isoKey = `${i}_${ok}${isoSuffix}`;
                    return !this.state.viewer3D.meshCache[isoKey] && !this.state.viewer3D.meshCache[`${i}_${ok}`];
                });
                if (need) frames.push(i);
            }
            if (frames.length > 0) {
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({
                            action_type: 'get_cached_meshes',
                            frames,
                            tracks: targets,
                            isovalue: iso
                        })
                    });
                    const d = await res.json();
                    if (this._preloadSeq !== seq) return;
                    if (d && d.success) {
                        if (d.meshes) {
                            for (const [k, v] of Object.entries(d.meshes)) {
                                this.state.viewer3D.meshCache[k] = v;
                                if (!k.endsWith(isoSuffix)) {
                                    this.state.viewer3D.meshCache[`${k}${isoSuffix}`] = v;
                                }
                            }
                        }
                        if (d.cubes) Object.assign(this.state.viewer3D.cubeCache, d.cubes);
                        if (frames.includes(this.state.currentFrame)) this.refreshViewer();
                    }
                } catch (e) {}
            }
        }

        if (this._preloadSeq === seq && force) {
            this.state.viewer3D?.setStatus(`✓ All ${total} frames remeshed (iso: ${iso.toFixed(3)})`);
            setTimeout(() => {
                if (this._preloadSeq === seq) this.state.viewer3D?.setStatus('Ready');
            }, 2000);
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

        const showEspSurface = this.state.showEspSurface || this.state.activeOrbitals.has('track_esp_charge');
        const showCharges = this.state.showCharges;
        const charges = (this.state.frameDataCache[f] && this.state.frameDataCache[f].charges) ||
                        (this.state.trajectoryCharges && this.state.trajectoryCharges[f]) ||
                        (this.state.megajson && this.state.megajson.trajectory_charges && this.state.megajson.trajectory_charges[f]) ||
                        (frame && frame.charges) || null;

        if (!charges && (showCharges || showEspSurface)) {
            if (!this._fetchingCharges) this._fetchingCharges = new Set();
            if (!this._fetchingCharges.has(f)) {
                this._fetchingCharges.add(f);
                fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_charges', frame_index: f })
                }).then(r => r.json()).then(data => {
                    this._fetchingCharges.delete(f);
                    if (data && data.charges) {
                        if (!this.state.frameDataCache[f]) this.state.frameDataCache[f] = {};
                        this.state.frameDataCache[f].charges = data.charges;
                        if (this.state.currentFrame === f) this.refreshViewer();
                    }
                }).catch(() => { this._fetchingCharges.delete(f); });
            }
        }
        if (this.state.svgManager) {
            this.state.svgManager.updatePushArrows(f, this.state.playDirection);
        }

        const activeTracks = this.state.tracks.filter(t => this.state.activeOrbitals.has(t.orbital_key));
        const moTracks = this.state.hideAllOrbitals ? [] : activeTracks.filter(t => t.orbital_key !== 'track_esp_charge');

        const doRender = () => {
            this.state.viewer3D.renderFastFrame(
                this.state.megajson.symbols, frame.coords, moTracks,
                this.state.isovalue, f, charges, showCharges, showEspSurface
            );
        };
        doRender();

        for (const orb of moTracks) {
            const ok = orb.orbital_key || orb.key || orb;
            const ck = orb.chem_key || ok;
            const isDefaultIso = Math.abs(this.state.isovalue - 0.03) < 0.003;
            const isoKey = `${f}_${ok}_${Math.round(this.state.isovalue * 1000)}`;
            const hasData = this.state.viewer3D.cubeCache[`${f}_${ok}`] ||
                this.state.viewer3D.meshCache[isoKey] ||
                (isDefaultIso && this.state.viewer3D.meshCache[`${f}_${ok}`]);

            if (!hasData) {
                this.fetchCube(f, ok, ck, this.state.isovalue).then(res => {
                    if (res && res.success) {
                        if (res.cube) {
                            this.state.viewer3D.cubeCache[`${f}_${ok}`] = res.cube;
                        }
                        if (res.pos || res.neg) {
                            const meshData = { pos: res.pos, neg: res.neg };
                            this.state.viewer3D.meshCache[`${f}_${ok}`] = meshData;
                            this.state.viewer3D.meshCache[isoKey] = meshData;
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
            const vk = norm(v.molecule_key || v.molecule_id || '');
            const vm = norm(v.molecule_name || v.name || v.title || '');
            return vk === target || vm === target || vk === idTarget || vm === idTarget || (v.directory && norm(v.directory).includes(target));
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
        if (!this.state.megajson.molecule_bundles) this.state.megajson.molecule_bundles = {};
        const b = this.state.megajson.molecule_bundles;
        b[molId] = bundleData;
        if (molName) b[molName] = bundleData;
        if (bundleData.molecule_key) b[bundleData.molecule_key] = bundleData;
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
