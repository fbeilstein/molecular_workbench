/**
 * Step 5: Master Multistep Reaction Studio & Viewer facade.
 * Modular architecture coordinates Bundle3DViewer, BundlePesChart, and bundle_export.
 * Completely dynamic: displays tracks configured in Step 4; full zoom and scaling.
 */
import { Bundle3DViewer } from './bundle/bundle_3d.js';
import { BundlePesChart } from './bundle/bundle_pes.js';
import { exportMultiModelXyz, exportBundleZip } from './bundle/bundle_export.js';

export const BundleViewerUI = {
    state: {
        stepId: 'step_5_bundle_viewer',
        megajson: null,
        currentFrame: 0,
        isPlaying: false,
        playTimer: null,
        viewer3D: null,
        pesChart: null,
        activeOrbitals: new Set(),
        showCharges: false,
        hideAllOrbitals: false,
        isovalue: 0.03,
        tracks: [],
        frameDataCache: {}
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

        const initialActive = new Set();
        this.state.tracks.forEach(trk => {
            if (trk.default_on && trk.orbital_key) initialActive.add(trk.orbital_key);
        });
        this.state.activeOrbitals = initialActive;

        const ep = (megajson.energy_profile && megajson.energy_profile.frame_energies_kcal) ||
                   megajson.frames.map(f => f.rel_energy_kcal || 0.0);

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 0.75rem; color: var(--text-main);">
                <!-- Header with Export Actions -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.75rem; background: rgba(0,0,0,0.25); padding: 0.6rem 0.8rem; border-radius: 6px; border: 1px solid rgba(255,255,255,0.08);">
                    <div>
                        <h2 style="font-size: 1.1rem; font-weight: 700; color: #f8fafc; margin: 0;">Master Reaction Studio Viewer</h2>
                        <div style="font-size: 0.78rem; color: #94a3b8; display: flex; gap: 0.5rem; margin-top: 0.25rem; flex-wrap: wrap;">
                            <span class="chip blue">${megajson.symbols.length} atoms</span>
                            <span class="chip">${totalFrames} frames</span>
                            <span class="chip amber">Barrier: ${(megajson.overall_barrier_kcal || 0).toFixed(1)} kcal/mol</span>
                            <span class="chip emerald">Net ΔG: ${(megajson.overall_reaction_energy_kcal || 0).toFixed(1)} kcal/mol</span>
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.5rem; flex-wrap: wrap;">
                        <button id="btn-export-bundle-zip" class="btn btn-primary btn-sm" style="background: linear-gradient(135deg, #10b981, #059669); border-color: #059669;">
                            Export Bundle (.zip)
                        </button>
                        <button id="btn-export-xyz" class="btn btn-secondary btn-sm">Multi-Model XYZ</button>
                    </div>
                </div>

                <!-- Main Work Area: Left Checkboxes, Center 3D Viewer -->
                <div style="display: flex; gap: 0.6rem; height: 480px; max-height: 480px;">
                    <!-- Orbital & Overlay Selection Sidebar -->
                    <div style="width: 280px; height: 100%; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.6rem; display: flex; flex-direction: column; gap: 0.6rem; overflow: hidden;">
                        <span style="font-size: 0.82rem; font-weight: 600; color: #f1f5f9;">Configured Reaction Tracks</span>
                        <div id="bundle-tracks-list" style="flex: 1; min-height: 0; overflow-y: auto; display: flex; flex-direction: column; gap: 0.35rem; padding-right: 4px;"></div>

                        <div style="border-top: 1px solid rgba(255,255,255,0.08); padding-top: 0.5rem; display: flex; flex-direction: column; gap: 0.4rem;">
                            <label style="display: flex; align-items: center; gap: 0.5rem; font-size: 0.75rem; color: #cbd5e1; cursor: pointer;">
                                <input type="checkbox" id="chk-bundle-esp" ${this.state.showCharges ? 'checked' : ''}>
                                ⚡ Electrostatic Surface & Charges
                            </label>
                            <label style="display: flex; align-items: center; gap: 0.5rem; font-size: 0.75rem; color: #cbd5e1; cursor: pointer;">
                                <input type="checkbox" id="chk-bundle-mute" ${this.state.hideAllOrbitals ? 'checked' : ''}>
                                🚫 Mute All Orbitals (Lecture Mode)
                            </label>
                        </div>
                    </div>

                    <!-- 3D Viewer Container with Zoom & Scale Controls -->
                    <div style="flex: 1; height: 100%; position: relative; border-radius: 6px; overflow: hidden; border: 1px solid rgba(255,255,255,0.08);">
                        <div id="bundle-3d-canvas" style="width: 100%; height: 100%;"></div>
                        <div style="position: absolute; top: 10px; right: 10px; display: flex; flex-direction: column; gap: 4px; z-index: 10;">
                            <button id="btn-bundle-zoom-in" class="btn btn-xs btn-secondary" title="Zoom In (+)" style="font-size: 0.85rem; padding: 2px 7px; background: rgba(15,23,42,0.85);">+</button>
                            <button id="btn-bundle-zoom-out" class="btn btn-xs btn-secondary" title="Zoom Out (-)" style="font-size: 0.85rem; padding: 2px 7px; background: rgba(15,23,42,0.85);">−</button>
                            <button id="btn-bundle-zoom-fit" class="btn btn-xs btn-secondary" title="Fit to View (0)" style="font-size: 0.72rem; padding: 2px 5px; background: rgba(15,23,42,0.85);">⟲</button>
                        </div>
                        <div id="bundle-status-bar" style="position: absolute; bottom: 8px; left: 8px; font-size: 0.72rem; color: #94a3b8; background: rgba(15,23,42,0.75); padding: 2px 6px; border-radius: 4px; pointer-events: none;">
                            Ready
                        </div>
                    </div>
                </div>

                <!-- Interactive 2D Canvas PES Energy Profile Chart & Scrubber -->
                <div style="background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.5rem;">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.25rem; gap: 0.6rem;">
                        <span style="font-size: 0.8rem; font-weight: 600; color: #e2e8f0; white-space: nowrap;">Energy Profile & Scrubber</span>
                        <div style="display: flex; gap: 0.5rem; align-items: center; flex: 1; max-width: 500px;">
                            <button id="btn-bundle-play" class="btn btn-xs btn-primary">▶ Play</button>
                            <input type="range" id="bundle-frame-scrubber" min="0" max="${totalFrames - 1}" value="0" style="flex: 1; cursor: pointer; accent-color: #38bdf8; height: 6px;">
                            <span id="bundle-frame-readout" style="font-size: 0.72rem; color: #94a3b8; min-width: 85px; text-align: right;">Frame 0 / ${totalFrames - 1}</span>
                        </div>
                    </div>
                    <canvas id="bundle-pes-canvas" style="width: 100%; height: 140px; display: block;"></canvas>
                </div>
            </div>
        `;

        const canvas3D = container.querySelector('#bundle-3d-canvas');
        const statusBar = container.querySelector('#bundle-status-bar');
        const canvasPes = container.querySelector('#bundle-pes-canvas');

        this.state.viewer3D = new Bundle3DViewer(canvas3D, statusBar);
        this.state.viewer3D.init();

        this.state.pesChart = new BundlePesChart(canvasPes, (f) => this.seekFrame(f));
        this.state.pesChart.setData(ep, uiSchema.landmarks || []);

        this.renderTracksList(container);
        this.bindEvents(container);

        await this.seekFrame(0);
    },

    renderTracksList(container) {
        const listEl = container.querySelector('#bundle-tracks-list');
        if (!listEl) return;

        if (this.state.tracks.length === 0) {
            listEl.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.5rem;">No tracks configured. Return to Step 4 to audition and add tracks.</div>`;
            return;
        }

        let html = '';
        this.state.tracks.forEach(trk => {
            if (trk.orbital_key === 'track_esp_charge') return;
            const isChecked = this.state.activeOrbitals.has(trk.orbital_key);
            html += `
                <label style="display: flex; align-items: center; justify-content: space-between; padding: 0.35rem 0.5rem; background: rgba(255,255,255,0.03); border: 1px solid rgba(255,255,255,0.06); border-radius: 4px; font-size: 0.75rem; cursor: pointer;">
                    <div style="display: flex; align-items: center; gap: 0.4rem; overflow: hidden;">
                        <input type="checkbox" class="chk-track-orb" data-key="${trk.orbital_key}" ${isChecked ? 'checked' : ''}>
                        <span style="color: ${trk.color || '#38bdf8'}; font-weight: 600;">●</span>
                        <span style="white-space: nowrap; overflow: hidden; text-overflow: ellipsis; color: #e2e8f0;" title="${trk.title}">${trk.title}</span>
                    </div>
                </label>
            `;
        });
        listEl.innerHTML = html;

        listEl.querySelectorAll('.chk-track-orb').forEach(chk => {
            chk.addEventListener('change', () => {
                const k = chk.dataset.key;
                if (chk.checked) this.state.activeOrbitals.add(k);
                else this.state.activeOrbitals.delete(k);
                this.refreshViewer();
            });
        });
    },

    bindEvents(container) {
        container.querySelector('#chk-bundle-esp').addEventListener('change', (e) => {
            this.state.showCharges = e.target.checked;
            this.refreshViewer();
        });
        container.querySelector('#chk-bundle-mute').addEventListener('change', (e) => {
            this.state.hideAllOrbitals = e.target.checked;
            this.refreshViewer();
        });
        container.querySelector('#btn-bundle-play').addEventListener('click', () => this.togglePlay());
        const scrubber = container.querySelector('#bundle-frame-scrubber');
        if (scrubber) {
            scrubber.addEventListener('input', (e) => this.seekFrame(parseInt(e.target.value, 10)));
        }
        container.querySelector('#btn-export-xyz').addEventListener('click', () => exportMultiModelXyz(this.state.megajson));
        container.querySelector('#btn-export-bundle-zip').addEventListener('click', () => exportBundleZip(this.state.megajson));

        container.querySelector('#btn-bundle-zoom-in')?.addEventListener('click', () => this.state.viewer3D?.zoomIn());
        container.querySelector('#btn-bundle-zoom-out')?.addEventListener('click', () => this.state.viewer3D?.zoomOut());
        container.querySelector('#btn-bundle-zoom-fit')?.addEventListener('click', () => this.state.viewer3D?.resetCamera());
    },

    async seekFrame(f) {
        this.state.currentFrame = f;
        if (this.state.pesChart) this.state.pesChart.setCurrentFrame(f);
        const scrubber = document.querySelector('#bundle-frame-scrubber');
        if (scrubber) scrubber.value = f;
        const readout = document.querySelector('#bundle-frame-readout');
        if (readout && this.state.megajson) readout.textContent = `Frame ${f} / ${this.state.megajson.frames.length - 1}`;
        await this.refreshViewer();
    },

    async fetchCube(frameIdx, orbKey, chemKey, isovalue) {
        const res = await fetch(`/api/step/step_4_orbitals/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'get_cube', frame_index: frameIdx, orbital: orbKey, chem_key: chemKey, isovalue: isovalue })
        });
        return await res.json();
    },

    async refreshViewer() {
        const f = this.state.currentFrame;
        const frame = this.state.megajson.frames[f];
        if (!frame) return;

        if (!this.state.frameDataCache[f]) {
            try {
                const res = await fetch(`/api/step/step_4_orbitals/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_orbitals', frame_index: f })
                });
                const d = await res.json();
                if (d.success) this.state.frameDataCache[f] = d;
            } catch (err) {}
        }

        const fd = this.state.frameDataCache[f];
        const activeTracks = this.state.tracks.filter(t => this.state.activeOrbitals.has(t.orbital_key));

        await this.state.viewer3D.renderScene({
            symbols: this.state.megajson.symbols,
            coords: frame.coords,
            charges: fd ? fd.charges : null,
            showCharges: this.state.showCharges,
            activeOrbitals: activeTracks,
            isovalue: this.state.isovalue,
            frameIdx: f,
            hideAllOrbitals: this.state.hideAllOrbitals,
            fetchCubeFn: (fi, ok, ck, iso) => this.fetchCube(fi, ok, ck, iso)
        });
    },

    togglePlay() {
        this.state.isPlaying = !this.state.isPlaying;
        const btn = document.querySelector('#btn-bundle-play');
        if (btn) btn.textContent = this.state.isPlaying ? '⏸ Pause' : '▶ Play';

        if (this.state.isPlaying) {
            const total = this.state.megajson.frames.length;
            this.state.playTimer = setInterval(() => {
                if (this.state.currentFrame >= total - 1) this.seekFrame(0);
                else this.seekFrame(this.state.currentFrame + 1);
            }, 120);
        } else {
            clearInterval(this.state.playTimer);
        }
    }
};
