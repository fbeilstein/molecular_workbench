/**
 * Step 4: Reaction Orbital Studio facade.
 * Modular architecture coordinates Orbitals3DViewer, OrbitalsPool, OrbitalsTimeline, and OrbitalsPrecompute.
 * Dynamic QM, interactive scrubber slider, bounded scrollable media pool, zoom & scaling.
 */
import { Orbitals3DViewer } from './orbitals/orbitals_3d.js';
import { OrbitalsPool } from './orbitals/orbitals_pool.js';
import { OrbitalsTimeline } from './orbitals/orbitals_timeline.js';
import { OrbitalsPrecompute } from './orbitals/orbitals_precompute.js';

export const OrbitalsUI = {
    state: {
        stepId: 'step_4_orbitals', uiSchema: null, totalFrames: 0, currentFrame: 0,
        isPlaying: false, playTimer: null, loadTimer: null, hideAllOrbitals: false,
        showCharges: false, isovalue: 0.03, viewer3D: null, pool: null,
        timeline: null, precompute: null, frameDataCache: {}, trajectory: []
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId; this.state.uiSchema = uiSchema;
        this.state.totalFrames = uiSchema.total_frames || 0;
        this.state.isovalue = uiSchema.default_isovalue || 0.03;
        this.state.trajectory = uiSchema.trajectory || [];
        this.state.trajectoryCharges = uiSchema.trajectory_charges || [];
        this.state.currentFrame = 0;
        const maxF = Math.max(0, this.state.totalFrames - 1);

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 0.5rem; background: var(--bg-dark); color: var(--text-main);">
                <div style="display: flex; justify-content: space-between; align-items: center; padding: 0.4rem 0.75rem; background: rgba(0,0,0,0.3); border-radius: 6px; border: 1px solid rgba(255,255,255,0.08); flex-wrap: wrap; gap: 0.5rem;">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <span style="font-weight: 700; font-size: 0.95rem; color: #38bdf8;">Reaction Orbital Studio</span>
                        <span class="badge" style="font-size: 0.7rem; background: rgba(56,189,248,0.15); color: #38bdf8;">Dynamic QM</span>
                        <span id="orbitals-status-label" style="font-size: 0.74rem; color: #94a3b8;">Initializing...</span>
                    </div>
                    <div style="display: flex; align-items: center; gap: 0.4rem; flex: 1; max-width: 480px; margin: 0 0.5rem;">
                        <button id="btn-play-pause" class="btn btn-sm btn-primary" style="padding: 0.2rem 0.6rem;">▶ Play</button>
                        <button id="btn-prev-frame" class="btn btn-sm btn-secondary">⏮</button>
                        <input type="range" id="transport-frame-scrubber" min="0" max="${maxF}" value="0" style="flex: 1; cursor: pointer; accent-color: #38bdf8; height: 6px;">
                        <button id="btn-next-frame" class="btn btn-sm btn-secondary">⏭</button>
                        <span id="transport-frame-counter" style="font-size: 0.74rem; font-weight: 600; color: #f1f5f9; min-width: 85px; text-align: right;">Frame 0 / ${maxF}</span>
                    </div>
                    <div style="display: flex; align-items: center; gap: 0.4rem;">
                        <button id="btn-toggle-lecture" class="btn btn-sm btn-secondary" title="Mute/Unmute Orbitals (M)">👁 Toggle MOs</button>
                        <button id="btn-toggle-esp" class="btn btn-sm btn-secondary" title="Toggle ESP / Charge Spheres (E)">⚡ ESP</button>
                        <button id="btn-precompute-all" class="btn btn-sm btn-outline" style="color: #10b981; border-color: #059669;">⚡ Precompute</button>
                    </div>
                </div>
                <div style="display: flex; gap: 0.5rem; height: 480px; max-height: 480px;">
                    <div style="width: 340px; height: 100%; display: flex; flex-direction: column; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.5rem; overflow: hidden;">
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
                            <span style="font-size: 0.8rem; font-weight: 600; color: #f1f5f9;">Orbital Media Pool</span>
                            <div style="display: flex; gap: 0.2rem;" id="pool-category-tabs">
                                <button class="btn btn-xs btn-primary tab-cat" data-cat="all">All</button>
                                <button class="btn btn-xs btn-secondary tab-cat" data-cat="frontier">Frontier</button>
                                <button class="btn btn-xs btn-secondary tab-cat" data-cat="bonding">Bonds</button>
                                <button class="btn btn-xs btn-secondary tab-cat" data-cat="lone_pairs">LP</button>
                                <button class="btn btn-xs btn-secondary tab-cat" data-cat="antibonding">π* / Antibonds</button>
                                <button class="btn btn-xs btn-secondary tab-cat" data-cat="esp">ESP</button>
                            </div>
                        </div>
                        <input type="text" id="pool-search-input" placeholder="Search atom, type, or label..." style="width: 100%; background: rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.1); border-radius: 4px; padding: 0.3rem 0.5rem; font-size: 0.74rem; color: #f8fafc; margin-bottom: 0.4rem;">
                        <div id="orbitals-pool-container" style="flex: 1; min-height: 0; overflow-y: auto; padding-right: 4px;"></div>
                    </div>
                    <div style="flex: 1; height: 100%; position: relative; border-radius: 6px; overflow: hidden; border: 1px solid rgba(255,255,255,0.08);">
                        <div id="orbitals-3d-canvas" style="width: 100%; height: 100%;"></div>
                        <div style="position: absolute; top: 10px; right: 10px; display: flex; flex-direction: column; gap: 4px; z-index: 10;">
                            <button id="btn-zoom-in" class="btn btn-xs btn-secondary" title="Zoom In (+)" style="font-size: 0.85rem; padding: 2px 7px; background: rgba(15,23,42,0.85);">+</button>
                            <button id="btn-zoom-out" class="btn btn-xs btn-secondary" title="Zoom Out (-)" style="font-size: 0.85rem; padding: 2px 7px; background: rgba(15,23,42,0.85);">−</button>
                            <button id="btn-zoom-fit" class="btn btn-xs btn-secondary" title="Fit to View (0)" style="font-size: 0.72rem; padding: 2px 5px; background: rgba(15,23,42,0.85);">⟲</button>
                        </div>
                        <div id="orbitals-viewer-overlay" style="position: absolute; bottom: 8px; left: 8px; font-size: 0.72rem; color: #94a3b8; background: rgba(15,23,42,0.75); padding: 2px 6px; border-radius: 4px; pointer-events: none;">Scroll wheel to zoom | Drag to rotate | Right-drag to pan</div>
                    </div>
                </div>
                <div id="orbitals-timeline-container" style="background: rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; overflow: hidden;"></div>
            </div>
        `;

        this.state.viewer3D = new Orbitals3DViewer(container.querySelector('#orbitals-3d-canvas'), container.querySelector('#orbitals-status-label'));
        this.state.viewer3D.init();
        this.state.pool = new OrbitalsPool(container.querySelector('#orbitals-pool-container'), (orb) => this.onAuditionOrbital(orb), (orb) => this.onAddToTimeline(orb));

        this.state.timeline = new OrbitalsTimeline(container.querySelector('#orbitals-timeline-container'), {
            totalFrames: this.state.totalFrames, currentFrame: this.state.currentFrame, isovalue: this.state.isovalue,
            onSeek: (f) => this.seekFrame(f), onVisibilityChange: () => this.refreshViewer(), onTimelineChange: (tl) => this.saveTimeline(tl)
        });
        this.state.timeline.setTimelineData(uiSchema.timeline, this.state.totalFrames);
        this.state.precompute = new OrbitalsPrecompute(this.state.stepId, () => this.getPrecomputeTargets(),
            (p) => {
                const b = container.querySelector('#btn-precompute-all');
                if (b) b.textContent = `Precomputing ${p.percent}%`;
                if (p.cubes) Object.assign(this.state.viewer3D.cubeCache, p.cubes);
            },
            () => { const b = container.querySelector('#btn-precompute-all'); if (b) b.textContent = '✓ Done (100%)'; this.refreshViewer(); }
        );
        this.bindEvents(container);
        await this.loadFrameData(this.state.currentFrame);
        this.preloadAllCachedCubes();
    },

    getPrecomputeTargets() {
        const targets = (this.state.timeline?.getTimelineObject().tracks || []).filter(t => t.orbital_key && t.orbital_key !== 'track_esp_charge').map(t => ({ orbital_key: t.orbital_key, chem_key: t.chem_key }));
        for (const a of (this.state.pool?.getAuditionedOrbitals() || [])) {
            const ok = a.key || a.orbital_key;
            if (!targets.some(t => t.orbital_key === ok)) targets.push({ orbital_key: ok, chem_key: a.chem_key || ok });
        }
        return targets.length ? targets : [{ orbital_key: 'HOMO', chem_key: 'HOMO' }, { orbital_key: 'LUMO', chem_key: 'LUMO' }];
    },

    bindEvents(container) {
        container.querySelector('#btn-play-pause').addEventListener('click', () => this.togglePlay());
        container.querySelector('#btn-prev-frame').addEventListener('click', () => this.stepFrame(-1));
        container.querySelector('#btn-next-frame').addEventListener('click', () => this.stepFrame(1));
        const scrubber = container.querySelector('#transport-frame-scrubber');
        if (scrubber) scrubber.addEventListener('input', (e) => this.seekFrame(parseInt(e.target.value, 10), false));

        container.querySelector('#btn-toggle-lecture').addEventListener('click', () => {
            this.state.hideAllOrbitals = !this.state.hideAllOrbitals;
            const b = container.querySelector('#btn-toggle-lecture');
            if (b) b.textContent = this.state.hideAllOrbitals ? '🚫 MOs: OFF' : '👁 MOs: ON';
            this.refreshViewer();
        });
        container.querySelector('#btn-toggle-esp').addEventListener('click', () => {
            this.state.showCharges = !this.state.showCharges;
            const b = container.querySelector('#btn-toggle-esp');
            if (b) {
                b.textContent = this.state.showCharges ? '⚡ ESP: ON' : '⚡ ESP: OFF';
                b.className = this.state.showCharges ? 'btn btn-sm btn-primary' : 'btn btn-sm btn-secondary';
            }
            if (!this.state.frameDataCache[this.state.currentFrame] && !this.state.trajectoryCharges[this.state.currentFrame]) {
                this.loadFrameData(this.state.currentFrame);
            } else {
                this.refreshViewer();
            }
        });
        container.querySelector('#btn-precompute-all').addEventListener('click', () => this.state.precompute.startPrecompute(this.state.totalFrames));

        const v3 = this.state.viewer3D;
        container.querySelector('#btn-zoom-in')?.addEventListener('click', () => v3?.zoomIn());
        container.querySelector('#btn-zoom-out')?.addEventListener('click', () => v3?.zoomOut());
        container.querySelector('#btn-zoom-fit')?.addEventListener('click', () => v3?.resetCamera());
        window.addEventListener('keydown', (e) => {
            if (e.target?.tagName === 'INPUT') return;
            if (e.key === '+' || e.key === '=') v3?.zoomIn(); else if (e.key === '-' || e.key === '_') v3?.zoomOut(); else if (e.key === '0') v3?.resetCamera();
        });

        container.querySelector('#pool-search-input').addEventListener('input', (e) => this.state.pool.setSearch(e.target.value));
        container.querySelectorAll('.tab-cat').forEach(btn => btn.addEventListener('click', () => {
            container.querySelectorAll('.tab-cat').forEach(b => { b.classList.remove('btn-primary'); b.classList.add('btn-secondary'); });
            btn.classList.add('btn-primary'); btn.classList.remove('btn-secondary');
            this.state.pool.setCategory(btn.dataset.cat);
        }));
    },

    async loadFrameData(frameIdx) {
        if (!this.state.frameDataCache[frameIdx]) {
            try {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST', headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_orbitals', frame_index: frameIdx })
                });
                const data = await res.json();
                if (data && data.success) this.state.frameDataCache[frameIdx] = data;
            } catch (err) { console.error('Failed to load frame data:', err); }
        }
        const fd = this.state.frameDataCache[frameIdx];
        if (fd && this.state.currentFrame === frameIdx) {
            this.state.pool.setOrbitalsFromFrameData(fd);
            if (this.state.timeline) this.state.timeline.updatePopulations(frameIdx, fd, this.state.uiSchema.symbols);
            this.refreshViewer();
        }
    },

    async fetchCube(frameIdx, orbKey, chemKey, isovalue) {
        const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'get_cube', frame_index: frameIdx, orbital: orbKey, chem_key: chemKey, isovalue: isovalue })
        });
        return await res.json();
    },

    _getActiveOrbitals(f) {
        const muted = new Set(this.state.timeline ? this.state.timeline.getMutedOrbitalKeys() : []);
        const activeTracks = (this.state.timeline ? this.state.timeline.getActiveOrbitalsAtFrame(f) : []).filter(t => !muted.has(t.orbital_key));
        const auditioned = (this.state.pool ? this.state.pool.getAuditionedOrbitals() : []).filter(a => !muted.has(a.key));
        const merged = [...activeTracks];
        for (const a of auditioned) {
            const aKey = a.orbital_key || a.key;
            if (!merged.some(t => (t.orbital_key || t.key) === aKey)) {
                merged.push({ id: 'audition_' + aKey, orbital_key: aKey, chem_key: a.chem_key || aKey, title: a.label || aKey, color: a.color || '#06b6d4' });
            }
        }
        return merged;
    },

    async refreshViewer() {
        const f = this.state.currentFrame, fd = this.state.frameDataCache[f];
        const coords = (fd && fd.coords) || (this.state.trajectory && this.state.trajectory[f]);
        if (!coords) return;
        const merged = this._getActiveOrbitals(f);
        const espActive = merged.some(t => t.orbital_key === 'track_esp_charge');
        const moTracks = merged.filter(t => t.orbital_key !== 'track_esp_charge');

        const charges = (fd && fd.charges) || this.state.trajectoryCharges[f] || null;
        await this.state.viewer3D.renderScene({
            symbols: this.state.uiSchema.symbols, coords, charges,
            showCharges: this.state.showCharges || espActive, activeOrbitals: moTracks,
            isovalue: this.state.isovalue, frameIdx: f, hideAllOrbitals: this.state.hideAllOrbitals,
            fetchCubeFn: (fi, ok, ck, iso) => this.fetchCube(fi, ok, ck, iso)
        });
    },

    seekFrame(frameIdx, isPlayback = false) {
        this.state.currentFrame = frameIdx;
        if (this.state.timeline) this.state.timeline.setCurrentFrame(frameIdx);
        const scrubber = document.querySelector('#transport-frame-scrubber');
        if (scrubber) scrubber.value = frameIdx;
        const counter = document.querySelector('#transport-frame-counter');
        if (counter) counter.textContent = `Frame ${frameIdx} / ${Math.max(0, this.state.totalFrames - 1)}`;

        const coords = (this.state.trajectory && this.state.trajectory[frameIdx]) || (this.state.frameDataCache[frameIdx] && this.state.frameDataCache[frameIdx].coords);
        if (coords) {
            const merged = this._getActiveOrbitals(frameIdx);
            const fd = this.state.frameDataCache[frameIdx];
            if (this.state.timeline && fd) this.state.timeline.updatePopulations(frameIdx, fd, this.state.uiSchema.symbols);
            const charges = (fd && fd.charges) || this.state.trajectoryCharges[frameIdx] || null;
            const espActive = merged.some(t => t.orbital_key === 'track_esp_charge');
            const moTracks = this.state.hideAllOrbitals ? [] : merged.filter(t => t.orbital_key !== 'track_esp_charge');
            this.state.viewer3D.renderFastFrame(
                this.state.uiSchema.symbols, coords, moTracks, this.state.isovalue,
                frameIdx, charges, this.state.showCharges || espActive
            );
            for (const orb of moTracks) {
                const ok = orb.orbital_key || orb.key || orb, ck = orb.chem_key || ok;
                if (!this.state.viewer3D.cubeCache[`${frameIdx}_${ok}`]) {
                    this.fetchCube(frameIdx, ok, ck, this.state.isovalue).then(res => {
                        if (res && res.success && res.cube) {
                            this.state.viewer3D.cubeCache[`${frameIdx}_${ok}`] = res.cube;
                            if (this.state.currentFrame === frameIdx) {
                                this.state.viewer3D.renderFastFrame(
                                    this.state.uiSchema.symbols, coords, moTracks, this.state.isovalue,
                                    frameIdx, charges, this.state.showCharges || espActive
                                );
                            }
                        }
                    }).catch(() => {});
                }
            }
        }

        if (!isPlayback) {
            clearTimeout(this.state.loadTimer);
            this.state.loadTimer = setTimeout(() => this.loadFrameData(frameIdx), 120);
        }
    },

    stepFrame(delta) {
        this.seekFrame(Math.max(0, Math.min(this.state.totalFrames - 1, this.state.currentFrame + delta)), this.state.isPlaying);
    },

    togglePlay() {
        this.state.isPlaying = !this.state.isPlaying;
        const btn = document.querySelector('#btn-play-pause');
        if (btn) btn.textContent = this.state.isPlaying ? '⏸ Pause' : '▶ Play';
        if (this.state.isPlaying) {
            this.preloadAllCachedCubes();
            if (this.state.currentFrame >= this.state.totalFrames - 1) this.seekFrame(0, true);
            this.state.playTimer = setInterval(() => {
                if (this.state.currentFrame >= this.state.totalFrames - 1) this.togglePlay();
                else this.stepFrame(1);
            }, 100);
        } else {
            clearInterval(this.state.playTimer);
            this.state.playTimer = null;
            this.loadFrameData(this.state.currentFrame);
        }
    },

    preloadAllCachedCubes() {
        if (this.state.precompute && this.state.viewer3D) {
            this.state.precompute.preloadCached(this.state.viewer3D.cubeCache, this.state.totalFrames, (f) => { if (f.includes(this.state.currentFrame)) this.refreshViewer(); });
        }
    },
    onAuditionOrbital() { this.refreshViewer(); this.preloadAllCachedCubes(); },
    onAddToTimeline(orbital) {
        if (this.state.timeline) {
            this.state.timeline.addTrack(orbital);
            if (this.state.pool && orbital) { this.state.pool.auditionKeys.delete(orbital.key); this.state.pool.render(); }
            this.refreshViewer(); this.preloadAllCachedCubes();
        }
    },
    async saveTimeline(timeline) {
        try {
            await fetch(`/api/step/${this.state.stepId}/custom`, { method: 'POST', headers: { 'Content-Type': 'application/json' }, body: JSON.stringify({ action_type: 'save_timeline', timeline }) });
        } catch (e) {}
    },
    getData() {
        return {
            timeline: this.state.timeline ? this.state.timeline.getTimelineObject() : {}
        };
    }
};
