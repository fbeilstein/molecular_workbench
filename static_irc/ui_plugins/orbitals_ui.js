/**
 * Step 6: Reaction Orbital Studio facade for IRC Pipeline.
 * Modular architecture coordinates Orbitals3DViewer, OrbitalsPool, OrbitalsTimeline, and OrbitalsPrecompute.
 * Integrates decoupled invariant track_ids for seamless cross-step tracking and 3D continuity.
 */
import { Orbitals3DViewer } from './orbitals/orbitals_3d.js';
import { OrbitalsPool } from './orbitals/orbitals_pool.js';
import { OrbitalsTimeline } from './orbitals/orbitals_timeline.js';
import { OrbitalsPrecompute } from './orbitals/orbitals_precompute.js';
import { isAntiBondingOrbital } from './orbitals/orbitals_types.js';
import { buildOrbitalsLayout } from './orbitals/orbitals_layout.js';
import { fetchFrameData, fetchOrbitalCube, postTimeline } from './orbitals/orbitals_api.js';

export const OrbitalsUI = {
    state: {
        stepId: 'step_6_orbitals', uiSchema: null, totalFrames: 0, currentFrame: 0, isPlaying: false,
        playTimer: null, loadTimer: null, hideAllOrbitals: false, showCharges: false, showLabels: false, isovalue: 0.03,
        viewer3D: null, pool: null, timeline: null, precompute: null, frameDataCache: {}, trajectory: [],
        trajectorySymbols: [], trajectoryAtomLabels: [], trajectoryBonds: [], trajectoryCharges: [], cubePromises: new Map(),
        visibleTrackIds: new Set()
    },

    toggleVisibility(tid) {
        if (!tid) return;
        const normId = (tid === 'track_esp_charge') ? 'trk_esp' : tid;
        if (this.state.visibleTrackIds.has(normId)) {
            this.state.visibleTrackIds.delete(normId);
        } else {
            this.state.visibleTrackIds.add(normId);
        }
        this.state.timeline?.syncVisibilitySet(this.state.visibleTrackIds);
        this.state.pool?.syncVisibilitySet(this.state.visibleTrackIds);
        this.refreshViewer();
        this.preloadAllCachedCubes();
    },

    async render(container, uiSchema, globalState, stepId) {
        Object.assign(this.state, {
            stepId: stepId || 'step_6_orbitals', uiSchema, totalFrames: uiSchema.total_frames || 0,
            isovalue: uiSchema.default_isovalue || 0.03, trajectory: uiSchema.trajectory || [],
            trajectorySymbols: uiSchema.trajectory_symbols || [], trajectoryAtomLabels: uiSchema.trajectory_atom_labels || [], trajectoryBonds: uiSchema.trajectory_bonds || [],
            trajectoryCharges: uiSchema.trajectory_charges || [], currentFrame: 0,
            showLabels: localStorage.getItem('orbitals_show_labels') === '1'
        });

        const initialVisible = (uiSchema.timeline?.tracks || [])
            .filter(t => t.default_on !== false)
            .map(t => (t.orbital_key === 'track_esp_charge' || t.track_id === 'trk_esp') ? 'trk_esp' : (t.track_id || t.orbital_key || t.id));
        this.state.visibleTrackIds = new Set(initialVisible);

        const maxF = Math.max(0, this.state.totalFrames - 1);
        buildOrbitalsLayout(container, maxF);

        this.state.viewer3D = new Orbitals3DViewer(
            container.querySelector('#orbitals-3d-canvas'),
            container.querySelector('#orbitals-status-label'),
            { showLabels: this.state.showLabels }
        );
        this.state.viewer3D.init();
        this.state.pool = new OrbitalsPool(
            container.querySelector('#orbitals-pool-container'),
            tid => this.toggleVisibility(tid),
            orb => this.onAddToTimeline(orb)
        );

        this.state.timeline = new OrbitalsTimeline(container.querySelector('#orbitals-timeline-container'), {
            totalFrames: this.state.totalFrames, currentFrame: this.state.currentFrame, isovalue: this.state.isovalue,
            onSeek: f => this.seekFrame(f),
            onToggleVisibility: tid => this.toggleVisibility(tid),
            onRemoveTrack: tid => {
                this.state.visibleTrackIds.delete(tid);
                this.state.pool?.syncVisibilitySet(this.state.visibleTrackIds);
                this.refreshViewer();
            },
            onTimelineChange: tl => this.saveTimeline(tl)
        });
        this.state.timeline.setTimelineData(uiSchema.timeline, this.state.totalFrames);
        this.state.timeline.syncVisibilitySet(this.state.visibleTrackIds);
        this.state.pool.syncVisibilitySet(this.state.visibleTrackIds);
        this.state.precompute = new OrbitalsPrecompute(this.state.stepId, () => this.getPrecomputeTargets(),
            p => { const b = container.querySelector('#btn-precompute-all'); if (b) b.textContent = `Precomputing ${p.percent}%`; if (p.cubes) Object.assign(this.state.viewer3D.cubeCache, p.cubes); },
            () => { const b = container.querySelector('#btn-precompute-all'); if (b) b.textContent = '✓ Done (100%)'; this.refreshViewer(); }
        );
        this.bindEvents(container);
        await this.loadFrameData(this.state.currentFrame);
        this.preloadAllCachedCubes();
    },

    getPrecomputeTargets() {
        const targets = (this.state.timeline?.getTimelineObject().tracks || []).filter(t => t.orbital_key && t.orbital_key !== 'track_esp_charge').map(t => ({ track_id: t.track_id || t.orbital_key, orbital_key: t.orbital_key, chem_key: t.chem_key }));
        for (const a of (this.state.pool?.getAuditionedOrbitals() || [])) {
            const tid = a.track_id || a.key || a.orbital_key;
            if (!targets.some(t => t.track_id === tid || t.orbital_key === tid)) targets.push({ track_id: tid, orbital_key: tid, chem_key: a.chem_key || tid });
        }
        return targets.length ? targets : [{ track_id: 'trk_homo', orbital_key: 'HOMO', chem_key: 'HOMO' }, { track_id: 'trk_lumo', orbital_key: 'LUMO', chem_key: 'LUMO' }];
    },

    bindEvents(container) {
        container.querySelector('#btn-play-pause')?.addEventListener('click', () => this.togglePlay());
        container.querySelector('#btn-prev-frame')?.addEventListener('click', () => this.stepFrame(-1));
        container.querySelector('#btn-next-frame')?.addEventListener('click', () => this.stepFrame(1));
        const sc = container.querySelector('#transport-frame-scrubber');
        if (sc) ['input', 'change'].forEach(ev => sc.addEventListener(ev, e => this.seekFrame(parseInt(e.target.value, 10), false)));

        container.querySelector('#btn-toggle-lecture')?.addEventListener('click', () => {
            this.state.hideAllOrbitals = !this.state.hideAllOrbitals;
            const b = container.querySelector('#btn-toggle-lecture');
            if (b) b.textContent = this.state.hideAllOrbitals ? '🚫 MOs: OFF' : '👁 MOs: ON';
            this.refreshViewer();
        });
        container.querySelector('#btn-toggle-esp')?.addEventListener('click', () => {
            this.state.showCharges = !this.state.showCharges;
            const b = container.querySelector('#btn-toggle-esp');
            if (b) { b.textContent = this.state.showCharges ? '⚡ ESP: ON' : '⚡ ESP: OFF'; b.className = this.state.showCharges ? 'btn btn-sm btn-primary' : 'btn btn-sm btn-secondary'; }
            (!this.state.frameDataCache[this.state.currentFrame] && !this.state.trajectoryCharges[this.state.currentFrame]) ? this.loadFrameData(this.state.currentFrame) : this.refreshViewer();
        });

        const updateLabelsBtn = () => {
            const b = container.querySelector('#btn-toggle-labels');
            if (b) {
                b.textContent = this.state.showLabels ? '🏷 Labels: ON' : '🏷 Labels';
                b.className = this.state.showLabels ? 'btn btn-sm btn-primary' : 'btn btn-sm btn-secondary';
            }
            const bf = container.querySelector('#btn-toggle-labels-float');
            if (bf) {
                bf.style.background = this.state.showLabels ? '#0284c7' : 'rgba(15,23,42,0.85)';
                bf.style.borderColor = this.state.showLabels ? '#38bdf8' : 'rgba(255,255,255,0.2)';
                bf.title = this.state.showLabels ? 'Hide Atom Labels (L)' : 'Show Atom Labels (L)';
            }
        };
        updateLabelsBtn();

        const toggleLabels = () => {
            this.state.showLabels = !this.state.showLabels;
            if (this.state.viewer3D) this.state.viewer3D.showLabels = this.state.showLabels;
            localStorage.setItem('orbitals_show_labels', this.state.showLabels ? '1' : '0');
            updateLabelsBtn();
            this.refreshViewer();
        };
        container.querySelector('#btn-toggle-labels')?.addEventListener('click', toggleLabels);
        container.querySelector('#btn-toggle-labels-float')?.addEventListener('click', toggleLabels);

        container.querySelector('#btn-precompute-all')?.addEventListener('click', () => this.state.precompute.startPrecompute(this.state.totalFrames, false));
        container.querySelector('#btn-purge-cache')?.addEventListener('click', async () => {
            if (!confirm('Are you sure you want to purge all memory and disk orbital cube caches?')) return;
            const btn = container.querySelector('#btn-purge-cache');
            const orig = btn.textContent;
            btn.disabled = true;
            btn.textContent = 'Purging...';
            try {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'clear_cache' })
                });
                const data = await res.json();
                if (data.success) {
                    if (this.state.viewer3D?.cubeCache) {
                        for (const k in this.state.viewer3D.cubeCache) delete this.state.viewer3D.cubeCache[k];
                    }
                    this.state.frameDataCache = {};
                    const preBtn = container.querySelector('#btn-precompute-all');
                    if (preBtn) preBtn.textContent = '⚡ Precompute';
                    btn.textContent = '✓ Cleared';
                    setTimeout(() => { btn.textContent = orig; btn.disabled = false; }, 1500);
                    this.refreshViewer();
                } else {
                    alert('Failed to clear cache: ' + (data.error || 'Unknown error'));
                    btn.textContent = orig;
                    btn.disabled = false;
                }
            } catch (e) {
                console.error('Error clearing cache:', e);
                btn.textContent = orig;
                btn.disabled = false;
            }
        });

        const v3 = this.state.viewer3D;
        container.querySelector('#btn-zoom-in')?.addEventListener('click', () => v3?.zoomIn());
        container.querySelector('#btn-zoom-out')?.addEventListener('click', () => v3?.zoomOut());
        container.querySelector('#btn-zoom-fit')?.addEventListener('click', () => v3?.resetCamera());
        window.addEventListener('keydown', e => {
            if (e.target?.tagName === 'INPUT') return;
            if (e.key === '+' || e.key === '=') v3?.zoomIn();
            else if (e.key === '-' || e.key === '_') v3?.zoomOut();
            else if (e.key === '0') v3?.resetCamera();
            else if (e.key === 'l' || e.key === 'L') toggleLabels();
        });
        container.querySelector('#pool-search-input').addEventListener('input', e => this.state.pool.setSearch(e.target.value));
        container.querySelectorAll('.tab-cat').forEach(btn => btn.addEventListener('click', () => {
            container.querySelectorAll('.tab-cat').forEach(b => { b.classList.remove('btn-primary'); b.classList.add('btn-secondary'); });
            btn.classList.add('btn-primary'); btn.classList.remove('btn-secondary');
            this.state.pool.setCategory(btn.dataset.cat);
        }));
    },

    async loadFrameData(frameIdx) {
        if (!this.state.frameDataCache[frameIdx]) {
            const data = await fetchFrameData(this.state.stepId, frameIdx);
            if (data && data.success) this.state.frameDataCache[frameIdx] = data;
        }
        const fd = this.state.frameDataCache[frameIdx];
        if (fd && this.state.currentFrame === frameIdx) {
            this.state.pool.setOrbitalsFromFrameData(fd);
            const syms = fd.symbols || (this.state.trajectorySymbols && this.state.trajectorySymbols[frameIdx]) || this.state.uiSchema.symbols;
            if (this.state.timeline) this.state.timeline.updatePopulations(frameIdx, fd, syms);
            this.refreshViewer();
        }
    },

    async fetchCube(frameIdx, orbKey, chemKey, isovalue, trackId) {
        return await fetchOrbitalCube(this.state.stepId, this.state.cubePromises, this.state.viewer3D?.cubeCache, frameIdx, orbKey, chemKey, isovalue, trackId);
    },

    _getActiveOrbitals(f) {
        const fd = this.state.frameDataCache[f];
        const allPool = fd ? [
            ...Object.values(fd.frontier_orbitals || {}),
            ...(fd.all_ibos || []),
            ...(fd.antibonds || []),
            ...(fd.canonical_pi || [])
        ] : [];
        const poolByTid = new Map(allPool.map(item => [item.track_id, item]));
        const tlTracks = this.state.timeline ? this.state.timeline.tracks : [];
        const tlByTid = new Map(tlTracks.map(t => [t.track_id || t.orbital_key || t.id, t]));

        const active = [];
        for (const tid of this.state.visibleTrackIds) {
            if (tid === 'trk_esp' || tid === 'track_esp_charge') {
                active.push({ id: 'trk_esp', track_id: 'trk_esp', orbital_key: 'track_esp_charge', title: 'ESP Surface' });
                continue;
            }
            const live = poolByTid.get(tid);
            const tl = tlByTid.get(tid);
            // If frame data is loaded and orbital does not exist at frame f, skip cleanly
            if (fd && !live) continue;
            if (!fd && tl?.frame_map && Array.isArray(tl.frame_map) && !tl.frame_map[f]) continue;
            const src = live || tl;
            if (!src) continue;
            const isAnti = isAntiBondingOrbital(src);
            active.push({
                id: tid,
                track_id: tid,
                orbital_key: tid,
                key: tid,
                orbital_class: src.orbital_class || (isAnti ? 'antibonding' : 'bonding'),
                type: src.type,
                palette: src.palette,
                color: (tl && tl.color) || src.color || (isAnti ? '#fb923c' : '#38bdf8'),
                title: (tl && tl.title) || src.label || tid
            });
        }
        return active;
    },

    async refreshViewer() {
        const f = this.state.currentFrame, fd = this.state.frameDataCache[f];
        const coords = (fd && fd.coords) || (this.state.trajectory && this.state.trajectory[f]);
        if (!coords) return;
        const symbols = (fd && fd.symbols) || (this.state.trajectorySymbols && this.state.trajectorySymbols[f]) || this.state.uiSchema.symbols;
        const bonds = (fd && fd.bonds) || (this.state.trajectoryBonds && this.state.trajectoryBonds[f]) || null;
        const atomLabels = (fd && fd.atom_labels) || (this.state.trajectoryAtomLabels && this.state.trajectoryAtomLabels[f]) || null;
        const merged = this._getActiveOrbitals(f);
        const espActive = merged.some(t => t.orbital_key === 'track_esp_charge' || t.track_id === 'trk_esp');
        const moTracks = merged.filter(t => t.orbital_key !== 'track_esp_charge' && t.track_id !== 'trk_esp');
        const charges = (fd && fd.charges) || this.state.trajectoryCharges[f] || null;

        await this.state.viewer3D.renderScene({
            symbols, coords, charges, bonds, atomLabels, showCharges: this.state.showCharges || espActive,
            showLabels: this.state.showLabels,
            activeOrbitals: moTracks, isovalue: this.state.isovalue, frameIdx: f,
            hideAllOrbitals: this.state.hideAllOrbitals, fetchCubeFn: (fi, ok, ck, iso) => {
                const trk = moTracks.find(t => (t.local_key && t.local_key === ok) || t.orbital_key === ok || t.track_id === ok);
                const tid = trk ? (trk.track_id || trk.orbital_key) : ok;
                return this.fetchCube(fi, ok, ck, iso, tid);
            }
        });
    },

    seekFrame(frameIdx, isPlayback = false) {
        this.state.currentFrame = frameIdx;
        if (this.state.viewer3D) this.state.viewer3D.currentFrameIdx = frameIdx;
        if (this.state.timeline) this.state.timeline.setCurrentFrame(frameIdx);
        const scrubber = document.querySelector('#transport-frame-scrubber');
        if (scrubber) scrubber.value = frameIdx;
        const counter = document.querySelector('#transport-frame-counter');
        if (counter) counter.textContent = `Frame ${frameIdx} / ${Math.max(0, this.state.totalFrames - 1)}`;

        const coords = (this.state.trajectory && this.state.trajectory[frameIdx]) || (this.state.frameDataCache[frameIdx] && this.state.frameDataCache[frameIdx].coords);
        if (coords) {
            const symbols = (this.state.frameDataCache[frameIdx] && this.state.frameDataCache[frameIdx].symbols) || (this.state.trajectorySymbols && this.state.trajectorySymbols[frameIdx]) || this.state.uiSchema.symbols;
            const bonds = (this.state.frameDataCache[frameIdx] && this.state.frameDataCache[frameIdx].bonds) || (this.state.trajectoryBonds && this.state.trajectoryBonds[frameIdx]) || null;
            const atomLabels = (this.state.frameDataCache[frameIdx] && this.state.frameDataCache[frameIdx].atom_labels) || (this.state.trajectoryAtomLabels && this.state.trajectoryAtomLabels[frameIdx]) || null;
            const merged = this._getActiveOrbitals(frameIdx);
            const fd = this.state.frameDataCache[frameIdx];
            if (this.state.timeline && fd) this.state.timeline.updatePopulations(frameIdx, fd, symbols);
            const charges = (fd && fd.charges) || this.state.trajectoryCharges[frameIdx] || null;
            const espActive = merged.some(t => t.orbital_key === 'track_esp_charge' || t.track_id === 'trk_esp');
            const moTracks = this.state.hideAllOrbitals ? [] : merged.filter(t => t.orbital_key !== 'track_esp_charge' && t.track_id !== 'trk_esp');

            this.state.viewer3D.renderFastFrame(symbols, coords, moTracks, this.state.isovalue, frameIdx, charges, this.state.showCharges || espActive, bonds, this.state.showLabels, atomLabels);

            for (const orb of moTracks) {
                const tid = orb.track_id || orb.orbital_key || orb.key || orb;
                const ok = orb.local_key || orb.orbital_key || tid, ck = orb.chem_key || ok;
                if (!this.state.viewer3D.cubeCache[`${frameIdx}_${tid}`]) {
                    this.state.viewer3D.setStatus(`Frame ${frameIdx} | Computing orbital ${orb.title || tid}...`);
                    this.fetchCube(frameIdx, ok, ck, this.state.isovalue, tid).then(res => {
                        if (res && res.success && res.cube && this.state.currentFrame === frameIdx) {
                            this.state.viewer3D.renderFastFrame(symbols, coords, moTracks, this.state.isovalue, frameIdx, charges, this.state.showCharges || espActive, bonds, this.state.showLabels, atomLabels);
                            this.state.viewer3D.setStatus(`Frame ${frameIdx} | Ready`);
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

    stepFrame(delta) { this.seekFrame(Math.max(0, Math.min(this.state.totalFrames - 1, this.state.currentFrame + delta)), this.state.isPlaying); },

    togglePlay() {
        this.state.isPlaying = !this.state.isPlaying;
        const btn = document.querySelector('#btn-play-pause');
        if (btn) btn.textContent = this.state.isPlaying ? '⏸ Pause' : '▶ Play';
        if (this.state.isPlaying) {
            this.preloadAllCachedCubes();
            if (this.state.currentFrame >= this.state.totalFrames - 1) this.seekFrame(0, true);
            this.state.playTimer = setInterval(() => { (this.state.currentFrame >= this.state.totalFrames - 1) ? this.togglePlay() : this.stepFrame(1); }, 100);
        } else {
            clearInterval(this.state.playTimer);
            this.state.playTimer = null;
            this.loadFrameData(this.state.currentFrame);
        }
    },

    preloadAllCachedCubes() {
        if (this.state.precompute && this.state.viewer3D) {
            this.state.precompute.preloadCached(this.state.viewer3D.cubeCache, this.state.totalFrames, f => { if (f.includes(this.state.currentFrame)) this.refreshViewer(); });
        }
    },
    onAuditionOrbital(tid) { this.toggleVisibility(tid); },
    onAddToTimeline(orbital) {
        if (!this.state.timeline || !orbital) return;
        this.state.timeline.addTrack(orbital);
        const tid = (orbital.track_id === 'trk_esp' || orbital.orbital_key === 'track_esp_charge') ? 'trk_esp' : (orbital.track_id || orbital.orbital_key || orbital.key);
        if (tid) {
            this.state.visibleTrackIds.add(tid);
            this.state.timeline.syncVisibilitySet(this.state.visibleTrackIds);
            this.state.pool?.syncVisibilitySet(this.state.visibleTrackIds);
        }
        this.refreshViewer();
        this.preloadAllCachedCubes();
    },
    async saveTimeline(timeline) {
        const res = await postTimeline(this.state.stepId, timeline);
        if (res?.success && res?.timeline && this.state.timeline) this.state.timeline.setTimelineData(res.timeline, this.state.totalFrames);
    },
    getData() { return { timeline: this.state.timeline ? this.state.timeline.getTimelineObject() : {} }; }
};
