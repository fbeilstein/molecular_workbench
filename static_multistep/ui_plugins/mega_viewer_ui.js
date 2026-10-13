export const MegaViewerUI = {
    state: {
        stepId: 'step_3_mega_viewer',
        megajson: null,
        currentFrame: 0,
        isPlaying: false,
        playTimer: null,
        viewer: null
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        const megajson = uiSchema.megajson;
        this.state.megajson = megajson;

        if (!megajson) {
            container.innerHTML = `<div class="empty-state">No stitched pathway available. Please complete Step 2.</div>`;
            return;
        }

        const totalFrames = megajson.frames ? megajson.frames.length : 0;
        const boundaries = megajson.step_boundaries || [];

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 1.25rem;">
                <!-- Header with Export Actions -->
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div>
                        <h2 style="font-size: 1.15rem; font-weight: 700; color: var(--text-main);">Master Reaction Pathway Playback</h2>
                        <div style="font-size: 0.8rem; color: var(--text-muted); display: flex; gap: 0.5rem; margin-top: 0.25rem;">
                            <span class="chip blue">${megajson.symbols.length} atoms</span>
                            <span class="chip">${totalFrames} total frames</span>
                            <span class="chip amber">Overall Barrier: ${megajson.overall_barrier_kcal.toFixed(2)} kcal/mol</span>
                            <span class="chip emerald">Net ΔE: ${megajson.overall_reaction_energy_kcal.toFixed(2)} kcal/mol</span>
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.6rem;">
                        <button id="btn-export-xyz" class="btn btn-secondary btn-sm">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                            </svg>
                            Export Multi-Model XYZ
                        </button>
                        <button id="btn-export-json" class="btn btn-primary btn-sm">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                            </svg>
                            Download MegaJSON
                        </button>
                    </div>
                </div>

                <!-- 3D Viewport -->
                <div class="viewer-wrapper" style="height: 480px;">
                    <div id="mega-viewport-3d" class="viewport-3d"></div>
                    <div id="mega-status-badge" class="viewer-overlay-badge">Frame 1 / ${totalFrames}</div>
                </div>

                <!-- Playback Control Bar -->
                <div class="playback-bar">
                    <button id="btn-mega-play" class="btn btn-sm btn-secondary" title="Play / Pause">
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                            <polygon points="5 3 19 12 5 21 5 3"/>
                        </svg>
                    </button>
                    <button id="btn-mega-prev" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                    <input type="range" id="mega-scrubber" class="scrub-slider" min="0" max="${Math.max(0, totalFrames - 1)}" value="0">
                    <button id="btn-mega-next" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                    <div id="mega-counter" class="frame-counter">1 / ${totalFrames}</div>
                </div>

                <!-- Step Jump Shortcuts -->
                <div style="display: flex; align-items: center; gap: 0.6rem; flex-wrap: wrap;">
                    <span style="font-size: 0.8rem; font-weight: 600; color: var(--text-muted); text-transform: uppercase;">Step Jumps:</span>
                    ${boundaries.map(b => `
                        <button class="btn btn-sm btn-secondary jump-step-btn" data-start="${b.start_frame}" data-ts="${b.ts_frame}">
                            Step ${b.step_index}: ${b.step_name} (TS: Frame ${b.ts_frame + 1})
                        </button>
                    `).join('')}
                </div>
            </div>
        `;

        this.initEventListeners();
        this.init3DViewer();
        this.setFrame(0, true);
    },

    initEventListeners() {
        const playBtn = document.getElementById('btn-mega-play');
        const prevBtn = document.getElementById('btn-mega-prev');
        const nextBtn = document.getElementById('btn-mega-next');
        const scrubber = document.getElementById('mega-scrubber');

        if (playBtn) playBtn.onclick = () => this.togglePlay();
        if (prevBtn) prevBtn.onclick = () => this.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepFrame(1);
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pause();
                this.setFrame(parseInt(e.target.value, 10));
            };
        }

        // Step Jump Buttons
        document.querySelectorAll('.jump-step-btn').forEach(btn => {
            btn.onclick = () => {
                this.pause();
                const tsFrame = parseInt(btn.dataset.ts, 10);
                this.setFrame(tsFrame);
            };
        });

        // Export Actions
        const expXyz = document.getElementById('btn-export-xyz');
        const expJson = document.getElementById('btn-export-json');

        if (expXyz) {
            expXyz.onclick = async () => {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_xyz' })
                });
                const data = await res.json();
                if (data.success && data.xyz) {
                    this._downloadFile(data.xyz, data.filename || 'pathway.xyz', 'text/plain');
                }
            };
        }

        if (expJson) {
            expJson.onclick = async () => {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_megajson' })
                });
                const data = await res.json();
                if (data.success && data.megajson) {
                    this._downloadFile(JSON.stringify(data.megajson, null, 2), data.filename || 'pathway.json', 'application/json');
                }
            };
        }
    },

    _downloadFile(content, filename, type) {
        const blob = new Blob([content], { type });
        const url = URL.createObjectURL(blob);
        const a = document.createElement('a');
        a.href = url;
        a.download = filename;
        document.body.appendChild(a);
        a.click();
        document.body.removeChild(a);
        URL.revokeObjectURL(url);
    },

    init3DViewer() {
        const container = document.getElementById('mega-viewport-3d');
        if (!container) return;

        container.innerHTML = '';
        this.state.viewer = $3Dmol.createViewer(container, {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: '0x0d121c'
        });
    },

    _buildXYZ(symbols, coords) {
        const n = symbols.length;
        let xyz = `${n}\nFrame\n`;
        for (let i = 0; i < n; i++) {
            const sym = symbols[i];
            const [x, y, z] = coords[i];
            xyz += `${sym} ${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)}\n`;
        }
        return xyz;
    },

    setFrame(index, resetCamera = false) {
        if (!this.state.megajson || !this.state.viewer) return;
        const frames = this.state.megajson.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.state.currentFrame = index;
        const frame = frames[index];
        const symbols = this.state.megajson.symbols;
        const coords = frame.coords;

        const xyz = this._buildXYZ(symbols, coords);
        this.state.viewer.clear();
        this.state.viewer.addModel(xyz, 'xyz');
        this.state.viewer.setStyle({}, {
            stick: { radius: 0.14 },
            sphere: { scale: 0.25 }
        });

        if (resetCamera) {
            this.state.viewer.zoomTo();
        }
        this.state.viewer.render();

        // Update UI
        const scrubber = document.getElementById('mega-scrubber');
        const counter = document.getElementById('mega-counter');
        const badge = document.getElementById('mega-status-badge');

        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;

        if (badge) {
            const stepIdx = frame.step_index !== undefined ? frame.step_index : '?';
            const energy = frame.energy_kcal !== undefined ? frame.energy_kcal.toFixed(2) : '0.00';
            const isTS = frame.is_ts ? ' ⭐ [TRANSITION STATE]' : '';
            badge.innerText = `Step ${stepIdx} | Frame ${index + 1} / ${frames.length} | Energy: ${energy} kcal/mol ${isTS}`;
        }
    },

    stepFrame(delta) {
        this.pause();
        const frames = this.state.megajson ? this.state.megajson.frames : [];
        if (!frames || frames.length === 0) return;
        let target = this.state.currentFrame + delta;
        target = Math.max(0, Math.min(frames.length - 1, target));
        this.setFrame(target);
    },

    togglePlay() {
        if (this.state.isPlaying) {
            this.pause();
        } else {
            this.play();
        }
    },

    play() {
        if (!this.state.megajson) return;
        const frames = this.state.megajson.frames || [];
        if (frames.length === 0) return;

        // If at the end, restart from the beginning
        if (this.state.currentFrame >= frames.length - 1) {
            this.setFrame(0);
        }

        this.state.isPlaying = true;
        const btn = document.getElementById('btn-mega-play');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.state.playTimer);
        this.state.playTimer = setInterval(() => {
            const next = this.state.currentFrame + 1;
            if (next >= frames.length) {
                this.pause();
                return;
            }
            this.setFrame(next);
        }, 100);
    },

    pause() {
        this.state.isPlaying = false;
        clearInterval(this.state.playTimer);
        const btn = document.getElementById('btn-mega-play');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    },


    async getData(uiSchema) {
        return { completed: true };
    }
};
