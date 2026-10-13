export const MegaViewerUI = {
    state: {
        stepId: 'step_4_mega_viewer',
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

        if (!megajson || !megajson.frames) {
            container.innerHTML = `<div class="empty-state">No stitched pathway available. Please complete earlier steps first.</div>`;
            return;
        }

        const totalFrames = megajson.frames ? megajson.frames.length : 0;
        const boundaries = megajson.step_boundaries || [];
        const ep = megajson.energy_profile || {};

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 1.25rem;">
                <!-- Header with Export Actions -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.75rem;">
                    <div>
                        <h2 style="font-size: 1.15rem; font-weight: 700; color: var(--text-main);">Master Multistep Pathway Playback</h2>
                        <div style="font-size: 0.8rem; color: var(--text-muted); display: flex; gap: 0.5rem; margin-top: 0.25rem; flex-wrap: wrap;">
                            <span class="chip blue">${megajson.symbols.length} atoms</span>
                            <span class="chip">${totalFrames} total frames</span>
                            <span class="chip amber">Overall Barrier ΔG‡: ${(megajson.overall_barrier_kcal || 0).toFixed(2)} kcal/mol</span>
                            <span class="chip emerald">Net Reaction ΔG°: ${(megajson.overall_reaction_energy_kcal || 0).toFixed(2)} kcal/mol</span>
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
                <div class="viewer-wrapper" style="height: 420px;">
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

                <!-- Synchronized PES Canvas Diagram -->
                <div style="background: #080c14; border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.75rem;">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
                        <span style="font-size: 0.75rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase;">Synchronized Potential Energy Profile:</span>
                        <span style="font-family: var(--font-mono); font-size: 0.75rem; color: #6ee7b7;" id="mega-pes-energy-label">0.00 kcal/mol</span>
                    </div>
                    <canvas id="mega-pes-canvas" width="1100" height="130" style="width: 100%; height: 130px; display: block; cursor: crosshair;"></canvas>
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
        this.drawPESChart();
        this.setFrame(0, true);
    },

    initEventListeners() {
        const playBtn = document.getElementById('btn-mega-play');
        const prevBtn = document.getElementById('btn-mega-prev');
        const nextBtn = document.getElementById('btn-mega-next');
        const scrubber = document.getElementById('mega-scrubber');
        const canvas = document.getElementById('mega-pes-canvas');

        if (playBtn) playBtn.onclick = () => this.togglePlay();
        if (prevBtn) prevBtn.onclick = () => this.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepFrame(1);
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pause();
                this.setFrame(parseInt(e.target.value, 10));
            };
        }

        if (canvas) {
            const handleCanvasClick = (e) => {
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const total = (this.state.megajson.frames || []).length;
                if (total <= 1) return;
                const padLeft = 45;
                const padRight = 15;
                const plotW = rect.width - padLeft - padRight;
                const relX = Math.max(0, Math.min(plotW, x - padLeft));
                const targetFrame = Math.round((relX / plotW) * (total - 1));
                this.pause();
                this.setFrame(targetFrame);
            };

            canvas.onmousedown = (e) => {
                handleCanvasClick(e);
                const onMove = (ev) => handleCanvasClick(ev);
                const onUp = () => {
                    window.removeEventListener('mousemove', onMove);
                    window.removeEventListener('mouseup', onUp);
                };
                window.addEventListener('mousemove', onMove);
                window.addEventListener('mouseup', onUp);
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
                    this._downloadFile(data.xyz, data.filename || 'multistep_pathway.xyz', 'text/plain');
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
                    this._downloadFile(JSON.stringify(data.megajson, null, 2), data.filename || 'multistep_pathway.json', 'application/json');
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

    drawPESChart() {
        const canvas = document.getElementById('mega-pes-canvas');
        if (!canvas || !this.state.megajson) return;

        const ctx = canvas.getContext('2d');
        const dpr = window.devicePixelRatio || 1;
        const rect = canvas.getBoundingClientRect();
        canvas.width = rect.width * dpr;
        canvas.height = rect.height * dpr;
        ctx.scale(dpr, dpr);

        const w = rect.width;
        const h = rect.height;
        const padLeft = 45;
        const padRight = 15;
        const padTop = 15;
        const padBottom = 20;
        const plotW = w - padLeft - padRight;
        const plotH = h - padTop - padBottom;

        const frames = this.state.megajson.frames || [];
        if (frames.length === 0) return;

        const energies = frames.map(f => f.energy_kcal !== undefined ? f.energy_kcal : 0.0);
        let minE = Math.min(...energies);
        let maxE = Math.max(...energies);
        if (minE === maxE) { minE -= 5; maxE += 5; }
        const margin = (maxE - minE) * 0.15 || 2;
        minE -= margin;
        maxE += margin;

        const toX = (idx) => padLeft + (idx / (energies.length - 1)) * plotW;
        const toY = (val) => padTop + plotH - ((val - minE) / (maxE - minE)) * plotH;

        ctx.clearRect(0, 0, w, h);

        // Grid
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.08)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let i = 0; i <= 3; i++) {
            const yVal = minE + (i / 3) * (maxE - minE);
            const yPos = toY(yVal);
            ctx.moveTo(padLeft, yPos);
            ctx.lineTo(w - padRight, yPos);
        }
        ctx.stroke();

        // Labels
        ctx.fillStyle = '#94a3b8';
        ctx.font = '9px monospace';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 3; i++) {
            const yVal = minE + (i / 3) * (maxE - minE);
            ctx.fillText(yVal.toFixed(1), padLeft - 6, toY(yVal) + 3);
        }

        // Shading
        const boundaries = this.state.megajson.step_boundaries || [];
        boundaries.forEach((b, k) => {
            const xStart = toX(b.start_frame);
            const xEnd = toX(b.end_frame);
            ctx.fillStyle = k % 2 === 0 ? 'rgba(59, 130, 246, 0.04)' : 'rgba(168, 85, 247, 0.04)';
            ctx.fillRect(xStart, padTop, xEnd - xStart, plotH);
        });

        // Curve
        ctx.lineWidth = 2.5;
        ctx.beginPath();
        energies.forEach((val, i) => {
            const x = toX(i);
            const y = toY(val);
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        });
        ctx.strokeStyle = '#38bdf8';
        ctx.stroke();

        // Fill
        const grad = ctx.createLinearGradient(0, padTop, 0, padTop + plotH);
        grad.addColorStop(0, 'rgba(56, 189, 248, 0.2)');
        grad.addColorStop(1, 'rgba(56, 189, 248, 0.0)');
        ctx.lineTo(toX(energies.length - 1), toY(minE));
        ctx.lineTo(toX(0), toY(minE));
        ctx.closePath();
        ctx.fillStyle = grad;
        ctx.fill();

        // TS markers
        boundaries.forEach(b => {
            const tsIdx = b.ts_frame;
            if (tsIdx !== undefined && tsIdx < energies.length) {
                const x = toX(tsIdx);
                const y = toY(energies[tsIdx]);
                ctx.fillStyle = '#f59e0b';
                ctx.beginPath();
                ctx.arc(x, y, 4, 0, 2 * Math.PI);
                ctx.fill();
            }
        });

        // Cursor
        const activeX = toX(this.state.currentFrame);
        const activeY = toY(energies[this.state.currentFrame]);
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(activeX, padTop);
        ctx.lineTo(activeX, padTop + plotH);
        ctx.stroke();

        ctx.fillStyle = '#ef4444';
        ctx.beginPath();
        ctx.arc(activeX, activeY, 4.5, 0, 2 * Math.PI);
        ctx.fill();
    },

    setFrame(index, resetCamera = false) {
        if (!this.state.megajson || !this.state.viewer) return;
        const frames = this.state.megajson.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.state.currentFrame = index;
        const frame = frames[index];
        const symbols = this.state.megajson.symbols;
        const coords = frame.coords || frame.positions;

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

        const scrubber = document.getElementById('mega-scrubber');
        const counter = document.getElementById('mega-counter');
        const badge = document.getElementById('mega-status-badge');
        const pesLabel = document.getElementById('mega-pes-energy-label');

        const eKcal = frame.energy_kcal !== undefined ? frame.energy_kcal : 0.0;
        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;
        if (pesLabel) pesLabel.innerText = `Frame ${index + 1}: ${eKcal.toFixed(2)} kcal/mol`;

        if (badge) {
            const stepIdx = frame.step_index !== undefined ? frame.step_index : '?';
            const isTS = frame.is_ts ? ' ⭐ [TRANSITION STATE]' : '';
            badge.innerText = `Step ${stepIdx} | Frame ${index + 1} / ${frames.length} | Energy: ${eKcal.toFixed(2)} kcal/mol ${isTS}`;
        }

        this.drawPESChart();
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
