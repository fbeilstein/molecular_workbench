export const BundleViewerUI = {
    state: {
        stepId: 'step_5_bundle_viewer',
        megajson: null,
        currentFrame: 0,
        isPlaying: false,
        playTimer: null,
        viewer: null,
        landmarkMeshes: {},
        activeOrbital: 'forming_bond',
        showOrbitals: true
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        const megajson = uiSchema.megajson;
        this.state.megajson = megajson;
        this.state.landmarkMeshes = megajson.landmark_orbital_meshes || {};

        if (!megajson || !megajson.frames) {
            container.innerHTML = `<div class="empty-state">No stitched pathway available. Please complete earlier steps first.</div>`;
            return;
        }

        const totalFrames = megajson.frames ? megajson.frames.length : 0;
        const boundaries = megajson.step_boundaries || [];
        const landmarks = uiSchema.landmarks || [];

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 1.25rem;">
                <!-- Header with Export Actions -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.75rem;">
                    <div>
                        <h2 style="font-size: 1.15rem; font-weight: 700; color: var(--text-main);">Master Multistep Reaction Bundle & Viewer</h2>
                        <div style="font-size: 0.8rem; color: var(--text-muted); display: flex; gap: 0.5rem; margin-top: 0.25rem; flex-wrap: wrap;">
                            <span class="chip blue">${megajson.symbols.length} atoms</span>
                            <span class="chip">${totalFrames} total frames</span>
                            <span class="chip amber">Overall Barrier ΔG‡: ${(megajson.overall_barrier_kcal || 0).toFixed(2)} kcal/mol</span>
                            <span class="chip emerald">Net Reaction ΔG°: ${(megajson.overall_reaction_energy_kcal || 0).toFixed(2)} kcal/mol</span>
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.6rem; flex-wrap: wrap;">
                        <button id="btn-export-bundle-zip" class="btn btn-primary btn-sm" style="background: linear-gradient(135deg, #10b981, #059669); border-color: #059669;">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                            </svg>
                            Export Bundle (.zip)
                        </button>
                        <button id="btn-export-xyz" class="btn btn-secondary btn-sm">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                            </svg>
                            Multi-Model XYZ
                        </button>
                        <button id="btn-export-json" class="btn btn-secondary btn-sm">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                            </svg>
                            MegaJSON
                        </button>
                    </div>
                </div>

                <!-- 3D Viewport with Orbital Overlay Controls -->
                <div class="viewer-wrapper" style="height: 420px; position: relative;">
                    <div id="bundle-viewport-3d" class="viewport-3d"></div>
                    <div id="bundle-status-badge" class="viewer-overlay-badge">Frame 1 / ${totalFrames}</div>
                    
                    <!-- Overlay orbital controls -->
                    <div style="position: absolute; top: 12px; right: 12px; background: rgba(11, 15, 23, 0.85); backdrop-filter: blur(8px); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.4rem 0.6rem; display: flex; align-items: center; gap: 0.6rem; z-index: 5; font-size: 0.75rem;">
                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer; color: var(--text-main); font-weight: 600;">
                            <input type="checkbox" id="bundle-toggle-orbitals" checked> Show Orbitals
                        </label>
                        <select id="bundle-select-orbital" style="background: rgba(0,0,0,0.5); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; padding: 0.2rem 0.4rem; font-size: 0.75rem;">
                            <option value="forming_bond">Forming C-O Bond (LP → σ)</option>
                            <option value="carbonyl_pi">Carbonyl π / Oxyanion LP</option>
                            <option value="leaving_group">Leaving Group C-Cl → Cl⁻</option>
                            <option value="LUMO">Acceptor LUMO (π*)</option>
                            <option value="HOMO">Canonical HOMO</option>
                        </select>
                    </div>
                </div>

                <!-- Playback Control Bar -->
                <div class="playback-bar">
                    <button id="btn-bundle-play" class="btn btn-sm btn-secondary" title="Play / Pause">
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                            <polygon points="5 3 19 12 5 21 5 3"/>
                        </svg>
                    </button>
                    <button id="btn-bundle-prev" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                    <input type="range" id="bundle-scrubber" class="scrub-slider" min="0" max="${Math.max(0, totalFrames - 1)}" value="0">
                    <button id="btn-bundle-next" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                    <div id="bundle-counter" class="frame-counter">1 / ${totalFrames}</div>
                </div>

                <!-- Synchronized Free Energy Canvas Diagram -->
                <div style="background: #080c14; border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.75rem;">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
                        <span style="font-size: 0.75rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase;">Synchronized Free Energy Profile (ΔG in Water):</span>
                        <span style="font-family: var(--font-mono); font-size: 0.75rem; color: #6ee7b7;" id="bundle-pes-energy-label">0.00 kcal/mol</span>
                    </div>
                    <canvas id="bundle-pes-canvas" width="1100" height="130" style="width: 100%; height: 130px; display: block; cursor: crosshair;"></canvas>
                </div>

                <!-- Step & Landmark Jump Shortcuts -->
                <div style="display: flex; align-items: center; gap: 0.6rem; flex-wrap: wrap;">
                    <span style="font-size: 0.8rem; font-weight: 600; color: var(--text-muted); text-transform: uppercase;">Landmarks:</span>
                    ${boundaries.map(b => `
                        <button class="btn btn-sm btn-secondary jump-step-btn" data-ts="${b.ts_frame}">
                            Step ${b.step_index} TS (f${b.ts_frame + 1})
                        </button>
                    `).join('')}
                    ${landmarks.map(lm => `
                        <button class="btn btn-sm btn-secondary jump-lm-btn" data-frame="${lm.frame_index}">
                            ${lm.name} (f${lm.frame_index + 1})
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
        const playBtn = document.getElementById('btn-bundle-play');
        const prevBtn = document.getElementById('btn-bundle-prev');
        const nextBtn = document.getElementById('btn-bundle-next');
        const scrubber = document.getElementById('bundle-scrubber');
        const canvas = document.getElementById('bundle-pes-canvas');
        const toggleOrb = document.getElementById('bundle-toggle-orbitals');
        const selectOrb = document.getElementById('bundle-select-orbital');

        if (playBtn) playBtn.onclick = () => this.togglePlay();
        if (prevBtn) prevBtn.onclick = () => this.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepFrame(1);
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pause();
                this.setFrame(parseInt(e.target.value, 10));
            };
        }

        if (toggleOrb) {
            toggleOrb.onchange = (e) => {
                this.state.showOrbitals = e.target.checked;
                this.setFrame(this.state.currentFrame);
            };
        }

        if (selectOrb) {
            selectOrb.onchange = (e) => {
                this.state.activeOrbital = e.target.value;
                this.setFrame(this.state.currentFrame);
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

        // Jump Buttons
        document.querySelectorAll('.jump-step-btn').forEach(btn => {
            btn.onclick = () => {
                this.pause();
                const tsFrame = parseInt(btn.dataset.ts, 10);
                this.setFrame(tsFrame);
            };
        });

        document.querySelectorAll('.jump-lm-btn').forEach(btn => {
            btn.onclick = () => {
                this.pause();
                const f = parseInt(btn.dataset.frame, 10);
                this.setFrame(f);
            };
        });

        // Export Standalone Bundle (.zip)
        const expBundle = document.getElementById('btn-export-bundle-zip');
        if (expBundle) {
            expBundle.onclick = async () => {
                expBundle.disabled = true;
                expBundle.innerText = "Packaging Bundle...";
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ action_type: 'export_bundle' })
                    });
                    const data = await res.json();
                    if (data.success && data.download_url) {
                        const a = document.createElement('a');
                        a.href = data.download_url;
                        a.download = data.filename || 'multistep_reaction_bundle.zip';
                        document.body.appendChild(a);
                        a.click();
                        document.body.removeChild(a);
                    }
                } catch (e) {
                    alert("Bundle export error: " + e.message);
                } finally {
                    expBundle.disabled = false;
                    expBundle.innerHTML = `
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                            <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                        </svg>
                        Export Bundle (.zip)
                    `;
                }
            };
        }

        // Export XYZ & JSON
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
        const container = document.getElementById('bundle-viewport-3d');
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
        const canvas = document.getElementById('bundle-pes-canvas');
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
        if (!frames.length) return;

        const energies = frames.map(f => f.energy_kcal !== undefined ? f.energy_kcal : (f.rel_energy_kcal || 0));
        const minE = Math.min(...energies);
        const maxE = Math.max(...energies);
        const rangeE = (maxE - minE) > 1e-3 ? (maxE - minE) : 1.0;

        const toX = (idx) => padLeft + (idx / (energies.length - 1)) * plotW;
        const toY = (val) => padTop + plotH - ((val - minE) / rangeE) * plotH;

        ctx.clearRect(0, 0, w, h);

        // Grid lines
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.06)';
        ctx.lineWidth = 1;
        for (let i = 0; i <= 4; i++) {
            const y = padTop + (plotH / 4) * i;
            ctx.beginPath();
            ctx.moveTo(padLeft, y);
            ctx.lineTo(padLeft + plotW, y);
            ctx.stroke();

            const eVal = maxE - (rangeE / 4) * i;
            ctx.fillStyle = '#6b7280';
            ctx.font = '10px Inter, sans-serif';
            ctx.textAlign = 'right';
            ctx.fillText(eVal.toFixed(1), padLeft - 6, y + 3);
        }

        // Reaction steps tint
        const boundaries = this.state.megajson.step_boundaries || [];
        boundaries.forEach((b, k) => {
            const xStart = toX(b.start_frame);
            const xEnd = toX(b.end_frame);
            ctx.fillStyle = k % 2 === 0 ? 'rgba(59, 130, 246, 0.04)' : 'rgba(168, 85, 247, 0.04)';
            ctx.fillRect(xStart, padTop, xEnd - xStart, plotH);
        });

        // PES curve
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

        // Fill gradient
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

        // Active frame cursor
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
            stick: { radius: 0.14, colorscheme: 'rasmol' },
            sphere: { scale: 0.28, colorscheme: 'rasmol' }
        });

        // If orbitals enabled, overlay cube VolumeData or mesh for landmark frames
        if (this.state.showOrbitals) {
            const frameStr = String(index);
            const orbName = this.state.activeOrbital;
            const landmarkCubes = this.state.megajson.landmark_cubes || {};
            const cubeKey = `${index}_${orbName}`;
            const cubeText = landmarkCubes[cubeKey];

            if (cubeText && typeof $3Dmol.VolumeData === 'function') {
                try {
                    const palettes = {
                        'forming_bond': { pos: '#06b6d4', neg: '#f97316' },
                        'carbonyl_pi': { pos: '#a855f7', neg: '#ec4899' },
                        'leaving_group': { pos: '#10b981', neg: '#84cc16' },
                        'LUMO': { pos: '#f59e0b', neg: '#fbbf24' },
                        'HOMO': { pos: '#38bdf8', neg: '#60a5fa' }
                    };
                    const pal = palettes[orbName] || { pos: '#2979ff', neg: '#ff9100' };
                    const vol = new $3Dmol.VolumeData(cubeText, 'cube');
                    this.state.viewer.addIsosurface(vol, { isoval: 0.03, color: pal.pos, alpha: 0.75, smoothness: 1 });
                    this.state.viewer.addIsosurface(vol, { isoval: -0.03, color: pal.neg, alpha: 0.75, smoothness: 1 });
                } catch (e) {
                    console.warn("VolumeData render error:", e);
                }
            } else {
                const meshes = this.state.landmarkMeshes[frameStr];
                if (meshes && meshes[orbName]) {
                    const m = meshes[orbName];
                    if (m.pos && (m.pos.vertexArr || m.pos.vertices)) {
                        this.state.viewer.addCustom({
                            vertexArr: m.pos.vertexArr || m.pos.vertices,
                            faceArr: m.pos.faceArr || m.pos.faces,
                            normalArr: m.pos.normalArr || undefined,
                            color: $3Dmol.CC.color('#2979ff'),
                            alpha: 0.75
                        });
                    }
                    if (m.neg && (m.neg.vertexArr || m.neg.vertices)) {
                        this.state.viewer.addCustom({
                            vertexArr: m.neg.vertexArr || m.neg.vertices,
                            faceArr: m.neg.faceArr || m.neg.faces,
                            normalArr: m.neg.normalArr || undefined,
                            color: $3Dmol.CC.color('#ff9100'),
                            alpha: 0.75
                        });
                    }
                }
            }
        }

        if (resetCamera) {
            this.state.viewer.zoomTo();
        }
        this.state.viewer.render();

        // Update UI
        const scrubber = document.getElementById('bundle-scrubber');
        if (scrubber) scrubber.value = index;
        const counter = document.getElementById('bundle-counter');
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;
        const badge = document.getElementById('bundle-status-badge');
        if (badge) {
            const stepNum = frame.step_index || 1;
            const tsBadge = frame.is_ts ? ' [Transition State]' : '';
            badge.innerText = `Frame ${index + 1} / ${frames.length} | Step ${stepNum}${tsBadge}`;
        }
        const eLabel = document.getElementById('bundle-pes-energy-label');
        if (eLabel) {
            const eVal = frame.energy_kcal !== undefined ? frame.energy_kcal : (frame.rel_energy_kcal || 0);
            eLabel.innerText = `${eVal >= 0 ? '+' : ''}${eVal.toFixed(2)} kcal/mol`;
        }

        this.drawPESChart();
    },

    stepFrame(delta) {
        this.pause();
        const frames = this.state.megajson.frames || [];
        const next = Math.max(0, Math.min(frames.length - 1, this.state.currentFrame + delta));
        this.setFrame(next);
    },

    togglePlay() {
        if (this.state.isPlaying) this.pause();
        else this.play();
    },

    play() {
        const frames = this.state.megajson.frames || [];
        if (this.state.currentFrame >= frames.length - 1) {
            this.setFrame(0);
        }

        this.state.isPlaying = true;
        const btn = document.getElementById('btn-bundle-play');
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
        const btn = document.getElementById('btn-bundle-play');
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
