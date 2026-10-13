export const OrbitalsUI = {
    state: {
        stepId: 'step_4_orbitals',
        uiSchema: null,
        totalFrames: 0,
        currentFrame: 0,
        landmarks: [],
        selectedOrbital: 'HOMO',
        isovalue: 0.03,
        isPlaying: false,
        playTimer: null,
        viewer: null,
        meshCache: {},
        cubeCache: {},
        orbitalsMetaCache: {},
        schedule: null
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.uiSchema = uiSchema;
        this.state.totalFrames = uiSchema.total_frames || 0;
        this.state.landmarks = uiSchema.landmarks || [];
        this.state.schedule = JSON.parse(JSON.stringify(uiSchema.schedule || {}));
        this.state.isovalue = this.state.schedule.default_isovalue || 0.03;

        if (!uiSchema.has_megajson || this.state.totalFrames === 0) {
            container.innerHTML = `<div class="empty-state">No reaction pathway available. Please complete earlier steps first.</div>`;
            return;
        }

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 1.25rem;">
                <!-- Header with Title and Badges -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.75rem;">
                    <div>
                        <h2 style="font-size: 1.15rem; font-weight: 700; color: var(--text-main);">Molecular Orbital Studio & Resonance π Systems</h2>
                        <div style="font-size: 0.8rem; color: var(--text-muted); display: flex; gap: 0.5rem; margin-top: 0.25rem; flex-wrap: wrap;">
                            <span class="chip blue">${uiSchema.symbols.length} atoms</span>
                            <span class="chip">${this.state.totalFrames} frames</span>
                            <span class="chip emerald">Canonical π Resonance Active</span>
                            <span class="chip purple">Phase Overlap Continuous</span>
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.6rem;">
                        <button id="btn-precompute-cubes" class="btn btn-secondary btn-sm" title="Precompute cubes for all landmark states for instant rendering">
                            ⚡ Precompute All Cubes
                        </button>
                        <button id="btn-dl-cube" class="btn btn-secondary btn-sm" title="Download Gaussian .cube for current frame and orbital">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/>
                            </svg>
                            Export .cube
                        </button>
                    </div>
                </div>

                <!-- Main Layout: 3D Viewport on Left, Control Panel on Right -->
                <div style="display: grid; grid-template-columns: 1fr 370px; gap: 1.25rem; align-items: start;">
                    <!-- 3D Viewport Card -->
                    <div style="display: flex; flex-direction: column; gap: 0.75rem;">
                        <div class="viewer-wrapper" style="height: 440px; position: relative;">
                            <div id="orb-viewport-3d" class="viewport-3d"></div>
                            <div id="orb-status-badge" class="viewer-overlay-badge">Frame 1 / ${this.state.totalFrames}</div>
                            <div id="orb-loading-spinner" style="display: none; position: absolute; inset: 0; background: rgba(11,15,23,0.65); backdrop-filter: blur(4px); justify-content: center; align-items: center; z-index: 10; color: #38bdf8; font-size: 0.9rem; font-weight: 600;">
                                <div style="display: flex; flex-direction: column; align-items: center; gap: 0.5rem;">
                                    <div class="spinner" style="width: 24px; height: 24px; border: 3px solid rgba(56,189,248,0.2); border-top-color: #38bdf8; border-radius: 50%; animation: spin 0.8s linear infinite;"></div>
                                    <span id="spinner-msg">Evaluating Wavefunction & Volume...</span>
                                </div>
                            </div>
                        </div>

                        <!-- Playback Bar -->
                        <div class="playback-bar">
                            <button id="btn-orb-play" class="btn btn-sm btn-secondary" title="Play / Pause">
                                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                                    <polygon points="5 3 19 12 5 21 5 3"/>
                                </svg>
                            </button>
                            <button id="btn-orb-prev" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                            <input type="range" id="orb-scrubber" class="scrub-slider" min="0" max="${Math.max(0, this.state.totalFrames - 1)}" value="0">
                            <button id="btn-orb-next" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                            <div id="orb-counter" class="frame-counter">1 / ${this.state.totalFrames}</div>
                        </div>

                        <!-- Landmark Jump Bar -->
                        <div style="background: rgba(17, 24, 39, 0.7); border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.6rem 0.8rem; display: flex; align-items: center; gap: 0.5rem; flex-wrap: wrap;">
                            <span style="font-size: 0.75rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase;">Landmarks:</span>
                            ${this.state.landmarks.map(lm => `
                                <button class="btn btn-sm btn-secondary orb-landmark-btn" data-frame="${lm.frame_index}" style="padding: 0.2rem 0.6rem; font-size: 0.75rem;">
                                    ${lm.name} (f${lm.frame_index + 1})
                                </button>
                            `).join('')}
                        </div>
                    </div>

                    <!-- Right Controls Panel -->
                    <div style="display: flex; flex-direction: column; gap: 1rem;">
                        <!-- Active Orbital Selection with Chemical Resonance Labels -->
                        <div class="card" style="padding: 1rem;">
                            <div style="font-size: 0.8rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase; margin-bottom: 0.6rem;">Chemical Orbitals & Resonance π</div>
                            <div style="display: flex; flex-direction: column; gap: 0.45rem;" id="orbital-chips-container">
                                <button class="btn btn-sm btn-primary orb-chip-btn" data-orb="HOMO" style="justify-content: space-between; padding: 0.45rem 0.65rem; text-align: left;">
                                    <div>
                                        <div style="font-weight: 700;" id="orb-label-HOMO">HOMO: Formate π (Nucleophile lone pair)</div>
                                    </div>
                                    <span style="font-size: 0.75rem; font-family: var(--font-mono); opacity: 0.9;" id="orb-energy-HOMO">-- eV</span>
                                </button>

                                <button class="btn btn-sm btn-secondary orb-chip-btn" data-orb="LUMO" style="justify-content: space-between; padding: 0.45rem 0.65rem; text-align: left;">
                                    <div>
                                        <div style="font-weight: 700;" id="orb-label-LUMO">LUMO: Carbonyl π* (Acceptor antibond)</div>
                                    </div>
                                    <span style="font-size: 0.75rem; font-family: var(--font-mono); opacity: 0.9;" id="orb-energy-LUMO">-- eV</span>
                                </button>

                                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.45rem;">
                                    <button class="btn btn-sm btn-secondary orb-chip-btn" data-orb="HOMO-1" style="justify-content: space-between; padding: 0.4rem 0.5rem;">
                                        <span style="font-weight: 600;">HOMO-1</span>
                                        <span style="font-size: 0.7rem; font-family: var(--font-mono);" id="orb-energy-HOMO-1">-- eV</span>
                                    </button>

                                    <button class="btn btn-sm btn-secondary orb-chip-btn" data-orb="LUMO+1" style="justify-content: space-between; padding: 0.4rem 0.5rem;">
                                        <span style="font-weight: 600;">LUMO+1</span>
                                        <span style="font-size: 0.7rem; font-family: var(--font-mono);" id="orb-energy-LUMO+1">-- eV</span>
                                    </button>
                                </div>

                                <button class="btn btn-sm btn-secondary orb-chip-btn" data-orb="PI_RES" style="justify-content: space-between; padding: 0.4rem 0.65rem; text-align: left;">
                                    <div>
                                        <span style="font-weight: 600;">Formate π₁ (Resonance backbone)</span>
                                    </div>
                                    <span style="font-size: 0.75rem; font-family: var(--font-mono); opacity: 0.9;" id="orb-energy-PI_RES">-- eV</span>
                                </button>
                            </div>

                            <div style="display: flex; align-items: center; justify-content: space-between; margin-top: 0.8rem; font-size: 0.75rem; color: var(--text-muted); background: rgba(0,0,0,0.25); padding: 0.35rem 0.6rem; border-radius: 4px;">
                                <div style="display: flex; align-items: center; gap: 0.35rem;">
                                    <span style="display: inline-block; width: 10px; height: 10px; border-radius: 2px; background: #2979ff;"></span>
                                    <span>Positive lobe (+ψ)</span>
                                </div>
                                <div style="display: flex; align-items: center; gap: 0.35rem;">
                                    <span style="display: inline-block; width: 10px; height: 10px; border-radius: 2px; background: #ff9100;"></span>
                                    <span>Negative lobe (-ψ)</span>
                                </div>
                            </div>
                        </div>

                        <!-- Isovalue Threshold Slider -->
                        <div class="card" style="padding: 1rem;">
                            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem;">
                                <label style="font-size: 0.8rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase;">Isovalue Threshold</label>
                                <span id="isovalue-val-display" style="font-family: var(--font-mono); font-size: 0.8rem; font-weight: 700; color: #38bdf8;">${this.state.isovalue.toFixed(3)}</span>
                            </div>
                            <input type="range" id="isovalue-slider" min="0.010" max="0.080" step="0.005" value="${this.state.isovalue}" style="width: 100%; cursor: pointer;">
                            <div style="display: flex; justify-content: space-between; font-size: 0.7rem; color: var(--text-dim); margin-top: 0.25rem;">
                                <span>0.010 (Diffuse)</span>
                                <span>0.030 (Standard)</span>
                                <span>0.080 (Dense)</span>
                            </div>
                        </div>

                        <!-- Schedule Configuration -->
                        <div class="card" style="padding: 1rem;">
                            <div style="font-size: 0.8rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase; margin-bottom: 0.6rem;">Pathway Visibility Schedule</div>
                            <p style="font-size: 0.75rem; color: var(--text-muted); margin-bottom: 0.8rem;">
                                Active orbitals featured during each reaction stage:
                            </p>
                            
                            <div style="display: flex; flex-direction: column; gap: 0.6rem; font-size: 0.8rem;">
                                <div style="background: rgba(0,0,0,0.25); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.6rem;">
                                    <div style="font-weight: 600; color: #93c5fd; margin-bottom: 0.3rem;">Step 1: Addition (Reactants → TS₁ → Int)</div>
                                    <div style="display: flex; gap: 1rem;">
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s1-homo" checked> HOMO (Formate π)
                                        </label>
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s1-lumo" checked> LUMO (Carbonyl π*)
                                        </label>
                                    </div>
                                </div>

                                <div style="background: rgba(0,0,0,0.25); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.6rem;">
                                    <div style="font-weight: 600; color: #c084fc; margin-bottom: 0.3rem;">Step 2: Elimination (Int → TS₂ → Prod)</div>
                                    <div style="display: flex; gap: 1rem;">
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s2-homo" checked> HOMO (Oxyanion lone pair)
                                        </label>
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s2-lumo" checked> LUMO (C-Cl σ*)
                                        </label>
                                    </div>
                                </div>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        this.initEventListeners();
        this.init3DViewer();
        await this.loadFrame(0, true);
    },

    initEventListeners() {
        const playBtn = document.getElementById('btn-orb-play');
        const prevBtn = document.getElementById('btn-orb-prev');
        const nextBtn = document.getElementById('btn-orb-next');
        const scrubber = document.getElementById('orb-scrubber');
        const isoSlider = document.getElementById('isovalue-slider');
        const isoDisplay = document.getElementById('isovalue-val-display');

        if (playBtn) playBtn.onclick = () => this.togglePlay();
        if (prevBtn) prevBtn.onclick = () => this.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepFrame(1);
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pause();
                this.loadFrame(parseInt(e.target.value, 10));
            };
        }

        if (isoSlider) {
            isoSlider.oninput = (e) => {
                const val = parseFloat(e.target.value);
                this.state.isovalue = val;
                if (isoDisplay) isoDisplay.innerText = val.toFixed(3);
            };
            isoSlider.onchange = () => {
                this.renderCurrentOrbitalMesh();
            };
        }

        // Landmark Jump Buttons
        document.querySelectorAll('.orb-landmark-btn').forEach(btn => {
            btn.onclick = () => {
                this.pause();
                const f = parseInt(btn.dataset.frame, 10);
                this.loadFrame(f);
            };
        });

        // Orbital Chip Selector Buttons
        document.querySelectorAll('.orb-chip-btn').forEach(btn => {
            btn.onclick = () => {
                const orb = btn.dataset.orb;
                this.state.selectedOrbital = orb;
                document.querySelectorAll('.orb-chip-btn').forEach(b => {
                    b.className = `btn btn-sm ${b.dataset.orb === orb ? 'btn-primary' : 'btn-secondary'} orb-chip-btn`;
                });
                this.renderCurrentOrbitalMesh();
            };
        });

        // Precompute All Cubes Button
        const precomputeBtn = document.getElementById('btn-precompute-cubes');
        if (precomputeBtn) {
            precomputeBtn.onclick = async () => {
                precomputeBtn.disabled = true;
                precomputeBtn.innerText = "⏳ Computing Cubes...";
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ action_type: 'pregenerate_all_cubes' })
                    });
                    const d = await res.json();
                    if (d.success) {
                        precomputeBtn.innerText = `✓ ${d.cached_total} Cubes Ready`;
                    }
                } catch (e) {
                    precomputeBtn.innerText = "⚡ Precompute Cubes";
                } finally {
                    precomputeBtn.disabled = false;
                }
            };
        }

        // Download Gaussian Cube
        const cubeBtn = document.getElementById('btn-dl-cube');
        if (cubeBtn) {
            cubeBtn.onclick = async () => {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'get_cube',
                        frame_index: this.state.currentFrame,
                        orbital: this.state.selectedOrbital
                    })
                });
                const data = await res.json();
                if (data.success && data.cube) {
                    this._downloadFile(data.cube, data.filename || 'orbital.cube', 'text/plain');
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
        const container = document.getElementById('orb-viewport-3d');
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

    async loadFrame(frameIndex, resetCamera = false) {
        if (frameIndex < 0 || frameIndex >= this.state.totalFrames) return;
        this.state.currentFrame = frameIndex;

        // Update UI controls
        const scrubber = document.getElementById('orb-scrubber');
        if (scrubber) scrubber.value = frameIndex;
        const counter = document.getElementById('orb-counter');
        if (counter) counter.innerText = `${frameIndex + 1} / ${this.state.totalFrames}`;
        const badge = document.getElementById('orb-status-badge');
        if (badge) {
            const lm = this.state.landmarks.find(l => l.frame_index === frameIndex);
            badge.innerText = lm ? `${lm.name} (Frame ${frameIndex + 1})` : `Frame ${frameIndex + 1} / ${this.state.totalFrames}`;
        }

        // Fetch frame orbitals metadata if not cached
        if (!this.state.orbitalsMetaCache[frameIndex]) {
            try {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_orbitals', frame_index: frameIndex })
                });
                const data = await res.json();
                if (data.success && data.orbitals) {
                    this.state.orbitalsMetaCache[frameIndex] = data;
                }
            } catch (e) {
                console.error("Failed to fetch frame orbitals:", e);
            }
        }

        // Update orbital energy labels on buttons
        const frameData = this.state.orbitalsMetaCache[frameIndex];
        if (frameData && frameData.orbitals) {
            ['HOMO-1', 'HOMO', 'LUMO', 'LUMO+1', 'PI_RES'].forEach(orb => {
                const span = document.getElementById(`orb-energy-${orb}`);
                if (span && frameData.orbitals[orb]) {
                    span.innerText = `${frameData.orbitals[orb].energy_ev.toFixed(2)} eV`;
                }
                const labelDiv = document.getElementById(`orb-label-${orb}`);
                if (labelDiv && frameData.orbitals[orb] && frameData.orbitals[orb].label) {
                    labelDiv.innerText = frameData.orbitals[orb].label;
                }
            });
        }

        await this.renderCurrentOrbitalMesh(resetCamera);
    },

    async renderCurrentOrbitalMesh(resetCamera = false) {
        const frameIndex = this.state.currentFrame;
        const orb = this.state.selectedOrbital;
        const iso = this.state.isovalue;
        const cacheKey = `${frameIndex}_${orb}_${iso.toFixed(3)}`;

        let cubeText = this.state.cubeCache[cacheKey];
        let mesh = this.state.meshCache[cacheKey];

        if (!cubeText && !mesh) {
            const spinner = document.getElementById('orb-loading-spinner');
            if (spinner) spinner.style.display = 'flex';
            try {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'get_orbital_mesh',
                        frame_index: frameIndex,
                        orbital: orb,
                        isovalue: iso
                    })
                });
                const data = await res.json();
                if (data.success) {
                    if (data.cube) {
                        cubeText = data.cube;
                        this.state.cubeCache[cacheKey] = cubeText;
                    }
                    if (data.mesh) {
                        mesh = data.mesh;
                        this.state.meshCache[cacheKey] = mesh;
                    }
                }
            } catch (e) {
                console.error("Failed to generate orbital mesh:", e);
            } finally {
                if (spinner) spinner.style.display = 'none';
            }
        }

        if (!this.state.viewer) return;

        const symbols = this.state.uiSchema.symbols;
        const frameData = this.state.orbitalsMetaCache[frameIndex];
        const coords = frameData ? frameData.coords : null;
        if (!coords) return;

        this.state.viewer.clear();
        const xyz = this._buildXYZ(symbols, coords);
        this.state.viewer.addModel(xyz, 'xyz');
        this.state.viewer.setStyle({}, {
            stick: { radius: 0.14, colorscheme: 'rasmol' },
            sphere: { scale: 0.28, colorscheme: 'rasmol' }
        });

        // 1. Primary rendering via 3Dmol VolumeData (identical to Molecular Explorer)
        if (cubeText && typeof $3Dmol.VolumeData === 'function') {
            try {
                const vol = new $3Dmol.VolumeData(cubeText, 'cube');
                this.state.viewer.addIsosurface(vol, {
                    isoval: iso,
                    color: '#2979ff',
                    alpha: 0.75,
                    smoothness: 1
                });
                this.state.viewer.addIsosurface(vol, {
                    isoval: -iso,
                    color: '#ff9100',
                    alpha: 0.75,
                    smoothness: 1
                });
            } catch (volErr) {
                console.warn("VolumeData render fallback to mesh:", volErr);
            }
        }
        // 2. Fallback rendering via addCustom if mesh is present
        else if (mesh) {
            if (mesh.pos && mesh.pos.vertexArr && mesh.pos.vertexArr.length > 0) {
                this.state.viewer.addCustom({
                    vertexArr: mesh.pos.vertexArr,
                    faceArr: mesh.pos.faceArr,
                    normalArr: mesh.pos.normalArr || undefined,
                    color: $3Dmol.CC.color('#2979ff'),
                    alpha: 0.75
                });
            }
            if (mesh.neg && mesh.neg.vertexArr && mesh.neg.vertexArr.length > 0) {
                this.state.viewer.addCustom({
                    vertexArr: mesh.neg.vertexArr,
                    faceArr: mesh.neg.faceArr,
                    normalArr: mesh.neg.normalArr || undefined,
                    color: $3Dmol.CC.color('#ff9100'),
                    alpha: 0.75
                });
            }
        }

        if (resetCamera) {
            this.state.viewer.zoomTo();
        }
        this.state.viewer.render();
    },

    stepFrame(delta) {
        this.pause();
        const next = Math.max(0, Math.min(this.state.totalFrames - 1, this.state.currentFrame + delta));
        this.loadFrame(next);
    },

    togglePlay() {
        if (this.state.isPlaying) this.pause();
        else this.play();
    },

    play() {
        if (this.state.currentFrame >= this.state.totalFrames - 1) {
            this.loadFrame(0);
        }
        this.state.isPlaying = true;
        const btn = document.getElementById('btn-orb-play');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.state.playTimer);
        this.state.playTimer = setInterval(() => {
            const next = this.state.currentFrame + 1;
            if (next >= this.state.totalFrames) {
                this.pause();
                return;
            }
            this.loadFrame(next);
        }, 120);
    },

    pause() {
        this.state.isPlaying = false;
        clearInterval(this.state.playTimer);
        const btn = document.getElementById('btn-orb-play');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    },

    async getData(uiSchema) {
        const s1Homo = document.getElementById('sched-s1-homo');
        const s1Lumo = document.getElementById('sched-s1-lumo');
        const s2Homo = document.getElementById('sched-s2-homo');
        const s2Lumo = document.getElementById('sched-s2-lumo');

        const s1Orbs = [];
        if (s1Homo && s1Homo.checked) s1Orbs.push('HOMO');
        if (s1Lumo && s1Lumo.checked) s1Orbs.push('LUMO');

        const s2Orbs = [];
        if (s2Homo && s2Homo.checked) s2Orbs.push('HOMO');
        if (s2Lumo && s2Lumo.checked) s2Orbs.push('LUMO');

        const schedule = {
            default_isovalue: this.state.isovalue,
            stage_assignments: {
                "0": { orbitals: s1Orbs.length ? s1Orbs : ['HOMO', 'LUMO'], isovalue: this.state.isovalue },
                "1": { orbitals: s2Orbs.length ? s2Orbs : ['HOMO', 'LUMO'], isovalue: this.state.isovalue }
            },
            landmarks: this.state.landmarks
        };

        return {
            schedule: schedule
        };
    }
};
