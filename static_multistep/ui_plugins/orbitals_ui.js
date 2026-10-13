export const OrbitalsUI = {
    state: {
        stepId: 'step_4_orbitals',
        uiSchema: null,
        totalFrames: 0,
        currentFrame: 0,
        landmarks: [],
        selectedOrbital: 'IBO_30',
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
        if (this.state.schedule.selected_orbital) {
            this.state.selectedOrbital = this.state.schedule.selected_orbital;
        }

        if (!uiSchema.has_megajson || this.state.totalFrames === 0) {
            container.innerHTML = `<div class="empty-state">No reaction pathway available. Please complete earlier steps first.</div>`;
            return;
        }

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 1.25rem;">
                <!-- Header with Title and Badges -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.75rem;">
                    <div>
                        <h2 style="font-size: 1.15rem; font-weight: 700; color: var(--text-main);">Reaction Probability Flow (IBO Movie)</h2>
                        <div style="font-size: 0.8rem; color: var(--text-muted); display: flex; gap: 0.5rem; margin-top: 0.25rem; flex-wrap: wrap;">
                            <span class="chip emerald">Intrinsic Bond Orbitals (IBO)</span>
                            <span class="chip blue">${uiSchema.symbols.length} atoms</span>
                            <span class="chip">${this.state.totalFrames} frames</span>
                            <span class="chip purple">Continuous Curly-Arrow Tracking</span>
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.6rem;">
                        <button id="btn-precompute-cubes" class="btn btn-secondary btn-sm" title="Precompute cubes across landmarks for instant 3D rendering">
                            ⚡ Precompute Landmark Cubes
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
                <div style="display: grid; grid-template-columns: 1fr 380px; gap: 1.25rem; align-items: start;">
                    <!-- 3D Viewport Card -->
                    <div style="display: flex; flex-direction: column; gap: 0.75rem;">
                        <div class="viewer-wrapper" style="height: 440px; position: relative;">
                            <div id="orb-viewport-3d" class="viewport-3d"></div>
                            <div id="orb-status-badge" class="viewer-overlay-badge">Frame 1 / ${this.state.totalFrames}</div>
                            <div id="orb-active-flow-badge" style="position: absolute; bottom: 12px; left: 12px; background: rgba(11, 15, 23, 0.85); backdrop-filter: blur(8px); border: 1px solid var(--border-color); border-radius: 6px; padding: 0.4rem 0.75rem; font-size: 0.78rem; font-weight: 600; color: #38bdf8; z-index: 5;">
                                Active: Loading...
                            </div>
                            <div id="orb-loading-spinner" style="display: none; position: absolute; inset: 0; background: rgba(11,15,23,0.65); backdrop-filter: blur(4px); justify-content: center; align-items: center; z-index: 10; color: #38bdf8; font-size: 0.9rem; font-weight: 600;">
                                <div style="display: flex; flex-direction: column; align-items: center; gap: 0.5rem;">
                                    <div class="spinner" style="width: 24px; height: 24px; border: 3px solid rgba(56,189,248,0.2); border-top-color: #38bdf8; border-radius: 50%; animation: spin 0.8s linear infinite;"></div>
                                    <span id="spinner-msg">Localizing IBO & Evaluating Volume...</span>
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
                        <!-- Primary Mechanism Flow Orbitals (The Curly Arrows!) -->
                        <div class="card" style="padding: 1rem;">
                            <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.6rem;">
                                <div style="font-size: 0.8rem; font-weight: 700; color: var(--text-muted); text-transform: uppercase;">Reaction Flow (Curly Arrows)</div>
                                <span class="chip emerald" style="font-size: 0.65rem;">Knizia IBO</span>
                            </div>
                            <div style="display: flex; flex-direction: column; gap: 0.5rem;" id="reaction-flow-buttons">
                                <button class="btn btn-sm ${this.state.selectedOrbital === 'IBO_30' ? 'btn-primary' : 'btn-secondary'} flow-btn" data-key="IBO_30" style="justify-content: space-between; text-align: left; padding: 0.45rem 0.65rem;">
                                    <div>
                                        <div style="font-weight: 700;">⚡ Forming C-O Bond</div>
                                        <div style="font-size: 0.72rem; opacity: 0.85;">LP(O9) → C1-O9 σ bond</div>
                                    </div>
                                    <span class="chip blue" style="font-size: 0.65rem;">Nucleophile</span>
                                </button>

                                <button class="btn btn-sm ${this.state.selectedOrbital === 'IBO_19' ? 'btn-primary' : 'btn-secondary'} flow-btn" data-key="IBO_19" style="justify-content: space-between; text-align: left; padding: 0.45rem 0.65rem;">
                                    <div>
                                        <div style="font-weight: 700;">⚡ Carbonyl Opening</div>
                                        <div style="font-size: 0.72rem; opacity: 0.85;">C1=O3 (π) → LP(O3⁻) oxyanion</div>
                                    </div>
                                    <span class="chip purple" style="font-size: 0.65rem;">Electrophile</span>
                                </button>

                                <button class="btn btn-sm ${this.state.selectedOrbital === 'IBO_20' ? 'btn-primary' : 'btn-secondary'} flow-btn" data-key="IBO_20" style="justify-content: space-between; text-align: left; padding: 0.45rem 0.65rem;">
                                    <div>
                                        <div style="font-weight: 700;">⚡ Leaving Group Departure</div>
                                        <div style="font-size: 0.72rem; opacity: 0.85;">C1-Cl2 σ bond → LP(Cl⁻)</div>
                                    </div>
                                    <span class="chip amber" style="font-size: 0.65rem;">Leaving Group</span>
                                </button>

                                <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.4rem; margin-top: 0.2rem;">
                                    <button class="btn btn-sm ${this.state.selectedOrbital === 'LUMO' ? 'btn-primary' : 'btn-secondary'} flow-btn" data-key="LUMO" style="padding: 0.35rem 0.5rem; justify-content: center;">
                                        <span style="font-weight: 600;">LUMO (Acceptor)</span>
                                    </button>
                                    <button class="btn btn-sm ${this.state.selectedOrbital === 'HOMO' ? 'btn-primary' : 'btn-secondary'} flow-btn" data-key="HOMO" style="padding: 0.35rem 0.5rem; justify-content: center;">
                                        <span style="font-weight: 600;">Canonical HOMO</span>
                                    </button>
                                </div>
                            </div>

                            <!-- Manual Localized Valence Orbital Chooser -->
                            <div style="margin-top: 0.8rem; border-top: 1px solid var(--border-color); padding-top: 0.75rem;">
                                <label style="font-size: 0.75rem; font-weight: 600; color: var(--text-muted); display: block; margin-bottom: 0.35rem;">Or Choose Any Localized Valence IBO Manually:</label>
                                <select id="manual-ibo-select" style="width: 100%; background: rgba(0,0,0,0.45); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; padding: 0.35rem 0.5rem; font-size: 0.75rem; font-family: var(--font-mono);">
                                    <option value="">Loading all localized IBOs...</option>
                                </select>
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
                            <div style="display: flex; flex-direction: column; gap: 0.6rem; font-size: 0.8rem;">
                                <div style="background: rgba(0,0,0,0.25); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.6rem;">
                                    <div style="font-weight: 600; color: #93c5fd; margin-bottom: 0.3rem;">Step 1: Addition (Reactants → TS₁ → Int)</div>
                                    <div style="display: flex; flex-direction: column; gap: 0.3rem;">
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s1-forming" checked> Forming C-O Bond (LP(O9) → C-O)
                                        </label>
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s1-opening" checked> Carbonyl Opening (C=O → LP)
                                        </label>
                                    </div>
                                </div>

                                <div style="background: rgba(0,0,0,0.25); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.6rem;">
                                    <div style="font-weight: 600; color: #c084fc; margin-bottom: 0.3rem;">Step 2: Elimination (Int → TS₂ → Prod)</div>
                                    <div style="display: flex; flex-direction: column; gap: 0.3rem;">
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s2-breaking" checked> Leaving Group (C-Cl → Cl⁻)
                                        </label>
                                        <label style="display: flex; align-items: center; gap: 0.35rem; cursor: pointer;">
                                            <input type="checkbox" id="sched-s2-carbonyl" checked> Carbonyl Reform (LP → C=O)
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
        const manualSelect = document.getElementById('manual-ibo-select');

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

        // Flow Buttons
        document.querySelectorAll('.flow-btn').forEach(btn => {
            btn.onclick = () => {
                const key = btn.dataset.key;
                this.state.selectedOrbital = key;
                document.querySelectorAll('.flow-btn').forEach(b => {
                    b.className = `btn btn-sm ${b.dataset.key === key ? 'btn-primary' : 'btn-secondary'} flow-btn`;
                });
                if (manualSelect) manualSelect.value = key;
                this.renderCurrentOrbitalMesh();
            };
        });

        // Manual IBO select dropdown
        if (manualSelect) {
            manualSelect.onchange = (e) => {
                const val = e.target.value;
                if (!val) return;
                this.state.selectedOrbital = val;
                document.querySelectorAll('.flow-btn').forEach(b => {
                    b.className = `btn btn-sm ${b.dataset.key === val ? 'btn-primary' : 'btn-secondary'} flow-btn`;
                });
                this.renderCurrentOrbitalMesh();
            };
        }

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
                    precomputeBtn.innerText = "⚡ Precompute Landmark Cubes";
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
                if (data.success) {
                    this.state.orbitalsMetaCache[frameIndex] = data;
                }
            } catch (e) {
                console.error("Failed to fetch frame orbitals:", e);
            }
        }

        const frameData = this.state.orbitalsMetaCache[frameIndex];
        if (frameData && frameData.all_ibos) {
            const manualSelect = document.getElementById('manual-ibo-select');
            if (manualSelect && manualSelect.options.length <= 1) {
                manualSelect.innerHTML = '<option value="">-- Choose any localized IBO --</option>';
                frameData.all_ibos.forEach(ib => {
                    const opt = document.createElement('option');
                    opt.value = ib.name;
                    opt.textContent = `${ib.name}: ${ib.label} (${ib.type.toUpperCase()})`;
                    if (ib.name === this.state.selectedOrbital) opt.selected = true;
                    manualSelect.appendChild(opt);
                });
            }
        }

        await this.renderCurrentOrbitalMesh(resetCamera);
    },

    async renderCurrentOrbitalMesh(resetCamera = false) {
        const frameIndex = this.state.currentFrame;
        const orb = this.state.selectedOrbital;
        const iso = this.state.isovalue;
        const cacheKey = `${frameIndex}_${orb}_${iso.toFixed(3)}`;

        // Update badge text
        const flowBadge = document.getElementById('orb-active-flow-badge');
        if (flowBadge) {
            const flowNames = {
                'IBO_30': 'Forming C-O Bond: LP(O9) → C-O',
                'IBO_19': 'Carbonyl Opening: C1=O3 (π) → LP(O3⁻)',
                'IBO_20': 'Leaving Group: C1-Cl2 → LP(Cl⁻)',
                'LUMO': 'Acceptor Virtual Space (LUMO)',
                'HOMO': 'Canonical HOMO'
            };
            flowBadge.innerText = `Active: ${flowNames[orb] || orb}`;
        }

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

        // 1. Primary rendering via 3Dmol VolumeData
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
        const schedule = {
            default_isovalue: this.state.isovalue,
            selected_orbital: this.state.selectedOrbital,
            stage_assignments: {
                "0": { orbitals: ["IBO_30", "IBO_19", "LUMO"], isovalue: this.state.isovalue },
                "1": { orbitals: ["IBO_30", "IBO_20", "LUMO"], isovalue: this.state.isovalue }
            },
            landmarks: this.state.landmarks
        };

        return {
            schedule: schedule
        };
    }
};
