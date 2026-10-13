export const StitchAlignUI = {
    state: {
        stepId: 'step_2_stitch_align',
        megajson: null,
        junctionPreviews: [],
        activeJunctionIdx: 0,
        junctionChoices: {},
        nBridgeFrames: 8,
        viewMode: 'transition', // 'transition', 'overlay', 'blink', 'prev_only', 'next_only', 'unaligned'
        transitionFrameIdx: 0,
        isTransitionPlaying: false,
        transitionTimer: null,
        blinkTimer: null,
        blinkState: 0,
        viewer: null
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.activeJunctionIdx = 0;
        this.stopBlink();
        this.pauseTransition();

        container.innerHTML = `
            <div class="choose-container">
                <!-- Left Panel: Junctions, Bridge Interpolation & Choices -->
                <div class="left-panel">
                    <div class="queue-card" style="flex: 1;">
                        <div class="section-title">
                            <span>Pathway Junctions</span>
                            <span class="chip blue" id="junctions-count-badge">Analyzing...</span>
                        </div>
                        <div id="stitch-status-banner" style="background: rgba(59, 130, 246, 0.1); border: 1px solid rgba(59, 130, 246, 0.25); border-radius: var(--radius-sm); padding: 0.65rem 0.85rem; font-size: 0.8rem; margin-bottom: 0.5rem;">
                            Detecting intermediate molecule overlaps and calculating Hungarian Kabsch alignment...
                        </div>
                        <div id="junctions-list" class="queue-list" style="max-height: 280px; overflow-y: auto;"></div>
                    </div>

                    <!-- Bridge Frame Interpolation Settings Card -->
                    <div class="catalog-card" style="flex: 0 0 auto;">
                        <div class="section-title">
                            <span>Intermittent Bridge Frames</span>
                            <span class="chip emerald" id="n-bridge-val">${this.state.nBridgeFrames} Frames</span>
                        </div>
                        <p style="font-size: 0.75rem; color: var(--text-muted); margin-bottom: 0.6rem;">
                            Generates smooth cubic smoothstep transition frames connecting Step K product to Step K+1 reactant to stitch pathways without geometric discontinuity.
                        </p>
                        <div style="display: flex; align-items: center; gap: 0.75rem;">
                            <span style="font-size: 0.75rem; font-family: var(--font-mono); color: var(--text-muted);">0</span>
                            <input type="range" id="n-bridge-slider" class="scrub-slider" min="0" max="20" step="1" value="${this.state.nBridgeFrames}" style="flex: 1;">
                            <span style="font-size: 0.75rem; font-family: var(--font-mono); color: var(--text-muted);">20</span>
                        </div>
                    </div>

                    <!-- Global Metrics Summary -->
                    <div class="catalog-card" style="flex: 0 0 auto;">
                        <div class="section-title">Global Master Pathway Metrics</div>
                        <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; font-size: 0.8rem;">
                            <div class="detail-box">
                                <div class="detail-box-title">Total Atoms</div>
                                <div id="metric-atoms" style="font-weight: 700; color: #93c5fd; font-family: var(--font-mono);">-</div>
                            </div>
                            <div class="detail-box">
                                <div class="detail-box-title">Total Frames</div>
                                <div id="metric-frames" style="font-weight: 700; color: #93c5fd; font-family: var(--font-mono);">-</div>
                            </div>
                            <div class="detail-box">
                                <div class="detail-box-title">Overall Barrier</div>
                                <div id="metric-barrier" style="font-weight: 700; color: #fcd34d; font-family: var(--font-mono);">-</div>
                            </div>
                            <div class="detail-box">
                                <div class="detail-box-title">Net Reaction Energy</div>
                                <div id="metric-rxn-energy" style="font-weight: 700; color: #6ee7b7; font-family: var(--font-mono);">-</div>
                            </div>
                        </div>
                    </div>
                </div>

                <!-- Right Panel: Interactive 3D Boundary & Transition Viewer -->
                <div class="preview-panel">
                    <div class="preview-header">
                        <div>
                            <div class="preview-rxn-title" id="junction-preview-title">Junction Alignment & Transition</div>
                            <div class="preview-meta-chips" id="junction-preview-chips"></div>
                        </div>
                        <!-- Legend -->
                        <div style="display: flex; gap: 0.6rem; align-items: center; font-size: 0.75rem; font-family: var(--font-mono);">
                            <span style="display: flex; align-items: center; gap: 0.3rem;">
                                <span style="display: inline-block; width: 10px; height: 10px; background: #38bdf8; border-radius: 50%;"></span>
                                Step K Product
                            </span>
                            <span style="display: flex; align-items: center; gap: 0.3rem;">
                                <span style="display: inline-block; width: 10px; height: 10px; background: #a855f7; border-radius: 50%;"></span>
                                Bridge
                            </span>
                            <span style="display: flex; align-items: center; gap: 0.3rem;">
                                <span style="display: inline-block; width: 10px; height: 10px; background: #f43f5e; border-radius: 50%;"></span>
                                Step K+1 Reactant
                            </span>
                        </div>
                    </div>

                    <!-- 3D Viewport -->
                    <div class="viewer-wrapper">
                        <div id="junction-viewport-3d" class="viewport-3d"></div>
                        <div id="junction-status-badge" class="viewer-overlay-badge">Loading Alignment...</div>
                    </div>

                    <!-- Alignment View Mode Toolbar -->
                    <div class="playback-bar" style="justify-content: space-between; flex-wrap: wrap; gap: 0.5rem;">
                        <div style="display: flex; gap: 0.4rem; align-items: center; flex-wrap: wrap;">
                            <span style="font-size: 0.75rem; font-weight: 600; color: var(--text-muted); text-transform: uppercase;">Mode:</span>
                            <button class="btn btn-sm btn-primary btn-view-mode" data-mode="transition" title="Smoothly animate across boundary bridge frames">
                                🎬 Play Transition
                            </button>
                            <button class="btn btn-sm btn-secondary btn-view-mode" data-mode="overlay" title="Superimpose both molecules to inspect fit">
                                🔄 Overlay Aligned
                            </button>
                            <button class="btn btn-sm btn-secondary btn-view-mode" data-mode="blink" title="Toggle back and forth to inspect atomic jump">
                                ⚡ Blink Toggle
                            </button>
                            <button class="btn btn-sm btn-secondary btn-view-mode" data-mode="prev_only" title="Show Step K final frame only">
                                🔵 Step K Product
                            </button>
                            <button class="btn btn-sm btn-secondary btn-view-mode" data-mode="next_only" title="Show Step K+1 initial frame only">
                                🔴 Step K+1 Reactant
                            </button>
                            <button class="btn btn-sm btn-secondary btn-view-mode" data-mode="unaligned" title="View before Kabsch alignment">
                                ⚠️ Raw Unaligned
                            </button>
                        </div>
                        <div id="alignment-rmsd-display" style="font-family: var(--font-mono); font-size: 0.8rem; font-weight: 700; color: #6ee7b7;">
                            RMSD: 0.0000 Å
                        </div>
                    </div>

                    <!-- Dedicated Transition Bridge Player Bar -->
                    <div class="playback-bar" id="transition-controls-bar" style="background: rgba(15, 23, 42, 0.75); border: 1px solid rgba(148, 163, 184, 0.2);">
                        <button id="btn-trans-play" class="btn btn-sm btn-secondary" title="Play / Pause Transition">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                                <polygon points="5 3 19 12 5 21 5 3"/>
                            </svg>
                        </button>
                        <button id="btn-trans-prev" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                        <input type="range" id="trans-scrubber" class="scrub-slider" min="0" max="9" value="0">
                        <button id="btn-trans-next" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                        <div id="trans-counter" class="frame-counter" style="min-width: 140px; text-align: right; color: #a855f7;">Frame 1 / 10</div>
                    </div>

                    <!-- Intermediate Scaffold Details -->
                    <div class="preview-details">
                        <div class="detail-box">
                            <div class="detail-box-title">Step K Final Product Scaffold</div>
                            <div id="j-prev-info" style="font-family: var(--font-mono); color: #93c5fd; font-size: 0.8rem;">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Step K+1 Initial Reactant Scaffold</div>
                            <div id="j-next-info" style="font-family: var(--font-mono); color: #fda4af; font-size: 0.8rem;">-</div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        this.initEventListeners();
        await this.runPreviewStitch();
    },

    initEventListeners() {
        // View mode buttons
        document.querySelectorAll('.btn-view-mode').forEach(btn => {
            btn.onclick = () => {
                document.querySelectorAll('.btn-view-mode').forEach(b => {
                    b.classList.remove('btn-primary');
                    b.classList.add('btn-secondary');
                });
                btn.classList.remove('btn-secondary');
                btn.classList.add('btn-primary');
                this.setViewMode(btn.dataset.mode);
            };
        });

        // Intermittent bridge frames slider
        const bridgeSlider = document.getElementById('n-bridge-slider');
        const bridgeValBadge = document.getElementById('n-bridge-val');
        if (bridgeSlider) {
            bridgeSlider.oninput = (e) => {
                const val = parseInt(e.target.value, 10);
                this.state.nBridgeFrames = val;
                if (bridgeValBadge) bridgeValBadge.innerText = `${val} Frames`;
            };
            bridgeSlider.onchange = (e) => {
                const val = parseInt(e.target.value, 10);
                this.state.nBridgeFrames = val;
                this.runPreviewStitch();
            };
        }

        // Transition player controls
        const playBtn = document.getElementById('btn-trans-play');
        const prevBtn = document.getElementById('btn-trans-prev');
        const nextBtn = document.getElementById('btn-trans-next');
        const scrubber = document.getElementById('trans-scrubber');

        if (playBtn) playBtn.onclick = () => this.toggleTransitionPlay();
        if (prevBtn) prevBtn.onclick = () => this.stepTransitionFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepTransitionFrame(1);
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pauseTransition();
                this.setTransitionFrame(parseInt(e.target.value, 10));
            };
        }
    },

    async runPreviewStitch(customChoices = null) {
        if (customChoices) {
            Object.assign(this.state.junctionChoices, customChoices);
        }

        const banner = document.getElementById('stitch-status-banner');
        try {
            const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'preview_stitch',
                    junction_choices: this.state.junctionChoices,
                    n_bridge_frames: this.state.nBridgeFrames
                })
            });

            const data = await res.json();
            if (!data.success || !data.megajson) {
                if (banner) {
                    banner.className = 'alert-banner alert-error';
                    banner.innerText = 'Alignment error: ' + (data.error || 'Unknown error');
                }
                return;
            }

            this.state.megajson = data.megajson;
            this.state.junctionPreviews = data.junction_previews || [];
            const previews = this.state.junctionPreviews;

            const badge = document.getElementById('junctions-count-badge');
            if (badge) badge.innerText = `${previews.length} Junction${previews.length > 1 ? 's' : ''}`;

            if (banner) {
                const maxRMSD = previews.reduce((acc, p) => Math.max(acc, p.rmsd || 0), 0);
                banner.style.background = 'rgba(16, 185, 129, 0.12)';
                banner.style.borderColor = 'rgba(16, 185, 129, 0.3)';
                banner.innerHTML = `<span style="color: #6ee7b7; font-weight: 600;">✓ Optimal Hungarian Alignment:</span> Maximum junction RMSD is <strong>${maxRMSD.toFixed(4)} Å</strong> with ${this.state.nBridgeFrames} transition frames.`;
            }

            // Global Metrics
            const mAtoms = document.getElementById('metric-atoms');
            const mFrames = document.getElementById('metric-frames');
            const mBarrier = document.getElementById('metric-barrier');
            const mRxn = document.getElementById('metric-rxn-energy');
            if (mAtoms) mAtoms.innerText = `${data.megajson.symbols.length} atoms`;
            if (mFrames) mFrames.innerText = `${data.megajson.frames.length} frames`;
            if (mBarrier) mBarrier.innerText = `${data.megajson.overall_barrier_kcal.toFixed(2)} kcal/mol`;
            if (mRxn) mRxn.innerText = `${data.megajson.overall_reaction_energy_kcal.toFixed(2)} kcal/mol`;

            this.renderJunctionsList();
            this.renderActiveJunction3D(this.state.activeJunctionIdx);

        } catch (e) {
            if (banner) {
                banner.className = 'alert-banner alert-error';
                banner.innerText = 'Network error: ' + e.message;
            }
        }
    },

    renderJunctionsList() {
        const list = document.getElementById('junctions-list');
        if (!list) return;

        const previews = this.state.junctionPreviews;
        if (previews.length === 0) {
            list.innerHTML = `<div class="empty-state">No consecutive junctions to align.</div>`;
            return;
        }

        list.innerHTML = previews.map((j, idx) => {
            const isSelected = idx === this.state.activeJunctionIdx;
            const candidates = j.all_candidates || [];
            const hasMultiple = candidates.length > 1;
            const isPerfect = (j.rmsd || 0) < 0.1;

            return `
                <div class="queue-item junction-card ${isSelected ? 'selected-junction' : ''}" data-idx="${idx}" style="flex-direction: column; align-items: stretch; gap: 0.5rem; cursor: pointer; border-left: 3px solid ${isSelected ? 'var(--primary)' : 'transparent'};">
                    <div style="display: flex; justify-content: space-between; align-items: center;">
                        <span class="queue-step-badge">Junction ${idx + 1}</span>
                        <span class="chip ${isPerfect ? 'emerald' : 'amber'}" style="font-family: var(--font-mono);">
                            RMSD: ${(j.rmsd || 0).toFixed(4)} Å
                        </span>
                    </div>

                    <div style="font-size: 0.8rem; font-weight: 600; color: var(--text-main); font-family: var(--font-mono);">
                        Step ${j.step_prev_index} ➔ Step ${j.step_next_index}
                    </div>

                    <!-- Coinciding Molecules Selector if 2 or more coincide -->
                    ${hasMultiple ? `
                        <div style="margin-top: 0.25rem;" onclick="event.stopPropagation();">
                            <div style="font-size: 0.7rem; color: #94a3b8; font-weight: 600; margin-bottom: 0.2rem;">
                                Coinciding Intermediate To Align:
                            </div>
                            <select class="junction-candidate-select" data-jkey="${j.junction_index}" style="width: 100%; padding: 0.35rem 0.5rem; background: var(--bg-surface); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; font-size: 0.775rem;">
                                ${candidates.map(c => `
                                    <option value="${c.id}" ${c.id === j.selected_candidate_id ? 'selected' : ''}>
                                        ${c.label} ${c.is_joint ? '⭐' : ''}
                                    </option>
                                `).join('')}
                            </select>
                        </div>
                    ` : `
                        <div style="font-size: 0.75rem; color: var(--text-muted); display: flex; gap: 0.5rem;">
                            <span>Hero: <strong style="color: #93c5fd;">${j.hero_formula || 'Scaffold'}</strong></span>
                            <span>(${j.n_matched_atoms} atoms)</span>
                        </div>
                    `}
                </div>
            `;
        }).join('');

        // Card click -> switch active junction
        list.querySelectorAll('.junction-card').forEach(el => {
            el.onclick = () => {
                const idx = parseInt(el.dataset.idx, 10);
                this.state.activeJunctionIdx = idx;
                this.renderJunctionsList();
                this.renderActiveJunction3D(idx);
            };
        });

        // Dropdown change -> re-stitch with new candidate selection
        list.querySelectorAll('.junction-candidate-select').forEach(sel => {
            sel.onchange = (e) => {
                const jkey = sel.dataset.jkey;
                const chosenId = sel.value;
                this.runPreviewStitch({ [jkey]: chosenId });
            };
        });
    },

    init3DViewer() {
        const container = document.getElementById('junction-viewport-3d');
        if (!container) return;

        if (!this.state.viewer) {
            container.innerHTML = '';
            this.state.viewer = $3Dmol.createViewer(container, {
                defaultcolors: $3Dmol.rasmolElementColors,
                backgroundColor: '0x0d121c'
            });
        }
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

    getTransitionSequence(j) {
        if (!j) return [];
        const seq = [
            {
                coords: j.prev_coords,
                symbols: j.prev_symbols,
                type: 'start',
                label: `Step ${j.step_prev_index} Product (Arrival Frame)`,
                color: '#38bdf8'
            }
        ];

        const bridges = j.bridge_frames || [];
        bridges.forEach((b, i) => {
            seq.push({
                coords: b.positions || b.coords,
                symbols: j.prev_symbols,
                type: 'bridge',
                label: `Bridge Frame ${i + 1} / ${bridges.length}`,
                energy: b.energy_kcal,
                color: '#a855f7'
            });
        });

        seq.push({
            coords: j.next_coords_aligned,
            symbols: j.next_symbols,
            type: 'end',
            label: `Step ${j.step_next_index} Reactant (Aligned Start)`,
            color: '#f43f5e'
        });

        return seq;
    },

    renderActiveJunction3D(idx, resetCamera = false) {
        this.init3DViewer();
        if (!this.state.viewer) return;

        const previews = this.state.junctionPreviews;
        if (!previews || idx >= previews.length) return;

        const j = previews[idx];

        // Update header details
        const title = document.getElementById('junction-preview-title');
        const chips = document.getElementById('junction-preview-chips');
        const rmsdLabel = document.getElementById('alignment-rmsd-display');
        const pInfo = document.getElementById('j-prev-info');
        const nInfo = document.getElementById('j-next-info');

        if (title) title.innerText = `Junction ${idx + 1}: Step ${j.step_prev_index} ➔ Step ${j.step_next_index}`;
        if (chips) {
            chips.innerHTML = `
                <span class="chip blue">${j.n_matched_atoms} Matched Atoms</span>
                <span class="chip emerald">RMSD: ${(j.rmsd || 0).toFixed(4)} Å</span>
                <span class="chip purple">${(j.bridge_frames || []).length} Bridge Frames</span>
                <span class="chip">${j.hero_formula || 'Intermediate Scaffold'}</span>
            `;
        }
        if (rmsdLabel) rmsdLabel.innerText = `RMSD: ${(j.rmsd || 0).toFixed(4)} Å`;
        if (pInfo) pInfo.innerText = `${j.step_prev_name} (${j.prev_symbols.length} atoms)`;
        if (nInfo) nInfo.innerText = `${j.step_next_name} (${j.next_symbols.length} atoms)`;

        // Update transition scrubber limits
        const seq = this.getTransitionSequence(j);
        const scrubber = document.getElementById('trans-scrubber');
        if (scrubber) {
            scrubber.max = Math.max(0, seq.length - 1);
            if (this.state.transitionFrameIdx >= seq.length) {
                this.state.transitionFrameIdx = 0;
            }
            scrubber.value = this.state.transitionFrameIdx;
        }

        this.apply3DScene(j, this.state.viewMode, resetCamera);
    },

    apply3DScene(j, mode, resetCamera = false) {
        if (!this.state.viewer) return;
        const badge = document.getElementById('junction-status-badge');
        const transControls = document.getElementById('transition-controls-bar');

        if (mode === 'transition') {
            if (transControls) transControls.style.display = 'flex';
            this.renderTransitionFrame(j, this.state.transitionFrameIdx, resetCamera);
            return;
        }

        if (transControls) transControls.style.display = 'none';
        this.pauseTransition();
        this.state.viewer.clear();

        const xyzPrev = this._buildXYZ(j.prev_symbols, j.prev_coords);
        const xyzNextAligned = this._buildXYZ(j.next_symbols, j.next_coords_aligned);
        const xyzNextUnaligned = this._buildXYZ(j.next_symbols, j.next_coords_unaligned);

        if (mode === 'overlay') {
            const m1 = this.state.viewer.addModel(xyzPrev, 'xyz');
            const m2 = this.state.viewer.addModel(xyzNextAligned, 'xyz');

            this.state.viewer.setStyle({ model: m1 }, {
                stick: { color: '#38bdf8', radius: 0.12 },
                sphere: { color: '#38bdf8', scale: 0.22 }
            });
            this.state.viewer.setStyle({ model: m2 }, {
                stick: { color: '#f43f5e', radius: 0.12 },
                sphere: { color: '#f43f5e', scale: 0.22 }
            });

            if (badge) badge.innerText = `Overlay Aligned | RMSD: ${(j.rmsd || 0).toFixed(4)} Å (Cyan = Step K, Magenta = Step K+1)`;

        } else if (mode === 'prev_only') {
            const m1 = this.state.viewer.addModel(xyzPrev, 'xyz');
            this.state.viewer.setStyle({ model: m1 }, {
                stick: { color: '#38bdf8', radius: 0.14 },
                sphere: { color: '#38bdf8', scale: 0.26 }
            });
            if (badge) badge.innerText = `Step ${j.step_prev_index} Final Product Frame`;

        } else if (mode === 'next_only') {
            const m2 = this.state.viewer.addModel(xyzNextAligned, 'xyz');
            this.state.viewer.setStyle({ model: m2 }, {
                stick: { color: '#f43f5e', radius: 0.14 },
                sphere: { color: '#f43f5e', scale: 0.26 }
            });
            if (badge) badge.innerText = `Step ${j.step_next_index} Initial Reactant Frame (Aligned)`;

        } else if (mode === 'unaligned') {
            const m1 = this.state.viewer.addModel(xyzPrev, 'xyz');
            const m2 = this.state.viewer.addModel(xyzNextUnaligned, 'xyz');
            this.state.viewer.setStyle({ model: m1 }, { stick: { color: '#38bdf8', radius: 0.12 }, sphere: { color: '#38bdf8', scale: 0.22 } });
            this.state.viewer.setStyle({ model: m2 }, { stick: { color: '#fbbf24', radius: 0.12 }, sphere: { color: '#fbbf24', scale: 0.22 } });
            if (badge) badge.innerText = `Raw Unaligned Coordinates (Yellow = Step K+1 before fit)`;

        } else if (mode === 'blink') {
            const isPrev = this.state.blinkState === 0;
            const xyzActive = isPrev ? xyzPrev : xyzNextAligned;
            const col = isPrev ? '#38bdf8' : '#f43f5e';
            const m = this.state.viewer.addModel(xyzActive, 'xyz');
            this.state.viewer.setStyle({ model: m }, {
                stick: { color: col, radius: 0.14 },
                sphere: { color: col, scale: 0.26 }
            });
            if (badge) badge.innerText = isPrev ? `[Blink] Step ${j.step_prev_index} End` : `[Blink] Step ${j.step_next_index} Start (RMSD: ${(j.rmsd || 0).toFixed(4)} Å)`;
        }

        if (resetCamera) {
            this.state.viewer.zoomTo();
        }
        this.state.viewer.render();
    },

    renderTransitionFrame(j, frameIdx, resetCamera = false) {
        if (!this.state.viewer || !j) return;
        const seq = this.getTransitionSequence(j);
        if (seq.length === 0) return;

        const idx = Math.max(0, Math.min(seq.length - 1, frameIdx));
        this.state.transitionFrameIdx = idx;
        const item = seq[idx];

        this.state.viewer.clear();
        const xyz = this._buildXYZ(item.symbols, item.coords);
        const m = this.state.viewer.addModel(xyz, 'xyz');

        if (item.type === 'start') {
            this.state.viewer.setStyle({ model: m }, {
                stick: { color: '#38bdf8', radius: 0.14 },
                sphere: { color: '#38bdf8', scale: 0.25 }
            });
        } else if (item.type === 'end') {
            this.state.viewer.setStyle({ model: m }, {
                stick: { color: '#f43f5e', radius: 0.14 },
                sphere: { color: '#f43f5e', scale: 0.25 }
            });
        } else {
            // Bridge morph: highlight with purple/indigo stick
            this.state.viewer.setStyle({ model: m }, {
                stick: { color: '#c084fc', radius: 0.13 },
                sphere: { scale: 0.25 }
            });
        }

        if (resetCamera) {
            this.state.viewer.zoomTo();
        }
        this.state.viewer.render();

        // Update transition UI controls
        const scrubber = document.getElementById('trans-scrubber');
        const counter = document.getElementById('trans-counter');
        const badge = document.getElementById('junction-status-badge');

        if (scrubber) scrubber.value = idx;
        if (counter) counter.innerText = `${item.label} (${idx + 1} / ${seq.length})`;
        if (badge) {
            badge.innerText = `🎬 [Transition Player] ${item.label} | ${idx + 1} / ${seq.length}`;
        }
    },

    setTransitionFrame(idx) {
        const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
        if (j) this.renderTransitionFrame(j, idx, false);
    },

    stepTransitionFrame(delta) {
        this.pauseTransition();
        const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
        if (!j) return;
        const seq = this.getTransitionSequence(j);
        let target = this.state.transitionFrameIdx + delta;
        target = Math.max(0, Math.min(seq.length - 1, target));
        this.renderTransitionFrame(j, target, false);
    },

    toggleTransitionPlay() {
        if (this.state.isTransitionPlaying) {
            this.pauseTransition();
        } else {
            this.playTransition();
        }
    },

    playTransition() {
        const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
        if (!j) return;
        const seq = this.getTransitionSequence(j);
        if (seq.length <= 1) return;

        if (this.state.transitionFrameIdx >= seq.length - 1) {
            this.state.transitionFrameIdx = 0;
            this.renderTransitionFrame(j, 0, false);
        }

        this.state.isTransitionPlaying = true;
        const btn = document.getElementById('btn-trans-play');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.state.transitionTimer);
        this.state.transitionTimer = setInterval(() => {
            const next = this.state.transitionFrameIdx + 1;
            if (next >= seq.length) {
                this.pauseTransition();
                return;
            }
            this.renderTransitionFrame(j, next, false);
        }, 120);
    },

    pauseTransition() {
        this.state.isTransitionPlaying = false;
        clearInterval(this.state.transitionTimer);
        const btn = document.getElementById('btn-trans-play');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    },

    setViewMode(mode) {
        this.state.viewMode = mode;
        this.stopBlink();
        this.pauseTransition();

        if (mode === 'blink') {
            this.startBlink();
        } else {
            const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
            if (j) this.apply3DScene(j, mode, false);
        }
    },

    startBlink() {
        this.stopBlink();
        this.state.blinkState = 0;
        const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
        if (j) this.apply3DScene(j, 'blink', false);

        this.state.blinkTimer = setInterval(() => {
            this.state.blinkState = 1 - this.state.blinkState;
            const currJ = this.state.junctionPreviews[this.state.activeJunctionIdx];
            if (currJ) this.apply3DScene(currJ, 'blink', false);
        }, 500);
    },

    stopBlink() {
        if (this.state.blinkTimer) {
            clearInterval(this.state.blinkTimer);
            this.state.blinkTimer = null;
        }
    },

    async getData(uiSchema) {
        this.stopBlink();
        this.pauseTransition();
        return {
            confirmed: true,
            junction_choices: this.state.junctionChoices,
            n_bridge_frames: this.state.nBridgeFrames
        };
    }
};
