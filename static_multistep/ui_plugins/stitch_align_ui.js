import { Stitch3DViewer } from './stitch/stitch_viewer_3d.js';
import { renderJunctionsQueue, updateGlobalMetrics } from './stitch/stitch_junctions_dom.js';

export const StitchAlignUI = {
    state: {
        stepId: 'step_2_stitch_align',
        megajson: null,
        junctionPreviews: [],
        activeJunctionIdx: 0,
        junctionChoices: {},
        nBridgeFrames: 8,
        viewer: null
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.activeJunctionIdx = 0;
        this.state.viewer = new Stitch3DViewer();

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

                    <!-- View Mode Selection Buttons -->
                    <div style="display: flex; gap: 0.4rem; background: var(--bg-surface); padding: 0.4rem 0.6rem; border-radius: var(--radius-sm); border: 1px solid var(--border-color); margin-bottom: 0.6rem; flex-wrap: wrap;">
                        <button class="btn btn-xs btn-mode ${this.state.viewer.viewMode === 'transition' ? 'btn-primary' : 'btn-secondary'}" data-mode="transition">
                            🎬 Smoothstep Transition
                        </button>
                        <button class="btn btn-xs btn-mode ${this.state.viewer.viewMode === 'overlay' ? 'btn-primary' : 'btn-secondary'}" data-mode="overlay">
                            🔍 Superposition Overlay
                        </button>
                        <button class="btn btn-xs btn-mode ${this.state.viewer.viewMode === 'blink' ? 'btn-primary' : 'btn-secondary'}" data-mode="blink">
                            💡 Alternating Blink
                        </button>
                        <button class="btn btn-xs btn-mode ${this.state.viewer.viewMode === 'prev_only' ? 'btn-primary' : 'btn-secondary'}" data-mode="prev_only">
                            ⬅ Step K Final Only
                        </button>
                        <button class="btn btn-xs btn-mode ${this.state.viewer.viewMode === 'next_only' ? 'btn-primary' : 'btn-secondary'}" data-mode="next_only">
                            ➡ Step K+1 Start Only
                        </button>
                        <button class="btn btn-xs btn-mode ${this.state.viewer.viewMode === 'unaligned' ? 'btn-primary' : 'btn-secondary'}" data-mode="unaligned">
                            ⚠️ Raw Unaligned
                        </button>
                    </div>

                    <!-- 3Dmol Boundary Viewport -->
                    <div class="viewer-wrapper">
                        <div id="junction-viewport-3d" class="viewport-3d"></div>
                        <div id="junction-status-badge" class="viewer-overlay-badge">Computing Boundary Overlap...</div>
                    </div>

                    <!-- Transition Playback Controller -->
                    <div class="playback-bar" id="transition-controls-bar" style="display: flex;">
                        <button id="btn-trans-play" class="btn btn-sm btn-secondary">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                                <polygon points="5 3 19 12 5 21 5 3"/>
                            </svg>
                        </button>
                        <button id="btn-trans-prev" class="btn btn-sm btn-secondary">⏮</button>
                        <input type="range" id="trans-scrubber" class="scrub-slider" min="0" max="10" value="0">
                        <button id="btn-trans-next" class="btn btn-sm btn-secondary">⏭</button>
                        <div id="trans-counter" class="frame-counter">1 / 10</div>
                    </div>

                    <!-- Details Inspector Grid -->
                    <div class="preview-details" style="display: grid; grid-template-columns: 1fr 1fr 1fr; margin-top: 0.6rem;">
                        <div class="detail-box">
                            <div class="detail-box-title">Step K Product</div>
                            <div id="j-prev-info" style="font-family: var(--font-mono); font-size: 0.75rem; color: #93c5fd;">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Aligned Boundary Fit</div>
                            <div id="alignment-rmsd-display" style="font-weight: 700; color: #6ee7b7; font-family: var(--font-mono);">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Step K+1 Reactant</div>
                            <div id="j-next-info" style="font-family: var(--font-mono); font-size: 0.75rem; color: #f43f5e;">-</div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        this.attachEvents();
        this.state.viewer.init(document.getElementById('junction-viewport-3d'));
        await this.computeAndRenderAlignment();
    },

    attachEvents() {
        const slider = document.getElementById('n-bridge-slider');
        const valBadge = document.getElementById('n-bridge-val');
        if (slider) {
            slider.oninput = (e) => {
                this.state.nBridgeFrames = parseInt(e.target.value, 10);
                if (valBadge) valBadge.innerText = `${this.state.nBridgeFrames} Frames`;
            };
            slider.onchange = () => this.computeAndRenderAlignment();
        }

        document.querySelectorAll('.btn-mode').forEach(btn => {
            btn.onclick = () => {
                const mode = btn.dataset.mode;
                document.querySelectorAll('.btn-mode').forEach(b => {
                    b.classList.remove('btn-primary');
                    b.classList.add('btn-secondary');
                });
                btn.classList.remove('btn-secondary');
                btn.classList.add('btn-primary');
                this.setMode(mode);
            };
        });

        const v = this.state.viewer;
        const playBtn = document.getElementById('btn-trans-play');
        const prevBtn = document.getElementById('btn-trans-prev');
        const nextBtn = document.getElementById('btn-trans-next');
        const scrubber = document.getElementById('trans-scrubber');

        if (playBtn) playBtn.onclick = () => {
            const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
            if (j) v.toggleTransitionPlay(j);
        };
        if (prevBtn) prevBtn.onclick = () => {
            const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
            if (j) v.stepTransition(j, -1);
        };
        if (nextBtn) nextBtn.onclick = () => {
            const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
            if (j) v.stepTransition(j, 1);
        };
        if (scrubber) scrubber.oninput = (e) => {
            const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
            if (j) {
                v.pauseTransition();
                v.renderTransitionFrame(j, parseInt(e.target.value, 10));
            }
        };
    },

    setMode(mode) {
        this.state.viewer.stopBlink();
        this.state.viewer.viewMode = mode;
        const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
        if (!j) return;

        if (mode === 'blink') {
            this.state.viewer.startBlink(j);
        } else {
            this.state.viewer.apply3DScene(j, mode);
        }
    },

    async computeAndRenderAlignment() {
        const banner = document.getElementById('stitch-status-banner');
        if (banner) {
            banner.innerText = 'Calculating Kabsch alignment and bridge interpolation...';
            banner.style.display = 'block';
        }

        try {
            const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'preview_alignment',
                    n_bridge_frames: this.state.nBridgeFrames,
                    junction_choices: this.state.junctionChoices
                })
            });

            const data = await res.json();
            if (!data.success) {
                if (banner) banner.innerText = data.error || 'Failed to align steps';
                return;
            }

            this.state.megajson = data.megajson;
            this.state.junctionPreviews = data.junction_previews || [];
            if (banner) banner.style.display = 'none';

            this.renderJunctions();
            updateGlobalMetrics(data.megajson);

            if (this.state.junctionPreviews.length > 0) {
                this.selectJunction(this.state.activeJunctionIdx, true);
            }
        } catch (err) {
            if (banner) banner.innerText = 'Error: ' + err.message;
        }
    },

    renderJunctions() {
        const container = document.getElementById('junctions-list');
        const badge = document.getElementById('junctions-count-badge');
        renderJunctionsQueue(container, badge, this.state, {
            onSelect: (idx) => this.selectJunction(idx),
            onStrategyChange: (idx, val) => {
                this.state.junctionChoices[idx] = val;
                this.computeAndRenderAlignment();
            }
        });
    },

    selectJunction(idx, resetCamera = false) {
        this.state.activeJunctionIdx = idx;
        this.renderJunctions();

        const previews = this.state.junctionPreviews;
        if (!previews || idx >= previews.length) return;
        const j = previews[idx];

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

        const seq = this.state.viewer.getTransitionSequence(j);
        const scrubber = document.getElementById('trans-scrubber');
        if (scrubber) {
            scrubber.max = Math.max(0, seq.length - 1);
            if (this.state.viewer.transitionFrameIdx >= seq.length) {
                this.state.viewer.transitionFrameIdx = 0;
            }
            scrubber.value = this.state.viewer.transitionFrameIdx;
        }

        this.setMode(this.state.viewer.viewMode);
    },

    async getData(uiSchema) {
        return {
            n_bridge_frames: this.state.nBridgeFrames,
            junction_choices: this.state.junctionChoices
        };
    }
};
