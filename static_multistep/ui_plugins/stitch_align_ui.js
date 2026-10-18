import { Stitch3DViewer } from './stitch/stitch_viewer_3d.js?v=2.0.8';
import { buildStitchAlignHTML, renderJunctionsQueue, updateGlobalMetrics } from './stitch/stitch_junctions_dom.js?v=2.0.8';

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

        if (uiSchema.n_bridge_frames !== undefined) {
            this.state.nBridgeFrames = parseInt(uiSchema.n_bridge_frames, 10);
        } else if (globalState && globalState[stepId] && globalState[stepId].n_bridge_frames !== undefined) {
            this.state.nBridgeFrames = parseInt(globalState[stepId].n_bridge_frames, 10);
        }

        container.innerHTML = buildStitchAlignHTML(this.state);

        this.attachEvents();
        this.state.viewer.init(document.getElementById('junction-viewport-3d'));

        if (uiSchema.megajson && uiSchema.junction_previews && uiSchema.junction_previews.length > 0) {
            this.state.megajson = uiSchema.megajson;
            this.state.junctionPreviews = uiSchema.junction_previews;
            const banner = document.getElementById('stitch-status-banner');
            if (banner) banner.style.display = 'none';
            this.renderJunctions();
            updateGlobalMetrics(uiSchema.megajson);
            this.selectJunction(0, true);
        } else if (uiSchema.megajson) {
            this.state.megajson = uiSchema.megajson;
            this.state.junctionPreviews = uiSchema.junction_previews || [];
            this.renderJunctions();
            updateGlobalMetrics(uiSchema.megajson);
            this.renderFullTrajectoryView(uiSchema.megajson);
        } else {
            await this.computeAndRenderAlignment();
        }
    },

    attachEvents() {
        const recalcBtn = document.getElementById('btn-recalculate-align');
        if (recalcBtn) {
            recalcBtn.onclick = async () => {
                recalcBtn.innerHTML = '⚡ Stitching...';
                recalcBtn.disabled = true;
                await this.computeAndRenderAlignment();
                recalcBtn.innerHTML = '⚡ Re-Stitch Pathway';
                recalcBtn.disabled = false;
            };
        }

        const applyBridgeBtn = document.getElementById('btn-apply-bridge');
        if (applyBridgeBtn) {
            applyBridgeBtn.onclick = async () => {
                applyBridgeBtn.innerHTML = '⚡...';
                applyBridgeBtn.disabled = true;
                await this.computeAndRenderAlignment();
                applyBridgeBtn.innerHTML = '⚡ Re-Stitch';
                applyBridgeBtn.disabled = false;
            };
        }

        const slider = document.getElementById('n-bridge-slider');
        const valBadge = document.getElementById('n-bridge-val');
        let debounceTimer = null;
        if (slider) {
            slider.oninput = (e) => {
                this.state.nBridgeFrames = parseInt(e.target.value, 10);
                if (valBadge) valBadge.innerText = `${this.state.nBridgeFrames} Frames`;
                const banner = document.getElementById('stitch-status-banner');
                if (banner) {
                    banner.innerText = `Re-stitching with ${this.state.nBridgeFrames} bridge frames...`;
                    banner.style.display = 'block';
                }
                clearTimeout(debounceTimer);
                debounceTimer = setTimeout(() => {
                    this.computeAndRenderAlignment();
                }, 250);
            };
            slider.onchange = () => {
                clearTimeout(debounceTimer);
                this.computeAndRenderAlignment();
            };
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
            if (v.viewMode === 'full_trajectory' && this.state.megajson) {
                v.toggleTrajectoryPlay(this.state.megajson);
            } else {
                const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
                if (j) v.toggleTransitionPlay(j);
            }
        };
        if (prevBtn) prevBtn.onclick = () => {
            if (v.viewMode === 'full_trajectory' && this.state.megajson) {
                v.stepTrajectory(this.state.megajson, -1);
            } else {
                const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
                if (j) v.stepTransition(j, -1);
            }
        };
        if (nextBtn) nextBtn.onclick = () => {
            if (v.viewMode === 'full_trajectory' && this.state.megajson) {
                v.stepTrajectory(this.state.megajson, 1);
            } else {
                const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
                if (j) v.stepTransition(j, 1);
            }
        };
        if (scrubber) scrubber.oninput = (e) => {
            const val = parseInt(e.target.value, 10);
            if (v.viewMode === 'full_trajectory' && this.state.megajson) {
                v.pauseTransition();
                v.renderTrajectoryFrame(this.state.megajson, val);
            } else {
                const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
                if (j) {
                    v.pauseTransition();
                    v.renderTransitionFrame(j, val);
                }
            }
        };
    },

    setMode(mode) {
        this.state.viewer.stopBlink();
        this.state.viewer.viewMode = mode;
        if (mode === 'full_trajectory') {
            this.renderFullTrajectoryView(this.state.megajson);
            return;
        }

        const j = this.state.junctionPreviews[this.state.activeJunctionIdx];
        if (!j) return;

        if (mode === 'blink') {
            this.state.viewer.startBlink(j);
        } else {
            this.state.viewer.apply3DScene(j, mode, false);
        }
    },

    renderFullTrajectoryView(megajson) {
        if (!megajson || !megajson.frames) return;
        const v = this.state.viewer;
        v.stopBlink();
        v.viewMode = 'full_trajectory';
        const frames = megajson.frames || [];
        const scrubber = document.getElementById('trans-scrubber');
        if (scrubber) {
            scrubber.max = Math.max(0, frames.length - 1);
            scrubber.value = 0;
        }
        const transControls = document.getElementById('transition-controls-bar');
        if (transControls) transControls.style.display = 'flex';
        v.renderTrajectoryFrame(megajson, 0, true);
    },

    async computeAndRenderAlignment() {
        const banner = document.getElementById('stitch-status-banner');
        if (banner) {
            banner.innerText = `Calculating alignment & generating ${this.state.nBridgeFrames} bridge frames...`;
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

            if (this.state.viewer && this.state.viewer.viewMode === 'full_trajectory') {
                this.renderFullTrajectoryView(data.megajson);
            } else if (this.state.junctionPreviews.length > 0) {
                this.selectJunction(this.state.activeJunctionIdx, false);
            } else if (data.megajson) {
                this.renderFullTrajectoryView(data.megajson);
            }
        } catch (err) {
            if (banner) banner.innerText = 'Error: ' + err.message;
        }
    },

    renderJunctions() {
        const container = document.getElementById('junctions-list');
        const badge = document.getElementById('junctions-count-badge');
        renderJunctionsQueue(container, badge, this.state, {
            onSelect: (idx) => this.selectJunction(idx, true),
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
        if (pInfo) pInfo.innerText = `${j.step_prev_name} (${(j.prev_symbols || []).length} atoms)`;
        if (nInfo) nInfo.innerText = `${j.step_next_name} (${(j.next_symbols || []).length} atoms)`;

        const seq = this.state.viewer.getTransitionSequence(j);
        const scrubber = document.getElementById('trans-scrubber');
        if (scrubber) {
            scrubber.max = Math.max(0, seq.length - 1);
            if (this.state.viewer.transitionFrameIdx >= seq.length) {
                this.state.viewer.transitionFrameIdx = 0;
            }
            scrubber.value = this.state.viewer.transitionFrameIdx;
        }

        this.setMode(this.state.viewer.viewMode || 'transition');
    },

    async getData(uiSchema) {
        return {
            n_bridge_frames: this.state.nBridgeFrames,
            junction_choices: this.state.junctionChoices
        };
    }
};
