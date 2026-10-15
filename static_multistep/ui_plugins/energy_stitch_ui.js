import { drawPESProfile } from './stitch/energy_pes_canvas.js';
import { renderFragmentInstrumentHTML } from './stitch/energy_fragment_dom.js';
import { Energy3DPlayer } from './stitch/energy_player_3d.js';

export const EnergyStitchUI = {
    state: {
        stepId: 'step_3_energy_stitch',
        megajson: null,
        energyProfile: null,
        stoichiometry: null,
        boundaries: [],
        player: null,
        solvent: 'water',
        fragmentRoles: {}
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.megajson = uiSchema.megajson;
        this.state.energyProfile = uiSchema.energy_profile || {};
        this.state.stoichiometry = uiSchema.stoichiometry_summary || {};
        this.state.boundaries = uiSchema.step_boundaries || [];
        this.state.solvent = uiSchema.solvent || 'water';
        this.state.fragmentRoles = Object.assign({}, uiSchema.fragment_roles || {});
        this.state.player = new Energy3DPlayer();
        this.state.player.megajson = this.state.megajson;

        const ep = this.state.energyProfile;
        const stoich = this.state.stoichiometry;
        const initReacs = (ep.initial_reactants || []).join(' + ') || 'Reactants';
        const finalProds = (ep.final_products || []).join(' + ') || 'Products';

        container.innerHTML = `
            <div class="choose-container">
                <!-- Left Panel: Interactive Fragment Instrument & Stoichiometry -->
                <div class="left-panel">
                    <div class="queue-card" style="flex: 0 0 auto;">
                        <div class="section-title">
                            <span>General Reaction Stoichiometry</span>
                            <span class="chip emerald">Free Energy (ΔG) Balanced</span>
                        </div>
                        <div style="background: rgba(16, 185, 129, 0.08); border: 1px solid rgba(16, 185, 129, 0.2); border-radius: var(--radius-sm); padding: 0.6rem 0.8rem; font-size: 0.8rem;">
                            <div style="color: var(--text-muted); font-size: 0.7rem; text-transform: uppercase; font-weight: 600;">Overall Reaction Equation:</div>
                            <div style="font-family: var(--font-mono); font-size: 0.85rem; font-weight: 700; color: #6ee7b7; margin-top: 0.2rem;">
                                ${initReacs} ➔ ${finalProds}
                            </div>
                        </div>
                    </div>

                    <div class="catalog-card" style="flex: 1; max-height: 380px; overflow-y: auto;">
                        <div class="section-title">
                            <span>Fragment Exclusion & Roles</span>
                            <span class="chip purple" id="active-solvent-badge">${this.state.solvent.toUpperCase()}</span>
                        </div>
                        <p style="font-size: 0.725rem; color: var(--text-muted); margin-bottom: 0.6rem;">
                            Auto-detected from trajectory graphs. Manually adjust which fragments are <strong>Hero Core</strong>, <strong>Incoming Reagents</strong>, or <strong>Excluded Leaving Groups</strong>.
                        </p>
                        <div id="fragment-instrument-list" style="display: flex; flex-direction: column; gap: 0.6rem;">
                            ${renderFragmentInstrumentHTML(stoich, this.state.fragmentRoles)}
                        </div>
                    </div>

                    <div class="catalog-card" style="flex: 0 0 auto;">
                        <div class="section-title">Solvation & Recalculation</div>
                        <div style="display: flex; flex-direction: column; gap: 0.5rem; margin-bottom: 0.75rem;">
                            <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.8rem;">
                                <label for="select-solvent" style="color: var(--text-muted); font-weight: 600;">Solvent Model:</label>
                                <select id="select-solvent" style="background: var(--bg-surface); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; padding: 0.25rem 0.5rem; font-size: 0.75rem;">
                                    <option value="water" ${this.state.solvent === 'water' ? 'selected' : ''}>💧 Water (ALPB implicit)</option>
                                    <option value="none" ${this.state.solvent === 'none' ? 'selected' : ''}>🌫️ Gas Phase / Vacuum</option>
                                    <option value="chcl3" ${this.state.solvent === 'chcl3' ? 'selected' : ''}>🧪 Chloroform</option>
                                    <option value="toluene" ${this.state.solvent === 'toluene' ? 'selected' : ''}>🛢️ Toluene</option>
                                </select>
                            </div>
                            <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; font-size: 0.8rem;">
                                <div class="detail-box">
                                    <div class="detail-box-title">Activation Barrier (ΔG‡)</div>
                                    <div id="pes-barrier-val" style="font-weight: 700; color: #fcd34d; font-family: var(--font-mono);">${(ep.overall_barrier_kcal || 0).toFixed(2)} kcal/mol</div>
                                </div>
                                <div class="detail-box">
                                    <div class="detail-box-title">Net Reaction (ΔG°)</div>
                                    <div id="pes-rxn-val" style="font-weight: 700; color: #6ee7b7; font-family: var(--font-mono);">${(ep.overall_reaction_energy_kcal || 0).toFixed(2)} kcal/mol</div>
                                </div>
                            </div>
                        </div>

                        <button id="btn-recalc-pes" class="btn btn-sm btn-primary" style="width: 100%;">
                            ⚡ Apply Roles & Recalculate Free Energy
                        </button>
                    </div>
                </div>

                <!-- Right Panel: Synchronized 2D PES & 3D Preview -->
                <div class="preview-panel">
                    <div class="preview-header">
                        <div>
                            <div class="preview-rxn-title">Hess's Law Free Energy Surface (ΔG at 298.15 K)</div>
                            <div class="preview-meta-chips">
                                <span class="chip blue">${this.state.megajson ? this.state.megajson.frames.length : 0} Stitched Frames</span>
                                <span class="chip emerald">ALPB Thermochemistry</span>
                                <span class="chip amber">${this.state.boundaries.length} IRC Boundaries</span>
                            </div>
                        </div>
                    </div>

                    <div style="background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.75rem; margin-bottom: 0.75rem; position: relative;">
                        <canvas id="pes-canvas-chart" style="width: 100%; height: 160px; display: block; cursor: crosshair;"></canvas>
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 0.4rem; font-size: 0.7rem; color: var(--text-muted); font-family: var(--font-mono);">
                            <span>Reaction Coordinate (Step 1 ➔ Step ${this.state.boundaries.length})</span>
                            <span id="pes-hover-tooltip" style="color: var(--primary);">Click or drag on chart to scrub timeline</span>
                        </div>
                    </div>

                    <div class="viewer-wrapper" style="height: 280px;">
                        <div id="energy-viewport-3d" class="viewport-3d"></div>
                        <div id="energy-status-badge" class="viewer-overlay-badge">Loading Stitched Structure...</div>
                    </div>

                    <div class="playback-bar">
                        <button id="btn-pes-play" class="btn btn-sm btn-secondary">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                                <polygon points="5 3 19 12 5 21 5 3"/>
                            </svg>
                        </button>
                        <button id="btn-pes-prev" class="btn btn-sm btn-secondary">⏮</button>
                        <input type="range" id="pes-scrubber" class="scrub-slider" min="0" max="${this.state.megajson ? this.state.megajson.frames.length - 1 : 0}" value="0">
                        <button id="btn-pes-next" class="btn btn-sm btn-secondary">⏭</button>
                        <div id="pes-counter" class="frame-counter">0 / 0</div>
                    </div>
                </div>
            </div>
        `;

        this.attachEvents();
        this.state.player.init(document.getElementById('energy-viewport-3d'));
        this.state.player.setFrame(0, true, () => this.drawPESChart());
        this.drawPESChart();
    },

    attachEvents() {
        const selectSolvent = document.getElementById('select-solvent');
        if (selectSolvent) {
            selectSolvent.onchange = (e) => {
                this.state.solvent = e.target.value;
                const badge = document.getElementById('active-solvent-badge');
                if (badge) badge.innerText = this.state.solvent.toUpperCase();
            };
        }

        const recalcBtn = document.getElementById('btn-recalc-pes');
        if (recalcBtn) {
            recalcBtn.onclick = () => this.recalculatePES();
        }

        const container = document.getElementById('fragment-instrument-list');
        if (container) {
            container.onchange = (e) => {
                if (e.target.classList.contains('select-frag-role')) {
                    const fKey = e.target.dataset.fragKey;
                    this.state.fragmentRoles[fKey] = e.target.value;
                }
            };
        }

        const p = this.state.player;
        const playBtn = document.getElementById('btn-pes-play');
        const prevBtn = document.getElementById('btn-pes-prev');
        const nextBtn = document.getElementById('btn-pes-next');
        const scrubber = document.getElementById('pes-scrubber');

        if (playBtn) playBtn.onclick = () => p.togglePlay(() => this.drawPESChart());
        if (prevBtn) prevBtn.onclick = () => p.stepFrame(-1, () => this.drawPESChart());
        if (nextBtn) nextBtn.onclick = () => p.stepFrame(1, () => this.drawPESChart());
        if (scrubber) {
            scrubber.oninput = (e) => {
                p.pause();
                p.setFrame(parseInt(e.target.value, 10), false, () => this.drawPESChart());
            };
        }

        const canvas = document.getElementById('pes-canvas-chart');
        if (canvas) {
            let isDragging = false;
            const handleCanvasPos = (e) => {
                const rect = canvas.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const padLeft = 45;
                const padRight = 15;
                const plotW = rect.width - padLeft - padRight;
                const frac = Math.max(0, Math.min(1, (x - padLeft) / plotW));
                const totalFrames = this.state.megajson ? this.state.megajson.frames.length : 1;
                const targetIdx = Math.round(frac * (totalFrames - 1));
                p.pause();
                p.setFrame(targetIdx, false, () => this.drawPESChart());
            };

            canvas.onmousedown = (e) => { isDragging = true; handleCanvasPos(e); };
            window.addEventListener('mousemove', (e) => { if (isDragging) handleCanvasPos(e); });
            window.addEventListener('mouseup', () => { isDragging = false; });
        }
    },

    drawPESChart() {
        const canvas = document.getElementById('pes-canvas-chart');
        const activeF = this.state.player ? this.state.player.activeFrame : 0;
        drawPESProfile(canvas, this.state.energyProfile, this.state.boundaries, activeF);
    },

    async recalculatePES() {
        const btn = document.getElementById('btn-recalc-pes');
        if (btn) {
            btn.disabled = true;
            btn.innerText = '⏳ Recalculating Hess Free Energies...';
        }

        try {
            const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'recalculate',
                    solvent: this.state.solvent,
                    fragment_roles: this.state.fragmentRoles
                })
            });

            const data = await res.json();
            if (data.success && data.energy_profile) {
                this.state.energyProfile = data.energy_profile;
                this.state.megajson = data.megajson;
                this.state.player.megajson = data.megajson;
                const bVal = document.getElementById('pes-barrier-val');
                const rVal = document.getElementById('pes-rxn-val');
                if (bVal) bVal.innerText = `${data.energy_profile.overall_barrier_kcal.toFixed(2)} kcal/mol`;
                if (rVal) rVal.innerText = `${data.energy_profile.overall_reaction_energy_kcal.toFixed(2)} kcal/mol`;
                this.drawPESChart();
                this.state.player.setFrame(this.state.player.activeFrame, false, () => this.drawPESChart());
            } else {
                alert('Recalculation error: ' + (data.error || 'Unknown'));
            }
        } catch (err) {
            alert('Server error during recalculation: ' + err.message);
        } finally {
            if (btn) {
                btn.disabled = false;
                btn.innerText = '⚡ Apply Roles & Recalculate Free Energy';
            }
        }
    },

    async getData(uiSchema) {
        return {
            solvent: this.state.solvent,
            fragment_roles: this.state.fragmentRoles
        };
    }
};
