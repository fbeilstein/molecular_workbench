export const EnergyStitchUI = {
    state: {
        stepId: 'step_3_energy_stitch',
        megajson: null,
        energyProfile: null,
        stoichiometry: null,
        boundaries: [],
        activeFrame: 0,
        isPlaying: false,
        playTimer: null,
        viewer: null,
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
        this.state.activeFrame = 0;
        this.pause();

        const ep = this.state.energyProfile;
        const stoich = this.state.stoichiometry;
        const initReacs = (ep.initial_reactants || []).join(' + ') || 'Reactants';
        const finalProds = (ep.final_products || []).join(' + ') || 'Products';

        container.innerHTML = `
            <div class="choose-container">
                <!-- Left Panel: Interactive Fragment Instrument & Stoichiometry -->
                <div class="left-panel">
                    <!-- General Stoichiometry Banner -->
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

                    <!-- Fragment Role & Exclusion Instrument -->
                    <div class="catalog-card" style="flex: 1; max-height: 380px; overflow-y: auto;">
                        <div class="section-title">
                            <span>Fragment Exclusion & Roles</span>
                            <span class="chip purple" id="active-solvent-badge">${this.state.solvent.toUpperCase()}</span>
                        </div>
                        <p style="font-size: 0.725rem; color: var(--text-muted); margin-bottom: 0.6rem;">
                            Auto-detected from trajectory graphs. Manually adjust which fragments are <strong>Hero Core</strong>, <strong>Incoming Reagents</strong>, or <strong>Excluded Leaving Groups</strong>.
                        </p>
                        <div id="fragment-instrument-list" style="display: flex; flex-direction: column; gap: 0.6rem;">
                            ${this.renderFragmentInstrumentHTML(stoich)}
                        </div>
                    </div>

                    <!-- Environment & Calculation Controls Card -->
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
                                <span class="chip blue">${(this.state.megajson.frames || []).length} Total Frames</span>
                                <span class="chip emerald">Net ΔG°: ${(ep.overall_reaction_energy_kcal || 0).toFixed(2)} kcal/mol</span>
                                <span class="chip amber">Barrier ΔG‡: ${(ep.overall_barrier_kcal || 0).toFixed(2)} kcal/mol</span>
                            </div>
                        </div>
                        <div style="font-size: 0.75rem; color: var(--text-muted); font-family: var(--font-mono);">
                            Frame <span id="active-frame-num" style="color: #93c5fd; font-weight: 700;">1</span>: <span id="active-frame-energy" style="color: #6ee7b7; font-weight: 700;">0.00</span> kcal/mol
                        </div>
                    </div>

                    <!-- 2D Interactive PES Canvas Chart -->
                    <div style="background: #080c14; border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.75rem; position: relative;">
                        <canvas id="pes-canvas" width="760" height="200" style="width: 100%; height: 200px; display: block; cursor: crosshair;"></canvas>
                    </div>

                    <!-- 3D Preview Viewport Synchronized with PES -->
                    <div class="viewer-wrapper" style="height: 300px;">
                        <div id="energy-viewport-3d" class="viewport-3d"></div>
                        <div id="energy-status-badge" class="viewer-overlay-badge">Frame 1 | Energy: 0.00 kcal/mol</div>
                    </div>

                    <!-- Playback Controls -->
                    <div class="playback-bar">
                        <button id="btn-pes-play" class="btn btn-sm btn-secondary" title="Play / Pause">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                                <polygon points="5 3 19 12 5 21 5 3"/>
                            </svg>
                        </button>
                        <button id="btn-pes-prev" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                        <input type="range" id="pes-scrubber" class="scrub-slider" min="0" max="${Math.max(0, (this.state.megajson.frames || []).length - 1)}" value="0">
                        <button id="btn-pes-next" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                        <div id="pes-counter" class="frame-counter">1 / ${(this.state.megajson.frames || []).length}</div>
                    </div>
                </div>
            </div>
        `;

        this.initEventListeners();
        this.init3DViewer();
        this.drawPESChart();
        this.setFrame(0, true);
    },

    renderFragmentInstrumentHTML(stoich) {
        if (!stoich || !stoich.steps) {
            return `<div style="font-size: 0.75rem; color: var(--text-dim);">No fragment data found.</div>`;
        }

        let html = '';
        stoich.steps.forEach(s => {
            html += `<div style="font-size: 0.75rem; font-weight: 700; color: #93c5fd; margin-top: 0.25rem;">Step ${s.step_index} Species</div>`;

            // Reactants
            s.reactants.forEach(r => {
                const currentRole = this.state.fragmentRoles[r.id] || r.role || 'core';
                html += `
                    <div style="background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.5rem; display: flex; justify-content: space-between; align-items: center; gap: 0.5rem;">
                        <div>
                            <div style="font-size: 0.775rem; font-weight: 600; color: #cbd5e1;">${r.formula || r.smiles} <span style="font-size: 0.675rem; color: var(--text-muted);">(Reactant, ${r.num_atoms} atoms)</span></div>
                            <div style="font-size: 0.675rem; color: var(--text-dim); font-family: var(--font-mono);">${r.smiles} [q=${r.charge}]</div>
                        </div>
                        <select class="frag-role-select" data-frag-id="${r.id}" style="padding: 0.25rem 0.4rem; background: var(--bg-body); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; font-size: 0.725rem;">
                            <option value="core" ${currentRole === 'core' ? 'selected' : ''}>🌟 Core Backbone</option>
                            <option value="incoming" ${currentRole === 'incoming' ? 'selected' : ''}>📥 Incoming Reagent</option>
                            <option value="leaving" ${currentRole === 'leaving' ? 'selected' : ''}>📤 Leaving Group (EXCLUDED)</option>
                            <option value="spectator" ${currentRole === 'spectator' ? 'selected' : ''}>🚫 Spectator (EXCLUDED)</option>
                        </select>
                    </div>
                `;
            });

            // Products
            s.products.forEach(p => {
                const currentRole = this.state.fragmentRoles[p.id] || p.role || (p.num_atoms <= 2 ? 'leaving' : 'core');
                html += `
                    <div style="background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: var(--radius-sm); padding: 0.5rem; display: flex; justify-content: space-between; align-items: center; gap: 0.5rem;">
                        <div>
                            <div style="font-size: 0.775rem; font-weight: 600; color: #cbd5e1;">${p.formula || p.smiles} <span style="font-size: 0.675rem; color: var(--text-muted);">(Product, ${p.num_atoms} atoms)</span></div>
                            <div style="font-size: 0.675rem; color: var(--text-dim); font-family: var(--font-mono);">${p.smiles} [q=${p.charge}]</div>
                        </div>
                        <select class="frag-role-select" data-frag-id="${p.id}" style="padding: 0.25rem 0.4rem; background: var(--bg-body); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; font-size: 0.725rem;">
                            <option value="core" ${currentRole === 'core' ? 'selected' : ''}>🌟 Core Product</option>
                            <option value="leaving" ${currentRole === 'leaving' ? 'selected' : ''}>📤 Leaving Group (EXCLUDED)</option>
                            <option value="incoming" ${currentRole === 'incoming' ? 'selected' : ''}>📥 Incoming Reagent</option>
                            <option value="spectator" ${currentRole === 'spectator' ? 'selected' : ''}>🚫 Spectator (EXCLUDED)</option>
                        </select>
                    </div>
                `;
            });
        });

        return html;
    },

    initEventListeners() {
        const playBtn = document.getElementById('btn-pes-play');
        const prevBtn = document.getElementById('btn-pes-prev');
        const nextBtn = document.getElementById('btn-pes-next');
        const scrubber = document.getElementById('pes-scrubber');
        const recalcBtn = document.getElementById('btn-recalc-pes');
        const solventSelect = document.getElementById('select-solvent');
        const canvas = document.getElementById('pes-canvas');

        if (playBtn) playBtn.onclick = () => this.togglePlay();
        if (prevBtn) prevBtn.onclick = () => this.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepFrame(1);
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pause();
                this.setFrame(parseInt(e.target.value, 10));
            };
        }

        if (solventSelect) {
            solventSelect.onchange = (e) => {
                this.state.solvent = e.target.value;
                const badge = document.getElementById('active-solvent-badge');
                if (badge) badge.innerText = e.target.value.toUpperCase();
            };
        }

        // Fragment role change handlers
        document.querySelectorAll('.frag-role-select').forEach(sel => {
            sel.onchange = (e) => {
                const fid = sel.dataset.fragId;
                this.state.fragmentRoles[fid] = sel.value;
            };
        });

        if (canvas) {
            const handleCanvasPos = (e) => {
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
                handleCanvasPos(e);
                const onMove = (ev) => handleCanvasPos(ev);
                const onUp = () => {
                    window.removeEventListener('mousemove', onMove);
                    window.removeEventListener('mouseup', onUp);
                };
                window.addEventListener('mousemove', onMove);
                window.addEventListener('mouseup', onUp);
            };
        }

        if (recalcBtn) {
            recalcBtn.onclick = async () => {
                recalcBtn.disabled = true;
                recalcBtn.innerText = 'Calculating xTB in ' + this.state.solvent.toUpperCase() + '...';
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({
                            action_type: 'recalculate_energies',
                            recalculate_with_xtb: true,
                            solvent: this.state.solvent,
                            fragment_roles: this.state.fragmentRoles
                        })
                    });
                    const data = await res.json();
                    if (data.success && data.megajson) {
                        this.state.megajson = data.megajson;
                        this.state.energyProfile = data.energy_profile;
                        this.state.stoichiometry = data.stoichiometry_summary;

                        const bVal = document.getElementById('pes-barrier-val');
                        const rVal = document.getElementById('pes-rxn-val');
                        if (bVal) bVal.innerText = `${(data.energy_profile.overall_barrier_kcal || 0).toFixed(2)} kcal/mol`;
                        if (rVal) rVal.innerText = `${(data.energy_profile.overall_reaction_energy_kcal || 0).toFixed(2)} kcal/mol`;

                        this.drawPESChart();
                        this.setFrame(this.state.activeFrame);
                    }
                } catch (e) {
                    console.error('Recalc error:', e);
                } finally {
                    recalcBtn.disabled = false;
                    recalcBtn.innerText = '⚡ Apply Roles & Recalculate PES';
                }
            };
        }
    },

    init3DViewer() {
        const container = document.getElementById('energy-viewport-3d');
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
        const canvas = document.getElementById('pes-canvas');
        if (!canvas) return;

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
        const padTop = 20;
        const padBottom = 30;
        const plotW = w - padLeft - padRight;
        const plotH = h - padTop - padBottom;

        const ep = this.state.energyProfile || {};
        const energies = ep.frame_energies_kcal || [];
        if (energies.length === 0) return;

        let minE = Math.min(...energies);
        let maxE = Math.max(...energies);
        if (minE === maxE) { minE -= 5; maxE += 5; }
        const margin = (maxE - minE) * 0.15 || 2;
        minE -= margin;
        maxE += margin;

        const toX = (idx) => padLeft + (idx / (energies.length - 1)) * plotW;
        const toY = (val) => padTop + plotH - ((val - minE) / (maxE - minE)) * plotH;

        ctx.clearRect(0, 0, w, h);

        // Grid lines
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.08)';
        ctx.lineWidth = 1;
        ctx.beginPath();
        for (let i = 0; i <= 4; i++) {
            const yVal = minE + (i / 4) * (maxE - minE);
            const yPos = toY(yVal);
            ctx.moveTo(padLeft, yPos);
            ctx.lineTo(w - padRight, yPos);
        }
        ctx.stroke();

        // Y-axis labels
        ctx.fillStyle = '#94a3b8';
        ctx.font = '10px monospace';
        ctx.textAlign = 'right';
        for (let i = 0; i <= 4; i++) {
            const yVal = minE + (i / 4) * (maxE - minE);
            ctx.fillText(yVal.toFixed(1), padLeft - 6, toY(yVal) + 3);
        }

        // Zero reference line
        if (minE <= 0 && maxE >= 0) {
            const zeroY = toY(0);
            ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)';
            ctx.setLineDash([4, 4]);
            ctx.beginPath();
            ctx.moveTo(padLeft, zeroY);
            ctx.lineTo(w - padRight, zeroY);
            ctx.stroke();
            ctx.setLineDash([]);
        }

        // Shading
        const boundaries = this.state.boundaries || [];
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
        grad.addColorStop(0, 'rgba(56, 189, 248, 0.25)');
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
                ctx.arc(x, y, 4.5, 0, 2 * Math.PI);
                ctx.fill();
                ctx.fillStyle = '#fde68a';
                ctx.font = 'bold 9px monospace';
                ctx.textAlign = 'center';
                ctx.fillText(`TS${b.step_index}`, x, y - 8);
            }
        });

        // Active frame line
        const activeX = toX(this.state.activeFrame);
        const activeY = toY(energies[this.state.activeFrame]);
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(activeX, padTop);
        ctx.lineTo(activeX, padTop + plotH);
        ctx.stroke();

        // Active point
        ctx.fillStyle = '#ef4444';
        ctx.beginPath();
        ctx.arc(activeX, activeY, 5, 0, 2 * Math.PI);
        ctx.fill();
    },

    setFrame(index, resetCamera = false) {
        if (!this.state.megajson || !this.state.viewer) return;
        const frames = this.state.megajson.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.state.activeFrame = index;
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

        const scrubber = document.getElementById('pes-scrubber');
        const counter = document.getElementById('pes-counter');
        const badge = document.getElementById('energy-status-badge');
        const numLabel = document.getElementById('active-frame-num');
        const energyLabel = document.getElementById('active-frame-energy');

        const eKcal = frame.energy_kcal !== undefined ? frame.energy_kcal : 0.0;
        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;
        if (numLabel) numLabel.innerText = `${index + 1}`;
        if (energyLabel) energyLabel.innerText = `${eKcal.toFixed(2)}`;

        if (badge) {
            const stepIdx = frame.step_index !== undefined ? frame.step_index : '?';
            const isTS = frame.is_ts ? ' ⭐ [TS]' : '';
            badge.innerText = `Step ${stepIdx} | Frame ${index + 1} / ${frames.length} | Energy: ${eKcal.toFixed(2)} kcal/mol ${isTS}`;
        }

        this.drawPESChart();
    },

    stepFrame(delta) {
        this.pause();
        const frames = this.state.megajson ? this.state.megajson.frames : [];
        if (!frames || frames.length === 0) return;
        let target = this.state.activeFrame + delta;
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

        if (this.state.activeFrame >= frames.length - 1) {
            this.setFrame(0);
        }

        this.state.isPlaying = true;
        const btn = document.getElementById('btn-pes-play');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.state.playTimer);
        this.state.playTimer = setInterval(() => {
            const next = this.state.activeFrame + 1;
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
        const btn = document.getElementById('btn-pes-play');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    },

    async getData(uiSchema) {
        this.pause();
        return {
            confirmed: true,
            recalculate_with_xtb: true,
            solvent: this.state.solvent,
            fragment_roles: this.state.fragmentRoles
        };
    }
};
