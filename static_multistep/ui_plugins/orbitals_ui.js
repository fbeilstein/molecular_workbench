export const OrbitalsUI = {
    state: {
        stepId: 'step_4_orbitals',
        uiSchema: null,
        totalFrames: 0,
        currentFrame: 0,
        landmarks: [],
        energyProfile: [],
        symbols: [],
        timeline: {
            default_isovalue: 0.03,
            tracks: []
        },
        // In/Out punching
        inFrame: null,
        outFrame: null,
        // Auditioning orbital
        auditionOrbital: null,
        auditionLabel: '',
        // Active 3D rendered orbitals
        activeOrbitals: new Set(),
        // Manual override set during lecture (overrides track default_on)
        manualOverrides: new Map(), // key -> boolean
        // Charge Distribution overlay mode
        showCharges: false,
        // Isovalue
        isovalue: 0.03,
        // Playback
        isPlaying: false,
        playTimer: null,
        // Caches & Viewer
        viewer: null,
        cubeCache: {},
        frameMetaCache: {},
        densityFlux: null,
        isPrecomputing: false,
        precomputeAbort: false,
        activeCategory: 'all',
        searchQuery: '',
        testedOrbitals: new Set(),
        sortMode: 'chemical',
        reviewFilter: 'all'
    },

    colorPalette: [
        '#06b6d4', '#a855f7', '#10b981', '#f59e0b', '#ec4899',
        '#3b82f6', '#84cc16', '#6366f1', '#14b8a6', '#f97316'
    ],

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.uiSchema = uiSchema;
        this.state.totalFrames = uiSchema.total_frames || 0;
        this.state.landmarks = uiSchema.landmarks || [];
        this.state.energyProfile = uiSchema.energy_profile || [];
        this.state.symbols = uiSchema.symbols || [];
        this.state.timeline = JSON.parse(JSON.stringify(uiSchema.timeline || { default_isovalue: 0.03, tracks: [] }));
        this.state.isovalue = this.state.timeline.default_isovalue || 0.03;

        try {
            const saved = localStorage.getItem('mw_tested_orbitals');
            if (saved) {
                this.state.testedOrbitals = new Set(JSON.parse(saved));
            }
        } catch (e) {}

        if (!uiSchema.has_megajson || this.state.totalFrames === 0) {
            container.innerHTML = `<div class="empty-state">No reaction pathway available. Please complete earlier steps first.</div>`;
            return;
        }

        container.innerHTML = `
            <div class="kdenlive-studio" style="display: flex; flex-direction: column; gap: 0.85rem; user-select: none;">
                <!-- Header Toolbar: Status, Precompute & Lecture Controls -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.75rem; background: rgba(15, 23, 42, 0.7); border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.65rem 1rem;">
                    <div>
                        <div style="display: flex; align-items: center; gap: 0.6rem;">
                            <span style="font-size: 1.1rem; font-weight: 700; color: var(--text-main);">Reaction Orbital Studio</span>
                            <span class="chip purple" style="font-size: 0.68rem; font-weight: 600;">Kdenlive Timeline</span>
                            <span class="chip blue" style="font-size: 0.68rem;">${this.state.totalFrames} frames</span>
                            <span class="chip emerald" style="font-size: 0.68rem;">${this.state.symbols.length} atoms</span>
                        </div>
                        <div style="font-size: 0.74rem; color: var(--text-muted); margin-top: 0.2rem;">
                            Audition any bonding or lone-pair orbital, punch In/Out clip ranges, and direct your lecture movie.
                        </div>
                    </div>
                    
                    <div style="display: flex; align-items: center; gap: 0.6rem; flex-wrap: wrap;">
                        <!-- Charge Distribution Mode Toggle -->
                        <button id="btn-toggle-charges" class="btn btn-sm btn-secondary" style="background: rgba(30, 41, 59, 0.9); font-weight: 600; font-size: 0.78rem; display: flex; align-items: center; gap: 0.4rem; border: 1px solid rgba(255,255,255,0.15);" title="Toggle atomic partial charges and electrostatic heat map (Hotkey: C)">
                            <span style="display: inline-block; width: 9px; height: 9px; border-radius: 50%; background: #64748b;" id="charge-toggle-indicator"></span>
                            <span>⚡ Partial Charges (ESP)</span>
                        </button>

                        <!-- Precompute Timeline Button -->
                        <button id="btn-precompute-timeline" class="btn btn-sm btn-primary" style="background: linear-gradient(135deg, #0284c7, #0369a1); font-size: 0.78rem; font-weight: 600; padding: 0.35rem 0.75rem;" title="Precalculate volumetric cubes for all timeline tracks for butter-smooth 60 FPS playback">
                            ▶ Precompute Timeline
                        </button>
                        <button id="btn-precompute-stop" class="btn btn-sm btn-secondary" style="display: none; font-size: 0.78rem;">⏸ Pause</button>

                        <!-- Cube Exporter -->
                        <button id="btn-dl-cube" class="btn btn-secondary btn-sm" title="Export Gaussian .cube for current active orbital">
                            <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>
                            Export .cube
                        </button>
                    </div>
                </div>

                <!-- Precompute Progress Bar -->
                <div id="timeline-progress-wrap" style="display: none; flex-direction: column; gap: 0.3rem; background: rgba(11, 15, 23, 0.9); border: 1px solid var(--border-color); border-radius: 6px; padding: 0.5rem 0.8rem;">
                    <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.74rem;">
                        <span id="timeline-progress-status" style="color: #38bdf8; font-weight: 600;">Precomputing timeline tracks...</span>
                        <span id="timeline-progress-pct" style="color: #6ee7b7; font-weight: 700; font-family: var(--font-mono);">0%</span>
                    </div>
                    <div style="width: 100%; height: 6px; background: rgba(0,0,0,0.5); border-radius: 3px; overflow: hidden;">
                        <div id="timeline-progress-fill" style="width: 0%; height: 100%; background: linear-gradient(90deg, #38bdf8, #10b981); transition: width 0.25s ease;"></div>
                    </div>
                </div>

                <!-- Main Middle Area: 3D Viewport on Left, Orbital Pool / Audition Inspector on Right -->
                <div style="display: grid; grid-template-columns: 1fr 390px; gap: 1rem; align-items: stretch;">
                    
                    <!-- 3D Viewport Column -->
                    <div style="display: flex; flex-direction: column; gap: 0.6rem;">
                        <div class="viewer-wrapper" style="height: 450px; position: relative; border-radius: var(--radius-md); overflow: hidden; border: 1px solid var(--border-color); background: #0b0f17;">
                            <div id="orb-viewport-3d" class="viewport-3d" style="width: 100%; height: 100%;"></div>
                            
                            <!-- Top Overlay: Current Frame, Energy & Landmark -->
                            <div style="position: absolute; top: 12px; left: 12px; display: flex; align-items: center; gap: 0.5rem; z-index: 5;">
                                <div id="orb-status-badge" class="viewer-overlay-badge" style="background: rgba(11,15,23,0.85); font-family: var(--font-mono); font-size: 0.75rem; font-weight: 700; border: 1px solid var(--border-color); padding: 0.35rem 0.65rem; border-radius: 6px; color: #f8fafc;">
                                    Frame 1 / ${this.state.totalFrames}
                                </div>
                                <div id="orb-energy-badge" class="chip blue" style="font-family: var(--font-mono); font-size: 0.72rem; font-weight: 600; padding: 0.3rem 0.55rem;">
                                    E: 0.00 kcal/mol
                                </div>
                                <div id="orb-landmark-chip" class="chip purple" style="display: none; font-size: 0.72rem; font-weight: 600;"></div>
                            </div>

                            <!-- Bottom Left: Active 3D Lobes Pill Badges with 1-click toggles and Hotkeys -->
                            <div id="active-lobes-panel" style="position: absolute; bottom: 12px; left: 12px; max-width: 80%; background: rgba(11, 15, 23, 0.88); backdrop-filter: blur(8px); border: 1px solid var(--border-color); border-radius: 6px; padding: 0.45rem 0.75rem; font-size: 0.75rem; z-index: 5; display: flex; flex-direction: column; gap: 0.3rem;">
                                <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.65rem; color: var(--text-dim); text-transform: uppercase; font-weight: 700;">
                                    <span>Active 3D Lobes (Click or keys 1..5 to toggle):</span>
                                </div>
                                <div id="active-lobes-chips" style="display: flex; gap: 0.4rem; flex-wrap: wrap; align-items: center;"></div>
                            </div>

                            <!-- Bottom Right: Charge Distribution Legend (When active) -->
                            <div id="charges-legend" style="display: none; position: absolute; bottom: 12px; right: 12px; background: rgba(11, 15, 23, 0.88); backdrop-filter: blur(8px); border: 1px solid var(--border-color); border-radius: 6px; padding: 0.4rem 0.65rem; font-size: 0.68rem; z-index: 5; flex-direction: column; gap: 0.25rem;">
                                <div style="font-weight: 700; color: #94a3b8; text-transform: uppercase;">Partial Charge (e⁻)</div>
                                <div style="display: flex; align-items: center; gap: 0.4rem; font-family: var(--font-mono);">
                                    <span style="color: #ef4444; font-weight: 700;">-1.0 (δ⁻)</span>
                                    <div style="width: 60px; height: 6px; border-radius: 3px; background: linear-gradient(90deg, #ef4444, #f8fafc, #3b82f6);"></div>
                                    <span style="color: #3b82f6; font-weight: 700;">+1.0 (δ⁺)</span>
                                </div>
                            </div>

                            <!-- Loading Spinner -->
                            <div id="orb-loading-spinner" style="display: none; position: absolute; inset: 0; background: rgba(11,15,23,0.65); backdrop-filter: blur(4px); justify-content: center; align-items: center; z-index: 10; color: #38bdf8; font-size: 0.85rem; font-weight: 600;">
                                <div style="display: flex; flex-direction: column; align-items: center; gap: 0.45rem;">
                                    <div class="spinner" style="width: 22px; height: 22px; border: 3px solid rgba(56,189,248,0.2); border-top-color: #38bdf8; border-radius: 50%; animation: spin 0.8s linear infinite;"></div>
                                    <span id="spinner-msg">Evaluating Volume...</span>
                                </div>
                            </div>
                        </div>

                        <!-- Transport & Mark In/Out Puncher Bar -->
                        <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.6rem; background: rgba(17, 24, 39, 0.7); border: 1px solid var(--border-color); border-radius: var(--radius-md); padding: 0.5rem 0.85rem;">
                            <!-- Transport Buttons -->
                            <div style="display: flex; align-items: center; gap: 0.35rem;">
                                <button id="btn-transport-start" class="btn btn-sm btn-secondary" title="Jump to First Frame (Home)">⏮</button>
                                <button id="btn-transport-prev" class="btn btn-sm btn-secondary" title="Previous Frame (Left Arrow)">⏴</button>
                                <button id="btn-transport-play" class="btn btn-sm btn-primary" style="padding: 0.3rem 0.8rem; font-weight: 700;" title="Play / Pause (Spacebar)">▶</button>
                                <button id="btn-transport-next" class="btn btn-sm btn-secondary" title="Next Frame (Right Arrow)">⏵</button>
                                <button id="btn-transport-end" class="btn btn-sm btn-secondary" title="Jump to Last Frame (End)">⏭</button>
                            </div>

                            <!-- In / Out Punch Controls -->
                            <div style="display: flex; align-items: center; gap: 0.4rem; background: rgba(0,0,0,0.35); padding: 0.25rem 0.5rem; border-radius: 6px; border: 1px solid var(--border-color);">
                                <button id="btn-mark-in" class="btn btn-sm btn-secondary" style="font-size: 0.75rem; padding: 0.25rem 0.55rem;" title="Start adding range at current frame (Hotkey: [)">
                                    [ Set In <span id="val-mark-in" style="font-family: var(--font-mono); color: #38bdf8; font-weight: 700; margin-left: 0.2rem;">--</span>
                                </button>
                                <button id="btn-mark-out" class="btn btn-sm btn-secondary" style="font-size: 0.75rem; padding: 0.25rem 0.55rem;" title="Stop adding range at current frame (Hotkey: ])">
                                    ] Set Out <span id="val-mark-out" style="font-family: var(--font-mono); color: #38bdf8; font-weight: 700; margin-left: 0.2rem;">--</span>
                                </button>
                                <button id="btn-add-clip" class="btn btn-sm btn-primary" style="background: linear-gradient(135deg, #10b981, #059669); font-weight: 600; font-size: 0.75rem; padding: 0.25rem 0.65rem;" title="Add selected orbital as a clip on timeline spanning [In..Out]">
                                    ➕ Add Clip
                                </button>
                                <button id="btn-clear-inout" class="btn btn-sm btn-ghost" style="padding: 0.25rem 0.4rem; font-size: 0.75rem; color: var(--text-dim);" title="Clear In/Out pins">✕</button>
                            </div>

                            <!-- Isovalue threshold -->
                            <div style="display: flex; align-items: center; gap: 0.45rem; font-size: 0.75rem;">
                                <span style="color: var(--text-muted); font-weight: 600;">Iso:</span>
                                <input type="range" id="isovalue-slider" min="0.010" max="0.100" step="0.005" value="${this.state.isovalue}" style="width: 85px; cursor: pointer;">
                                <span id="isovalue-val-display" style="font-family: var(--font-mono); color: #38bdf8; font-weight: 700; font-size: 0.75rem;">${this.state.isovalue.toFixed(3)}</span>
                            </div>
                        </div>
                    </div>

                    <!-- Right Column: Orbital Media Pool & Audition Inspector -->
                    <div class="card" style="display: flex; flex-direction: column; height: 505px; padding: 0.85rem; background: rgba(17, 24, 39, 0.85); border: 1px solid var(--border-color); border-radius: var(--radius-md); overflow: hidden;">
                        
                        <!-- Pool Header & Search -->
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.45rem;">
                            <div>
                                <span style="font-size: 0.84rem; font-weight: 700; color: var(--text-main);">Orbital Media Pool</span>
                                <span id="pool-count-badge" class="chip blue" style="font-size: 0.62rem; margin-left: 0.35rem;">32 orbitals</span>
                            </div>
                            <button id="btn-calc-flux" class="btn btn-sm btn-secondary" style="font-size: 0.7rem; padding: 0.2rem 0.5rem;" title="Compute net electron redistribution between current In/Out frames to rank reacting orbitals">
                                ⚡ Rank Density Flux
                            </button>
                        </div>

                        <!-- Filter Tabs -->
                        <div style="display: flex; flex-wrap: wrap; gap: 0.25rem; margin-bottom: 0.45rem; border-bottom: 1px solid var(--border-color); padding-bottom: 0.4rem;">
                            <button class="filter-tab active" data-cat="all" style="background: none; border: none; color: #38bdf8; font-size: 0.72rem; font-weight: 700; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">All</button>
                            <button class="filter-tab" data-cat="reacting" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">⚡ Reacting Center</button>
                            <button class="filter-tab" data-cat="canon_pi" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">✨ Conjugated π</button>
                            <button class="filter-tab" data-cat="formyl" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">🔷 Formyl</button>
                            <button class="filter-tab" data-cat="methyl" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">⚪ Methyl</button>
                            <button class="filter-tab" data-cat="lp" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">🔋 Lone Pairs</button>
                            <button class="filter-tab" data-cat="bond" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">Bonds (σ/π)</button>
                            <button class="filter-tab" data-cat="antibond" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">Antibonds (σ*/π*)</button>
                            <button class="filter-tab" data-cat="virt" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">🎯 Frontier</button>
                            <button class="filter-tab" data-cat="flux" style="background: none; border: none; color: #94a3b8; font-size: 0.72rem; font-weight: 600; cursor: pointer; padding: 0.15rem 0.35rem; border-radius: 4px;">🔥 High Flux</button>
                        </div>

                        <!-- Review Progress & Checklist Toolbar -->
                        <div style="display: flex; flex-direction: column; gap: 0.35rem; margin-bottom: 0.45rem; background: rgba(0,0,0,0.3); border: 1px solid var(--border-color); border-radius: 6px; padding: 0.4rem 0.55rem;">
                            <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.7rem;">
                                <div style="display: flex; align-items: center; gap: 0.4rem;">
                                    <span style="color: #94a3b8; font-weight: 600;">Audit Checklist:</span>
                                    <span id="review-stat-badge" style="color: #34d399; font-weight: 700; font-family: var(--font-mono);">0 / 27 Tested (0%)</span>
                                </div>
                                <div style="display: flex; align-items: center; gap: 0.45rem;">
                                    <button id="btn-mark-all-reviewed" style="background: none; border: none; color: #38bdf8; font-size: 0.65rem; cursor: pointer; text-decoration: underline;" title="Mark all orbitals as tested">✓ Mark All</button>
                                    <button id="btn-reset-reviewed" style="background: none; border: none; color: #64748b; font-size: 0.65rem; cursor: pointer; text-decoration: underline;" title="Clear tested status checkmarks">↺ Reset</button>
                                </div>
                            </div>
                            <div style="width: 100%; height: 5px; background: rgba(255,255,255,0.08); border-radius: 3px; overflow: hidden;">
                                <div id="review-progress-fill" style="width: 0%; height: 100%; background: linear-gradient(90deg, #10b981, #34d399); transition: width 0.25s ease;"></div>
                            </div>
                            <div style="display: flex; justify-content: space-between; align-items: center; gap: 0.3rem; margin-top: 0.15rem; flex-wrap: wrap;">
                                <!-- Review status filter pills -->
                                <div style="display: flex; gap: 0.2rem;" id="review-filter-group">
                                    <button class="review-filter-btn active" data-filter="all" style="background: rgba(255,255,255,0.1); border: 1px solid rgba(255,255,255,0.15); color: #f8fafc; font-size: 0.65rem; padding: 0.1rem 0.35rem; border-radius: 3px; cursor: pointer;">All</button>
                                    <button class="review-filter-btn" data-filter="unreviewed" style="background: none; border: 1px solid transparent; color: #94a3b8; font-size: 0.65rem; padding: 0.1rem 0.35rem; border-radius: 3px; cursor: pointer;" title="Show only orbitals you have not tested yet">○ Untested Only</button>
                                    <button class="review-filter-btn" data-filter="reviewed" style="background: none; border: 1px solid transparent; color: #94a3b8; font-size: 0.65rem; padding: 0.1rem 0.35rem; border-radius: 3px; cursor: pointer;">✓ Tested</button>
                                    <button class="review-filter-btn" data-filter="timeline" style="background: none; border: 1px solid transparent; color: #94a3b8; font-size: 0.65rem; padding: 0.1rem 0.35rem; border-radius: 3px; cursor: pointer;">🎬 On Timeline</button>
                                </div>
                                <!-- Sort mode select -->
                                <div style="display: flex; align-items: center; gap: 0.25rem;">
                                    <span style="font-size: 0.65rem; color: #64748b;">Sort:</span>
                                    <select id="select-pool-sort" style="background: rgba(15,23,42,0.9); border: 1px solid var(--border-color); color: #38bdf8; font-size: 0.68rem; font-weight: 600; padding: 0.1rem 0.25rem; border-radius: 3px; cursor: pointer;">
                                        <option value="chemical">📍 Fixed Chemical Order</option>
                                        <option value="untested">○ Untested First</option>
                                        <option value="flux">⚡ Reactivity Flux</option>
                                        <option value="energy">📈 Energy</option>
                                    </select>
                                </div>
                            </div>
                        </div>

                        <!-- Search Box -->
                        <div style="margin-bottom: 0.45rem;">
                            <input type="text" id="orbital-search-input" placeholder="Search by name (e.g. C1-Cl2, LUMO, π1) or atom..." style="width: 100%; background: rgba(0,0,0,0.4); border: 1px solid var(--border-color); border-radius: 4px; padding: 0.3rem 0.5rem; font-size: 0.74rem; color: #f8fafc;">
                        </div>

                        <!-- Current Audition Pill Banner -->
                        <div id="audition-banner" style="background: rgba(6, 182, 212, 0.12); border: 1px solid rgba(6, 182, 212, 0.35); border-radius: 6px; padding: 0.35rem 0.6rem; margin-bottom: 0.5rem; display: flex; justify-content: space-between; align-items: center; font-size: 0.72rem; gap: 0.5rem;">
                            <div style="display: flex; align-items: center; gap: 0.4rem; flex: 1; min-width: 0;">
                                <span style="color: #06b6d4; font-weight: 800; white-space: nowrap;">Audition:</span>
                                <span id="audition-title" style="font-weight: 700; color: #f8fafc; font-family: var(--font-mono); overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">None</span>
                            </div>
                            <div style="display: flex; align-items: center; gap: 0.25rem;">
                                <button id="btn-prev-orbital" class="btn btn-sm btn-secondary" style="padding: 0.15rem 0.45rem; font-size: 0.65rem;" title="Previous orbital (P)">◀</button>
                                <button id="btn-next-orbital" class="btn btn-sm btn-secondary" style="padding: 0.15rem 0.45rem; font-size: 0.65rem;" title="Next orbital">▶</button>
                                <button id="btn-next-untested" class="btn btn-sm btn-primary" style="padding: 0.15rem 0.55rem; font-size: 0.65rem; background: #059669; border-color: #10b981; font-weight: 700; display: flex; align-items: center; gap: 0.25rem;" title="Advance to next untested orbital (Hotkey: N)">
                                    <span>▶ Next Untested</span>
                                </button>
                            </div>
                        </div>

                        <!-- Scrollable Orbital Cards List -->
                        <div id="orbital-pool-list" style="flex: 1; overflow-y: auto; display: flex; flex-direction: column; gap: 0.35rem; padding-right: 0.2rem;">
                            <div style="text-align: center; color: var(--text-dim); padding: 2rem 0; font-size: 0.78rem;">Loading orbitals for frame 1...</div>
                        </div>
                    </div>
                </div>

                <!-- Bottom Area: Kdenlive Multi-Track Timeline -->
                <div class="card" style="padding: 0.75rem 1rem; background: #0c111c; border: 1px solid var(--border-color); border-radius: var(--radius-md); display: flex; flex-direction: column; gap: 0.4rem;">
                    
                    <!-- Timeline Header & Controls -->
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.15rem;">
                        <div style="display: flex; align-items: center; gap: 0.6rem;">
                            <span style="font-weight: 700; font-size: 0.82rem; color: var(--text-main); text-transform: uppercase;">Timeline Tracks</span>
                            <span id="track-count-badge" class="chip blue" style="font-size: 0.65rem;">${this.state.timeline.tracks.length} tracks</span>
                            <span style="font-size: 0.72rem; color: var(--text-dim);">Click or drag anywhere on ruler to seek</span>
                        </div>
                        <div style="display: flex; gap: 0.5rem; align-items: center; font-size: 0.72rem; color: var(--text-dim);">
                            <span>👁️ Default ON</span>
                            <span>•</span>
                            <span>✕ Delete Clip</span>
                        </div>
                    </div>

                    <!-- Timeline Workspace Container -->
                    <div id="timeline-container" style="position: relative; width: 100%; border: 1px solid var(--border-color); border-radius: 6px; background: #080c14; overflow-x: hidden;">
                        
                        <!-- Track 0: Ruler & Energy Profile Track -->
                        <div style="display: grid; grid-template-columns: 200px 1fr; border-bottom: 1px solid rgba(255,255,255,0.1); background: #0f172a; height: 52px;">
                            <div style="padding: 0.4rem 0.6rem; border-right: 1px solid rgba(255,255,255,0.08); display: flex; flex-direction: column; justify-content: center;">
                                <div style="font-weight: 700; font-size: 0.74rem; color: #38bdf8;">Energy Profile</div>
                                <div style="font-size: 0.65rem; color: var(--text-dim);">Erel (kcal/mol)</div>
                            </div>
                            <div id="timeline-ruler" style="position: relative; width: 100%; height: 100%; cursor: pointer;">
                                <svg id="ruler-energy-svg" style="width: 100%; height: 100%; display: block;"></svg>
                                
                                <!-- In/Out Shaded Band Overlay -->
                                <div id="inout-band-overlay" style="display: none; position: absolute; top: 0; bottom: 0; background: rgba(56, 189, 248, 0.15); border-left: 2px solid #38bdf8; border-right: 2px solid #38bdf8; pointer-events: none;"></div>
                            </div>
                        </div>

                        <!-- Scrollable Track Lanes Container -->
                        <div id="timeline-tracks-body" style="display: flex; flex-direction: column; max-height: 230px; overflow-y: auto;">
                            <!-- Populated dynamically by renderTimelineTracks -->
                        </div>

                        <!-- Global Needle Playhead (Sweeps across ruler and all tracks) -->
                        <div id="timeline-playhead-needle" style="position: absolute; top: 0; bottom: 0; left: 200px; width: 2px; background: #38bdf8; box-shadow: 0 0 8px #38bdf8; z-index: 20; pointer-events: none; transition: left 0.05s linear;">
                            <div style="width: 10px; height: 10px; background: #38bdf8; transform: translate(-4px, 0) rotate(45deg); border-radius: 2px;"></div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        this.initEventListeners();
        this.init3DViewer();
        await this.loadFrame(0, true);
        this.renderEnergyProfileSVG();
        this.renderTimelineTracks();
    },

    initEventListeners() {
        // Transport
        const btnPlay = document.getElementById('btn-transport-play');
        const btnPrev = document.getElementById('btn-transport-prev');
        const btnNext = document.getElementById('btn-transport-next');
        const btnStart = document.getElementById('btn-transport-start');
        const btnEnd = document.getElementById('btn-transport-end');

        if (btnPlay) btnPlay.onclick = () => this.togglePlay();
        if (btnPrev) btnPrev.onclick = () => this.stepFrame(-1);
        if (btnNext) btnNext.onclick = () => this.stepFrame(1);
        if (btnStart) btnStart.onclick = () => this.loadFrame(0);
        if (btnEnd) btnEnd.onclick = () => this.loadFrame(this.state.totalFrames - 1);

        // In/Out punching
        const btnMarkIn = document.getElementById('btn-mark-in');
        const btnMarkOut = document.getElementById('btn-mark-out');
        const btnAddClip = document.getElementById('btn-add-clip');
        const btnClearInOut = document.getElementById('btn-clear-inout');

        if (btnMarkIn) btnMarkIn.onclick = () => this.setMarkIn();
        if (btnMarkOut) btnMarkOut.onclick = () => this.setMarkOut();
        if (btnAddClip) btnAddClip.onclick = () => this.addCurrentClipToTimeline();
        if (btnClearInOut) btnClearInOut.onclick = () => this.clearInOut();

        // Charge Distribution Toggle
        const btnCharges = document.getElementById('btn-toggle-charges');
        if (btnCharges) btnCharges.onclick = () => this.toggleChargeDistribution();

        // Isovalue slider
        const isoSlider = document.getElementById('isovalue-slider');
        const isoDisplay = document.getElementById('isovalue-val-display');
        if (isoSlider) {
            isoSlider.oninput = (e) => {
                const val = parseFloat(e.target.value);
                this.state.isovalue = val;
                if (isoDisplay) isoDisplay.innerText = val.toFixed(3);
            };
            isoSlider.onchange = () => {
                this.renderCurrent3DScene();
            };
        }

        // Density Flux helper
        const btnFlux = document.getElementById('btn-calc-flux');
        if (btnFlux) btnFlux.onclick = () => this.fetchDensityFlux();

        // Search and filters
        const searchInput = document.getElementById('orbital-search-input');
        if (searchInput) {
            searchInput.oninput = (e) => {
                this.state.searchQuery = e.target.value.toLowerCase().trim();
                this.renderOrbitalPoolList();
            };
        }

        document.querySelectorAll('.filter-tab').forEach(tab => {
            tab.onclick = () => {
                document.querySelectorAll('.filter-tab').forEach(t => {
                    t.classList.remove('active');
                    t.style.color = '#94a3b8';
                });
                tab.classList.add('active');
                tab.style.color = '#38bdf8';
                this.state.activeCategory = tab.dataset.cat;
                this.renderOrbitalPoolList();
            };
        });

        // Review filter buttons
        document.querySelectorAll('.review-filter-btn').forEach(btn => {
            btn.onclick = () => {
                document.querySelectorAll('.review-filter-btn').forEach(b => {
                    b.classList.remove('active');
                    b.style.background = 'none';
                    b.style.borderColor = 'transparent';
                    b.style.color = '#94a3b8';
                });
                btn.classList.add('active');
                btn.style.background = 'rgba(255,255,255,0.1)';
                btn.style.borderColor = 'rgba(255,255,255,0.15)';
                btn.style.color = '#f8fafc';
                this.state.reviewFilter = btn.dataset.filter;
                this.renderOrbitalPoolList();
            };
        });

        // Sort mode selector
        const sortSelect = document.getElementById('select-pool-sort');
        if (sortSelect) {
            sortSelect.value = this.state.sortMode || 'chemical';
            sortSelect.onchange = (e) => {
                this.state.sortMode = e.target.value;
                this.renderOrbitalPoolList();
            };
        }

        // Mark all as reviewed
        const btnMarkAll = document.getElementById('btn-mark-all-reviewed');
        if (btnMarkAll) {
            btnMarkAll.onclick = () => {
                const catalog = this.getCatalogList();
                catalog.forEach(item => {
                    this.state.testedOrbitals.add(item.id);
                    if (item.chem_key) this.state.testedOrbitals.add(item.chem_key);
                    if (item.name) this.state.testedOrbitals.add(item.name);
                });
                try { localStorage.setItem('mw_tested_orbitals', JSON.stringify([...this.state.testedOrbitals])); } catch (e) {}
                this.renderOrbitalPoolList();
            };
        }

        // Reset reviewed checkmarks
        const btnResetRev = document.getElementById('btn-reset-reviewed');
        if (btnResetRev) {
            btnResetRev.onclick = () => {
                if (confirm('Reset review checkmarks?')) {
                    this.state.testedOrbitals.clear();
                    try { localStorage.removeItem('mw_tested_orbitals'); } catch (e) {}
                    this.renderOrbitalPoolList();
                }
            };
        }

        // Audition Banner Next Untested / Prev / Next buttons
        const btnNextUntested = document.getElementById('btn-next-untested');
        if (btnNextUntested) {
            btnNextUntested.onclick = () => this.auditionNextUntested();
        }

        const btnPrevOrb = document.getElementById('btn-prev-orbital');
        if (btnPrevOrb) {
            btnPrevOrb.onclick = () => this.auditionStep(-1);
        }

        const btnNextOrb = document.getElementById('btn-next-orbital');
        if (btnNextOrb) {
            btnNextOrb.onclick = () => this.auditionStep(1);
        }

        // Precompute All
        const btnPrecompute = document.getElementById('btn-precompute-timeline');
        const btnPreStop = document.getElementById('btn-precompute-stop');
        if (btnPrecompute) btnPrecompute.onclick = () => this.runPrecomputeTimeline();
        if (btnPreStop) btnPreStop.onclick = () => { this.state.precomputeAbort = true; };

        // Export cube
        const btnDl = document.getElementById('btn-dl-cube');
        if (btnDl) btnDl.onclick = () => this.downloadCurrentCube();

        // Timeline ruler scrubbing
        const ruler = document.getElementById('timeline-ruler');
        if (ruler) {
            const handleRulerClick = (e) => {
                const rect = ruler.getBoundingClientRect();
                const frac = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width));
                const targetFrame = Math.round(frac * (this.state.totalFrames - 1));
                this.pause();
                this.loadFrame(targetFrame);
            };
            ruler.onmousedown = (e) => {
                handleRulerClick(e);
                const onMouseMove = (ev) => handleRulerClick(ev);
                const onMouseUp = () => {
                    window.removeEventListener('mousemove', onMouseMove);
                    window.removeEventListener('mouseup', onMouseUp);
                };
                window.addEventListener('mousemove', onMouseMove);
                window.addEventListener('mouseup', onMouseUp);
            };
        }

        // Keyboard Shortcuts
        window.onkeydown = (e) => {
            if (e.target.tagName === 'INPUT') return;
            if (e.code === 'Space') {
                e.preventDefault();
                this.togglePlay();
            } else if (e.code === 'ArrowLeft') {
                e.preventDefault();
                this.stepFrame(-1);
            } else if (e.code === 'ArrowRight') {
                e.preventDefault();
                this.stepFrame(1);
            } else if (e.key === 'n' || e.key === 'N') {
                e.preventDefault();
                this.auditionNextUntested();
            } else if (e.key === 'p' || e.key === 'P') {
                e.preventDefault();
                this.auditionStep(-1);
            } else if (e.key === '[') {
                this.setMarkIn();
            } else if (e.key === ']') {
                this.setMarkOut();
            } else if (e.key === 'c' || e.key === 'C') {
                this.toggleChargeDistribution();
            } else if (['1', '2', '3', '4', '5'].includes(e.key)) {
                this.toggleActiveLobeByHotkey(parseInt(e.key, 10));
            }
        };
    },

    setMarkIn() {
        this.state.inFrame = this.state.currentFrame;
        const el = document.getElementById('val-mark-in');
        if (el) el.innerText = `f${this.state.inFrame + 1}`;
        this.updateInOutOverlay();
    },

    setMarkOut() {
        this.state.outFrame = this.state.currentFrame;
        const el = document.getElementById('val-mark-out');
        if (el) el.innerText = `f${this.state.outFrame + 1}`;
        this.updateInOutOverlay();
    },

    clearInOut() {
        this.state.inFrame = null;
        this.state.outFrame = null;
        const elIn = document.getElementById('val-mark-in');
        const elOut = document.getElementById('val-mark-out');
        if (elIn) elIn.innerText = '--';
        if (elOut) elOut.innerText = '--';
        this.updateInOutOverlay();
    },

    updateInOutOverlay() {
        const band = document.getElementById('inout-band-overlay');
        if (!band) return;
        const total = this.state.totalFrames;
        if (this.state.inFrame === null || this.state.outFrame === null || total <= 1) {
            band.style.display = 'none';
            return;
        }
        const start = Math.min(this.state.inFrame, this.state.outFrame);
        const end = Math.max(this.state.inFrame, this.state.outFrame);
        const leftPct = (start / (total - 1)) * 100;
        const widthPct = Math.max(0.5, ((end - start) / (total - 1)) * 100);
        band.style.display = 'block';
        band.style.left = `${leftPct}%`;
        band.style.width = `${widthPct}%`;
    },

    toggleChargeDistribution() {
        this.state.showCharges = !this.state.showCharges;
        const ind = document.getElementById('charge-toggle-indicator');
        const leg = document.getElementById('charges-legend');
        const btn = document.getElementById('btn-toggle-charges');
        if (ind) ind.style.background = this.state.showCharges ? '#38bdf8' : '#64748b';
        if (btn) btn.style.borderColor = this.state.showCharges ? '#38bdf8' : 'rgba(255,255,255,0.15)';
        if (leg) leg.style.display = this.state.showCharges ? 'flex' : 'none';
        this.renderCurrent3DScene();
    },

    toggleActiveLobeByHotkey(index) {
        const activeList = Array.from(this.getActiveOrbitalsAtFrame(this.state.currentFrame));
        if (index >= 1 && index <= activeList.length) {
            const key = activeList[index - 1];
            const currentVis = this.isOrbitalVisibleNow(key);
            this.state.manualOverrides.set(key, !currentVis);
            this.renderCurrent3DScene();
            this.updateActiveLobesPills();
        }
    },

    init3DViewer() {
        const container = document.getElementById('orb-viewport-3d');
        if (!container) return;
        container.innerHTML = '';
        this.state.viewer = $3Dmol.createViewer(container, {
            defaultcolors: $3Dmol.rasmolElementColors,
            backgroundColor: '0x0b0f17'
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

        // Update Top Badges
        const statusBadge = document.getElementById('orb-status-badge');
        if (statusBadge) statusBadge.innerText = `Frame ${frameIndex + 1} / ${this.state.totalFrames}`;

        const energyBadge = document.getElementById('orb-energy-badge');
        if (energyBadge && this.state.energyProfile.length > frameIndex) {
            const val = this.state.energyProfile[frameIndex];
            energyBadge.innerText = `E: ${val >= 0 ? '+' : ''}${val.toFixed(2)} kcal/mol`;
        }

        const landmarkChip = document.getElementById('orb-landmark-chip');
        if (landmarkChip) {
            const lm = this.state.landmarks.find(l => l.frame_index === frameIndex);
            if (lm) {
                landmarkChip.style.display = 'inline-flex';
                landmarkChip.innerText = `⭐ ${lm.name}`;
            } else {
                landmarkChip.style.display = 'none';
            }
        }

        // Update timeline playhead needle
        this.updatePlayheadNeedle();

        // Fetch frame metadata if not cached
        if (!this.state.frameMetaCache[frameIndex]) {
            try {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_orbitals', frame_index: frameIndex })
                });
                const data = await res.json();
                if (data.success) {
                    this.state.frameMetaCache[frameIndex] = data;
                }
            } catch (err) {
                console.error("Frame orbitals fetch failed:", err);
            }
        }

        // Keep active audition tracking this catalog item across frames
        if (this.state.auditionCatalogId) {
            const catItem = this.getCatalogItemById(this.state.auditionCatalogId);
            if (catItem) {
                const resolved = this.resolveOrbitalKeyAtFrame(catItem, frameIndex);
                if (resolved) {
                    this.state.auditionOrbital = resolved;
                }
            }
        }

        this.renderOrbitalPoolList();
        this.updateActiveLobesPills();
        await this.renderCurrent3DScene(resetCamera);
    },

    updatePlayheadNeedle() {
        const needle = document.getElementById('timeline-playhead-needle');
        const ruler = document.getElementById('timeline-ruler');
        if (!needle || !ruler || this.state.totalFrames <= 1) return;

        const rulerWidth = ruler.offsetWidth;
        const frac = this.state.currentFrame / (this.state.totalFrames - 1);
        // Header left column width is 200px
        const leftPx = 200 + frac * rulerWidth;
        needle.style.left = `${leftPx}px`;
    },

    getActiveOrbitalsAtFrame(frameIndex) {
        const active = new Set();
        // Check timeline tracks
        this.state.timeline.tracks.forEach(trk => {
            const hasClip = trk.clips.some(c => frameIndex >= c.start && frameIndex <= c.end);
            if (hasClip) {
                active.add(trk.orbital_key);
            }
        });
        // If auditioning an orbital, include it
        if (this.state.auditionOrbital) {
            active.add(this.state.auditionOrbital);
        }
        return active;
    },

    isOrbitalVisibleNow(orbitalKey) {
        // If user manually clicked override in lecture
        if (this.state.manualOverrides.has(orbitalKey)) {
            return this.state.manualOverrides.get(orbitalKey);
        }
        // If it's the audition orbital, it's visible
        if (this.state.auditionOrbital === orbitalKey) {
            return true;
        }
        // Look up track default_on
        const trk = this.state.timeline.tracks.find(t => t.orbital_key === orbitalKey);
        if (trk) {
            const hasClip = trk.clips.some(c => this.state.currentFrame >= c.start && this.state.currentFrame <= c.end);
            return hasClip && trk.default_on;
        }
        return false;
    },

    updateActiveLobesPills() {
        const container = document.getElementById('active-lobes-chips');
        if (!container) return;

        const activeKeys = Array.from(this.getActiveOrbitalsAtFrame(this.state.currentFrame));
        if (activeKeys.length === 0) {
            container.innerHTML = `<span style="color: var(--text-dim); font-size: 0.7rem;">None on timeline at frame ${this.state.currentFrame + 1} (Audition an orbital or add a clip)</span>`;
            return;
        }

        container.innerHTML = activeKeys.map((key, idx) => {
            const trk = this.state.timeline.tracks.find(t => t.orbital_key === key);
            const isAudition = (this.state.auditionOrbital === key);
            const isVisible = this.isOrbitalVisibleNow(key);
            const title = trk ? trk.title : (isAudition ? `${key} (${this.state.auditionLabel})` : key);
            const color = trk ? trk.color : '#06b6d4';
            const hotkey = idx + 1 <= 5 ? `[${idx + 1}]` : '';

            return `
                <div class="active-lobe-chip" data-key="${key}" style="display: inline-flex; align-items: center; gap: 0.35rem; padding: 0.2rem 0.5rem; border-radius: 4px; font-size: 0.72rem; cursor: pointer; background: ${isVisible ? 'rgba(0,0,0,0.5)' : 'rgba(0,0,0,0.2)'}; border: 1px solid ${isVisible ? color : 'rgba(255,255,255,0.1)'}; opacity: ${isVisible ? 1.0 : 0.45}; transition: all 0.2s;" title="Click or press ${idx + 1} to toggle visibility">
                    <span style="display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: ${color};"></span>
                    <span style="font-family: var(--font-mono); font-size: 0.65rem; color: #38bdf8; font-weight: 700;">${hotkey}</span>
                    <span style="color: #f8fafc; font-weight: 600;">${title}</span>
                    ${isAudition ? '<span class="chip blue" style="font-size: 0.55rem; padding: 0.05rem 0.3rem;">AUDITION</span>' : ''}
                </div>
            `;
        }).join('');

        container.querySelectorAll('.active-lobe-chip').forEach(chip => {
            chip.onclick = () => {
                const k = chip.dataset.key;
                const cur = this.isOrbitalVisibleNow(k);
                this.state.manualOverrides.set(k, !cur);
                this.renderCurrent3DScene();
                this.updateActiveLobesPills();
            };
        });
    },

    async renderCurrent3DScene(resetCamera = false) {
        if (!this.state.viewer) return;
        const frameIndex = this.state.currentFrame;
        const frameData = this.state.frameMetaCache[frameIndex];
        if (!frameData || !frameData.coords) return;

        const coords = frameData.coords;
        const symbols = this.state.symbols;
        const charges = frameData.charges || [];

        this.state.viewer.clear();
        const xyz = this._buildXYZ(symbols, coords);
        this.state.viewer.addModel(xyz, 'xyz');

        // Style atoms & bonds
        if (this.state.showCharges && charges.length === symbols.length) {
            // Apply partial charge coloring and labels
            for (let i = 0; i < symbols.length; i++) {
                const q = charges[i];
                // Blue for positive, White for 0, Red for negative
                let colorHex = '#f8fafc';
                if (q < 0) {
                    const norm = Math.min(1.0, Math.abs(q));
                    const r = 255;
                    const g = Math.round(255 * (1 - norm));
                    const b = Math.round(255 * (1 - norm));
                    colorHex = `rgb(${r},${g},${b})`;
                } else if (q > 0) {
                    const norm = Math.min(1.0, q);
                    const r = Math.round(255 * (1 - norm));
                    const g = Math.round(255 * (1 - norm));
                    const b = 255;
                    colorHex = `rgb(${r},${g},${b})`;
                }

                this.state.viewer.setStyle({ index: i }, {
                    stick: { radius: 0.14, color: colorHex },
                    sphere: { scale: 0.30, color: colorHex }
                });

                // Add 3D text label for significant charges
                if (Math.abs(q) >= 0.15) {
                    const qStr = `${q >= 0 ? '+' : ''}${q.toFixed(2)}`;
                    this.state.viewer.addLabel(qStr, {
                        position: { x: coords[i][0], y: coords[i][1] + 0.35, z: coords[i][2] },
                        backgroundColor: 'rgba(15, 23, 42, 0.85)',
                        fontColor: q < 0 ? '#fca5a5' : '#93c5fd',
                        fontSize: 10,
                        font: 'monospace',
                        showBackground: true,
                        inFront: true
                    });
                }
            }
        } else {
            this.state.viewer.setStyle({}, {
                stick: { radius: 0.14, colorscheme: 'rasmol' },
                sphere: { scale: 0.28, colorscheme: 'rasmol' }
            });
        }

        // Render Active 3D Lobes
        const activeKeys = Array.from(this.getActiveOrbitalsAtFrame(frameIndex));
        const iso = this.state.isovalue;

        for (const orbKey of activeKeys) {
            if (!this.isOrbitalVisibleNow(orbKey)) continue;

            const cacheKey = `${frameIndex}_${orbKey}`;
            let cubeText = this.state.cubeCache[cacheKey];

            if (!cubeText && !this.state.isPrecomputing) {
                const spinner = document.getElementById('orb-loading-spinner');
                if (spinner) spinner.style.display = 'flex';
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({
                            action_type: 'get_cube',
                            frame_index: frameIndex,
                            orbital: orbKey,
                            isovalue: iso
                        })
                    });
                    const d = await res.json();
                    if (d.success && d.cube) {
                        cubeText = d.cube;
                        this.state.cubeCache[cacheKey] = cubeText;
                    }
                } catch (e) {
                    console.error("Cube fetch error:", e);
                } finally {
                    if (spinner) spinner.style.display = 'none';
                }
            }

            if (cubeText && typeof $3Dmol.VolumeData === 'function') {
                const trk = this.state.timeline.tracks.find(t => t.orbital_key === orbKey);
                const baseColor = trk ? trk.color : '#06b6d4';
                try {
                    const vol = new $3Dmol.VolumeData(cubeText, 'cube');
                    this.state.viewer.addIsosurface(vol, {
                        isoval: iso,
                        color: baseColor,
                        alpha: 0.78,
                        smoothness: 1
                    });
                    this.state.viewer.addIsosurface(vol, {
                        isoval: -iso,
                        color: '#f97316',
                        alpha: 0.78,
                        smoothness: 1
                    });
                } catch (err) {
                    console.warn("VolumeData render error:", err);
                }
            }
        }

        if (resetCamera) {
            this.state.viewer.zoomTo();
        }
        this.state.viewer.render();
    },

    CANONICAL_ORBITAL_CATALOG: [
        // ⚡ Reacting Center (C1, Cl2, O3, O9)
        {
            id: 'bond_C1_Cl2',
            chem_key: 'sigma_C1-Cl2',
            group: 'reacting',
            groupTitle: '⚡ Reacting Center (C1, Cl2, O3, O9)',
            name: 'σ(C1–Cl2)',
            desc: 'Cleaving C–Cl bond → Chloride (Cl⁻) leaving group',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C1', 'Cl2'],
            matchKeys: ['sigma_C1-Cl2', 'lp_Cl2', 'IBO_20']
        },
        {
            id: 'bond_C1_O3_sigma',
            chem_key: 'sigma_C1-O3',
            group: 'reacting',
            name: 'σ(C1–O3)',
            desc: 'Carbonyl C1=O3 σ framework',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C1', 'O3'],
            matchKeys: ['sigma_C1-O3', 'IBO_19']
        },
        {
            id: 'bond_C1_O3_pi',
            chem_key: 'pi_C1-O3',
            group: 'reacting',
            name: 'π(C1–O3)',
            desc: 'Carbonyl C1=O3 π bond → tetrahedral oxyanion polarization',
            category: 'bonding',
            type: 'pi',
            atoms: ['C1', 'O3'],
            matchKeys: ['pi_C1-O3', 'IBO_21']
        },
        {
            id: 'bond_C1_O9',
            chem_key: 'sigma_C1-O9',
            group: 'reacting',
            name: 'σ(C1–O9) / Nucleophile LP',
            desc: 'Attacking formate O9 lone pair → forming C1–O9 covalent bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C1', 'O9'],
            matchKeys: ['sigma_C1-O9', 'lp_O9', 'IBO_12', 'IBO_04']
        },
        {
            id: 'lp_O3_1',
            chem_key: 'lp_O3',
            group: 'reacting',
            name: 'LP 1 (O3)',
            desc: 'In-plane carbonyl lone pair on acetyl oxygen',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['O3'],
            matchKeys: ['lp_O3', 'IBO_11']
        },
        {
            id: 'lp_O3_2',
            chem_key: 'lp_O3',
            group: 'reacting',
            name: 'LP 2 (O3)',
            desc: 'Out-of-plane carbonyl lone pair on acetyl oxygen',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['O3'],
            matchKeys: ['lp_O3', 'IBO_25']
        },
        {
            id: 'lp_Cl2_1',
            chem_key: 'lp_Cl2',
            group: 'reacting',
            name: 'LP 1 (Cl2)',
            desc: 'Valence lone pair on chlorine atom',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['Cl2'],
            matchKeys: ['lp_Cl2', 'IBO_14']
        },
        {
            id: 'lp_Cl2_2',
            chem_key: 'lp_Cl2',
            group: 'reacting',
            name: 'LP 2 (Cl2)',
            desc: 'Valence lone pair on chlorine atom',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['Cl2'],
            matchKeys: ['lp_Cl2', 'IBO_23']
        },
        {
            id: 'lp_Cl2_3',
            chem_key: 'lp_Cl2',
            group: 'reacting',
            name: 'LP 3 (Cl2)',
            desc: 'Valence lone pair on chlorine atom',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['Cl2'],
            matchKeys: ['lp_Cl2', 'IBO_24']
        },

        // ✨ Conjugated π System (O=C–O–C=O)
        {
            id: 'deloc_pi_1',
            chem_key: 'deloc_pi_1',
            group: 'conjugated',
            groupTitle: '✨ Conjugated π System (O=C–O–C=O)',
            name: 'Delocalized π₁ (0 nodes)',
            desc: 'Completely bonding π orbital across O=C–O–C=O (0.00% methyl contamination)',
            category: 'conjugated_pi',
            type: 'delocalized_pi',
            matchKeys: ['DELOC_PI_01', 'deloc_pi_1', 'CANONICAL_PI_01']
        },
        {
            id: 'deloc_pi_2',
            chem_key: 'deloc_pi_2',
            group: 'conjugated',
            name: 'Delocalized π₂ (1 node)',
            desc: 'Conjugated π orbital with central nodal surface across anhydride core',
            category: 'conjugated_pi',
            type: 'delocalized_pi',
            matchKeys: ['DELOC_PI_02', 'deloc_pi_2', 'CANONICAL_PI_02']
        },
        {
            id: 'deloc_pi_3',
            chem_key: 'deloc_pi_3',
            group: 'conjugated',
            name: 'Delocalized π₃ (2 nodes)',
            desc: 'Highest occupied conjugated π orbital',
            category: 'conjugated_pi',
            type: 'delocalized_pi',
            matchKeys: ['DELOC_PI_03', 'deloc_pi_3', 'CANONICAL_PI_03']
        },
        {
            id: 'deloc_pistar_1',
            chem_key: 'deloc_pistar_1',
            group: 'conjugated',
            name: 'Delocalized π*₁ (LUMO)',
            desc: 'Lowest unoccupied conjugated π* antibonding acceptor orbital',
            category: 'conjugated_pi',
            type: 'delocalized_pistar',
            matchKeys: ['DELOC_PISTAR_01', 'deloc_pistar_1', 'CANONICAL_PISTAR_01']
        },
        {
            id: 'deloc_pistar_2',
            chem_key: 'deloc_pistar_2',
            group: 'conjugated',
            name: 'Delocalized π*₂ (LUMO+1)',
            desc: 'Second unoccupied conjugated π* antibonding level',
            category: 'conjugated_pi',
            type: 'delocalized_pistar',
            matchKeys: ['DELOC_PISTAR_02', 'deloc_pistar_2', 'CANONICAL_PISTAR_02']
        },

        // 🔷 Formyl Partner (C8, O7, H10, O9)
        {
            id: 'bond_C8_O7_sigma',
            chem_key: 'sigma_C8-O7',
            group: 'formyl',
            groupTitle: '🔷 Formyl Partner (C8, O7, H10, O9)',
            name: 'σ(C8–O7)',
            desc: 'Formyl carbonyl C8=O7 σ framework',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C8', 'O7'],
            matchKeys: ['sigma_C8-O7', 'IBO_28']
        },
        {
            id: 'bond_C8_O7_pi',
            chem_key: 'pi_C8-O7',
            group: 'formyl',
            name: 'π(C8–O7)',
            desc: 'Formyl carbonyl C8=O7 π bond',
            category: 'bonding',
            type: 'pi',
            atoms: ['C8', 'O7'],
            matchKeys: ['pi_C8-O7', 'IBO_27']
        },
        {
            id: 'bond_C8_O9',
            chem_key: 'sigma_C8-O9',
            group: 'formyl',
            name: 'σ(C8–O9)',
            desc: 'Formate C8–O9 single bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C8', 'O9'],
            matchKeys: ['sigma_C8-O9', 'IBO_26']
        },
        {
            id: 'bond_C8_H10',
            chem_key: 'sigma_C8-H10',
            group: 'formyl',
            name: 'σ(C8–H10)',
            desc: 'Formyl C–H bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C8', 'H10'],
            matchKeys: ['sigma_C8-H10', 'IBO_22']
        },
        {
            id: 'lp_O7_1',
            chem_key: 'lp_O7',
            group: 'formyl',
            name: 'LP 1 (O7)',
            desc: 'Formyl carbonyl oxygen lone pair',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['O7'],
            matchKeys: ['lp_O7', 'IBO_15']
        },
        {
            id: 'lp_O7_2',
            chem_key: 'lp_O7',
            group: 'formyl',
            name: 'LP 2 (O7)',
            desc: 'Formyl carbonyl oxygen lone pair',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['O7'],
            matchKeys: ['lp_O7', 'IBO_29']
        },
        {
            id: 'lp_O9_1',
            chem_key: 'lp_O9',
            group: 'formyl',
            name: 'LP 1 (O9)',
            desc: 'Non-bonding pair on bridging formate oxygen O9',
            category: 'non_bonding',
            type: 'lone_pair',
            atoms: ['O9'],
            matchKeys: ['lp_O9', 'IBO_30']
        },

        // ⚪ Spectator Methyl Group (-CH3)
        {
            id: 'bond_C0_C1',
            chem_key: 'sigma_C0-C1',
            group: 'methyl',
            groupTitle: '⚪ Spectator Methyl (-CH3)',
            name: 'σ(C0–C1)',
            desc: 'Methyl-to-carbonyl C–C single bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C0', 'C1'],
            matchKeys: ['sigma_C0-C1', 'IBO_16']
        },
        {
            id: 'bond_C0_H4',
            chem_key: 'sigma_C0-H4',
            group: 'methyl',
            name: 'σ(C0–H4)',
            desc: 'Methyl C–H bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C0', 'H4'],
            matchKeys: ['sigma_C0-H4', 'IBO_18']
        },
        {
            id: 'bond_C0_H5',
            chem_key: 'sigma_C0-H5',
            group: 'methyl',
            name: 'σ(C0–H5)',
            desc: 'Methyl C–H bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C0', 'H5'],
            matchKeys: ['sigma_C0-H5', 'IBO_13']
        },
        {
            id: 'bond_C0_H6',
            chem_key: 'sigma_C0-H6',
            group: 'methyl',
            name: 'σ(C0–H6)',
            desc: 'Methyl C–H bond',
            category: 'bonding',
            type: 'sigma',
            atoms: ['C0', 'H6'],
            matchKeys: ['sigma_C0-H6', 'IBO_17']
        },

        // 🎯 Frontier Orbitals & Antibonds
        {
            id: 'mo_lumo',
            chem_key: 'LUMO',
            group: 'frontier',
            groupTitle: '🎯 Frontier Orbitals & Antibonds',
            name: 'LUMO',
            desc: 'Lowest Unoccupied MO; primary electrophilic acceptor lobe',
            category: 'frontier',
            type: 'frontier',
            matchKeys: ['LUMO']
        },
        {
            id: 'mo_homo',
            chem_key: 'HOMO',
            group: 'frontier',
            name: 'HOMO',
            desc: 'Highest Occupied MO; primary nucleophilic donor density',
            category: 'frontier',
            type: 'frontier',
            matchKeys: ['HOMO']
        },
        {
            id: 'antibond_C1_Cl2',
            chem_key: 'sigma_star_C1-Cl2',
            group: 'frontier',
            name: 'σ*(C1–Cl2)',
            desc: 'Antibonding σ* orbital of cleaving C–Cl bond (SN2 backside attack)',
            category: 'antibonding',
            type: 'sigma_star',
            atoms: ['C1', 'Cl2'],
            matchKeys: ['sigma_star_C1-Cl2', 'STAR_20']
        },
        {
            id: 'antibond_C1_O3',
            chem_key: 'pi_star_C1-O3',
            group: 'frontier',
            name: 'π*(C1–O3)',
            desc: 'Carbonyl C1=O3 π* antibonding orbital',
            category: 'antibonding',
            type: 'pi_star',
            atoms: ['C1', 'O3'],
            matchKeys: ['pi_star_C1-O3', 'STAR_21']
        }
    ],

    getCatalogItemById(id) {
        return this.CANONICAL_ORBITAL_CATALOG.find(i => i.id === id) || null;
    },

    getCatalogList() {
        return this.CANONICAL_ORBITAL_CATALOG;
    },

    isOrbitalTested(key, chemKey, label) {
        if (!key) return false;
        return this.state.testedOrbitals.has(key) ||
               (chemKey && this.state.testedOrbitals.has(chemKey)) ||
               (label && this.state.testedOrbitals.has(label));
    },

    markOrbitalTested(key, chemKey, label) {
        if (!key) return;
        this.state.testedOrbitals.add(key);
        if (chemKey) this.state.testedOrbitals.add(chemKey);
        if (label) this.state.testedOrbitals.add(label);
        try {
            localStorage.setItem('mw_tested_orbitals', JSON.stringify([...this.state.testedOrbitals]));
        } catch (e) {}
    },

    toggleOrbitalTested(key, chemKey, label) {
        const isTested = this.isOrbitalTested(key, chemKey, label);
        if (isTested) {
            this.state.testedOrbitals.delete(key);
            if (chemKey) this.state.testedOrbitals.delete(chemKey);
            if (label) this.state.testedOrbitals.delete(label);
        } else {
            this.state.testedOrbitals.add(key);
            if (chemKey) this.state.testedOrbitals.add(chemKey);
            if (label) this.state.testedOrbitals.add(label);
        }
        try {
            localStorage.setItem('mw_tested_orbitals', JSON.stringify([...this.state.testedOrbitals]));
        } catch (e) {}
        this.renderOrbitalPoolList();
    },

    auditionNextUntested() {
        const catalog = this.getCatalogList();
        const untested = catalog.filter(o => !this.isOrbitalTested(o.id, o.chem_key, o.name));
        if (untested.length === 0) {
            alert('🎉 Outstanding! You have reviewed all orbitals in this reaction!');
            return;
        }
        this.setAuditionCatalogItem(untested[0]);
    },

    auditionStep(delta) {
        const catalog = this.getCatalogList();
        let curIdx = catalog.findIndex(o => o.id === this.state.auditionCatalogId);
        if (curIdx < 0) curIdx = 0;
        const nextIdx = (curIdx + delta + catalog.length) % catalog.length;
        this.setAuditionCatalogItem(catalog[nextIdx]);
    },

    setAuditionCatalogItem(item) {
        if (!item) return;
        if (this.state.auditionCatalogId === item.id) {
            // Toggle off
            this.state.auditionCatalogId = null;
            this.state.auditionOrbital = null;
            this.state.auditionLabel = '';
            this.state.auditionChemKey = null;
        } else {
            this.state.auditionCatalogId = item.id;
            this.state.auditionLabel = item.name;
            this.state.auditionChemKey = item.chem_key;
            const resolved = this.resolveOrbitalKeyAtFrame(item, this.state.currentFrame);
            this.state.auditionOrbital = resolved || item.id;
            this.markOrbitalTested(item.id, item.chem_key, item.name);
        }

        const bannerTitle = document.getElementById('audition-title');
        if (bannerTitle) {
            bannerTitle.innerText = this.state.auditionCatalogId ? `${item.name}: ${item.desc || ''}` : 'None';
        }

        this.renderOrbitalPoolList();
        this.updateActiveLobesPills();
        this.renderCurrent3DScene();
    },

    setAuditionOrbital(key, label, chemKey) {
        const catItem = this.getCatalogItemById(key) || this.CANONICAL_ORBITAL_CATALOG.find(i => i.chem_key === chemKey || i.chem_key === key);
        if (catItem) {
            this.setAuditionCatalogItem(catItem);
            return;
        }

        if (this.state.auditionOrbital === key) {
            this.state.auditionOrbital = null;
            this.state.auditionCatalogId = null;
            this.state.auditionLabel = '';
            this.state.auditionChemKey = null;
        } else {
            this.state.auditionOrbital = key;
            this.state.auditionCatalogId = null;
            this.state.auditionLabel = label || key;
            this.state.auditionChemKey = chemKey || null;
            this.markOrbitalTested(key, chemKey, label);
        }

        const bannerTitle = document.getElementById('audition-title');
        if (bannerTitle) {
            bannerTitle.innerText = this.state.auditionOrbital ? `${this.state.auditionOrbital}: ${this.state.auditionLabel}` : 'None';
        }

        this.renderOrbitalPoolList();
        this.updateActiveLobesPills();
        this.renderCurrent3DScene();
    },

    resolveOrbitalKeyAtFrame(item, frameIndex) {
        const frameData = this.state.frameMetaCache[frameIndex];
        if (!frameData) return item.matchKeys ? item.matchKeys[0] : item.id;

        const allCandidates = [
            ...(frameData.canonical_pi || []),
            ...(frameData.all_ibos || []),
            ...(frameData.antibonds || []),
            ...(frameData.frontier_orbitals ? Object.values(frameData.frontier_orbitals) : [])
        ];

        // 1. Exact match by chem_key
        if (item.chem_key) {
            const m = allCandidates.find(o => o.chem_key === item.chem_key);
            if (m) return m.name;
        }

        // 2. Match by matchKeys
        if (item.matchKeys) {
            for (const mk of item.matchKeys) {
                const m = allCandidates.find(o => o.name === mk || o.chem_key === mk);
                if (m) return m.name;
            }
        }

        // 3. Match by atoms & type
        if (item.atoms && item.atoms.length > 0) {
            const m = allCandidates.find(o => {
                if (!o.atoms) return false;
                return item.atoms.every(a => o.atoms.includes(a));
            });
            if (m) return m.name;
        }

        return item.matchKeys ? item.matchKeys[0] : item.id;
    },

    resolveOrbitalDataAtFrame(item, frameIndex) {
        const frameData = this.state.frameMetaCache[frameIndex];
        if (!frameData) return null;

        const allCandidates = [
            ...(frameData.canonical_pi || []),
            ...(frameData.all_ibos || []),
            ...(frameData.antibonds || []),
            ...(frameData.frontier_orbitals ? Object.values(frameData.frontier_orbitals) : [])
        ];

        // 1. Exact chem_key match
        if (item.chem_key) {
            const m = allCandidates.find(o => o.chem_key === item.chem_key);
            if (m) return m;
        }

        // 2. matchKeys match
        if (item.matchKeys) {
            for (const mk of item.matchKeys) {
                const m = allCandidates.find(o => o.name === mk || o.chem_key === mk);
                if (m) return m;
            }
        }

        // 3. atoms match
        if (item.atoms && item.atoms.length > 0) {
            const m = allCandidates.find(o => {
                if (!o.atoms) return false;
                return item.atoms.every(a => o.atoms.includes(a));
            });
            if (m) return m;
        }

        return null;
    },

    isOrbitalOnTimeline(key, label, chemKey) {
        if (!this.state.timeline || !this.state.timeline.tracks) return null;
        for (let i = 0; i < this.state.timeline.tracks.length; i++) {
            const trk = this.state.timeline.tracks[i];
            if (trk.orbital_key === key) return i + 1;
            if (chemKey && trk.chem_key === chemKey) return i + 1;
            if (label && trk.title && trk.title.includes(label)) return i + 1;
            if (key && trk.title && trk.title.includes(key)) return i + 1;
        }
        return null;
    },

    renderOrbitalPoolList() {
        const container = document.getElementById('orbital-pool-list');
        if (!container) return;

        let catalog = [...this.getCatalogList()];

        // Tag each item with current state
        catalog.forEach((item, index) => {
            item.catalogIndex = index;
            item.isTested = this.isOrbitalTested(item.id, item.chem_key, item.name);
            item.onTrack = this.isOrbitalOnTimeline(item.id, item.name, item.chem_key);
            
            const liveData = this.resolveOrbitalDataAtFrame(item, this.state.currentFrame);
            if (liveData) {
                item.currentName = liveData.name;
                item.currentEnergy = liveData.energy_ev;
                item.currentLabel = liveData.label;
                item.currentPop = liveData.pop;
            } else {
                item.currentName = item.matchKeys ? item.matchKeys[0] : item.id;
                item.currentEnergy = null;
                item.currentLabel = item.name;
            }

            // Density flux
            if (this.state.densityFlux) {
                const fluxEntry = this.state.densityFlux.find(f => f.key === item.currentName || f.key === item.chem_key);
                item.delta_pop = fluxEntry ? fluxEntry.delta_pop : 0.0;
            } else {
                item.delta_pop = 0.0;
            }
        });

        // Update Review Progress Stats in Toolbar
        const totalCatalogCount = catalog.length;
        const testedCount = catalog.filter(o => o.isTested).length;
        const pctTested = totalCatalogCount > 0 ? Math.round((testedCount / totalCatalogCount) * 100) : 0;

        const statBadge = document.getElementById('review-stat-badge');
        if (statBadge) {
            statBadge.innerText = `${testedCount} / ${totalCatalogCount} Tested (${pctTested}%)`;
        }
        const progFill = document.getElementById('review-progress-fill');
        if (progFill) {
            progFill.style.width = `${pctTested}%`;
        }

        let filteredList = [...catalog];

        // Filter by Category Tabs
        if (this.state.activeCategory === 'reacting') {
            filteredList = filteredList.filter(o => o.group === 'reacting');
        } else if (this.state.activeCategory === 'canon_pi') {
            filteredList = filteredList.filter(o => o.group === 'conjugated');
        } else if (this.state.activeCategory === 'formyl') {
            filteredList = filteredList.filter(o => o.group === 'formyl');
        } else if (this.state.activeCategory === 'methyl') {
            filteredList = filteredList.filter(o => o.group === 'methyl');
        } else if (this.state.activeCategory === 'lp') {
            filteredList = filteredList.filter(o => o.category === 'non_bonding' || o.type === 'lone_pair');
        } else if (this.state.activeCategory === 'bond') {
            filteredList = filteredList.filter(o => o.category === 'bonding' || o.type === 'sigma' || o.type === 'pi');
        } else if (this.state.activeCategory === 'antibond') {
            filteredList = filteredList.filter(o => o.category === 'antibonding' || o.type === 'sigma_star' || o.type === 'pi_star');
        } else if (this.state.activeCategory === 'virt') {
            filteredList = filteredList.filter(o => o.group === 'frontier' || o.category === 'frontier');
        } else if (this.state.activeCategory === 'flux') {
            filteredList = filteredList.filter(o => (o.delta_pop || 0) >= 0.35);
        }

        // Filter by Review Status Filter Pills
        if (this.state.reviewFilter === 'unreviewed') {
            filteredList = filteredList.filter(o => !o.isTested);
        } else if (this.state.reviewFilter === 'reviewed') {
            filteredList = filteredList.filter(o => o.isTested);
        } else if (this.state.reviewFilter === 'timeline') {
            filteredList = filteredList.filter(o => o.onTrack !== null);
        }

        // Filter by Search Query
        if (this.state.searchQuery) {
            const q = this.state.searchQuery.toLowerCase();
            filteredList = filteredList.filter(o => 
                o.name.toLowerCase().includes(q) || 
                (o.desc && o.desc.toLowerCase().includes(q)) || 
                (o.currentLabel && o.currentLabel.toLowerCase().includes(q)) ||
                (o.chem_key && o.chem_key.toLowerCase().includes(q)) ||
                (o.atoms && o.atoms.some(a => a.toLowerCase().includes(q)))
            );
        }

        // Sorting
        const mode = this.state.sortMode || 'chemical';
        if (mode === 'chemical') {
            // STRICT IMMUTABLE CATALOG ORDER (Zero row jumping!)
            filteredList.sort((a, b) => a.catalogIndex - b.catalogIndex);
        } else if (mode === 'untested') {
            filteredList.sort((a, b) => {
                if (a.isTested !== b.isTested) return (a.isTested ? 1 : 0) - (b.isTested ? 1 : 0);
                return a.catalogIndex - b.catalogIndex;
            });
        } else if (mode === 'flux') {
            filteredList.sort((a, b) => (b.delta_pop || 0) - (a.delta_pop || 0));
        } else if (mode === 'energy') {
            filteredList.sort((a, b) => (a.currentEnergy || 0) - (b.currentEnergy || 0));
        }

        const countBadge = document.getElementById('pool-count-badge');
        if (countBadge) countBadge.innerText = `${filteredList.length} orbitals`;

        if (filteredList.length === 0) {
            container.innerHTML = `
                <div style="text-align: center; color: var(--text-dim); padding: 2.5rem 1rem; font-size: 0.78rem;">
                    No orbitals match the active checklist filter.
                </div>
            `;
            return;
        }

        let currentSectionGroup = '';
        let htmlParts = [];

        filteredList.forEach(item => {
            const isAudition = (this.state.auditionCatalogId === item.id) || 
                               (this.state.auditionOrbital === item.id) || 
                               (item.currentName && this.state.auditionOrbital === item.currentName);
            const hasFlux = (item.delta_pop && item.delta_pop >= 0.35);

            // Group divider in chemical mode
            if (mode === 'chemical' && !this.state.searchQuery && this.state.reviewFilter === 'all' && this.state.activeCategory === 'all') {
                if (item.groupTitle && item.groupTitle !== currentSectionGroup) {
                    currentSectionGroup = item.groupTitle;
                    // Count tested within group
                    const groupItems = catalog.filter(o => o.group === item.group);
                    const grpTested = groupItems.filter(o => o.isTested).length;
                    const grpTotal = groupItems.length;

                    htmlParts.push(`
                        <div style="font-size: 0.68rem; font-weight: 800; color: #94a3b8; text-transform: uppercase; letter-spacing: 0.04em; padding: 0.45rem 0.2rem 0.2rem 0.2rem; margin-top: 0.35rem; border-top: 1px solid rgba(255,255,255,0.08); display: flex; align-items: center; justify-content: space-between;">
                            <span>${currentSectionGroup}</span>
                            <span style="font-size: 0.62rem; color: ${grpTested === grpTotal ? '#34d399' : '#64748b'}; font-weight: 700; font-family: var(--font-mono);">${grpTested}/${grpTotal} Tested</span>
                        </div>
                    `);
                }
            }

            let badgeHtml = '';
            if (item.category === 'non_bonding' || item.type === 'lone_pair') {
                badgeHtml = '<span class="chip purple" style="font-size: 0.6rem; padding: 0.05rem 0.35rem; font-weight: 700;">LONE PAIR</span>';
            } else if (item.category === 'conjugated_pi') {
                badgeHtml = '<span class="chip" style="background: rgba(34, 197, 94, 0.2); color: #4ade80; border: 1px solid rgba(34, 197, 94, 0.4); font-size: 0.6rem; padding: 0.05rem 0.35rem; font-weight: 700;">CONJUGATED π</span>';
            } else if (item.category === 'antibonding') {
                badgeHtml = '<span class="chip" style="background: rgba(244, 63, 94, 0.2); color: #fb7185; border: 1px solid rgba(244, 63, 94, 0.4); font-size: 0.6rem; padding: 0.05rem 0.35rem; font-weight: 700;">ANTIBOND</span>';
            } else if (item.category === 'bonding') {
                badgeHtml = '<span class="chip blue" style="font-size: 0.6rem; padding: 0.05rem 0.35rem; font-weight: 700;">BOND</span>';
            } else {
                badgeHtml = `<span class="chip amber" style="font-size: 0.6rem; padding: 0.05rem 0.35rem; font-weight: 700;">${(item.category || item.type).toUpperCase()}</span>`;
            }

            // Timeline Status Chip
            let timelineChip = '';
            if (item.onTrack) {
                timelineChip = `<span class="chip" style="background: rgba(6, 182, 212, 0.2); color: #38bdf8; border: 1px solid rgba(6, 182, 212, 0.4); font-size: 0.58rem; padding: 0.05rem 0.35rem; font-weight: 700;">🎬 Track ${item.onTrack}</span>`;
            }

            htmlParts.push(`
                <div class="orbital-pool-card" data-id="${item.id}" data-name="${item.name}" data-chemkey="${item.chem_key}" style="display: flex; justify-content: space-between; align-items: center; padding: 0.45rem 0.65rem; background: ${isAudition ? 'rgba(6, 182, 212, 0.16)' : 'rgba(0,0,0,0.3)'}; border: 1px solid ${isAudition ? '#06b6d4' : (item.isTested ? 'rgba(16, 185, 129, 0.35)' : 'var(--border-color)')}; border-radius: 6px; cursor: pointer; transition: all 0.15s; gap: 0.5rem;">
                    
                    <!-- Left: Interactive Checkbox Toggle -->
                    <button class="btn-toggle-test" data-id="${item.id}" data-chemkey="${item.chem_key}" title="Click to mark tested/untested" style="background: ${item.isTested ? 'rgba(16, 185, 129, 0.2)' : 'rgba(255,255,255,0.04)'}; border: 1px solid ${item.isTested ? 'rgba(16, 185, 129, 0.45)' : 'rgba(255,255,255,0.15)'}; color: ${item.isTested ? '#34d399' : '#64748b'}; width: 24px; height: 24px; border-radius: 4px; display: flex; align-items: center; justify-content: center; cursor: pointer; font-size: 0.8rem; font-weight: 800; flex-shrink: 0; transition: all 0.15s;">
                        ${item.isTested ? '✓' : '○'}
                    </button>

                    <!-- Middle: Orbital Chemical Info -->
                    <div style="flex: 1; min-width: 0;">
                        <div style="display: flex; align-items: center; gap: 0.35rem; flex-wrap: wrap;">
                            <span style="font-family: var(--font-mono); font-weight: 700; font-size: 0.78rem; color: #f8fafc;">${item.name}</span>
                            ${badgeHtml}
                            ${item.isTested ? '<span class="chip" style="background: rgba(16, 185, 129, 0.15); color: #34d399; border: 1px solid rgba(16, 185, 129, 0.3); font-size: 0.58rem; padding: 0.05rem 0.35rem; font-weight: 700;">Tested</span>' : ''}
                            ${timelineChip}
                            ${hasFlux ? `<span class="chip" style="background: rgba(239, 68, 68, 0.2); color: #f87171; border-color: rgba(239,68,68,0.4); font-size: 0.58rem; font-weight: 700;">🔥 Δ${item.delta_pop}</span>` : ''}
                        </div>
                        <div style="font-size: 0.72rem; color: #cbd5e1; font-weight: 500; margin-top: 0.15rem; line-height: 1.25;">
                            ${item.desc}
                        </div>
                        <div style="font-size: 0.65rem; color: var(--text-dim); margin-top: 0.15rem; font-family: var(--font-mono); display: flex; gap: 0.6rem; flex-wrap: wrap;">
                            <span>Frame ${this.state.currentFrame + 1}: <b style="color: #94a3b8;">${item.currentLabel}</b></span>
                            ${item.currentEnergy !== null ? `<span>Energy: <b style="color: #94a3b8;">${item.currentEnergy} eV</b></span>` : ''}
                        </div>
                    </div>

                    <!-- Right: Actions -->
                    <div style="display: flex; gap: 0.3rem; flex-shrink: 0;">
                        <button class="btn btn-sm btn-secondary btn-audition" data-id="${item.id}" style="padding: 0.25rem 0.5rem; font-size: 0.68rem; ${isAudition ? 'background: rgba(6, 182, 212, 0.25); border-color: #06b6d4; color: #38bdf8;' : ''}" title="Audition this orbital in 3D (Space / scrubs)">
                            ${isAudition ? '✓ Auditioning' : '👁️ Audition'}
                        </button>
                        <button class="btn btn-sm btn-primary btn-add-pool-clip" data-id="${item.id}" style="padding: 0.25rem 0.5rem; font-size: 0.68rem; background: #0284c7;" title="Add as clip on lecture timeline">
                            ➕ Clip
                        </button>
                    </div>
                </div>
            `);
        });

        container.innerHTML = htmlParts.join('');

        // Wire click on Checkbox button
        container.querySelectorAll('.btn-toggle-test').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const id = btn.dataset.id;
                const catItem = this.getCatalogItemById(id);
                if (catItem) {
                    this.toggleOrbitalTested(catItem.id, catItem.chem_key, catItem.name);
                }
            };
        });

        // Wire click on Audition button
        container.querySelectorAll('.btn-audition').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const id = btn.dataset.id;
                const catItem = this.getCatalogItemById(id);
                if (catItem) {
                    this.setAuditionCatalogItem(catItem);
                }
            };
        });

        // Wire click on entire card
        container.querySelectorAll('.orbital-pool-card').forEach(card => {
            card.onclick = () => {
                const id = card.dataset.id;
                const catItem = this.getCatalogItemById(id);
                if (catItem) {
                    this.setAuditionCatalogItem(catItem);
                }
            };
        });

        // Wire click on Add Clip button
        container.querySelectorAll('.btn-add-pool-clip').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const id = btn.dataset.id;
                const catItem = this.getCatalogItemById(id);
                if (catItem) {
                    this.setAuditionCatalogItem(catItem);
                    this.addCurrentClipToTimeline();
                }
            };
        });
    },

    async fetchDensityFlux() {
        const start = this.state.inFrame !== null ? this.state.inFrame : 0;
        const end = this.state.outFrame !== null ? this.state.outFrame : (this.state.totalFrames - 1);
        const btn = document.getElementById('btn-calc-flux');
        if (btn) btn.innerText = 'Calculating...';

        try {
            const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'get_density_flux',
                    start_frame: start,
                    end_frame: end
                })
            });
            const d = await res.json();
            if (d.success && d.fluxes) {
                this.state.densityFlux = d.fluxes;
                // Switch tab to high flux
                this.state.activeCategory = 'flux';
                document.querySelectorAll('.filter-tab').forEach(t => {
                    t.classList.remove('active');
                    t.style.color = '#94a3b8';
                    if (t.dataset.cat === 'flux') {
                        t.classList.add('active');
                        t.style.color = '#38bdf8';
                    }
                });
                this.renderOrbitalPoolList();
            }
        } catch (e) {
            console.error("Density flux error:", e);
        } finally {
            if (btn) btn.innerText = '⚡ Rank Density Flux';
        }
    },

    addCurrentClipToTimeline() {
        const catItem = this.state.auditionCatalogId ? this.getCatalogItemById(this.state.auditionCatalogId) : null;
        const orbKey = this.state.auditionOrbital || (catItem ? catItem.id : 'IBO_23');
        const orbLabel = catItem ? catItem.name : (this.state.auditionLabel || orbKey);
        const chemKey = catItem ? catItem.chem_key : (this.state.auditionChemKey || null);
        const start = this.state.inFrame !== null ? Math.min(this.state.inFrame, this.state.outFrame || this.state.inFrame) : 0;
        const end = this.state.outFrame !== null ? Math.max(this.state.inFrame || this.state.outFrame, this.state.outFrame) : Math.min(this.state.totalFrames - 1, start + 30);

        let trk = this.state.timeline.tracks.find(t => t.orbital_key === orbKey || (chemKey && t.chem_key === chemKey));
        if (!trk) {
            const colorIdx = this.state.timeline.tracks.length % this.colorPalette.length;
            trk = {
                id: `trk_${(chemKey || orbKey).toLowerCase()}_${Date.now()}`,
                orbital_key: orbKey,
                chem_key: chemKey,
                title: `${orbLabel}`,
                color: this.colorPalette[colorIdx],
                default_on: true,
                clips: []
            };
            this.state.timeline.tracks.push(trk);
        }

        const clipId = `clip_${Date.now()}`;
        trk.clips.push({
            id: clipId,
            start: start,
            end: end,
            title: `${orbLabel} (f${start + 1}–f${end + 1})`
        });

        // Mark as tested
        this.markOrbitalTested(orbKey, chemKey, orbLabel);
        if (catItem) this.markOrbitalTested(catItem.id, catItem.chem_key, catItem.name);

        this.renderTimelineTracks();
        this.renderOrbitalPoolList();
        this.updateActiveLobesPills();
        this.saveTimeline();
    },

    renderTimelineTracks() {
        const container = document.getElementById('timeline-tracks-body');
        if (!container) return;

        const countBadge = document.getElementById('track-count-badge');
        if (countBadge) countBadge.innerText = `${this.state.timeline.tracks.length} tracks`;

        if (this.state.timeline.tracks.length === 0) {
            container.innerHTML = `
                <div style="padding: 2rem; text-align: center; color: var(--text-dim); font-size: 0.78rem;">
                    No tracks yet. Choose an orbital above, set In/Out frames, and click <b>[➕ Add Clip]</b>.
                </div>
            `;
            return;
        }

        const total = this.state.totalFrames;

        container.innerHTML = this.state.timeline.tracks.map((trk, trkIdx) => {
            const clipsHTML = trk.clips.map(clip => {
                const leftPct = (clip.start / (total - 1)) * 100;
                const widthPct = Math.max(1.5, ((clip.end - clip.start) / (total - 1)) * 100);

                return `
                    <div class="timeline-clip-block" data-trkid="${trk.id}" data-clipid="${clip.id}" style="position: absolute; left: ${leftPct}%; width: ${widthPct}%; top: 4px; bottom: 4px; background: rgba(30, 41, 59, 0.9); border: 1px solid ${trk.color}; border-left: 3px solid ${trk.color}; border-radius: 4px; padding: 0.2rem 0.4rem; display: flex; justify-content: space-between; align-items: center; overflow: hidden; cursor: grab; box-shadow: 0 2px 6px rgba(0,0,0,0.4);" title="${clip.title} (Frames ${clip.start + 1} to ${clip.end + 1})">
                        <span style="font-size: 0.65rem; font-weight: 700; color: #f8fafc; text-overflow: ellipsis; white-space: nowrap; overflow: hidden;">
                            ${clip.title}
                        </span>
                        <button class="btn-del-clip" data-trkid="${trk.id}" data-clipid="${clip.id}" style="background: none; border: none; color: #94a3b8; font-size: 0.65rem; cursor: pointer; padding: 0 0.2rem;" title="Delete this clip">✕</button>
                    </div>
                `;
            }).join('');

            return `
                <div class="timeline-track-row" data-trkid="${trk.id}" style="display: grid; grid-template-columns: 200px 1fr; border-bottom: 1px solid rgba(255,255,255,0.06); height: 42px; background: ${trkIdx % 2 === 0 ? 'rgba(0,0,0,0.15)' : 'rgba(0,0,0,0.25)'};">
                    
                    <!-- Track Header Column -->
                    <div style="padding: 0.3rem 0.6rem; border-right: 1px solid rgba(255,255,255,0.08); display: flex; justify-content: space-between; align-items: center;">
                        <div style="display: flex; align-items: center; gap: 0.4rem; overflow: hidden;">
                            <span style="display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: ${trk.color}; flex-shrink: 0;"></span>
                            <div style="overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">
                                <div style="font-weight: 700; font-size: 0.74rem; color: #f8fafc; font-family: var(--font-mono);">${trk.orbital_key}</div>
                                <div style="font-size: 0.65rem; color: var(--text-muted); overflow: hidden; text-overflow: ellipsis;">${trk.title}</div>
                            </div>
                        </div>

                        <div style="display: flex; align-items: center; gap: 0.25rem;">
                            <!-- Default ON Toggle (Eye) -->
                            <button class="btn-toggle-default-on" data-trkid="${trk.id}" style="background: none; border: none; font-size: 0.8rem; cursor: pointer; opacity: ${trk.default_on ? 1.0 : 0.35}; color: ${trk.default_on ? '#38bdf8' : '#64748b'};" title="Default ON: Automatically visible in 3D when playback enters this clip">
                                👁️
                            </button>
                            <!-- Delete Track -->
                            <button class="btn-del-track" data-trkid="${trk.id}" style="background: none; border: none; font-size: 0.72rem; cursor: pointer; color: #ef4444; opacity: 0.7;" title="Delete Track">
                                🗑️
                            </button>
                        </div>
                    </div>

                    <!-- Track Lane Body -->
                    <div class="timeline-track-lane" data-trkid="${trk.id}" style="position: relative; width: 100%; height: 100%;">
                        ${clipsHTML}
                    </div>
                </div>
            `;
        }).join('');

        // Attach event handlers to clips and tracks
        container.querySelectorAll('.btn-del-clip').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const trkId = btn.dataset.trkid;
                const clipId = btn.dataset.clipid;
                const trk = this.state.timeline.tracks.find(t => t.id === trkId);
                if (trk) {
                    trk.clips = trk.clips.filter(c => c.id !== clipId);
                    this.renderTimelineTracks();
                    this.updateActiveLobesPills();
                    this.saveTimeline();
                }
            };
        });

        container.querySelectorAll('.btn-del-track').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const trkId = btn.dataset.trkid;
                this.state.timeline.tracks = this.state.timeline.tracks.filter(t => t.id !== trkId);
                this.renderTimelineTracks();
                this.updateActiveLobesPills();
                this.saveTimeline();
            };
        });

        container.querySelectorAll('.btn-toggle-default-on').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const trkId = btn.dataset.trkid;
                const trk = this.state.timeline.tracks.find(t => t.id === trkId);
                if (trk) {
                    trk.default_on = !trk.default_on;
                    this.renderTimelineTracks();
                    this.updateActiveLobesPills();
                    this.renderCurrent3DScene();
                    this.saveTimeline();
                }
            };
        });

        // Click on track lane seeks to that frame
        container.querySelectorAll('.timeline-track-lane').forEach(lane => {
            lane.onclick = (e) => {
                const rect = lane.getBoundingClientRect();
                const frac = Math.max(0, Math.min(1, (e.clientX - rect.left) / rect.width));
                const targetFrame = Math.round(frac * (this.state.totalFrames - 1));
                this.pause();
                this.loadFrame(targetFrame);
            };
        });
    },

    renderEnergyProfileSVG() {
        const svg = document.getElementById('ruler-energy-svg');
        if (!svg) return;
        const total = this.state.totalFrames;
        const ep = this.state.energyProfile;
        if (!ep || ep.length === 0 || total <= 1) return;

        const width = 1000;
        const height = 52;
        svg.setAttribute('viewBox', `0 0 ${width} ${height}`);
        svg.setAttribute('preserveAspectRatio', 'none');

        const minE = Math.min(...ep);
        const maxE = Math.max(...ep);
        const rangeE = Math.max(1.0, maxE - minE);

        // Build points
        const points = ep.map((e, idx) => {
            const x = (idx / (total - 1)) * width;
            // Invert y because SVG y goes downwards
            const y = height - 8 - ((e - minE) / rangeE) * (height - 18);
            return { x, y };
        });

        const pathD = `M 0,${height} ` + points.map(p => `L ${p.x.toFixed(1)},${p.y.toFixed(1)}`).join(' ') + ` L ${width},${height} Z`;
        const lineD = points.map((p, i) => `${i === 0 ? 'M' : 'L'} ${p.x.toFixed(1)},${p.y.toFixed(1)}`).join(' ');

        // Landmarks markers
        const landmarksHTML = this.state.landmarks.map(lm => {
            const x = (lm.frame_index / (total - 1)) * width;
            return `
                <line x1="${x}" y1="0" x2="${x}" y2="${height}" stroke="rgba(192, 132, 252, 0.45)" stroke-dasharray="3,3" stroke-width="1.5" />
                <text x="${x + 3}" y="12" fill="#c084fc" font-size="9" font-family="monospace" font-weight="bold">${lm.name}</text>
            `;
        }).join('');

        svg.innerHTML = `
            <defs>
                <linearGradient id="energyGrad" x1="0" y1="0" x2="0" y2="1">
                    <stop offset="0%" stop-color="#38bdf8" stop-opacity="0.35" />
                    <stop offset="100%" stop-color="#38bdf8" stop-opacity="0.02" />
                </linearGradient>
            </defs>
            <path d="${pathD}" fill="url(#energyGrad)" />
            <path d="${lineD}" fill="none" stroke="#38bdf8" stroke-width="2" />
            ${landmarksHTML}
        `;
    },

    async saveTimeline() {
        try {
            await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'save_timeline',
                    timeline: this.state.timeline
                })
            });
        } catch (err) {
            console.error("Save timeline error:", err);
        }
    },

    async runPrecomputeTimeline() {
        if (this.state.isPrecomputing) return;
        this.state.isPrecomputing = true;
        this.state.precomputeAbort = false;

        const btnPre = document.getElementById('btn-precompute-timeline');
        const btnStop = document.getElementById('btn-precompute-stop');
        const progWrap = document.getElementById('timeline-progress-wrap');
        const progFill = document.getElementById('timeline-progress-fill');
        const progStatus = document.getElementById('timeline-progress-status');
        const progPct = document.getElementById('timeline-progress-pct');

        if (btnPre) btnPre.disabled = true;
        if (btnStop) btnStop.style.display = 'inline-block';
        if (progWrap) progWrap.style.display = 'flex';

        // Collect all orbital keys from timeline
        const allKeys = Array.from(new Set(this.state.timeline.tracks.map(t => t.orbital_key)));
        if (allKeys.length === 0) allKeys.push('IBO_23', 'IBO_21', 'IBO_20', 'LUMO');

        const total = this.state.totalFrames;
        const batchSz = 6;
        let cur = 0;

        while (cur < total && !this.state.precomputeAbort) {
            const endF = Math.min(total - 1, cur + batchSz - 1);
            if (progStatus) progStatus.innerText = `Precomputing frames ${cur + 1}..${endF + 1} of ${total}...`;

            try {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'batch_precompute_cubes',
                        start_frame: cur,
                        batch_size: batchSz,
                        selected_tracks: allKeys
                    })
                });
                const d = await res.json();
                if (!d.success) break;

                if (d.cubes) {
                    Object.assign(this.state.cubeCache, d.cubes);
                }

                const pct = d.percent || Math.round(((endF + 1) / total) * 100);
                if (progFill) progFill.style.width = `${pct}%`;
                if (progPct) progPct.innerText = `${pct}%`;

                cur = endF + 1;
            } catch (e) {
                console.error("Batch error:", e);
                break;
            }
        }

        this.state.isPrecomputing = false;
        if (btnPre) btnPre.disabled = false;
        if (btnStop) btnStop.style.display = 'none';

        if (progStatus) {
            if (this.state.precomputeAbort) {
                progStatus.innerText = 'Precomputation paused.';
                progStatus.style.color = '#fbbf24';
            } else {
                progStatus.innerText = '✓ Timeline precomputed! Butter-smooth 60 FPS playback active.';
                progStatus.style.color = '#10b981';
            }
        }

        await this.renderCurrent3DScene();
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
        const btn = document.getElementById('btn-transport-play');
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
        const btn = document.getElementById('btn-transport-play');
        if (btn) btn.innerHTML = '▶';
    },

    async downloadCurrentCube() {
        const key = this.state.auditionOrbital || (this.state.timeline.tracks[0] ? this.state.timeline.tracks[0].orbital_key : 'IBO_23');
        const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                action_type: 'get_cube',
                frame_index: this.state.currentFrame,
                orbital: key,
                isovalue: this.state.isovalue
            })
        });
        const data = await res.json();
        if (data.success && data.cube) {
            const blob = new Blob([data.cube], { type: 'text/plain' });
            const url = URL.createObjectURL(blob);
            const a = document.createElement('a');
            a.href = url;
            a.download = data.filename || `${key}_frame${this.state.currentFrame + 1}.cube`;
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);
            URL.revokeObjectURL(url);
        }
    },

    async getData(uiSchema) {
        return {
            timeline: this.state.timeline
        };
    }
};
