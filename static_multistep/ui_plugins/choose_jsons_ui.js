import { Choose3DPreviewer } from './choose/choose_preview_3d.js';
import { renderCatalogList, renderSequenceQueue } from './choose/choose_queue_dom.js';

export const ChooseJSONsUI = {
    state: {
        available: [],
        sequence: [],
        activeFilename: null,
        player: null,
        stepId: 'step_1_choose_jsons'
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.available = uiSchema.available_trajectories || [];
        this.state.sequence = [...(uiSchema.selected_sequence || [])];
        this.state.player = new Choose3DPreviewer();

        container.innerHTML = `
            <div class="choose-container">
                <!-- Left Column: Pathway Sequence Queue + Catalog -->
                <div class="left-panel">
                    <!-- Sequence Assembly Queue -->
                    <div class="queue-card">
                        <div class="section-title">
                            <span>Pathway Sequence</span>
                            <div style="display: flex; gap: 0.5rem; align-items: center;">
                                <span class="chip blue" id="queue-count-badge">${this.state.sequence.length} Steps</span>
                                <button id="btn-clear-queue" class="btn btn-sm btn-ghost danger-hover" style="padding: 0.15rem 0.5rem; font-size: 0.725rem;" title="Clear entire sequence">Clear</button>
                            </div>
                        </div>
                        <div id="sequence-queue-list" class="queue-list"></div>
                    </div>

                    <!-- Available Catalog -->
                    <div class="catalog-card">
                        <div class="section-title">
                            <span>Available Trajectories</span>
                            <button id="btn-rescan-traj" class="btn-icon" title="Refresh files">
                                <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2">
                                    <path d="M21.5 2v6h-6M2.5 22v-6h6M2 11.5a10 10 0 0 1 18.8-4.3M22 12.5a10 10 0 0 1-18.8 4.2"/>
                                </svg>
                            </button>
                        </div>
                        <div id="catalog-list" class="catalog-list"></div>
                        
                        <!-- Upload Dropzone -->
                        <div id="upload-zone" class="upload-dropzone">
                            <input type="file" id="file-input-upload" accept=".json" style="display: none;">
                            <svg width="20" height="20" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" style="margin-bottom: 4px;">
                                <path d="M21 15v4a2 2 0 0 1-2 2H5a2 2 0 0 1-2-2v-4M17 8l-5-5-5 5M12 3v12"/>
                            </svg>
                            <div>Drop or <strong>browse</strong> custom IRC trajectory JSON</div>
                        </div>
                    </div>
                </div>

                <!-- Right Column: Rich Interactive 3D Preview -->
                <div class="preview-panel" id="preview-panel">
                    <div class="preview-header">
                        <div>
                            <div class="preview-rxn-title" id="preview-title">Select a trajectory to preview</div>
                            <div class="preview-meta-chips" id="preview-chips"></div>
                        </div>
                        <div id="preview-quick-actions">
                            <button id="btn-add-to-sequence" class="btn btn-sm btn-primary" style="display: none;">
                                + Add to Sequence
                            </button>
                        </div>
                    </div>

                    <!-- 3Dmol Viewport -->
                    <div class="viewer-wrapper">
                        <div id="preview-viewport-3d" class="viewport-3d"></div>
                        <div id="preview-status-badge" class="viewer-overlay-badge">No Trajectory Selected</div>
                    </div>

                    <!-- Playback Controls -->
                    <div class="playback-bar" id="playback-bar" style="opacity: 0.5; pointer-events: none;">
                        <button id="btn-play-pause" class="btn btn-sm btn-secondary" title="Play / Pause">
                            <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                                <polygon points="5 3 19 12 5 21 5 3"/>
                            </svg>
                        </button>
                        <button id="btn-prev-frame" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                        <input type="range" id="frame-scrubber" class="scrub-slider" min="0" max="0" value="0">
                        <button id="btn-next-frame" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                        <button id="btn-jump-ts" class="btn btn-sm btn-secondary" title="Jump to Transition State">TS</button>
                        <div id="frame-label" class="frame-counter">0 / 0</div>
                    </div>

                    <!-- Details Accordion -->
                    <div class="preview-details" id="preview-details" style="display: none;">
                        <div class="detail-box">
                            <div class="detail-box-title">Reactants</div>
                            <div id="preview-reactants-text" style="font-family: var(--font-mono); color: #cbd5e1;">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Products</div>
                            <div id="preview-products-text" style="font-family: var(--font-mono); color: #cbd5e1;">-</div>
                        </div>
                    </div>
                </div>
            </div>
        `;

        this.attachGlobalEvents();
        this.renderQueue();
        this.renderCatalog();

        if (this.state.available.length > 0) {
            const first = this.state.sequence[0] || this.state.available[0].filename;
            this.loadAndPreview(first);
        }
    },

    attachGlobalEvents() {
        const clearBtn = document.getElementById('btn-clear-queue');
        if (clearBtn) {
            clearBtn.onclick = () => {
                this.state.sequence = [];
                this.renderQueue();
                this.renderCatalog();
            };
        }

        const rescanBtn = document.getElementById('btn-rescan-traj');
        if (rescanBtn) {
            rescanBtn.onclick = async () => {
                rescanBtn.classList.add('rotating');
                try {
                    const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ action_type: 'rescan' })
                    });
                    const data = await res.json();
                    if (data.success) {
                        this.state.available = data.available;
                        this.renderCatalog();
                    }
                } finally {
                    rescanBtn.classList.remove('rotating');
                }
            };
        }

        const uploadZone = document.getElementById('upload-zone');
        const fileInput = document.getElementById('file-input-upload');
        if (uploadZone && fileInput) {
            uploadZone.onclick = () => fileInput.click();
            fileInput.onchange = (e) => {
                if (e.target.files.length > 0) this.handleFileUpload(e.target.files[0]);
            };
            uploadZone.ondragover = (e) => { e.preventDefault(); uploadZone.classList.add('dragover'); };
            uploadZone.ondragleave = () => uploadZone.classList.remove('dragover');
            uploadZone.ondrop = (e) => {
                e.preventDefault();
                uploadZone.classList.remove('dragover');
                if (e.dataTransfer.files.length > 0) this.handleFileUpload(e.dataTransfer.files[0]);
            };
        }

        const addBtn = document.getElementById('btn-add-to-sequence');
        if (addBtn) {
            addBtn.onclick = () => {
                if (this.state.activeFilename && !this.state.sequence.includes(this.state.activeFilename)) {
                    this.state.sequence.push(this.state.activeFilename);
                    this.renderQueue();
                    this.renderCatalog();
                    addBtn.style.display = 'none';
                }
            };
        }

        const p = this.state.player;
        const playBtn = document.getElementById('btn-play-pause');
        const prevBtn = document.getElementById('btn-prev-frame');
        const nextBtn = document.getElementById('btn-next-frame');
        const tsBtn = document.getElementById('btn-jump-ts');
        const scrubber = document.getElementById('frame-scrubber');

        if (playBtn) playBtn.onclick = () => p.togglePlay();
        if (prevBtn) prevBtn.onclick = () => p.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => p.stepFrame(1);
        if (tsBtn) tsBtn.onclick = () => {
            if (p.activeTraj && p.activeTraj.ts_frame_index !== undefined) {
                p.pause();
                p.setFrame(p.activeTraj.ts_frame_index);
            }
        };
        if (scrubber) {
            scrubber.oninput = (e) => {
                p.pause();
                p.setFrame(parseInt(e.target.value, 10));
            };
        }
    },

    renderQueue() {
        const queueList = document.getElementById('sequence-queue-list');
        const badge = document.getElementById('queue-count-badge');
        renderSequenceQueue(queueList, badge, this.state, {
            onSelect: (fn) => this.loadAndPreview(fn),
            onChanged: () => {
                this.renderQueue();
                this.renderCatalog();
            }
        });
    },

    renderCatalog() {
        const catalogList = document.getElementById('catalog-list');
        renderCatalogList(catalogList, this.state, {
            onSelect: (fn) => this.loadAndPreview(fn),
            onChanged: () => {
                this.renderQueue();
                this.renderCatalog();
            }
        });
    },

    async loadAndPreview(filename) {
        this.state.player.pause();
        this.state.activeFilename = filename;
        this.renderCatalog();

        const title = document.getElementById('preview-title');
        const chips = document.getElementById('preview-chips');
        const addBtn = document.getElementById('btn-add-to-sequence');
        const statusBadge = document.getElementById('preview-status-badge');

        if (title) title.innerText = `Loading ${filename}...`;
        if (statusBadge) statusBadge.innerText = 'Fetching trajectory coordinates...';

        try {
            const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ action_type: 'get_trajectory', filename: filename })
            });
            const data = await res.json();
            if (!data.success || !data.trajectory) {
                if (title) title.innerText = `Error loading ${filename}`;
                if (statusBadge) statusBadge.innerText = data.error || 'Failed to load';
                return;
            }

            const traj = data.trajectory;
            const frames = traj.frames || [];
            const barriers = traj.barriers || {};

            if (title) title.innerText = traj.reaction_smiles || filename;
            if (chips) {
                chips.innerHTML = `
                    <span class="chip blue">${traj.n_atoms || (traj.symbols ? traj.symbols.length : 0)} atoms</span>
                    <span class="chip">${frames.length} frames</span>
                    <span class="chip amber">Forward Barrier: ${(barriers.forward_barrier_kcal || 0).toFixed(2)} kcal/mol</span>
                    <span class="chip emerald">Reaction Energy: ${(barriers.reaction_energy_kcal || 0).toFixed(2)} kcal/mol</span>
                `;
            }

            if (addBtn) {
                const inQueue = this.state.sequence.includes(filename);
                addBtn.style.display = inQueue ? 'none' : 'inline-flex';
            }

            const details = document.getElementById('preview-details');
            const rText = document.getElementById('preview-reactants-text');
            const pText = document.getElementById('preview-products-text');
            if (details && traj.fragments) {
                details.style.display = 'grid';
                if (rText) rText.innerText = (traj.fragments.reactants || []).map(f => f.formula || f.smiles).join(' + ') || 'Reactants';
                if (pText) pText.innerText = (traj.fragments.products || []).map(f => f.formula || f.smiles).join(' + ') || 'Products';
            }

            const pbar = document.getElementById('playback-bar');
            const scrubber = document.getElementById('frame-scrubber');
            if (pbar) {
                pbar.style.opacity = '1';
                pbar.style.pointerEvents = 'auto';
            }
            if (scrubber) {
                scrubber.max = Math.max(0, frames.length - 1);
                scrubber.value = 0;
            }

            this.state.player.init(document.getElementById('preview-viewport-3d'));
            this.state.player.setTrajectory(traj);

        } catch (err) {
            if (title) title.innerText = `Error: ${err.message}`;
        }
    },

    async handleFileUpload(file) {
        try {
            const text = await file.text();
            const data = JSON.parse(text);

            const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'upload',
                    filename: file.name,
                    content: data
                })
            });

            const respData = await res.json();
            if (respData.success) {
                this.state.available = respData.available;
                this.renderCatalog();
                this.loadAndPreview(respData.summary.filename);
            } else {
                alert("Upload failed: " + (respData.error || "Invalid file"));
            }
        } catch (e) {
            alert("Error parsing JSON: " + e.message);
        }
    },

    async getData(uiSchema) {
        if (this.state.sequence.length < 2) {
            throw new Error("Please select at least 2 consecutive reaction trajectories in the Pathway Sequence before clicking Next.");
        }
        return {
            sequence: this.state.sequence
        };
    }
};
