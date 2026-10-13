export const ChooseJSONsUI = {
    state: {
        available: [],
        sequence: [],
        activeFilename: null,
        activeTraj: null,
        currentFrame: 0,
        isPlaying: false,
        playTimer: null,
        viewer: null,
        stepId: 'step_1_choose_jsons'
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        this.state.available = uiSchema.available_trajectories || [];
        this.state.sequence = [...(uiSchema.selected_sequence || [])];

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

        this.initEventListeners();
        this.renderQueue();
        this.renderCatalog();

        // Default preview to the first file in sequence or catalog
        const initialToPreview = this.state.sequence[0] || (this.state.available[0] ? this.state.available[0].filename : null);
        if (initialToPreview) {
            this.loadAndPreview(initialToPreview);
        }
    },

    initEventListeners() {
        // Clear Queue button
        const clearBtn = document.getElementById('btn-clear-queue');
        if (clearBtn) {
            clearBtn.onclick = () => {
                this.state.sequence = [];
                this.renderQueue();
                this.renderCatalog();
            };
        }

        // Upload Dropzone
        const dropzone = document.getElementById('upload-zone');
        const fileInput = document.getElementById('file-input-upload');
        if (dropzone && fileInput) {
            dropzone.onclick = () => fileInput.click();
            fileInput.onchange = (e) => {
                if (e.target.files && e.target.files[0]) {
                    this.handleFileUpload(e.target.files[0]);
                }
            };
            dropzone.ondragover = (e) => {
                e.preventDefault();
                dropzone.classList.add('dragover');
            };
            dropzone.ondragleave = () => dropzone.classList.remove('dragover');
            dropzone.ondrop = (e) => {
                e.preventDefault();
                dropzone.classList.remove('dragover');
                if (e.dataTransfer.files && e.dataTransfer.files[0]) {
                    this.handleFileUpload(e.dataTransfer.files[0]);
                }
            };
        }

        // Rescan button
        const rescanBtn = document.getElementById('btn-rescan-traj');
        if (rescanBtn) {
            rescanBtn.onclick = async () => {
                const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'rescan' })
                });
                const data = await res.json();
                if (data.available) {
                    this.state.available = data.available;
                    this.renderCatalog();
                }
            };
        }

        // Add to sequence button in preview header
        const addBtn = document.getElementById('btn-add-to-sequence');
        if (addBtn) {
            addBtn.onclick = () => {
                if (this.state.activeFilename && !this.state.sequence.includes(this.state.activeFilename)) {
                    this.state.sequence.push(this.state.activeFilename);
                    this.renderQueue();
                    this.renderCatalog();
                }
            };
        }

        // Playback Controls
        const playBtn = document.getElementById('btn-play-pause');
        const scrubber = document.getElementById('frame-scrubber');
        const prevBtn = document.getElementById('btn-prev-frame');
        const nextBtn = document.getElementById('btn-next-frame');
        const tsBtn = document.getElementById('btn-jump-ts');

        if (playBtn) playBtn.onclick = () => this.togglePlay();
        if (prevBtn) prevBtn.onclick = () => this.stepFrame(-1);
        if (nextBtn) nextBtn.onclick = () => this.stepFrame(1);
        if (tsBtn) {
            tsBtn.onclick = () => {
                if (this.state.activeTraj) {
                    const tsIdx = this.state.activeTraj.ts_frame_index || 0;
                    this.setFrame(tsIdx);
                }
            };
        }
        if (scrubber) {
            scrubber.oninput = (e) => {
                this.pause();
                this.setFrame(parseInt(e.target.value, 10));
            };
        }
    },

    renderQueue() {
        const queueList = document.getElementById('sequence-queue-list');
        const badge = document.getElementById('queue-count-badge');
        if (badge) badge.innerText = `${this.state.sequence.length} Steps`;

        if (!queueList) return;

        if (this.state.sequence.length === 0) {
            queueList.innerHTML = `<div class="empty-state" style="padding: 1.5rem; font-size: 0.8rem; color: var(--text-muted);">No trajectories selected yet. Click "+ Add" on any trajectory from the catalog below.</div>`;
            return;
        }

        queueList.innerHTML = this.state.sequence.map((filename, idx) => {
            const summary = this.state.available.find(a => a.filename === filename);
            const rxnTitle = summary ? summary.reaction_smiles : filename;
            return `
                <div class="queue-item" data-filename="${filename}">
                    <div style="display: flex; align-items: center; overflow: hidden; cursor: pointer;" class="queue-item-click">
                        <span class="queue-step-badge">Step ${idx + 1}</span>
                        <div class="queue-file-title" title="${rxnTitle}">${filename}</div>
                    </div>
                    <div class="queue-actions">
                        <button class="btn-icon" data-action="up" data-idx="${idx}" title="Move up" ${idx === 0 ? 'disabled' : ''}>▲</button>
                        <button class="btn-icon" data-action="down" data-idx="${idx}" title="Move down" ${idx === this.state.sequence.length - 1 ? 'disabled' : ''}>▼</button>
                        <button class="btn-icon danger" data-action="remove" data-idx="${idx}" title="Remove">✕</button>
                    </div>
                </div>
            `;
        }).join('');

        // Attach action handlers
        queueList.querySelectorAll('.queue-item-click').forEach(el => {
            el.onclick = () => {
                const fn = el.closest('.queue-item').dataset.filename;
                this.loadAndPreview(fn);
            };
        });

        queueList.querySelectorAll('button[data-action]').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const action = btn.dataset.action;
                const idx = parseInt(btn.dataset.idx, 10);
                if (action === 'up' && idx > 0) {
                    const temp = this.state.sequence[idx - 1];
                    this.state.sequence[idx - 1] = this.state.sequence[idx];
                    this.state.sequence[idx] = temp;
                } else if (action === 'down' && idx < this.state.sequence.length - 1) {
                    const temp = this.state.sequence[idx + 1];
                    this.state.sequence[idx + 1] = this.state.sequence[idx];
                    this.state.sequence[idx] = temp;
                } else if (action === 'remove') {
                    this.state.sequence.splice(idx, 1);
                }
                this.renderQueue();
                this.renderCatalog();
            };
        });
    },

    renderCatalog() {
        const catalogList = document.getElementById('catalog-list');
        if (!catalogList) return;

        if (this.state.available.length === 0) {
            catalogList.innerHTML = `<div class="empty-state" style="padding: 1.5rem; font-size: 0.8rem;">No trajectories found. Upload one to get started.</div>`;
            return;
        }

        catalogList.innerHTML = this.state.available.map(item => {
            const isSelected = item.filename === this.state.activeFilename;
            const inQueue = this.state.sequence.includes(item.filename);
            const queueIdx = inQueue ? this.state.sequence.indexOf(item.filename) + 1 : null;

            return `
                <div class="catalog-item ${isSelected ? 'selected-for-preview' : ''}" data-filename="${item.filename}">
                    <div class="catalog-header">
                        <span class="catalog-filename">${item.filename}</span>
                        <div style="display: flex; gap: 0.4rem; align-items: center;">
                            ${inQueue ? `<span class="chip blue" style="font-size: 0.7rem;">Step ${queueIdx}</span>` : ''}
                            <button class="btn btn-sm ${inQueue ? 'btn-ghost danger-hover' : 'btn-primary'} btn-toggle-queue" data-filename="${item.filename}" style="padding: 0.2rem 0.55rem; font-size: 0.725rem;">
                                ${inQueue ? '✕ Remove' : '+ Add'}
                            </button>
                        </div>
                    </div>
                    <div style="font-size: 0.775rem; color: #e2e8f0; font-family: var(--font-mono); overflow: hidden; text-overflow: ellipsis; white-space: nowrap;">
                        ${item.reaction_smiles}
                    </div>
                    <div class="catalog-chips">
                        <span class="chip">${item.n_atoms} atoms</span>
                        <span class="chip">${item.total_frames} frames</span>
                        <span class="chip amber">ΔE‡ fwd: ${item.forward_barrier_kcal.toFixed(1)} kcal/mol</span>
                        <span class="chip emerald">ΔE rxn: ${item.reaction_energy_kcal.toFixed(1)} kcal/mol</span>
                    </div>
                </div>
            `;
        }).join('');

        // Clicking card previews
        catalogList.querySelectorAll('.catalog-item').forEach(el => {
            el.onclick = () => {
                const fn = el.dataset.filename;
                this.loadAndPreview(fn);
            };
        });

        // Add/remove button handler
        catalogList.querySelectorAll('.btn-toggle-queue').forEach(btn => {
            btn.onclick = (e) => {
                e.stopPropagation();
                const fn = btn.dataset.filename;
                const idx = this.state.sequence.indexOf(fn);
                if (idx >= 0) {
                    this.state.sequence.splice(idx, 1);
                } else {
                    this.state.sequence.push(fn);
                }
                this.renderQueue();
                this.renderCatalog();
            };
        });
    },


    async loadAndPreview(filename) {
        this.pause();
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

            this.state.activeTraj = data.trajectory;
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

            // Details
            const details = document.getElementById('preview-details');
            const rText = document.getElementById('preview-reactants-text');
            const pText = document.getElementById('preview-products-text');
            if (details && traj.fragments) {
                details.style.display = 'grid';
                if (rText) rText.innerText = (traj.fragments.reactants || []).map(f => f.formula || f.smiles).join(' + ') || 'Reactants';
                if (pText) pText.innerText = (traj.fragments.products || []).map(f => f.formula || f.smiles).join(' + ') || 'Products';
            }

            // Enable playback bar
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

            this.init3DViewer();
            this.setFrame(0, true);

        } catch (err) {
            if (title) title.innerText = `Error: ${err.message}`;
        }
    },

    init3DViewer() {
        const container = document.getElementById('preview-viewport-3d');
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

    setFrame(index, resetCamera = false) {
        if (!this.state.activeTraj || !this.state.viewer) return;
        const frames = this.state.activeTraj.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.state.currentFrame = index;
        const frame = frames[index];
        const symbols = this.state.activeTraj.symbols;
        const coords = frame.coords;

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

        // Update UI counters
        const scrubber = document.getElementById('frame-scrubber');
        const counter = document.getElementById('frame-label');
        const badge = document.getElementById('preview-status-badge');

        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;

        if (badge) {
            const isTS = index === this.state.activeTraj.ts_frame_index;
            const energy = frame.energy_kcal !== undefined ? frame.energy_kcal.toFixed(2) : '0.00';
            badge.innerText = `Frame ${index + 1} / ${frames.length} | Energy: ${energy} kcal/mol ${isTS ? ' ⭐ [TRANSITION STATE]' : ''}`;
        }
    },

    stepFrame(delta) {
        this.pause();
        const frames = this.state.activeTraj ? this.state.activeTraj.frames : [];
        if (!frames || frames.length === 0) return;
        let target = this.state.currentFrame + delta;
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
        if (!this.state.activeTraj) return;
        const frames = this.state.activeTraj.frames || [];
        if (frames.length === 0) return;

        // If at the end, restart from the beginning
        if (this.state.currentFrame >= frames.length - 1) {
            this.setFrame(0);
        }

        this.state.isPlaying = true;
        const btn = document.getElementById('btn-play-pause');
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
        const btn = document.getElementById('btn-play-pause');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
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
