/**
 * High-Resolution Orbital Recomputation Manager & Progress Controller.
 * Allows recomputing chosen reaction tracks with finer 3D grids (48³, 56³, 64³).
 * Supports chunked batch processing, live progress bars, and instant 3D viewer cache update.
 */

export class BundleHighresManager {
    constructor(options = {}) {
        this.stepId = options.stepId || 'step_7_bundle';
        this.tracks = options.tracks || [];
        this.totalFrames = options.totalFrames || 147;
        this.landmarks = options.landmarks || [];
        this.onComplete = options.onComplete || (() => {});
        this.isRunning = false;
        this.abortRequested = false;
    }

    renderModal(container, tracks = null) {
        if (tracks && tracks.length > 0) {
            this.tracks = tracks;
        }
        let existingModal = document.getElementById('modal-highres-recompute');
        if (existingModal) existingModal.remove();

        let quantumTracks = this.tracks.filter(t => t.orbital_key && t.orbital_key !== 'track_esp_charge');
        if (quantumTracks.length === 0) {
            quantumTracks = [
                { orbital_key: 'HOMO', chem_key: 'HOMO', title: 'HOMO', color: '#06b6d4' },
                { orbital_key: 'LUMO', chem_key: 'LUMO', title: 'LUMO', color: '#a855f7' }
            ];
        }

        const modal = document.createElement('div');
        modal.id = 'modal-highres-recompute';
        modal.className = 'custom-modal-backdrop';
        modal.innerHTML = `
            <div class="custom-modal-dialog" style="max-width: 540px;">
                <div class="modal-header">
                    <div>
                        <h3 style="font-size: 1.05rem; font-weight: 700; color: #f8fafc; margin: 0;">Recompute Orbitals (High Resolution)</h3>
                        <p style="font-size: 0.78rem; color: #94a3b8; margin-top: 2px;">Evaluate wavefunctions on a dense 3D volumetric grid for smooth, publication-grade isosurfaces.</p>
                    </div>
                    <button class="btn-modal-close" id="btn-close-highres-modal">&times;</button>
                </div>

                <div class="modal-body" style="display: flex; flex-direction: column; gap: 1rem; padding: 1.25rem;">
                    <!-- Resolution Preset Picker -->
                    <div>
                        <label style="font-size: 0.8rem; font-weight: 600; color: #cbd5e1; display: block; margin-bottom: 0.4rem;">Target Grid Resolution</label>
                        <div style="display: grid; grid-template-columns: 1fr 1fr 1fr; gap: 0.5rem;">
                            <label class="res-preset-card">
                                <input type="radio" name="res-preset" value="48">
                                <div class="res-info">
                                    <span class="res-title">48³ Grid</span>
                                    <span class="res-sub">~110k pts (Fast)</span>
                                </div>
                            </label>
                            <label class="res-preset-card active">
                                <input type="radio" name="res-preset" value="56" checked>
                                <div class="res-info">
                                    <span class="res-title">56³ Grid</span>
                                    <span class="res-sub">~175k pts (Optimal)</span>
                                </div>
                            </label>
                            <label class="res-preset-card">
                                <input type="radio" name="res-preset" value="64">
                                <div class="res-info">
                                    <span class="res-title">64³ Grid</span>
                                    <span class="res-sub">~262k pts (Ultra)</span>
                                </div>
                            </label>
                        </div>
                    </div>

                    <!-- Target Tracks Selection -->
                    <div>
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
                            <label style="font-size: 0.8rem; font-weight: 600; color: #cbd5e1;">Target Reaction Tracks</label>
                            <span style="font-size: 0.72rem; color: #38bdf8; cursor: pointer;" id="btn-select-all-tracks">Select All</span>
                        </div>
                        <div style="max-height: 140px; overflow-y: auto; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.5rem; display: flex; flex-direction: column; gap: 0.35rem;">
                            ${quantumTracks.map(trk => `
                                <label style="display: flex; align-items: center; gap: 0.5rem; font-size: 0.76rem; color: #e2e8f0; cursor: pointer;">
                                    <input type="checkbox" class="chk-highres-track" value="${trk.orbital_key}" data-chem="${trk.chem_key || ''}" checked>
                                    <span style="color: ${trk.color || '#38bdf8'};">●</span>
                                    <span>${trk.title}</span>
                                </label>
                            `).join('')}
                        </div>
                    </div>

                    <!-- Frame Scope Selection -->
                    <div>
                        <label style="font-size: 0.8rem; font-weight: 600; color: #cbd5e1; display: block; margin-bottom: 0.4rem;">Frame Scope</label>
                        <select id="sel-highres-scope" class="input-select" style="width: 100%; font-size: 0.8rem;">
                            <option value="all" selected>All ${this.totalFrames} Pathway Frames (Full High-Res Video)</option>
                            <option value="landmarks">Key Landmark Frames Only (Minima & TS states)</option>
                            <option value="stride2">Every 2nd Frame (Half size / 2x speed)</option>
                        </select>
                    </div>

                    <!-- Progress Bar (Initially hidden) -->
                    <div id="highres-progress-container" style="display: none; background: rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.75rem;">
                        <div style="display: flex; justify-content: space-between; font-size: 0.75rem; margin-bottom: 0.35rem;">
                            <span id="highres-status-label" style="color: #38bdf8; font-weight: 600;">Evaluating wavefunction grid...</span>
                            <span id="highres-percent-label" style="color: #f8fafc; font-family: monospace;">0%</span>
                        </div>
                        <div style="width: 100%; height: 8px; background: rgba(255,255,255,0.1); border-radius: 9999px; overflow: hidden;">
                            <div id="highres-progress-bar" style="width: 0%; height: 100%; background: linear-gradient(90deg, #38bdf8, #2563eb); transition: width 0.2s ease;"></div>
                        </div>
                        <div id="highres-detail-label" style="font-size: 0.72rem; color: #94a3b8; margin-top: 0.4rem;">Processing frame 0...</div>
                    </div>
                </div>

                <div class="modal-footer" style="padding: 0.8rem 1.25rem; display: flex; justify-content: flex-end; gap: 0.6rem; background: rgba(0,0,0,0.2); border-top: 1px solid rgba(255,255,255,0.08);">
                    <button class="btn btn-secondary btn-sm" id="btn-cancel-highres">Cancel</button>
                    <button class="btn btn-primary btn-sm" id="btn-start-highres" style="background: linear-gradient(135deg, #0284c7, #2563eb);">Start Recomputation</button>
                </div>
            </div>
        `;

        document.body.appendChild(modal);

        // Bind preset radio cards
        modal.querySelectorAll('.res-preset-card').forEach(card => {
            card.addEventListener('click', () => {
                modal.querySelectorAll('.res-preset-card').forEach(c => c.classList.remove('active'));
                card.classList.add('active');
            });
        });

        const selectAllBtn = modal.querySelector('#btn-select-all-tracks');
        if (selectAllBtn) {
            selectAllBtn.onclick = () => {
                const checkboxes = modal.querySelectorAll('.chk-highres-track');
                const anyUnchecked = Array.from(checkboxes).some(c => !c.checked);
                checkboxes.forEach(c => c.checked = anyUnchecked);
                selectAllBtn.textContent = anyUnchecked ? 'Deselect All' : 'Select All';
            };
        }

        modal.querySelector('#btn-close-highres-modal').onclick = () => this.closeModal();
        modal.querySelector('#btn-cancel-highres').onclick = () => {
            if (this.isRunning) {
                this.abortRequested = true;
            } else {
                this.closeModal();
            }
        };

        modal.querySelector('#btn-start-highres').onclick = () => this.startRecomputation(modal);
    }

    closeModal() {
        const modal = document.getElementById('modal-highres-recompute');
        if (modal) modal.remove();
        this.isRunning = false;
        this.abortRequested = false;
    }

    async startRecomputation(modal) {
        const selectedRes = parseInt(modal.querySelector('input[name="res-preset"]:checked').value, 10);
        const selectedTracks = [];
        modal.querySelectorAll('.chk-highres-track:checked').forEach(chk => {
            selectedTracks.push({
                orbital_key: chk.value,
                chem_key: chk.dataset.chem || null
            });
        });

        if (selectedTracks.length === 0) {
            alert('Please select at least one track to recompute.');
            return;
        }

        const scope = modal.querySelector('#sel-highres-scope').value;
        const progressBox = modal.querySelector('#highres-progress-container');
        const progressBar = modal.querySelector('#highres-progress-bar');
        const statusLabel = modal.querySelector('#highres-status-label');
        const percentLabel = modal.querySelector('#highres-percent-label');
        const detailLabel = modal.querySelector('#highres-detail-label');
        const startBtn = modal.querySelector('#btn-start-highres');
        const cancelBtn = modal.querySelector('#btn-cancel-highres');

        progressBox.style.display = 'block';
        startBtn.disabled = true;
        cancelBtn.textContent = 'Stop Process';

        this.isRunning = true;
        this.abortRequested = false;

        const batchSize = 10;
        let startFrame = 0;
        const total = this.totalFrames;

        try {
            while (startFrame < total && !this.abortRequested) {
                const endFrame = Math.min(startFrame + batchSize - 1, total - 1);
                detailLabel.textContent = `Computing grid points for frames ${startFrame} to ${endFrame}...`;

                const res = await fetch(`/api/step/${this.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'recompute_highres_cubes',
                        resolution: selectedRes,
                        selected_tracks: selectedTracks,
                        start_frame: startFrame,
                        end_frame: endFrame
                    })
                });

                const data = await res.json();
                if (!data.success) {
                    throw new Error(data.error || 'Recomputation failed.');
                }

                const pct = Math.min(100, Math.round(((endFrame + 1) / total) * 100));
                progressBar.style.width = `${pct}%`;
                percentLabel.textContent = `${pct}%`;

                if (data.is_complete || endFrame >= total - 1) {
                    break;
                }
                startFrame = endFrame + 1;
            }

            if (this.abortRequested) {
                statusLabel.textContent = 'Process stopped by user.';
                statusLabel.style.color = '#f59e0b';
            } else {
                statusLabel.textContent = `✓ High-Res Grid (${selectedRes}³) Complete!`;
                statusLabel.style.color = '#10b981';
                progressBar.style.width = '100%';
                percentLabel.textContent = '100%';
                detailLabel.textContent = `Successfully evaluated all selected orbitals at ${selectedRes}³ resolution.`;

                setTimeout(() => {
                    this.closeModal();
                    this.onComplete(selectedRes);
                }, 1200);
            }
        } catch (err) {
            statusLabel.textContent = `Error: ${err.message}`;
            statusLabel.style.color = '#ef4444';
            startBtn.disabled = false;
        } finally {
            this.isRunning = false;
        }
    }
}
