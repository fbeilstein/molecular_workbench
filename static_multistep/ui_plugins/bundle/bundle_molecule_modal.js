/**
 * Molecule Orbital Studio Modal & New Window Viewer.
 * Triggered by right-clicking a molecule in the interactive SVG mechanism.
 * Allows uploading and unpacking Molecular Workbench ZIP bundles,
 * listing all discovery orbitals, and rendering 3D isosurfaces.
 */

export class BundleMoleculeModal {
    constructor(options = {}) {
        this.stepId = options.stepId || 'step_5_bundle_viewer';
        this.moleculeId = options.moleculeId || 'molecule';
        this.moleculeName = options.moleculeName || 'Molecule';
        this.onBundleUploaded = options.onBundleUploaded || (() => {});
        this.viewer = null;
        this.activeOrbital = null;
        this.isovalue = 0.03;
        this.bundleData = null;
    }

    async show(container, moleculeId, moleculeName, existingBundle = null) {
        this.moleculeId = moleculeId;
        this.moleculeName = moleculeName || moleculeId;
        this.bundleData = existingBundle;

        let existingModal = document.getElementById('modal-molecule-orbitals');
        if (existingModal) existingModal.remove();

        const modal = document.createElement('div');
        modal.id = 'modal-molecule-orbitals';
        modal.className = 'custom-modal-backdrop';
        modal.innerHTML = `
            <div class="custom-modal-dialog" style="max-width: 880px; width: 90vw; height: 560px; display: flex; flex-direction: column;">
                <div class="modal-header">
                    <div>
                        <h3 style="font-size: 1.05rem; font-weight: 700; color: #f8fafc; margin: 0; display: flex; align-items: center; gap: 0.5rem;">
                            <span>⚛️ ${this.moleculeName}</span>
                            <span class="chip blue" style="font-size: 0.72rem;">Molecular Workbench Bundle</span>
                        </h3>
                        <p style="font-size: 0.78rem; color: #94a3b8; margin-top: 2px;">Inspect all quantum molecular orbitals from uploaded Molecular Workbench archive.</p>
                    </div>
                    <div style="display: flex; gap: 0.5rem; align-items: center;">
                        <button class="btn btn-xs btn-secondary" id="btn-popout-window" title="Open in new window">⧉ New Window</button>
                        <button class="btn-modal-close" id="btn-close-mol-modal">&times;</button>
                    </div>
                </div>

                <div class="modal-body" style="flex: 1; min-height: 0; display: flex; gap: 0.75rem; padding: 0.8rem;">
                    <!-- Left Sidebar: Orbital List or Upload Zone -->
                    <div style="width: 280px; height: 100%; display: flex; flex-direction: column; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.6rem; overflow: hidden;">
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem;">
                            <span style="font-size: 0.82rem; font-weight: 600; color: #f1f5f9;">Molecular Orbitals</span>
                            <label class="btn btn-xs btn-secondary" style="cursor: pointer; font-size: 0.72rem; padding: 2px 6px;">
                                ⬆ Upload .zip
                                <input type="file" id="input-mol-zip" accept=".zip" style="display: none;">
                            </label>
                        </div>

                        <!-- Orbitals List or Upload Prompt -->
                        <div id="mol-orbitals-list" style="flex: 1; min-height: 0; overflow-y: auto; display: flex; flex-direction: column; gap: 0.35rem; padding-right: 4px;">
                            <div class="empty-state" style="padding: 1.5rem 0.5rem; font-size: 0.75rem;">
                                No bundle uploaded yet for this molecule.<br>Upload a Molecular Workbench .zip bundle above.
                            </div>
                        </div>

                        <!-- Isovalue Controller -->
                        <div style="border-top: 1px solid rgba(255,255,255,0.08); padding-top: 0.5rem; display: flex; flex-direction: column; gap: 0.3rem;">
                            <div style="display: flex; justify-content: space-between; font-size: 0.72rem; color: #94a3b8;">
                                <span>Isovalue</span>
                                <span id="mol-isoval-readout">0.030</span>
                            </div>
                            <input type="range" id="mol-isoval-slider" min="0.005" max="0.10" step="0.005" value="0.03" style="width: 100%; cursor: pointer; accent-color: #38bdf8;">
                        </div>
                    </div>

                    <!-- Right Side: 3Dmol Interactive Molecule & Orbital Viewer -->
                    <div style="flex: 1; height: 100%; position: relative; background: #090d16; border-radius: 6px; border: 1px solid rgba(255,255,255,0.08); overflow: hidden;">
                        <div id="mol-3d-canvas" style="width: 100%; height: 100%;"></div>
                        <div id="mol-status-badge" style="position: absolute; bottom: 8px; left: 8px; font-size: 0.72rem; color: #94a3b8; background: rgba(0,0,0,0.7); padding: 2px 6px; border-radius: 4px;">
                            Ready
                        </div>
                        <div style="position: absolute; top: 8px; right: 8px; display: flex; gap: 4px;">
                            <button id="btn-mol-zoom-in" class="btn btn-xs btn-secondary">+</button>
                            <button id="btn-mol-zoom-out" class="btn btn-xs btn-secondary">−</button>
                            <button id="btn-mol-fit" class="btn btn-xs btn-secondary">⟲</button>
                        </div>
                    </div>
                </div>
            </div>
        `;

        document.body.appendChild(modal);

        modal.querySelector('#btn-close-mol-modal').onclick = () => modal.remove();
        modal.querySelector('#btn-popout-window').onclick = () => this.openPopoutWindow();

        const fileInput = modal.querySelector('#input-mol-zip');
        fileInput.addEventListener('change', (e) => this.handleFileUpload(e.target.files[0]));

        const slider = modal.querySelector('#mol-isoval-slider');
        slider.addEventListener('input', (e) => {
            this.isovalue = parseFloat(e.target.value);
            modal.querySelector('#mol-isoval-readout').textContent = this.isovalue.toFixed(3);
            if (this.activeOrbital) this.renderActiveOrbital();
        });

        modal.querySelector('#btn-mol-zoom-in').onclick = () => this.viewer?.zoom(1.2);
        modal.querySelector('#btn-mol-zoom-out').onclick = () => this.viewer?.zoom(0.83);
        modal.querySelector('#btn-mol-fit').onclick = () => { this.viewer?.zoomTo(); this.viewer?.render(); };

        this.init3DViewer(modal.querySelector('#mol-3d-canvas'));

        // If existing bundle was provided or can be fetched from backend
        if (this.bundleData && this.bundleData.orbitals) {
            this.populateOrbitals(this.bundleData.orbitals);
        } else {
            this.fetchExistingBundle();
        }
    }

    init3DViewer(containerEl) {
        if (!containerEl || typeof $3Dmol === 'undefined') return;
        this.viewer = $3Dmol.createViewer(containerEl, {
            backgroundColor: '#090d16',
            id: 'mol-3dmol-inner'
        });
    }

    async fetchExistingBundle() {
        try {
            const res = await fetch(`/api/step/${this.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'get_molecule_bundle',
                    molecule_id: this.moleculeId
                })
            });
            const data = await res.json();
            if (data.success && data.orbitals && data.orbitals.length > 0) {
                this.bundleData = data;
                this.populateOrbitals(data.orbitals);
            }
        } catch (err) {}
    }

    async handleFileUpload(file) {
        if (!file) return;
        const reader = new FileReader();
        reader.onload = async () => {
            const base64Str = reader.result.split(',')[1];
            const modal = document.getElementById('modal-molecule-orbitals');
            const statusBadge = modal ? modal.querySelector('#mol-status-badge') : null;
            if (statusBadge) statusBadge.textContent = 'Uploading & parsing bundle...';

            try {
                const res = await fetch(`/api/step/${this.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'upload_molecule_bundle',
                        molecule_id: this.moleculeId,
                        filename: file.name,
                        data_b64: base64Str
                    })
                });
                const data = await res.json();
                if (data.success) {
                    this.bundleData = data;
                    this.populateOrbitals(data.orbitals);
                    this.onBundleUploaded(this.moleculeId, data);
                    if (statusBadge) statusBadge.textContent = `Parsed ${data.count} orbitals from ${file.name}`;
                } else {
                    alert('Error: ' + (data.error || 'Failed to parse bundle.'));
                }
            } catch (err) {
                alert('Upload failed: ' + err.message);
            }
        };
        reader.readAsDataURL(file);
    }

    populateOrbitals(orbitalNames) {
        const modal = document.getElementById('modal-molecule-orbitals');
        if (!modal) return;
        const listEl = modal.querySelector('#mol-orbitals-list');
        if (!listEl) return;

        if (!orbitalNames || orbitalNames.length === 0) {
            listEl.innerHTML = `<div class="empty-state">No orbital cubes found in this bundle.</div>`;
            return;
        }

        let html = '';
        orbitalNames.forEach((name, idx) => {
            html += `
                <div class="mol-orb-item ${idx === 0 ? 'active' : ''}" data-name="${name}" style="padding: 0.4rem 0.6rem; background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.08); border-radius: 4px; font-size: 0.75rem; cursor: pointer; display: flex; justify-content: space-between; align-items: center; transition: all 0.15s ease;">
                    <div style="display: flex; align-items: center; gap: 0.4rem;">
                        <span style="color: #38bdf8;">●</span>
                        <span style="font-weight: 600; color: #f1f5f9;">${name}</span>
                    </div>
                    <span style="font-size: 0.68rem; color: #64748b;">.cube</span>
                </div>
            `;
        });
        listEl.innerHTML = html;

        listEl.querySelectorAll('.mol-orb-item').forEach(item => {
            item.addEventListener('click', () => {
                listEl.querySelectorAll('.mol-orb-item').forEach(i => i.classList.remove('active'));
                item.classList.add('active');
                this.selectOrbital(item.dataset.name);
            });
        });

        // Select first orbital automatically
        if (orbitalNames.length > 0) {
            this.selectOrbital(orbitalNames[0]);
        }
    }

    async selectOrbital(orbName) {
        this.activeOrbital = orbName;
        this.renderActiveOrbital();
    }

    async renderActiveOrbital() {
        if (!this.viewer || !this.activeOrbital || !this.bundleData) return;
        const cubeText = this.bundleData.files ? this.bundleData.files[this.activeOrbital] : null;
        if (!cubeText) return;

        this.viewer.clear();

        try {
            // Volume data from cube
            const vdo = new $3Dmol.VolumeData(cubeText, 'cube');

            // Render atoms from cube
            this.viewer.addModel(cubeText, 'cube');
            this.viewer.setStyle({}, {
                stick: { radius: 0.14, colorscheme: 'Jmol' },
                sphere: { scale: 0.28, colorscheme: 'Jmol' }
            });

            // Positive Lobe (Cyan)
            this.viewer.addIsosurface(vdo, {
                isoval: Math.abs(this.isovalue),
                color: '#38bdf8',
                opacity: 0.75,
                wireframe: false
            });

            // Negative Lobe (Amber)
            this.viewer.addIsosurface(vdo, {
                isoval: -Math.abs(this.isovalue),
                color: '#f59e0b',
                opacity: 0.75,
                wireframe: false
            });

            this.viewer.zoomTo();
            this.viewer.render();

            const modal = document.getElementById('modal-molecule-orbitals');
            const statusBadge = modal ? modal.querySelector('#mol-status-badge') : null;
            if (statusBadge) statusBadge.textContent = `${this.activeOrbital} (isovalue: ±${this.isovalue.toFixed(3)})`;
        } catch (err) {
            console.error('Error rendering molecule orbital:', err);
        }
    }

    openPopoutWindow() {
        const title = `Orbitals - ${this.moleculeName}`;
        const win = window.open('', '_blank', 'width=900,height=600');
        if (!win) {
            alert('Popup was blocked by your browser. Please allow popups for this site.');
            return;
        }

        const orbNames = (this.bundleData && this.bundleData.orbitals) ? JSON.stringify(this.bundleData.orbitals) : '[]';
        const filesData = (this.bundleData && this.bundleData.files) ? JSON.stringify(this.bundleData.files) : '{}';

        win.document.write(`<!DOCTYPE html>
<html>
<head>
    <title>${title}</title>
    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <style>
        body { margin: 0; background: #0b0f17; color: #fff; font-family: sans-serif; display: flex; height: 100vh; }
        #sidebar { width: 260px; background: #111827; border-right: 1px solid rgba(255,255,255,0.1); padding: 12px; display: flex; flex-direction: column; gap: 8px; }
        #viewer { flex: 1; position: relative; }
        .item { padding: 6px 10px; background: rgba(255,255,255,0.05); border-radius: 4px; cursor: pointer; font-size: 13px; }
        .item.active { background: #2563eb; font-weight: bold; }
    </style>
</head>
<body>
    <div id="sidebar">
        <h3 style="margin: 0; font-size: 14px;">${this.moleculeName} Orbitals</h3>
        <div id="list" style="flex:1; overflow-y:auto; display:flex; flex-direction:column; gap:4px;"></div>
    </div>
    <div id="viewer"></div>
    <script>
        const orbs = ${orbNames};
        const files = ${filesData};
        const viewer = $3Dmol.createViewer(document.getElementById('viewer'), { backgroundColor: '#090d16' });
        const listEl = document.getElementById('list');

        function renderOrb(name) {
            const cube = files[name];
            if (!cube) return;
            viewer.clear();
            const vdo = new $3Dmol.VolumeData(cube, 'cube');
            viewer.addModel(cube, 'cube');
            viewer.setStyle({}, { stick: { radius: 0.14 }, sphere: { scale: 0.28 } });
            viewer.addIsosurface(vdo, { isoval: 0.03, color: '#38bdf8', opacity: 0.75 });
            viewer.addIsosurface(vdo, { isoval: -0.03, color: '#f59e0b', opacity: 0.75 });
            viewer.zoomTo();
            viewer.render();
        }

        orbs.forEach((name, i) => {
            const div = document.createElement('div');
            div.className = 'item ' + (i === 0 ? 'active' : '');
            div.textContent = name;
            div.onclick = () => {
                document.querySelectorAll('.item').forEach(d => d.classList.remove('active'));
                div.classList.add('active');
                renderOrb(name);
            };
            listEl.appendChild(div);
        });

        if (orbs.length > 0) renderOrb(orbs[0]);
    </script>
</body>
</html>`);
    }
}
