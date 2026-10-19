/**
 * Molecule Orbital Studio Modal & New Window Viewer.
 * Triggered by right-clicking a molecule in the interactive SVG mechanism.
 * Allows uploading and unpacking Molecular Workbench ZIP bundles,
 * listing all discovery orbitals, and rendering 3D isosurfaces & ESP meshes.
 */

import { openMoleculePopout } from './bundle_molecule_popout.js';

export class BundleMoleculeModal {
    constructor(options = {}) {
        this.stepId = options.stepId || 'step_7_bundle';
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

        document.getElementById('modal-molecule-orbitals')?.remove();

        const modal = document.createElement('div');
        modal.id = 'modal-molecule-orbitals';
        modal.className = 'custom-modal-backdrop';
        modal.innerHTML = `
            <div class="custom-modal-dialog" style="max-width: 900px; width: 92vw; height: 580px; display: flex; flex-direction: column;">
                <div class="modal-header">
                    <div>
                        <h3 style="font-size: 1.05rem; font-weight: 700; color: #f8fafc; margin: 0; display: flex; align-items: center; gap: 0.5rem;">
                            <span>⚛️ ${this.moleculeName}</span>
                            <span class="chip blue" style="font-size: 0.72rem;">Molecular Workbench Bundle</span>
                        </h3>
                        <p style="font-size: 0.78rem; color: #94a3b8; margin-top: 2px;">Quantum orbitals & electrostatic surfaces from bundle.</p>
                    </div>
                    <div style="display: flex; gap: 0.5rem; align-items: center;">
                        <button class="btn btn-xs btn-secondary" id="btn-popout-window" title="Open in new window">⧉ New Window</button>
                        <button class="btn-modal-close" id="btn-close-mol-modal">&times;</button>
                    </div>
                </div>

                <div class="modal-body" style="flex: 1; min-height: 0; display: flex; gap: 0.75rem; padding: 0.8rem;">
                    <div style="width: 290px; height: 100%; display: flex; flex-direction: column; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.6rem; overflow: hidden;">
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.5rem;">
                            <span style="font-size: 0.82rem; font-weight: 600; color: #f1f5f9;">Orbital Catalog</span>
                            <label class="btn btn-xs btn-secondary" style="cursor: pointer; font-size: 0.72rem; padding: 2px 6px;">
                                ⬆ Upload .zip
                                <input type="file" id="input-mol-zip" accept=".zip" style="display: none;">
                            </label>
                        </div>

                        <div id="mol-orbitals-list" style="flex: 1; min-height: 0; overflow-y: auto; display: flex; flex-direction: column; gap: 0.35rem; padding-right: 4px;">
                            <div class="empty-state" style="padding: 1.5rem 0.5rem; font-size: 0.75rem;">
                                No bundle uploaded yet.<br>Upload a Molecular Workbench .zip bundle above.
                            </div>
                        </div>

                        <div style="border-top: 1px solid rgba(255,255,255,0.08); padding-top: 0.5rem; display: flex; flex-direction: column; gap: 0.3rem;">
                            <div style="display: flex; justify-content: space-between; font-size: 0.72rem; color: #94a3b8;">
                                <span>Isovalue</span>
                                <span id="mol-isoval-readout">0.030</span>
                            </div>
                            <input type="range" id="mol-isoval-slider" min="0.005" max="0.10" step="0.005" value="0.03" style="width: 100%; cursor: pointer; accent-color: #38bdf8;">
                        </div>
                    </div>

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
        modal.querySelector('#btn-popout-window').onclick = () => openMoleculePopout(this.moleculeName, this.bundleData);

        modal.querySelector('#input-mol-zip').addEventListener('change', (e) => this.handleFileUpload(e.target.files[0]));

        modal.querySelector('#mol-isoval-slider').addEventListener('input', (e) => {
            this.isovalue = parseFloat(e.target.value);
            modal.querySelector('#mol-isoval-readout').textContent = this.isovalue.toFixed(3);
            if (this.activeOrbital) this.renderActiveOrbital();
        });

        modal.querySelector('#btn-mol-zoom-in').onclick = () => this.viewer?.zoom(1.2);
        modal.querySelector('#btn-mol-zoom-out').onclick = () => this.viewer?.zoom(0.83);
        modal.querySelector('#btn-mol-fit').onclick = () => { this.viewer?.zoomTo(); this.viewer?.render(); };

        this.init3DViewer(modal.querySelector('#mol-3d-canvas'));

        if (this.bundleData && (this.bundleData.orbitals || this.bundleData.xyz || this.bundleData.esp_surface)) {
            this.populateOrbitals(this.bundleData);
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
                    molecule_id: this.moleculeId,
                    molecule_name: this.moleculeName
                })
            });
            const data = await res.json();
            if (data.success) {
                this.bundleData = data;
                this.populateOrbitals(data);
            }
        } catch (err) {}
    }

    async handleFileUpload(file) {
        if (!file) return;
        const reader = new FileReader();
        reader.onload = async () => {
            const base64Str = reader.result.split(',')[1];
            const statusBadge = document.querySelector('#modal-molecule-orbitals #mol-status-badge');
            if (statusBadge) statusBadge.textContent = 'Uploading & parsing bundle...';

            try {
                const res = await fetch(`/api/step/${this.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'upload_molecule_bundle',
                        molecule_id: this.moleculeId,
                        molecule_name: this.moleculeName,
                        filename: file.name,
                        data_b64: base64Str
                    })
                });
                const data = await res.json();
                if (data.success) {
                    this.bundleData = data;
                    this.populateOrbitals(data);
                    this.onBundleUploaded(this.moleculeId, data, this.moleculeName);
                    if (statusBadge) statusBadge.textContent = `Loaded ${data.count || 0} orbitals from ${file.name}`;
                } else {
                    alert('Error: ' + (data.error || 'Failed to parse bundle.'));
                }
            } catch (err) {
                alert('Upload failed: ' + err.message);
            }
        };
        reader.readAsDataURL(file);
    }

    populateOrbitals(data) {
        const listEl = document.querySelector('#modal-molecule-orbitals #mol-orbitals-list');
        if (!listEl) return;

        const tree = data.orbitals_tree || [];
        const flatList = data.orbitals || [];
        const hasEsp = !!data.esp_surface;
        let html = '';

        if (hasEsp) {
            html += `
                <div class="mol-orb-item" id="btn-mol-esp" style="padding: 0.4rem 0.6rem; background: rgba(56,189,248,0.12); border: 1px solid rgba(56,189,248,0.3); border-radius: 4px; font-size: 0.75rem; cursor: pointer; display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem; color: #38bdf8; font-weight: 600;">
                    <div style="display: flex; align-items: center; gap: 0.4rem;">
                        <span>⚡</span>
                        <span>ESP Surface</span>
                    </div>
                    <span style="font-size: 0.68rem; color: #94a3b8;">±kcal/mol</span>
                </div>
            `;
        }

        if (tree.length > 0) {
            tree.forEach(cat => {
                html += `<div style="font-size: 0.68rem; font-weight: 700; color: #94a3b8; text-transform: uppercase; margin: 0.4rem 0 0.2rem 0; letter-spacing: 0.5px;">${cat.name}</div>`;
                (cat.items || []).forEach(item => {
                    const base = item.file ? item.file.replace(/^molecules\//, '').replace(/\.json$/, '').replace(/\.cube$/, '') : item.label;
                    const evStr = (item.energy_ev !== undefined && item.energy_ev !== null) ? `${item.energy_ev} eV` : '';
                    html += `
                        <div class="mol-orb-item" data-name="${base}" data-label="${item.label}" style="padding: 0.35rem 0.55rem; background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.08); border-radius: 4px; font-size: 0.74rem; cursor: pointer; display: flex; justify-content: space-between; align-items: center;">
                            <div style="display: flex; align-items: center; gap: 0.4rem;">
                                <span style="color: ${item.color || '#38bdf8'}; font-size: 0.8rem;">●</span>
                                <span style="font-weight: 600; color: #f1f5f9;">${item.label}</span>
                            </div>
                            <span style="font-size: 0.66rem; color: #94a3b8;">${evStr}</span>
                        </div>
                    `;
                });
            });
        } else if (flatList.length > 0) {
            flatList.forEach((name, idx) => {
                html += `
                    <div class="mol-orb-item ${idx === 0 ? 'active' : ''}" data-name="${name}" style="padding: 0.4rem 0.6rem; background: rgba(255,255,255,0.04); border: 1px solid rgba(255,255,255,0.08); border-radius: 4px; font-size: 0.75rem; cursor: pointer; display: flex; justify-content: space-between; align-items: center;">
                        <div style="display: flex; align-items: center; gap: 0.4rem;">
                            <span style="color: #38bdf8;">●</span>
                            <span style="font-weight: 600; color: #f1f5f9;">${name}</span>
                        </div>
                        <span style="font-size: 0.68rem; color: #64748b;">orbital</span>
                    </div>
                `;
            });
        } else if (!hasEsp) {
            listEl.innerHTML = `<div class="empty-state" style="padding: 1.5rem 0.5rem; font-size: 0.75rem;">No orbitals or surfaces found in bundle.</div>`;
            return;
        }

        listEl.innerHTML = html;

        listEl.querySelectorAll('.mol-orb-item').forEach(item => {
            item.addEventListener('click', () => {
                listEl.querySelectorAll('.mol-orb-item').forEach(i => i.classList.remove('active'));
                item.classList.add('active');
                if (item.id === 'btn-mol-esp') {
                    this.renderEsp();
                } else {
                    this.selectOrbital(item.dataset.name);
                }
            });
        });

        const firstOrb = listEl.querySelector('.mol-orb-item[data-name]');
        if (firstOrb) {
            firstOrb.classList.add('active');
            this.selectOrbital(firstOrb.dataset.name);
        } else if (hasEsp) {
            const espBtn = listEl.querySelector('#btn-mol-esp');
            if (espBtn) {
                espBtn.classList.add('active');
                this.renderEsp();
            }
        } else {
            this.renderMoleculeModel();
        }
    }

    renderMoleculeModel() {
        if (!this.viewer || !this.bundleData) return false;
        if (this.bundleData.xyz) {
            this.viewer.addModel(this.bundleData.xyz, 'xyz');
            this.viewer.setStyle({}, {
                stick: { radius: 0.12, colorscheme: 'Jmol' },
                sphere: { scale: 0.25, colorscheme: 'Jmol' }
            });
            return true;
        }
        return false;
    }

    async selectOrbital(orbName) {
        this.activeOrbital = orbName;
        this.renderActiveOrbital();
    }

    async renderActiveOrbital() {
        if (!this.viewer || !this.bundleData) return;
        this.viewer.removeAllModels();
        this.viewer.removeAllShapes();
        this.viewer.removeAllSurfaces();

        const hasAtoms = this.renderMoleculeModel();
        const orbName = this.activeOrbital;
        const files = this.bundleData.files || {};

        const posKey = Object.keys(files).find(k => k === `${orbName}_pos.json` || k.endsWith(`/${orbName}_pos.json`));
        const negKey = Object.keys(files).find(k => k === `${orbName}_neg.json` || k.endsWith(`/${orbName}_neg.json`));

        if (posKey && negKey) {
            try {
                const posMesh = typeof files[posKey] === 'string' ? JSON.parse(files[posKey]) : files[posKey];
                const negMesh = typeof files[negKey] === 'string' ? JSON.parse(files[negKey]) : files[negKey];
                if (posMesh.vertices && posMesh.faces) {
                    this.viewer.addCustom({
                        vertexArr: posMesh.vertices, faceArr: posMesh.faces,
                        normalArr: posMesh.normals || undefined, color: $3Dmol.CC.color('#38bdf8'), opacity: 0.80
                    });
                }
                if (negMesh.vertices && negMesh.faces) {
                    this.viewer.addCustom({
                        vertexArr: negMesh.vertices, faceArr: negMesh.faces,
                        normalArr: negMesh.normals || undefined, color: $3Dmol.CC.color('#f59e0b'), opacity: 0.80
                    });
                }
            } catch (err) {
                console.error('Error parsing orbital mesh:', err);
            }
        } else if (files[orbName] || files[`${orbName}.cube`]) {
            const cubeText = files[orbName] || files[`${orbName}.cube`];
            const vdo = new $3Dmol.VolumeData(cubeText, 'cube');
            if (!hasAtoms) {
                this.viewer.addModel(cubeText, 'cube');
                this.viewer.setStyle({}, { stick: { radius: 0.12, colorscheme: 'Jmol' }, sphere: { scale: 0.25, colorscheme: 'Jmol' } });
            }
            this.viewer.addIsosurface(vdo, { isoval: Math.abs(this.isovalue), color: '#38bdf8', opacity: 0.75 });
            this.viewer.addIsosurface(vdo, { isoval: -Math.abs(this.isovalue), color: '#f59e0b', opacity: 0.75 });
        }

        this.viewer.zoomTo();
        this.viewer.render();

        const statusBadge = document.querySelector('#modal-molecule-orbitals #mol-status-badge');
        if (statusBadge) statusBadge.textContent = `${this.activeOrbital} (isovalue: ±${this.isovalue.toFixed(3)})`;
    }

    async renderEsp() {
        if (!this.viewer || !this.bundleData) return;
        this.viewer.removeAllModels();
        this.viewer.removeAllShapes();
        this.viewer.removeAllSurfaces();

        this.renderMoleculeModel();

        const files = this.bundleData.files || {};
        const espKey = Object.keys(files).find(k => k.endsWith('_esp.json') || k === 'esp.json');
        if (espKey && files[espKey]) {
            try {
                const espData = typeof files[espKey] === 'string' ? JSON.parse(files[espKey]) : files[espKey];
                if (espData.vertices && espData.faces && espData.esp_values) {
                    const espMin = espData.esp_min || -30;
                    const espMax = espData.esp_max || 30;
                    const absMax = Math.max(Math.abs(espMin), Math.abs(espMax), 20);
                    const colors = espData.esp_values.map(v => {
                        const t = Math.max(-1, Math.min(1, v / absMax));
                        return t < 0 ? { r: 1.0, g: 1.0 + t, b: 1.0 + t } : { r: 1.0 - t, g: 1.0 - t, b: 1.0 };
                    });
                    this.viewer.addCustom({
                        vertexArr: espData.vertices, faceArr: espData.faces,
                        normalArr: espData.normals || undefined, color: colors, opacity: 0.85
                    });
                    const statusBadge = document.querySelector('#modal-molecule-orbitals #mol-status-badge');
                    if (statusBadge) statusBadge.textContent = `ESP Surface [${espMin.toFixed(1)} to +${espMax.toFixed(1)} kcal/mol]`;
                }
            } catch (err) {
                console.error('Error rendering ESP:', err);
            }
        }
        this.viewer.zoomTo();
        this.viewer.render();
    }
}
