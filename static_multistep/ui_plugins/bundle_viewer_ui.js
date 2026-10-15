/**
 * Step 5: Master Multistep Reaction Studio & Viewer facade.
 * Coordinates:
 * - Bundle3DViewer (3Dmol molecular coordinates, high-res orbital volumetric isosurfaces, ESP surface)
 * - BundlePesChart (Interactive 2D Free Energy Profile)
 * - SvgMechanismManager (Interactive reaction mechanism SVG with nested group hierarchy wiring,
 *   direction-aware electron-pushing arrow dynamic highlighting, reaction triggers, and state jumps)
 * - BundleHighresManager (Wavefunction grid recomputation at 48³, 56³, 64³)
 * - BundleMoleculeModal (Right-click molecule orbital inspector & Molecular Workbench ZIP bundle unpacker)
 * - Maximal-Compression Repacker (Level 9 Deflate standalone bundle)
 */

import { Bundle3DViewer } from './bundle/bundle_3d.js';
import { BundlePesChart } from './bundle/bundle_pes.js';
import { SvgMechanismManager } from './bundle/bundle_svg.js';
import { BundleHighresManager } from './bundle/bundle_highres.js';
import { BundleMoleculeModal } from './bundle/bundle_molecule_modal.js';
import { exportMultiModelXyz, exportMaxCompBundleZip } from './bundle/bundle_export.js';

export const BundleViewerUI = {
    state: {
        stepId: 'step_5_bundle_viewer',
        megajson: null,
        currentFrame: 0,
        isPlaying: false,
        playDirection: 'forward', // 'forward' or 'reverse'
        playbackSpeed: 1.0,
        playTimer: null,
        viewer3D: null,
        pesChart: null,
        svgManager: null,
        highresManager: null,
        moleculeModal: null,
        activeOrbitals: new Set(),
        showCharges: false,
        hideAllOrbitals: false,
        isovalue: 0.03,
        orbitalResolution: 32,
        tracks: [],
        frameDataCache: {},
        activeTab: 'interactive' // 'interactive' or 'wiring'
    },

    async render(container, uiSchema, globalState, stepId) {
        this.state.stepId = stepId;
        const megajson = uiSchema.megajson;
        this.state.megajson = megajson;

        if (!megajson || !megajson.frames) {
            container.innerHTML = `<div class="empty-state" style="padding: 2rem; text-align: center;">No stitched reaction pathway found. Complete earlier steps first.</div>`;
            return;
        }

        const totalFrames = megajson.frames.length;
        const timeline = uiSchema.orbital_timeline || megajson.orbital_timeline || { tracks: [] };
        this.state.tracks = timeline.tracks || [];
        this.state.isovalue = timeline.default_isovalue || 0.03;
        this.state.showCharges = this.state.tracks.some(t => t.orbital_key === 'track_esp_charge');
        this.state.orbitalResolution = uiSchema.orbital_resolution || megajson.orbital_resolution || 32;

        const initialActive = new Set();
        this.state.tracks.forEach(trk => {
            if (trk.default_on && trk.orbital_key) initialActive.add(trk.orbital_key);
        });
        this.state.activeOrbitals = initialActive;

        const ep = (megajson.energy_profile && megajson.energy_profile.frame_energies_kcal) ||
                   megajson.frames.map(f => f.rel_energy_kcal || 0.0);

        container.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 0.65rem; color: var(--text-main);">
                <!-- Top Master Bar with Action Buttons -->
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.6rem; background: rgba(0,0,0,0.3); padding: 0.6rem 0.9rem; border-radius: 8px; border: 1px solid rgba(255,255,255,0.08);">
                    <div>
                        <div style="display: flex; align-items: center; gap: 0.6rem;">
                            <h2 style="font-size: 1.1rem; font-weight: 700; color: #f8fafc; margin: 0;">Master Reaction Studio Viewer</h2>
                            <span id="badge-orbital-res" class="chip blue" style="font-weight: 600;">Grid: ${this.state.orbitalResolution}³</span>
                        </div>
                        <div style="font-size: 0.76rem; color: #94a3b8; display: flex; gap: 0.4rem; margin-top: 0.25rem; flex-wrap: wrap;">
                            <span class="chip">${megajson.symbols.length} atoms</span>
                            <span class="chip">${totalFrames} frames</span>
                            <span class="chip amber">Barrier: ${(megajson.overall_barrier_kcal || 0).toFixed(1)} kcal/mol</span>
                            <span class="chip emerald">Net ΔG: ${(megajson.overall_reaction_energy_kcal || 0).toFixed(1)} kcal/mol</span>
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.5rem; flex-wrap: wrap; align-items: center;">
                        <button id="btn-recompute-highres" class="btn btn-secondary btn-sm" style="display: flex; align-items: center; gap: 0.35rem; border-color: rgba(56,189,248,0.4);">
                            <span>⚡ Recompute Orbitals (Hi-Res)</span>
                        </button>
                        <label class="btn btn-secondary btn-sm" style="cursor: pointer; display: flex; align-items: center; gap: 0.35rem;">
                            <span>📁 Upload Mechanism SVG</span>
                            <input type="file" id="input-upload-svg" accept=".svg" style="display: none;">
                        </label>
                        <button id="btn-export-bundle-zip" class="btn btn-primary btn-sm" style="background: linear-gradient(135deg, #10b981, #059669); border-color: #059669; font-weight: 700;">
                            Export Bundle (.zip)
                        </button>
                        <button id="btn-export-xyz" class="btn btn-secondary btn-sm">Multi-Model XYZ</button>
                    </div>
                </div>

                <!-- Dual Directional Playback Transport Bar -->
                <div style="display: flex; align-items: center; justify-content: space-between; flex-wrap: wrap; gap: 0.6rem; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 8px; padding: 0.45rem 0.8rem;">
                    <!-- Playback Controls -->
                    <div style="display: flex; align-items: center; gap: 0.4rem;">
                        <button id="btn-play-fwd" class="btn btn-xs btn-primary" title="Play Forward (0 → End)">▶ Play Fwd</button>
                        <button id="btn-play-rev" class="btn btn-xs btn-amber" title="Play Reverse (End → 0)">◀ Play Rev</button>
                        <button id="btn-playback-pause" class="btn btn-xs btn-secondary" title="Pause">⏸ Pause</button>
                        <span id="badge-playback-dir" class="dir-badge dir-fwd">Forward ⏩</span>
                    </div>

                    <!-- Scrubber & Readout -->
                    <div style="display: flex; align-items: center; gap: 0.6rem; flex: 1; min-width: 280px; max-width: 650px;">
                        <input type="range" id="bundle-frame-scrubber" min="0" max="${totalFrames - 1}" value="0" style="flex: 1; cursor: pointer; accent-color: #38bdf8; height: 6px;">
                        <span id="bundle-frame-readout" style="font-size: 0.75rem; color: #94a3b8; font-family: monospace; min-width: 105px; text-align: right;">Frame 0 / ${totalFrames - 1}</span>
                    </div>

                    <!-- Speed Selector -->
                    <div style="display: flex; align-items: center; gap: 0.35rem; font-size: 0.75rem; color: #cbd5e1;">
                        <span>Speed:</span>
                        <select id="sel-playback-speed" class="input-select" style="font-size: 0.72rem; padding: 2px 5px;">
                            <option value="0.5">0.5x</option>
                            <option value="1.0" selected>1.0x</option>
                            <option value="1.5">1.5x</option>
                            <option value="2.0">2.0x</option>
                        </select>
                    </div>
                </div>

                <!-- Main Work Area: Left (Tracks + 3D Viewer), Right (Mechanism SVG & Group Wiring) -->
                <div style="display: flex; gap: 0.65rem; height: 500px; max-height: 500px;">
                    <!-- Left: Orbital Tracks Sidebar & 3D Viewer Container -->
                    <div style="flex: 1.1; display: flex; gap: 0.5rem; height: 100%;">
                        <!-- Tracks Sidebar -->
                        <div style="width: 220px; height: 100%; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 8px; padding: 0.5rem; display: flex; flex-direction: column; gap: 0.5rem; overflow: hidden;">
                            <span style="font-size: 0.78rem; font-weight: 600; color: #f1f5f9;">Configured Reaction Tracks</span>
                            <div id="bundle-tracks-list" style="flex: 1; min-height: 0; overflow-y: auto; display: flex; flex-direction: column; gap: 0.3rem; padding-right: 3px;"></div>

                            <div style="border-top: 1px solid rgba(255,255,255,0.08); padding-top: 0.4rem; display: flex; flex-direction: column; gap: 0.35rem;">
                                <label style="display: flex; align-items: center; gap: 0.4rem; font-size: 0.73rem; color: #cbd5e1; cursor: pointer;">
                                    <input type="checkbox" id="chk-bundle-esp" ${this.state.showCharges ? 'checked' : ''}>
                                    ⚡ ESP Surface & Charges
                                </label>
                                <label style="display: flex; align-items: center; gap: 0.4rem; font-size: 0.73rem; color: #cbd5e1; cursor: pointer;">
                                    <input type="checkbox" id="chk-bundle-mute" ${this.state.hideAllOrbitals ? 'checked' : ''}>
                                    🚫 Mute Orbitals (Lecture Mode)
                                </label>
                            </div>
                        </div>

                        <!-- 3D Canvas Box -->
                        <div style="flex: 1; height: 100%; position: relative; border-radius: 8px; overflow: hidden; border: 1px solid rgba(255,255,255,0.08); background: #090d16;">
                            <div id="bundle-3d-canvas" style="width: 100%; height: 100%;"></div>
                            <div style="position: absolute; top: 8px; right: 8px; display: flex; flex-direction: column; gap: 4px; z-index: 10;">
                                <button id="btn-bundle-zoom-in" class="btn btn-xs btn-secondary" title="Zoom In (+)">+</button>
                                <button id="btn-bundle-zoom-out" class="btn btn-xs btn-secondary" title="Zoom Out (-)">−</button>
                                <button id="btn-bundle-zoom-fit" class="btn btn-xs btn-secondary" title="Fit to View">⟲</button>
                            </div>
                            <div id="bundle-status-bar" style="position: absolute; bottom: 8px; left: 8px; font-size: 0.72rem; color: #94a3b8; background: rgba(15,23,42,0.85); padding: 2px 6px; border-radius: 4px; pointer-events: none;">
                                Ready
                            </div>
                        </div>
                    </div>

                    <!-- Right: Reaction Mechanism Studio & Group Wiring Panel -->
                    <div style="flex: 1.1; height: 100%; display: flex; flex-direction: column; background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 8px; overflow: hidden;">
                        <!-- Mode Tabs & Mechanism Toolbar -->
                        <div style="display: flex; justify-content: space-between; align-items: center; padding: 0.45rem 0.75rem; background: rgba(0,0,0,0.3); border-bottom: 1px solid rgba(255,255,255,0.08);">
                            <div style="display: flex; gap: 0.35rem;">
                                <button id="tab-mode-interactive" class="btn btn-xs btn-primary active">Interactive Mechanism</button>
                                <button id="tab-mode-wiring" class="btn btn-xs btn-secondary">Group Wiring Studio</button>
                            </div>
                            <div style="display: flex; gap: 0.4rem; align-items: center;">
                                <span id="label-wiring-hint" style="font-size: 0.7rem; color: #38bdf8;">Click arrows to play • Right-click molecules</span>
                            </div>
                        </div>

                        <!-- Mechanism SVG Canvas / Preview Area -->
                        <div id="svg-canvas-container" style="flex: 1; min-height: 0; position: relative; padding: 0.5rem; display: flex; align-items: center; justify-content: center; background: rgba(15,23,42,0.6); overflow: auto;">
                            <!-- Injected SVG -->
                        </div>

                        <!-- Bottom Wiring Inspector / Hierarchy Drawer (Visible in Wiring Mode) -->
                        <div id="wiring-inspector-drawer" style="display: none; height: 190px; border-top: 1px solid rgba(255,255,255,0.1); background: rgba(15,23,42,0.95); padding: 0.5rem; display: none; gap: 0.6rem;">
                            <!-- Groups Tree Sidebar -->
                            <div style="width: 220px; height: 100%; display: flex; flex-direction: column; background: rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.4rem; overflow: hidden;">
                                <span style="font-size: 0.74rem; font-weight: 600; color: #cbd5e1; margin-bottom: 0.3rem;">SVG Group Hierarchy</span>
                                <div id="svg-groups-tree" style="flex: 1; min-height: 0; overflow-y: auto; display: flex; flex-direction: column; gap: 0.25rem;"></div>
                            </div>

                            <!-- Selected Group Wiring Form -->
                            <div id="wiring-form-panel" style="flex: 1; height: 100%; background: rgba(0,0,0,0.3); border: 1px solid rgba(255,255,255,0.08); border-radius: 6px; padding: 0.5rem; display: flex; flex-direction: column; gap: 0.4rem; overflow-y: auto;">
                                <div style="font-size: 0.72rem; color: #94a3b8;">Select a group from the hierarchy tree or SVG canvas to configure actions.</div>
                            </div>
                        </div>
                    </div>
                </div>

                <!-- Interactive 2D Canvas PES Energy Profile Chart -->
                <div style="background: rgba(0,0,0,0.25); border: 1px solid rgba(255,255,255,0.08); border-radius: 8px; padding: 0.5rem;">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.25rem;">
                        <span style="font-size: 0.78rem; font-weight: 600; color: #e2e8f0;">Reaction Energy Profile & Free Energy Coordinates</span>
                        <span id="pes-legend-readout" style="font-size: 0.72rem; color: #94a3b8;">Click graph to jump to frame</span>
                    </div>
                    <canvas id="bundle-pes-canvas" style="width: 100%; height: 130px; display: block;"></canvas>
                </div>
            </div>
        `;

        // Initialize 3Dmol Viewer
        const canvas3D = container.querySelector('#bundle-3d-canvas');
        const statusBar = container.querySelector('#bundle-status-bar');
        this.state.viewer3D = new Bundle3DViewer(canvas3D, statusBar);
        this.state.viewer3D.init();

        // Initialize PES Chart
        const canvasPes = container.querySelector('#bundle-pes-canvas');
        this.state.pesChart = new BundlePesChart(canvasPes, (f) => this.seekFrame(f));
        this.state.pesChart.setData(ep, uiSchema.landmarks || []);

        // Initialize SVG Mechanism Manager
        const svgContainer = container.querySelector('#svg-canvas-container');
        this.state.svgManager = new SvgMechanismManager(svgContainer, {
            initialWiring: uiSchema.svg_wiring || megajson.svg_wiring || {},
            onPlaySegment: (sf, ef, dir) => this.playSegment(sf, ef, dir),
            onJumpFrame: (tf) => this.seekFrame(tf),
            onOpenMolecule: (molId, cfg) => this.openMoleculeOrbitals(molId, cfg),
            onSaveWiring: (w) => this.saveWiringToBackend(w),
            onGroupSelected: (grpId, cfg) => this.renderWiringForm(grpId, cfg)
        });

        const initialSvgText = uiSchema.mechanism_svg || megajson.mechanism_svg || '';
        this.state.svgManager.loadSvg(initialSvgText);

        // Initialize Highres & Molecule Managers
        this.state.highresManager = new BundleHighresManager({
            stepId: this.state.stepId,
            tracks: this.state.tracks,
            totalFrames: totalFrames,
            landmarks: uiSchema.landmarks || [],
            onComplete: (newRes) => this.onHighresCompleted(newRes)
        });

        this.state.moleculeModal = new BundleMoleculeModal({
            stepId: this.state.stepId,
            onBundleUploaded: (mId, data) => this.onMoleculeBundleUploaded(mId, data)
        });

        this.renderTracksList(container);
        this.renderGroupsTree();
        this.bindEvents(container);

        await this.seekFrame(0);
    },

    renderTracksList(container) {
        const listEl = container.querySelector('#bundle-tracks-list');
        if (!listEl) return;

        if (this.state.tracks.length === 0) {
            listEl.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.4rem;">No tracks configured. Return to Step 4 to audition tracks.</div>`;
            return;
        }

        let html = '';
        this.state.tracks.forEach(trk => {
            if (trk.orbital_key === 'track_esp_charge') {
                const isChecked = this.state.showCharges;
                html += `
                    <label style="display: flex; align-items: center; justify-content: space-between; padding: 0.35rem 0.5rem; background: rgba(56,189,248,0.06); border: 1px solid rgba(56,189,248,0.25); border-radius: 4px; font-size: 0.74rem; cursor: pointer;">
                        <div style="display: flex; align-items: center; gap: 0.4rem; overflow: hidden;">
                            <input type="checkbox" class="chk-track-esp" ${isChecked ? 'checked' : ''}>
                            <span style="color: #38bdf8; font-weight: 700;">⚡</span>
                            <span style="white-space: nowrap; overflow: hidden; text-overflow: ellipsis; color: #38bdf8; font-weight: 600;" title="${trk.title}">${trk.title}</span>
                        </div>
                    </label>
                `;
                return;
            }
            const isChecked = this.state.activeOrbitals.has(trk.orbital_key);
            html += `
                <label style="display: flex; align-items: center; justify-content: space-between; padding: 0.35rem 0.5rem; background: rgba(255,255,255,0.03); border: 1px solid rgba(255,255,255,0.06); border-radius: 4px; font-size: 0.74rem; cursor: pointer;">
                    <div style="display: flex; align-items: center; gap: 0.4rem; overflow: hidden;">
                        <input type="checkbox" class="chk-track-orb" data-key="${trk.orbital_key}" ${isChecked ? 'checked' : ''}>
                        <span style="color: ${trk.color || '#38bdf8'}; font-weight: 600;">●</span>
                        <span style="white-space: nowrap; overflow: hidden; text-overflow: ellipsis; color: #e2e8f0;" title="${trk.title}">${trk.title}</span>
                    </div>
                </label>
            `;
        });
        listEl.innerHTML = html;

        listEl.querySelectorAll('.chk-track-esp').forEach(chk => {
            chk.addEventListener('change', () => {
                this.state.showCharges = chk.checked;
                const bottomChk = container.querySelector('#chk-bundle-esp');
                if (bottomChk) bottomChk.checked = chk.checked;
                this.refreshViewer();
            });
        });

        listEl.querySelectorAll('.chk-track-orb').forEach(chk => {
            chk.addEventListener('change', () => {
                const k = chk.dataset.key;
                if (chk.checked) this.state.activeOrbitals.add(k);
                else this.state.activeOrbitals.delete(k);
                this.refreshViewer();
            });
        });
    },

    renderGroupsTree() {
        const treeEl = document.querySelector('#svg-groups-tree');
        if (!treeEl || !this.state.svgManager) return;

        const groups = this.state.svgManager.groupsTree || [];
        if (groups.length === 0) {
            treeEl.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.4rem;">No labeled groups found in SVG.</div>`;
            return;
        }

        let html = '';
        const renderNode = (grp) => {
            const cfg = this.state.svgManager.wiring[grp.id] || {};
            const role = cfg.role || 'none';
            let roleBadge = '';
            if (role === 'arrow') roleBadge = `<span class="chip-sm blue">➔ ${cfg.startFrame || 0}→${cfg.endFrame || 0}</span>`;
            else if (role === 'state') roleBadge = `<span class="chip-sm emerald">📍 F${cfg.targetFrame || 0}</span>`;
            else if (role === 'molecule') roleBadge = `<span class="chip-sm purple">⚛️ Mol</span>`;
            else if (role === 'push_arrow') roleBadge = `<span class="chip-sm amber">⚡ ${cfg.startFrame || 0}→${cfg.endFrame || 0} ${cfg.direction === 'reverse' ? '◀' : '▶'}</span>`;

            const indent = grp.depth * 14;
            html += `
                <div class="tree-group-node ${this.state.svgManager.selectedGroupId === grp.id ? 'active' : ''}" data-id="${grp.id}" style="padding-left: ${indent + 6}px; padding-top: 3px; padding-bottom: 3px; font-size: 0.73rem; cursor: pointer; display: flex; align-items: center; justify-content: space-between; border-radius: 4px; gap: 4px;">
                    <div style="display: flex; align-items: center; gap: 4px; overflow: hidden;">
                        <span style="color: #64748b;">${grp.hasChildren ? '▾' : '•'}</span>
                        <span style="font-weight: 600; color: #f1f5f9; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="${grp.label}">${grp.label}</span>
                    </div>
                    <div>${roleBadge}</div>
                </div>
            `;
            if (grp.children && grp.children.length > 0) {
                grp.children.forEach(c => renderNode(c));
            }
        };

        groups.forEach(g => renderNode(g));
        treeEl.innerHTML = html;

        treeEl.querySelectorAll('.tree-group-node').forEach(node => {
            node.addEventListener('click', () => {
                const grpId = node.dataset.id;
                this.state.svgManager.selectGroup(grpId);
            });
        });
    },

    renderWiringForm(grpId, cfg = {}) {
        const panel = document.querySelector('#wiring-form-panel');
        if (!panel) return;

        const role = cfg.role || 'none';
        const startFrame = cfg.startFrame !== undefined ? cfg.startFrame : 0;
        const endFrame = cfg.endFrame !== undefined ? cfg.endFrame : this.state.megajson.frames.length - 1;
        const targetFrame = cfg.targetFrame !== undefined ? cfg.targetFrame : this.state.currentFrame;
        const direction = cfg.direction || 'forward';
        const molName = cfg.moleculeName || grpId;

        panel.innerHTML = `
            <div style="display: flex; justify-content: space-between; align-items: center; border-bottom: 1px solid rgba(255,255,255,0.08); padding-bottom: 0.35rem;">
                <div style="display: flex; align-items: center; gap: 0.4rem;">
                    <span style="font-weight: 700; color: #f8fafc; font-size: 0.8rem;">Group: ${grpId}</span>
                </div>
                <button id="btn-save-group-wiring" class="btn btn-xs btn-primary">Save Wiring</button>
            </div>

            <!-- Role Selector -->
            <div style="display: flex; gap: 0.35rem; align-items: center; margin-top: 0.2rem;">
                <span style="font-size: 0.72rem; color: #94a3b8;">Role:</span>
                <button class="btn btn-xs btn-role ${role === 'arrow' ? 'active' : ''}" data-role="arrow">➔ Reaction Arrow</button>
                <button class="btn btn-xs btn-role ${role === 'state' ? 'active' : ''}" data-role="state">📍 State / Intermediate</button>
                <button class="btn btn-xs btn-role ${role === 'molecule' ? 'active' : ''}" data-role="molecule">⚛️ Molecule</button>
                <button class="btn btn-xs btn-role ${role === 'push_arrow' ? 'active' : ''}" data-role="push_arrow">⚡ Push-Arrow</button>
                <button class="btn btn-xs btn-role ${role === 'none' ? 'active' : ''}" data-role="none">⚪ None</button>
            </div>

            <!-- Role Details Editor -->
            <div id="role-fields-container" style="flex: 1; display: flex; flex-direction: column; gap: 0.35rem; padding-top: 0.2rem;">
                ${this.renderRoleFieldsHtml(role, { startFrame, endFrame, targetFrame, direction, molName })}
            </div>
        `;

        panel.querySelectorAll('.btn-role').forEach(btn => {
            btn.addEventListener('click', () => {
                panel.querySelectorAll('.btn-role').forEach(b => b.classList.remove('active'));
                btn.classList.add('active');
                const newRole = btn.dataset.role;
                const fieldsBox = panel.querySelector('#role-fields-container');
                fieldsBox.innerHTML = this.renderRoleFieldsHtml(newRole, { startFrame, endFrame, targetFrame, direction, molName });
                this.bindRoleFieldEvents(panel, grpId);
            });
        });

        this.bindRoleFieldEvents(panel, grpId);

        panel.querySelector('#btn-save-group-wiring').addEventListener('click', () => {
            const roleBtn = panel.querySelector('.btn-role.active');
            const activeRole = roleBtn ? roleBtn.dataset.role : 'none';
            const newConfig = this.gatherWiringForm(panel, activeRole);
            this.state.svgManager.setGroupWiring(grpId, newConfig);
            this.renderGroupsTree();
        });
    },

    renderRoleFieldsHtml(role, data) {
        if (role === 'arrow') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Start Frame:</label>
                        <input type="number" id="inp-arrow-start" value="${data.startFrame}" class="input-xs" style="width: 60px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-start">Cur</button>
                    </div>
                    <div>
                        <label>End Frame:</label>
                        <input type="number" id="inp-arrow-end" value="${data.endFrame}" class="input-xs" style="width: 60px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-end">Cur</button>
                    </div>
                    <div>
                        <label>Direction:</label>
                        <select id="sel-arrow-dir" class="input-xs">
                            <option value="forward" ${data.direction === 'forward' ? 'selected' : ''}>Forward ▶</option>
                            <option value="reverse" ${data.direction === 'reverse' ? 'selected' : ''}>Reverse ◀</option>
                        </select>
                    </div>
                </div>
                <div style="font-size: 0.7rem; color: #94a3b8;">Clicking this arrow in interactive mode will play the reaction between these two frames.</div>
            `;
        }
        if (role === 'state') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Target Frame:</label>
                        <input type="number" id="inp-state-target" value="${data.targetFrame}" class="input-xs" style="width: 60px;">
                        <button class="btn btn-xxs btn-secondary" id="btn-use-cur-target">Cur</button>
                    </div>
                </div>
                <div style="font-size: 0.7rem; color: #94a3b8;">Clicking this state will instantly jump the 3D viewer and scrubber to this frame.</div>
            `;
        }
        if (role === 'molecule') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Molecule Name:</label>
                        <input type="text" id="inp-mol-name" value="${data.molName}" class="input-xs" style="width: 140px;">
                    </div>
                    <button class="btn btn-xs btn-secondary" id="btn-inspect-mol">Open Orbitals</button>
                </div>
                <div style="font-size: 0.7rem; color: #94a3b8;">Right-clicking this molecule in interactive mode opens the full orbital studio window.</div>
            `;
        }
        if (role === 'push_arrow') {
            return `
                <div style="display: flex; gap: 0.8rem; align-items: center; font-size: 0.74rem;">
                    <div>
                        <label>Active Start:</label>
                        <input type="number" id="inp-push-start" value="${data.startFrame}" class="input-xs" style="width: 60px;">
                    </div>
                    <div>
                        <label>Active End:</label>
                        <input type="number" id="inp-push-end" value="${data.endFrame}" class="input-xs" style="width: 60px;">
                    </div>
                    <div>
                        <label>Direction:</label>
                        <select id="sel-push-dir" class="input-xs">
                            <option value="forward" ${data.direction === 'forward' ? 'selected' : ''}>Forward Only ▶</option>
                            <option value="reverse" ${data.direction === 'reverse' ? 'selected' : ''}>Reverse Only ◀</option>
                            <option value="both" ${data.direction === 'both' ? 'selected' : ''}>Both Directions ◀▶</option>
                        </select>
                    </div>
                </div>
                <div style="font-size: 0.7rem; color: #94a3b8;">This arrow will automatically glow/highlight when the reaction is in this frame range and direction.</div>
            `;
        }
        return `<div style="font-size: 0.72rem; color: #64748b;">No interactive action assigned.</div>`;
    },

    bindRoleFieldEvents(panel, grpId) {
        panel.querySelector('#btn-use-cur-start')?.addEventListener('click', () => {
            const inp = panel.querySelector('#inp-arrow-start');
            if (inp) inp.value = this.state.currentFrame;
        });
        panel.querySelector('#btn-use-cur-end')?.addEventListener('click', () => {
            const inp = panel.querySelector('#inp-arrow-end');
            if (inp) inp.value = this.state.currentFrame;
        });
        panel.querySelector('#btn-use-cur-target')?.addEventListener('click', () => {
            const inp = panel.querySelector('#inp-state-target');
            if (inp) inp.value = this.state.currentFrame;
        });
        panel.querySelector('#btn-inspect-mol')?.addEventListener('click', () => {
            const molName = panel.querySelector('#inp-mol-name')?.value || grpId;
            this.openMoleculeOrbitals(grpId, { moleculeName: molName });
        });
    },

    gatherWiringForm(panel, role) {
        if (role === 'none') return { role: 'none' };
        if (role === 'arrow') {
            return {
                role: 'arrow',
                startFrame: parseInt(panel.querySelector('#inp-arrow-start')?.value || 0, 10),
                endFrame: parseInt(panel.querySelector('#inp-arrow-end')?.value || 100, 10),
                direction: panel.querySelector('#sel-arrow-dir')?.value || 'forward'
            };
        }
        if (role === 'state') {
            return {
                role: 'state',
                targetFrame: parseInt(panel.querySelector('#inp-state-target')?.value || 0, 10)
            };
        }
        if (role === 'molecule') {
            return {
                role: 'molecule',
                moleculeName: panel.querySelector('#inp-mol-name')?.value || 'Molecule'
            };
        }
        if (role === 'push_arrow') {
            return {
                role: 'push_arrow',
                startFrame: parseInt(panel.querySelector('#inp-push-start')?.value || 0, 10),
                endFrame: parseInt(panel.querySelector('#inp-push-end')?.value || 100, 10),
                direction: panel.querySelector('#sel-push-dir')?.value || 'both'
            };
        }
        return { role: 'none' };
    },

    bindEvents(container) {
        // SVG Upload
        const svgInput = container.querySelector('#input-upload-svg');
        svgInput.addEventListener('change', (e) => {
            const file = e.target.files[0];
            if (!file) return;
            const reader = new FileReader();
            reader.onload = async () => {
                const svgText = reader.result;
                this.state.svgManager.loadSvg(svgText);
                this.renderGroupsTree();
                await this.saveSvgToBackend(svgText);
            };
            reader.readAsText(file);
        });

        // Mode Toggles (Interactive vs Wiring)
        const tabInteractive = container.querySelector('#tab-mode-interactive');
        const tabWiring = container.querySelector('#tab-mode-wiring');
        const drawer = container.querySelector('#wiring-inspector-drawer');
        const hint = container.querySelector('#label-wiring-hint');

        tabInteractive.addEventListener('click', () => {
            tabInteractive.classList.add('active', 'btn-primary');
            tabInteractive.classList.remove('btn-secondary');
            tabWiring.classList.remove('active', 'btn-primary');
            tabWiring.classList.add('btn-secondary');
            drawer.style.display = 'none';
            hint.textContent = 'Click arrows to play • Right-click molecules';
            this.state.svgManager.setMode('interactive');
        });

        tabWiring.addEventListener('click', () => {
            tabWiring.classList.add('active', 'btn-primary');
            tabWiring.classList.remove('btn-secondary');
            tabInteractive.classList.remove('active', 'btn-primary');
            tabInteractive.classList.add('btn-secondary');
            drawer.style.display = 'flex';
            hint.textContent = 'Click group in tree/canvas to wire actions';
            this.state.svgManager.setMode('edit');
            this.renderGroupsTree();
        });

        // Dual Direction Playback Buttons
        container.querySelector('#btn-play-fwd').addEventListener('click', () => this.playForward());
        container.querySelector('#btn-play-rev').addEventListener('click', () => this.playReverse());
        container.querySelector('#btn-playback-pause').addEventListener('click', () => this.pause());

        // Playback Speed
        container.querySelector('#sel-playback-speed').addEventListener('change', (e) => {
            this.state.playbackSpeed = parseFloat(e.target.value);
            if (this.state.isPlaying) {
                this.pause();
                if (this.state.playDirection === 'forward') this.playForward();
                else this.playReverse();
            }
        });

        // Scrubber
        const scrubber = container.querySelector('#bundle-frame-scrubber');
        if (scrubber) {
            scrubber.addEventListener('input', (e) => {
                this.pause();
                const newF = parseInt(e.target.value, 10);
                // Infer direction from scrub vector
                if (newF > this.state.currentFrame) this.setDirection('forward');
                else if (newF < this.state.currentFrame) this.setDirection('reverse');
                this.seekFrame(newF);
            });
        }

        // Recompute High-Res Orbitals
        container.querySelector('#btn-recompute-highres').addEventListener('click', () => {
            this.state.highresManager.renderModal(container);
        });

        // Export Actions
        container.querySelector('#btn-export-bundle-zip').addEventListener('click', () => exportMaxCompBundleZip(this.state.stepId));
        container.querySelector('#btn-export-xyz').addEventListener('click', () => exportMultiModelXyz(this.state.megajson));

        // 3D Viewer Zoom Buttons
        container.querySelector('#btn-bundle-zoom-in')?.addEventListener('click', () => this.state.viewer3D?.zoomIn());
        container.querySelector('#btn-bundle-zoom-out')?.addEventListener('click', () => this.state.viewer3D?.zoomOut());
        container.querySelector('#btn-bundle-zoom-fit')?.addEventListener('click', () => this.state.viewer3D?.resetCamera());

        // Mute & Charges
        container.querySelector('#chk-bundle-esp').addEventListener('change', (e) => {
            this.state.showCharges = e.target.checked;
            const topChk = container.querySelector('.chk-track-esp');
            if (topChk) topChk.checked = e.target.checked;
            this.refreshViewer();
        });
        container.querySelector('#chk-bundle-mute').addEventListener('change', (e) => {
            this.state.hideAllOrbitals = e.target.checked;
            this.refreshViewer();
        });
    },

    setDirection(dir) {
        this.state.playDirection = dir;
        const badge = document.querySelector('#badge-playback-dir');
        if (badge) {
            badge.className = `dir-badge ${dir === 'forward' ? 'dir-fwd' : 'dir-rev'}`;
            badge.textContent = dir === 'forward' ? 'Forward ⏩' : 'Reverse ⏪';
        }
    },

    playForward() {
        this.pause();
        this.setDirection('forward');
        const total = this.state.megajson.frames.length;
        if (this.state.currentFrame >= total - 1) this.seekFrame(0);

        this.state.isPlaying = true;
        const interval = Math.max(25, Math.round(110 / this.state.playbackSpeed));
        this.state.playTimer = setInterval(() => {
            if (this.state.currentFrame >= total - 1) {
                this.pause();
            } else {
                this.seekFrame(this.state.currentFrame + 1);
            }
        }, interval);
    },

    playReverse() {
        this.pause();
        this.setDirection('reverse');
        const total = this.state.megajson.frames.length;
        if (this.state.currentFrame <= 0) this.seekFrame(total - 1);

        this.state.isPlaying = true;
        const interval = Math.max(25, Math.round(110 / this.state.playbackSpeed));
        this.state.playTimer = setInterval(() => {
            if (this.state.currentFrame <= 0) {
                this.pause();
            } else {
                this.seekFrame(this.state.currentFrame - 1);
            }
        }, interval);
    },

    pause() {
        this.state.isPlaying = false;
        if (this.state.playTimer) {
            clearInterval(this.state.playTimer);
            this.state.playTimer = null;
        }
    },

    playSegment(startFrame, endFrame, direction = 'forward') {
        this.pause();
        this.setDirection(direction);

        const sf = Math.max(0, startFrame);
        const ef = Math.min(this.state.megajson.frames.length - 1, endFrame);

        this.seekFrame(direction === 'forward' ? sf : ef);
        this.state.isPlaying = true;

        const interval = Math.max(25, Math.round(95 / this.state.playbackSpeed));
        this.state.playTimer = setInterval(() => {
            if (direction === 'forward') {
                if (this.state.currentFrame >= ef) {
                    this.pause();
                } else {
                    this.seekFrame(this.state.currentFrame + 1);
                }
            } else {
                if (this.state.currentFrame <= sf) {
                    this.pause();
                } else {
                    this.seekFrame(this.state.currentFrame - 1);
                }
            }
        }, interval);
    },

    async seekFrame(f) {
        this.state.currentFrame = f;
        if (this.state.pesChart) this.state.pesChart.setCurrentFrame(f);

        const scrubber = document.querySelector('#bundle-frame-scrubber');
        if (scrubber) scrubber.value = f;

        const readout = document.querySelector('#bundle-frame-readout');
        if (readout && this.state.megajson) {
            readout.textContent = `Frame ${f} / ${this.state.megajson.frames.length - 1}`;
        }

        // Dynamically update electron-pushing arrows in mechanism SVG
        if (this.state.svgManager) {
            this.state.svgManager.updatePushArrows(f, this.state.playDirection);
        }

        await this.refreshViewer();
    },

    async fetchCube(frameIdx, orbKey, chemKey, isovalue) {
        const res = await fetch(`/api/step/${this.state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({
                action_type: 'get_cube',
                frame_index: frameIdx,
                orbital: orbKey,
                chem_key: chemKey,
                isovalue: isovalue
            })
        });
        return await res.json();
    },

    async refreshViewer() {
        const f = this.state.currentFrame;
        const frame = this.state.megajson.frames[f];
        if (!frame) return;

        if (!this.state.frameDataCache[f]) {
            try {
                const res = await fetch(`/api/step/step_4_orbitals/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_frame_orbitals', frame_index: f })
                });
                const d = await res.json();
                if (d.success) this.state.frameDataCache[f] = d;
            } catch (err) {}
        }

        const fd = this.state.frameDataCache[f];
        const activeTracks = this.state.tracks.filter(t => this.state.activeOrbitals.has(t.orbital_key));

        await this.state.viewer3D.renderScene({
            symbols: this.state.megajson.symbols,
            coords: frame.coords,
            charges: fd ? fd.charges : null,
            showCharges: this.state.showCharges,
            activeOrbitals: activeTracks,
            isovalue: this.state.isovalue,
            frameIdx: f,
            hideAllOrbitals: this.state.hideAllOrbitals,
            fetchCubeFn: (fi, ok, ck, iso) => this.fetchCube(fi, ok, ck, iso)
        });
    },

    openMoleculeOrbitals(molId, cfg) {
        const molName = (cfg && cfg.moleculeName) || molId;
        const existingBundle = (this.state.megajson.molecule_bundles && this.state.megajson.molecule_bundles[molId]) || null;
        this.state.moleculeModal.show(document.body, molId, molName, existingBundle);
    },

    onHighresCompleted(newResolution) {
        this.state.orbitalResolution = newResolution;
        const badge = document.querySelector('#badge-orbital-res');
        if (badge) badge.textContent = `Grid: ${newResolution}³ (High-Res)`;
        // Flush viewer cache to force reload of high-res cubes
        if (this.state.viewer3D) {
            this.state.viewer3D.cubeCache = {};
        }
        this.refreshViewer();
    },

    onMoleculeBundleUploaded(molId, bundleData) {
        if (!this.state.megajson.molecule_bundles) {
            this.state.megajson.molecule_bundles = {};
        }
        this.state.megajson.molecule_bundles[molId] = bundleData;
    },

    async saveSvgToBackend(svgText) {
        try {
            await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'save_mechanism_svg',
                    svg_text: svgText
                })
            });
        } catch (err) {}
    },

    async saveWiringToBackend(wiring) {
        try {
            await fetch(`/api/step/${this.state.stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({
                    action_type: 'save_svg_wiring',
                    wiring: wiring
                })
            });
        } catch (err) {}
    }
};
