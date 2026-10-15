/**
 * HTML Template Builder for Step 5 Bundle Viewer.
 * SVG Reaction Mechanism is positioned on top spanning full horizontal width.
 * 3D Coordinates, Orbital Tracks, and SVG Hierarchy/Wiring are arranged side-by-side underneath.
 */

export function getBundleViewerTemplate(megajson, totalFrames, state) {
    return `
        <div style="display: flex; flex-direction: column; gap: 0.65rem; color: var(--text-main);">
            <!-- Top Master Bar with Action Buttons -->
            <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 0.6rem; background: rgba(0,0,0,0.3); padding: 0.6rem 0.9rem; border-radius: 8px; border: 1px solid rgba(255,255,255,0.08);">
                <div>
                    <div style="display: flex; align-items: center; gap: 0.6rem;">
                        <h2 style="font-size: 1.1rem; font-weight: 700; color: #f8fafc; margin: 0;">Master Reaction Studio Viewer</h2>
                        <span id="badge-orbital-res" class="chip blue" style="font-weight: 600;">Grid: ${state.orbitalResolution}³</span>
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
                    <button id="btn-purge-cubes" class="btn btn-secondary btn-sm" style="color: #f87171; border-color: rgba(248,113,113,0.4);" title="Delete temporary .cube files to free up disk space">
                        <span>🗑 Purge Cubes</span>
                    </button>
                </div>
            </div>

            <!-- Dual Direction Playback Transport Bar -->
            <div style="display: flex; align-items: center; gap: 0.6rem; background: rgba(17,24,39,0.9); padding: 0.5rem 0.8rem; border-radius: 8px; border: 1px solid rgba(255,255,255,0.08);">
                <div style="display: flex; gap: 0.3rem;">
                    <button id="btn-play-fwd" class="btn btn-primary btn-sm" title="Play Forward">▶ Play Fwd</button>
                    <button id="btn-play-rev" class="btn btn-amber btn-sm" title="Play Reverse">◀ Play Rev</button>
                    <button id="btn-playback-pause" class="btn btn-secondary btn-sm" title="Pause">⏸ Pause</button>
                </div>
                <span id="badge-playback-dir" class="dir-badge dir-fwd">Forward ⏩</span>
                <input type="range" id="bundle-frame-scrubber" min="0" max="${totalFrames - 1}" value="0" class="scrub-slider" style="flex: 1;">
                <div style="display: flex; align-items: center; gap: 0.6rem;">
                    <select id="sel-playback-speed" class="input-xs" style="background: rgba(0,0,0,0.4); color: #fff; border: 1px solid rgba(255,255,255,0.1); border-radius: 4px; padding: 2px 4px;">
                        <option value="0.5">0.5x</option>
                        <option value="1.0" selected>1.0x</option>
                        <option value="2.0">2.0x</option>
                    </select>
                    <span id="bundle-frame-readout" style="font-family: monospace; font-size: 0.8rem; color: #94a3b8; min-width: 105px; text-align: right;">
                        Frame 0 / ${totalFrames - 1}
                    </span>
                </div>
            </div>

            <!-- Top Card: Interactive SVG Reaction Mechanism (Full-Width Panorama) -->
            <div style="display: flex; flex-direction: column; background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: 8px; overflow: hidden;">
                <div style="padding: 0.4rem 0.75rem; border-bottom: 1px solid var(--border-color); display: flex; justify-content: space-between; align-items: center; background: rgba(0,0,0,0.25);">
                    <div style="display: flex; gap: 0.6rem; align-items: center;">
                        <span style="font-size: 0.82rem; font-weight: 700; color: #f8fafc;">Interactive Reaction Mechanism</span>
                        <span id="label-wiring-hint" style="font-size: 0.72rem; color: #94a3b8;">Click arrows to play • Right-click molecules</span>
                    </div>
                    <div style="display: flex; gap: 0.4rem; align-items: center;">
                        <button id="tab-mode-interactive" class="btn btn-xxs btn-primary active">Interactive</button>
                        <button id="tab-mode-wiring" class="btn btn-xxs btn-secondary">Element Wiring</button>
                    </div>
                </div>
                <div style="height: 270px; width: 100%; position: relative; overflow: hidden; background: #090d16;">
                    <div id="bundle-svg-view" style="width: 100%; height: 100%;"></div>
                </div>
            </div>

            <!-- Middle Row: 3D Viewport + Orbitals Sidebar + SVG Hierarchy & Wiring -->
            <div style="display: flex; gap: 0.75rem; min-height: 420px;">
                <!-- Left: 3D Molecular Coordinates & Volumetric Isosurfaces Viewport -->
                <div style="flex: 1.2; display: flex; flex-direction: column; background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: 8px; overflow: hidden; position: relative;">
                    <div style="padding: 0.4rem 0.75rem; border-bottom: 1px solid var(--border-color); display: flex; justify-content: space-between; align-items: center; font-size: 0.78rem;">
                        <span style="font-weight: 600;">Synchronized 3D Coordinate Space</span>
                        <div style="display: flex; gap: 0.4rem;">
                            <button id="btn-bundle-zoom-in" class="btn btn-xxs btn-secondary">+</button>
                            <button id="btn-bundle-zoom-out" class="btn btn-xxs btn-secondary">−</button>
                            <button id="btn-bundle-zoom-fit" class="btn btn-xxs btn-secondary">⟲</button>
                        </div>
                    </div>
                    <div id="bundle-viewport-3d" style="flex: 1; min-height: 380px; position: relative;"></div>
                </div>

                <!-- Middle: Tracks Sidebar -->
                <div style="width: 250px; display: flex; flex-direction: column; background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: 8px; padding: 0.5rem; gap: 0.4rem;">
                    <div style="font-size: 0.8rem; font-weight: 600; display: flex; justify-content: space-between; align-items: center;">
                        <span>Orbital Tracks</span>
                        <label style="font-size: 0.7rem; color: #94a3b8; display: flex; align-items: center; gap: 4px; cursor: pointer;">
                            <input type="checkbox" id="chk-bundle-mute" style="accent-color: #ef4444;"> Mute
                        </label>
                    </div>
                    <div style="display: flex; flex-direction: column; gap: 0.3rem;">
                        <div style="display: flex; justify-content: space-between; font-size: 0.7rem; color: #94a3b8;">
                            <span>Isovalue:</span>
                            <span id="label-bundle-isovalue">${state.isovalue}</span>
                        </div>
                        <input type="range" id="slider-bundle-isovalue" min="0.01" max="0.08" step="0.005" value="${state.isovalue}" class="scrub-slider">
                    </div>
                    <div style="display: flex; gap: 0.75rem; align-items: center; margin-top: 2px;">
                        <label style="font-size: 0.72rem; color: #cbd5e1; display: flex; align-items: center; gap: 5px; cursor: pointer;">
                            <input type="checkbox" id="chk-bundle-charges" ${state.showCharges ? 'checked' : ''} style="accent-color: #a855f7;"> 🏷 Atom Charges
                        </label>
                        <label style="font-size: 0.72rem; color: #cbd5e1; display: flex; align-items: center; gap: 5px; cursor: pointer;">
                            <input type="checkbox" id="chk-bundle-esp-surface" ${state.showEspSurface ? 'checked' : ''} style="accent-color: #38bdf8;"> ⚡ ESP Surface
                        </label>
                    </div>
                    <div id="bundle-tracks-list" style="flex: 1; overflow-y: auto; display: flex; flex-direction: column; gap: 0.25rem;"></div>
                </div>

                <!-- Right: SVG Hierarchy & Wiring Panel -->
                <div id="wiring-inspector-drawer" style="width: 280px; display: flex; flex-direction: column; background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: 8px; overflow: hidden;">
                    <div style="padding: 0.4rem 0.6rem; border-bottom: 1px solid var(--border-color); font-size: 0.75rem; font-weight: 600; display: flex; justify-content: space-between; align-items: center;">
                        <span>SVG Hierarchy & Wiring</span>
                        <span id="label-wiring-status" class="chip blue" style="font-size: 0.65rem;">Hierarchy</span>
                    </div>
                    <div id="group-tree-container" style="flex: 1; overflow-y: auto; padding: 0.35rem; min-height: 160px;"></div>
                    <div id="group-wiring-form" style="border-top: 1px solid var(--border-color); padding: 0.5rem; background: rgba(0,0,0,0.2);"></div>
                </div>
            </div>

            <!-- Bottom: 2D PES Chart -->
            <div style="background: var(--bg-surface); border: 1px solid var(--border-color); border-radius: 8px; padding: 0.6rem; height: 140px;">
                <canvas id="bundle-pes-canvas" style="width: 100%; height: 100%; display: block;"></canvas>
            </div>
        </div>
    `;
}
