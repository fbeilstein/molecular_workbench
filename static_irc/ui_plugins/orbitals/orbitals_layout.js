/**
 * DOM Layout and toolbar scaffolding for Step 6: Reaction Orbital Studio.
 */
export function buildOrbitalsLayout(container, maxF) {
    container.innerHTML = `
        <div style="display: flex; flex-direction: column; gap: 0.5rem; color: #1e293b;">
            <div style="display: flex; justify-content: space-between; align-items: center; padding: 0.4rem 0.75rem; background: #ffffff; border-radius: 8px; border: 1px solid #e2e8f0; gap: 0.5rem; box-shadow: 0 1px 2px rgba(0,0,0,0.04);">
                <div style="display: flex; align-items: center; gap: 0.5rem; flex-shrink: 0;">
                    <span style="font-weight: 700; font-size: 0.95rem; color: #0f172a;">Reaction Orbital Studio</span>
                    <span class="badge" style="font-size: 0.7rem; background: #e0f2fe; color: #0284c7;">Step 6</span>
                </div>
                <div style="display: flex; align-items: center; gap: 0.4rem; flex: 1; max-width: 520px; margin: 0 0.75rem;">
                    <button id="btn-play-pause" class="btn btn-sm btn-primary" style="padding: 0.2rem 0.6rem; flex-shrink: 0;">▶ Play</button>
                    <button id="btn-prev-frame" class="btn btn-sm btn-secondary" style="flex-shrink: 0;">⏮</button>
                    <input type="range" id="transport-frame-scrubber" min="0" max="${maxF}" value="0" style="flex: 1; cursor: pointer; accent-color: #0284c7; height: 6px;">
                    <button id="btn-next-frame" class="btn btn-sm btn-secondary" style="flex-shrink: 0;">⏭</button>
                    <span id="transport-frame-counter" style="font-size: 0.74rem; font-weight: 600; color: #334155; min-width: 105px; text-align: right; flex-shrink: 0;">Frame 0 / ${maxF}</span>
                </div>
                <div style="display: flex; align-items: center; gap: 0.4rem; flex-shrink: 0;">
                    <button id="btn-toggle-lecture" class="btn btn-sm btn-secondary" title="Mute/Unmute Orbitals (M)">👁 Toggle MOs</button>
                    <button id="btn-toggle-esp" class="btn btn-sm btn-secondary" title="Toggle ESP / Charge Spheres (E)">⚡ ESP</button>
                    <button id="btn-precompute-all" class="btn btn-sm btn-outline" style="color: #059669; border-color: #059669; background: #ecfdf5;">⚡ Precompute</button>
                </div>
            </div>
            <div style="display: flex; gap: 0.5rem; height: 480px; max-height: 480px;">
                <div style="width: 340px; height: 100%; display: flex; flex-direction: column; background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 0.5rem; overflow: hidden; box-shadow: 0 1px 2px rgba(0,0,0,0.04);">
                    <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 0.4rem;">
                        <span style="font-size: 0.8rem; font-weight: 700; color: #1e293b;">Orbital Media Pool</span>
                        <div style="display: flex; gap: 0.2rem;" id="pool-category-tabs">
                            <button class="btn btn-xs btn-primary tab-cat" data-cat="all">All</button>
                            <button class="btn btn-xs btn-secondary tab-cat" data-cat="frontier">Frontier</button>
                            <button class="btn btn-xs btn-secondary tab-cat" data-cat="bonding">Bonds</button>
                            <button class="btn btn-xs btn-secondary tab-cat" data-cat="lone_pairs">LP</button>
                            <button class="btn btn-xs btn-secondary tab-cat" data-cat="antibonding">π* / Antibonds</button>
                            <button class="btn btn-xs btn-secondary tab-cat" data-cat="esp">ESP</button>
                        </div>
                    </div>
                    <input type="text" id="pool-search-input" placeholder="Search atom, type, or label..." style="width: 100%; background: #f8fafc; border: 1px solid #cbd5e1; border-radius: 4px; padding: 0.3rem 0.5rem; font-size: 0.74rem; color: #0f172a; margin-bottom: 0.4rem;">
                    <div id="orbitals-pool-container" style="flex: 1; min-height: 0; overflow-y: auto; padding-right: 4px;"></div>
                </div>
                <div style="flex: 1; height: 100%; position: relative; border-radius: 8px; overflow: hidden; border: 1px solid #cbd5e1; background: #0f172a;">
                    <div id="orbitals-3d-canvas" style="width: 100%; height: 100%;"></div>
                    <div id="orbitals-status-label" style="position: absolute; top: 10px; left: 10px; font-size: 0.72rem; color: #f8fafc; background: rgba(15,23,42,0.85); padding: 3px 8px; border-radius: 4px; border: 1px solid rgba(255,255,255,0.1); pointer-events: none; z-index: 10; max-width: 70%; white-space: nowrap; overflow: hidden; text-overflow: ellipsis; backdrop-filter: blur(4px);">Ready</div>
                    <div style="position: absolute; top: 10px; right: 10px; display: flex; flex-direction: column; gap: 4px; z-index: 10;">
                        <button id="btn-zoom-in" class="btn btn-xs btn-secondary" title="Zoom In (+)" style="font-size: 0.85rem; padding: 2px 7px; background: rgba(15,23,42,0.85); color: #fff;">+</button>
                        <button id="btn-zoom-out" class="btn btn-xs btn-secondary" title="Zoom Out (-)" style="font-size: 0.85rem; padding: 2px 7px; background: rgba(15,23,42,0.85); color: #fff;">−</button>
                        <button id="btn-zoom-fit" class="btn btn-xs btn-secondary" title="Fit to View (0)" style="font-size: 0.72rem; padding: 2px 5px; background: rgba(15,23,42,0.85); color: #fff;">⟲</button>
                    </div>
                    <div id="orbitals-viewer-overlay" style="position: absolute; bottom: 8px; left: 8px; font-size: 0.72rem; color: #94a3b8; background: rgba(15,23,42,0.75); padding: 2px 6px; border-radius: 4px; pointer-events: none;">Scroll wheel to zoom | Drag to rotate | Right-drag to pan</div>
                </div>
            </div>
            <div id="orbitals-timeline-container" style="background: #ffffff; border: 1px solid #e2e8f0; border-radius: 8px; overflow: hidden; box-shadow: 0 1px 2px rgba(0,0,0,0.04);"></div>
        </div>
    `;
}
