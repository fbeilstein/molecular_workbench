/**
 * Junctions Queue DOM Builder & Metrics for Step 2 Stitch Align.
 */

export function buildStitchAlignHTML(state) {
    const nBridge = state.nBridgeFrames || 8;
    const viewMode = (state.viewer && state.viewer.viewMode) || 'transition';

    return `
        <!-- Top High-Visibility Action Bar -->
        <div class="stitch-header-banner" style="display: flex; justify-content: space-between; align-items: center; background: linear-gradient(135deg, rgba(30, 41, 59, 0.85), rgba(15, 23, 42, 0.95)); border: 1px solid rgba(59, 130, 246, 0.35); border-radius: var(--radius-md); padding: 0.85rem 1.25rem; margin-bottom: 1.25rem; box-shadow: 0 4px 16px rgba(0,0,0,0.35); flex-wrap: wrap; gap: 0.75rem;">
            <div style="display: flex; align-items: center; gap: 0.85rem;">
                <div style="width: 38px; height: 38px; border-radius: 8px; background: linear-gradient(135deg, #2563eb, #7c3aed); display: flex; align-items: center; justify-content: center; font-size: 1.25rem; box-shadow: 0 4px 12px rgba(37, 99, 235, 0.4);">
                    🔗
                </div>
                <div>
                    <div style="font-weight: 700; font-size: 0.95rem; color: #f8fafc;">Step 2: Coordinate Alignment & Bridge Stitching</div>
                    <div id="stitch-banner-subtitle" style="font-size: 0.775rem; color: #94a3b8;">
                        Verify intermediate Kabsch overlap and bridge continuity across reaction steps.
                    </div>
                </div>
            </div>
            <div style="display: flex; gap: 0.75rem; align-items: center;">
                <button id="btn-recalculate-align" class="btn btn-secondary" style="font-size: 0.825rem; padding: 0.5rem 1rem;">
                    ⚡ Re-Stitch Pathway
                </button>
            </div>
        </div>

        <div class="choose-container">
            <!-- Left Panel: Junctions, Bridge Interpolation & Choices -->
            <div class="left-panel">
                <div class="queue-card" style="flex: 1;">
                    <div class="section-title">
                        <span>Pathway Junctions</span>
                        <span class="chip blue" id="junctions-count-badge">Analyzing...</span>
                    </div>
                    <div id="stitch-status-banner" style="background: rgba(59, 130, 246, 0.1); border: 1px solid rgba(59, 130, 246, 0.25); border-radius: var(--radius-sm); padding: 0.65rem 0.85rem; font-size: 0.8rem; margin-bottom: 0.5rem;">
                        Detecting intermediate molecule overlaps and calculating Hungarian Kabsch alignment...
                    </div>
                    <div id="junctions-list" class="queue-list" style="max-height: 280px; overflow-y: auto;"></div>
                </div>

                <!-- Bridge Frame Interpolation Settings Card -->
                <div class="catalog-card" style="flex: 0 0 auto;">
                    <div class="section-title">
                        <span>Intermittent Bridge Frames</span>
                        <div style="display: flex; align-items: center; gap: 0.5rem;">
                            <span class="chip emerald" id="n-bridge-val">${nBridge} Frames</span>
                            <button id="btn-apply-bridge" class="btn btn-xs btn-primary" style="font-size: 0.725rem; padding: 0.2rem 0.55rem; font-weight: 600;">
                                ⚡ Re-Stitch
                            </button>
                        </div>
                    </div>
                    <p style="font-size: 0.75rem; color: var(--text-muted); margin-bottom: 0.6rem;">
                        Generates smooth cubic smoothstep transition frames connecting Step K product to Step K+1 reactant to stitch pathways without geometric discontinuity.
                    </p>
                    <div style="display: flex; align-items: center; gap: 0.75rem;">
                        <span style="font-size: 0.75rem; font-family: var(--font-mono); color: var(--text-muted);">0</span>
                        <input type="range" id="n-bridge-slider" class="scrub-slider" min="0" max="20" step="1" value="${nBridge}" style="flex: 1;">
                        <span style="font-size: 0.75rem; font-family: var(--font-mono); color: var(--text-muted);">20</span>
                    </div>
                </div>

                <!-- Boundary Scissors & Frame Trimmer Card -->
                <div class="catalog-card" style="flex: 0 0 auto;">
                    <div class="section-title">
                        <span>✂️ Boundary Scissors & Trimmer</span>
                        <div style="display: flex; align-items: center; gap: 0.4rem;">
                            <span class="chip purple" id="trim-active-badge" style="font-size: 0.7rem;">Junction 1</span>
                            <button id="btn-auto-trim" class="btn btn-xs btn-primary" title="Snap to minimal RMSD overlap" style="font-size: 0.7rem; padding: 0.15rem 0.5rem; font-weight: 600;">
                                ⚡ Best Overlap
                            </button>
                            <button id="btn-reset-trim" class="btn btn-xs btn-secondary" title="Reset cuts to full trajectory" style="font-size: 0.7rem; padding: 0.15rem 0.4rem;">
                                ↺ Reset
                            </button>
                        </div>
                    </div>
                    <p style="font-size: 0.75rem; color: var(--text-muted); margin-bottom: 0.5rem;">
                        Cut overshooting product tail or reactant head frames to stitch at the exact minimal-strain intermediate conformer.
                    </p>

                    <!-- Step K Product Tail Cut -->
                    <div style="margin-bottom: 0.5rem;">
                        <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.725rem; margin-bottom: 0.2rem;">
                            <span style="color: #93c5fd; font-weight: 600;">Step K Product Tail Cut:</span>
                            <span class="chip blue" id="trim-prev-val" style="font-size: 0.7rem; padding: 0.1rem 0.35rem;">Frame -</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 0.5rem;">
                            <span style="font-size: 0.7rem; font-family: var(--font-mono); color: var(--text-muted);" id="trim-prev-min">0</span>
                            <input type="range" id="trim-prev-slider" class="scrub-slider" min="0" max="10" value="10" style="flex: 1;">
                            <span style="font-size: 0.7rem; font-family: var(--font-mono); color: var(--text-muted);" id="trim-prev-max">10</span>
                        </div>
                    </div>

                    <!-- Step K+1 Reactant Head Cut -->
                    <div style="margin-bottom: 0.5rem;">
                        <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.725rem; margin-bottom: 0.2rem;">
                            <span style="color: #f43f5e; font-weight: 600;">Step K+1 Reactant Head Cut:</span>
                            <span class="chip rose" id="trim-next-val" style="font-size: 0.7rem; padding: 0.1rem 0.35rem;">Frame -</span>
                        </div>
                        <div style="display: flex; align-items: center; gap: 0.5rem;">
                            <span style="font-size: 0.7rem; font-family: var(--font-mono); color: var(--text-muted);" id="trim-next-min">0</span>
                            <input type="range" id="trim-next-slider" class="scrub-slider" min="0" max="10" value="0" style="flex: 1;">
                            <span style="font-size: 0.7rem; font-family: var(--font-mono); color: var(--text-muted);" id="trim-next-max">10</span>
                        </div>
                    </div>

                    <!-- Overlap Metrics Bar -->
                    <div id="trim-stats-box" style="background: rgba(15, 23, 42, 0.6); border: 1px solid rgba(255, 255, 255, 0.08); border-radius: 4px; padding: 0.35rem 0.5rem; font-size: 0.725rem; display: flex; justify-content: space-between; align-items: center;">
                        <span style="color: var(--text-muted);">Optimal Window Match:</span>
                        <span id="trim-recommended-label" style="font-weight: 600; color: #6ee7b7; font-family: var(--font-mono);">-</span>
                    </div>
                </div>

                <!-- Global Metrics Summary -->
                <div class="catalog-card" style="flex: 0 0 auto;">
                    <div class="section-title">Global Master Pathway Metrics</div>
                    <div style="display: grid; grid-template-columns: 1fr 1fr; gap: 0.5rem; font-size: 0.8rem; margin-bottom: 0.65rem;">
                        <div class="detail-box">
                            <div class="detail-box-title">Total Atoms</div>
                            <div id="metric-atoms" style="font-weight: 700; color: #93c5fd; font-family: var(--font-mono);">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Total Frames</div>
                            <div id="metric-frames" style="font-weight: 700; color: #93c5fd; font-family: var(--font-mono);">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Overall Barrier</div>
                            <div id="metric-barrier" style="font-weight: 700; color: #fcd34d; font-family: var(--font-mono);">-</div>
                        </div>
                        <div class="detail-box">
                            <div class="detail-box-title">Net Reaction Energy</div>
                            <div id="metric-rxn-energy" style="font-weight: 700; color: #6ee7b7; font-family: var(--font-mono);">-</div>
                        </div>
                    </div>
                </div>
            </div>

            <!-- Right Panel: Interactive 3D Boundary & Transition Viewer -->
            <div class="preview-panel">
                <div class="preview-header">
                    <div>
                        <div class="preview-rxn-title" id="junction-preview-title">Junction Alignment & Transition</div>
                        <div class="preview-meta-chips" id="junction-preview-chips"></div>
                    </div>
                    <div style="display: flex; gap: 0.6rem; align-items: center; font-size: 0.75rem; font-family: var(--font-mono);">
                        <span style="display: flex; align-items: center; gap: 0.3rem;">
                            <span style="display: inline-block; width: 10px; height: 10px; background: #38bdf8; border-radius: 50%;"></span>
                            Step K Product
                        </span>
                        <span style="display: flex; align-items: center; gap: 0.3rem;">
                            <span style="display: inline-block; width: 10px; height: 10px; background: #a855f7; border-radius: 50%;"></span>
                            Bridge
                        </span>
                        <span style="display: flex; align-items: center; gap: 0.3rem;">
                            <span style="display: inline-block; width: 10px; height: 10px; background: #f43f5e; border-radius: 50%;"></span>
                            Step K+1 Reactant
                        </span>
                    </div>
                </div>

                <!-- View Mode Selection Buttons -->
                <div style="display: flex; gap: 0.4rem; background: var(--bg-surface); padding: 0.4rem 0.6rem; border-radius: var(--radius-sm); border: 1px solid var(--border-color); margin-bottom: 0.6rem; flex-wrap: wrap;">
                    <button class="btn btn-xs btn-mode ${viewMode === 'transition' ? 'btn-primary' : 'btn-secondary'}" data-mode="transition">
                        🎬 Smoothstep Transition
                    </button>
                    <button class="btn btn-xs btn-mode ${viewMode === 'overlay' ? 'btn-primary' : 'btn-secondary'}" data-mode="overlay">
                        🔍 Superposition Overlay
                    </button>
                    <button class="btn btn-xs btn-mode ${viewMode === 'blink' ? 'btn-primary' : 'btn-secondary'}" data-mode="blink">
                        💡 Alternating Blink
                    </button>
                    <button class="btn btn-xs btn-mode ${viewMode === 'prev_only' ? 'btn-primary' : 'btn-secondary'}" data-mode="prev_only">
                        ⬅ Step K Final Only
                    </button>
                    <button class="btn btn-xs btn-mode ${viewMode === 'next_only' ? 'btn-primary' : 'btn-secondary'}" data-mode="next_only">
                        ➡ Step K+1 Start Only
                    </button>
                    <button class="btn btn-xs btn-mode ${viewMode === 'unaligned' ? 'btn-primary' : 'btn-secondary'}" data-mode="unaligned">
                        ⚠️ Raw Unaligned
                    </button>
                    <button class="btn btn-xs btn-mode ${viewMode === 'full_trajectory' ? 'btn-primary' : 'btn-secondary'}" data-mode="full_trajectory">
                        🎞 Entire Pathway
                    </button>
                </div>

                <!-- 3Dmol Boundary Viewport -->
                <div class="viewer-wrapper">
                    <div id="junction-viewport-3d" class="viewport-3d"></div>
                    <div id="junction-status-badge" class="viewer-overlay-badge">Computing Boundary Overlap...</div>
                </div>

                <!-- Transition Playback Controller -->
                <div class="playback-bar" id="transition-controls-bar" style="display: flex;">
                    <button id="btn-trans-play" class="btn btn-sm btn-secondary" title="Play / Pause">
                        <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                            <polygon points="5 3 19 12 5 21 5 3"/>
                        </svg>
                    </button>
                    <button id="btn-trans-prev" class="btn btn-sm btn-secondary" title="Previous Frame">⏮</button>
                    <input type="range" id="trans-scrubber" class="scrub-slider" min="0" max="10" value="0">
                    <button id="btn-trans-next" class="btn btn-sm btn-secondary" title="Next Frame">⏭</button>
                    <div id="trans-counter" class="frame-counter">1 / 10</div>
                </div>

                <!-- Details Inspector Grid -->
                <div class="preview-details" style="display: grid; grid-template-columns: 1fr 1fr 1fr; margin-top: 0.6rem;">
                    <div class="detail-box">
                        <div class="detail-box-title">Step K Product</div>
                        <div id="j-prev-info" style="font-family: var(--font-mono); font-size: 0.75rem; color: #93c5fd;">-</div>
                    </div>
                    <div class="detail-box">
                        <div class="detail-box-title">Aligned Boundary Fit</div>
                        <div id="alignment-rmsd-display" style="font-weight: 700; color: #6ee7b7; font-family: var(--font-mono);">-</div>
                    </div>
                    <div class="detail-box">
                        <div class="detail-box-title">Step K+1 Reactant</div>
                        <div id="j-next-info" style="font-family: var(--font-mono); font-size: 0.75rem; color: #f43f5e;">-</div>
                    </div>
                </div>
            </div>
        </div>
    `;
}

export function renderJunctionsQueue(containerEl, badgeEl, state, callbacks) {
    const previews = state.junctionPreviews || [];
    if (badgeEl) badgeEl.innerText = `${previews.length} Junctions`;
    if (!containerEl) return;

    if (previews.length === 0) {
        containerEl.innerHTML = `
            <div class="empty-state" style="padding: 1.25rem 0.75rem; font-size: 0.8rem; color: #94a3b8; background: rgba(59, 130, 246, 0.08); border: 1px dashed rgba(59, 130, 246, 0.25); border-radius: var(--radius-sm); text-align: center;">
                ✓ Continuous Pathway Loaded.<br>
                <span style="color: #6ee7b7; font-weight: 600;">No unaligned boundaries detected.</span>
            </div>
        `;
        return;
    }

    containerEl.innerHTML = previews.map((j, idx) => {
        const isSelected = idx === state.activeJunctionIdx;
        const currentChoice = state.junctionChoices[idx] || 'kabsch_bridge';
        const isTrimmed = j.is_trimmed;
        const pCut = j.prev_cut !== undefined ? j.prev_cut : (j.prev_total_frames ? j.prev_total_frames - 1 : '-');
        const nCut = j.next_cut !== undefined ? j.next_cut : 0;
        const hasRec = (!isTrimmed && j.rmsd_improvement > 0.02);

        return `
            <div class="queue-item ${isSelected ? 'selected-for-preview' : ''}" style="flex-direction: column; align-items: stretch; gap: 0.4rem; padding: 0.6rem 0.75rem;" data-j-idx="${idx}">
                <div style="display: flex; justify-content: space-between; align-items: center; cursor: pointer;" class="j-header-click">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <span class="queue-step-badge">Junction ${idx + 1}</span>
                        <span style="font-size: 0.775rem; font-weight: 600; color: #f8fafc;">
                            Step ${j.step_prev_index} ➔ Step ${j.step_next_index}
                        </span>
                    </div>
                    <div style="display: flex; align-items: center; gap: 0.35rem;">
                        ${isTrimmed ? '<span class="chip purple" style="font-family: var(--font-mono); font-size: 0.65rem; padding: 0.1rem 0.35rem;">✂️ Trimmed</span>' : ''}
                        ${hasRec ? `<span class="chip emerald" style="font-family: var(--font-mono); font-size: 0.65rem; padding: 0.1rem 0.35rem;" title="Optimal window overlap improves RMSD by ${j.rmsd_improvement.toFixed(3)} Å">⚡ +${j.rmsd_improvement.toFixed(2)} Å</span>` : ''}
                        <span class="chip ${j.rmsd < 0.25 ? 'emerald' : 'amber'}" style="font-family: var(--font-mono); font-size: 0.7rem;">
                            RMSD: ${(j.rmsd || 0).toFixed(3)} Å
                        </span>
                    </div>
                </div>

                <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.725rem; color: var(--text-muted);">
                    <span>${j.hero_formula || 'Matched Core'}: ${j.n_matched_atoms} atoms</span>
                    <span style="color: #c084fc;">${(j.bridge_frames || []).length} Bridge Frames</span>
                </div>

                <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.7rem; color: #94a3b8;">
                    <span>Active Frames: [0..${pCut}] ➔ [${nCut}..${(j.next_total_frames || 1) - 1}]</span>
                    <span style="color: #64748b;">Tot: ${j.prev_total_frames || '-'} | ${j.next_total_frames || '-'}</span>
                </div>

                <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 0.2rem; font-size: 0.725rem;">
                    <span style="color: var(--text-muted);">Stitch Strategy:</span>
                    <select class="select-j-strategy" data-j-idx="${idx}" style="background: var(--bg-surface); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; padding: 0.2rem 0.4rem; font-size: 0.725rem;">
                        <option value="kabsch_bridge" ${currentChoice === 'kabsch_bridge' ? 'selected' : ''}>Kabsch Fit + Bridge Frames</option>
                        <option value="kabsch_direct" ${currentChoice === 'kabsch_direct' ? 'selected' : ''}>Direct Kabsch (No Bridge)</option>
                        <option value="as_is" ${currentChoice === 'as_is' ? 'selected' : ''}>Raw Unaligned Join</option>
                    </select>
                </div>
            </div>
        `;
    }).join('');

    containerEl.querySelectorAll('.j-header-click').forEach(el => {
        el.onclick = () => {
            const idx = parseInt(el.closest('.queue-item').dataset.jIdx, 10);
            callbacks.onSelect(idx);
        };
    });

    containerEl.querySelectorAll('.select-j-strategy').forEach(sel => {
        sel.onchange = (e) => {
            const idx = parseInt(sel.dataset.jIdx, 10);
            callbacks.onStrategyChange(idx, e.target.value);
        };
    });
}

export function updateGlobalMetrics(megajson) {
    if (!megajson) return;
    const atomsEl = document.getElementById('metric-atoms');
    const framesEl = document.getElementById('metric-frames');
    const barrierEl = document.getElementById('metric-barrier');
    const rxnEl = document.getElementById('metric-rxn-energy');

    if (atomsEl) atomsEl.innerText = `${(megajson.symbols || []).length} atoms`;
    if (framesEl) framesEl.innerText = `${(megajson.frames || []).length} frames`;

    const ep = megajson.energy_profile;
    if (ep) {
        if (barrierEl) barrierEl.innerText = `${(ep.overall_barrier_kcal || 0).toFixed(2)} kcal/mol`;
        if (rxnEl) rxnEl.innerText = `${(ep.overall_reaction_energy_kcal || 0).toFixed(2)} kcal/mol`;
    }
}

export function updateTrimmerControls(j, idx) {
    if (!j) return;
    const badgeTrim = document.getElementById('trim-active-badge');
    if (badgeTrim) badgeTrim.innerText = `Junction ${idx + 1}`;

    const pSlider = document.getElementById('trim-prev-slider');
    const pVal = document.getElementById('trim-prev-val');
    const pMin = document.getElementById('trim-prev-min');
    const pMax = document.getElementById('trim-prev-max');

    const nSlider = document.getElementById('trim-next-slider');
    const nVal = document.getElementById('trim-next-val');
    const nMin = document.getElementById('trim-next-min');
    const nMax = document.getElementById('trim-next-max');

    const recLabel = document.getElementById('trim-recommended-label');

    const pTot = j.prev_total_frames || 1;
    const nTot = j.next_total_frames || 1;
    const pCut = (j.prev_cut !== undefined) ? j.prev_cut : pTot - 1;
    const nCut = (j.next_cut !== undefined) ? j.next_cut : 0;

    if (pSlider) {
        pSlider.min = 0;
        pSlider.max = Math.max(0, pTot - 1);
        pSlider.value = pCut;
    }
    if (pMin) pMin.innerText = '0';
    if (pMax) pMax.innerText = String(Math.max(0, pTot - 1));
    if (pVal) pVal.innerText = `Frame ${pCut} / ${Math.max(0, pTot - 1)}`;

    if (nSlider) {
        nSlider.min = 0;
        nSlider.max = Math.max(0, nTot - 1);
        nSlider.value = nCut;
    }
    if (nMin) nMin.innerText = '0';
    if (nMax) nMax.innerText = String(Math.max(0, nTot - 1));
    if (nVal) nVal.innerText = `Frame ${nCut} / ${Math.max(0, nTot - 1)}`;

    if (recLabel) {
        if (j.recommended_prev_cut !== undefined && j.recommended_next_cut !== undefined) {
            const imp = (j.rmsd_improvement || 0).toFixed(3);
            recLabel.innerText = `Frame ${j.recommended_prev_cut} ➔ ${j.recommended_next_cut} (RMSD: ${(j.recommended_rmsd || 0).toFixed(3)} Å, Δ: -${imp} Å)`;
        } else {
            recLabel.innerText = 'Optimal Overlap Analysis Ready';
        }
    }
}

