/**
 * Junctions Queue DOM Builder & Metrics for Step 2 Stitch Align.
 */

export function renderJunctionsQueue(containerEl, badgeEl, state, callbacks) {
    const previews = state.junctionPreviews || [];
    if (badgeEl) badgeEl.innerText = `${previews.length} Junctions`;
    if (!containerEl) return;

    if (previews.length === 0) {
        containerEl.innerHTML = `<div class="empty-state" style="padding: 1.5rem; font-size: 0.8rem; color: var(--text-muted);">No multi-step junctions to align (single-step pathway).</div>`;
        return;
    }

    containerEl.innerHTML = previews.map((j, idx) => {
        const isSelected = idx === state.activeJunctionIdx;
        const currentChoice = state.junctionChoices[idx] || 'kabsch_bridge';

        return `
            <div class="queue-item ${isSelected ? 'selected-for-preview' : ''}" style="flex-direction: column; align-items: stretch; gap: 0.4rem; padding: 0.6rem 0.75rem;" data-j-idx="${idx}">
                <div style="display: flex; justify-content: space-between; align-items: center; cursor: pointer;" class="j-header-click">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <span class="queue-step-badge">Junction ${idx + 1}</span>
                        <span style="font-size: 0.775rem; font-weight: 600; color: #f8fafc;">
                            Step ${j.step_prev_index} ➔ Step ${j.step_next_index}
                        </span>
                    </div>
                    <span class="chip ${j.rmsd < 0.25 ? 'emerald' : 'amber'}" style="font-family: var(--font-mono); font-size: 0.7rem;">
                        RMSD: ${(j.rmsd || 0).toFixed(3)} Å
                    </span>
                </div>

                <div style="display: flex; justify-content: space-between; align-items: center; font-size: 0.725rem; color: var(--text-muted);">
                    <span>${j.hero_formula || 'Matched Core'}: ${j.n_matched_atoms} atoms</span>
                    <span style="color: #c084fc;">${(j.bridge_frames || []).length} Bridge Frames</span>
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
