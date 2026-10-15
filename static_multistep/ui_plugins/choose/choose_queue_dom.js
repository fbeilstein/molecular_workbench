/**
 * Queue and Catalog DOM Rendering for Step 1 Choose JSONs.
 */

export function renderSequenceQueue(queueListEl, badgeEl, state, callbacks) {
    if (badgeEl) badgeEl.innerText = `${state.sequence.length} Steps`;
    if (!queueListEl) return;

    if (state.sequence.length === 0) {
        queueListEl.innerHTML = `<div class="empty-state" style="padding: 1.5rem; font-size: 0.8rem; color: var(--text-muted);">No trajectories selected yet. Click "+ Add" on any trajectory from the catalog below.</div>`;
        return;
    }

    queueListEl.innerHTML = state.sequence.map((filename, idx) => {
        const summary = state.available.find(a => a.filename === filename);
        const rxnTitle = summary ? summary.reaction_smiles : filename;
        return `
            <div class="queue-item" data-filename="${filename}">
                <div style="display: flex; align-items: center; overflow: hidden; cursor: pointer;" class="queue-item-click">
                    <span class="queue-step-badge">Step ${idx + 1}</span>
                    <div class="queue-file-title" title="${rxnTitle}">${filename}</div>
                </div>
                <div class="queue-actions">
                    <button class="btn-icon" data-action="up" data-idx="${idx}" title="Move up" ${idx === 0 ? 'disabled' : ''}>▲</button>
                    <button class="btn-icon" data-action="down" data-idx="${idx}" title="Move down" ${idx === state.sequence.length - 1 ? 'disabled' : ''}>▼</button>
                    <button class="btn-icon danger" data-action="remove" data-idx="${idx}" title="Remove">✕</button>
                </div>
            </div>
        `;
    }).join('');

    queueListEl.querySelectorAll('.queue-item-click').forEach(el => {
        el.onclick = () => {
            const fn = el.closest('.queue-item').dataset.filename;
            callbacks.onSelect(fn);
        };
    });

    queueListEl.querySelectorAll('button[data-action]').forEach(btn => {
        btn.onclick = (e) => {
            e.stopPropagation();
            const action = btn.dataset.action;
            const idx = parseInt(btn.dataset.idx, 10);
            if (action === 'up' && idx > 0) {
                const temp = state.sequence[idx - 1];
                state.sequence[idx - 1] = state.sequence[idx];
                state.sequence[idx] = temp;
            } else if (action === 'down' && idx < state.sequence.length - 1) {
                const temp = state.sequence[idx + 1];
                state.sequence[idx + 1] = state.sequence[idx];
                state.sequence[idx] = temp;
            } else if (action === 'remove') {
                state.sequence.splice(idx, 1);
            }
            callbacks.onChanged();
        };
    });
}

export function renderCatalogList(catalogListEl, state, callbacks) {
    if (!catalogListEl) return;

    if (state.available.length === 0) {
        catalogListEl.innerHTML = `<div class="empty-state" style="padding: 1.5rem; font-size: 0.8rem;">No trajectories found. Upload one to get started.</div>`;
        return;
    }

    catalogListEl.innerHTML = state.available.map(item => {
        const isSelected = item.filename === state.activeFilename;
        const inQueue = state.sequence.includes(item.filename);
        const queueIdx = inQueue ? state.sequence.indexOf(item.filename) + 1 : null;

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

    catalogListEl.querySelectorAll('.catalog-item').forEach(el => {
        el.onclick = () => {
            const fn = el.dataset.filename;
            callbacks.onSelect(fn);
        };
    });

    catalogListEl.querySelectorAll('.btn-toggle-queue').forEach(btn => {
        btn.onclick = (e) => {
            e.stopPropagation();
            const fn = btn.dataset.filename;
            const idx = state.sequence.indexOf(fn);
            if (idx >= 0) {
                state.sequence.splice(idx, 1);
            } else {
                state.sequence.push(fn);
            }
            callbacks.onChanged();
        };
    });
}
