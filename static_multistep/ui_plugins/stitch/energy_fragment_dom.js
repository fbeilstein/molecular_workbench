/**
 * Stoichiometry & Fragment Role DOM Builders for Step 3 Energy Stitch.
 */

export function renderFragmentInstrumentHTML(stoichiometry, fragmentRoles) {
    if (!stoichiometry || !stoichiometry.all_fragments || stoichiometry.all_fragments.length === 0) {
        return `<div class="empty-state" style="padding: 1rem; font-size: 0.75rem;">No independent fragments detected across the reaction steps.</div>`;
    }

    return stoichiometry.all_fragments.map(frag => {
        const fKey = frag.key;
        const currentRole = fragmentRoles[fKey] || frag.suggested_role || 'core';
        const roleLabel = {
            'core': '🛡️ Hero Core Molecule',
            'incoming': '📥 Incoming Reactant',
            'leaving': '📤 Excluded Leaving Group',
            'bystander': '👀 Bystander / Spectator'
        }[currentRole] || currentRole;

        return `
            <div class="queue-item" style="flex-direction: column; align-items: stretch; gap: 0.4rem; padding: 0.5rem 0.6rem;" data-frag-key="${fKey}">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div style="display: flex; align-items: center; gap: 0.4rem;">
                        <span class="queue-step-badge" style="background: rgba(255,255,255,0.06); color: #e2e8f0; font-family: var(--font-mono);">${frag.formula || fKey}</span>
                        <span style="font-size: 0.75rem; color: var(--text-muted); font-family: var(--font-mono);">${frag.smiles || ''}</span>
                    </div>
                    <span class="chip ${currentRole === 'core' ? 'blue' : (currentRole === 'incoming' ? 'emerald' : 'amber')}" style="font-size: 0.65rem;">
                        ${roleLabel}
                    </span>
                </div>
                <div style="display: flex; justify-content: space-between; align-items: center; margin-top: 0.2rem; font-size: 0.725rem;">
                    <span style="color: var(--text-muted);">Role Assignment:</span>
                    <select class="select-frag-role" data-frag-key="${fKey}" style="background: var(--bg-surface); color: var(--text-main); border: 1px solid var(--border-color); border-radius: 4px; padding: 0.2rem 0.4rem; font-size: 0.725rem;">
                        <option value="core" ${currentRole === 'core' ? 'selected' : ''}>Hero Core Molecule</option>
                        <option value="incoming" ${currentRole === 'incoming' ? 'selected' : ''}>Incoming Reactant</option>
                        <option value="leaving" ${currentRole === 'leaving' ? 'selected' : ''}>Excluded Leaving Group</option>
                        <option value="bystander" ${currentRole === 'bystander' ? 'selected' : ''}>Spectator</option>
                    </select>
                </div>
            </div>
        `;
    }).join('');
}
