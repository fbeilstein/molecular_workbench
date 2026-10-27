/**
 * Navigation, tab rendering, and step determination for the wizard.
 */

export function updateProjectUI(reactionStepSelect, deleteStepBtn, proj) {
    if (!proj || !reactionStepSelect) return;
    reactionStepSelect.innerHTML = '';

    const steps = proj.steps || [];
    steps.forEach(s => {
        const opt = document.createElement('option');
        opt.value = s.id;
        
        let statusTag = '[New]';
        if (s.has_irc) {
            const bFwd = s.barrier_fwd != null ? `${s.barrier_fwd.toFixed(1)}k` : '';
            statusTag = bFwd ? `[IRC ✓ ${bFwd}]` : '[IRC ✓]';
        } else if (s.has_ts) {
            statusTag = '[TS ✓]';
        } else if (s.has_mapping) {
            statusTag = '[Mapped ✓]';
        } else if (s.has_input) {
            statusTag = '[Input ✓]';
        }

        const cleanTitle = s.title || `Step ${s.index}: Reaction`;
        const shortTitle = cleanTitle.length > 34 ? cleanTitle.slice(0, 31) + '...' : cleanTitle;
        opt.textContent = `${shortTitle} ${statusTag}`;
        if (s.id === proj.active_step_id) {
            opt.selected = true;
        }
        reactionStepSelect.appendChild(opt);
    });

    if (deleteStepBtn) {
        deleteStepBtn.disabled = (steps.length <= 1 && !steps[0]?.has_input);
    }
}

export function getMaxAvailableIndex(plugins, state, project) {
    let maxAvailable = 0;
    for (let i = 0; i < plugins.length; i++) {
        if (!state[plugins[i].step_id]) break;
        maxAvailable = (i + 1 < plugins.length) ? i + 1 : i;
    }
    const completedIRC = (project?.steps || []).filter(s => s.has_irc).length;
    if (completedIRC >= 1 && plugins.length > 4) {
        maxAvailable = Math.max(maxAvailable, 4);
    }
    if (state['step_5_cascade'] && plugins.length > 5) {
        maxAvailable = Math.max(maxAvailable, 5);
    }
    const hasMegajson = plugins.find(p => p.step_id === 'step_7_bundle')?.ui_schema?.has_megajson;
    if ((state['step_6_orbitals'] || state['step_5_cascade'] || hasMegajson) && plugins.length > 6) {
        maxAvailable = Math.max(maxAvailable, 6);
    }
    return maxAvailable;
}

export function determineActiveStep(plugins, state, project, currentActive, preferredIndex = null) {
    const maxAvailable = getMaxAvailableIndex(plugins, state, project);
    const hash = (window.location.hash || '').replace('#', '') || new URLSearchParams(window.location.search).get('step');
    if (hash) {
        const idx = plugins.findIndex((p, i) => p.step_id === hash || String(i + 1) === hash || p.step_name.toLowerCase().includes(hash.toLowerCase()));
        if (idx !== -1) {
            try { localStorage.setItem('irc_active_step', plugins[idx].step_id); } catch (e) {}
            return idx;
        }
    }
    if (preferredIndex !== null) {
        const idx = Math.max(0, Math.min(preferredIndex, maxAvailable));
        try { localStorage.setItem('irc_active_step', plugins[idx].step_id); } catch (e) {}
        return idx;
    }
    if (currentActive !== null && currentActive >= 0 && currentActive <= maxAvailable) {
        return currentActive;
    }
    try {
        const savedStep = localStorage.getItem('irc_active_step');
        if (savedStep) {
            const sIdx = plugins.findIndex(p => p.step_id === savedStep);
            if (sIdx !== -1 && sIdx <= maxAvailable) {
                return sIdx;
            }
        }
    } catch (e) {}
    return maxAvailable;
}

export function renderTabs(tabsContainer, plugins, state, project, activeStepIndex, onSelectStep) {
    if (!tabsContainer) return;
    tabsContainer.innerHTML = '';
    const maxAvailable = getMaxAvailableIndex(plugins, state, project);
    const completedIRC = (project?.steps || []).filter(s => s.has_irc).length;
    const hasCascade = !!state['step_5_cascade'];

    plugins.forEach((plugin, index) => {
        const isCompleted = !!state[plugin.step_id] || (plugin.step_id === 'step_5_cascade' && hasCascade);
        let isUnlocked = index <= maxAvailable || isCompleted;
        if (plugin.step_id === 'step_5_cascade' && completedIRC >= 1) isUnlocked = true;
        if (plugin.step_id === 'step_6_orbitals' && hasCascade) isUnlocked = true;
        if (plugin.step_id === 'step_7_bundle' && (plugin.ui_schema?.has_megajson || hasCascade || !!state['step_6_orbitals'])) isUnlocked = true;
        const isLocked = !isUnlocked;

        const btn = document.createElement('button');
        btn.className = `tab ${index === activeStepIndex ? 'active' : ''} ${isCompleted ? 'completed' : ''} ${isLocked ? 'locked' : ''}`;
        btn.innerText = plugin.step_name;
        btn.disabled = isLocked;

        btn.onclick = () => {
            if (isLocked) return;
            try {
                window.location.hash = plugin.step_id;
                localStorage.setItem('irc_active_step', plugin.step_id);
            } catch (e) {}
            onSelectStep(index, plugin);
        };

        tabsContainer.appendChild(btn);
    });
}
