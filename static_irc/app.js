import { KetcherUI } from './ui_plugins/ketcher_ui.js';
import { FieldsUI } from './ui_plugins/fields_ui.js';
import { MappingUI } from './ui_plugins/mapping_ui.js';
import { TSEditorUI } from './ui_plugins/ts_editor_ui.js';
import { IRCViewerUI } from './ui_plugins/irc_viewer_ui.js';
import { CascadeUI } from './ui_plugins/cascade_ui.js';
import { OrbitalsUI } from './ui_plugins/orbitals_ui.js';
import { BundleViewerUI } from './ui_plugins/bundle_viewer_ui.js';

const UIRegistry = {
    'ketcher': KetcherUI,
    'mapping_ui': MappingUI,
    'fields': FieldsUI,
    'ts_editor_ui': TSEditorUI,
    'irc_viewer_ui': IRCViewerUI,
    'cascade_ui': CascadeUI,
    'orbitals_ui': OrbitalsUI,
    'bundle_viewer_ui': BundleViewerUI
};

async function executeStepRun(stepId, userInput) {
    const res = await fetch(`/api/step/${stepId}/run`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(userInput)
    });
    return await res.json();
}

document.addEventListener('DOMContentLoaded', () => {
    let plugins = [];
    let state = {};
    let project = {};
    let activeStepIndex = 0;

    const tabsContainer = document.getElementById('wizard-tabs');
    const schemaContainer = document.getElementById('schema-container');
    const errorContainer = document.getElementById('error-container');
    const runBtn = document.getElementById('btn-run');
    const resetBtn = document.getElementById('btn-reset');

    const reactionStepSelect = document.getElementById('reaction-step-select');
    const addStepBtn = document.getElementById('btn-add-step');
    const deleteStepBtn = document.getElementById('btn-delete-step');

    function updateProjectUI(proj) {
        if (!proj || !reactionStepSelect) return;
        project = proj;
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

    async function loadPipeline(preferredIndex = null) {
        const [pluginsRes, stateRes, projectRes] = await Promise.all([
            fetch('/api/plugins'),
            fetch('/api/state'),
            fetch('/api/project')
        ]);

        plugins = await pluginsRes.json();
        state = await stateRes.json();
        const projData = await projectRes.json();
        updateProjectUI(projData);

        const globalSolv = document.getElementById('global-solvent-select');
        if (globalSolv) {
            globalSolv.value = state.solvent || projData.solvent || 'water';
            globalSolv.onchange = async () => {
                const newSolv = globalSolv.value;
                state.solvent = newSolv;
                await fetch('/api/environment', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ solvent: newSolv })
                });
                renderActiveStep();
            };
        }

        determineActiveStep(preferredIndex);
        renderTabs();
        renderActiveStep();
    }

    function getMaxAvailableIndex() {
        let maxAvailable = 0;
        for (let i = 0; i < plugins.length; i++) {
            if (!state[plugins[i].step_id]) break;
            maxAvailable = (i + 1 < plugins.length) ? i + 1 : i;
        }
        const completedIRC = (project.steps || []).filter(s => s.has_irc).length;
        if (completedIRC >= 1 && plugins.length > 4) {
            maxAvailable = Math.max(maxAvailable, 4);
        }
        if (state['step_5_cascade'] && plugins.length > 5) {
            maxAvailable = Math.max(maxAvailable, 5);
        }
        if ((state['step_6_orbitals'] || state['step_5_cascade'] || plugins.find(p => p.step_id === 'step_7_bundle')?.ui_schema?.has_megajson) && plugins.length > 6) {
            maxAvailable = Math.max(maxAvailable, 6);
        }
        return maxAvailable;
    }

    function determineActiveStep(preferredIndex = null) {
        const maxAvailable = getMaxAvailableIndex();
        const hash = (window.location.hash || '').replace('#', '') || new URLSearchParams(window.location.search).get('step');
        if (hash) {
            const idx = plugins.findIndex((p, i) => p.step_id === hash || String(i + 1) === hash || p.step_name.toLowerCase().includes(hash.toLowerCase()));
            if (idx !== -1) {
                activeStepIndex = idx;
                return;
            }
        }
        if (preferredIndex !== null) {
            activeStepIndex = Math.max(0, Math.min(preferredIndex, maxAvailable));
        } else if (activeStepIndex !== null && activeStepIndex >= 0 && activeStepIndex <= maxAvailable) {
            // Keep current user tab selection without forceful auto-redirects
        } else {
            activeStepIndex = maxAvailable;
        }
    }

    function renderTabs() {
        tabsContainer.innerHTML = '';
        const maxAvailable = getMaxAvailableIndex();
        const completedIRC = (project.steps || []).filter(s => s.has_irc).length;
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
                activeStepIndex = index;
                try { window.location.hash = plugin.step_id; } catch (e) {}
                renderTabs();
                renderActiveStep();
            };

            tabsContainer.appendChild(btn);
        });
    }

    function renderActiveStep() {
        errorContainer.style.display = 'none';
        schemaContainer.innerHTML = '';

        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema;

        const envBar = document.querySelector('.env-bar');
        if (envBar) {
            const isCalcStep = plugin.step_id === 'step_3_ts' || plugin.step_id === 'step_4_irc';
            envBar.style.display = isCalcStep ? 'flex' : 'none';
        }

        const msg = document.createElement('p');
        msg.innerText = ui.message || '';
        schemaContainer.appendChild(msg);

        let handlerType = ui.type;
        if (ui.fields && handlerType !== 'ketcher' && handlerType !== 'mapping_ui') {
            handlerType = 'fields';
        }

        const handler = UIRegistry[handlerType];
        if (handler) {
            handler.render(schemaContainer, ui, state, plugin.step_id);
            if (activeStepIndex < plugins.length - 1) {
                runBtn.style.display = 'inline-block';
                runBtn.innerText = "Next";
            } else {
                runBtn.style.display = 'none';
            }
        } else {
            const infoMsg = ui.message || (handlerType ? `Unknown UI schema type: ${handlerType}` : 'Prerequisite step required. Please complete previous steps.');
            schemaContainer.innerHTML = `
                <div style="padding: 40px 24px; text-align: center; background: rgba(255,255,255,0.02); border: 1px solid rgba(255,255,255,0.08); border-radius: 12px; margin: 24px 0;">
                    <div style="font-size: 2.2rem; margin-bottom: 12px;">ℹ️</div>
                    <div style="font-size: 1.15rem; font-weight: 600; color: var(--text-main, #fff); margin-bottom: 8px;">Action Required</div>
                    <p style="color: var(--text-dim, #999); font-size: 0.95rem; max-width: 520px; margin: 0 auto; line-height: 1.5;">${infoMsg}</p>
                </div>
            `;
            runBtn.style.display = 'none';
        }
    }

    // Reaction Step Selector Event
    if (reactionStepSelect) {
        reactionStepSelect.onchange = async () => {
            const targetId = reactionStepSelect.value;
            const curTabIndex = activeStepIndex;
            try {
                const res = await fetch('/api/project/step/select', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ step_id: targetId })
                });
                const r = await res.json();
                if (r.success) {
                    await loadPipeline(curTabIndex);
                } else {
                    alert('Failed to switch step: ' + (r.error || 'Unknown error'));
                }
            } catch (err) {
                alert('Request failed: ' + err.message);
            }
        };
    }

    // Add Step Event (seeds product of current step)
    if (addStepBtn) {
        addStepBtn.onclick = async () => {
            addStepBtn.disabled = true;
            try {
                const res = await fetch('/api/project/step/add', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ seed: true })
                });
                const r = await res.json();
                if (r.success) {
                    await loadPipeline();
                    activeStepIndex = 0; // Go to Tab 1 (Reaction Input) of the new step
                    renderTabs();
                    renderActiveStep();
                } else {
                    alert('Failed to add step: ' + (r.error || 'Unknown error'));
                }
            } catch (err) {
                alert('Add step request failed: ' + err.message);
            } finally {
                addStepBtn.disabled = false;
            }
        };
    }

    // Delete Step Event
    if (deleteStepBtn) {
        deleteStepBtn.onclick = async () => {
            const curStep = reactionStepSelect ? reactionStepSelect.value : 'step_1';
            if (!confirm(`Are you sure you want to delete ${curStep}?`)) return;
            deleteStepBtn.disabled = true;
            try {
                const res = await fetch('/api/project/step/delete', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ step_id: curStep })
                });
                const r = await res.json();
                if (r.success) {
                    await loadPipeline();
                } else {
                    alert('Failed to delete step: ' + (r.error || 'Unknown error'));
                }
            } catch (err) {
                alert('Delete step request failed: ' + err.message);
            } finally {
                deleteStepBtn.disabled = false;
            }
        };
    }


    runBtn.onclick = async () => {
        errorContainer.style.display = 'none';
        runBtn.disabled = true;

        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema;

        let handlerType = ui.type;
        if (ui.fields && handlerType !== 'ketcher' && handlerType !== 'mapping_ui') {
            handlerType = 'fields';
        }

        const handler = UIRegistry[handlerType];
        if (!handler) {
            runBtn.disabled = false;
            return;
        }

        let userInput = {};
        try {
            userInput = await handler.getData(ui);
        } catch(e) {
            errorContainer.innerText = e.message;
            errorContainer.style.display = 'block';
            runBtn.disabled = false;
            return;
        }

        try {
            const data = await executeStepRun(plugin.step_id, userInput);
            if (!data.success) {
                errorContainer.innerText = data.error || "Unknown error";
                errorContainer.style.display = 'block';
                return;
            }
            const nextIndex = Math.min(activeStepIndex + 1, plugins.length - 1);
            await loadPipeline(nextIndex);
        } catch (err) {
            errorContainer.innerText = err.message;
            errorContainer.style.display = 'block';
        } finally {
            runBtn.disabled = false;
        }
    };

    resetBtn.onclick = async () => {
        const totalSteps = project.steps?.length || 1;
        let resetAll = false;
        if (totalSteps > 1) {
            resetAll = confirm("Click OK to reset the ENTIRE project (all reaction steps), or Cancel to reset ONLY the currently active step.");
            if (!resetAll) {
                if (!confirm(`Reset only the active step (${project.active_step_id})?`)) return;
            }
        } else {
            if (!confirm("Reset the current reaction step?")) return;
        }

        await fetch('/api/reset', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ all_steps: resetAll })
        });
        await loadPipeline();
    };

    const saveProjectBtn = document.getElementById('btn-save-project');
    const openProjectBtn = document.getElementById('btn-open-project');
    const fileOpenProjectInput = document.getElementById('file-open-project');

    if (saveProjectBtn) {
        saveProjectBtn.onclick = () => {
            window.location.href = '/api/project/archive/export';
        };
    }

    if (openProjectBtn && fileOpenProjectInput) {
        openProjectBtn.onclick = () => {
            fileOpenProjectInput.value = '';
            fileOpenProjectInput.click();
        };

        fileOpenProjectInput.onchange = async () => {
            const file = fileOpenProjectInput.files[0];
            if (!file) return;

            openProjectBtn.disabled = true;
            openProjectBtn.textContent = 'Opening...';
            try {
                const res = await fetch('/api/project/archive/import', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/octet-stream' },
                    body: file
                });
                const r = await res.json();
                if (r.success) {
                    alert(r.message || 'Project successfully loaded!');
                    await loadPipeline(0);
                } else {
                    alert('Failed to load project: ' + (r.error || 'Unknown error'));
                }
            } catch (err) {
                alert('Project import failed: ' + err.message);
            } finally {
                openProjectBtn.disabled = false;
                openProjectBtn.textContent = '📂 Open Project';
            }
        };
    }

    loadPipeline();
});
