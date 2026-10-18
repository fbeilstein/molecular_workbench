
import { KetcherUI } from './ui_plugins/ketcher_ui.js';
import { FieldsUI } from './ui_plugins/fields_ui.js';
import { MappingUI } from './ui_plugins/mapping_ui.js';
import { TSEditorUI } from './ui_plugins/ts_editor_ui.js';
import { IRCViewerUI } from './ui_plugins/irc_viewer_ui.js';
import { CascadeUI } from './ui_plugins/cascade_ui.js';

const UIRegistry = {
    'ketcher': KetcherUI,
    'mapping_ui': MappingUI,
    'fields': FieldsUI,
    'ts_editor_ui': TSEditorUI,
    'irc_viewer_ui': IRCViewerUI,
    'cascade_ui': CascadeUI
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
    let activeStepIndex = 0;

    const tabsContainer = document.getElementById('wizard-tabs');
    const schemaContainer = document.getElementById('schema-container');
    const errorContainer = document.getElementById('error-container');
    const runBtn = document.getElementById('btn-run');
    const resetBtn = document.getElementById('btn-reset');

    async function loadPipeline() {
        const [pluginsRes, stateRes] = await Promise.all([
            fetch('/api/plugins'),
            fetch('/api/state')
        ]);
        
        plugins = await pluginsRes.json();
        state = await stateRes.json();
        
        const globalSolv = document.getElementById('global-solvent-select');
        if (globalSolv) {
            globalSolv.value = state.solvent || 'water';
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

        determineActiveStep();
        renderTabs();
        renderActiveStep();
    }

    function determineActiveStep() {
        activeStepIndex = 0;
        for (let i = 0; i < plugins.length; i++) {
            if (!state[plugins[i].step_id]) break;
            activeStepIndex = (i + 1 < plugins.length) ? i + 1 : i;
        }
    }

    function renderTabs() {
        tabsContainer.innerHTML = '';
        plugins.forEach((plugin, index) => {
            const isCompleted = !!state[plugin.step_id];
            const isLocked = !isCompleted && index > activeStepIndex;
            
            const btn = document.createElement('button');
            btn.className = `tab ${index === activeStepIndex ? 'active' : ''} ${isCompleted ? 'completed' : ''} ${isLocked ? 'locked' : ''}`;
            btn.innerText = plugin.step_name;
            btn.disabled = isLocked;
            
            btn.onclick = () => {
                if (isLocked) return;
                activeStepIndex = index;
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
        
        // Use UI Plugin Registry
        let handlerType = ui.type;
        if (ui.fields && handlerType !== 'ketcher' && handlerType !== 'mapping_ui') {
            handlerType = 'fields'; // fallback
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
            await loadPipeline();
        } catch (err) {
            errorContainer.innerText = err.message;
            errorContainer.style.display = 'block';
        } finally {
            runBtn.disabled = false;
        }
    };
    
    resetBtn.onclick = async () => {
        if (!confirm("Reset entire pipeline?")) return;
        await fetch('/api/reset', { method: 'POST' });
        await loadPipeline();
    };

    loadPipeline();
});
