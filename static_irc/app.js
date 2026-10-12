
import { KetcherUI } from './ui_plugins/ketcher_ui.js';
import { FieldsUI } from './ui_plugins/fields_ui.js';
import { MappingUI } from './ui_plugins/mapping_ui.js';
import { TSEditorUI } from './ui_plugins/ts_editor_ui.js';
import { IRCViewerUI } from './ui_plugins/irc_viewer_ui.js';

const UIRegistry = {
    'ketcher': KetcherUI,
    'mapping_ui': MappingUI,
    'fields': FieldsUI,
    'ts_editor_ui': TSEditorUI,
    'irc_viewer_ui': IRCViewerUI
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
        } else {
            schemaContainer.innerText = "Unknown UI schema type: " + handlerType;
        }
        
        if (activeStepIndex < plugins.length - 1) {
            runBtn.style.display = 'inline-block';
            runBtn.innerText = "Next";
        } else {
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
        
        let userInput = {};
        const handler = UIRegistry[handlerType];
        
        try {
            if (handler) {
                userInput = await handler.getData(ui);
            }
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
