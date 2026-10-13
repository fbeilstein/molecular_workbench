import { ChooseJSONsUI } from './ui_plugins/choose_jsons_ui.js';
import { StitchAlignUI } from './ui_plugins/stitch_align_ui.js';
import { EnergyStitchUI } from './ui_plugins/energy_stitch_ui.js';
import { MegaViewerUI } from './ui_plugins/mega_viewer_ui.js';

const UIRegistry = {
    'choose_jsons_ui': ChooseJSONsUI,
    'stitch_align_ui': StitchAlignUI,
    'energy_stitch_ui': EnergyStitchUI,
    'mega_viewer_ui': MegaViewerUI
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
    const statusText = document.getElementById('footer-status-text');
    const runBtn = document.getElementById('btn-run');
    const resetBtn = document.getElementById('btn-reset');

    async function loadPipeline() {
        try {
            const [pluginsRes, stateRes] = await Promise.all([
                fetch('/api/plugins'),
                fetch('/api/state')
            ]);

            plugins = await pluginsRes.json();
            state = await stateRes.json();

            determineActiveStep();
            renderTabs();
            renderActiveStep();
        } catch (err) {
            showError("Failed to connect to backend server: " + err.message);
        }
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
            btn.className = `wizard-tab ${index === activeStepIndex ? 'active' : ''} ${isCompleted ? 'completed' : ''} ${isLocked ? 'locked' : ''}`;
            btn.disabled = isLocked;

            btn.innerHTML = `
                <span class="tab-num">${index + 1}</span>
                <span>${plugin.step_name.replace(/^\d+\.\s*/, '')}</span>
            `;

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
        hideError();
        schemaContainer.innerHTML = '';

        if (!plugins || plugins.length === 0) {
            schemaContainer.innerHTML = `<div class="empty-state">No pipeline steps found.</div>`;
            return;
        }

        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema || {};

        if (statusText) {
            statusText.innerText = ui.message || '';
        }

        const handler = UIRegistry[ui.type];
        if (handler) {
            handler.render(schemaContainer, ui, state, plugin.step_id);
        } else {
            schemaContainer.innerHTML = `<div class="alert-banner alert-error">Unknown UI type: ${ui.type}</div>`;
        }

        // Advance button configuration
        if (activeStepIndex < plugins.length - 1) {
            runBtn.style.display = 'inline-flex';
            runBtn.innerText = "Next";
        } else {
            runBtn.style.display = 'none';
        }
    }

    function showError(msg) {
        if (errorContainer) {
            errorContainer.innerText = msg;
            errorContainer.style.display = 'block';
        }
    }

    function hideError() {
        if (errorContainer) {
            errorContainer.style.display = 'none';
        }
    }

    runBtn.onclick = async () => {
        hideError();
        runBtn.disabled = true;

        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema || {};
        const handler = UIRegistry[ui.type];

        let userInput = {};
        try {
            if (handler && handler.getData) {
                userInput = await handler.getData(ui);
            }
        } catch (e) {
            showError(e.message);
            runBtn.disabled = false;
            return;
        }

        try {
            const data = await executeStepRun(plugin.step_id, userInput);
            if (!data.success) {
                showError(data.error || "Step execution failed");
                return;
            }
            await loadPipeline();
        } catch (err) {
            showError(err.message);
        } finally {
            runBtn.disabled = false;
        }
    };

    resetBtn.onclick = async () => {
        if (!confirm("Reset entire pathway pipeline state?")) return;
        try {
            await fetch('/api/reset', { method: 'POST' });
            await loadPipeline();
        } catch (err) {
            showError("Failed to reset: " + err.message);
        }
    };

    loadPipeline();
});
