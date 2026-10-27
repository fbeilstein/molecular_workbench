import { KetcherUI } from './ui_plugins/ketcher_ui.js';
import { FieldsUI } from './ui_plugins/fields_ui.js';
import { MappingUI } from './ui_plugins/mapping_ui.js';
import { TSEditorUI } from './ui_plugins/ts_editor_ui.js';
import { IRCViewerUI } from './ui_plugins/irc_viewer_ui.js';
import { CascadeUI } from './ui_plugins/cascade_ui.js';
import { OrbitalsUI } from './ui_plugins/orbitals_ui.js';
import { BundleViewerUI } from './ui_plugins/bundle_viewer_ui.js?v=2.1.4';
import { fetchPipelineData, executeStepRun, setGlobalSolvent } from './api_client.js';
import { updateProjectUI, determineActiveStep, renderTabs } from './wizard_nav.js';
import { setupProjectControls } from './project_controls.js';

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

function boot() {
    let plugins = [];
    let state = {};
    let project = {};
    let activeStepIndex = null;

    const tabsContainer = document.getElementById('wizard-tabs');
    const schemaContainer = document.getElementById('schema-container');
    const errorContainer = document.getElementById('error-container');
    const runBtn = document.getElementById('btn-run');
    const resetBtn = document.getElementById('btn-reset');
    const reactionStepSelect = document.getElementById('reaction-step-select');
    const addStepBtn = document.getElementById('btn-add-step');
    const deleteStepBtn = document.getElementById('btn-delete-step');
    const saveProjectBtn = document.getElementById('btn-save-project');
    const openProjectBtn = document.getElementById('btn-open-project');
    const fileOpenProjectInput = document.getElementById('file-open-project');

    async function loadPipeline(preferredIndex = null) {
        const data = await fetchPipelineData();
        plugins = data.plugins;
        state = data.state;
        project = data.project;

        updateProjectUI(reactionStepSelect, deleteStepBtn, project);
        setupSolvent();
        activeStepIndex = determineActiveStep(plugins, state, project, activeStepIndex, preferredIndex);
        refreshTabs();
        renderActiveStep();
    }

    function refreshTabs() {
        renderTabs(tabsContainer, plugins, state, project, activeStepIndex, (idx) => {
            activeStepIndex = idx;
            refreshTabs();
            renderActiveStep();
        });
    }

    function setupSolvent() {
        const globalSolv = document.getElementById('global-solvent-select');
        if (globalSolv) {
            globalSolv.value = state.solvent || project.solvent || 'water';
            globalSolv.onchange = async () => {
                state.solvent = globalSolv.value;
                await setGlobalSolvent(globalSolv.value);
                renderActiveStep();
            };
        }
    }

    function renderActiveStep() {
        errorContainer.style.display = 'none';
        schemaContainer.innerHTML = '';
        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema;

        const envBar = document.querySelector('.env-bar');
        if (envBar) {
            const isCalc = plugin.step_id === 'step_3_ts' || plugin.step_id === 'step_4_irc';
            envBar.style.display = isCalc ? 'flex' : 'none';
        }

        const msg = document.createElement('p');
        msg.innerText = ui.message || '';
        schemaContainer.appendChild(msg);

        let handlerType = ui.type;
        if (ui.fields && handlerType !== 'ketcher' && handlerType !== 'mapping_ui') {
            handlerType = 'fields';
        }

        const applyFeedback = document.getElementById('apply-feedback');
        if (applyFeedback) applyFeedback.style.display = 'none';

        const handler = UIRegistry[handlerType];
        if (handler) {
            handler.render(schemaContainer, ui, state, plugin.step_id);
            runBtn.style.display = 'inline-flex';
            runBtn.innerHTML = '💾 Apply';
            runBtn.title = `Apply and save ${plugin.step_name}`;
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
        const applyFeedback = document.getElementById('apply-feedback');
        if (applyFeedback) applyFeedback.style.display = 'none';
        runBtn.disabled = true;
        runBtn.innerHTML = '⏳ Applying...';

        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema;
        let handlerType = ui.type;
        if (ui.fields && handlerType !== 'ketcher' && handlerType !== 'mapping_ui') handlerType = 'fields';

        const handler = UIRegistry[handlerType];
        if (!handler) {
            runBtn.disabled = false;
            runBtn.innerHTML = '💾 Apply';
            return;
        }

        let userInput = {};
        try {
            userInput = await handler.getData(ui);
        } catch(e) {
            errorContainer.innerText = e.message;
            errorContainer.style.display = 'block';
            runBtn.disabled = false;
            runBtn.innerHTML = '💾 Apply';
            return;
        }

        try {
            const data = await executeStepRun(plugin.step_id, userInput);
            if (!data.success) {
                errorContainer.innerText = data.error || "Unknown error";
                errorContainer.style.display = 'block';
                return;
            }

            state[plugin.step_id] = data.data;
            const resData = await fetchPipelineData();
            project = resData.project;
            plugins = resData.plugins;
            updateProjectUI(reactionStepSelect, deleteStepBtn, project);
            refreshTabs();

            if (applyFeedback) {
                applyFeedback.textContent = `✓ ${plugin.step_name} Applied & Saved!`;
                applyFeedback.style.display = 'inline-block';
                applyFeedback.style.opacity = '1';
                setTimeout(() => {
                    applyFeedback.style.opacity = '0';
                    setTimeout(() => { applyFeedback.style.display = 'none'; }, 300);
                }, 3000);
            }
        } catch (err) {
            errorContainer.innerText = err.message;
            errorContainer.style.display = 'block';
        } finally {
            runBtn.disabled = false;
            runBtn.innerHTML = '💾 Apply';
        }
    };

    setupProjectControls({
        reactionStepSelect,
        addStepBtn,
        deleteStepBtn,
        resetBtn,
        saveProjectBtn,
        openProjectBtn,
        fileOpenProjectInput,
        getActiveStepIndex: () => activeStepIndex,
        setActiveStepIndex: (idx) => { activeStepIndex = idx; },
        getProject: () => project,
        loadPipeline,
        refreshTabs,
        renderActiveStep
    });

    loadPipeline();
}

if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', boot);
} else {
    boot();
}
