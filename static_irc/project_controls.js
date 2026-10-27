/**
 * Event handlers and bindings for multi-step project actions.
 */

import {
    selectProjectStep,
    addProjectStep,
    deleteProjectStep,
    resetPipeline,
    importProjectArchive,
    exportProjectArchive
} from './api_client.js';

export function setupProjectControls({
    reactionStepSelect,
    addStepBtn,
    deleteStepBtn,
    resetBtn,
    saveProjectBtn,
    openProjectBtn,
    fileOpenProjectInput,
    getActiveStepIndex,
    setActiveStepIndex,
    getProject,
    loadPipeline,
    refreshTabs,
    renderActiveStep
}) {
    if (reactionStepSelect) {
        reactionStepSelect.onchange = async () => {
            const targetId = reactionStepSelect.value;
            const curTabIndex = getActiveStepIndex();
            try {
                const r = await selectProjectStep(targetId);
                if (r.success) await loadPipeline(curTabIndex);
                else alert('Failed to switch step: ' + (r.error || 'Unknown error'));
            } catch (err) {
                alert('Request failed: ' + err.message);
            }
        };
    }

    if (addStepBtn) {
        addStepBtn.onclick = async () => {
            addStepBtn.disabled = true;
            try {
                const r = await addProjectStep(true);
                if (r.success) {
                    await loadPipeline();
                    setActiveStepIndex(0);
                    refreshTabs();
                    renderActiveStep();
                } else alert('Failed to add step: ' + (r.error || 'Unknown error'));
            } catch (err) {
                alert('Add step request failed: ' + err.message);
            } finally {
                addStepBtn.disabled = false;
            }
        };
    }

    if (deleteStepBtn) {
        deleteStepBtn.onclick = async () => {
            const curStep = reactionStepSelect ? reactionStepSelect.value : 'step_1';
            if (!confirm(`Are you sure you want to delete ${curStep}?`)) return;
            deleteStepBtn.disabled = true;
            try {
                const r = await deleteProjectStep(curStep);
                if (r.success) await loadPipeline();
                else alert('Failed to delete step: ' + (r.error || 'Unknown error'));
            } catch (err) {
                alert('Delete step request failed: ' + err.message);
            } finally {
                deleteStepBtn.disabled = false;
            }
        };
    }

    if (resetBtn) {
        resetBtn.onclick = async () => {
            const proj = getProject();
            const totalSteps = proj.steps?.length || 1;
            let resetAll = false;
            if (totalSteps > 1) {
                resetAll = confirm("Click OK to reset the ENTIRE project (all reaction steps), or Cancel to reset ONLY the currently active step.");
                if (!resetAll && !confirm(`Reset only the active step (${proj.active_step_id})?`)) return;
            } else if (!confirm("Reset the current reaction step?")) return;

            await resetPipeline(resetAll);
            await loadPipeline();
        };
    }

    if (saveProjectBtn) saveProjectBtn.onclick = exportProjectArchive;

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
                const r = await importProjectArchive(file);
                if (r.success) {
                    alert(r.message || 'Project successfully loaded!');
                    await loadPipeline(0);
                } else alert('Failed to load project: ' + (r.error || 'Unknown error'));
            } catch (err) {
                alert('Project import failed: ' + err.message);
            } finally {
                openProjectBtn.disabled = false;
                openProjectBtn.textContent = '📂 Open Project';
            }
        };
    }
}
