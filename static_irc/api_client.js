/**
 * API client for Molecular Workbench IRC pipeline.
 */

export async function fetchPipelineData() {
    const [pluginsRes, stateRes, projectRes] = await Promise.all([
        fetch('/api/plugins'),
        fetch('/api/state'),
        fetch('/api/project')
    ]);
    const plugins = await pluginsRes.json();
    const state = await stateRes.json();
    const project = await projectRes.json();
    return { plugins, state, project };
}

export async function executeStepRun(stepId, userInput) {
    const res = await fetch(`/api/step/${stepId}/run`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify(userInput)
    });
    return await res.json();
}

export async function setGlobalSolvent(newSolv) {
    return await fetch('/api/environment', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ solvent: newSolv })
    });
}

export async function selectProjectStep(stepId) {
    const res = await fetch('/api/project/step/select', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ step_id: stepId })
    });
    return await res.json();
}

export async function addProjectStep(seed = true) {
    const res = await fetch('/api/project/step/add', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ seed })
    });
    return await res.json();
}

export async function deleteProjectStep(stepId) {
    const res = await fetch('/api/project/step/delete', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ step_id: stepId })
    });
    return await res.json();
}

export async function resetPipeline(resetAll) {
    const res = await fetch('/api/reset', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ all_steps: resetAll })
    });
    return await res.json();
}

export async function importProjectArchive(file) {
    const res = await fetch('/api/project/archive/import', {
        method: 'POST',
        headers: { 'Content-Type': 'application/octet-stream' },
        body: file
    });
    return await res.json();
}

export function exportProjectArchive() {
    window.location.href = '/api/project/archive/export';
}
