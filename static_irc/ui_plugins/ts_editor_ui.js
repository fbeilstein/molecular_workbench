import { TSEditorViewer } from './ts_editor_viewer.js';
import { TSEditorActions } from './ts_editor_actions.js';
import { updateMolBlockPositions } from './ts_sdf_parser.js';
import { createReactionCard, createRefinementCard, createOptimizationCard } from './ts_editor_cards.js';

function el(tag, style = {}, text = null) {
    const e = document.createElement(tag);
    Object.assign(e.style, style);
    if (text !== null) e.textContent = text;
    return e;
}

function styledButton(text, primary = false, small = false) {
    return el('button', {
        padding: small ? '5px 9px' : '9px 14px',
        margin: '4px 6px 4px 0',
        borderRadius: '6px',
        cursor: 'pointer',
        fontSize: small ? '11.5px' : '13px',
        border: primary ? 'none' : '1px solid #c7ccd6',
        fontWeight: primary ? '600' : '500',
        background: primary ? 'linear-gradient(135deg,#4f6ef7,#7b4ff7)' : '#fff',
        color: primary ? '#fff' : '#333',
        boxShadow: primary ? '0 2px 8px rgba(79,110,247,.35)' : 'none'
    }, text);
}

function card(title) {
    const c = el('div', {
        background: '#fff', border: '1px solid #e1e4ea', borderRadius: '10px',
        padding: '12px 14px', boxShadow: '0 2px 8px rgba(0,0,0,.03)',
        flexShrink: '0'
    });
    if (title) {
        c.appendChild(el('div', { fontWeight: '700', fontSize: '13px', color: '#1a1f36', marginBottom: '8px' }, title));
    }
    return c;
}

let activeViewer = null;
let activeActions = null;
let state = { molblocks: [], rxnBonds: { forming: [], breaking: [] }, overridden: false, picked: [] };

export function renderTSEditor(arg1, arg2) {
    let container, ui;
    if (arg1 && (arg1.nodeType || arg1 instanceof HTMLElement)) {
        container = arg1;
        ui = arg2 || {};
    } else {
        ui = arg1 || {};
        container = arg2;
    }
    if (!container || !container.appendChild) {
        console.error("Invalid container provided to renderTSEditor:", arg1, arg2);
        return;
    }
    container.innerHTML = '';
    state.molblocks = [...(ui.reactants_3d || [])];
    const info = ui.reaction_info || {};
    state.rxnBonds = { forming: [...(info.forming || [])], breaking: [...(info.breaking || [])] };
    state.overridden = false;
    state.picked = [];

    const wrapper = el('div', {
        display: 'flex', flexDirection: 'row', width: '100%', height: 'calc(100vh - 120px)',
        minHeight: '620px', gap: '14px', boxSizing: 'border-box'
    });
    const viewerContainer = el('div', {
        flex: '1 1 65%', minWidth: '450px', height: '100%', position: 'relative',
        background: '#0f111a', borderRadius: '12px', overflow: 'hidden'
    });
    const controlsContainer = el('div', {
        flex: '1 1 35%', minWidth: '360px', height: '100%', display: 'flex',
        flexDirection: 'column', gap: '10px', overflowY: 'auto', paddingRight: '4px'
    });

    const distanceHUD = el('div', {
        position: 'absolute', top: '12px', left: '12px', background: 'rgba(15,17,26,0.85)',
        color: '#e2e8f0', padding: '6px 10px', borderRadius: '6px', fontSize: '12px',
        fontFamily: 'monospace', zIndex: '10', display: 'none', border: '1px solid #2d3748'
    });
    viewerContainer.appendChild(distanceHUD);

    wrapper.appendChild(viewerContainer);
    wrapper.appendChild(controlsContainer);
    container.appendChild(wrapper);

    function updateLiveHUD() {
        if (state.picked.length < 2) { distanceHUD.style.display = 'none'; return; }
        const p1 = activeViewer.getAtomWorldPos(state.picked[0].fi, state.picked[0].ai);
        const p2 = activeViewer.getAtomWorldPos(state.picked[1].fi, state.picked[1].ai);
        if (!p1 || !p2) { distanceHUD.style.display = 'none'; return; }
        distanceHUD.textContent = `Distance: ${p1.distanceTo(p2).toFixed(3)} Å (${state.picked[0].sym}···${state.picked[1].sym})`;
        distanceHUD.style.display = 'block';
    }

    activeViewer = new TSEditorViewer(viewerContainer, (fi, ai, sym, mesh) => {
        const existIdx = state.picked.findIndex(p => p.fi === fi && p.ai === ai);
        if (existIdx >= 0) {
            mesh.material.emissive.setHex(0x000000);
            state.picked.splice(existIdx, 1);
        } else {
            if (state.picked.length >= 2) {
                const old = state.picked.shift();
                const oldMesh = activeViewer.moleculeGroups[old.fi]?.userData.atomMeshes[old.ai];
                if (oldMesh) oldMesh.material.emissive.setHex(0x000000);
            }
            state.picked.push({ fi, ai, sym });
            mesh.material.emissive.setHex(0x3355ff);
        }
        updateLiveHUD();
        rxnCard.refreshBondList();
    }, () => {
        state.overridden = true;
        activeViewer.updateDynamicBonds(state.rxnBonds);
        updateLiveHUD();
    });

    const rxnCard = createReactionCard({
        rxnBonds: state.rxnBonds,
        getPicked: () => state.picked,
        onBondsChanged: () => {
            state.overridden = true;
            activeViewer.updateDynamicBonds(state.rxnBonds);
        },
        el, styledButton, card
    });
    controlsContainer.appendChild(rxnCard.card);

    let refineCard, optCard;
    refineCard = createRefinementCard({
        numReactants: state.molblocks.length,
        onSelect: fi => activeViewer.selectGroup(fi),
        onRefine: (fi, m) => activeActions.refineFragment(fi, m),
        onGizmoMode: m => activeViewer.setGizmoMode(m),
        onUndo: () => activeActions.popHistory(() => rxnCard.refreshBondList()),
        onReset: () => activeActions.resetToClean(ui.clean_reactants_3d),
        onSeed: () => activeActions.seedFromPrevious(),
        el, styledButton, card
    });
    controlsContainer.appendChild(refineCard.card);

    const strats = [
        { id: 'genetic_swarm', name: '🧬 Genetic Swarm: Multi-Cycle SE(3) Rigid-Body Evolution (Recommended)' },
        { id: 'lego', name: '🧱 Dual-Endpoint IDPP: Universal Reaction Interpolation → Sella' },
        { id: 'neb', name: '⛰️ Climbing Image NEB (CI-NEB with IDPP)' },
        { id: 'manual_valley', name: 'Quantum Valley: Refine Crude Manual 3D Guess' },
        { id: 'auto_dock', name: '🎯 Auto-Dock: 32-Start Geometric Align → Sella' },
        { id: 'orca_dft', name: '🔬 ORCA 6.1: High-Level r2SCAN-3c DFT Refinement' },
        { id: 'drive', name: 'Constrained Scan (1D Drive) → Sella' },
        { id: 'dimer', name: 'Dimer Method (ase.mep.dimer)' },
        { id: 'rits', name: '⚡ RitS: Flow-Matching ML (Requires Checkpoint)' }
    ];

    optCard = createOptimizationCard({
        strats,
        currentSolvent: state.solvent || ui.solvent || 'water',
        onSolventChange: async (s) => {
            state.solvent = s;
            const g = document.getElementById('global-solvent-select');
            if (g) g.value = s;
            try {
                await fetch('/api/environment', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ solvent: s })
                });
            } catch (e) {}
        },
        onRun: (strat, solv) => activeActions.runOptimization(strat, solv),
        el, styledButton, card
    });
    controlsContainer.appendChild(optCard.card);
    controlsContainer.appendChild(optCard.profCanvas);
    if (optCard.logHeader) controlsContainer.appendChild(optCard.logHeader);
    controlsContainer.appendChild(optCard.logBox);

    activeActions = new TSEditorActions(state, optCard, refineCard, activeViewer);
    activeActions.onGeometryUpdated = () => updateLiveHUD();
    activeViewer.buildMolecules(state.molblocks);
    activeViewer.updateDynamicBonds(state.rxnBonds);

    window.addEventListener('keydown', e => {
        if (['INPUT', 'TEXTAREA'].includes(document.activeElement?.tagName)) return;
        if (e.key.toLowerCase() === 't') activeViewer.setGizmoMode('translate');
        if (e.key.toLowerCase() === 'r') activeViewer.setGizmoMode('rotate');
        if (e.key === 'Escape') activeViewer.setGizmoMode('detach');
    });
}

export const TSEditorUI = {
    render: renderTSEditor,
    getData: async () => ({
        docked_molblocks: TSEditorUI.extractData(),
        reaction_override: state.overridden ? state.rxnBonds : null
    }),
    extractData: () => {
        if (!activeViewer) return state.molblocks;
        return state.molblocks.map((mb, fi) => {
            const g = activeViewer.moleculeGroups[fi];
            return g ? updateMolBlockPositions(mb, g) : mb;
        });
    }
};
