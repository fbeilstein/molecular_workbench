/**
 * UI control cards for the 3D Transition State Editor:
 * - Reaction Center Bonds Card (forming/breaking bonds)
 * - Fragment Positioning & In-place Refinement Card (MMFF / xTB)
 * - Strategy Selection & Optimization Runner Card
 */

export function createReactionCard({ rxnBonds, getPicked, onBondsChanged, el, styledButton, card }) {
    const rxnCard = card('Reaction Center Bonds');
    const bondList = el('div', { fontSize: '12px', color: '#444', marginBottom: '8px' });
    rxnCard.appendChild(bondList);

    function addBondRow(pair, type, color) {
        const [[f1, a1], [f2, a2]] = pair;
        const row = el('div', { display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '3px 0' });
        row.innerHTML = `<span style="color:${color};font-weight:600">${type}</span> R${f1 + 1}:${a1} ··· R${f2 + 1}:${a2}`;
        const del = styledButton('✕', false, true);
        del.onclick = () => {
            const key = type.toLowerCase();
            rxnBonds[key] = rxnBonds[key].filter(p => p !== pair);
            onBondsChanged();
            refreshBondList();
        };
        row.appendChild(del);
        bondList.appendChild(row);
    }

    function refreshBondList() {
        bondList.innerHTML = '';
        rxnBonds.forming.forEach(p => addBondRow(p, 'Forming', '#18a558'));
        rxnBonds.breaking.forEach(p => addBondRow(p, 'Breaking', '#e0245e'));

        if (!rxnBonds.forming.length && !rxnBonds.breaking.length) {
            bondList.textContent = 'No reacting bonds defined. Pick two atoms to assign.';
        }
    }
    refreshBondList();

    const pickRow = el('div', { display: 'flex', gap: '6px', flexWrap: 'wrap' });
    const formBtn = styledButton('+ Mark Forming', false, true);
    formBtn.onclick = () => {
        const picked = getPicked();
        if (picked.length !== 2) return alert('Select 2 atoms first by clicking them.');
        rxnBonds.forming.push([[picked[0].fi, picked[0].ai], [picked[1].fi, picked[1].ai]]);
        onBondsChanged();
        refreshBondList();
    };

    const breakBtn = styledButton('+ Mark Breaking', false, true);
    breakBtn.onclick = () => {
        const picked = getPicked();
        if (picked.length !== 2) return alert('Select 2 atoms first by clicking them.');
        rxnBonds.breaking.push([[picked[0].fi, picked[0].ai], [picked[1].fi, picked[1].ai]]);
        onBondsChanged();
        refreshBondList();
    };

    pickRow.appendChild(formBtn);
    pickRow.appendChild(breakBtn);
    rxnCard.appendChild(pickRow);

    return { card: rxnCard, refreshBondList };
}

export function createRefinementCard({ numReactants, onSelect, onRefine, onGizmoMode, onUndo, onReset, el, styledButton, card }) {
    const tuneCard = card('Positioning & Shape Refinement');

    for (let idx = 0; idx < numReactants; idx++) {
        const row = el('div', { display: 'flex', alignItems: 'center', gap: '6px', marginBottom: '6px', flexWrap: 'wrap' });
        const selBtn = styledButton(`Select Reactant ${idx + 1}`, false, true);
        selBtn.onclick = () => onSelect(idx);

        const mmffBtn = styledButton('Clean (MMFF)', false, true);
        mmffBtn.title = 'Relax bond lengths and angles with MMFF94 in-place without losing position';
        mmffBtn.onclick = () => onRefine(idx, 'mmff');

        const xtbBtn = styledButton('Relax (xTB)', false, true);
        xtbBtn.title = 'Quantum geometry relaxation with GFN2-xTB in-place without losing position';
        xtbBtn.onclick = () => onRefine(idx, 'xtb');

        row.appendChild(selBtn);
        row.appendChild(mmffBtn);
        row.appendChild(xtbBtn);
        tuneCard.appendChild(row);
    }

    const gizmoRow = el('div', { marginTop: '6px' });
    const transBtn = styledButton('Move (T)', false, true);
    transBtn.onclick = () => onGizmoMode('translate');
    const rotBtn = styledButton('Rotate (R)', false, true);
    rotBtn.onclick = () => onGizmoMode('rotate');
    const hideBtn = styledButton('Hide gizmo', false, true);
    hideBtn.onclick = () => onGizmoMode('detach');

    gizmoRow.appendChild(transBtn);
    gizmoRow.appendChild(rotBtn);
    gizmoRow.appendChild(hideBtn);
    tuneCard.appendChild(gizmoRow);

    const recoveryRow = el('div', { marginTop: '8px', paddingTop: '8px', borderTop: '1px dashed #e1e4ea', display: 'flex', gap: '6px', flexWrap: 'wrap' });
    const undoBtn = styledButton('↩ Undo Last Change', false, true);
    undoBtn.disabled = true;
    undoBtn.onclick = () => onUndo();

    const resetCleanBtn = styledButton('↺ Reset to Clean Reactants', false, true);
    resetCleanBtn.onclick = () => onReset();

    recoveryRow.appendChild(undoBtn);
    recoveryRow.appendChild(resetCleanBtn);
    tuneCard.appendChild(recoveryRow);

    return { card: tuneCard, undoBtn, resetCleanBtn };
}

export function createOptimizationCard({ strats, showNextPose, onRun, onNextPose, el, styledButton, card }) {
    const actCard = card('Transition State Optimization');
    const stratLabel = el('div', { fontSize: '12px', fontWeight: '600', color: '#444', marginBottom: '4px' }, 'Optimization Strategy:');
    const stratSelect = el('select', {
        padding: '7px 10px', borderRadius: '6px', border: '1px solid #c7ccd6',
        width: '100%', marginBottom: '10px', fontSize: '12.5px', background: '#fff', color: '#222'
    });

    strats.forEach(s => {
        const opt = el('option', {}, s.name);
        opt.value = s.id;
        stratSelect.appendChild(opt);
    });

    const runBtn = styledButton('⚡ Run / Refine TS', true);
    const nextPoseBtn = styledButton('↻ Next docking pose');
    const optStatus = el('div', { fontWeight: '600', color: '#4f6ef7', marginTop: '8px', fontSize: '13px' });

    actCard.appendChild(stratLabel);
    actCard.appendChild(stratSelect);
    actCard.appendChild(runBtn);
    if (showNextPose) actCard.appendChild(nextPoseBtn);

    actCard.appendChild(el('div', { fontSize: '11.5px', color: '#777', marginTop: '6px' },
        'RitS: Generative flow-matching TS predictor. ' +
        'Lego: R+P dual-endpoint scaffold superposition. ' +
        'Quantum Valley: relaxes crude manual placements. ' +
        'ORCA: high-level DFT (r2SCAN-3c) refinement.'));
    actCard.appendChild(optStatus);

    const profCanvas = el('canvas', { width: '100%', border: '1px solid #e1e4ea', borderRadius: '8px', background: '#fff', display: 'none', flexShrink: '0' });
    profCanvas.width = 380;
    profCanvas.height = 150;

    const logHeader = el('div', {
        display: 'flex', justifyContent: 'space-between', alignItems: 'center',
        marginTop: '6px', marginBottom: '2px', flexShrink: '0'
    });
    const logTitle = el('span', { fontWeight: '700', fontSize: '12px', color: '#2d3748' }, 'Optimization Console Log:');
    const logClearBtn = styledButton('Clear Log', false, true);
    logClearBtn.onclick = () => { logBox.value = ''; };
    logHeader.appendChild(logTitle);
    logHeader.appendChild(logClearBtn);

    const logBox = el('textarea', {
        width: '100%', height: '240px', minHeight: '200px', flexShrink: '0',
        fontFamily: 'monospace', fontSize: '11.5px', lineHeight: '1.45',
        backgroundColor: '#1e2130', color: '#d6e2ff', border: '1px solid #2d3348',
        borderRadius: '8px', padding: '10px', boxSizing: 'border-box',
        whiteSpace: 'pre-wrap', wordBreak: 'break-word', overflowX: 'hidden', overflowY: 'auto',
        resize: 'vertical'
    });
    logBox.readOnly = true;
    logBox.value = 'Ready.';

    runBtn.onclick = () => onRun(stratSelect.value);
    nextPoseBtn.onclick = () => onNextPose();

    return { card: actCard, profCanvas, logBox, logHeader, runBtn, nextPoseBtn, optStatus };
}
