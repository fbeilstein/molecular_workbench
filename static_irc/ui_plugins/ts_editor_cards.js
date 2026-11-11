/**
 * UI control cards for the 3D Transition State Editor:
 * - Reaction Center Bonds Card (forming/breaking bonds)
 * - Fragment Positioning & In-place Refinement Card (MMFF / xTB)
 * - Strategy Selection & Optimization Runner Card
 */

export function createReactionCard({ rxnBonds, getPicked, onBondsChanged, el, styledButton, card, viewer }) {
    const rxnCard = card('Reaction Center Bonds');
    const bondList = el('div', { fontSize: '12px', color: '#444', marginBottom: '8px' });
    rxnCard.appendChild(bondList);

    function addBondRow(pair, type, color) {
        const [[f1, a1], [f2, a2]] = pair;
        const curD = viewer ? viewer.getAtomDistance(f1, a1, f2, a2) : 1.85;

        const row = el('div', { display: 'flex', alignItems: 'center', justifyContent: 'space-between', padding: '3px 0' });
        const titleSpan = el('div', { display: 'flex', alignItems: 'center', gap: '6px' });
        titleSpan.innerHTML = `<span style="color:${color};font-weight:600">${type}</span> R${f1 + 1}:${a1} ··· R${f2 + 1}:${a2}`;

        const dSpan = el('span', { fontSize: '11.5px', color: '#0284c7', fontWeight: 'bold' }, `${curD.toFixed(3)} Å`);
        titleSpan.appendChild(dSpan);

        const del = styledButton('✕', false, true);
        del.onclick = () => {
            const key = type.toLowerCase();
            rxnBonds[key] = rxnBonds[key].filter(p => p !== pair);
            onBondsChanged();
            refreshBondList();
        };
        row.appendChild(titleSpan);
        row.appendChild(del);
        bondList.appendChild(row);

        // Interactive Distance Slider & Numeric Input Control
        const ctrlRow = el('div', { display: 'flex', alignItems: 'center', gap: '8px', marginTop: '2px', marginBottom: '8px', paddingLeft: '4px' });
        const slider = el('input', { type: 'range', min: '1.20', max: '3.60', step: '0.01', value: curD.toFixed(2) });
        slider.style.width = '120px';
        slider.style.cursor = 'pointer';

        const numInput = el('input', { type: 'number', min: '1.0', max: '4.5', step: '0.02', value: curD.toFixed(2) });
        numInput.style.width = '55px';
        numInput.style.fontSize = '11.5px';
        numInput.style.padding = '2px 4px';
        numInput.style.border = '1px solid #cbd5e1';
        numInput.style.borderRadius = '4px';

        const updateDist = (val) => {
            const clamped = Math.max(1.0, Math.min(4.5, val));
            slider.value = clamped.toFixed(2);
            numInput.value = clamped.toFixed(2);
            dSpan.textContent = `${clamped.toFixed(3)} Å`;
            if (viewer) viewer.setAtomDistance(f1, a1, f2, a2, clamped);
            if (onBondsChanged) onBondsChanged();
        };

        slider.oninput = () => updateDist(parseFloat(slider.value));
        numInput.onchange = () => updateDist(parseFloat(numInput.value));

        const minusBtn = styledButton('-0.05', false, true);
        minusBtn.style.padding = '2px 6px';
        minusBtn.style.fontSize = '10.5px';
        minusBtn.onclick = () => updateDist(parseFloat(slider.value) - 0.05);

        const plusBtn = styledButton('+0.05', false, true);
        plusBtn.style.padding = '2px 6px';
        plusBtn.style.fontSize = '10.5px';
        plusBtn.onclick = () => updateDist(parseFloat(slider.value) + 0.05);

        ctrlRow.appendChild(slider);
        ctrlRow.appendChild(numInput);
        ctrlRow.appendChild(minusBtn);
        ctrlRow.appendChild(plusBtn);
        bondList.appendChild(ctrlRow);
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

export function createRefinementCard({ numReactants, onSelect, onRefine, onGizmoMode, onUndo, onReset, onSeed, el, styledButton, card }) {
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

    const resetCleanBtn = styledButton('↺ Reset Clean', false, true);
    resetCleanBtn.onclick = () => onReset();

    const seedBtn = styledButton('🌱 Seed from Prev Product', false, true);
    seedBtn.title = 'Seed exact 3D geometry of reactants from previous IRC product';
    seedBtn.onclick = () => onSeed && onSeed();

    recoveryRow.appendChild(undoBtn);
    recoveryRow.appendChild(resetCleanBtn);
    recoveryRow.appendChild(seedBtn);
    tuneCard.appendChild(recoveryRow);

    return { card: tuneCard, undoBtn, resetCleanBtn, seedBtn };
}

export function createOptimizationCard({
    strats, currentSolvent, onSolventChange, onRun, el, styledButton, card,
    reactionInfo, onComputeDefaults
}) {
    const actCard = card('Transition State Optimization');
    const solvLabel = el('div', { fontSize: '12px', fontWeight: '600', color: '#444', marginBottom: '4px' }, 'Solvent / Environment:');
    const solventSelect = el('select', {
        padding: '7px 10px', borderRadius: '6px', border: '1px solid #c7ccd6',
        width: '100%', marginBottom: '10px', fontSize: '12.5px', background: '#fff', color: '#222'
    });
    const solvOptions = [
        { val: 'water', label: 'Water (ALPB, ε=80.4)' },
        { val: 'gas', label: 'Gas Phase (Vacuum)' },
        { val: 'acetonitrile', label: 'Acetonitrile (ε=37.5)' },
        { val: 'ch2cl2', label: 'DCM (ε=8.9)' }
    ];
    solvOptions.forEach(opt => {
        const o = el('option', {}, opt.label);
        o.value = opt.val;
        solventSelect.appendChild(o);
    });
    solventSelect.value = currentSolvent || 'water';
    if (onSolventChange) {
        solventSelect.onchange = () => onSolventChange(solventSelect.value);
    }

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

    // Drive-specific options (Scan Start, End, Points)
    const driveBox = el('div', {
        display: stratSelect.value === 'drive' ? 'flex' : 'none',
        flexDirection: 'column',
        gap: '6px',
        background: '#f8fafc',
        border: '1px solid #e2e8f0',
        borderRadius: '6px',
        padding: '8px 10px',
        marginBottom: '10px'
    });

    const rec = reactionInfo?.recommended_scan || null;

    const driveHeader = el('div', { display: 'flex', justifyContent: 'space-between', alignItems: 'center' });
    const driveTitle = el('div', { fontSize: '11.5px', fontWeight: '700', color: '#475569' }, '1D Coordinate Scan Settings:');
    const btnResetDefaults = el('button', {
        background: 'none', border: 'none', color: '#3b82f6', cursor: 'pointer',
        fontSize: '11px', textDecoration: 'underline', padding: '0'
    }, '↺ Prefill Defaults');
    driveHeader.appendChild(driveTitle);
    driveHeader.appendChild(btnResetDefaults);

    const driveSubtitle = el('div', { fontSize: '10.5px', color: '#64748b', fontStyle: 'italic', marginBottom: '2px' });
    if (rec?.label) {
        driveSubtitle.textContent = `Target: ${rec.label}` + (rec.d0 ? ` (d0 = ${rec.d0} Å)` : '');
    }

    const driveRow = el('div', { display: 'flex', alignItems: 'center', gap: '8px' });

    const scanStartBox = el('div', { flex: '1' });
    const lblStart = el('div', { fontSize: '11px', color: '#64748b' }, 'Start (Å):');
    const inpStart = document.createElement('input');
    inpStart.type = 'number';
    inpStart.step = '0.02';
    inpStart.min = '0.2';
    inpStart.max = '10.0';
    inpStart.placeholder = rec?.start != null ? String(rec.start) : 'auto';
    inpStart.value = rec?.start != null ? String(rec.start) : '';
    inpStart.style.width = '100%';
    inpStart.style.padding = '4px 6px';
    inpStart.style.fontSize = '12px';
    inpStart.style.border = '1px solid #cbd5e1';
    inpStart.style.borderRadius = '4px';
    inpStart.style.boxSizing = 'border-box';
    scanStartBox.appendChild(lblStart);
    scanStartBox.appendChild(inpStart);

    const scanEndBox = el('div', { flex: '1' });
    const lblEnd = el('div', { fontSize: '11px', color: '#64748b' }, 'End (Å):');
    const inpEnd = document.createElement('input');
    inpEnd.type = 'number';
    inpEnd.step = '0.02';
    inpEnd.min = '0.2';
    inpEnd.max = '10.0';
    inpEnd.placeholder = rec?.end != null ? String(rec.end) : 'auto';
    inpEnd.value = rec?.end != null ? String(rec.end) : '';
    inpEnd.style.width = '100%';
    inpEnd.style.padding = '4px 6px';
    inpEnd.style.fontSize = '12px';
    inpEnd.style.border = '1px solid #cbd5e1';
    inpEnd.style.borderRadius = '4px';
    inpEnd.style.boxSizing = 'border-box';
    scanEndBox.appendChild(lblEnd);
    scanEndBox.appendChild(inpEnd);

    const scanPtsBox = el('div', { width: '65px' });
    const lblPts = el('div', { fontSize: '11px', color: '#64748b' }, 'Points:');
    const inpPts = document.createElement('input');
    inpPts.type = 'number';
    inpPts.step = '1';
    inpPts.min = '5';
    inpPts.max = '60';
    inpPts.placeholder = rec?.points != null ? String(rec.points) : '21';
    inpPts.value = rec?.points != null ? String(rec.points) : '21';
    inpPts.style.width = '100%';
    inpPts.style.padding = '4px 6px';
    inpPts.style.fontSize = '12px';
    inpPts.style.border = '1px solid #cbd5e1';
    inpPts.style.borderRadius = '4px';
    inpPts.style.boxSizing = 'border-box';
    scanPtsBox.appendChild(lblPts);
    scanPtsBox.appendChild(inpPts);

    driveRow.appendChild(scanStartBox);
    driveRow.appendChild(scanEndBox);
    driveRow.appendChild(scanPtsBox);
    driveBox.appendChild(driveHeader);
    if (rec?.label) driveBox.appendChild(driveSubtitle);
    driveBox.appendChild(driveRow);

    const applyDefaults = (defaults) => {
        if (!defaults) return;
        if (defaults.start != null) {
            inpStart.value = String(defaults.start);
            inpStart.placeholder = String(defaults.start);
        }
        if (defaults.end != null) {
            inpEnd.value = String(defaults.end);
            inpEnd.placeholder = String(defaults.end);
        }
        if (defaults.points != null) {
            inpPts.value = String(defaults.points);
            inpPts.placeholder = String(defaults.points);
        }
        if (defaults.label) {
            driveSubtitle.textContent = `Target: ${defaults.label}` + (defaults.d0 ? ` (d0 = ${defaults.d0} Å)` : '');
            if (!driveSubtitle.parentNode) {
                driveBox.insertBefore(driveSubtitle, driveRow);
            }
        }
    };

    btnResetDefaults.onclick = async (e) => {
        e.preventDefault();
        if (onComputeDefaults) {
            btnResetDefaults.textContent = '⏳ Calculating...';
            try {
                const defaults = await onComputeDefaults();
                applyDefaults(defaults);
            } catch (err) {
                console.error("Scan default error:", err);
            } finally {
                btnResetDefaults.textContent = '↺ Prefill Defaults';
            }
        }
    };

    stratSelect.addEventListener('change', async () => {
        const isDrive = stratSelect.value === 'drive';
        driveBox.style.display = isDrive ? 'flex' : 'none';
        if (isDrive && (!inpStart.value || !inpEnd.value) && onComputeDefaults) {
            try {
                const defaults = await onComputeDefaults();
                applyDefaults(defaults);
            } catch (err) {}
        }
    });

    const runBtn = styledButton('⚡ Run / Refine TS', true);
    const optStatus = el('div', { fontWeight: '600', color: '#4f6ef7', marginTop: '8px', fontSize: '13px' });

    actCard.appendChild(solvLabel);
    actCard.appendChild(solventSelect);
    actCard.appendChild(stratLabel);
    actCard.appendChild(stratSelect);
    actCard.appendChild(driveBox);
    actCard.appendChild(runBtn);

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

    runBtn.onclick = () => onRun(stratSelect.value, solventSelect.value);

    const prefillDefaults = (startVal, endVal, ptsVal) => {
        if (startVal != null) inpStart.value = String(startVal);
        if (endVal != null) inpEnd.value = String(endVal);
        if (ptsVal != null) inpPts.value = String(ptsVal);
    };

    return {
        card: actCard, profCanvas, logBox, logHeader, runBtn, optStatus, solventSelect, stratSelect,
        getScanStart: () => {
            const v = parseFloat(inpStart.value);
            return isNaN(v) ? null : v;
        },
        getScanEnd: () => {
            const v = parseFloat(inpEnd.value);
            return isNaN(v) ? null : v;
        },
        getScanPoints: () => {
            const v = parseInt(inpPts.value, 10);
            return isNaN(v) ? null : v;
        },
        prefillDefaults,
        applyDefaults
    };
}
