export function buildCascadeLayout(container, ui, options = {}) {
    const root = document.createElement('div');
    root.style.display = 'grid';
    root.style.gridTemplateColumns = '370px 1fr';
    root.style.gap = '16px';
    root.style.minHeight = '620px';

    // ============ LEFT PANEL ============
    const leftPanel = document.createElement('div');
    leftPanel.style.display = 'flex';
    leftPanel.style.flexDirection = 'column';
    leftPanel.style.gap = '12px';

    // 1. Sequence Card
    const seqCard = document.createElement('div');
    seqCard.style.padding = '12px';
    seqCard.style.background = '#f8fafc';
    seqCard.style.border = '1px solid #e2e8f0';
    seqCard.style.borderRadius = '8px';

    const seqTitle = document.createElement('div');
    seqTitle.style.fontWeight = '700';
    seqTitle.style.fontSize = '13.5px';
    seqTitle.style.color = '#1e293b';
    seqTitle.style.marginBottom = '8px';
    seqTitle.textContent = 'Reaction Cascade Steps';
    seqCard.appendChild(seqTitle);

    const seqList = document.createElement('div');
    seqList.style.display = 'flex';
    seqList.style.flexDirection = 'column';
    seqList.style.gap = '6px';
    seqCard.appendChild(seqList);
    leftPanel.appendChild(seqCard);

    // 2. Junction Card
    const juncCard = document.createElement('div');
    juncCard.style.padding = '12px';
    juncCard.style.background = '#f8fafc';
    juncCard.style.border = '1px solid #e2e8f0';
    juncCard.style.borderRadius = '8px';

    const juncTitle = document.createElement('div');
    juncTitle.style.fontWeight = '700';
    juncTitle.style.fontSize = '13.5px';
    juncTitle.style.color = '#1e293b';
    juncTitle.style.marginBottom = '4px';
    juncTitle.textContent = 'Junction Intermediate Alignment';
    juncCard.appendChild(juncTitle);

    const juncDesc = document.createElement('div');
    juncDesc.style.fontSize = '11.5px';
    juncDesc.style.color = '#64748b';
    juncDesc.style.marginBottom = '8px';
    juncDesc.textContent = 'Toggle Rotamer-Aware Scaffold (preserves flexible dihedrals) or Rigid All-Atom.';
    juncCard.appendChild(juncDesc);

    const juncList = document.createElement('div');
    juncList.style.display = 'flex';
    juncList.style.flexDirection = 'column';
    juncList.style.gap = '8px';
    juncCard.appendChild(juncList);
    leftPanel.appendChild(juncCard);

    // 3. Bridge Frames Card
    const bridgeCard = document.createElement('div');
    bridgeCard.style.padding = '12px';
    bridgeCard.style.background = '#f8fafc';
    bridgeCard.style.border = '1px solid #e2e8f0';
    bridgeCard.style.borderRadius = '8px';

    const bTitleRow = document.createElement('div');
    bTitleRow.style.display = 'flex';
    bTitleRow.style.justifyContent = 'space-between';
    bTitleRow.style.marginBottom = '6px';

    const bTitle = document.createElement('span');
    bTitle.style.fontWeight = '700';
    bTitle.style.fontSize = '13.5px';
    bTitle.textContent = 'Bridge Transition Frames';

    const bVal = document.createElement('span');
    bVal.style.fontSize = '12px';
    bVal.style.fontWeight = '600';
    bVal.style.color = '#0284c7';
    bVal.textContent = `${options.initialBridgeFrames || 8} Frames`;

    bTitleRow.appendChild(bTitle);
    bTitleRow.appendChild(bVal);
    bridgeCard.appendChild(bTitleRow);

    const bSlider = document.createElement('input');
    bSlider.type = 'range';
    bSlider.min = '0';
    bSlider.max = '20';
    bSlider.step = '1';
    bSlider.value = (options.initialBridgeFrames || 8).toString();
    bSlider.style.width = '100%';
    bridgeCard.appendChild(bSlider);
    leftPanel.appendChild(bridgeCard);

    // 4. Actions
    const actionRow = document.createElement('div');
    actionRow.style.display = 'flex';
    actionRow.style.gap = '8px';

    const restitchBtn = document.createElement('button');
    restitchBtn.className = 'primary';
    restitchBtn.style.flex = '1';
    restitchBtn.style.padding = '8px 12px';
    restitchBtn.style.fontWeight = '600';
    restitchBtn.textContent = '⚡ Re-Stitch Pathway';

    actionRow.appendChild(restitchBtn);
    leftPanel.appendChild(actionRow);

    // 5. Metrics Box
    const metricsBox = document.createElement('div');
    metricsBox.style.display = 'grid';
    metricsBox.style.gridTemplateColumns = '1fr 1fr';
    metricsBox.style.gap = '8px';
    metricsBox.style.padding = '10px';
    metricsBox.style.background = '#ffffff';
    metricsBox.style.border = '1px solid #e2e8f0';
    metricsBox.style.borderRadius = '8px';

    function createMetric(title) {
        const m = document.createElement('div');
        const mt = document.createElement('div');
        mt.style.fontSize = '11px';
        mt.style.color = '#64748b';
        mt.textContent = title;
        const mv = document.createElement('div');
        mv.style.fontSize = '14px';
        mv.style.fontWeight = '700';
        mv.style.color = '#1e293b';
        mv.textContent = '--';
        m.appendChild(mt);
        m.appendChild(mv);
        return { el: m, val: mv };
    }

    const mFrames = createMetric('Total Frames');
    const mSteps = createMetric('Total Steps');
    const mBar = createMetric('Max Barrier ΔE‡');
    const mRxn = createMetric('Net ΔE rxn');

    metricsBox.appendChild(mFrames.el);
    metricsBox.appendChild(mSteps.el);
    metricsBox.appendChild(mBar.el);
    metricsBox.appendChild(mRxn.el);
    leftPanel.appendChild(metricsBox);

    root.appendChild(leftPanel);

    // ============ RIGHT PANEL ============
    const rightPanel = document.createElement('div');
    rightPanel.style.display = 'flex';
    rightPanel.style.flexDirection = 'column';
    rightPanel.style.gap = '8px';

    // Top View Mode Bar
    const viewModeBar = document.createElement('div');
    viewModeBar.style.display = 'flex';
    viewModeBar.style.justifyContent = 'space-between';
    viewModeBar.style.alignItems = 'center';
    viewModeBar.style.background = '#ffffff';
    viewModeBar.style.padding = '6px 10px';
    viewModeBar.style.border = '1px solid #e2e8f0';
    viewModeBar.style.borderRadius = '8px';

    const modeToggleBox = document.createElement('div');
    modeToggleBox.style.display = 'flex';
    modeToggleBox.style.gap = '6px';

    const btnModeCascade = document.createElement('button');
    btnModeCascade.style.padding = '5px 12px';
    btnModeCascade.style.fontSize = '12px';
    btnModeCascade.style.fontWeight = '600';
    btnModeCascade.style.borderRadius = '5px';
    btnModeCascade.style.background = '#0284c7';
    btnModeCascade.style.color = '#ffffff';
    btnModeCascade.style.border = 'none';
    btnModeCascade.textContent = '🎬 Full Cascade';

    const btnModeOverlap = document.createElement('button');
    btnModeOverlap.style.padding = '5px 12px';
    btnModeOverlap.style.fontSize = '12px';
    btnModeOverlap.style.fontWeight = '600';
    btnModeOverlap.style.borderRadius = '5px';
    btnModeOverlap.style.background = '#f1f5f9';
    btnModeOverlap.style.color = '#334155';
    btnModeOverlap.style.border = '1px solid #cbd5e1';
    btnModeOverlap.textContent = '🔍 Boundary Overlap';

    modeToggleBox.appendChild(btnModeCascade);
    modeToggleBox.appendChild(btnModeOverlap);
    viewModeBar.appendChild(modeToggleBox);

    // Overlap Sub-Controls (shown when in Boundary Overlap mode)
    const overlapSubControls = document.createElement('div');
    overlapSubControls.style.display = 'none';
    overlapSubControls.style.alignItems = 'center';
    overlapSubControls.style.gap = '8px';

    const juncSelect = document.createElement('select');
    juncSelect.style.padding = '4px 8px';
    juncSelect.style.fontSize = '11.5px';
    juncSelect.style.borderRadius = '4px';
    juncSelect.style.border = '1px solid #cbd5e1';

    const alignTypeBtn = document.createElement('button');
    alignTypeBtn.style.padding = '4px 8px';
    alignTypeBtn.style.fontSize = '11.5px';
    alignTypeBtn.style.borderRadius = '4px';
    alignTypeBtn.style.border = '1px solid #cbd5e1';
    alignTypeBtn.style.background = '#ffffff';
    alignTypeBtn.textContent = 'Aligned (Fit)';

    const rmsdBadge = document.createElement('span');
    rmsdBadge.style.fontSize = '11.5px';
    rmsdBadge.style.fontWeight = '600';
    rmsdBadge.style.color = '#0284c7';
    rmsdBadge.style.background = '#e0f2fe';
    rmsdBadge.style.padding = '3px 8px';
    rmsdBadge.style.borderRadius = '4px';
    rmsdBadge.textContent = 'RMSD: --';

    overlapSubControls.appendChild(juncSelect);
    overlapSubControls.appendChild(alignTypeBtn);
    overlapSubControls.appendChild(rmsdBadge);
    viewModeBar.appendChild(overlapSubControls);

    rightPanel.appendChild(viewModeBar);

    const viewerBox = document.createElement('div');
    viewerBox.style.width = '100%';
    viewerBox.style.height = '390px';
    viewerBox.style.background = '#0f172a';
    viewerBox.style.borderRadius = '8px';
    viewerBox.style.position = 'relative';

    const hudOverlay = document.createElement('div');
    hudOverlay.style.position = 'absolute';
    hudOverlay.style.top = '10px';
    hudOverlay.style.left = '12px';
    hudOverlay.style.background = 'rgba(15, 23, 42, 0.75)';
    hudOverlay.style.color = '#f8fafc';
    hudOverlay.style.padding = '4px 10px';
    hudOverlay.style.borderRadius = '6px';
    hudOverlay.style.fontSize = '12px';
    hudOverlay.style.fontFamily = 'monospace';
    hudOverlay.style.zIndex = '10';
    hudOverlay.textContent = 'Multi-Step Pathway Viewer';
    viewerBox.appendChild(hudOverlay);
    rightPanel.appendChild(viewerBox);

    const ctrlRow = document.createElement('div');
    ctrlRow.style.display = 'flex';
    ctrlRow.style.alignItems = 'center';
    ctrlRow.style.gap = '8px';

    const playBtn = document.createElement('button');
    playBtn.style.padding = '6px 12px';
    playBtn.style.fontWeight = '600';
    playBtn.textContent = '▶ Play';

    const prevBtn = document.createElement('button');
    prevBtn.textContent = '⏮';
    const nextBtn = document.createElement('button');
    nextBtn.textContent = '⏭';

    const frameSlider = document.createElement('input');
    frameSlider.type = 'range';
    frameSlider.min = '0';
    frameSlider.max = '0';
    frameSlider.value = '0';
    frameSlider.style.flex = '1';

    const frameLabel = document.createElement('span');
    frameLabel.style.fontSize = '12px';
    frameLabel.style.fontFamily = 'monospace';
    frameLabel.style.color = '#475569';
    frameLabel.textContent = 'Frame 0/0';

    ctrlRow.appendChild(playBtn);
    ctrlRow.appendChild(prevBtn);
    ctrlRow.appendChild(nextBtn);
    ctrlRow.appendChild(frameSlider);
    ctrlRow.appendChild(frameLabel);
    rightPanel.appendChild(ctrlRow);

    const energyCanvas = document.createElement('canvas');
    energyCanvas.width = 680;
    energyCanvas.height = 140;
    energyCanvas.style.width = '100%';
    energyCanvas.style.height = '140px';
    energyCanvas.style.background = '#ffffff';
    energyCanvas.style.border = '1px solid #cbd5e1';
    energyCanvas.style.borderRadius = '8px';
    rightPanel.appendChild(energyCanvas);

    root.appendChild(rightPanel);
    container.appendChild(root);

    return {
        root, seqList, juncList, bSlider, bVal,
        restitchBtn, viewerBox, hudOverlay,
        btnModeCascade, btnModeOverlap, overlapSubControls,
        juncSelect, alignTypeBtn, rmsdBadge,
        playBtn, prevBtn, nextBtn, frameSlider, frameLabel,
        energyCanvas, ctrlRow,
        mFrames: mFrames.val,
        mSteps: mSteps.val,
        mBar: mBar.val,
        mRxn: mRxn.val
    };
}
