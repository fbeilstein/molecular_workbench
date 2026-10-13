export function buildIRCLayout(container) {
    const root = document.createElement('div');
    root.style.display = 'flex';
    root.style.gap = '16px';
    root.style.flexWrap = 'wrap';
    root.style.width = '100%';

    // --- Left: 3D Viewport ---
    const viewerBox = document.createElement('div');
    viewerBox.style.flex = '1 1 540px';
    viewerBox.style.height = '480px';
    viewerBox.style.background = '#0f0f1a';
    viewerBox.style.border = '1px solid #2a2a45';
    viewerBox.style.borderRadius = '10px';
    viewerBox.style.position = 'relative';
    viewerBox.style.overflow = 'hidden';

    const hudOverlay = document.createElement('div');
    hudOverlay.style.position = 'absolute';
    hudOverlay.style.top = '10px';
    hudOverlay.style.left = '12px';
    hudOverlay.style.padding = '6px 12px';
    hudOverlay.style.background = 'rgba(22, 22, 37, 0.85)';
    hudOverlay.style.borderRadius = '6px';
    hudOverlay.style.border = '1px solid #2a2a45';
    hudOverlay.style.fontSize = '12px';
    hudOverlay.style.fontFamily = 'monospace';
    hudOverlay.style.color = '#e0e0f0';
    hudOverlay.style.zIndex = '5';
    hudOverlay.textContent = 'Ready for IRC trajectory.';
    viewerBox.appendChild(hudOverlay);
    root.appendChild(viewerBox);

    // --- Right: Side Panel ---
    const sidePanel = document.createElement('div');
    sidePanel.style.flex = '1 1 380px';
    sidePanel.style.display = 'flex';
    sidePanel.style.flexDirection = 'column';
    sidePanel.style.gap = '12px';

    // 1. Header & Metrics Card
    const metricsCard = document.createElement('div');
    metricsCard.style.padding = '12px 16px';
    metricsCard.style.background = '#ffffff';
    metricsCard.style.border = '1px solid #e2e8f0';
    metricsCard.style.borderRadius = '8px';

    const title = document.createElement('h3');
    title.style.margin = '0 0 6px 0';
    title.style.fontSize = '15px';
    title.style.color = '#0f172a';
    title.textContent = 'Reaction Path (IRC) Trajectory';
    metricsCard.appendChild(title);

    const badgesRow = document.createElement('div');
    badgesRow.style.display = 'flex';
    badgesRow.style.gap = '8px';
    badgesRow.style.flexWrap = 'wrap';
    badgesRow.style.fontSize = '12px';

    const makeBadge = (bg, col) => {
        const span = document.createElement('span');
        span.style.padding = '3px 8px';
        span.style.borderRadius = '4px';
        span.style.background = bg;
        span.style.color = col;
        span.style.fontWeight = '600';
        return span;
    };

    const fwdBadge = makeBadge('#eff6ff', '#2563eb');
    fwdBadge.textContent = 'ΔE‡ fwd: --';
    const revBadge = makeBadge('#fef2f2', '#dc2626');
    revBadge.textContent = 'ΔE‡ rev: --';
    const rxnBadge = makeBadge('#ecfdf5', '#059669');
    rxnBadge.textContent = 'ΔE rxn: --';

    badgesRow.appendChild(fwdBadge);
    badgesRow.appendChild(revBadge);
    badgesRow.appendChild(rxnBadge);
    metricsCard.appendChild(badgesRow);
    sidePanel.appendChild(metricsCard);

    // 2. Profile Canvas
    const chartBox = document.createElement('div');
    chartBox.style.padding = '8px 12px';
    chartBox.style.background = '#ffffff';
    chartBox.style.border = '1px solid #e2e8f0';
    chartBox.style.borderRadius = '8px';

    const chartTitle = document.createElement('div');
    chartTitle.style.fontSize = '12px';
    chartTitle.style.fontWeight = '600';
    chartTitle.style.color = '#334155';
    chartTitle.style.marginBottom = '4px';
    chartTitle.textContent = 'Energy Profile ΔE (kcal/mol) vs Reaction Coordinate s:';
    chartBox.appendChild(chartTitle);

    const energyCanvas = document.createElement('canvas');
    energyCanvas.width = 380;
    energyCanvas.height = 140;
    energyCanvas.style.width = '100%';
    energyCanvas.style.borderRadius = '6px';
    energyCanvas.style.cursor = 'pointer';
    chartBox.appendChild(energyCanvas);
    sidePanel.appendChild(chartBox);

    // 3. Player Card
    const playerCard = document.createElement('div');
    playerCard.style.padding = '12px 16px';
    playerCard.style.background = '#ffffff';
    playerCard.style.border = '1px solid #e2e8f0';
    playerCard.style.borderRadius = '8px';
    playerCard.style.display = 'flex';
    playerCard.style.flexDirection = 'column';
    playerCard.style.gap = '8px';

    const sliderRow = document.createElement('div');
    sliderRow.style.display = 'flex';
    sliderRow.style.alignItems = 'center';
    sliderRow.style.gap = '10px';

    const frameSlider = document.createElement('input');
    frameSlider.type = 'range';
    frameSlider.min = '0';
    frameSlider.max = '0';
    frameSlider.value = '0';
    frameSlider.style.flex = '1';
    frameSlider.disabled = true;

    const frameCountLabel = document.createElement('span');
    frameCountLabel.style.fontSize = '12px';
    frameCountLabel.style.fontFamily = 'monospace';
    frameCountLabel.style.color = '#475569';
    frameCountLabel.textContent = 'Frame 0/0';

    sliderRow.appendChild(frameSlider);
    sliderRow.appendChild(frameCountLabel);
    playerCard.appendChild(sliderRow);

    const btnRow = document.createElement('div');
    btnRow.style.display = 'flex';
    btnRow.style.gap = '8px';
    btnRow.style.alignItems = 'center';

    const playBtn = document.createElement('button');
    playBtn.className = 'primary';
    playBtn.style.padding = '6px 14px';
    playBtn.style.fontSize = '13px';
    playBtn.disabled = true;
    playBtn.textContent = '▶ Play';

    const prevBtn = document.createElement('button');
    prevBtn.style.padding = '6px 10px';
    prevBtn.disabled = true;
    prevBtn.textContent = '⏮';

    const nextBtn = document.createElement('button');
    nextBtn.style.padding = '6px 10px';
    nextBtn.disabled = true;
    nextBtn.textContent = '⏭';

    const loopCheckLabel = document.createElement('label');
    loopCheckLabel.style.fontSize = '12px';
    loopCheckLabel.style.color = '#475569';
    loopCheckLabel.style.display = 'flex';
    loopCheckLabel.style.alignItems = 'center';
    loopCheckLabel.style.gap = '4px';
    const loopCheckbox = document.createElement('input');
    loopCheckbox.type = 'checkbox';
    loopCheckbox.checked = true;
    loopCheckLabel.appendChild(loopCheckbox);
    loopCheckLabel.appendChild(document.createTextNode('Loop'));

    btnRow.appendChild(playBtn);
    btnRow.appendChild(prevBtn);
    btnRow.appendChild(nextBtn);
    btnRow.appendChild(loopCheckLabel);
    playerCard.appendChild(btnRow);

    // Solvent & Endpoint Controls
    const optRow = document.createElement('div');
    optRow.style.display = 'flex';
    optRow.style.gap = '8px';
    optRow.style.alignItems = 'center';
    optRow.style.flexWrap = 'wrap';
    optRow.style.padding = '4px 0';
    optRow.style.fontSize = '12px';

    const solvLabel = document.createElement('span');
    solvLabel.style.color = '#475569';
    solvLabel.style.fontWeight = '500';
    solvLabel.textContent = 'Solvent:';

    const solventSelect = document.createElement('select');
    solventSelect.style.padding = '3px 6px';
    solventSelect.style.fontSize = '12px';
    solventSelect.style.borderRadius = '4px';
    solventSelect.style.border = '1px solid #cbd5e1';
    solventSelect.style.background = '#f8fafc';

    const solvOptions = [
        { val: 'water', label: 'Water (ALPB, ε=80.4)' },
        { val: 'gas', label: 'Gas Phase (Vacuum)' },
        { val: 'acetonitrile', label: 'Acetonitrile (ε=37.5)' },
        { val: 'ch2cl2', label: 'DCM (ε=8.9)' }
    ];
    for (const opt of solvOptions) {
        const el = document.createElement('option');
        el.value = opt.val;
        el.textContent = opt.label;
        solventSelect.appendChild(el);
    }

    const relaxLabel = document.createElement('label');
    relaxLabel.style.display = 'flex';
    relaxLabel.style.alignItems = 'center';
    relaxLabel.style.gap = '4px';
    relaxLabel.style.color = '#475569';
    relaxLabel.style.cursor = 'pointer';

    const relaxCheckbox = document.createElement('input');
    relaxCheckbox.type = 'checkbox';
    relaxCheckbox.checked = true;
    relaxLabel.appendChild(relaxCheckbox);
    relaxLabel.appendChild(document.createTextNode('Relax to isolated products'));

    optRow.appendChild(solvLabel);
    optRow.appendChild(solventSelect);
    optRow.appendChild(relaxLabel);
    playerCard.appendChild(optRow);

    const actionRow = document.createElement('div');
    actionRow.style.display = 'flex';
    actionRow.style.gap = '8px';
    actionRow.style.marginTop = '4px';

    const calcBtn = document.createElement('button');
    calcBtn.className = 'primary';
    calcBtn.style.flex = '1';
    calcBtn.style.padding = '8px 12px';
    calcBtn.style.fontSize = '13px';
    calcBtn.style.fontWeight = '600';
    calcBtn.textContent = '⚡ Calculate IRC Path (xTB)';

    const downloadBtn = document.createElement('button');
    downloadBtn.style.padding = '8px 12px';
    downloadBtn.style.fontSize = '13px';
    downloadBtn.disabled = true;
    downloadBtn.textContent = '💾 Export JSON';

    actionRow.appendChild(calcBtn);
    actionRow.appendChild(downloadBtn);
    playerCard.appendChild(actionRow);
    sidePanel.appendChild(playerCard);

    // 4. Log Box
    const logBox = document.createElement('textarea');
    logBox.style.width = '100%';
    logBox.style.height = '90px';
    logBox.style.padding = '8px';
    logBox.style.background = '#f1f5f9';
    logBox.style.border = '1px solid #cbd5e1';
    logBox.style.borderRadius = '6px';
    logBox.style.fontFamily = 'monospace';
    logBox.style.fontSize = '11.5px';
    logBox.style.color = '#334155';
    logBox.readOnly = true;
    logBox.value = 'Ready to compute or view IRC trajectory.';
    sidePanel.appendChild(logBox);

    root.appendChild(sidePanel);
    container.appendChild(root);

    return {
        viewerBox, hudOverlay,
        fwdBadge, revBadge, rxnBadge,
        energyCanvas, frameSlider, frameCountLabel,
        playBtn, prevBtn, nextBtn, loopCheckbox,
        solventSelect, relaxCheckbox,
        calcBtn, downloadBtn, logBox
    };
}
