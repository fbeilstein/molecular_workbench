/**
 * 2D Canvas Chart for Reaction Free Energy Profile (PES).
 */

export function drawPESProfile(canvas, energyProfile, boundaries, activeFrame) {
    if (!canvas) return;

    const ctx = canvas.getContext('2d');
    const dpr = window.devicePixelRatio || 1;
    const rect = canvas.getBoundingClientRect();
    canvas.width = rect.width * dpr;
    canvas.height = rect.height * dpr;
    ctx.scale(dpr, dpr);

    const w = rect.width;
    const h = rect.height;
    const padLeft = 45;
    const padRight = 15;
    const padTop = 20;
    const padBottom = 30;
    const plotW = w - padLeft - padRight;
    const plotH = h - padTop - padBottom;

    const energies = (energyProfile && energyProfile.frame_energies_kcal) || [];
    if (energies.length === 0) return;

    let minE = Math.min(...energies);
    let maxE = Math.max(...energies);
    if (minE === maxE) { minE -= 5; maxE += 5; }
    const margin = (maxE - minE) * 0.15 || 2;
    minE -= margin;
    maxE += margin;

    const toX = (idx) => padLeft + (idx / (energies.length - 1)) * plotW;
    const toY = (val) => padTop + plotH - ((val - minE) / (maxE - minE)) * plotH;

    ctx.clearRect(0, 0, w, h);

    // Grid lines
    ctx.strokeStyle = 'rgba(255, 255, 255, 0.08)';
    ctx.lineWidth = 1;
    ctx.beginPath();
    for (let i = 0; i <= 4; i++) {
        const yVal = minE + (i / 4) * (maxE - minE);
        const yPos = toY(yVal);
        ctx.moveTo(padLeft, yPos);
        ctx.lineTo(w - padRight, yPos);
    }
    ctx.stroke();

    // Y-axis labels
    ctx.fillStyle = '#94a3b8';
    ctx.font = '10px monospace';
    ctx.textAlign = 'right';
    for (let i = 0; i <= 4; i++) {
        const yVal = minE + (i / 4) * (maxE - minE);
        ctx.fillText(yVal.toFixed(1), padLeft - 6, toY(yVal) + 3);
    }

    // Zero reference line
    if (minE <= 0 && maxE >= 0) {
        const zeroY = toY(0);
        ctx.strokeStyle = 'rgba(255, 255, 255, 0.2)';
        ctx.setLineDash([4, 4]);
        ctx.beginPath();
        ctx.moveTo(padLeft, zeroY);
        ctx.lineTo(w - padRight, zeroY);
        ctx.stroke();
        ctx.setLineDash([]);
    }

    // Step boundary shading
    (boundaries || []).forEach((b, k) => {
        const xStart = toX(b.start_frame);
        const xEnd = toX(b.end_frame);
        ctx.fillStyle = k % 2 === 0 ? 'rgba(59, 130, 246, 0.04)' : 'rgba(168, 85, 247, 0.04)';
        ctx.fillRect(xStart, padTop, xEnd - xStart, plotH);
    });

    // Energy line curve
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    energies.forEach((val, i) => {
        const x = toX(i);
        const y = toY(val);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    });
    ctx.strokeStyle = '#38bdf8';
    ctx.stroke();

    // Gradient fill under curve
    const grad = ctx.createLinearGradient(0, padTop, 0, padTop + plotH);
    grad.addColorStop(0, 'rgba(56, 189, 248, 0.25)');
    grad.addColorStop(1, 'rgba(56, 189, 248, 0.0)');
    ctx.lineTo(toX(energies.length - 1), toY(minE));
    ctx.lineTo(toX(0), toY(minE));
    ctx.closePath();
    ctx.fillStyle = grad;
    ctx.fill();

    // Transition state markers
    (boundaries || []).forEach(b => {
        const tsIdx = b.ts_frame;
        if (tsIdx !== undefined && tsIdx < energies.length) {
            const x = toX(tsIdx);
            const y = toY(energies[tsIdx]);
            ctx.fillStyle = '#f59e0b';
            ctx.beginPath();
            ctx.arc(x, y, 4.5, 0, 2 * Math.PI);
            ctx.fill();
            ctx.fillStyle = '#fde68a';
            ctx.font = 'bold 9px monospace';
            ctx.textAlign = 'center';
            ctx.fillText(`TS${b.step_index}`, x, y - 8);
        }
    });

    // Active frame marker
    if (activeFrame >= 0 && activeFrame < energies.length) {
        const activeX = toX(activeFrame);
        const activeY = toY(energies[activeFrame]);
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 1.5;
        ctx.beginPath();
        ctx.moveTo(activeX, padTop);
        ctx.lineTo(activeX, padTop + plotH);
        ctx.stroke();

        ctx.fillStyle = '#ef4444';
        ctx.beginPath();
        ctx.arc(activeX, activeY, 5, 0, 2 * Math.PI);
        ctx.fill();
    }
}
