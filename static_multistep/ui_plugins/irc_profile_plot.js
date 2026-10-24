export function renderEnergyProfile(canvas, trajectory, currentIdx, onPointClick) {
    if (!canvas || !trajectory || !trajectory.frames || trajectory.frames.length === 0) return;
    if (canvas.parentElement && canvas.parentElement.clientWidth > 0) {
        const parentW = canvas.parentElement.clientWidth;
        if (canvas.width !== parentW) {
            canvas.width = parentW;
        }
    }
    const ctx = canvas.getContext('2d');
    const w = canvas.width;
    const h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    const padLeft = 45;
    const padRight = 20;
    const padTop = 20;
    const padBottom = 30;
    const plotW = w - padLeft - padRight;
    const plotH = h - padTop - padBottom;

    const frames = trajectory.frames;
    const s_vals = frames.map((f, i) => f.s_arc !== undefined ? f.s_arc : (f.step !== undefined ? f.step : i));
    const e_vals = frames.map(f => f.rel_energy_kcal !== undefined ? f.rel_energy_kcal : (f.energy_kcal !== undefined ? f.energy_kcal : 0.0));

    const minS = Math.min(...s_vals);
    const maxS = Math.max(...s_vals);
    const minE = Math.min(...e_vals, 0);
    const maxE = Math.max(...e_vals, 0);
    const rangeS = (maxS - minS) || 1.0;
    const rangeE = (maxE - minE) || 1.0;

    const toX = s => padLeft + ((s - minS) / rangeS) * plotW;
    const toY = e => padTop + plotH - ((e - minE) / rangeE) * plotH;

    // Grid lines & Axes
    ctx.strokeStyle = '#e2e8f0';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padLeft, padTop);
    ctx.lineTo(padLeft, padTop + plotH);
    ctx.lineTo(padLeft + plotW, padTop + plotH);
    ctx.stroke();

    // Zero energy line
    const y0 = toY(0);
    ctx.strokeStyle = '#cbd5e1';
    ctx.setLineDash([4, 4]);
    ctx.beginPath();
    ctx.moveTo(padLeft, y0);
    ctx.lineTo(padLeft + plotW, y0);
    ctx.stroke();
    ctx.setLineDash([]);

    // Asymptotic Approach boundary if present
    const asympCount = trajectory.asymp_r_count || 0;
    if (asympCount > 0 && asympCount < frames.length) {
        const asympS = frames[asympCount].s_arc !== undefined ? frames[asympCount].s_arc : asympCount;
        const ax = toX(asympS);
        ctx.strokeStyle = 'rgba(148, 163, 184, 0.4)';
        ctx.setLineDash([2, 2]);
        ctx.beginPath();
        ctx.moveTo(ax, padTop);
        ctx.lineTo(ax, padTop + plotH);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = '#94a3b8';
        ctx.font = '9px sans-serif';
        ctx.fillText('Approach', padLeft + 6, padTop + 12);
        ctx.fillText('Chemical IRC', ax + 6, padTop + 12);
    }

    // Curve
    ctx.strokeStyle = '#4f6ef7';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    frames.forEach((f, idx) => {
        const s = f.s_arc !== undefined ? f.s_arc : (f.step !== undefined ? f.step : idx);
        const e = f.rel_energy_kcal !== undefined ? f.rel_energy_kcal : (f.energy_kcal !== undefined ? f.energy_kcal : 0.0);
        const x = toX(s);
        const y = toY(e);
        if (idx === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    });
    ctx.stroke();

    const tsIndices = new Set(
        trajectory.step_boundaries && trajectory.step_boundaries.length > 0
            ? trajectory.step_boundaries.map(b => b.ts_frame)
            : [trajectory.ts_frame_index]
    );

    // Points on curve
    frames.forEach((f, idx) => {
        const s = f.s_arc !== undefined ? f.s_arc : (f.step !== undefined ? f.step : idx);
        const e = f.rel_energy_kcal !== undefined ? f.rel_energy_kcal : (f.energy_kcal !== undefined ? f.energy_kcal : 0.0);
        const x = toX(s);
        const y = toY(e);
        if (tsIndices.has(idx)) {
            ctx.fillStyle = '#e11d48';
        } else if (idx < asympCount) {
            ctx.fillStyle = '#93c5fd';
        } else {
            ctx.fillStyle = '#3b82f6';
        }
        ctx.beginPath();
        ctx.arc(x, y, tsIndices.has(idx) ? 4.5 : 2.5, 0, Math.PI * 2);
        ctx.fill();
    });

    // Active Frame Highlight
    if (frames[currentIdx]) {
        const curF = frames[currentIdx];
        const curS = curF.s_arc !== undefined ? curF.s_arc : (curF.step !== undefined ? curF.step : currentIdx);
        const curE = curF.rel_energy_kcal !== undefined ? curF.rel_energy_kcal : (curF.energy_kcal !== undefined ? curF.energy_kcal : 0.0);
        const curX = toX(curS);
        const curY = toY(curE);

        // Vertical highlight bar
        ctx.strokeStyle = 'rgba(239, 68, 68, 0.4)';
        ctx.lineWidth = 1.5;
        ctx.setLineDash([3, 3]);
        ctx.beginPath();
        ctx.moveTo(curX, padTop);
        ctx.lineTo(curX, padTop + plotH);
        ctx.stroke();
        ctx.setLineDash([]);

        // Active Ring
        ctx.fillStyle = '#ef4444';
        ctx.beginPath();
        ctx.arc(curX, curY, 6, 0, Math.PI * 2);
        ctx.fill();
        ctx.strokeStyle = '#fff';
        ctx.lineWidth = 2;
        ctx.stroke();

        // Active Energy Label (Dark text with subtle background for pristine contrast)
        const eStr = `${curE.toFixed(1)} kcal/mol`;
        ctx.font = 'bold 11px Inter, monospace, sans-serif';
        const txtW = ctx.measureText(eStr).width;
        const tagX = curX > w - txtW - 20 ? curX - txtW - 10 : curX + 10;
        const tagY = Math.max(padTop + 14, Math.min(curY - 4, padTop + plotH - 6));

        ctx.fillStyle = 'rgba(255, 255, 255, 0.9)';
        ctx.fillRect(tagX - 4, tagY - 11, txtW + 8, 15);
        ctx.strokeStyle = '#cbd5e1';
        ctx.lineWidth = 1;
        ctx.strokeRect(tagX - 4, tagY - 11, txtW + 8, 15);

        ctx.fillStyle = '#1e293b';
        ctx.textAlign = 'left';
        ctx.fillText(eStr, tagX, tagY);
    }

    // Text Labels
    ctx.fillStyle = '#64748b';
    ctx.font = '10.5px monospace';
    ctx.textAlign = 'left';
    ctx.fillText(`${maxE.toFixed(1)} kcal`, 4, padTop + 10);
    ctx.fillText(`${minE.toFixed(1)} kcal`, 4, padTop + plotH);
    ctx.fillText(`s: ${minS.toFixed(1)}`, padLeft, h - 8);

    const endX = padLeft + plotW;
    ctx.fillText(`s: +${maxS.toFixed(1)}`, endX - 45, h - 8);

    if (trajectory.step_boundaries && trajectory.step_boundaries.length > 1) {
        trajectory.step_boundaries.forEach(sb => {
            const tsF = frames[sb.ts_frame];
            if (tsF) {
                const tsS = tsF.s_arc !== undefined ? tsF.s_arc : tsF.step;
                const tsX = toX(tsS);
                ctx.fillText(`TS ${sb.step_index}`, tsX - 14, h - 8);
            }
        });
    } else {
        const tsX = toX(0);
        ctx.fillText(`s: 0 (TS)`, Math.min(tsX - 22, endX - 55), h - 8);
    }

    // Interactive Click
    canvas.onclick = e => {
        const rect = canvas.getBoundingClientRect();
        const clickX = e.clientX - rect.left;
        let closestIdx = 0;
        let minDist = 9999;
        frames.forEach((f, idx) => {
            const s = f.s_arc !== undefined ? f.s_arc : (f.step !== undefined ? f.step : idx);
            const x = toX(s);
            const dist = Math.abs(x - clickX);
            if (dist < minDist) {
                minDist = dist;
                closestIdx = idx;
            }
        });
        if (onPointClick) onPointClick(closestIdx);
    };
}
