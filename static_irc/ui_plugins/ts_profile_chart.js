/**
 * 2D Canvas chart renderer for reaction path energy profiles.
 */

export function drawProfile(canvas, profile) {
    if (!profile || profile.length < 2) return;
    const ctx = canvas.getContext('2d');
    const w = canvas.width, h = canvas.height;
    ctx.clearRect(0, 0, w, h);

    const padL = 40, padR = 20, padT = 20, padB = 30;
    const pw = w - padL - padR, ph = h - padT - padB;

    const xs = profile.map(p => p[0]);
    const ys = profile.map(p => p[1]);
    const minX = Math.min(...xs), maxX = Math.max(...xs);
    const minY = Math.min(...ys), maxY = Math.max(...ys);
    const rangeX = maxX - minX || 1.0;
    const rangeY = maxY - minY || 1.0;

    const toX = x => padL + ((x - minX) / rangeX) * pw;
    const toY = y => padT + ph - ((y - minY) / rangeY) * ph;

    // Draw grid and zero line
    ctx.strokeStyle = '#e1e4ea';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, padT + ph);
    ctx.lineTo(padL + pw, padT + ph);
    ctx.moveTo(padL, padT);
    ctx.lineTo(padL, padT + ph);
    ctx.stroke();

    // Energy curve
    ctx.strokeStyle = '#4f6ef7';
    ctx.lineWidth = 2.5;
    ctx.beginPath();
    profile.forEach((pt, i) => {
        const x = toX(pt[0]), y = toY(pt[1]);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    });
    ctx.stroke();

    // Data points & barrier peak marker
    let maxIdx = 0;
    ys.forEach((y, i) => { if (y > ys[maxIdx]) maxIdx = i; });

    profile.forEach((pt, i) => {
        const x = toX(pt[0]), y = toY(pt[1]);
        ctx.fillStyle = (i === maxIdx) ? '#e0245e' : '#4f6ef7';
        ctx.beginPath();
        ctx.arc(x, y, (i === maxIdx) ? 5 : 3.5, 0, 2 * Math.PI);
        ctx.fill();
    });

    // Labels
    ctx.fillStyle = '#666';
    ctx.font = '10px sans-serif';
    ctx.fillText(`${minY.toFixed(1)} kcal/mol`, 5, toY(minY));
    ctx.fillText(`${maxY.toFixed(1)} kcal/mol`, 5, toY(maxY) + 10);
    ctx.fillText('Reaction coordinate (λ)', padL + pw / 2 - 50, h - 10);
}
