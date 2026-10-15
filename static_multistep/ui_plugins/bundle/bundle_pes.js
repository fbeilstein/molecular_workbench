/**
 * Interactive 2D Canvas Potential Energy Surface (PES) chart for Step 5 Bundle Viewer.
 */

export class BundlePesChart {
    constructor(canvasEl, optionsOrSeek) {
        this.canvas = canvasEl;
        this.ctx = canvasEl ? canvasEl.getContext('2d') : null;
        if (typeof optionsOrSeek === 'function') {
            this.onSeek = optionsOrSeek;
            this.energies = [];
            this.landmarks = [];
        } else {
            const opts = optionsOrSeek || {};
            this.onSeek = opts.onSeekFrame || opts.onSeek || null;
            this.energies = opts.energies || [];
            this.landmarks = opts.boundaries || opts.landmarks || [];
        }
        this.currentFrame = 0;
        this.isDragging = false;
        this.initEvents();
        if (this.energies.length > 0) {
            this.draw();
        }
    }

    setData(energies, landmarks) {
        this.energies = energies || [];
        this.landmarks = landmarks || [];
        this.draw();
    }

    setCurrentFrame(frameIdx) {
        this.currentFrame = frameIdx;
        this.draw();
    }

    initEvents() {
        if (!this.canvas) return;
        const handleEvent = (e) => {
            const rect = this.canvas.getBoundingClientRect();
            const x = e.clientX - rect.left;
            const pct = Math.max(0, Math.min(1, x / rect.width));
            const total = this.energies.length || 1;
            const targetFrame = Math.round(pct * (total - 1));
            if (this.onSeek) this.onSeek(targetFrame);
        };

        this.canvas.addEventListener('mousedown', (e) => {
            this.isDragging = true;
            handleEvent(e);
        });
        window.addEventListener('mousemove', (e) => {
            if (this.isDragging) handleEvent(e);
        });
        window.addEventListener('mouseup', () => {
            this.isDragging = false;
        });
    }

    draw() {
        if (!this.ctx || !this.canvas) return;
        const w = this.canvas.width = this.canvas.parentElement ? this.canvas.parentElement.clientWidth : 800;
        const h = this.canvas.height = 140;
        const ctx = this.ctx;
        ctx.clearRect(0, 0, w, h);

        const n = this.energies.length;
        if (n < 2) {
            ctx.fillStyle = '#64748b';
            ctx.font = '12px Inter, sans-serif';
            ctx.textAlign = 'center';
            ctx.fillText('Energy profile data not available', w / 2, h / 2);
            return;
        }

        const minE = Math.min(...this.energies);
        const maxE = Math.max(...this.energies);
        const rangeE = Math.max(0.1, maxE - minE);
        const padX = 24, padY = 24;
        const chartW = w - padX * 2;
        const chartH = h - padY * 2;

        const getX = (i) => padX + (i / (n - 1)) * chartW;
        const getY = (val) => padY + chartH - ((val - minE) / rangeE) * chartH;

        // Draw background gradient area
        ctx.beginPath();
        ctx.moveTo(getX(0), h - padY);
        for (let i = 0; i < n; i++) {
            ctx.lineTo(getX(i), getY(this.energies[i]));
        }
        ctx.lineTo(getX(n - 1), h - padY);
        ctx.closePath();
        const grad = ctx.createLinearGradient(0, padY, 0, h - padY);
        grad.addColorStop(0, 'rgba(56, 189, 248, 0.35)');
        grad.addColorStop(1, 'rgba(56, 189, 248, 0.0)');
        ctx.fillStyle = grad;
        ctx.fill();

        // Draw energy curve line
        ctx.beginPath();
        for (let i = 0; i < n; i++) {
            const x = getX(i), y = getY(this.energies[i]);
            if (i === 0) ctx.moveTo(x, y);
            else ctx.lineTo(x, y);
        }
        ctx.strokeStyle = '#38bdf8';
        ctx.lineWidth = 2.5;
        ctx.stroke();

        // Energy curve drawn cleanly without arbitrary TS/Products landmarks

        // Draw current frame cursor line
        const curX = getX(this.currentFrame);
        const curY = getY(this.energies[this.currentFrame] || minE);
        ctx.beginPath();
        ctx.moveTo(curX, padY / 2);
        ctx.lineTo(curX, h - padY / 2);
        ctx.strokeStyle = '#ef4444';
        ctx.lineWidth = 2;
        ctx.stroke();

        ctx.beginPath();
        ctx.arc(curX, curY, 5, 0, Math.PI * 2);
        ctx.fillStyle = '#ef4444';
        ctx.fill();
        ctx.strokeStyle = '#ffffff';
        ctx.lineWidth = 2;
        ctx.stroke();

        // Energy text at cursor
        ctx.fillStyle = '#ffffff';
        ctx.font = 'bold 11px Inter, sans-serif';
        ctx.textAlign = curX > w - 80 ? 'right' : 'left';
        const eVal = (this.energies[this.currentFrame] || 0).toFixed(1);
        ctx.fillText(`${eVal} kcal/mol`, curX > w - 80 ? curX - 8 : curX + 8, curY);
    }
}
