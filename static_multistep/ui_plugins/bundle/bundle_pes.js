/**
 * Interactive 2D Canvas Potential Energy Surface (PES) chart for Step 7 Bundle Viewer.
 * Uses the proven renderEnergyProfile engine from Step 5 for full visual consistency and clear contrast.
 */
import { renderEnergyProfile } from '../irc_profile_plot.js';

export class BundlePesChart {
    constructor(canvasEl, optionsOrSeek) {
        this.canvas = canvasEl;
        if (typeof optionsOrSeek === 'function') {
            this.onSeek = optionsOrSeek;
            this.trajectory = null;
            this.energies = [];
            this.boundaries = [];
        } else {
            const opts = optionsOrSeek || {};
            this.onSeek = opts.onSeekFrame || opts.onSeek || null;
            this.trajectory = opts.trajectory || null;
            this.energies = opts.energies || [];
            this.boundaries = opts.boundaries || [];
        }
        this.currentFrame = 0;
        this.draw();
    }

    setData(trajectory, boundaries) {
        if (trajectory && trajectory.frames) {
            this.trajectory = trajectory;
        } else if (Array.isArray(trajectory)) {
            this.energies = trajectory;
        }
        if (boundaries) this.boundaries = boundaries;
        this.draw();
    }

    setCurrentFrame(frameIdx) {
        this.currentFrame = frameIdx;
        this.draw();
    }

    draw() {
        if (!this.canvas) return;
        const traj = this.trajectory || (this.energies && this.energies.length ? {
            frames: this.energies.map((e, idx) => ({ s_arc: idx, rel_energy_kcal: e })),
            step_boundaries: this.boundaries || []
        } : null);

        if (!traj) return;

        renderEnergyProfile(this.canvas, traj, this.currentFrame, (clickedIdx) => {
            this.currentFrame = clickedIdx;
            if (this.onSeek) this.onSeek(clickedIdx);
        });
    }
}
