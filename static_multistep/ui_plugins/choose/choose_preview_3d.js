/**
 * 3D Trajectory Preview & Playback for Step 1 Choose JSONs.
 */

export class Choose3DPreviewer {
    constructor() {
        this.viewer = null;
        this.activeTraj = null;
        this.currentFrame = 0;
        this.isPlaying = false;
        this.playTimer = null;
    }

    init(containerEl) {
        if (!containerEl) return;
        if (!this.viewer) {
            containerEl.innerHTML = '';
            this.viewer = $3Dmol.createViewer(containerEl, {
                defaultcolors: $3Dmol.rasmolElementColors,
                backgroundColor: '0x0d121c'
            });
        }
    }

    setTrajectory(traj) {
        this.pause();
        this.activeTraj = traj;
        this.currentFrame = 0;
        this.setFrame(0, true);
    }

    _buildXYZ(symbols, coords) {
        const n = symbols.length;
        let xyz = `${n}\nFrame\n`;
        for (let i = 0; i < n; i++) {
            const sym = symbols[i];
            const [x, y, z] = coords[i];
            xyz += `${sym} ${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)}\n`;
        }
        return xyz;
    }

    setFrame(index, resetCamera = false) {
        if (!this.activeTraj || !this.viewer) return;
        const frames = this.activeTraj.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.currentFrame = index;
        const frame = frames[index];
        const symbols = this.activeTraj.symbols;
        const coords = frame.coords;

        const xyz = this._buildXYZ(symbols, coords);
        this.viewer.clear();
        this.viewer.addModel(xyz, 'xyz');
        this.viewer.setStyle({}, {
            stick: { radius: 0.14 },
            sphere: { scale: 0.25 }
        });

        if (resetCamera) {
            this.viewer.zoomTo();
        }
        this.viewer.render();

        const scrubber = document.getElementById('frame-scrubber');
        const counter = document.getElementById('frame-label');
        const badge = document.getElementById('preview-status-badge');

        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;

        if (badge) {
            const isTS = index === this.activeTraj.ts_frame_index;
            const energy = frame.energy_kcal !== undefined ? frame.energy_kcal.toFixed(2) : '0.00';
            badge.innerText = `Frame ${index + 1} / ${frames.length} | Energy: ${energy} kcal/mol ${isTS ? ' ⭐ [TRANSITION STATE]' : ''}`;
        }
    }

    stepFrame(delta) {
        this.pause();
        const frames = this.activeTraj ? this.activeTraj.frames : [];
        if (!frames || frames.length === 0) return;
        let target = this.currentFrame + delta;
        target = Math.max(0, Math.min(frames.length - 1, target));
        this.setFrame(target);
    }

    togglePlay() {
        if (this.isPlaying) {
            this.pause();
        } else {
            this.play();
        }
    }

    play() {
        if (!this.activeTraj) return;
        const frames = this.activeTraj.frames || [];
        if (frames.length === 0) return;

        if (this.currentFrame >= frames.length - 1) {
            this.setFrame(0);
        }

        this.isPlaying = true;
        const btn = document.getElementById('btn-play-pause');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.playTimer);
        this.playTimer = setInterval(() => {
            const next = this.currentFrame + 1;
            if (next >= frames.length) {
                this.pause();
                return;
            }
            this.setFrame(next);
        }, 100);
    }

    pause() {
        this.isPlaying = false;
        clearInterval(this.playTimer);
        const btn = document.getElementById('btn-play-pause');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    }
}
