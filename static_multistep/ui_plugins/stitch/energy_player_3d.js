/**
 * 3D Playback Player for Step 3 Energy Stitch.
 */

export class Energy3DPlayer {
    constructor() {
        this.viewer = null;
        this.megajson = null;
        this.activeFrame = 0;
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

    _buildXYZ(symbols, coords) {
        const n = symbols.length;
        let xyz = `${n}\nStitched Frame\n`;
        for (let i = 0; i < n; i++) {
            const sym = symbols[i];
            const [x, y, z] = coords[i];
            xyz += `${sym} ${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)}\n`;
        }
        return xyz;
    }

    setFrame(index, resetCamera = false, onFrameChange = null) {
        if (!this.megajson || !this.viewer) return;
        const frames = this.megajson.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.activeFrame = index;
        const frame = frames[index];
        const symbols = this.megajson.symbols;
        const coords = frame.coords || frame.positions;

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

        const scrubber = document.getElementById('pes-scrubber');
        const counter = document.getElementById('pes-counter');
        const badge = document.getElementById('energy-status-badge');

        const eKcal = frame.energy_kcal !== undefined ? frame.energy_kcal : 0.0;
        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;

        if (badge) {
            const stepIdx = frame.step_index !== undefined ? frame.step_index : '?';
            const isTS = frame.is_ts ? ' ⭐ [TS]' : '';
            badge.innerText = `Step ${stepIdx} | Frame ${index + 1} / ${frames.length} | Energy: ${eKcal.toFixed(2)} kcal/mol ${isTS}`;
        }

        if (onFrameChange) onFrameChange(index);
    }

    stepFrame(delta, onFrameChange = null) {
        this.pause();
        const frames = this.megajson ? this.megajson.frames : [];
        if (!frames || frames.length === 0) return;
        let target = this.activeFrame + delta;
        target = Math.max(0, Math.min(frames.length - 1, target));
        this.setFrame(target, false, onFrameChange);
    }

    togglePlay(onFrameChange = null) {
        if (this.isPlaying) this.pause();
        else this.play(onFrameChange);
    }

    play(onFrameChange = null) {
        if (!this.megajson) return;
        const frames = this.megajson.frames || [];
        if (frames.length === 0) return;

        if (this.activeFrame >= frames.length - 1) {
            this.setFrame(0, false, onFrameChange);
        }

        this.isPlaying = true;
        const btn = document.getElementById('btn-pes-play');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.playTimer);
        this.playTimer = setInterval(() => {
            const next = this.activeFrame + 1;
            if (next >= frames.length) {
                this.pause();
                return;
            }
            this.setFrame(next, false, onFrameChange);
        }, 100);
    }

    pause() {
        this.isPlaying = false;
        clearInterval(this.playTimer);
        const btn = document.getElementById('btn-pes-play');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    }
}
