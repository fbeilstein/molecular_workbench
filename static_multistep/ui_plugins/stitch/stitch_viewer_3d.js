/**
 * 3Dmol Boundary Alignment & Transition Viewer for Step 2 Stitch Align.
 */

export class Stitch3DViewer {
    constructor() {
        this.viewer = null;
        this.viewMode = 'transition';
        this.transitionFrameIdx = 0;
        this.isTransitionPlaying = false;
        this.transitionTimer = null;
        this.blinkTimer = null;
        this.blinkState = 0;
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
        let xyz = `${n}\nFrame\n`;
        for (let i = 0; i < n; i++) {
            const sym = symbols[i];
            const [x, y, z] = coords[i];
            xyz += `${sym} ${x.toFixed(6)} ${y.toFixed(6)} ${z.toFixed(6)}\n`;
        }
        return xyz;
    }

    getTransitionSequence(j) {
        if (!j) return [];
        const seq = [
            {
                coords: j.prev_coords,
                symbols: j.prev_symbols,
                type: 'start',
                label: `Step ${j.step_prev_index} Product (Arrival Frame)`,
                color: '#38bdf8'
            }
        ];

        const bridges = j.bridge_frames || [];
        bridges.forEach((b, i) => {
            seq.push({
                coords: b.positions || b.coords,
                symbols: j.prev_symbols,
                type: 'bridge',
                label: `Bridge Frame ${i + 1} / ${bridges.length}`,
                energy: b.energy_kcal,
                color: '#a855f7'
            });
        });

        seq.push({
            coords: j.next_coords_aligned,
            symbols: j.next_symbols,
            type: 'end',
            label: `Step ${j.step_next_index} Reactant (Aligned Start)`,
            color: '#f43f5e'
        });

        return seq;
    }

    apply3DScene(j, mode, resetCamera = false) {
        if (!this.viewer || !j) return;
        this.viewMode = mode;
        const badge = document.getElementById('junction-status-badge');
        const transControls = document.getElementById('transition-controls-bar');

        if (mode === 'transition') {
            if (transControls) transControls.style.display = 'flex';
            this.renderTransitionFrame(j, this.transitionFrameIdx, resetCamera);
            return;
        }

        if (transControls) transControls.style.display = 'none';
        this.pauseTransition();
        this.viewer.clear();

        const xyzPrev = this._buildXYZ(j.prev_symbols, j.prev_coords);
        const xyzNextAligned = this._buildXYZ(j.next_symbols, j.next_coords_aligned);
        const xyzNextUnaligned = this._buildXYZ(j.next_symbols, j.next_coords_unaligned);

        if (mode === 'overlay') {
            const m1 = this.viewer.addModel(xyzPrev, 'xyz');
            const m2 = this.viewer.addModel(xyzNextAligned, 'xyz');

            this.viewer.setStyle({ model: m1 }, {
                stick: { color: '#38bdf8', radius: 0.12 },
                sphere: { color: '#38bdf8', scale: 0.22 }
            });
            this.viewer.setStyle({ model: m2 }, {
                stick: { color: '#f43f5e', radius: 0.12 },
                sphere: { color: '#f43f5e', scale: 0.22 }
            });

            if (badge) badge.innerText = `Overlay Aligned | RMSD: ${(j.rmsd || 0).toFixed(4)} Å (Cyan = Step K, Magenta = Step K+1)`;

        } else if (mode === 'prev_only') {
            const m1 = this.viewer.addModel(xyzPrev, 'xyz');
            this.viewer.setStyle({ model: m1 }, {
                stick: { color: '#38bdf8', radius: 0.14 },
                sphere: { color: '#38bdf8', scale: 0.26 }
            });
            if (badge) badge.innerText = `Step ${j.step_prev_index} Final Product Frame`;

        } else if (mode === 'next_only') {
            const m2 = this.viewer.addModel(xyzNextAligned, 'xyz');
            this.viewer.setStyle({ model: m2 }, {
                stick: { color: '#f43f5e', radius: 0.14 },
                sphere: { color: '#f43f5e', scale: 0.26 }
            });
            if (badge) badge.innerText = `Step ${j.step_next_index} Initial Reactant Frame (Aligned)`;

        } else if (mode === 'unaligned') {
            const m1 = this.viewer.addModel(xyzPrev, 'xyz');
            const m2 = this.viewer.addModel(xyzNextUnaligned, 'xyz');
            this.viewer.setStyle({ model: m1 }, { stick: { color: '#38bdf8', radius: 0.12 }, sphere: { color: '#38bdf8', scale: 0.22 } });
            this.viewer.setStyle({ model: m2 }, { stick: { color: '#fbbf24', radius: 0.12 }, sphere: { color: '#fbbf24', scale: 0.22 } });
            if (badge) badge.innerText = `Raw Unaligned Coordinates (Yellow = Step K+1 before fit)`;

        } else if (mode === 'blink') {
            const isPrev = this.blinkState === 0;
            const xyzActive = isPrev ? xyzPrev : xyzNextAligned;
            const col = isPrev ? '#38bdf8' : '#f43f5e';
            const m = this.viewer.addModel(xyzActive, 'xyz');
            this.viewer.setStyle({ model: m }, {
                stick: { color: col, radius: 0.14 },
                sphere: { color: col, scale: 0.26 }
            });
            if (badge) badge.innerText = isPrev ? `[Blink] Step ${j.step_prev_index} End` : `[Blink] Step ${j.step_next_index} Start (RMSD: ${(j.rmsd || 0).toFixed(4)} Å)`;
        }

        if (resetCamera) {
            this.viewer.zoomTo();
        }
        this.viewer.render();
    }

    renderTransitionFrame(j, frameIdx, resetCamera = false) {
        if (!this.viewer || !j) return;
        const seq = this.getTransitionSequence(j);
        if (seq.length === 0) return;

        const idx = Math.max(0, Math.min(seq.length - 1, frameIdx));
        this.transitionFrameIdx = idx;
        const item = seq[idx];

        this.viewer.clear();
        const xyz = this._buildXYZ(item.symbols, item.coords);
        const m = this.viewer.addModel(xyz, 'xyz');

        if (item.type === 'start') {
            this.viewer.setStyle({ model: m }, {
                stick: { color: '#38bdf8', radius: 0.14 },
                sphere: { color: '#38bdf8', scale: 0.25 }
            });
        } else if (item.type === 'end') {
            this.viewer.setStyle({ model: m }, {
                stick: { color: '#f43f5e', radius: 0.14 },
                sphere: { color: '#f43f5e', scale: 0.25 }
            });
        } else {
            this.viewer.setStyle({ model: m }, {
                stick: { color: '#a855f7', radius: 0.13 },
                sphere: { color: '#c084fc', scale: 0.23 }
            });
        }

        if (resetCamera) {
            this.viewer.zoomTo();
        }
        this.viewer.render();

        const scrubber = document.getElementById('trans-scrubber');
        const counter = document.getElementById('trans-counter');
        const badge = document.getElementById('junction-status-badge');

        if (scrubber) scrubber.value = idx;
        if (counter) counter.innerText = `${idx + 1} / ${seq.length}`;
        if (badge) badge.innerText = `${item.label} (RMSD: ${(j.rmsd || 0).toFixed(4)} Å)`;
    }

    stepTransition(j, delta) {
        this.pauseTransition();
        const seq = this.getTransitionSequence(j);
        if (seq.length === 0) return;
        let target = this.transitionFrameIdx + delta;
        target = Math.max(0, Math.min(seq.length - 1, target));
        this.renderTransitionFrame(j, target);
    }

    toggleTransitionPlay(j) {
        if (this.isTransitionPlaying) this.pauseTransition();
        else this.playTransition(j);
    }

    playTransition(j) {
        const seq = this.getTransitionSequence(j);
        if (seq.length === 0) return;

        if (this.transitionFrameIdx >= seq.length - 1) {
            this.transitionFrameIdx = 0;
        }

        this.isTransitionPlaying = true;
        const btn = document.getElementById('btn-trans-play');
        if (btn) btn.innerHTML = '⏸';

        clearInterval(this.transitionTimer);
        this.transitionTimer = setInterval(() => {
            const next = this.transitionFrameIdx + 1;
            if (next >= seq.length) {
                this.pauseTransition();
                return;
            }
            this.renderTransitionFrame(j, next);
        }, 110);
    }

    pauseTransition() {
        this.isTransitionPlaying = false;
        clearInterval(this.transitionTimer);
        const btn = document.getElementById('btn-trans-play');
        if (btn) {
            btn.innerHTML = `
                <svg width="14" height="14" viewBox="0 0 24 24" fill="currentColor">
                    <polygon points="5 3 19 12 5 21 5 3"/>
                </svg>
            `;
        }
    }

    startBlink(j) {
        this.stopBlink();
        this.blinkState = 0;
        this.blinkTimer = setInterval(() => {
            this.blinkState = 1 - this.blinkState;
            this.apply3DScene(j, 'blink');
        }, 400);
    }

    stopBlink() {
        clearInterval(this.blinkTimer);
        this.blinkTimer = null;
    }
}
