/**
 * 3D Trajectory Preview & Playback for Step 1 Choose JSONs.
 * Uses 3Dmol SDF models with explicit bond order and dynamic distance filtering,
 * matching Molecular Workbench server_irc aesthetics.
 */

const COVALENT_RADII = {
    'H': 0.31, 'He': 0.28, 'Li': 1.28, 'Be': 0.96, 'B': 0.84, 'C': 0.76, 'N': 0.71, 'O': 0.66,
    'F': 0.57, 'Ne': 0.58, 'Na': 1.66, 'Mg': 1.41, 'Al': 1.21, 'Si': 1.11, 'P': 1.07, 'S': 1.05,
    'Cl': 1.02, 'Ar': 1.06, 'K': 2.03, 'Ca': 1.76, 'Br': 1.20, 'I': 1.39
};

function filterBondsByDistance(symbols, positions, bonds, frameIdx = 0, tsIdx = 0) {
    if (!bonds || !bonds.length || !positions) return [];
    const active = [];
    bonds.forEach(b => {
        const u = b.from !== undefined ? b.from : b[0];
        const v = b.to !== undefined ? b.to : b[1];
        const p1 = positions[u], p2 = positions[v];
        if (!p1 || !p2) return;
        const dx = p1[0] - p2[0], dy = p1[1] - p2[1], dz = p1[2] - p2[2];
        const dist = Math.sqrt(dx * dx + dy * dy + dz * dz);
        const r1 = COVALENT_RADII[symbols[u]] || 0.77;
        const r2 = COVALENT_RADII[symbols[v]] || 0.77;
        const maxDist = r1 + r2 + 0.45;
        if (dist <= maxDist) {
            let ord = 1;
            if (b.r_order !== undefined && b.p_order !== undefined) {
                ord = frameIdx <= tsIdx ? b.r_order : b.p_order;
            } else if (b.order !== undefined) {
                ord = b.order;
            }
            active.push({ from: u, to: v, order: Math.max(1, Math.min(3, Math.round(ord))) });
        }
    });
    return active;
}

export class Choose3DPreviewer {
    constructor() {
        this.viewer = null;
        this.activeTraj = null;
        this.currentFrame = 0;
        this.isPlaying = false;
        this.playTimer = null;
        this.hasFrames = false;
    }

    init(containerEl) {
        if (!containerEl) return;
        if (!this.viewer) {
            containerEl.innerHTML = '';
            this.viewer = window.$3Dmol.createViewer(containerEl, {
                backgroundColor: '#0f0f1a'
            });
        }
    }

    _formatXyz(symbols, positions, comment = 'Frame') {
        let out = `${symbols.length}\n${comment}\n`;
        positions.forEach((pos, i) => {
            const sym = symbols[i] || 'C';
            out += `${sym} ${pos[0].toFixed(5)} ${pos[1].toFixed(5)} ${pos[2].toFixed(5)}\n`;
        });
        return out;
    }

    _formatSdfRecord(symbols, positions, activeBonds = [], title = 'Frame') {
        const na = String(symbols.length).padStart(3, ' ');
        const nb = String(activeBonds.length).padStart(3, ' ');
        let out = `${title}\n  3Dmol\n\n${na}${nb}  0  0  0  0  0  0  0  0999 V2000\n`;

        positions.forEach((pos, i) => {
            const x = pos[0].toFixed(4).padStart(10, ' ');
            const y = pos[1].toFixed(4).padStart(10, ' ');
            const z = pos[2].toFixed(4).padStart(10, ' ');
            const sym = (symbols[i] || 'C').padEnd(3, ' ');
            out += `${x}${y}${z} ${sym} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
        });

        activeBonds.forEach(b => {
            const u = String(b.from + 1).padStart(3, ' ');
            const v = String(b.to + 1).padStart(3, ' ');
            const t = String(b.order || 1).padStart(3, ' ');
            out += `${u}${v}${t}  0  0  0  0\n`;
        });
        out += `M  END\n$$$$\n`;
        return out;
    }

    setTrajectory(traj) {
        this.pause();
        this.activeTraj = traj;
        this.currentFrame = 0;
        if (!this.viewer || !traj || !traj.frames || !traj.frames.length) return;

        const symbols = traj.symbols || [];
        const hasBonds = traj.bonds && traj.bonds.length > 0;
        const tsIdx = traj.ts_frame_index || 0;

        try {
            this.viewer.removeAllModels();
            this.viewer.removeAllShapes();

            if (hasBonds) {
                let multiSdf = '';
                traj.frames.forEach((frame, idx) => {
                    const pos = frame.positions || frame.coords;
                    const activeBonds = filterBondsByDistance(symbols, pos, traj.bonds, idx, tsIdx);
                    multiSdf += this._formatSdfRecord(symbols, pos, activeBonds, `Frame ${idx}`);
                });
                this.viewer.addModelsAsFrames(multiSdf, 'sdf');
            } else {
                let multiXyz = '';
                const n = symbols.length;
                traj.frames.forEach((frame, idx) => {
                    multiXyz += `${n}\nFrame ${idx}\n`;
                    const pos = frame.positions || frame.coords;
                    pos.forEach((p, i) => {
                        const sym = symbols[i] || 'C';
                        multiXyz += `${sym} ${p[0].toFixed(5)} ${p[1].toFixed(5)} ${p[2].toFixed(5)}\n`;
                    });
                });
                this.viewer.addModelsAsFrames(multiXyz, 'xyz');
            }

            this.viewer.setStyle({}, {
                stick: { radius: 0.12, doubleBondScaling: 0.5, tripleBondScaling: 0.3 },
                sphere: { scale: 0.25 }
            });
            this.viewer.zoomTo();
            this.hasFrames = true;
            this.setFrame(0, true);
        } catch (err) {
            console.warn('Multimodel load fallback:', err);
            this.hasFrames = false;
            this.setFrame(0, true);
        }
    }

    async setFrame(index, resetCamera = false) {
        if (!this.activeTraj || !this.viewer) return;
        const frames = this.activeTraj.frames || [];
        if (index < 0 || index >= frames.length) return;

        this.currentFrame = index;
        const frame = frames[index];
        const symbols = this.activeTraj.symbols;
        const pos = frame.positions || frame.coords;

        if (this.hasFrames) {
            try {
                await this.viewer.setFrame(index);
                if (resetCamera) this.viewer.zoomTo();
                this.viewer.render();
            } catch (err) {
                this.hasFrames = false;
            }
        }

        if (!this.hasFrames) {
            this.viewer.removeAllModels();
            this.viewer.removeAllShapes();
            const tsIdx = this.activeTraj.ts_frame_index || 0;
            if (this.activeTraj.bonds && this.activeTraj.bonds.length > 0) {
                const activeBonds = filterBondsByDistance(symbols, pos, this.activeTraj.bonds, index, tsIdx);
                const sdf = this._formatSdfRecord(symbols, pos, activeBonds, `Frame ${index}`);
                this.viewer.addModel(sdf, 'sdf');
            } else {
                const xyz = this._formatXyz(symbols, pos, `Frame ${index}`);
                this.viewer.addModel(xyz, 'xyz');
            }
            this.viewer.setStyle({}, {
                stick: { radius: 0.12, doubleBondScaling: 0.5, tripleBondScaling: 0.3 },
                sphere: { scale: 0.25 }
            });
            if (resetCamera) this.viewer.zoomTo();
            this.viewer.render();
        }

        const scrubber = document.getElementById('frame-scrubber');
        const counter = document.getElementById('frame-label');
        const badge = document.getElementById('preview-status-badge');

        if (scrubber) scrubber.value = index;
        if (counter) counter.innerText = `${index + 1} / ${frames.length}`;

        if (badge) {
            const isTS = index === this.activeTraj.ts_frame_index;
            const energy = frame.energy_kcal !== undefined ? frame.energy_kcal.toFixed(2) : (frame.rel_energy_kcal !== undefined ? frame.rel_energy_kcal.toFixed(2) : '0.00');
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
