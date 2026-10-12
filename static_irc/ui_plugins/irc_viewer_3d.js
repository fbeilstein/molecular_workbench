// 3Dmol.js viewer for Stage 4 IRC Trajectory playback
// Replaces Three.js with 3Dmol matching Molecular Workbench

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

export class IRC3DViewer {
    constructor(container) {
        this.container = container;
        this.viewer = null;
        this.symbols = [];
        this.currentTraj = null;
        this.hasFrames = false;
        this._initViewer();
    }

    _initViewer() {
        if (!window.$3Dmol) {
            console.error('3Dmol library is not loaded on window.$3Dmol');
            return;
        }

        // Molecular Workbench signature dark aesthetic
        this.viewer = window.$3Dmol.createViewer(this.container, {
            backgroundColor: '#0f0f1a'
        });

        window.addEventListener('resize', () => this.resize());
    }

    resize() {
        if (this.viewer) {
            this.viewer.resize();
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

    buildStructure(symbols, positions, bonds = []) {
        if (!this.viewer) return;
        this.symbols = symbols;
        this.viewer.removeAllModels();
        this.viewer.removeAllShapes();

        const activeBonds = filterBondsByDistance(symbols, positions, bonds, 0, 0);
        if (activeBonds.length > 0) {
            const sdf = this._formatSdfRecord(symbols, positions, activeBonds, 'Structure');
            this.viewer.addModel(sdf, 'sdf');
        } else {
            const xyz = this._formatXyz(symbols, positions, 'Structure');
            this.viewer.addModel(xyz, 'xyz');
        }

        // Molecular Workbench ball-and-stick styling with multiple bond support
        this.viewer.setStyle({}, {
            stick: { radius: 0.12, doubleBondScaling: 0.5, tripleBondScaling: 0.3 },
            sphere: { scale: 0.25 }
        });
        this.viewer.zoomTo();
        this.viewer.render();
    }

    loadTrajectory(traj) {
        if (!this.viewer || !traj || !traj.frames || !traj.frames.length) return;
        this.currentTraj = traj;
        this.symbols = traj.symbols;
        const tsIdx = traj.ts_frame_index || 0;
        const hasBonds = traj.bonds && traj.bonds.length > 0;

        try {
            this.viewer.removeAllModels();
            this.viewer.removeAllShapes();

            if (hasBonds) {
                let multiSdf = '';
                traj.frames.forEach((frame, idx) => {
                    const activeBonds = filterBondsByDistance(
                        traj.symbols,
                        frame.positions,
                        traj.bonds,
                        idx,
                        tsIdx
                    );
                    multiSdf += this._formatSdfRecord(
                        traj.symbols,
                        frame.positions,
                        activeBonds,
                        `Frame ${idx} s=${frame.s_arc.toFixed(3)}`
                    );
                });
                this.viewer.addModelsAsFrames(multiSdf, 'sdf');
            } else {
                let multiXyz = '';
                const n = traj.symbols.length;
                traj.frames.forEach((frame, idx) => {
                    multiXyz += `${n}\nFrame ${idx} s=${frame.s_arc.toFixed(3)}\n`;
                    frame.positions.forEach((pos, i) => {
                        const sym = traj.symbols[i] || 'C';
                        multiXyz += `${sym} ${pos[0].toFixed(5)} ${pos[1].toFixed(5)} ${pos[2].toFixed(5)}\n`;
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
            this.setFrame(tsIdx);
        } catch (err) {
            console.warn('Multimodel frame load fallback:', err);
            this.hasFrames = false;
            this.updatePositions(traj.frames[tsIdx].positions, true, tsIdx);
        }
    }

    async setFrame(idx) {
        if (!this.viewer) return;
        if (this.hasFrames) {
            try {
                await this.viewer.setFrame(idx);
                this.viewer.render();
                return;
            } catch (err) {
                this.hasFrames = false;
            }
        }
        // Fallback: update model geometry directly
        if (this.currentTraj && this.currentTraj.frames && this.currentTraj.frames[idx]) {
            this.updatePositions(this.currentTraj.frames[idx].positions, false, idx);
        }
    }

    updatePositions(positions, autoZoom = false, frameIdx = 0) {
        if (!this.viewer || !this.symbols || !positions) return;
        this.viewer.removeAllModels();
        this.viewer.removeAllShapes();

        if (this.currentTraj && this.currentTraj.bonds && this.currentTraj.bonds.length > 0) {
            const tsIdx = this.currentTraj.ts_frame_index || 0;
            const activeBonds = filterBondsByDistance(
                this.symbols,
                positions,
                this.currentTraj.bonds,
                frameIdx,
                tsIdx
            );
            const sdf = this._formatSdfRecord(this.symbols, positions, activeBonds);
            this.viewer.addModel(sdf, 'sdf');
        } else {
            const xyz = this._formatXyz(this.symbols, positions);
            this.viewer.addModel(xyz, 'xyz');
        }

        this.viewer.setStyle({}, {
            stick: { radius: 0.12, doubleBondScaling: 0.5, tripleBondScaling: 0.3 },
            sphere: { scale: 0.25 }
        });
        if (autoZoom) {
            this.viewer.zoomTo();
        }
        this.viewer.render();
    }
}
