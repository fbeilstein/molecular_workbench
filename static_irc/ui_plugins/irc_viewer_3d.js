// 3Dmol.js viewer for Stage 4 IRC Trajectory playback
// Replaces Three.js with 3Dmol matching Molecular Workbench

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

    buildStructure(symbols, positions, bonds = []) {
        if (!this.viewer) return;
        this.symbols = symbols;
        const xyz = this._formatXyz(symbols, positions, 'Transition State');

        this.viewer.removeAllModels();
        this.viewer.removeAllShapes();
        this.viewer.addModel(xyz, 'xyz');
        // Molecular Workbench ball-and-stick styling
        this.viewer.setStyle({}, {
            stick: { radius: 0.12, color: 'spectrum' },
            sphere: { scale: 0.25 }
        });
        this.viewer.zoomTo();
        this.viewer.render();
    }

    loadTrajectory(traj) {
        if (!this.viewer || !traj || !traj.frames || !traj.frames.length) return;
        this.currentTraj = traj;
        this.symbols = traj.symbols;
        const n = traj.symbols.length;

        // Build multi-model XYZ text for native 3Dmol frame animation
        let multiXyz = '';
        traj.frames.forEach((frame, idx) => {
            multiXyz += `${n}\nFrame ${idx} s=${frame.s_arc.toFixed(3)} E=${frame.rel_energy_kcal.toFixed(2)}\n`;
            frame.positions.forEach((pos, i) => {
                const sym = traj.symbols[i] || 'C';
                multiXyz += `${sym} ${pos[0].toFixed(5)} ${pos[1].toFixed(5)} ${pos[2].toFixed(5)}\n`;
            });
        });

        try {
            this.viewer.removeAllModels();
            this.viewer.removeAllShapes();
            this.viewer.addModelsAsFrames(multiXyz, 'xyz');
            this.viewer.setStyle({}, {
                stick: { radius: 0.12 },
                sphere: { scale: 0.25 }
            });
            this.viewer.zoomTo();
            this.hasFrames = true;
            this.setFrame(traj.ts_frame_index || 0);
        } catch (err) {
            console.warn('Multimodel frame load fallback:', err);
            this.hasFrames = false;
            const tsIdx = traj.ts_frame_index || 0;
            this.updatePositions(traj.frames[tsIdx].positions, true);
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
            this.updatePositions(this.currentTraj.frames[idx].positions, false);
        }
    }

    updatePositions(positions, autoZoom = false) {
        if (!this.viewer || !this.symbols || !positions) return;
        const xyz = this._formatXyz(this.symbols, positions);

        this.viewer.removeAllModels();
        this.viewer.removeAllShapes();
        this.viewer.addModel(xyz, 'xyz');
        this.viewer.setStyle({}, {
            stick: { radius: 0.12 },
            sphere: { scale: 0.25 }
        });
        if (autoZoom) {
            this.viewer.zoomTo();
        }
        this.viewer.render();
    }
}
