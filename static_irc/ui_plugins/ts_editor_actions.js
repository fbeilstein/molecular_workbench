import { updateMolBlockPositions } from './ts_sdf_parser.js';
import { drawProfile } from './ts_profile_chart.js';

export class TSEditorActions {
    constructor(state, optCard, refineCard, viewer) {
        this.state = state;
        this.optCard = optCard;
        this.refineCard = refineCard;
        this.viewer = viewer;
        this.historyStack = [];
        this.dockVariant = 0;
    }

    pushHistory() {
        const snapshot = this.state.molblocks.map((mb, fi) => {
            const g = this.viewer.moleculeGroups[fi];
            return g ? updateMolBlockPositions(mb, g) : mb;
        });
        this.historyStack.push({
            molblocks: snapshot,
            rxnBonds: JSON.parse(JSON.stringify(this.state.rxnBonds)),
            overridden: this.state.overridden
        });
        if (this.refineCard.undoBtn) this.refineCard.undoBtn.disabled = false;
    }

    popHistory(onRestored) {
        if (!this.historyStack.length) return;
        const prev = this.historyStack.pop();
        this.state.molblocks = prev.molblocks;
        this.state.rxnBonds = prev.rxnBonds;
        this.state.overridden = prev.overridden;
        this.viewer.buildMolecules(this.state.molblocks);
        if (this.refineCard.undoBtn) this.refineCard.undoBtn.disabled = this.historyStack.length === 0;
        this.setStatus('Reverted to previous state.', '#4f6ef7');
        if (onRestored) onRestored();
    }

    setStatus(text, color = '#4f6ef7') {
        if (this.optCard.optStatus) {
            this.optCard.optStatus.style.color = color;
            this.optCard.optStatus.textContent = text;
        }
    }

    async fetchLog() {
        try {
            const r = await fetch('/api/log');
            const t = await r.text();
            if (t && this.optCard.logBox) {
                this.optCard.logBox.value = t;
                this.optCard.logBox.scrollTop = this.optCard.logBox.scrollHeight;
            }
        } catch (e) { }
    }

    async refineFragment(fi, method) {
        this.pushHistory();
        this.setStatus(`Refining Reactant ${fi + 1} with ${method.toUpperCase()}...`);
        const g = this.viewer.moleculeGroups[fi];
        const currentBlock = updateMolBlockPositions(this.state.molblocks[fi], g);
        try {
            const payload = { action_type: 'refine_fragment', fragment_index: fi, method, molblock: currentBlock };
            const res = await fetch('/api/step/step_3_ts/custom', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });
            const r = await res.json();
            if (r.success) {
                this.state.molblocks[fi] = r.molblock;
                this.viewer.buildMolecules(this.state.molblocks);
                this.setStatus(r.message, '#18a558');
            } else {
                this.setStatus(r.error, '#e0245e');
            }
        } catch (err) {
            this.setStatus('Refine error: ' + err.message, '#e0245e');
        }
    }

    async resetToClean(cleanReactants) {
        this.pushHistory();
        if (cleanReactants && cleanReactants.length) {
            this.state.molblocks = [...cleanReactants];
            this.viewer.buildMolecules(this.state.molblocks);
            this.setStatus('Reset scene to clean reactant conformers.', '#4f6ef7');
            return;
        }
        try {
            const res = await fetch('/api/step/step_3_ts/custom', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ action_type: 'reset_clean' })
            });
            const r = await res.json();
            if (r.success) {
                this.state.molblocks = r.molblocks;
                this.viewer.buildMolecules(this.state.molblocks);
                this.setStatus(r.message, '#4f6ef7');
            }
        } catch (e) {
            this.setStatus('Reset failed: ' + e.message, '#e0245e');
        }
    }

    async nextDockingPose() {
        this.pushHistory();
        this.dockVariant += 1;
        if (this.optCard.nextPoseBtn) this.optCard.nextPoseBtn.disabled = true;
        try {
            const payload = {
                action_type: 'auto_dock',
                variant: this.dockVariant,
                reaction_override: this.state.overridden ? this.state.rxnBonds : null
            };
            const res = await fetch('/api/step/step_3_ts/custom', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });
            const r = await res.json();
            if (r.success) {
                this.state.molblocks = r.molblocks;
                this.viewer.buildMolecules(r.molblocks);
                this.viewer.updateDynamicBonds(this.state.rxnBonds);
                this.viewer.flashAtoms(0x4f6ef7);
                if (this.onGeometryUpdated) this.onGeometryUpdated();
                this.setStatus(r.message + ' (Reactant docking input — click "Run / Refine TS" to optimize)', '#4f6ef7');
            } else {
                this.setStatus(r.error, '#e0245e');
            }
        } finally {
            if (this.optCard.nextPoseBtn) this.optCard.nextPoseBtn.disabled = false;
        }
    }

    async runOptimization(strategy, onDone) {
        this.pushHistory();
        const runBtn = this.optCard.runBtn;
        const origLabel = runBtn ? runBtn.textContent : '';
        if (runBtn) { runBtn.disabled = true; runBtn.textContent = 'Searching…'; }
        this.setStatus(`Running ${strategy.toUpperCase()} TS optimization...`);
        if (this.optCard.logBox) this.optCard.logBox.value = `Starting ${strategy.toUpperCase()} TS optimization…\n`;

        const logTimer = setInterval(() => this.fetchLog(), 700);
        const currentData = this.state.molblocks.map((mb, fi) => {
            const g = this.viewer.moleculeGroups[fi];
            return g ? updateMolBlockPositions(mb, g) : mb;
        });

        try {
            const payload = {
                action_type: 'optimize_ts',
                strategy,
                docked_molblocks: currentData,
                reaction_override: this.state.overridden ? this.state.rxnBonds : null
            };
            const res = await fetch('/api/step/step_3_ts/custom', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(payload)
            });
            const result = await res.json();
            await this.fetchLog();
            this._handleOptResult(result);
        } catch (e) {
            this.setStatus('Request failed: ' + e.message, '#e0245e');
        } finally {
            clearInterval(logTimer);
            if (runBtn) { runBtn.disabled = false; runBtn.textContent = origLabel; }
            if (onDone) onDone();
        }
    }

    _handleOptResult(result) {
        if (result.success) {
            this.setStatus(result.message, '#18a558');
            this.state.molblocks = result.optimized_molblocks;
            this.viewer.buildMolecules(this.state.molblocks);
            this.viewer.updateDynamicBonds(this.state.rxnBonds);
            this.viewer.flashAtoms(0x22c55e);
            if (this.onGeometryUpdated) this.onGeometryUpdated();
            if (result.profile && result.profile.length > 1 && this.optCard.profCanvas) {
                this.optCard.profCanvas.style.display = 'block';
                drawProfile(this.optCard.profCanvas, result.profile);
            }
            return;
        }
        if (result.diverged) {
            this.setStatus('⚠️ ' + result.error, '#d97706');
            if (this.viewer.flashAtoms) this.viewer.flashAtoms(0xd97706);
            return;
        }
        this.setStatus('❌ Error: ' + (result.error || 'Optimization failed'), '#e0245e');
    }
}
