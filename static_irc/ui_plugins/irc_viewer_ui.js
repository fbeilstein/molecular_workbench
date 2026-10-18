import { IRC3DViewer } from './irc_viewer_3d.js';
import { renderEnergyProfile } from './irc_profile_plot.js';
import { buildIRCLayout } from './irc_dom_builder.js';

let activeViewer = null;
let currentTrajectory = null;
let playbackTimer = null;
let currentFrameIndex = 0;
let isPlaying = false;

export const IRCViewerUI = {
    render: (container, ui, state) => {
        container.innerHTML = '';
        currentTrajectory = ui.trajectory || null;

        const dom = buildIRCLayout(container);
        activeViewer = new IRC3DViewer(dom.viewerBox);

        function setFrame(idx) {
            if (!currentTrajectory || !currentTrajectory.frames || currentTrajectory.frames.length === 0) return;
            const n = currentTrajectory.frames.length;
            currentFrameIndex = Math.max(0, Math.min(n - 1, idx));
            const frame = currentTrajectory.frames[currentFrameIndex];

            activeViewer.setFrame(currentFrameIndex);
            dom.frameSlider.value = currentFrameIndex.toString();
            dom.frameCountLabel.textContent = `Frame ${currentFrameIndex + 1}/${n}`;

            const dirTag = currentFrameIndex === currentTrajectory.ts_frame_index ? '🔴 Transition State' :
                           currentFrameIndex < currentTrajectory.ts_frame_index ? '🔵 Reactant Basin' : '🟢 Product Basin';
            dom.hudOverlay.textContent = `${dirTag} | s = ${frame.s_arc.toFixed(2)} amu½·Å | ΔE = ${frame.rel_energy_kcal >= 0 ? '+' : ''}${frame.rel_energy_kcal.toFixed(2)} kcal/mol`;

            renderEnergyProfile(dom.energyCanvas, currentTrajectory, currentFrameIndex, clickedIdx => setFrame(clickedIdx));
        }

        function loadTrajectory(traj) {
            if (!traj || !traj.frames || traj.frames.length === 0) return;
            currentTrajectory = traj;
            const n = traj.frames.length;

            dom.frameSlider.max = (n - 1).toString();
            dom.frameSlider.disabled = false;
            dom.playBtn.disabled = false;
            dom.prevBtn.disabled = false;
            dom.nextBtn.disabled = false;
            dom.downloadBtn.disabled = false;
            if (dom.prolongBtn) dom.prolongBtn.disabled = false;

            const b = traj.barriers || {};
            dom.fwdBadge.textContent = `ΔE‡ fwd: ${b.forward_barrier_kcal != null ? b.forward_barrier_kcal.toFixed(1) : '--'} kcal`;
            dom.revBadge.textContent = `ΔE‡ rev: ${b.reverse_barrier_kcal != null ? b.reverse_barrier_kcal.toFixed(1) : '--'} kcal`;
            dom.rxnBadge.textContent = `ΔE rxn: ${b.reaction_energy_kcal != null ? (b.reaction_energy_kcal > 0 ? '+' : '') + b.reaction_energy_kcal.toFixed(1) : '--'} kcal`;

            const tsIdx = traj.ts_frame_index || 0;
            activeViewer.loadTrajectory(traj);
            setFrame(tsIdx);
        }

        // Timeline controls
        dom.frameSlider.oninput = () => setFrame(parseInt(dom.frameSlider.value));
        dom.prevBtn.onclick = () => setFrame(currentFrameIndex - 1);
        dom.nextBtn.onclick = () => setFrame(currentFrameIndex + 1);

        dom.playBtn.onclick = () => {
            if (isPlaying) {
                clearInterval(playbackTimer);
                isPlaying = false;
                dom.playBtn.textContent = '▶ Play';
            } else {
                isPlaying = true;
                dom.playBtn.textContent = '⏸ Pause';
                playbackTimer = setInterval(() => {
                    if (!currentTrajectory) return;
                    let next = currentFrameIndex + 1;
                    if (next >= currentTrajectory.frames.length) {
                        if (dom.loopCheckbox.checked) next = 0;
                        else {
                            clearInterval(playbackTimer);
                            isPlaying = false;
                            dom.playBtn.textContent = '▶ Play';
                            return;
                        }
                    }
                    setFrame(next);
                }, 80);
            }
        };

        // Download JSON Trajectory
        dom.downloadBtn.onclick = () => {
            if (!currentTrajectory) return;
            const dataStr = 'data:text/json;charset=utf-8,' + encodeURIComponent(JSON.stringify(currentTrajectory, null, 2));
            const a = document.createElement('a');
            a.setAttribute('href', dataStr);
            a.setAttribute('download', 'irc_trajectory.json');
            document.body.appendChild(a);
            a.click();
            a.remove();
        };

        if (dom.prolongBtn) {
            dom.prolongBtn.onclick = async () => {
                if (!confirm('Prolong trajectory: Stitch previous reaction trajectory with this current IRC trajectory?')) return;
                dom.prolongBtn.disabled = true;
                dom.prolongBtn.textContent = 'Stitching...';
                try {
                    const res = await fetch('/api/step/step_4_irc/custom', {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ action_type: 'prolong_trajectory' })
                    });
                    const r = await res.json();
                    if (r.success && r.trajectory) {
                        loadTrajectory(r.trajectory);
                        dom.logBox.value += `\n✓ ${r.message}`;
                    } else {
                        alert('Prolong Error: ' + (r.error || 'Failed to prolong trajectory'));
                    }
                } catch (err) {
                    alert('Prolong Request failed: ' + err.message);
                } finally {
                    dom.prolongBtn.disabled = false;
                    dom.prolongBtn.textContent = '🔗 Prolong Pathway';
                }
            };
        }

        const curSolv = state.solvent || ui.solvent || 'water';
        if (dom.solventSelect) {
            dom.solventSelect.value = curSolv;
            dom.solventSelect.onchange = async () => {
                const s = dom.solventSelect.value;
                state.solvent = s;
                const g = document.getElementById('global-solvent-select');
                if (g) g.value = s;
                await fetch('/api/environment', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ solvent: s })
                });
            };
        }

        // Calculate IRC Action
        dom.calcBtn.onclick = async () => {
            dom.calcBtn.disabled = true;
            dom.calcBtn.textContent = 'Integrating IRC...';
            dom.logBox.value = 'Starting Intrinsic Reaction Coordinate integration...\n';

            const pollId = setInterval(async () => {
                try {
                    const lRes = await fetch('/api/log');
                    const text = await lRes.text();
                    dom.logBox.value = text;
                    dom.logBox.scrollTop = dom.logBox.scrollHeight;
                } catch (e) {}
            }, 800);

            try {
                const res = await fetch('/api/step/step_4_irc/custom', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'calculate_irc',
                        max_steps: 30,
                        step_size: 0.08,
                        ts_molblocks: ui.ts_molblocks,
                        engine: dom.engineSelect ? dom.engineSelect.value : 'xtb',
                        solvent: dom.solventSelect.value,
                        relax_fragments: dom.relaxCheckbox.checked
                    })
                });
                const r = await res.json();
                clearInterval(pollId);
                const finalLogRes = await fetch('/api/log');
                dom.logBox.value = await finalLogRes.text();
                dom.logBox.scrollTop = dom.logBox.scrollHeight;

                if (r.success && r.trajectory) {
                    loadTrajectory(r.trajectory);
                    dom.logBox.value += `\n✓ ${r.message}`;
                } else {
                    alert('IRC Error: ' + (r.error || 'Failed to compute trajectory'));
                }
            } catch (err) {
                clearInterval(pollId);
                alert('IRC Request failed: ' + err.message);
            } finally {
                dom.calcBtn.disabled = false;
                const engName = (dom.engineSelect && dom.engineSelect.value === 'orca') ? 'ORCA 6.1 DFT' : 'xTB';
                dom.calcBtn.textContent = `⚡ Re-Calculate IRC Path (${engName})`;
            }
        };

        if (dom.engineSelect) {
            dom.engineSelect.onchange = () => {
                const engName = dom.engineSelect.value === 'orca' ? 'ORCA 6.1 DFT' : 'xTB';
                dom.calcBtn.textContent = `⚡ Calculate IRC Path (${engName})`;
            };
        }

        if (currentTrajectory) {
            loadTrajectory(currentTrajectory);
        }
    },

    getData: async () => ({
        trajectory: currentTrajectory
    })
};
