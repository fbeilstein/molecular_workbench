import { MoleculeViewer3D } from './common/molecule_viewer_3d.js';
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
        activeViewer = new MoleculeViewer3D(dom.viewerBox, {
            statusEl: dom.hudOverlay
        });
        activeViewer.init();

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
            const isGibbs = currentTrajectory.barriers && currentTrajectory.barriers.energy_type === 'gibbs';
            const eUnit = isGibbs ? 'ΔG' : 'ΔE';
            dom.hudOverlay.textContent = `${dirTag} | s = ${frame.s_arc.toFixed(2)} amu½·Å | ${eUnit} = ${frame.rel_energy_kcal >= 0 ? '+' : ''}${frame.rel_energy_kcal.toFixed(2)} kcal/mol`;

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

            const b = traj.barriers || {};
            const isGibbs = b.energy_type === 'gibbs';
            const unitTag = isGibbs ? 'ΔG' : 'ΔE';

            dom.fwdBadge.textContent = `${unitTag}‡ fwd: ${b.forward_barrier_kcal != null ? b.forward_barrier_kcal.toFixed(1) : '--'} kcal`;
            dom.revBadge.textContent = `${unitTag}‡ rev: ${b.reverse_barrier_kcal != null ? b.reverse_barrier_kcal.toFixed(1) : '--'} kcal`;
            dom.rxnBadge.textContent = `${unitTag} rxn: ${b.reaction_energy_kcal != null ? (b.reaction_energy_kcal > 0 ? '+' : '') + b.reaction_energy_kcal.toFixed(1) : '--'} kcal`;

            if (isGibbs) {
                dom.fwdBadge.title = 'Standard-State Gibbs Free Energy Barrier (298.15 K, 1 M in solvent)';
                dom.chartTitle.textContent = 'Gibbs Free Energy Profile ΔG (kcal/mol) vs Reaction Coordinate s:';
            } else {
                dom.fwdBadge.title = 'Electronic Energy Barrier';
                dom.chartTitle.textContent = 'Energy Profile ΔE (kcal/mol) vs Reaction Coordinate s:';
            }

            if (dom.freqBadge) {
                const imFreqs = b.imaginary_freq_cm || [];
                if (imFreqs.length > 0) {
                    dom.freqBadge.style.display = 'inline-block';
                    dom.freqBadge.textContent = `ν‡: ${imFreqs[0].toFixed(1)} cm⁻¹`;
                    dom.freqBadge.title = `Imaginary vibrational mode (${imFreqs[0].toFixed(1)} cm⁻¹) confirming TS saddle point`;
                } else {
                    dom.freqBadge.style.display = 'none';
                }
            }

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

        // ORCA Mode descriptions & button texts
        const orcaDescriptions = {
            sp: 'Calculates r2SCAN-3c single points across trajectory anchors with monotonic PCHIP interpolation (~2 min).',
            thermo: 'Computes exact Gibbs free energy G (ZPE, thermal entropy -TΔS, standard-state solution volume) & verifies 1 imaginary frequency (~10 min).',
            irc: 'Full ab initio DFT intrinsic reaction coordinate integration walking downhill directly on the r2SCAN-3c surface (~30-45 min).',
            neb: 'Climbing Image Nudged Elastic Band on DFT surface. Unattended overnight search for exact saddle point and MEP (2–6 hrs).'
        };

        const calcBtnLabels = {
            xtb: '⚡ Calculate IRC Path (xTB)',
            sp: '⚡ Calculate IRC Path (ORCA DFT Preview)',
            thermo: '🌡️ Calculate Thermochemical ΔG Profile (ORCA)',
            irc: '🎯 Run Native ORCA IRC (DFT)',
            neb: '🌙 Launch Overnight NEB-CI (ORCA)'
        };

        const runningBtnLabels = {
            xtb: 'Integrating IRC (xTB)...',
            sp: 'Calculating ORCA DFT Preview...',
            thermo: 'Calculating Thermochemical ΔG (ORCA Freq)...',
            irc: 'Integrating Native ORCA IRC (DFT)...',
            neb: 'Running Overnight NEB-CI (ORCA)...'
        };

        function updateEngineAndModeUI() {
            const isOrca = dom.engineSelect && dom.engineSelect.value === 'orca';
            if (dom.orcaModeRow) {
                dom.orcaModeRow.style.display = isOrca ? 'flex' : 'none';
            }
            const mode = isOrca && dom.orcaModeSelect ? dom.orcaModeSelect.value : 'xtb';
            if (dom.orcaDescBox && isOrca && orcaDescriptions[mode]) {
                dom.orcaDescBox.textContent = orcaDescriptions[mode];
            }
            if (dom.calcBtn) {
                dom.calcBtn.textContent = calcBtnLabels[mode] || calcBtnLabels.xtb;
            }
        }

        if (dom.engineSelect) {
            dom.engineSelect.onchange = () => updateEngineAndModeUI();
        }
        if (dom.orcaModeSelect) {
            dom.orcaModeSelect.onchange = () => updateEngineAndModeUI();
        }
        updateEngineAndModeUI();

        // Calculate IRC Action
        dom.calcBtn.onclick = async () => {
            const isOrca = dom.engineSelect && dom.engineSelect.value === 'orca';
            const modeKey = isOrca && dom.orcaModeSelect ? dom.orcaModeSelect.value : 'xtb';

            dom.calcBtn.disabled = true;
            dom.calcBtn.textContent = runningBtnLabels[modeKey] || 'Calculating...';
            dom.logBox.value = `Starting calculation: Engine=${dom.engineSelect ? dom.engineSelect.value : 'xtb'}, Mode=${modeKey}...\n`;

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
                        orca_mode: dom.orcaModeSelect ? dom.orcaModeSelect.value : 'sp',
                        solvent: dom.solventSelect ? dom.solventSelect.value : 'water',
                        relax_fragments: dom.relaxCheckbox ? dom.relaxCheckbox.checked : true
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
                    alert('Calculation Error: ' + (r.error || 'Failed to compute trajectory'));
                }
            } catch (err) {
                clearInterval(pollId);
                alert('Calculation Request failed: ' + err.message);
            } finally {
                dom.calcBtn.disabled = false;
                const activeMode = dom.engineSelect && dom.engineSelect.value === 'orca' && dom.orcaModeSelect ? dom.orcaModeSelect.value : 'xtb';
                dom.calcBtn.textContent = calcBtnLabels[activeMode] || calcBtnLabels.xtb;
            }
        };

        if (currentTrajectory) {
            loadTrajectory(currentTrajectory);
        }
    },

    getData: async () => ({
        trajectory: currentTrajectory
    })
};
