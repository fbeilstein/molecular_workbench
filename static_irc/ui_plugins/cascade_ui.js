import { IRC3DViewer } from './irc_viewer_3d.js';
import { renderEnergyProfile } from './irc_profile_plot.js';
import { buildCascadeLayout } from './cascade_dom_builder.js';

let activeViewer = null;
let currentMegajson = null;
let playbackTimer = null;
let currentFrameIndex = 0;
let isPlaying = false;
let currentJunctionChoices = {};
let currentBridgeFrames = 8;

export const CascadeUI = {
    render: (container, ui, state, stepId) => {
        container.innerHTML = '';
        currentMegajson = ui.megajson || null;
        currentBridgeFrames = ui.n_bridge_frames || 8;
        currentJunctionChoices = Object.assign({}, ui.junction_choices || {});

        const dom = buildCascadeLayout(container, ui, { initialBridgeFrames: currentBridgeFrames });

        // Populate Sequence List
        const trajs = ui.trajectories || [];
        if (trajs.length === 0) {
            const noSeq = document.createElement('div');
            noSeq.style.fontSize = '12px';
            noSeq.style.color = '#64748b';
            noSeq.textContent = 'No previous trajectories found. Complete Step 4 first.';
            dom.seqList.appendChild(noSeq);
        } else {
            trajs.forEach((t, i) => {
                const sRow = document.createElement('div');
                sRow.style.display = 'flex';
                sRow.style.alignItems = 'center';
                sRow.style.gap = '8px';
                sRow.style.fontSize = '12px';
                sRow.style.padding = '6px 8px';
                sRow.style.background = '#ffffff';
                sRow.style.border = '1px solid #e2e8f0';
                sRow.style.borderRadius = '6px';

                const badge = document.createElement('span');
                badge.style.padding = '2px 6px';
                badge.style.borderRadius = '4px';
                badge.style.fontSize = '11px';
                badge.style.fontWeight = '600';
                badge.style.background = '#dbeafe';
                badge.style.color = '#1e40af';
                badge.textContent = `Step ${i + 1}`;

                const sName = document.createElement('span');
                sName.style.flex = '1';
                sName.style.fontFamily = 'monospace';
                sName.style.overflow = 'hidden';
                sName.style.textOverflow = 'ellipsis';
                sName.style.whiteSpace = 'nowrap';
                sName.title = t.reaction_smiles || t.filename;
                sName.textContent = t.reaction_smiles || t.filename;

                const fCount = document.createElement('span');
                fCount.style.color = '#64748b';
                fCount.style.fontSize = '11px';
                fCount.textContent = `${t.n_frames}f`;

                sRow.appendChild(badge);
                sRow.appendChild(sName);
                sRow.appendChild(fCount);
                dom.seqList.appendChild(sRow);
            });
        }

        // Populate Junctions
        function populateJunctions(reports) {
            dom.juncList.innerHTML = '';
            if (!reports || reports.length === 0) {
                const noJ = document.createElement('div');
                noJ.style.fontSize = '12px';
                noJ.style.color = '#94a3b8';
                noJ.textContent = 'No intermediate junctions detected.';
                dom.juncList.appendChild(noJ);
                return;
            }
            reports.forEach((rep, jIdx) => {
                const jBox = document.createElement('div');
                jBox.style.padding = '8px';
                jBox.style.background = '#ffffff';
                jBox.style.border = '1px solid #cbd5e1';
                jBox.style.borderRadius = '6px';

                const jHdr = document.createElement('div');
                jHdr.style.display = 'flex';
                jHdr.style.justifyContent = 'space-between';
                jHdr.style.fontSize = '12px';
                jHdr.style.fontWeight = '600';
                jHdr.style.marginBottom = '4px';
                jHdr.textContent = `Junction ${jIdx + 1}➔${jIdx + 2}: ${rep.hero_formula || 'Intermediate'}`;

                const sel = document.createElement('select');
                sel.style.width = '100%';
                sel.style.padding = '5px 8px';
                sel.style.fontSize = '11.5px';
                sel.style.borderRadius = '4px';
                sel.style.border = '1px solid #cbd5e1';

                (rep.all_candidates || []).forEach(cand => {
                    const opt = document.createElement('option');
                    opt.value = cand.id;
                    opt.textContent = `${cand.label} (RMSD: ${cand.rmsd != null ? cand.rmsd.toFixed(4) : '--'} Å)`;
                    if (currentJunctionChoices[jIdx] === cand.id || (!currentJunctionChoices[jIdx] && cand.id.includes('scaffold'))) {
                        opt.selected = true;
                        currentJunctionChoices[jIdx] = cand.id;
                    }
                    sel.appendChild(opt);
                });

                sel.onchange = () => {
                    currentJunctionChoices[jIdx] = sel.value;
                };

                jBox.appendChild(jHdr);
                jBox.appendChild(sel);
                dom.juncList.appendChild(jBox);
            });
        }

        populateJunctions((currentMegajson && currentMegajson.stitching_alignment_report) || ui.stitching_report || []);

        dom.bSlider.oninput = () => {
            currentBridgeFrames = parseInt(dom.bSlider.value);
            dom.bVal.textContent = `${currentBridgeFrames} Frames`;
        };

        // Wire 3D Viewer & Controls
        activeViewer = new IRC3DViewer(dom.viewerBox);

        function setFrame(idx) {
            if (!currentMegajson || !currentMegajson.frames || currentMegajson.frames.length === 0) return;
            const n = currentMegajson.frames.length;
            currentFrameIndex = Math.max(0, Math.min(n - 1, idx));
            const frame = currentMegajson.frames[currentFrameIndex];

            activeViewer.setFrame(currentFrameIndex);
            dom.frameSlider.value = currentFrameIndex.toString();
            dom.frameLabel.textContent = `Frame ${currentFrameIndex + 1}/${n}`;

            const stepTag = frame.step_index ? (typeof frame.step_index === 'number' ? `Step ${frame.step_index}` : frame.step_index) : 'Path';
            const tsTag = frame.is_ts ? ' 🔴 TS' : '';
            const eVal = frame.rel_energy_kcal != null ? frame.rel_energy_kcal : (frame.energy_kcal || 0);
            dom.hudOverlay.textContent = `${stepTag}${tsTag} | Frame ${currentFrameIndex + 1}/${n} | ΔE = ${eVal >= 0 ? '+' : ''}${eVal.toFixed(2)} kcal/mol`;

            renderEnergyProfile(dom.energyCanvas, currentMegajson, currentFrameIndex, clicked => setFrame(clicked));
        }

        function loadMegajson(mj) {
            if (!mj || !mj.frames || mj.frames.length === 0) return;
            currentMegajson = mj;
            const n = mj.frames.length;

            dom.frameSlider.max = (n - 1).toString();
            dom.frameSlider.disabled = false;
            dom.playBtn.disabled = false;
            dom.prevBtn.disabled = false;
            dom.nextBtn.disabled = false;

            dom.mFrames.textContent = `${n}`;
            dom.mSteps.textContent = `${mj.n_steps || 2}`;
            if (mj.overall_barrier_kcal != null) dom.mBar.textContent = `+${mj.overall_barrier_kcal} kcal`;
            if (mj.overall_reaction_energy_kcal != null) dom.mRxn.textContent = `${mj.overall_reaction_energy_kcal} kcal`;

            activeViewer.loadTrajectory(mj);
            setFrame(0);
        }

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
                    if (!currentMegajson) return;
                    let next = currentFrameIndex + 1;
                    if (next >= currentMegajson.frames.length) next = 0;
                    setFrame(next);
                }, 75);
            }
        };

        dom.restitchBtn.onclick = async () => {
            dom.restitchBtn.disabled = true;
            dom.restitchBtn.textContent = 'Stitching...';
            try {
                const res = await fetch('/api/step/step_5_cascade/custom', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'stitch_cascade',
                        junction_choices: currentJunctionChoices,
                        n_bridge_frames: currentBridgeFrames
                    })
                });
                const r = await res.json();
                if (r.success && r.megajson) {
                    loadMegajson(r.megajson);
                    if (r.stitching_report) populateJunctions(r.stitching_report);
                } else {
                    alert('Stitching Error: ' + (r.error || 'Failed to stitch'));
                }
            } catch (err) {
                alert('Request failed: ' + err.message);
            } finally {
                dom.restitchBtn.disabled = false;
                dom.restitchBtn.textContent = '⚡ Re-Stitch Pathway';
            }
        };

        dom.exportBtn.onclick = () => {
            if (!currentMegajson) return;
            const dataStr = 'data:text/json;charset=utf-8,' + encodeURIComponent(JSON.stringify(currentMegajson, null, 2));
            const a = document.createElement('a');
            a.setAttribute('href', dataStr);
            a.setAttribute('download', 'master_reaction_cascade.json');
            document.body.appendChild(a);
            a.click();
            a.remove();
        };

        if (currentMegajson) {
            loadMegajson(currentMegajson);
        }
    },

    getData: async () => ({
        megajson: currentMegajson
    })
};
