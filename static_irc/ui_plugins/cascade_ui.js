/**
 * Step 5: Master Cascade Stitcher and Boundary Overlap Inspector UI.
 * Sequence list and junction builders delegated to cascade_preview.js.
 */

import { MoleculeViewer3D } from './common/molecule_viewer_3d.js';
import { renderEnergyProfile } from './irc_profile_plot.js';
import { buildCascadeLayout } from './cascade_dom_builder.js';
import { renderSequenceList, populateJunctions } from './cascade_preview.js';

let activeViewer = null;
let currentMegajson = null;
let playbackTimer = null;
let currentFrameIndex = 0;
let isPlaying = false;
let currentJunctionChoices = {};
let currentBridgeFrames = 8;
let currentViewMode = 'cascade';
let activeJunctionIdx = 0;
let isOverlapAligned = true;
let isHeroOnly = true;

export const CascadeUI = {
    render: (container, ui, state, stepId) => {
        container.innerHTML = '';
        currentMegajson = ui.megajson || null;
        currentBridgeFrames = ui.n_bridge_frames || 8;
        currentJunctionChoices = Object.assign({}, ui.junction_choices || {});
        currentViewMode = 'cascade';
        activeJunctionIdx = 0;
        isOverlapAligned = true;

        const dom = buildCascadeLayout(container, ui, { initialBridgeFrames: currentBridgeFrames });

        const seqState = {
            currentTrajs: [...(ui.trajectories || [])],
            originalTrajs: [...(ui.trajectories || [])],
            orderModified: false
        };

        const refreshSeq = () => renderSequenceList(dom, seqState);
        refreshSeq();

        function getJunctionReports() {
            return (currentMegajson && (currentMegajson.stitching_alignment_report || currentMegajson.junction_previews)) || ui.stitching_report || [];
        }

        populateJunctions(dom, getJunctionReports(), currentJunctionChoices);

        dom.bSlider.oninput = () => {
            currentBridgeFrames = parseInt(dom.bSlider.value);
            dom.bVal.textContent = `${currentBridgeFrames} Frames`;
        };

        activeViewer = new MoleculeViewer3D(dom.viewerBox, { statusEl: dom.hudOverlay });
        activeViewer.init();

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

        function renderCurrentOverlap() {
            const reports = getJunctionReports();
            if (!reports || reports.length === 0) {
                dom.hudOverlay.textContent = 'No intermediate junctions available for overlap inspection.';
                return;
            }
            const j = reports[activeJunctionIdx] || reports[0];
            if (!j || !j.prev_symbols || !j.prev_coords) return;

            const nextCoords = isOverlapAligned ? (j.next_coords_aligned || j.next_coords) : (j.next_coords_unaligned || j.next_coords);
            const rmsd = j.rmsd || j.joint_rmsd_angstrom || 0.0;

            let pSyms = j.prev_symbols, pCoords = j.prev_coords;
            let nSyms = j.next_symbols || j.prev_symbols, nCoords = nextCoords;

            if (isHeroOnly && j.matched_prev_indices?.length && j.matched_next_indices?.length) {
                pSyms = j.matched_prev_indices.map(i => j.prev_symbols[i]);
                pCoords = j.matched_prev_indices.map(i => j.prev_coords[i]);
                nSyms = j.matched_next_indices.map(i => (j.next_symbols || j.prev_symbols)[i]);
                nCoords = j.matched_next_indices.map(i => nextCoords[i]);
            }

            activeViewer.renderOverlap(pSyms, pCoords, nSyms, nCoords, rmsd, isOverlapAligned);
            dom.rmsdBadge.textContent = `RMSD: ${rmsd.toFixed(4)} Å`;
            const focusTag = isHeroOnly ? `[Hero Only: ${j.hero_formula || 'Intermediate'}]` : '[Full System]';
            dom.hudOverlay.textContent = `Boundary Overlap | Junction ${activeJunctionIdx + 1}➔${activeJunctionIdx + 2}: ${j.hero_formula || 'Intermediate'} | RMSD: ${rmsd.toFixed(4)} Å ${focusTag} (${isOverlapAligned ? 'Cyan: Step K, Magenta: Step K+1' : 'Cyan: Step K, Yellow: Step K+1 Unaligned'})`;
        }

        function setViewMode(mode) {
            currentViewMode = mode;
            if (mode === 'cascade') {
                dom.btnModeCascade.style.cssText = 'background:#0284c7;color:#ffffff;border:none;';
                dom.btnModeOverlap.style.cssText = 'background:#f1f5f9;color:#334155;border:1px solid #cbd5e1;';
                dom.overlapSubControls.style.display = 'none';
                dom.ctrlRow.style.display = 'flex';
                dom.energyCanvas.style.display = 'block';
                if (currentMegajson) { activeViewer.loadTrajectory(currentMegajson); setFrame(currentFrameIndex); }
            } else {
                if (isPlaying) { clearInterval(playbackTimer); isPlaying = false; dom.playBtn.textContent = '▶ Play'; }
                dom.btnModeCascade.style.cssText = 'background:#f1f5f9;color:#334155;border:1px solid #cbd5e1;';
                dom.btnModeOverlap.style.cssText = 'background:#0284c7;color:#ffffff;border:none;';
                dom.overlapSubControls.style.display = 'flex';
                dom.ctrlRow.style.display = 'none';
                dom.energyCanvas.style.display = 'none';
                renderCurrentOverlap();
            }
        }

        dom.btnModeCascade.onclick = () => setViewMode('cascade');
        dom.btnModeOverlap.onclick = () => setViewMode('overlap');

        dom.juncSelect.onchange = () => {
            activeJunctionIdx = parseInt(dom.juncSelect.value) || 0;
            renderCurrentOverlap();
        };

        dom.alignTypeBtn.onclick = () => {
            isOverlapAligned = !isOverlapAligned;
            dom.alignTypeBtn.textContent = isOverlapAligned ? 'Aligned (Fit)' : 'Raw Unaligned';
            dom.alignTypeBtn.style.color = isOverlapAligned ? '#0284c7' : '#d97706';
            renderCurrentOverlap();
        };

        dom.heroFilterBtn.onclick = () => {
            isHeroOnly = !isHeroOnly;
            dom.heroFilterBtn.textContent = isHeroOnly ? 'Hero Only' : 'Show All';
            dom.heroFilterBtn.style.color = isHeroOnly ? '#059669' : '#64748b';
            renderCurrentOverlap();
        };

        function loadMegajson(mj) {
            if (!mj || !mj.frames || mj.frames.length === 0) return;
            currentMegajson = mj;
            const n = mj.frames.length;

            dom.frameSlider.max = (n - 1).toString();
            dom.frameSlider.disabled = dom.playBtn.disabled = dom.prevBtn.disabled = dom.nextBtn.disabled = false;
            dom.mFrames.textContent = `${n}`;
            dom.mSteps.textContent = `${mj.n_steps || 2}`;
            if (mj.overall_barrier_kcal != null) dom.mBar.textContent = `+${mj.overall_barrier_kcal} kcal`;
            if (mj.overall_reaction_energy_kcal != null) dom.mRxn.textContent = `${mj.overall_reaction_energy_kcal} kcal`;

            if (currentViewMode === 'cascade') { activeViewer.loadTrajectory(mj); setFrame(0); }
            else { renderCurrentOverlap(); }
        }

        dom.frameSlider.oninput = () => setFrame(parseInt(dom.frameSlider.value));
        dom.frameSlider.onchange = () => setFrame(parseInt(dom.frameSlider.value));
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
                        n_bridge_frames: currentBridgeFrames,
                        ordered_files: seqState.currentTrajs.map(t => t.step_id || t.filename),
                        reversed_steps: seqState.currentTrajs.filter(t => t.is_reversed).map(t => t.step_id || t.filename)
                    })
                });
                const r = await res.json();
                if (r.success && r.megajson) {
                    seqState.orderModified = false;
                    refreshSeq();
                    loadMegajson(r.megajson);
                    if (r.stitching_report) populateJunctions(dom, r.stitching_report, currentJunctionChoices);
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

        if (currentMegajson) {
            loadMegajson(currentMegajson);
        }
    },

    getData: async () => ({
        megajson: currentMegajson
    })
};
