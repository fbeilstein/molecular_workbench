/**
 * DOM Event Wiring for Step 5 & Step 7 Bundle Viewer.
 */

import { exportMaxCompBundleZip } from './bundle_export.js';
import { saveSvgToBackend } from './bundle_preload_sync.js';

export function bindBundleEvents(container, ui) {
    const s = ui.state;
    const t = s.transport;

    const svgInput = container.querySelector('#input-upload-svg');
    svgInput?.addEventListener('change', (e) => {
        const file = e.target.files[0];
        if (!file) return;
        const reader = new FileReader();
        reader.onload = async () => {
            try {
                const svgText = reader.result;
                if (s.svgManager.loadSvg(svgText)) {
                    if (s.megajson) s.megajson.mechanism_svg = svgText;
                    ui.renderGroupsTree();
                    await saveSvgToBackend(s, svgText);
                    const hint = container.querySelector('#label-wiring-hint');
                    if (hint) { hint.textContent = `✓ Loaded ${file.name}`; hint.style.color = '#10b981'; }
                }
            } catch (err) { alert('Error reading SVG: ' + err.message); }
            finally { svgInput.value = ''; }
        };
        reader.readAsText(file);
    });

    const tabInteractive = container.querySelector('#tab-mode-interactive'), tabWiring = container.querySelector('#tab-mode-wiring');
    const hint = container.querySelector('#label-wiring-hint'), wiringStatus = container.querySelector('#label-wiring-status');
    const setMode = (isWiring) => {
        tabWiring.className = isWiring ? 'btn-primary active' : 'btn-secondary';
        tabInteractive.className = !isWiring ? 'btn-primary active' : 'btn-secondary';
        hint.textContent = isWiring ? 'Click element in tree/canvas to wire' : 'Click arrows to play • Right-click molecules';
        if (wiringStatus) { wiringStatus.textContent = isWiring ? 'Wiring Mode' : 'Hierarchy'; wiringStatus.className = isWiring ? 'chip amber' : 'chip blue'; }
        s.svgManager.setMode(isWiring ? 'edit' : 'interactive');
        if (isWiring) ui.renderGroupsTree();
    };
    tabInteractive?.addEventListener('click', () => setMode(false));
    tabWiring?.addEventListener('click', () => setMode(true));

    container.querySelector('#btn-step-prev')?.addEventListener('click', () => t.stepBackward());
    container.querySelector('#btn-play-rev')?.addEventListener('click', () => t.playReverse());
    container.querySelector('#btn-playback-pause')?.addEventListener('click', () => t.pause());
    container.querySelector('#btn-play-fwd')?.addEventListener('click', () => t.playForward());
    container.querySelector('#btn-step-next')?.addEventListener('click', () => t.stepForward());

    container.querySelector('#sel-playback-speed')?.addEventListener('change', (e) => {
        s.playbackSpeed = parseFloat(e.target.value);
        if (s.isPlaying) { t.pause(); if (s.playDirection === 'forward') t.playForward(); else t.playReverse(); }
    });

    container.querySelector('#bundle-frame-scrubber')?.addEventListener('input', (e) => {
        t.pause();
        const newF = parseInt(e.target.value, 10);
        t.setDirection(newF > s.currentFrame ? 'forward' : 'reverse');
        t.seekFrame(newF);
    });

    container.querySelector('#btn-recompute-highres')?.addEventListener('click', () => s.highresManager.renderModal(container, s.tracks));
    container.querySelector('#btn-export-bundle-zip')?.addEventListener('click', () => exportMaxCompBundleZip(s.stepId, { isovalue: s.isovalue, esp_surfaces: s.viewer3D?.espMeshCache || {} }));
    container.querySelector('#btn-bundle-zoom-in')?.addEventListener('click', () => s.viewer3D?.zoomIn());
    container.querySelector('#btn-bundle-zoom-out')?.addEventListener('click', () => s.viewer3D?.zoomOut());
    container.querySelector('#btn-bundle-zoom-fit')?.addEventListener('click', () => s.viewer3D?.resetCamera());
    container.querySelector('#chk-bundle-charges')?.addEventListener('change', (e) => { s.showCharges = e.target.checked; ui.refreshViewer(); });
    container.querySelector('#chk-bundle-esp-surface')?.addEventListener('change', (e) => {
        s.showEspSurface = e.target.checked;
        if (e.target.checked) {
            s.activeOrbitals.add('track_esp_charge');
            s.viewer3D?.prewarmEspCache(s.trajectorySymbols || s.megajson.symbols, s.megajson.frames, s.megajson.trajectory_charges);
        } else s.activeOrbitals.delete('track_esp_charge');
        ui.renderTracksList(); ui.refreshViewer();
    });
    container.querySelector('#chk-bundle-mute')?.addEventListener('change', (e) => { s.hideAllOrbitals = e.target.checked; ui.refreshViewer(); });
}
