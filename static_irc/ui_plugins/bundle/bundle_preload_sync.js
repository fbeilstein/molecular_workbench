/**
 * Cube and Mesh Preloading & Backend Sync Engine for Bundle Viewer.
 */

export async function preloadAllCachedCubes(state, force = false, refreshCallback) {
    if (!state.megajson || !state.viewer3D) return;
    const total = state.megajson.frames.length;
    const targets = state.tracks.filter(t => t.orbital_key && t.orbital_key !== 'track_esp_charge');
    if (!targets.length) return;

    const iso = state.isovalue;
    const isoSuffix = `_${Math.round(iso * 1000)}`;
    const seq = ++state._preloadSeq || (state._preloadSeq = 1);

    if (force) {
        state.viewer3D?.setStatus(`Remeshing all frames for isovalue ${iso.toFixed(3)}...`);
        fetch(`/api/step/${state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'set_isovalue', isovalue: iso })
        }).catch(() => {});
    }

    const batchSize = 35;
    for (let start = 0; start < total; start += batchSize) {
        if (state._preloadSeq !== seq) return;
        const end = Math.min(total, start + batchSize);
        const frames = [];
        for (let i = start; i < end; i++) {
            const need = force || targets.some(t => {
                const ok = t.orbital_key;
                const isoKey = `${i}_${ok}${isoSuffix}`;
                return !state.viewer3D.meshCache[isoKey] && !state.viewer3D.meshCache[`${i}_${ok}`];
            });
            if (need) frames.push(i);
        }
        if (frames.length > 0) {
            try {
                const res = await fetch(`/api/step/${state.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({ action_type: 'get_cached_meshes', frames, tracks: targets, isovalue: iso })
                });
                const d = await res.json();
                if (state._preloadSeq !== seq) return;
                if (d && d.success) {
                    if (d.meshes) {
                        for (const [k, v] of Object.entries(d.meshes)) {
                            state.viewer3D.meshCache[k] = v;
                            if (!k.endsWith(isoSuffix)) state.viewer3D.meshCache[`${k}${isoSuffix}`] = v;
                        }
                    }
                    if (d.cubes) Object.assign(state.viewer3D.cubeCache, d.cubes);
                    if (frames.includes(state.currentFrame) && refreshCallback) refreshCallback();
                }
            } catch (e) {}
        }
    }

    if (state._preloadSeq === seq && force) {
        state.viewer3D?.setStatus(`✓ All ${total} frames remeshed (iso: ${iso.toFixed(3)})`);
        setTimeout(() => { if (state._preloadSeq === seq) state.viewer3D?.setStatus('Ready'); }, 2000);
    }
}

export async function fetchCube(state, frameIdx, orbKey, chemKey, isovalue) {
    const res = await fetch(`/api/step/${state.stepId}/custom`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({
            action_type: 'get_orbital_mesh', frame_index: frameIdx,
            orbital: orbKey, chem_key: chemKey, isovalue: isovalue
        })
    });
    return await res.json();
}

export async function saveSvgToBackend(state, svgText) {
    try {
        await fetch(`/api/step/${state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'save_mechanism_svg', svg_text: svgText })
        });
    } catch (err) {}
}

export async function saveWiringToBackend(state, wiring) {
    try {
        await fetch(`/api/step/${state.stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'save_svg_wiring', wiring })
        });
    } catch (err) {}
}
