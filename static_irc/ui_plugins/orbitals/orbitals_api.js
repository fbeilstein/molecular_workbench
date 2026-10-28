/**
 * Step 6: API and network communication for Reaction Orbital Studio.
 */
export async function fetchFrameData(stepId, frameIdx) {
    try {
        const res = await fetch(`/api/step/${stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'get_frame_orbitals', frame_index: frameIdx })
        });
        return await res.json();
    } catch (err) {
        console.error('Failed to load frame data:', err);
        return { success: false, error: err.message };
    }
}

export async function fetchOrbitalCube(stepId, cubePromises, cubeCache, frameIdx, orbKey, chemKey, isovalue, trackId) {
    const tid = trackId || orbKey;
    const k = `${frameIdx}_${tid}`;
    if (cubeCache && cubeCache[k]) return { success: true, cube: cubeCache[k] };
    if (cubePromises.has(k)) return await cubePromises.get(k);

    const promise = (async () => {
        try {
            const res = await fetch(`/api/step/${stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ action_type: 'get_cube', frame_index: frameIdx, track_id: tid, orbital: orbKey, chem_key: chemKey, isovalue })
            });
            const data = await res.json();
            if (data && data.success && data.cube && cubeCache) {
                cubeCache[k] = data.cube;
                if (data.resolved_orbital && data.resolved_orbital !== tid) {
                    cubeCache[`${frameIdx}_${data.resolved_orbital}`] = data.cube;
                }
            }
            return data;
        } catch (err) {
            console.error('fetchCube network error:', err);
            return { success: false, error: err.message };
        } finally {
            cubePromises.delete(k);
        }
    })();

    cubePromises.set(k, promise);
    return await promise;
}

export async function postTimeline(stepId, timeline) {
    try {
        const res = await fetch(`/api/step/${stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'save_timeline', timeline })
        });
        return await res.json();
    } catch (err) {
        console.error('Failed to save timeline:', err);
        return { success: false, error: err.message };
    }
}
