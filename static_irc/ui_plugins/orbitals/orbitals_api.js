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

export async function fetchOrbitalCube(stepId, cubePromises, cubeCache, frameIdx, orbKey, chemKey, isovalue) {
    const k = `${frameIdx}_${orbKey}`;
    if (cubeCache && cubeCache[k]) return { success: true, cube: cubeCache[k] };
    if (cubePromises.has(k)) return await cubePromises.get(k);

    const promise = (async () => {
        try {
            const res = await fetch(`/api/step/${stepId}/custom`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ action_type: 'get_cube', frame_index: frameIdx, orbital: orbKey, chem_key: chemKey, isovalue })
            });
            const data = await res.json();
            if (data && data.success && data.cube && cubeCache) {
                cubeCache[k] = data.cube;
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
        await fetch(`/api/step/${stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'save_timeline', timeline })
        });
    } catch (e) {}
}
