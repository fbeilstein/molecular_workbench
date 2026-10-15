/**
 * Batch precomputation controller and progress bar for Reaction Orbital Studio.
 */

export class OrbitalsPrecompute {
    constructor(stepId, getTargetsFn, onProgress, onComplete) {
        this.stepId = stepId;
        this.getTargetsFn = getTargetsFn;
        this.onProgress = onProgress;
        this.onComplete = onComplete;
        this.isPrecomputing = false;
        this.abortFlag = false;
    }

    async startPrecompute(totalFrames) {
        if (this.isPrecomputing) return;
        this.isPrecomputing = true;
        this.abortFlag = false;

        let targets = this.getTargetsFn ? this.getTargetsFn() : [];
        if (!targets || targets.length === 0) {
            targets = [{ orbital_key: 'HOMO', chem_key: 'HOMO' }, { orbital_key: 'LUMO', chem_key: 'LUMO' }];
        }

        let currentFrame = 0;
        const batchSize = 2;

        if (this.onProgress) {
            this.onProgress({ percent: 1, currentFrame: 0, totalFrames, isComplete: false });
        }

        while (currentFrame < totalFrames && !this.abortFlag) {
            try {
                const res = await fetch(`/api/step/${this.stepId}/custom`, {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'batch_precompute_cubes',
                        start_frame: currentFrame,
                        batch_size: batchSize,
                        selected_tracks: targets
                    })
                });
                const data = await res.json();
                if (!data.success) {
                    console.error('Batch precompute returned error:', data);
                    break;
                }

                currentFrame = data.end_frame + 1;
                const percent = Math.min(100, Math.round((currentFrame / totalFrames) * 100));

                if (this.onProgress) {
                    this.onProgress({
                        percent,
                        currentFrame,
                        totalFrames,
                        cubes: data.cubes || {},
                        isComplete: currentFrame >= totalFrames
                    });
                }

                if (data.is_complete || currentFrame >= totalFrames) {
                    break;
                }
            } catch (err) {
                console.error('Batch precompute network error:', err);
                break;
            }
        }

        this.isPrecomputing = false;
        if (this.onComplete) {
            this.onComplete(this.abortFlag);
        }
    }

    async preloadCached(cubeCache, totalFrames, onFrameLoaded) {
        const targets = this.getTargetsFn ? this.getTargetsFn() : [];
        if (!targets || targets.length === 0) return;
        const batchSize = 35;
        for (let start = 0; start < totalFrames; start += batchSize) {
            const end = Math.min(totalFrames, start + batchSize), frames = [];
            for (let i = start; i < end; i++) {
                const need = targets.some(t => {
                    const ok = t.orbital_key || t.key;
                    return !cubeCache[`${i}_${ok}`];
                });
                if (need) frames.push(i);
            }
            if (frames.length > 0) {
                try {
                    const res = await fetch(`/api/step/${this.stepId}/custom`, {
                        method: 'POST',
                        headers: { 'Content-Type': 'application/json' },
                        body: JSON.stringify({ action_type: 'get_cached_cubes', frames, tracks: targets })
                    });
                    const data = await res.json();
                    if (data && data.success && data.cubes) {
                        Object.assign(cubeCache, data.cubes);
                        if (onFrameLoaded) onFrameLoaded(frames);
                    }
                } catch (e) {}
            }
        }
    }

    abort() {
        if (this.isPrecomputing) {
            this.abortFlag = true;
        }
    }
}
