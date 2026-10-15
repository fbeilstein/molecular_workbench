/**
 * Dual Direction Transport & Playback Controller for Step 5 Bundle Viewer.
 */

export class BundleTransportController {
    constructor(state, callbacks) {
        this.state = state;
        this.callbacks = callbacks || {};
    }

    setDirection(dir) {
        this.state.playDirection = dir;
        const badge = document.querySelector('#badge-playback-dir');
        if (badge) {
            badge.className = `dir-badge ${dir === 'forward' ? 'dir-fwd' : 'dir-rev'}`;
            badge.textContent = dir === 'forward' ? 'Forward ⏩' : 'Reverse ⏪';
        }
    }

    playForward() {
        this.pause();
        this.setDirection('forward');
        const total = this.state.megajson.frames.length;
        if (this.state.currentFrame >= total - 1) this.seekFrame(0);

        this.state.isPlaying = true;
        const interval = Math.max(25, Math.round(110 / this.state.playbackSpeed));
        this.state.playTimer = setInterval(() => {
            if (this.state.currentFrame >= total - 1) {
                this.pause();
            } else {
                this.seekFrame(this.state.currentFrame + 1);
            }
        }, interval);
    }

    playReverse() {
        this.pause();
        this.setDirection('reverse');
        const total = this.state.megajson.frames.length;
        if (this.state.currentFrame <= 0) this.seekFrame(total - 1);

        this.state.isPlaying = true;
        const interval = Math.max(25, Math.round(110 / this.state.playbackSpeed));
        this.state.playTimer = setInterval(() => {
            if (this.state.currentFrame <= 0) {
                this.pause();
            } else {
                this.seekFrame(this.state.currentFrame - 1);
            }
        }, interval);
    }

    pause() {
        this.state.isPlaying = false;
        if (this.state.playTimer) {
            clearInterval(this.state.playTimer);
            this.state.playTimer = null;
        }
    }

    playSegment(startFrame, endFrame, direction = 'forward') {
        this.pause();
        this.setDirection(direction);

        const total = (this.state.megajson && this.state.megajson.frames) ? this.state.megajson.frames.length : 147;
        const f1 = Math.max(0, Math.min(total - 1, parseInt(startFrame, 10) || 0));
        const f2 = Math.max(0, Math.min(total - 1, parseInt(endFrame, 10) || 0));
        const minF = Math.min(f1, f2);
        const maxF = Math.max(f1, f2);

        const fromF = (direction === 'reverse') ? maxF : minF;
        const toF = (direction === 'reverse') ? minF : maxF;

        this.seekFrame(fromF);
        this.state.isPlaying = true;

        const interval = Math.max(25, Math.round(95 / this.state.playbackSpeed));
        this.state.playTimer = setInterval(() => {
            if (direction === 'forward') {
                if (this.state.currentFrame >= toF) {
                    this.pause();
                } else {
                    this.seekFrame(this.state.currentFrame + 1);
                }
            } else {
                if (this.state.currentFrame <= toF) {
                    this.pause();
                } else {
                    this.seekFrame(this.state.currentFrame - 1);
                }
            }
        }, interval);
    }

    async seekFrame(f) {
        this.state.currentFrame = f;
        if (this.state.pesChart) this.state.pesChart.setCurrentFrame(f);

        const scrubber = document.querySelector('#bundle-frame-scrubber');
        if (scrubber) scrubber.value = f;

        const readout = document.querySelector('#bundle-frame-readout');
        if (readout && this.state.megajson) {
            readout.textContent = `Frame ${f} / ${this.state.megajson.frames.length - 1}`;
        }

        if (this.state.svgManager) {
            this.state.svgManager.updatePushArrows(f, this.state.playDirection);
        }

        if (this.callbacks.onFrameSeek) {
            await this.callbacks.onFrameSeek(f);
        }
    }
}
