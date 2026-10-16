/**
 * viewer_app.js
 * Main coordinator orchestrating loader, 3Dmol viewer, static PES plot, and playback loop.
 */

import { loadReactionBundle } from './viewer_loader.js';
import { Viewer3D } from './viewer_3d.js';
import { ViewerPES } from './viewer_pes.js';
import { ViewerUI } from './viewer_ui.js';

class ReactionViewerApp {
  constructor() {
    this.bundle = null;
    this.currentFrame = 0;
    this.totalFrames = 0;
    this.isPlaying = false;
    this.playSpeed = 1.0;
    this.loop = true;
    this.playTimer = null;
    this.baseFps = 15;

    // Elements
    this.overlay = document.getElementById('overlay');
    this.overlayStatus = document.getElementById('overlay-status');
    this.progressBarContainer = document.getElementById('progress-bar-container');
    this.progressBar = document.getElementById('progress-bar');
    this.fileInput = document.getElementById('file-input');

    // Submodules
    this.viewer3d = new Viewer3D(document.getElementById('gl-canvas-container'));
    this.viewerPes = new ViewerPES(document.getElementById('pes-canvas'));
    this.ui = new ViewerUI(this);

    this.bindDropEvents();
    this.bindKeyboard();
  }

  init() {
    const params = new URLSearchParams(window.location.search);
    const bundleUrl = params.get('bundle');
    const initialFrame = parseInt(params.get('frame') || '0', 10);
    const autoplay = params.get('autoplay') === '1' || params.get('autoplay') === 'true';

    if (bundleUrl) {
      this.loadUrl(bundleUrl, initialFrame, autoplay);
    } else {
      this.showDropzone();
    }
  }

  showDropzone() {
    if (this.overlay) {
      this.overlay.style.display = 'flex';
      if (this.progressBarContainer) this.progressBarContainer.style.display = 'none';
      if (this.overlayStatus) this.overlayStatus.textContent = 'Drag and drop a .rxnbundle.zip here, or click to choose a file';
    }
  }

  showProgress(pct, statusText) {
    if (this.overlay) {
      this.overlay.style.display = 'flex';
      if (this.progressBarContainer) {
        this.progressBarContainer.style.display = 'block';
        if (this.progressBar) this.progressBar.style.width = `${pct}%`;
      }
      if (this.overlayStatus) this.overlayStatus.textContent = statusText;
    }
  }

  hideOverlay() {
    if (this.overlay) {
      this.overlay.style.display = 'none';
    }
  }

  async loadUrl(url, initialFrame = 0, autoplay = false) {
    try {
      this.showProgress(10, 'Fetching bundle...');
      const bundle = await loadReactionBundle(url, (pct, msg) => this.showProgress(pct, msg));
      this.onBundleLoaded(bundle, initialFrame, autoplay);
    } catch (err) {
      console.error(err);
      this.showDropzone();
      if (this.overlayStatus) {
        const dlBtn = (url.startsWith('http://') || url.startsWith('https://'))
          ? `<br><a href="${url}" target="_blank" rel="noopener" class="btn btn-primary" style="display:inline-block;margin:12px 0 6px 0;text-decoration:none;">⬇ Download Bundle Archive</a><br><span style="font-size:12px;color:var(--text-secondary);">After downloading, drop it here</span>`
          : '<br><br>Drag & drop the .zip file manually:';
        this.overlayStatus.innerHTML = `<span style="color:var(--danger-color);font-weight:600;">Error: ${err.message}</span>${dlBtn}`;
      }
    }
  }

  async loadFile(file) {
    try {
      this.showProgress(10, 'Reading file...');
      const bundle = await loadReactionBundle(file, (pct, msg) => this.showProgress(pct, msg));
      this.onBundleLoaded(bundle, 0, false);
    } catch (err) {
      console.error(err);
      this.showDropzone();
      if (this.overlayStatus) {
        this.overlayStatus.innerHTML = `<span style="color:var(--danger-color);font-weight:600;">Error: ${err.message}</span><br><br>Please select a valid bundle .zip`;
      }
    }
  }

  onBundleLoaded(bundle, initialFrame = 0, autoplay = false) {
    this.bundle = bundle;
    const frames = (bundle.trajectory && bundle.trajectory.frames) || [];
    this.totalFrames = frames.length;

    this.viewer3d.setBundle(bundle);
    this.viewerPes.setTrajectory(bundle.trajectory);
    this.ui.setBundle(bundle);

    this.goToFrame(initialFrame >= 0 && initialFrame < this.totalFrames ? initialFrame : 0);
    this.hideOverlay();

    if (autoplay) {
      this.play();
    }
  }

  goToFrame(frameIdx) {
    if (!this.bundle || this.totalFrames === 0) return;
    this.currentFrame = Math.max(0, Math.min(frameIdx, this.totalFrames - 1));

    this.viewer3d.renderFrame(this.currentFrame, false);
    this.viewerPes.setCurrentFrame(this.currentFrame);

    const frame = this.bundle.trajectory.frames[this.currentFrame] || {};
    const e = frame.rel_energy_kcal !== undefined ? frame.rel_energy_kcal : (frame.energy_kcal || 0.0);
    this.ui.updateFrameDisplay(this.currentFrame, this.totalFrames, e);
  }

  stepFrame(delta) {
    let next = this.currentFrame + delta;
    if (next >= this.totalFrames) {
      next = this.loop ? 0 : this.totalFrames - 1;
    } else if (next < 0) {
      next = this.loop ? this.totalFrames - 1 : 0;
    }
    this.goToFrame(next);
  }

  togglePlay() {
    if (this.isPlaying) this.pause();
    else this.play();
  }

  play() {
    if (this.isPlaying || this.totalFrames <= 1) return;
    this.isPlaying = true;
    this.ui.updatePlayState(true);

    const intervalMs = Math.max(25, Math.floor(1000 / (this.baseFps * this.playSpeed)));
    this.playTimer = setInterval(() => {
      let next = this.currentFrame + 1;
      if (next >= this.totalFrames) {
        if (this.loop) {
          next = 0;
        } else {
          this.pause();
          return;
        }
      }
      this.goToFrame(next);
    }, intervalMs);
  }

  pause() {
    this.isPlaying = false;
    if (this.playTimer) {
      clearInterval(this.playTimer);
      this.playTimer = null;
    }
    this.ui.updatePlayState(false);
  }

  setSpeed(speed) {
    this.playSpeed = speed;
    if (this.isPlaying) {
      this.pause();
      this.play();
    }
  }

  bindDropEvents() {
    const box = document.getElementById('drop-box');
    const container = document.getElementById('app-container');

    ['dragenter', 'dragover', 'dragleave', 'drop'].forEach(name => {
      container.addEventListener(name, (e) => e.preventDefault());
    });

    ['dragenter', 'dragover'].forEach(name => {
      container.addEventListener(name, () => box?.classList.add('dragover'));
    });

    ['dragleave', 'drop'].forEach(name => {
      container.addEventListener(name, () => box?.classList.remove('dragover'));
    });

    container.addEventListener('drop', (e) => {
      const files = e.dataTransfer?.files;
      if (files && files.length > 0) {
        this.loadFile(files[0]);
      }
    });

    if (this.fileInput) {
      this.fileInput.addEventListener('change', (e) => {
        if (e.target.files && e.target.files.length > 0) {
          this.loadFile(e.target.files[0]);
        }
      });
    }

    document.getElementById('btn-browse')?.addEventListener('click', () => {
      this.fileInput?.click();
    });
  }

  bindKeyboard() {
    window.addEventListener('keydown', (e) => {
      if (['input', 'textarea', 'select'].includes(document.activeElement?.tagName?.toLowerCase())) return;
      if (e.code === 'Space') {
        e.preventDefault();
        this.togglePlay();
      } else if (e.code === 'ArrowRight') {
        e.preventDefault();
        this.stepFrame(1);
      } else if (e.code === 'ArrowLeft') {
        e.preventDefault();
        this.stepFrame(-1);
      } else if (e.code === 'Home') {
        e.preventDefault();
        this.goToFrame(0);
      } else if (e.code === 'End') {
        e.preventDefault();
        this.goToFrame(this.totalFrames - 1);
      }
    });
  }
}

// Instantiate and start app on page load
window.addEventListener('DOMContentLoaded', () => {
  const app = new ReactionViewerApp();
  app.init();
});
