/**
 * viewer_app.js
 * Main coordinator orchestrating loader, 3Dmol viewer, static PES plot, and playback loop.
 * Supports multi-bundle dropdown switching, molecule inspection modal, and self-invoking popups.
 */

import { loadReactionBundle, buildMoleculeBundle } from './viewer_loader.js';
import { Viewer3D } from './viewer_3d.js';
import { ViewerPES } from './viewer_pes.js';
import { ViewerUI } from './viewer_ui.js';
import { ViewerMoleculeModal } from './viewer_modal.js';

class ReactionViewerApp {
  constructor() {
    this.bundle = null;
    this.currentBundleUrl = '';
    this.bundleList = [];
    this.currentFrame = 0;
    this.totalFrames = 0;
    this.isPlaying = false;
    this.playDirection = 'forward';
    this.playSpeed = 1.0;
    this.playTimer = null;
    this.baseFps = 15;

    // Elements
    this.overlay = document.getElementById('overlay');
    this.overlayStatus = document.getElementById('overlay-status');
    this.progressBarContainer = document.getElementById('progress-bar-container');
    this.progressBar = document.getElementById('progress-bar');
    this.fileInput = document.getElementById('file-input');
    this.bundleSelector = document.getElementById('sel-bundle-switcher');

    // Submodules
    this.viewer3d = new Viewer3D(document.getElementById('gl-canvas-container'));
    this.viewerPes = new ViewerPES(document.getElementById('pes-canvas'));
    this.ui = new ViewerUI(this);
    this.moleculeModal = new ViewerMoleculeModal(this);

    this.bindDropEvents();
    this.bindKeyboard();
    this.bindBundleSelector();
  }

  init() {
    const params = new URLSearchParams(window.location.search);
    const bundleParams = [];

    // Parse ?bundle=... and ?bundles=...
    params.getAll('bundle').forEach(b => {
      b.split(',').forEach(item => {
        const tr = item.trim();
        if (tr && !bundleParams.includes(tr)) bundleParams.push(tr);
      });
    });
    if (params.get('bundles')) {
      params.get('bundles').split(',').forEach(item => {
        const tr = item.trim();
        if (tr && !bundleParams.includes(tr)) bundleParams.push(tr);
      });
    }

    this.bundleList = bundleParams;
    this.populateBundleSelector(this.bundleList);

    const initialFrame = parseInt(params.get('frame') || '0', 10);
    const autoplay = params.get('autoplay') === '1' || params.get('autoplay') === 'true';
    const targetMolecule = params.get('molecule');

    if (this.bundleList.length > 0) {
      this.loadUrl(this.bundleList[0], initialFrame, autoplay, targetMolecule);
    } else {
      this.showDropzone();
    }
  }

  bindBundleSelector() {
    if (!this.bundleSelector) return;
    this.bundleSelector.addEventListener('change', (e) => {
      const val = e.target.value;
      if (val && val !== this.currentBundleUrl) {
        this.loadUrl(val, 0, false);
      }
    });
  }

  populateBundleSelector(bundleUrls) {
    if (!this.bundleSelector) return;
    if (!bundleUrls || bundleUrls.length <= 1) {
      this.bundleSelector.style.display = 'none';
      return;
    }

    this.bundleSelector.style.display = 'inline-block';
    this.bundleSelector.innerHTML = '';

    bundleUrls.forEach((url) => {
      const opt = document.createElement('option');
      opt.value = url;
      const fname = url.split('/').pop().replace(/\.zip$/i, '').replace(/_/g, ' ');
      opt.textContent = `📦 ${fname}`;
      if (url === this.currentBundleUrl) opt.selected = true;
      this.bundleSelector.appendChild(opt);
    });
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

  async loadUrl(url, initialFrame = 0, autoplay = false, targetMolecule = null) {
    try {
      this.currentBundleUrl = url;
      this.showProgress(10, 'Fetching bundle...');
      const bundle = await loadReactionBundle(url, (pct, msg) => this.showProgress(pct, msg));
      this.onBundleLoaded(bundle, initialFrame, autoplay, targetMolecule);
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
      this.currentBundleUrl = file.name;
      this.onBundleLoaded(bundle, 0, false, null);
    } catch (err) {
      console.error(err);
      this.showDropzone();
      if (this.overlayStatus) {
        this.overlayStatus.innerHTML = `<span style="color:var(--danger-color);font-weight:600;">Error: ${err.message}</span><br><br>Please select a valid bundle .zip`;
      }
    }
  }

  onBundleLoaded(bundle, initialFrame = 0, autoplay = false, targetMolecule = null) {
    if (targetMolecule) {
      const molBundle = buildMoleculeBundle(targetMolecule, bundle);
      if (molBundle) {
        this.loadMoleculeAsMain(molBundle);
        return;
      }
    }

    this.bundle = bundle;
    const frames = (bundle.trajectory && bundle.trajectory.frames) || [];
    this.totalFrames = frames.length;

    this.viewer3d.setBundle(bundle);
    this.viewerPes.setTrajectory(bundle.trajectory);
    this.ui.setBundle(bundle);

    this.populateBundleSelector(this.bundleList);

    this.goToFrame(initialFrame >= 0 && initialFrame < this.totalFrames ? initialFrame : 0);
    this.hideOverlay();

    if (autoplay) {
      this.playForward();
    }
  }

  loadMoleculeAsMain(molBundle) {
    this.bundle = molBundle;
    this.totalFrames = 1;
    this.currentFrame = 0;
    this.viewer3d.setBundle(molBundle);
    this.ui.setMoleculeMode(molBundle);
    this.hideOverlay();
  }

  openMoleculeModal(molId, cfg) {
    this.moleculeModal.show(molId, cfg);
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

  setDirection(dir) {
    this.playDirection = dir;
  }

  stepFrame(delta) {
    this.pause();
    this.setDirection(delta >= 0 ? 'forward' : 'reverse');
    const next = Math.max(0, Math.min(this.totalFrames - 1, this.currentFrame + delta));
    this.goToFrame(next);
  }

  togglePlay() {
    this.togglePlayForward();
  }

  togglePlayForward() {
    if (this.isPlaying) {
      this.pause();
    } else {
      this.playForward();
    }
  }

  togglePlayReverse() {
    if (this.isPlaying) {
      this.pause();
    } else {
      this.playReverse();
    }
  }

  playForward() {
    this.pause();
    this.setDirection('forward');
    if (this.totalFrames <= 1) return;
    if (this.currentFrame >= this.totalFrames - 1) this.goToFrame(0);

    this.isPlaying = true;
    this.ui.updatePlayState(true, 'forward');

    const intervalMs = Math.max(25, Math.floor(1000 / (this.baseFps * this.playSpeed)));
    this.playTimer = setInterval(() => {
      if (this.currentFrame >= this.totalFrames - 1) {
        this.pause();
      } else {
        this.goToFrame(this.currentFrame + 1);
      }
    }, intervalMs);
  }

  playReverse() {
    this.pause();
    this.setDirection('reverse');
    if (this.totalFrames <= 1) return;
    if (this.currentFrame <= 0) this.goToFrame(this.totalFrames - 1);

    this.isPlaying = true;
    this.ui.updatePlayState(true, 'reverse');

    const intervalMs = Math.max(25, Math.floor(1000 / (this.baseFps * this.playSpeed)));
    this.playTimer = setInterval(() => {
      if (this.currentFrame <= 0) {
        this.pause();
      } else {
        this.goToFrame(this.currentFrame - 1);
      }
    }, intervalMs);
  }

  playSegment(startFrame, endFrame, direction = 'forward') {
    this.pause();
    if (this.totalFrames <= 1) return;

    this.setDirection(direction);

    const f1 = Math.max(0, Math.min(this.totalFrames - 1, parseInt(startFrame, 10) || 0));
    const f2 = Math.max(0, Math.min(this.totalFrames - 1, parseInt(endFrame, 10) || 0));
    const minF = Math.min(f1, f2);
    const maxF = Math.max(f1, f2);

    const fromF = (direction === 'reverse') ? maxF : minF;
    const toF = (direction === 'reverse') ? minF : maxF;

    this.goToFrame(fromF);
    this.isPlaying = true;
    this.ui.updatePlayState(true, direction);

    const intervalMs = Math.max(25, Math.floor(1000 / (this.baseFps * this.playSpeed)));
    this.playTimer = setInterval(() => {
      if (direction === 'forward') {
        if (this.currentFrame >= toF) {
          this.pause();
        } else {
          this.goToFrame(this.currentFrame + 1);
        }
      } else {
        if (this.currentFrame <= toF) {
          this.pause();
        } else {
          this.goToFrame(this.currentFrame - 1);
        }
      }
    }, intervalMs);
  }

  pause() {
    this.isPlaying = false;
    if (this.playTimer) {
      clearInterval(this.playTimer);
      this.playTimer = null;
    }
    this.ui.updatePlayState(false, this.playDirection);
  }

  setSpeed(speed) {
    this.playSpeed = speed;
    if (this.isPlaying) {
      const dir = this.playDirection;
      this.pause();
      if (dir === 'reverse') this.playReverse();
      else this.playForward();
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

// Start app on DOMContentLoaded
window.addEventListener('DOMContentLoaded', () => {
  const app = new ReactionViewerApp();
  window.app = app;
  app.init();
});
