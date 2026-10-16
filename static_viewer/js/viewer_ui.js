/**
 * viewer_ui.js
 * Static UI manager for mechanism header, side orbital controls, and playback transport dock.
 */

export class ViewerUI {
  constructor(app) {
    this.app = app;

    // DOM references
    this.elTitle = document.getElementById('reaction-title');
    this.elBarrier = document.getElementById('stat-barrier');
    this.elSolvent = document.getElementById('stat-solvent');
    this.elMechanism = document.getElementById('mechanism-strip');
    this.elSlider = document.getElementById('frame-slider');
    this.elFrameInfo = document.getElementById('frame-info');
    this.elBtnPlay = document.getElementById('btn-play');
    this.elOrbitalsList = document.getElementById('orbitals-list');
    this.elEspCheck = document.getElementById('chk-esp');

    this.bindEvents();
  }

  bindEvents() {
    // Play / Pause
    if (this.elBtnPlay) {
      this.elBtnPlay.addEventListener('click', () => this.app.togglePlay());
    }

    // Step Buttons
    document.getElementById('btn-first')?.addEventListener('click', () => this.app.goToFrame(0));
    document.getElementById('btn-last')?.addEventListener('click', () => this.app.goToFrame(this.app.totalFrames - 1));
    document.getElementById('btn-prev')?.addEventListener('click', () => this.app.stepFrame(-1));
    document.getElementById('btn-next')?.addEventListener('click', () => this.app.stepFrame(1));

    // Slider
    if (this.elSlider) {
      this.elSlider.addEventListener('input', (e) => {
        this.app.goToFrame(parseInt(e.target.value, 10));
      });
    }

    // Speed selector
    const speedSelect = document.getElementById('select-speed');
    if (speedSelect) {
      speedSelect.addEventListener('change', (e) => {
        this.app.setSpeed(parseFloat(e.target.value));
      });
    }

    // Loop
    const chkLoop = document.getElementById('chk-loop');
    if (chkLoop) {
      chkLoop.addEventListener('change', (e) => {
        this.app.loop = e.target.checked;
      });
    }

    // Reset Camera
    document.getElementById('btn-reset-cam')?.addEventListener('click', () => {
      if (this.app.viewer3d && this.app.viewer3d.glViewer) {
        this.app.viewer3d.glViewer.zoomTo();
        this.app.viewer3d.glViewer.render();
      }
    });

    // ESP Toggle
    if (this.elEspCheck) {
      this.elEspCheck.addEventListener('change', (e) => {
        if (this.app.viewer3d) {
          this.app.viewer3d.setEspActive(e.target.checked);
        }
      });
    }
  }

  setBundle(bundle) {
    const manifest = bundle.manifest || {};
    const traj = bundle.trajectory || {};

    // 1. Meta / Title
    if (this.elTitle) {
      this.elTitle.textContent = manifest.title || traj.title || 'Multistep Reaction Pathway';
    }
    if (this.elBarrier) {
      const b = manifest.overall_barrier_kcal !== undefined ? manifest.overall_barrier_kcal : traj.overall_barrier_kcal;
      this.elBarrier.textContent = b !== undefined ? `${b} kcal/mol` : '—';
    }
    if (this.elSolvent) {
      this.elSolvent.textContent = manifest.solvent || traj.solvent || 'Gas Phase';
    }

    // 2. Mechanism SVG
    if (this.elMechanism && bundle.mechanismSvg) {
      this.elMechanism.innerHTML = bundle.mechanismSvg;
    }

    // 3. Slider configuration
    const nFrames = (traj.frames || []).length;
    if (this.elSlider) {
      this.elSlider.min = 0;
      this.elSlider.max = Math.max(0, nFrames - 1);
      this.elSlider.value = 0;
    }

    // 4. Populate Orbitals Checklist
    this.populateOrbitals(bundle);
  }

  populateOrbitals(bundle) {
    if (!this.elOrbitalsList) return;
    this.elOrbitalsList.innerHTML = '';

    const meshes = bundle.orbitalMeshes || {};
    const keys = Object.keys(meshes);
    if (!keys.length) {
      this.elOrbitalsList.innerHTML = '<span style="color:var(--text-secondary);font-size:11px;">No orbital surfaces in bundle</span>';
      return;
    }

    // Extract unique orbital names across frames
    const orbNames = new Set();
    for (const k of keys) {
      const parts = k.split('_');
      if (parts.length > 1 && !isNaN(parts[0])) {
        orbNames.add(parts.slice(1).join('_'));
      } else {
        orbNames.add(k);
      }
    }

    // "None" option
    const noneItem = document.createElement('label');
    noneItem.innerHTML = '<input type="radio" name="orbital-choice" value="" checked> None';
    noneItem.querySelector('input').addEventListener('change', () => {
      this.app.viewer3d.setOrbital(null);
    });
    this.elOrbitalsList.appendChild(noneItem);

    for (const name of Array.from(orbNames).sort()) {
      const item = document.createElement('label');
      item.innerHTML = `<input type="radio" name="orbital-choice" value="${name}"> ${name}`;
      item.querySelector('input').addEventListener('change', () => {
        this.app.viewer3d.setOrbital(name);
      });
      this.elOrbitalsList.appendChild(item);
    }
  }

  updateFrameDisplay(frameIdx, totalFrames, energy) {
    if (this.elSlider) {
      this.elSlider.value = frameIdx;
    }
    if (this.elFrameInfo) {
      const eStr = energy !== undefined ? ` | ${energy >= 0 ? '+' : ''}${energy.toFixed(2)} kcal/mol` : '';
      this.elFrameInfo.innerHTML = `Frame <strong>${frameIdx + 1}</strong> / ${totalFrames}${eStr}`;
    }
  }

  updatePlayState(isPlaying) {
    if (this.elBtnPlay) {
      this.elBtnPlay.textContent = isPlaying ? '❚❚ Pause' : '▶ Play';
    }
  }
}
