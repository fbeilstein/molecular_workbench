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

    // 2. Mechanism SVG & Wiring
    this.wiring = bundle.svgWiring || (bundle.manifest && bundle.manifest.svg_wiring) || (bundle.trajectory && bundle.trajectory.svg_wiring) || {};
    if (this.elMechanism && bundle.mechanismSvg) {
      this.elMechanism.innerHTML = bundle.mechanismSvg;
      this.initInteractiveSvg();
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

  initInteractiveSvg() {
    if (!this.elMechanism || !this.wiring) return;
    const svg = this.elMechanism.querySelector('svg');
    if (!svg) return;

    for (const [key, cfg] of Object.entries(this.wiring)) {
      const el = svg.getElementById(key) ||
                 svg.querySelector(`[inkscape\\:label="${key}"]`) ||
                 svg.querySelector(`[label="${key}"]`);
      if (!el) continue;

      if (cfg.role === 'arrow' || cfg.role === 'push_arrow') {
        el.classList.add('svg-clickable-arrow');
        const dir = cfg.direction === 'reverse' ? '⟵' : '⟶';
        el.title = `Reaction Step: ${cfg.startFrame ?? 0} ${dir} ${cfg.endFrame ?? 0} (Click to play)`;
      } else if (cfg.role === 'state' || cfg.role === 'state_molecule') {
        el.classList.add('svg-clickable-state');
        el.title = `State (Frame ${cfg.targetFrame ?? 0}): Click to jump`;
      } else if (cfg.role === 'molecule') {
        el.classList.add('svg-clickable-molecule');
        el.title = `${cfg.moleculeName || 'Molecule'}`;
      }
    }

    svg.addEventListener('click', (e) => {
      let cur = e.target;
      while (cur && cur !== svg && cur.tagName && cur.tagName.toLowerCase() !== 'svg') {
        const id = cur.getAttribute('id');
        const lbl = cur.getAttribute('inkscape:label') || cur.getAttribute('label');
        const cfg = this.wiring[id] || this.wiring[lbl];
        if (cfg) {
          e.stopPropagation();
          if (cfg.role === 'arrow' || cfg.role === 'push_arrow') {
            const sf = cfg.startFrame !== undefined ? parseInt(cfg.startFrame, 10) : 0;
            const ef = cfg.endFrame !== undefined ? parseInt(cfg.endFrame, 10) : this.app.totalFrames - 1;
            this.app.playSegment(sf, ef, cfg.direction || 'forward');
          } else if (cfg.role === 'state' || cfg.role === 'state_molecule') {
            const tf = cfg.targetFrame !== undefined ? parseInt(cfg.targetFrame, 10) : 0;
            this.app.goToFrame(tf);
          }
          return;
        }
        cur = cur.parentElement;
      }
    });
  }

  updateSvgActive(frameIdx) {
    if (!this.elMechanism || !this.wiring) return;
    const svg = this.elMechanism.querySelector('svg');
    if (!svg) return;

    svg.querySelectorAll('.svg-active-element').forEach(el => el.classList.remove('svg-active-element'));

    let closestState = null;
    let minDiff = 12;
    for (const [key, cfg] of Object.entries(this.wiring)) {
      if ((cfg.role === 'state' || cfg.role === 'state_molecule') && cfg.targetFrame !== undefined) {
        const diff = Math.abs(cfg.targetFrame - frameIdx);
        if (diff < minDiff) {
          minDiff = diff;
          closestState = key;
        }
      }
    }
    if (closestState) {
      const el = svg.getElementById(closestState) || svg.querySelector(`[inkscape\\:label="${closestState}"]`);
      if (el) el.classList.add('svg-active-element');
    }
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
    this.updateSvgActive(frameIdx);
  }

  updatePlayState(isPlaying) {
    if (this.elBtnPlay) {
      this.elBtnPlay.textContent = isPlaying ? '❚❚ Pause' : '▶ Play';
    }
  }
}
