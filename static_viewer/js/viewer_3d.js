/**
 * viewer_3d.js
 * 3Dmol.js molecular trajectory and surface mesh renderer.
 * Operates entirely from pre-unpacked in-memory meshes with 0 ms extraction latency.
 * Chemistry bond perception and surface rendering delegated to viewer_bonds.js and viewer_surface.js.
 */

import { buildMolBlock } from './viewer_bonds.js';
import {
  renderActiveOrbitals,
  renderEspLayer,
  renderAtomCharges,
  findOrbitalMesh,
  updateEspColorbar,
} from './viewer_surface.js';

export class Viewer3D {
  constructor(containerElement) {
    this.container = containerElement;
    this.glViewer = null;
    this.bundle = null;
    this.currentFrameIdx = 0;
    this.activeOrbitals = new Set();
    this.activeEsp = false;
    this.showCharges = false;
    this.orbitalIsovalue = 0.035;

    this.initGL();
  }

  initGL() {
    if (!window.$3Dmol) {
      console.error('$3Dmol library is not loaded');
      return;
    }
    this.glViewer = window.$3Dmol.createViewer(this.container, {
      backgroundColor: '0x0f141c',
      id: 'reactionViewer'
    });
    window.addEventListener('resize', () => {
      if (this.glViewer) this.glViewer.resize();
    });
  }

  setBundle(bundle) {
    this.bundle = bundle;
    this.currentFrameIdx = 0;
    this.activeOrbitals.clear();
    this.activeEsp = false;
    this.showCharges = false;
    this.updateColorbar(false);
    if (this.glViewer) {
      this.glViewer.clear();
      this.renderFrame(0, true);
    }
  }

  renderFrame(frameIdx, autoZoom = false) {
    if (!this.glViewer || !this.bundle || !this.bundle.trajectory) return;
    const traj = this.bundle.trajectory;
    const frames = traj.frames || [];
    if (frameIdx < 0 || frameIdx >= frames.length) return;

    this.currentFrameIdx = frameIdx;
    const frame = frames[frameIdx];
    const coords = frame.coords || frame.positions || [];
    const symbols = frame.symbols || traj.symbols || [];
    const natoms = Math.min(symbols.length, coords.length);

    // 1. Build authentic SDF format with double/triple bonds
    const sdf = buildMolBlock(symbols, coords);

    // 2. Render molecular structure
    this.glViewer.removeAllSurfaces();
    this.glViewer.removeAllShapes();
    this.glViewer.clear();

    const m = this.glViewer.addModel(sdf, 'sdf');
    m.setStyle({}, {
      stick: { radius: 0.13, colorscheme: 'Jmol' },
      sphere: { scale: 0.28, colorscheme: 'Jmol' }
    });

    // 3. Render all active orbitals concurrently
    if (this.activeOrbitals.size > 0) {
      renderActiveOrbitals(this.glViewer, this.bundle, this.activeOrbitals, frameIdx, this.orbitalIsovalue);
    }

    // 4. Render ESP Surface if active
    if (this.activeEsp) {
      renderEspLayer(this.glViewer, this.bundle, frame, symbols, coords, natoms, frameIdx);
      this.updateColorbar(true);
    }

    // 5. Render Atom Charge labels if showCharges is active
    if (this.showCharges) {
      renderAtomCharges(this.glViewer, this.bundle, frame, symbols, coords, natoms, frameIdx);
    }

    if (autoZoom) {
      this.glViewer.zoomTo();
    }
    this.glViewer.render();
  }

  findOrbitalMesh(frameIdx, orbKey) {
    return findOrbitalMesh(this.bundle, frameIdx, orbKey);
  }

  toggleOrbital(orbKey, isChecked) {
    if (isChecked) this.activeOrbitals.add(orbKey);
    else this.activeOrbitals.delete(orbKey);
    this.renderFrame(this.currentFrameIdx, false);
  }

  clearOrbitals() {
    this.activeOrbitals.clear();
    this.renderFrame(this.currentFrameIdx, false);
  }

  setEspActive(active) {
    this.activeEsp = !!active;
    this.updateColorbar(this.activeEsp);
    this.renderFrame(this.currentFrameIdx, false);
  }

  setChargesActive(active) {
    this.showCharges = !!active;
    this.renderFrame(this.currentFrameIdx, false);
  }

  updateColorbar(active) {
    updateEspColorbar(this.container, active, this.bundle, this.currentFrameIdx);
  }
}
