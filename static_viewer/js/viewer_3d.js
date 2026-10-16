/**
 * viewer_3d.js
 * 3Dmol.js molecular trajectory and surface mesh renderer.
 * Operates entirely from pre-unpacked in-memory meshes with 0 ms extraction latency.
 */

export class Viewer3D {
  constructor(containerElement) {
    this.container = containerElement;
    this.glViewer = null;
    this.bundle = null;
    this.currentFrameIdx = 0;
    this.activeOrbital = null;
    this.activeEsp = false;
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
    const symbols = traj.symbols || [];

    // 1. Build SDF format in memory
    const natoms = Math.min(symbols.length, coords.length);
    let sdf = `\n  Reaction Pathway Viewer\n\n${String(natoms).padStart(3, ' ')}  0  0  0  0  0  0  0  0  0999 V2000\n`;
    for (let i = 0; i < natoms; i++) {
      const sym = symbols[i];
      const pos = coords[i];
      const x = Number(pos[0]).toFixed(4).padStart(10, ' ');
      const y = Number(pos[1]).toFixed(4).padStart(10, ' ');
      const z = Number(pos[2]).toFixed(4).padStart(10, ' ');
      sdf += `${x}${y}${z} ${sym.padEnd(3, ' ')} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
    }
    sdf += 'M  END\n$$$$\n';

    // 2. Render molecular structure
    this.glViewer.removeAllSurfaces();
    this.glViewer.removeAllShapes();
    this.glViewer.clear();

    const m = this.glViewer.addModel(sdf, 'sdf');
    m.setStyle({}, {
      stick: { radius: 0.15, colorscheme: 'Jmol' },
      sphere: { scale: 0.28, colorscheme: 'Jmol' }
    });

    // 3. Render Orbital Mesh if active
    if (this.activeOrbital && this.bundle.orbitalMeshes) {
      const meshPair = this.findOrbitalMesh(frameIdx, this.activeOrbital);
      if (meshPair) {
        if (meshPair.pos) {
          this.glViewer.addCustom({
            vertexArr: meshPair.pos.vertexArr,
            normalArr: meshPair.pos.normalArr,
            faceArr: meshPair.pos.faceArr,
            color: 0x38bdf8,
            alpha: 0.85
          });
        }
        if (meshPair.neg) {
          this.glViewer.addCustom({
            vertexArr: meshPair.neg.vertexArr,
            normalArr: meshPair.neg.normalArr,
            faceArr: meshPair.neg.faceArr,
            color: 0xf97316,
            alpha: 0.85
          });
        }
      }
    }

    // 4. Render ESP Surface if active
    if (this.activeEsp && this.bundle.espMeshes) {
      const espMesh = this.bundle.espMeshes[frameIdx] || this.bundle.espMeshes[String(frameIdx)];
      if (espMesh && espMesh.vertexArr && espMesh.faceArr) {
        this.glViewer.addCustom({
          vertexArr: espMesh.vertexArr,
          normalArr: espMesh.normalArr || espMesh.vertexArr,
          faceArr: espMesh.faceArr,
          color: espMesh.color || 0x10b981,
          alpha: 0.8
        });
      }
    }

    if (autoZoom) {
      this.glViewer.zoomTo();
    }
    this.glViewer.render();
  }

  findOrbitalMesh(frameIdx, orbKey) {
    if (!this.bundle || !this.bundle.orbitalMeshes) return null;
    const meshes = this.bundle.orbitalMeshes;

    // Exact key: e.g. "0_IBO_21" or "28_HOMO"
    const k1 = `${frameIdx}_${orbKey}`;
    if (meshes[k1]) return meshes[k1];

    // Fallback search
    for (const k in meshes) {
      if (k.startsWith(`${frameIdx}_`) && (k.includes(orbKey) || orbKey.includes(k.replace(`${frameIdx}_`, '')))) {
        return meshes[k];
      }
    }

    // Direct key without frame prefix (for single molecule bundles)
    if (meshes[orbKey]) return meshes[orbKey];
    return null;
  }

  setOrbital(orbKey) {
    this.activeOrbital = orbKey;
    this.renderFrame(this.currentFrameIdx, false);
  }

  setEspActive(active) {
    this.activeEsp = !!active;
    this.renderFrame(this.currentFrameIdx, false);
  }
}
