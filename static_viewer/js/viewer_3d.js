/**
 * viewer_3d.js
 * 3Dmol.js molecular trajectory and surface mesh renderer.
 * Operates entirely from pre-unpacked in-memory meshes with 0 ms extraction latency.
 */

const COV_RADII = {
  H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
  P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39
};
const MAX_VALENCE = { H: 1, F: 1, Cl: 1, Br: 1, I: 1, O: 2, N: 4, C: 4, S: 6, P: 5 };

function computeBonds(symbols, coords) {
  const bonds = [];
  const n = Math.min(symbols.length, coords.length);
  const valences = new Array(n).fill(0);
  const candidates = [];

  for (let i = 0; i < n; i++) {
    const s1 = symbols[i];
    const r1 = COV_RADII[s1] || 0.75;
    for (let j = i + 1; j < n; j++) {
      const s2 = symbols[j];
      const r2 = COV_RADII[s2] || 0.75;
      if (s1 === 'H' && s2 === 'H' && n > 2) continue;

      const dx = coords[i][0] - coords[j][0];
      const dy = coords[i][1] - coords[j][1];
      const dz = coords[i][2] - coords[j][2];
      const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
      const maxD = 1.28 * (r1 + r2);

      if (d >= 0.4 && d <= maxD) {
        const pair = [s1, s2].sort().join('-');
        let targetOrder = 1;
        if (pair === 'C-O') {
          if (d < 1.16) targetOrder = 3;
          else if (d < 1.265) targetOrder = 2;
        } else if (pair === 'C-C') {
          if (d < 1.24) targetOrder = 3;
          else if (d < 1.39) targetOrder = 2;
        } else if (pair === 'C-N') {
          if (d < 1.18) targetOrder = 3;
          else if (d < 1.32) targetOrder = 2;
        } else if (pair === 'N-N') {
          if (d < 1.14) targetOrder = 3;
          else if (d < 1.28) targetOrder = 2;
        } else if (pair === 'N-O' || pair === 'O-O') {
          if (d < 1.25) targetOrder = 2;
        }
        candidates.push({ d, i, j, s1, s2, targetOrder });
      }
    }
  }

  // Shortest bonds assigned first
  candidates.sort((a, b) => a.d - b.d);

  for (const c of candidates) {
    bonds.push({ a1: c.i, a2: c.j, order: 1, target: c.targetOrder, s1: c.s1, s2: c.s2, d: c.d });
    valences[c.i]++;
    valences[c.j]++;
  }

  // Promote to double/triple bond if atom valences allow
  for (const b of bonds) {
    if (b.target > 1) {
      const v1Cap = MAX_VALENCE[b.s1] || 4;
      const v2Cap = MAX_VALENCE[b.s2] || 4;
      const desiredBoost = b.target - 1;
      const maxBoost = Math.min(v1Cap - valences[b.a1], v2Cap - valences[b.a2]);
      const boost = Math.max(0, Math.min(desiredBoost, maxBoost));
      if (boost > 0) {
        b.order += boost;
        valences[b.a1] += boost;
        valences[b.a2] += boost;
      }
    }
  }

  return bonds;
}

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

    // 1. Build SDF format with dynamic bond orders (single, double, triple)
    const natoms = Math.min(symbols.length, coords.length);
    const bonds = computeBonds(symbols, coords);
    let sdf = `\n  Reaction Pathway Viewer\n\n${String(natoms).padStart(3, ' ')}${String(bonds.length).padStart(3, ' ')}  0  0  0  0  0  0  0  0999 V2000\n`;
    for (let i = 0; i < natoms; i++) {
      const sym = symbols[i];
      const pos = coords[i];
      const x = Number(pos[0]).toFixed(4).padStart(10, ' ');
      const y = Number(pos[1]).toFixed(4).padStart(10, ' ');
      const z = Number(pos[2]).toFixed(4).padStart(10, ' ');
      sdf += `${x}${y}${z} ${sym.padEnd(3, ' ')} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
    }
    for (let i = 0; i < bonds.length; i++) {
      const b = bonds[i];
      sdf += `${String(b.a1 + 1).padStart(3, ' ')}${String(b.a2 + 1).padStart(3, ' ')}${String(b.order).padStart(3, ' ')}  0  0  0  0\n`;
    }
    sdf += 'M  END\n$$$$\n';

    // 2. Render molecular structure
    this.glViewer.removeAllSurfaces();
    this.glViewer.removeAllShapes();
    this.glViewer.clear();

    const m = this.glViewer.addModel(sdf, 'sdf');
    m.setStyle({}, {
      stick: { radius: 0.13, colorscheme: 'Jmol' },
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
