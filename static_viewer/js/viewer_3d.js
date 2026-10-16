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
    this.activeOrbitals = new Set();
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

    // 3. Render all Active Orbitals concurrently
    if (this.activeOrbitals.size > 0 && this.bundle.orbitalMeshes) {
      const CC = window.$3Dmol ? window.$3Dmol.CC : null;
      for (const orbKey of this.activeOrbitals) {
        const meshPair = this.findOrbitalMesh(frameIdx, orbKey);
        if (meshPair) {
          const pV = meshPair.pos?.vertices || meshPair.pos?.vertexArr;
          const pF = meshPair.pos?.faces || meshPair.pos?.faceArr;
          const pN = meshPair.pos?.normals || meshPair.pos?.normalArr;
          if (pV && pF) {
            this.glViewer.addCustom({
              vertexArr: pV,
              faceArr: pF,
              normalArr: pN,
              color: CC ? CC.color('#38bdf8') : 0x38bdf8,
              opacity: 0.75
            });
          }
          const nV = meshPair.neg?.vertices || meshPair.neg?.vertexArr;
          const nF = meshPair.neg?.faces || meshPair.neg?.faceArr;
          const nN = meshPair.neg?.normals || meshPair.neg?.normalArr;
          if (nV && nF) {
            this.glViewer.addCustom({
              vertexArr: nV,
              faceArr: nF,
              normalArr: nN,
              color: CC ? CC.color('#f97316') : 0xf97316,
              opacity: 0.75
            });
          }
        } else if (this.bundle.cubeCache) {
          const cubeData = this.bundle.cubeCache[`${frameIdx}_${orbKey}`] || this.bundle.cubeCache[orbKey];
          if (cubeData && window.$3Dmol) {
            try {
              const vdo = new window.$3Dmol.VolumeData(cubeData, 'cube');
              const iso = this.orbitalIsovalue || 0.035;
              this.glViewer.addIsosurface(vdo, { isoval: iso, color: '#38bdf8', opacity: 0.72 });
              this.glViewer.addIsosurface(vdo, { isoval: -iso, color: '#f97316', opacity: 0.72 });
            } catch (_) {}
          }
        }
      }
    }

    // 4. Render ESP Surface if active
    if (this.activeEsp) {
      const espMesh = this.bundle.espMeshes ? (
        this.bundle.espMeshes[frameIdx] ||
        this.bundle.espMeshes[String(frameIdx)] ||
        (this.bundle.isSingleMolecule ? (this.bundle.espMeshes['0'] || Object.values(this.bundle.espMeshes)[0]) : null)
      ) : null;

      const eV = espMesh?.vertices || espMesh?.vertexArr;
      const eF = espMesh?.faces || espMesh?.faceArr;
      const eN = espMesh?.normals || espMesh?.normalArr;

      if (eV && eF) {
        let colorParam = 0x10b981;
        if (espMesh.esp_values && espMesh.esp_values.length === eV.length) {
          const absMax = Math.max(Math.abs(espMesh.esp_min || -30), Math.abs(espMesh.esp_max || 30), 20);
          colorParam = espMesh.esp_values.map(v => {
            const t = Math.max(-1, Math.min(1, v / absMax));
            return t < 0 ? { r: 1.0, g: 1.0 + t, b: 1.0 + t } : { r: 1.0 - t, g: 1.0 - t, b: 1.0 };
          });
        }
        const customObj = {
          vertexArr: eV,
          faceArr: eF,
          color: colorParam,
          opacity: 0.82
        };
        if (eN && eN.length === eV.length) {
          customObj.normalArr = eN;
        }
        this.glViewer.addCustom(customObj);
      } else if (window.$3Dmol) {
        // Fallback for reaction frames without precomputed mesh: dynamic molecular surface
        try {
          this.glViewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
            opacity: 0.70,
            colorscheme: 'Jmol'
          });
        } catch (e) {
          console.warn('Dynamic surface rendering failed:', e);
        }
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

    // 1. Exact frame match: e.g. "0_IBO_21"
    const k1 = `${frameIdx}_${orbKey}`;
    if (meshes[k1]) return meshes[k1];

    // 2. Direct key without frame prefix (only for single molecule mode)
    if (this.bundle.isSingleMolecule && meshes[orbKey]) {
      return meshes[orbKey];
    }

    // 3. Search within current frame
    for (const k in meshes) {
      if (k.startsWith(`${frameIdx}_`)) {
        const baseKey = k.replace(`${frameIdx}_`, '');
        if (baseKey === orbKey || k.includes(orbKey) || orbKey.includes(baseKey)) {
          return meshes[k];
        }
      }
    }

    // 4. In single-molecule mode only: allow matching regardless of prefix
    if (this.bundle.isSingleMolecule) {
      for (const k in meshes) {
        const stripped = k.replace(/^\d+_/, '');
        if (stripped === orbKey || k.includes(orbKey) || orbKey.includes(stripped)) {
          return meshes[k];
        }
      }
    }

    // In reaction/trajectory mode: if no orbital mesh was computed for this exact frameIdx,
    // return null. Orbitals are tied to atomic coordinates; showing another frame's mesh
    // leaves a frozen ghost orbital floating in space while the molecule moves away.
    return null;
  }

  toggleOrbital(orbKey, isChecked) {
    if (isChecked) {
      this.activeOrbitals.add(orbKey);
    } else {
      this.activeOrbitals.delete(orbKey);
    }
    this.renderFrame(this.currentFrameIdx, false);
  }

  clearOrbitals() {
    this.activeOrbitals.clear();
    this.renderFrame(this.currentFrameIdx, false);
  }

  setEspActive(active) {
    this.activeEsp = !!active;
    this.renderFrame(this.currentFrameIdx, false);
  }
}
