/**
 * Surface and Isosurface Rendering Engine for Standalone 3D Viewer.
 * Manages volumetric orbitals, dual-phase coloring, and ESP electrostatic surfaces.
 */

export function findOrbitalMesh(bundle, frameIdx, orbKey) {
  if (!bundle || !bundle.orbitalMeshes) return null;
  const meshes = bundle.orbitalMeshes;
  const exactKey = `${frameIdx}_${orbKey}`;
  if (meshes[exactKey]) return meshes[exactKey];
  if (bundle.isSingleMolecule) {
    if (meshes[orbKey]) return meshes[orbKey];
    if (meshes[`0_${orbKey}`]) return meshes[`0_${orbKey}`];
  }
  return null;
}

export function renderActiveOrbitals(glViewer, bundle, activeOrbitals, frameIdx, isovalue) {
  if (!activeOrbitals || activeOrbitals.size === 0 || !bundle) return;
  const CC = window.$3Dmol ? window.$3Dmol.CC : null;

  for (const orbKey of activeOrbitals) {
    const sKey = String(orbKey).toLowerCase();
    const isAnti = sKey.includes('star') || sKey.includes('*') || sKey.includes('lumo') || sKey.includes('antibond') || sKey.includes('virtual');
    const posColorHex = isAnti ? '#fb923c' : '#38bdf8';
    const negColorHex = isAnti ? '#db2777' : '#1d4ed8';
    const posColorInt = isAnti ? 0xfb923c : 0x38bdf8;
    const negColorInt = isAnti ? 0xdb2777 : 0x1d4ed8;

    const meshPair = findOrbitalMesh(bundle, frameIdx, orbKey);
    if (meshPair) {
      const pV = meshPair.pos?.vertices || meshPair.pos?.vertexArr;
      const pF = meshPair.pos?.faces || meshPair.pos?.faceArr;
      const pN = meshPair.pos?.normals || meshPair.pos?.normalArr;
      if (pV && pF) {
        glViewer.addCustom({
          vertexArr: pV, faceArr: pF, normalArr: pN,
          color: CC ? CC.color(posColorHex) : posColorInt, opacity: 0.75
        });
      }
      const nV = meshPair.neg?.vertices || meshPair.neg?.vertexArr;
      const nF = meshPair.neg?.faces || meshPair.neg?.faceArr;
      const nN = meshPair.neg?.normals || meshPair.neg?.normalArr;
      if (nV && nF) {
        glViewer.addCustom({
          vertexArr: nV, faceArr: nF, normalArr: nN,
          color: CC ? CC.color(negColorHex) : negColorInt, opacity: 0.75
        });
      }
    } else if (bundle.cubeCache) {
      const cubeData = bundle.cubeCache[`${frameIdx}_${orbKey}`] || bundle.cubeCache[orbKey];
      if (cubeData && window.$3Dmol) {
        try {
          const vdo = new window.$3Dmol.VolumeData(cubeData, 'cube');
          const iso = isovalue || 0.035;
          glViewer.addIsosurface(vdo, { isoval: iso, color: posColorHex, opacity: 0.72 });
          glViewer.addIsosurface(vdo, { isoval: -iso, color: negColorHex, opacity: 0.72 });
        } catch (_) {}
      }
    }
  }
}

export function renderEspLayer(glViewer, bundle, frame, symbols, coords, natoms, frameIdx) {
  let espMesh = bundle.espMeshes ? (
    bundle.espMeshes[frameIdx] ||
    bundle.espMeshes[String(frameIdx)] ||
    (bundle.isSingleMolecule ? (bundle.espMeshes['0'] || Object.values(bundle.espMeshes)[0]) : null)
  ) : null;

  // Fallback during initial prewarming: display closest cached frame so surface never disappears
  if (!espMesh && bundle.espMeshes) {
    const keys = Object.keys(bundle.espMeshes).map(k => parseInt(k, 10)).filter(n => !isNaN(n));
    if (keys.length > 0) {
      let closest = keys[0];
      for (const k of keys) {
        if (Math.abs(k - frameIdx) < Math.abs(closest - frameIdx)) closest = k;
      }
      espMesh = bundle.espMeshes[closest] || bundle.espMeshes[String(closest)];
    }
  }

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
    } else if (espMesh.color && Array.isArray(espMesh.color)) {
      colorParam = espMesh.color;
    }
    const customObj = { vertexArr: eV, faceArr: eF, color: colorParam, opacity: 0.65 };
    if (eN && eN.length === eV.length) customObj.normalArr = eN;
    glViewer.addCustom(customObj);
  } else if (window.$3Dmol) {
    try {
      const traj = bundle.trajectory || {};
      const trajCharges = traj.trajectory_charges ? traj.trajectory_charges[frameIdx] : null;
      const frameCharges = frame.charges || trajCharges || bundle.charges || [];
      const atoms = glViewer.selectedAtoms({});
      for (let i = 0; i < atoms.length && i < natoms; i++) {
        let q = (frameCharges && frameCharges[i] !== undefined) ? frameCharges[i] : null;
        if (q === null) {
          const s = symbols[i];
          q = s === 'H' ? 0.08 : (s === 'O' ? -0.50 : (s === 'Cl' ? -0.45 : 0.0));
        }
        atoms[i].partialCharge = q;
        atoms[i].properties = Object.assign(atoms[i].properties || {}, { partialCharge: q });
      }
      const p = glViewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
        opacity: 0.65,
        colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
      });
      if (p && typeof p.then === 'function') {
        p.then(sId => {
          const surf = glViewer.surfaces?.[sId];
          const meshObj = surf?.[0] || surf;
          const grp = meshObj?.geo?.geometryGroups?.[0];
          if (grp && grp.vertexArray && grp.faceArray) {
            const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
            const fCount = grp.faceidx || grp.faceArray.length || 0;
            const verts = [], faces = Array.from(grp.faceArray.slice(0, fCount)), cols = [];
            for (let vi = 0; vi < vCount; vi++) {
              verts.push({ x: grp.vertexArray[vi * 3], y: grp.vertexArray[vi * 3 + 1], z: grp.vertexArray[vi * 3 + 2] });
            }
            if (grp.colorArray && grp.colorArray.length >= vCount * 3) {
              for (let ci = 0; ci < vCount; ci++) {
                cols.push({ r: grp.colorArray[ci * 3], g: grp.colorArray[ci * 3 + 1], b: grp.colorArray[ci * 3 + 2] });
              }
            }
            bundle.espMeshes = bundle.espMeshes || {};
            bundle.espMeshes[frameIdx] = { vertexArr: verts, faceArr: faces, color: cols.length ? cols : undefined };
          }
        }).catch(() => {});
      }
    } catch (e) {
      console.warn('Dynamic surface rendering failed:', e);
    }
  }
}

export function renderAtomCharges(glViewer, bundle, frame, symbols, coords, natoms, frameIdx) {
  const traj = bundle.trajectory || {};
  const trajCharges = traj.trajectory_charges ? traj.trajectory_charges[frameIdx] : null;
  const charges = trajCharges || frame.charges || traj.charges || bundle.charges || bundle.manifest?.charges || [];
  const labels = frame.atom_labels || (traj.trajectory_atom_labels ? traj.trajectory_atom_labels[frameIdx] : null) || null;
  for (let i = 0; i < natoms; i++) {
    let q = (charges && charges[i] !== undefined) ? charges[i] : null;
    if (q === null) {
      const s = symbols[i];
      q = s === 'H' ? 0.06 : (s === 'O' ? -0.50 : (s === 'N' ? -0.40 : 0.0));
    }
    const sign = q > 0 ? '+' : '';
    const tag = (labels && labels[i]) ? labels[i] : `${symbols[i]}`;
    glViewer.addLabel(`${tag}: ${sign}${Number(q).toFixed(2)}`, {
      position: { x: coords[i][0], y: coords[i][1] + 0.35, z: coords[i][2] },
      fontSize: 10, fontColor: '#ffffff',
      backgroundColor: 'rgba(15, 23, 42, 0.75)', backgroundOpacity: 0.75
    });
  }
}

export function updateEspColorbar(container, active, bundle, frameIdx) {
  const colorbar = container.parentElement?.querySelector('#esp-colorbar') || document.getElementById('esp-colorbar');
  if (!colorbar) return;
  colorbar.style.display = active ? 'flex' : 'none';
  if (active && bundle) {
    const espMesh = bundle.espMeshes ? (
      bundle.espMeshes[frameIdx] ||
      bundle.espMeshes[String(frameIdx)] ||
      (bundle.isSingleMolecule ? (bundle.espMeshes['0'] || Object.values(bundle.espMeshes)[0]) : null)
    ) : null;
    const absMax = Math.max(Math.abs(espMesh?.esp_min || -30), Math.abs(espMesh?.esp_max || 30), 20);
    const minEl = colorbar.querySelector('#esp-colorbar-min') || document.getElementById('esp-colorbar-min');
    const maxEl = colorbar.querySelector('#esp-colorbar-max') || document.getElementById('esp-colorbar-max');
    if (minEl) minEl.textContent = `${(-absMax).toFixed(1)}`;
    if (maxEl) maxEl.textContent = `+${absMax.toFixed(1)}`;
  }
}

let _bgViewerDiv = null;
let _bgViewerInstance = null;

export async function prewarmViewerEspCache(bundle, onFramePrewarmed) {
  if (!bundle || !bundle.trajectory || !window.$3Dmol) return;
  const traj = bundle.trajectory;
  const frames = traj.frames || [];
  const total = frames.length;
  if (total === 0) return;

  const trajectoryCharges = traj.trajectory_charges || [];
  const symbols = traj.symbols || [];
  bundle.espMeshes = bundle.espMeshes || {};

  try {
    if (!_bgViewerDiv) {
      _bgViewerDiv = document.createElement('div');
      _bgViewerDiv.style.cssText = 'position:fixed;left:-9999px;width:120px;height:120px;opacity:0;pointer-events:none;';
      document.body.appendChild(_bgViewerDiv);
      _bgViewerInstance = window.$3Dmol.createViewer(_bgViewerDiv);
    }
    const bg = _bgViewerInstance;

    for (let f = 0; f < total; f++) {
      if (bundle.espMeshes[f] || bundle.espMeshes[String(f)]) continue;
      const coords = frames[f]?.coords || frames[f]?.positions;
      const charges = trajectoryCharges?.[f] || frames[f]?.charges;
      const syms = Array.isArray(symbols?.[f]) ? symbols[f] : (frames[f]?.symbols || symbols);
      if (!coords || !syms) continue;
      const n = Math.min(syms.length, coords.length);

      bg.clear();
      let mdl = `${n}\nFrame Model\n`;
      for (let i = 0; i < n; i++) {
        const [x, y, z] = coords[i];
        mdl += `${syms[i]} ${x.toFixed(4)} ${y.toFixed(4)} ${z.toFixed(4)}\n`;
      }
      bg.addModel(mdl, 'xyz');
      const atoms = bg.selectedAtoms({});
      for (let i = 0; i < atoms.length && i < n; i++) {
        let qVal = (charges && charges[i] !== undefined) ? charges[i] : null;
        if (qVal === null) {
          const s = syms[i];
          qVal = s === 'H' ? 0.08 : (s === 'O' ? -0.50 : (s === 'Cl' ? -0.45 : 0.0));
        }
        atoms[i].partialCharge = qVal;
        atoms[i].properties = { partialCharge: qVal };
      }
      const p = bg.addSurface(window.$3Dmol.SurfaceType.VDW, {
        opacity: 0.65,
        colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
      });
      if (p && typeof p.then === 'function') {
        const sId = await p;
        const surf = bg.surfaces?.[sId];
        const meshObj = surf?.[0] || surf;
        const grp = meshObj?.geo?.geometryGroups?.[0];
        if (grp && grp.vertexArray && grp.faceArray) {
          const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
          const fCount = grp.faceidx || grp.faceArray.length || 0;
          const verts = [];
          for (let vi = 0; vi < vCount; vi++) {
            verts.push({ x: grp.vertexArray[vi * 3], y: grp.vertexArray[vi * 3 + 1], z: grp.vertexArray[vi * 3 + 2] });
          }
          const faces = Array.from(grp.faceArray.slice(0, fCount));
          const cols = [];
          if (grp.colorArray && grp.colorArray.length >= vCount * 3) {
            for (let ci = 0; ci < vCount; ci++) {
              cols.push({ r: grp.colorArray[ci * 3], g: grp.colorArray[ci * 3 + 1], b: grp.colorArray[ci * 3 + 2] });
            }
          }
          bundle.espMeshes[f] = { vertexArr: verts, faceArr: faces, color: cols.length ? cols : undefined };
        }
        try { bg.removeSurface(sId); } catch (e) {}
      }
      if (onFramePrewarmed) onFramePrewarmed(f);
    }
  } catch (err) {
    console.warn('[ViewerSurface] ESP prewarm error:', err);
  }
}
