/**
 * Surface and Isosurface Rendering Engine for Standalone 3D Viewer.
 * Manages volumetric orbitals, dual-phase coloring, and ESP electrostatic surfaces.
 */

export function findOrbitalMesh(bundle, frameIdx, orbKey) {
  if (!bundle || !bundle.orbitalMeshes) return null;
  const meshes = bundle.orbitalMeshes;

  const k1 = `${frameIdx}_${orbKey}`;
  if (meshes[k1]) return meshes[k1];

  if (bundle.isSingleMolecule && meshes[orbKey]) {
    return meshes[orbKey];
  }

  for (const k in meshes) {
    if (k.startsWith(`${frameIdx}_`)) {
      const baseKey = k.replace(`${frameIdx}_`, '');
      if (baseKey === orbKey || k.includes(orbKey) || orbKey.includes(baseKey)) {
        return meshes[k];
      }
    }
  }

  if (bundle.isSingleMolecule) {
    for (const k in meshes) {
      const stripped = k.replace(/^\d+_/, '');
      if (stripped === orbKey || k.includes(orbKey) || orbKey.includes(stripped)) {
        return meshes[k];
      }
    }
  }
  return null;
}

export function renderActiveOrbitals(glViewer, bundle, activeOrbitals, frameIdx, isovalue) {
  if (!activeOrbitals || activeOrbitals.size === 0 || !bundle) return;
  const CC = window.$3Dmol ? window.$3Dmol.CC : null;

  for (const orbKey of activeOrbitals) {
    const sKey = String(orbKey).toLowerCase();
    const isAnti = sKey.includes('star') || sKey.includes('*') || sKey.includes('lumo') || sKey.includes('antibond') || sKey.includes('virtual');
    const posColorHex = isAnti ? '#d946ef' : '#38bdf8';
    const negColorHex = isAnti ? '#22c55e' : '#f97316';
    const posColorInt = isAnti ? 0xd946ef : 0x38bdf8;
    const negColorInt = isAnti ? 0x22c55e : 0xf97316;

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
  const espMesh = bundle.espMeshes ? (
    bundle.espMeshes[frameIdx] ||
    bundle.espMeshes[String(frameIdx)] ||
    (bundle.isSingleMolecule ? (bundle.espMeshes['0'] || Object.values(bundle.espMeshes)[0]) : null)
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
      glViewer.addSurface(window.$3Dmol.SurfaceType.VDW, {
        opacity: 0.65,
        colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
      });
    } catch (e) {
      console.warn('Dynamic surface rendering failed:', e);
    }
  }
}

export function renderAtomCharges(glViewer, bundle, frame, symbols, coords, natoms, frameIdx) {
  const traj = bundle.trajectory || {};
  const trajCharges = traj.trajectory_charges ? traj.trajectory_charges[frameIdx] : null;
  const charges = trajCharges || frame.charges || traj.charges || bundle.charges || bundle.manifest?.charges || [];
  for (let i = 0; i < natoms; i++) {
    let q = (charges && charges[i] !== undefined) ? charges[i] : null;
    if (q === null) {
      const s = symbols[i];
      q = s === 'H' ? 0.06 : (s === 'O' ? -0.50 : (s === 'N' ? -0.40 : 0.0));
    }
    const sign = q > 0 ? '+' : '';
    glViewer.addLabel(`${symbols[i]}: ${sign}${Number(q).toFixed(2)}`, {
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
