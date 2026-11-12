/**
 * Surface and Isosurface Rendering Engine for Standalone 3D Viewer.
 * Manages volumetric orbitals, dual-phase coloring, and ESP electrostatic surfaces.
 */

function findTrack(bundle, orbKey) {
  if (!bundle || !orbKey) return null;
  const tracks = (bundle.manifest && bundle.manifest.orbital_timeline && bundle.manifest.orbital_timeline.tracks) ||
                 (bundle.trajectory && bundle.trajectory.orbital_timeline && bundle.trajectory.orbital_timeline.tracks) ||
                 (bundle.orbital_timeline && bundle.orbital_timeline.tracks) || [];
  const sTarget = String(orbKey).toLowerCase();
  for (const trk of tracks) {
    if (!trk) continue;
    if (String(trk.track_id || '').toLowerCase() === sTarget ||
        String(trk.id || '').toLowerCase() === sTarget ||
        String(trk.orbital_key || '').toLowerCase() === sTarget ||
        String(trk.local_key || '').toLowerCase() === sTarget ||
        String(trk.chem_key || '').toLowerCase() === sTarget) {
      return trk;
    }
  }
  return null;
}

function findOrbitalMesh(bundle, frameIdx, orbKey) {
  if (!bundle || !bundle.orbitalMeshes) return null;
  const meshes = bundle.orbitalMeshes;
  const exactKey = `${frameIdx}_${orbKey}`;
  if (meshes[exactKey]) return meshes[exactKey];

  const track = findTrack(bundle, orbKey);
  if (track) {
    if (track.frame_map && track.frame_map[frameIdx]) {
      const mapKey = `${frameIdx}_${track.frame_map[frameIdx]}`;
      if (meshes[mapKey]) return meshes[mapKey];
    }
    if (track.orbital_key) {
      const altKey = `${frameIdx}_${track.orbital_key}`;
      if (meshes[altKey]) return meshes[altKey];
    }
    if (track.chem_key) {
      const altKey = `${frameIdx}_${track.chem_key}`;
      if (meshes[altKey]) return meshes[altKey];
    }
    if (track.local_key) {
      const altKey = `${frameIdx}_${track.local_key}`;
      if (meshes[altKey]) return meshes[altKey];
    }
  }

  if (bundle.isSingleMolecule) {
    if (meshes[orbKey]) return meshes[orbKey];
    if (meshes[`0_${orbKey}`]) return meshes[`0_${orbKey}`];
  }
  return null;
}

function renderActiveOrbitals(glViewer, bundle, activeOrbitals, frameIdx, isovalue) {
  if (!activeOrbitals || activeOrbitals.size === 0 || !bundle) return;
  const CC = window.$3Dmol ? window.$3Dmol.CC : null;

  for (const orbKey of activeOrbitals) {
    const sKey = String(orbKey).toLowerCase();
    const track = findTrack(bundle, orbKey);

    const trackClass = String(track?.orbital_class || track?.category || track?.type || '').toLowerCase();
    const trackLabel = String(track?.label || track?.title || track?.short_label || track?.role_transition || '').toLowerCase();
    const isAnti = trackClass.includes('antibond') ||
                   trackClass.includes('star') ||
                   trackClass.includes('*') ||
                   trackClass.includes('lumo') ||
                   trackClass.includes('virtual') ||
                   trackLabel.includes('star') ||
                   trackLabel.includes('*') ||
                   trackLabel.includes('lumo') ||
                   trackLabel.includes('antibond') ||
                   trackLabel.includes('virtual') ||
                   sKey.includes('star') ||
                   sKey.includes('*') ||
                   sKey.includes('lumo') ||
                   sKey.includes('antibond') ||
                   sKey.includes('virtual');

    // Dual-phase coloring strictly matches separate molecule inspection:
    // Antibonding: Orange (#fb923c) pos, Magenta (#db2777) neg
    // Bonding:     Cyan   (#38bdf8) pos, Deep Blue (#1d4ed8) neg
    let posColorHex = isAnti ? '#fb923c' : '#38bdf8';
    let negColorHex = isAnti ? '#db2777' : '#1d4ed8';

    if (track && track.palette && track.palette.pos && track.palette.neg) {
      if (isAnti) {
        const pPos = track.palette.pos.toLowerCase();
        const pNeg = track.palette.neg.toLowerCase();
        if (pPos !== '#38bdf8' && pPos !== '#1d4ed8') posColorHex = track.palette.pos;
        if (pNeg !== '#1d4ed8' && pNeg !== '#38bdf8') negColorHex = track.palette.neg;
      } else {
        const pPos = track.palette.pos.toLowerCase();
        const pNeg = track.palette.neg.toLowerCase();
        if (pPos !== '#fb923c' && pPos !== '#db2777') posColorHex = track.palette.pos;
        if (pNeg !== '#db2777' && pNeg !== '#fb923c') negColorHex = track.palette.neg;
      }
    }

    const posColorInt = parseInt(posColorHex.replace('#', ''), 16);
    const negColorInt = parseInt(negColorHex.replace('#', ''), 16);

    const meshPair = findOrbitalMesh(bundle, frameIdx, orbKey);
    if (meshPair) {
      const pV = meshPair.pos?.vertices || meshPair.pos?.vertexArr;
      const pF = meshPair.pos?.faces || meshPair.pos?.faceArr;
      const pN = meshPair.pos?.normals || meshPair.pos?.normalArr;
      if (pV && pF) {
        glViewer.addCustom({
          vertexArr: pV, faceArr: pF, normalArr: pN,
          color: CC ? CC.color(posColorHex) : posColorInt, opacity: 0.90
        });
      }
      const nV = meshPair.neg?.vertices || meshPair.neg?.vertexArr;
      const nF = meshPair.neg?.faces || meshPair.neg?.faceArr;
      const nN = meshPair.neg?.normals || meshPair.neg?.normalArr;
      if (nV && nF) {
        glViewer.addCustom({
          vertexArr: nV, faceArr: nF, normalArr: nN,
          color: CC ? CC.color(negColorHex) : negColorInt, opacity: 0.90
        });
      }
    } else if (bundle.cubeCache) {
      const cubeData = bundle.cubeCache[`${frameIdx}_${orbKey}`] || bundle.cubeCache[orbKey];
      if (cubeData && window.$3Dmol) {
        try {
          const vdo = new window.$3Dmol.VolumeData(cubeData, 'cube');
          const iso = isovalue || 0.035;
          glViewer.addIsosurface(vdo, { isoval: iso, color: posColorHex, opacity: 0.90 });
          glViewer.addIsosurface(vdo, { isoval: -iso, color: negColorHex, opacity: 0.90 });
        } catch (_) {}
      }
    }
  }
}

// Patch 3Dmol once to support ambient lighting on ESP semi-transparent surfaces
function ensureEspLightingSupport(glViewer) {
  if (window.$3Dmol && window.$3Dmol.GLShape && !window.$3Dmol.GLShape._espPatchApplied) {
    window.$3Dmol.GLShape._espPatchApplied = true;
    const origGlobj = window.$3Dmol.GLShape.prototype.globj;
    window.$3Dmol.GLShape.prototype.globj = function(e) {
      const ret = origGlobj.call(this, e);
      if (this.stylespec && this.stylespec.isEsp && this.shapeObj) {
        const mesh = this.shapeObj.children && this.shapeObj.children[0];
        if (mesh && mesh.material) {
          mesh.material.isEsp = true;
        }
      }
      return ret;
    };
  }

  if (glViewer && glViewer.renderer && !glViewer.renderer._espPatchApplied) {
    glViewer.renderer._espPatchApplied = true;
    const origInit = glViewer.renderer.initMaterial;
    glViewer.renderer.initMaterial = function(mat, t, i, r) {
      origInit.call(this, mat, t, i, r);
      if (mat.isEsp && mat.vertexShader) {
        mat.vertexShader = mat.vertexShader
          .replace(
            'vec3 directionalLightWeighting = vec3( max( dotProduct, 0.0 ) );',
            'vec3 directionalLightWeighting = vec3( 0.65 + 0.35 * max( dotProduct, 0.0 ) );'
          )
          .replace(
            'vec3 directionalLightWeightingBack = vec3( max( -dotProduct, 0.0 ) );',
            'vec3 directionalLightWeightingBack = vec3( 0.65 + 0.35 * max( -dotProduct, 0.0 ) );'
          );
        const n = { wireframe: mat.wireframe, fragdepth: mat.imposter, volumetric: mat.volumetric, shaded: mat.shaded };
        mat.program = this.buildProgram(mat.fragmentShader, mat.vertexShader, mat.uniforms, n);
      }
    };
  }
}

function getGlobalEspRange(bundle) {
  if (!bundle) return { min: -30, max: 30, isCharged: false };
  if (bundle._globalEspRange) return bundle._globalEspRange;

  let gMin = Infinity;
  let gMax = -Infinity;
  let hasValues = false;

  if (bundle.espMeshes) {
    const keys = Object.keys(bundle.espMeshes);
    for (const key of keys) {
      const mesh = bundle.espMeshes[key];
      if (!mesh) continue;
      if (mesh.esp_min !== undefined && mesh.esp_max !== undefined) {
        gMin = Math.min(gMin, mesh.esp_min);
        gMax = Math.max(gMax, mesh.esp_max);
        hasValues = true;
      } else if (mesh.esp_values && mesh.esp_values.length > 0) {
        for (let i = 0; i < mesh.esp_values.length; i++) {
          const v = mesh.esp_values[i];
          if (v < gMin) gMin = v;
          if (v > gMax) gMax = v;
        }
        hasValues = true;
      }
    }
  }

  if (!hasValues || !isFinite(gMin) || !isFinite(gMax)) {
    gMin = -30.0;
    gMax = 30.0;
  }

  const isCharged = (gMax < -10) || (gMin > 10);
  const result = { min: gMin, max: gMax, isCharged };
  bundle._globalEspRange = result;
  return result;
}

function renderEspLayer(glViewer, bundle, frame, symbols, coords, natoms, frameIdx) {
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
    ensureEspLightingSupport(glViewer);

    if (!espMesh._orientedMesh) {
      // 1. Calculate centroid
      let cx = 0, cy = 0, cz = 0;
      const vLen = eV.length;
      for (let i = 0; i < vLen; i++) {
        cx += eV[i].x; cy += eV[i].y; cz += eV[i].z;
      }
      cx /= vLen; cy /= vLen; cz /= vLen;

      // 2. Check face winding: if inward (CW from outside), reverse so CCW outward
      let faceDotSum = 0;
      const fLen = eF.length;
      const sampleStepF = Math.max(3, Math.floor(fLen / 300) * 3);
      for (let i = 0; i < fLen; i += sampleStepF) {
        const v0 = eV[eF[i]], v1 = eV[eF[i+1]], v2 = eV[eF[i+2]];
        if (!v0 || !v1 || !v2) continue;
        const ax = v0.x - v1.x, ay = v0.y - v1.y, az = v0.z - v1.z;
        const hx = v2.x - v1.x, hy = v2.y - v1.y, hz = v2.z - v1.z;
        const nx = hy * az - hz * ay, ny = hz * ax - hx * az, nz = hx * ay - hy * ax;
        const mx = (v0.x + v1.x + v2.x) / 3 - cx;
        const my = (v0.y + v1.y + v2.y) / 3 - cy;
        const mz = (v0.z + v1.z + v2.z) / 3 - cz;
        faceDotSum += nx * mx + ny * my + nz * mz;
      }

      let finalFaces = eF;
      if (faceDotSum < 0) {
        finalFaces = new (eF.constructor)(fLen);
        for (let i = 0; i < fLen; i += 3) {
          finalFaces[i] = eF[i];
          finalFaces[i+1] = eF[i+2];
          finalFaces[i+2] = eF[i+1];
        }
      }

      // 3. Check vertex normals: if inward, negate to point outward
      let finalNormals = eN;
      if (eN && eN.length === vLen) {
        let normDotSum = 0;
        const sampleStepN = Math.max(1, Math.floor(vLen / 100));
        for (let i = 0; i < vLen; i += sampleStepN) {
          normDotSum += (eV[i].x - cx) * eN[i].x + (eV[i].y - cy) * eN[i].y + (eV[i].z - cz) * eN[i].z;
        }
        if (normDotSum < 0) {
          finalNormals = eN.map(n => ({ x: -n.x, y: -n.y, z: -n.z }));
        }
      }

      espMesh._orientedMesh = { faces: finalFaces, normals: finalNormals };
    }

    const finalFaces = espMesh._orientedMesh.faces;
    const finalNormals = espMesh._orientedMesh.normals;

    // 4. Adaptive colormapping
    let colorParam = 0x10b981;
    if (espMesh.esp_values && espMesh.esp_values.length === eV.length) {
      let minVal, maxVal, isCharged;
      if (!bundle.isSingleMolecule) {
        const gRange = getGlobalEspRange(bundle);
        minVal = gRange.min;
        maxVal = gRange.max;
        isCharged = gRange.isCharged;
      } else {
        minVal = espMesh.esp_min !== undefined ? espMesh.esp_min : Math.min(...espMesh.esp_values);
        maxVal = espMesh.esp_max !== undefined ? espMesh.esp_max : Math.max(...espMesh.esp_values);
        isCharged = (maxVal < -10) || (minVal > 10);
      }
      const span = Math.max(1e-5, maxVal - minVal);
      const absMax = Math.max(Math.abs(minVal || -30), Math.abs(maxVal || 30), 20);

      colorParam = espMesh.esp_values.map(v => {
        let t;
        if (isCharged) {
          t = -1.0 + 2.0 * ((v - minVal) / span);
        } else {
          t = Math.max(-1, Math.min(1, v / absMax));
        }
        if (t < 0) {
          const f = 1.0 + t;
          return {
            r: 1.0,
            g: Math.min(1, Math.max(0.04, 0.04 + 0.96 * f)),
            b: Math.min(1, Math.max(0.04, 0.04 + 0.96 * f))
          };
        } else {
          const f = 1.0 - t;
          return {
            r: Math.min(1, Math.max(0.04, 0.04 + 0.96 * f)),
            g: Math.min(1, Math.max(0.20, 0.20 + 0.80 * f)),
            b: 1.0
          };
        }
      });
    } else if (espMesh.color && Array.isArray(espMesh.color)) {
      // Enhance vividness and saturation of precomputed mesh colors
      colorParam = espMesh.color.map(c => {
        const r = typeof c.r === 'number' ? c.r : (c[0] !== undefined ? c[0] : 1);
        const g = typeof c.g === 'number' ? c.g : (c[1] !== undefined ? c[1] : 1);
        const b = typeof c.b === 'number' ? c.b : (c[2] !== undefined ? c[2] : 1);
        const minC = Math.min(r, g, b);
        const maxC = Math.max(r, g, b);
        const span = maxC - minC;
        if (span < 0.005 || maxC < 0.01) return c;
        const s = span / maxC;
        const sNew = Math.min(1.0, Math.pow(s, 0.60));
        const f = sNew / s;
        return {
          r: Math.max(0, Math.min(1, maxC - (maxC - r) * f)),
          g: Math.max(0, Math.min(1, maxC - (maxC - g) * f)),
          b: Math.max(0, Math.min(1, maxC - (maxC - b) * f))
        };
      });
    }

    const customObj = {
      vertexArr: eV,
      faceArr: finalFaces,
      color: colorParam,
      opacity: 0.85,
      isEsp: true
    };
    if (finalNormals && finalNormals.length === eV.length) customObj.normalArr = finalNormals;
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

function renderAtomCharges(glViewer, bundle, frame, symbols, coords, natoms, frameIdx) {
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

function updateEspColorbar(container, active, bundle, frameIdx) {
  const colorbar = container?.parentElement?.querySelector('.esp-colorbar') ||
                   container?.parentElement?.querySelector('#esp-colorbar') ||
                   document.querySelector('#modal-esp-colorbar') ||
                   document.getElementById('esp-colorbar');
  if (!colorbar) return;
  colorbar.style.display = active ? 'flex' : 'none';
  if (active && bundle) {
    const minEl = colorbar.querySelector('#esp-colorbar-min') || document.getElementById('esp-colorbar-min');
    const midEl = colorbar.querySelector('#esp-colorbar-mid') || document.getElementById('esp-colorbar-mid');
    const maxEl = colorbar.querySelector('#esp-colorbar-max') || document.getElementById('esp-colorbar-max');

    let minVal, maxVal, isCharged;
    if (!bundle.isSingleMolecule) {
      const gRange = getGlobalEspRange(bundle);
      minVal = gRange.min;
      maxVal = gRange.max;
      isCharged = gRange.isCharged;
    } else {
      const espMesh = bundle.espMeshes ? (
        bundle.espMeshes[frameIdx] ||
        bundle.espMeshes[String(frameIdx)] ||
        (bundle.espMeshes['0'] || Object.values(bundle.espMeshes)[0])
      ) : null;
      if (espMesh && espMesh.esp_values && espMesh.esp_values.length > 0) {
        minVal = espMesh.esp_min !== undefined ? espMesh.esp_min : Math.min(...espMesh.esp_values);
        maxVal = espMesh.esp_max !== undefined ? espMesh.esp_max : Math.max(...espMesh.esp_values);
        isCharged = (maxVal < -10) || (minVal > 10);
      } else {
        minVal = -30.0;
        maxVal = 30.0;
        isCharged = false;
      }
    }

    if (isCharged) {
      if (minEl) minEl.textContent = `${minVal.toFixed(1)}`;
      if (midEl) midEl.textContent = `${((minVal + maxVal) / 2).toFixed(1)}`;
      if (maxEl) maxEl.textContent = `${maxVal > 0 ? '+' : ''}${maxVal.toFixed(1)}`;
    } else {
      const absMax = Math.max(Math.abs(minVal || -30), Math.abs(maxVal || 30), 20);
      if (minEl) minEl.textContent = `${(-absMax).toFixed(1)}`;
      if (midEl) midEl.textContent = '0';
      if (maxEl) maxEl.textContent = `+${absMax.toFixed(1)}`;
    }
  }
}

let _bgViewerDiv = null;
let _bgViewerInstance = null;

async function prewarmViewerEspCache(bundle, onFramePrewarmed) {
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

window.findTrack = findTrack;
window.findOrbitalMesh = findOrbitalMesh;
window.renderActiveOrbitals = renderActiveOrbitals;
window.renderEspLayer = renderEspLayer;
window.renderAtomCharges = renderAtomCharges;
window.getGlobalEspRange = getGlobalEspRange;
window.updateEspColorbar = updateEspColorbar;
window.prewarmViewerEspCache = prewarmViewerEspCache;
