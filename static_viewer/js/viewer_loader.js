/**
 * viewer_loader.js
 * In-memory bundle fetcher and full zip unpacker using JSZip.
 * Unpacks all trajectory coords, SVG, orbital meshes, and ESP surfaces upfront for zero-lag playback.
 */

export async function loadReactionBundle(source, onProgress) {
  let arrayBuffer;

  if (typeof source === 'string') {
    let url = source.trim();

    // 1. If pointing to a release zip from this repo, check if it's hosted locally in ./bundles/
    if (url.includes('/releases/download/') && url.endsWith('.zip')) {
      const fname = url.split('/').pop().split('?')[0];
      const localCandidate = `./bundles/${fname}`;
      try {
        if (onProgress) onProgress(10, `Checking local repository bundles for ${fname}...`);
        const check = await fetch(localCandidate, { method: 'HEAD' });
        if (check.ok) {
          url = localCandidate;
        }
      } catch (_) {}
    }

    if (onProgress) onProgress(15, 'Fetching bundle from ' + url + '...');

    // 2. Direct fetch
    try {
      const resp = await fetch(url, { mode: 'cors' });
      if (!resp.ok) {
        throw new Error(`Failed to download bundle (${resp.status} ${resp.statusText})`);
      }
      arrayBuffer = await resp.arrayBuffer();
    } catch (directErr) {
      // 3. Fallback proxy cascade for cross-origin URLs blocked by CORS
      if (url.startsWith('http://') || url.startsWith('https://')) {
        const proxies = [
          (u) => `https://api.allorigins.win/raw?url=${encodeURIComponent(u)}`,
          (u) => `https://corsproxy.org/?${encodeURIComponent(u)}`,
          (u) => `https://cors.eu.org/${u}`
        ];
        let fetched = false;
        for (let i = 0; i < proxies.length; i++) {
          try {
            if (onProgress) onProgress(20 + i * 5, `Direct fetch blocked (CORS). Retrying with proxy ${i + 1}/${proxies.length}...`);
            const pResp = await fetch(proxies[i](url));
            if (pResp.ok) {
              arrayBuffer = await pResp.arrayBuffer();
              fetched = true;
              break;
            }
          } catch (_) {}
        }
        if (!fetched) {
          throw new Error(`Remote server blocked direct download due to CORS policy.`);
        }
      } else {
        throw directErr;
      }
    }
  } else if (source instanceof Blob || source instanceof File) {
    if (onProgress) onProgress(10, 'Reading local bundle archive...');
    arrayBuffer = await source.arrayBuffer();
  } else if (source instanceof ArrayBuffer) {
    arrayBuffer = source;
  } else {
    throw new Error('Unsupported bundle source type');
  }

  if (onProgress) onProgress(30, 'Unpacking archive entries into memory...');
  if (!window.JSZip) {
    throw new Error('JSZip library is not loaded');
  }

  const zip = await window.JSZip.loadAsync(arrayBuffer);
  const bundle = {
    manifest: null,
    trajectory: null,
    mechanismSvg: '',
    svgWiring: null,
    orbitalMeshes: {},
    espMeshes: {},
    molecules: {}
  };

  const fileKeys = Object.keys(zip.files);
  const totalFiles = fileKeys.length;
  let processed = 0;

  for (const name of fileKeys) {
    const file = zip.files[name];
    if (file.dir) continue;

    processed++;
    if (onProgress && processed % 15 === 0) {
      const pct = 30 + Math.floor((processed / totalFiles) * 65);
      onProgress(pct, `Unpacking assets (${processed}/${totalFiles})...`);
    }

    if (name === 'bundle.json') {
      const text = await file.async('string');
      try { bundle.manifest = JSON.parse(text); } catch (e) { console.warn(e); }
    } else if (name === 'multistep_pathway.json') {
      const text = await file.async('string');
      try { bundle.trajectory = JSON.parse(text); } catch (e) { console.warn(e); }
    } else if (name === 'mechanism.svg') {
      bundle.mechanismSvg = await file.async('string');
    } else if (name === 'svg_wiring.json') {
      const text = await file.async('string');
      try { bundle.svgWiring = JSON.parse(text); } catch (e) { console.warn(e); }
    } else if (name.startsWith('orbitals/') && name.endsWith('.json')) {
      const base = name.replace('orbitals/', '');
      const text = await file.async('string');
      try {
        const meshData = JSON.parse(text);
        if (base.endsWith('_pos.json')) {
          const key = base.replace('_pos.json', '');
          bundle.orbitalMeshes[key] = bundle.orbitalMeshes[key] || {};
          bundle.orbitalMeshes[key].pos = meshData;
        } else if (base.endsWith('_neg.json')) {
          const key = base.replace('_neg.json', '');
          bundle.orbitalMeshes[key] = bundle.orbitalMeshes[key] || {};
          bundle.orbitalMeshes[key].neg = meshData;
        }
      } catch (e) {
        console.warn('Failed parsing orbital mesh:', name, e);
      }
    } else if (name.startsWith('esp/') && name.endsWith('.json')) {
      const base = name.replace('esp/', '').replace('.json', '');
      const frameKey = base.replace('frame_', '');
      const text = await file.async('string');
      try {
        bundle.espMeshes[frameKey] = JSON.parse(text);
      } catch (e) {
        console.warn('Failed parsing ESP mesh:', name, e);
      }
    } else if (name.startsWith('molecules/')) {
      const parts = name.split('/');
      if (parts.length > 2) {
        const molKey = parts[1];
        const sub = parts.slice(2).join('/');
        bundle.molecules[molKey] = bundle.molecules[molKey] || {};
        if (sub.endsWith('.json')) {
          const text = await file.async('string');
          try { bundle.molecules[molKey][sub] = JSON.parse(text); } catch (_) { bundle.molecules[molKey][sub] = text; }
        }
      }
    }
  }

  if (!bundle.trajectory && bundle.manifest && bundle.manifest.files && bundle.manifest.files.trajectory_json) {
    const trajFile = zip.file(bundle.manifest.files.trajectory_json);
    if (trajFile) {
      const text = await trajFile.async('string');
      bundle.trajectory = JSON.parse(text);
    }
  }

  if (!bundle.trajectory) {
    throw new Error('Archive does not contain a valid multistep_pathway.json');
  }

  if (onProgress) onProgress(100, 'Unpack complete. Initializing 3D viewer...');
  return bundle;
}

/**
 * Adapt a molecule sub-directory inside a reaction bundle into a standalone 1-state bundle.
 */
export function buildMoleculeBundle(molIdentifier, parentBundle) {
  if (!parentBundle || !parentBundle.molecules) return null;
  const mols = parentBundle.molecules;

  let molKey = null;
  let molDict = null;
  if (mols[molIdentifier]) {
    molKey = molIdentifier;
    molDict = mols[molIdentifier];
  } else {
    for (const k in mols) {
      if (k.toLowerCase() === String(molIdentifier).toLowerCase() ||
          k.toLowerCase().includes(String(molIdentifier).toLowerCase()) ||
          String(molIdentifier).toLowerCase().includes(k.toLowerCase())) {
        molKey = k;
        molDict = mols[k];
        break;
      }
    }
  }

  if (!molDict && parentBundle.trajectory && parentBundle.trajectory.molecule_bundles) {
    const mb = parentBundle.trajectory.molecule_bundles;
    const entry = mb[molIdentifier] || Object.values(mb).find(v => v.key === molIdentifier || v.directory?.includes(molIdentifier));
    if (entry && entry.directory) {
      const dirClean = entry.directory.replace(/^molecules\//, '').replace(/\/$/, '');
      if (mols[dirClean]) {
        molKey = dirClean;
        molDict = mols[dirClean];
      }
    }
  }

  if (!molDict) return null;

  const manifest = molDict['manifest.json'] || {};
  const molInfo = (manifest.molecules && manifest.molecules[0]) || molDict[`${molKey}.json`] || {};

  const xyzText = molInfo.xyz || '';
  const lines = xyzText.trim().split('\n');
  const symbols = [];
  const coords = [];
  if (lines.length >= 3) {
    for (let i = 2; i < lines.length; i++) {
      const parts = lines[i].trim().split(/\s+/);
      if (parts.length >= 4) {
        symbols.push(parts[0]);
        coords.push([parseFloat(parts[1]), parseFloat(parts[2]), parseFloat(parts[3])]);
      }
    }
  }

  const orbitalMeshes = {};
  const espMeshes = {};
  for (const [fname, content] of Object.entries(molDict)) {
    if (fname.endsWith('_pos.json')) {
      const base = fname.replace(/_pos\.json$/, '');
      orbitalMeshes[base] = orbitalMeshes[base] || {};
      orbitalMeshes[base].pos = content;
    } else if (fname.endsWith('_neg.json')) {
      const base = fname.replace(/_neg\.json$/, '');
      orbitalMeshes[base] = orbitalMeshes[base] || {};
      orbitalMeshes[base].neg = content;
    } else if (fname.endsWith('_esp.json')) {
      espMeshes['0'] = content;
    }
  }

  const rawCategories = molInfo.orbitals || manifest.orbitals || [];
  const categories = rawCategories.map(cat => ({
    name: cat.name,
    items: (cat.items || []).map(it => {
      const rawKey = it.file ? it.file.replace(/^molecules\//, '').replace(/\.json$/, '') : (it.key || it.label);
      return {
        key: rawKey,
        label: it.label || rawKey,
        color: it.color || '#38bdf8',
        energy_ev: it.energy_ev
      };
    })
  }));

  const title = molInfo.name || molInfo.key || molKey.replace(/_/g, ' ');

  return {
    manifest: {
      title,
      molKey,
      smiles: molInfo.smiles || '',
      isSingleMolecule: true
    },
    trajectory: {
      symbols,
      frames: [{ coords, energy_kcal: 0 }]
    },
    orbitalMeshes,
    espMeshes,
    categories,
    hasEsp: true,
    isSingleMolecule: true
  };
}
