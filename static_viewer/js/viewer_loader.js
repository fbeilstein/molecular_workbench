/**
 * viewer_loader.js
 * In-memory bundle fetcher and full zip unpacker using JSZip.
 * Unpacks all trajectory coords, SVG, orbital meshes, and ESP surfaces upfront for zero-lag playback.
 */

export async function loadReactionBundle(source, onProgress) {
  let arrayBuffer;

  if (typeof source === 'string') {
    if (onProgress) onProgress(10, 'Fetching bundle from ' + source + '...');
    const resp = await fetch(source, { mode: 'cors' });
    if (!resp.ok) {
      throw new Error(`Failed to download bundle (${resp.status} ${resp.statusText})`);
    }
    arrayBuffer = await resp.arrayBuffer();
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
