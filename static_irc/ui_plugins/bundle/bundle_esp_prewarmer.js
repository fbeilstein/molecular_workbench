/**
 * Background ESP Surface Cache Prewarmer for Step 5 / Step 7 Bundle Viewer.
 */

let _bgDiv = null;
let _bgViewer = null;

export async function prewarmEspCache(symbols, frames, trajectoryCharges, espMeshCache, onProgress) {
    if (!symbols || !frames || !frames.length || typeof $3Dmol === 'undefined') return;

    try {
        if (!_bgDiv) {
            _bgDiv = document.createElement('div');
            _bgDiv.style.cssText = 'position:fixed;left:-9999px;width:120px;height:120px;opacity:0;pointer-events:none;';
            document.body.appendChild(_bgDiv);
            _bgViewer = $3Dmol.createViewer(_bgDiv);
        }
        const bg = _bgViewer;
        const total = frames.length, n = symbols.length;

        for (let f = 0; f < total; f++) {
            if (espMeshCache[f]) continue;
            const coords = frames[f]?.coords || frames[f]?.positions;
            const charges = trajectoryCharges?.[f] || frames[f]?.charges;
            if (!coords || !charges || charges.length < n) continue;

            bg.clear();
            let mdl = `${n}\nFrame Model\n`;
            for (let i = 0; i < n; i++) {
                const [x, y, z] = coords[i];
                mdl += `${symbols[i]} ${x.toFixed(4)} ${y.toFixed(4)} ${z.toFixed(4)}\n`;
            }
            bg.addModel(mdl, 'xyz');
            const atoms = bg.selectedAtoms({});
            for (let i = 0; i < atoms.length && i < charges.length; i++) {
                atoms[i].partialCharge = charges[i] || 0.0;
                atoms[i].properties = { partialCharge: charges[i] || 0.0 };
            }
            const p = bg.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.65,
                colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
            });
            if (p && typeof p.then === 'function') {
                const sId = await p;
                const grp = bg.surfaces?.[sId]?.geo?.geometryGroups?.[0];
                if (grp && grp.vertexArray && grp.faceArray) {
                    const vCount = grp.vertices || Math.floor(grp.vertexArray.length / 3);
                    const verts = [];
                    for (let vi = 0; vi < vCount; vi++) {
                        verts.push({ x: grp.vertexArray[vi * 3], y: grp.vertexArray[vi * 3 + 1], z: grp.vertexArray[vi * 3 + 2] });
                    }
                    const faces = Array.from(grp.faceArray.slice(0, grp.faces * 3));
                    const cols = [];
                    if (grp.colorArray && grp.colorArray.length >= vCount * 3) {
                        for (let ci = 0; ci < vCount; ci++) {
                            cols.push({ r: grp.colorArray[ci * 3], g: grp.colorArray[ci * 3 + 1], b: grp.colorArray[ci * 3 + 2] });
                        }
                    }
                    espMeshCache[f] = { vertexArr: verts, faceArr: faces, color: cols.length ? cols : undefined };
                }
                try { bg.removeSurface(sId); } catch (e) {}
            }
            if (onProgress) onProgress(f + 1, total);
        }
    } catch (e) {
        console.warn('[BundleESPPrewarmer] ESP prewarm error:', e);
    }
}
