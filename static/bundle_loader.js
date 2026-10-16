/* ═══════════════════════════════════════════════════════════════════
   Molecular Workbench — Bundle Item & Orbital Sidebar Loader
   ═══════════════════════════════════════════════════════════════════ */

window.WB = window.WB || {};

async function loadEspSurfaceData(base, espFile) {
    if (window._currentEspCache) return window._currentEspCache;
    try {
        const res = await fetch(`${base}/${espFile}`);
        window._currentEspCache = await res.json();
        return window._currentEspCache;
    } catch (e) {
        console.error('ESP load failed', e);
        return null;
    }
}

function renderOrbitalSidebar(mol, base) {
    const sidebar = document.getElementById('orbital-toggles');
    if (!sidebar) return;
    sidebar.innerHTML = '';

    if (mol.orbitals && mol.orbitals.length) {
        window._currentOrbitalCache = {};
        window._currentEspCache = null;

        mol.orbitals.forEach(group => {
            if (!group.items || !group.items.length) return;

            const groupHeader = document.createElement('div');
            groupHeader.style.cssText = 'font-weight:bold;color:#fff;margin-top:10px;margin-bottom:4px;padding-left:4px;font-size:13px;border-bottom:1px solid #335;';
            groupHeader.textContent = group.name;
            sidebar.appendChild(groupHeader);

            group.items.forEach(orb => {
                const row = document.createElement('label');
                row.className = 'orbital-toggle';
                row.style.cssText = 'display:flex;align-items:center;gap:6px;padding:4px 8px;cursor:pointer;';
                const cb = document.createElement('input');
                cb.type = 'checkbox';
                cb.onchange = window.WB.onChangeHandler;
                cb.setAttribute('data-orb-file', orb.file);
                cb.setAttribute('data-orb-color', orb.color);
                const dot = document.createElement('span');
                dot.style.cssText = `display:inline-block;width:10px;height:10px;border-radius:50%;background:${orb.color};`;
                const lbl = document.createElement('span');
                lbl.style.cssText = 'font-size:12px;color:#ccc;';
                lbl.textContent = orb.energy_ev !== undefined ? `${orb.label} (${orb.energy_ev} eV)` : orb.label;
                row.appendChild(cb);
                row.appendChild(dot);
                row.appendChild(lbl);
                sidebar.appendChild(row);
            });
        });
    }

    if (mol.esp_surface) {
        const espHeader = document.createElement('div');
        espHeader.style.cssText = 'font-weight:bold;color:#fff;margin-top:10px;margin-bottom:4px;padding-left:4px;font-size:13px;border-bottom:1px solid #335;';
        espHeader.textContent = 'Electrostatic Potential';
        sidebar.appendChild(espHeader);

        const espRow = document.createElement('label');
        espRow.className = 'orbital-toggle';
        espRow.style.cssText = 'display:flex;align-items:center;gap:6px;padding:4px 8px;cursor:pointer;';
        const espCb = document.createElement('input');
        espCb.type = 'checkbox';
        espCb.setAttribute('data-esp', 'true');
        const espDot = document.createElement('span');
        espDot.style.cssText = 'display:inline-block;width:10px;height:10px;border-radius:50%;background:linear-gradient(90deg, #ff3333, #ffffff, #3333ff);';
        const espLbl = document.createElement('span');
        espLbl.style.cssText = 'font-size:12px;color:#ccc;';
        espLbl.textContent = `ESP Surface (${mol.esp_surface.esp_min} to ${mol.esp_surface.esp_max} kcal/mol)`;
        espRow.appendChild(espCb);
        espRow.appendChild(espDot);
        espRow.appendChild(espLbl);
        sidebar.appendChild(espRow);

        espCb.onchange = async () => {
            if (!espCb.checked) {
                if (window.WB.onChangeHandler) window.WB.onChangeHandler();
                return;
            }
            const cache = await loadEspSurfaceData(base, mol.esp_surface.file);
            if (cache && cache.vertices && window.WB.onChangeHandler) {
                window.WB.onChangeHandler();
            }
        };
    }

    const chargesRow = document.createElement('label');
    chargesRow.className = 'orbital-toggle';
    chargesRow.style.cssText = 'display:flex;align-items:center;gap:6px;padding:4px 8px;cursor:pointer;';
    const chargesCb = document.createElement('input');
    chargesCb.type = 'checkbox';
    chargesCb.setAttribute('data-charges', 'true');
    const chargesDot = document.createElement('span');
    chargesDot.style.cssText = 'display:inline-block;width:10px;height:10px;border-radius:50%;background:#a855f7;';
    const chargesLbl = document.createElement('span');
    chargesLbl.style.cssText = 'font-size:12px;color:#ccc;';
    chargesLbl.textContent = '🏷 Atom Charges';
    chargesRow.appendChild(chargesCb);
    chargesRow.appendChild(chargesDot);
    chargesRow.appendChild(chargesLbl);
    sidebar.appendChild(chargesRow);

    chargesCb.onchange = () => {
        if (window.WB.onChangeHandler) window.WB.onChangeHandler();
    };
}

window.WB.loadBundleItem = async function(value) {
    if (!value || !value.startsWith('bundle:')) return;
    const parts = value.split(':');
    const bundleName = parts[1];
    const type = parts[2];
    const base = `/api/bundle/${bundleName}`;
    const cacheBuster = `?t=${Date.now()}`;

    if (type !== 'mol') return;

    const key = parts[3];
    const mol = await (await fetch(`${base}/molecules/${key}.json${cacheBuster}`)).json();

    window.WB.viewer.removeAllModels();
    window.WB.viewer.removeAllShapes();
    window.WB.viewer.removeAllSurfaces();

    window._currentMolCharges = mol.charges || (window._currentBundleManifest && window._currentBundleManifest.charges) || null;

    if (mol.xyz) {
        window._currentMolXyz = mol.xyz;
        window.WB.viewer.addModel(mol.xyz, 'xyz');
        window.WB.viewer.setStyle({}, { stick: { radius: 0.12 }, sphere: { scale: 0.25 } });
        window.WB.viewer.zoomTo();
        window.WB.viewer.render();
    }

    renderOrbitalSidebar(mol, base);

    let engineInfo = '';
    if (window._currentBundleManifest) {
        const eng = window._currentBundleManifest.engine || 'unknown';
        const meth = window._currentBundleManifest.method || 'unknown';
        engineInfo = `<br><span style="font-size:11px;color:#aaa;">Engine: ${eng} (${meth})</span>`;
    }
    const infoEl = document.getElementById('viewer-info');
    if (infoEl) {
        infoEl.innerHTML = `${mol.smiles} &mdash; ${mol.role}${engineInfo}`;
    }
};
