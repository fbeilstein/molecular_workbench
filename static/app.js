/* ═══════════════════════════════════════════════════════════════════
   Organic Reactions Workbench — Main Application Shell
   ═══════════════════════════════════════════════════════════════════ */

window.WB = window.WB || {};

// ── Init ────────────────────────────────────────────────────────

window.WB.init = function() {
    window.WB.initResizer();
    window.WB.viewer = $3Dmol.createViewer('viewer-3d', {
        backgroundColor: '#0f0f1a',
    });
    window.WB.status('Ready — draw a reaction, then click 📦 Compute Script');

    window.WB.refreshDirs();
};

window.WB.refreshDirs = function() {
    fetch('/api/dirs').then(r => r.json()).then(dirs => {
        const chooser = document.getElementById('run-chooser');
        const currentVal = chooser.value;
        chooser.innerHTML = '<option value="">-- Select Run --</option>';
        
        dirs.forEach(d => {
            const opt = document.createElement('option');
            opt.value = opt.textContent = d;
            chooser.appendChild(opt);
        });
        if (currentVal && dirs.includes(currentVal)) {
            chooser.value = currentVal;
        }
    }).catch(console.error);
};

// ── SMILES helpers ──────────────────────────────────────────────

window.WB.getSmiles = async function() {
    try {
        const frame = document.getElementById('ketcher-frame');
        const ketcher = frame.contentWindow.ketcher;
        return await ketcher.getSmiles();
    } catch (e) {
        window.WB.status('Error reading SMILES from Ketcher');
        return null;
    }
};

// ── Compute Script ──────────────────────────────────────────────

let pollInterval = null;
let currentJobName = null;

window.WB.stopCompute = async function() {
    if (!currentJobName) return;
    await fetch('/api/job/stop', {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name: currentJobName }),
    });
    if (pollInterval) clearInterval(pollInterval);
    document.getElementById('btn-stop').style.display = 'none';
    window.WB.status('Computation stopped.');
};

window.WB.generateComputeScript = async function() {
    const name = document.getElementById('job-name').value || 'molecule';
    const charge = parseInt(document.getElementById('mol-charge').value) || 0;
    const engineSelect = document.getElementById('qm-engine');
    const engine = engineSelect ? engineSelect.value : 'xtb';
    const levelShiftCb = document.getElementById('pyscf-levelshift');
    const levelshift = levelShiftCb ? levelShiftCb.checked : false;
    const smiles = await window.WB.getSmiles();
    currentJobName = name;

    if (!smiles) {
        window.WB.status('Error: Draw a molecule in Ketcher first');
        return;
    }

    window.WB.status('Generating compute script...');
    try {
        const frame = document.getElementById('ketcher-frame');
        const ketcher = frame.contentWindow.ketcher;
        let rxn = '';
        let svg = '';
        let ket = '';
        let mol = '';
        try {
            mol = await ketcher.getMolfile();
            svg = await ketcher.generateImage(mol, { outputFormat: 'svg' });
            ket = await ketcher.getKet();
        } catch (err) {
            console.warn("Could not get MOL/SVG/KET from Ketcher", err);
        }

        const res = await fetch('/api/rxn/compute-script', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ smiles, name, charge, engine, levelshift, rxn, svg, ket, mol }),
        });
        const data = await res.json();
        if (data.error) { 
            window.WB.status(`Error: ${data.error}`);
            return; 
        }

        window.WB.status('Computation started...');
        document.getElementById('btn-stop').style.display = 'inline-block';
        window.WB.openLog();
        
        // Poll for status
        if (pollInterval) clearInterval(pollInterval);
        pollInterval = setInterval(async () => {
            const sres = await fetch(`/api/job/status?name=${name}`);
            const sdata = await sres.json();
            
            window.WB.clearLog();
            window.WB.appendLog('═══ Computation Log ═══\n');
            window.WB.appendLog(sdata.log || 'Starting up...');
            
            if (!sdata.running) {
                clearInterval(pollInterval);
                document.getElementById('btn-stop').style.display = 'none';
                if (sdata.ready) {
                    window.WB.status('✓ Computation finished!');
                    await window.WB.refreshDirs();
                    const chooser = document.getElementById('run-chooser');
                    if (chooser) chooser.value = name;
                    const jobNameInput = document.getElementById('job-name');
                    if (jobNameInput) jobNameInput.value = name;
                    await window.WB.verifyBundle();
                } else {
                    window.WB.status('Computation failed (bundle not ready).');
                }
            }
        }, 1000);

    } catch (e) {
        window.WB.status(`Error: ${e.message}`);
    }
};

// ── Verify Bundle ───────────────────────────────────────────────

window.WB.verifyBundle = async function() {
    const name = document.getElementById('job-name').value ||
                 document.getElementById('run-chooser').value;
    if (!name) {
        window.WB.status('Error: Enter a job name or select a job from Load Job');
        return;
    }

    window.WB.status(`Loading bundle: ${name}...`);
    try {
        const res = await fetch(`/api/bundle/${name}/manifest.json?t=${Date.now()}`);
        if (!res.ok) {
            window.WB.status(`Error: No bundle found for "${name}". Run the compute script first.`);
            return;
        }
        const manifest = await res.json();
        window._currentBundleManifest = manifest;

        // Clear viewer
        window.WB.viewer.removeAllModels();
        window.WB.viewer.removeAllShapes();
        window.WB.viewer.removeAllSurfaces();

        // Load KET layout back into Ketcher if available
        const frame = document.getElementById('ketcher-frame');
        const ketcherObj = (frame && frame.contentWindow) ? frame.contentWindow.ketcher : null;
        
        if (manifest.ket_file && ketcherObj) {
            try {
                const ketRes = await fetch(`/api/bundle/${name}/${manifest.ket_file}?t=${Date.now()}`);
                if (ketRes.ok) {
                    const ketData = await ketRes.text();
                    await ketcherObj.setMolecule(ketData);
                }
            } catch (err) {
                console.warn("Could not load KET file into Ketcher:", err);
                if (manifest.smiles) await ketcherObj.setMolecule(manifest.smiles);
            }
        } else if (ketcherObj && manifest.smiles) {
            // Fallback to SMILES for older bundles that lack a .ket file
            await ketcherObj.setMolecule(manifest.smiles);
        }

        // Build chemical chooser from bundle manifest
        const chooser = document.getElementById('chemical-chooser');
        chooser.innerHTML = '';
        document.getElementById('chemical-chooser-wrap').style.display = '';

        manifest.molecules.forEach(m => {
            const opt = document.createElement('option');
            opt.value = `bundle:${name}:mol:${m.key}`;
            opt.textContent = `${m.role === 'reactant' ? '⬅' : '➡'} ${m.smiles} (${m.role})`;
            chooser.appendChild(opt);
        });

        // Summary in log
        window.WB.appendLog('\n\n\n' + '═'.repeat(50) + '\n');
        window.WB.appendLog(`Job Result: ${manifest.title || name}\n`);
        window.WB.appendLog('═'.repeat(50));
        window.WB.appendLog(`SMILES: ${manifest.smiles}`);
        window.WB.appendLog(`Method: ${manifest.method}`);
        window.WB.appendLog(`Molecules: ${manifest.molecules.length}`);
        if (manifest.barrier_forward_kcal != null)
            window.WB.appendLog(`Barrier: ${manifest.barrier_forward_kcal.toFixed(1)} kcal/mol`);
        window.WB.appendLog(`Orbitals: ${manifest.has_orbitals ? manifest.n_orbital_cubes + ' cubes' : 'none'}`);
        manifest.molecules.forEach(m => {
            const nOrb = m.orbitals ? m.orbitals.length : 0;
            window.WB.appendLog(`  ${m.role} ${m.key}: ${m.smiles} (${nOrb} orbitals)`);
        });
        window.WB.appendLog('═'.repeat(50));
        window.WB.appendLog('Select items from the Chemical dropdown to verify 3D/orbitals.');
        window.WB.openLog();

        // Enable slider / reset button states
        const btnPurge = document.getElementById('btn-purge-cubes');
        if (btnPurge) {
            btnPurge.textContent = '🗑 Purge Cube Files';
            btnPurge.disabled = false;
        }
        const slider = document.getElementById('isovalue-slider');
        if (slider) {
            slider.disabled = false;
            slider.title = '';
        }

        // Auto-load first molecule
        window.WB.loadBundleItem(chooser.value);

        window.WB.status(`✓ Bundle loaded: ${manifest.molecules.length} molecules`);
    } catch(e) {
        window.WB.status(`Error loading bundle: ${e.message}`);
    }
};

// ── Bundle Item Loading ─────────────────────────────────────────

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

    if (mol.xyz) {
        window._currentMolXyz = mol.xyz;
        window.WB.viewer.addModel(mol.xyz, 'xyz');
        window.WB.viewer.setStyle({}, {stick:{radius:0.12}, sphere:{scale:0.25}});
        window.WB.viewer.zoomTo();
        window.WB.viewer.render();
    }

    // Build orbital toggles for hierarchical groups
    const sidebar = document.getElementById('orbital-toggles');
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
                row.appendChild(cb); row.appendChild(dot); row.appendChild(lbl);
                sidebar.appendChild(row);
            });
        });
    }

    // ESP Surface toggle
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
        espRow.appendChild(espCb); espRow.appendChild(espDot); espRow.appendChild(espLbl);
        sidebar.appendChild(espRow);

        espCb.onchange = async () => {
            if (espCb.checked) {
                if (!window._currentEspCache) {
                    try {
                        window._currentEspCache = await (await fetch(`${base}/${mol.esp_surface.file}`)).json();
                    } catch(e) { console.error('ESP load failed', e); return; }
                }
                if (window._currentEspCache && window._currentEspCache.vertices) {
                    window.WB.onChangeHandler();
                }
            } else {
                if (window.WB.onChangeHandler) window.WB.onChangeHandler();
            }
        };
    }

    let engineInfo = '';
    if (window._currentBundleManifest) {
        const eng = window._currentBundleManifest.engine || 'unknown';
        const meth = window._currentBundleManifest.method || 'unknown';
        engineInfo = `<br><span style="font-size:11px;color:#aaa;">Engine: ${eng} (${meth})</span>`;
    }
    document.getElementById('viewer-info').innerHTML = `${mol.smiles} &mdash; ${mol.role}${engineInfo}`;
};

// ── Load Chemical (dropdown handler) ─────────────────────────────

window.WB.loadChemical = function() {
    const chooser = document.getElementById('chemical-chooser');
    const val = chooser.value;
    if (!val) return;
    if (val.startsWith('bundle:')) {
        window.WB.loadBundleItem(val);
    }
};

// ── Load Job (dropdown) ──────────────────────────────────────────

window.WB.loadRun = function() {
    const dir = document.getElementById('run-chooser').value;
    if (!dir) return;
    document.getElementById('job-name').value = dir;
    window.WB.status(`Job "${dir}" selected. Loading...`);
    window.WB.verifyBundle();
};

// ── Log Panel ────────────────────────────────────────────────────

window.WB.toggleLog = function() {
    const d = document.getElementById('log-drawer');
    d.classList.toggle('open');
};
window.WB.openLog = function() {
    document.getElementById('log-drawer').classList.add('open');
};
window.WB.clearLog = function() {
    document.getElementById('log-content').textContent = '';
};
window.WB.appendLog = function(text) {
    const el = document.getElementById('log-content');
    el.textContent += text + '\n';
    el.scrollTop = el.scrollHeight;
};

window.WB.status = function(msg) {
    document.getElementById('status-bar').textContent = msg;
};

// Attach boot sequence
window.onload = window.WB.init;
