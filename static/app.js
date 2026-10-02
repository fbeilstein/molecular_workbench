/* ═══════════════════════════════════════════════════════════════════
   Organic Reactions Workbench — UI Logic
   ═══════════════════════════════════════════════════════════════════ */

const WB = (() => {
    let viewer = null;

    // ── Init ────────────────────────────────────────────────────────

    
    // ── Resizer ─────────────────────────────────────────────────────────────
    function initResizer() {
        const resizer = document.getElementById('dragMe');
        const leftSide = document.querySelector('.ketcher-panel');
        const rightSide = document.querySelector('.viewer-panel');

        if(!resizer || !leftSide || !rightSide) return;

        let x = 0;
        let leftWidth = 0;

        const mouseDownHandler = function(e) {
            x = e.clientX;
            leftWidth = leftSide.getBoundingClientRect().width;
            document.addEventListener('mousemove', mouseMoveHandler);
            document.addEventListener('mouseup', mouseUpHandler);
            
            // disable pointer events on iframes so they don't swallow mouse events
            const iframes = document.querySelectorAll('iframe');
            iframes.forEach(ifr => ifr.style.pointerEvents = 'none');
            document.getElementById('viewer-3d').style.pointerEvents = 'none';
        };

        const mouseMoveHandler = function(e) {
            const dx = e.clientX - x;
            const containerWidth = resizer.parentNode.getBoundingClientRect().width;
            const newLeftWidth = ((leftWidth + dx) * 100) / containerWidth;
            leftSide.style.flex = `0 0 ${newLeftWidth}%`;
            rightSide.style.flex = '1 1 0%';
        };

        const mouseUpHandler = function() {
            document.removeEventListener('mousemove', mouseMoveHandler);
            document.removeEventListener('mouseup', mouseUpHandler);
            
            // re-enable pointer events
            const iframes = document.querySelectorAll('iframe');
            iframes.forEach(ifr => ifr.style.pointerEvents = '');
            document.getElementById('viewer-3d').style.pointerEvents = '';
        };

        resizer.addEventListener('mousedown', mouseDownHandler);

        // Log Resizer
        const logResizer = document.getElementById('logDragMe');
        const logDrawer = document.getElementById('log-drawer');
        if (logResizer && logDrawer) {
            let startY = 0;
            let startHeight = 0;

            const logMouseDownHandler = function(e) {
                if (!logDrawer.classList.contains('open')) return;
                startY = e.clientY;
                startHeight = logDrawer.getBoundingClientRect().height;
                logDrawer.classList.add('dragging');
                document.addEventListener('mousemove', logMouseMoveHandler);
                document.addEventListener('mouseup', logMouseUpHandler);
                
                const iframes = document.querySelectorAll('iframe');
                iframes.forEach(ifr => ifr.style.pointerEvents = 'none');
                const viewer = document.getElementById('viewer-3d');
                if(viewer) viewer.style.pointerEvents = 'none';
            };

            const logMouseMoveHandler = function(e) {
                const dy = startY - e.clientY; // moving up increases height
                let newHeight = startHeight + dy;
                if (newHeight < 100) newHeight = 100;
                logDrawer.style.setProperty('--log-height', newHeight + 'px');
            };

            const logMouseUpHandler = function() {
                logDrawer.classList.remove('dragging');
                document.removeEventListener('mousemove', logMouseMoveHandler);
                document.removeEventListener('mouseup', logMouseUpHandler);
                
                const iframes = document.querySelectorAll('iframe');
                iframes.forEach(ifr => ifr.style.pointerEvents = '');
                const viewer = document.getElementById('viewer-3d');
                if(viewer) viewer.style.pointerEvents = '';
            };

            logResizer.addEventListener('mousedown', logMouseDownHandler);

            logResizer.addEventListener('mousedown', logMouseDownHandler);
        }

    }

    
    function toggleExportSelectAll(checkbox) {
        const listDiv = document.getElementById('export-orbital-list');
        const chks = listDiv.querySelectorAll('input[type="checkbox"]');
        chks.forEach(chk => chk.checked = checkbox.checked);
    }

    async function downloadCustomBundle() {
        const name = document.getElementById('job-name').value;
        if (!name) {
            status('No job loaded to download.');
            return;
        }
        
        const listDiv = document.getElementById('export-orbital-list');
        listDiv.innerHTML = '';
        
        const sidebar = document.getElementById('orbital-toggles');
        if (sidebar) {
            const inputs = sidebar.querySelectorAll('input[data-orb-file], input[data-orb-key]');
            inputs.forEach(inp => {
                const fileAttr = inp.getAttribute('data-orb-file');
                const keyAttr = inp.getAttribute('data-orb-key');
                const text = inp.closest('label').innerText.trim();
                
                const div = document.createElement('div');
                div.style.marginBottom = '5px';
                
                const cb = document.createElement('input');
                cb.type = 'checkbox';
                if (fileAttr) {
                    cb.value = fileAttr;
                } else if (keyAttr) {
                    // For trajectory orbitals, they are stored in orbitals/ folder
                    cb.value = 'orb_key:' + keyAttr;
                }
                cb.checked = document.getElementById('export-select-all').checked;
                cb.id = 'export_cb_' + Math.random().toString(36).substring(7);
                
                const cl = document.createElement('label');
                cl.htmlFor = cb.id;
                cl.innerText = text;
                cl.style.marginLeft = '8px';
                cl.style.fontSize = '12px';
                
                div.appendChild(cb);
                div.appendChild(cl);
                listDiv.appendChild(div);
            });
            
            // Handle ESP separately if present
            const espCb = sidebar.querySelector('input[data-esp]');
            if (espCb) {
                const text = espCb.closest('label').innerText.trim();
                const div = document.createElement('div');
                div.style.marginBottom = '5px';
                
                const cb = document.createElement('input');
                cb.type = 'checkbox';
                cb.value = `molecules/${name}_esp.json`;
                cb.checked = document.getElementById('export-select-all').checked;
                cb.id = 'export_cb_esp';
                
                const cl = document.createElement('label');
                cl.htmlFor = cb.id;
                cl.innerText = text || "ESP Surface";
                cl.style.marginLeft = '8px';
                cl.style.fontSize = '12px';
                
                div.appendChild(cb);
                div.appendChild(cl);
                listDiv.appendChild(div);
            }
        }
        
        if (listDiv.children.length === 0) {
            listDiv.innerHTML = '<div style="color:#888; font-size:12px;">No orbitals found to export.</div>';
        }
        
        document.getElementById('export-modal').style.display = 'block';
    }

    async function confirmDownloadBundle() {
        document.getElementById('export-modal').style.display = 'none';
        const name = document.getElementById('job-name').value;
        
        const keepFiles = [];
        const chks = document.getElementById('export-orbital-list').querySelectorAll('input[type="checkbox"]:checked');
        chks.forEach(chk => {
            const fileAttr = chk.value;
            if (fileAttr) {
                if (fileAttr.startsWith('orb_key:')) {
                    const key = fileAttr.replace('orb_key:', '');
                    // For trajectory orbitals, they are stored per-frame in the bundle, e.g. orbitals/name_mep_key_pos.json
                    // Wait, we don't know the exact filenames here for trajectories, but bundle_exporter put them in orbitals/
                    // Let's just pass the key prefix and handle it in server.py, or we can just pass the generic pattern.
                    keepFiles.push('trajectory_key:' + key);
                } else if (fileAttr.endsWith('.json') && !fileAttr.endsWith('_esp.json')) {
                    keepFiles.push(fileAttr.replace('.json', '_pos.json'));
                    keepFiles.push(fileAttr.replace('.json', '_neg.json'));
                    keepFiles.push(fileAttr);
                } else {
                    keepFiles.push(fileAttr);
                }
            }
        });
        
        try {
            status('Packaging bundle for download...');
            const res = await fetch('/api/bundle/package', {
                method: 'POST',
                headers: {'Content-Type': 'application/json'},
                body: JSON.stringify({name: name, keep_files: keepFiles})
            });
            const data = await res.json();
            if (data.download_url) {
                window.location.href = data.download_url;
                status('Download started.');
            } else {
                status('Error packaging bundle.');
            }
        } catch(e) {
            console.error('Failed to package bundle', e);
            status('Error packaging bundle: ' + e.message);
        }
    }

    function init() {
        initResizer();
        viewer = $3Dmol.createViewer('viewer-3d', {
            backgroundColor: '#0f0f1a',
        });
        status('Ready — draw a reaction, then click 📦 Compute Script');

        refreshDirs();
    }

    function refreshDirs() {
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
    }

    // ── SMILES helpers ──────────────────────────────────────────────

    async function getSmiles() {
        try {
            const frame = document.getElementById('ketcher-frame');
            const ketcher = frame.contentWindow.ketcher;
            return await ketcher.getSmiles();
        } catch (e) {
            status('Error reading SMILES from Ketcher');
            return null;
        }
    }

    // ── Compute Script ──────────────────────────────────────────────

    let pollInterval = null;
    let currentJobName = null;

    async function stopCompute() {
        if (!currentJobName) return;
        await fetch('/api/job/stop', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ name: currentJobName }),
        });
        if (pollInterval) clearInterval(pollInterval);
        document.getElementById('btn-stop').style.display = 'none';
        status('Computation stopped.');
    }

    async function generateComputeScript() {
        const name = document.getElementById('job-name').value || 'molecule';
        const charge = parseInt(document.getElementById('mol-charge').value) || 0;
        const engineSelect = document.getElementById('qm-engine');
        const engine = engineSelect ? engineSelect.value : 'xtb';
        const smiles = await getSmiles();
        currentJobName = name;

        if (!smiles) {
            status('Error: Draw a molecule in Ketcher first');
            return;
        }

        status('Generating compute script...');
        try {
            const frame = document.getElementById('ketcher-frame');
            const ketcher = frame.contentWindow.ketcher;
            let rxn = '';
            let svg = '';
            let ket = '';
            try {
                const mol = await ketcher.getMolfile();
                svg = await ketcher.generateImage(mol, { outputFormat: 'svg' });
                ket = await ketcher.getKet();
            } catch (err) {
                console.warn("Could not get MOL/SVG/KET from Ketcher", err);
            }

            const res = await fetch('/api/rxn/compute-script', {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify({ smiles, name, charge, engine, rxn, svg, ket }),
            });
            const data = await res.json();
            if (data.error) { 
                status(`Error: ${data.error}`);
                return; 
            }

            status('Computation started...');
            document.getElementById('btn-stop').style.display = 'inline-block';
            openLog();
            
            // Poll for status
            if (pollInterval) clearInterval(pollInterval);
            pollInterval = setInterval(async () => {
                const sres = await fetch(`/api/job/status?name=${name}`);
                const sdata = await sres.json();
                
                clearLog();
                appendLog('═══ Computation Log ═══\n');
                appendLog(sdata.log || 'Starting up...');
                
                if (!sdata.running) {
                    clearInterval(pollInterval);
                    document.getElementById('btn-stop').style.display = 'none';
                    if (sdata.ready) {
                        status('✓ Computation finished!');
                        await refreshDirs();
                        const chooser = document.getElementById('run-chooser');
                        if (chooser) chooser.value = name;
                        const jobNameInput = document.getElementById('job-name');
                        if (jobNameInput) jobNameInput.value = name;
                        await verifyBundle();
                    } else {
                        status('Computation failed (bundle not ready).');
                    }
                }
            }, 1000);

        } catch (e) {
            status(`Error: ${e.message}`);
        }
    }
    // ── Verify Bundle ───────────────────────────────────────────────

    async function verifyBundle() {
        const name = document.getElementById('job-name').value ||
                     document.getElementById('run-chooser').value;
        if (!name) {
            status('Error: Enter a job name or select a job from Load Job');
            return;
        }

        status(`Loading bundle: ${name}...`);
        try {
            // Add cache buster to bypass stale browser cache on old bundles
            const res = await fetch(`/api/bundle/${name}/manifest.json?t=${Date.now()}`);
            if (!res.ok) {
                status(`Error: No bundle found for "${name}". Run the compute script first.`);
                return;
            }
            const manifest = await res.json();

            // Clear viewer
            viewer.removeAllModels();
            viewer.removeAllShapes();
            viewer.removeAllSurfaces();

            // Load SMILES back into Ketcher
            if (manifest.reaction_smiles) {
                try {
                    const frame = document.getElementById('ketcher-frame');
                    if (frame && frame.contentWindow && frame.contentWindow.ketcher) {
                        await frame.contentWindow.ketcher.setMolecule(manifest.reaction_smiles);
                    }
                } catch (err) {
                    console.warn("Could not load SMILES into Ketcher:", err);
                }
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

            if (manifest.has_trajectory) {
                const opt = document.createElement('option');
                opt.value = `bundle:${name}:trajectory`;
                opt.textContent = `▶ IRC Trajectory (${manifest.n_frames} frames)`;
                chooser.appendChild(opt);
            }

            // Summary in log
            appendLog('\n\n\n' + '═'.repeat(50) + '\n');
            appendLog(`Job Result: ${manifest.title || name}\n`);
            appendLog('═'.repeat(50));
            appendLog(`SMILES: ${manifest.smiles}`);
            appendLog(`Method: ${manifest.method}`);
            appendLog(`Molecules: ${manifest.molecules.length}`);
            if (manifest.barrier_forward_kcal != null)
                appendLog(`Barrier: ${manifest.barrier_forward_kcal.toFixed(1)} kcal/mol`);
            appendLog(`Trajectory: ${manifest.has_trajectory ? manifest.n_frames + ' frames' : 'none'}`);
            appendLog(`Orbitals: ${manifest.has_orbitals ? manifest.n_orbital_cubes + ' cubes' : 'none'}`);
            manifest.molecules.forEach(m => {
                const nOrb = m.orbitals ? m.orbitals.length : 0;
                appendLog(`  ${m.role} ${m.key}: ${m.smiles} (${nOrb} orbitals)`);
            });
            appendLog('═'.repeat(50));
            appendLog('Select items from the Chemical dropdown to verify 3D/orbitals.');
            openLog();

            // Auto-load first molecule
            loadBundleItem(chooser.value);

            status(`✓ Bundle loaded: ${manifest.molecules.length} molecules, ${manifest.n_frames || 0} frames`);
        } catch(e) {
            status(`Error loading bundle: ${e.message}`);
        }
    }

    // ── Bundle Item Loading ─────────────────────────────────────────

    async function loadBundleItem(value) {
        if (!value || !value.startsWith('bundle:')) return;
        const parts = value.split(':');
        const bundleName = parts[1];
        const type = parts[2];
        const base = `/api/bundle/${bundleName}`;
        const cacheBuster = `?t=${Date.now()}`;

        // Clear bundle trajectory state when switching to molecule
        if (type === 'mol') window._bundleTraj = null;

        if (type === 'mol') {
            const key = parts[3];
            const mol = await (await fetch(`${base}/molecules/${key}.json${cacheBuster}`)).json();

            viewer.removeAllModels();
            viewer.removeAllShapes();
            viewer.removeAllSurfaces();

            if (mol.xyz) {
                viewer.addModel(mol.xyz, 'xyz');
                viewer.setStyle({}, {stick:{radius:0.12}, sphere:{scale:0.25}});
                viewer.zoomTo();
                viewer.render();
            }

            // Build orbital toggles for hierarchical groups
            const sidebar = document.getElementById('orbital-toggles');
            sidebar.innerHTML = '';
            

            if (mol.orbitals && mol.orbitals.length) {
                // Flatten items for the change handler
                const allItems = mol.orbitals.reduce((acc, group) => acc.concat(group.items), []);
                const cache = {};
                let currentRenderId = 0;

                window.onChangeHandler = async () => {
                    const renderId = ++currentRenderId;
                    
                    const shapesToAdd = [];
                    // Iterate checked checkboxes directly (not allItems) to avoid
                    // querySelector collisions when multiple items share the same file
                    const checkedBoxes = sidebar.querySelectorAll('input[data-orb-file]:checked');
                    for (const chk of checkedBoxes) {
                        const file = chk.getAttribute('data-orb-file');
                        const color = chk.getAttribute('data-orb-color');
                        if (!file) continue;
                        
                        if (file.endsWith('.json')) {
                            for (let sign of ['pos', 'neg']) {
                                try {
                                    const actualFile = file.replace('.json', `_${sign}.json`);
                                    let mesh = cache[actualFile];
                                    if (!mesh) {
                                        mesh = await (await fetch(`${base}/${actualFile}`)).json();
                                        cache[actualFile] = mesh;
                                    }
                                    if (mesh && mesh.vertices) {
                                        const finalColor = (sign === 'neg') ? shiftColor(color) : color;
                                        const alpha = (sign === 'neg') ? 0.5 : 0.7;
                                        shapesToAdd.push({vertexArr: mesh.vertices, faceArr: mesh.faces, color: $3Dmol.CC.color(finalColor), opacity: alpha});
                                    }
                                } catch(e) {}
                            }
                        } else if (file.endsWith('.cube')) {
                            let txt = cache[file];
                            if (!txt) {
                                try {
                                    const r = await fetch(`${base}/${file}`);
                                    const buf = await r.arrayBuffer();
                                    txt = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'))).text();
                                    cache[file] = txt;
                                } catch(e) { continue; }
                            }
                            const vol = new $3Dmol.VolumeData(txt, 'cube');
                            shapesToAdd.push({vol: vol, isoval: 0.03, color: color, opacity: 0.55, smoothness: 1});
                            shapesToAdd.push({vol: vol, isoval: -0.03, color: lighten(color), opacity: 0.35, smoothness: 1});
                        }
                    }

                    if (renderId !== currentRenderId) return; // Abort if a newer render started
                    
                    viewer.removeAllModels(); viewer.removeAllShapes();
                    viewer.addModel(mol.xyz, 'xyz');
                    viewer.setStyle({}, {stick:{radius:0.12}, sphere:{scale:0.25}});
                    
                    for (const shape of shapesToAdd) {
                        if (shape.vol) {
                            viewer.addIsosurface(shape.vol, {isoval: shape.isoval, color: shape.color, opacity: shape.opacity, smoothness: shape.smoothness});
                        } else {
                            viewer.addCustom(shape);
                        }
                    }
                    viewer.render();
                };

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
                        cb.onchange = window.onChangeHandler;
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

                let espCache = null;
                espCb.onchange = async () => {
                    if (espCb.checked) {
                        if (!espCache) {
                            try {
                                espCache = await (await fetch(`${base}/${mol.esp_surface.file}`)).json();
                            } catch(e) { console.error('ESP load failed', e); return; }
                        }
                        if (espCache && espCache.vertices) {
                            // Color each vertex by ESP value using red↔blue diverging colormap
                            const espMin = espCache.esp_min;
                            const espMax = espCache.esp_max;
                            const absMax = Math.max(Math.abs(espMin), Math.abs(espMax), 20); // at least ±20 kcal/mol
                            const colors = espCache.esp_values.map(v => {
                                const t = Math.max(-1, Math.min(1, v / absMax));
                                let r, g, b;
                                if (t < 0) {
                                    r = 1.0; g = 1.0 + t; b = 1.0 + t;
                                } else {
                                    r = 1.0 - t; g = 1.0 - t; b = 1.0;
                                }
                                return {r: r, g: g, b: b};
                            });

                            // Build per-face color array for addCustom
                            const verts = espCache.vertices;
                            const faces = espCache.faces;
                            const normals = espCache.normals || [];
                            // addCustom with per-vertex coloring
                            console.log("Colors:", colors[0]);
                            viewer.addCustom({
                                vertexArr: verts,
                                faceArr: faces,
                                normalArr: normals.length ? normals : undefined,
                                                                color: colors.length ? colors : $3Dmol.CC.color('#aaaaaa'),
                                opacity: 0.85,
                            });
                            viewer.render();
                        }
                    } else {
                        // Remove ESP by re-triggering orbital render
                        if (window.onChangeHandler) window.onChangeHandler();
                    }
                };
            }

            document.getElementById('viewer-info').textContent = `${mol.smiles} — ${mol.role}`;

        } else if (type === 'trajectory') {
            const frames = await (await fetch(`${base}/trajectory/frames.json`)).json();
            let energy = null;
            try { energy = await (await fetch(`${base}/trajectory/energy.json`)).json(); } catch(e) {}
            let orbData = null;
            try { orbData = await (await fetch(`${base}/trajectory/orbitals.json`)).json(); } catch(e) {}

            window._bundleTraj = { frames, energy, orbData, base, orbCache: {}, orbActive: new Set() };

            // Setup slider
            const controls = document.getElementById('playback-controls');
            controls.style.display = 'flex';
            const slider = document.getElementById('frame-slider');
            slider.max = frames.length - 1;
            slider.value = 0;

            // Build orbital toggles
            const sidebar = document.getElementById('orbital-toggles');
            sidebar.innerHTML = '';
            if (orbData && orbData.ibos) {
                orbData.ibos.forEach((ibo, i) => {
                    if (i >= orbData.orbital_keys.length) return;
                    const key = orbData.orbital_keys[i];
                    const row = document.createElement('label');
                    row.className = 'orbital-toggle';
                    row.style.cssText = 'display:flex;align-items:center;gap:6px;padding:4px 8px;cursor:pointer;';
                    const cb = document.createElement('input');
                    cb.type = 'checkbox';
                    cb.checked = ibo.significant;
                    cb.setAttribute('data-orb-key', key);
                    cb.onchange = () => {
                        const t = window._bundleTraj;
                        if (cb.checked) t.orbActive.add(key);
                        else t.orbActive.delete(key);
                        renderBundleFrame(parseInt(slider.value));
                    };
                    if (ibo.significant) window._bundleTraj.orbActive.add(key);
                    const dot = document.createElement('span');
                    dot.style.cssText = `display:inline-block;width:10px;height:10px;border-radius:50%;background:${ibo.color || '#aaa'};`;
                    const lbl = document.createElement('span');
                    lbl.style.cssText = 'font-size:12px;color:#ccc;';
                    lbl.textContent = `${ibo.significant ? '★ ' : ''}${ibo.label}`;
                    row.appendChild(cb); row.appendChild(dot); row.appendChild(lbl);
                    sidebar.appendChild(row);
                });
            }

            document.getElementById('viewer-info').textContent = 'IRC Trajectory';
            renderBundleFrame(0);
        }
    }

    // ── Bundle Frame Rendering ──────────────────────────────────────
    let trajRenderId = 0;

    async function renderBundleFrame(idx) {
        const t = window._bundleTraj;
        if (!t) return;
        const fr = t.frames[idx];
        if (!fr) return;
        
        const currentRenderId = ++trajRenderId;

        const n = fr.atoms.length;
        let xyz = n + '\n' + (fr.comment || '') + '\n';
        fr.atoms.forEach(a => { xyz += a.sym + ' ' + a.x + ' ' + a.y + ' ' + a.z + '\n'; });

        const shapesToAdd = [];
        
        // Fetch and prepare all shapes
        if (t.orbData && t.orbActive.size > 0) {
            const frame = t.orbData.frames[idx];
            if (frame) {
                for (const key of t.orbActive) {
                    const info = frame[key];
                    if (!info) continue;
                    const ii = t.orbData.orbital_keys.indexOf(key);
                    const col = (ii >= 0 && t.orbData.ibos[ii]) ? t.orbData.ibos[ii].color : '#aaa';
                    try {
                        if (info.file.endsWith('.json')) {
                            for (let sign of ['pos', 'neg']) {
                                try {
                                    const actualFile = info.file.replace('.json', `_${sign}.json`);
                                    let mesh = t.orbCache[actualFile];
                                    if (!mesh) {
                                        mesh = await (await fetch(`${t.base}/${actualFile}`)).json();
                                        t.orbCache[actualFile] = mesh;
                                    }
                                    if (mesh && mesh.vertices) {
                                        const finalColor = (sign === 'neg') ? shiftColor(col) : col;
                                        const alpha = (sign === 'neg') ? 0.5 : 0.7;
                                        shapesToAdd.push({vertexArr: mesh.vertices, faceArr: mesh.faces, color: $3Dmol.CC.color(finalColor), opacity: alpha});
                                    }
                                } catch(e) {}
                            }
                        } else if (info.file.endsWith('.cube')) {
                            let txt = t.orbCache[info.file];
                            if (!txt) {
                                try {
                                    const r = await fetch(`${t.base}/${info.file}`);
                                    const buf = await r.arrayBuffer();
                                    txt = await new Response(new Blob([buf]).stream().pipeThrough(new DecompressionStream('gzip'))).text();
                                    t.orbCache[info.file] = txt;
                                } catch(e) { continue; }
                            }
                            const vol = new $3Dmol.VolumeData(txt, 'cube');
                            shapesToAdd.push({vol: vol, isoval: 0.03, color: col, opacity: 0.55, smoothness: 1});
                            shapesToAdd.push({vol: vol, isoval: -0.03, color: lighten(col), opacity: 0.35, smoothness: 1});
                        }
                    } catch(e) {}
                }
            }
        }
        
        if (currentRenderId !== trajRenderId) return; // Abort if a newer render started
        
        viewer.removeAllModels(); viewer.removeAllShapes();
        viewer.addModel(xyz, 'xyz');
        viewer.setStyle({}, {stick:{radius:0.12}, sphere:{scale:0.25}});
        if (idx === 0) viewer.zoomTo();
        
        for (const shape of shapesToAdd) {
            if (shape.vol) {
                viewer.addIsosurface(shape.vol, {isoval: shape.isoval, color: shape.color, opacity: shape.opacity, smoothness: shape.smoothness});
            } else {
                viewer.addCustom(shape);
            }
        }
        
        viewer.render();

        // Update counter
        let info = `${idx}/${t.frames.length - 1}`;
        if (t.energy && t.energy.frames && t.energy.frames[idx]) {
            const e = t.energy.frames[idx].energy_kcal;
            if (e != null) info += ` · ${e.toFixed(1)} kcal/mol`;
        }
        document.getElementById('frame-counter').textContent = info;
    }

    // ── Frame Slider ────────────────────────────────────────────────

    function setFrame(idx) {
        idx = parseInt(idx);
        if (window._bundleTraj) {
            renderBundleFrame(idx);
        }
    }

    // ── Load Chemical (dropdown handler) ─────────────────────────────

    function loadChemical() {
        const chooser = document.getElementById('chemical-chooser');
        const val = chooser.value;
        if (!val) return;
        if (val.startsWith('bundle:')) {
            loadBundleItem(val);
        }
    }

    // ── Load Job (dropdown) ──────────────────────────────────────────

    function loadRun() {
        const dir = document.getElementById('run-chooser').value;
        if (!dir) return;
        document.getElementById('job-name').value = dir;
        status(`Job "${dir}" selected. Loading...`);
        verifyBundle();
    }


    function shiftColor(hex) {
        if (!hex || hex.length < 7) return '#cccccc';
        const r = parseInt(hex.slice(1, 3), 16);
        const g = parseInt(hex.slice(3, 5), 16);
        const b = parseInt(hex.slice(5, 7), 16);
        const nr = Math.min(255, 255 - r + 40);
        const ng = Math.min(255, 255 - g + 40);
        const nb = Math.min(255, 255 - b + 40);
        return `#${nr.toString(16).padStart(2, '0')}${ng.toString(16).padStart(2, '0')}${nb.toString(16).padStart(2, '0')}`;
    }

    // ── Log Panel ────────────────────────────────────────────────────

    function toggleLog() {
        const d = document.getElementById('log-drawer');
        d.classList.toggle('open');
    }
    function openLog() {
        document.getElementById('log-drawer').classList.add('open');
    }
    function clearLog() {
        document.getElementById('log-content').textContent = '';
    }
    function appendLog(text) {
        const el = document.getElementById('log-content');
        el.textContent += text + '\n';
        el.scrollTop = el.scrollHeight;
    }

    function status(msg) {
        document.getElementById('status-bar').textContent = msg;
    }

    // ── Boot ─────────────────────────────────────────────────────────

    // ── Editor ─────────────────────────────────────────────────────────

    function openEditor() {
        document.getElementById('editor-modal').style.display = 'flex';
        switchEditorTab('trim');
        
        // Populate merge dropdown
        const chooser = document.getElementById('run-chooser');
        
        for (let i = 1; i < chooser.options.length; i++) {
            const val = chooser.options[i].value;
            const opt = document.createElement('option');
            opt.value = val;
            opt.text = val;
            mergeSelect.appendChild(opt);
        }
        
        // Populate trim range based on current bundle
        if (state.currentRun && state.manifest && state.manifest.n_frames) {
            document.getElementById('trim-end').value = state.manifest.n_frames - 1;
        }
    }

    function switchEditorTab(tab) {
        document.getElementById('panel-trim').style.display = (tab === 'trim') ? 'block' : 'none';
        document.getElementById('panel-merge').style.display = (tab === 'merge') ? 'block' : 'none';
        document.getElementById('tab-trim').style.background = (tab === 'trim') ? '#4CAF50' : '#333355';
        document.getElementById('tab-merge').style.background = (tab === 'merge') ? '#2196F3' : '#333355';
        document.getElementById('editor-status').textContent = '';
    }

    async function submitTrim() {
        const name = state.currentRun;
        if (!name) {
            document.getElementById('editor-status').textContent = 'Error: Load a bundle first.';
            return;
        }
        const start = parseInt(document.getElementById('trim-start').value);
        const end = parseInt(document.getElementById('trim-end').value);
        document.getElementById('editor-status').textContent = 'Processing trim...';
        try {
            const res = await fetch('/api/editor/trim', {
                method: 'POST',
                body: JSON.stringify({ name: name, start: start, end: end })
            });
            const data = await res.json();
            if (data.error) throw new Error(data.error);
            refreshDirs();
            document.getElementById('editor-status').textContent = `Success! Saved as ${data.new_bundle}. You can load it from the dropdown.`;
        } catch(e) {
            document.getElementById('editor-status').textContent = `Error: ${e.message}`;
        }
    }

    async function submitMerge() {
        const sel = document.getElementById('merge-select');
        const bundles = Array.from(sel.selectedOptions).map(o => o.value);
        if (bundles.length < 2) {
            document.getElementById('editor-status').textContent = 'Error: Select at least 2 bundles.';
            return;
        }
        document.getElementById('editor-status').textContent = 'Processing merge...';
        try {
            const res = await fetch('/api/editor/merge', {
                method: 'POST',
                body: JSON.stringify({ bundles: bundles })
            });
            const data = await res.json();
            if (data.error) throw new Error(data.error);
            refreshDirs();
            
            document.getElementById('editor-status').textContent = `Success! Saved as ${data.new_bundle}. Loading...`;
            document.getElementById('job-name').value = data.new_bundle;
            setTimeout(() => {
                document.getElementById('editor-modal').style.display = 'none';
                verifyBundle();
            }, 1000);
            
        } catch(e) {
            document.getElementById('editor-status').textContent = `Error: ${e.message}`;
        }
    }

    window.addEventListener('DOMContentLoaded', init);

    return {
        generateComputeScript, verifyBundle, stopCompute,
        loadChemical, loadRun, toggleLog,
        setFrame,
        openEditor, switchEditorTab, submitTrim, submitMerge, downloadCustomBundle,
        toggleExportSelectAll, confirmDownloadBundle
    };
})();
