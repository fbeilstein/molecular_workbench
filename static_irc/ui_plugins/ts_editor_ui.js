import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';

let currentDockedMolblocks = [];
let animationId = null;
let keyHandlerInstalled = false;
let activeTransformControl = null;

// Simple CPK color mapping
const elementColors = {
    'C': 0x909090, 'H': 0xFFFFFF, 'O': 0xFF0D0D, 'N': 0x2250FF, 'S': 0xFFFF30,
    'P': 0xFF8000, 'F': 0x90E050, 'Cl': 0x1FF01F, 'Br': 0xA62929, 'I': 0x940094,
    'Mg': 0x8AFF00, 'Na': 0xAB5CF2, 'K': 0x8F40D4, 'Li': 0xCC80FF, 'B': 0xFFB5B5, 'Si': 0xF0C8A0
};
const elementRadii = {
    'H': 0.3, 'C': 0.7, 'O': 0.6, 'N': 0.65, 'Br': 0.9, 'Cl': 0.8, 'I': 1.0, 'S': 0.8, 'Mg': 0.8
};

function parseSDF(sdfString) {
    let lines = sdfString.split('\n');
    if (lines.length < 4) return null;
    let countsLine = lines[3];
    let numAtoms = parseInt(countsLine.substring(0, 3));
    let numBonds = parseInt(countsLine.substring(3, 6));

    let atoms = [];
    for (let i = 0; i < numAtoms; i++) {
        let line = lines[4 + i];
        let x = parseFloat(line.substring(0, 10));
        let y = parseFloat(line.substring(10, 20));
        let z = parseFloat(line.substring(20, 30));
        let symbol = line.substring(31, 34).trim();
        atoms.push({ x, y, z, symbol });
    }

    let bonds = [];
    for (let i = 0; i < numBonds; i++) {
        let line = lines[4 + numAtoms + i];
        let a1 = parseInt(line.substring(0, 3)) - 1;
        let a2 = parseInt(line.substring(3, 6)) - 1;
        let order = parseInt(line.substring(6, 9));
        bonds.push({ a1, a2, order });
    }

    return { atoms, bonds, originalLines: lines, numAtoms, numBonds };
}

function updateSDFCoordinates(parsedSDF, group) {
    // Transform the (centroid-centred) atom coordinates by the group's world matrix
    let matrix = group.matrixWorld;
    let newLines = [...parsedSDF.originalLines];

    for (let i = 0; i < parsedSDF.numAtoms; i++) {
        let atom = parsedSDF.atoms[i];
        let vec = new THREE.Vector3(atom.x, atom.y, atom.z);
        vec.applyMatrix4(matrix);

        let lineIdx = 4 + i;
        let line = newLines[lineIdx];
        let newXStr = vec.x.toFixed(4).padStart(10, ' ');
        let newYStr = vec.y.toFixed(4).padStart(10, ' ');
        let newZStr = vec.z.toFixed(4).padStart(10, ' ');
        newLines[lineIdx] = newXStr + newYStr + newZStr + line.substring(30);
    }
    return newLines.join('\n');
}

function el(tag, style = {}, text = null) {
    const e = document.createElement(tag);
    Object.assign(e.style, style);
    if (text !== null) e.textContent = text;
    return e;
}

function styledButton(text, primary = false) {
    const b = el('button', {
        padding: '9px 14px', margin: '4px 6px 4px 0', borderRadius: '6px', cursor: 'pointer',
        border: primary ? 'none' : '1px solid #c7ccd6', fontWeight: primary ? '600' : '500',
        background: primary ? 'linear-gradient(135deg,#4f6ef7,#7b4ff7)' : '#fff',
        color: primary ? '#fff' : '#333', boxShadow: primary ? '0 2px 8px rgba(79,110,247,.35)' : 'none'
    }, text);
    return b;
}

function drawProfile(canvas, profile) {
    const ctx = canvas.getContext('2d');
    const W = canvas.width, H = canvas.height, pad = 30;
    ctx.clearRect(0, 0, W, H);
    if (!profile || profile.length < 2) return;
    const xs = profile.map(p => p[0]), ys = profile.map(p => p[1]);
    const ymin = Math.min(0, ...ys), ymax = Math.max(...ys) + 1;
    const X = x => pad + (x - Math.min(...xs)) / (Math.max(...xs) - Math.min(...xs) || 1) * (W - pad - 10);
    const Y = y => H - pad + 5 - (y - ymin) / (ymax - ymin || 1) * (H - pad - 15);
    ctx.strokeStyle = '#bbb'; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(pad, 8); ctx.lineTo(pad, H - pad + 5); ctx.lineTo(W - 5, H - pad + 5); ctx.stroke();
    ctx.fillStyle = '#666'; ctx.font = '11px sans-serif';
    ctx.fillText('ΔE kcal/mol', pad + 4, 14);
    ctx.fillText('λ (R → P)', W - 70, H - 8);
    ctx.fillText(ymax.toFixed(0), 2, Y(ymax) + 4);
    ctx.fillText('0', 14, Y(0) + 4);
    ctx.strokeStyle = '#7b4ff7'; ctx.lineWidth = 2; ctx.beginPath();
    profile.forEach((p, i) => i ? ctx.lineTo(X(p[0]), Y(p[1])) : ctx.moveTo(X(p[0]), Y(p[1])));
    ctx.stroke();
    let iMax = ys.indexOf(Math.max(...ys));
    profile.forEach((p, i) => {
        ctx.fillStyle = i === iMax ? '#e0245e' : '#4f6ef7';
        ctx.beginPath(); ctx.arc(X(p[0]), Y(p[1]), i === iMax ? 5 : 3, 0, 2 * Math.PI); ctx.fill();
    });
}

export const TSEditorUI = {
    render: (container, ui, state, stepId) => {
        // NB: app.js already prints ui.message above us
        currentDockedMolblocks = [...ui.reactants_3d];
        if (animationId) cancelAnimationFrame(animationId);

        const info = ui.reaction_info || { forming: [], breaking: [], summary: [] };
        let dockVariant = 0;
        // Editable copy of the reacting bonds (mapping-derived unless the user overrides)
        const clonePairs = (p) => (p || []).map(([[a, b], [c, d]]) => [[a, b], [c, d]]);
        let rxnBonds = { forming: clonePairs(info.forming), breaking: clonePairs(info.breaking) };
        let overridden = false;
        let pickMode = false;
        let picked = [];   // [{fi, ai}] max 2

        const editorLayout = el('div', { display: 'flex', gap: '20px', marginTop: '12px', flexWrap: 'wrap' });
        const viewerContainer = el('div', {
            width: '700px', height: '600px', position: 'relative', border: '1px solid #d5d9e2',
            borderRadius: '10px', overflow: 'hidden', boxShadow: '0 4px 18px rgba(0,0,0,.06)'
        });
        editorLayout.appendChild(viewerContainer);
        const controlsContainer = el('div', { flex: '1', minWidth: '340px', display: 'flex', flexDirection: 'column', gap: '12px' });
        editorLayout.appendChild(controlsContainer);
        container.appendChild(editorLayout);

        // --- THREE.JS SETUP ---
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0xf5f6fa);
        const camera = new THREE.PerspectiveCamera(45, 700 / 600, 0.1, 200);
        camera.position.set(0, 0, 15);
        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(700, 600);
        viewerContainer.appendChild(renderer.domElement);
        scene.add(new THREE.AmbientLight(0x404040, 2));
        const dirLight = new THREE.DirectionalLight(0xffffff, 1.5);
        dirLight.position.set(10, 10, 10);
        scene.add(dirLight);

        const orbit = new OrbitControls(camera, renderer.domElement);
        orbit.enableDamping = true;
        const transformControl = new TransformControls(camera, renderer.domElement);
        transformControl.addEventListener('dragging-changed', e => { orbit.enabled = !e.value; });
        scene.add(transformControl);
        activeTransformControl = transformControl;

        let moleculeGroups = [];
        let parsedSDFs = [];
        let reactionLines = [];   // dashed forming-bond indicators

        const sphereGeo = new THREE.SphereGeometry(1, 32, 32);
        const cylGeo = new THREE.CylinderGeometry(0.15, 0.15, 1, 16);

        const samePair = ([[f1, x1], [f2, x2]], fa, a, fb, b) =>
            (f1 === fa && x1 === a && f2 === fb && x2 === b) || (f1 === fb && x1 === b && f2 === fa && x2 === a);

        function isBreaking(fi, a1, a2) {
            return rxnBonds.breaking.some(p => samePair(p, fi, a1, fi, a2));
        }

        function buildMolecules(molblocks) {
            if (transformControl.object) transformControl.detach();
            moleculeGroups.forEach(g => scene.remove(g));
            reactionLines.forEach(l => scene.remove(l.line));
            moleculeGroups = []; parsedSDFs = []; reactionLines = [];

            molblocks.forEach((molBlock, fi) => {
                const parsed = molBlock ? parseSDF(molBlock) : null;
                if (!parsed) return;
                parsedSDFs.push(parsed);
                const group = new THREE.Group();

                // Centre on the centroid so the gizmo sits on the molecule
                let cx = 0, cy = 0, cz = 0;
                parsed.atoms.forEach(a => { cx += a.x; cy += a.y; cz += a.z; });
                const n = parsed.atoms.length || 1;
                cx /= n; cy /= n; cz /= n;
                parsed.atoms.forEach(a => { a.x -= cx; a.y -= cy; a.z -= cz; });
                group.position.set(cx, cy, cz);

                group.userData.atomMeshes = parsed.atoms.map((atom, ai) => {
                    const mat = new THREE.MeshPhongMaterial({ color: elementColors[atom.symbol] || 0xcccccc, shininess: 50 });
                    const mesh = new THREE.Mesh(sphereGeo, mat);
                    const r = elementRadii[atom.symbol] || 0.5;
                    mesh.scale.set(r, r, r);
                    mesh.position.set(atom.x, atom.y, atom.z);
                    mesh.userData = { fi, ai, isAtom: true };
                    if (picked.some(p => p.fi === fi && p.ai === ai)) mat.emissive.setHex(0x3355ff);
                    group.add(mesh);
                    return mesh;
                });

                parsed.bonds.forEach(bond => {
                    const a1 = parsed.atoms[bond.a1], a2 = parsed.atoms[bond.a2];
                    const p1 = new THREE.Vector3(a1.x, a1.y, a1.z);
                    const p2 = new THREE.Vector3(a2.x, a2.y, a2.z);
                    const breaking = isBreaking(fi, bond.a1, bond.a2);
                    const mat = new THREE.MeshPhongMaterial({ color: breaking ? 0xe0245e : 0x999999 });
                    const mesh = new THREE.Mesh(cylGeo, mat);
                    mesh.position.copy(p1).lerp(p2, 0.5);
                    mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), p2.clone().sub(p1).normalize());
                    mesh.scale.set(breaking ? 1.3 : 1, p1.distanceTo(p2), breaking ? 1.3 : 1);
                    group.add(mesh);
                });

                scene.add(group);
                moleculeGroups.push(group);
            });

            // Reacting bonds as dashed lines: forming = green, breaking = red
            const addLine = ([[f1, a1], [f2, a2]], color) => {
                const g1 = moleculeGroups[f1], g2 = moleculeGroups[f2];
                if (!g1 || !g2) return;
                const geom = new THREE.BufferGeometry().setFromPoints([new THREE.Vector3(), new THREE.Vector3()]);
                const mat = new THREE.LineDashedMaterial({ color, dashSize: 0.18, gapSize: 0.12, linewidth: 2 });
                const line = new THREE.Line(geom, mat);
                scene.add(line);
                reactionLines.push({ line, m1: g1.userData.atomMeshes[a1], m2: g2.userData.atomMeshes[a2] });
            };
            rxnBonds.forming.forEach(p => addLine(p, 0x18a558));
            rxnBonds.breaking.forEach(p => addLine(p, 0xe0245e));
        }

        function updateReactionLines() {
            const v1 = new THREE.Vector3(), v2 = new THREE.Vector3();
            reactionLines.forEach(({ line, m1, m2 }) => {
                if (!m1 || !m2) return;
                m1.getWorldPosition(v1); m2.getWorldPosition(v2);
                const pos = line.geometry.attributes.position;
                pos.setXYZ(0, v1.x, v1.y, v1.z); pos.setXYZ(1, v2.x, v2.y, v2.z);
                pos.needsUpdate = true;
                line.computeLineDistances();
            });
        }

        // Rebuild in place, keeping whatever the user has moved
        function rebuildKeepingPose() {
            const attachedIdx = moleculeGroups.indexOf(transformControl.object);
            const blocks = TSEditorUI.extractData ? TSEditorUI.extractData() : currentDockedMolblocks;
            currentDockedMolblocks = blocks;
            buildMolecules(blocks);
            if (attachedIdx >= 0 && moleculeGroups[attachedIdx]) transformControl.attach(moleculeGroups[attachedIdx]);
        }

        buildMolecules(ui.reactants_3d);

        // --- UI CONTROLS ---
        const card = (title) => {
            const c = el('div', { background: '#fff', border: '1px solid #e1e4ea', borderRadius: '10px', padding: '12px 14px' });
            if (title) c.appendChild(el('div', { fontWeight: '700', marginBottom: '6px', color: '#2b2f3a' }, title));
            return c;
        };

        // Reaction analysis
        const rxnCard = card('Reaction analysis (from mapping)');
        const pre = el('pre', { margin: 0, fontSize: '12px', whiteSpace: 'pre-wrap', color: '#333' });
        pre.textContent = (info.summary || []).join('\n') || 'No information';
        rxnCard.appendChild(pre);
        controlsContainer.appendChild(rxnCard);

        // Main actions
        const actCard = card('Transition state');
        const optBtn = styledButton('⚡ Find TS (automatic scan)', true);
        const refineBtn = styledButton('🎯 Refine my guess as TS');
        const nextPoseBtn = styledButton('↻ Next docking pose');
        const optStatus = el('div', { fontWeight: '600', color: '#4f6ef7', marginTop: '6px', fontSize: '13px' });
        actCard.appendChild(optBtn);
        actCard.appendChild(refineBtn);
        if (moleculeGroups.length > 1) actCard.appendChild(nextPoseBtn);
        actCard.appendChild(el('div', { fontSize: '11.5px', color: '#777', marginTop: '4px' },
            'Scan: starts from reactant-like distances of the current pose, drives the reacting bonds, ' +
            'takes the energy maximum. Refine: your current geometry is the TS guess (reacting bonds ' +
            'frozen while the rest relaxes, then saddle search) - use after placing fragments by hand.'));
        actCard.appendChild(optStatus);
        controlsContainer.appendChild(actCard);

        // Reacting bonds editor (semi-automatic fallback)
        const bondCard = card('Reacting bonds');
        const bondList = el('pre', { margin: '0 0 6px 0', fontSize: '12px', whiteSpace: 'pre-wrap', color: '#333' });
        const pickBtn = styledButton('Pick atoms: OFF');
        const formBtn = styledButton('+ Forming');
        const breakBtn = styledButton('+ Breaking');
        const removeBtn = styledButton('Remove pair');
        const resetBondsBtn = styledButton('Reset to mapping');
        const pickInfo = el('div', { fontSize: '11.5px', color: '#777' },
            'Turn picking on, click two atoms, then mark them. Overrides the mapping for scan, refine and re-docking.');
        [bondList, pickBtn, formBtn, breakBtn, removeBtn, resetBondsBtn, pickInfo].forEach(x => bondCard.appendChild(x));
        controlsContainer.appendChild(bondCard);

        const atomLabel = (fi, ai) => {
            const s = parsedSDFs[fi]?.atoms[ai]?.symbol || '?';
            return `${s}(R${fi + 1}:${ai})`;
        };
        function refreshBondList() {
            const lines = [];
            rxnBonds.forming.forEach(([[a, b], [c, d]]) => lines.push(`  FORMING  ${atomLabel(a, b)} ... ${atomLabel(c, d)}`));
            rxnBonds.breaking.forEach(([[a, b], [c, d]]) => lines.push(`  BREAKING ${atomLabel(a, b)} -- ${atomLabel(c, d)}`));
            bondList.textContent = (overridden ? 'Source: MANUAL\n' : 'Source: mapping\n') + (lines.join('\n') || '  (none)');
            const sel = picked.map(p => atomLabel(p.fi, p.ai)).join(' + ');
            pickInfo.textContent = picked.length ? `Selected: ${sel}` :
                'Turn picking on, click two atoms, then mark them. Overrides the mapping for scan, refine and re-docking.';
            formBtn.disabled = breakBtn.disabled = removeBtn.disabled = picked.length !== 2;
            resetBondsBtn.disabled = !overridden;
        }
        function editPair(kind) {
            if (picked.length !== 2) return;
            const [p, q] = picked;
            const notThis = pr => !samePair(pr, p.fi, p.ai, q.fi, q.ai);
            rxnBonds.forming = rxnBonds.forming.filter(notThis);
            rxnBonds.breaking = rxnBonds.breaking.filter(notThis);
            if (kind) rxnBonds[kind].push([[p.fi, p.ai], [q.fi, q.ai]]);
            overridden = true;
            picked = [];
            rebuildKeepingPose();
            refreshBondList();
        }
        formBtn.onclick = () => editPair('forming');
        breakBtn.onclick = () => editPair('breaking');
        removeBtn.onclick = () => editPair(null);
        resetBondsBtn.onclick = () => {
            rxnBonds = { forming: clonePairs(info.forming), breaking: clonePairs(info.breaking) };
            overridden = false; picked = [];
            rebuildKeepingPose(); refreshBondList();
        };
        pickBtn.onclick = () => {
            pickMode = !pickMode;
            pickBtn.textContent = pickMode ? 'Pick atoms: ON' : 'Pick atoms: OFF';
            pickBtn.style.background = pickMode ? '#e8edff' : '#fff';
            if (!pickMode) { picked = []; rebuildKeepingPose(); }
            refreshBondList();
        };
        refreshBondList();

        // Atom picking (click without dragging)
        const raycaster = new THREE.Raycaster();
        let downXY = null;
        renderer.domElement.addEventListener('pointerdown', e => { downXY = [e.clientX, e.clientY]; });
        renderer.domElement.addEventListener('pointerup', e => {
            if (!pickMode || !downXY || transformControl.dragging) return;
            if (Math.hypot(e.clientX - downXY[0], e.clientY - downXY[1]) > 4) return;
            const rect = renderer.domElement.getBoundingClientRect();
            const ndc = new THREE.Vector2(((e.clientX - rect.left) / rect.width) * 2 - 1,
                -((e.clientY - rect.top) / rect.height) * 2 + 1);
            raycaster.setFromCamera(ndc, camera);
            const meshes = moleculeGroups.flatMap(g => g.userData.atomMeshes);
            const hit = raycaster.intersectObjects(meshes, false)[0];
            if (!hit) return;
            const { fi, ai } = hit.object.userData;
            const idx = picked.findIndex(p => p.fi === fi && p.ai === ai);
            if (idx >= 0) picked.splice(idx, 1);
            else { picked.push({ fi, ai }); if (picked.length > 2) picked.shift(); }
            moleculeGroups.forEach((g, gi) => g.userData.atomMeshes.forEach((m, mi) =>
                m.material.emissive.setHex(picked.some(p => p.fi === gi && p.ai === mi) ? 0x3355ff : 0x000000)));
            refreshBondList();
        });

        // Optional fine-tuning
        const tuneCard = card('Optional fine-tuning');
        tuneCard.appendChild(el('div', { fontSize: '12px', color: '#666', marginBottom: '6px' },
            'Select a reactant, move/rotate it with the gizmo. Left-drag: rotate view · right-drag: pan · scroll: zoom.'));
        ui.reactants_3d.forEach((_, idx) => {
            const btn = styledButton(`Reactant ${idx + 1}`);
            btn.onclick = () => { if (moleculeGroups[idx]) transformControl.attach(moleculeGroups[idx]); };
            tuneCard.appendChild(btn);
        });
        const transBtn = styledButton('Move (T)');
        transBtn.onclick = () => transformControl.setMode('translate');
        const rotBtn = styledButton('Rotate (R)');
        rotBtn.onclick = () => transformControl.setMode('rotate');
        const hideBtn = styledButton('Hide gizmo');
        hideBtn.onclick = () => transformControl.detach();
        tuneCard.appendChild(el('br'));
        [transBtn, rotBtn, hideBtn].forEach(b => tuneCard.appendChild(b));
        controlsContainer.appendChild(tuneCard);

        // Profile + log
        const profCanvas = el('canvas', { width: '100%', border: '1px solid #e1e4ea', borderRadius: '8px', background: '#fff', display: 'none' });
        profCanvas.width = 380; profCanvas.height = 150;
        controlsContainer.appendChild(profCanvas);

        const logBox = el('textarea', {
            width: '100%', height: '220px', fontFamily: 'monospace', fontSize: '11.5px',
            backgroundColor: '#1e2130', color: '#d6e2ff', border: 'none', borderRadius: '8px', padding: '8px', boxSizing: 'border-box'
        });
        logBox.readOnly = true;
        logBox.value = 'Ready.';
        controlsContainer.appendChild(logBox);

        nextPoseBtn.onclick = async () => {
            dockVariant += 1;
            nextPoseBtn.disabled = true;
            try {
                const res = await fetch('/api/step/step_3_ts/custom', {
                    method: 'POST', headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'auto_dock', variant: dockVariant,
                        reaction_override: overridden ? rxnBonds : null
                    })
                });
                const result = await res.json();
                if (result.success) {
                    currentDockedMolblocks = result.molblocks;
                    buildMolecules(result.molblocks);
                    optStatus.style.color = '#4f6ef7';
                    optStatus.textContent = result.message;
                } else {
                    optStatus.style.color = '#e0245e';
                    optStatus.textContent = result.error;
                }
            } finally { nextPoseBtn.disabled = false; }
        };

        async function fetchLog() {
            try {
                const r = await fetch('/api/log');
                const t = await r.text();
                if (t) { logBox.value = t; logBox.scrollTop = logBox.scrollHeight; }
            } catch (e) { }
        }

        async function runOptimization(mode) {
            const buttons = [optBtn, refineBtn, nextPoseBtn];
            const activeBtn = mode === 'refine' ? refineBtn : optBtn;
            const label = activeBtn.textContent;
            buttons.forEach(b => b.disabled = true);
            activeBtn.textContent = 'Searching…';
            optStatus.style.color = '#4f6ef7';
            optStatus.textContent = mode === 'refine'
                ? 'Constrained relax → Sella → frequency check'
                : 'Bond drive → Sella → frequency check';
            logBox.value = 'Starting GFN2-xTB…\n';
            const logInterval = setInterval(fetchLog, 700);

            const data = TSEditorUI.extractData();
            try {
                const res = await fetch('/api/step/step_3_ts/custom', {
                    method: 'POST', headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'optimize_ts', mode, docked_molblocks: data,
                        reaction_override: overridden ? rxnBonds : null
                    })
                });
                const result = await res.json();
                await fetchLog();
                if (result.success) {
                    optStatus.textContent = result.message;
                    currentDockedMolblocks = result.optimized_molblocks;
                    picked = [];
                    buildMolecules(currentDockedMolblocks);
                    if (result.profile && result.profile.length > 1) {
                        profCanvas.style.display = 'block';
                        drawProfile(profCanvas, result.profile);
                    } else {
                        profCanvas.style.display = 'none';
                    }
                } else {
                    optStatus.style.color = '#e0245e';
                    optStatus.textContent = 'Error: ' + result.error;
                }
            } catch (e) {
                optStatus.style.color = '#e0245e';
                optStatus.textContent = 'Request failed: ' + e.message;
            } finally {
                clearInterval(logInterval);
                buttons.forEach(b => b.disabled = false);
                activeBtn.textContent = label;
                refreshBondList();
            }
        }
        optBtn.onclick = () => runOptimization('auto');
        refineBtn.onclick = () => runOptimization('refine');

        if (!keyHandlerInstalled) {
            window.addEventListener('keydown', (event) => {
                if (!activeTransformControl || ['INPUT', 'TEXTAREA'].includes(document.activeElement?.tagName)) return;
                switch (event.key.toLowerCase()) {
                    case 't': activeTransformControl.setMode('translate'); break;
                    case 'r': activeTransformControl.setMode('rotate'); break;
                }
            });
            keyHandlerInstalled = true;
        }

        // --- ANIMATION LOOP ---
        function animate() {
            animationId = requestAnimationFrame(animate);
            orbit.update();
            updateReactionLines();
            renderer.render(scene, camera);
        }
        animate();

        // Read back the (possibly user-modified) coordinates from the scene
        TSEditorUI.extractData = () => moleculeGroups.map((group, i) => {
            group.updateMatrixWorld(true);
            return updateSDFCoordinates(parsedSDFs[i], group);
        });
    },

    getData: async (ui) => {
        if (TSEditorUI.extractData) {
            currentDockedMolblocks = TSEditorUI.extractData();
        }
        return { docked_molblocks: currentDockedMolblocks };
    }
};
