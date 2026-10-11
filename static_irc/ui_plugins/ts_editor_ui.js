/**
 * 3D Transition State Editor UI:
 * Three.js scene, interactive raycast atom picking, live distance HUD,
 * and orchestration of TS optimization strategies.
 */

import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';
import { parseMolBlock, updateMolBlockPositions } from './ts_sdf_parser.js';
import { drawProfile } from './ts_profile_chart.js';
import { createReactionCard, createRefinementCard, createOptimizationCard } from './ts_editor_cards.js';

let currentDockedMolblocks = [];
let moleculeGroups = [];
let picked = [];
let historyStack = [];
let rxnBonds = { forming: [], breaking: [] };
let overridden = false;
let dockVariant = 0;
let animationId = null;
let keyHandlerInstalled = false;
let activeTransformControl = null;
let initialCameraSet = false;

let scene, camera, renderer, controls, transformControl, dynamicBondGroup, distanceHUD;
let reactionCardObj, refineCardObj, optCardObj;
const raycaster = new THREE.Raycaster();
const mouse = new THREE.Vector2();

const ELEMENT_COLORS = {
    'C': 0x909090, 'H': 0xFFFFFF, 'O': 0xFF0D0D, 'N': 0x2250FF, 'S': 0xFFFF30,
    'P': 0xFF8000, 'F': 0x90E050, 'Cl': 0x1FF01F, 'Br': 0xA62929, 'I': 0x940094,
    'Mg': 0x8AFF00, 'Na': 0xAB5CF2, 'K': 0x8F40D4, 'Li': 0xCC80FF, 'B': 0xFFB5B5, 'Si': 0xF0C8A0
};
const ELEMENT_RADII = {
    'H': 0.3, 'C': 0.7, 'O': 0.6, 'N': 0.65, 'Br': 0.9, 'Cl': 0.8, 'I': 1.0, 'S': 0.8, 'Mg': 0.8
};

const sphereGeo = new THREE.SphereGeometry(1, 24, 24);
const cylGeo = new THREE.CylinderGeometry(0.12, 0.12, 1, 16);
cylGeo.translate(0, 0.5, 0);
cylGeo.rotateX(Math.PI / 2);
const bondMat = new THREE.MeshPhongMaterial({ color: 0x888888 });

function el(tag, style = {}, text = null) {
    const e = document.createElement(tag);
    Object.assign(e.style, style);
    if (text !== null) e.textContent = text;
    return e;
}

function styledButton(text, primary = false, small = false) {
    return el('button', {
        padding: small ? '5px 9px' : '9px 14px',
        margin: '4px 6px 4px 0',
        borderRadius: '6px',
        cursor: 'pointer',
        fontSize: small ? '11.5px' : '13px',
        border: primary ? 'none' : '1px solid #c7ccd6',
        fontWeight: primary ? '600' : '500',
        background: primary ? 'linear-gradient(135deg,#4f6ef7,#7b4ff7)' : '#fff',
        color: primary ? '#fff' : '#333',
        boxShadow: primary ? '0 2px 8px rgba(79,110,247,.35)' : 'none'
    }, text);
}

function card(title) {
    const c = el('div', {
        background: '#fff', border: '1px solid #e1e4ea', borderRadius: '10px',
        padding: '12px 14px', boxShadow: '0 2px 8px rgba(0,0,0,.03)',
        flexShrink: '0'
    });
    if (title) {
        c.appendChild(el('div', { fontWeight: '700', fontSize: '13px', color: '#1a1f36', marginBottom: '8px' }, title));
    }
    return c;
}

function computeCentroid(atoms) {
    if (!atoms || !atoms.length) return { cx: 0, cy: 0, cz: 0 };
    let sx = 0, sy = 0, sz = 0;
    atoms.forEach(a => { sx += a.x; sy += a.y; sz += a.z; });
    const n = atoms.length;
    return { cx: sx / n, cy: sy / n, cz: sz / n };
}

function buildAtomMeshes(atoms, fi, centroid) {
    const { cx, cy, cz } = centroid;
    return atoms.map((atom, ai) => {
        const mat = new THREE.MeshPhongMaterial({ color: ELEMENT_COLORS[atom.symbol] || 0xcccccc });
        const mesh = new THREE.Mesh(sphereGeo, mat);
        const r = ELEMENT_RADII[atom.symbol] || 0.5;
        mesh.scale.set(r, r, r);
        mesh.position.set(atom.x - cx, atom.y - cy, atom.z - cz);
        mesh.userData = { fi, ai, isAtom: true };
        if (picked.some(p => p.fi === fi && p.ai === ai)) {
            mat.emissive.setHex(0x3355ff);
        }
        return mesh;
    });
}

function buildBondMeshes(bonds, atomMeshes) {
    const bondMeshes = [];
    bonds.forEach(b => {
        const m1 = atomMeshes[b.from];
        const m2 = atomMeshes[b.to];
        if (!m1 || !m2) return;
        const mesh = new THREE.Mesh(cylGeo, bondMat);
        mesh.matrixAutoUpdate = false;
        mesh.userData = { m1, m2 };
        bondMeshes.push(mesh);
    });
    return bondMeshes;
}

function updateCylinder(mesh, p1, p2) {
    const dir = new THREE.Vector3().subVectors(p2, p1);
    const len = dir.length();
    mesh.position.copy(p1);
    mesh.scale.set(1, 1, Math.max(len, 0.001));
    mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), dir.clone().normalize());
    mesh.updateMatrix();
}

function createDashedLine(p1, p2, color) {
    const geo = new THREE.BufferGeometry().setFromPoints([p1, p2]);
    const mat = new THREE.LineDashedMaterial({ color, dashSize: 0.15, gapSize: 0.1, linewidth: 2 });
    const line = new THREE.Line(geo, mat);
    line.computeLineDistances();
    return line;
}

function getAtomWorldPos(groups, fi, ai) {
    const g = groups[fi];
    const atomMesh = g?.userData?.atomMeshes?.[ai];
    if (!atomMesh) return null;
    g.updateMatrixWorld(true);
    return atomMesh.position.clone().applyMatrix4(g.matrixWorld);
}

function centerCameraOnMolecules() {
    const box = new THREE.Box3();
    let hasAtoms = false;
    moleculeGroups.forEach(g => {
        g.updateMatrixWorld(true);
        (g.userData.atomMeshes || []).forEach(m => {
            const wp = m.position.clone().applyMatrix4(g.matrixWorld);
            box.expandByPoint(wp);
            hasAtoms = true;
        });
    });
    if (!hasAtoms) return;
    const center = new THREE.Vector3();
    box.getCenter(center);
    const size = new THREE.Vector3();
    box.getSize(size);
    const maxDim = Math.max(size.x, size.y, size.z, 4);
    controls.target.copy(center);
    camera.position.set(center.x, center.y, center.z + maxDim * 2.8);
    camera.lookAt(center);
    controls.update();
}

function pushHistory() {
    const snapshot = currentDockedMolblocks.map((mb, fi) => {
        const g = moleculeGroups[fi];
        return g ? updateMolBlockPositions(mb, g) : mb;
    });
    historyStack.push({
        molblocks: snapshot,
        rxnBonds: JSON.parse(JSON.stringify(rxnBonds)),
        overridden
    });
}

function popHistory() {
    if (!historyStack.length) return;
    const prev = historyStack.pop();
    currentDockedMolblocks = prev.molblocks;
    rxnBonds = prev.rxnBonds;
    overridden = prev.overridden;
    picked = [];
    buildMolecules(currentDockedMolblocks);
    reactionCardObj.refreshBondList();
    refineCardObj.undoBtn.disabled = historyStack.length === 0;
    optCardObj.optStatus.style.color = '#4f6ef7';
    optCardObj.optStatus.textContent = 'Reverted to previous state.';
}

function buildMolecules(molblocks) {
    transformControl.detach();
    moleculeGroups.forEach(g => scene.remove(g));
    moleculeGroups = [];

    molblocks.forEach((mb, fi) => {
        const parsed = parseMolBlock(mb);
        const centroid = computeCentroid(parsed.atoms);
        const group = new THREE.Group();
        group.position.set(centroid.cx, centroid.cy, centroid.cz);
        group.userData = { fi, parsed, centroid };

        const atomMeshes = buildAtomMeshes(parsed.atoms, fi, centroid);
        atomMeshes.forEach(m => group.add(m));
        group.userData.atomMeshes = atomMeshes;

        const bondMeshes = buildBondMeshes(parsed.bonds, atomMeshes);
        bondMeshes.forEach(bm => group.add(bm));
        group.userData.bondMeshes = bondMeshes;

        scene.add(group);
        moleculeGroups.push(group);
    });

    if (moleculeGroups.length > 1) transformControl.attach(moleculeGroups[1]);
    else if (moleculeGroups.length === 1) transformControl.attach(moleculeGroups[0]);

    updateBondsAll();
    updateDynamicBonds();
    updateLiveDistanceHUD();

    if (!initialCameraSet && moleculeGroups.length > 0) {
        centerCameraOnMolecules();
        initialCameraSet = true;
    }
}

function updateBondsAll() {
    moleculeGroups.forEach(g => {
        (g.userData.bondMeshes || []).forEach(bm => {
            updateCylinder(bm, bm.userData.m1.position, bm.userData.m2.position);
        });
    });
}

function updateDynamicBonds() {
    while (dynamicBondGroup.children.length) {
        dynamicBondGroup.remove(dynamicBondGroup.children[0]);
    }
    const drawDashed = (pair, color) => {
        const [[f1, a1], [f2, a2]] = pair;
        const p1 = getAtomWorldPos(moleculeGroups, f1, a1);
        const p2 = getAtomWorldPos(moleculeGroups, f2, a2);
        if (p1 && p2) dynamicBondGroup.add(createDashedLine(p1, p2, color));
    };
    rxnBonds.forming.forEach(p => drawDashed(p, 0x18a558));
    rxnBonds.breaking.forEach(p => drawDashed(p, 0xe0245e));
}

function updateLiveDistanceHUD() {
    if (picked.length < 2) {
        distanceHUD.style.display = 'none';
        return;
    }
    const p1 = getAtomWorldPos(moleculeGroups, picked[0].fi, picked[0].ai);
    const p2 = getAtomWorldPos(moleculeGroups, picked[1].fi, picked[1].ai);
    if (!p1 || !p2) {
        distanceHUD.style.display = 'none';
        return;
    }
    const dist = p1.distanceTo(p2);
    distanceHUD.textContent = `Distance: ${dist.toFixed(3)} Å (${picked[0].sym}···${picked[1].sym})`;
    distanceHUD.style.display = 'block';
}

function handleRaycastPick(event) {
    if (transformControl.dragging) return;
    const rect = renderer.domElement.getBoundingClientRect();
    mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
    mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
    raycaster.setFromCamera(mouse, camera);

    const allAtoms = moleculeGroups.flatMap(g => g.userData.atomMeshes || []);
    const intersects = raycaster.intersectObjects(allAtoms, false);
    if (!intersects.length) return;

    const mesh = intersects[0].object;
    const { fi, ai } = mesh.userData;
    const sym = moleculeGroups[fi].userData.parsed.atoms[ai].symbol;

    const existingIdx = picked.findIndex(p => p.fi === fi && p.ai === ai);
    if (existingIdx >= 0) {
        mesh.material.emissive.setHex(0x000000);
        picked.splice(existingIdx, 1);
    } else {
        if (picked.length >= 2) {
            const old = picked.shift();
            const oldMesh = moleculeGroups[old.fi]?.userData.atomMeshes[old.ai];
            if (oldMesh) oldMesh.material.emissive.setHex(0x000000);
        }
        picked.push({ fi, ai, sym });
        mesh.material.emissive.setHex(0x3355ff);
    }
    updateLiveDistanceHUD();
    reactionCardObj.refreshBondList();
}

async function refineFragment(fi, method) {
    pushHistory();
    optCardObj.optStatus.style.color = '#4f6ef7';
    optCardObj.optStatus.textContent = `Refining Reactant ${fi + 1} with ${method.toUpperCase()}...`;
    const g = moleculeGroups[fi];
    const currentBlock = updateMolBlockPositions(currentDockedMolblocks[fi], g);

    try {
        const payload = { action_type: 'refine_fragment', fragment_index: fi, method, molblock: currentBlock };
        const res = await fetch('/api/step/step_3_ts/custom', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        const r = await res.json();
        if (r.success) {
            currentDockedMolblocks[fi] = r.molblock;
            buildMolecules(currentDockedMolblocks);
            optCardObj.optStatus.style.color = '#18a558';
            optCardObj.optStatus.textContent = r.message;
        } else {
            optCardObj.optStatus.style.color = '#e0245e';
            optCardObj.optStatus.textContent = r.error;
        }
    } catch (err) {
        optCardObj.optStatus.style.color = '#e0245e';
        optCardObj.optStatus.textContent = 'Refine error: ' + err.message;
    }
}

async function resetToClean(cleanReactants) {
    pushHistory();
    if (cleanReactants && cleanReactants.length) {
        currentDockedMolblocks = [...cleanReactants];
        picked = [];
        buildMolecules(currentDockedMolblocks);
        optCardObj.optStatus.style.color = '#4f6ef7';
        optCardObj.optStatus.textContent = 'Reset scene to clean reactant conformers.';
        return;
    }
    try {
        const res = await fetch('/api/step/step_3_ts/custom', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'reset_clean' })
        });
        const r = await res.json();
        if (r.success) {
            currentDockedMolblocks = r.molblocks;
            picked = [];
            buildMolecules(currentDockedMolblocks);
            optCardObj.optStatus.style.color = '#4f6ef7';
            optCardObj.optStatus.textContent = r.message;
        }
    } catch (e) {
        optCardObj.optStatus.style.color = '#e0245e';
        optCardObj.optStatus.textContent = 'Reset failed: ' + e.message;
    }
}

async function nextDockingPose() {
    pushHistory();
    dockVariant += 1;
    optCardObj.nextPoseBtn.disabled = true;
    try {
        const payload = { action_type: 'auto_dock', variant: dockVariant, reaction_override: overridden ? rxnBonds : null };
        const res = await fetch('/api/step/step_3_ts/custom', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        const result = await res.json();
        if (result.success) {
            currentDockedMolblocks = result.molblocks;
            buildMolecules(result.molblocks);
            optCardObj.optStatus.style.color = '#4f6ef7';
            optCardObj.optStatus.textContent = result.message;
        } else {
            optCardObj.optStatus.style.color = '#e0245e';
            optCardObj.optStatus.textContent = result.error;
        }
    } finally {
        optCardObj.nextPoseBtn.disabled = false;
    }
}

async function fetchLog() {
    try {
        const r = await fetch('/api/log');
        const t = await r.text();
        if (t) {
            optCardObj.logBox.value = t;
            optCardObj.logBox.scrollTop = optCardObj.logBox.scrollHeight;
        }
    } catch (e) { }
}

function handleOptimizationResult(result) {
    if (result.success) {
        optCardObj.optStatus.style.color = '#18a558';
        optCardObj.optStatus.textContent = result.message;
        currentDockedMolblocks = result.optimized_molblocks;
        picked = [];
        buildMolecules(currentDockedMolblocks);
        if (result.profile && result.profile.length > 1) {
            optCardObj.profCanvas.style.display = 'block';
            drawProfile(optCardObj.profCanvas, result.profile);
        } else {
            optCardObj.profCanvas.style.display = 'none';
        }
        return;
    }
    if (result.diverged) {
        optCardObj.optStatus.style.color = '#d97706';
        optCardObj.optStatus.textContent = '⚠️ ' + result.error;
        if (result.profile && result.profile.length > 1) {
            optCardObj.profCanvas.style.display = 'block';
            drawProfile(optCardObj.profCanvas, result.profile);
        }
        return;
    }
    optCardObj.optStatus.style.color = '#e0245e';
    optCardObj.optStatus.textContent = '❌ Error: ' + (result.error || 'Optimization failed');
}

async function runOptimization(strategy) {
    pushHistory();
    const buttons = [optCardObj.runBtn, optCardObj.nextPoseBtn, refineCardObj.undoBtn, refineCardObj.resetCleanBtn];
    const label = optCardObj.runBtn.textContent;
    buttons.forEach(b => { if (b) b.disabled = true; });
    optCardObj.runBtn.textContent = 'Searching…';
    optCardObj.optStatus.style.color = '#4f6ef7';
    optCardObj.optStatus.textContent = `Running ${strategy.toUpperCase()} TS optimization...`;
    optCardObj.logBox.value = `Starting ${strategy.toUpperCase()} TS optimization…\n`;
    optCardObj.logBox.scrollIntoView({ behavior: 'smooth', block: 'nearest' });
    const logInterval = setInterval(fetchLog, 700);

    const data = TSEditorUI.extractData();
    try {
        const payload = {
            action_type: 'optimize_ts',
            strategy,
            docked_molblocks: data,
            reaction_override: overridden ? rxnBonds : null
        };
        const res = await fetch('/api/step/step_3_ts/custom', {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify(payload)
        });
        const result = await res.json();
        await fetchLog();
        handleOptimizationResult(result);
    } catch (e) {
        optCardObj.optStatus.style.color = '#e0245e';
        optCardObj.optStatus.textContent = 'Request failed: ' + e.message;
    } finally {
        clearInterval(logInterval);
        buttons.forEach(b => { if (b) b.disabled = false; });
        refineCardObj.undoBtn.disabled = historyStack.length === 0;
        optCardObj.runBtn.textContent = label;
        reactionCardObj.refreshBondList();
    }
}

function initSceneViewer(viewerContainer) {
    scene = new THREE.Scene();
    scene.background = new THREE.Color(0x0f111a);
    
    const w = viewerContainer.clientWidth || 650;
    const h = viewerContainer.clientHeight || 620;
    camera = new THREE.PerspectiveCamera(45, w / h, 0.1, 1000);
    camera.position.set(0, 0, 15);

    renderer = new THREE.WebGLRenderer({ antialias: true });
    renderer.setSize(w, h);
    viewerContainer.appendChild(renderer.domElement);

    scene.add(new THREE.AmbientLight(0xffffff, 0.8));
    const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
    dirLight.position.set(10, 20, 15);
    scene.add(dirLight);

    controls = new OrbitControls(camera, renderer.domElement);
    controls.enableDamping = true;

    transformControl = new TransformControls(camera, renderer.domElement);
    transformControl.size = 0.75;
    scene.add(transformControl);
    activeTransformControl = transformControl;

    transformControl.addEventListener('dragging-changed', event => {
        controls.enabled = !event.value;
        if (event.value) pushHistory();
    });
    transformControl.addEventListener('change', () => {
        updateDynamicBonds();
        updateLiveDistanceHUD();
    });

    dynamicBondGroup = new THREE.Group();
    scene.add(dynamicBondGroup);

    distanceHUD = el('div', {
        position: 'absolute', bottom: '10px', left: '10px', background: 'rgba(20,25,40,0.85)',
        color: '#fff', padding: '6px 12px', borderRadius: '6px', fontSize: '12px',
        fontFamily: 'monospace', pointerEvents: 'none', display: 'none'
    });
    viewerContainer.appendChild(distanceHUD);

    let pointerDownPos = null;
    renderer.domElement.addEventListener('pointerdown', (e) => {
        pointerDownPos = { x: e.clientX, y: e.clientY };
    });
    renderer.domElement.addEventListener('pointerup', (e) => {
        if (!pointerDownPos) return;
        const dist = Math.hypot(e.clientX - pointerDownPos.x, e.clientY - pointerDownPos.y);
        pointerDownPos = null;
        if (dist > 5) return;
        handleRaycastPick(e);
    });

    const resizeObserver = new ResizeObserver(() => {
        const cw = viewerContainer.clientWidth;
        const ch = viewerContainer.clientHeight;
        if (!cw || !ch) return;
        camera.aspect = cw / ch;
        camera.updateProjectionMatrix();
        renderer.setSize(cw, ch);
    });
    resizeObserver.observe(viewerContainer);
}

function renderTSEditor(schemaContainer, ui, state, stepId) {
    if (animationId) {
        cancelAnimationFrame(animationId);
        animationId = null;
    }

    const stepState = state?.[stepId] || {};
    currentDockedMolblocks = [
        ...(stepState.docked_reactants || stepState.optimized_molblocks || ui.reactants_3d || [])
    ];
    picked = [];
    historyStack = [];
    overridden = !!stepState.reaction_override;
    initialCameraSet = false;

    const rInfo = stepState.reaction_override || ui.reaction_info || {};
    rxnBonds = {
        forming: (rInfo.forming || []).map(p => [[p[0][0], p[0][1]], [p[1][0], p[1][1]]]),
        breaking: (rInfo.breaking || []).map(p => [[p[0][0], p[0][1]], [p[1][0], p[1][1]]])
    };

    schemaContainer.innerHTML = '';
    const layout = el('div', { display: 'flex', gap: '16px', height: '680px', width: '100%' });
    const viewerContainer = el('div', { flex: '1', position: 'relative', borderRadius: '10px', overflow: 'hidden', background: '#0f111a' });
    const controlsContainer = el('div', { width: '380px', display: 'flex', flexDirection: 'column', gap: '12px', overflowY: 'auto', paddingRight: '4px' });

    layout.appendChild(viewerContainer);
    layout.appendChild(controlsContainer);
    schemaContainer.appendChild(layout);

    initSceneViewer(viewerContainer);

    reactionCardObj = createReactionCard({
        rxnBonds,
        getPicked: () => picked,
        onBondsChanged: () => {
            overridden = true;
            updateDynamicBonds();
        },
        el, styledButton, card
    });
    controlsContainer.appendChild(reactionCardObj.card);

    refineCardObj = createRefinementCard({
        numReactants: currentDockedMolblocks.length || (ui.reactants_3d || []).length,
        onSelect: (idx) => { if (moleculeGroups[idx]) transformControl.attach(moleculeGroups[idx]); },
        onRefine: refineFragment,
        onGizmoMode: (mode) => {
            if (mode === 'detach') transformControl.detach();
            else transformControl.setMode(mode);
        },
        onUndo: popHistory,
        onReset: () => resetToClean(ui.clean_reactants_3d),
        el, styledButton, card
    });
    controlsContainer.appendChild(refineCardObj.card);

    const strats = [
        { id: 'drive', name: 'Constrained Scan (Drive R → P) → Sella (Recommended)' },
        { id: 'neb', name: 'Climbing Image NEB (CI-NEB with IDPP)' },
        { id: 'dimer', name: 'Dimer Method (ase.mep.dimer)' },
        { id: 'sella', name: 'Pre-Relaxation → Seeded Sella (v₀)' }
    ];

    optCardObj = createOptimizationCard({
        strats,
        showNextPose: (currentDockedMolblocks.length || (ui.reactants_3d || []).length) > 1,
        onRun: runOptimization,
        onNextPose: nextDockingPose,
        el, styledButton, card
    });
    controlsContainer.appendChild(optCardObj.card);
    controlsContainer.appendChild(optCardObj.profCanvas);
    if (optCardObj.logHeader) controlsContainer.appendChild(optCardObj.logHeader);
    controlsContainer.appendChild(optCardObj.logBox);

    if (!keyHandlerInstalled) {
        window.addEventListener('keydown', event => {
            if (!activeTransformControl || ['INPUT', 'TEXTAREA'].includes(document.activeElement?.tagName)) return;
            switch (event.key.toLowerCase()) {
                case 't': activeTransformControl.setMode('translate'); break;
                case 'r': activeTransformControl.setMode('rotate'); break;
                case 'escape': activeTransformControl.detach(); break;
            }
        });
        keyHandlerInstalled = true;
    }

    function animate() {
        animationId = requestAnimationFrame(animate);
        controls.update();
        renderer.render(scene, camera);
    }
    animate();

    buildMolecules(currentDockedMolblocks);
}

export const TSEditorUI = {
    render: renderTSEditor,
    getData: async (ui) => ({
        docked_molblocks: TSEditorUI.extractData(),
        reaction_override: overridden ? rxnBonds : null
    }),
    extractData: () => {
        return currentDockedMolblocks.map((mb, fi) => {
            const g = moleculeGroups[fi];
            return g ? updateMolBlockPositions(mb, g) : mb;
        });
    }
};
