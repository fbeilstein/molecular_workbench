import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';
import { parseMolBlock } from './ts_sdf_parser.js';

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

export class TSEditorViewer {
    constructor(container, onAtomPicked, onTransformChange) {
        this.container = container;
        this.onAtomPicked = onAtomPicked;
        this.onTransformChange = onTransformChange;
        this.moleculeGroups = [];
        this.initialCameraSet = false;
        this.picked = [];
        this._initScene();
    }

    _initScene() {
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x0f111a);
        const w = this.container.clientWidth || 650;
        const h = this.container.clientHeight || 620;
        this.camera = new THREE.PerspectiveCamera(45, w / h, 0.1, 1000);
        this.camera.position.set(0, 0, 15);

        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(w, h);
        this.container.appendChild(this.renderer.domElement);

        this.scene.add(new THREE.AmbientLight(0xffffff, 0.8));
        const dirLight = new THREE.DirectionalLight(0xffffff, 0.6);
        dirLight.position.set(10, 20, 15);
        this.scene.add(dirLight);

        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;

        this.transformControl = new TransformControls(this.camera, this.renderer.domElement);
        this.transformControl.size = 0.75;
        this.transformControl.addEventListener('dragging-changed', e => {
            this.controls.enabled = !e.value;
        });
        this.transformControl.addEventListener('change', () => {
            this._updateBondsAll();
            if (this.onTransformChange) this.onTransformChange();
        });
        this.scene.add(this.transformControl);

        this.dynamicBondGroup = new THREE.Group();
        this.scene.add(this.dynamicBondGroup);

        this.raycaster = new THREE.Raycaster();
        this.mouse = new THREE.Vector2();

        this.renderer.domElement.addEventListener('pointerdown', e => this._onPointerDown(e));
        this._animate();
    }

    _animate() {
        requestAnimationFrame(() => this._animate());
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    _onPointerDown(event) {
        if (this.transformControl.dragging) return;
        const rect = this.renderer.domElement.getBoundingClientRect();
        this.mouse.x = ((event.clientX - rect.left) / rect.width) * 2 - 1;
        this.mouse.y = -((event.clientY - rect.top) / rect.height) * 2 + 1;
        this.raycaster.setFromCamera(this.mouse, this.camera);

        const allAtoms = this.moleculeGroups.flatMap(g => g.userData.atomMeshes || []);
        const hits = this.raycaster.intersectObjects(allAtoms, false);
        if (!hits.length) return;
        const mesh = hits[0].object;
        const { fi, ai } = mesh.userData;
        const sym = this.moleculeGroups[fi].userData.parsed.atoms[ai].symbol;
        if (this.onAtomPicked) this.onAtomPicked(fi, ai, sym, mesh);
    }

    buildMolecules(molblocks) {
        this.transformControl.detach();
        this.moleculeGroups.forEach(g => this.scene.remove(g));
        this.moleculeGroups = [];

        molblocks.forEach((mb, fi) => {
            const parsed = parseMolBlock(mb);
            const centroid = this._computeCentroid(parsed.atoms);
            const group = new THREE.Group();
            group.position.set(centroid.cx, centroid.cy, centroid.cz);
            group.userData = { fi, parsed, centroid };

            const atomMeshes = parsed.atoms.map((atom, ai) => {
                const mat = new THREE.MeshPhongMaterial({ color: ELEMENT_COLORS[atom.symbol] || 0xcccccc });
                const mesh = new THREE.Mesh(sphereGeo, mat);
                const r = ELEMENT_RADII[atom.symbol] || 0.5;
                mesh.scale.set(r, r, r);
                mesh.position.set(atom.x - centroid.cx, atom.y - centroid.cy, atom.z - centroid.cz);
                mesh.userData = { fi, ai, isAtom: true };
                return mesh;
            });
            atomMeshes.forEach(m => group.add(m));
            group.userData.atomMeshes = atomMeshes;

            const bondMeshes = [];
            parsed.bonds.forEach(b => {
                const m1 = atomMeshes[b.from];
                const m2 = atomMeshes[b.to];
                if (!m1 || !m2) return;
                const bm = new THREE.Mesh(cylGeo, bondMat);
                bm.matrixAutoUpdate = false;
                bm.userData = { m1, m2 };
                group.add(bm);
                bondMeshes.push(bm);
            });
            group.userData.bondMeshes = bondMeshes;

            this.scene.add(group);
            this.moleculeGroups.push(group);
        });

        if (this.moleculeGroups.length > 1) this.transformControl.attach(this.moleculeGroups[1]);
        else if (this.moleculeGroups.length === 1) this.transformControl.attach(this.moleculeGroups[0]);

        this._updateBondsAll();
        if (!this.initialCameraSet && this.moleculeGroups.length > 0) {
            this.centerCamera();
            this.initialCameraSet = true;
        }
    }

    _computeCentroid(atoms) {
        if (!atoms || !atoms.length) return { cx: 0, cy: 0, cz: 0 };
        let sx = 0, sy = 0, sz = 0;
        atoms.forEach(a => { sx += a.x; sy += a.y; sz += a.z; });
        return { cx: sx / atoms.length, cy: sy / atoms.length, cz: sz / atoms.length };
    }

    _updateBondsAll() {
        this.moleculeGroups.forEach(g => {
            (g.userData.bondMeshes || []).forEach(bm => {
                const p1 = bm.userData.m1.position;
                const p2 = bm.userData.m2.position;
                const dir = new THREE.Vector3().subVectors(p2, p1);
                bm.position.copy(p1);
                bm.scale.set(1, 1, Math.max(dir.length(), 0.001));
                bm.quaternion.setFromUnitVectors(new THREE.Vector3(0, 0, 1), dir.clone().normalize());
                bm.updateMatrix();
            });
        });
    }

    getAtomWorldPos(fi, ai) {
        const g = this.moleculeGroups[fi];
        const mesh = g?.userData?.atomMeshes?.[ai];
        if (!mesh) return null;
        g.updateMatrixWorld(true);
        return mesh.position.clone().applyMatrix4(g.matrixWorld);
    }

    updateDynamicBonds(rxnBonds) {
        while (this.dynamicBondGroup.children.length) {
            this.dynamicBondGroup.remove(this.dynamicBondGroup.children[0]);
        }
        const drawDashed = (pair, color) => {
            const [[f1, a1], [f2, a2]] = pair;
            const p1 = this.getAtomWorldPos(f1, a1);
            const p2 = this.getAtomWorldPos(f2, a2);
            if (!p1 || !p2) return;
            const geo = new THREE.BufferGeometry().setFromPoints([p1, p2]);
            const mat = new THREE.LineDashedMaterial({ color, dashSize: 0.15, gapSize: 0.1, linewidth: 2 });
            const line = new THREE.Line(geo, mat);
            line.computeLineDistances();
            this.dynamicBondGroup.add(line);
        };
        (rxnBonds.forming || []).forEach(p => drawDashed(p, 0x18a558));
        (rxnBonds.breaking || []).forEach(p => drawDashed(p, 0xe0245e));
    }

    centerCamera() {
        const box = new THREE.Box3();
        let hasAtoms = false;
        this.moleculeGroups.forEach(g => {
            g.updateMatrixWorld(true);
            (g.userData.atomMeshes || []).forEach(m => {
                box.expandByPoint(m.position.clone().applyMatrix4(g.matrixWorld));
                hasAtoms = true;
            });
        });
        if (!hasAtoms) return;
        const center = new THREE.Vector3(), size = new THREE.Vector3();
        box.getCenter(center);
        box.getSize(size);
        const maxDim = Math.max(size.x, size.y, size.z, 4);
        this.controls.target.copy(center);
        this.camera.position.set(center.x, center.y, center.z + maxDim * 2.8);
        this.camera.lookAt(center);
        this.controls.update();
    }

    setGizmoMode(mode) {
        if (mode === 'detach') this.transformControl.detach();
        else this.transformControl.setMode(mode);
    }

    selectGroup(fi) {
        const g = this.moleculeGroups[fi];
        if (g) this.transformControl.attach(g);
    }
}
