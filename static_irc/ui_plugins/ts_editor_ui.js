import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { TransformControls } from 'three/addons/controls/TransformControls.js';

let currentDockedMolblocks = [];

// Simple CPK color mapping
const elementColors = {
    'C': 0x909090,
    'H': 0xFFFFFF,
    'O': 0xFF0D0D,
    'N': 0x2250FF,
    'S': 0xFFFF30,
    'P': 0xFF8000,
    'F': 0x90E050,
    'Cl': 0x1FF01F,
    'Br': 0xA62929,
    'I': 0x940094,
    'Mg': 0x8AFF00
};
const elementRadii = {
    'H': 0.3,
    'C': 0.7,
    'O': 0.6,
    'N': 0.65,
    'Br': 0.9,
    'Mg': 0.8
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
    // group.children contains Meshes, but we can just use the global transform of the group
    // to transform the original atomic coordinates.
    let matrix = group.matrixWorld;
    let newLines = [...parsedSDF.originalLines];
    
    for (let i = 0; i < parsedSDF.numAtoms; i++) {
        let atom = parsedSDF.atoms[i];
        let vec = new THREE.Vector3(atom.x, atom.y, atom.z);
        vec.applyMatrix4(matrix); // applies translation, rotation, scale
        
        let lineIdx = 4 + i;
        let line = newLines[lineIdx];
        
        let newXStr = vec.x.toFixed(4).padStart(10, ' ');
        let newYStr = vec.y.toFixed(4).padStart(10, ' ');
        let newZStr = vec.z.toFixed(4).padStart(10, ' ');
        
        newLines[lineIdx] = newXStr + newYStr + newZStr + line.substring(30);
    }
    return newLines.join('\n');
}

export const TSEditorUI = {
    render: (container, ui, state, stepId) => {
        container.innerHTML = '';
        currentDockedMolblocks = [...ui.reactants_3d]; 

        if (ui.message) {
            const msg = document.createElement('p');
            msg.textContent = ui.message;
            container.appendChild(msg);
        }

        const editorLayout = document.createElement('div');
        editorLayout.style.display = 'flex';
        editorLayout.style.gap = '20px';
        editorLayout.style.marginTop = '20px';

        const viewerContainer = document.createElement('div');
        viewerContainer.style.width = '700px';
        viewerContainer.style.height = '600px';
        viewerContainer.style.position = 'relative';
        viewerContainer.style.border = '1px solid #ccc';
        viewerContainer.style.borderRadius = '8px';
        editorLayout.appendChild(viewerContainer);

        const controlsContainer = document.createElement('div');
        controlsContainer.style.flex = '1';
        controlsContainer.style.display = 'flex';
        controlsContainer.style.flexDirection = 'column';
        controlsContainer.style.gap = '15px';
        editorLayout.appendChild(controlsContainer);

        container.appendChild(editorLayout);

        // --- THREE.JS SETUP ---
        const scene = new THREE.Scene();
        scene.background = new THREE.Color(0xf5f5f5);

        const camera = new THREE.PerspectiveCamera(45, 700 / 600, 0.1, 100);
        camera.position.set(0, 0, 15);

        const renderer = new THREE.WebGLRenderer({ antialias: true });
        renderer.setSize(700, 600);
        viewerContainer.appendChild(renderer.domElement);

        const ambientLight = new THREE.AmbientLight(0x404040, 2); // soft white light
        scene.add(ambientLight);
        const dirLight = new THREE.DirectionalLight(0xffffff, 1.5);
        dirLight.position.set(10, 10, 10);
        scene.add(dirLight);

        const orbit = new OrbitControls(camera, renderer.domElement);
        orbit.enableDamping = true;

        const transformControl = new TransformControls(camera, renderer.domElement);
        transformControl.addEventListener('dragging-changed', function (event) {
            orbit.enabled = !event.value;
        });
        scene.add(transformControl);

        let moleculeGroups = [];
        let parsedSDFs = [];

        // Build 3D Models
        ui.reactants_3d.forEach((molBlock, idx) => {
            if (!molBlock) return;
            const parsed = parseSDF(molBlock);
            if (!parsed) return;
            
            parsedSDFs.push(parsed);
            
            const group = new THREE.Group();
            
            const sphereGeo = new THREE.SphereGeometry(1, 32, 32);
            
            // Atoms
            parsed.atoms.forEach(atom => {
                const color = elementColors[atom.symbol] || 0xcccccc;
                const radius = elementRadii[atom.symbol] || 0.5;
                const mat = new THREE.MeshPhongMaterial({ color: color, shininess: 50 });
                const mesh = new THREE.Mesh(sphereGeo, mat);
                mesh.scale.set(radius, radius, radius);
                mesh.position.set(atom.x, atom.y, atom.z);
                group.add(mesh);
            });
            
            // Bonds
            const cylGeo = new THREE.CylinderGeometry(0.15, 0.15, 1, 16);
            parsed.bonds.forEach(bond => {
                const a1 = parsed.atoms[bond.a1];
                const a2 = parsed.atoms[bond.a2];
                const p1 = new THREE.Vector3(a1.x, a1.y, a1.z);
                const p2 = new THREE.Vector3(a2.x, a2.y, a2.z);
                
                const distance = p1.distanceTo(p2);
                const mat = new THREE.MeshPhongMaterial({ color: 0x999999 });
                const mesh = new THREE.Mesh(cylGeo, mat);
                
                // Position and orient cylinder
                mesh.position.copy(p1).lerp(p2, 0.5);
                mesh.quaternion.setFromUnitVectors(new THREE.Vector3(0, 1, 0), p2.clone().sub(p1).normalize());
                mesh.scale.set(1, distance, 1);
                
                group.add(mesh);
            });
            
            scene.add(group);
            moleculeGroups.push(group);
        });

        // --- UI CONTROLS ---
        const helpText = document.createElement('div');
        helpText.innerHTML = `
            <h3>Blender-style 3D Editor</h3>
            <p><strong>Left Click & Drag:</strong> Rotate Camera</p>
            <p><strong>Right Click & Drag:</strong> Pan Camera</p>
            <p><strong>Scroll:</strong> Zoom</p>
            <p>Select a molecule below to activate its transform gizmo. Use the arrows to Translate and the rings to Rotate.</p>
        `;
        controlsContainer.appendChild(helpText);

        const activeGroupSelect = document.createElement('div');
        moleculeGroups.forEach((group, idx) => {
            const btn = document.createElement('button');
            btn.textContent = `Select Reactant ${idx + 1}`;
            btn.style.padding = '10px';
            btn.style.margin = '5px';
            btn.onclick = () => {
                transformControl.attach(group);
            };
            activeGroupSelect.appendChild(btn);
        });
        controlsContainer.appendChild(activeGroupSelect);

        const modeSelect = document.createElement('div');
        modeSelect.innerHTML = `<p><strong>Gizmo Mode:</strong></p>`;
        
        const transBtn = document.createElement('button');
        transBtn.textContent = 'Translate (T)';
        transBtn.onclick = () => transformControl.setMode('translate');
        
        const rotBtn = document.createElement('button');
        rotBtn.textContent = 'Rotate (R)';
        rotBtn.onclick = () => transformControl.setMode('rotate');
        
        transBtn.style.margin = '5px';
        rotBtn.style.margin = '5px';
        modeSelect.appendChild(transBtn);
        modeSelect.appendChild(rotBtn);
        controlsContainer.appendChild(modeSelect);
        
        // Add Optimize button
        const optBtn = document.createElement('button');
        optBtn.className = 'primary';
        optBtn.style.marginTop = '20px';
        optBtn.textContent = 'Optimize TS Guess';
        
        const optStatus = document.createElement('p');
        optStatus.style.fontWeight = 'bold';
        optStatus.style.color = '#007bff';
        
        optBtn.onclick = async () => {
            optBtn.disabled = true;
            optBtn.textContent = 'Optimizing... (Simulated)';
            optStatus.textContent = '';
            
            let data = TSEditorUI.extractData();
            try {
                let res = await fetch('/api/step/step_3_ts/custom', {
                    method: 'POST',
                    headers: { 'Content-Type': 'application/json' },
                    body: JSON.stringify({
                        action_type: 'optimize_ts',
                        docked_molblocks: data
                    })
                });
                let result = await res.json();
                if (result.success) {
                    optStatus.textContent = result.message;
                    currentDockedMolblocks = result.optimized_molblocks;
                    // In a real app we'd re-parse the optimized molblocks and render
                    // For now, they are identical, so we just show the success message
                } else {
                    optStatus.style.color = 'red';
                    optStatus.textContent = 'Error: ' + result.error;
                }
            } catch (e) {
                optStatus.style.color = 'red';
                optStatus.textContent = 'Request failed: ' + e.message;
            } finally {
                optBtn.disabled = false;
                optBtn.textContent = 'Optimize TS Guess';
            }
        };
        
        controlsContainer.appendChild(optBtn);
        controlsContainer.appendChild(optStatus);
        
        // Key bindings for transform modes
        window.addEventListener('keydown', function (event) {
            switch (event.key.toLowerCase()) {
                case 't': transformControl.setMode('translate'); break;
                case 'r': transformControl.setMode('rotate'); break;
            }
        });

        if (moleculeGroups.length > 0) {
            transformControl.attach(moleculeGroups[0]);
        }

        // --- ANIMATION LOOP ---
        function animate() {
            requestAnimationFrame(animate);
            orbit.update();
            renderer.render(scene, camera);
        }
        animate();

        // Save method attached to window/UI for data retrieval
        // We will compute the modified SDFs directly from the scene groups!
        TSEditorUI.extractData = () => {
            return moleculeGroups.map((group, i) => {
                group.updateMatrixWorld(true);
                return updateSDFCoordinates(parsedSDFs[i], group);
            });
        };
    },

    getData: async (ui) => {
        if (TSEditorUI.extractData) {
            currentDockedMolblocks = TSEditorUI.extractData();
        }
        return { docked_molblocks: currentDockedMolblocks };
    }
};
