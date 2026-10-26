/**
 * Unified 3D Molecular & Orbital Viewer for IRC Pipeline.
 * Consolidates all 3Dmol rendering into a single high-performance component:
 * 1. True bond order rendering (MolBlock V2000: single, double, triple bonds).
 * 2. Trajectory playback & frame scrubbing with camera orientation preservation.
 * 3. Volumetric orbital isosurface rendering with dual-phase coloring and VolumeData caching.
 * 4. Electrostatic potential (ESP) charge surface mapping.
 * 5. Boundary overlap inspection (Cyan/Magenta dual-model superposition).
 */

export class MoleculeViewer3D {
    constructor(containerEl, options = {}) {
        this.containerEl = containerEl;
        this.options = options;
        this.viewer = null;
        this.statusEl = options.statusEl || null;
        this.cubeCache = options.cubeCache || {};
        this.vdoCache = {};
        this.hasInitializedCamera = false;
        this.moleculeScale = options.scale || 1.0;
        this.currentFrameIdx = 0;
        this.bondTable = options.bonds || null; // [[i, j, order], ...]
    }

    init() {
        if (!this.containerEl) return;
        this.containerEl.innerHTML = '';
        if (typeof $3Dmol === 'undefined') {
            this.containerEl.innerHTML = '<div style="padding: 1rem; color: #dc2626; font-size: 0.85rem;">3Dmol.js library not loaded.</div>';
            return;
        }

        this.viewer = $3Dmol.createViewer(this.containerEl, {
            backgroundColor: this.options.backgroundColor || '#0f172a',
            id: 'unified-3dmol-canvas'
        });

        // Mouse wheel scales 3D molecule without interfering with page scroll
        this.containerEl.addEventListener('wheel', (e) => {
            e.preventDefault();
            const factor = e.deltaY < 0 ? 1.12 : 0.89;
            this.zoom(factor);
        }, { passive: false });

        // Ensure WebGL viewport resizes properly after DOM layout reflow
        requestAnimationFrame(() => {
            if (this.viewer) {
                this.viewer.resize();
                this.viewer.render();
            }
        });
        setTimeout(() => {
            if (this.viewer) {
                this.viewer.resize();
                this.viewer.render();
            }
        }, 150);

        this._resizeHandler = () => this.resize();
        window.addEventListener('resize', this._resizeHandler);
    }

    zoom(factor) {
        if (!this.viewer) return;
        this.viewer.zoom(factor);
        this.viewer.render();
    }

    zoomIn() {
        this.zoom(1.2);
    }

    zoomOut() {
        this.zoom(0.83);
    }

    resetCamera() {
        if (!this.viewer) return;
        this.viewer.zoomTo();
        this.viewer.render();
    }

    setMoleculeScale(scale) {
        this.moleculeScale = Math.max(0.4, Math.min(3.0, scale));
        if (this.viewer) {
            this.viewer.setStyle({}, {
                stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
                sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
            });
            this.viewer.render();
        }
    }

    setStatus(text, isError = false) {
        if (!this.statusEl) return;
        this.statusEl.textContent = text;
        this.statusEl.style.color = isError ? '#ef4444' : '#94a3b8';
    }

    setBondTable(bonds) {
        this.bondTable = bonds;
    }

    /**
     * Builds a standard V2000 MolBlock string to preserve authentic single,
     * double, and triple bond rendering in 3Dmol.
     */
    buildMolBlock(symbols, coords) {
        const nAtoms = Math.min(symbols.length, coords.length);
        const bonds = this._inferBonds(symbols, coords);
        const nBonds = bonds.length;

        // V2000 Header
        let mol = `Molecule\n  IRC-Pipeline\n\n`;
        // Counts line: 3 chars for nAtoms, 3 chars for nBonds, V2000 tag
        mol += `${String(nAtoms).padStart(3, ' ')}${String(nBonds).padStart(3, ' ')}  0  0  0  0  0  0  0  0999 V2000\n`;

        // Atom block: x, y, z (10.4f), symbol (3s)
        for (let i = 0; i < nAtoms; i++) {
            const sym = symbols[i];
            const [x, y, z] = coords[i];
            const xs = x.toFixed(4).padStart(10, ' ');
            const ys = y.toFixed(4).padStart(10, ' ');
            const zs = z.toFixed(4).padStart(10, ' ');
            const s = sym.padEnd(3, ' ');
            mol += `${xs}${ys}${zs} ${s} 0  0  0  0  0  0  0  0  0  0  0  0\n`;
        }

        // Bond block: atom1 (3d), atom2 (3d), bond_type (3d), bond_stereo (3d)
        for (const [a1, a2, order] of bonds) {
            const b1 = String(a1 + 1).padStart(3, ' ');
            const b2 = String(a2 + 1).padStart(3, ' ');
            const bo = String(order).padStart(3, ' ');
            mol += `${b1}${b2}${bo}  0  0  0  0\n`;
        }

        mol += `M  END\n`;
        return mol;
    }

    _inferBonds(symbols, coords) {
        // Physical distance-based bond computation with 6-ring Kekulization and trivalent oxygen
        const bonds = [];
        const n = Math.min(symbols.length, coords.length);
        const valences = new Array(n).fill(0);
        const candidates = [];
        const adj = Array.from({ length: n }, () => []);
        const distMat = {};

        const covRadii = {
            H: 0.31, C: 0.76, N: 0.71, O: 0.66, F: 0.57,
            P: 1.07, S: 1.05, Cl: 1.02, Br: 1.20, I: 1.39
        };
        const maxValence = { H: 1, F: 1, Cl: 1, Br: 1, I: 1, O: 3, N: 4, C: 4, S: 6, P: 5 };

        for (let i = 0; i < n; i++) {
            const s1 = symbols[i];
            const r1 = covRadii[s1] || 0.75;
            for (let j = i + 1; j < n; j++) {
                const s2 = symbols[j];
                const r2 = covRadii[s2] || 0.75;
                if (s1 === 'H' && s2 === 'H' && n > 2) continue;

                const dx = coords[i][0] - coords[j][0];
                const dy = coords[i][1] - coords[j][1];
                const dz = coords[i][2] - coords[j][2];
                const d = Math.sqrt(dx * dx + dy * dy + dz * dz);
                const maxD = 1.28 * (r1 + r2);

                if (d >= 0.4 && d <= maxD) {
                    adj[i].push(j);
                    adj[j].push(i);
                    const key = `${i}_${j}`;
                    distMat[key] = d;
                    candidates.push({ i, j, s1, s2, d });
                }
            }
        }

        const getDist = (u, v) => {
            const key = u < v ? `${u}_${v}` : `${v}_${u}`;
            return distMat[key] || 1.4;
        };

        // 1. Detect 6-membered aromatic/conjugated rings (primarily C, N)
        const ringDoubleBonds = new Set();
        const ringBondsAll = new Set();
        const ringAtoms = [];
        for (let i = 0; i < n; i++) {
            if (symbols[i] === 'C' || symbols[i] === 'N') {
                ringAtoms.push(i);
            }
        }

        const cycles = [];
        const visited = new Set();
        const path = [];

        const findCycles = (v, start, depth) => {
            if (depth === 6) {
                if (adj[v].includes(start)) {
                    cycles.push([...path]);
                }
                return;
            }
            for (const u of adj[v]) {
                if ((symbols[u] === 'C' || symbols[u] === 'N') && !visited.has(u)) {
                    visited.add(u);
                    path.push(u);
                    findCycles(u, start, depth + 1);
                    path.pop();
                    visited.delete(u);
                }
            }
        };

        for (const startAtom of ringAtoms) {
            visited.clear();
            visited.add(startAtom);
            path.length = 0;
            path.push(startAtom);
            findCycles(startAtom, startAtom, 1);
        }

        const uniqueRings = [];
        const seenRings = new Set();
        for (const c of cycles) {
            const sorted = [...c].sort((a, b) => a - b).join(',');
            if (!seenRings.has(sorted)) {
                seenRings.add(sorted);
                uniqueRings.push(c);
            }
        }

        // Kekulize each 6-membered ring:
        // A 6-cycle has two alternating matchings: edges0 and edges1.
        // Pick matching with smaller total bond length sum.
        for (const ring of uniqueRings) {
            let isConjugated = true;
            for (let k = 0; k < 6; k++) {
                const u = ring[k];
                const v = ring[(k + 1) % 6];
                if (getDist(u, v) > 1.48) {
                    isConjugated = false;
                    break;
                }
            }
            if (!isConjugated) continue;

            const edges0 = [[ring[0], ring[1]], [ring[2], ring[3]], [ring[4], ring[5]]];
            const edges1 = [[ring[1], ring[2]], [ring[3], ring[4]], [ring[5], ring[0]]];
            const sum0 = edges0.reduce((acc, [u, v]) => acc + getDist(u, v), 0);
            const sum1 = edges1.reduce((acc, [u, v]) => acc + getDist(u, v), 0);

            const chosen = sum0 <= sum1 ? edges0 : edges1;
            for (const [u, v] of chosen) {
                const key = u < v ? `${u}_${v}` : `${v}_${u}`;
                ringDoubleBonds.add(key);
            }
            for (let k = 0; k < 6; k++) {
                const u = ring[k];
                const v = ring[(k + 1) % 6];
                const key = u < v ? `${u}_${v}` : `${v}_${u}`;
                ringBondsAll.add(key);
            }
        }

        for (const c of candidates) {
            const pairKey = c.i < c.j ? `${c.i}_${c.j}` : `${c.j}_${c.i}`;
            const pair = [c.s1, c.s2].sort().join('-');
            let targetOrder = 1;

            if (ringDoubleBonds.has(pairKey)) {
                targetOrder = 2;
            } else if (!ringBondsAll.has(pairKey)) {
                if (pair === 'C-O') {
                    if (c.d < 1.16) targetOrder = 3;
                    else if (c.d < 1.265) targetOrder = 2;
                } else if (pair === 'C-C') {
                    if (c.d < 1.24) targetOrder = 3;
                    else if (c.d < 1.39) targetOrder = 2;
                } else if (pair === 'C-N') {
                    if (c.d < 1.18) targetOrder = 3;
                    else if (c.d < 1.32) targetOrder = 2;
                } else if (pair === 'N-N') {
                    if (c.d < 1.14) targetOrder = 3;
                    else if (c.d < 1.28) targetOrder = 2;
                } else if (pair === 'N-O' || pair === 'O-O') {
                    if (c.d < 1.25) targetOrder = 2;
                }
            }

            bonds.push({
                a1: c.i, a2: c.j, order: 1, target: targetOrder,
                s1: c.s1, s2: c.s2, d: c.d, pairKey
            });
            valences[c.i]++;
            valences[c.j]++;
        }

        // Priority sort: triple bonds first (target 3), then ring double bonds, then other double bonds, then by length
        bonds.sort((a, b) => {
            const prioA = a.target === 3 ? 0 : (ringDoubleBonds.has(a.pairKey) ? 1 : (a.target === 2 ? 2 : 3));
            const prioB = b.target === 3 ? 0 : (ringDoubleBonds.has(b.pairKey) ? 1 : (b.target === 2 ? 2 : 3));
            if (prioA !== prioB) return prioA - prioB;
            return a.d - b.d;
        });

        // Promote to double/triple bond if atom valences allow
        for (const b of bonds) {
            if (b.target > 1) {
                const v1Cap = maxValence[b.s1] || 4;
                const v2Cap = maxValence[b.s2] || 4;
                const desiredBoost = b.target - 1;
                const maxBoost = Math.min(v1Cap - valences[b.a1], v2Cap - valences[b.a2]);
                const boost = Math.max(0, Math.min(desiredBoost, maxBoost));
                if (boost > 0) {
                    b.order += boost;
                    valences[b.a1] += boost;
                    valences[b.a2] += boost;
                }
            }
        }

        return bonds.map(b => [b.a1, b.a2, b.order]);
    }

    /**
     * Renders a molecular structure with authentic double bonds.
     */
    renderMolecule(symbols, coords, charges = null, showCharges = false) {
        if (!this.viewer || !symbols || !coords) return;
        this.viewer.clear();

        const molblock = this.buildMolBlock(symbols, coords);
        this.viewer.addModel(molblock, 'mol');

        this.viewer.setStyle({}, {
            stick: { radius: 0.14 * this.moleculeScale, colorscheme: 'Jmol' },
            sphere: { scale: 0.28 * this.moleculeScale, colorscheme: 'Jmol' }
        });

        if (showCharges && charges && charges.length >= symbols.length) {
            const atoms = this.viewer.selectedAtoms({});
            for (let i = 0; i < atoms.length && i < charges.length; i++) {
                atoms[i].partialCharge = charges[i] || 0.0;
                atoms[i].properties = atoms[i].properties || {};
                atoms[i].properties.partialCharge = charges[i] || 0.0;
            }
            this.viewer.addSurface($3Dmol.SurfaceType.VDW, {
                opacity: 0.65,
                colorscheme: { prop: 'partialCharge', gradient: 'rwb', min: -0.6, max: 0.6 }
            });
            for (let i = 0; i < symbols.length; i++) {
                const q = charges[i];
                const [x, y, z] = coords[i];
                this.viewer.addLabel(`${symbols[i]}: ${q > 0 ? '+' : ''}${q.toFixed(2)}`, {
                    position: { x, y: y + 0.35, z },
                    fontSize: 10, fontColor: '#ffffff', backgroundColor: 'rgba(15, 23, 42, 0.75)', backgroundOpacity: 0.75
                });
            }
        }

        this.viewer.render();
    }

    /**
     * Renders active volumetric orbital isosurfaces with phase coloring.
     */
    renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx) {
        if (!this.viewer || !activeOrbitals) return;

        for (const orb of activeOrbitals) {
            const orbKey = orb.orbital_key || orb.key || orb;
            const cubeData = this.cubeCache[`${frameIdx}_${orbKey}`];
            if (cubeData) {
                const vKey = `${frameIdx}_${orbKey}`;
                let vdo = this.vdoCache[vKey];
                if (!vdo) {
                    try {
                        vdo = new $3Dmol.VolumeData(cubeData, 'cube');
                        this.vdoCache[vKey] = vdo;
                    } catch (e) {
                        console.error('VolumeData error:', e);
                    }
                }
                if (vdo) {
                    try {
                        const isAnti = (orb.orbital_class === 'antibonding') || (orb.type && orb.type.includes('star'));
                        // Bonding: Cyan (+), Blue (-). Anti-bonding: Orange (+), Magenta (-)
                        const posColor = isAnti ? '#fb923c' : '#38bdf8';
                        const negColor = isAnti ? '#db2777' : '#1d4ed8';

                        this.viewer.addIsosurface(vdo, { isoval: Math.abs(isovalue), color: posColor, opacity: 0.72 });
                        this.viewer.addIsosurface(vdo, { isoval: -Math.abs(isovalue), color: negColor, opacity: 0.72 });
                    } catch (e) {
                        console.error('addIsosurface error:', e);
                    }
                }
            }
        }
        this.viewer.render();
    }

    resize() {
        if (this.viewer) {
            this.viewer.resize();
            this.viewer.render();
        }
    }

    loadTrajectory(traj) {
        if (!this.viewer || !traj || !traj.frames || !traj.frames.length) return;
        this.currentTraj = traj;
        this.symbols = traj.symbols || [];
        const tsIdx = traj.ts_frame_index || 0;

        try {
            this.viewer.removeAllModels();
            this.viewer.removeAllShapes();

            let multiSdf = '';
            traj.frames.forEach((frame, idx) => {
                const coords = frame.coords || frame.positions || [];
                const syms = frame.symbols || this.symbols;
                multiSdf += this.buildMolBlock(syms, coords);
                multiSdf += '$$$$\n';
            });

            this.viewer.addModelsAsFrames(multiSdf, 'sdf');
            this.viewer.setStyle({}, {
                stick: { radius: 0.12 * this.moleculeScale, doubleBondScaling: 0.5, tripleBondScaling: 0.3, colorscheme: 'Jmol' },
                sphere: { scale: 0.24 * this.moleculeScale, colorscheme: 'Jmol' }
            });
            this.viewer.zoomTo();
            this.hasFrames = true;
            this.setFrame(tsIdx);
        } catch (err) {
            console.warn('Multimodel frame load fallback:', err);
            this.hasFrames = false;
            this.setFrame(tsIdx);
        }
    }

    async setFrame(idx) {
        if (!this.viewer) return;
        this.currentFrameIdx = idx;
        if (this.hasFrames) {
            try {
                await this.viewer.setFrame(idx);
                this.viewer.render();
                return;
            } catch (err) {
                this.hasFrames = false;
            }
        }
        if (this.currentTraj && this.currentTraj.frames && this.currentTraj.frames[idx]) {
            const frame = this.currentTraj.frames[idx];
            const coords = frame.coords || frame.positions || [];
            const syms = frame.symbols || this.symbols;
            this.renderMolecule(syms, coords, frame.charges || null, false);
        }
    }

    /**
     * Superimposes Step K Product and Step K+1 Reactant
     * for Boundary Overlap Inspection.
     */
    renderOverlap(prevSymbols, prevCoords, nextSymbols, nextCoords, rmsd = 0, isAligned = true) {
        if (!this.viewer || !prevSymbols || !nextSymbols) return;
        this.viewer.clear();

        const molPrev = this.buildMolBlock(prevSymbols, prevCoords);
        const molNext = this.buildMolBlock(nextSymbols, nextCoords);

        const m1 = this.viewer.addModel(molPrev, 'mol');
        const m2 = this.viewer.addModel(molNext, 'mol');

        // Step K Product in Cyan (#38bdf8)
        this.viewer.setStyle({ model: m1 }, {
            stick: { color: '#38bdf8', radius: 0.12, doubleBondScaling: 0.5, tripleBondScaling: 0.3 },
            sphere: { color: '#38bdf8', scale: 0.22 }
        });

        // Step K+1 Reactant in Magenta (#f43f5e) if aligned, Yellow (#fbbf24) if unaligned
        const nextColor = isAligned ? '#f43f5e' : '#fbbf24';
        this.viewer.setStyle({ model: m2 }, {
            stick: { color: nextColor, radius: 0.12, doubleBondScaling: 0.5, tripleBondScaling: 0.3 },
            sphere: { color: nextColor, scale: 0.22 }
        });

        this.viewer.zoomTo();
        this.viewer.render();
        const modeLabel = isAligned ? 'Aligned' : 'Raw Unaligned';
        this.setStatus(`Boundary Overlap (${modeLabel}) | RMSD: ${rmsd.toFixed(4)} Å (Cyan: Step K Product, ${isAligned ? 'Magenta' : 'Yellow'}: Step K+1 Reactant)`);
    }

    renderFastFrame(symbols, coords, activeOrbitals, isovalue, frameIdx, charges = null, showCharges = false) {
        this.currentFrameIdx = frameIdx;
        this.renderMolecule(symbols, coords, charges, showCharges);

        if (!this.hasInitializedCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (activeOrbitals && activeOrbitals.length > 0) {
            this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
        }
    }

    async renderScene(options) {
        const {
            symbols, coords, charges, showCharges,
            activeOrbitals, isovalue, frameIdx,
            hideAllOrbitals, fetchCubeFn, resetCamera
        } = options;

        if (!this.viewer) return;
        this.currentFrameIdx = frameIdx;

        // 1. Immediately render molecular coordinates with authentic double bonds
        this.renderMolecule(symbols, coords, charges, showCharges);

        if (!this.hasInitializedCamera || resetCamera) {
            this.viewer.zoomTo();
            this.hasInitializedCamera = true;
        }

        if (hideAllOrbitals || !activeOrbitals || activeOrbitals.length === 0) {
            this.setStatus(`Frame ${frameIdx} | Ready`);
            return;
        }

        // 2. Render any already-cached orbital surfaces immediately
        this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);

        // 3. Find missing orbital cubes
        const missing = [];
        if (fetchCubeFn) {
            for (const orb of activeOrbitals) {
                const orbKey = orb.orbital_key || orb.key || orb;
                const cacheKey = `${frameIdx}_${orbKey}`;
                if (!this.cubeCache[cacheKey]) {
                    missing.push({ orb, orbKey, cacheKey });
                }
            }
        }

        if (missing.length > 0) {
            this.setStatus(`Frame ${frameIdx} | Computing orbital (${missing.map(m => m.orbKey).join(', ')})...`);
            const fetchErrors = [];
            try {
                await Promise.all(missing.map(async (m) => {
                    const res = await fetchCubeFn(frameIdx, m.orbKey, m.orb.chem_key || null, isovalue);
                    if (res && res.success && res.cube) {
                        this.cubeCache[m.cacheKey] = res.cube;
                    } else if (res && res.error) {
                        fetchErrors.push(res.error);
                    }
                }));
            } catch (err) {
                fetchErrors.push(err.message);
                console.error('fetchCube error:', err);
            }

            // 4. Redraw scene once missing cubes are fetched
            if (this.currentFrameIdx === frameIdx) {
                this.renderMolecule(symbols, coords, charges, showCharges);
                this.renderActiveIsosurfaces(activeOrbitals, isovalue, frameIdx);
                if (fetchErrors.length > 0) {
                    this.setStatus(`Frame ${frameIdx} | ${fetchErrors.join(' · ')}`, true);
                } else {
                    this.setStatus(`Frame ${frameIdx} | Ready`);
                }
            }
        } else {
            this.setStatus(`Frame ${frameIdx} | Ready`);
        }
    }

    clearCache() {
        this.cubeCache = {};
        this.vdoCache = {};
    }
}
