/**
 * Interactive SVG Reaction Mechanism & Hierarchical Group Wiring Manager.
 * Supports nested Inkscape groups (<g inkscape:label="..." id="...">),
 * dynamic electron-pushing arrow highlighting (direction-aware: forward vs reverse),
 * reaction arrow triggers, state jumps, and molecule right-click orbital bundles.
 */

export class SvgMechanismManager {
    constructor(containerEl, options = {}) {
        this.containerEl = containerEl;
        this.options = options;
        this.svgDoc = null;
        this.svgElement = null;
        this.wiring = options.initialWiring || {};
        this.groupsTree = [];
        this.selectedGroupId = null;
        this.mode = 'interactive'; // 'interactive' or 'edit'
        this.currentFrame = 0;
        this.playDirection = 'forward'; // 'forward' or 'reverse'

        // Callbacks
        this.onPlaySegment = options.onPlaySegment || (() => {});
        this.onJumpFrame = options.onJumpFrame || (() => {});
        this.onOpenMolecule = options.onOpenMolecule || (() => {});
        this.onSaveWiring = options.onSaveWiring || (() => {});
        this.onSaveSvg = options.onSaveSvg || (() => {});
    }

    getDefaultSvg() {
        return `<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" viewBox="0 0 920 280" width="100%" height="100%">
            <defs>
                <marker id="arrow-head-fwd" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
                    <polygon points="0 0, 10 3.5, 0 7" fill="#38bdf8"/>
                </marker>
                <marker id="arrow-head-rev" markerWidth="10" markerHeight="7" refX="9" refY="3.5" orient="auto">
                    <polygon points="0 0, 10 3.5, 0 7" fill="#f59e0b"/>
                </marker>
                <marker id="push-head-fwd" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto">
                    <polygon points="0 0, 8 3, 0 6" fill="#38bdf8"/>
                </marker>
                <marker id="push-head-rev" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto">
                    <polygon points="0 0, 8 3, 0 6" fill="#f59e0b"/>
                </marker>
                <filter id="glow-cyan" x="-20%" y="-20%" width="140%" height="140%">
                    <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#38bdf8" flood-opacity="0.9"/>
                </filter>
                <filter id="glow-amber" x="-20%" y="-20%" width="140%" height="140%">
                    <feDropShadow dx="0" dy="0" stdDeviation="4" flood-color="#f59e0b" flood-opacity="0.9"/>
                </filter>
            </defs>
            <style>
                text { font-family: 'Inter', -apple-system, sans-serif; }
                .state-rect { fill: rgba(30, 41, 59, 0.7); stroke: #475569; stroke-width: 1.5; rx: 8; transition: all 0.2s ease; }
                .state-rect:hover { stroke: #10b981; fill: rgba(16, 185, 129, 0.15); cursor: pointer; }
                .rxn-arrow { stroke: #38bdf8; stroke-width: 3; fill: none; transition: stroke-width 0.2s ease; cursor: pointer; }
                .rxn-arrow:hover { stroke-width: 5; filter: url(#glow-cyan); }
                .push-arrow-path { fill: none; stroke-width: 2.5; stroke-dasharray: 4,3; transition: all 0.25s ease; }
            </style>

            <!-- Reaction Step 1 Hierarchy Group -->
            <g id="group_step_1" inkscape:label="Step 1: Nucleophilic Addition">
                <!-- Reactants State -->
                <g id="group_reactants" inkscape:label="Reactants: CH3COCl + HCOO-">
                    <rect class="state-rect" x="30" y="70" width="220" height="140"/>
                    <text x="140" y="105" text-anchor="middle" fill="#f8fafc" font-size="14" font-weight="700">Reactants (Frame 0)</text>
                    <text x="140" y="135" text-anchor="middle" fill="#cbd5e1" font-size="13">CH₃C(=O)Cl + HCOO⁻</text>
                    <text x="140" y="165" text-anchor="middle" fill="#94a3b8" font-size="11">Acetyl chloride & Formate</text>
                    <text x="140" y="195" text-anchor="middle" fill="#10b981" font-size="11">ΔG = 0.0 kcal/mol</text>
                </g>

                <!-- Reaction Arrow Step 1 (Addition) -->
                <g id="group_arrow_step1" inkscape:label="Reaction Arrow: Addition (0 → 68)">
                    <path class="rxn-arrow" d="M 265 140 L 335 140" marker-end="url(#arrow-head-fwd)"/>
                    <text x="300" y="125" text-anchor="middle" fill="#38bdf8" font-size="11" font-weight="600">k₁ (TS 39)</text>
                </g>

                <!-- Forward Electron Pushing Arrow Step 1 -->
                <g id="group_push_fwd_step1" inkscape:label="Push Arrow: Formate O- attack on Carbonyl (Fwd)">
                    <path class="push-arrow-path" stroke="#38bdf8" marker-end="url(#push-head-fwd)" d="M 190 145 C 230 180, 270 170, 290 150"/>
                    <text x="245" y="185" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">O⁻ ➔ C=O</text>
                </g>

                <!-- Reverse Electron Pushing Arrow Step 1 -->
                <g id="group_push_rev_step1" inkscape:label="Push Arrow: Tetrahedral carbonyl reform / expel Formate (Rev)">
                    <path class="push-arrow-path" stroke="#f59e0b" marker-end="url(#push-head-rev)" d="M 330 155 C 310 185, 250 185, 200 155"/>
                    <text x="260" y="198" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="600">O⁻ ➔ C-O expel</text>
                </g>
            </g>

            <!-- Intermediate State Hierarchy Group -->
            <g id="group_intermediate" inkscape:label="Tetrahedral Intermediate">
                <rect class="state-rect" x="350" y="70" width="220" height="140"/>
                <text x="460" y="105" text-anchor="middle" fill="#f8fafc" font-size="14" font-weight="700">Tetrahedral Intermediate</text>
                <text x="460" y="135" text-anchor="middle" fill="#e2e8f0" font-size="13">CH₃C(O⁻)(Cl)(OCH=O)</text>
                <text x="460" y="165" text-anchor="middle" fill="#94a3b8" font-size="11">Frame 68-76 (Min 1)</text>
                <text x="460" y="195" text-anchor="middle" fill="#38bdf8" font-size="11">ΔG = -1.5 kcal/mol</text>
            </g>

            <!-- Reaction Step 2 Hierarchy Group -->
            <g id="group_step_2" inkscape:label="Step 2: Elimination">
                <!-- Reaction Arrow Step 2 (Elimination) -->
                <g id="group_arrow_step2" inkscape:label="Reaction Arrow: Elimination (77 → 146)">
                    <path class="rxn-arrow" d="M 585 140 L 655 140" marker-end="url(#arrow-head-fwd)"/>
                    <text x="620" y="125" text-anchor="middle" fill="#38bdf8" font-size="11" font-weight="600">k₂ (TS 107)</text>
                </g>

                <!-- Forward Electron Pushing Arrow Step 2 -->
                <g id="group_push_fwd_step2" inkscape:label="Push Arrow: Oxyanion collapse & Cl- expulsion (Fwd)">
                    <path class="push-arrow-path" stroke="#38bdf8" marker-end="url(#push-head-fwd)" d="M 520 145 C 560 180, 600 170, 620 150"/>
                    <text x="575" y="185" text-anchor="middle" fill="#38bdf8" font-size="10" font-weight="600">O⁻ ➔ C=O & Cl⁻</text>
                </g>

                <!-- Reverse Electron Pushing Arrow Step 2 -->
                <g id="group_push_rev_step2" inkscape:label="Push Arrow: Chloride re-attack on Anhydride (Rev)">
                    <path class="push-arrow-path" stroke="#f59e0b" marker-end="url(#push-head-rev)" d="M 660 155 C 640 185, 580 185, 530 155"/>
                    <text x="595" y="198" text-anchor="middle" fill="#f59e0b" font-size="10" font-weight="600">Cl⁻ ➔ Anhydride</text>
                </g>

                <!-- Products State -->
                <g id="group_products" inkscape:label="Products: Acetic-formic anhydride + Cl-">
                    <rect class="state-rect" x="670" y="70" width="220" height="140"/>
                    <text x="780" y="105" text-anchor="middle" fill="#f8fafc" font-size="14" font-weight="700">Products (Frame 146)</text>
                    <text x="780" y="135" text-anchor="middle" fill="#cbd5e1" font-size="13">CH₃COOCHO + Cl⁻</text>
                    <text x="780" y="165" text-anchor="middle" fill="#94a3b8" font-size="11">Anhydride + Chloride</text>
                    <text x="780" y="195" text-anchor="middle" fill="#10b981" font-size="11">ΔG = -12.7 kcal/mol</text>
                </g>
            </g>
        </svg>`;
    }

    getDefaultWiring() {
        return {
            "group_arrow_step1": { role: "arrow", startFrame: 0, endFrame: 68, direction: "forward" },
            "group_arrow_step2": { role: "arrow", startFrame: 77, endFrame: 146, direction: "forward" },
            "group_reactants": { role: "state", targetFrame: 0, moleculeName: "Reactants" },
            "group_intermediate": { role: "state", targetFrame: 68, moleculeName: "Tetrahedral Intermediate" },
            "group_products": { role: "state", targetFrame: 146, moleculeName: "Products" },
            "group_push_fwd_step1": { role: "push_arrow", startFrame: 0, endFrame: 45, direction: "forward" },
            "group_push_rev_step1": { role: "push_arrow", startFrame: 20, endFrame: 68, direction: "reverse" },
            "group_push_fwd_step2": { role: "push_arrow", startFrame: 77, endFrame: 115, direction: "forward" },
            "group_push_rev_step2": { role: "push_arrow", startFrame: 100, endFrame: 146, direction: "reverse" }
        };
    }

    loadSvg(svgText) {
        if (!svgText || !svgText.trim()) {
            svgText = this.getDefaultSvg();
            if (Object.keys(this.wiring).length === 0) {
                this.wiring = this.getDefaultWiring();
            }
        }

        const parser = new DOMParser();
        const doc = parser.parseFromString(svgText, "image/svg+xml");
        const svgEl = doc.querySelector("svg");
        if (!svgEl) {
            console.error("Invalid SVG document parsed.");
            return false;
        }

        this.svgDoc = doc;
        this.svgElement = svgEl;
        this.groupsTree = this.parseGroupTree(svgEl);

        this.renderSvgToContainer();
        this.wireInteractiveElements();
        this.updatePushArrows(this.currentFrame, this.playDirection);
        return true;
    }

    parseGroupTree(parentEl, depth = 0) {
        const groups = [];
        const children = Array.from(parentEl.children);

        children.forEach((child, idx) => {
            if (child.tagName.toLowerCase() === 'g') {
                const id = child.getAttribute('id') || `grp_${depth}_${idx}`;
                if (!child.getAttribute('id')) child.setAttribute('id', id);

                const label = child.getAttribute('inkscape:label') ||
                              child.getAttribute('label') ||
                              child.getAttribute('id') ||
                              `Group ${idx + 1}`;

                const subChildren = this.parseGroupTree(child, depth + 1);
                const pathCount = child.querySelectorAll('path, rect, circle, polygon, text').length;

                groups.push({
                    id: id,
                    label: label,
                    depth: depth,
                    pathCount: pathCount,
                    hasChildren: subChildren.length > 0,
                    children: subChildren,
                    element: child
                });
            }
        });
        return groups;
    }

    renderSvgToContainer() {
        if (!this.containerEl || !this.svgElement) return;
        this.containerEl.innerHTML = '';
        const cloned = this.svgElement.cloneNode(true);
        cloned.style.width = '100%';
        cloned.style.height = '100%';
        cloned.style.maxHeight = '100%';
        cloned.style.display = 'block';
        cloned.id = 'active-mechanism-svg';
        this.containerEl.appendChild(cloned);
    }

    setMode(newMode) {
        this.mode = newMode;
        this.wireInteractiveElements();
    }

    wireInteractiveElements() {
        const svg = this.containerEl.querySelector('svg');
        if (!svg) return;

        // Clear existing classes
        svg.querySelectorAll('.svg-clickable-arrow, .svg-clickable-state, .svg-clickable-molecule, .svg-edit-selectable').forEach(el => {
            el.classList.remove('svg-clickable-arrow', 'svg-clickable-state', 'svg-clickable-molecule', 'svg-edit-selectable');
        });

        const allGroups = svg.querySelectorAll('g');
        allGroups.forEach(g => {
            const grpId = g.getAttribute('id');
            const cfg = this.wiring[grpId];

            if (this.mode === 'interactive') {
                if (cfg) {
                    if (cfg.role === 'arrow') {
                        g.classList.add('svg-clickable-arrow');
                        g.style.cursor = 'pointer';
                        g.title = `Reaction Arrow: Play frames ${cfg.startFrame || 0} → ${cfg.endFrame || 0} (${cfg.direction || 'forward'})`;
                    } else if (cfg.role === 'state') {
                        g.classList.add('svg-clickable-state');
                        g.style.cursor = 'pointer';
                        g.title = `State: Jump to frame ${cfg.targetFrame || 0}`;
                    } else if (cfg.role === 'molecule') {
                        g.classList.add('svg-clickable-molecule');
                        g.style.cursor = 'pointer';
                        g.title = `Molecule: Right-click to inspect orbitals`;
                    }
                }
            } else {
                // Edit mode
                g.classList.add('svg-edit-selectable');
                g.style.cursor = 'crosshair';
            }

            // Remove prior listeners by cloning or replaceWith
            const newG = g.cloneNode(true);
            g.parentNode.replaceChild(newG, g);

            // Re-bind click
            newG.addEventListener('click', (e) => {
                e.stopPropagation();
                this.handleGroupClick(grpId, newG);
            });

            // Re-bind right-click (context menu)
            newG.addEventListener('contextmenu', (e) => {
                if (this.mode === 'interactive') {
                    const c = this.wiring[grpId];
                    if (c && c.role === 'molecule') {
                        e.preventDefault();
                        e.stopPropagation();
                        this.onOpenMolecule(grpId, c);
                    }
                }
            });
        });
    }

    handleGroupClick(grpId, el) {
        if (this.mode === 'interactive') {
            const cfg = this.wiring[grpId];
            if (!cfg) return;

            if (cfg.role === 'arrow') {
                const sf = cfg.startFrame !== undefined ? parseInt(cfg.startFrame, 10) : 0;
                const ef = cfg.endFrame !== undefined ? parseInt(cfg.endFrame, 10) : 100;
                const dir = cfg.direction || 'forward';
                this.onPlaySegment(sf, ef, dir);
            } else if (cfg.role === 'state') {
                const tf = cfg.targetFrame !== undefined ? parseInt(cfg.targetFrame, 10) : 0;
                this.onJumpFrame(tf);
            } else if (cfg.role === 'molecule') {
                this.onOpenMolecule(grpId, cfg);
            }
        } else {
            // Edit wiring mode: select group
            this.selectGroup(grpId);
        }
    }

    selectGroup(grpId) {
        this.selectedGroupId = grpId;
        const svg = this.containerEl.querySelector('svg');
        if (svg) {
            svg.querySelectorAll('.svg-group-selected').forEach(el => el.classList.remove('svg-group-selected'));
            const el = svg.getElementById(grpId);
            if (el) el.classList.add('svg-group-selected');
        }
        if (this.options.onGroupSelected) {
            this.options.onGroupSelected(grpId, this.wiring[grpId]);
        }
    }

    updatePushArrows(frameIdx, direction) {
        this.currentFrame = frameIdx;
        this.playDirection = direction || 'forward';

        const svg = this.containerEl.querySelector('svg');
        if (!svg) return;

        Object.entries(this.wiring).forEach(([grpId, cfg]) => {
            if (cfg.role !== 'push_arrow') return;

            const el = svg.getElementById(grpId);
            if (!el) return;

            const sf = cfg.startFrame !== undefined ? parseInt(cfg.startFrame, 10) : 0;
            const ef = cfg.endFrame !== undefined ? parseInt(cfg.endFrame, 10) : 9999;
            const arrowDir = cfg.direction || 'both';

            const dirMatches = (arrowDir === 'both') || (arrowDir === this.playDirection);
            const frameMatches = (frameIdx >= sf && frameIdx <= ef);

            if (dirMatches && frameMatches) {
                el.classList.add('svg-push-active');
                el.classList.remove('svg-push-inactive');
            } else {
                el.classList.add('svg-push-inactive');
                el.classList.remove('svg-push-active');
            }
        });
    }

    setGroupWiring(grpId, roleConfig) {
        if (!roleConfig || roleConfig.role === 'none') {
            delete this.wiring[grpId];
        } else {
            this.wiring[grpId] = { ...roleConfig };
        }
        this.wireInteractiveElements();
        this.updatePushArrows(this.currentFrame, this.playDirection);
        this.onSaveWiring(this.wiring);
    }

    getExportSvg() {
        const svg = this.containerEl.querySelector('svg');
        return svg ? svg.outerHTML : '';
    }
}
