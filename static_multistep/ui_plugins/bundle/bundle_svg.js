/**
 * Interactive SVG Reaction Mechanism & Hierarchical Group Wiring Manager.
 */

import { parseLabeledElements, resolveEventTarget } from './bundle_svg_parser.js?v=2.0.3';
import { renderSvgPaperCard } from './bundle_svg_view.js?v=2.0.3';

export class SvgMechanismManager {
    constructor(containerEl, options = {}) {
        this.containerEl = containerEl;
        this.options = options;
        this.rawSvgText = '';
        this.wiring = options.initialWiring || {};
        this.groupsTree = [];
        this.selectedGroupId = null;
        this.mode = 'interactive'; // 'interactive' or 'edit'
        this.currentFrame = 0;
        this.playDirection = 'forward'; // 'forward' or 'reverse'
        this.zoomLevel = 1.0;

        // Callbacks
        this.onPlaySegment = options.onPlaySegment || (() => {});
        this.onJumpFrame = options.onJumpFrame || (() => {});
        this.onOpenMolecule = options.onOpenMolecule || (() => {});
        this.onSaveWiring = options.onSaveWiring || (() => {});
        this.onSaveSvg = options.onSaveSvg || (() => {});
    }

    loadSvg(svgText) {
        if (!svgText || !svgText.trim()) {
            return false;
        }

        // Ensure default SVG namespace is present so DOMParser treats it as an SVGSVGElement
        if (!/<svg[^>]+xmlns\s*=\s*['"]http:\/\/www\.w3\.org\/2000\/svg['"]/i.test(svgText)) {
            if (/<svg[^>]+xmlns\s*=/i.test(svgText)) {
                svgText = svgText.replace(/xmlns\s*=\s*['"][^'"]*['"]/i, 'xmlns="http://www.w3.org/2000/svg"');
            } else {
                svgText = svgText.replace(/<svg\b/i, '<svg xmlns="http://www.w3.org/2000/svg"');
            }
        }

        this.rawSvgText = svgText;

        const parser = new DOMParser();
        let doc = parser.parseFromString(svgText, "image/svg+xml");

        if (doc.querySelector("parsererror")) {
            console.warn("XML parser error, falling back to HTML parser for SVG extraction...");
            const htmlDoc = parser.parseFromString(svgText, "text/html");
            const svgInHtml = htmlDoc.querySelector("svg");
            if (svgInHtml) {
                doc = htmlDoc;
            } else {
                console.error("Failed to parse SVG document.");
                return false;
            }
        }

        const svgEl = doc.querySelector("svg");
        if (!svgEl) {
            console.error("No <svg> root element found.");
            return false;
        }

        // Normalize raw SVG text elements to match Inkscape rendering and prevent bold/overlapping text
        svgEl.querySelectorAll('text').forEach(t => {
            const st = t.getAttribute('style') || '';
            let newSt = st.replace(/stroke:\s*#[0-9a-fA-F]+;?/gi, 'stroke:none;');
            if (!/font-size/i.test(newSt) && !t.getAttribute('font-size')) {
                newSt += ';font-size:10.5px;';
            }
            if (!/font-family/i.test(newSt) && !t.getAttribute('font-family')) {
                newSt += ';font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;';
            }
            t.setAttribute('style', newSt);
        });

        this.groupsTree = parseLabeledElements(svgEl);
        this.renderSvgToContainer(svgEl);
        this.attachEventDelegation();
        this.updatePushArrows(this.currentFrame, this.playDirection);
        return true;
    }

    renderSvgToContainer(svgEl) {
        renderSvgPaperCard(this.containerEl, svgEl, {
            zoom: (factor) => this.zoom(factor),
            resetZoom: () => this.resetZoom()
        });
        this.applyClassesToWiredElements();
    }

    zoom(factor) {
        this.zoomLevel = Math.max(0.4, Math.min(4.0, this.zoomLevel * factor));
        const paper = this.containerEl.querySelector('#svg-paper-card');
        if (paper) {
            paper.style.transform = `scale(${this.zoomLevel})`;
        }
    }

    resetZoom() {
        this.zoomLevel = 1.0;
        const paper = this.containerEl.querySelector('#svg-paper-card');
        if (paper) {
            paper.style.transform = `scale(1.0)`;
        }
    }

    setMode(newMode) {
        this.mode = newMode;
        const paper = this.containerEl.querySelector('#svg-paper-card') || this.containerEl;
        if (paper) {
            paper.classList.toggle('mode-wiring', newMode === 'edit');
            paper.classList.toggle('mode-interactive', newMode === 'interactive');
        }
        this.applyClassesToWiredElements();
        this.updatePushArrows(this.currentFrame, this.playDirection);
    }

    applyClassesToWiredElements() {
        const svg = this.containerEl.querySelector('#active-mechanism-svg');
        if (!svg) return;

        svg.querySelectorAll('.svg-clickable-arrow, .svg-clickable-state, .svg-clickable-molecule, .svg-edit-selectable, .svg-wired-push-arrow').forEach(el => {
            el.classList.remove('svg-clickable-arrow', 'svg-clickable-state', 'svg-clickable-molecule', 'svg-edit-selectable', 'svg-wired-push-arrow');
        });

        Object.entries(this.wiring).forEach(([key, cfg]) => {
            const el = svg.getElementById(key) ||
                       svg.querySelector(`[inkscape\\:label="${key}"]`) ||
                       svg.querySelector(`[label="${key}"]`);
            if (!el) return;

            if (this.mode === 'interactive') {
                if (cfg.role === 'arrow') {
                    el.classList.add('svg-clickable-arrow');
                    el.title = `Reaction Arrow: Play frames ${cfg.startFrame || 0} → ${cfg.endFrame || 0} (${cfg.direction || 'forward'})`;
                } else if (cfg.role === 'state' || cfg.role === 'state_molecule') {
                    el.classList.add('svg-clickable-state');
                    if (cfg.role === 'state_molecule' || cfg.moleculeName) {
                        el.classList.add('svg-clickable-molecule');
                        el.title = `${cfg.moleculeName || 'Intermediate'} (Frame ${cfg.targetFrame || 0}): Click to jump, Right-click for orbitals`;
                    } else {
                        el.title = `State: Jump to frame ${cfg.targetFrame || 0}`;
                    }
                } else if (cfg.role === 'molecule') {
                    el.classList.add('svg-clickable-molecule');
                    el.title = `${cfg.moleculeName || 'Molecule'}: Right-click for orbitals`;
                } else if (cfg.role === 'push_arrow') {
                    el.classList.add('svg-wired-push-arrow');
                    el.title = `Push Arrow: Active frames ${cfg.startFrame || 0} → ${cfg.endFrame || 0} (${cfg.direction || 'both'})`;
                }
            } else {
                el.classList.add('svg-edit-selectable');
            }
        });
    }

    attachEventDelegation() {
        const svg = this.containerEl.querySelector('#active-mechanism-svg');
        if (!svg) return;

        svg.addEventListener('click', (e) => {
            const target = resolveEventTarget(e.target, this.wiring);
            if (!target) return;

            const grpId = target.getAttribute('id');
            const lbl = target.getAttribute('inkscape:label') || target.getAttribute('label') || grpId;
            const cfg = this.wiring[grpId] || this.wiring[lbl];

            if (this.mode === 'interactive') {
                if (!cfg) return;
                e.stopPropagation();

                if (cfg.role === 'arrow') {
                    const sf = cfg.startFrame !== undefined ? parseInt(cfg.startFrame, 10) : 0;
                    const ef = cfg.endFrame !== undefined ? parseInt(cfg.endFrame, 10) : 100;
                    const dir = cfg.direction || 'forward';
                    this.onPlaySegment(sf, ef, dir);
                } else if (cfg.role === 'state' || cfg.role === 'state_molecule') {
                    const tf = cfg.targetFrame !== undefined ? parseInt(cfg.targetFrame, 10) : 0;
                    this.onJumpFrame(tf);
                } else if (cfg.role === 'molecule') {
                    this.onOpenMolecule(grpId, cfg);
                }
            } else {
                e.stopPropagation();
                this.selectGroup(grpId || lbl);
            }
        });

        svg.addEventListener('contextmenu', (e) => {
            const target = resolveEventTarget(e.target, this.wiring);
            if (!target) return;

            const grpId = target.getAttribute('id');
            const lbl = target.getAttribute('inkscape:label') || target.getAttribute('label') || grpId;
            const cfg = this.wiring[grpId] || this.wiring[lbl];

            if (cfg && (cfg.role === 'molecule' || cfg.role === 'state_molecule' || (cfg.role === 'state' && cfg.moleculeName))) {
                e.preventDefault();
                e.stopPropagation();
                this.onOpenMolecule(grpId, cfg);
            }
        });
    }

    selectGroup(grpId) {
        this.selectedGroupId = grpId;
        const svg = this.containerEl.querySelector('#active-mechanism-svg');
        let lbl = grpId;
        if (svg) {
            svg.querySelectorAll('.svg-group-selected').forEach(el => el.classList.remove('svg-group-selected'));
            let el = svg.getElementById(grpId);
            if (!el) {
                for (const elem of svg.querySelectorAll('*')) {
                    if (elem.getAttribute('inkscape:label') === grpId || elem.getAttribute('label') === grpId) {
                        el = elem;
                        break;
                    }
                }
            }
            if (el) {
                el.classList.add('svg-group-selected');
                lbl = el.getAttribute('inkscape:label') || el.getAttribute('label') || grpId;
            }
        }
        const cfg = this.wiring[grpId] || this.wiring[lbl] || {};
        if (this.options.onGroupSelected) {
            this.options.onGroupSelected(grpId, cfg);
        }
    }

    updatePushArrows(frameIdx, direction) {
        this.currentFrame = frameIdx;
        this.playDirection = direction || 'forward';

        const svg = this.containerEl.querySelector('#active-mechanism-svg');
        // In wiring / edit mode, all arrows remain 100% visible and selectable
        if (this.mode === 'edit') {
            svg.querySelectorAll('.svg-push-inactive, .svg-push-active').forEach(el => {
                el.classList.remove('svg-push-inactive', 'svg-push-active');
            });
            return;
        }

        Object.entries(this.wiring).forEach(([key, cfg]) => {
            if (cfg.role !== 'push_arrow') return;

            let el = svg.getElementById(key);
            if (!el) {
                try {
                    el = svg.querySelector(`[inkscape\\:label="${key}"], [label="${key}"]`);
                } catch (e) {}
            }
            if (!el) {
                for (const elem of svg.querySelectorAll('*')) {
                    if (elem.getAttribute('inkscape:label') === key || elem.getAttribute('label') === key) {
                        el = elem;
                        break;
                    }
                }
            }
            if (!el) return;

            const sf = !isNaN(parseInt(cfg.startFrame, 10)) ? parseInt(cfg.startFrame, 10) : 0;
            const ef = !isNaN(parseInt(cfg.endFrame, 10)) ? parseInt(cfg.endFrame, 10) : 9999;
            const arrowDir = cfg.direction || 'both';

            const dirMatches = (arrowDir === 'both') || (arrowDir === this.playDirection);
            const minF = Math.min(sf, ef);
            const maxF = Math.max(sf, ef);
            const frameMatches = (frameIdx >= minF && frameIdx <= maxF);

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
        this.applyClassesToWiredElements();
        this.updatePushArrows(this.currentFrame, this.playDirection);
        this.onSaveWiring(this.wiring);
    }

    getExportSvg() {
        const svg = this.containerEl.querySelector('#active-mechanism-svg');
        return svg ? svg.outerHTML : '';
    }
}
