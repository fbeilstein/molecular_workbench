/**
 * viewer_svg.js
 * Interactive SVG Reaction Mechanism Manager.
 * Step 5 mechanism interactions:
 * - Click forward/reverse reaction arrows to play segments
 * - Click states/intermediates to jump to frame
 * - Click or right-click molecules to open molecule inspection
 * - Dynamic electron-pushing arrow highlighting & pulsing on playback
 */

class ViewerSvgMechanism {
  constructor(containerElement, options = {}) {
    this.container = containerElement;
    this.wiring = options.wiring || {};
    this.onPlaySegment = options.onPlaySegment || (() => {});
    this.onJumpFrame = options.onJumpFrame || (() => {});
    this.onOpenMolecule = options.onOpenMolecule || (() => {});
    this.currentFrame = 0;
    this.playDirection = 'forward';
  }

  setMechanism(svgText, wiring = {}) {
    this.wiring = wiring || {};
    if (!this.container || !svgText) return;

    // Ensure proper SVG namespace and bottom-anchored aspect ratio
    let text = svgText;
    if (!/<svg[^>]+xmlns\s*=\s*['"]http:\/\/www\.w3\.org\/2000\/svg['"]/i.test(text)) {
      text = text.replace(/<svg\b/i, '<svg xmlns="http://www.w3.org/2000/svg"');
    }
    // Anchor SVG to bottom so wide or constrained mechanism images sit directly above the playback controls
    if (!/preserveAspectRatio=/i.test(text)) {
      text = text.replace(/<svg\b/i, '<svg preserveAspectRatio="xMidYMax meet"');
    } else {
      text = text.replace(/preserveAspectRatio=["'][^"']*["']/i, 'preserveAspectRatio="xMidYMax meet"');
    }

    // Convert black strokes/fills to crisp white #f1f5f9 directly in SVG markup.
    // This preserves colored arrows/heteroatoms while allowing vibrant, unclashed CSS hover effects.
    text = text
      .replace(/stroke:\s*(#000000|#000|black|rgb\(0,\s*0,\s*0\))/gi, 'stroke:#f1f5f9')
      .replace(/fill:\s*(#000000|#000|black|rgb\(0,\s*0,\s*0\))/gi, 'fill:#f1f5f9')
      .replace(/stroke=["'](#000000|#000|black|rgb\(0,\s*0,\s*0\))["']/gi, 'stroke="#f1f5f9"')
      .replace(/fill=["'](#000000|#000|black|rgb\(0,\s*0,\s*0\))["']/gi, 'fill="#f1f5f9"');

    // Fix Inkscape text glyph groups: ensure g[aria-label] and g[id^="text"] have solid white fill and no stroke outline
    text = text.replace(/<g([^>]*aria-label=[^>]*)>/gi, (match) => {
      if (/style=["']/i.test(match)) {
        return match.replace(/style=["']([^"']*)["']/i, (m, style) => {
          let s = style.replace(/stroke:[^;]+;?/gi, '').replace(/stroke-width:[^;]+;?/gi, '').replace(/fill:[^;]+;?/gi, '');
          return `style="${s};fill:#f1f5f9;stroke:none;"`;
        });
      } else {
        return match.replace(/<g\b/i, '<g style="fill:#f1f5f9;stroke:none;" ');
      }
    });

    // Strip stroke and ensure fill on any <text> and <tspan> elements
    text = text.replace(/<(text|tspan)([^>]*)>/gi, (match, tag, attrs) => {
      if (/style=["']/i.test(attrs)) {
        const newAttrs = attrs.replace(/style=["']([^"']*)["']/i, (m, style) => {
          let s = style.replace(/stroke:[^;]+;?/gi, '').replace(/stroke-width:[^;]+;?/gi, '').replace(/fill:[^;]+;?/gi, '');
          return `style="${s};fill:#f1f5f9;stroke:none;"`;
        });
        return `<${tag}${newAttrs}>`;
      } else {
        return `<${tag} style="fill:#f1f5f9;stroke:none;"${attrs}>`;
      }
    });

    this.container.innerHTML = text;
    this.applyClasses();
    this.bindEvents();
    this.updatePushArrows(this.currentFrame, this.playDirection);
  }

  applyClasses() {
    const svg = this.container.querySelector('svg');
    if (!svg) return;

    for (const [key, cfg] of Object.entries(this.wiring)) {
      let el = this.findElement(svg, key);
      if (!el) continue;

      // Ensure standalone paths are wrapped in a group so hit cushions can be nested cleanly
      if (el.tagName && el.tagName.toLowerCase() === 'path') {
        const g = document.createElementNS('http://www.w3.org/2000/svg', 'g');
        const id = el.getAttribute('id');
        if (id) g.setAttribute('id', id);
        const lbl = el.getAttribute('inkscape:label');
        if (lbl) g.setAttribute('inkscape:label', lbl);
        el.removeAttribute('id');
        el.removeAttribute('inkscape:label');
        el.parentNode.insertBefore(g, el);
        g.appendChild(el);
        el = g;
      }

      if (cfg.role === 'arrow') {
        el.classList.add('svg-clickable-arrow');
        const dir = cfg.direction === 'reverse' ? '⟵' : '⟶';
        el.setAttribute('title', `Reaction Arrow: Play frames ${cfg.startFrame || 0} ${dir} ${cfg.endFrame || 0}`);
        this.addArrowHitCushion(el);
      } else if (cfg.role === 'state' || cfg.role === 'state_molecule') {
        // If this is a pure 'state' container group whose children are individually wired as molecules,
        // do not add an enclosing bounding cushion that would cause overlapping hover rectangles.
        let isContainer = false;
        if (cfg.role === 'state') {
          isContainer = Object.keys(this.wiring).some(childKey => {
            if (childKey === key) return false;
            const childCfg = this.wiring[childKey];
            if (!childCfg || (childCfg.role !== 'molecule' && childCfg.role !== 'state_molecule')) return false;
            const childEl = this.findElement(el, childKey);
            return childEl && childEl !== el;
          });
        }

        if (!isContainer) {
          el.classList.add('svg-clickable-state');
          if (cfg.role === 'state_molecule' || cfg.moleculeName) {
            el.classList.add('svg-clickable-molecule');
            el.setAttribute('title', `${cfg.moleculeName || 'Intermediate'} (Frame ${cfg.targetFrame || 0}): Click to jump, Double-click for molecule`);
          } else {
            el.setAttribute('title', `State: Jump to frame ${cfg.targetFrame || 0}`);
          }
          this.addBoundingHitCushion(el);
        }
      } else if (cfg.role === 'molecule') {
        el.classList.add('svg-clickable-molecule');
        const stateTarget = this.findStateTarget(el.parentElement, svg);
        const stateInfo = (stateTarget && stateTarget.cfg && stateTarget.cfg.targetFrame !== undefined)
          ? ` (Frame ${stateTarget.cfg.targetFrame}): Click to jump, `
          : ': ';
        el.setAttribute('title', `${cfg.moleculeName || 'Molecule'}${stateInfo}Double-click to inspect`);
        this.addBoundingHitCushion(el);
      }
    }
  }

  addArrowHitCushion(el) {
    if (el.querySelector('.svg-hit-cushion')) return;

    // 1. Add generous rectangular hit cushion centered on the reaction arrow shaft
    const paths = Array.from(el.querySelectorAll('path'));
    let shaftPath = paths.find(p => {
      try {
        const b = p.getBBox();
        return b.width > 20;
      } catch (_) {
        return false;
      }
    });

    if (!shaftPath && paths.length > 0) {
      shaftPath = paths[0];
    }

    if (shaftPath) {
      try {
        const b = shaftPath.getBBox();
        if (b && b.width > 0) {
          const padX = Math.max(2.5, b.width * 0.03);
          const padY = Math.max(4.0, (b.height || 0) * 0.6);
          const rectH = Math.max(8.0, (b.height || 0) + padY * 2);
          const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
          rect.setAttribute('x', (b.x - padX).toFixed(2));
          rect.setAttribute('y', (b.y - padY).toFixed(2));
          rect.setAttribute('width', (b.width + padX * 2).toFixed(2));
          rect.setAttribute('height', rectH.toFixed(2));
          rect.setAttribute('rx', '3.5');
          rect.setAttribute('ry', '3.5');
          rect.setAttribute('fill', 'rgba(0,0,0,0.001)');
          rect.setAttribute('stroke', 'none');
          rect.setAttribute('pointer-events', 'all');
          rect.classList.add('svg-hit-cushion');
          el.insertBefore(rect, el.firstChild);
        }
      } catch (_) {}
    }

    // 2. Also add stroke hit cushions along paths as robust fallback
    paths.forEach(path => {
      const ghost = path.cloneNode(false);
      ghost.removeAttribute('id');
      ghost.removeAttribute('inkscape:label');
      ghost.removeAttribute('marker-end');
      ghost.removeAttribute('marker-start');
      ghost.removeAttribute('style'); // Strip inline style so stroke-width actually applies!
      ghost.setAttribute('style', 'fill:none !important; stroke:rgba(0,0,0,0.001) !important; stroke-width:7.0 !important; stroke-linecap:round !important; stroke-linejoin:round !important; pointer-events:stroke !important;');
      ghost.classList.add('svg-hit-cushion');
      el.insertBefore(ghost, el.firstChild);
    });
  }

  addBoundingHitCushion(el) {
    if (el.querySelector('.svg-hit-cushion')) return;
    try {
      const bbox = el.getBBox();
      if (bbox && bbox.width > 0 && bbox.height > 0) {
        const padX = 3.5; // ~16px padding in SVG mm units
        const padY = 3.5;
        const rect = document.createElementNS('http://www.w3.org/2000/svg', 'rect');
        rect.setAttribute('x', (bbox.x - padX).toFixed(2));
        rect.setAttribute('y', (bbox.y - padY).toFixed(2));
        rect.setAttribute('width', (bbox.width + padX * 2).toFixed(2));
        rect.setAttribute('height', (bbox.height + padY * 2).toFixed(2));
        rect.setAttribute('rx', '4.0');
        rect.setAttribute('ry', '4.0');
        rect.setAttribute('fill', 'rgba(0,0,0,0.001)');
        rect.setAttribute('stroke', 'none');
        rect.setAttribute('pointer-events', 'all');
        rect.classList.add('svg-hit-cushion');
        el.insertBefore(rect, el.firstChild);
        return;
      }
    } catch (_) {}

    // Fallback if getBBox() unavailable: add ghost stroke cushions to all child paths
    el.querySelectorAll('path').forEach(path => {
      const ghost = path.cloneNode(false);
      ghost.removeAttribute('id');
      ghost.removeAttribute('marker-end');
      ghost.removeAttribute('style');
      ghost.setAttribute('style', 'fill:none !important; stroke:rgba(0,0,0,0.001) !important; stroke-width:6.0 !important; pointer-events:stroke !important;');
      ghost.classList.add('svg-hit-cushion');
      el.insertBefore(ghost, el.firstChild);
    });
  }

  findStateTarget(el, rootSvg) {
    let cur = el;
    while (cur && cur !== rootSvg && cur.tagName && cur.tagName.toLowerCase() !== 'svg') {
      const id = cur.getAttribute('id');
      const lbl = cur.getAttribute('inkscape:label') || cur.getAttribute('label');
      const cfg = (id && this.wiring[id]) || (lbl && this.wiring[lbl]);
      if (cfg && (cfg.role === 'state' || cfg.role === 'state_molecule' || cfg.targetFrame !== undefined)) {
        return { el: cur, id: id || lbl, cfg };
      }
      cur = cur.parentElement;
    }
    return null;
  }

  bindEvents() {
    const svg = this.container.querySelector('svg');
    if (!svg) return;

    // Single click handler: Arrows (play segment) & States (jump frame)
    svg.addEventListener('click', (e) => {
      const target = this.resolveTarget(e.target, svg);
      if (!target) return;

      const { cfg } = target;
      if (!cfg) return;

      if (cfg.role === 'arrow') {
        e.stopPropagation();
        const sf = cfg.startFrame !== undefined ? parseInt(cfg.startFrame, 10) : 0;
        const ef = cfg.endFrame !== undefined ? parseInt(cfg.endFrame, 10) : 100;
        const dir = cfg.direction || 'forward';
        this.onPlaySegment(sf, ef, dir);
      } else if (cfg.role === 'state' || cfg.role === 'state_molecule' || cfg.targetFrame !== undefined) {
        e.stopPropagation();
        const tf = cfg.targetFrame !== undefined ? parseInt(cfg.targetFrame, 10) : 0;
        this.onJumpFrame(tf);
      } else {
        // e.g. clicked on a molecule (reactants or products) belonging to a state
        const stateTarget = this.findStateTarget(e.target, svg);
        if (stateTarget && stateTarget.cfg && stateTarget.cfg.targetFrame !== undefined) {
          e.stopPropagation();
          const tf = parseInt(stateTarget.cfg.targetFrame, 10) || 0;
          this.onJumpFrame(tf);
        }
      }
    });

    // Double-click handler: Open molecule inspection modal
    svg.addEventListener('dblclick', (e) => {
      const target = this.resolveTarget(e.target, svg);
      if (!target) return;

      const { id, cfg } = target;
      if (cfg && (cfg.role === 'molecule' || cfg.role === 'state_molecule' || cfg.moleculeName)) {
        e.stopPropagation();
        e.preventDefault();
        this.onOpenMolecule(id, cfg);
      }
    });

    // Right-click contextmenu handler
    svg.addEventListener('contextmenu', (e) => {
      const target = this.resolveTarget(e.target, svg);
      if (!target) return;

      const { id, cfg } = target;
      if (cfg && (cfg.role === 'molecule' || cfg.role === 'state_molecule' || cfg.moleculeName)) {
        e.preventDefault();
        e.stopPropagation();
        this.onOpenMolecule(id, cfg);
      }
    });
  }

  updatePushArrows(frameIdx, direction = 'forward') {
    this.currentFrame = frameIdx;
    this.playDirection = direction || 'forward';

    const svg = this.container ? (this.container.querySelector('#active-mechanism-svg') || this.container.querySelector('svg')) : null;
    if (!svg) return;

    Object.entries(this.wiring).forEach(([key, cfg]) => {
      if (cfg.role !== 'push_arrow') return;

      const el = this.findElement(svg, key);
      if (!el) return;

      const sf = cfg.startFrame !== undefined ? parseInt(cfg.startFrame, 10) : 0;
      const ef = cfg.endFrame !== undefined ? parseInt(cfg.endFrame, 10) : 9999;
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

  highlightCurrentState(frameIdx) {
    // User requested NO persistent highlights on states/molecules when mouse is not hovering.
    // Reaction steps can contain multiple molecules and states; only show interactive highlights on hover.
    const svg = this.container ? this.container.querySelector('svg') : null;
    if (!svg) return;
    svg.querySelectorAll('.svg-active-element').forEach(el => el.classList.remove('svg-active-element'));
  }

  findElement(svg, key) {
    if (!svg) return null;
    let el = null;
    if (typeof svg.getElementById === 'function') {
      el = svg.getElementById(key);
    }
    if (!el) {
      try {
        el = svg.querySelector(`[id="${key}"], [inkscape\\:label="${key}"], [label="${key}"], [data-id="${key}"]`);
      } catch (e) {}
    }
    if (!el) {
      for (const elem of svg.querySelectorAll('*')) {
        if (elem.getAttribute('id') === key || elem.getAttribute('inkscape:label') === key || elem.getAttribute('label') === key) {
          return elem;
        }
      }
    }
    return el;
  }

  resolveTarget(el, rootSvg) {
    let cur = el;
    while (cur && cur !== rootSvg && cur.tagName && cur.tagName.toLowerCase() !== 'svg') {
      const id = cur.getAttribute('id');
      const lbl = cur.getAttribute('inkscape:label') || cur.getAttribute('label');
      const cfg = (id && this.wiring[id]) || (lbl && this.wiring[lbl]);
      if (cfg) {
        return { el: cur, id: id || lbl, cfg };
      }
      cur = cur.parentElement;
    }
    return null;
  }
}

window.ViewerSvgMechanism = ViewerSvgMechanism;
