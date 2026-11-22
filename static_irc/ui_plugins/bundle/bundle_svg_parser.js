/**
 * SVG Parser & Element Resolution for Interactive Mechanisms.
 */

export function getElementLabel(el) {
    return el.getAttribute('inkscape:label') ||
           el.getAttributeNS('http://www.inkscape.org/namespaces/inkscape', 'label') ||
           el.getAttribute('label') ||
           el.getAttribute('data-label') ||
           null;
}

export function parseLabeledElements(rootSvg) {
    const tree = [];
    const seenIds = new Set();

    const traverse = (parentEl, depth = 0) => {
        const result = [];
        const children = Array.from(parentEl.children);

        children.forEach((child, idx) => {
            const tag = child.tagName.toLowerCase();
            if (tag === 'defs' || tag === 'sodipodi:namedview' || tag === 'metadata') {
                return;
            }

            const label = getElementLabel(child);
            const hasLabel = Boolean(label);
            const isGroup = (tag === 'g');
            let id = child.getAttribute('id');

            if (!id) {
                id = `elem_${depth}_${idx}`;
                child.setAttribute('id', id);
            }

            const subChildren = isGroup ? traverse(child, depth + 1) : [];

            if (hasLabel || (isGroup && subChildren.length > 0)) {
                if (!seenIds.has(id)) {
                    seenIds.add(id);
                    result.push({
                        id: id,
                        label: label || id,
                        tag: tag,
                        depth: depth,
                        hasChildren: subChildren.length > 0,
                        children: subChildren,
                        isExplicitLabel: hasLabel
                    });
                }
            } else if (subChildren.length > 0) {
                result.push(...subChildren);
            }
        });
        return result;
    };

    const topNodes = traverse(rootSvg, 0);

    const allLabeled = rootSvg.querySelectorAll('[inkscape\\:label], [label]');
    allLabeled.forEach(el => {
        if (el.closest('defs')) return;
        const id = el.getAttribute('id');
        const lbl = getElementLabel(el);
        if (id && lbl && !seenIds.has(id)) {
            seenIds.add(id);
            topNodes.push({
                id: id,
                label: lbl,
                tag: el.tagName.toLowerCase(),
                depth: 0,
                hasChildren: false,
                children: [],
                isExplicitLabel: true
            });
        }
    });

    return topNodes;
}

export function resolveEventTarget(el, wiring, mode = 'interactive') {
    if (!el || !el.tagName || el.tagName.toLowerCase() === 'svg') return null;

    // Resolve external hit-cushion rects attached for standalone <path> elements
    const hitFor = el.getAttribute && el.getAttribute('data-hit-for');
    if (hitFor) {
        const rootSvg = el.ownerSVGElement || el.closest('svg');
        if (rootSvg) {
            const ref = rootSvg.getElementById(hitFor) ||
                        rootSvg.querySelector(`[inkscape\\:label="${hitFor}"], [label="${hitFor}"]`);
            if (ref) el = ref;
        }
    }

    if (mode === 'interactive') {
        let cur = el;
        while (cur && cur.tagName && cur.tagName.toLowerCase() !== 'svg') {
            const id = cur.getAttribute('id');
            const lbl = getElementLabel(cur);
            const cfg = (id && wiring[id]) || (lbl && wiring[lbl]);
            if (cfg && cfg.role && cfg.role !== 'none' && cfg.role !== 'push_arrow') {
                return cur;
            }
            cur = cur.parentElement;
        }
        return null;
    }

    // Edit / Wiring mode: prefer an ancestor that is wired or explicitly labeled over unnamed inner Inkscape <g> sub-groups
    let cur = el;
    let fallbackGroup = null;
    while (cur && cur.tagName && cur.tagName.toLowerCase() !== 'svg') {
        const id = cur.getAttribute('id');
        const lbl = getElementLabel(cur);
        if ((id && wiring[id]) || (lbl && wiring[lbl]) || lbl) {
            return cur;
        }
        if (!fallbackGroup && id && cur.tagName.toLowerCase() === 'g') {
            fallbackGroup = cur;
        }
        cur = cur.parentElement;
    }
    return fallbackGroup || el;
}
