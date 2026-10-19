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

export function resolveEventTarget(el, wiring) {
    if (!el || el.tagName.toLowerCase() === 'svg') return null;

    let cur = el;
    while (cur && cur.tagName && cur.tagName.toLowerCase() !== 'svg') {
        const id = cur.getAttribute('id');
        const lbl = cur.getAttribute('inkscape:label') || cur.getAttribute('label');
        if ((id && wiring[id]) || (lbl && wiring[lbl]) ||
            cur.classList.contains('svg-clickable-arrow') ||
            cur.classList.contains('svg-clickable-state') ||
            cur.classList.contains('svg-clickable-molecule')) {
            return cur;
        }
        if (id && cur.tagName.toLowerCase() === 'g') {
            return cur;
        }
        cur = cur.parentElement;
    }
    return el;
}
