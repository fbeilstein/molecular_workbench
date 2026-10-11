/**
 * Mapping UI Plugin: Interactive 2D atom-to-atom reaction mapping canvas.
 * Allows visual alignment, drag-and-drop fragment repositioning, and click-to-wire atom pairing.
 */

function createAtomCircle(side, fragIdx, atomIdx, pt, onAtomClick) {
    const circle = document.createElementNS('http://www.w3.org/2000/svg', 'circle');
    circle.setAttribute('cx', pt.x);
    circle.setAttribute('cy', pt.y);
    circle.setAttribute('r', '15');
    circle.setAttribute('fill', 'transparent');
    circle.setAttribute('stroke', 'transparent');
    circle.style.cursor = 'pointer';
    circle.style.pointerEvents = 'all';
    circle.id = `atom_${side}_${fragIdx}_${atomIdx}`;
    circle.dataset.side = side;
    circle.dataset.key = `${side}_${fragIdx}_${atomIdx}`;

    circle.addEventListener('mouseover', () => circle.setAttribute('stroke', '#3b82f6'));
    circle.addEventListener('mouseout', () => circle.setAttribute('stroke', 'transparent'));
    circle.addEventListener('click', () => onAtomClick(circle, circle.dataset.key));

    return circle;
}

function createFragCard(side, svgStr, fragIdx, coords, startX, getScale, onAtomClick, onDrag) {
    const wrapper = document.createElement('div');
    wrapper.style.position = 'absolute';
    wrapper.style.top = '0';
    wrapper.style.left = '0';
    wrapper.style.cursor = 'grab';
    wrapper.style.borderRadius = '8px';
    wrapper.style.padding = '10px';
    wrapper.style.boxShadow = '0 2px 8px rgba(0,0,0,0.15)';
    wrapper.style.zIndex = '10';

    if (side === 'R') {
        wrapper.style.backgroundColor = 'rgba(235, 245, 255, 0.92)';
        wrapper.style.border = '2px solid #a3c2f0';
    } else {
        wrapper.style.backgroundColor = 'rgba(255, 240, 235, 0.92)';
        wrapper.style.border = '2px solid #f0bda3';
    }

    wrapper.innerHTML = svgStr;

    let fx = startX;
    let fy = 50 + fragIdx * 300;
    wrapper.style.transform = `translate(${fx}px, ${fy}px)`;
    let isDragging = false;

    wrapper.addEventListener('mousedown', (e) => {
        if (e.target.tagName === 'circle') return;
        e.stopPropagation();
        isDragging = true;
        wrapper.style.cursor = 'grabbing';
    });

    window.addEventListener('mousemove', (e) => {
        if (!isDragging) return;
        const scale = getScale();
        fx += e.movementX / scale;
        fy += e.movementY / scale;
        wrapper.style.transform = `translate(${fx}px, ${fy}px)`;
        onDrag();
    });

    window.addEventListener('mouseup', () => {
        if (!isDragging) return;
        isDragging = false;
        wrapper.style.cursor = 'grab';
    });

    const svgEl = wrapper.querySelector('svg');
    if (svgEl && coords) {
        Object.entries(coords).forEach(([atomIdx, pt]) => {
            svgEl.appendChild(createAtomCircle(side, fragIdx, atomIdx, pt, onAtomClick));
        });
    }

    return wrapper;
}

function createWirePath(rKey, pKey, rEl, pEl, zoomRect, scale, isApproved, onToggleApproved, onDeleteMapping) {
    const rRect = rEl.getBoundingClientRect();
    const pRect = pEl.getBoundingClientRect();

    const x1 = (rRect.left + rRect.width / 2 - zoomRect.left) / scale;
    const y1 = (rRect.top + rRect.height / 2 - zoomRect.top) / scale;
    const x2 = (pRect.left + pRect.width / 2 - zoomRect.left) / scale;
    const y2 = (pRect.top + pRect.height / 2 - zoomRect.top) / scale;

    const line = document.createElementNS('http://www.w3.org/2000/svg', 'path');
    line.setAttribute('d', `M ${x1} ${y1} C ${x1 + 100} ${y1} ${x2 - 100} ${y2} ${x2} ${y2}`);
    line.setAttribute('stroke', isApproved ? 'rgba(0,180,0,0.5)' : 'rgba(255,165,0,0.7)');
    line.setAttribute('stroke-width', '4');
    line.setAttribute('fill', 'none');
    line.style.pointerEvents = 'all';
    line.style.cursor = 'pointer';

    line.addEventListener('mouseover', () => {
        line.setAttribute('stroke-width', '6');
        rEl.setAttribute('stroke', 'red');
        rEl.setAttribute('stroke-width', '3');
        pEl.setAttribute('stroke', 'red');
        pEl.setAttribute('stroke-width', '3');
    });

    line.addEventListener('mouseout', () => {
        line.setAttribute('stroke-width', '4');
        rEl.setAttribute('stroke', 'transparent');
        pEl.setAttribute('stroke', 'transparent');
    });

    line.addEventListener('click', (e) => {
        if (e.shiftKey) {
            onDeleteMapping(rKey);
            return;
        }
        onToggleApproved(rKey);
    });

    line.addEventListener('contextmenu', (e) => {
        e.preventDefault();
        onDeleteMapping(rKey);
    });

    return line;
}

function setupPanZoom(container, zoomContainer, onPanZoom) {
    let scale = 1, panX = 0, panY = 0, isPanning = false;

    container.addEventListener('wheel', (e) => {
        e.preventDefault();
        const factor = e.deltaY > 0 ? 0.9 : 1.1;
        const rect = container.getBoundingClientRect();
        const mx = e.clientX - rect.left, my = e.clientY - rect.top;

        panX = mx - (mx - panX) * factor;
        panY = my - (my - panY) * factor;
        scale *= factor;
        zoomContainer.style.transform = `translate(${panX}px, ${panY}px) scale(${scale})`;
        onPanZoom(scale);
    });

    container.addEventListener('mousedown', (e) => {
        if (e.target.tagName === 'circle' || e.target.tagName === 'path') return;
        isPanning = true;
        container.style.cursor = 'move';
    });

    window.addEventListener('mousemove', (e) => {
        if (!isPanning) return;
        panX += e.movementX;
        panY += e.movementY;
        zoomContainer.style.transform = `translate(${panX}px, ${panY}px) scale(${scale})`;
        onPanZoom(scale);
    });

    window.addEventListener('mouseup', () => {
        if (!isPanning) return;
        isPanning = false;
        container.style.cursor = 'default';
    });

    return () => scale;
}

function connectAtoms(rk, pk) {
    for (let k in window.currentMapping) {
        if (k !== rk && window.currentMapping[k] !== pk) continue;
        delete window.currentMapping[k];
        delete window.approvedMappings[k];
    }
    window.currentMapping[rk] = pk;
}

function updateAtomHighlight(rEl, pEl, isApproved) {
    const color = isApproved ? 'rgba(0,255,0,0.1)' : 'transparent';
    if (rEl) rEl.setAttribute('fill', color);
    if (pEl) pEl.setAttribute('fill', color);
}

function renderMappingUI(schemaContainer, ui, state, step_id) {
    window.currentMapping = ui.initial_mapping ? Object.assign({}, ui.initial_mapping) : {};
    window.approvedMappings = {};
    let selectedNode = null;
    let hideApproved = false;

    const btnBar = document.createElement('div');
    btnBar.style.marginBottom = '10px';

    const toggleBtn = document.createElement('button');
    toggleBtn.innerText = "Hide Approved Wires";
    toggleBtn.onclick = () => {
        hideApproved = !hideApproved;
        toggleBtn.innerText = hideApproved ? "Show Approved Wires" : "Hide Approved Wires";
        drawWires();
    };

    const approveAllBtn = document.createElement('button');
    approveAllBtn.innerText = "Approve All Remaining";
    approveAllBtn.style.marginLeft = '10px';
    approveAllBtn.onclick = () => {
        for (let k in window.currentMapping) window.approvedMappings[k] = true;
        drawWires();
    };

    btnBar.appendChild(toggleBtn);
    btnBar.appendChild(approveAllBtn);
    schemaContainer.appendChild(btnBar);

    const container = document.createElement('div');
    container.style.position = 'relative';
    container.style.width = '100%';
    container.style.minHeight = '600px';
    container.style.overflow = 'hidden';
    container.style.border = '1px solid #ccc';
    container.style.backgroundColor = '#f4f4f4';

    const zoomContainer = document.createElement('div');
    zoomContainer.style.position = 'absolute';
    zoomContainer.style.top = '0';
    zoomContainer.style.left = '0';
    zoomContainer.style.width = '100%';
    zoomContainer.style.height = '100%';
    zoomContainer.style.transformOrigin = '0 0';

    const overlay = document.createElementNS('http://www.w3.org/2000/svg', 'svg');
    overlay.style.position = 'absolute';
    overlay.style.top = '0';
    overlay.style.left = '0';
    overlay.style.width = '100%';
    overlay.style.height = '100%';
    overlay.style.overflow = 'visible';
    overlay.style.pointerEvents = 'none';
    overlay.style.zIndex = '100';
    zoomContainer.appendChild(overlay);

    const getScale = setupPanZoom(container, zoomContainer, () => drawWires());

    function handleAtomClick(circle, key) {
        if (!selectedNode) {
            selectedNode = { key, el: circle };
            circle.setAttribute('fill', 'rgba(0,0,255,0.2)');
            return;
        }

        if (selectedNode.key.startsWith('R') && key.startsWith('P')) {
            connectAtoms(selectedNode.key, key);
        } else if (selectedNode.key.startsWith('P') && key.startsWith('R')) {
            connectAtoms(key, selectedNode.key);
        }

        selectedNode.el.setAttribute('fill', 'transparent');
        selectedNode = null;
        drawWires();
    }

    const buildSideCards = (side, data, startX) => {
        data.svgs.forEach((svgStr, fragIdx) => {
            const card = createFragCard(side, svgStr, fragIdx, data.coords[fragIdx], startX, getScale, handleAtomClick, drawWires);
            zoomContainer.appendChild(card);
        });
    };

    buildSideCards('R', ui.reactants, 50);
    buildSideCards('P', ui.products, 600);

    container.appendChild(zoomContainer);
    schemaContainer.appendChild(container);

    function drawWires() {
        overlay.innerHTML = '';
        const zoomRect = zoomContainer.getBoundingClientRect();
        const scale = getScale();

        for (let rKey in window.currentMapping) {
            const pKey = window.currentMapping[rKey];
            const rEl = document.getElementById(`atom_${rKey}`);
            const pEl = document.getElementById(`atom_${pKey}`);
            if (!rEl || !pEl) continue;

            const isApproved = !!window.approvedMappings[rKey];
            updateAtomHighlight(rEl, pEl, isApproved);
            if (isApproved && hideApproved) continue;

            const path = createWirePath(
                rKey, pKey, rEl, pEl, zoomRect, scale, isApproved,
                (rk) => {
                    window.approvedMappings[rk] = !window.approvedMappings[rk];
                    drawWires();
                },
                (rk) => {
                    delete window.currentMapping[rk];
                    delete window.approvedMappings[rk];
                    drawWires();
                }
            );
            overlay.appendChild(path);
        }
    }

    setTimeout(drawWires, 100);
    window.addEventListener('resize', drawWires);
}

export const MappingUI = {
    render: renderMappingUI,
    getData: async () => ({ mapping: window.currentMapping })
};
