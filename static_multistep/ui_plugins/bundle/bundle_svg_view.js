/**
 * High-contrast Paper Card & Zoom Viewport for SVG Mechanisms.
 */

export function renderSvgPaperCard(containerEl, svgEl, zoomCallbacks) {
    if (!containerEl) return null;
    containerEl.innerHTML = '';

    const wrapper = document.createElement('div');
    wrapper.id = 'svg-paper-wrapper';
    wrapper.className = 'svg-paper-wrapper';
    wrapper.style.width = '100%';
    wrapper.style.height = '100%';
    wrapper.style.display = 'flex';
    wrapper.style.flexDirection = 'column';
    wrapper.style.position = 'relative';
    wrapper.style.overflow = 'hidden';

    // Zoom toolbar overlay
    const toolbar = document.createElement('div');
    toolbar.style.position = 'absolute';
    toolbar.style.top = '8px';
    toolbar.style.right = '8px';
    toolbar.style.display = 'flex';
    toolbar.style.gap = '4px';
    toolbar.style.zIndex = '20';
    toolbar.innerHTML = `
        <button id="btn-svg-zoom-in" class="btn btn-xs btn-secondary" title="Zoom In">+</button>
        <button id="btn-svg-zoom-out" class="btn btn-xs btn-secondary" title="Zoom Out">−</button>
        <button id="btn-svg-zoom-fit" class="btn btn-xs btn-secondary" title="Fit Canvas">⟲</button>
    `;
    wrapper.appendChild(toolbar);

    // Scrollable viewport
    const viewport = document.createElement('div');
    viewport.id = 'svg-paper-viewport';
    viewport.style.flex = '1';
    viewport.style.overflow = 'auto';
    viewport.style.display = 'flex';
    viewport.style.alignItems = 'center';
    viewport.style.justifyContent = 'center';
    viewport.style.padding = '1rem';
    viewport.style.background = 'radial-gradient(circle at center, rgba(30,41,59,0.5), rgba(15,23,42,0.9))';

    // White paper card backing
    const paperCard = document.createElement('div');
    paperCard.id = 'svg-paper-card';
    paperCard.style.background = '#ffffff';
    paperCard.style.borderRadius = '8px';
    paperCard.style.padding = '16px 24px';
    paperCard.style.boxShadow = '0 10px 35px -5px rgba(0, 0, 0, 0.6), 0 0 0 1px rgba(255, 255, 255, 0.15)';
    paperCard.style.transformOrigin = 'center center';
    paperCard.style.transition = 'transform 0.15s ease-out';
    paperCard.style.minWidth = 'min(100%, 880px)';
    paperCard.style.maxWidth = '1100px';
    paperCard.style.height = 'fit-content';
    paperCard.style.margin = 'auto';

    const clonedSvg = svgEl.cloneNode(true);
    clonedSvg.setAttribute('id', 'active-mechanism-svg');

    // Ensure valid viewBox
    let w = parseFloat(clonedSvg.getAttribute('width'));
    let h = parseFloat(clonedSvg.getAttribute('height'));
    const vb = clonedSvg.getAttribute('viewBox');
    if (vb) {
        const parts = vb.trim().split(/[\s,]+/);
        if (parts.length === 4) {
            w = parseFloat(parts[2]) || w;
            h = parseFloat(parts[3]) || h;
        }
    } else if (w && h) {
        clonedSvg.setAttribute('viewBox', `0 0 ${w} ${h}`);
    }

    clonedSvg.removeAttribute('width');
    clonedSvg.removeAttribute('height');
    clonedSvg.setAttribute('style', 'display: block; width: 100%; height: auto;');

    paperCard.appendChild(clonedSvg);
    viewport.appendChild(paperCard);
    wrapper.appendChild(viewport);
    containerEl.appendChild(wrapper);

    toolbar.querySelector('#btn-svg-zoom-in').onclick = () => zoomCallbacks.zoom(1.25);
    toolbar.querySelector('#btn-svg-zoom-out').onclick = () => zoomCallbacks.zoom(0.8);
    toolbar.querySelector('#btn-svg-zoom-fit').onclick = () => zoomCallbacks.resetZoom();

    return { wrapper, paperCard, clonedSvg };
}
