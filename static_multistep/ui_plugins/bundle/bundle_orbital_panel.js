/**
 * Orbital Tracks Sidebar & Toggles Manager for Step 5 Bundle Viewer.
 */

export function renderOrbitalTracksSidebar(container, tracks, activeOrbitals, callbacks) {
    if (!container) return;
    container.innerHTML = '';

    if (!tracks || tracks.length === 0) {
        container.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.5rem;">No tracked orbital continuous series available.</div>`;
        return;
    }

    tracks.forEach(trk => {
        const row = document.createElement('label');
        row.style.cssText = 'display: flex; align-items: center; justify-content: space-between; gap: 0.5rem; padding: 0.35rem 0.6rem; border-radius: 4px; background: rgba(255,255,255,0.03); border: 1px solid rgba(255,255,255,0.06); cursor: pointer; font-size: 0.75rem;';

        const left = document.createElement('div');
        left.style.cssText = 'display: flex; align-items: center; gap: 0.45rem; overflow: hidden;';

        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.checked = activeOrbitals.has(trk.orbital_key);
        cb.style.accentColor = trk.color || '#38bdf8';

        cb.onchange = (e) => {
            if (e.target.checked) activeOrbitals.add(trk.orbital_key);
            else activeOrbitals.delete(trk.orbital_key);
            if (callbacks.onChanged) callbacks.onChanged();
        };

        const dot = document.createElement('span');
        dot.style.cssText = `display: inline-block; width: 9px; height: 9px; border-radius: 50%; background: ${trk.color || '#38bdf8'}; flex-shrink: 0;`;

        const title = document.createElement('span');
        title.style.cssText = 'overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-weight: 500;';
        title.textContent = trk.label || trk.orbital_key;
        title.title = `${trk.label} (${trk.category})`;

        left.appendChild(cb);
        left.appendChild(dot);
        left.appendChild(title);

        const badge = document.createElement('span');
        badge.className = 'chip-mini';
        badge.style.fontSize = '0.62rem';
        badge.textContent = trk.category || 'val';

        row.appendChild(left);
        row.appendChild(badge);
        container.appendChild(row);
    });
}
