/**
 * Orbital Tracks Sidebar & Toggles Manager for Step 5 Bundle Viewer.
 */
import { isAntiBondingOrbital, getOrbitalTypeBadge } from '../orbitals/orbitals_types.js';

export function renderOrbitalTracksSidebar(container, tracks, activeOrbitals, callbacks) {
    if (!container) return;
    container.innerHTML = '';

    if (!tracks || tracks.length === 0) {
        container.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.5rem;">No tracked orbital continuous series available.</div>`;
        return;
    }

    // Filter out ESP from the orbital tracks list to avoid duplicate ESP checkboxes (ESP is cleanly controlled in the toolbar)
    const moTracks = tracks.filter(trk => {
        const tid = (trk.orbital_key === 'track_esp_charge' || trk.track_id === 'trk_esp') ? 'track_esp_charge' : (trk.track_id || trk.orbital_key || trk.id);
        return tid !== 'track_esp_charge' && tid !== 'trk_esp';
    });

    if (moTracks.length === 0) {
        container.innerHTML = `<div style="font-size: 0.72rem; color: #64748b; padding: 0.5rem;">No tracked orbital continuous series available.</div>`;
        return;
    }

    const moveTrack = (fromIdx, toIdx) => {
        if (fromIdx < 0 || fromIdx >= moTracks.length || toIdx < 0 || toIdx >= moTracks.length) return;
        const item = moTracks[fromIdx];
        moTracks.splice(fromIdx, 1);
        moTracks.splice(toIdx, 0, item);

        // Keep non-MO tracks (like trk_esp) preserved at their original position
        const nonMoTracks = tracks.filter(t => {
            const tid = (t.orbital_key === 'track_esp_charge' || t.track_id === 'trk_esp') ? 'track_esp_charge' : (t.track_id || t.orbital_key || t.id);
            return tid === 'track_esp_charge' || tid === 'trk_esp';
        });
        const newTracks = [...nonMoTracks, ...moTracks];
        if (callbacks && callbacks.onReorder) callbacks.onReorder(newTracks);
    };

    moTracks.forEach((trk, idx) => {
        const display = formatOrbitalDisplay(trk);
        const row = document.createElement('div');
        row.className = 'orbital-track-row';
        row.draggable = true;
        row.dataset.idx = idx;
        row.style.cssText = 'display: flex; align-items: center; justify-content: space-between; gap: 0.4rem; padding: 0.35rem 0.5rem; border-radius: 4px; background: var(--bg-card, #f8fafc); border: 1px solid var(--border-main, #e2e8f0); color: var(--text-main, #0f172a); font-size: 0.75rem; user-select: none; transition: background 0.15s, border-color 0.15s; margin-bottom: 2px;';

        row.addEventListener('mouseenter', () => {
            if (!row.dataset.dragging) row.style.background = 'var(--bg-card-hover, #f1f5f9)';
        });
        row.addEventListener('mouseleave', () => {
            if (!row.dataset.dragging) row.style.background = 'var(--bg-card, #f8fafc)';
        });

        // Drag and drop support
        row.addEventListener('dragstart', (e) => {
            row.dataset.dragging = 'true';
            e.dataTransfer.setData('text/plain', String(idx));
            e.dataTransfer.effectAllowed = 'move';
            row.style.opacity = '0.4';
        });
        row.addEventListener('dragend', () => {
            delete row.dataset.dragging;
            row.style.opacity = '1';
            container.querySelectorAll('.orbital-track-row').forEach(r => {
                r.style.borderTop = '1px solid var(--border-main, #e2e8f0)';
                r.style.background = 'var(--bg-card, #f8fafc)';
            });
        });
        row.addEventListener('dragover', (e) => {
            e.preventDefault();
            e.dataTransfer.dropEffect = 'move';
            row.style.borderTop = '2px solid #0284c7';
            row.style.background = 'rgba(2,132,199,0.08)';
        });
        row.addEventListener('dragleave', () => {
            row.style.borderTop = '1px solid var(--border-main, #e2e8f0)';
            row.style.background = 'var(--bg-card, #f8fafc)';
        });
        row.addEventListener('drop', (e) => {
            e.preventDefault();
            row.style.borderTop = '1px solid var(--border-main, #e2e8f0)';
            row.style.background = 'var(--bg-card, #f8fafc)';
            const fromIdx = parseInt(e.dataTransfer.getData('text/plain'), 10);
            if (!isNaN(fromIdx) && fromIdx !== idx) {
                moveTrack(fromIdx, idx);
            }
        });

        const left = document.createElement('div');
        left.style.cssText = 'display: flex; align-items: center; gap: 0.4rem; overflow: hidden; min-width: 0; flex: 1;';

        // Move Up / Move Down buttons
        const btnBox = document.createElement('div');
        btnBox.style.cssText = 'display: flex; flex-direction: column; gap: 1px; flex-shrink: 0;';

        const btnUp = document.createElement('button');
        btnUp.type = 'button';
        btnUp.innerHTML = '▲';
        btnUp.title = 'Move orbital up';
        btnUp.disabled = (idx === 0);
        btnUp.style.cssText = `background: none; border: none; padding: 0; font-size: 0.55rem; line-height: 1; color: ${idx === 0 ? '#cbd5e1' : '#64748b'}; cursor: ${idx === 0 ? 'default' : 'pointer'};`;
        if (idx > 0) {
            btnUp.onmouseenter = () => { btnUp.style.color = '#0284c7'; };
            btnUp.onmouseleave = () => { btnUp.style.color = '#64748b'; };
            btnUp.onclick = (e) => { e.stopPropagation(); moveTrack(idx, idx - 1); };
        }

        const btnDown = document.createElement('button');
        btnDown.type = 'button';
        btnDown.innerHTML = '▼';
        btnDown.title = 'Move orbital down';
        btnDown.disabled = (idx === moTracks.length - 1);
        btnDown.style.cssText = `background: none; border: none; padding: 0; font-size: 0.55rem; line-height: 1; color: ${idx === moTracks.length - 1 ? '#cbd5e1' : '#64748b'}; cursor: ${idx === moTracks.length - 1 ? 'default' : 'pointer'};`;
        if (idx < moTracks.length - 1) {
            btnDown.onmouseenter = () => { btnDown.style.color = '#0284c7'; };
            btnDown.onmouseleave = () => { btnDown.style.color = '#64748b'; };
            btnDown.onclick = (e) => { e.stopPropagation(); moveTrack(idx, idx + 1); };
        }

        btnBox.appendChild(btnUp);
        btnBox.appendChild(btnDown);
        left.appendChild(btnBox);

        const tid = (trk.orbital_key === 'track_esp_charge' || trk.track_id === 'trk_esp') ? 'track_esp_charge' : (trk.track_id || trk.orbital_key || trk.id);
        const isAnti = isAntiBondingOrbital(trk);
        const baseColor = isAnti ? '#fb923c' : '#38bdf8';
        const accentCol = isAnti ? '#ea580c' : '#0284c7';

        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.checked = activeOrbitals.has(tid);
        cb.style.cssText = `accent-color: ${trk.color || accentCol}; cursor: pointer; flex-shrink: 0; margin: 0; width: 14px; height: 14px;`;

        cb.onchange = (e) => {
            if (e.target.checked) {
                if (tid) activeOrbitals.add(tid);
            } else {
                if (tid) activeOrbitals.delete(tid);
            }
            if (callbacks && callbacks.onChanged) callbacks.onChanged();
        };

        const dot = document.createElement('span');
        dot.style.cssText = `display: inline-block; width: 8px; height: 8px; border-radius: 50%; background: ${trk.color || baseColor}; flex-shrink: 0; cursor: pointer;`;
        dot.onclick = () => { cb.checked = !cb.checked; cb.dispatchEvent(new Event('change')); };

        const title = document.createElement('span');
        title.style.cssText = 'overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-weight: 500; font-size: 0.74rem; cursor: pointer; color: var(--text-main, #0f172a);';
        title.textContent = display.name;
        title.title = `${display.name} (drag or use ▲/▼ to reorder)`;
        title.onclick = () => { cb.checked = !cb.checked; cb.dispatchEvent(new Event('change')); };

        left.appendChild(cb);
        left.appendChild(dot);
        left.appendChild(title);

        const badgeInfo = getOrbitalTypeBadge(trk);
        const badge = document.createElement('span');
        badge.className = `chip-mini ${badgeInfo.cls}`;
        badge.style.cssText = `font-size: 0.62rem; padding: 1px 6px; border-radius: 3px; font-weight: 700; white-space: nowrap; flex-shrink: 0; ${badgeInfo.style}`;
        badge.textContent = badgeInfo.label;
        badge.title = display.tooltip;

        row.appendChild(left);
        row.appendChild(badge);
        container.appendChild(row);
    });
}

export function formatOrbitalDisplay(trk) {
    if (!trk) return { name: '', badge: 'val', tooltip: '' };

    if (trk.orbital_key === 'track_esp_charge' || trk.id === 'trk_esp_charge') {
        return {
            name: 'ESP Surface',
            badge: 'ESP',
            tooltip: 'Electrostatic Potential Surface & Partial Charges'
        };
    }

    let name = trk.label || trk.title || '';
    let badge = trk.role_transition || trk.category || '';

    if (trk.role_transition && trk.start_role && trk.end_role && trk.start_role !== trk.end_role && !name.includes('→') && !name.includes('->')) {
        name = `${trk.start_role} → ${trk.end_role}`;
    }

    if (!name || name.startsWith('IBO_') || name.startsWith('STAR_') || name.startsWith('trk_')) {
        name = cleanTechnicalKey(trk.chem_key || trk.orbital_key || '', trk.type);
    }

    name = name.replace(/Non-Bonding\s*/gi, '');
    name = name.replace(/\[Antibond\]/i, '').trim();

    if (badge.includes('→') && badge.includes('(')) {
        badge = badge.split('→').map(p => p.trim().split('(')[0]).join(' → ');
    } else if (!badge || badge === 'val' || badge === 'bonding' || badge === 'antibonding') {
        if (name.includes('→')) {
            const parts = name.split('→').map(p => p.trim().split('(')[0]);
            badge = parts.join(' → ');
        } else if (name.includes('π*') || name.toLowerCase().includes('pistar')) {
            badge = 'π*';
        } else if (name.includes('σ*') || name.toLowerCase().includes('sigmastar')) {
            badge = 'σ*';
        } else if (name.includes('π') || name.toLowerCase().includes('pi')) {
            badge = 'π';
        } else if (name.includes('σ') || name.toLowerCase().includes('sigma')) {
            badge = 'σ';
        } else if (name.includes('LP') || name.toLowerCase().includes('lp')) {
            badge = 'LP';
        } else {
            badge = 'val';
        }
    }

    const tooltip = `${name} [${badge}] (${trk.orbital_key || ''})`;
    return { name, badge, tooltip };
}

function cleanTechnicalKey(ck, typ) {
    if (!ck) return 'Orbital';
    if (ck.startsWith('lp_')) {
        const at = ck.replace('lp_', '');
        return `LP(${at})`;
    }
    if (ck.startsWith('sigma_')) {
        const at = ck.replace('sigma_', '');
        return `σ(${at})`;
    }
    if (ck.startsWith('pi_')) {
        const at = ck.replace('pi_', '');
        return `π(${at})`;
    }
    if (ck.startsWith('pistar_')) {
        const at = ck.replace('pistar_', '').replace(/_\d+$/, '');
        return `π*(${at})`;
    }
    if (ck.startsWith('sigmastar_')) {
        const at = ck.replace('sigmastar_', '').replace(/_\d+$/, '');
        return `σ*(${at})`;
    }
    return ck.replace(/_/g, ' ');
}
