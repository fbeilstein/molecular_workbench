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

    moTracks.forEach(trk => {
        const display = formatOrbitalDisplay(trk);
        const row = document.createElement('label');
        row.style.cssText = 'display: flex; align-items: center; justify-content: space-between; gap: 0.5rem; padding: 0.35rem 0.6rem; border-radius: 4px; background: rgba(255,255,255,0.03); border: 1px solid rgba(255,255,255,0.06); cursor: pointer; font-size: 0.75rem;';

        const left = document.createElement('div');
        left.style.cssText = 'display: flex; align-items: center; gap: 0.45rem; overflow: hidden; min-width: 0; flex: 1;';

        const tid = (trk.orbital_key === 'track_esp_charge' || trk.track_id === 'trk_esp') ? 'track_esp_charge' : (trk.track_id || trk.orbital_key || trk.id);
        const isAnti = isAntiBondingOrbital(trk);
        const baseColor = isAnti ? '#fb923c' : '#38bdf8';
        const accentCol = isAnti ? '#ea580c' : '#0284c7';

        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.checked = activeOrbitals.has(tid);
        cb.style.accentColor = trk.color || accentCol;

        cb.onchange = (e) => {
            if (e.target.checked) {
                if (tid) activeOrbitals.add(tid);
            } else {
                if (tid) activeOrbitals.delete(tid);
            }
            if (callbacks.onChanged) callbacks.onChanged();
        };

        const dot = document.createElement('span');
        dot.style.cssText = `display: inline-block; width: 9px; height: 9px; border-radius: 50%; background: ${trk.color || baseColor}; flex-shrink: 0;`;

        const title = document.createElement('span');
        title.style.cssText = 'overflow: hidden; text-overflow: ellipsis; white-space: nowrap; font-weight: 500; font-size: 0.76rem;';
        title.textContent = display.name;
        title.title = display.tooltip;

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
