/**
 * Orbital Media Pool for Step 6 Orbitals.
 * Categorizes and displays available molecular orbitals and localized IBOs
 * with distinctive styling for bonding (cyan) vs anti-bonding (orange/magenta).
 */
import { getOrbitalTypeBadge, isAntiBondingOrbital } from './orbitals_types.js';

export class OrbitalsPool {
    constructor(containerEl, onAudition, onAddToTimeline) {
        this.containerEl = containerEl;
        this.onAudition = onAudition;
        this.onAddToTimeline = onAddToTimeline;
        this.orbitals = [];
        this.auditionKeys = new Set();
        this.currentCategory = 'all';
        this.searchQuery = '';
    }

    setOrbitalsFromFrameData(frameData) {
        this.setFrameData(frameData);
    }

    setFrameData(frameData) {
        if (!frameData) return;
        const pool = [];
        const stIdx = frameData.step_index;

        if (frameData.frontier_orbitals) {
            for (const [k, fmo] of Object.entries(frameData.frontier_orbitals)) {
                const isLumo = k.toUpperCase().includes('LUMO');
                pool.push({
                    key: k, audition_id: k, step_index: null, chem_key: fmo.chem_key || k,
                    label: fmo.label || k, type: isLumo ? 'LUMO' : 'HOMO', category: 'frontier',
                    orbital_class: isLumo ? 'antibonding' : 'bonding',
                    desc: fmo.description || `Frontier MO (${fmo.energy_ev != null ? fmo.energy_ev.toFixed(2) + ' eV' : ''})`,
                    atoms: fmo.atoms || [], color: isLumo ? '#fb923c' : '#38bdf8', palette: fmo.palette || null, energy: fmo.energy_ev
                });
            }
        }
        if (frameData.all_ibos && Array.isArray(frameData.all_ibos)) {
            for (const ibo of frameData.all_ibos) {
                const isCore = ibo.type === 'core', isLp = ibo.type === 'lone_pair';
                const cat = isCore ? 'core' : (isLp ? 'lone_pairs' : 'bonding');
                const key = ibo.key || ibo.name;
                pool.push({
                    key, audition_id: `${stIdx}::${key}`, step_index: stIdx, chem_key: ibo.chem_key || key,
                    label: ibo.label || key, type: ibo.type || 'ibo', category: cat, orbital_class: 'bonding',
                    desc: ibo.description || (isLp ? `Lone pair on ${ibo.atoms ? ibo.atoms[0] : 'atom'}` : `Bond ${ibo.atoms ? ibo.atoms.join('-') : ''}`),
                    atoms: ibo.atoms || [], color: '#38bdf8', palette: ibo.palette || null, energy: ibo.energy_ev
                });
            }
        }
        if (frameData.canonical_pi && Array.isArray(frameData.canonical_pi)) {
            for (const cpi of frameData.canonical_pi) {
                const isStar = (cpi.type && cpi.type.includes('star')) || (cpi.label && cpi.label.includes('π*'));
                const key = cpi.key || cpi.name;
                pool.push({
                    key, audition_id: `${stIdx}::${key}`, step_index: stIdx, chem_key: cpi.chem_key || key,
                    label: cpi.label || key, type: cpi.type || (isStar ? 'delocalized_pistar' : 'delocalized_pi'),
                    category: isStar ? 'antibonding' : 'delocalized', orbital_class: isStar ? 'antibonding' : 'bonding',
                    desc: cpi.description || (isStar ? `Delocalized π* on ${cpi.atoms ? cpi.atoms.join('-') : 'fragment'}` : 'Delocalized conjugated π system'),
                    atoms: cpi.atoms || [], color: isStar ? '#fb923c' : '#38bdf8', palette: cpi.palette || null, energy: cpi.energy_ev
                });
            }
        }
        if (frameData.antibonds && Array.isArray(frameData.antibonds)) {
            for (const ab of frameData.antibonds) {
                const isPiStar = (ab.type && ab.type.includes('pi')) || (ab.label && ab.label.includes('π*'));
                const key = ab.key || ab.name;
                pool.push({
                    key, audition_id: `${stIdx}::${key}`, step_index: stIdx, chem_key: ab.chem_key || key,
                    label: ab.label || key, type: ab.type || 'antibonding', category: 'antibonding', orbital_class: 'antibonding',
                    desc: ab.description || `Antibond ${ab.label}`, atoms: ab.atoms || [],
                    color: isPiStar ? '#db2777' : '#fb923c', palette: ab.palette || null, energy: ab.energy_ev
                });
            }
        }

        const espItem = {
            key: 'track_esp_charge', audition_id: 'track_esp_charge', step_index: null, chem_key: 'esp_charge',
            label: '⚡ Electrostatic Surface & Charges (ESP)', type: 'esp', category: 'esp', orbital_class: 'esp',
            desc: 'Molecular electrostatic surface mapped with atomic partial charges', atoms: [], energy: null, color: '#c084fc'
        };

        const order = { esp: -1, frontier: 0, bonding: 1, lone_pairs: 2, delocalized: 3, antibonding: 4, core: 5 };
        pool.unshift(espItem);
        pool.sort((a, b) => (order[a.category] ?? 9) - (order[b.category] ?? 9));
        this.orbitals = pool;
        this.render();
    }

    setCategory(cat) { this.currentCategory = cat; this.render(); }
    setSearch(query) { this.searchQuery = (query || '').toLowerCase().trim(); this.render(); }

    toggleAudition(auditionId) {
        if (this.auditionKeys.has(auditionId)) this.auditionKeys.delete(auditionId);
        else this.auditionKeys.add(auditionId);
        this.render();
    }

    getAuditionedOrbitals() {
        return this.orbitals.filter(o => this.auditionKeys.has(o.audition_id));
    }

    render() {
        if (!this.containerEl) return;
        const filtered = this.orbitals.filter(o => {
            if (this.currentCategory !== 'all') {
                if (this.currentCategory === 'esp' && o.key === 'track_esp_charge') { /* pass */ }
                else if (o.category !== this.currentCategory) return false;
            }
            if (this.searchQuery) {
                const q = this.searchQuery;
                return o.key.toLowerCase().includes(q) || o.label.toLowerCase().includes(q) ||
                    (o.desc || '').toLowerCase().includes(q) || (o.chem_key || '').toLowerCase().includes(q) ||
                    (o.type || '').toLowerCase().includes(q) || (o.atoms && o.atoms.some(a => String(a).toLowerCase().includes(q)));
            }
            return true;
        });

        if (filtered.length === 0) {
            this.containerEl.innerHTML = `<div class="empty-state" style="padding: 1.5rem; text-align: center; color: var(--text-muted); font-size: 0.75rem;">
                ${this.orbitals.length === 0 ? 'Computing molecular orbitals for current frame...' : 'No orbitals match current filter.'}
            </div>`;
            return;
        }

        let html = '<div style="display: flex; flex-direction: column; gap: 0.35rem;">';
        for (const item of filtered) {
            const isAuditioning = this.auditionKeys.has(item.audition_id);
            const badge = getOrbitalTypeBadge(item), isAnti = isAntiBondingOrbital(item);
            const atomStr = item.atoms && item.atoms.length > 0 ? `[${item.atoms.join(', ')}]` : '';
            const isEsp = item.key === 'track_esp_charge';

            let cardBorder = isAuditioning ? (isAnti ? '#ea580c' : '#0284c7') : (isAnti ? '#fed7aa' : (isEsp ? '#e9d5ff' : '#e2e8f0'));
            let cardBg = isAuditioning ? (isAnti ? '#fff7ed' : '#f0f9ff') : (isAnti ? '#fffbeb' : '#ffffff');

            html += `
                <div class="orbital-pool-card" data-aid="${item.audition_id}" style="background: ${cardBg}; border: 1px solid ${cardBorder}; border-radius: 6px; padding: 0.45rem 0.6rem; display: flex; justify-content: space-between; align-items: center; gap: 0.4rem; box-shadow: 0 1px 2px rgba(0,0,0,0.03);">
                    <div style="flex: 1; min-width: 0;">
                        <div style="display: flex; align-items: center; gap: 0.35rem; flex-wrap: wrap;">
                            <span style="font-weight: 600; font-size: 0.82rem; color: #0f172a;">${item.label}</span>
                            <span class="badge ${badge.cls}" style="font-size: 0.65rem; padding: 0.05rem 0.35rem; border-radius: 3px; ${badge.style}">${badge.label}</span>
                            <span style="font-size: 0.7rem; color: #64748b;">${item.key}</span>
                        </div>
                        <div style="font-size: 0.7rem; color: #475569; margin-top: 0.15rem; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
                            ${item.desc} ${atomStr}
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.3rem; align-items: center;">
                        <button class="btn btn-xs ${isAuditioning ? 'btn-primary' : 'btn-secondary'} btn-audition" data-aid="${item.audition_id}" data-key="${item.key}" title="Toggle 3D visibility" style="font-size: 0.68rem; padding: 0.15rem 0.45rem;">
                            ${isAuditioning ? '👁 ON' : '👁 OFF'}
                        </button>
                        <button class="btn btn-xs btn-outline btn-add-track" data-aid="${item.audition_id}" data-key="${item.key}" title="Add track to timeline" style="font-size: 0.68rem; padding: 0.15rem 0.45rem;">
                            + Track
                        </button>
                    </div>
                </div>`;
        }
        html += '</div>';
        this.containerEl.innerHTML = html;

        this.containerEl.querySelectorAll('.btn-audition').forEach(btn => {
            btn.addEventListener('click', e => {
                e.stopPropagation();
                const aid = btn.dataset.aid || btn.dataset.key;
                this.toggleAudition(aid);
                const orb = this.orbitals.find(o => o.audition_id === aid || o.key === btn.dataset.key);
                if (this.onAudition) this.onAudition(orb, this.auditionKeys.has(aid));
            });
        });

        this.containerEl.querySelectorAll('.btn-add-track').forEach(btn => {
            btn.addEventListener('click', e => {
                e.stopPropagation();
                const aid = btn.dataset.aid || btn.dataset.key;
                const orb = this.orbitals.find(o => o.audition_id === aid || o.key === btn.dataset.key);
                if (this.onAddToTimeline && orb) this.onAddToTimeline(orb);
            });
        });
    }
}
