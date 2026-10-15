/**
 * Dynamic Media Pool: discovers and organizes all molecular orbitals dynamically from QM calculations.
 * Supports toggleable audition (multi-orbital ON/OFF in 3D view). Zero hardcoded molecules.
 */
import { getOrbitalTypeBadge } from './orbitals_types.js';

export class OrbitalsPool {
    constructor(containerEl, onAudition, onAddToTimeline) {
        this.containerEl = containerEl;
        this.onAudition = onAudition;
        this.onAddToTimeline = onAddToTimeline;
        this.orbitals = [];
        this.currentCategory = 'all';
        this.searchQuery = '';
        this.auditionKeys = new Set();
    }

    setOrbitalsFromFrameData(frameData) {
        const pool = [];
        if (frameData.frontier_orbitals) {
            for (const [k, meta] of Object.entries(frameData.frontier_orbitals)) {
                pool.push({
                    key: k, chem_key: meta.chem_key || k, label: meta.label || k,
                    type: 'frontier', category: 'frontier', desc: meta.description || `Frontier MO ${k}`,
                    atoms: meta.atoms || [], energy: meta.energy_ev,
                    color: k.includes('HOMO') ? '#06b6d4' : '#f59e0b'
                });
            }
        }
        if (frameData.all_ibos && Array.isArray(frameData.all_ibos)) {
            for (const ibo of frameData.all_ibos) {
                const typ = ibo.type || 'sigma';
                let cat = 'bonding';
                if (typ.includes('lone') || typ.includes('lp')) cat = 'lone_pairs';
                else if (typ.includes('deloc')) cat = 'delocalized';
                else if (typ === 'core') cat = 'core';

                pool.push({
                    key: ibo.key || ibo.name, chem_key: ibo.chem_key || ibo.key || ibo.name,
                    label: ibo.label || ibo.key || ibo.name, type: typ, category: cat,
                    desc: ibo.description || `${typ.toUpperCase()} on ${ibo.atoms ? ibo.atoms.join('-') : 'atoms'}`,
                    atoms: ibo.atoms || [], energy: ibo.energy_ev,
                    color: typ === 'pi' ? '#ec4899' : (cat === 'lone_pairs' ? '#a855f7' : (cat === 'core' ? '#64748b' : '#10b981'))
                });
            }
        }
        if (frameData.canonical_pi && Array.isArray(frameData.canonical_pi)) {
            for (const cpi of frameData.canonical_pi) {
                pool.push({
                    key: cpi.key || cpi.name, chem_key: cpi.chem_key || cpi.key,
                    label: cpi.label || cpi.key, type: 'delocalized_pi', category: 'delocalized',
                    desc: cpi.description || 'Delocalized conjugated π system', atoms: cpi.atoms || [], color: '#3b82f6'
                });
            }
        }
        if (frameData.antibonds && Array.isArray(frameData.antibonds)) {
            for (const ab of frameData.antibonds) {
                pool.push({
                    key: ab.key || ab.name, chem_key: ab.chem_key || ab.key,
                    label: ab.label || ab.key, type: ab.type || 'antibonding', category: 'antibonding',
                    desc: ab.description || `Antibond ${ab.label}`, atoms: ab.atoms || [], color: '#f43f5e'
                });
            }
        }
        const order = { frontier: 0, bonding: 1, lone_pairs: 2, delocalized: 3, antibonding: 4, core: 5 };
        pool.sort((a, b) => (order[a.category] ?? 9) - (order[b.category] ?? 9));
        this.orbitals = pool;
        this.render();
    }

    setCategory(cat) {
        this.currentCategory = cat;
        this.render();
    }

    setSearch(query) {
        this.searchQuery = (query || '').toLowerCase().trim();
        this.render();
    }

    toggleAudition(key) {
        if (this.auditionKeys.has(key)) {
            this.auditionKeys.delete(key);
        } else {
            this.auditionKeys.add(key);
        }
        this.render();
    }

    getAuditionedOrbitals() {
        return this.orbitals.filter(o => this.auditionKeys.has(o.key));
    }

    render() {
        if (!this.containerEl) return;
        const filtered = this.orbitals.filter(o => {
            if (this.currentCategory !== 'all' && o.category !== this.currentCategory) return false;
            if (this.searchQuery) {
                const q = this.searchQuery;
                const matchKey = o.key.toLowerCase().includes(q);
                const matchLbl = o.label.toLowerCase().includes(q);
                const matchAtoms = o.atoms && o.atoms.some(a => String(a).toLowerCase().includes(q));
                return matchKey || matchLbl || matchAtoms;
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
            const isAuditioning = this.auditionKeys.has(item.key);
            const badge = getOrbitalTypeBadge(item.type);
            const atomStr = item.atoms && item.atoms.length > 0 ? `[${item.atoms.join(', ')}]` : '';

            html += `
                <div class="orbital-pool-card" data-key="${item.key}" style="
                    background: ${isAuditioning ? 'rgba(6, 182, 212, 0.15)' : 'rgba(255, 255, 255, 0.03)'};
                    border: 1px solid ${isAuditioning ? '#06b6d4' : 'rgba(255, 255, 255, 0.08)'};
                    border-radius: 5px; padding: 0.4rem 0.6rem; display: flex; justify-content: space-between; align-items: center; gap: 0.4rem;
                ">
                    <div style="flex: 1; min-width: 0;">
                        <div style="display: flex; align-items: center; gap: 0.35rem; flex-wrap: wrap;">
                            <span style="font-weight: 600; font-size: 0.82rem; color: #f8fafc;">${item.label}</span>
                            <span class="badge ${badge.cls}" style="font-size: 0.65rem; padding: 0.05rem 0.3rem; border-radius: 3px; background: rgba(255,255,255,0.1);">${badge.label}</span>
                            <span style="font-size: 0.7rem; color: #94a3b8;">${item.key}</span>
                        </div>
                        <div style="font-size: 0.7rem; color: #64748b; margin-top: 0.15rem; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
                            ${item.desc} ${atomStr}
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.3rem; align-items: center;">
                        <button class="btn btn-xs ${isAuditioning ? 'btn-primary' : 'btn-secondary'} btn-audition" data-key="${item.key}" title="Toggle 3D visibility" style="font-size: 0.68rem; padding: 0.15rem 0.45rem;">
                            ${isAuditioning ? '👁 ON' : '👁 OFF'}
                        </button>
                        <button class="btn btn-xs btn-outline btn-add-track" data-key="${item.key}" title="Add track across all frames" style="font-size: 0.68rem; padding: 0.15rem 0.45rem;">
                            + Track
                        </button>
                    </div>
                </div>
            `;
        }
        html += '</div>';
        this.containerEl.innerHTML = html;

        this.containerEl.querySelectorAll('.btn-audition').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.stopPropagation();
                const key = btn.dataset.key;
                this.toggleAudition(key);
                const orb = this.orbitals.find(o => o.key === key);
                if (this.onAudition) this.onAudition(orb, this.auditionKeys.has(key));
            });
        });

        this.containerEl.querySelectorAll('.btn-add-track').forEach(btn => {
            btn.addEventListener('click', (e) => {
                e.stopPropagation();
                const key = btn.dataset.key;
                const orb = this.orbitals.find(o => o.key === key);
                if (this.onAddToTimeline) this.onAddToTimeline(orb);
            });
        });
    }
}
