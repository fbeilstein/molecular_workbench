/**
 * Orbital Media Pool for Step 6 Orbitals.
 * Categorizes and displays available molecular orbitals and localized IBOs
 * with distinctive styling for bonding (cyan) vs anti-bonding (orange/magenta).
 * Uses decoupled invariant track_ids for rock-solid cross-frame tracking.
 */
import { getOrbitalTypeBadge, isAntiBondingOrbital } from './orbitals_types.js';

export class OrbitalsPool {
    constructor(containerEl, onAudition, onAddToTimeline) {
        this.containerEl = containerEl;
        this.onAudition = onAudition;
        this.onAddToTimeline = onAddToTimeline;
        this.orbitals = [];
        this.visibleSet = new Set();
        this.currentCategory = 'all';
        this.searchQuery = '';
    }

    syncVisibilitySet(set) {
        this.visibleSet = new Set(set);
        this.render();
    }

    setCalculating(isCalc, progress = null) {
        this.isCalculating = isCalc;
        this.calcProgress = progress || (isCalc ? { percent: 0, completed_frames: 0, total_frames: 0, message: 'Initializing parallel PySCF SCF + IBO calculations...' } : null);
        this.render();
    }

    updateProgress(progress) {
        if (!progress) return;
        this.calcProgress = progress;
        if (!this.isCalculating) {
            this.isCalculating = true;
            this.render();
            return;
        }
        const bar = this.containerEl.querySelector('#pool-calc-bar');
        const pctEl = this.containerEl.querySelector('#pool-calc-pct');
        const stageEl = this.containerEl.querySelector('#pool-calc-stage');
        const msgEl = this.containerEl.querySelector('#pool-calc-msg');
        const pct = Math.min(100, Math.max(0, progress.percent || 0));

        if (bar) bar.style.width = `${pct}%`;
        if (pctEl) pctEl.textContent = `${pct}%`;
        if (stageEl) {
            if (progress.completed_frames && progress.total_frames) {
                stageEl.textContent = `Frame ${progress.completed_frames} / ${progress.total_frames}`;
            } else if (progress.stage === 'tracking') {
                stageEl.textContent = 'Propagating Tracks';
            } else if (progress.stage === 'saving') {
                stageEl.textContent = 'Finalizing';
            } else {
                stageEl.textContent = 'Computing';
            }
        }
        if (msgEl && progress.message) msgEl.textContent = progress.message;
    }

    setNotComputed(notComputed) {
        this.notComputed = notComputed;
        this.render();
    }

    setOrbitalsFromFrameData(frameData) { this.setFrameData(frameData); }

    setFrameData(frameData) {
        if (!frameData) return;
        this.frameData = frameData;
        const pool = [];
        const stIdx = frameData.step_index;

        if (frameData.frontier_orbitals) {
            for (const [k, fmo] of Object.entries(frameData.frontier_orbitals)) {
                const isLumo = k.toUpperCase().includes('LUMO');
                const tid = fmo.track_id || (isLumo ? 'trk_lumo' : 'trk_homo');
                pool.push({
                    track_id: tid, orbital_key: tid, key: tid, local_key: k, audition_id: tid, step_index: null,
                    chem_key: fmo.chem_key || k, label: fmo.label || k, type: isLumo ? 'LUMO' : 'HOMO',
                    category: 'frontier', orbital_class: isLumo ? 'antibonding' : 'bonding',
                    desc: fmo.description || `Frontier MO (${fmo.energy_ev != null ? (fmo.energy_ev > 0 ? '+' : '') + fmo.energy_ev.toFixed(2) + ' eV' : ''})`,
                    atoms: fmo.atoms || [], color: isLumo ? '#fb923c' : '#38bdf8', palette: fmo.palette || null,
                    energy: fmo.energy_ev, energy_ev: fmo.energy_ev, pop: fmo.pop || []
                });
            }
        }
        if (frameData.all_ibos && Array.isArray(frameData.all_ibos)) {
            for (const ibo of frameData.all_ibos) {
                const isCore = ibo.type === 'core', isLp = ibo.type === 'lone_pair';
                const cat = isCore ? 'core' : (isLp ? 'lone_pairs' : 'bonding');
                const key = ibo.key || ibo.name, tid = ibo.track_id || key;
                const trkLbl = ibo.track_label || ibo.decoupled_label || ibo.label || key;
                pool.push({
                    track_id: tid, orbital_key: tid, key: tid, local_key: key, audition_id: tid, step_index: stIdx,
                    chem_key: ibo.chem_key || key, label: ibo.label || key, frame_label: ibo.frame_label || ibo.label || key,
                    track_label: trkLbl, type: ibo.type || 'ibo', category: cat, orbital_class: isCore ? 'core' : 'bonding',
                    desc: ibo.description || (isLp ? `Lone pair on ${ibo.atoms ? ibo.atoms[0] : 'atom'}` : `Bond ${ibo.atoms ? ibo.atoms.join('-') : ''}`),
                    atoms: ibo.atoms || [], track_atoms: ibo.track_atoms || ibo.atoms || [],
                    color: isCore ? '#94a3b8' : '#38bdf8', palette: ibo.palette || null,
                    energy: ibo.energy_ev, energy_ev: ibo.energy_ev, pop: ibo.pop || []
                });
            }
        }
        if (frameData.canonical_pi && Array.isArray(frameData.canonical_pi)) {
            for (const cpi of frameData.canonical_pi) {
                const isStar = (cpi.type && cpi.type.includes('star')) || (cpi.label && cpi.label.includes('π*'));
                const key = cpi.key || cpi.name, tid = cpi.track_id || key;
                const trkLbl = cpi.track_label || cpi.decoupled_label || cpi.label || key;
                pool.push({
                    track_id: tid, orbital_key: tid, key: tid, local_key: key, audition_id: tid, step_index: stIdx,
                    chem_key: cpi.chem_key || key, label: cpi.label || key, frame_label: cpi.frame_label || cpi.label || key,
                    track_label: trkLbl, type: cpi.type || (isStar ? 'delocalized_pistar' : 'delocalized_pi'),
                    category: isStar ? 'antibonding' : 'delocalized', orbital_class: isStar ? 'antibonding' : 'bonding',
                    desc: cpi.description || (isStar ? `Delocalized π* on ${cpi.atoms ? cpi.atoms.join('-') : 'fragment'}` : 'Delocalized conjugated π system'),
                    atoms: cpi.atoms || [], track_atoms: cpi.track_atoms || cpi.atoms || [],
                    color: isStar ? '#fb923c' : '#38bdf8', palette: cpi.palette || null,
                    energy: cpi.energy_ev, energy_ev: cpi.energy_ev, pop: cpi.pop || []
                });
            }
        }
        if (frameData.antibonds && Array.isArray(frameData.antibonds)) {
            for (const ab of frameData.antibonds) {
                const isPiStar = (ab.type && ab.type.includes('pi')) || (ab.label && ab.label.includes('π*'));
                const key = ab.key || ab.name, tid = ab.track_id || key;
                const trkLbl = ab.track_label || ab.decoupled_label || ab.label || key;
                pool.push({
                    track_id: tid, orbital_key: tid, key: tid, local_key: key, audition_id: tid, step_index: stIdx,
                    chem_key: ab.chem_key || key, label: ab.label || key, frame_label: ab.frame_label || ab.label || key,
                    track_label: trkLbl, type: ab.type || 'antibonding', category: 'antibonding', orbital_class: 'antibonding',
                    desc: ab.description || `Antibond ${ab.label}`, atoms: ab.atoms || [], track_atoms: ab.track_atoms || ab.atoms || [],
                    color: isPiStar ? '#db2777' : '#fb923c', palette: ab.palette || null,
                    energy: ab.energy_ev, energy_ev: ab.energy_ev, pop: ab.pop || []
                });
            }
        }

        const espItem = {
            track_id: 'trk_esp', orbital_key: 'track_esp_charge', key: 'trk_esp', local_key: 'track_esp_charge',
            audition_id: 'trk_esp', step_index: null, chem_key: 'esp_charge',
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

    isAuditioning(item) {
        if (!item) return false;
        const tid = item.track_id || item.audition_id || item.key;
        return Boolean(this.visibleSet.has(tid) || (item.track_id && this.visibleSet.has(item.track_id)));
    }

    toggleAudition(tid) {
        if (this.visibleSet.has(tid)) this.visibleSet.delete(tid);
        else this.visibleSet.add(tid);
        this.render();
    }

    getAuditionedOrbitals() { return this.orbitals.filter(o => this.isAuditioning(o)); }

    render() {
        if (!this.containerEl) return;
        const filtered = this.orbitals.filter(o => {
            if (this.currentCategory !== 'all') {
                if (this.currentCategory === 'esp' && (o.track_id === 'trk_esp' || o.key === 'track_esp_charge')) { /* pass */ }
                else if (o.category !== this.currentCategory) return false;
            }
            if (this.searchQuery) {
                const q = this.searchQuery;
                const eVal = (typeof o.energy_ev === 'number') ? o.energy_ev : (typeof o.energy === 'number' ? o.energy : null);
                const eStr = eVal !== null ? `${eVal > 0 ? '+' : ''}${eVal.toFixed(2)}` : '';
                return (o.track_id || '').toLowerCase().includes(q) || (o.key || '').toLowerCase().includes(q) ||
                    (o.label || '').toLowerCase().includes(q) || (o.track_label || '').toLowerCase().includes(q) ||
                    (o.desc || '').toLowerCase().includes(q) || (o.type || '').toLowerCase().includes(q) ||
                    eStr.includes(q) ||
                    (o.atoms && o.atoms.some(a => String(a).toLowerCase().includes(q))) ||
                    (o.track_atoms && o.track_atoms.some(a => String(a).toLowerCase().includes(q)));
            }
            return true;
        });

        if (this.isCalculating) {
            const p = this.calcProgress || { percent: 0, completed_frames: 0, total_frames: 0, message: 'Starting parallel PySCF calculations...' };
            const pct = Math.min(100, Math.max(0, p.percent || 0));
            const frameTxt = (p.completed_frames && p.total_frames) ? `Frame ${p.completed_frames} / ${p.total_frames}` : (p.stage === 'tracking' ? 'Propagating Tracks' : 'Initializing');
            this.containerEl.innerHTML = `<div class="empty-state" style="padding: 2rem 1rem; text-align: center; color: #0284c7; font-size: 0.85rem;">
                <div style="font-size: 1.6rem; margin-bottom: 0.5rem; animation: pulse 1.5s infinite;">⚡</div>
                <strong style="color: #0369a1; font-size: 0.9rem;">Calculating Solvated Orbitals (PySCF ddCOSMO)...</strong>
                <div style="margin: 0.9rem 1.2rem 0.4rem 1.2rem; background: #e2e8f0; border-radius: 9999px; height: 10px; overflow: hidden; position: relative;">
                    <div id="pool-calc-bar" style="background: linear-gradient(90deg, #0284c7, #38bdf8); height: 100%; width: ${pct}%; transition: width 0.3s ease;"></div>
                </div>
                <div style="display: flex; justify-content: space-between; font-size: 0.74rem; font-weight: 600; color: #0369a1; padding: 0 1.3rem;">
                    <span id="pool-calc-stage">${frameTxt}</span>
                    <span id="pool-calc-pct">${pct}%</span>
                </div>
                <div id="pool-calc-msg" style="font-size: 0.72rem; color: #64748b; margin-top: 0.5rem; line-height: 1.4; padding: 0 0.5rem;">${p.message || 'Computing localized IBOs and tracking across frames. Please wait...'}</div>
            </div>`;
            return;
        }

        if (this.notComputed || this.orbitals.length === 0) {
            this.containerEl.innerHTML = `<div class="empty-state" style="padding: 1.8rem 1rem; text-align: center; color: #475569; font-size: 0.82rem; background: #f8fafc; border: 1px dashed #cbd5e1; border-radius: 8px; margin: 0.4rem;">
                <div style="font-size: 1.6rem; margin-bottom: 0.4rem;">⚛️</div>
                <strong style="color: #0f172a; font-size: 0.88rem;">Orbital Tracking Not Computed Yet</strong>
                <div style="font-size: 0.72rem; color: #64748b; margin: 0.4rem 0 0.8rem 0; line-height: 1.4;">
                    Run PySCF ddCOSMO implicit solvation wavefunctions and Knizia IBO tracking across all frames.
                </div>
                <button id="pool-trigger-calc" class="btn btn-sm btn-primary" style="padding: 0.35rem 0.9rem; font-weight: 600;">⚡ Calculate Orbitals</button>
            </div>`;
            this.containerEl.querySelector('#pool-trigger-calc')?.addEventListener('click', () => {
                const btn = document.querySelector('#btn-calc-orbitals');
                if (btn) btn.click();
            });
            return;
        }

        if (filtered.length === 0) {
            this.containerEl.innerHTML = `<div class="empty-state" style="padding: 1.5rem; text-align: center; color: var(--text-muted); font-size: 0.75rem;">
                No orbitals match current filter.
            </div>`;
            return;
        }

        let html = '<div style="display: flex; flex-direction: column; gap: 0.35rem;">';
        const atomLabels = this.frameData?.atom_labels || [];
        const symbols = this.frameData?.symbols || [];

        for (const item of filtered) {
            const isAuditioning = this.isAuditioning(item);
            const badge = getOrbitalTypeBadge(item), isAnti = isAntiBondingOrbital(item);
            const atomStr = item.atoms && item.atoms.length > 0 ? `[${item.atoms.join(', ')}]` : '';
            const isEsp = item.track_id === 'trk_esp' || item.key === 'track_esp_charge';
            const tid = item.track_id || item.key;
            const transPill = (item.track_label && item.track_label.includes('→'))
                ? `<span style="font-size: 0.65rem; padding: 0.05rem 0.35rem; border-radius: 3px; background: #e0f2fe; color: #0369a1; font-weight: 500;" title="Reaction trajectory">${item.track_label}</span>`
                : '';

            const energyVal = (typeof item.energy_ev === 'number') ? item.energy_ev : (typeof item.energy === 'number' ? item.energy : null);
            const energyBadge = (energyVal !== null && !isNaN(energyVal))
                ? `<span class="badge-energy" style="font-size: 0.68rem; font-family: monospace; font-weight: 700; padding: 0.05rem 0.35rem; border-radius: 3px; background: ${isAnti ? 'rgba(251,146,60,0.18)' : 'rgba(14,165,233,0.15)'}; color: ${isAnti ? '#c2410c' : '#0369a1'}; border: 1px solid ${isAnti ? 'rgba(251,146,60,0.4)' : 'rgba(14,165,233,0.35)'};" title="Orbital Energy at current frame">${energyVal > 0 ? '+' : ''}${energyVal.toFixed(2)} eV</span>`
                : '';

            let cardBorder = isAuditioning ? (isAnti ? '#ea580c' : '#0284c7') : (isAnti ? '#fed7aa' : (isEsp ? '#e9d5ff' : '#e2e8f0'));
            let cardBg = isAuditioning ? (isAnti ? '#fff7ed' : '#f0f9ff') : (isAnti ? '#fffbeb' : '#ffffff');

            // Compute top atomic populations if available (e.g. for delocalized & antibonding orbitals)
            let popDetail = '';
            if (item.pop && Array.isArray(item.pop) && item.pop.length > 0) {
                const scored = [];
                const n = Math.min(item.pop.length, symbols.length || item.pop.length);
                for (let i = 0; i < n; i++) {
                    const p = item.pop[i] || 0;
                    if (p > 0.08) {
                        const lbl = atomLabels[i] || `${symbols[i] || 'Atom'}${i}`;
                        scored.push({ lbl, pct: (p * 50).toFixed(0), val: p });
                    }
                }
                if (scored.length > 0) {
                    scored.sort((a, b) => b.val - a.val);
                    popDetail = scored.slice(0, 4).map(s => `${s.lbl}: ${s.pct}%`).join(' · ');
                }
            }

            html += `
                <div class="orbital-pool-card" data-tid="${tid}" style="background: ${cardBg}; border: 1px solid ${cardBorder}; border-radius: 6px; padding: 0.45rem 0.6rem; display: flex; justify-content: space-between; align-items: center; gap: 0.4rem; box-shadow: 0 1px 2px rgba(0,0,0,0.03);">
                    <div style="flex: 1; min-width: 0;">
                        <div style="display: flex; align-items: center; gap: 0.35rem; flex-wrap: wrap;">
                            <span style="font-weight: 600; font-size: 0.82rem; color: #0f172a;">${item.label}</span>
                            <span class="badge ${badge.cls}" style="font-size: 0.65rem; padding: 0.05rem 0.35rem; border-radius: 3px; ${badge.style}">${badge.label}</span>
                            ${energyBadge}
                            ${transPill}
                            <span style="font-size: 0.7rem; color: #64748b;">${tid}</span>
                        </div>
                        <div style="font-size: 0.7rem; color: #475569; margin-top: 0.15rem; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="${popDetail ? `${item.desc} (${popDetail})` : item.desc}">
                            ${popDetail ? `<strong style="color: #0f172a;">${popDetail}</strong> · ` : ''}${item.desc} ${atomStr}
                        </div>
                    </div>
                    <div style="display: flex; gap: 0.3rem; align-items: center;">
                        <button class="btn btn-xs ${isAuditioning ? 'btn-primary' : 'btn-secondary'} btn-audition" data-tid="${tid}" title="Toggle 3D visibility" style="font-size: 0.68rem; padding: 0.15rem 0.45rem;">
                            ${isAuditioning ? '👁 ON' : '👁 OFF'}
                        </button>
                        <button class="btn btn-xs btn-outline btn-add-track" data-tid="${tid}" title="Add track to timeline" style="font-size: 0.68rem; padding: 0.15rem 0.45rem;">
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
                const tid = btn.dataset.tid;
                if (this.onAudition) this.onAudition(tid);
            });
        });

        this.containerEl.querySelectorAll('.btn-add-track').forEach(btn => {
            btn.addEventListener('click', e => {
                e.stopPropagation();
                const tid = btn.dataset.tid;
                const orb = this.orbitals.find(o => (o.track_id && o.track_id === tid) || o.audition_id === tid || o.key === tid);
                if (this.onAddToTimeline && orb) this.onAddToTimeline(orb);
            });
        });
    }
}
