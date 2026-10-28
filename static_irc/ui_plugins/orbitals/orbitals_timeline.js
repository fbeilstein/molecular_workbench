/**
 * Reaction Orbital Selection & Real-Time Electron Tracking Panel.
 * Uses decoupled invariant track_ids to reliably follow electron pairs across cascade steps.
 * Real-time electron tracking displays electron population distribution across atoms per frame.
 * Visually distinguishes bonding (Cyan/Blue) vs anti-bonding (Orange/Magenta) tracks.
 */
import { COLOR_PALETTE, isAntiBondingOrbital, getOrbitalTypeBadge } from './orbitals_types.js';

export class OrbitalsTimeline {
    constructor(containerEl, options) {
        this.containerEl = containerEl;
        this.options = options || {};
        this.tracks = [];
        this.totalFrames = options.totalFrames || 0;
        this.currentFrame = options.currentFrame || 0;
        this.manualOverrides = new Map();
        this.lastFrameData = null;
        this.symbols = [];
    }

    setTimelineData(timeline, totalFrames) {
        this.tracks = (timeline && timeline.tracks) ? JSON.parse(JSON.stringify(timeline.tracks)) : [];
        if (totalFrames) this.totalFrames = totalFrames;
        this.render();
    }

    setCurrentFrame(frameIdx) {
        this.currentFrame = Math.max(0, Math.min(this.totalFrames - 1, frameIdx));
        if (this.lastFrameData) this.updatePopulations(this.currentFrame, this.lastFrameData, this.symbols);
    }

    addTrack(orbital) {
        if (!orbital) return;
        const tid = orbital.track_id || orbital.key || orbital.orbital_key;
        if (tid === 'trk_esp' || tid === 'track_esp_charge') { this.addEspTrack(); return; }
        const exists = this.tracks.find(t => (t.track_id && t.track_id === tid) || t.id === tid || t.orbital_key === tid);
        if (exists) {
            this.manualOverrides.set(exists.id, true);
            this.render();
            if (this.options.onVisibilityChange) this.options.onVisibilityChange();
            return;
        }
        const isAnti = isAntiBondingOrbital(orbital);
        const color = orbital.color || (isAnti ? '#fb923c' : (COLOR_PALETTE[this.tracks.length % COLOR_PALETTE.length]));
        const isFrontier = tid.startsWith('HOMO') || tid.startsWith('LUMO') || tid.startsWith('trk_');
        const newTrack = {
            id: tid, track_id: tid, orbital_key: tid, local_key: orbital.local_key || tid,
            chem_key: orbital.chem_key || tid, ref_frame: this.currentFrame, step_index: isFrontier ? null : orbital.step_index,
            orbital_class: orbital.orbital_class || (isAnti ? 'antibonding' : 'bonding'),
            type: orbital.type, palette: orbital.palette, title: orbital.label || orbital.name || tid,
            color: color, default_on: true,
            clips: [{ id: `c_${Date.now()}`, start: 0, end: Math.max(0, this.totalFrames - 1), title: orbital.label || orbital.name || tid }]
        };
        this.tracks.push(newTrack);
        this.manualOverrides.set(newTrack.id, true);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
        if (this.options.onVisibilityChange) this.options.onVisibilityChange();
    }

    addEspTrack() {
        if (this.tracks.some(t => t.track_id === 'trk_esp' || t.orbital_key === 'track_esp_charge')) return;
        const espTrack = {
            id: 'trk_esp', track_id: 'trk_esp', orbital_key: 'track_esp_charge', chem_key: 'esp_charge',
            orbital_class: 'esp', title: '⚡ Electrostatic Surface & Charges', color: '#38bdf8', default_on: true,
            clips: [{ id: 'c_esp', start: 0, end: Math.max(0, this.totalFrames - 1), title: 'ESP Surface' }]
        };
        this.tracks.unshift(espTrack);
        this.manualOverrides.set(espTrack.id, true);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
        if (this.options.onVisibilityChange) this.options.onVisibilityChange();
    }

    removeTrack(trackId) {
        this.tracks = this.tracks.filter(t => t.id !== trackId && t.track_id !== trackId);
        this.manualOverrides.delete(trackId);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
        if (this.options.onVisibilityChange) this.options.onVisibilityChange();
    }

    toggleTrackVisibility(trackId) {
        const trk = this.tracks.find(t => t.id === trackId || t.track_id === trackId);
        if (trk) {
            const cur = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            this.manualOverrides.set(trk.id, !cur);
            this.render();
            if (this.options.onVisibilityChange) this.options.onVisibilityChange();
        }
    }

    getActiveOrbitalsAtFrame(frameIdx) {
        const active = [];
        for (const trk of this.tracks) {
            const isVisible = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            if (!isVisible) continue;
            const tid = trk.track_id || trk.orbital_key || trk.id;
            if (tid === 'trk_esp' || tid === 'track_esp_charge') {
                active.push({ ...trk, track_id: 'trk_esp', orbital_key: 'track_esp_charge' });
                continue;
            }
            let activeKey = trk.orbital_key;
            if (trk.frame_map && Array.isArray(trk.frame_map)) {
                activeKey = trk.frame_map[frameIdx];
                if (!activeKey) continue;
            }
            active.push({ ...trk, track_id: tid, orbital_key: tid, local_key: activeKey });
        }
        return active;
    }

    getMutedOrbitalKeys() {
        const muted = [];
        for (const trk of this.tracks) {
            const isVisible = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            if (!isVisible) muted.push(trk.track_id || trk.orbital_key || trk.id);
        }
        return muted;
    }

    getTimelineObject() { return { default_isovalue: this.options.isovalue || 0.03, tracks: this.tracks }; }

    updatePopulations(frameIdx, frameData, symbols) {
        this.currentFrame = frameIdx;
        this.lastFrameData = frameData;
        if (symbols && symbols.length) this.symbols = symbols;
        if (!frameData || !this.containerEl) return;

        const allIbos = frameData.all_ibos || [];
        for (const trk of this.tracks) {
            const tid = trk.track_id || trk.orbital_key || trk.id;
            const el = this.containerEl.querySelector(`#pop-readout-${tid}`) || this.containerEl.querySelector(`#pop-readout-${trk.orbital_key}`);
            if (!el) continue;
            if (tid === 'trk_esp' || trk.orbital_key === 'track_esp_charge') {
                el.innerHTML = `<span style="color: #64748b;">F${frameIdx}:</span> <strong style="color: #0284c7;">Electrostatic Surface mapped with partial charges</strong>`;
                continue;
            }
            let activeKey = trk.orbital_key;
            if (trk.frame_map && Array.isArray(trk.frame_map)) {
                activeKey = trk.frame_map[frameIdx];
                if (!activeKey) {
                    el.innerHTML = `<span style="color: #94a3b8; font-style: italic;">Not present at frame ${frameIdx}</span>`;
                    continue;
                }
            }
            const ibo = allIbos.find(i => (i.track_id && i.track_id === tid) || i.key === activeKey || i.local_key === activeKey);
            if (ibo && ibo.pop && this.symbols.length) {
                const major = [];
                for (let i = 0; i < this.symbols.length; i++) {
                    const p = ibo.pop[i] || 0;
                    if (p > 0.15) major.push(`${this.symbols[i]}${i}: ${(p * 50).toFixed(0)}%`);
                }
                const roleTag = (ibo.label && ibo.label !== trk.title) ? ` <span style="font-size:0.68rem;color:#0284c7;">[${ibo.label}]</span>` : '';
                el.innerHTML = `<span style="color: #64748b;">Frame ${frameIdx}:</span> <strong style="color: #0f172a;">${major.join(' · ') || 'delocalized'}</strong>${roleTag}`;
            } else if (tid.startsWith('trk_homo') || tid.startsWith('trk_lumo') || trk.orbital_key.includes('STAR')) {
                el.innerHTML = `<span style="color: #64748b;">F${frameIdx}:</span> <span style="color: #334155; font-weight: 500;">Active across reaction</span>`;
            } else {
                el.innerHTML = `<span style="color: #64748b;">F${frameIdx}:</span> <span style="color: #64748b; font-style: italic;">Bond not formed at this frame</span>`;
            }
        }
    }

    render() {
        if (!this.containerEl) return;
        const activeCount = this.getActiveOrbitalsAtFrame(this.currentFrame).length;
        let cardsHtml = '';
        if (this.tracks.length === 0) {
            cardsHtml = `<div style="padding: 1rem; text-align: center; color: #94a3b8; font-size: 0.78rem;">No reaction orbitals selected. Click <strong>[+ Track]</strong> on any orbital from the pool to track it.</div>`;
        } else {
            cardsHtml = `<div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 0.5rem; padding: 0.2rem;">`;
            for (const trk of this.tracks) {
                const isChecked = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
                const tid = trk.track_id || trk.orbital_key || trk.id;
                const isEsp = tid === 'trk_esp' || trk.orbital_key === 'track_esp_charge', isAnti = isAntiBondingOrbital(trk), badgeInfo = getOrbitalTypeBadge(trk);
                const badgeIndicator = isEsp ? `
                    <div style="display: flex; align-items: center; gap: 2px; font-size: 0.65rem; font-weight: 700; flex-shrink: 0;" title="ESP Surface">
                        <span style="background: #ef4444; color: #fff; padding: 1px 4px; border-radius: 3px 0 0 3px;">δ−</span>
                        <span style="background: #3b82f6; color: #fff; padding: 1px 4px; border-radius: 0 3px 3px 0;">δ+</span>
                    </div>` : `
                    <div style="display: flex; align-items: center; gap: 2px; font-size: 0.65rem; font-weight: 700; flex-shrink: 0;" title="${isAnti ? 'Anti-Bonding' : 'Bonding'}">
                        <span style="background: ${isAnti ? '#ea580c' : '#0284c7'}; color: #fff; padding: 1px 4px; border-radius: 3px 0 0 3px; font-weight: 800;">+</span>
                        <span style="background: ${isAnti ? '#db2777' : '#1d4ed8'}; color: #fff; padding: 1px 4px; border-radius: 0 3px 3px 0; font-weight: 800;">−</span>
                    </div>`;
                const borderCol = isChecked ? (isEsp ? '#9333ea' : (isAnti ? '#ea580c' : '#0284c7')) : '#e2e8f0';
                const bgCol = isChecked ? (isEsp ? '#faf5ff' : (isAnti ? '#fff7ed' : '#f0f9ff')) : '#ffffff';
                const stepTag = trk.step_index ? `<span style="font-size:0.65rem;color:#64748b;background:#f1f5f9;padding:1px 4px;border-radius:3px;">Step ${trk.step_index}</span>` : '';
                cardsHtml += `
                    <div style="display: flex; align-items: center; justify-content: space-between; gap: 0.5rem; padding: 0.45rem 0.65rem; background: ${bgCol}; border: 1px solid ${borderCol}; border-radius: 6px; box-shadow: 0 1px 2px rgba(0,0,0,0.03);">
                        <div style="display: flex; align-items: center; gap: 0.5rem; flex: 1; min-width: 0;">
                            <input type="checkbox" class="track-vis-checkbox" data-id="${trk.id}" ${isChecked ? 'checked' : ''} style="cursor: pointer; width: 16px; height: 16px; accent-color: ${isAnti ? '#ea580c' : '#0284c7'}; margin: 0;">
                            ${badgeIndicator}
                            <div style="display: flex; flex-direction: column; min-width: 0; flex: 1;">
                                <div style="display: flex; align-items: center; gap: 0.35rem;">
                                    <span style="font-weight: 700; font-size: 0.82rem; color: #0f172a; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">${trk.title}</span>
                                    <span class="badge ${badgeInfo.cls}" style="font-size: 0.65rem; padding: 1px 4px; border-radius: 3px; ${badgeInfo.style}">${badgeInfo.label}</span>
                                    ${stepTag}
                                </div>
                                <div id="pop-readout-${tid}" style="font-size: 0.72rem; color: ${isAnti ? '#c2410c' : '#0284c7'}; margin-top: 2px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
                                    ${isEsp ? 'Electrostatic Surface mapped with atomic charges' : (isAnti ? 'Tracking antibonding acceptor orbital...' : 'Tracking bonding electron pair...')}
                                </div>
                            </div>
                        </div>
                        <button class="btn-track-remove" data-id="${trk.id}" title="Remove Orbital" style="background: none; border: none; color: #ef4444; cursor: pointer; font-size: 1rem; padding: 2px 5px; opacity: 0.8;">✕</button>
                    </div>`;
            }
            cardsHtml += `</div>`;
        }

        this.containerEl.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 0.4rem; padding: 0.5rem 0.75rem;">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <span style="font-weight: 700; font-size: 0.85rem; color: #0f172a;">Selected Reaction Orbitals</span>
                        <span class="badge" style="font-size: 0.7rem; background: #e0f2fe; color: #0284c7;">${activeCount} Active</span>
                    </div>
                    <span style="font-size: 0.72rem; color: #64748b;">Turn ON/OFF with checkbox · Live electron tracking across reaction</span>
                </div>
                ${cardsHtml}
            </div>`;

        this.containerEl.querySelectorAll('.track-vis-checkbox').forEach(cb => cb.addEventListener('change', e => this.toggleTrackVisibility(e.target.dataset.id)));
        this.containerEl.querySelectorAll('.btn-track-remove').forEach(btn => btn.addEventListener('click', e => this.removeTrack(e.currentTarget.dataset.id)));
        if (this.lastFrameData) this.updatePopulations(this.currentFrame, this.lastFrameData, this.symbols);
    }
}
