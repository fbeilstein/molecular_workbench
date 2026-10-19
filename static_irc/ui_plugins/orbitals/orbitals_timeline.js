/**
 * Reaction Orbital Selection & Real-Time Electron Tracking Panel.
 * Streamlined orbital choice: toggles orbitals ON/OFF via simple checkboxes with zero clip dragging.
 * Real-time electron tracking displays electron population distribution across atoms per frame.
 * Visually distinguishes bonding (Cyan/Blue) vs anti-bonding (Orange/Magenta) tracks.
 */
import { COLOR_PALETTE, getPhaseColors, isAntiBondingOrbital, getOrbitalTypeBadge } from './orbitals_types.js';

export class OrbitalsTimeline {
    constructor(containerEl, options) {
        this.containerEl = containerEl;
        this.options = options || {};
        this.tracks = [];
        this.totalFrames = options.totalFrames || 0;
        this.currentFrame = options.currentFrame || 0;
        this.manualOverrides = new Map(); // trackId -> boolean
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
        if (this.lastFrameData) {
            this.updatePopulations(this.currentFrame, this.lastFrameData, this.symbols);
        }
    }

    addTrack(orbital) {
        if (!orbital) return;
        if (orbital.key === 'track_esp_charge') {
            this.addEspTrack();
            return;
        }
        const exists = this.tracks.find(t => t.orbital_key === orbital.key);
        if (exists) {
            this.manualOverrides.set(exists.id, true);
            this.render();
            if (this.options.onVisibilityChange) this.options.onVisibilityChange();
            return;
        }
        const isAnti = isAntiBondingOrbital(orbital);
        const color = orbital.color || (isAnti ? '#fb923c' : (COLOR_PALETTE[this.tracks.length % COLOR_PALETTE.length]));
        const newTrack = {
            id: `trk_${Date.now()}_${Math.floor(Math.random() * 1000)}`,
            orbital_key: orbital.key,
            chem_key: orbital.chem_key || orbital.key,
            orbital_class: orbital.orbital_class || (isAnti ? 'antibonding' : 'bonding'),
            type: orbital.type,
            palette: orbital.palette,
            title: orbital.label || orbital.key,
            color: color,
            default_on: true,
            clips: [{
                id: `c_${Date.now()}`,
                start: 0,
                end: Math.max(0, this.totalFrames - 1),
                title: orbital.label || orbital.key
            }]
        };
        this.tracks.push(newTrack);
        this.manualOverrides.set(newTrack.id, true);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
        if (this.options.onVisibilityChange) this.options.onVisibilityChange();
    }

    addEspTrack() {
        if (this.tracks.some(t => t.orbital_key === 'track_esp_charge')) return;
        const espTrack = {
            id: 'trk_esp_charge',
            orbital_key: 'track_esp_charge',
            chem_key: 'esp_charge',
            orbital_class: 'esp',
            title: '⚡ Electrostatic Surface & Charges',
            color: '#38bdf8',
            default_on: true,
            clips: [{ id: 'c_esp', start: 0, end: Math.max(0, this.totalFrames - 1), title: 'ESP Surface' }]
        };
        this.tracks.unshift(espTrack);
        this.manualOverrides.set(espTrack.id, true);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
        if (this.options.onVisibilityChange) this.options.onVisibilityChange();
    }

    removeTrack(trackId) {
        this.tracks = this.tracks.filter(t => t.id !== trackId);
        this.manualOverrides.delete(trackId);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
        if (this.options.onVisibilityChange) this.options.onVisibilityChange();
    }

    toggleTrackVisibility(trackId) {
        const trk = this.tracks.find(t => t.id === trackId);
        if (trk) {
            const cur = this.manualOverrides.has(trackId) ? this.manualOverrides.get(trackId) : trk.default_on;
            this.manualOverrides.set(trackId, !cur);
            this.render();
            if (this.options.onVisibilityChange) this.options.onVisibilityChange();
        }
    }

    getActiveOrbitalsAtFrame(frameIdx) {
        const active = [];
        for (const trk of this.tracks) {
            const isVisible = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            if (isVisible) active.push(trk);
        }
        return active;
    }

    getMutedOrbitalKeys() {
        const muted = [];
        for (const trk of this.tracks) {
            const isVisible = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            if (!isVisible) {
                if (trk.orbital_key) muted.push(trk.orbital_key);
            }
        }
        return muted;
    }

    getTimelineObject() {
        return {
            default_isovalue: this.options.isovalue || 0.03,
            tracks: this.tracks
        };
    }

    updatePopulations(frameIdx, frameData, symbols) {
        this.currentFrame = frameIdx;
        this.lastFrameData = frameData;
        if (symbols && symbols.length) this.symbols = symbols;
        if (!frameData || !this.containerEl) return;

        const iboMap = new Map();
        if (frameData.all_ibos) {
            for (const ibo of frameData.all_ibos) {
                iboMap.set(ibo.name || ibo.key, ibo);
            }
        }

        for (const trk of this.tracks) {
            const el = this.containerEl.querySelector(`#pop-readout-${trk.orbital_key}`);
            if (!el) continue;
            if (trk.orbital_key === 'track_esp_charge') {
                el.innerHTML = `<span style="color: #94a3b8;">F${frameIdx}:</span> <strong style="color: #38bdf8;">Electrostatic Surface mapped with partial charges</strong>`;
                continue;
            }
            const ibo = iboMap.get(trk.orbital_key);
            if (ibo && ibo.pop && this.symbols.length) {
                const major = [];
                for (let i = 0; i < this.symbols.length; i++) {
                    const p = ibo.pop[i] || 0;
                    if (p > 0.15) major.push(`${this.symbols[i]}${i}: ${(p * 50).toFixed(0)}%`);
                }
                const curLabel = ibo.label || trk.title;
                el.innerHTML = `<span style="color: #94a3b8;">F${frameIdx} [${curLabel}]:</span> <strong style="color: #e2e8f0;">${major.join(' · ') || 'delocalized'}</strong>`;
            } else if (trk.orbital_key.includes('HOMO') || trk.orbital_key.includes('LUMO') || trk.orbital_key.includes('STAR')) {
                el.innerHTML = `<span style="color: #94a3b8;">F${frameIdx}:</span> <span style="color: #cbd5e1;">Active across reaction</span>`;
            }
        }
    }

    render() {
        if (!this.containerEl) return;
        const activeCount = this.getActiveOrbitalsAtFrame(this.currentFrame).length;

        let cardsHtml = '';
        if (this.tracks.length === 0) {
            cardsHtml = `
                <div style="padding: 1rem; text-align: center; color: #94a3b8; font-size: 0.78rem;">
                    No reaction orbitals selected. Click <strong>[+ Track]</strong> on any orbital or ESP from the pool on the left to track it.
                </div>
            `;
        } else {
            cardsHtml = `<div style="display: grid; grid-template-columns: repeat(auto-fill, minmax(280px, 1fr)); gap: 0.5rem; padding: 0.2rem;">`;
            for (const trk of this.tracks) {
                const isChecked = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
                const isEsp = trk.orbital_key === 'track_esp_charge';
                const isAnti = isAntiBondingOrbital(trk);
                const pc = getPhaseColors(trk);
                const badgeInfo = getOrbitalTypeBadge(trk);

                const badgeIndicator = isEsp ? `
                    <div style="display: flex; align-items: center; gap: 2px; font-size: 0.65rem; font-weight: 700; flex-shrink: 0;" title="Electrostatic Potential (δ- red, δ+ blue)">
                        <span style="background: #ef4444; color: #fff; padding: 1px 4px; border-radius: 3px 0 0 3px;">δ−</span>
                        <span style="background: #3b82f6; color: #fff; padding: 1px 4px; border-radius: 0 3px 3px 0;">δ+</span>
                    </div>
                ` : `
                    <div style="display: flex; align-items: center; gap: 2px; font-size: 0.65rem; font-weight: 700; flex-shrink: 0;" title="${isAnti ? 'Anti-Bonding: Orange (+) / Magenta (-)' : 'Bonding: Cyan (+) / Royal Blue (-)'}">
                        <span style="background: ${pc.pos}; color: #000; padding: 1px 4px; border-radius: 3px 0 0 3px; font-weight: 800;">+</span>
                        <span style="background: ${pc.neg}; color: #fff; padding: 1px 4px; border-radius: 0 3px 3px 0; font-weight: 800;">−</span>
                    </div>
                `;

                const borderCol = isChecked ? (isEsp ? '#38bdf8' : (isAnti ? '#fb923c' : '#38bdf8')) : 'rgba(255,255,255,0.06)';
                const bgCol = isChecked ? (isEsp ? 'rgba(56,189,248,0.12)' : (isAnti ? 'rgba(251,146,60,0.12)' : 'rgba(56,189,248,0.08)')) : 'rgba(15,23,42,0.4)';
                const subtitle = isEsp ? 'Electrostatic Surface mapped with atomic charges' : (isAnti ? 'Tracking empty/antibonding acceptor orbital...' : 'Tracking bonding electron pair across pathway...');

                cardsHtml += `
                    <div style="display: flex; align-items: center; justify-content: space-between; gap: 0.5rem; padding: 0.45rem 0.65rem; background: ${bgCol}; border: 1px solid ${borderCol}; border-radius: 6px; transition: border-color 0.15s ease;">
                        <div style="display: flex; align-items: center; gap: 0.5rem; flex: 1; min-width: 0;">
                            <input type="checkbox" class="track-vis-checkbox" data-id="${trk.id}" ${isChecked ? 'checked' : ''} style="cursor: pointer; width: 16px; height: 16px; accent-color: ${isAnti ? '#fb923c' : '#38bdf8'}; margin: 0;">
                            ${badgeIndicator}
                            <div style="display: flex; flex-direction: column; min-width: 0; flex: 1;">
                                <div style="display: flex; align-items: center; gap: 0.35rem;">
                                    <span style="font-weight: 700; font-size: 0.82rem; color: #f8fafc; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">${trk.title}</span>
                                    <span style="font-size: 0.65rem; padding: 1px 4px; border-radius: 3px; ${badgeInfo.style}">${badgeInfo.label}</span>
                                </div>
                                <div id="pop-readout-${trk.orbital_key}" style="font-size: 0.72rem; color: ${isAnti ? '#fb923c' : '#38bdf8'}; margin-top: 2px; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;">
                                    ${subtitle}
                                </div>
                            </div>
                        </div>
                        <button class="btn-track-remove" data-id="${trk.id}" title="Remove Orbital" style="background: none; border: none; color: #ef4444; cursor: pointer; font-size: 1rem; padding: 2px 5px; opacity: 0.8;">✕</button>
                    </div>
                `;
            }
            cardsHtml += `</div>`;
        }

        this.containerEl.innerHTML = `
            <div style="display: flex; flex-direction: column; gap: 0.4rem; padding: 0.5rem 0.75rem;">
                <div style="display: flex; justify-content: space-between; align-items: center;">
                    <div style="display: flex; align-items: center; gap: 0.5rem;">
                        <span style="font-weight: 700; font-size: 0.85rem; color: #38bdf8;">Selected Reaction Orbitals</span>
                        <span class="badge" style="font-size: 0.7rem; background: rgba(56,189,248,0.15); color: #38bdf8;">${activeCount} Active</span>
                    </div>
                    <span style="font-size: 0.72rem; color: #94a3b8;">Turn ON/OFF with checkbox · Live electron tracking across reaction</span>
                </div>
                ${cardsHtml}
            </div>
        `;

        this.containerEl.querySelectorAll('.track-vis-checkbox').forEach(cb => {
            cb.addEventListener('change', (e) => this.toggleTrackVisibility(e.target.dataset.id));
        });
        this.containerEl.querySelectorAll('.btn-track-remove').forEach(btn => {
            btn.addEventListener('click', (e) => this.removeTrack(e.currentTarget.dataset.id));
        });

        if (this.lastFrameData) {
            this.updatePopulations(this.currentFrame, this.lastFrameData, this.symbols);
        }
    }
}
