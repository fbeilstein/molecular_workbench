/**
 * Kdenlive-style multi-track timeline and clip sequencer for Reaction Orbital Studio.
 */
import { COLOR_PALETTE } from './orbitals_types.js';

export class OrbitalsTimeline {
    constructor(containerEl, options) {
        this.containerEl = containerEl;
        this.options = options || {};
        this.tracks = [];
        this.totalFrames = options.totalFrames || 0;
        this.currentFrame = options.currentFrame || 0;
        this.manualOverrides = new Map(); // trackId -> boolean
    }

    setTimelineData(timeline, totalFrames) {
        this.tracks = (timeline && timeline.tracks) ? JSON.parse(JSON.stringify(timeline.tracks)) : [];
        if (totalFrames) this.totalFrames = totalFrames;
        this.render();
    }

    setCurrentFrame(frameIdx) {
        this.currentFrame = Math.max(0, Math.min(this.totalFrames - 1, frameIdx));
        this.updatePlayhead();
    }

    addTrack(orbital) {
        if (!orbital) return;
        const color = COLOR_PALETTE[this.tracks.length % COLOR_PALETTE.length];
        const newTrack = {
            id: `trk_${Date.now()}_${Math.floor(Math.random() * 1000)}`,
            orbital_key: orbital.key,
            chem_key: orbital.chem_key || orbital.key,
            title: orbital.label || orbital.key,
            color: orbital.color || color,
            default_on: true,
            clips: [{
                id: `c_${Date.now()}`,
                start: 0,
                end: Math.max(0, this.totalFrames - 1),
                title: orbital.label || orbital.key
            }]
        };
        this.tracks.push(newTrack);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
    }

    addEspTrack() {
        if (this.tracks.some(t => t.orbital_key === 'track_esp_charge')) return;
        const espTrack = {
            id: 'trk_esp_charge',
            orbital_key: 'track_esp_charge',
            chem_key: 'esp_charge',
            title: '⚡ Electrostatic Surface & Charges',
            color: '#38bdf8',
            default_on: false,
            clips: [{
                id: 'c_esp',
                start: 0,
                end: Math.max(0, this.totalFrames - 1),
                title: 'ESP Surface'
            }]
        };
        this.tracks.unshift(espTrack);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
    }

    removeTrack(trackId) {
        this.tracks = this.tracks.filter(t => t.id !== trackId);
        this.render();
        if (this.options.onTimelineChange) this.options.onTimelineChange(this.getTimelineObject());
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
            if (!isVisible) continue;
            // Check if frame falls within any clip
            const inClip = trk.clips && trk.clips.some(c => frameIdx >= c.start && frameIdx <= c.end);
            if (inClip) {
                active.push(trk);
            }
        }
        return active;
    }

    getMutedOrbitalKeys() {
        const muted = [];
        for (const trk of this.tracks) {
            const isVisible = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            if (!isVisible) {
                if (trk.orbital_key) muted.push(trk.orbital_key);
                if (trk.chem_key) muted.push(trk.chem_key);
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

    updatePlayhead() {
        const needle = this.containerEl ? this.containerEl.querySelector('#timeline-playhead-needle') : null;
        if (!needle || this.totalFrames <= 1) return;
        const frac = this.currentFrame / (this.totalFrames - 1);
        needle.style.left = `calc(190px + (100% - 190px) * ${frac.toFixed(4)})`;
        const timeReadout = this.containerEl.querySelector('#timeline-frame-indicator');
        if (timeReadout) timeReadout.textContent = `Frame ${this.currentFrame} / ${this.totalFrames - 1}`;
    }

    render() {
        if (!this.containerEl) return;
        const total = Math.max(1, this.totalFrames);

        let rulerHtml = '';
        const step = Math.max(1, Math.floor(total / 10));
        for (let i = 0; i < total; i += step) {
            const pct = ((i / (total - 1 || 1)) * 100).toFixed(1);
            rulerHtml += `<div style="position: absolute; left: ${pct}%; transform: translateX(-50%); font-size: 0.65rem; color: #64748b;">${i}</div>`;
        }

        let tracksHtml = '';
        for (const trk of this.tracks) {
            const isVisible = this.manualOverrides.has(trk.id) ? this.manualOverrides.get(trk.id) : trk.default_on;
            let clipsHtml = '';
            if (trk.clips) {
                for (const clip of trk.clips) {
                    const lPct = ((clip.start / (total - 1 || 1)) * 100).toFixed(1);
                    const wPct = (((clip.end - clip.start + 1) / (total - 1 || 1)) * 100).toFixed(1);
                    clipsHtml += `
                        <div style="position: absolute; left: ${lPct}%; width: ${wPct}%; top: 2px; bottom: 2px; background: ${trk.color}44; border: 1px solid ${trk.color}; border-radius: 3px; display: flex; align-items: center; padding: 0 4px; overflow: hidden;">
                            <span style="font-size: 0.65rem; color: #fff; white-space: nowrap; text-overflow: ellipsis; overflow: hidden;">${clip.title}</span>
                        </div>
                    `;
                }
            }

            tracksHtml += `
                <div class="timeline-track-row" style="display: flex; height: 32px; border-bottom: 1px solid rgba(255,255,255,0.06); align-items: center;">
                    <div style="width: 190px; padding: 0 0.5rem; display: flex; align-items: center; justify-content: space-between; background: rgba(0,0,0,0.2); height: 100%;">
                        <div style="display: flex; align-items: center; gap: 0.4rem; overflow: hidden;">
                            <button class="btn-track-vis" data-id="${trk.id}" style="background: none; border: none; cursor: pointer; color: ${isVisible ? '#38bdf8' : '#64748b'}; font-size: 0.8rem; padding: 0;">
                                ${isVisible ? '👁' : '🚫'}
                            </button>
                            <span style="font-size: 0.75rem; font-weight: 500; color: #f1f5f9; white-space: nowrap; overflow: hidden; text-overflow: ellipsis;" title="${trk.title}">${trk.title}</span>
                        </div>
                        <button class="btn-track-del" data-id="${trk.id}" style="background: none; border: none; cursor: pointer; color: #ef4444; font-size: 0.75rem; padding: 0;">✕</button>
                    </div>
                    <div class="track-canvas" data-id="${trk.id}" style="flex: 1; height: 100%; position: relative; background: rgba(15,23,42,0.6);">
                        ${clipsHtml}
                    </div>
                </div>
            `;
        }

        this.containerEl.innerHTML = `
            <div style="display: flex; justify-content: space-between; align-items: center; padding: 0.4rem 0.6rem; background: rgba(0,0,0,0.3); border-bottom: 1px solid rgba(255,255,255,0.08);">
                <div style="display: flex; gap: 0.5rem; align-items: center;">
                    <span style="font-size: 0.8rem; font-weight: 600; color: #e2e8f0;">Timeline Tracks</span>
                    <button id="btn-add-esp-track" class="btn btn-xs btn-outline" style="font-size: 0.7rem; color: #38bdf8; border-color: #0284c7;">+ Add ESP Surface</button>
                </div>
                <div id="timeline-frame-indicator" style="font-size: 0.72rem; color: #94a3b8;">Frame ${this.currentFrame} / ${total - 1}</div>
            </div>
            <div style="display: flex; height: 20px; background: rgba(0,0,0,0.4); border-bottom: 1px solid rgba(255,255,255,0.08);">
                <div style="width: 190px; border-right: 1px solid rgba(255,255,255,0.06); font-size: 0.65rem; color: #64748b; padding-left: 0.5rem; line-height: 20px;">Ruler</div>
                <div id="timeline-ruler" style="flex: 1; position: relative; cursor: pointer;">
                    ${rulerHtml}
                </div>
            </div>
            <div id="timeline-tracks-body" style="position: relative; max-height: 180px; overflow-y: auto;">
                <div id="timeline-playhead-needle" style="position: absolute; top: 0; bottom: 0; width: 2px; background: #ef4444; z-index: 10; pointer-events: none; transform: translateX(-50%);">
                    <div style="width: 10px; height: 10px; background: #ef4444; transform: rotate(45deg); position: absolute; top: -5px; left: -4px;"></div>
                </div>
                ${tracksHtml || '<div class="empty-state" style="padding: 1rem; text-align: center; color: #64748b;">No tracks added. Select an orbital from the pool or add ESP surface.</div>'}
            </div>
        `;

        this.updatePlayhead();
        this.bindEvents();
    }

    bindEvents() {
        const ruler = this.containerEl.querySelector('#timeline-ruler');
        if (ruler) {
            ruler.addEventListener('click', (e) => {
                const rect = ruler.getBoundingClientRect();
                const x = e.clientX - rect.left;
                const pct = Math.max(0, Math.min(1, x / rect.width));
                const targetF = Math.round(pct * (this.totalFrames - 1));
                if (this.options.onSeek) this.options.onSeek(targetF);
            });
        }

        const btnEsp = this.containerEl.querySelector('#btn-add-esp-track');
        if (btnEsp) {
            btnEsp.addEventListener('click', () => this.addEspTrack());
        }

        this.containerEl.querySelectorAll('.btn-track-vis').forEach(b => {
            b.addEventListener('click', (e) => {
                e.stopPropagation();
                this.toggleTrackVisibility(b.dataset.id);
            });
        });

        this.containerEl.querySelectorAll('.btn-track-del').forEach(b => {
            b.addEventListener('click', (e) => {
                e.stopPropagation();
                this.removeTrack(b.dataset.id);
            });
        });
    }
}
