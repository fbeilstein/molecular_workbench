/**
 * Shared constants and utility functions for Reaction Orbital Studio.
 */

export const COLOR_PALETTE = [
    '#06b6d4', '#a855f7', '#10b981', '#f59e0b', '#ec4899',
    '#3b82f6', '#84cc16', '#6366f1', '#14b8a6', '#f97316'
];

export const ELEMENT_COLORS = {
    H: 0xffffff,
    C: 0x909090,
    N: 0x3050f8,
    O: 0xff0d0d,
    F: 0x90e050,
    Cl: 0x1ff01f,
    Br: 0xa62929,
    I: 0x940094,
    S: 0xffff30,
    P: 0xff8000
};

export const ELEMENT_RADII = {
    H: 0.32, C: 0.77, N: 0.75, O: 0.73, F: 0.71,
    Cl: 0.99, Br: 1.14, I: 1.33, S: 1.02, P: 1.06
};

export function getOrbitalTypeBadge(type) {
    const t = String(type || '').toLowerCase();
    if (t.includes('esp')) return { label: 'ESP', cls: 'badge-esp' };
    if (t.includes('pi') && t.includes('star')) return { label: 'π*', cls: 'badge-pistar' };
    if (t.includes('sigma') && t.includes('star')) return { label: 'σ*', cls: 'badge-sigmastar' };
    if (t.includes('pi')) return { label: 'π', cls: 'badge-pi' };
    if (t.includes('sigma')) return { label: 'σ', cls: 'badge-sigma' };
    if (t.includes('lone') || t.includes('lp')) return { label: 'LP', cls: 'badge-lp' };
    if (t.includes('frontier') || t.includes('homo') || t.includes('lumo')) return { label: 'FMO', cls: 'badge-fmo' };
    return { label: 'MO', cls: 'badge-mo' };
}

export function getEspColor(charge) {
    const clamped = Math.max(-1.0, Math.min(1.0, charge));
    if (clamped < 0) {
        const intensity = Math.abs(clamped);
        return `rgb(${Math.round(255 * intensity)}, ${Math.round(50 * (1 - intensity))}, ${Math.round(50 * (1 - intensity))})`;
    } else {
        const intensity = clamped;
        return `rgb(${Math.round(50 * (1 - intensity))}, ${Math.round(50 * (1 - intensity))}, ${Math.round(255 * intensity)})`;
    }
}

export function getPhaseColors(baseColor) {
    const map = {
        '#ec4899': { pos: '#ec4899', neg: '#06b6d4' },
        '#f43f5e': { pos: '#f43f5e', neg: '#0284c7' },
        '#10b981': { pos: '#10b981', neg: '#a855f7' },
        '#06b6d4': { pos: '#06b6d4', neg: '#f97316' },
        '#3b82f6': { pos: '#3b82f6', neg: '#ef4444' },
        '#a855f7': { pos: '#a855f7', neg: '#10b981' },
        '#f59e0b': { pos: '#f59e0b', neg: '#3b82f6' },
        '#6366f1': { pos: '#6366f1', neg: '#eab308' },
        '#84cc16': { pos: '#84cc16', neg: '#ec4899' },
        '#14b8a6': { pos: '#14b8a6', neg: '#f43f5e' },
        '#f97316': { pos: '#f97316', neg: '#0ea5e9' }
    };
    return map[baseColor] || { pos: baseColor || '#0284c7', neg: '#f97316' };
}


