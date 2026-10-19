/**
 * Shared constants and utility functions for Reaction Orbital Studio in IRC Pipeline.
 * Features dedicated distinct color schemes for Bonding vs Anti-Bonding orbitals.
 */

export const COLOR_PALETTE = [
    '#38bdf8', '#fb923c', '#10b981', '#ec4899', '#6366f1',
    '#f59e0b', '#3b82f6', '#84cc16', '#14b8a6', '#db2777'
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

export function isAntiBondingOrbital(orb) {
    if (!orb) return false;
    if (typeof orb === 'string') {
        const s = orb.toUpperCase();
        return s.includes('*') || s.includes('LUMO') || s.includes('STAR') || s.includes('ANTIBOND');
    }
    if (orb.orbital_class === 'antibonding') return true;
    const name = String(orb.name || orb.label || orb.orbital_key || orb.chem_key || orb.title || '').toUpperCase();
    const type = String(orb.type || '').toLowerCase();
    if (name.includes('*') || name.includes('LUMO') || type.includes('star') || type.includes('antibond')) {
        return true;
    }
    return false;
}

export function getOrbitalTypeBadge(orb) {
    const isAnti = isAntiBondingOrbital(orb);
    const type = String(orb.type || orb || '').toLowerCase();
    const name = String(orb.name || orb.label || '').toUpperCase();

    if (type.includes('esp')) {
        return { label: 'ESP', cls: 'badge-esp', style: 'background: rgba(168,85,247,0.2); color: #c084fc; border: 1px solid rgba(168,85,247,0.4);' };
    }
    if (isAnti) {
        const lbl = name.includes('LUMO') ? name : (type.includes('pi') ? 'π* Antibond' : 'σ* Antibond');
        return { label: lbl, cls: 'badge-antibonding', style: 'background: rgba(251,146,60,0.2); color: #fb923c; border: 1px solid rgba(251,146,60,0.5); font-weight: 700;' };
    }
    if (type.includes('core')) {
        return { label: 'Core', cls: 'badge-core', style: 'background: rgba(148,163,184,0.15); color: #94a3b8; border: 1px solid rgba(148,163,184,0.3);' };
    }
    if (type.includes('lone') || type.includes('lp')) {
        return { label: 'Lone Pair', cls: 'badge-bonding', style: 'background: rgba(56,189,248,0.2); color: #38bdf8; border: 1px solid rgba(56,189,248,0.4); font-weight: 600;' };
    }
    if (name.includes('HOMO')) {
        return { label: name, cls: 'badge-bonding', style: 'background: rgba(56,189,248,0.25); color: #38bdf8; border: 1px solid rgba(56,189,248,0.5); font-weight: 700;' };
    }
    return { label: 'Bonding', cls: 'badge-bonding', style: 'background: rgba(56,189,248,0.15); color: #38bdf8; border: 1px solid rgba(56,189,248,0.35); font-weight: 600;' };
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

/**
 * Returns distinct phase lobe colors for 3D isosurfaces.
 * Bonding: Cyan (+) / Royal Blue (-)
 * Anti-Bonding: Bright Orange (+) / Hot Magenta-Red (-)
 */
export function getPhaseColors(orb) {
    // If orb is an object with explicit palette from backend
    if (orb && orb.palette && orb.palette.pos && orb.palette.neg) {
        return orb.palette;
    }

    if (isAntiBondingOrbital(orb)) {
        return {
            pos: '#fb923c', // Bright Orange (+ phase)
            neg: '#db2777'  // Hot Magenta / Deep Red (- phase)
        };
    }

    // Default: Bonding Orbital Palette
    return {
        pos: '#38bdf8', // Light Cyan (+ phase)
        neg: '#1d4ed8'  // Deep Royal Blue (- phase)
    };
}
