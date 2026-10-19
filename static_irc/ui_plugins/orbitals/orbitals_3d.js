/**
 * 3Dmol scene and volumetric isosurface renderer for Step 6 Orbitals.
 * Extends unified MoleculeViewer3D with authentic double bond rendering,
 * camera preservation, and dual-phase isosurface coloring.
 */
import { MoleculeViewer3D } from '../common/molecule_viewer_3d.js';
import { getEspColor, getPhaseColors, isAntiBondingOrbital } from './orbitals_types.js';

export class Orbitals3DViewer extends MoleculeViewer3D {
    constructor(containerEl, statusEl, options = {}) {
        super(containerEl, { statusEl, ...options });
    }
}
