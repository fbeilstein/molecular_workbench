/**
 * 3Dmol viewer for Stage 4 IRC Trajectory playback.
 * Extends unified MoleculeViewer3D with true double bond rendering.
 */
import { MoleculeViewer3D } from './common/molecule_viewer_3d.js';

export class IRC3DViewer extends MoleculeViewer3D {
    constructor(container, options = {}) {
        super(container, options);
        this.init();
    }
}
