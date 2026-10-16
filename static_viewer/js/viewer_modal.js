/**
 * viewer_modal.js
 * In-page modal dialog for inspecting individual molecules.
 * Reuses Viewer3D and ViewerOrbitalsPanel with identical collapsible drawer UI.
 */

import { buildMoleculeBundle } from './viewer_loader.js';
import { Viewer3D } from './viewer_3d.js';
import { ViewerOrbitalsPanel } from './viewer_orbitals.js';

export class ViewerMoleculeModal {
  constructor(app) {
    this.app = app;
    this.modalEl = null;
    this.viewer3d = null;
    this.orbitalsPanel = null;
    this.isPanelCollapsed = false;
  }

  show(molId, cfg = {}) {
    this.close();

    const parentBundle = this.app.bundle;
    const molBundle = buildMoleculeBundle(molId, parentBundle) ||
                      buildMoleculeBundle(cfg.moleculeName, parentBundle);

    if (!molBundle) {
      console.warn(`Molecule '${molId}' could not be resolved in bundle.`);
      return;
    }

    const title = molBundle.manifest.title || cfg.moleculeName || molId;

    const modal = document.createElement('div');
    modal.id = 'molecule-modal-backdrop';
    modal.className = 'custom-modal-backdrop';
    modal.innerHTML = `
      <div class="custom-modal-dialog">
        <div class="modal-header">
          <div class="modal-title-group">
            <h3>⚛️ ${title}</h3>
            <span class="chip-blue">Molecule Inspection</span>
          </div>
          <div class="modal-actions">
            <button class="btn btn-sm" id="btn-modal-popout" title="Open in dedicated pop-up window">↗ Pop-up Window</button>
            <button class="btn-modal-close" id="btn-modal-close">&times;</button>
          </div>
        </div>
        <div class="modal-viewport-container">
          <div id="modal-gl-canvas" class="modal-gl-canvas"></div>
          <div class="viewport-overlay-controls">
            <button class="btn-panel-toggle" id="modal-btn-toggle-orbitals">Orbitals ▾</button>
          </div>
          <aside id="modal-side-panel" class="side-panel">
            <div class="panel-header">
              <span>Orbitals & Surfaces</span>
            </div>
            <div class="panel-body">
              <div id="modal-orbitals-list" class="orbitals-list-container"></div>
            </div>
          </aside>
        </div>
      </div>
    `;

    document.body.appendChild(modal);
    this.modalEl = modal;

    // Close button and backdrop click
    modal.querySelector('#btn-modal-close').addEventListener('click', () => this.close());
    modal.addEventListener('click', (e) => {
      if (e.target === modal) this.close();
    });

    // Pop-up window button: Opens standalone single-molecule viewer
    modal.querySelector('#btn-modal-popout').addEventListener('click', () => {
      const bundleUrl = this.app.currentBundleUrl || './bundles/tetrahedral_intermediate_CO.zip';
      const targetMol = molBundle.manifest.molKey || molId;
      const popUrl = `viewer.html?bundle=${encodeURIComponent(bundleUrl)}&molecule=${encodeURIComponent(targetMol)}`;
      window.open(popUrl, '_blank', 'width=1100,height=750,menubar=no,toolbar=no');
    });

    // Collapsible drawer toggle
    const btnToggle = modal.querySelector('#modal-btn-toggle-orbitals');
    const sidePanel = modal.querySelector('#modal-side-panel');
    this.isPanelCollapsed = false;
    btnToggle.addEventListener('click', () => {
      this.isPanelCollapsed = !this.isPanelCollapsed;
      sidePanel.classList.toggle('panel-hidden', this.isPanelCollapsed);
      btnToggle.innerHTML = `Orbitals ${this.isPanelCollapsed ? '▸' : '▾'}`;
    });

    // Reusable 3D Viewer
    const canvasContainer = modal.querySelector('#modal-gl-canvas');
    this.viewer3d = new Viewer3D(canvasContainer);
    this.viewer3d.setBundle(molBundle);

    // Reusable Orbitals Panel
    const orbContainer = modal.querySelector('#modal-orbitals-list');
    this.orbitalsPanel = new ViewerOrbitalsPanel(
      orbContainer,
      (orbKey, isChecked) => this.viewer3d?.toggleOrbital(orbKey, isChecked),
      (isEsp) => this.viewer3d?.setEspActive(isEsp)
    );
    this.orbitalsPanel.renderGroups(molBundle.categories, molBundle.hasEsp);
  }

  close() {
    if (this.modalEl) {
      this.modalEl.remove();
      this.modalEl = null;
    }
    this.viewer3d = null;
    this.orbitalsPanel = null;
  }
}
