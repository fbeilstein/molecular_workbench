/**
 * viewer_modal.js
 * In-page modal dialog for inspecting individual molecules.
 * Features 3Dmol WebGL rendering, collapsible orbital checklist, ESP surface toggle,
 * and self-invoking "Open in Pop-up Window" (viewer.html?bundle=...&molecule=...).
 */

import { ViewerOrbitalsPanel } from './viewer_orbitals.js';

export class ViewerMoleculeModal {
  constructor(app) {
    this.app = app;
    this.modalEl = null;
    this.glViewer = null;
    this.orbitalsPanel = null;
    this.activeOrbitals = new Set();
    this.isEspActive = false;
    this.currentMoleculeData = null;
  }

  show(molId, cfg = {}) {
    this.close();

    const molName = cfg.moleculeName || molId;
    const bundle = this.app.bundle;
    const molData = this.resolveMoleculeData(molId, molName, bundle);
    this.currentMoleculeData = molData;

    const modal = document.createElement('div');
    modal.id = 'molecule-modal-backdrop';
    modal.className = 'custom-modal-backdrop';
    modal.innerHTML = `
      <div class="custom-modal-dialog">
        <div class="modal-header">
          <div class="modal-title-group">
            <h3>⚛️ ${molName}</h3>
            <span class="chip-blue">Molecule Inspection</span>
          </div>
          <div class="modal-actions">
            <button class="btn btn-sm" id="btn-modal-popout" title="Open in dedicated pop-up window">↗ Pop-up Window</button>
            <button class="btn-modal-close" id="btn-modal-close">&times;</button>
          </div>
        </div>
        <div class="modal-body">
          <div class="modal-viewport">
            <div id="modal-3d-canvas"></div>
            <div class="modal-view-tools">
              <button class="btn btn-sm" id="btn-modal-zoom-fit">⟲ Reset</button>
            </div>
          </div>
          <div class="modal-sidebar">
            <div class="modal-sidebar-header">
              <span>Orbitals & Surfaces</span>
            </div>
            <div id="modal-orbitals-container" class="modal-orbitals-list"></div>
          </div>
        </div>
      </div>
    `;

    document.body.appendChild(modal);
    this.modalEl = modal;

    // Close button
    modal.querySelector('#btn-modal-close').addEventListener('click', () => this.close());
    modal.addEventListener('click', (e) => {
      if (e.target === modal) this.close();
    });

    // Pop-up button: Calls viewer.html with bundle and molecule params!
    modal.querySelector('#btn-modal-popout').addEventListener('click', () => {
      const bundleUrl = this.app.currentBundleUrl || './bundles/tetrahedral_intermediate_CO.zip';
      const targetMol = molData?.key || molId;
      const popUrl = `viewer.html?bundle=${encodeURIComponent(bundleUrl)}&molecule=${encodeURIComponent(targetMol)}`;
      window.open(popUrl, '_blank', 'width=1000,height=750,menubar=no,toolbar=no');
    });

    // 3D Canvas
    const canvasContainer = modal.querySelector('#modal-3d-canvas');
    if (window.$3Dmol) {
      this.glViewer = window.$3Dmol.createViewer(canvasContainer, {
        backgroundColor: '0x090d16',
        id: 'modalGlViewer'
      });
    }

    modal.querySelector('#btn-modal-zoom-fit').addEventListener('click', () => {
      if (this.glViewer) {
        this.glViewer.zoomTo();
        this.glViewer.render();
      }
    });

    // Orbitals Panel
    const orbContainer = modal.querySelector('#modal-orbitals-container');
    this.orbitalsPanel = new ViewerOrbitalsPanel(
      orbContainer,
      (orbKey, isChecked) => this.toggleOrbital(orbKey, isChecked),
      (isEsp) => this.toggleEsp(isEsp)
    );

    this.renderMolecule(molData);
  }

  resolveMoleculeData(molId, molName, bundle) {
    if (!bundle) return null;
    const mb = (bundle.trajectory && bundle.trajectory.molecule_bundles) || {};
    const mols = bundle.molecules || {};

    let candidate = mb[molId] || mb[molName] || mols[molId] || mols[molName];
    if (!candidate) {
      for (const k in mb) {
        if (k.toLowerCase() === molId.toLowerCase() || k.toLowerCase() === molName.toLowerCase()) {
          candidate = mb[k];
          break;
        }
      }
    }
    if (!candidate) {
      for (const k in mols) {
        if (k.toLowerCase() === molId.toLowerCase() || k.toLowerCase() === molName.toLowerCase()) {
          candidate = mols[k];
          break;
        }
      }
    }

    if (candidate && !candidate.key) {
      candidate.key = candidate.molecule_id || candidate.molecule_key || molId;
    }
    return candidate || { key: molId, name: molName };
  }

  renderMolecule(molData) {
    if (!this.glViewer || !molData) return;
    this.activeOrbitals.clear();
    this.isEspActive = false;

    this.glViewer.clear();
    this.glViewer.removeAllSurfaces();
    this.glViewer.removeAllShapes();

    // 1. Render Structure
    if (molData.xyz) {
      this.glViewer.addModel(molData.xyz, 'xyz');
      this.glViewer.setStyle({}, {
        stick: { radius: 0.13, colorscheme: 'Jmol' },
        sphere: { scale: 0.28, colorscheme: 'Jmol' }
      });
      this.glViewer.zoomTo();
      this.glViewer.render();
    }

    // 2. Populate Orbitals
    const tree = molData.orbitals_tree || molData.orbitals || [];
    let categories = [];
    if (Array.isArray(tree) && tree.length > 0 && tree[0].items) {
      categories = tree;
    } else if (Array.isArray(tree)) {
      categories = ViewerOrbitalsPanel.parseRawOrbitalKeys(tree.map(t => t.label || t.file || t));
    }
    const hasEsp = !!(molData.esp_surface || molData[`${molData.key}_esp.json`]);
    this.orbitalsPanel.renderGroups(categories, hasEsp);
  }

  toggleOrbital(orbKey, isChecked) {
    if (isChecked) this.activeOrbitals.add(orbKey);
    else this.activeOrbitals.delete(orbKey);
    this.refreshSurfaces();
  }

  toggleEsp(isEsp) {
    this.isEspActive = isEsp;
    this.refreshSurfaces();
  }

  refreshSurfaces() {
    if (!this.glViewer || !this.currentMoleculeData) return;
    this.glViewer.removeAllSurfaces();
    this.glViewer.removeAllShapes();

    const molData = this.currentMoleculeData;
    const bundle = this.app.bundle;
    const molFiles = (bundle && bundle.molecules && bundle.molecules[molData.key]) || molData;

    // Render active orbitals
    for (const orbKey of this.activeOrbitals) {
      const cleanKey = orbKey.replace(/^molecules\//, '').replace(/\.json$/, '');
      const posObj = molFiles[`${cleanKey}_pos.json`] || molFiles[`${cleanKey}_pos`] || bundle?.orbitalMeshes?.[`${cleanKey}_pos`];
      const negObj = molFiles[`${cleanKey}_neg.json`] || molFiles[`${cleanKey}_neg`] || bundle?.orbitalMeshes?.[`${cleanKey}_neg`];

      if (posObj && posObj.vertexArr) {
        this.glViewer.addCustom({
          vertexArr: posObj.vertexArr, normalArr: posObj.normalArr, faceArr: posObj.faceArr,
          color: 0x38bdf8, alpha: 0.78
        });
      }
      if (negObj && negObj.vertexArr) {
        this.glViewer.addCustom({
          vertexArr: negObj.vertexArr, normalArr: negObj.normalArr, faceArr: negObj.faceArr,
          color: 0xf97316, alpha: 0.78
        });
      }
    }

    // Render ESP
    if (this.isEspActive) {
      const espObj = molData.esp_surface || molFiles[`${molData.key}_esp.json`] || molFiles['esp.json'];
      if (espObj && espObj.vertexArr && espObj.faceArr) {
        this.glViewer.addCustom({
          vertexArr: espObj.vertexArr, normalArr: espObj.normalArr || espObj.vertexArr, faceArr: espObj.faceArr,
          color: espObj.color || 0x10b981, alpha: 0.8
        });
      }
    }

    this.glViewer.render();
  }

  close() {
    if (this.modalEl) {
      this.modalEl.remove();
      this.modalEl = null;
    }
    this.glViewer = null;
  }
}
