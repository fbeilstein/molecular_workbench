/**
 * viewer_orbitals.js
 * Grouped, collapsible orbital checklist matching Molecular Explorer design.
 * Supports multiple active orbital checkboxes and embedded ESP toggle.
 */

class ViewerOrbitalsPanel {
  constructor(containerElement, onOrbitalToggle, onEspToggle, onChargesToggle) {
    this.container = containerElement;
    this.onOrbitalToggle = onOrbitalToggle || (() => {});
    this.onEspToggle = onEspToggle || (() => {});
    this.onChargesToggle = onChargesToggle || (() => {});
    this.activeOrbitals = new Set();
    this.isEspActive = false;
    this.isChargesActive = false;
    this.collapsedGroups = new Set();
  }

  clear() {
    this.container.innerHTML = '';
    this.activeOrbitals.clear();
    this.isEspActive = false;
    this.isChargesActive = false;
  }

  /**
   * Render grouped categories from molecule JSON or reaction pathway tracks
   */
  renderGroups(categories, hasEsp = false) {
    this.container.innerHTML = '';
    if (!categories || categories.length === 0) {
      this.container.innerHTML = '<span class="empty-hint">No orbitals available</span>';
      return;
    }

    // 1. Render orbital category groups
    categories.forEach((cat, catIdx) => {
      const groupEl = document.createElement('div');
      groupEl.className = 'orbital-group';
      const isCollapsed = this.collapsedGroups.has(cat.name);

      const header = document.createElement('div');
      header.className = 'orbital-group-header';
      header.innerHTML = `
        <span class="group-title">${cat.name}</span>
        <span class="group-chevron">${isCollapsed ? '▸' : '▾'}</span>
      `;
      header.addEventListener('click', () => {
        if (this.collapsedGroups.has(cat.name)) {
          this.collapsedGroups.delete(cat.name);
        } else {
          this.collapsedGroups.add(cat.name);
        }
        this.renderGroups(categories, hasEsp);
      });
      groupEl.appendChild(header);

      const body = document.createElement('div');
      body.className = 'orbital-group-body';
      if (isCollapsed) body.style.display = 'none';

      (cat.items || []).forEach(item => {
        if (item.isEsp) {
          const row = document.createElement('div');
          row.className = 'orbital-item-row';
          row.style.cssText = 'display: flex; align-items: center; justify-content: space-between; gap: 0.35rem; padding: 0.25rem 0.4rem;';

          const left = document.createElement('div');
          left.style.cssText = 'display: flex; align-items: center; gap: 0.35rem; overflow: hidden; min-width: 0;';

          // Surface checkbox & label
          const surfBox = document.createElement('label');
          surfBox.style.cssText = 'display: flex; align-items: center; gap: 3px; cursor: pointer; margin: 0; user-select: none;';
          surfBox.title = 'Toggle 3D ESP continuous surface';
          const surfChk = document.createElement('input');
          surfChk.type = 'checkbox';
          surfChk.className = 'orb-chk';
          surfChk.checked = this.isEspActive;
          surfChk.style.accentColor = '#06b6d4';
          surfChk.addEventListener('change', (e) => {
            this.isEspActive = e.target.checked;
            this.onEspToggle(this.isEspActive);
          });
          const surfLbl = document.createElement('span');
          surfLbl.style.cssText = 'font-size: 0.68rem; font-weight: 600; color: #0284c7;';
          surfLbl.textContent = 'Surface';
          surfBox.appendChild(surfChk);
          surfBox.appendChild(surfLbl);

          // Divider
          const div = document.createElement('span');
          div.style.cssText = 'width: 1px; height: 12px; background: #cbd5e1; margin: 0 1px;';

          // Charges checkbox & label
          const chgBox = document.createElement('label');
          chgBox.style.cssText = 'display: flex; align-items: center; gap: 3px; cursor: pointer; margin: 0; user-select: none;';
          chgBox.title = 'Toggle atom partial charges (±q)';
          const chgChk = document.createElement('input');
          chgChk.type = 'checkbox';
          chgChk.className = 'orb-chk';
          chgChk.checked = this.isChargesActive;
          chgChk.style.accentColor = '#a855f7';
          chgChk.addEventListener('change', (e) => {
            this.isChargesActive = e.target.checked;
            this.onChargesToggle(this.isChargesActive);
          });
          const chgLbl = document.createElement('span');
          chgLbl.style.cssText = 'font-size: 0.68rem; font-weight: 600; color: #9333ea;';
          chgLbl.textContent = '±q';
          chgBox.appendChild(chgChk);
          chgBox.appendChild(chgLbl);

          const title = document.createElement('span');
          title.style.cssText = 'font-size: 0.74rem; font-weight: 600; color: var(--text-main, #0f172a); white-space: nowrap; margin-left: 2px;';
          title.textContent = 'ESP & Charges';

          left.appendChild(surfBox);
          left.appendChild(div);
          left.appendChild(chgBox);
          left.appendChild(title);

          const badge = document.createElement('span');
          badge.className = 'chip-mini badge-esp';
          badge.style.cssText = 'font-size: 0.62rem; padding: 1px 6px; border-radius: 3px; font-weight: 700; background: rgba(168,85,247,0.2); color: #c084fc; border: 1px solid rgba(168,85,247,0.4);';
          badge.textContent = 'ESP';

          row.appendChild(left);
          row.appendChild(badge);
          body.appendChild(row);
          return;
        }

        const row = document.createElement('label');
        row.className = 'orbital-item-row';

        const chk = document.createElement('input');
        chk.type = 'checkbox';
        chk.className = 'orb-chk';
        chk.checked = this.activeOrbitals.has(item.key || item.file || item.label);
        chk.addEventListener('change', (e) => {
          const k = item.key || item.file || item.label;
          if (e.target.checked) this.activeOrbitals.add(k);
          else this.activeOrbitals.delete(k);
          this.onOrbitalToggle(k, e.target.checked, item);
        });

        const dot = document.createElement('span');
        dot.className = 'orb-dot';
        const fKeyLower = String(item.key || item.file || '').toLowerCase();
        const isAntiItem = item.isAntiBonding || fKeyLower.includes('star') || fKeyLower.includes('*') || fKeyLower.includes('lumo') || fKeyLower.includes('antibond') || fKeyLower.includes('virtual');
        dot.style.backgroundColor = isAntiItem ? '#fb923c' : '#38bdf8';

        const labelSpan = document.createElement('span');
        labelSpan.className = 'orb-label-text';
        let txt = item.label || item.title || item.key;
        if (item.energy_ev !== undefined) {
          txt += ` (${item.energy_ev > 0 ? '+' : ''}${item.energy_ev.toFixed(2)} eV)`;
        }
        labelSpan.textContent = txt;

        row.appendChild(chk);
        row.appendChild(dot);
        row.appendChild(labelSpan);
        body.appendChild(row);
      });

      groupEl.appendChild(body);
      this.container.appendChild(groupEl);
    });

    // 2. Render ESP & Charges section (for single-molecule bundles, or if not embedded in tracks)
    const hasEspInTracks = categories.some(cat => (cat.items || []).some(it => it.isEsp));
    if (hasEspInTracks) return;

    const groupTitle = hasEsp ? 'Surfaces & Charges' : 'Atom Charges';
    const espGroup = document.createElement('div');
    espGroup.className = 'orbital-group esp-group';
    const isEspCollapsed = this.collapsedGroups.has(groupTitle);

    const espHeader = document.createElement('div');
    espHeader.className = 'orbital-group-header';
    espHeader.innerHTML = `
      <span class="group-title">${groupTitle}</span>
      <span class="group-chevron">${isEspCollapsed ? '▸' : '▾'}</span>
    `;
    espHeader.addEventListener('click', () => {
      if (this.collapsedGroups.has(groupTitle)) {
        this.collapsedGroups.delete(groupTitle);
      } else {
        this.collapsedGroups.add(groupTitle);
      }
      this.renderGroups(categories, hasEsp);
    });
    espGroup.appendChild(espHeader);

    const espBody = document.createElement('div');
    espBody.className = 'orbital-group-body';
    if (isEspCollapsed) espBody.style.display = 'none';

      if (hasEsp) {
        const espRow = document.createElement('label');
        espRow.className = 'orbital-item-row';

        const espChk = document.createElement('input');
        espChk.type = 'checkbox';
        espChk.className = 'orb-chk';
        espChk.checked = this.isEspActive;
        espChk.addEventListener('change', (e) => {
          this.isEspActive = e.target.checked;
          this.onEspToggle(this.isEspActive);
        });

        const espDot = document.createElement('span');
        espDot.className = 'orb-dot';
        espDot.style.backgroundColor = '#06b6d4';

        const espText = document.createElement('span');
        espText.className = 'orb-label-text';
        espText.textContent = '⚡ ESP Surface';

        espRow.appendChild(espChk);
        espRow.appendChild(espDot);
        espRow.appendChild(espText);
        espBody.appendChild(espRow);
      }

      // Atom Charges checkbox
      const chargesRow = document.createElement('label');
      chargesRow.className = 'orbital-item-row';

      const chargesChk = document.createElement('input');
      chargesChk.type = 'checkbox';
      chargesChk.className = 'orb-chk';
      chargesChk.checked = this.isChargesActive;
      chargesChk.addEventListener('change', (e) => {
        this.isChargesActive = e.target.checked;
        this.onChargesToggle(this.isChargesActive);
      });

      const chargesDot = document.createElement('span');
      chargesDot.className = 'orb-dot';
      chargesDot.style.backgroundColor = '#a855f7';

      const chargesText = document.createElement('span');
      chargesText.className = 'orb-label-text';
      chargesText.textContent = '🏷 Atom Charges';

      chargesRow.appendChild(chargesChk);
      chargesRow.appendChild(chargesDot);
      chargesRow.appendChild(chargesText);
      espBody.appendChild(chargesRow);

      espGroup.appendChild(espBody);
      this.container.appendChild(espGroup);
  }

  /**
   * Helper to format reaction pathway tracks into Molecular Explorer style categories
   */
  static parsePathwayTracks(tracks = []) {
    const items = [];
    tracks.forEach(trk => {
      const isEsp = (trk.role === 'esp' || trk.orbital_key === 'track_esp_charge' || trk.track_id === 'trk_esp');
      if (isEsp) {
        items.push({
          key: 'track_esp_charge',
          track_id: trk.track_id || 'trk_esp',
          label: 'ESP & Charges',
          color: '#c084fc',
          isEsp: true,
          badge: 'ESP'
        });
        return;
      }
      const key = trk.track_id || trk.orbital_key || trk.id;
      let title = trk.label || trk.title || key;
      title = title.replace(/Non-Bonding\s*/gi, '').replace(/\[Antibond\]/gi, '').trim();
      const sKey = `${trk.track_id || ''} ${trk.orbital_key || ''} ${trk.id || ''} ${title} ${trk.type || ''} ${trk.orbital_class || ''} ${trk.category || ''}`.toLowerCase();
      const isAnti = (trk.orbital_class === 'antibonding' || trk.category === 'antibonding' || (trk.category && trk.category.includes('*'))) ||
                     sKey.includes('star') || sKey.includes('*') || sKey.includes('lumo') || sKey.includes('antibond') || sKey.includes('virtual');
      const color = isAnti ? '#fb923c' : '#38bdf8';
      items.push({
        key,
        label: title,
        color,
        isAntiBonding: isAnti,
        energy_ev: trk.energy_ev
      });
    });

    if (items.length === 0) return [];
    return [{ name: 'Reaction Orbitals & Surfaces', items }];
  }

  /**
   * Helper to format raw orbital keys (e.g. from bundle.orbitalMeshes) into Molecular Explorer categories
   */
  static parseRawOrbitalKeys(keys = []) {
    const uniqueBases = new Set();
    keys.forEach(k => {
      const parts = k.split('_');
      if (parts.length > 1 && !isNaN(parts[0])) {
        uniqueBases.add(parts.slice(1).join('_'));
      } else {
        uniqueBases.add(k);
      }
    });

    const categoriesMap = {
      'Frontier': [],
      'σ Bonds': [],
      'π Bonds': [],
      'σ* Anti-Bonds': [],
      'π* Anti-Bonds': [],
      'Lone Pairs': [],
      'Localized Orbitals': []
    };

    const sorted = Array.from(uniqueBases).sort();
    sorted.forEach(key => {
      const upper = key.toUpperCase();
      let label = key.replace(/_/g, ' ');
      let cat = 'Localized Orbitals';
      let color = '#38bdf8';

      if (upper.includes('HOMO')) {
        cat = 'Frontier';
        color = '#38bdf8';
        label = key;
      } else if (upper.includes('LUMO')) {
        cat = 'Frontier';
        color = '#fb923c';
        label = key;
      } else if (upper.includes('SIGMASTAR') || upper.includes('SIGMA*') || key.includes('σ*')) {
        cat = 'σ* Anti-Bonds';
        color = '#fb923c';
        label = key.replace(/sigmastar[_\s]*/i, 'σ*(').replace(/$/, ')');
      } else if (upper.includes('PISTAR') || upper.includes('PI*') || key.includes('π*')) {
        cat = 'π* Anti-Bonds';
        color = '#fb923c';
        label = key.replace(/pistar[_\s]*/i, 'π*(').replace(/$/, ')');
      } else if (upper.includes('PI') || key.includes('π')) {
        cat = 'π Bonds';
        color = '#38bdf8';
        label = key.replace(/pi[_\s]*/i, 'π(').replace(/$/, ')');
      } else if (upper.includes('SIGMA') || key.includes('σ')) {
        cat = 'σ Bonds';
        color = '#38bdf8';
        label = key.replace(/sigma[_\s]*/i, 'σ(').replace(/$/, ')');
      } else if (upper.includes('LP') || upper.includes('LONE')) {
        cat = 'Lone Pairs';
        color = '#38bdf8';
        label = key.replace(/lp[_\s]*/i, 'LP(').replace(/$/, ')');
      } else if (upper.startsWith('IBO')) {
        label = `IBO ${key.replace('IBO_', '')}`;
        color = '#38bdf8';
      }

      categoriesMap[cat].push({ key, label, color });
    });

    return Object.entries(categoriesMap)
      .filter(([_, items]) => items.length > 0)
      .map(([name, items]) => ({ name, items }));
  }
}

window.ViewerOrbitalsPanel = ViewerOrbitalsPanel;
