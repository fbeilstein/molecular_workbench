// Molecular Gallery Interactive Script

function inflateTable() {
  const table = document.getElementById('molTable');
  if (!table) return;
  const rows = table.querySelectorAll('tbody tr');

  rows.forEach(row => {
    const cells = row.children;
    if (cells.length < 4) return;

    // Support 4-column variant: if column 3 (Open) was omitted, insert an empty cell
    if (cells.length === 4) {
      const openTd = document.createElement('td');
      row.insertBefore(openTd, cells[2]);
    }

    const tdImg = row.children[0];
    const tdName = row.children[1];
    const tdOpen = row.children[2];
    const tdNotes = row.children[3];
    const tdSmiles = row.children[4];

    // 1. Structure Image
    let imgSrc = '';
    const imgEl = tdImg.querySelector('img');
    if (imgEl) {
      imgEl.classList.add('mol-thumb');
      imgSrc = imgEl.getAttribute('src') || '';
    } else {
      const text = tdImg.textContent.trim();
      if (text) {
        imgSrc = text.includes('/') ? text : `bundles/molecules/${text}${text.endsWith('.svg') ? '' : '.svg'}`;
        const alt = imgSrc.split('/').pop().replace(/\.svg$/i, '');
        tdImg.innerHTML = `<img src="${imgSrc}" alt="${alt}" class="mol-thumb">`;
      }
    }

    // 2. Name: parse **bold** -> <strong>bold</strong>
    tdName.innerHTML = tdName.innerHTML.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');

    // 3. Open button
    const aEl = tdOpen.querySelector('a');
    if (aEl) {
      aEl.classList.add('view-btn');
      aEl.setAttribute('target', '_blank');
      if (!aEl.textContent.trim()) aEl.textContent = 'Open';
    } else {
      let zipPath = tdOpen.textContent.trim();
      if (!zipPath && imgSrc) {
        zipPath = imgSrc.replace(/\.svg$/i, '.zip');
      } else if (zipPath && !zipPath.includes('/')) {
        zipPath = `bundles/molecules/${zipPath}${zipPath.endsWith('.zip') ? '' : '.zip'}`;
      }
      if (zipPath) {
        tdOpen.innerHTML = `<a href="viewer.html?bundle=${zipPath}&download" target="_blank" class="view-btn">Open</a>`;
      }
    }

    // 4. Notes: parse **bold** and add class notes-cell
    tdNotes.classList.add('notes-cell');
    tdNotes.innerHTML = tdNotes.innerHTML.replace(/\*\*(.+?)\*\*/g, '<strong>$1</strong>');

    // 5. SMILES: wrap plain text in <code class="smiles-copy" title="Click to copy">
    const codeEl = tdSmiles.querySelector('.smiles-copy');
    if (!codeEl) {
      const smiles = tdSmiles.textContent.trim();
      if (smiles) {
        tdSmiles.innerHTML = `<code class="smiles-copy" title="Click to copy">${smiles}</code>`;
      }
    }
  });

  // Dynamic counter in header
  const rowCount = document.getElementById('rowCount');
  if (rowCount) {
    rowCount.textContent = `${rows.length} molecules`;
  }
}

document.addEventListener('DOMContentLoaded', () => {
  // 1. Inflate bare-minimum table cells into interactive components
  inflateTable();

  const searchInput = document.getElementById('searchInput');
  const rowCount = document.getElementById('rowCount');
  const tableRows = document.querySelectorAll('#molTable tbody tr');
  const totalCount = tableRows.length;
  const popover = document.getElementById('imgPreview');
  const popoverImg = popover ? popover.querySelector('img') : null;

  // 2. Search Filter
  if (searchInput) {
    searchInput.addEventListener('input', () => {
      const q = searchInput.value.toLowerCase().trim();
      let visibleCount = 0;
      tableRows.forEach(row => {
        const text = row.textContent.toLowerCase();
        const match = !q || text.includes(q);
        row.style.display = match ? '' : 'none';
        if (match) visibleCount++;
      });
      if (rowCount) {
        rowCount.textContent = q ? `Showing ${visibleCount} of ${totalCount} molecules` : `${totalCount} molecules`;
      }
    });
  }

  // 3. Click to Copy SMILES
  document.addEventListener('click', (e) => {
    const btn = e.target.closest('.smiles-copy');
    if (!btn) return;
    const smiles = btn.getAttribute('data-smiles') || btn.textContent.trim();
    navigator.clipboard.writeText(smiles).then(() => {
      const origText = btn.textContent;
      btn.textContent = 'Copied!';
      btn.classList.add('copied');
      setTimeout(() => {
        btn.textContent = origText;
        btn.classList.remove('copied');
      }, 1200);
    }).catch(err => {
      console.warn('Clipboard copy failed:', err);
    });
  });

  // 4. Hover Formula Image Enlargement
  document.addEventListener('mouseover', (e) => {
    const thumb = e.target.closest('.mol-thumb');
    if (!thumb || !popover || !popoverImg) return;
    popoverImg.src = thumb.src;
    popover.style.display = 'block';
  });

  document.addEventListener('mousemove', (e) => {
    if (!popover || popover.style.display !== 'block') return;
    const padding = 20;
    let left = e.clientX + padding;
    let top = e.clientY - 40;

    const popoverWidth = 340;
    const popoverHeight = 260;
    if (left + popoverWidth > window.innerWidth) {
      left = e.clientX - popoverWidth - padding;
    }
    if (top + popoverHeight > window.innerHeight) {
      top = window.innerHeight - popoverHeight - 10;
    }
    if (top < 10) top = 10;

    popover.style.left = `${left}px`;
    popover.style.top = `${top}px`;
  });

  document.addEventListener('mouseout', (e) => {
    const thumb = e.target.closest('.mol-thumb');
    if (thumb && popover) {
      popover.style.display = 'none';
    }
  });
});
