// Molecular Gallery Interactive Script
document.addEventListener('DOMContentLoaded', () => {
  const searchInput = document.getElementById('searchInput');
  const rowCount = document.getElementById('rowCount');
  const tableRows = document.querySelectorAll('#molTable tbody tr');
  const totalCount = tableRows.length;
  const popover = document.getElementById('imgPreview');
  const popoverImg = popover ? popover.querySelector('img') : null;

  // 1. Search Filter
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

  // 2. Click to Copy SMILES
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

  // 3. Hover Formula Image Enlargement
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

