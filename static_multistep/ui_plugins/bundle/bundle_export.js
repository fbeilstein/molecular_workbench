/**
 * Export utilities for multi-model XYZ and Maximal-Compression ZIP Bundle Repacker.
 */

export function exportMultiModelXyz(megajson) {
    if (!megajson || !megajson.frames) return;
    const symbols = megajson.symbols || [];
    const nAtoms = symbols.length;
    let xyzContent = '';

    megajson.frames.forEach((frame, idx) => {
        const e = frame.rel_energy_kcal !== undefined ? frame.rel_energy_kcal.toFixed(2) : '0.00';
        xyzContent += `${nAtoms}\nFrame ${idx} | E = ${e} kcal/mol\n`;
        frame.coords.forEach((coord, atomIdx) => {
            const sym = symbols[atomIdx] || 'X';
            xyzContent += `${sym.padEnd(2)} ${coord[0].toFixed(5).padStart(12)} ${coord[1].toFixed(5).padStart(12)} ${coord[2].toFixed(5).padStart(12)}\n`;
        });
    });

    const blob = new Blob([xyzContent], { type: 'chemical/x-xyz' });
    const url = URL.createObjectURL(blob);
    const a = document.createElement('a');
    a.href = url;
    a.download = `reaction_pathway_${Date.now()}.xyz`;
    document.body.appendChild(a);
    a.click();
    document.body.removeChild(a);
    URL.revokeObjectURL(url);
}

export async function exportMaxCompBundleZip(stepId = 'step_5_bundle_viewer') {
    // Show loading indicator
    let toast = document.createElement('div');
    toast.id = 'toast-repack-progress';
    toast.className = 'custom-toast';
    toast.innerHTML = `
        <div style="display: flex; align-items: center; gap: 0.6rem;">
            <div class="spinner-sm"></div>
            <span>Repacking project into maximal-compression archive (Level 9 Deflate)...</span>
        </div>
    `;
    document.body.appendChild(toast);

    try {
        const res = await fetch(`/api/step/${stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'repack_bundle' })
        });
        const data = await res.json();
        toast.remove();

        if (data.success && data.download_url) {
            const sizeMb = (data.size_bytes / (1024 * 1024)).toFixed(2);
            const origMb = (data.uncompressed_bytes / (1024 * 1024)).toFixed(2);

            // Trigger browser download
            const a = document.createElement('a');
            a.href = data.download_url;
            a.download = data.filename || 'multistep_reaction_bundle_maxcomp.zip';
            document.body.appendChild(a);
            a.click();
            document.body.removeChild(a);

            // Display success modal
            showRepackSuccessModal(data.filename, sizeMb, origMb, data.compression_ratio_pct);
        } else {
            alert('Repack error: ' + (data.error || 'Unknown error occurred.'));
        }
    } catch (err) {
        toast.remove();
        alert('Repack failed: ' + err.message);
    }
}

function showRepackSuccessModal(filename, sizeMb, origMb, ratioPct) {
    let existing = document.getElementById('modal-repack-success');
    if (existing) existing.remove();

    const modal = document.createElement('div');
    modal.id = 'modal-repack-success';
    modal.className = 'custom-modal-backdrop';
    modal.innerHTML = `
        <div class="custom-modal-dialog" style="max-width: 480px;">
            <div class="modal-header">
                <div>
                    <h3 style="font-size: 1.05rem; font-weight: 700; color: #10b981; margin: 0; display: flex; align-items: center; gap: 0.5rem;">
                        <span>✓ Maximal-Compression Bundle Created</span>
                    </h3>
                    <p style="font-size: 0.78rem; color: #94a3b8; margin-top: 2px;">All project data, high-res orbitals, SVG mechanism, and standalone viewer repacked.</p>
                </div>
                <button class="btn-modal-close" onclick="document.getElementById('modal-repack-success').remove()">&times;</button>
            </div>
            <div class="modal-body" style="padding: 1.25rem; display: flex; flex-direction: column; gap: 0.75rem;">
                <div style="background: rgba(16, 185, 129, 0.08); border: 1px solid rgba(16, 185, 129, 0.25); border-radius: 6px; padding: 0.8rem;">
                    <div style="display: flex; justify-content: space-between; font-size: 0.82rem; margin-bottom: 0.4rem;">
                        <span style="color: #cbd5e1;">Archive Name:</span>
                        <span style="color: #f8fafc; font-weight: 600; font-family: monospace;">${filename}</span>
                    </div>
                    <div style="display: flex; justify-content: space-between; font-size: 0.82rem; margin-bottom: 0.4rem;">
                        <span style="color: #cbd5e1;">Compressed Size:</span>
                        <span style="color: #38bdf8; font-weight: 700;">${sizeMb} MB</span>
                    </div>
                    <div style="display: flex; justify-content: space-between; font-size: 0.82rem; margin-bottom: 0.4rem;">
                        <span style="color: #cbd5e1;">Uncompressed Volume:</span>
                        <span style="color: #94a3b8;">${origMb} MB</span>
                    </div>
                    <div style="display: flex; justify-content: space-between; font-size: 0.82rem;">
                        <span style="color: #cbd5e1;">Space Reduction:</span>
                        <span style="color: #10b981; font-weight: 700;">${ratioPct}% compressed</span>
                    </div>
                </div>

                <div style="font-size: 0.75rem; color: #94a3b8; line-height: 1.4;">
                    📦 Includes:
                    <ul style="margin-left: 1.2rem; margin-top: 0.3rem;">
                        <li>Standalone zero-dependency HTML viewer (<code>index.html</code>)</li>
                        <li>Lightweight 3D orbital surface meshes (<code>orbitals/*_pos.json</code>, <code>*_neg.json</code>)</li>
                        <li>Reaction mechanism SVG & dynamic wiring (<code>mechanism.svg</code>, <code>svg_wiring.json</code>)</li>
                        <li>Trajectory JSON & Multi-model XYZ</li>
                    </ul>
                </div>
            </div>
            <div class="modal-footer" style="padding: 0.75rem 1.25rem; display: flex; justify-content: flex-end;">
                <button class="btn btn-primary btn-sm" onclick="document.getElementById('modal-repack-success').remove()">Done</button>
            </div>
        </div>
    `;
    document.body.appendChild(modal);
}

export async function purgeTemporaryCubes(stepId = 'step_5_bundle_viewer') {
    const btn = document.getElementById('btn-purge-cubes');
    if (btn) btn.disabled = true;
    try {
        const res = await fetch(`/api/step/${stepId}/custom`, {
            method: 'POST',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ action_type: 'purge_cubes' })
        });
        const data = await res.json();
        if (data.success) {
            if (btn) {
                btn.textContent = `✓ Purged (${data.purged})`;
                btn.style.color = '#10b981';
                btn.style.borderColor = 'rgba(16,185,129,0.4)';
            }
            alert(`Successfully purged ${data.purged} temporary cube files. Disk space freed!`);
        } else {
            alert('Purge error: ' + (data.error || 'Unknown error'));
            if (btn) btn.disabled = false;
        }
    } catch (e) {
        alert('Purge failed: ' + e.message);
        if (btn) btn.disabled = false;
    }
}

