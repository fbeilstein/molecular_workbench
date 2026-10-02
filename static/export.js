window.WB = window.WB || {};

window.WB.toggleExportSelectAll = function(checkbox) {
    const listDiv = document.getElementById('export-orbital-list');
    const chks = listDiv.querySelectorAll('input[type="checkbox"]');
    chks.forEach(chk => chk.checked = checkbox.checked);
};

window.WB.downloadCustomBundle = async function() {
    const name = document.getElementById('job-name').value;
    if (!name) {
        window.WB.status('No job loaded to download.');
        return;
    }
    
    const listDiv = document.getElementById('export-orbital-list');
    listDiv.innerHTML = '';
    
    const sidebar = document.getElementById('orbital-toggles');
    if (sidebar) {
        const inputs = sidebar.querySelectorAll('input[data-orb-file]');
        inputs.forEach(inp => {
            const fileAttr = inp.getAttribute('data-orb-file');
            const text = inp.closest('label').innerText.trim();
            
            const div = document.createElement('div');
            div.style.marginBottom = '5px';
            
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            if (fileAttr) {
                cb.value = fileAttr;
            }
            cb.checked = document.getElementById('export-select-all').checked;
            cb.id = 'export_cb_' + Math.random().toString(36).substring(7);
            
            const cl = document.createElement('label');
            cl.htmlFor = cb.id;
            cl.innerText = text;
            cl.style.marginLeft = '8px';
            cl.style.fontSize = '12px';
            
            div.appendChild(cb);
            div.appendChild(cl);
            listDiv.appendChild(div);
        });
        
        // Handle ESP separately if present
        const espCb = sidebar.querySelector('input[data-esp]');
        if (espCb) {
            const text = espCb.closest('label').innerText.trim();
            const div = document.createElement('div');
            div.style.marginBottom = '5px';
            
            const cb = document.createElement('input');
            cb.type = 'checkbox';
            cb.value = `molecules/${name}_esp.json`;
            cb.checked = document.getElementById('export-select-all').checked;
            cb.id = 'export_cb_esp';
            
            const cl = document.createElement('label');
            cl.htmlFor = cb.id;
            cl.innerText = text || "ESP Surface";
            cl.style.marginLeft = '8px';
            cl.style.fontSize = '12px';
            
            div.appendChild(cb);
            div.appendChild(cl);
            listDiv.appendChild(div);
        }
    }
    
    if (listDiv.children.length === 0) {
        listDiv.innerHTML = '<div style="color:#888; font-size:12px;">No orbitals found to export.</div>';
    }
    
    document.getElementById('export-modal').style.display = 'block';
};

window.WB.confirmDownloadBundle = async function() {
    document.getElementById('export-modal').style.display = 'none';
    const name = document.getElementById('job-name').value;
    
    const keepFiles = [];
    const chks = document.getElementById('export-orbital-list').querySelectorAll('input[type="checkbox"]:checked');
    chks.forEach(chk => {
        const fileAttr = chk.value;
        if (fileAttr) {
            if (fileAttr.endsWith('.json') && !fileAttr.endsWith('_esp.json')) {
                keepFiles.push(fileAttr.replace('.json', '_pos.json'));
                keepFiles.push(fileAttr.replace('.json', '_neg.json'));
                keepFiles.push(fileAttr);
            } else {
                keepFiles.push(fileAttr);
            }
        }
    });
    
    const payload = { name: name, keep_files: keepFiles };
    if (window._customIsovalue) {
        payload.isovalue = window._customIsovalue;
    }

    try {
        window.WB.status('Packaging bundle for download...');
        const res = await fetch('/api/bundle/package', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify(payload)
        });
        const data = await res.json();
        if (data.download_url) {
            window.location.href = data.download_url;
            window.WB.status('Download started.');
        } else {
            window.WB.status('Error packaging bundle.');
        }
    } catch(e) {
        console.error('Failed to package bundle', e);
        window.WB.status('Error packaging bundle: ' + e.message);
    }
};

window.WB.purgeCubes = async function() {
    const name = document.getElementById('job-name').value;
    if (!name) return;
    
    try {
        window.WB.status('Purging cube files...');
        const res = await fetch('/api/bundle/purge-cubes', {
            method: 'POST',
            headers: {'Content-Type': 'application/json'},
            body: JSON.stringify({ name })
        });
        const data = await res.json();
        
        const btn = document.getElementById('btn-purge-cubes');
        if (btn) {
            btn.textContent = `Purged (${data.purged})`;
            btn.disabled = true;
        }
        
        const slider = document.getElementById('isovalue-slider');
        if (slider) {
            slider.disabled = true;
            slider.title = 'Cube files purged — remeshing unavailable';
        }
        
        window.WB.status(`Purged ${data.purged} cube files.`);
    } catch(e) {
        console.error('Failed to purge cubes', e);
        window.WB.status('Error purging cubes: ' + e.message);
    }
};
