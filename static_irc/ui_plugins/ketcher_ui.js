function loadPreviousKetcherData(iframe, rawRxnfile, rawSmiles) {
    if (!rawRxnfile && !rawSmiles) return;
    let attempts = 0;
    const maxAttempts = 100; // 100 * 150ms = 15 seconds
    let applied = false;

    const tryApply = () => {
        if (applied) return true;
        const ketcher = iframe.contentWindow?.ketcher;
        if (ketcher && typeof ketcher.setMolecule === 'function') {
            try {
                if (rawRxnfile) {
                    const p = ketcher.setMolecule(rawRxnfile);
                    if (p && typeof p.catch === 'function') {
                        p.catch(() => {
                            if (rawSmiles) {
                                try { ketcher.setMolecule(rawSmiles); } catch (e) {}
                            }
                        });
                    }
                } else if (rawSmiles) {
                    ketcher.setMolecule(rawSmiles);
                }
                applied = true;
                return true;
            } catch (e) {
                if (rawSmiles) {
                    try {
                        ketcher.setMolecule(rawSmiles);
                        applied = true;
                        return true;
                    } catch (err) {}
                }
            }
        }
        return false;
    };

    const interval = setInterval(() => {
        attempts++;
        if (tryApply() || attempts >= maxAttempts) {
            clearInterval(interval);
        }
    }, 150);

    iframe.onload = () => {
        setTimeout(tryApply, 200);
    };
}

function renderKetcherUI(schemaContainer, ui, state, step_id) {
    const wrapper = document.createElement('div');
    wrapper.style.height = '500px';
    wrapper.style.border = '1px solid #ccc';
    wrapper.style.marginBottom = '1rem';

    const iframe = document.createElement('iframe');
    iframe.id = 'ketcher-frame';
    iframe.src = '/ketcher/index.html';
    iframe.style.width = '100%';
    iframe.style.height = '100%';
    iframe.style.border = 'none';
    wrapper.appendChild(iframe);
    schemaContainer.appendChild(wrapper);

    const prevData = state[step_id];
    const rawRxn = prevData?.raw_rxnfile || ui?.raw_rxnfile || '';
    const rawSmi = prevData?.raw_smiles || ui?.raw_smiles || '';
    loadPreviousKetcherData(iframe, rawRxn, rawSmi);
}

async function getKetcherData() {
    const iframe = document.getElementById('ketcher-frame');
    const ketcher = iframe?.contentWindow?.ketcher;
    if (!ketcher) throw new Error("Ketcher frame not found or not loaded.");

    const smiles = await ketcher.getSmiles();
    let rxnfile = "";
    try {
        rxnfile = await ketcher.getRxn();
    } catch (e) {
        // Multi-arrow or query reactions may not export to rxnfile, SMILES will be used
    }
    return { smiles, rxnfile };
}

export const KetcherUI = {
    render: renderKetcherUI,
    getData: getKetcherData
};
