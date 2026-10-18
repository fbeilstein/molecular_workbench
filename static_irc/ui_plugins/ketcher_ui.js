function loadPreviousKetcherData(iframe, rawRxnfile, rawSmiles) {
    if (!rawRxnfile && !rawSmiles) return;
    iframe.onload = () => {
        setTimeout(() => {
            const ketcher = iframe.contentWindow?.ketcher;
            if (ketcher) {
                if (rawRxnfile) ketcher.setMolecule(rawRxnfile);
                else if (rawSmiles) ketcher.setMolecule(rawSmiles);
            }
        }, 500);
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
    loadPreviousKetcherData(iframe, prevData?.raw_rxnfile, prevData?.raw_smiles);
}

async function getKetcherData() {
    const iframe = document.getElementById('ketcher-frame');
    const ketcher = iframe?.contentWindow?.ketcher;
    if (!ketcher) throw new Error("Ketcher frame not found or not loaded.");

    const smiles = await ketcher.getSmiles();
    const rxnfile = await ketcher.getRxn();
    return { smiles, rxnfile };
}

export const KetcherUI = {
    render: renderKetcherUI,
    getData: getKetcherData
};
