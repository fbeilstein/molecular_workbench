
export const KetcherUI = {
    render: (schemaContainer, ui, state, step_id) => {
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

        // Give Ketcher a moment to load, then push any previous data
        const prevData = state[step_id];
        if (prevData && prevData.raw_rxnfile) {
            iframe.onload = () => {
                setTimeout(() => {
                    if(iframe.contentWindow.ketcher) {
                        iframe.contentWindow.ketcher.setMolecule(prevData.raw_rxnfile);
                    }
                }, 500);
            };
        }
    },
    
    getData: async () => {
        const iframe = document.getElementById('ketcher-frame');
        if (iframe && iframe.contentWindow.ketcher) {
            const smiles = await iframe.contentWindow.ketcher.getSmiles();
            const rxnfile = await iframe.contentWindow.ketcher.getRxn();
            return { smiles, rxnfile };
        }
        throw new Error("Ketcher frame not found or not loaded.");
    }
};
