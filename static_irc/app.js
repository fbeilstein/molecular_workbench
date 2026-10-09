document.addEventListener('DOMContentLoaded', () => {
    let plugins = [];
    let state = {};
    let activeStepIndex = 0;

    const tabsContainer = document.getElementById('wizard-tabs');
    const schemaContainer = document.getElementById('schema-container');
    const errorContainer = document.getElementById('error-container');
    const runBtn = document.getElementById('btn-run');
    const resetBtn = document.getElementById('btn-reset');

    async function loadPipeline() {
        const [pluginsRes, stateRes] = await Promise.all([
            fetch('/api/plugins'),
            fetch('/api/state')
        ]);
        
        plugins = await pluginsRes.json();
        state = await stateRes.json();
        
        determineActiveStep();
        renderTabs();
        renderActiveStep();
    }

    function determineActiveStep() {
        activeStepIndex = 0;
        for (let i = 0; i < plugins.length; i++) {
            if (state[plugins[i].step_id]) {
                activeStepIndex = i; // Move to next or stay if it's the last
                if (i + 1 < plugins.length) {
                    activeStepIndex = i + 1;
                }
            } else {
                break;
            }
        }
    }

    function renderTabs() {
        tabsContainer.innerHTML = '';
        plugins.forEach((plugin, index) => {
            const isCompleted = !!state[plugin.step_id];
            const isLocked = !isCompleted && index > activeStepIndex;
            
            const btn = document.createElement('button');
            btn.className = `tab ${index === activeStepIndex ? 'active' : ''} ${isCompleted ? 'completed' : ''} ${isLocked ? 'locked' : ''}`;
            btn.innerText = plugin.step_name;
            btn.disabled = isLocked;
            
            btn.onclick = () => {
                if (!isLocked) {
                    activeStepIndex = index;
                    renderTabs();
                    renderActiveStep();
                }
            };
            
            tabsContainer.appendChild(btn);
        });
    }

    function renderActiveStep() {
        errorContainer.style.display = 'none';
        schemaContainer.innerHTML = '';
        
        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema;
        
        const msg = document.createElement('p');
        msg.innerText = ui.message || '';
        schemaContainer.appendChild(msg);
        
        if (ui.fields) {
            ui.fields.forEach(field => {
                const group = document.createElement('div');
                group.className = 'form-group';
                
                const label = document.createElement('label');
                label.innerText = field.label;
                
                const input = document.createElement('input');
                input.type = field.type;
                input.id = `field_${field.name}`;
                
                // pre-fill if correcting
                const prevData = state[plugin.step_id];
                if (prevData && prevData[field.name]) {
                    input.value = prevData[field.name];
                }
                
                group.appendChild(label);
                group.appendChild(input);
                schemaContainer.appendChild(group);
            });
        }
        
        // Show status of current step data if any
        if (state[plugin.step_id]) {
            const pre = document.createElement('pre');
            pre.innerText = "Current State: " + JSON.stringify(state[plugin.step_id], null, 2);
            schemaContainer.appendChild(pre);
            runBtn.innerText = "Re-Run / Correct";
        } else {
            runBtn.innerText = "Run & Continue";
        }
    }

    runBtn.onclick = async () => {
        errorContainer.style.display = 'none';
        runBtn.disabled = true;
        
        const plugin = plugins[activeStepIndex];
        const ui = plugin.ui_schema;
        const userInput = {};
        
        if (ui.fields) {
            ui.fields.forEach(field => {
                const input = document.getElementById(`field_${field.name}`);
                if (input) userInput[field.name] = input.value;
            });
        }
        
        try {
            const res = await fetch(`/api/step/${plugin.step_id}/run`, {
                method: 'POST',
                headers: { 'Content-Type': 'application/json' },
                body: JSON.stringify(userInput)
            });
            const data = await res.json();
            
            if (data.success) {
                // reload pipeline to get new UI schemas (since they depend on state)
                await loadPipeline(); 
            } else {
                errorContainer.innerText = data.error || "Unknown error";
                errorContainer.style.display = 'block';
            }
        } catch (err) {
            errorContainer.innerText = err.message;
            errorContainer.style.display = 'block';
        } finally {
            runBtn.disabled = false;
        }
    };
    
    resetBtn.onclick = async () => {
        if(confirm("Reset entire pipeline?")) {
            await fetch('/api/reset', { method: 'POST' });
            await loadPipeline();
        }
    };

    // Init
    loadPipeline();
});
