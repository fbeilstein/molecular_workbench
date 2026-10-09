
export const FieldsUI = {
    render: (schemaContainer, ui, state, step_id) => {
        ui.fields.forEach(field => {
            const group = document.createElement('div');
            group.className = 'form-group';
            
            const label = document.createElement('label');
            label.innerText = field.label;
            
            const input = document.createElement('input');
            input.type = field.type;
            input.id = `field_${field.name}`;
            
            const prevData = state[step_id];
            if (prevData && prevData[field.name]) {
                input.value = prevData[field.name];
            }
            
            group.appendChild(label);
            group.appendChild(input);
            schemaContainer.appendChild(group);
        });
    },
    
    getData: async (ui) => {
        const userInput = {};
        ui.fields.forEach(field => {
            const input = document.getElementById(`field_${field.name}`);
            if (input) userInput[field.name] = input.value;
        });
        return userInput;
    }
};
