import os
import json

def detect_charge_from_smiles(smiles):
    """Detect total charge from SMILES string by calculating formal charges."""
    try:
        from rdkit import Chem
        mol = Chem.MolFromSmiles(smiles, sanitize=False)
        if mol:
            mol.UpdatePropertyCache(strict=False)
            return sum(atom.GetFormalCharge() for atom in mol.GetAtoms())
    except ImportError:
        pass
        
    import re
    total = 0
    # Match patterns like [C+], [NH4+], [Br-], [O-2], [Fe+3], [O+:6]
    for match in re.finditer(r'\[([^\]]*?)([+-])(\d*)(?::\d+)?\]', smiles):
        sign = 1 if match.group(2) == '+' else -1
        magnitude = int(match.group(3)) if match.group(3) else 1
        total += sign * magnitude
    return total

def build_compute_script(body, output_dir, venv_python, tools_dir):
    rxn_smiles = body.get('smiles', '')
    name = body.get('name', 'molecule')
    charge = body.get('charge', 0)
    engine = body.get('engine', 'xtb')

    smiles_charge = detect_charge_from_smiles(rxn_smiles)
    if charge == 0 and smiles_charge != 0:
        charge = smiles_charge

    job_dir = os.path.join(output_dir, name)
    os.makedirs(job_dir, exist_ok=True)

    ket_data = body.get('ket', '')
    if ket_data:
        ket_path = os.path.join(job_dir, f"{name}.ket")
        with open(ket_path, 'w') as f:
            f.write(ket_data)

    svg_data = body.get('svg', '')
    if svg_data:
        svg_path = os.path.join(job_dir, f"{name}.svg")
        with open(svg_path, 'w') as f:
            f.write(svg_data)

    script_path = os.path.join(job_dir, f'compute_{name}.sh')

    safe_smiles_json = json.dumps(rxn_smiles)
    
    lines = [
        '#!/bin/bash',
        'set -e',
        'export PYTHONUNBUFFERED=1',
        f'PYTHON="{venv_python}"',
        f'TOOLS="{tools_dir}"',
        f'OUT="{job_dir}"',
        f'echo "═══ Computing: {name} ═══"',
        '',
        '# Step 1: Prepare molecule (3D structure, SVG, MOL)',
        f'$PYTHON $TOOLS/mol_prep.py {safe_smiles_json} --name {name} -o $OUT --charge {charge} --engine {engine}',
        f'obabel $OUT/{name}.mol -O $OUT/{name}.cdxml 2>/dev/null || true',
        '',
        '# Step 2: Create bundle manifest',
        f'cat <<\'EOF\' > $OUT/{name}_bundle.json',
        '{',
        f'  "name": "{name}",',
        f'  "smiles": {safe_smiles_json},',
        f'  "charge": {charge},',
        f'  "engine": "{engine}",',
        f'  "method": "b3lyp",',
        '  "molecules": [',
        '    {',
        f'      "key": "{name}",',
        '      "role": "molecule",',
        f'      "smiles": {safe_smiles_json},',
        f'      "xyz": "{name}.xyz"',
        '    }',
        '  ]',
        '}',
        'EOF',
        '',
        '# Export data-only archive + clean intermediates',
        f'$PYTHON $TOOLS/bundle_exporter.py $OUT --clean',
        '',
        'echo "DONE"'
    ]

    with open(script_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    os.chmod(script_path, 0o755)

    return script_path, name, job_dir
