import os
import json
import shlex

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

    mol_data = body.get('mol', '')
    molfile_arg = ''
    if mol_data and engine == 'pyscf-flat':
        mol_path = os.path.join(job_dir, "input_2d.mol")
        with open(mol_path, 'w') as f:
            f.write(mol_data)
        molfile_arg = ' --molfile "$OUT/input_2d.mol"'

    levelshift_arg = ''
    if body.get('levelshift'):
        levelshift_arg = ' --levelshift'

    script_path = os.path.join(job_dir, f'compute_{name}.sh')

    manifest = {
        "name": name,
        "smiles": rxn_smiles,
        "charge": charge,
        "engine": engine,
        "method": "b3lyp",
    }
    if ket_data:
        manifest["ket_file"] = f"{name}.ket"
    manifest["molecules"] = [
        {
            "key": name,
            "role": "molecule",
            "smiles": rxn_smiles,
            "xyz": f"{name}.xyz"
        }
    ]
    manifest_json = json.dumps(manifest, indent=2)

    lines = [
        '#!/bin/bash',
        'set -e',
        'export PYTHONUNBUFFERED=1',
        f'PYTHON={shlex.quote(venv_python)}',
        f'TOOLS={shlex.quote(tools_dir)}',
        f'OUT={shlex.quote(job_dir)}',
        f'NAME={shlex.quote(name)}',
        'echo "═══ Computing: $NAME ═══"',
        '',
        '# Step 1: Prepare molecule (3D structure, SVG, MOL)',
        f'"$PYTHON" "$TOOLS/mol_prep.py" {shlex.quote(rxn_smiles)} --name "$NAME" -o "$OUT" --charge {charge} --engine {engine}{molfile_arg}{levelshift_arg}',
        'obabel "$OUT/$NAME.mol" -O "$OUT/$NAME.cdxml" 2>/dev/null || true',
        '',
        '# Step 2: Create bundle manifest',
        'cat <<\'__BUNDLE_MANIFEST_EOF__\' > "$OUT/${NAME}_bundle.json"',
        manifest_json,
        '__BUNDLE_MANIFEST_EOF__',
        '',
        '# Export data-only archive + clean intermediates',
        '"$PYTHON" "$TOOLS/bundle_exporter.py" "$OUT" --clean',
        '',
        'echo "DONE"'
    ]

    with open(script_path, 'w') as f:
        f.write('\n'.join(lines) + '\n')
    os.chmod(script_path, 0o755)

    return script_path, name, job_dir
