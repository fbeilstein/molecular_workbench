import numpy as np

def parse_xyz(filepath):
    """Parse XYZ file → PySCF atom string + atom labels list."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    n_atoms = int(lines[0].strip())
    atoms = []
    labels = []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if len(parts) >= 4:
            sym = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append(f'{sym} {x} {y} {z}')
            labels.append(sym)
    return '; '.join(atoms), labels


def make_atom_ids(atom_labels):
    """Create unique atom IDs: ['C','C','O','H'] → ['C1','C2','O1','H1']."""
    ids = []
    count = {}
    for sym in atom_labels:
        count[sym] = count.get(sym, 0) + 1
        ids.append(f"{sym}{count[sym]}")
    return ids


