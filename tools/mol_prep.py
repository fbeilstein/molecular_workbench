#!/usr/bin/env python3
"""
mol_prep — SMILES to optimized 3D geometry

Converts a SMILES string into an optimized 3D structure with multiple
output formats suitable for slides and further computation.

Outputs:
    {name}.xyz   — XYZ coordinates (for PySCF, ORCA, viewers)
    {name}.sdf   — SDF with 3D coords (for OpenBabel/ChemDraw)
    {name}.mol   — MDL Molfile (for ChemDraw)
    {name}.svg   — 2D depiction (for slides)

Usage:
    python mol_prep.py "CCO" --name ethanol -o output/
    python mol_prep.py "CC(=O)Oc1ccccc1C(=O)O" --name aspirin
"""

import argparse
import json
import os
import subprocess
import sys

# ── RDKit imports ────────────────────────────────────────────────────────────

def _get_rdkit():
    """Import RDKit lazily so CLI --help works without it."""
    from rdkit import Chem
    from rdkit.Chem import AllChem, Draw, rdmolfiles
    return Chem, AllChem, Draw, rdmolfiles


# ── xTB path discovery ──────────────────────────────────────────────────────

def _find_xtb():
    """Find xTB binary."""
    # Primary: xtb bundled alongside this module in workbench/tools/xtb/
    _here = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(_here, 'xtb', 'bin', 'xtb'),
        'xtb',  # system PATH
    ]
    for c in candidates:
        c = os.path.abspath(c) if c != 'xtb' else c
        if c != 'xtb' and os.path.isfile(c) and os.access(c, os.X_OK):
            return c
        elif c == 'xtb':
            import shutil
            if shutil.which('xtb'):
                return 'xtb'
    return None


# ── Core pipeline ────────────────────────────────────────────────────────────

def smiles_to_3d(smiles, name, out_dir, charge=0, engine="xtb", molfile=None, levelshift=False):
    """
    Full pipeline: SMILES → RDKit 3D → geometry optimization → output files.

    Returns dict with paths to generated files and metadata.
    """
    Chem, AllChem, Draw, rdmolfiles = _get_rdkit()

    os.makedirs(out_dir, exist_ok=True)

    # 1. Parse SMILES
    params = Chem.SmilesParserParams()
    params.removeHs = False
    mol = Chem.MolFromSmiles(smiles, params)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # Store canonical SMILES (with explicit Hs if mapped)
    canon = Chem.MolToSmiles(mol)

    # 2. Add hydrogens and generate coordinates
    if engine == "pyscf-flat" and molfile is not None:
        import numpy as np
        with open(molfile, 'r') as f:
            mol_2d = Chem.MolFromMolBlock(f.read(), sanitize=False)
        mol_2d.UpdatePropertyCache(strict=False)
        mol_h = Chem.AddHs(mol_2d, addCoords=True)
        
        conf = mol_h.GetConformer()
        pos = np.array(conf.GetPositions())
        bond_lens = []
        for bond in mol_2d.GetBonds():
            a, b = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
            bond_lens.append(np.linalg.norm(pos[a] - pos[b]))
        
        if bond_lens:
            scale = 1.40 / np.mean(bond_lens)
            pos *= scale
            
        pos[:, 2] = 0.0
        
        for i in range(mol_h.GetNumAtoms()):
            conf.SetAtomPosition(i, pos[i])
            
        print(f"  Generated flat 2D starting geometry from molfile for {name}")
    else:
        mol_h = Chem.AddHs(mol)
        
        result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3())
        if result == -1:
            # Retry with random coords
            result = AllChem.EmbedMolecule(mol_h, AllChem.ETKDGv3(),
                                           useRandomCoords=True)
            if result == -1:
                raise RuntimeError(f"Failed to embed 3D coordinates for {smiles}")

        # 3. MMFF94 force-field optimization
        try:
            AllChem.MMFFOptimizeMolecule(mol_h, maxIters=500)
            print(f"  Optimized {name} with MMFF94")
        except Exception:
            print(f"  Warning: MMFF94 optimization failed for {name}, using raw embed")
            
        if engine == "pyscf-flat":
            import numpy as np
            conf = mol_h.GetConformer()
            coords = np.array(conf.GetPositions())
            coords -= coords.mean(axis=0)
            
            # SVD to find best-fit plane
            U, S, Vt = np.linalg.svd(coords)
            normal = Vt[2, :]
            
            # Rotate so normal is aligned with Z-axis
            z_axis = np.array([0, 0, 1])
            v = np.cross(normal, z_axis)
            s = np.linalg.norm(v)
            c = np.dot(normal, z_axis)
            if s > 1e-6:
                vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
                R = np.eye(3) + vx + (vx @ vx) * ((1 - c) / (s ** 2))
            else:
                R = np.eye(3) if c > 0 else -np.eye(3)
                
            flat_coords = coords @ R.T
            flat_coords[:, 2] = 0.0  # Force perfect flatness
            
            for i in range(mol_h.GetNumAtoms()):
                conf.SetAtomPosition(i, flat_coords[i])
                
            print(f"  Flattened geometry to perfect Z=0 plane for {name}")

    # 4. Compute formula
    formula = Chem.rdMolDescriptors.CalcMolFormula(mol_h)
    n_atoms = mol_h.GetNumAtoms()
    n_heavy = mol_h.GetNumHeavyAtoms()

    # 5. Write SDF (has 3D coords, good for OpenBabel)
    sdf_path = os.path.join(out_dir, f"{name}.sdf")
    writer = rdmolfiles.SDWriter(sdf_path)
    writer.write(mol_h)
    writer.close()

    # 6. Write MOL file (for ChemDraw conversion)
    mol_path = os.path.join(out_dir, f"{name}.mol")
    mol_block = Chem.MolToMolBlock(mol_h)
    with open(mol_path, 'w') as f:
        f.write(mol_block)

    # 7. Write initial XYZ
    xyz_path = os.path.join(out_dir, f"{name}.xyz")
    _write_xyz(mol_h, xyz_path)

    # 8. Geometry optimization
    opt = None
    if engine == "xtb":
        from qm.optimizers.xtb_opt import XtbOptimizer
        opt = XtbOptimizer()
    elif engine == "pyscf" or engine == "pyscf-flat":
        from qm.optimizers.pyscf_opt import PyscfOptimizer
        opt = PyscfOptimizer(levelshift=levelshift)
    
    if opt:
        xyz_path = opt.optimize(xyz_path, charge)

    # 9. Generate 2D SVG depiction
    svg_path = os.path.join(out_dir, f"{name}.svg")
    _write_svg(smiles, svg_path)
    print(f"  Generated {name}.svg")

    # 10. Save SMILES file
    smi_path = os.path.join(out_dir, f"{name}.smi")
    with open(smi_path, 'w') as f:
        f.write(canon)

    result = {
        'name': name,
        'smiles': canon,
        'formula': formula,
        'n_atoms': n_atoms,
        'n_heavy': n_heavy,
        'charge': charge,
        'files': {
            'xyz': f"{name}.xyz",
            'sdf': f"{name}.sdf",
            'mol': f"{name}.mol",
            'svg': f"{name}.svg",
            'smi': f"{name}.smi",
        }
    }

    print(f"  {name}: {formula}, {n_atoms} atoms ({n_heavy} heavy)")
    print(f"  Written: {xyz_path}")

    return result


# ── Helper functions ─────────────────────────────────────────────────────────

def _write_xyz(mol, filepath):
    """Write RDKit molecule to XYZ format."""
    Chem, _, _, _ = _get_rdkit()
    conf = mol.GetConformer()
    n = mol.GetNumAtoms()
    lines = [str(n), f"Generated by mol_prep"]
    for i in range(n):
        atom = mol.GetAtomWithIdx(i)
        pos = conf.GetAtomPosition(i)
        lines.append(f"{atom.GetSymbol():2s}  {pos.x:12.6f}  {pos.y:12.6f}  {pos.z:12.6f}")
    with open(filepath, 'w') as f:
        f.write('\n'.join(lines) + '\n')


def _write_svg(smiles, filepath, width=400, height=300):
    """Generate a 2D SVG depiction without background."""
    Chem, _, Draw, _ = _get_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return
    try:
        drawer = Draw.MolDraw2DSVG(width, height)
        opts = drawer.drawOptions()
        opts.addAtomIndices = False
        opts.clearBackground = False  # No white rectangle
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        # Also strip any remaining rect fills just in case
        import re
        svg = re.sub(r"<rect[^>]*style='[^']*fill:[^']*white[^']*'[^/]*/?>", '', svg)
        svg = re.sub(r"<rect[^>]*fill='[^']*white[^']*'[^/]*/?>", '', svg)
        with open(filepath, 'w') as f:
            f.write(svg)
    except Exception as e:
        print(f"  Warning: SVG generation failed: {e}")





# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Convert SMILES to optimized 3D geometry with multiple output formats',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s "CCO" --name ethanol -o output/
  %(prog)s "CC(=O)Oc1ccccc1C(=O)O" --name aspirin
  %(prog)s "[Cl-]" --name chloride --charge -1
        """
    )
    parser.add_argument('smiles', help='SMILES string')
    parser.add_argument('--name', '-n', required=True, help='Base name for output files')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory')
    parser.add_argument('--charge', type=int, default=0, help='Molecular charge')
    parser.add_argument('--engine', choices=['xtb', 'pyscf', 'pyscf-flat', 'none'], default='xtb', help='Geometry optimization engine (default: xtb)')
    parser.add_argument('--molfile', default=None, help='Input 2D molfile for flat starting geometry')
    parser.add_argument('--levelshift', action='store_true', help='Use SCF level shift for difficult convergence in PySCF')
    parser.add_argument('--json', action='store_true', help='Print result as JSON')

    args = parser.parse_args()

    result = smiles_to_3d(args.smiles, args.name, args.output_dir,
                          charge=args.charge, engine=args.engine, molfile=args.molfile, levelshift=args.levelshift)

    if args.json:
        print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
