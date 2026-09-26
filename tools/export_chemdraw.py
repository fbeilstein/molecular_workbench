#!/usr/bin/env python3
"""
export_chemdraw — Convert molecule files to ChemDraw-editable + SVG formats

Produces:
    {name}.cdxml — ChemDraw XML format (editable in ChemDraw)
    {name}.svg   — SVG depiction (for slides)

Uses OpenBabel for .mol/.sdf → .cdxml conversion.

Usage:
    python export_chemdraw.py aspirin.mol -o output/
    python export_chemdraw.py aspirin.sdf --name aspirin_pretty
"""

import argparse
import os
import subprocess
import sys


def export(mol_file, name=None, out_dir=None):
    """Convert .mol/.sdf to .cdxml + .svg. Returns dict of output paths."""
    if name is None:
        name = os.path.splitext(os.path.basename(mol_file))[0]
    if out_dir is None:
        out_dir = os.path.dirname(mol_file) or '.'
    os.makedirs(out_dir, exist_ok=True)

    results = {'name': name, 'files': {}}

    # 1. Convert to CDXML via OpenBabel
    cdxml_path = os.path.join(out_dir, f"{name}.cdxml")
    try:
        cmd = ['obabel', mol_file, '-O', cdxml_path]
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
        if os.path.exists(cdxml_path) and os.path.getsize(cdxml_path) > 0:
            results['files']['cdxml'] = f"{name}.cdxml"
            print(f"  ✓ {name}.cdxml (ChemDraw XML)")
        else:
            print(f"  Warning: OpenBabel produced empty .cdxml")
            if res.stderr:
                print(f"    {res.stderr.strip()}")
    except FileNotFoundError:
        print(f"  Warning: OpenBabel (obabel) not found, skipping .cdxml")
    except Exception as e:
        print(f"  Warning: .cdxml conversion failed: {e}")

    # 2. Generate SVG via RDKit (nicer than OpenBabel SVG)
    svg_path = os.path.join(out_dir, f"{name}.svg")
    try:
        from rdkit import Chem
        from rdkit.Chem import Draw

        # Try reading from .mol file
        mol = Chem.MolFromMolFile(mol_file, removeHs=True)
        if mol is None:
            # Try as SDF
            suppl = Chem.SDMolSupplier(mol_file, removeHs=True)
            for m in suppl:
                if m is not None:
                    mol = m
                    break

        if mol is not None:
            drawer = Draw.MolDraw2DSVG(500, 400)
            opts = drawer.drawOptions()
            opts.addAtomIndices = False
            opts.clearBackground = False
            drawer.DrawMolecule(mol)
            drawer.FinishDrawing()
            with open(svg_path, 'w') as f:
                f.write(drawer.GetDrawingText())
            results['files']['svg'] = f"{name}.svg"
            print(f"  ✓ {name}.svg")
        else:
            print(f"  Warning: Could not parse {mol_file} for SVG")
    except ImportError:
        # Fallback to OpenBabel SVG
        try:
            cmd = ['obabel', mol_file, '-O', svg_path]
            subprocess.run(cmd, capture_output=True, timeout=30)
            if os.path.exists(svg_path):
                results['files']['svg'] = f"{name}.svg"
                print(f"  ✓ {name}.svg (via OpenBabel)")
        except Exception:
            print(f"  Warning: SVG generation failed")

    return results


def main():
    parser = argparse.ArgumentParser(
        description='Export molecule to ChemDraw (.cdxml) and SVG formats',
        epilog="""
Examples:
  %(prog)s aspirin.mol -o output/
  %(prog)s reaction_ts.sdf --name transition_state
        """
    )
    parser.add_argument('mol_file', help='Input .mol or .sdf file')
    parser.add_argument('--name', '-n', default=None, help='Base name')
    parser.add_argument('-o', '--output-dir', default=None, help='Output directory')

    args = parser.parse_args()
    if not os.path.exists(args.mol_file):
        print(f"ERROR: {args.mol_file} not found", file=sys.stderr)
        sys.exit(1)

    export(args.mol_file, args.name, args.output_dir)


if __name__ == '__main__':
    main()
