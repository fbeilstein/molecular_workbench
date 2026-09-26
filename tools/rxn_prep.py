#!/usr/bin/env python3
"""
rxn_prep — Prepare reaction geometries with atom mapping
"""
import argparse
import json
import os
import sys

from rxn_mapping import parse_reaction_smiles, detect_charge, auto_map_atoms, _validate_pt_mapping
from rxn_mechanism import detect_driving_coordinates
from rxn_complex import _build_neb_complexes, _read_xyz

def prepare_reaction(rxn_smiles, name='reaction', out_dir='.', auto_confirm=False):
    """
    Full reaction preparation pipeline.
    Returns dict with reactant/product info, mapping, and per-fragment charges.
    """
    os.makedirs(out_dir, exist_ok=True)

    # 1. Parse
    reactants, products = parse_reaction_smiles(rxn_smiles)
    print(f"Reaction: {' + '.join(reactants)} → {' + '.join(products)}")

    # 2. Detect per-fragment charges
    r_charges = [detect_charge(s) for s in reactants]
    p_charges = [detect_charge(s) for s in products]
    total_charge = sum(r_charges)
    print(f"Total charge: {total_charge}")
    for i, (s, c) in enumerate(zip(reactants, r_charges)):
        print(f"  R{i+1}: {s}  charge={c}")
    for i, (s, c) in enumerate(zip(products, p_charges)):
        print(f"  P{i+1}: {s}  charge={c}")

    # 3. Auto-map
    mapping, summary = auto_map_atoms(reactants, products)
    print(f"\nAtom mapping:\n{summary}")

    # 3b. Validate mapping for proton transfer consistency
    mapping = _validate_pt_mapping(mapping, reactants, products)

    # 4. Interactive confirmation (CLI mode)
    if not auto_confirm and sys.stdin.isatty():
        print(f"\n  Mapped {len(mapping)} atom pairs.")
        response = input("  Accept mapping? [Y/n/edit]: ").strip().lower()
        if response == 'n':
            print("  Mapping rejected. Please provide manual mapping.")
            print("  Format: R<frag>:<atom> P<frag>:<atom> (e.g. R1:0 P1:2)")
            print("  Enter mappings one per line, empty line to finish:")
            mapping = {}
            while True:
                line = input("  > ").strip()
                if not line:
                    break
                try:
                    r_part, p_part = line.split()
                    r_frag, r_atom = r_part[1:].split(':')
                    p_frag, p_atom = p_part[1:].split(':')
                    mapping[('R', int(r_frag)-1, int(r_atom))] = ('P', int(p_frag)-1, int(p_atom))
                except Exception as e:
                    print(f"    Parse error: {e}. Try again.")

    # 5. Generate 3D for each fragment
    tools_dir = os.path.dirname(os.path.abspath(__file__))
    if tools_dir not in sys.path:
        sys.path.insert(0, tools_dir)
    from mol_prep import smiles_to_3d

    print(f"\nGenerating 3D structures...")
    fragments = {}
    r_xyz_paths = []
    p_xyz_paths = []

    for side, frags, charges in [('reactant', reactants, r_charges),
                                  ('product', products, p_charges)]:
        for i, (smi, charge) in enumerate(zip(frags, charges)):
            frag_name = f"{name}_{side[0]}{i+1}"
            print(f"\n  [{frag_name}] {smi} (charge={charge})")
            try:
                result = smiles_to_3d(smi, frag_name, out_dir, charge=charge, skip_xtb=True)
                fragments[frag_name] = result
                fragments[frag_name]['charge'] = charge
                xyz_path = os.path.join(out_dir, f"{frag_name}.xyz")
                if os.path.exists(xyz_path):
                    if side == 'reactant':
                        r_xyz_paths.append(xyz_path)
                    else:
                        p_xyz_paths.append(xyz_path)
            except Exception as e:
                print(f"  ERROR generating 3D for {smi}: {e}")
                fragments[frag_name] = {'error': str(e), 'charge': charge}

    # 6. Generate matched complex .xyz files (for ORCA NEB)
    combined_r = None
    combined_p = None
    complex_driving = None
    if len(r_xyz_paths) >= 1 and len(p_xyz_paths) >= 1:
        print(f"\nGenerating matched NEB complexes...")
        combined_r, combined_p, complex_driving = _build_neb_complexes(
            r_xyz_paths, p_xyz_paths, mapping,
            reactants, products, out_dir, name)
        print(f"  ✓ Reactant complex: {combined_r}")
        print(f"  ✓ Product complex:  {combined_p}")

        # Verify atom counts match
        r_atoms = _read_xyz(combined_r)
        p_atoms = _read_xyz(combined_p)
        if len(r_atoms) == len(p_atoms):
            print(f"  ✓ Atom counts match: {len(r_atoms)} atoms")
            types_match = all(r[0] == p[0] for r, p in zip(r_atoms, p_atoms))
            if types_match:
                print(f"  ✓ Atom types match at all positions")
            else:
                print(f"  ⚠ Atom type mismatch detected!")
        else:
            print(f"  ⚠ Atom counts differ: R={len(r_atoms)}, P={len(p_atoms)}")

    # Detect driving coordinates
    ket_path = os.path.join(out_dir, f"{name}_reaction.ket")
    if not os.path.exists(ket_path):
        ket_path = os.path.join(out_dir, f"{name}.ket")
    if not os.path.exists(ket_path):
        ket_path = None
        
    driving = detect_driving_coordinates(reactants, products, mapping, ket_path=ket_path)
    print(f"\nDriving coordinates:")
    if driving['formed_labels']:
        print(f"  Bonds forming: {', '.join(driving['formed_labels'])}")
    if driving['broken_labels']:
        print(f"  Bonds breaking: {', '.join(driving['broken_labels'])}")
    if driving.get('pt_labels'):
        print(f"  Proton transfers: {', '.join(driving['pt_labels'])}")
    if not driving['formed_labels'] and not driving['broken_labels'] and not driving.get('pt_labels'):
        print(f"  (none detected — check atom mapping)")

    # 7. Save reaction metadata
    # Serialize driving coordinates for JSON
    driving_json = {
        'formed': [{'atoms': [f"R{b[0][0]+1}:{b[0][1]}", f"R{b[1][0]+1}:{b[1][1]}"],
                    'label': driving['formed_labels'][i]}
                   for i, b in enumerate(driving['formed_bonds'])],
        'broken': [{'atoms': [f"R{b[0][0]+1}:{b[0][1]}", f"R{b[1][0]+1}:{b[1][1]}"],
                    'label': driving['broken_labels'][i]}
                   for i, b in enumerate(driving['broken_bonds'])],
        'proton_transfers': [{'atoms': [f"R{b[0][0]+1}:{b[0][1]}", f"R{b[1][0]+1}:{b[1][1]}"],
                               'label': driving['pt_labels'][i]}
                             for i, b in enumerate(driving.get('proton_transfers', []))],
    }
    # Add complex-global indices if available
    if complex_driving:
        driving_json['complex_formed_1based'] = complex_driving['formed_1based']
        driving_json['complex_broken_1based'] = complex_driving['broken_1based']
        driving_json['complex_pt_1based'] = complex_driving.get('pt_1based', [])

    meta = {
        'reaction_smiles': rxn_smiles,
        'reactants': reactants,
        'products': products,
        'reactant_charges': r_charges,
        'product_charges': p_charges,
        'total_charge': total_charge,
        'mapping': {f"R{k[1]+1}:{k[2]}": f"P{v[1]+1}:{v[2]}"
                    for k, v in mapping.items()},
        'driving_coordinates': driving_json,
        'fragments': fragments,
        'combined_reactant_xyz': combined_r,
        'combined_product_xyz': combined_p,
    }

    meta_path = os.path.join(out_dir, f"{name}_reaction.json")
    with open(meta_path, 'w') as f:
        json.dump(meta, f, indent=2, default=str)
    print(f"\n✓ Reaction metadata: {meta_path}")

    return meta

def main():
    parser = argparse.ArgumentParser(
        description='Prepare reaction geometries with atom mapping',
        epilog="""
Examples:
  %(prog)s "[Cl-].CBr>>CCl.[Br-]" -o output/sn2/
  %(prog)s "CC(=O)O.OCC>>CC(=O)OCC.O" --name ester --auto
        """
    )
    parser.add_argument('reaction', help='Reaction SMILES (reactants>>products)')
    parser.add_argument('--name', '-n', default='reaction', help='Base name')
    parser.add_argument('-o', '--output-dir', default='.', help='Output directory')
    parser.add_argument('--auto', action='store_true',
                        help='Auto-confirm mapping without prompting')

    args = parser.parse_args()
    prepare_reaction(args.reaction, args.name, args.output_dir, args.auto)

if __name__ == '__main__':
    main()
