import os
import sys
import json
import argparse
from qm.geometry import parse_xyz
from qm.scf import run_scf
from qm.localization import compute_canonical, compute_localized
from qm.esp import compute_esp_surface

def compute_all_orbitals(xyz_file, name=None, out_dir=None, charge=0, spin=0,
                         basis='6-31g*', method='b3lyp', grid_points=50,
                         chkfile=None, read_chk=False):
    """Full pipeline: SCF → canonical + localized orbitals → cube files + JSON."""
    if name is None:
        name = os.path.splitext(os.path.basename(xyz_file))[0]
    if out_dir is None:
        out_dir = os.path.dirname(xyz_file) or '.'
    os.makedirs(out_dir, exist_ok=True)

    # Run SCF
    mol, mf, labels = run_scf(xyz_file, charge, spin, basis, method,
                               chkfile, read_chk)

    # Canonical orbitals
    canonical = compute_canonical(mol, mf, labels, name, out_dir, grid_points)

    # Localized orbitals
    localized = compute_localized(mol, mf, labels, name, out_dir, grid_points)

    # ESP surface
    esp_data = compute_esp_surface(mol, mf, name, out_dir, grid_points)

    # Write manifest
    manifest = {
        'name': name,
        'method': method,
        'basis': basis,
        'charge': charge,
        'spin': spin,
        'n_electrons': mol.nelectron,
        'energy_hartree': float(mf.e_tot),
        'canonical': canonical,
        'localized': localized,
    }
    if esp_data:
        manifest['esp_surface'] = esp_data

    manifest_path = os.path.join(out_dir, f"{name}_orbitals.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"  Manifest: {manifest_path}")

    return manifest


def main():
    parser = argparse.ArgumentParser(
        description='Compute molecular orbitals (HOMO/LUMO/σ/π/LP) and generate cube files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s water.xyz --name water
  %(prog)s aspirin.xyz --charge 0 --basis 6-31g* -o output/
  %(prog)s chloride.xyz --charge -1 --name Cl
        """
    )
    parser.add_argument('xyz_file', help='Input .xyz file')
    parser.add_argument('--name', '-n', default=None, help='Base name (default: from filename)')
    parser.add_argument('-o', '--output-dir', default=None, help='Output directory')
    parser.add_argument('--charge', type=int, default=0, help='Molecular charge')
    parser.add_argument('--spin', type=int, default=0, help='Spin (2S, e.g. 0=singlet)')
    parser.add_argument('--basis', default='6-31g*', help='Basis set')
    parser.add_argument('--method', default='b3lyp', help='DFT functional or hf')
    parser.add_argument('--grid', type=int, default=50, help='Cube grid points per dimension')
    parser.add_argument('--chkfile', default=None, help='PySCF checkpoint file')
    parser.add_argument('--read-chk', action='store_true', help='Read guess from chkfile')

    args = parser.parse_args()
    if not os.path.exists(args.xyz_file):
        print(f"ERROR: {args.xyz_file} not found", file=sys.stderr)
        sys.exit(1)

    compute_all_orbitals(args.xyz_file, args.name, args.output_dir,
                         args.charge, args.spin, args.basis, args.method,
                         args.grid, args.chkfile, args.read_chk)



if __name__ == '__main__':
    main()
