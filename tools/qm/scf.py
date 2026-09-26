from pyscf import gto, scf
import sys
import os
from .geometry import parse_xyz

def run_scf(xyz_file, charge=0, spin=0, basis='6-31g*', method='b3lyp',
            chkfile=None, read_chk=False):
    """Run DFT/HF calculation, return (mol, mf, atom_labels)."""
    from pyscf import gto, scf, dft

    atom_str, labels = parse_xyz(xyz_file)
    n_atoms = len(labels)
    print(f"  SCF: {n_atoms} atoms, basis={basis}")

    # Auto-detect odd electron count before building molecule
    _Z = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
          'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
          'K':19,'Ca':20,'Br':35,'I':53}
    total_z = sum(_Z.get(s, 0) for s in labels)
    n_elec_check = total_z - charge
    if n_elec_check % 2 != 0 and spin == 0:
        spin = 1  # doublet for odd-electron systems
        print(f"  ⚠ Odd electron count ({n_elec_check}), auto-setting spin=1 (doublet)")

    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=0)
    n_elec = mol.nelectron
    print(f"  {n_elec} electrons (charge={charge}, spin={spin})")

    if spin > 0:
        # Open-shell: use unrestricted methods
        if method.lower() in ('hf', 'rhf', 'uhf'):
            mf = scf.UHF(mol)
        else:
            mf = dft.UKS(mol)
            mf.xc = method
        print(f"  Using unrestricted (open-shell) method")
    else:
        if method.lower() in ('hf', 'rhf'):
            mf = scf.RHF(mol)
        else:
            mf = dft.RKS(mol)
            mf.xc = method

    if chkfile:
        mf.chkfile = chkfile
        if read_chk and os.path.exists(chkfile):
            mf.init_guess = 'chkfile'

    mf.verbose = 0
    print(f"  Running {method.upper()} calculation...")
    energy = mf.kernel()

    if not mf.converged:
        print(f"  WARNING: SCF did not converge!", file=sys.stderr)

    print(f"  Total energy: {energy:.6f} Hartree")
    return mol, mf, labels


