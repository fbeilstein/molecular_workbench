from rdkit import Chem
from rdkit.Chem import AllChem
from .common import align_to_reference, set_positions


def _refine_mmff(mol):
    AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
    return mol.GetConformer().GetPositions()


def _refine_xtb(mol, p_orig, charge, mult):
    from ase import Atoms
    from ase.optimize import LBFGS
    from tblite.ase import TBLite
    symbols = [a.GetSymbol() for a in mol.GetAtoms()]
    frag_atoms = Atoms(symbols=symbols, positions=p_orig.copy())
    frag_atoms.calc = TBLite(method="GFN2-xTB", charge=charge, multiplicity=mult, verbosity=0)
    opt = LBFGS(frag_atoms, logfile=None, maxstep=0.2)
    opt.run(fmax=0.04, steps=100)
    return frag_atoms.get_positions()


def refine_fragment_shape(mol_block: str, method="mmff", charge=0, mult=1) -> str:
    """Refine isolated molecule geometry in-place using MMFF94 or GFN2-xTB."""
    mol = Chem.MolFromMolBlock(mol_block, removeHs=False)
    if mol is None:
        raise ValueError("Invalid molblock")

    p_orig = mol.GetConformer().GetPositions().copy()
    m = method.lower()
    if m == "mmff":
        p_new = _refine_mmff(mol)
    elif m == "xtb":
        p_new = _refine_xtb(mol, p_orig, charge, mult)
    else:
        raise ValueError(f"Unknown refinement method: {method}")

    # Align back to original pose so docking orientation is preserved
    p_aligned = align_to_reference(p_orig, p_new)
    return Chem.MolToMolBlock(set_positions(mol, p_aligned))
