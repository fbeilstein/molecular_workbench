import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation

_PT = Chem.GetPeriodicTable()

try:
    from ..constants import EV_TO_KCAL, DOCK_EXTRA, BREAK_EXTRA
except (ImportError, ValueError):
    from irc_plugins.constants import EV_TO_KCAL, DOCK_EXTRA, BREAK_EXTRA


def rcov(sym: str) -> float:
    return _PT.GetRcovalent(_PT.GetAtomicNumber(sym))


def rvdw(sym: str) -> float:
    return _PT.GetRvdw(_PT.GetAtomicNumber(sym))


def _unit(v: np.ndarray, fallback: np.ndarray = None) -> np.ndarray:
    n = np.linalg.norm(v)
    if n > 1e-9:
        return v / n
    if fallback is not None:
        fb_n = np.linalg.norm(fallback)
        if fb_n > 1e-9:
            return fallback / fb_n
    if v.ndim == 1 and len(v) == 3:
        return np.array([1.0, 0.0, 0.0])
    return v * 0.0


def dir_pen(mode, R, u):
    """Penalty (0 = ideal) for approaching along unit vector u."""
    if mode is None:
        return 0.0
    kind, v = mode
    c = float(np.dot(R @ v, u))
    if kind == "vec":
        return 1.0 - c
    if kind == "axis":
        return 1.0 - c * c
    return c * c  # perp


def align_to_reference(p_orig: np.ndarray, p_new: np.ndarray) -> np.ndarray:
    """Rigid Kabsch alignment: rotate and translate p_new to best fit p_orig."""
    c_orig = np.mean(p_orig, axis=0)
    c_new = np.mean(p_new, axis=0)
    p0 = p_orig - c_orig
    p1 = p_new - c_new
    if len(p0) < 2:
        return p_new - c_new + c_orig
    try:
        rot, _ = Rotation.align_vectors(p0, p1)
        return rot.apply(p1) + c_orig
    except Exception:
        return p_new - c_new + c_orig


def set_positions(mol: Chem.Mol, xyz: np.ndarray) -> Chem.Mol:
    mol = Chem.Mol(mol)
    if mol.GetNumConformers() == 0:
        mol.AddConformer(Chem.Conformer(mol.GetNumAtoms()), assignId=True)
    conf = mol.GetConformer()
    for i, p in enumerate(xyz):
        conf.SetAtomPosition(i, Point3D(float(p[0]), float(p[1]), float(p[2])))
    return mol


def positions_to_molblocks(mols, positions: np.ndarray):
    out, k = [], 0
    for m in mols:
        n = m.GetNumAtoms()
        out.append(Chem.MolToMolBlock(set_positions(m, positions[k:k + n])))
        k += n
    return out


def molblocks_to_positions(molblocks):
    mols = []
    for mb in molblocks:
        m = Chem.MolFromMolBlock(mb, removeHs=False) if mb else None
        if m is None and mb:
            m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
        mols.append(m)
    if any(m is None for m in mols):
        raise ValueError("Failed to parse molblocks")
    return mols, [m.GetConformer().GetPositions() for m in mols]


def clean_smiles(smiles: str) -> str:
    """Normalize SMILES without stereochemistry for robust intermediate matching."""
    if not smiles:
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Chem.RemoveStereochemistry(mol)
            return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        pass
    try:
        return Chem.CanonSmiles(smiles)
    except Exception:
        return smiles


def get_formal_charge(smiles: str) -> int:
    """Extract molecular formal charge from SMILES using RDKit."""
    if not smiles:
        return 0
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.GetFormalCharge(mol) if mol else 0
    except Exception:
        return 0


def get_spin_multiplicity(smiles: str) -> int:
    """Extract spin multiplicity (2S + 1) from SMILES radical count using RDKit."""
    if not smiles:
        return 1
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            return int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())) + 1
    except Exception:
        pass
    return 1


def compute_isolated_fragment_energies(
    fragments,
    solvent: str = "water",
    relax: bool = True,
    fmax: float = 0.03,
    steps: int = 40,
):
    """Compute isolated ground-state xTB electronic energies (eV) for solvated fragments."""
    from ase import Atoms
    from ase.optimize import BFGS
    from rdkit.Chem import AllChem
    from .calculators import make_xtb

    cache = {}
    for fr in fragments:
        smi = fr.get("smiles", "")
        if not smi:
            continue
        c_smi = clean_smiles(smi)
        if c_smi in cache:
            continue

        syms = list(fr.get("symbols", []))
        pos = list(fr.get("positions", []))
        chg = fr.get("charge", get_formal_charge(smi))
        mult = fr.get("multiplicity", get_spin_multiplicity(smi))

        if (not syms or not pos or len(syms) != len(pos)) and c_smi:
            try:
                mol = Chem.AddHs(Chem.MolFromSmiles(c_smi))
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                syms = [a.GetSymbol() for a in mol.GetAtoms()]
                pos = mol.GetConformer().GetPositions().tolist()
                chg = Chem.GetFormalCharge(mol)
                mult = int(sum(a.GetNumRadicalElectrons() for a in mol.GetAtoms())) + 1
            except Exception:
                pass

        if not syms or not pos or len(syms) != len(pos):
            cache[c_smi] = 0.0
            continue

        try:
            atoms = Atoms(symbols=syms, positions=np.array(pos, dtype=float))
            atoms.calc = make_xtb(charge=chg, multiplicity=mult, solvent=solvent)
            if relax and len(syms) > 1:
                try:
                    opt = BFGS(atoms, logfile=None)
                    opt.run(fmax=fmax, steps=steps)
                except Exception:
                    pass
            cache[c_smi] = float(atoms.get_potential_energy())
        except Exception:
            try:
                atoms = Atoms(symbols=syms, positions=np.array(pos, dtype=float))
                atoms.calc = make_xtb(charge=chg, multiplicity=mult, solvent=solvent)
                cache[c_smi] = float(atoms.get_potential_energy())
            except Exception:
                cache[c_smi] = 0.0
    return cache
