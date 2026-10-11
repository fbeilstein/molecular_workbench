import numpy as np
from rdkit import Chem
from rdkit.Geometry import Point3D
from scipy.spatial.transform import Rotation

_PT = Chem.GetPeriodicTable()
EV_TO_KCAL = 23.060548

DOCK_EXTRA = 1.0      # pre-reaction complex: forming bond = d0 + 1.0 A
BREAK_EXTRA = 1.4     # end of drive: breaking bond = d0 + 1.4 A


def rcov(sym: str) -> float:
    return _PT.GetRcovalent(_PT.GetAtomicNumber(sym))


def rvdw(sym: str) -> float:
    return _PT.GetRvdw(_PT.GetAtomicNumber(sym))


def _unit(v: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(v)
    return v / n if n > 1e-9 else v * 0.0


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
    mols = [Chem.MolFromMolBlock(mb, removeHs=False) for mb in molblocks]
    if any(m is None for m in mols):
        raise ValueError("Failed to parse molblocks")
    return mols, [m.GetConformer().GetPositions() for m in mols]
