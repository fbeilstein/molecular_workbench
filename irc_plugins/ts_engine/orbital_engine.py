import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from pyscf import gto, scf, lo
from skimage import measure
from scipy.optimize import linear_sum_assignment
import trimesh

BOHR_TO_ANG = 0.52917721092
ANG_TO_BOHR = 1.8897259886


def build_global_grid(
    all_positions: List[List[List[float]]],
    grid_points: int = 32,
    margin_angstrom: float = 2.5
) -> Tuple[np.ndarray, np.ndarray, Tuple[int, int, int]]:
    """Compute invariant static 3D bounding box across all frames."""
    all_coords = np.vstack(all_positions)
    lo = all_coords.min(axis=0) - margin_angstrom
    hi = all_coords.max(axis=0) + margin_angstrom
    n = grid_points
    step = (hi - lo) / (n - 1)
    return lo, step, (n, n, n)


def compute_frame_scf(
    symbols: List[str],
    positions: List[List[float]],
    charge: int = -1,
    spin: int = 0,
    basis: str = "sto-3g"
) -> Tuple[Any, Any]:
    """Run RHF/UHF on a frame geometry and return (mol, mf)."""
    atom_str = "; ".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(symbols, positions))
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol) if spin > 0 else scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    return mol, mf


def _classify_ibo(symbols: List[str], pop: np.ndarray) -> Tuple[str, str]:
    ranked = sorted(range(len(symbols)), key=lambda i: -pop[i])
    top = [(i, pop[i]) for i in ranked if pop[i] > 0.25]
    if len(top) == 1 or (len(top) >= 1 and top[0][1] > 1.75):
        ai = top[0][0]
        return "lp", f"LP({symbols[ai]}{ai})"
    elif len(top) >= 2:
        a1, a2 = sorted([top[0][0], top[1][0]])
        return "bond", f"{symbols[a1]}{a1}-{symbols[a2]}{a2}"
    return "core", "core"


def compute_frame_ibos(
    mol: Any,
    mf: Any,
    symbols: List[str],
    prev_data: Optional[Dict[str, Any]] = None
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, np.ndarray], np.ndarray]:
    """Compute Intrinsic Bond Orbitals (IBO) with cross-frame tracking and phase continuity."""
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    occ_mo = mf.mo_coeff[:, :nocc] if mol.spin == 0 else mf.mo_coeff[0][:, :nocc]
    ibo_coeff = lo.ibo.ibo(mol, occ_mo, verbose=0)

    ovlp = mol.intor_symmetric("int1e_ovlp")
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]

    n_ibo = ibo_coeff.shape[1]
    pops = np.zeros((n_ibo, mol.natm))
    vecs = [ibo_coeff[:, k].copy() for k in range(n_ibo)]
    labels, types = [], []

    for k in range(n_ibo):
        v = vecs[k]
        ps = np.outer(v, v) @ ovlp
        for mu in range(mol.nao):
            pops[k, ao_atom[mu]] += ps[mu, mu]
        pops[k] *= 2.0
        typ, lbl = _classify_ibo(symbols, pops[k])
        types.append(typ)
        labels.append(lbl)

    # Match to previous frame by population fingerprint (Hungarian algorithm)
    if prev_data and "pops" in prev_data and "vecs" in prev_data:
        cost = np.sum((prev_data["pops"][:, None, :] - pops[None, :, :]) ** 2, axis=2)
        r, c = linear_sum_assignment(cost)
        pops = pops[c]
        vecs = [vecs[i] for i in c]
        labels = [labels[i] for i in c]
        types = [types[i] for i in c]

        # Fix phase overlap with previous frame
        prev_vec_list = prev_data["vec_list"]
        for k in range(len(vecs)):
            if k < len(prev_vec_list):
                s = float(prev_vec_list[k] @ ovlp @ vecs[k])
                if s < 0.0:
                    vecs[k] = -vecs[k]

    ibo_dict = {}
    vec_dict = {}
    for k in range(n_ibo):
        key = f"IBO_{k:02d}"
        ibo_dict[key] = {
            "name": key,
            "index": k,
            "type": types[k],
            "label": labels[k],
            "is_core": bool(types[k] == "core"),
            "pop": pops[k].round(2).tolist()
        }
        vec_dict[key] = vecs[k]

    return ibo_dict, vec_dict, pops


def extract_frontier_orbitals(
    mol: Any,
    mf: Any,
    prev_mol: Optional[Any] = None,
    prev_mo_coeffs: Optional[Dict[str, np.ndarray]] = None
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]:
    """Extract canonical frontier MOs with smooth phase sign continuity."""
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    mo_coeff = mf.mo_coeff[0] if mol.spin > 0 else mf.mo_coeff
    mo_energy = mf.mo_energy[0] if mol.spin > 0 else mf.mo_energy

    s_cross = gto.intor_cross("int1e_ovlp", prev_mol, mol) if (prev_mol and prev_mo_coeffs) else None

    target_keys = [
        ("HOMO-1", nocc - 2, "HOMO-1"),
        ("HOMO", nocc - 1, "HOMO: Formate π (Nucleophile lone pair)"),
        ("LUMO", nocc, "LUMO: Carbonyl π* (Acceptor antibond)"),
        ("LUMO+1", nocc + 1, "LUMO+1"),
    ]

    orbitals_meta = {}
    current_vectors = {}

    for name, mo_idx, label in target_keys:
        if 0 <= mo_idx < mo_coeff.shape[1]:
            vec = mo_coeff[:, mo_idx].copy()
            if s_cross is not None and name in prev_mo_coeffs:
                prev_vec = prev_mo_coeffs[name]
                overlap = float(prev_vec @ s_cross @ vec)
                if overlap < 0.0:
                    vec = -vec

            current_vectors[name] = vec
            e_ev = float(mo_energy[mo_idx] * 27.2114)
            orbitals_meta[name] = {
                "name": name,
                "label": label,
                "mo_index": int(mo_idx),
                "energy_ev": round(e_ev, 3),
                "is_occupied": bool(mo_idx < nocc)
            }

    return orbitals_meta, current_vectors


def evaluate_mo_grid(
    mol: Any,
    mo_vector: np.ndarray,
    origin_ang: np.ndarray,
    step_ang: np.ndarray,
    npts: Tuple[int, int, int]
) -> np.ndarray:
    """Evaluate MO on a 3D grid and return volumetric 3D numpy array."""
    nx, ny, nz = npts
    xs = origin_ang[0] + np.arange(nx) * step_ang[0]
    ys = origin_ang[1] + np.arange(ny) * step_ang[1]
    zs = origin_ang[2] + np.arange(nz) * step_ang[2]
    gx, gy, gz = np.meshgrid(xs, ys, zs, indexing="ij")
    grid_bohr = np.column_stack([gx.ravel(), gy.ravel(), gz.ravel()]) * ANG_TO_BOHR
    ao_vals = mol.eval_gto("GTOval", grid_bohr)
    return (ao_vals @ mo_vector).reshape(npts)


def extract_isosurface_mesh(
    vol: np.ndarray,
    origin_ang: np.ndarray,
    step_ang: np.ndarray,
    isovalue: float = 0.03
) -> Dict[str, Any]:
    """Generate 3Dmol-compatible surface meshes for positive and negative lobes."""
    mesh_data = {"pos": None, "neg": None, "isovalue": isovalue}
    level = abs(isovalue)

    try:
        vp, fp, np_arr, _ = measure.marching_cubes(vol, level, spacing=step_ang)
        vp += origin_ang
        tp = trimesh.Trimesh(vertices=vp, faces=fp, vertex_normals=np_arr, process=True)
        v_list = [{"x": float(v[0]), "y": float(v[1]), "z": float(v[2])} for v in tp.vertices.round(4)]
        f_list = [int(idx) for face in tp.faces for idx in face]
        mesh_data["pos"] = {
            "vertexArr": v_list,
            "faceArr": f_list,
            "normalArr": [{"x": float(n[0]), "y": float(n[1]), "z": float(n[2])} for n in tp.vertex_normals.round(4)]
        }
    except Exception:
        mesh_data["pos"] = None

    try:
        vn, fn, nn_arr, _ = measure.marching_cubes(vol, -level, spacing=step_ang)
        vn += origin_ang
        fn = fn[:, ::-1]
        tn = trimesh.Trimesh(vertices=vn, faces=fn, vertex_normals=-nn_arr, process=True)
        v_list = [{"x": float(v[0]), "y": float(v[1]), "z": float(v[2])} for v in tn.vertices.round(4)]
        f_list = [int(idx) for face in tn.faces for idx in face]
        mesh_data["neg"] = {
            "vertexArr": v_list,
            "faceArr": f_list,
            "normalArr": [{"x": float(n[0]), "y": float(n[1]), "z": float(n[2])} for n in tn.vertex_normals.round(4)]
        }
    except Exception:
        mesh_data["neg"] = None

    return mesh_data


def generate_cube_text(
    comment: str,
    symbols: List[str],
    positions: List[List[float]],
    vol: np.ndarray,
    origin_ang: np.ndarray,
    step_ang: np.ndarray,
    npts: Tuple[int, int, int]
) -> str:
    """Generate standard Gaussian .cube format text for 3Dmol VolumeData."""
    nx, ny, nz = npts
    orig_bohr = origin_ang * ANG_TO_BOHR
    step_bohr = step_ang * ANG_TO_BOHR
    lines = [
        f"{comment}\n",
        "Generated by Molecular Workbench Orbital Engine\n",
        f"{len(symbols):5d}{orig_bohr[0]:12.6f}{orig_bohr[1]:12.6f}{orig_bohr[2]:12.6f}\n",
        f"{nx:5d}{step_bohr[0]:12.6f}{0.0:12.6f}{0.0:12.6f}\n",
        f"{ny:5d}{0.0:12.6f}{step_bohr[1]:12.6f}{0.0:12.6f}\n",
        f"{nz:5d}{0.0:12.6f}{0.0:12.6f}{step_bohr[2]:12.6f}\n"
    ]
    atom_z = {"H": 1, "C": 6, "N": 7, "O": 8, "F": 9, "Cl": 17, "Br": 35, "I": 53}
    for s, p in zip(symbols, positions):
        pb = np.array(p) * ANG_TO_BOHR
        lines.append(f"{atom_z.get(s, 0):5d}{0.0:12.6f}{pb[0]:12.6f}{pb[1]:12.6f}{pb[2]:12.6f}\n")

    flat = vol.ravel()
    for i in range(0, len(flat), 6):
        chunk = flat[i : i + 6]
        lines.append("".join(f"{v:13.5E}" for v in chunk) + "\n")

    return "".join(lines)
