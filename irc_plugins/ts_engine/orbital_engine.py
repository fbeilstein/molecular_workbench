import os
import sys
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from pyscf import gto, scf, lo
from skimage import measure
from scipy.optimize import linear_sum_assignment
import trimesh

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

from qm.analyzers.canonical import CanonicalAnalyzer
from qm.classification import _atom_populations, _classify_orbital, _classify_sigma_pi

BOHR_TO_ANG = 0.52917721092
ANG_TO_BOHR = 1.8897259886


def build_global_grid(
    all_positions: List[List[List[float]]],
    grid_points: int = 50,
    margin_angstrom: float = 2.5
) -> Tuple[np.ndarray, np.ndarray, Tuple[int, int, int]]:
    """Compute invariant static 3D bounding box across all frames."""
    all_coords = np.vstack(all_positions)
    lo_corner = all_coords.min(axis=0) - margin_angstrom
    hi_corner = all_coords.max(axis=0) + margin_angstrom
    n = grid_points
    step = (hi_corner - lo_corner) / (n - 1)
    return lo_corner, step, (n, n, n)


def compute_frame_scf(
    symbols: List[str],
    positions: List[List[float]],
    charge: int = -1,
    spin: int = 0,
    basis: str = "sto-3g"
) -> Tuple[Any, Any, List[float]]:
    """Run RHF/UHF on a frame geometry and return (mol, mf, mulliken_charges)."""
    atom_str = "; ".join(f"{s} {p[0]} {p[1]} {p[2]}" for s, p in zip(symbols, positions))
    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=0)
    mf = scf.UHF(mol) if spin > 0 else scf.RHF(mol)
    mf.verbose = 0
    mf.kernel()
    try:
        mull = mf.mulliken_pop(verbose=0)[1].round(3).tolist()
    except Exception:
        mull = [0.0] * len(symbols)
    return mol, mf, mull


def _classify_ibo_quick(symbols: List[str], pop: np.ndarray) -> Tuple[str, str]:
    ranked = sorted(range(len(symbols)), key=lambda i: -pop[i])
    top = [(i, pop[i]) for i in ranked if pop[i] > 0.25]
    if len(top) == 1 or (len(top) >= 1 and top[0][1] > 1.75):
        ai = top[0][0]
        return "lp", f"LP({symbols[ai]}{ai})"
    elif len(top) >= 2:
        a1, a2 = sorted([top[0][0], top[1][0]])
        return "bond", f"{symbols[a1]}{a1}-{symbols[a2]}{a2}"
    return "core", "core"


def compute_density_flux(
    pops_start: np.ndarray,
    pops_end: np.ndarray,
    symbols: List[str]
) -> List[Dict[str, Any]]:
    """Rank orbitals by net population redistribution between two frames."""
    nocc = pops_start.shape[0]
    diff = np.abs(pops_end - pops_start)
    delta_sums = diff.sum(axis=1)
    fluxes = []
    for k in range(nocc):
        typ_s, lbl_s = _classify_ibo_quick(symbols, pops_start[k])
        typ_e, lbl_e = _classify_ibo_quick(symbols, pops_end[k])
        fluxes.append({
            "key": f"IBO_{k:02d}",
            "index": k,
            "delta_pop": round(float(delta_sums[k]), 3),
            "label_start": lbl_s,
            "label_end": lbl_e,
            "type": typ_s
        })
    fluxes.sort(key=lambda x: -x["delta_pop"])
    return fluxes


def compute_frame_ibos(
    mol: Any,
    mf: Any,
    symbols: List[str],
    prev_data: Optional[Dict[str, Any]] = None
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """Compute all IBOs with cross-frame orbital overlap continuity and textbook antibonds."""
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    occ_mo = mf.mo_coeff[:, :nocc] if mol.spin == 0 else mf.mo_coeff[0][:, :nocc]
    mo_e = (mf.mo_energy[:nocc] if mol.spin == 0 else mf.mo_energy[0][:nocc]) * 27.2114
    core_mask = mo_e < -50.0
    core_mo = occ_mo[:, core_mask]
    val_mo = occ_mo[:, ~core_mask]
    ibo_val = lo.ibo.ibo(mol, val_mo, verbose=0)

    # Match and orient core and valence orbitals using cross-frame overlap matrix
    if prev_data and "mol" in prev_data and "ibo_coeff" in prev_data:
        n_prev_core = prev_data.get("n_core", core_mo.shape[1])
        s_cross = gto.intor_cross("int1e_ovlp", prev_data["mol"], mol)
        if core_mo.shape[1] == n_prev_core:
            prev_core = prev_data["ibo_coeff"][:, :n_prev_core]
            ovlp_core = prev_core.T @ s_cross @ core_mo
            rc, cc = linear_sum_assignment(-np.abs(ovlp_core))
            core_mo = core_mo[:, cc]
            for k in range(core_mo.shape[1]):
                if ovlp_core[k, cc[k]] < 0.0:
                    core_mo[:, k] = -core_mo[:, k]

        prev_val = prev_data["ibo_coeff"][:, n_prev_core:]
        ovlp_mat = prev_val.T @ s_cross @ ibo_val
        r, c = linear_sum_assignment(-np.abs(ovlp_mat))
        ibo_val = ibo_val[:, c]
        for k in range(ibo_val.shape[1]):
            if ovlp_mat[k, c[k]] < 0.0:
                ibo_val[:, k] = -ibo_val[:, k]

    ibo_coeff = np.hstack([core_mo, ibo_val])

    ovlp = mol.intor_symmetric("int1e_ovlp")
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]

    pops = np.zeros((nocc, mol.natm))
    for k in range(nocc):
        v = ibo_coeff[:, k]
        ps = np.outer(v, v) @ ovlp
        for mu in range(mol.nao):
            pops[k, ao_atom[mu]] += ps[mu, mu]
        pops[k] *= 2.0

    fock = mf.get_fock()
    ibo_fock = ibo_coeff.T @ fock @ ibo_coeff
    energies_ev = np.diag(ibo_fock) * 27.2114

    atom_labels = [mol.atom_symbol(i) for i in range(mol.natm)]
    atom_ids = [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
    ao_labels = mol.ao_labels(fmt=False)
    virt_coeff = mf.mo_coeff[:, nocc:] if mol.spin == 0 else mf.mo_coeff[0][:, nocc:]

    ibo_dict = {}
    star_dict = {}
    vec_dict = {}

    for k in range(nocc):
        key = f"IBO_{k:02d}"
        vec = ibo_coeff[:, k]
        e_ev = float(energies_ev[k])

        if e_ev < -50.0:
            ai = int(np.argmax(pops[k]))
            typ = "core"
            cat = "core"
            lbl = f"Core({symbols[ai]}{ai})"
            atoms = [f"{symbols[ai]}{ai}"]
            c_key = f"core_{symbols[ai]}{ai}_{k}"
        else:
            info = _classify_orbital(mol, vec, pops[k], atom_labels, atom_ids, ovlp)
            typ = info.get("type", "unknown")
            atoms = info.get("atoms", [info.get("atom", "?")])

            if len(atoms) >= 2:
                def _atom_sort_key(at_str):
                    sym = "".join(filter(str.isalpha, at_str))
                    num = int("".join(filter(str.isdigit, at_str)) or 0)
                    prio = 0 if sym == "C" else (2 if sym == "H" else 1)
                    return (prio, num, sym)
                atoms = sorted(atoms, key=_atom_sort_key)

            if typ == "lone_pair":
                cat = "non_bonding"
                lbl = f"Non-Bonding LP({atoms[0]})"
                c_key = f"lp_{atoms[0]}"
            elif typ in ("sigma", "delocalized_sigma"):
                cat = "bonding"
                lbl = f"σ({'-'.join(atoms)})"
                c_key = f"sigma_{'-'.join(atoms)}"
            elif typ in ("pi", "delocalized_pi"):
                cat = "bonding"
                lbl = f"π({'-'.join(atoms)})"
                c_key = f"pi_{'-'.join(atoms)}"
            else:
                cat = "other"
                lbl = f"{typ}({'-'.join(atoms)})"
                c_key = f"{typ}_{'-'.join(atoms)}"

        ibo_dict[key] = {
            "name": key, "index": k, "type": typ, "category": cat,
            "label": lbl, "atoms": atoms,
            "chem_key": c_key,
            "energy_ev": round(e_ev, 2),
            "pop": pops[k].round(2).tolist()
        }
        vec_dict[key] = vec

        # Textbook antisymmetrized virtual projected antibond for 2-center bonds
        if virt_coeff is not None and typ in ("sigma", "delocalized_sigma", "pi", "delocalized_pi") and len(atoms) >= 2:
            a1, a2 = atoms[0], atoms[1]
            if a1 in atom_ids and a2 in atom_ids:
                idx1 = atom_ids.index(a1)
                idx2 = atom_ids.index(a2)
                cA = np.zeros_like(vec)
                cB = np.zeros_like(vec)
                for mu, (ai, *_) in enumerate(ao_labels):
                    if ai == idx1:
                        cA[mu] = vec[mu]
                    elif ai == idx2:
                        cB[mu] = vec[mu]
                c_anti = cA - cB
                c_virt = virt_coeff.T @ ovlp @ c_anti
                c_star = virt_coeff @ c_virt
                norm = np.sqrt(float(c_star.T @ ovlp @ c_star))
                if norm > 1e-3:
                    c_star /= norm
                    star_key = f"STAR_{k:02d}"
                    star_type = "sigma_star" if "sigma" in typ else "pi_star"
                    star_sym = "σ*" if "sigma" in typ else "π*"
                    e_star = float(c_star.T @ fock @ c_star) * 27.2114
                    star_dict[star_key] = {
                        "name": star_key, "index": 100 + k, "type": star_type,
                        "category": "antibonding", "label": f"{star_sym}({a1}-{a2}) [Antibond]",
                        "atoms": [a1, a2], "chem_key": f"{star_type}_{a1}-{a2}",
                        "energy_ev": round(e_star, 2), "pop": []
                    }
                    vec_dict[star_key] = c_star

    return ibo_dict, star_dict, vec_dict, pops, ibo_coeff


def extract_frontier_orbitals(
    mol: Any,
    mf: Any,
    prev_mol: Optional[Any] = None,
    prev_mo_coeffs: Optional[Dict[str, np.ndarray]] = None,
    ibos_m: Optional[Dict[str, Any]] = None,
    star_m: Optional[Dict[str, Any]] = None,
    ibos_v: Optional[Dict[str, np.ndarray]] = None
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]:
    """Extract canonical frontier MOs and clean fragment delocalized pi systems with smooth phase continuity."""
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    mo_coeff = mf.mo_coeff[0] if mol.spin > 0 else mf.mo_coeff
    mo_energy = mf.mo_energy[0] if mol.spin > 0 else mf.mo_energy

    s_cross = gto.intor_cross("int1e_ovlp", prev_mol, mol) if (prev_mol and prev_mo_coeffs) else None

    target_keys = [
        ("HOMO-3", nocc - 4, "HOMO-3"),
        ("HOMO-2", nocc - 3, "HOMO-2"),
        ("HOMO-1", nocc - 2, "HOMO-1"),
        ("HOMO", nocc - 1, "HOMO (Frontier Donor)"),
        ("LUMO", nocc, "LUMO (Frontier Acceptor)"),
        ("LUMO+1", nocc + 1, "LUMO+1"),
        ("LUMO+2", nocc + 2, "LUMO+2"),
    ]

    orbitals_meta = {}
    frontier_vectors = {}

    for name, mo_idx, label in target_keys:
        if 0 <= mo_idx < mo_coeff.shape[1]:
            vec = mo_coeff[:, mo_idx].copy()
            if s_cross is not None and name in prev_mo_coeffs:
                prev_vec = prev_mo_coeffs[name]
                if float(prev_vec @ s_cross @ vec) < 0.0:
                    vec = -vec

            frontier_vectors[name] = vec
            e_ev = float(mo_energy[mo_idx] * 27.2114)
            orbitals_meta[name] = {
                "name": name,
                "label": label,
                "category": "frontier",
                "type": "frontier",
                "mo_index": int(mo_idx),
                "energy_ev": round(e_ev, 3),
                "is_occupied": bool(mo_idx < nocc)
            }

    # Clean Fragment Quasi-Canonical Delocalized π Systems
    canonical_pi_meta = {}
    canonical_pi_vecs = {}

    coords = mol.atom_coords() * BOHR_TO_ANG
    # Exclude saturated aliphatic carbons (>= 2 bonded hydrogens, e.g. -CH3, -CH2-)
    conjugated_candidates = []
    for i in range(mol.natm):
        sym = mol.atom_symbol(i)
        if sym == "H":
            continue
        if sym == "C":
            n_h = sum(1 for j in range(mol.natm) if mol.atom_symbol(j) == "H" and np.linalg.norm(coords[i] - coords[j]) < 1.25)
            if n_h >= 2:
                continue
        conjugated_candidates.append(i)

    adj = {i: [] for i in conjugated_candidates}
    for i in range(len(conjugated_candidates)):
        for j in range(i + 1, len(conjugated_candidates)):
            ai, aj = conjugated_candidates[i], conjugated_candidates[j]
            cutoff = 2.15 if (mol.atom_symbol(ai) == "Cl" or mol.atom_symbol(aj) == "Cl") else 1.65
            if np.linalg.norm(coords[ai] - coords[aj]) < cutoff:
                adj[ai].append(aj)
                adj[aj].append(ai)

    visited = set()
    frags = []
    for node in conjugated_candidates:
        if node not in visited:
            comp = []
            q = [node]
            visited.add(node)
            while q:
                curr = q.pop(0)
                comp.append(curr)
                for nb in adj[curr]:
                    if nb not in visited:
                        visited.add(nb)
                        q.append(nb)
            if len(comp) >= 3:
                frags.append(sorted(comp))

    fock = mf.get_fock()
    ovlp_mat = mol.intor_symmetric("int1e_ovlp")

    for f in frags:
        atom_ids = [f"{mol.atom_symbol(i)}{i}" for i in f]
        frag_lbl = "-".join(atom_ids)
        coords_f = coords[f]
        c = np.mean(coords_f, axis=0)
        cov = (coords_f - c).T @ (coords_f - c)
        evals, evecs = np.linalg.eigh(cov)
        normal = evecs[:, 0]

        # 1. Identify localized pi bonds and out-of-plane p-lone pairs on this fragment
        basis_keys = []
        basis_vecs = []

        if ibos_m and ibos_v:
            for k, d in ibos_m.items():
                typ = d.get("type", "")
                ats = d.get("atoms", [])
                if typ in ("pi", "delocalized_pi") and all(at in atom_ids for at in ats):
                    basis_keys.append(k)
                    basis_vecs.append(ibos_v[k])
                elif typ == "lone_pair" and ats and (ats[0] in atom_ids):
                    idx = int("".join(filter(str.isdigit, ats[0])))
                    pt_up = (coords[idx] + 0.5 * normal).reshape(1, 3) * ANG_TO_BOHR
                    pt_dn = (coords[idx] - 0.5 * normal).reshape(1, 3) * ANG_TO_BOHR
                    v = ibos_v[k]
                    v_up = float((mol.eval_gto("GTOval", pt_up) @ v)[0])
                    v_dn = float((mol.eval_gto("GTOval", pt_dn) @ v)[0])
                    # An out-of-plane conjugated p-lone pair has anti-symmetric lobes across the plane
                    if v_up * v_dn < 0.0 and abs(v_up) > 0.04:
                        basis_keys.append(k)
                        basis_vecs.append(v)

        # 2. If we have localized pi basis orbitals on this fragment, diagonalize the subspace Fock matrix
        if len(basis_keys) >= 2:
            aligned_vecs = []
            for k, v in zip(basis_keys, basis_vecs):
                v_copy = v.copy()
                at = ibos_m[k]["atoms"][0]
                idx = int("".join(filter(str.isdigit, at)))
                pt_up = (coords[idx] + 0.5 * normal).reshape(1, 3) * ANG_TO_BOHR
                val = float((mol.eval_gto("GTOval", pt_up) @ v_copy)[0])
                if val < 0.0:
                    v_copy = -v_copy
                aligned_vecs.append(v_copy)

            C_pi = np.column_stack(aligned_vecs)
            F_sub = C_pi.T @ fock @ C_pi
            e_vals, e_vecs = np.linalg.eigh(F_sub)
            C_deloc = C_pi @ e_vecs

            for i in range(len(e_vals)):
                vec = C_deloc[:, i].copy()
                norm = np.sqrt(float(vec.T @ ovlp_mat @ vec))
                if norm > 1e-4:
                    vec /= norm

                key = f"DELOC_PI_{i+1:02d}"
                if s_cross is not None and key in prev_mo_coeffs:
                    if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                        vec = -vec

                canonical_pi_vecs[key] = vec
                node_str = f"{i} node{'s' if i != 1 else ''}" if i > 0 else "all-bonding, 0 nodes"
                canonical_pi_meta[key] = {
                    "name": key,
                    "label": f"Delocalized π ({frag_lbl}) [π{i+1}, {node_str}]",
                    "category": "conjugated_pi",
                    "type": "delocalized_pi",
                    "chem_key": f"deloc_pi_{i+1}",
                    "mo_index": i + 1,
                    "energy_ev": round(float(e_vals[i] * 27.2114), 2),
                    "atoms": atom_ids
                }
                # Alias for compatibility with canonical naming
                canonical_pi_vecs[f"CANONICAL_PI_{i+1:02d}"] = vec

            # 3. Virtual antibonding subspace (pi*)
            virt_keys = []
            virt_vecs = []
            if star_m and ibos_v:
                for k, d in star_m.items():
                    if d.get("type") == "pi_star" and all(at in atom_ids for at in d.get("atoms", [])):
                        virt_keys.append(k)
                        virt_vecs.append(ibos_v[k])

            if len(virt_keys) >= 2:
                aligned_virt = []
                for k, v in zip(virt_keys, virt_vecs):
                    v_copy = v.copy()
                    at = star_m[k]["atoms"][0]
                    idx = int("".join(filter(str.isdigit, at)))
                    pt_up = (coords[idx] + 0.5 * normal).reshape(1, 3) * ANG_TO_BOHR
                    val = float((mol.eval_gto("GTOval", pt_up) @ v_copy)[0])
                    if val < 0.0:
                        v_copy = -v_copy
                    aligned_virt.append(v_copy)

                C_pistar = np.column_stack(aligned_virt)
                F_virt_sub = C_pistar.T @ fock @ C_pistar
                e_vals_star, e_vecs_star = np.linalg.eigh(F_virt_sub)
                C_deloc_star = C_pistar @ e_vecs_star

                for i in range(len(e_vals_star)):
                    vec = C_deloc_star[:, i].copy()
                    norm = np.sqrt(float(vec.T @ ovlp_mat @ vec))
                    if norm > 1e-4:
                        vec /= norm

                    key = f"DELOC_PISTAR_{i+1:02d}"
                    if s_cross is not None and key in prev_mo_coeffs:
                        if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                            vec = -vec

                    canonical_pi_vecs[key] = vec
                    canonical_pi_meta[key] = {
                        "name": key,
                        "label": f"Delocalized π* ({frag_lbl}) [π*{i+1}]",
                        "category": "conjugated_pi",
                        "type": "delocalized_pistar",
                        "chem_key": f"deloc_pistar_{i+1}",
                        "mo_index": 100 + i + 1,
                        "energy_ev": round(float(e_vals_star[i] * 27.2114), 2),
                        "atoms": atom_ids
                    }
                    canonical_pi_vecs[f"CANONICAL_PISTAR_{i+1:02d}"] = vec

        else:
            # Fallback to CanonicalAnalyzer if no localized pi basis found
            try:
                ca = CanonicalAnalyzer(mol, mo_coeff, mo_energy)
                pi_occ, pi_virt = ca.get_pi_system(f)
                for idx in pi_occ:
                    depth = nocc - 1 - idx
                    hlbl = f"HOMO-{depth}" if depth > 0 else "HOMO"
                    key = f"CANONICAL_PI_{idx:02d}"
                    vec = mo_coeff[:, idx].copy()
                    if s_cross is not None and key in prev_mo_coeffs:
                        if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                            vec = -vec
                    canonical_pi_vecs[key] = vec
                    canonical_pi_meta[key] = {
                        "name": key,
                        "label": f"Delocalized π ({frag_lbl}) [{hlbl}]",
                        "category": "conjugated_pi",
                        "type": "delocalized_pi",
                        "mo_index": int(idx),
                        "energy_ev": round(float(mo_energy[idx] * 27.2114), 2),
                        "atoms": atom_ids
                    }
                for idx in pi_virt[:2]:
                    depth = idx - nocc
                    llbl = f"LUMO+{depth}" if depth > 0 else "LUMO"
                    key = f"CANONICAL_PISTAR_{idx:02d}"
                    vec = mo_coeff[:, idx].copy()
                    if s_cross is not None and key in prev_mo_coeffs:
                        if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                            vec = -vec
                    canonical_pi_vecs[key] = vec
                    canonical_pi_meta[key] = {
                        "name": key,
                        "label": f"Delocalized π* ({frag_lbl}) [{llbl}]",
                        "category": "conjugated_pi",
                        "type": "delocalized_pistar",
                        "mo_index": int(idx),
                        "energy_ev": round(float(mo_energy[idx] * 27.2114), 2),
                        "atoms": atom_ids
                    }
            except Exception:
                pass

    all_vecs = {**frontier_vectors, **canonical_pi_vecs}
    return orbitals_meta, canonical_pi_meta, all_vecs


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
