import os
import sys
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from pyscf import gto

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

from qm.analyzers.canonical import CanonicalAnalyzer

BOHR_TO_ANG = 0.52917721092
ANG_TO_BOHR = 1.8897259886


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
                    if v_up * v_dn < 0.0 and abs(v_up) > 0.04:
                        basis_keys.append(k)
                        basis_vecs.append(v)

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

                key = f"DELOC_PI_{frag_lbl}_{i+1:02d}"
                if s_cross is not None and key in prev_mo_coeffs:
                    if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                        vec = -vec

                canonical_pi_vecs[key] = vec
                canonical_pi_vecs[f"DELOC_PI_{i+1:02d}"] = vec
                canonical_pi_vecs[f"CANONICAL_PI_{i+1:02d}"] = vec
                node_str = f"{i} node{'s' if i != 1 else ''}" if i > 0 else "all-bonding, 0 nodes"
                canonical_pi_meta[key] = {
                    "name": key,
                    "label": f"Delocalized π ({frag_lbl}) [π{i+1}, {node_str}]",
                    "category": "conjugated_pi",
                    "type": "delocalized_pi",
                    "chem_key": f"deloc_pi_{frag_lbl}_{i+1}",
                    "mo_index": i + 1,
                    "energy_ev": round(float(e_vals[i] * 27.2114), 2),
                    "atoms": atom_ids
                }

            virt_keys = []
            virt_vecs = []
            if star_m and ibos_v:
                for k, d in star_m.items():
                    if d.get("type") == "pi_star" and any(at in atom_ids for at in d.get("atoms", [])):
                        virt_keys.append(k)
                        virt_vecs.append(ibos_v[k])

            if len(virt_keys) == 1:
                k = virt_keys[0]
                vec = virt_vecs[0].copy()
                key = f"PISTAR_{frag_lbl}_01"
                canonical_pi_vecs[key] = vec
                canonical_pi_meta[key] = {
                    "name": key,
                    "label": f"π* ({frag_lbl}) [Antibond]",
                    "category": "antibonding",
                    "type": "delocalized_pistar",
                    "chem_key": f"pistar_{frag_lbl}_1",
                    "mo_index": 101,
                    "energy_ev": star_m[k].get("energy_ev", 0.0),
                    "atoms": atom_ids
                }
            elif len(virt_keys) >= 2:
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

                    key = f"DELOC_PISTAR_{frag_lbl}_{i+1:02d}"
                    if s_cross is not None and key in prev_mo_coeffs:
                        if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                            vec = -vec

                    canonical_pi_vecs[key] = vec
                    canonical_pi_meta[key] = {
                        "name": key,
                        "label": f"Delocalized π* ({frag_lbl}) [π*{i+1}]",
                        "category": "antibonding",
                        "type": "delocalized_pistar",
                        "chem_key": f"deloc_pistar_{frag_lbl}_{i+1}",
                        "mo_index": 100 + i + 1,
                        "energy_ev": round(float(e_vals_star[i] * 27.2114), 2),
                        "atoms": atom_ids
                    }
                    canonical_pi_vecs[f"CANONICAL_PISTAR_{i+1:02d}"] = vec

            try:
                ca = CanonicalAnalyzer(mol, mo_coeff, mo_energy)
                _, pi_virt = ca.get_pi_system(f)
                for idx in pi_virt:
                    depth = idx - nocc
                    llbl = f"LUMO+{depth}" if depth > 0 else "LUMO"
                    key = f"CANONICAL_PISTAR_{frag_lbl}_{idx:02d}"
                    if key not in canonical_pi_vecs and f"CANONICAL_PISTAR_{idx:02d}" not in canonical_pi_vecs:
                        vec = mo_coeff[:, idx].copy()
                        if s_cross is not None and key in prev_mo_coeffs:
                            if float(prev_mo_coeffs[key] @ s_cross @ vec) < 0.0:
                                vec = -vec
                        canonical_pi_vecs[key] = vec
                        canonical_pi_vecs[f"CANONICAL_PISTAR_{idx:02d}"] = vec
                        canonical_pi_meta[key] = {
                            "name": key,
                            "label": f"Canonical π* ({frag_lbl}) [{llbl}]",
                            "category": "antibonding",
                            "type": "delocalized_pistar",
                            "chem_key": f"pistar_{frag_lbl}_{idx}",
                            "mo_index": int(idx),
                            "energy_ev": round(float(mo_energy[idx] * 27.2114), 2),
                            "atoms": atom_ids
                        }
            except Exception:
                pass

        else:
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
