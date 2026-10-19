import os
import sys
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from pyscf import gto, lo
from scipy.optimize import linear_sum_assignment

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

from qm.classification import _classify_orbital


def classify_ibo_quick(symbols: List[str], pop: np.ndarray) -> Tuple[str, str]:
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
        typ_s, lbl_s = classify_ibo_quick(symbols, pops_start[k])
        typ_e, lbl_e = classify_ibo_quick(symbols, pops_end[k])
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
    ovlp = mol.intor_symmetric("int1e_ovlp")
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]

    n_val = ibo_val.shape[1]
    val_pops = np.zeros((n_val, mol.natm))
    for k in range(n_val):
        v = ibo_val[:, k]
        ps = np.outer(v, v) @ ovlp
        for mu in range(mol.nao):
            val_pops[k, ao_atom[mu]] += ps[mu, mu]
        val_pops[k] *= 2.0

    # Match and orient core and valence orbitals using cross-frame population continuity and overlap
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

        # 2. Valence matching via 3D absolute cross-overlap (retains 3D spatial & nodal symmetry)
        prev_val = prev_data["ibo_coeff"][:, n_prev_core:]
        ovlp_val = prev_val.T @ s_cross @ ibo_val
        r, c = linear_sum_assignment(-np.abs(ovlp_val))
        ibo_val = ibo_val[:, c]
        val_pops = val_pops[c]

        # Phase alignment along the diagonal
        for k in range(n_val):
            s_val = float(prev_val[:, k].T @ s_cross @ ibo_val[:, k])
            if s_val < 0.0:
                ibo_val[:, k] = -ibo_val[:, k]

        # 3. Generalized Subspace Procrustes SVD Alignment
        # Only allow unitary rotation within near-degenerate manifolds (|Delta eps| < 1.5 eV).
        # Mixing non-degenerate orbitals (e.g. sigma-bond at -25 eV with LP at -6 eV) is unphysical
        # and creates artificial chimeric orbitals due to geometry translation cross-overlap.
        F = mf.get_fock()
        eps_val = np.array([float(ibo_val[:, k].T @ F @ ibo_val[:, k]) * 27.2114 for k in range(n_val)])

        O_val = prev_val.T @ s_cross @ ibo_val
        adj = {i: set() for i in range(n_val)}
        for i in range(n_val):
            for j in range(i + 1, n_val):
                if max(abs(O_val[i, j]), abs(O_val[j, i])) > 0.20 and abs(eps_val[i] - eps_val[j]) < 1.5:
                    adj[i].add(j)
                    adj[j].add(i)

        visited = set()
        for i in range(n_val):
            if i not in visited and adj[i]:
                comp = []
                q = [i]
                visited.add(i)
                while q:
                    curr = q.pop(0)
                    comp.append(curr)
                    for nb in adj[curr]:
                        if nb not in visited:
                            visited.add(nb)
                            q.append(nb)

                if len(comp) >= 2:
                    sub_prev = prev_val[:, comp]
                    sub_curr = ibo_val[:, comp]
                    M_sub = sub_prev.T @ s_cross @ sub_curr
                    V, s, Wt = np.linalg.svd(M_sub)
                    U = Wt.T @ V.T
                    if np.linalg.det(U) < 0:
                        Wt_fixed = Wt.copy()
                        Wt_fixed[-1, :] *= -1
                        U = Wt_fixed.T @ V.T

                    ibo_val[:, comp] = sub_curr @ U

        for k in range(n_val):
            s_val = float(prev_val[:, k].T @ s_cross @ ibo_val[:, k])
            if s_val < 0.0:
                ibo_val[:, k] = -ibo_val[:, k]

        # Recalculate atomic populations from final rotated coefficients
        for idx in range(n_val):
            v = ibo_val[:, idx]
            ps = np.outer(v, v) @ ovlp
            val_pops[idx] = 0.0
            for mu in range(mol.nao):
                val_pops[idx, ao_atom[mu]] += ps[mu, mu]
            val_pops[idx] *= 2.0

    ibo_coeff = np.hstack([core_mo, ibo_val])

    pops = np.zeros((nocc, mol.natm))
    n_c = core_mo.shape[1]
    for k in range(n_c):
        v = core_mo[:, k]
        ps = np.outer(v, v) @ ovlp
        for mu in range(mol.nao):
            pops[k, ao_atom[mu]] += ps[mu, mu]
        pops[k] *= 2.0
    pops[n_c:] = val_pops

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
            "key": key, "name": key, "index": k, "type": typ, "category": cat,
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
                        "key": star_key, "name": star_key, "index": 100 + k, "type": star_type,
                        "category": "antibonding", "label": f"{star_sym}({a1}-{a2}) [Antibond]",
                        "atoms": [a1, a2], "chem_key": f"{star_type}_{a1}-{a2}",
                        "energy_ev": round(e_star, 2), "pop": []
                    }
                    vec_dict[star_key] = c_star

    return ibo_dict, star_dict, vec_dict, pops, ibo_coeff

