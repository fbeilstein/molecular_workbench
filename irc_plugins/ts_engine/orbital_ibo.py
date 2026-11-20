import os
import sys
import numpy as np
from typing import Any, Dict, List, Optional, Tuple
from pyscf import gto
from scipy.optimize import linear_sum_assignment

BASE_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
TOOLS_DIR = os.path.join(BASE_DIR, "tools")
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)

from qm.analyzers.valence import ValenceLocalizer
from qm.classification import count_core_orbitals
from qm.constants import HARTREE_TO_EV


def classify_ibo_quick(symbols: List[str], pop: np.ndarray, atom_labels: Optional[List[str]] = None) -> Tuple[str, str]:
    ranked = sorted(range(len(symbols)), key=lambda i: -pop[i])
    top = [(i, pop[i]) for i in ranked if pop[i] > 0.25]
    tags = atom_labels if (atom_labels and len(atom_labels) == len(symbols)) else [f"{symbols[i]}{i}" for i in range(len(symbols))]
    if len(top) == 1 or (len(top) >= 1 and top[0][1] > 1.75):
        ai = top[0][0]
        return "lp", f"LP({tags[ai]})"
    elif len(top) >= 2:
        a1, a2 = sorted([top[0][0], top[1][0]])
        return "bond", f"{tags[a1]}-{tags[a2]}"
    return "core", "core"


def compute_density_flux(
    pops_start: np.ndarray, pops_end: np.ndarray, symbols: List[str], atom_labels: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """Rank orbitals by net population redistribution between two frames."""
    nocc = pops_start.shape[0]
    delta_sums = np.abs(pops_end - pops_start).sum(axis=1)
    fluxes = []
    for k in range(nocc):
        typ_s, lbl_s = classify_ibo_quick(symbols, pops_start[k], atom_labels)
        typ_e, lbl_e = classify_ibo_quick(symbols, pops_end[k], atom_labels)
        fluxes.append({
            "key": f"IBO_{k:02d}", "index": k, "delta_pop": round(float(delta_sums[k]), 3),
            "label_start": lbl_s, "label_end": lbl_e, "type": typ_s
        })
    fluxes.sort(key=lambda x: -x["delta_pop"])
    return fluxes


def compute_frame_ibos(
    mol: Any,
    mf: Any,
    symbols: Optional[List[str]] = None,
    prev_data: Optional[Dict[str, Any]] = None,
    atom_labels: Optional[List[str]] = None,
    non_pi_occ_coeff: Optional[np.ndarray] = None,
    non_pi_virt_coeff: Optional[np.ndarray] = None,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, np.ndarray], np.ndarray, np.ndarray]:
    """
    Thin adapter delegating occupied + valence virtual (VVO) localization to the
    unified ValenceLocalizer shared with server.py (tools/qm/analyzers/valence.py).
    """
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    full_mo_coeff = mf.mo_coeff if mol.spin == 0 else mf.mo_coeff[0]
    full_mo_energy = mf.mo_energy if mol.spin == 0 else mf.mo_energy[0]
    occ_mo = full_mo_coeff[:, :nocc]

    n_core = count_core_orbitals(mol, full_mo_energy, nocc)
    core_mo = occ_mo[:, :n_core].copy()
    val_mo = non_pi_occ_coeff.copy() if non_pi_occ_coeff is not None else occ_mo[:, n_core:].copy()
    virt_mo = (
        non_pi_virt_coeff.copy()
        if non_pi_virt_coeff is not None
        else full_mo_coeff[:, nocc:].copy()
    )

    ovlp = mol.intor_symmetric("int1e_ovlp")
    raw_fock = mf.get_fock()
    fock = (
        raw_fock[0]
        if (mol.spin > 0 and isinstance(raw_fock, (tuple, list, np.ndarray)) and np.ndim(raw_fock) == 3)
        else raw_fock
    )
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]
    atom_ids = (
        list(atom_labels)
        if (atom_labels and len(atom_labels) == mol.natm)
        else [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
    )

    ev_s, U_s = np.linalg.eigh(ovlp)
    S_half = U_s @ np.diag(np.sqrt(np.maximum(ev_s, 1e-12))) @ U_s.T

    def _lowdin_pop(vec: np.ndarray, scale: float) -> np.ndarray:
        v_low = S_half @ vec
        sq = v_low ** 2
        ap = np.zeros(mol.natm)
        for mu in range(mol.nao):
            ap[ao_atom[mu]] += sq[mu]
        s_tot = float(ap.sum())
        if s_tot > 1e-10:
            ap = (ap / s_tot) * scale
        return ap

    # Delegate localization to unified ValenceLocalizer (shared with server.py)
    localizer = ValenceLocalizer(
        mol,
        full_mo_coeff,
        fock=fock,
        atom_labels=atom_ids,
        occ_all_coeff=occ_mo,
    )
    localized_orbitals = localizer.localize(val_mo, virt_mo)

    occ_orbs = []
    virt_orbs = []
    for orb in localized_orbitals:
        if orb.type in ("sigma_star", "pi_star", "p_star") or "star" in orb.type:
            virt_orbs.append(orb)
        else:
            occ_orbs.append(orb)

    ibo_val = (
        np.column_stack([orb.mo_coeff for orb in occ_orbs])
        if occ_orbs
        else np.zeros((mol.nao, 0))
    )
    n_val = ibo_val.shape[1]

    # Optional backward-compatible ordering against prev_data if directly supplied
    if prev_data and "mol" in prev_data and "ibo_coeff" in prev_data:
        n_prev_core = prev_data.get("n_core", core_mo.shape[1])
        if (
            prev_data["mol"].natm == mol.natm
            and core_mo.shape[1] == n_prev_core
            and prev_data["ibo_coeff"].shape[1] - n_prev_core == n_val
        ):
            s_cross = gto.intor_cross("int1e_ovlp", prev_data["mol"], mol)
            if n_prev_core > 0:
                prev_core = prev_data["ibo_coeff"][:, :n_prev_core]
                ovlp_core = prev_core.T @ s_cross @ core_mo
                _, cc = linear_sum_assignment(-np.abs(ovlp_core))
                core_mo = core_mo[:, cc]
                for k in range(core_mo.shape[1]):
                    if ovlp_core[k, cc[k]] < 0.0:
                        core_mo[:, k] = -core_mo[:, k]
            if n_val > 0:
                prev_val = prev_data["ibo_coeff"][:, n_prev_core:]
                ovlp_val = prev_val.T @ s_cross @ ibo_val
                _, c = linear_sum_assignment(-np.abs(ovlp_val))
                ibo_val = ibo_val[:, c]
                occ_orbs = [occ_orbs[idx] for idx in c]
                for k in range(n_val):
                    if float(prev_val[:, k].T @ s_cross @ ibo_val[:, k]) < 0.0:
                        ibo_val[:, k] = -ibo_val[:, k]
                        occ_orbs[k].mo_coeff = ibo_val[:, k]

    ibo_coeff = np.hstack([core_mo, ibo_val]) if (core_mo.shape[1] > 0 or n_val > 0) else np.zeros((mol.nao, 0))
    n_c = core_mo.shape[1]
    n_ibo_tot = ibo_coeff.shape[1]
    pops = np.zeros((n_ibo_tot, mol.natm))

    ibo_dict: Dict[str, Dict[str, Any]] = {}
    star_dict: Dict[str, Dict[str, Any]] = {}
    vec_dict: Dict[str, np.ndarray] = {}

    # 1. Core orbitals
    for k in range(n_c):
        key = f"IBO_{k:02d}"
        vec = core_mo[:, k].copy()
        pop_k = _lowdin_pop(vec, scale=2.0)
        pops[k] = pop_k
        ai = int(np.argmax(pop_k))
        e_ev = float(vec.T @ fock @ vec) * HARTREE_TO_EV
        ibo_dict[key] = {
            "key": key, "name": key, "index": k, "type": "core", "category": "core",
            "label": f"Core({atom_ids[ai]})", "atoms": [atom_ids[ai]],
            "chem_key": f"core_{atom_ids[ai]}_{k}",
            "energy_ev": round(e_ev, 2), "pop": pop_k.round(2).tolist(),
        }
        vec_dict[key] = vec

    # 2. Occupied valence IBOs
    for idx_v, orb in enumerate(occ_orbs):
        k = n_c + idx_v
        key = f"IBO_{k:02d}"
        vec = orb.mo_coeff.copy()
        pop_k = _lowdin_pop(vec, scale=2.0)
        pops[k] = pop_k
        typ = orb.type
        atoms = list(orb.atoms)
        if typ == "lone_pair":
            cat, lbl = "non_bonding", f"LP({atoms[0]})"
            c_key = f"lp_{atoms[0]}"
        elif typ in ("sigma", "delocalized_sigma"):
            cat, lbl = "bonding", f"σ({'-'.join(atoms)})"
            c_key = f"sigma_{'-'.join(atoms)}"
        elif typ in ("pi", "delocalized_pi"):
            cat, lbl = "bonding", f"π({'-'.join(atoms)})"
            c_key = f"pi_{'-'.join(atoms)}"
        else:
            cat, lbl = "other", f"{typ}({'-'.join(atoms)})"
            c_key = f"{typ}_{'-'.join(atoms)}"

        e_ev = orb.energy_ev if orb.energy_ev is not None else round(float(vec.T @ fock @ vec) * HARTREE_TO_EV, 2)
        ibo_dict[key] = {
            "key": key, "name": key, "index": k, "type": typ, "category": cat,
            "label": lbl, "atoms": atoms, "chem_key": c_key,
            "energy_ev": round(float(e_ev), 2), "pop": pop_k.round(2).tolist(),
        }
        vec_dict[key] = vec

    # 3. Orthonormal Valence Virtual Orbitals (VVOs)
    for idx_s, orb in enumerate(virt_orbs):
        star_key = f"STAR_{idx_s:02d}"
        vec = orb.mo_coeff.copy()
        pop_s = _lowdin_pop(vec, scale=1.0)
        v_type = orb.type
        v_atoms = list(orb.atoms)
        if v_type == "p_star":
            lbl = f"p*({v_atoms[0]})" if v_atoms else "p*"
            c_key = f"p_star_{v_atoms[0]}" if v_atoms else f"p_star_{idx_s}"
        elif v_type == "pi_star":
            lbl = f"π*({'-'.join(v_atoms)})"
            c_key = f"pi_star_{'-'.join(v_atoms)}"
        else:
            lbl = f"σ*({'-'.join(v_atoms)})"
            c_key = f"sigma_star_{'-'.join(v_atoms)}"

        e_star = orb.energy_ev if orb.energy_ev is not None else round(float(vec.T @ fock @ vec) * HARTREE_TO_EV, 2)
        star_dict[star_key] = {
            "key": star_key, "name": star_key, "index": 100 + idx_s, "type": v_type,
            "category": "antibonding", "label": lbl, "atoms": v_atoms,
            "chem_key": c_key, "energy_ev": round(float(e_star), 2),
            "pop": pop_s.round(2).tolist(),
        }
        vec_dict[star_key] = vec

    return ibo_dict, star_dict, vec_dict, pops, ibo_coeff
