"""
Frontier MOs & Step-0 Interface Module.
"""

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from pyscf import gto

from qm.constants import HARTREE_TO_EV
from .subspace import extract_step0_planar_pi


def extract_frontier_orbitals(
    mol: Any,
    mf: Any,
    prev_mol: Optional[Any] = None,
    prev_mo_coeffs: Optional[Dict[str, np.ndarray]] = None,
    ibos_m: Optional[Dict[str, Any]] = None,
    star_m: Optional[Dict[str, Any]] = None,
    ibos_v: Optional[Dict[str, np.ndarray]] = None,
    atom_labels: Optional[List[str]] = None,
    step0_pi: Optional[Tuple[Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]] = None,
) -> Tuple[Dict[str, Dict[str, Any]], Dict[str, Dict[str, Any]], Dict[str, np.ndarray]]:
    """Extract canonical frontier MOs and Step-0 planar conjugated pi systems."""
    nocc = mol.nelectron // 2 if mol.spin == 0 else mol.nelec[0]
    mo_coeff = mf.mo_coeff[0] if mol.spin > 0 else mf.mo_coeff
    mo_energy = mf.mo_energy[0] if mol.spin > 0 else mf.mo_energy

    all_atom_ids = (
        list(atom_labels)
        if (atom_labels and len(atom_labels) == mol.natm)
        else [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]
    )

    ovlp = mol.intor_symmetric("int1e_ovlp")
    ao_atom = [int(l.split()[0]) for l in mol.ao_labels()]
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

    orbitals_meta: Dict[str, Dict[str, Any]] = {}
    frontier_vectors: Dict[str, np.ndarray] = {}

    for name, mo_idx, label in target_keys:
        if 0 <= mo_idx < mo_coeff.shape[1]:
            vec = mo_coeff[:, mo_idx].copy()
            if s_cross is not None and prev_mo_coeffs and name in prev_mo_coeffs:
                if float(prev_mo_coeffs[name] @ s_cross @ vec) < 0.0:
                    vec = -vec

            ps = np.outer(vec, vec) @ ovlp
            ap = np.zeros(mol.natm)
            for mu in range(mol.nao):
                ap[ao_atom[mu]] += ps[mu, mu]
            dom_atoms = [all_atom_ids[i] for i in np.argsort(-ap) if ap[i] >= 0.12][:4]

            frontier_vectors[name] = vec
            e_ev = float(mo_energy[mo_idx] * HARTREE_TO_EV)
            orbitals_meta[name] = {
                "name": name,
                "label": label,
                "category": "frontier",
                "type": "frontier",
                "mo_index": int(mo_idx),
                "energy_ev": round(e_ev, 3),
                "is_occupied": bool(mo_idx < nocc),
                "atoms": dom_atoms,
            }

    if step0_pi is not None:
        canonical_pi_meta, canonical_pi_vecs = step0_pi
    else:
        _, _, canonical_pi_meta, canonical_pi_vecs = extract_step0_planar_pi(mol, mf, atom_labels=atom_labels)

    if s_cross is not None and prev_mo_coeffs:
        for k, v in list(canonical_pi_vecs.items()):
            if k in prev_mo_coeffs and float(prev_mo_coeffs[k] @ s_cross @ v) < 0.0:
                canonical_pi_vecs[k] = -v

    all_vecs = {**frontier_vectors, **canonical_pi_vecs}
    return orbitals_meta, canonical_pi_meta, all_vecs

