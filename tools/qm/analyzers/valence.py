from typing import Optional
import numpy as np
from pyscf import gto, lo
from qm.classification import _atom_populations, _classify_orbital
from qm.models import Orbital

HARTREE_TO_EV = 27.2114


class ValenceLocalizer:
    def __init__(self, mol: gto.M, mo_coeff: np.ndarray, fock: Optional[np.ndarray] = None):
        self.mol = mol
        self.mo_coeff = mo_coeff
        self.fock = fock
        self.ovlp = mol.intor_symmetric('int1e_ovlp')
        self.atom_labels = [mol.atom_symbol(i) for i in range(mol.natm)]
        # For simplicity, we assume atom_ids are just the symbol + index (e.g. C0, H1)
        self.atom_ids = [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]

    def _eval_energy_ev(self, vec: np.ndarray) -> Optional[float]:
        if self.fock is None:
            return None
        denom = float(vec.T @ self.ovlp @ vec)
        if denom < 1e-8:
            return None
        return round(float((vec.T @ self.fock @ vec) / denom) * HARTREE_TO_EV, 2)

    def localize(self, cols_to_localize: list, virt_cols: list):
        """
        Localize the given columns (occupied space) and return a list of Orbital objects.
        Automatically constructs textbook sigma* orbitals by antisymmetrizing the localized sigma bonds 
        and projecting them onto the virtual space.
        """
        orbitals = []
        if cols_to_localize is None:
            return orbitals
        if isinstance(cols_to_localize, (list, tuple)) and len(cols_to_localize) == 0:
            return orbitals
        if isinstance(cols_to_localize, np.ndarray) and cols_to_localize.size == 0:
            return orbitals

        if isinstance(cols_to_localize, np.ndarray) and cols_to_localize.ndim == 2:
            occ_coeff = cols_to_localize
        else:
            occ_coeff = self.mo_coeff[:, cols_to_localize]

        if virt_cols is not None:
            if isinstance(virt_cols, np.ndarray) and virt_cols.ndim == 2:
                virt_coeff = virt_cols
            elif isinstance(virt_cols, (list, tuple)) and len(virt_cols) > 0:
                virt_coeff = self.mo_coeff[:, virt_cols]
            else:
                virt_coeff = None
        else:
            virt_coeff = None
        
        ao_labels = self.mol.ao_labels(fmt=False)
        
        try:
            # IAO/IBO localization
            iaos = lo.iao.iao(self.mol, occ_coeff)
            loc_occ = lo.ibo.ibo(self.mol, occ_coeff, iaos=iaos)
            
            # Track pi bonds to prevent "triple bonds" to terminal oxygens
            seen_pi = set()
            
            for i in range(loc_occ.shape[1]):
                mo = loc_occ[:, i]
                pop = _atom_populations(self.mol, mo, self.ovlp)
                info = _classify_orbital(self.mol, mo, pop, self.atom_labels, self.atom_ids, self.ovlp)
                
                orb_type = info.get('type', 'unknown')
                atoms = info.get('atoms', [info.get('atom', '?')])
                
                # Heuristic: Limit hypervalent M=O bonds to exactly 1 pi bond.
                # If a second polarized lone pair crosses the threshold, demote it to a lone pair.
                if orb_type == 'pi' and len(atoms) == 2:
                    pair = tuple(sorted(atoms))
                    s1, s2 = ''.join(filter(str.isalpha, pair[0])), ''.join(filter(str.isalpha, pair[1]))
                    hyper = {'S', 'P', 'Cl', 'Br', 'I', 'N'}
                    if (s1 in hyper and s2 == 'O') or (s2 in hyper and s1 == 'O'):
                        if pair in seen_pi:
                            orb_type = 'lone_pair'
                            atoms = [pair[0] if s1 == 'O' else pair[1]]
                        else:
                            seen_pi.add(pair)
                
                orbitals.append(Orbital(
                    type=orb_type,
                    atoms=atoms,
                    mo_coeff=mo,
                    energy_ev=self._eval_energy_ev(mo),
                ))
                
                # If this is a 2-center bond and we have virtual space, construct the antibond
                if virt_coeff is not None and orb_type in ('sigma', 'delocalized_sigma', 'pi', 'delocalized_pi') and len(atoms) >= 2:
                    a1 = atoms[0]
                    a2 = atoms[1]
                    idx1 = self.atom_ids.index(a1)
                    idx2 = self.atom_ids.index(a2)
                    
                    c_A = np.zeros_like(mo)
                    c_B = np.zeros_like(mo)
                    
                    for mu, (atom_idx, *_) in enumerate(ao_labels):
                        if atom_idx == idx1:
                            c_A[mu] = mo[mu]
                        elif atom_idx == idx2:
                            c_B[mu] = mo[mu]
                            
                    # True textbook sigma* is the antisymmetrized combination
                    c_anti = c_A - c_B
                    
                    # Project exactly onto the virtual space
                    c_virt_basis = virt_coeff.T @ self.ovlp @ c_anti
                    c_star = virt_coeff @ c_virt_basis
                    
                    norm = np.sqrt(c_star.T @ self.ovlp @ c_star)
                    if norm > 1e-3:
                        c_star /= norm
                        star_type = 'sigma_star' if 'sigma' in orb_type else 'pi_star'
                        orbitals.append(Orbital(
                            type=star_type, 
                            atoms=[a1, a2], 
                            mo_coeff=c_star, 
                            energy_ev=self._eval_energy_ev(c_star),
                            sym_score=1.0  # Perfect by definition
                        ))

            # Intra-atomic Lone-Pair Subspace Fock Purification ("No Rabbit Ears" / NBO-like natural hybrids):
            # Diagonalize the projected Fock operator strictly inside each single atom's localized LP subspace
            # (m >= 2) to eliminate arbitrary s/p mixing without delocalizing onto neighboring bonds.
            if self.fock is not None:
                lp_groups = {}
                for idx_orb, orb in enumerate(orbitals):
                    if orb.type == 'lone_pair' and orb.atoms:
                        lp_groups.setdefault(orb.atoms[0], []).append(idx_orb)

                for _, orb_indices in lp_groups.items():
                    if len(orb_indices) < 2:
                        continue
                    C_sub = np.column_stack([orbitals[k].mo_coeff for k in orb_indices])
                    S_sub = C_sub.T @ self.ovlp @ C_sub
                    ev_s, U_s = np.linalg.eigh(S_sub)
                    C_orth = C_sub @ (U_s @ np.diag(1.0 / np.sqrt(np.maximum(ev_s, 1e-10))) @ U_s.T)
                    F_sub = C_orth.T @ self.fock @ C_orth
                    ev_f, U_f = np.linalg.eigh(F_sub)
                    C_pur = C_orth @ U_f
                    for j, k in enumerate(orb_indices):
                        orbitals[k].mo_coeff = C_pur[:, j].copy()
                        orbitals[k].energy_ev = round(float(ev_f[j] * HARTREE_TO_EV), 2)
                
        except Exception as e:
            print(f"Valence localization failed: {e}")

        return orbitals
