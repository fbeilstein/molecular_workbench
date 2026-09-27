import numpy as np
from pyscf import gto, lo
from qm.classification import _atom_populations, _classify_orbital
from qm.models import Orbital

class ValenceLocalizer:
    def __init__(self, mol: gto.M, mo_coeff: np.ndarray):
        self.mol = mol
        self.mo_coeff = mo_coeff
        self.ovlp = mol.intor_symmetric('int1e_ovlp')
        self.atom_labels = [mol.atom_symbol(i) for i in range(mol.natm)]
        # For simplicity, we assume atom_ids are just the symbol + index (e.g. C0, H1)
        self.atom_ids = [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]

    def localize(self, cols_to_localize: list, virt_cols: list):
        """
        Localize the given columns (occupied space) and return a list of Orbital objects.
        Automatically constructs textbook sigma* orbitals by antisymmetrizing the localized sigma bonds 
        and projecting them onto the virtual space.
        """
        orbitals = []
        if not cols_to_localize:
            return orbitals

        occ_coeff = self.mo_coeff[:, cols_to_localize]
        virt_coeff = self.mo_coeff[:, virt_cols] if virt_cols else None
        
        ao_labels = self.mol.ao_labels(fmt=False)
        
        try:
            # IAO/IBO localization
            iaos = lo.iao.iao(self.mol, occ_coeff)
            loc_occ = lo.ibo.ibo(self.mol, occ_coeff, iaos=iaos)
            
            for i in range(loc_occ.shape[1]):
                mo = loc_occ[:, i]
                pop = _atom_populations(self.mol, mo, self.ovlp)
                info = _classify_orbital(self.mol, mo, pop, self.atom_labels, self.atom_ids, self.ovlp)
                
                orb_type = info.get('type', 'unknown')
                atoms = info.get('atoms', [info.get('atom', '?')])
                
                orbitals.append(Orbital(type=orb_type, atoms=atoms, mo_coeff=mo))
                
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
                            sym_score=1.0  # Perfect by definition
                        ))
                
        except Exception as e:
            print(f"Valence localization failed: {e}")

        return orbitals
