import numpy as np
from pyscf import gto, lo
from qm.classification import _atom_populations, _classify_orbital
from qm.models import Orbital

class VirtualLocalizer:
    def __init__(self, mol: gto.M, mo_coeff: np.ndarray):
        self.mol = mol
        self.mo_coeff = mo_coeff
        self.ovlp = mol.intor_symmetric('int1e_ovlp')
        self.atom_labels = [mol.atom_symbol(i) for i in range(mol.natm)]
        self.atom_ids = [f"{mol.atom_symbol(i)}{i}" for i in range(mol.natm)]

    def localize(self, cols_to_localize: list):
        """
        Localize the virtual space and return filtered sigma_star orbitals.
        Uses Pipek-Mezey localization to prevent diffuse function mixing.
        """
        orbitals = []
        if not cols_to_localize:
            return orbitals

        virt_coeff = self.mo_coeff[:, cols_to_localize]
        
        try:
            pm = lo.PipekMezey(self.mol, virt_coeff)
            pm.init_guess = 'random'
            virt_loc = pm.kernel()
            
            best_sigmastar = {}
            for i in range(virt_loc.shape[1]):
                mo = virt_loc[:, i]
                pop = _atom_populations(self.mol, mo, self.ovlp)
                info = _classify_orbital(self.mol, mo, pop, self.atom_labels, self.atom_ids, self.ovlp)
                
                atoms = info.get('atoms', [])
                if len(atoms) == 2:
                    a1, a2 = atoms[0], atoms[1]
                    pair = tuple(sorted([a1, a2]))
                    
                    # Calculate symmetry score (how close to 50/50 population)
                    abs_pop = np.abs(pop)
                    frac = abs_pop / abs_pop.sum()
                    idx1 = self.atom_ids.index(a1)
                    idx2 = self.atom_ids.index(a2)
                    sym_score = abs(frac[idx1] - frac[idx2])

                    if pair not in best_sigmastar or sym_score < best_sigmastar[pair]['sym_score']:
                        best_sigmastar[pair] = {
                            'mo': mo,
                            'sym_score': sym_score,
                            'atoms': [a1, a2]
                        }
                        
            for data in best_sigmastar.values():
                orbitals.append(Orbital(
                    type="sigma_star",
                    atoms=data['atoms'],
                    mo_coeff=data['mo'],
                    sym_score=data['sym_score']
                ))
                
        except Exception as e:
            print(f"Virtual localization failed: {e}")

        return orbitals
