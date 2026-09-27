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

    def localize(self, cols_to_localize: list):
        """
        Localize the given columns (occupied space) and return a list of Orbital objects.
        Uses Knizia's IBO method.
        """
        orbitals = []
        if not cols_to_localize:
            return orbitals

        occ_coeff = self.mo_coeff[:, cols_to_localize]
        
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
                
        except Exception as e:
            print(f"Valence localization failed: {e}")

        return orbitals
