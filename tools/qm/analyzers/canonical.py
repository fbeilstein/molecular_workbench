import numpy as np
from pyscf import gto

class CanonicalAnalyzer:
    def __init__(self, mol: gto.M, mo_coeff: np.ndarray, mo_energy: np.ndarray):
        self.mol = mol
        self.mo_coeff = mo_coeff
        self.mo_energy = mo_energy
        self.n_occ = mol.nelectron // 2  # Assuming closed shell for simplicity for now
        self.ovlp = mol.intor_symmetric('int1e_ovlp')

    def get_frontier_orbitals(self):
        """Returns the HOMO and LUMO indices."""
        homo_idx = self.n_occ - 1
        lumo_idx = self.n_occ
        return homo_idx, lumo_idx

    def get_pi_system(self, aromatic_atom_indices):
        """
        Identify canonical pi and pi* orbitals for the given aromatic ring atoms.
        """
        if not aromatic_atom_indices:
            return [], []

        # Find the geometric plane of the aromatic atoms
        coords = self.mol.atom_coords()[aromatic_atom_indices]
        center = np.mean(coords, axis=0)
        centered = coords - center
        cov = np.dot(centered.T, centered)
        evals, evecs = np.linalg.eigh(cov)
        normal = evecs[:, 0]  # The eigenvector corresponding to the smallest eigenvalue is the normal
        
        pi_occ = []
        pi_virt = []
        
        # Analyze each MO
        for idx in range(self.mo_coeff.shape[1]):
            mo = self.mo_coeff[:, idx]
            is_pi = self._is_pi_orbital(mo, aromatic_atom_indices, normal)
            if is_pi:
                if idx < self.n_occ:
                    pi_occ.append(idx)
                else:
                    pi_virt.append(idx)
                    
        return pi_occ, pi_virt

    def _is_pi_orbital(self, mo, ring_indices, normal):
        """
        Check if an MO is a pi orbital centered on the ring.
        Uses a rotation-invariant nodal plane check: a pi orbital must have 
        nodes exactly at the atomic centers of the ring atoms.
        """
        # 1. Check if the orbital has significant density on the ring
        dm = np.outer(mo, mo)
        ps = dm * self.ovlp
        
        total_ring = 0.0
        ao_labels = self.mol.ao_labels(fmt=False)
        for mu, (atom_idx, _, ao_type, axis) in enumerate(ao_labels):
            if atom_idx in ring_indices:
                total_ring += ps[mu, :].sum()
                
        if total_ring < 0.2:
            return False
            
        # 2. Check if the MO has a nodal plane at the ring atoms
        coords = self.mol.atom_coords()[ring_indices]
        ao_values = self.mol.eval_gto("GTOval", coords)
        mo_values = np.dot(ao_values, mo)
        
        # Max absolute value of the MO at any ring atom's nucleus
        max_val_at_nuclei = np.max(np.abs(mo_values))
        
        # Pi orbitals should be extremely close to 0 at the nuclei (<0.05)
        # Sigma orbitals usually have values > 0.1 at the nuclei
        if max_val_at_nuclei < 0.05:
            return True
            
        return False
