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
        """Check if an MO is a pi orbital centered on the ring."""
        dm = np.outer(mo, mo)
        ps = dm * self.ovlp
        
        total_ring = 0.0
        p_perp = 0.0
        
        ao_labels = self.mol.ao_labels(fmt=False)
        
        for mu, (atom_idx, _, ao_type, axis) in enumerate(ao_labels):
            if atom_idx in ring_indices:
                pop = ps[mu, :].sum()
                total_ring += pop
                
                # If it's a p-orbital, project onto the normal vector
                if 'p' in ao_type:
                    vec = np.zeros(3)
                    if axis == 'x': vec[0] = 1.0
                    elif axis == 'y': vec[1] = 1.0
                    elif axis == 'z': vec[2] = 1.0
                    
                    # Contribution perpendicular to the ring plane
                    proj = abs(np.dot(vec, normal))
                    p_perp += pop * proj

        # To be a ring pi orbital, it must have significant density on the ring
        if total_ring < 0.2:
            return False
            
        # And that density must be predominantly perpendicular p-character
        if p_perp > total_ring * 0.75:
            return True
            
        return False
