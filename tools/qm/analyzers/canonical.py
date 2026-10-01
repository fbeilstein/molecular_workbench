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
        
        # Check planarity: if the ring is significantly puckered (e.g. tub-shaped cyclooctatetraene),
        # a single global plane cannot define the pi system.
        distances = np.abs(np.dot(centered, normal))
        if np.max(distances) > 0.5:  # 0.5 Bohr ~ 0.26 Angstroms
            return [], []
        
        from pyscf.data.elements import chemcore
        n_core = chemcore(self.mol)
        
        pi_occ = []
        pi_virt = []
        
        # Analyze each MO (valence only, skip core)
        for idx in range(n_core, self.mo_coeff.shape[1]):
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
        Uses an anti-symmetry test: a true pi orbital must have opposite phases
        above and below the molecular plane.
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
            
        # 2. Check for anti-symmetry across the ring plane
        coords = self.mol.atom_coords()[ring_indices]
        
        # Sample points above and below the plane (e.g. 0.5 Bohr)
        d = 0.5
        coords_up = coords + d * normal
        coords_dn = coords - d * normal
        
        ao_up = self.mol.eval_gto("GTOval", coords_up)
        ao_dn = self.mol.eval_gto("GTOval", coords_dn)
        
        mo_up = np.dot(ao_up, mo)
        mo_dn = np.dot(ao_dn, mo)
        
        # We need significant amplitude above/below the plane
        max_abs = max(np.max(np.abs(mo_up)), np.max(np.abs(mo_dn)))
        if max_abs < 0.05:
            return False
            
        # For significant samples, check if they are anti-symmetric
        n_significant = 0
        n_anti_symmetric = 0
        
        for u, d in zip(mo_up, mo_dn):
            mag_u, mag_d = abs(u), abs(d)
            max_mag = max(mag_u, mag_d)
            if max_mag > max_abs * 0.1:
                n_significant += 1
                
                # Must be reasonably symmetric in magnitude to be a true canonical pi orbital.
                # If highly asymmetric, it's just an arbitrary node in an asymmetric molecule.
                if u * d < 0 and abs(mag_u - mag_d) / max_mag < 0.45:
                    n_anti_symmetric += 1
                    
        if n_significant > 0 and (n_anti_symmetric / n_significant) > 0.8:
            return True
            
        return False
