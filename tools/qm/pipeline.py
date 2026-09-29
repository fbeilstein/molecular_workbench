from pyscf import gto
from qm.analyzers.topology import TopologyDetector
from qm.analyzers.canonical import CanonicalAnalyzer
from qm.analyzers.valence import ValenceLocalizer
from qm.renderer import Renderer
from qm.models import Orbital

class QuantumPipeline:
    def __init__(self, mol: gto.M, mf, smiles: str, name: str, out_dir: str, grid_points: int = 50):
        self.mol = mol
        self.mf = mf
        self.smiles = smiles
        self.name = name
        self.out_dir = out_dir
        self.grid_points = grid_points
        self.orbitals = []

    def run(self):
        print(f"--- Running Refactored Quantum Pipeline for {self.name} ---")
        
        # 1. Topology (find aromatic rings and conjugated systems)
        topology = TopologyDetector(self.smiles)
        conjugated_fragments = topology.get_conjugated_fragments(min_size=3)
        
        # 2. Canonical Analyzer
        canonical = CanonicalAnalyzer(self.mol, self.mf.mo_coeff, self.mf.mo_energy)
        homo_idx, lumo_idx = canonical.get_frontier_orbitals()
        
        self.orbitals.append(Orbital(type='homo', atoms=[], mo_coeff=self.mf.mo_coeff[:, homo_idx], energy_ev=self.mf.mo_energy[homo_idx]*27.211))
        self.orbitals.append(Orbital(type='lumo', atoms=[], mo_coeff=self.mf.mo_coeff[:, lumo_idx], energy_ev=self.mf.mo_energy[lumo_idx]*27.211))
        
        all_pi_occ = []
        all_pi_virt = []
        
        for fragment_atoms in conjugated_fragments:
            print(f"  Detected conjugated system: {len(fragment_atoms)} atoms")
            pi_occ, pi_virt = canonical.get_pi_system(fragment_atoms)
            
            # In larger basis sets, many high-energy virtual pi orbitals appear.
            valence_pi_virt = pi_virt[:len(pi_occ)]

            self._snap_degenerate_pairs(pi_occ, fragment_atoms)
            self._snap_degenerate_pairs(valence_pi_virt, fragment_atoms)
            
            atom_ids = [f"{self.mol.atom_symbol(i)}{i}" for i in fragment_atoms]
            
            for idx in pi_occ:
                energy_ev = self.mf.mo_energy[idx] * 27.211
                depth = homo_idx - idx
                lbl = f"HOMO-{depth}" if depth > 0 else "HOMO"
                self.orbitals.append(Orbital(
                    type='delocalized_pi', atoms=atom_ids, mo_coeff=self.mf.mo_coeff[:, idx], 
                    energy_ev=energy_ev, canonical_label=f"π ({lbl})"
                ))
                all_pi_occ.append(idx)
                
            for idx in valence_pi_virt:
                energy_ev = self.mf.mo_energy[idx] * 27.211
                depth = idx - lumo_idx
                lbl = f"LUMO+{depth}" if depth > 0 else "LUMO"
                self.orbitals.append(Orbital(
                    type='delocalized_pistar', atoms=atom_ids, mo_coeff=self.mf.mo_coeff[:, idx], 
                    energy_ev=energy_ev, canonical_label=f"π* ({lbl})"
                ))
                all_pi_virt.append(idx)

        # 3. Valence Localizer (and Sigma Star projection)
        from pyscf.data.elements import chemcore
        n_core = chemcore(self.mol)
        n_occ = self.mol.nelectron // 2
        non_pi_occ = [i for i in range(n_core, n_occ) if i not in all_pi_occ]
        non_pi_virt = [i for i in range(n_occ, self.mf.mo_coeff.shape[1]) if i not in all_pi_virt]
        
        valence_loc = ValenceLocalizer(self.mol, self.mf.mo_coeff)
        self.orbitals.extend(valence_loc.localize(non_pi_occ, non_pi_virt))

        # 5. Render
        renderer = Renderer(self.mol, self.name, self.out_dir, self.grid_points)
        manifest = renderer.render(self.orbitals)
        
        return manifest

    def _snap_degenerate_pairs(self, indices, aromatic_atoms):
        import numpy as np
        if not indices or not aromatic_atoms:
            return
        ovlp = self.mol.intor_symmetric('int1e_ovlp')
        ao_labels = self.mol.ao_labels(fmt=False)
        target_atom = aromatic_atoms[0]
        
        # Build population operator for the target atom
        P_A = np.zeros_like(ovlp)
        for mu, (atom_idx, *_) in enumerate(ao_labels):
            if atom_idx == target_atom:
                P_A[mu, :] += 0.5 * ovlp[mu, :]
                P_A[:, mu] += 0.5 * ovlp[:, mu]
                
        # Group indices into degenerate sets (within 0.01 Hartree)
        i = 0
        while i < len(indices) - 1:
            if abs(self.mf.mo_energy[indices[i]] - self.mf.mo_energy[indices[i+1]]) < 1e-4:
                idx1, idx2 = indices[i], indices[i+1]
                C = self.mf.mo_coeff[:, [idx1, idx2]]
                P_red = C.T @ P_A @ C
                evals, evecs = np.linalg.eigh(P_red)
                self.mf.mo_coeff[:, [idx1, idx2]] = C @ evecs
                i += 2
            else:
                i += 1
