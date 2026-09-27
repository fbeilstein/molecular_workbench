from pyscf import gto
from qm.analyzers.topology import TopologyDetector
from qm.analyzers.canonical import CanonicalAnalyzer
from qm.analyzers.valence import ValenceLocalizer
from qm.analyzers.virtual import VirtualLocalizer
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
        
        # 1. Topology (find aromatic rings)
        topology = TopologyDetector(self.smiles)
        aromatic_atoms = topology.get_aromatic_atoms()
        
        # 2. Canonical Analyzer
        canonical = CanonicalAnalyzer(self.mol, self.mf.mo_coeff, self.mf.mo_energy)
        homo_idx, lumo_idx = canonical.get_frontier_orbitals()
        
        self.orbitals.append(Orbital(type='homo', atoms=[], mo_coeff=self.mf.mo_coeff[:, homo_idx], energy_ev=self.mf.mo_energy[homo_idx]*27.211))
        self.orbitals.append(Orbital(type='lumo', atoms=[], mo_coeff=self.mf.mo_coeff[:, lumo_idx], energy_ev=self.mf.mo_energy[lumo_idx]*27.211))
        
        pi_occ, pi_virt = canonical.get_pi_system(aromatic_atoms)
        
        for idx in pi_occ:
            energy_ev = self.mf.mo_energy[idx] * 27.211
            depth = homo_idx - idx
            lbl = f"HOMO-{depth}" if depth > 0 else "HOMO"
            self.orbitals.append(Orbital(
                type='pi_canonical', atoms=['ring'], mo_coeff=self.mf.mo_coeff[:, idx], 
                energy_ev=energy_ev, canonical_label=f"π ({lbl})"
            ))
            
        # In larger basis sets, many high-energy virtual pi orbitals appear.
        # We only want the core valence pi* orbitals, so we take a number equal to pi_occ
        valence_pi_virt = pi_virt[:len(pi_occ)]
        
        for idx in valence_pi_virt:
            energy_ev = self.mf.mo_energy[idx] * 27.211
            depth = idx - lumo_idx
            lbl = f"LUMO+{depth}" if depth > 0 else "LUMO"
            self.orbitals.append(Orbital(
                type='pi_star_canonical', atoms=['ring'], mo_coeff=self.mf.mo_coeff[:, idx], 
                energy_ev=energy_ev, canonical_label=f"π* ({lbl})"
            ))

        # 3. Valence Localizer
        n_occ = self.mol.nelectron // 2
        non_pi_occ = [i for i in range(n_occ) if i not in pi_occ]
        valence_loc = ValenceLocalizer(self.mol, self.mf.mo_coeff)
        self.orbitals.extend(valence_loc.localize(non_pi_occ))

        # 4. Virtual Localizer
        non_pi_virt = [i for i in range(n_occ, self.mf.mo_coeff.shape[1]) if i not in pi_virt]
        virtual_loc = VirtualLocalizer(self.mol, self.mf.mo_coeff)
        self.orbitals.extend(virtual_loc.localize(non_pi_virt))

        # 5. Render
        renderer = Renderer(self.mol, self.name, self.out_dir, self.grid_points)
        manifest = renderer.render(self.orbitals)
        
        return manifest
