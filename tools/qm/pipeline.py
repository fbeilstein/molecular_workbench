from pyscf import gto
from qm.analyzers.planar_pi import extract_step0_planar_pi
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
        print(f"--- Running Unified Quantum Pipeline for {self.name} ---")
        
        all_atom_ids = [f"{self.mol.atom_symbol(i)}{i}" for i in range(self.mol.natm)]
        
        # 1. Step-0 Planar Conjugated π Subspace Separation (unified engine shared with server_irc.py)
        C_non_pi_occ, C_non_pi_virt, pi_meta, pi_vecs = extract_step0_planar_pi(
            self.mol, self.mf, atom_labels=all_atom_ids, smiles=self.smiles
        )
        
        # 2. Frontier Orbitals (HOMO, LUMO)
        nocc = self.mol.nelectron // 2 if self.mol.spin == 0 else self.mol.nelec[0]
        homo_idx = nocc - 1
        lumo_idx = nocc
        self.orbitals.append(Orbital(type='homo', atoms=[], mo_coeff=self.mf.mo_coeff[:, homo_idx], energy_ev=self.mf.mo_energy[homo_idx]*27.211))
        self.orbitals.append(Orbital(type='lumo', atoms=[], mo_coeff=self.mf.mo_coeff[:, lumo_idx], energy_ev=self.mf.mo_energy[lumo_idx]*27.211))
        
        # 3. Add extracted delocalized π and π* orbitals
        for key, meta in pi_meta.items():
            vec = pi_vecs[key]
            orb_type = 'delocalized_pi' if meta.get('type') == 'delocalized_pi' else 'delocalized_pistar'
            self.orbitals.append(Orbital(
                type=orb_type,
                atoms=meta.get('atoms', []),
                mo_coeff=vec,
                energy_ev=meta.get('energy_ev'),
                canonical_label=meta.get('label')
            ))

        # 4. Valence Localizer (IBO on non-pi subspaces)
        valence_loc = ValenceLocalizer(self.mol, self.mf.mo_coeff)
        self.orbitals.extend(valence_loc.localize(C_non_pi_occ, C_non_pi_virt))

        # 5. Render
        renderer = Renderer(self.mol, self.name, self.out_dir, self.grid_points)
        manifest = renderer.render(self.orbitals)
        
        return manifest
