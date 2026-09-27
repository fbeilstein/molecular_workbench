import os
import json
from pyscf.tools import cubegen
from qm.models import Orbital
from typing import List, Dict

class Renderer:
    def __init__(self, mol, name: str, out_dir: str, grid_points: int = 50):
        self.mol = mol
        self.name = name
        self.out_dir = out_dir
        self.grid_points = grid_points
        self.manifest = {
            'canonical': {},
            'localized': {
                'sigma': [],
                'pi': [],
                'lone_pairs': [],
                'sigma_star': [],
                'pi_star': []
            }
        }

    def render(self, orbitals: List[Orbital]):
        for idx, orb in enumerate(orbitals):
            if orb.type in ['core']:
                continue
                
            atoms_str = "_".join(orb.atoms)
            
            if orb.type in ('homo', 'lumo'):
                cube_name = f"{self.name}_{orb.type}.cube"
                self.manifest['canonical'][orb.type] = orb.to_dict(cube_name)
                
            elif orb.type == 'pi_canonical':
                cube_name = f"{self.name}_pi_canonical_{idx}.cube"
                if 'pi_system' not in self.manifest['canonical']:
                    self.manifest['canonical']['pi_system'] = []
                d = orb.to_dict(cube_name)
                d['homo_label'] = orb.canonical_label
                self.manifest['canonical']['pi_system'].append(d)
                
            elif orb.type == 'pi_star_canonical':
                cube_name = f"{self.name}_pistar_canonical_{idx}.cube"
                self.manifest['localized']['pi_star'].append(orb.to_dict(cube_name))
                
            elif orb.type in ('sigma', 'delocalized_sigma'):
                cube_name = f"{self.name}_sigma_{atoms_str}_{idx}.cube"
                self.manifest['localized']['sigma'].append(orb.to_dict(cube_name))
                
            elif orb.type == 'sigma_star':
                cube_name = f"{self.name}_sigmastar_{atoms_str}_{idx}.cube"
                self.manifest['localized']['sigma_star'].append(orb.to_dict(cube_name))
                
            elif orb.type == 'lone_pair':
                cube_name = f"{self.name}_lp_{atoms_str}_{idx}.cube"
                self.manifest['localized']['lone_pairs'].append(orb.to_dict(cube_name))
                
            elif orb.type in ('pi', 'delocalized_pi'):
                cube_name = f"{self.name}_pi_{atoms_str}_{idx}.cube"
                self.manifest['localized']['pi'].append(orb.to_dict(cube_name))
                
            elif orb.type == 'pi_star':
                cube_name = f"{self.name}_pistar_{atoms_str}_{idx}.cube"
                self.manifest['localized']['pi_star'].append(orb.to_dict(cube_name))
            
            else:
                cube_name = f"{self.name}_{orb.type}_{atoms_str}_{idx}.cube"
                
            cube_path = os.path.join(self.out_dir, cube_name)
            try:
                cubegen.orbital(self.mol, cube_path, orb.mo_coeff, 
                                nx=self.grid_points, ny=self.grid_points, nz=self.grid_points, margin=5.0)
                print(f"  ✓ Rendered {cube_name}")
            except Exception as e:
                print(f"  ⚠ Failed to render {cube_name}: {e}")

        return self.manifest
