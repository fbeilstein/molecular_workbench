import os
import shutil
import numpy as np
from pyscf.tools import cubegen
from qm.models import Orbital
from typing import List, Tuple


class Renderer:
    def __init__(self, mol, name: str, out_dir: str, grid_points: int = 50):
        self.mol = mol
        self.name = name
        self.out_dir = out_dir
        self.grid_points = grid_points
        self.manifest = {
            'canonical': {
                'frontier': [],
            },
            'localized': {
                'sigma': [],
                'pi': [],
                'lone_pairs': [],
                'sigma_star': [],
                'pi_star': []
            }
        }

    def render(self, orbitals: List[Orbital]):
        rendered_cache: List[Tuple[np.ndarray, str]] = []

        for idx, orb in enumerate(orbitals):
            if orb.type in ['core']:
                continue
                
            atoms_str = "_".join(orb.atoms)
            
            if orb.type in ('homo', 'lumo', 'frontier'):
                fkey = orb.frontier_key or orb.type
                cube_name = f"{self.name}_{fkey}.cube"
                entry = orb.to_dict(cube_name)
                if 'canonical_label' not in entry:
                    entry['canonical_label'] = fkey.upper()
                self.manifest['canonical']['frontier'].append(entry)
                if orb.type in ('homo', 'lumo'):
                    self.manifest['canonical'][orb.type] = entry
                
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
                
            elif orb.type in ('pi_star', 'delocalized_pistar'):
                cube_name = f"{self.name}_pistar_{atoms_str}_{idx}.cube"
                self.manifest['localized']['pi_star'].append(orb.to_dict(cube_name))
            
            else:
                cube_name = f"{self.name}_{orb.type}_{atoms_str}_{idx}.cube"
                
            cube_path = os.path.join(self.out_dir, cube_name)

            # Reuse already-rendered cube when orbital snapped identically to a canonical MO
            reused = False
            for prev_vec, prev_path in rendered_cache:
                if os.path.exists(prev_path) and np.allclose(orb.mo_coeff, prev_vec, atol=1e-7):
                    if os.path.abspath(prev_path) != os.path.abspath(cube_path):
                        shutil.copyfile(prev_path, cube_path)
                    reused = True
                    print(f"  ✓ Reused identical cube for {cube_name}")
                    break

            if not reused:
                try:
                    cubegen.orbital(
                        self.mol,
                        cube_path,
                        orb.mo_coeff,
                        nx=self.grid_points,
                        ny=self.grid_points,
                        nz=self.grid_points,
                        margin=5.0,
                    )
                    rendered_cache.append((orb.mo_coeff.copy(), cube_path))
                    print(f"  ✓ Rendered {cube_name}")
                except Exception as e:
                    print(f"  ⚠ Failed to render {cube_name}: {e}")

        return self.manifest
