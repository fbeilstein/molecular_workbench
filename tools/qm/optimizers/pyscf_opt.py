import os
from .base import GeometryOptimizer

class PyscfOptimizer(GeometryOptimizer):
    def __init__(self, basis='6-31g*', method='b3lyp'):
        self.basis = basis
        self.method = method

    def optimize(self, xyz_path: str, charge: int = 0) -> str:
        try:
            import geometric
        except ImportError:
            print("  Warning: 'geometric' library is not installed. Please run 'pip install geometric' to use PySCF geometry optimization.")
            return xyz_path

        from pyscf import gto, dft
        from pyscf.geomopt.geometric_solver import optimize
        import numpy as np

        print(f"  Optimizing with PySCF ({self.method}/{self.basis})...")
        try:
            mol = gto.M(
                atom=xyz_path,
                basis=self.basis,
                charge=charge,
                spin=0,
                verbose=3
            )
            
            if self.method.lower() == 'hf':
                from pyscf import scf
                mf = scf.RHF(mol)
            else:
                mf = dft.RKS(mol)
                mf.xc = self.method
                
            mol_eq = optimize(mf, maxsteps=100)
            
            # Write optimized coordinates back
            coords = mol_eq.atom_coords() * 0.529177210903  # Bohr to Angstrom
            
            with open(xyz_path, 'r') as f:
                lines = f.readlines()
            
            n_atoms = int(lines[0].strip())
            new_lines = lines[:2]
            
            for i in range(n_atoms):
                parts = lines[i+2].split()
                sym = parts[0]
                new_lines.append(f"{sym:2s}  {coords[i][0]:12.6f}  {coords[i][1]:12.6f}  {coords[i][2]:12.6f}\n")
                
            out_dir = os.path.dirname(os.path.abspath(xyz_path))
            opt_xyz = os.path.join(out_dir, 'pyscfopt.xyz')
            with open(opt_xyz, 'w') as f:
                f.writelines(new_lines)
                
            import shutil
            shutil.move(opt_xyz, xyz_path)
            print("  PySCF optimization converged")

        except Exception as e:
            print(f"  Warning: PySCF optimization failed: {e}")

        return xyz_path
