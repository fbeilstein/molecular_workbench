import os
import subprocess
from .base import GeometryOptimizer

class XtbOptimizer(GeometryOptimizer):
    def _find_xtb(self):
        """Find xTB binary."""
        _here = os.path.dirname(os.path.abspath(__file__))
        _tools = os.path.abspath(os.path.join(_here, '..', '..'))
        candidates = [
            os.path.join(_tools, 'xtb', 'bin', 'xtb'),
            'xtb',
        ]
        for c in candidates:
            c = os.path.abspath(c) if c != 'xtb' else c
            if c != 'xtb' and os.path.isfile(c) and os.access(c, os.X_OK):
                return c
            elif c == 'xtb':
                import shutil
                if shutil.which('xtb'):
                    return 'xtb'
        return None

    def optimize(self, xyz_path: str, charge: int = 0) -> str:
        xtb_bin = self._find_xtb()
        if not xtb_bin:
            print("  Warning: xTB binary not found. Skipping xTB optimization.")
            return xyz_path
            
        out_dir = os.path.dirname(os.path.abspath(xyz_path))
        print("  Optimizing with xTB...")
        try:
            cmd = [xtb_bin, xyz_path, '--opt', '--chrg', str(charge)]
            result = subprocess.run(cmd, capture_output=True, text=True,
                                    cwd=out_dir, timeout=120)

            opt_xyz = os.path.join(out_dir, 'xtbopt.xyz')
            if os.path.exists(opt_xyz):
                import shutil
                shutil.move(opt_xyz, xyz_path)
                print("  xTB optimization converged")
            else:
                print("  Warning: xTB did not produce optimized geometry")

        except subprocess.TimeoutExpired:
            print("  Warning: xTB timed out after 120s, using MMFF geometry")
        except Exception as e:
            print(f"  Warning: xTB failed: {e}")

        return xyz_path
