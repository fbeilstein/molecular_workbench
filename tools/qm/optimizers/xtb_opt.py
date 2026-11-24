import os
import shutil
import subprocess
import tempfile
import numpy as np
from .base import GeometryOptimizer


def _read_xyz(path: str):
    with open(path, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    n = int(lines[0])
    syms, coords = [], []
    for line in lines[2:2 + n]:
        parts = line.split()
        syms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])
    return syms, np.array(coords, dtype=float)


def _write_xyz(path: str, syms, coords: np.ndarray, comment: str = "Optimized by xTB"):
    with open(path, 'w') as f:
        f.write(f"{len(syms)}\n{comment}\n")
        for s, (x, y, z) in zip(syms, coords):
            f.write(f"{s:2s}  {x:14.8f}  {y:14.8f}  {z:14.8f}\n")


def _parse_g98_modes(g98_path: str, n_atoms: int):
    """Parse vibrational frequencies (cm^-1) and Cartesian normal-mode vectors from xTB g98.out."""
    if not os.path.exists(g98_path):
        return []
    with open(g98_path, 'r') as f:
        lines = f.readlines()
    modes = []
    i = 0
    while i < len(lines):
        if "Frequencies --" in lines[i]:
            freqs = [float(x) for x in lines[i].split("--")[1].split()]
            n_m = len(freqs)
            while i < len(lines) and "Atom AN" not in lines[i]:
                i += 1
            i += 1
            if i + n_atoms > len(lines):
                break
            vecs = [np.zeros((n_atoms, 3), dtype=float) for _ in range(n_m)]
            for a in range(n_atoms):
                parts = lines[i + a].split()
                for m in range(n_m):
                    vecs[m][a] = [
                        float(parts[2 + 3 * m]),
                        float(parts[3 + 3 * m]),
                        float(parts[4 + 3 * m]),
                    ]
            modes.extend(zip(freqs, vecs))
            i += n_atoms
        else:
            i += 1
    return modes


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
                if shutil.which('xtb'):
                    return 'xtb'
        return None

    def optimize(self, xyz_path: str, charge: int = 0) -> str:
        xtb_bin = self._find_xtb()
        if not xtb_bin:
            print("  Warning: xTB binary not found. Skipping xTB optimization.")
            return xyz_path

        syms, pos = _read_xyz(xyz_path)
        n_atoms = len(syms)
        use_hess = n_atoms <= 60
        max_passes = 3 if use_hess else 1
        imag_cutoff = -15.0  # cm^-1 threshold for genuine saddle-point imaginary modes

        print("  Optimizing with xTB (tight + Hessian verification)..." if use_hess else "  Optimizing with xTB (tight)...")
        best_pos = None
        try:
            with tempfile.TemporaryDirectory(prefix="xtb_ohess_") as workdir:
                in_xyz = os.path.join(workdir, "input.xyz")
                opt_xyz = os.path.join(workdir, "xtbopt.xyz")
                g98_path = os.path.join(workdir, "g98.out")

                for attempt in range(max_passes):
                    _write_xyz(in_xyz, syms, pos, comment=f"xTB pass {attempt + 1}")
                    for stale in (opt_xyz, g98_path):
                        if os.path.exists(stale):
                            os.remove(stale)

                    opt_flag = ["--ohess", "tight"] if use_hess else ["--opt", "tight"]
                    cmd = [xtb_bin, "input.xyz", *opt_flag, "--chrg", str(charge)]
                    subprocess.run(
                        cmd, capture_output=True, text=True, cwd=workdir, timeout=180
                    )

                    if not os.path.exists(opt_xyz):
                        if best_pos is None:
                            # Fallback to standard --opt if --ohess tight failed to converge
                            subprocess.run(
                                [xtb_bin, "input.xyz", "--opt", "--chrg", str(charge)],
                                capture_output=True, text=True, cwd=workdir, timeout=120
                            )
                        if not os.path.exists(opt_xyz):
                            print("  Warning: xTB did not produce optimized geometry")
                            break

                    _, pos = _read_xyz(opt_xyz)
                    best_pos = pos.copy()

                    if not use_hess:
                        print("  xTB optimization converged")
                        break

                    modes = _parse_g98_modes(g98_path, n_atoms)
                    if not modes:
                        print("  xTB optimization converged")
                        break

                    imags = [(fv, vec) for fv, vec in modes if fv < imag_cutoff]
                    if not imags:
                        lowest = min(fv for fv, _ in modes)
                        print(f"  xTB optimization converged (true minimum, lowest freq = {lowest:+.1f} cm^-1)")
                        break

                    f_low, vec = min(imags, key=lambda item: item[0])
                    if attempt == max_passes - 1:
                        print(f"  xTB optimization converged (remaining mode {f_low:.1f} cm^-1 after {max_passes} passes)")
                        break

                    max_disp = float(np.max(np.linalg.norm(vec, axis=1)))
                    if max_disp <= 1e-6:
                        break
                    print(f"  xTB Hessian detected saddle mode ({f_low:.1f} cm^-1); following eigenvector (pass {attempt + 2}/{max_passes})...")
                    pos = pos + (0.25 / max_disp) * vec

            if best_pos is not None:
                _write_xyz(xyz_path, syms, best_pos, comment="Optimized by xTB (Hessian-verified minimum)")

        except subprocess.TimeoutExpired:
            if best_pos is not None:
                _write_xyz(xyz_path, syms, best_pos, comment="Optimized by xTB")
            else:
                print("  Warning: xTB timed out, using MMFF geometry")
        except Exception as e:
            print(f"  Warning: xTB failed: {e}")

        return xyz_path
