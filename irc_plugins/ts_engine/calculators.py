import os
import subprocess
import tempfile
import numpy as np
from ase.calculators.calculator import Calculator, all_changes
from .common import _unit

_HERE = os.path.dirname(os.path.abspath(__file__))
XTB_BIN = os.path.abspath(os.path.join(_HERE, "..", "..", "tools", "xtb", "bin", "xtb"))


class XtbCalculator(Calculator):
    """ASE Calculator wrapper for xTB binary supporting ALPB implicit solvation."""
    implemented_properties = ["energy", "forces"]

    def __init__(self, charge=0, multiplicity=1, solvent="water", **kwargs):
        super().__init__(**kwargs)
        self.charge = int(charge)
        self.multiplicity = int(multiplicity)
        self.solvent = solvent.lower() if solvent and solvent.lower() not in ("gas", "none", "vacuum", "") else None
        self._tempdir = tempfile.TemporaryDirectory(prefix="xtb_calc_")
        self._workdir = self._tempdir.name

    def calculate(self, atoms=None, properties=("energy", "forces"), system_changes=all_changes):
        super().calculate(atoms, properties, system_changes)
        xyz_file = os.path.join(self._workdir, "input.xyz")
        symbols = self.atoms.get_chemical_symbols()
        positions = self.atoms.get_positions()
        with open(xyz_file, "w") as f:
            f.write(f"{len(symbols)}\n\n")
            for s, (x, y, z) in zip(symbols, positions):
                f.write(f"{s} {x:.8f} {y:.8f} {z:.8f}\n")

        cmd = [XTB_BIN, "input.xyz", "--grad", "--chrg", str(self.charge), "--uhf", str(self.multiplicity - 1)]
        if self.solvent:
            cmd += ["--alpb", self.solvent]

        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "4")
        env.setdefault("MKL_NUM_THREADS", "1")
        res = subprocess.run(cmd, cwd=self._workdir, env=env, capture_output=True, text=True)
        e = None
        for line in res.stdout.splitlines():
            if "TOTAL ENERGY" in line:
                e = float(line.split()[3]) * 27.211386245988
                break

        grad_file = os.path.join(self._workdir, "gradient")
        forces = np.zeros((len(symbols), 3))
        if os.path.exists(grad_file):
            with open(grad_file) as f:
                lines = f.readlines()
            n = len(symbols)
            grad_lines = lines[2 + n : 2 + 2 * n]
            grads = [[float(x.replace("D", "e")) for x in gl.split()[:3]] for gl in grad_lines]
            forces = - np.array(grads) * 51.4220674763

        self.results = {"energy": e, "forces": forces}


def make_xtb(charge: int, multiplicity: int, solvent: str = None):
    """Factory for GFN2-xTB calculator with optional implicit solvent (e.g. 'water', 'alpb')."""
    if solvent and solvent.lower() not in ("gas", "none", "vacuum", ""):
        return XtbCalculator(charge=charge, multiplicity=multiplicity, solvent=solvent)
    from tblite.ase import TBLite
    return TBLite(method="GFN2-xTB", charge=charge, multiplicity=multiplicity, verbosity=0,
                  electronic_temperature=1000.0, max_iterations=600, mixer_damping=0.2)


def _apply_reaction_restraints(pos, f, restraints, k, model=None, k_angle=8.0):
    er = 0.0
    for i, j, t in restraints:
        r = pos[j] - pos[i]
        d = np.linalg.norm(r)
        u = r / d
        dd = d - t
        er += 0.5 * k * dd * dd
        f[j] -= k * dd * u
        f[i] += k * dd * u
    if model is not None:
        for f1, f2 in model.formed:
            for b1, b2 in model.broken:
                shared = None
                if f1 == b1: shared, a, b = f1, f2, b2
                elif f1 == b2: shared, a, b = f1, f2, b1
                elif f2 == b1: shared, a, b = f2, f1, b2
                elif f2 == b2: shared, a, b = f2, f1, b1
                if shared is not None:
                    va = pos[a] - pos[shared]
                    vb = pos[b] - pos[shared]
                    na, nb = np.linalg.norm(va), np.linalg.norm(vb)
                    ua, ub = va / (na + 1e-9), vb / (nb + 1e-9)
                    cos_ang = float(np.dot(ua, ub))
                    if cos_ang > -0.90:
                        diff = cos_ang - (-0.90)
                        er += 0.5 * k_angle * diff * diff
                        grad_a = (ub - cos_ang * ua) / (na + 1e-9)
                        grad_b = (ua - cos_ang * ub) / (nb + 1e-9)
                        f[a] -= k_angle * diff * grad_a
                        f[b] -= k_angle * diff * grad_b
                        f[shared] += k_angle * diff * (grad_a + grad_b)
    return er



def _apply_spectator_restraints(pos, f, spectator_restraints, k_spec):
    if not spectator_restraints:
        return 0.0
    er = 0.0
    for i, j, t0 in spectator_restraints:
        r = pos[j] - pos[i]
        d = np.linalg.norm(r)
        if d <= t0 + 0.15:
            continue
        u = r / d
        dd = d - (t0 + 0.15)
        er += 0.5 * k_spec * dd * dd
        f[j] -= k_spec * dd * u
        f[i] += k_spec * dd * u
    return er


def _make_restrained_calc(base_atoms, restraints, k, spectator_restraints=None, k_spec=50.0, model=None):
    """Calculator adding harmonic restraints for reaction coordinates and spectator protection."""
    from ase.calculators.calculator import Calculator, all_changes

    class RestrainedCalculator(Calculator):
        implemented_properties = ["energy", "forces"]

        def __init__(self):
            super().__init__()
            self.e_base = None

        def calculate(self, atoms=None, properties=("energy",), system_changes=all_changes):
            super().calculate(atoms, properties, system_changes)
            pos = self.atoms.get_positions()
            base_atoms.set_positions(pos)
            e = base_atoms.get_potential_energy()
            f = base_atoms.get_forces().copy()

            e_restr = _apply_reaction_restraints(pos, f, restraints, k, model=model)
            e_spec = _apply_spectator_restraints(pos, f, spectator_restraints, k_spec)
            self.e_base = e
            self.results = {"energy": e + e_restr + e_spec, "forces": f}

    return RestrainedCalculator()


def check_geometry_sanity(atoms, model, log, tolerance=0.35) -> bool:
    """Verify that non-reacting spectator bonds were not stretched into unphysical states."""
    pos = atoms.get_positions()
    bad_bonds = []
    for i, j in model.spectator_bonds:
        d = float(np.linalg.norm(pos[i] - pos[j]))
        d_ideal = model.d0(i, j)
        if d <= d_ideal + tolerance:
            continue
        bad_bonds.append((i, j, d, d_ideal))

    if not bad_bonds:
        return True

    log("\n⚠️ SANITY CHECK FAILED: The following spectator bond(s) were severely distorted:")
    for i, j, d, di in bad_bonds:
        log(f"   {model.label(i)}-{model.label(j)}: length {d:.2f} A (expected ~{di:.2f} A)")
    return False


def build_reaction_eigenvector(atoms, model):
    """Construct an explicit reaction vector (compressing formed, stretching broken/order-changed)."""
    pos = atoms.get_positions()
    v = np.zeros_like(pos)
    for i, j in model.formed:
        u = _unit(pos[j] - pos[i])
        v[i] += u
        v[j] -= u
    for i, j in model.broken:
        u = _unit(pos[j] - pos[i])
        v[i] -= u
        v[j] += u

    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-6 else None


def _xtb_hess_imaginary_modes(atoms, charge, multiplicity, solvent=None, imag_cutoff=25.0):
    """Run fast native xTB Hessian and return imaginary modes [(abs_freq_cm, vec_nx3), ...] > imag_cutoff."""
    if not os.path.isfile(XTB_BIN):
        return []
    syms = atoms.get_chemical_symbols()
    pos = atoms.get_positions()
    n = len(syms)
    try:
        with tempfile.TemporaryDirectory(prefix="xtb_ts_hess_") as td:
            xyz = os.path.join(td, "ts.xyz")
            with open(xyz, "w") as f:
                f.write(f"{n}\n\n")
                for s, (x, y, z) in zip(syms, pos):
                    f.write(f"{s} {x:.8f} {y:.8f} {z:.8f}\n")
            cmd = [XTB_BIN, "ts.xyz", "--hess", "--chrg", str(int(charge)), "--uhf", str(max(0, int(multiplicity) - 1))]
            if solvent and str(solvent).lower() not in ("gas", "none", "vacuum", ""):
                cmd += ["--alpb", str(solvent).lower()]
            env = os.environ.copy()
            env.setdefault("OMP_NUM_THREADS", "4")
            env.setdefault("MKL_NUM_THREADS", "1")
            subprocess.run(cmd, cwd=td, env=env, capture_output=True, text=True, timeout=120)
            g98 = os.path.join(td, "g98.out")
            if not os.path.exists(g98):
                return []
            with open(g98, "r") as f:
                lines = f.readlines()
    except Exception:
        return []

    imags = []
    i = 0
    while i < len(lines):
        if "Frequencies --" in lines[i]:
            freqs = [float(x) for x in lines[i].split("--")[1].split()]
            n_m = len(freqs)
            while i < len(lines) and "Atom AN" not in lines[i]:
                i += 1
            i += 1
            if i + n > len(lines):
                break
            vecs = [np.zeros((n, 3), dtype=float) for _ in range(n_m)]
            for a in range(n):
                parts = lines[i + a].split()
                for m in range(n_m):
                    vecs[m][a] = [
                        float(parts[2 + 3 * m]),
                        float(parts[3 + 3 * m]),
                        float(parts[4 + 3 * m]),
                    ]
            for fv, v in zip(freqs, vecs):
                if fv < -float(imag_cutoff):
                    imags.append((abs(fv), v))
            i += n
        else:
            i += 1
    imags.sort(key=lambda item: -item[0])
    return imags


def _pick_spurious_secondary_mode(atoms, model, imag_modes):
    """Given >=2 imaginary modes [(abs_freq, vec), ...], identify the primary reaction mode and return the strongest spurious secondary mode."""
    pos = atoms.get_positions()
    reacting = list(getattr(model, "formed", []) or []) + list(getattr(model, "broken", []) or [])
    if reacting:
        def _react_score(item):
            freq, vec = item
            stretch = sum(abs(float(np.dot(_unit(pos[j] - pos[i]), vec[j] - vec[i]))) for i, j in reacting)
            return (stretch, freq)
        primary_idx = int(np.argmax([_react_score(m)[0] for m in imag_modes]))
        # Guard: if stretch scores are degenerate/near-zero, fall back to highest imaginary frequency
        if _react_score(imag_modes[primary_idx])[0] < 1e-4:
            primary_idx = 0
    else:
        primary_idx = 0

    secondary = [m for idx, m in enumerate(imag_modes) if idx != primary_idx]
    secondary.sort(key=lambda item: -item[0])
    return secondary[0]


def _run_sella_pass(atoms, model, log, target_dists, max_steps=50, max_drift=0.35):
    """Run a single guarded Sella (order=1) polish pass and restore best in-trust-zone geometry."""
    best_pos = atoms.positions.copy()
    atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=getattr(model, "solvent", None))
    try:
        best_fmax = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
    except Exception:
        best_fmax = 1.0

    profile = []
    drifted = False

    try:
        from sella import Sella
        v0 = build_reaction_eigenvector(atoms, model)
        sella_kwargs = {"internal": True, "order": 1, "logfile": None}
        if v0 is not None and not sella_kwargs.get("internal"):
            sella_kwargs["v0"] = v0
        opt = Sella(atoms, **sella_kwargs)
        for step in range(max_steps):
            opt.step()
            cur_pos = atoms.positions.copy()
            cur_fmax = float(np.max(np.linalg.norm(atoms.get_forces(), axis=1)))
            profile.append({"step": step + 1, "energy_ev": float(atoms.get_potential_energy()), "fmax": cur_fmax})

            for (i, j), t in target_dists.items():
                d = np.linalg.norm(cur_pos[i] - cur_pos[j])
                if abs(d - t) > max_drift:
                    drifted = True
                    break

            if drifted:
                log(f"  Sella step {step + 1}: Reaction coordinate drifted (> {max_drift:.2f} A). Stopping Sella to preserve TS saddle.")
                break

            if cur_fmax < best_fmax:
                best_fmax = cur_fmax
                best_pos = cur_pos.copy()

            if opt.converged():
                best_pos = cur_pos.copy()
                log(f"  Sella saddle polish converged in {step + 1} steps (fmax = {cur_fmax:.3f} eV/A).")
                break
    except Exception as e:
        log(f"  Sella polish note: {e}")

    atoms.set_positions(best_pos)
    return profile


def guard_saddle_polish(atoms, model, log, max_steps=50, max_drift=0.35):
    """Safely polish a candidate saddle point with Sella DLC while guarding the reaction coordinate,
    and deterministically escape higher-order saddle points via native xTB Hessian check.
    """
    target_dists = {(i, j): float(np.linalg.norm(atoms.positions[i] - atoms.positions[j]))
                    for i, j in model.formed + model.broken}
    profile = _run_sella_pass(atoms, model, log, target_dists, max_steps=max_steps, max_drift=max_drift)

    # Deterministic Hessian check: only act if >=2 imaginary modes (> 25i cm^-1) reveal a higher-order saddle
    if len(atoms) <= 60:
        solvent = getattr(model, "solvent", None)
        for esc_pass in range(2):
            imags = _xtb_hess_imaginary_modes(atoms, model.charge, model.multiplicity, solvent=solvent, imag_cutoff=25.0)
            if len(imags) <= 1:
                break
            freq_str = ", ".join(f"{f:.0f}i" for f, _ in imags)
            sec_f, sec_vec = _pick_spurious_secondary_mode(atoms, model, imags)
            max_disp = float(np.max(np.linalg.norm(sec_vec, axis=1)))
            if max_disp <= 1e-6:
                break
            log(f"  Higher-order saddle detected ({freq_str} cm-1); escaping along spurious {sec_f:.0f}i mode (pass {esc_pass + 1}/2)...")
            pre_pos = atoms.positions.copy()
            atoms.set_positions(pre_pos + (0.25 / max_disp) * sec_vec)
            esc_targets = {(i, j): float(np.linalg.norm(pre_pos[i] - pre_pos[j]))
                           for i, j in model.formed + model.broken}
            esc_profile = _run_sella_pass(
                atoms, model, log, esc_targets,
                max_steps=max(max_steps, 50),
                max_drift=max(max_drift, 0.35),
            )
            if not check_geometry_sanity(atoms, model, log):
                log("  Secondary saddle escape distorted spectator bonds; restoring prior TS geometry.")
                atoms.set_positions(pre_pos)
                break
            profile.extend(esc_profile)

    return profile

