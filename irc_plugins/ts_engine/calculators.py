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

        res = subprocess.run(cmd, cwd=self._workdir, capture_output=True, text=True)
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


def guard_saddle_polish(atoms, model, log, max_steps=50, max_drift=0.35):
    """Safely polish a candidate saddle point with Sella DLC while guarding the reaction coordinate.

    If Sella drifts away from the transition state (> max_drift A on reacting bonds),
    Sella is immediately aborted and the best geometry within the TS trust zone is restored.
    """
    target_dists = {(i, j): float(np.linalg.norm(atoms.positions[i] - atoms.positions[j]))
                    for i, j in model.formed + model.broken}
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

