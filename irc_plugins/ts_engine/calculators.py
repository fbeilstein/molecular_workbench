import numpy as np
from .common import _unit


def make_xtb(charge: int, multiplicity: int):
    """Factory for GFN2-xTB calculator with convergence parameters for charged/anion systems."""
    from tblite.ase import TBLite
    return TBLite(method="GFN2-xTB", charge=charge, multiplicity=multiplicity, verbosity=0,
                  electronic_temperature=1000.0, max_iterations=600, mixer_damping=0.2)


def _apply_reaction_restraints(pos, f, restraints, k):
    er = 0.0
    for i, j, t in restraints:
        r = pos[j] - pos[i]
        d = np.linalg.norm(r)
        u = r / d
        dd = d - t
        er += 0.5 * k * dd * dd
        f[j] -= k * dd * u
        f[i] += k * dd * u
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


def _make_restrained_calc(base_atoms, restraints, k, spectator_restraints=None, k_spec=50.0):
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

            e_restr = _apply_reaction_restraints(pos, f, restraints, k)
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
    for i, j, _, _ in model.order_changed:
        u = _unit(pos[j] - pos[i])
        v[i] -= u * 0.5
        v[j] += u * 0.5

    norm = np.linalg.norm(v)
    return v / norm if norm > 1e-6 else None
