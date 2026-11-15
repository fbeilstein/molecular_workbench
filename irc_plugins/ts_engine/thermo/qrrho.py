"""
Grimme Quasi-Rigid-Rotor Harmonic-Oscillator (qRRHO) Solution-Phase Thermochemistry Engine.

Computes standard-state Gibbs free energy in solution:
    G°_soln(T) = E_el,solv + G_corr^qRRHO(T) + ΔG°->*(T)

References:
    - S. Grimme, Chem. Eur. J. 2012, 18, 9955-9964.
    - C. A. P. P. et al., J. Phys. Chem. B 2019, 123, 30, 6601-6609.
"""

import os
import subprocess
import tempfile
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from ase import Atoms
from ase.data import atomic_masses, covalent_radii

from ..calculators import XTB_BIN

# Universal Physical Constants
HARTREE_TO_KCAL = 627.5094740631
EV_TO_KCAL = 23.06054885
R_CAL = 1.98720425864083          # cal / (mol * K)
R_KCAL = R_CAL / 1000.0           # kcal / (mol * K)
R_J = 8.314462618                 # J / (mol * K)
KB_J = 1.380649e-23               # J / K
PLANCK_H = 6.62607015e-34         # J * s
SPEED_OF_LIGHT = 2.99792458e10    # cm / s
AVOGADRO_NA = 6.02214076e23       # mol^-1
AMU_TO_KG = 1.66053906660e-27     # kg
ATM_TO_PA = 101325.0              # Pa

# Standard-state solution phase compression: 1 atm (24.46 L/mol) -> 1 M (1 mol/L) at 298.15 K
# ΔG°->* = R * T * ln(24.465) = +1.894 kcal/mol per independent translating solute molecule
STANDARD_STATE_COMPRESSION_KCAL = 1.894

PROTIC_SOLVENTS = {"water", "h2o", "methanol", "meoh", "ethanol", "etoh", "alcohol"}
PROTIC_HALIDE_HB_CORR_KCAL = {
    "F": -7.50,
    "Cl": -5.50,
    "Br": -4.50,
    "I": -3.50,
}


def sackur_tetrode_single_atom(
    symbol: str,
    e_el_kcal: float = 0.0,
    T: float = 298.15,
    pressure_pa: float = ATM_TO_PA,
    charge: int = 0,
    solvent: Optional[str] = "water"
) -> Dict[str, Any]:
    """
    Exact analytical Sackur-Tetrode translational thermochemistry for a single atom/ion.
    For N_atoms = 1, rotational and vibrational degrees of freedom are strictly zero.
    Includes first-shell protic hydrogen-bond solvation correction for bare halide anions.
    """
    from ase.data import atomic_numbers
    z = atomic_numbers.get(symbol, 1)
    mass_amu = float(atomic_masses[z])
    mass_kg = mass_amu * AMU_TO_KG

    # Translational partition function (ideal gas at T, P)
    # q_trans = (2 * pi * m * k_B * T / h^2)^(3/2) * (k_B * T / P)
    q_trans = ((2.0 * np.pi * mass_kg * KB_J * T) / (PLANCK_H ** 2)) ** 1.5 * (KB_J * T / pressure_pa)
    s_trans_cal = R_CAL * (np.log(max(q_trans, 1e-12)) + 2.5)
    ts_trans_kcal = T * s_trans_cal / 1000.0
    h_trans_kcal = 2.5 * R_KCAL * T

    # First-shell hydrogen-bonding hydration for bare halide anions in protic solvents
    hb_corr_kcal = 0.0
    if int(charge) == -1 and solvent and solvent.lower().strip() in PROTIC_SOLVENTS:
        hb_corr_kcal = PROTIC_HALIDE_HB_CORR_KCAL.get(symbol, 0.0)

    # G_corr in solution includes 1 atm -> 1 M standard state correction and protic HB hydration
    g_corr_kcal = h_trans_kcal - ts_trans_kcal + STANDARD_STATE_COMPRESSION_KCAL + hb_corr_kcal

    return {
        "G_soln_kcal": round(float(e_el_kcal + g_corr_kcal), 4),
        "E_el_kcal": round(float(e_el_kcal), 4),
        "G_corr_kcal": round(float(g_corr_kcal), 4),
        "ZPE_kcal": 0.0,
        "H_corr_kcal": round(float(h_trans_kcal), 4),
        "TS_qrrho_kcal": round(float(ts_trans_kcal), 4),
        "n_imag": 0,
        "imag_freq_cm": 0.0,
        "frequencies": [],
        "method": "sackur_tetrode_analytical"
    }


def _run_xtb_hessian_thermo(
    atoms: Atoms,
    charge: int = 0,
    mult: int = 1,
    solvent: Optional[str] = "water",
    T: float = 298.15,
) -> Optional[Dict[str, Any]]:
    """Run xtb --hess to obtain native Grimme qRRHO thermochemistry."""
    if not os.path.exists(XTB_BIN):
        return None

    symbols = atoms.get_chemical_symbols()
    positions = atoms.get_positions()
    n_atoms = len(symbols)
    if n_atoms < 2:
        return None

    with tempfile.TemporaryDirectory(prefix="xtb_hess_") as td:
        xyz_file = os.path.join(td, "mol.xyz")
        with open(xyz_file, "w") as f:
            f.write(f"{n_atoms}\n\n")
            for s, (x, y, z) in zip(symbols, positions):
                f.write(f"{s:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")

        cmd = [XTB_BIN, "mol.xyz", "--hess", "--chrg", str(int(charge)), "--uhf", str(int(mult) - 1)]
        if solvent and solvent.lower() not in ("gas", "none", "vacuum", ""):
            cmd += ["--alpb", solvent.lower()]

        env = os.environ.copy()
        env.setdefault("OMP_NUM_THREADS", "4")
        env.setdefault("MKL_NUM_THREADS", "1")
        res = subprocess.run(cmd, cwd=td, env=env, capture_output=True, text=True)

        e_el_eh = None
        zpe_eh = None
        h_tot_eh = None
        g_tot_eh = None
        n_imag = 0
        imag_freqs = []

        def _get_eh(line: str) -> Optional[float]:
            tokens = line.replace("|", " ").replace("::", " ").replace(":", " ").split()
            for t in tokens:
                try:
                    return float(t)
                except ValueError:
                    continue
            return None

        lines = res.stdout.splitlines()
        for i, line in enumerate(lines):
            if "TOTAL ENERGY" in line and e_el_eh is None:
                e_el_eh = _get_eh(line)
            elif "zero point energy" in line and zpe_eh is None:
                zpe_eh = _get_eh(line)
            elif "TOTAL ENTHALPY" in line and h_tot_eh is None:
                h_tot_eh = _get_eh(line)
            elif "TOTAL FREE ENERGY" in line and g_tot_eh is None:
                g_tot_eh = _get_eh(line)
            elif "# imaginary freq." in line:
                val = _get_eh(line)
                if val is not None:
                    n_imag = int(val)

        # Parse frequencies from g98.out or vibspectrum if present
        vib_file = os.path.join(td, "vibspectrum")
        all_freqs = []
        if os.path.exists(vib_file):
            try:
                with open(vib_file) as vf:
                    for vline in vf:
                        vparts = vline.split()
                        if len(vparts) >= 2:
                            try:
                                freq_val = float(vparts[0])
                                all_freqs.append(freq_val)
                                if freq_val < -10.0:
                                    imag_freqs.append(abs(freq_val))
                            except ValueError:
                                pass
            except Exception:
                pass

        if e_el_eh is None or g_tot_eh is None:
            return None

        e_el_kcal = e_el_eh * HARTREE_TO_KCAL
        g_tot_kcal = g_tot_eh * HARTREE_TO_KCAL
        h_tot_kcal = (h_tot_eh * HARTREE_TO_KCAL) if h_tot_eh is not None else e_el_kcal
        zpe_kcal = (zpe_eh * HARTREE_TO_KCAL) if zpe_eh is not None else 0.0

        # Raw gas-phase thermal correction from xtb qRRHO
        raw_g_corr_kcal = g_tot_kcal - e_el_kcal
        raw_h_corr_kcal = h_tot_kcal - e_el_kcal
        ts_qrrho_kcal = h_tot_kcal - g_tot_kcal

        # Add 1 atm -> 1 M solution-phase standard-state compression
        g_corr_kcal = raw_g_corr_kcal + STANDARD_STATE_COMPRESSION_KCAL
        g_soln_kcal = e_el_kcal + g_corr_kcal

        primary_imag = imag_freqs[0] if imag_freqs else 0.0

        return {
            "G_soln_kcal": round(float(g_soln_kcal), 4),
            "E_el_kcal": round(float(e_el_kcal), 4),
            "G_corr_kcal": round(float(g_corr_kcal), 4),
            "ZPE_kcal": round(float(zpe_kcal), 4),
            "H_corr_kcal": round(float(raw_h_corr_kcal), 4),
            "TS_qrrho_kcal": round(float(ts_qrrho_kcal), 4),
            "n_imag": int(n_imag),
            "imag_freq_cm": round(float(primary_imag), 2),
            "frequencies": [round(f, 1) for f in all_freqs],
            "method": "xtb_qrrho_native"
        }


def _python_qrrho_from_frequencies(
    atoms: Atoms,
    freqs_cm: List[float],
    e_el_kcal: float,
    T: float = 298.15,
    nu0_cm: float = 100.0,
    is_ts: bool = False
) -> Dict[str, Any]:
    """Pure-Python Grimme qRRHO thermochemistry implementation."""
    symbols = atoms.get_chemical_symbols()
    masses = atoms.get_masses()
    n_atoms = len(symbols)

    if n_atoms == 1:
        return sackur_tetrode_single_atom(symbols[0], e_el_kcal=e_el_kcal, T=T)

    # Separate real and imaginary modes
    real_freqs = []
    imag_freqs = []
    for f in freqs_cm:
        if isinstance(f, complex):
            if abs(f.imag) > 10.0:
                imag_freqs.append(abs(f.imag))
            elif f.real > 10.0:
                real_freqs.append(f.real)
        else:
            if f < -10.0:
                imag_freqs.append(abs(f))
            elif f > 10.0:
                real_freqs.append(f)

    # For TS, project out the primary reaction mode
    if is_ts and imag_freqs:
        pass  # imaginary mode already not in real_freqs

    # Average moment of inertia for free rotor interpolation
    # B_av = h / (8 * pi^2 * I_av)
    positions = atoms.get_positions()
    center_of_mass = np.average(positions, axis=0, weights=masses)
    rel_pos = positions - center_of_mass
    inertia_tensor = np.zeros((3, 3))
    for m, r in zip(masses, rel_pos):
        inertia_tensor += m * (np.dot(r, r) * np.eye(3) - np.outer(r, r))
    principal_moments = np.sort(np.linalg.eigvalsh(inertia_tensor))
    # Filter positive moments of inertia
    pos_moments = principal_moments[principal_moments > 1e-4]
    b_av_cm = 1.0
    if len(pos_moments) > 0:
        i_av_kg_m2 = np.mean(pos_moments) * AMU_TO_KG * 1e-20
        b_av_cm = (PLANCK_H / (8.0 * (np.pi ** 2) * i_av_kg_m2 * SPEED_OF_LIGHT)) if i_av_kg_m2 > 0 else 1.0

    zpe_kcal = 0.0
    u_vib_kcal = 0.0
    s_vib_cal = 0.0

    beta = 1.0 / (KB_J * T)
    for nu in real_freqs:
        nu_si = nu * SPEED_OF_LIGHT
        h_nu = PLANCK_H * nu_si
        x = h_nu * beta

        # ZPE = 0.5 * h * c * nu
        zpe_kcal += (0.5 * h_nu * AVOGADRO_NA) / (4184.0)

        # Thermal vibrational energy: U_vib = h*nu / (exp(x) - 1)
        if x < 100.0:
            exp_x = np.exp(x)
            u_vib_kcal += (h_nu * AVOGADRO_NA) / (4184.0 * (exp_x - 1.0))
            s_ho_cal = R_CAL * (x / (exp_x - 1.0) - np.log(1.0 - np.exp(-x)))
        else:
            s_ho_cal = 0.0

        # Grimme rotor entropy
        # mu' = h / (8 * pi^2 * c * nu)
        mu_prime = PLANCK_H / (8.0 * (np.pi ** 2) * SPEED_OF_LIGHT * max(nu, 1.0))
        # S_free_rotor
        s_fr_cal = R_CAL * (0.5 + np.log(max(1e-12, np.sqrt((8.0 * (np.pi ** 3) * mu_prime * KB_J * T) / (PLANCK_H ** 2)))))

        # Damping function w(nu) = 1 / (1 + (nu0 / nu)^4)
        w = 1.0 / (1.0 + (nu0_cm / max(nu, 1e-6)) ** 4)
        s_vib_cal += w * s_ho_cal + (1.0 - w) * s_fr_cal

    # Translational & Rotational Enthalpy and Entropy
    total_mass_amu = float(np.sum(masses))
    total_mass_kg = total_mass_amu * AMU_TO_KG
    q_trans = ((2.0 * np.pi * total_mass_kg * KB_J * T) / (PLANCK_H ** 2)) ** 1.5 * (KB_J * T / ATM_TO_PA)
    s_trans_cal = R_CAL * (np.log(max(q_trans, 1e-12)) + 2.5)
    h_trans_kcal = 2.5 * R_KCAL * T

    is_linear = len(pos_moments) == 2 or (len(pos_moments) == 3 and pos_moments[0] < 1e-3)
    if is_linear:
        h_rot_kcal = 1.0 * R_KCAL * T
        ia_kg_m2 = pos_moments[-1] * AMU_TO_KG * 1e-20
        q_rot = (8.0 * (np.pi ** 2) * ia_kg_m2 * KB_J * T) / (PLANCK_H ** 2)
        s_rot_cal = R_CAL * (np.log(max(q_rot, 1e-12)) + 1.0)
    else:
        h_rot_kcal = 1.5 * R_KCAL * T
        ia, ib, ic = [(m * AMU_TO_KG * 1e-20) for m in pos_moments[-3:]]
        q_rot = (np.sqrt(np.pi) / 1.0) * ((8.0 * (np.pi ** 2) * KB_J * T / (PLANCK_H ** 2)) ** 1.5) * np.sqrt(max(1e-60, ia * ib * ic))
        s_rot_cal = R_CAL * (np.log(max(q_rot, 1e-12)) + 1.5)

    s_tot_cal = s_trans_cal + s_rot_cal + s_vib_cal
    ts_tot_kcal = T * s_tot_cal / 1000.0
    h_tot_kcal = zpe_kcal + u_vib_kcal + h_trans_kcal + h_rot_kcal

    g_corr_kcal = h_tot_kcal - ts_tot_kcal + STANDARD_STATE_COMPRESSION_KCAL
    g_soln_kcal = e_el_kcal + g_corr_kcal

    return {
        "G_soln_kcal": round(float(g_soln_kcal), 4),
        "E_el_kcal": round(float(e_el_kcal), 4),
        "G_corr_kcal": round(float(g_corr_kcal), 4),
        "ZPE_kcal": round(float(zpe_kcal), 4),
        "H_corr_kcal": round(float(h_tot_kcal), 4),
        "TS_qrrho_kcal": round(float(ts_tot_kcal), 4),
        "n_imag": len(imag_freqs),
        "imag_freq_cm": round(float(imag_freqs[0]), 2) if imag_freqs else 0.0,
        "frequencies": [round(f, 1) for f in real_freqs],
        "method": "python_qrrho"
    }


def compute_species_gibbs_free_energy(
    atoms: Atoms,
    model: Optional[Any] = None,
    charge: int = 0,
    mult: int = 1,
    solvent: Optional[str] = "water",
    T: float = 298.15,
    nu0_cm: float = 100.0,
    is_ts: bool = False,
    e_el_override_kcal: Optional[float] = None
) -> Dict[str, Any]:
    """
    Compute rigorous solution-phase Gibbs free energy G°_soln for a molecular species or TS.
    """
    symbols = atoms.get_chemical_symbols()
    n_atoms = len(symbols)

    if n_atoms == 1:
        if e_el_override_kcal is not None:
            e_el = float(e_el_override_kcal)
        else:
            # Single-point xTB electronic energy with solvent
            try:
                from ..calculators import make_xtb
                at_sp = atoms.copy()
                at_sp.calc = make_xtb(charge, mult, solvent=solvent)
                e_el = float(at_sp.get_potential_energy()) * EV_TO_KCAL
            except Exception:
                e_el = 0.0
        return sackur_tetrode_single_atom(symbols[0], e_el_kcal=e_el, T=T, charge=charge, solvent=solvent)

    # Try fast native xtb --hess
    thermo = _run_xtb_hessian_thermo(atoms, charge=charge, mult=mult, solvent=solvent, T=T)
    if thermo is not None:
        if e_el_override_kcal is not None:
            thermo["E_el_kcal"] = round(float(e_el_override_kcal), 4)
            thermo["G_soln_kcal"] = round(float(e_el_override_kcal + thermo["G_corr_kcal"]), 4)
        return thermo

    # Fallback to pure-Python qRRHO using ASE vibrations
    from ..vibrations import _compute_vibrations
    try:
        from ..calculators import make_xtb
        at_vib = atoms.copy()
        at_vib.calc = make_xtb(charge, mult, solvent=solvent)
        freqs, _ = _compute_vibrations(at_vib)
        if e_el_override_kcal is not None:
            e_el = float(e_el_override_kcal)
        else:
            e_el = float(at_vib.get_potential_energy()) * EV_TO_KCAL
        return _python_qrrho_from_frequencies(atoms, freqs, e_el_kcal=e_el, T=T, nu0_cm=nu0_cm, is_ts=is_ts)
    except Exception as err:
        # Emergency safety fallback
        e_el = float(e_el_override_kcal if e_el_override_kcal is not None else 0.0)
        return {
            "G_soln_kcal": e_el,
            "E_el_kcal": e_el,
            "G_corr_kcal": 0.0,
            "ZPE_kcal": 0.0,
            "H_corr_kcal": 0.0,
            "TS_qrrho_kcal": 0.0,
            "n_imag": 0,
            "imag_freq_cm": 0.0,
            "frequencies": [],
            "method": f"fallback_error: {err}"
        }


def split_into_independent_fragments(atoms: Atoms, scale: float = 1.30) -> List[List[int]]:
    """
    Find disconnected molecular fragments based on covalent bonding threshold.
    Returns list of atom index lists (e.g. [[0, 1, 2], [3]]).
    """
    n_atoms = len(atoms)
    if n_atoms <= 1:
        return [list(range(n_atoms))]

    pos = atoms.get_positions()
    atomic_nums = atoms.get_atomic_numbers()
    radii = np.array([covalent_radii[z] for z in atomic_nums])

    # Build adjacency
    adj = [[] for _ in range(n_atoms)]
    for i in range(n_atoms):
        for j in range(i + 1, n_atoms):
            d = np.linalg.norm(pos[i] - pos[j])
            thresh = scale * (radii[i] + radii[j])
            if d <= thresh:
                adj[i].append(j)
                adj[j].append(i)

    # Connected components via BFS
    visited = [False] * n_atoms
    fragments = []
    for i in range(n_atoms):
        if not visited[i]:
            comp = []
            queue = [i]
            visited[i] = True
            while queue:
                curr = queue.pop(0)
                comp.append(curr)
                for neighbor in adj[curr]:
                    if not visited[neighbor]:
                        visited[neighbor] = True
                        queue.append(neighbor)
            fragments.append(sorted(comp))

    return fragments


def compute_fragment_charges_and_multiplicities(
    atoms: Atoms,
    fragment_indices: List[List[int]],
    total_charge: int = 0,
    total_mult: int = 1,
    solvent: str = "water"
) -> List[Tuple[int, int]]:
    """
    Compute integer formal charges and spin multiplicities for disconnected fragments.
    """
    n_frags = len(fragment_indices)
    if n_frags <= 1:
        return [(int(total_charge), int(total_mult))]

    # Get atomic partial charges via xTB
    symbols = atoms.get_chemical_symbols()
    atomic_nums = atoms.get_atomic_numbers()
    partial_charges = None

    if os.path.exists(XTB_BIN):
        with tempfile.TemporaryDirectory(prefix="xtb_chrg_") as td:
            xyz_file = os.path.join(td, "mol.xyz")
            with open(xyz_file, "w") as f:
                f.write(f"{len(symbols)}\n\n")
                for s, (x, y, z) in zip(symbols, atoms.get_positions()):
                    f.write(f"{s:<2} {x:12.6f} {y:12.6f} {z:12.6f}\n")
            cmd = [XTB_BIN, "mol.xyz", "--grad", "--chrg", str(total_charge), "--uhf", str(total_mult - 1)]
            if solvent:
                cmd += ["--alpb", solvent]
            res = subprocess.run(cmd, cwd=td, capture_output=True, text=True)
            chrg_file = os.path.join(td, "charges")
            if os.path.exists(chrg_file):
                try:
                    with open(chrg_file) as cf:
                        partial_charges = [float(line.strip()) for line in cf if line.strip()]
                except Exception:
                    pass

    if partial_charges is None or len(partial_charges) != len(symbols):
        # Fallback to equal or electronegativity distribution
        partial_charges = [float(total_charge) / len(symbols)] * len(symbols)

    # Sum partial charges per fragment
    frag_raw_charges = [sum(partial_charges[i] for i in f_idx) for f_idx in fragment_indices]
    
    # Round to nearest integer conserving total charge
    rounded_charges = [int(round(q)) for q in frag_raw_charges]
    diff = total_charge - sum(rounded_charges)
    if diff != 0:
        # Adjust the fragment with largest rounding error
        errors = [abs(frag_raw_charges[k] - rounded_charges[k]) for k in range(n_frags)]
        k_adjust = int(np.argmax(errors))
        rounded_charges[k_adjust] += diff

    # Determine multiplicities from electron counts
    results = []
    for f_idx, q in zip(fragment_indices, rounded_charges):
        n_elec = sum(atomic_nums[i] for i in f_idx) - q
        mult = 1 if (n_elec % 2 == 0) else 2
        results.append((int(q), int(mult)))

    return results
