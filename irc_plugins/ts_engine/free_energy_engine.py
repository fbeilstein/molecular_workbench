"""
Free Energy and Thermodynamics Engine for Chemical Reactions.

Provides mathematically rigorous, physically continuous Free Energy (G) profiles
for single-step IRC pathways and multistep reaction cascades in solution.

Physical Theory:
----------------
Along a reaction path on the Born-Oppenheimer potential energy surface, electronic structure
methods (xTB, DFT, etc.) compute the raw electronic potential energy E_elec.
In solution (standard state: 298.15 K, 1 M, 1 atm):
  G°_solv = E_elec + ΔG_solv + G_vib - T·ΔS_trans/rot

For an elementary step R -> [TS]‡ -> P:
1. Association (N_R > 1 -> N_P = 1, e.g. A + B -> C):
   N_R independent molecules condense into 1 transition state complex.
   Translational & rotational entropy loss:
     ΔG_assoc° = (N_R - 1) × 13.76 kcal/mol  (Whiting-Pacey / Amzel-Eckstein model)
   The resulting adduct (C) is ALSO 1 bound molecule and maintains this identical
   entropic baseline. Hence, no discontinuity exists between TS and bound intermediates.

2. Dissociation (N_R = 1 -> N_P > 1, e.g. A -> B + C or A -> B + C + D):
   1 bound reactant forms 1 bound TS, which subsequently fragments into N_P free particles.
   Entropy gain upon asymptotic separation:
     ΔG_dissoc° = -(N_P - 1) × 13.76 kcal/mol
   As fragment separation d increases beyond contact distance, free energy smoothly decreases
   into the fully dissociated asymptotic limit.

3. Spectators / Catalysts (A + B -> C + A):
   Spectator species present on both reactant and product sides cancel identically in net ΔG_rxn°
   via Hess's Law. If the spectator coordinates to the TS, it pays association entropy at TS
   which is smoothly restored upon product departure.

Continuity Guarantee:
---------------------
Energy transitions are modulated by an analytic switching function λ(d) with C¹ smoothness:
  λ(d) = sin²((π / 2) · clip((d - d_contact) / (d_inf - d_contact), 0, 1))
This guarantees zero cliffs, zero artificial jumps, and strict derivative continuity
across the entire reaction coordinate s.
"""

from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from ase import Atoms
from rdkit import Chem
from rdkit.Chem import AllChem

from .common import EV_TO_KCAL
from .calculators import make_xtb

# Standard-state association free energy penalty in solution (kcal/mol at 298.15 K, 1 M)
# References: Whiting & Pacey (1974); Amzel & Eckstein (1979); J. Am. Chem. Soc. 1998, 120, 24, 6032-6040
DELTA_G_ASSOC_KCAL: float = 13.76
DELTA_G_ASSOC_EV: float = DELTA_G_ASSOC_KCAL / EV_TO_KCAL

# Interaction distance thresholds (Angstrom)
DEFAULT_D_CONTACT_MARGIN: float = 0.60   # Contact radius beyond covalent equilibrium d0
DEFAULT_D_INF_MARGIN: float = 3.20       # Asymptotic separation distance threshold


def smooth_switching_function(d: float, d_contact: float, d_inf: float) -> float:
    """
    Compute C¹ smooth switching parameter λ(d) ∈ [0, 1].
    λ = 0 when d <= d_contact (bound complex)
    λ = 1 when d >= d_inf (isolated free molecules)
    """
    if d_inf <= d_contact:
        return 0.0 if d <= d_contact else 1.0
    u = float(np.clip((d - d_contact) / (d_inf - d_contact), 0.0, 1.0))
    return float(np.sin(0.5 * np.pi * u) ** 2)


def clean_smiles(smiles: str) -> str:
    """Normalize SMILES representation removing arbitrary stereochemical flags."""
    if not smiles:
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Chem.RemoveStereochemistry(mol)
            return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        pass
    try:
        return Chem.CanonSmiles(smiles)
    except Exception:
        return smiles


def get_formal_charge(smiles: str) -> int:
    """Extract molecular formal charge from SMILES."""
    if not smiles:
        return 0
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.GetFormalCharge(mol) if mol else 0
    except Exception:
        return 0


def compute_isolated_fragment_energies(
    fragments: List[Dict[str, Any]],
    solvent: str = "water",
    relax: bool = True
) -> Dict[str, float]:
    """
    Calculate electronic ground-state potential energy (eV) of isolated fragments in solvent.
    Uses geometry relaxation (BFGS) with xTB.
    """
    from ase.optimize import BFGS

    cache: Dict[str, float] = {}
    for fr in fragments:
        smi = fr.get("smiles", "")
        if not smi:
            continue
        c_smi = clean_smiles(smi)
        if c_smi in cache:
            continue

        syms = list(fr.get("symbols", []))
        pos = list(fr.get("positions", []))
        chg = fr.get("charge", get_formal_charge(smi))

        # Fallback 3D embed if positions are absent or mismatched
        if (not syms or not pos or len(syms) != len(pos)) and c_smi:
            try:
                mol = Chem.AddHs(Chem.MolFromSmiles(c_smi))
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                syms = [a.GetSymbol() for a in mol.GetAtoms()]
                pos = mol.GetConformer().GetPositions().tolist()
                chg = Chem.GetFormalCharge(mol)
            except Exception:
                pass

        if not syms or not pos or len(syms) != len(pos):
            cache[c_smi] = 0.0
            continue

        try:
            atoms = Atoms(symbols=syms, positions=np.array(pos, dtype=float))
            atoms.calc = make_xtb(charge=chg, multiplicity=1, solvent=solvent)
            if relax and len(syms) > 1:
                try:
                    opt = BFGS(atoms, logfile=None)
                    opt.run(fmax=0.03, steps=35)
                except Exception:
                    pass
            cache[c_smi] = float(atoms.get_potential_energy())
        except Exception:
            cache[c_smi] = 0.0

    return cache


def measure_interfragment_distance(
    pos: np.ndarray,
    frag_indices_a: List[int],
    frag_indices_b: List[int]
) -> float:
    """Measure minimum inter-atomic or center-of-mass distance between two fragments."""
    if not frag_indices_a or not frag_indices_b:
        return 0.0
    sub_a = pos[frag_indices_a]
    sub_b = pos[frag_indices_b]
    dists = np.linalg.norm(sub_a[:, None, :] - sub_b[None, :, :], axis=-1)
    return float(np.min(dists))


def compute_pathway_free_energy(
    trajectory: Dict[str, Any],
    model: Any,
    solvent: str = "water",
    log_fn: Optional[Any] = None
) -> Dict[str, Any]:
    """
    Equip an IRC trajectory with rigorous, smooth Free Energy (ΔG) metrics.
    
    Parameters:
    -----------
    trajectory : dict
        Raw IRC trajectory dict containing 'frames', 'ts_frame_index', 'symbols', etc.
    model : ReactionModel
        Model defining reactants, products, bonds formed/broken, fragments, and charges.
    solvent : str
        Implicit solvent name (e.g. 'water').
    log_fn : callable, optional
        Logging callback.
        
    Returns:
    --------
    dict
        Updated trajectory with continuous 'rel_free_energy_kcal', 'free_energy_ev',
        and thermodynamically consistent forward/reverse barriers.
    """
    log = log_fn or (lambda msg: None)
    frames = trajectory.get("frames", [])
    if not frames:
        return trajectory

    # If ab initio Gibbs free energies were already computed by ORCA frequency analysis, preserve them
    if trajectory.get("barriers", {}).get("energy_type") == "gibbs" and trajectory.get("thermo_data"):
        log("  [Free Energy Engine] Trajectory already possesses ab initio Gibbs free energies (from ORCA Freq). Preserving.")
        return trajectory

    ts_idx = trajectory.get("ts_frame_index", len(frames) // 2)
    n_frames = len(frames)
    symbols = trajectory.get("symbols", getattr(model, "symbols", []))

    # Identify reactant and product molecularity
    r_mols = getattr(model, "r_mols", [])
    p_mols = getattr(model, "p_mols", [])
    n_r = max(1, len(r_mols))
    n_p = max(1, len(p_mols))

    log(f"  [Free Energy Engine] Molecularity: {n_r} reactant(s) -> [TS]‡ -> {n_p} product(s)")

    # Build fragment descriptors for isolated energy evaluation
    r_frags = []
    for fi, pm in enumerate(r_mols):
        f_idxs = [k for k in range(model.n_atoms) if model.frag_of[k] == fi]
        f_syms = [symbols[k] for k in f_idxs]
        f_pos0 = [frames[0]["positions"][k] for k in f_idxs]
        smi = model.r_smiles[fi] if fi < len(model.r_smiles) else ""
        r_frags.append({
            "smiles": smi, "symbols": f_syms, "positions": f_pos0,
            "charge": get_formal_charge(smi), "indices": f_idxs
        })

    p_frags = []
    for fi, pm in enumerate(p_mols):
        f_idxs = [k for k in range(model.n_atoms) if model.r2p.get(k, (-1, -1))[0] == fi]
        f_syms = [symbols[k] for k in f_idxs]
        f_pos_end = [frames[-1]["positions"][k] for k in f_idxs]
        smi = model.p_smiles[fi] if fi < len(model.p_smiles) else ""
        p_frags.append({
            "smiles": smi, "symbols": f_syms, "positions": f_pos_end,
            "charge": get_formal_charge(smi), "indices": f_idxs
        })

    # Evaluate isolated ground-state energies
    all_frags = r_frags + p_frags
    iso_cache = compute_isolated_fragment_energies(all_frags, solvent=solvent, relax=True)

    e_iso_r = sum(iso_cache.get(clean_smiles(f["smiles"]), 0.0) for f in r_frags)
    e_iso_p = sum(iso_cache.get(clean_smiles(f["smiles"]), 0.0) for f in p_frags)

    raw_e_list = []
    for f in frames:
        if "energy_ev" in f:
            raw_e_list.append(float(f["energy_ev"]))
        elif "energy_kcal" in f:
            raw_e_list.append(float(f["energy_kcal"]) / EV_TO_KCAL)
        elif "rel_energy_kcal" in f:
            raw_e_list.append(float(f["rel_energy_kcal"]) / EV_TO_KCAL)
        else:
            raw_e_list.append(0.0)
    raw_e = np.array(raw_e_list)
    e_ts_raw = raw_e[ts_idx]
    e_r_raw = raw_e[0]
    e_p_raw = raw_e[-1]

    # Baseline thermodynamic delta G for isolated endpoints
    has_iso = (e_iso_r != 0.0 and e_iso_p != 0.0)
    if has_iso:
        # Standard-state Hess's Law delta G between infinitely separated products and reactants
        delta_g_rxn_kcal = (e_iso_p - e_iso_r) * EV_TO_KCAL + (n_r - n_p) * DELTA_G_ASSOC_KCAL
        e_ref_baseline = e_iso_r
    else:
        delta_g_rxn_kcal = (e_p_raw - e_r_raw) * EV_TO_KCAL + (n_r - n_p) * DELTA_G_ASSOC_KCAL
        e_ref_baseline = e_r_raw

    # Contact distance parameters based on formed/broken reacting bonds
    if model.formed:
        i_f, j_f = model.formed[0]
        d0_form = model.d0(i_f, j_f)
    else:
        i_f, j_f = 0, 1
        d0_form = 1.50

    if model.broken:
        i_b, j_b = model.broken[0]
        d0_brk = model.d0(i_b, j_b)
    else:
        i_b, j_b = 0, 1
        d0_brk = 1.50

    d_contact_r = d0_form + DEFAULT_D_CONTACT_MARGIN
    d_inf_r = max(5.40, d_contact_r + DEFAULT_D_INF_MARGIN)

    d_contact_p = d0_brk + DEFAULT_D_CONTACT_MARGIN
    d_inf_p = max(5.40, d_contact_p + DEFAULT_D_INF_MARGIN)

    # Compute continuous switching progress λ_R and λ_P across frames
    lambda_r = np.zeros(n_frames)
    lambda_p = np.zeros(n_frames)

    for k in range(n_frames):
        pos = np.array(frames[k]["positions"])
        # Reactant branch (k <= ts_idx)
        if k <= ts_idx:
            if n_r > 1 and len(r_frags) >= 2:
                d_val = float(np.linalg.norm(pos[i_f] - pos[j_f])) if model.formed else measure_interfragment_distance(pos, r_frags[0]["indices"], r_frags[1]["indices"])
                lambda_r[k] = smooth_switching_function(d_val, d_contact_r, d_inf_r)
            else:
                lambda_r[k] = 0.0
        # Product branch (k >= ts_idx)
        if k >= ts_idx:
            if n_p > 1 and len(p_frags) >= 2:
                d_val = float(np.linalg.norm(pos[i_b] - pos[j_b])) if model.broken else measure_interfragment_distance(pos, p_frags[0]["indices"], p_frags[1]["indices"])
                lambda_p[k] = smooth_switching_function(d_val, d_contact_p, d_inf_p)
            else:
                lambda_p[k] = 0.0

    # Ensure TS frame has exact zero switching (pure bound complex)
    lambda_r[ts_idx] = 0.0
    lambda_p[ts_idx] = 0.0

    # Assemble smooth continuous Free Energy curve G(k)
    # At any point k:
    # 1. Reverse branch (k <= ts_idx):
    #    ΔG(k) = (E_raw(k) - e_r_raw) * EV_TO_KCAL + (1.0 - lambda_r[k]) * (n_r - 1) * DELTA_G_ASSOC_KCAL
    # 2. Forward branch (k >= ts_idx):
    #    ΔG(k) = (E_raw(k) - e_r_raw) * EV_TO_KCAL + (n_r - 1) * DELTA_G_ASSOC_KCAL - lambda_p[k] * (n_p - 1) * DELTA_G_ASSOC_KCAL
    g_kcal = np.zeros(n_frames)

    for k in range(ts_idx + 1):
        lam_r = lambda_r[k]
        g_kcal[k] = (raw_e[k] - e_r_raw) * EV_TO_KCAL + (1.0 - lam_r) * (n_r - 1) * DELTA_G_ASSOC_KCAL

    for k in range(ts_idx, n_frames):
        lam_p = lambda_p[k]
        g_kcal[k] = (raw_e[k] - e_r_raw) * EV_TO_KCAL + (n_r - 1) * DELTA_G_ASSOC_KCAL - lam_p * (n_p - 1) * DELTA_G_ASSOC_KCAL

    # Exact value match at TS saddle point
    g_kcal[ts_idx] = (raw_e[ts_idx] - e_r_raw) * EV_TO_KCAL + (n_r - 1) * DELTA_G_ASSOC_KCAL

    # Normalize relative to initial reactant state G(0) = 0.0
    rel_g_kcal = g_kcal - g_kcal[0]

    # Compute thermodynamically consistent barriers
    min_reac_g = float(np.min(rel_g_kcal[:ts_idx + 1])) if ts_idx > 0 else 0.0
    barrier_fwd = max(0.0, float(rel_g_kcal[ts_idx] - min_reac_g))
    barrier_rev = max(0.0, float(rel_g_kcal[ts_idx] - rel_g_kcal[-1]))
    rxn_g = float(rel_g_kcal[-1])

    log(f"  [Free Energy Engine] ΔG‡_fwd = {barrier_fwd:.2f} kcal/mol, ΔG‡_rev = {barrier_rev:.2f} kcal/mol, ΔG_rxn = {rxn_g:+.2f} kcal/mol")

    # Update frame records
    for k, fr in enumerate(frames):
        fr["free_energy_kcal"] = round(float(rel_g_kcal[k]), 3)
        fr["rel_free_energy_kcal"] = round(float(rel_g_kcal[k]), 3)
        # Update rel_energy_kcal to display Free Energy by default in frontend
        fr["rel_energy_kcal"] = round(float(rel_g_kcal[k]), 3)

    trajectory["barriers"].update({
        "energy_type": "gibbs",
        "forward_barrier_kcal": round(barrier_fwd, 2),
        "reverse_barrier_kcal": round(barrier_rev, 2),
        "reaction_energy_kcal": round(rxn_g, 2),
        "delta_g_assoc_applied_kcal": round((n_r - 1) * DELTA_G_ASSOC_KCAL, 2),
        "delta_g_dissoc_applied_kcal": round((n_p - 1) * DELTA_G_ASSOC_KCAL, 2)
    })
    trajectory["frames"] = frames
    return trajectory
