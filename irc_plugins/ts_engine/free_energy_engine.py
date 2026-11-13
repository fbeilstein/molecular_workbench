"""
Free Energy and Thermodynamics Engine for Chemical Reactions.

Provides rigorous, physically sound solution-phase Gibbs free energy (ΔG°_soln)
profiles based on Grimme qRRHO thermochemistry, variational transition state
alignment, and barrierless reaction coordinate handling.
"""

from typing import Any, Dict, List, Optional
import numpy as np
from ase import Atoms

from .common import EV_TO_KCAL
from .thermo import (
    compute_species_gibbs_free_energy,
    split_into_independent_fragments,
    compute_fragment_charges_and_multiplicities,
    STANDARD_STATE_COMPRESSION_KCAL,
)

# Interaction distance margins for smooth fragment approach/departure (Angstrom)
DEFAULT_D_CONTACT_MARGIN: float = 0.60
DEFAULT_D_INF_MARGIN: float = 3.20


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


def measure_interfragment_distance(
    pos: np.ndarray,
    frag_indices_a: List[int],
    frag_indices_b: List[int]
) -> float:
    """Measure minimum inter-atomic distance between two fragments."""
    if not frag_indices_a or not frag_indices_b:
        return 0.0
    sub_a = pos[frag_indices_a]
    sub_b = pos[frag_indices_b]
    dists = np.linalg.norm(sub_a[:, None, :] - sub_b[None, :, :], axis=-1)
    return float(np.min(dists))


def compute_pathway_free_energy(
    trajectory: Dict[str, Any],
    model: Optional[Any] = None,
    solvent: str = "water",
    log_fn: Optional[Any] = None
) -> Dict[str, Any]:
    """
    Equip an IRC trajectory with rigorous solution-phase Gibbs Free Energy (ΔG°_soln) metrics,
    relocate ts_frame_index to the true variational ΔG maximum (or mark as barrierless),
    and guarantee continuity across the reaction coordinate s.
    """
    log = log_fn or (lambda msg: None)
    frames = trajectory.get("frames", [])
    if not frames:
        return trajectory

    n_frames = len(frames)
    symbols = trajectory.get("symbols") or (getattr(model, "symbols", None) if model else None)
    if not symbols and frames and "symbols" in frames[0]:
        symbols = frames[0]["symbols"]
    if not symbols and frames and "positions" in frames[0]:
        symbols = ["C"] * len(frames[0]["positions"])

    chem_start = int(trajectory.get("approach_frames") or trajectory.get("asymp_r_count") or 0)
    chem_end = n_frames - int(trajectory.get("departure_frames") or trajectory.get("asymp_p_count") or 0) - 1
    chem_start = max(0, min(chem_start, n_frames - 1))
    chem_end = max(chem_start, min(chem_end, n_frames - 1))

    # Extract electronic potential energies across trajectory
    raw_e_ev = []
    for f in frames:
        if "energy_ev" in f:
            raw_e_ev.append(float(f["energy_ev"]))
        elif "energy_kcal" in f:
            raw_e_ev.append(float(f["energy_kcal"]) / EV_TO_KCAL)
        elif "rel_energy_kcal" in f:
            raw_e_ev.append(float(f["rel_energy_kcal"]) / EV_TO_KCAL)
        else:
            raw_e_ev.append(0.0)
    raw_e = np.array(raw_e_ev)

    sella_idx = trajectory.get("sella_ts_frame_index")
    ts_idx = trajectory.get("ts_frame_index")
    if sella_idx is not None and 0 <= sella_idx < n_frames:
        orig_ts_idx = int(sella_idx)
    elif ts_idx is not None and 0 <= ts_idx < n_frames:
        orig_ts_idx = int(ts_idx)
    elif len(raw_e) > 0 and chem_end > chem_start:
        orig_ts_idx = int(chem_start + np.argmax(raw_e[chem_start:chem_end + 1]))
    else:
        orig_ts_idx = n_frames // 2

    # Guard against corrupted midpoint defaults when a clear electronic peak exists in chemical region
    if orig_ts_idx == n_frames // 2 and len(raw_e) > 0 and chem_end > chem_start:
        e_peak_idx = int(chem_start + np.argmax(raw_e[chem_start:chem_end + 1]))
        if (raw_e[e_peak_idx] - raw_e[orig_ts_idx]) * EV_TO_KCAL > 0.05:
            orig_ts_idx = e_peak_idx

    orig_ts_idx = int(np.clip(orig_ts_idx, 0, n_frames - 1))
    trajectory["sella_ts_frame_index"] = orig_ts_idx

    charge = int(getattr(model, "charge", 0) if model else trajectory.get("charge", 0))
    mult = int(getattr(model, "multiplicity", 1) if model else trajectory.get("multiplicity", 1))

    e_r_raw = raw_e[0]
    e_ts_raw = raw_e[orig_ts_idx]
    e_p_raw = raw_e[-1]

    # Identify reactant and product molecularity & fragments
    pos_r = np.array(frames[0]["positions"])
    pos_ts = np.array(frames[orig_ts_idx]["positions"])
    pos_p = np.array(frames[-1]["positions"])

    atoms_r = Atoms(symbols=symbols, positions=pos_r)
    atoms_ts = Atoms(symbols=symbols, positions=pos_ts)
    atoms_p = Atoms(symbols=symbols, positions=pos_p)

    r_mols = getattr(model, "r_mols", None) if model else None
    p_mols = getattr(model, "p_mols", None) if model else None

    # Fragment identification
    if r_mols and len(r_mols) > 1 and hasattr(model, "frag_of"):
        r_frags = [[k for k in range(len(symbols)) if model.frag_of[k] == fi] for fi in range(len(r_mols))]
    else:
        r_frags = split_into_independent_fragments(atoms_r)

    if p_mols and len(p_mols) > 1 and hasattr(model, "r2p"):
        p_frags = [[k for k in range(len(symbols)) if model.r2p.get(k, (-1, -1))[0] == fi] for fi in range(len(p_mols))]
    else:
        p_frags = split_into_independent_fragments(atoms_p)

    n_r = len(r_frags)
    n_p = len(p_frags)
    log(f"  [Free Energy Engine] Molecularity: {n_r} reactant fragment(s) -> supermolecule -> {n_p} product fragment(s)")

    # 1. Compute rigorous qRRHO free energy of TS supermolecule
    thermo_ts = compute_species_gibbs_free_energy(
        atoms_ts, charge=charge, mult=mult, solvent=solvent, is_ts=True,
        e_el_override_kcal=e_ts_raw * EV_TO_KCAL
    )
    g_corr_ts = thermo_ts.get("G_corr_kcal", 0.0)

    # 2. Compute qRRHO free energy of Reactant state
    if n_r > 1:
        frag_props_r = compute_fragment_charges_and_multiplicities(atoms_r, r_frags, charge, mult, solvent=solvent)
        g_corr_r_sum = 0.0
        for f_idx, (f_chg, f_mult) in zip(r_frags, frag_props_r):
            sub_atoms = Atoms(symbols=[symbols[i] for i in f_idx], positions=pos_r[f_idx])
            sub_thermo = compute_species_gibbs_free_energy(sub_atoms, charge=f_chg, mult=f_mult, solvent=solvent, is_ts=False)
            g_corr_r_sum += sub_thermo.get("G_corr_kcal", 0.0)
        g_corr_r = g_corr_r_sum
    else:
        thermo_r = compute_species_gibbs_free_energy(
            atoms_r, charge=charge, mult=mult, solvent=solvent, is_ts=False,
            e_el_override_kcal=e_r_raw * EV_TO_KCAL
        )
        g_corr_r = thermo_r.get("G_corr_kcal", 0.0)

    # 3. Compute qRRHO free energy of Product state
    if n_p > 1:
        frag_props_p = compute_fragment_charges_and_multiplicities(atoms_p, p_frags, charge, mult, solvent=solvent)
        g_corr_p_sum = 0.0
        for f_idx, (f_chg, f_mult) in zip(p_frags, frag_props_p):
            sub_atoms = Atoms(symbols=[symbols[i] for i in f_idx], positions=pos_p[f_idx])
            sub_thermo = compute_species_gibbs_free_energy(sub_atoms, charge=f_chg, mult=f_mult, solvent=solvent, is_ts=False)
            g_corr_p_sum += sub_thermo.get("G_corr_kcal", 0.0)
        g_corr_p = g_corr_p_sum
    else:
        thermo_p = compute_species_gibbs_free_energy(
            atoms_p, charge=charge, mult=mult, solvent=solvent, is_ts=False,
            e_el_override_kcal=e_p_raw * EV_TO_KCAL
        )
        g_corr_p = thermo_p.get("G_corr_kcal", 0.0)

    # 4. Construct continuous G_corr(k) along the path
    bonds = trajectory.get("bonds", [])
    forming_pairs = [(b["from"], b["to"]) for b in bonds if b.get("type") == "forming"]
    breaking_pairs = [(b["from"], b["to"]) for b in bonds if b.get("type") == "breaking"]
    if not forming_pairs and model and hasattr(model, "formed") and model.formed:
        forming_pairs = list(model.formed)
    if not breaking_pairs and model and hasattr(model, "broken") and model.broken:
        breaking_pairs = list(model.broken)

    if forming_pairs:
        active_pair = forming_pairs[0]
    elif breaking_pairs:
        active_pair = breaking_pairs[0]
    elif n_r > 1 and r_frags and len(r_frags) > 1:
        sub_a = pos_ts[r_frags[0]]
        sub_b = pos_ts[r_frags[1]]
        dists = np.linalg.norm(sub_a[:, None, :] - sub_b[None, :, :], axis=-1)
        min_idx = np.unravel_index(np.argmin(dists), dists.shape)
        active_pair = (r_frags[0][min_idx[0]], r_frags[1][min_idx[1]])
    elif n_p > 1 and p_frags and len(p_frags) > 1:
        sub_a = pos_ts[p_frags[0]]
        sub_b = pos_ts[p_frags[1]]
        dists = np.linalg.norm(sub_a[:, None, :] - sub_b[None, :, :], axis=-1)
        min_idx = np.unravel_index(np.argmin(dists), dists.shape)
        active_pair = (p_frags[0][min_idx[0]], p_frags[1][min_idx[1]])
    else:
        active_pair = (0, 1)

    i_act, j_act = active_pair
    d_ts = float(np.linalg.norm(pos_ts[i_act] - pos_ts[j_act]))
    d_p = float(np.linalg.norm(pos_p[i_act] - pos_p[j_act]))
    d_inf = max(5.00, d_ts + 2.50)

    # Build continuous G(k) = E_el(k) + G_corr(k)
    g_kcal = np.zeros(n_frames)
    for k in range(n_frames):
        pos_k = np.array(frames[k]["positions"])
        d_val = float(np.linalg.norm(pos_k[i_act] - pos_k[j_act]))
        e_el_k = (raw_e[k] - e_r_raw) * EV_TO_KCAL

        if k <= orig_ts_idx:
            if n_r > 1:
                u = np.clip((d_inf - d_val) / max(0.1, d_inf - d_ts), 0.0, 1.0)
                w = np.sin(0.5 * np.pi * u) ** 2
            else:
                u = k / max(1, orig_ts_idx)
                w = np.sin(0.5 * np.pi * u) ** 2
            g_corr_k = g_corr_r + w * (g_corr_ts - g_corr_r)
        else:
            if n_p > 1:
                u = np.clip((d_val - d_ts) / max(0.1, d_inf - d_ts), 0.0, 1.0)
                w = np.sin(0.5 * np.pi * u) ** 2
            else:
                denom = abs(d_ts - d_p)
                if denom > 0.1:
                    u = np.clip((d_ts - d_val) / denom, 0.0, 1.0)
                else:
                    u = (k - orig_ts_idx) / max(1, n_frames - 1 - orig_ts_idx)
                w = np.sin(0.5 * np.pi * u) ** 2
            g_corr_k = g_corr_ts + w * (g_corr_p - g_corr_ts)

        g_kcal[k] = e_el_k + (g_corr_k - g_corr_r)

    # Re-normalize so reactant baseline G(0) = 0.0
    rel_g_kcal = g_kcal - g_kcal[0]

    # 5. Variational TS Relocation & Barrierless Step Detection
    # Search for peak in chemical region
    k_star = int(chem_start + np.argmax(rel_g_kcal[chem_start:chem_end + 1]))
    g_peak = float(rel_g_kcal[k_star])
    rxn_g = float(rel_g_kcal[-1])

    min_reac_g = float(np.min(rel_g_kcal[:k_star + 1])) if k_star > 0 else 0.0
    barrier_fwd = max(0.0, float(g_peak - min_reac_g))
    barrier_rev = max(0.0, float(g_peak - rxn_g))

    # A step is only truly barrierless if it is downhill/thermoneutral AND has no significant forward barrier
    is_barrierless = bool(barrier_fwd <= 0.50 and rxn_g <= 0.50)

    # Reaction coordinate s_arc alignment
    orig_s = [float(f.get("s_arc", float(k - orig_ts_idx))) for k, f in enumerate(frames)]

    if not is_barrierless:
        # Activated step: relocate ts_frame_index to true variational maximum k_star
        ts_frame_idx = k_star
        trajectory["is_barrierless"] = False
        trajectory["ts_frame_index"] = ts_frame_idx
        trajectory["sella_ts_frame_index"] = orig_ts_idx

        # Re-zero s_arc at true TS peak
        s_offset = orig_s[ts_frame_idx]
        new_s = [round(s - s_offset, 2) for s in orig_s]

        min_reac_g = float(np.min(rel_g_kcal[:ts_frame_idx + 1])) if ts_frame_idx > 0 else 0.0
        barrier_fwd = max(0.0, float(rel_g_kcal[ts_frame_idx] - min_reac_g))
        barrier_rev = max(0.0, float(rel_g_kcal[ts_frame_idx] - rel_g_kcal[-1]))
        rxn_g = float(rel_g_kcal[-1])
        log(f"  [Free Energy Engine] Activated Step: Variational TS relocated to frame {ts_frame_idx} (was {orig_ts_idx}).")
        log(f"  [Free Energy Engine] ΔG‡_fwd = {barrier_fwd:.2f} kcal/mol, ΔG‡_rev = {barrier_rev:.2f} kcal/mol, ΔG_rxn = {rxn_g:+.2f} kcal/mol")
    else:
        # Barrierless / diffusion-limited step: no artificial saddle peak
        ts_frame_idx = None
        trajectory["is_barrierless"] = True
        trajectory["ts_frame_index"] = None
        trajectory["sella_ts_frame_index"] = orig_ts_idx

        # Re-zero s_arc at frame 0
        s_offset = orig_s[0]
        new_s = [round(s - s_offset, 2) for s in orig_s]

        barrier_fwd = max(0.0, float(rel_g_kcal[-1]))
        barrier_rev = max(0.0, float(-rel_g_kcal[-1]))
        rxn_g = float(rel_g_kcal[-1])
        log(f"  [Free Energy Engine] Step is BARRIERLESS (ΔG‡_fwd = {barrier_fwd:.2f} <= 0.5 kcal/mol, ΔG_rxn = {rxn_g:+.2f} <= 0.5 kcal/mol). Suppressing fake TS.")
        log(f"  [Free Energy Engine] ΔG‡_fwd = {barrier_fwd:.2f} kcal/mol, ΔG_rxn = {rxn_g:+.2f} kcal/mol")

    # Update frame records
    for k, fr in enumerate(frames):
        fr["free_energy_kcal"] = round(float(rel_g_kcal[k]), 3)
        fr["rel_free_energy_kcal"] = round(float(rel_g_kcal[k]), 3)
        fr["rel_energy_kcal"] = round(float(rel_g_kcal[k]), 3)
        fr["s_arc"] = new_s[k]
        if ts_frame_idx is not None:
            fr["is_ts"] = bool(k == ts_frame_idx)
        else:
            fr["is_ts"] = False

    trajectory["barriers"].update({
        "energy_type": "gibbs",
        "forward_barrier_kcal": round(barrier_fwd, 2),
        "reverse_barrier_kcal": round(barrier_rev, 2),
        "reaction_energy_kcal": round(rxn_g, 2),
        "is_barrierless": is_barrierless,
        "delta_g_corr_ts_kcal": round(g_corr_ts, 2),
        "delta_g_corr_r_kcal": round(g_corr_r, 2),
        "delta_g_corr_p_kcal": round(g_corr_p, 2),
    })
    trajectory["frames"] = frames
    return trajectory
