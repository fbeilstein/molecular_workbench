import numpy as np
from typing import Any, Dict, List, Optional
from rdkit import Chem
from ase import Atoms
from .calculators import make_xtb
from .common import rcov


def get_formal_charge(smiles: str) -> int:
    """Calculate formal charge from SMILES string using RDKit."""
    if not smiles:
        return 0
    try:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.GetFormalCharge(mol) if mol else 0
    except Exception:
        return 0


def analyze_pathway_stoichiometry(
    trajectories: List[Dict[str, Any]], fragment_roles: Optional[Dict[str, str]] = None
) -> Dict[str, Any]:
    """Map global reaction stoichiometry with automatic & user-configurable roles."""
    roles = fragment_roles or {}
    n_steps = len(trajectories)
    step_info, all_external_reactants, all_external_byproducts = [], [], []

    for k, traj in enumerate(trajectories):
        frags = traj.get("fragments", {})
        symbols = traj.get("symbols", [])
        f0_pos = traj.get("frames", [{}])[0].get("positions", [])
        flast_pos = traj.get("frames", [{}])[-1].get("positions", [])

        r_list = [
            {
                "id": f"s{k+1}_r{r_idx}", "step_index": k + 1, "formula": r.get("formula", ""),
                "smiles": Chem.CanonSmiles(r.get("smiles", "")) if r.get("smiles") else f"R_{k}_{r_idx}",
                "charge": get_formal_charge(r.get("smiles", "")),
                "num_atoms": r.get("num_atoms", len(r.get("atom_indices", []))),
                "atom_indices": r.get("atom_indices", []),
                "symbols": [symbols[i] for i in r.get("atom_indices", [])] if symbols else [],
                "positions": [f0_pos[i] for i in r.get("atom_indices", [])] if f0_pos else [],
                "role": roles.get(f"s{k+1}_r{r_idx}", roles.get(Chem.CanonSmiles(r.get("smiles", "")), "auto"))
            }
            for r_idx, r in enumerate(frags.get("reactants", []))
        ]

        p_list = [
            {
                "id": f"s{k+1}_p{p_idx}", "step_index": k + 1, "formula": p.get("formula", ""),
                "smiles": Chem.CanonSmiles(p.get("smiles", "")) if p.get("smiles") else f"P_{k}_{p_idx}",
                "charge": get_formal_charge(p.get("smiles", "")),
                "num_atoms": p.get("num_atoms", len(p.get("atom_indices", []))),
                "atom_indices": p.get("atom_indices", []),
                "symbols": [symbols[i] for i in p.get("atom_indices", [])] if symbols else [],
                "positions": [flast_pos[i] for i in p.get("atom_indices", [])] if flast_pos else [],
                "role": roles.get(f"s{k+1}_p{p_idx}", roles.get(Chem.CanonSmiles(p.get("smiles", "")), "auto"))
            }
            for p_idx, p in enumerate(frags.get("products", []))
        ]

        step_info.append({"step_index": k + 1, "reactants": r_list, "products": p_list})

    # Forward resolution: resolve core vs incoming
    produced = {}
    for k, s in enumerate(step_info):
        s["incoming_external"] = []
        max_r_atoms = max((r["num_atoms"] for r in s["reactants"]), default=0)
        for r in s["reactants"]:
            if r["role"] == "auto":
                if r["smiles"] in produced and produced[r["smiles"]] > 0:
                    produced[r["smiles"]] -= 1
                    r["role"] = "core"
                elif k == 0:
                    r["role"] = "core" if (r["num_atoms"] == max_r_atoms and len(s["reactants"]) > 1) else "incoming"
                else:
                    r["role"] = "incoming"
            if r["role"] == "incoming":
                s["incoming_external"].append(r)
                all_external_reactants.append(r)
        for p in s["products"]:
            produced[p["smiles"]] = produced.get(p["smiles"], 0) + 1

    # Backward resolution: resolve core vs leaving
    consumed = {}
    for k in reversed(range(n_steps)):
        s = step_info[k]
        s["departing_byproducts"] = []
        max_p_atoms = max((p["num_atoms"] for p in s["products"]), default=0)
        for p in s["products"]:
            if p["role"] == "auto":
                if p["smiles"] in consumed and consumed[p["smiles"]] > 0:
                    consumed[p["smiles"]] -= 1
                    p["role"] = "core"
                elif k == n_steps - 1:
                    p["role"] = "core" if (p["num_atoms"] == max_p_atoms and len(s["products"]) > 1) else "leaving"
                else:
                    p["role"] = "leaving"
            if p["role"] == "leaving":
                s["departing_byproducts"].append(p)
                all_external_byproducts.append(p)
        for r in s["reactants"]:
            consumed[r["smiles"]] = consumed.get(r["smiles"], 0) + 1

    for k in range(n_steps):
        step_info[k]["waiting_reagents"] = [r for m in range(k + 1, n_steps) for r in step_info[m]["incoming_external"]]
        step_info[k]["departed_byproducts"] = [p for m in range(0, k) for p in step_info[m]["departing_byproducts"]]

    return {"steps": step_info, "all_external_reactants": all_external_reactants, "all_external_byproducts": all_external_byproducts}


def compute_isolated_fragment_energies(fragments: List[Dict[str, Any]], solvent: str = "water") -> Dict[str, float]:
    """Compute isolated electronic energy (eV) for fragments in chosen solvent using xTB."""
    cache = {}
    for fr in fragments:
        smi = fr["smiles"]
        if smi in cache:
            continue
        syms, pos, chg = fr.get("symbols", []), fr.get("positions", []), fr.get("charge", 0)
        if not syms or not pos or len(syms) != len(pos):
            cache[smi] = 0.0
            continue
        try:
            atoms = Atoms(symbols=syms, positions=pos)
            atoms.calc = make_xtb(charge=chg, multiplicity=1, solvent=solvent)
            cache[smi] = float(atoms.get_potential_energy())
        except Exception:
            cache[smi] = 0.0
    return cache


def _smoothstep(x: float) -> float:
    x = float(np.clip(x, 0.0, 1.0))
    return x * x * (3.0 - 2.0 * x)


def _hermite_spline(t: np.ndarray, y0: float, y1: float, g0: float, g1: float) -> np.ndarray:
    h00 = 2.0 * t**3 - 3.0 * t**2 + 1.0
    h10 = t**3 - 2.0 * t**2 + t
    h01 = -2.0 * t**3 + 3.0 * t**2
    h11 = t**3 - t**2
    return h00 * y0 + h10 * g0 + h01 * y1 + h11 * g1


def balance_pathway_energies(
    megajson: Dict[str, Any],
    trajectories: List[Dict[str, Any]],
    recalculate_with_xtb: bool = True,
    solvent: str = "water",
    fragment_roles: Optional[Dict[str, str]] = None,
    truncate_at_minimum: bool = False
) -> Dict[str, Any]:
    """Calculate Hess's Law balanced Free Energy (G) profile with dynamic association entropy and leaving group exclusion."""
    stoich = analyze_pathway_stoichiometry(trajectories, fragment_roles)
    all_ext = stoich["all_external_reactants"] + stoich["all_external_byproducts"]
    frag_energies = compute_isolated_fragment_energies(all_ext, solvent=solvent)

    boundaries = megajson.get("step_boundaries", [])
    master_frames = megajson.get("frames", [])
    symbols = megajson.get("symbols", [])
    n_frames = len(master_frames)
    if n_frames == 0:
        return megajson

    # 1. Base electronic energy in chosen solvent
    active_e = np.zeros(n_frames)
    if recalculate_with_xtb and symbols:
        init_chg = sum(r.get("charge", 0) for r in stoich["steps"][0]["reactants"])
        calc = make_xtb(charge=init_chg, multiplicity=1, solvent=solvent)
        atoms = Atoms(symbols=symbols)
        atoms.calc = calc
        for i, fr in enumerate(master_frames):
            atoms.set_positions(fr.get("positions", fr.get("coords", [])))
            active_e[i] = float(atoms.get_potential_energy())
    else:
        for i, fr in enumerate(master_frames):
            active_e[i] = float(fr.get("energy_ev", fr.get("energy_kcal", 0.0) / 23.0605))

    delta_g_assoc = 11.87 / 23.0605  # Bimolecular association entropy penalty (~11.87 kcal/mol in eV)
    balanced_ev = np.copy(active_e)

    # 2. Process each elementary reaction step with dynamic free energy & fragment isolation
    for k, b in enumerate(boundaries):
        s_info = stoich["steps"][k]
        start, end = b["start_frame"], b["end_frame"]
        traj = trajectories[k] if k < len(trajectories) else {}
        bonds = traj.get("bonds", [])

        # Hess's Law offset for waiting reagents and past departed byproducts
        s_offset = sum(frag_energies.get(r["smiles"], 0.0) for r in s_info["waiting_reagents"])
        s_offset += sum(frag_energies.get(p["smiles"], 0.0) for p in s_info["departed_byproducts"])

        # Check for incoming fragment in this step
        inc_frags = s_info.get("incoming_external", [])
        forming_bonds = [bd for bd in bonds if bd.get("type") == "forming"]
        # Check for leaving fragment in this step
        leaving_frags = s_info.get("departing_byproducts", [])
        breaking_bonds = [bd for bd in bonds if bd.get("type") == "breaking"]

        # Case A: Associative Step (incoming reagent arrives)
        if inc_frags and forming_bonds:
            bi, bj = forming_bonds[0]["from"], forming_bonds[0]["to"]
            d0_form = rcov(symbols[bi]) + rcov(symbols[bj])
            d_contact, d_cov = d0_form + 1.25, d0_form * 1.15
            for idx in range(start, end + 1):
                pos = np.array(master_frames[idx].get("coords", master_frames[idx].get("positions", [])))
                d_val = float(np.linalg.norm(pos[bi] - pos[bj]))
                assoc_w = _smoothstep((d_contact - d_val) / max(0.2, (d_contact - d_cov)))
                balanced_ev[idx] = active_e[idx] + delta_g_assoc * assoc_w + s_offset

        # Case B: Dissociative Step (leaving group departs)
        elif leaving_frags and breaking_bonds:
            bi, bj = breaking_bonds[0]["from"], breaking_bonds[0]["to"]
            d0_break = rcov(symbols[bi]) + rcov(symbols[bj])
            d_tight, d_free = d0_break * 1.25, d0_break * 1.65

            leaving_atom_set = set(leaving_frags[0].get("atom_indices", []))
            if not leaving_atom_set:
                leaving_atom_set = {bj}
            hero_idxs = [i for i in range(len(symbols)) if i not in leaving_atom_set]

            q_leaving = leaving_frags[0].get("charge", 0)
            tot_step_chg = sum(r.get("charge", 0) for r in s_info["reactants"])
            q_hero = tot_step_chg - q_leaving

            hero_atoms = Atoms(symbols=[symbols[i] for i in hero_idxs])
            hero_calc = make_xtb(charge=q_hero, multiplicity=1, solvent=solvent)
            hero_atoms.calc = hero_calc
            e_leaving_iso = frag_energies.get(leaving_frags[0]["smiles"], 0.0)

            for idx in range(start, end + 1):
                pos = np.array(master_frames[idx].get("coords", master_frames[idx].get("positions", [])))
                d_val = float(np.linalg.norm(pos[bi] - pos[bj]))
                e_bound = active_e[idx] + delta_g_assoc
                hero_atoms.set_positions(pos[hero_idxs])
                e_hero = float(hero_atoms.get_potential_energy())
                e_decoupled = e_hero + e_leaving_iso

                if d_val <= d_tight:
                    balanced_ev[idx] = e_bound + s_offset
                elif d_val >= d_free:
                    balanced_ev[idx] = e_decoupled + s_offset
                else:
                    w = _smoothstep((d_val - d_tight) / max(0.1, (d_free - d_tight)))
                    balanced_ev[idx] = ((1.0 - w) * e_bound + w * e_decoupled) + s_offset

        # Case C: Isomerization / Rearrangement
        else:
            balanced_ev[start : end + 1] += s_offset

    # 3. Continuous inter-step bridge alignment
    for k in range(len(boundaries) - 1):
        b_curr, b_next = boundaries[k], boundaries[k + 1]
        shift = balanced_ev[b_curr["end_frame"]] - balanced_ev[b_next["start_frame"]]
        balanced_ev[b_next["start_frame"] : b_next["end_frame"] + 1] += shift

        b_start, b_end = b_curr["end_frame"] + 1, b_next["start_frame"] - 1
        if b_start <= b_end:
            n_b = b_end - b_start + 1
            y0, y1 = balanced_ev[b_curr["end_frame"]], balanced_ev[b_next["start_frame"]]
            g0 = balanced_ev[b_curr["end_frame"]] - balanced_ev[max(0, b_curr["end_frame"] - 1)]
            g1 = balanced_ev[min(n_frames - 1, b_next["start_frame"] + 1)] - balanced_ev[b_next["start_frame"]]
            tau = np.linspace(0, 1, n_b + 2)[1:-1]
            balanced_ev[b_start : b_end + 1] = _hermite_spline(tau, y0, y1, g0, g1)

    # 4. Guarantee product refinement: exit channel monotonically descends into product minimum
    last_b = boundaries[-1]
    s_last_start, s_last_end = last_b["start_frame"], last_b["end_frame"]
    min_idx = s_last_start + int(np.argmin(balanced_ev[s_last_start : s_last_end + 1]))
    min_val = balanced_ev[min_idx]
    for idx in range(min_idx, s_last_end + 1):
        balanced_ev[idx] = min(balanced_ev[idx], min_val)

    ref_ev = balanced_ev[0]
    rel_kcal = (balanced_ev - ref_ev) * 23.0605

    step_metrics = [
        {
            "step_index": k + 1, "step_name": b.get("step_name", f"Step {k+1}"),
            "start_frame": b["start_frame"], "ts_frame": b["ts_frame"], "end_frame": b["end_frame"],
            "barrier_kcal": round(float(max(0.0, np.max(rel_kcal[b["start_frame"]:b["end_frame"]+1]) - rel_kcal[b["start_frame"]])), 2),
            "reaction_energy_kcal": round(float(rel_kcal[b["end_frame"]] - rel_kcal[b["start_frame"]]), 2),
            "ts_energy_kcal": round(float(rel_kcal[b["ts_frame"]]), 2),
            "reactant_energy_kcal": round(float(rel_kcal[b["start_frame"]]), 2),
            "product_energy_kcal": round(float(rel_kcal[b["end_frame"]]), 2)
        }
        for k, b in enumerate(boundaries)
    ]
    for b, m in zip(boundaries, step_metrics):
        b.update(m)

    for i, fr in enumerate(master_frames):
        fr["energy_kcal"], fr["balanced_energy_ev"] = round(float(rel_kcal[i]), 3), round(float(balanced_ev[i]), 5)

    megajson["overall_barrier_kcal"] = round(float(max(np.max(rel_kcal), 0.0)), 2)
    megajson["overall_reaction_energy_kcal"] = round(float(rel_kcal[-1]), 2)
    megajson["step_boundaries"] = boundaries
    megajson["energy_profile"] = {
        "frame_energies_kcal": rel_kcal.round(3).tolist(),
        "step_metrics": step_metrics,
        "overall_barrier_kcal": megajson["overall_barrier_kcal"],
        "overall_reaction_energy_kcal": megajson["overall_reaction_energy_kcal"],
        "initial_reactants": [r["formula"] for r in stoich["all_external_reactants"]],
        "final_products": [p["formula"] for p in stoich["all_external_byproducts"]]
    }
    megajson["stoichiometry_summary"], megajson["solvent"] = stoich, solvent
    return megajson
