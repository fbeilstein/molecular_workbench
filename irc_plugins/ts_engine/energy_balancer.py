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


def _clean_smiles(smiles: str) -> str:
    """Normalize SMILES without stereochemistry for robust intermediate matching."""
    if not smiles:
        return ""
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol:
            Chem.RemoveStereochemistry(mol)
            return Chem.MolToSmiles(mol, isomericSmiles=False)
    except Exception:
        pass
    return Chem.CanonSmiles(smiles) if smiles else ""


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
                "smiles": _clean_smiles(r.get("smiles", "")) or f"R_{k}_{r_idx}",
                "charge": get_formal_charge(r.get("smiles", "")),
                "num_atoms": r.get("num_atoms", len(r.get("atom_indices", []))),
                "atom_indices": r.get("atom_indices", []),
                "symbols": [symbols[i] for i in r.get("atom_indices", [])] if symbols else [],
                "positions": [f0_pos[i] for i in r.get("atom_indices", [])] if f0_pos else [],
                "role": roles.get(f"s{k+1}_r{r_idx}", roles.get(_clean_smiles(r.get("smiles", "")), "auto"))
            }
            for r_idx, r in enumerate(frags.get("reactants", []))
        ]

        p_list = [
            {
                "id": f"s{k+1}_p{p_idx}", "step_index": k + 1, "formula": p.get("formula", ""),
                "smiles": _clean_smiles(p.get("smiles", "")) or f"P_{k}_{p_idx}",
                "charge": get_formal_charge(p.get("smiles", "")),
                "num_atoms": p.get("num_atoms", len(p.get("atom_indices", []))),
                "atom_indices": p.get("atom_indices", []),
                "symbols": [symbols[i] for i in p.get("atom_indices", [])] if symbols else [],
                "positions": [flast_pos[i] for i in p.get("atom_indices", [])] if flast_pos else [],
                "role": roles.get(f"s{k+1}_p{p_idx}", roles.get(_clean_smiles(p.get("smiles", "")), "auto"))
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


def compute_isolated_fragment_energies(
    fragments: List[Dict[str, Any]], solvent: str = "water", relax: bool = True
) -> Dict[str, float]:
    """Compute isolated electronic energy (eV) for fully relaxed ground-state fragments in chosen solvent using xTB."""
    from ase.optimize import BFGS
    from rdkit import Chem
    from rdkit.Chem import AllChem

    cache = {}
    for fr in fragments:
        smi = fr.get("smiles", "")
        if not smi:
            continue
        clean_s = _clean_smiles(smi)
        if clean_s in cache:
            continue

        syms = list(fr.get("symbols", []))
        pos = list(fr.get("positions", []))
        chg = fr.get("charge", 0)

        # If 3D coordinates are missing, embed a clean 3D conformer with RDKit
        if (not syms or not pos or len(syms) != len(pos)) and clean_s:
            try:
                mol = Chem.AddHs(Chem.MolFromSmiles(clean_s))
                AllChem.EmbedMolecule(mol, AllChem.ETKDGv3())
                syms = [a.GetSymbol() for a in mol.GetAtoms()]
                pos = mol.GetConformer().GetPositions().tolist()
                chg = Chem.GetFormalCharge(mol)
            except Exception:
                pass

        if not syms or not pos or len(syms) != len(pos):
            cache[clean_s] = 0.0
            continue

        try:
            atoms = Atoms(symbols=syms, positions=np.array(pos, dtype=float))
            atoms.calc = make_xtb(charge=chg, multiplicity=1, solvent=solvent)
            if relax and len(syms) > 1:
                try:
                    opt = BFGS(atoms, logfile=None)
                    opt.run(fmax=0.03, steps=40)
                except Exception:
                    pass
            cache[clean_s] = float(atoms.get_potential_energy())
        except Exception:
            try:
                atoms = Atoms(symbols=syms, positions=np.array(pos, dtype=float))
                atoms.calc = make_xtb(charge=chg, multiplicity=1, solvent=solvent)
                cache[clean_s] = float(atoms.get_potential_energy())
            except Exception:
                cache[clean_s] = 0.0
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
    init_reactants = stoich["steps"][0]["reactants"] if stoich["steps"] else []
    all_ext = stoich["all_external_reactants"] + stoich["all_external_byproducts"] + init_reactants
    frag_energies = compute_isolated_fragment_energies(all_ext, solvent=solvent, relax=True)

    boundaries = megajson.get("step_boundaries", [])
    master_frames = megajson.get("frames", [])
    symbols = megajson.get("symbols", [])
    n_frames = len(master_frames)
    if n_frames == 0:
        return megajson

    # Combined translational/rotational entropy loss + 1M solution standard-state correction (~13.76 kcal/mol)
    delta_g_assoc = 13.76 / 23.0605

    if not recalculate_with_xtb:
        rel_kcal = np.array([float(fr.get("rel_energy_kcal", fr.get("energy_kcal", 0.0))) for fr in master_frames])
        balanced_ev = rel_kcal / 23.0605
        # Apply Hess offsets and association penalty even if using pre-computed energies
        for k, b in enumerate(boundaries):
            s_info = stoich["steps"][k]
            start, end = b["start_frame"], b["end_frame"]
            s_offset = sum(frag_energies.get(_clean_smiles(r["smiles"]), 0.0) for r in s_info["waiting_reagents"])
            s_offset += sum(frag_energies.get(_clean_smiles(p["smiles"]), 0.0) for p in s_info["departed_byproducts"])
            balanced_ev[start : end + 1] += s_offset
            if s_info.get("incoming_external"):
                ts_f = b.get("ts_frame", (start + end) // 2)
                balanced_ev[ts_f] += delta_g_assoc
    else:
        # 1. Base electronic energy normalized from IRC trajectory frames
        active_e = np.zeros(n_frames)
        for k, b in enumerate(boundaries):
            start, end = b["start_frame"], b["end_frame"]
            traj = trajectories[k] if k < len(trajectories) else {}
            e_base = traj.get("barriers", {}).get("e_reactant_ev", 0.0)
            for idx in range(start, min(end + 1, n_frames)):
                fr = master_frames[idx]
                e_fr = float(fr.get("energy_ev", 0.0))
                if abs(e_fr) < 100.0 and abs(e_base) > 100.0:
                    active_e[idx] = e_base + e_fr
                else:
                    active_e[idx] = e_fr

        # Interpolate active_e for inter-step bridge frames
        for k in range(len(boundaries) - 1):
            b_curr, b_next = boundaries[k], boundaries[k + 1]
            b_start, b_end = b_curr["end_frame"] + 1, b_next["start_frame"] - 1
            if b_start <= b_end:
                n_b = b_end - b_start + 1
                e0 = active_e[b_curr["end_frame"]]
                e1 = active_e[b_next["start_frame"]]
                tau = np.linspace(0, 1, n_b + 2)[1:-1]
                active_e[b_start : b_end + 1] = (1.0 - tau) * e0 + tau * e1

        balanced_ev = np.copy(active_e)

        # 2. Process each elementary reaction step with dynamic free energy & fragment isolation
        for k, b in enumerate(boundaries):
            s_info = stoich["steps"][k]
            start, end = b["start_frame"], b["end_frame"]
            traj = trajectories[k] if k < len(trajectories) else {}
            bonds = traj.get("bonds", [])
            step_symbols = traj.get("symbols", [])

            # Hess's Law offset for waiting reagents and past departed byproducts
            s_offset = sum(frag_energies.get(_clean_smiles(r["smiles"]), 0.0) for r in s_info["waiting_reagents"])
            s_offset += sum(frag_energies.get(_clean_smiles(p["smiles"]), 0.0) for p in s_info["departed_byproducts"])

            # Check for incoming fragment in this step
            inc_frags = s_info.get("incoming_external", [])
            forming_bonds = [bd for bd in bonds if bd.get("type") == "forming"]
            # Check for leaving fragment in this step
            leaving_frags = s_info.get("departing_byproducts", [])
            breaking_bonds = [bd for bd in bonds if bd.get("type") == "breaking"]

            # Case A: Associative Step (incoming reagent arrives)
            if inc_frags and forming_bonds and step_symbols:
                bi, bj = forming_bonds[0]["from"], forming_bonds[0]["to"]
                d0_form = rcov(step_symbols[bi]) + rcov(step_symbols[bj])
                d_contact, d_cov = d0_form + 1.35, d0_form * 1.15

                inc_atom_set = set(inc_frags[0].get("atom_indices", []))
                if not inc_atom_set:
                    inc_atom_set = {bj}
                core_idxs = [i for i in range(len(step_symbols)) if i not in inc_atom_set]

                q_inc = inc_frags[0].get("charge", 0)
                tot_step_chg = sum(r.get("charge", 0) for r in s_info["reactants"])
                q_core = tot_step_chg - q_inc

                core_atoms = Atoms(symbols=[step_symbols[i] for i in core_idxs])
                core_calc = make_xtb(charge=q_core, multiplicity=1, solvent=solvent)
                core_atoms.calc = core_calc
                e_inc_iso = frag_energies.get(_clean_smiles(inc_frags[0]["smiles"]), 0.0)

                for idx in range(start, end + 1):
                    pos = np.array(master_frames[idx].get("coords", master_frames[idx].get("positions", [])))
                    d_val = float(np.linalg.norm(pos[bi] - pos[bj]))
                    e_bound = active_e[idx] + delta_g_assoc

                    try:
                        core_atoms.set_positions(pos[core_idxs])
                        e_core = float(core_atoms.get_potential_energy())
                    except Exception:
                        e_core = active_e[idx] - e_inc_iso
                    e_decoupled = e_core + e_inc_iso

                    if d_val >= d_contact:
                        balanced_ev[idx] = e_decoupled + s_offset
                    elif d_val <= d_cov:
                        balanced_ev[idx] = e_bound + s_offset
                    else:
                        assoc_w = _smoothstep((d_contact - d_val) / max(0.1, (d_contact - d_cov)))
                        balanced_ev[idx] = ((1.0 - assoc_w) * e_decoupled + assoc_w * e_bound) + s_offset

            # Case B: Dissociative Step (leaving group departs)
            elif leaving_frags and breaking_bonds and step_symbols:
                bi, bj = breaking_bonds[0]["from"], breaking_bonds[0]["to"]
                d0_break = rcov(step_symbols[bi]) + rcov(step_symbols[bj])
                d_tight, d_free = d0_break * 1.25, d0_break * 1.65

                leaving_atom_set = set(leaving_frags[0].get("atom_indices", []))
                if not leaving_atom_set:
                    leaving_atom_set = {bj}
                hero_idxs = [i for i in range(len(step_symbols)) if i not in leaving_atom_set]

                q_leaving = leaving_frags[0].get("charge", 0)
                tot_step_chg = sum(r.get("charge", 0) for r in s_info["reactants"])
                q_hero = tot_step_chg - q_leaving

                hero_atoms = Atoms(symbols=[step_symbols[i] for i in hero_idxs])
                hero_calc = make_xtb(charge=q_hero, multiplicity=1, solvent=solvent)
                hero_atoms.calc = hero_calc
                e_leaving_iso = frag_energies.get(_clean_smiles(leaving_frags[0]["smiles"]), 0.0)

                for idx in range(start, end + 1):
                    pos = np.array(master_frames[idx].get("coords", master_frames[idx].get("positions", [])))
                    d_val = float(np.linalg.norm(pos[bi] - pos[bj]))
                    e_bound = active_e[idx] + delta_g_assoc
                    try:
                        hero_atoms.set_positions(pos[hero_idxs])
                        e_hero = float(hero_atoms.get_potential_energy())
                    except Exception:
                        e_hero = active_e[idx] - e_leaving_iso
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
            balanced_ev[b_next["start_frame"] :] += shift

            b_start, b_end = b_curr["end_frame"] + 1, b_next["start_frame"] - 1
            if b_start <= b_end:
                n_b = b_end - b_start + 1
                y0, y1 = balanced_ev[b_curr["end_frame"]], balanced_ev[b_next["start_frame"]]
                g0 = balanced_ev[b_curr["end_frame"]] - balanced_ev[max(0, b_curr["end_frame"] - 1)]
                g1 = balanced_ev[min(n_frames - 1, b_next["start_frame"] + 1)] - balanced_ev[b_next["start_frame"]]
                tau = np.linspace(0, 1, n_b + 2)[1:-1]
                balanced_ev[b_start : b_end + 1] = _hermite_spline(tau, y0, y1, g0, g1)

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
        fr["energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["rel_energy_kcal"] = round(float(rel_kcal[i]), 3)
        fr["balanced_energy_ev"] = round(float(balanced_ev[i]), 5)

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
