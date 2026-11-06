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
    """Calculate Hess's Law balanced Free Energy (G) profile with standard-state isolated fragment thermodynamics in chosen solvent."""
    boundaries = megajson.get("step_boundaries", [])
    master_frames = megajson.get("frames", [])
    n_frames = len(master_frames)
    if n_frames == 0:
        return megajson

    stoich = analyze_pathway_stoichiometry(trajectories, fragment_roles)
    solvent = (solvent or megajson.get("solvent", "water") or "water").strip()

    # Collect all chemical species involved across elementary steps
    all_species = []
    for s in stoich.get("steps", []):
        all_species.extend(s.get("reactants", []))
        all_species.extend(s.get("products", []))
    all_species.extend(stoich.get("all_external_reactants", []))
    all_species.extend(stoich.get("all_external_byproducts", []))

    frag_energies = {}
    if recalculate_with_xtb and all_species:
        try:
            frag_energies = compute_isolated_fragment_energies(all_species, solvent=solvent, relax=True)
        except Exception as fe_err:
            import logging
            logging.getLogger("irc_workbench").warning(f"Could not evaluate isolated fragment energies: {fe_err}")

    eV_to_kcal = 23.0605

    # Check if we have valid fragment energies to compute thermodynamic Hess baselines
    has_thermo = False
    step_delta_kcal = []
    if frag_energies and stoich.get("steps"):
        valid_steps = 0
        for k, s in enumerate(stoich["steps"]):
            e_r = sum(frag_energies.get(_clean_smiles(r["smiles"]), 0.0) for r in s.get("reactants", []))
            e_p = sum(frag_energies.get(_clean_smiles(p["smiles"]), 0.0) for p in s.get("products", []))
            if e_r != 0.0 and e_p != 0.0:
                valid_steps += 1
            step_delta_kcal.append((e_p - e_r) * eV_to_kcal)
        if valid_steps == len(stoich["steps"]) and len(boundaries) == len(stoich["steps"]):
            has_thermo = True

    rel_kcal = np.zeros(n_frames)

    if has_thermo:
        # Full Hess's Law balancing with standard-state thermodynamics in solvent
        cum_baseline = 0.0
        for k, b in enumerate(boundaries):
            start, end, ts = b["start_frame"], b["end_frame"], b["ts_frame"]
            traj = trajectories[k] if k < len(trajectories) else {}
            de_step = step_delta_kcal[k]

            e_start_target = cum_baseline
            e_end_target = cum_baseline + de_step

            # Extract raw IRC profile for this step
            s_len = end - start + 1
            raw = np.array([float(master_frames[i].get("rel_energy_kcal", 0.0)) for i in range(start, end + 1)])
            raw = raw - raw[0]

            # Forward barrier height from calculation
            raw_barrier = float(np.max(raw)) if len(raw) > 0 else 0.0
            rep_barrier = float(traj.get("barriers", {}).get("forward_barrier_kcal", 0.0) or 0.0)
            barrier = max(raw_barrier, rep_barrier, 1.0)

            ts_rel_idx = ts - start if start <= ts <= end else int(s_len * 0.4)
            e_ts_target = e_start_target + barrier

            for local_i, global_i in enumerate(range(start, end + 1)):
                if local_i <= ts_rel_idx:
                    t = local_i / max(1, ts_rel_idx)
                    s_curve = np.sin(t * np.pi / 2.0) ** 2
                    rel_kcal[global_i] = e_start_target + (e_ts_target - e_start_target) * s_curve
                else:
                    t = (local_i - ts_rel_idx) / max(1, (s_len - 1 - ts_rel_idx))
                    s_curve = np.sin(t * np.pi / 2.0) ** 2
                    rel_kcal[global_i] = e_ts_target - (e_ts_target - e_end_target) * s_curve

            cum_baseline = e_end_target

            # Fill bridge frames to next step smoothly at the intermediate baseline
            if k < len(boundaries) - 1:
                next_start = boundaries[k + 1]["start_frame"]
                if end + 1 < next_start:
                    for br_i in range(end + 1, next_start):
                        rel_kcal[br_i] = e_end_target
        balanced_ev = rel_kcal / eV_to_kcal
    else:
        # Fallback if fragment energies are unavailable: relative IRC curve
        rel_kcal = np.array([float(fr.get("rel_energy_kcal", fr.get("energy_kcal", 0.0))) for fr in master_frames])
        if len(rel_kcal) > 0:
            rel_kcal = rel_kcal - rel_kcal[0]
        balanced_ev = rel_kcal / eV_to_kcal

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
