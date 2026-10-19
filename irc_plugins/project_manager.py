import os
import json
import shutil
from typing import Any, Dict, List, Optional, Tuple

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
STEPS_DIR = os.path.join(OUT_DIR, 'steps')
PROJECT_FILE = os.path.join(OUT_DIR, 'project.json')
LEGACY_STATE_FILE = os.path.join(OUT_DIR, 'reaction_state.json')
LEGACY_TRAJ_FILE = os.path.join(OUT_DIR, 'irc_trajectory.json')


def _ensure_dirs():
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(STEPS_DIR, exist_ok=True)


def get_step_dir(step_id: str) -> str:
    _ensure_dirs()
    step_dir = os.path.join(STEPS_DIR, step_id)
    os.makedirs(step_dir, exist_ok=True)
    return step_dir


def _step_state_path(step_id: str) -> str:
    return os.path.join(get_step_dir(step_id), 'reaction_state.json')


def _step_traj_path(step_id: str) -> str:
    return os.path.join(get_step_dir(step_id), 'irc_trajectory.json')


def init_project() -> Dict[str, Any]:
    """Initialize or load the multistep project structure."""
    _ensure_dirs()
    proj = None
    if os.path.exists(PROJECT_FILE):
        try:
            with open(PROJECT_FILE, 'r') as f:
                proj = json.load(f)
        except Exception:
            pass

    if not proj or not proj.get("steps") or not proj.get("active_step_id"):
        # Migrate from legacy single-reaction setup if present
        step1_dir = get_step_dir("step_1")
        step1_state_path = os.path.join(step1_dir, 'reaction_state.json')
        step1_traj_path = os.path.join(step1_dir, 'irc_trajectory.json')

        if os.path.exists(LEGACY_STATE_FILE) and not os.path.exists(step1_state_path):
            try:
                shutil.copy2(LEGACY_STATE_FILE, step1_state_path)
            except Exception:
                pass

        if os.path.exists(LEGACY_TRAJ_FILE) and not os.path.exists(step1_traj_path):
            try:
                shutil.copy2(LEGACY_TRAJ_FILE, step1_traj_path)
            except Exception:
                pass

        proj = {
            "active_step_id": "step_1",
            "solvent": "water",
            "steps": [
                {
                    "id": "step_1",
                    "index": 1,
                    "title": "Step 1: Reaction",
                    "smiles": ""
                }
            ]
        }

    # Self-healing check: if project currently has only 1 step, check if step 1 has an undecomposed multi-arrow reaction
    if len(proj.get("steps", [])) <= 1:
        step1_state_path = _step_state_path("step_1")
        if os.path.exists(step1_state_path):
            try:
                with open(step1_state_path, 'r') as f:
                    st1 = json.load(f)
                raw_smi = st1.get("step_1_input", {}).get("raw_smiles", "")
                raw_rxn = st1.get("step_1_input", {}).get("raw_rxnfile", "")
                if raw_smi and (raw_rxn or ">>" in raw_smi):
                    subreactions = decompose_multistep_reaction(raw_smi, raw_rxn)
                    if len(subreactions) > 1:
                        solvent = st1.get("solvent", proj.get("solvent", "water"))
                        new_steps = []
                        for i, sub in enumerate(subreactions, start=1):
                            sid = f"step_{i}"
                            s_dir = os.path.join(STEPS_DIR, sid)
                            os.makedirs(s_dir, exist_ok=True)
                            s_state_file = _step_state_path(sid)
                            cur_st = {
                                "solvent": solvent,
                                "step_1_input": {
                                    "raw_smiles": sub["smiles"],
                                    "raw_rxnfile": sub["rxnfile"]
                                }
                            }
                            with open(s_state_file, 'w') as f:
                                json.dump(cur_st, f, indent=2)
                            new_steps.append({
                                "id": sid,
                                "index": i,
                                "title": f"Step {i}: {sub['smiles']}",
                                "smiles": sub["smiles"]
                            })
                        proj["steps"] = new_steps
                        proj["active_step_id"] = "step_1"
                        shutil.copy2(_step_state_path("step_1"), LEGACY_STATE_FILE)
            except Exception:
                pass

    save_project_meta(proj)
    return proj


def save_project_meta(proj: Dict[str, Any]):
    """Persist project.json metadata."""
    _ensure_dirs()
    with open(PROJECT_FILE, 'w') as f:
        json.dump(proj, f, indent=2)


def get_project_summary() -> Dict[str, Any]:
    """Get high-level overview of all reaction steps and their progress."""
    proj = init_project()
    steps_meta = []
    
    for s in proj.get("steps", []):
        sid = s["id"]
        s_path = _step_state_path(sid)
        t_path = _step_traj_path(sid)
        
        state_data = {}
        if os.path.exists(s_path):
            try:
                with open(s_path, 'r') as f:
                    state_data = json.load(f)
            except Exception:
                pass
                
        traj_data = {}
        if os.path.exists(t_path):
            try:
                with open(t_path, 'r') as f:
                    traj_data = json.load(f)
            except Exception:
                pass

        raw_smiles = state_data.get("step_1_input", {}).get("raw_smiles", s.get("smiles", ""))
        has_input = bool(raw_smiles and ">>" in raw_smiles)
        has_mapping = bool(state_data.get("step_2_mapping", {}).get("final_mapping"))
        has_ts = bool(state_data.get("step_3_ts", {}).get("docked_reactants") or state_data.get("step_3_ts", {}).get("ts_geometry"))
        has_irc = bool(traj_data and traj_data.get("frames"))

        barriers = traj_data.get("barriers", {})
        title = f"Step {s.get('index', 1)}: {raw_smiles}" if raw_smiles else f"Step {s.get('index', 1)}: New Reaction"

        steps_meta.append({
            "id": sid,
            "index": s.get("index", 1),
            "title": title,
            "smiles": raw_smiles,
            "has_input": has_input,
            "has_mapping": has_mapping,
            "has_ts": has_ts,
            "has_irc": has_irc,
            "barrier_fwd": barriers.get("forward_barrier_kcal"),
            "barrier_rev": barriers.get("reverse_barrier_kcal"),
            "rxn_energy": barriers.get("reaction_energy_kcal"),
            "energy_type": barriers.get("energy_type", "electronic")
        })

    return {
        "active_step_id": proj.get("active_step_id", "step_1"),
        "solvent": proj.get("solvent", "water"),
        "total_steps": len(steps_meta),
        "steps": steps_meta
    }


def load_active_state() -> Dict[str, Any]:
    """Load the state of the currently active reaction step."""
    proj = init_project()
    active_id = proj.get("active_step_id", "step_1")
    s_path = _step_state_path(active_id)
    state = {}
    if os.path.exists(s_path):
        try:
            with open(s_path, 'r') as f:
                state = json.load(f)
        except Exception:
            pass

    # Ensure solvent is preserved
    if "solvent" not in state and proj.get("solvent"):
        state["solvent"] = proj.get("solvent", "water")

    # Load active trajectory if present and valid elementary IRC
    t_path = _step_traj_path(active_id)
    if os.path.exists(t_path):
        try:
            with open(t_path, 'r') as f:
                traj_data = json.load(f)
            if traj_data and traj_data.get("frames") and traj_data.get("n_steps", 1) == 1 and traj_data.get("title") != "Master Multistep IRC Trajectory":
                state.setdefault("step_4_irc", {})["trajectory"] = traj_data
            shutil.copy2(t_path, LEGACY_TRAJ_FILE)
        except Exception:
            pass
    else:
        if "step_4_irc" in state:
            state["step_4_irc"].pop("trajectory", None)
        if os.path.exists(LEGACY_TRAJ_FILE):
            try:
                os.remove(LEGACY_TRAJ_FILE)
            except Exception:
                pass

    # Sanity check: Ensure step_4_irc never contains a multistep master trajectory in an elementary step!
    cur_traj = state.get("step_4_irc", {}).get("trajectory")
    if cur_traj and (cur_traj.get("n_steps", 1) > 1 or cur_traj.get("title") == "Master Multistep IRC Trajectory"):
        if os.path.exists(t_path):
            try:
                with open(t_path, 'r') as f:
                    tj = json.load(f)
                if tj and tj.get("n_steps", 1) == 1 and tj.get("title") != "Master Multistep IRC Trajectory":
                    state["step_4_irc"]["trajectory"] = tj
                else:
                    state.pop("step_4_irc", None)
            except Exception:
                state.pop("step_4_irc", None)
        else:
            state.pop("step_4_irc", None)

    # Sync to legacy out/reaction_state.json for downstream tools
    try:
        shutil.copy2(s_path, LEGACY_STATE_FILE)
    except Exception:
        pass

    return state


def save_step_state(step_id: str, state: Dict[str, Any]):
    """Save state into a specific reaction step directory with cross-step pollution safeguards."""
    proj = init_project()
    s_path = _step_state_path(step_id)

    # Cross-step pollution safeguard:
    # If this step is registered in project.json with a specific smiles, and state has raw_smiles
    # belonging to a DIFFERENT step, redirect to the matching step instead of overwriting this step.
    registered_smiles = None
    for s in proj.get("steps", []):
        if s["id"] == step_id:
            registered_smiles = s.get("smiles")
            break

    state_to_save = dict(state)
    raw_smiles = state.get("step_1_input", {}).get("raw_smiles")
    if registered_smiles and raw_smiles and registered_smiles != raw_smiles:
        other_matches = [s["id"] for s in proj.get("steps", []) if s.get("smiles") == raw_smiles and s["id"] != step_id]
        if other_matches:
            step_id = other_matches[0]
            s_path = _step_state_path(step_id)

    # Strip master cascade trajectory from elementary step state
    traj = state_to_save.get("step_4_irc", {}).get("trajectory")
    if traj and (traj.get("n_steps", 1) > 1 or traj.get("title") == "Master Multistep IRC Trajectory"):
        s4 = dict(state_to_save.get("step_4_irc", {}))
        t_path = _step_traj_path(step_id)
        if os.path.exists(t_path):
            try:
                with open(t_path, 'r') as f:
                    tj = json.load(f)
                if tj.get("n_steps", 1) == 1 and tj.get("title") != "Master Multistep IRC Trajectory":
                    s4["trajectory"] = tj
                else:
                    s4.pop("trajectory", None)
            except Exception:
                s4.pop("trajectory", None)
        else:
            s4.pop("trajectory", None)
        state_to_save["step_4_irc"] = s4
        traj = s4.get("trajectory")

    with open(s_path, 'w') as f:
        json.dump(state_to_save, f, indent=2)

    # Sync elementary trajectory if valid
    if traj and traj.get("n_steps", 1) == 1 and traj.get("title") != "Master Multistep IRC Trajectory":
        t_path = _step_traj_path(step_id)
        with open(t_path, 'w') as f:
            json.dump(traj, f, indent=2)
        if step_id == proj.get("active_step_id"):
            try:
                shutil.copy2(t_path, LEGACY_TRAJ_FILE)
            except Exception:
                pass

    if step_id == proj.get("active_step_id"):
        try:
            shutil.copy2(s_path, LEGACY_STATE_FILE)
        except Exception:
            pass

    # Update SMILES and solvent in project.json if legitimately changed for this step
    solvent = state.get("solvent")
    updated = False
    if solvent and solvent != proj.get("solvent"):
        proj["solvent"] = solvent
        updated = True
    for s in proj.get("steps", []):
        if s["id"] == step_id and raw_smiles:
            if s.get("smiles") != raw_smiles:
                s["smiles"] = raw_smiles
                s["title"] = f"Step {s.get('index', 1)}: {raw_smiles}"
                updated = True
    if updated:
        save_project_meta(proj)


def save_active_state(state: Dict[str, Any]):
    """Save state into the currently active reaction step directory."""
    proj = init_project()
    active_id = proj.get("active_step_id", "step_1")
    save_step_state(active_id, state)


def select_step(target_step_id: str, current_state: Optional[Dict[str, Any]] = None) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Switch active reaction step, persisting current state and loading target state."""
    proj = init_project()
    active_id = proj.get("active_step_id", "step_1")
    if current_state:
        save_active_state(current_state)

    valid_ids = [s["id"] for s in proj.get("steps", [])]
    if target_step_id not in valid_ids:
        target_step_id = valid_ids[0] if valid_ids else "step_1"

    proj["active_step_id"] = target_step_id
    save_project_meta(proj)

    loaded_state = load_active_state()
    summary = get_project_summary()
    return loaded_state, summary


def add_step(seed_from_previous: bool = True, current_state: Optional[Dict[str, Any]] = None) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """Add a new reaction step to the project, optionally seeding product from the previous step."""
    proj = init_project()
    if current_state:
        save_active_state(current_state)

    steps = proj.get("steps", [])
    new_index = len(steps) + 1
    new_id = f"step_{new_index}"
    step_dir = get_step_dir(new_id)

    new_state: Dict[str, Any] = {"solvent": proj.get("solvent", "water")}

    # Seed product of previous step as reactant of this step
    if seed_from_previous and steps:
        prev_id = proj.get("active_step_id") or steps[-1]["id"]
        prev_state_path = _step_state_path(prev_id)
        prev_traj_path = _step_traj_path(prev_id)
        
        prod_smiles = ""
        seeded_blocks = []

        if os.path.exists(prev_traj_path):
            try:
                with open(prev_traj_path, 'r') as f:
                    traj = json.load(f)
                prods = traj.get("fragments", {}).get("products", [])
                if prods:
                    from .ts_engine.trajectory_stitcher import extract_fragment_molblock
                    p_smiles_list = [p.get("smiles") for p in prods if p.get("smiles")]
                    prod_smiles = ".".join(p_smiles_list)
                    seeded_blocks = [extract_fragment_molblock(traj, p) for p in prods]
            except Exception:
                pass

        if not prod_smiles and os.path.exists(prev_state_path):
            try:
                with open(prev_state_path, 'r') as f:
                    pst = json.load(f)
                raw = pst.get("step_1_input", {}).get("raw_smiles", "")
                if ">>" in raw:
                    prod_smiles = raw.split(">>")[1].strip()
            except Exception:
                pass

        if prod_smiles:
            new_state["step_1_input"] = {
                "raw_smiles": f"{prod_smiles}>>",
                "raw_rxnfile": ""
            }
        if seeded_blocks:
            new_state["step_3_ts"] = {
                "docked_reactants": seeded_blocks
            }

    # Save new state and register in project.json
    with open(_step_state_path(new_id), 'w') as f:
        json.dump(new_state, f, indent=2)

    steps.append({
        "id": new_id,
        "index": new_index,
        "title": f"Step {new_index}: {new_state.get('step_1_input', {}).get('raw_smiles', 'New Reaction')}",
        "smiles": new_state.get('step_1_input', {}).get('raw_smiles', '')
    })
    proj["steps"] = steps
    proj["active_step_id"] = new_id
    save_project_meta(proj)

    summary = get_project_summary()
    return new_id, new_state, summary


def delete_step(step_id: str) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """Delete a step from the project, cleanly re-indexing remaining steps."""
    proj = init_project()
    steps = proj.get("steps", [])
    
    # If only 1 step remains, clear its state instead of leaving 0 steps
    if len(steps) <= 1:
        s_path = _step_state_path("step_1")
        t_path = _step_traj_path("step_1")
        if os.path.exists(s_path):
            os.remove(s_path)
        if os.path.exists(t_path):
            os.remove(t_path)
        proj["steps"] = [{
            "id": "step_1",
            "index": 1,
            "title": "Step 1: New Reaction",
            "smiles": ""
        }]
        proj["active_step_id"] = "step_1"
        save_project_meta(proj)
        return "step_1", {}, get_project_summary()

    target_dir = os.path.join(STEPS_DIR, step_id)
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir, ignore_errors=True)

    remaining = [s for s in steps if s["id"] != step_id]

    # Re-index remaining steps sequentially
    new_steps = []
    old_to_new = {}
    for idx, s in enumerate(remaining, start=1):
        nid = f"step_{idx}"
        old_to_new[s["id"]] = nid
        old_dir = os.path.join(STEPS_DIR, s["id"])
        new_dir = os.path.join(STEPS_DIR, nid)
        if old_dir != new_dir and os.path.exists(old_dir):
            if os.path.exists(new_dir):
                shutil.rmtree(new_dir, ignore_errors=True)
            os.rename(old_dir, new_dir)
        new_steps.append({
            "id": nid,
            "index": idx,
            "title": s.get("title", f"Step {idx}"),
            "smiles": s.get("smiles", "")
        })

    active_id = proj.get("active_step_id", "step_1")
    new_active_id = old_to_new.get(active_id, new_steps[0]["id"])

    proj["steps"] = new_steps
    proj["active_step_id"] = new_active_id
    save_project_meta(proj)

    loaded_state = load_active_state()
    summary = get_project_summary()
    return new_active_id, loaded_state, summary


def decompose_multistep_reaction(smiles: str, rxnfile: str = "") -> List[Dict[str, str]]:
    """
    Decompose a multi-arrow reaction into sequential elementary reaction steps.
    Supports:
    1. Spatial decomposition of MDL RXN files from Ketcher (where multiple arrows
       are exported into a single RXN block with multiple products or stages).
    2. Multi-arrow SMILES (e.g., "A.B >> C >> D.E").
    3. Multi-line reaction SMILES.
    """
    # 1. Attempt spatial decomposition from rxnfile
    if rxnfile and "$MOL" in rxnfile:
        try:
            from rdkit import Chem
            raw_blocks = rxnfile.split("$MOL\n")[1:]
            mols_info = []
            for b in raw_blocks:
                clean_b = b
                if "M  END" in clean_b:
                    clean_b = clean_b[:clean_b.index("M  END") + 6]
                m = Chem.MolFromMolBlock(clean_b)
                if not m or m.GetNumAtoms() == 0:
                    continue
                smi = Chem.MolToSmiles(m)
                conf = m.GetConformer()
                xs = [conf.GetAtomPosition(a).x for a in range(m.GetNumAtoms())]
                ys = [conf.GetAtomPosition(a).y for a in range(m.GetNumAtoms())]
                mols_info.append({
                    "mol": m,
                    "smiles": smi,
                    "min_x": min(xs),
                    "max_x": max(xs),
                    "min_y": min(ys),
                    "max_y": max(ys),
                    "block": clean_b
                })

            if len(mols_info) >= 2:
                span_x = max(m["max_x"] for m in mols_info) - min(m["min_x"] for m in mols_info)
                span_y = max(m["max_y"] for m in mols_info) - min(m["min_y"] for m in mols_info)

                if span_y > span_x * 1.5:
                    # Vertical flow: top to bottom (decreasing y)
                    sorted_mols = sorted(mols_info, key=lambda m: -m["max_y"])
                    get_start = lambda m: -m["max_y"]
                    get_end = lambda m: -m["min_y"]
                else:
                    # Horizontal flow: left to right (increasing x)
                    sorted_mols = sorted(mols_info, key=lambda m: m["min_x"])
                    get_start = lambda m: m["min_x"]
                    get_end = lambda m: m["max_x"]

                # Group molecules into stages separated by reaction arrows
                stages = []
                current_stage = [sorted_mols[0]]
                current_stage_end = get_end(sorted_mols[0])
                gap_threshold = 2.8

                for m in sorted_mols[1:]:
                    start = get_start(m)
                    if start - current_stage_end > gap_threshold:
                        stages.append(current_stage)
                        current_stage = [m]
                        current_stage_end = get_end(m)
                    else:
                        current_stage.append(m)
                        current_stage_end = max(current_stage_end, get_end(m))
                if current_stage:
                    stages.append(current_stage)

                if len(stages) > 2:
                    steps = []
                    for i in range(len(stages) - 1):
                        r_mols = stages[i]
                        p_mols = stages[i + 1]

                        # Build valid MDL $RXN file for this subreaction
                        rxn_lines = ["$RXN\n\nGenerated by Molecular Workbench\n\n"]
                        rxn_lines.append(f"{len(r_mols):>3}{len(p_mols):>3}\n")
                        for rm in r_mols:
                            rxn_lines.append("$MOL\n" + Chem.MolToMolBlock(rm["mol"]))
                        for pm in p_mols:
                            rxn_lines.append("$MOL\n" + Chem.MolToMolBlock(pm["mol"]))
                        step_rxn = "".join(rxn_lines)

                        r_smis = ".".join(rm["smiles"] for rm in r_mols)
                        p_smis = ".".join(pm["smiles"] for pm in p_mols)
                        step_smi = f"{r_smis}>>{p_smis}"

                        steps.append({
                            "smiles": step_smi,
                            "rxnfile": step_rxn
                        })
                    return steps
        except Exception:
            pass

    # Fallback 1: Multi-arrow in single line: A + B >> C >> D
    clean_smi = smiles.strip() if smiles else ""
    if clean_smi.count(">>") > 1:
        parts = [p.strip() for p in clean_smi.split(">>") if p.strip()]
        steps = []
        for i in range(len(parts) - 1):
            steps.append({
                "smiles": f"{parts[i]}>>{parts[i+1]}",
                "rxnfile": rxnfile if i == 0 else ""
            })
        return steps

    # Fallback 2: Multi-line reactions: A >> B \n C >> D
    lines = [l.strip() for l in clean_smi.splitlines() if ">>" in l]
    if len(lines) > 1:
        return [{"smiles": l, "rxnfile": rxnfile if i == 0 else ""} for i, l in enumerate(lines)]

    return [{"smiles": smiles, "rxnfile": rxnfile}]


def parse_multistep_smiles(raw_text: str, rxnfile: str = "") -> List[str]:
    """Parse text or multi-arrow SMILES into elementary reaction steps."""
    subreactions = decompose_multistep_reaction(raw_text, rxnfile)
    return [s["smiles"] for s in subreactions]


def stitch_project_trajectories() -> Tuple[bool, Optional[Dict[str, Any]], str]:
    """Stitch all completed step trajectories into one continuous master reaction profile."""
    import glob
    proj = init_project()
    trajs = []
    
    for s in sorted(proj.get("steps", []), key=lambda x: x.get("index", 1)):
        sid = s["id"]
        s_idx = s.get("index", 1)
        t_path = _step_traj_path(sid)
        tj = None

        if os.path.exists(t_path):
            try:
                with open(t_path, 'r') as f:
                    tj = json.load(f)
                if not (tj.get("frames") and tj.get("n_steps", 1) == 1):
                    tj = None
            except Exception:
                tj = None

        # Self-healing: if step trajectory file is missing or invalid, search candidate files
        if not tj or not tj.get("frames"):
            target_smi = s.get("smiles", "")
            cands = []
            cand_indexed = os.path.join(ROOT_DIR, "trajectories", f"irc_trajectory{s_idx}.json")
            if os.path.exists(cand_indexed):
                cands.append(cand_indexed)
            t_dir = os.path.join(ROOT_DIR, "trajectories")
            if os.path.exists(t_dir):
                cands.extend(sorted(glob.glob(os.path.join(t_dir, "*.json"))))
            cands.extend([
                os.path.join(OUT_DIR, "irc_trajectory.json"),
                os.path.join(ROOT_DIR, "irc_trajectory.json")
            ])

            for cp in cands:
                if cp.endswith("prolonged_pathway.json") or cp.endswith("master_pathway.json"):
                    continue
                if not os.path.exists(cp):
                    continue
                try:
                    with open(cp, 'r') as f:
                        cand_d = json.load(f)
                    if cand_d.get("frames") and cand_d.get("n_steps", 1) == 1 and cand_d.get("title") != "Master Multistep IRC Trajectory":
                        cp_smi = cand_d.get("reaction_smiles", "")
                        if cp_smi and target_smi and ">>" in cp_smi and ">>" in target_smi:
                            from .ts_engine.trajectory_stitcher import get_canonical_smiles_set
                            r_cp, p_cp = cp_smi.split(">>")[:2]
                            r_tgt, p_tgt = target_smi.split(">>")[:2]
                            if get_canonical_smiles_set(r_cp) == get_canonical_smiles_set(r_tgt) and get_canonical_smiles_set(p_cp) == get_canonical_smiles_set(p_tgt):
                                tj = cand_d
                                try:
                                    with open(t_path, 'w') as out_f:
                                        json.dump(cand_d, out_f, indent=2)
                                except Exception:
                                    pass
                                break
                except Exception:
                    continue

        if tj and tj.get("frames"):
            trajs.append(tj)

    if len(trajs) < 2:
        return False, None, "At least 2 completed IRC step trajectories are required to stitch a multistep pathway."

    from .ts_engine.trajectory_stitcher import stitch_multistep_trajectories
    try:
        master = stitch_multistep_trajectories(trajs, n_bridge_frames=8, auto_trim=True, auto_sort=True)
        master_path = os.path.join(OUT_DIR, 'master_pathway.json')
        with open(master_path, 'w') as f:
            json.dump(master, f, indent=2)
        # Also sync to prolonged_pathway.json in trajectories
        prolong_path = os.path.join(ROOT_DIR, 'trajectories', 'prolonged_pathway.json')
        try:
            with open(prolong_path, 'w') as f:
                json.dump(master, f, indent=2)
        except Exception:
            pass
        # Also copy to legacy trajectory file so viewer can render it
        with open(LEGACY_TRAJ_FILE, 'w') as f:
            json.dump(master, f, indent=2)
        n_frames = master.get("total_frames", len(master["frames"]))
        n_steps = master.get("n_steps", len(trajs))
        return True, master, f"Master pathway stitched across {n_steps} reaction steps ({n_frames} total frames)!"
    except Exception as e:
        return False, None, f"Failed to stitch master pathway: {e}"
