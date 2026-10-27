"""
Project Manager Module
Manages multistep reaction cascade project lifecycle, step indexing, and master stitching.
Storage and path resolution delegated to project_store.py.
"""

import os
import shutil
import json
from typing import Any, Dict, List, Optional, Tuple

from .project_store import (
    ROOT_DIR,
    OUT_DIR,
    STEPS_DIR,
    PROJECT_FILE,
    LEGACY_STATE_FILE,
    LEGACY_TRAJ_FILE,
    ensure_dirs,
    get_step_dir,
    step_state_path,
    step_traj_path,
    read_json_file,
    write_json_file,
    sync_legacy_files,
    reindex_step_directories,
    seed_previous_product,
    init_project,
)

# Backward-compatibility aliases
_ensure_dirs = ensure_dirs
_step_state_path = step_state_path
_step_traj_path = step_traj_path
from .rxn_decomposer import decompose_multistep_reaction


def save_project_meta(proj: Dict[str, Any]):
    """Persist the project metadata."""
    write_json_file(PROJECT_FILE, proj)


def get_project_summary() -> Dict[str, Any]:
    """Return a high-level summary of all reaction steps and their completion status."""
    proj = init_project()
    active_id = proj.get("active_step_id", "step_1")
    steps_meta = []

    for s in proj.get("steps", []):
        sid = s["id"]
        st = read_json_file(step_state_path(sid), {})
        has_input = bool(st.get("step_1_input", {}).get("raw_smiles"))
        has_mapping = bool(st.get("step_2_mapping", {}).get("final_mapping"))
        has_ts = bool(st.get("step_3_ts", {}).get("ts_geometry") or st.get("step_3_ts", {}).get("ts_guess"))
        has_irc = os.path.exists(step_traj_path(sid))

        status = "Completed" if has_irc else "In Progress" if has_ts else "Configured" if has_input else "Empty"
        steps_meta.append({
            "id": sid, "index": s.get("index", 1), "title": s.get("title", f"Step {s.get('index', 1)}"),
            "smiles": s.get("smiles", st.get("step_1_input", {}).get("raw_smiles", "")),
            "status": status, "is_active": (sid == active_id),
            "completed_phases": {"input": has_input, "mapping": has_mapping, "ts": has_ts, "irc": has_irc}
        })

    return {
        "active_step_id": active_id, "solvent": proj.get("solvent", "water"),
        "total_steps": len(steps_meta), "steps": steps_meta
    }


def load_active_state() -> Dict[str, Any]:
    """Load the state of the currently active reaction step."""
    proj = init_project()
    active_id = proj.get("active_step_id", "step_1")
    s_path = step_state_path(active_id)
    t_path = step_traj_path(active_id)
    state = read_json_file(s_path, {})

    if "solvent" not in state and proj.get("solvent"):
        state["solvent"] = proj.get("solvent", "water")

    traj_data = read_json_file(t_path)
    if traj_data and traj_data.get("frames") and traj_data.get("n_steps", 1) == 1 and traj_data.get("title") != "Master Multistep IRC Trajectory":
        state.setdefault("step_4_irc", {})["trajectory"] = traj_data
        sync_legacy_files(active_id, s_path, t_path)
    else:
        if "step_4_irc" in state:
            state["step_4_irc"].pop("trajectory", None)
        sync_legacy_files(active_id, s_path, None)

    return state


def save_step_state(step_id: str, state: Dict[str, Any]):
    """Save state into a specific reaction step directory with cross-step pollution safeguards."""
    proj = init_project()
    s_path = step_state_path(step_id)
    state_to_save = dict(state)
    raw_smiles = state.get("step_1_input", {}).get("raw_smiles")

    registered_smiles = next((s.get("smiles") for s in proj.get("steps", []) if s["id"] == step_id), None)
    if registered_smiles and raw_smiles and registered_smiles != raw_smiles:
        other_matches = [s["id"] for s in proj.get("steps", []) if s.get("smiles") == raw_smiles and s["id"] != step_id]
        if other_matches:
            step_id = other_matches[0]
            s_path = step_state_path(step_id)

    traj = state_to_save.get("step_4_irc", {}).get("trajectory")
    if traj and (traj.get("n_steps", 1) > 1 or traj.get("title") == "Master Multistep IRC Trajectory"):
        s4 = dict(state_to_save.get("step_4_irc", {}))
        s4.pop("trajectory", None)
        state_to_save["step_4_irc"] = s4
        traj = None

    write_json_file(s_path, state_to_save)

    t_path = step_traj_path(step_id)
    if traj and traj.get("n_steps", 1) == 1 and traj.get("title") != "Master Multistep IRC Trajectory":
        write_json_file(t_path, traj)
        if step_id == proj.get("active_step_id"):
            sync_legacy_files(step_id, s_path, t_path)
    elif step_id == proj.get("active_step_id"):
        sync_legacy_files(step_id, s_path, None)

    solvent = state.get("solvent")
    updated = False
    if solvent and solvent != proj.get("solvent"):
        proj["solvent"] = solvent
        updated = True
    for s in proj.get("steps", []):
        if s["id"] == step_id and raw_smiles and s.get("smiles") != raw_smiles:
            s["smiles"] = raw_smiles
            s["title"] = f"Step {s.get('index', 1)}: {raw_smiles}"
            updated = True
    if updated:
        save_project_meta(proj)


def save_active_state(state: Dict[str, Any]):
    """Save state into the currently active reaction step directory."""
    proj = init_project()
    save_step_state(proj.get("active_step_id", "step_1"), state)


def select_step(target_step_id: str, current_state: Optional[Dict[str, Any]] = None) -> Tuple[Dict[str, Any], Dict[str, Any]]:
    """Switch active reaction step."""
    proj = init_project()
    if current_state:
        save_active_state(current_state)

    valid_ids = [s["id"] for s in proj.get("steps", [])]
    if target_step_id not in valid_ids:
        target_step_id = valid_ids[0] if valid_ids else "step_1"

    proj["active_step_id"] = target_step_id
    save_project_meta(proj)
    return load_active_state(), get_project_summary()


def add_step(seed_from_previous: bool = True, current_state: Optional[Dict[str, Any]] = None) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """Add a new reaction step to the project, optionally seeding product from the previous step."""
    proj = init_project()
    if current_state:
        save_active_state(current_state)

    steps = proj.get("steps", [])
    new_index = len(steps) + 1
    new_id = f"step_{new_index}"
    new_state: Dict[str, Any] = {"solvent": proj.get("solvent", "water")}

    if seed_from_previous and steps:
        prev_id = proj.get("active_step_id") or steps[-1]["id"]
        seeded = seed_previous_product(prev_id)
        new_state.update(seeded)

    write_json_file(step_state_path(new_id), new_state)
    steps.append({
        "id": new_id, "index": new_index,
        "title": f"Step {new_index}: {new_state.get('step_1_input', {}).get('raw_smiles', 'New Reaction')}",
        "smiles": new_state.get('step_1_input', {}).get('raw_smiles', '')
    })
    proj["steps"] = steps
    proj["active_step_id"] = new_id
    save_project_meta(proj)
    return new_id, new_state, get_project_summary()


def delete_step(step_id: str) -> Tuple[str, Dict[str, Any], Dict[str, Any]]:
    """Delete a step from the project, cleanly re-indexing remaining steps."""
    proj = init_project()
    steps = proj.get("steps", [])

    if len(steps) <= 1:
        s_path, t_path = step_state_path("step_1"), step_traj_path("step_1")
        if os.path.exists(s_path): os.remove(s_path)
        if os.path.exists(t_path): os.remove(t_path)
        proj["steps"] = [{"id": "step_1", "index": 1, "title": "Step 1: New Reaction", "smiles": ""}]
        proj["active_step_id"] = "step_1"
        save_project_meta(proj)
        return "step_1", {}, get_project_summary()

    remaining = [s for s in steps if s["id"] != step_id]
    new_steps, old_to_new = reindex_step_directories(step_id, remaining)

    active_id = proj.get("active_step_id", "step_1")
    new_active_id = old_to_new.get(active_id, new_steps[0]["id"])
    proj["steps"] = new_steps
    proj["active_step_id"] = new_active_id
    save_project_meta(proj)
    return new_active_id, load_active_state(), get_project_summary()


def stitch_project_trajectories() -> Tuple[bool, Optional[Dict[str, Any]], str]:
    """Stitch all completed step trajectories into one continuous master reaction profile in strict project order."""
    proj = init_project()
    steps = proj.get("steps", [])
    trajs = []

    for s in steps:
        sid = s["id"]
        tj = read_json_file(step_traj_path(sid))
        if not tj or not tj.get("frames"):
            st = read_json_file(step_state_path(sid), {})
            tj = st.get("step_4_irc", {}).get("trajectory")
            if tj and tj.get("frames") and tj.get("n_steps", 1) == 1:
                write_json_file(step_traj_path(sid), tj)
            else:
                tj = None
        if tj and tj.get("frames"):
            trajs.append(tj)

    if not trajs:
        return False, None, "No completed IRC step trajectories found to stitch. Complete Step 4 first."

    from .ts_engine.trajectory_stitcher import stitch_multistep_trajectories
    try:
        master = stitch_multistep_trajectories(trajs, n_bridge_frames=8, auto_trim=True, auto_sort=False)
        write_json_file(os.path.join(OUT_DIR, 'master_pathway.json'), master)
        write_json_file(os.path.join(ROOT_DIR, 'trajectories', 'prolonged_pathway.json'), master)
        write_json_file(LEGACY_TRAJ_FILE, master)
        n_frames = master.get("total_frames", len(master["frames"]))
        n_steps = master.get("n_steps", len(trajs))
        return True, master, f"Master pathway prepared across {n_steps} reaction step{'s' if n_steps != 1 else ''} ({n_frames} total frames)!"
    except Exception as e:
        return False, None, f"Failed to stitch master pathway: {e}"
