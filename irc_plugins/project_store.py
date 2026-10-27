"""
Project Store Module
Handles filesystem paths, step directory layout, JSON serialization, and legacy file syncing.
"""

import json
import os
import shutil
from typing import Any, Dict, Optional

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
STEPS_DIR = os.path.join(OUT_DIR, 'steps')
PROJECT_FILE = os.path.join(OUT_DIR, 'project.json')
LEGACY_STATE_FILE = os.path.join(OUT_DIR, 'reaction_state.json')
LEGACY_TRAJ_FILE = os.path.join(OUT_DIR, 'irc_trajectory.json')


def ensure_dirs():
    """Ensure base project and step output directories exist."""
    os.makedirs(OUT_DIR, exist_ok=True)
    os.makedirs(STEPS_DIR, exist_ok=True)


def get_step_dir(step_id: str) -> str:
    """Get absolute path to a specific step directory, creating it if needed."""
    ensure_dirs()
    step_dir = os.path.join(STEPS_DIR, step_id)
    os.makedirs(step_dir, exist_ok=True)
    return step_dir


def step_state_path(step_id: str) -> str:
    """Path to a step's reaction_state.json."""
    return os.path.join(get_step_dir(step_id), 'reaction_state.json')


def step_traj_path(step_id: str) -> str:
    """Path to a step's irc_trajectory.json."""
    return os.path.join(get_step_dir(step_id), 'irc_trajectory.json')


def init_project() -> Dict[str, Any]:
    """Initialize or load the multistep project structure."""
    ensure_dirs()
    proj = read_json_file(PROJECT_FILE)
    if not proj or not proj.get("steps") or not proj.get("active_step_id"):
        step1_dir = get_step_dir("step_1")
        s1_state = os.path.join(step1_dir, 'reaction_state.json')
        s1_traj = os.path.join(step1_dir, 'irc_trajectory.json')
        if os.path.exists(LEGACY_STATE_FILE) and not os.path.exists(s1_state):
            try: shutil.copy2(LEGACY_STATE_FILE, s1_state)
            except Exception: pass
        if os.path.exists(LEGACY_TRAJ_FILE) and not os.path.exists(s1_traj):
            try: shutil.copy2(LEGACY_TRAJ_FILE, s1_traj)
            except Exception: pass
        proj = {
            "active_step_id": "step_1",
            "solvent": "water",
            "steps": [{"id": "step_1", "index": 1, "title": "Step 1: Primary Reaction", "smiles": ""}]
        }
        write_json_file(PROJECT_FILE, proj)
    return proj


def read_json_file(filepath: str, default: Any = None) -> Any:
    """Safely read and parse a JSON file."""
    if os.path.exists(filepath):
        try:
            with open(filepath, 'r') as f:
                return json.load(f)
        except Exception:
            pass
    return default


def write_json_file(filepath: str, data: Any):
    """Safely write JSON data with indentation."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, 'w') as f:
        json.dump(data, f, indent=2)


def sync_legacy_files(step_id: str, state_path: str, traj_path: Optional[str] = None):
    """Sync active step state and trajectory to legacy out/ files for downstream tool compatibility."""
    if os.path.exists(state_path):
        try:
            shutil.copy2(state_path, LEGACY_STATE_FILE)
        except Exception:
            pass

    if traj_path and os.path.exists(traj_path):
        try:
            shutil.copy2(traj_path, LEGACY_TRAJ_FILE)
        except Exception:
            pass
    elif not traj_path and os.path.exists(LEGACY_TRAJ_FILE):
        try:
            os.remove(LEGACY_TRAJ_FILE)
        except Exception:
            pass


def reindex_step_directories(step_id_to_delete: str, remaining_steps: list) -> tuple:
    """Delete a step directory and renumber the remaining step directories sequentially."""
    target_dir = os.path.join(STEPS_DIR, step_id_to_delete)
    if os.path.exists(target_dir):
        shutil.rmtree(target_dir, ignore_errors=True)

    new_steps, old_to_new = [], {}
    for idx, s in enumerate(remaining_steps, start=1):
        nid = f"step_{idx}"
        old_to_new[s["id"]] = nid
        old_dir, new_dir = os.path.join(STEPS_DIR, s["id"]), os.path.join(STEPS_DIR, nid)
        if old_dir != new_dir and os.path.exists(old_dir):
            if os.path.exists(new_dir):
                shutil.rmtree(new_dir, ignore_errors=True)
            os.rename(old_dir, new_dir)
        new_steps.append({"id": nid, "index": idx, "title": s.get("title", f"Step {idx}"), "smiles": s.get("smiles", "")})
    return new_steps, old_to_new


def seed_previous_product(prev_id: str) -> dict:
    """Extract seeded product SMILES and 3D MolBlock from the preceding reaction step."""
    traj = read_json_file(step_traj_path(prev_id))
    prods = traj.get("fragments", {}).get("products", []) if traj else []
    prod_smiles = ""
    seeded_blocks = []
    if prods:
        from .ts_engine.alignment import extract_fragment_molblock
        prod_smiles = ".".join(p.get("smiles") for p in prods if p.get("smiles"))
        seeded_blocks = [extract_fragment_molblock(traj, p) for p in prods]

    if not prod_smiles:
        pst = read_json_file(step_state_path(prev_id), {})
        raw = pst.get("step_1_input", {}).get("raw_smiles", "")
        if ">>" in raw:
            prod_smiles = raw.split(">>")[1].strip()

    res = {}
    if prod_smiles:
        res["step_1_input"] = {"raw_smiles": f"{prod_smiles}>>", "raw_rxnfile": ""}
    if seeded_blocks:
        res["step_3_ts"] = {"docked_reactants": seeded_blocks}
    return res
