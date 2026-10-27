import glob
import json
import os
from typing import Any, Dict, List
from .base import PipelineStep, StepResult
from .ts_engine.trajectory_stitcher import stitch_multistep_trajectories

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, 'trajectories')
PROLONGED_PATH = os.path.join(TRAJECTORIES_DIR, 'prolonged_pathway.json')


def invert_trajectory_direction(traj: Dict[str, Any]) -> Dict[str, Any]:
    """Invert downhill IRC trajectory (Product -> TS -> Reactant becomes Reactant -> TS -> Product)."""
    inv = json.loads(json.dumps(traj))
    frames = inv.get("frames", [])
    if not frames:
        return inv
    inv_frames = frames[::-1]
    n = len(inv_frames)
    ts_idx = inv.get("ts_frame_index", 0)
    new_ts_idx = max(0, n - 1 - ts_idx)
    inv["ts_frame_index"] = new_ts_idx

    e_base = inv_frames[0].get("energy_ev", 0.0)
    for k, fr in enumerate(inv_frames):
        fr["step"] = k - new_ts_idx
        fr["is_ts"] = (k == new_ts_idx)
        if "energy_ev" in fr:
            fr["rel_energy_kcal"] = round((fr["energy_ev"] - e_base) * 23.0605, 3)

    smi = inv.get("reaction_smiles", "")
    if ">>" in smi:
        parts = smi.split(">>")
        inv["reaction_smiles"] = f"{parts[1]}>>{parts[0]}"

    inv["frames"] = inv_frames
    return inv


class Step5Cascade(PipelineStep):
    """Step 5: Inspect and tune multi-step reaction cascades with rotamer alignment & bridge frames."""

    @property
    def step_id(self) -> str:
        return "step_5_cascade"

    @property
    def step_name(self) -> str:
        return "5. Multistep Cascade"

    def _discover_trajectories(self) -> List[Dict[str, Any]]:
        """
        Discover trajectories in strict project index order [Step 1, Step 2, ...].
        No automated topological re-sorting or fuzzy multi-directory globbing.
        """
        from .project_manager import init_project, _step_traj_path, _step_state_path

        proj = init_project()
        steps = proj.get("steps", [])

        results = []
        for s in sorted(steps, key=lambda x: x.get("index", 1)):
            sid = s["id"]
            s_idx = s.get("index", 1)
            t_path = _step_traj_path(sid)
            traj_data = None

            # 1. Primary isolated source: step directory irc_trajectory.json
            if os.path.exists(t_path):
                try:
                    with open(t_path, "r") as f:
                        tj = json.load(f)
                    if tj.get("frames") and tj.get("n_steps", 1) == 1:
                        traj_data = tj
                except Exception:
                    traj_data = None

            # 2. Secondary isolated source: step directory reaction_state.json
            if not traj_data:
                s_state_path = _step_state_path(sid)
                if os.path.exists(s_state_path):
                    try:
                        with open(s_state_path, "r") as f:
                            st = json.load(f)
                        tj = st.get("step_4_irc", {}).get("trajectory")
                        if tj and tj.get("frames") and tj.get("n_steps", 1) == 1:
                            traj_data = tj
                            # Resync to step trajectory file
                            try:
                                with open(t_path, "w") as out_f:
                                    json.dump(tj, out_f, indent=2)
                            except Exception:
                                pass
                    except Exception:
                        pass

            if traj_data and traj_data.get("frames"):
                smi = traj_data.get("reaction_smiles") or s.get("smiles", "")
                results.append({
                    "step_id": sid, "step_index": s_idx, "title": s.get("title", f"Step {s_idx}"),
                    "path": t_path, "filename": f"{sid}_irc_trajectory.json", "reaction_smiles": smi,
                    "n_frames": len(traj_data.get("frames", [])), "n_atoms": len(traj_data.get("symbols", [])),
                    "n_steps": 1, "barriers": traj_data.get("barriers", {}), "data": traj_data
                })
        return results

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        available = self._discover_trajectories()
        step_state = current_state.get(self.step_id, {})
        megajson = step_state.get("megajson")

        # Validate cached prolonged_pathway.json against available sequence
        has_coords = bool(megajson and megajson.get("stitching_alignment_report") and "prev_symbols" in megajson["stitching_alignment_report"][0])
        if (not megajson or not has_coords) and os.path.exists(PROLONGED_PATH):
            try:
                with open(PROLONGED_PATH, "r") as f:
                    cand_mj = json.load(f)
                if cand_mj and cand_mj.get("step_boundaries") and len(available) >= 1:
                    boundaries = cand_mj.get("step_boundaries", [])
                    b_smiles = [b.get("reaction_smiles") for b in boundaries]
                    a_smiles = [a.get("reaction_smiles") for a in available]
                    if len(boundaries) == len(available) and b_smiles == a_smiles and b_smiles:
                        megajson = cand_mj
            except Exception:
                pass

        has_coords = bool(megajson and megajson.get("stitching_alignment_report") and "prev_symbols" in megajson["stitching_alignment_report"][0])
        if (not megajson or not has_coords) and len(available) >= 1:
            try:
                megajson = stitch_multistep_trajectories([a["data"] for a in available], n_bridge_frames=8, auto_sort=False, balance_energies=True)
                os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
                with open(PROLONGED_PATH, "w") as f:
                    json.dump(megajson, f, indent=2)
            except Exception:
                megajson = None

        stitching_report = megajson.get("stitching_alignment_report", []) if megajson else []
        step_boundaries = megajson.get("step_boundaries", []) if megajson else []
        traj_meta = [{
            "step_id": a.get("step_id", f"step_{i+1}"), "step_index": a.get("step_index", i + 1),
            "title": a.get("title", f"Step {i+1}"), "filename": a["filename"], "path": a["path"],
            "reaction_smiles": a["reaction_smiles"], "n_frames": a["n_frames"], "n_atoms": a["n_atoms"]
        } for i, a in enumerate(available)]

        return {
            "type": "cascade_ui", "trajectories": traj_meta, "megajson": megajson,
            "has_cascade": bool(megajson), "stitching_report": stitching_report,
            "step_boundaries": step_boundaries, "n_bridge_frames": step_state.get("n_bridge_frames", 8),
            "junction_choices": step_state.get("junction_choices", {}),
            "message": "Inspect aligned reaction sequence, select rotamer/scaffold modes, and tune bridge interpolation."
        }

    def _handle_stitch_cascade(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> Dict[str, Any]:
        user_input = user_input or {}
        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))
        ordered_files = user_input.get("ordered_files", [])
        reversed_steps = set(user_input.get("reversed_steps", []))

        available_list = self._discover_trajectories()
        if not available_list:
            return {"error": "No trajectories found to stitch. Complete Step 4 first."}

        to_stitch = []
        if ordered_files:
            id_map = {a.get("step_id"): a["data"] for a in available_list}
            for a in available_list:
                id_map[a.get("filename")] = a["data"]
                id_map[str(a.get("step_index"))] = a["data"]
                id_map[a.get("path")] = a["data"]

            for identifier in ordered_files:
                if identifier in id_map:
                    t_data = id_map[identifier]
                    if identifier in reversed_steps:
                        t_data = invert_trajectory_direction(t_data)
                    to_stitch.append(t_data)

        if not to_stitch:
            to_stitch = [a["data"] for a in available_list]
        if not to_stitch:
            return {"error": "No trajectories available to stitch."}

        try:
            megajson = stitch_multistep_trajectories(to_stitch, junction_choices=junction_choices, n_bridge_frames=n_bridge, auto_sort=False, balance_energies=True)
            os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
            with open(PROLONGED_PATH, "w") as f:
                json.dump(megajson, f, indent=2)

            try:
                with open(os.path.join(OUT_DIR, "irc_trajectory.json"), "w") as f:
                    json.dump(megajson, f, indent=2)
            except Exception:
                pass

            return {
                "success": True, "megajson": megajson,
                "stitching_report": megajson.get("stitching_alignment_report", []),
                "step_boundaries": megajson.get("step_boundaries", []),
                "message": f"Successfully stitched {megajson.get('n_steps', len(to_stitch))} step{'s' if megajson.get('n_steps', len(to_stitch)) != 1 else ''} ({megajson.get('total_frames', len(megajson['frames']))} frames)!"
            }
        except Exception as e:
            return {"error": f"Stitching failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        action = (user_input or {}).get("action_type")
        if action == "stitch_cascade":
            return self._handle_stitch_cascade(global_state or {}, user_input)
        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = user_input.get("megajson")
        if not megajson and os.path.exists(PROLONGED_PATH):
            try:
                with open(PROLONGED_PATH, "r") as f:
                    megajson = json.load(f)
            except Exception:
                pass
        return StepResult(success=True, data={"megajson": megajson})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("megajson"))
