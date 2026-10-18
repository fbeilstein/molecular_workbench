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


class Step5Cascade(PipelineStep):
    """Step 5: Inspect and tune multi-step reaction cascades with rotamer alignment & bridge frames."""

    @property
    def step_id(self) -> str:
        return "step_5_cascade"

    @property
    def step_name(self) -> str:
        return "5. Multistep Cascade"

    def _discover_trajectories(self) -> List[Dict[str, Any]]:
        from .project_manager import init_project, _step_traj_path, STEPS_DIR
        from .ts_engine.trajectory_stitcher import sort_trajectories_chemically, get_canonical_smiles_set

        proj = init_project()
        steps = proj.get("steps", [])

        results = []
        seen_smiles = set()
        seen_paths = set()

        # 1. Primary source: Discover step trajectories in strict project index order
        for s in sorted(steps, key=lambda x: x.get("index", 1)):
            sid = s["id"]
            s_idx = s.get("index", 1)
            t_path = _step_traj_path(sid)
            traj_data = None

            if os.path.exists(t_path):
                try:
                    with open(t_path, "r") as f:
                        tj = json.load(f)
                    if tj.get("frames") and tj.get("n_steps", 1) == 1:
                        traj_data = tj
                except Exception:
                    pass

            # Self-healing if step trajectory missing: look for matching SMILES in trajectories/ or out/
            if not traj_data:
                target_smi = s.get("smiles", "")
                cands = []
                cand_indexed = os.path.join(TRAJECTORIES_DIR, f"irc_trajectory{s_idx}.json")
                if os.path.exists(cand_indexed):
                    cands.append(cand_indexed)
                if os.path.exists(TRAJECTORIES_DIR):
                    cands.extend(sorted(glob.glob(os.path.join(TRAJECTORIES_DIR, "*.json"))))
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
                        with open(cp, "r") as f:
                            tj = json.load(f)
                        if not tj.get("frames") or tj.get("n_steps", 1) > 1 or tj.get("title") == "Master Multistep IRC Trajectory":
                            continue
                        cp_smi = tj.get("reaction_smiles", "")
                        if cp_smi and target_smi and ">>" in cp_smi and ">>" in target_smi:
                            r_cp, p_cp = cp_smi.split(">>")[:2]
                            r_tgt, p_tgt = target_smi.split(">>")[:2]
                            if get_canonical_smiles_set(r_cp) == get_canonical_smiles_set(r_tgt) and get_canonical_smiles_set(p_cp) == get_canonical_smiles_set(p_tgt):
                                traj_data = tj
                                try:
                                    with open(t_path, "w") as out_f:
                                        json.dump(tj, out_f, indent=2)
                                except Exception:
                                    pass
                                break
                    except Exception:
                        continue

            if traj_data and traj_data.get("frames"):
                smi = traj_data.get("reaction_smiles") or s.get("smiles", "")
                seen_smiles.add(smi)
                seen_paths.add(t_path)
                results.append({
                    "step_id": sid,
                    "step_index": s_idx,
                    "title": s.get("title", f"Step {s_idx}"),
                    "path": t_path,
                    "filename": f"{sid}_irc_trajectory.json",
                    "reaction_smiles": smi,
                    "n_frames": len(traj_data.get("frames", [])),
                    "n_atoms": len(traj_data.get("symbols", [])),
                    "n_steps": 1,
                    "barriers": traj_data.get("barriers", {}),
                    "data": traj_data
                })

        # 2. Check trajectories directory for additional standalone trajectories
        if os.path.exists(TRAJECTORIES_DIR):
            for p in sorted(glob.glob(os.path.join(TRAJECTORIES_DIR, "*.json"))):
                if p.endswith("prolonged_pathway.json") or p.endswith("master_pathway.json") or p in seen_paths:
                    continue
                try:
                    with open(p, "r") as f:
                        tj = json.load(f)
                    if not tj.get("frames") or tj.get("n_steps", 1) > 1 or tj.get("title") == "Master Multistep IRC Trajectory":
                        continue
                    smi = tj.get("reaction_smiles") or os.path.basename(p)
                    if smi in seen_smiles:
                        continue
                    seen_smiles.add(smi)
                    seen_paths.add(p)
                    idx = len(results) + 1
                    results.append({
                        "step_id": f"file_{idx}",
                        "step_index": idx,
                        "title": f"Step {idx}: {smi}",
                        "path": p,
                        "filename": os.path.basename(p),
                        "reaction_smiles": tj.get("reaction_smiles", ""),
                        "n_frames": len(tj.get("frames", [])),
                        "n_atoms": len(tj.get("symbols", [])),
                        "n_steps": 1,
                        "barriers": tj.get("barriers", {}),
                        "data": tj
                    })
                except Exception:
                    continue

        # 3. Chemical ordering validation: ensure sorted chemically
        raw_trajs = [r["data"] for r in results]
        sorted_trajs = sort_trajectories_chemically(raw_trajs)
        if sorted_trajs and len(sorted_trajs) == len(results):
            data_to_result = {id(r["data"]): r for r in results}
            ordered_results = [data_to_result[id(t)] for t in sorted_trajs if id(t) in data_to_result]
            if len(ordered_results) == len(results):
                for idx, r in enumerate(ordered_results, start=1):
                    r["step_index"] = idx
                return ordered_results

        return results

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        available = self._discover_trajectories()
        step_state = current_state.get(self.step_id, {})
        megajson = step_state.get("megajson")

        # Validate cached prolonged_pathway.json against available sequence
        if not megajson and os.path.exists(PROLONGED_PATH):
            try:
                with open(PROLONGED_PATH, "r") as f:
                    cand_mj = json.load(f)
                if cand_mj and cand_mj.get("step_boundaries") and len(available) >= 2:
                    boundaries = cand_mj.get("step_boundaries", [])
                    b_smiles = [b.get("reaction_smiles") for b in boundaries]
                    a_smiles = [a.get("reaction_smiles") for a in available[:len(boundaries)]]
                    # Only accept cached cascade if its step boundaries match available sequence
                    if b_smiles == a_smiles and b_smiles:
                        megajson = cand_mj
            except Exception:
                pass

        if not megajson and len(available) >= 2:
            try:
                megajson = stitch_multistep_trajectories([a["data"] for a in available[:2]], n_bridge_frames=8, auto_sort=True)
                with open(PROLONGED_PATH, "w") as f:
                    json.dump(megajson, f, indent=2)
            except Exception:
                megajson = None

        stitching_report = megajson.get("stitching_alignment_report", []) if megajson else []
        step_boundaries = megajson.get("step_boundaries", []) if megajson else []

        traj_meta = [{
            "step_id": a.get("step_id", f"step_{i+1}"),
            "step_index": a.get("step_index", i + 1),
            "title": a.get("title", f"Step {i+1}"),
            "filename": a["filename"],
            "path": a["path"],
            "reaction_smiles": a["reaction_smiles"],
            "n_frames": a["n_frames"],
            "n_atoms": a["n_atoms"]
        } for i, a in enumerate(available)]

        return {
            "type": "cascade_ui",
            "trajectories": traj_meta,
            "megajson": megajson,
            "has_cascade": bool(megajson),
            "stitching_report": stitching_report,
            "step_boundaries": step_boundaries,
            "n_bridge_frames": step_state.get("n_bridge_frames", 8),
            "junction_choices": step_state.get("junction_choices", {}),
            "message": "Inspect aligned reaction sequence, select rotamer/scaffold modes, and tune bridge interpolation."
        }

    def _handle_stitch_cascade(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> Dict[str, Any]:
        user_input = user_input or {}
        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))
        ordered_files = user_input.get("ordered_files", [])

        available_list = self._discover_trajectories()
        if not available_list:
            return {"error": "No trajectories found to stitch."}

        to_stitch = []
        if ordered_files:
            # Map filenames, step_ids, or paths given by user
            for identifier in ordered_files:
                for a in available_list:
                    if identifier in (a.get("filename"), a.get("step_id"), a.get("path"), str(a.get("step_index"))):
                        if a["data"] not in to_stitch:
                            to_stitch.append(a["data"])
                        break

        if len(to_stitch) < 2:
            to_stitch = [a["data"] for a in available_list][:2]

        if len(to_stitch) < 2:
            return {"error": "Need at least 2 trajectories to stitch a cascade."}

        try:
            # If user provided explicit manual ordering, honor it (auto_sort=False), otherwise auto_sort=True
            megajson = stitch_multistep_trajectories(
                to_stitch,
                junction_choices=junction_choices,
                n_bridge_frames=n_bridge,
                auto_sort=(not ordered_files)
            )
            os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
            with open(PROLONGED_PATH, "w") as f:
                json.dump(megajson, f, indent=2)

            # Also update out/irc_trajectory.json master
            try:
                out_traj = os.path.join(OUT_DIR, "irc_trajectory.json")
                with open(out_traj, "w") as f:
                    json.dump(megajson, f, indent=2)
            except Exception:
                pass

            return {
                "success": True,
                "megajson": megajson,
                "stitching_report": megajson.get("stitching_alignment_report", []),
                "step_boundaries": megajson.get("step_boundaries", []),
                "message": f"Successfully stitched {megajson.get('n_steps', len(to_stitch))} steps ({megajson.get('total_frames', len(megajson['frames']))} frames)!"
            }
        except Exception as e:
            return {"error": f"Stitching failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action = user_input.get("action_type")
        if action == "stitch_cascade":
            return self._handle_stitch_cascade(global_state, user_input)
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
