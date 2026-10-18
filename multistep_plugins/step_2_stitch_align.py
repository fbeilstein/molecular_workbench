from typing import Any, Dict, List
from .base import PipelineStep, StepResult
from irc_plugins.ts_engine.trajectory_stitcher import stitch_multistep_trajectories


def _ensure_junction_previews(mega: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Synthesize junction preview dictionaries from step_boundaries if missing."""
    if not mega:
        return []
    previews = mega.get("junction_previews", [])
    if previews:
        return previews
    boundaries = mega.get("step_boundaries", [])
    frames = mega.get("frames", [])
    symbols = mega.get("symbols", [])
    if len(boundaries) < 2 or not frames or not symbols:
        return []

    synthesized = []
    for i in range(len(boundaries) - 1):
        b_prev = boundaries[i]
        b_next = boundaries[i + 1]
        prev_idx = min(b_prev.get("end_frame", 0), len(frames) - 1)
        next_idx = min(b_next.get("start_frame", 0), len(frames) - 1)
        prev_coords = frames[prev_idx].get("coords") or frames[prev_idx].get("positions", [])
        next_coords = frames[next_idx].get("coords") or frames[next_idx].get("positions", [])

        bridge_frames = []
        if next_idx > prev_idx + 1:
            for bf in frames[prev_idx + 1 : next_idx]:
                b_c = bf.get("coords") or bf.get("positions", [])
                bridge_frames.append({
                    "coords": b_c,
                    "positions": b_c,
                    "energy_kcal": bf.get("energy_kcal", 0.0),
                    "is_bridge": True
                })

        synthesized.append({
            "junction_index": i,
            "step_prev_index": b_prev.get("step_index", i + 1),
            "step_next_index": b_next.get("step_index", i + 2),
            "step_prev_name": b_prev.get("step_name", f"Step {i + 1}"),
            "step_next_name": b_next.get("step_name", f"Step {i + 2}"),
            "hero_formula": "Intermediate Scaffold",
            "hero_smiles": "",
            "n_matched_atoms": len(symbols),
            "joint_rmsd_angstrom": 0.0,
            "rmsd": 0.0,
            "is_reflected": False,
            "prev_symbols": symbols,
            "prev_coords": prev_coords,
            "next_symbols": symbols,
            "next_coords_unaligned": next_coords,
            "next_coords_aligned": next_coords,
            "bridge_frames": bridge_frames,
            "n_bridge_frames": len(bridge_frames),
            "matched_prev_indices": list(range(len(symbols))),
            "matched_next_indices": list(range(len(symbols)))
        })
    return synthesized


def _extract_elementary_trajectories_from_megajson(mega: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Extract individual elementary reaction steps from a pre-stitched cascade megajson."""
    boundaries = mega.get("step_boundaries", [])
    frames = mega.get("frames", [])
    symbols = mega.get("symbols", [])
    if not boundaries or not frames:
        return []

    elementary = []
    for idx, b in enumerate(boundaries):
        start = b.get("start_frame", 0)
        end = b.get("end_frame", len(frames) - 1)
        sub_frames = [f for f in frames[start : end + 1] if not f.get("is_bridge")]
        if not sub_frames:
            continue
        ts_frame = b.get("ts_frame", 0) - start
        elem = {
            "name": b.get("step_name", f"Step_{idx+1}"),
            "reaction_smiles": b.get("reaction_smiles") or b.get("step_name", f"Step_{idx+1}"),
            "symbols": symbols,
            "bonds": mega.get("bonds", []),
            "n_atoms": len(symbols),
            "frames": sub_frames,
            "ts_frame_index": max(0, min(len(sub_frames) - 1, ts_frame)),
            "barriers": {
                "forward_barrier_kcal": b.get("barrier_kcal", 0.0),
                "reaction_energy_kcal": b.get("reaction_energy_kcal", 0.0)
            }
        }
        elementary.append(elem)
    return elementary


class Step2StitchAlign(PipelineStep):
    """Step 2: Detect Hero intermediate overlap, align coordinates with Kabsch fit, and stitch."""

    @property
    def step_id(self) -> str:
        return "step_2_stitch_align"

    @property
    def step_name(self) -> str:
        return "2. Align & Stitch"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        step_1_data = current_state.get("step_1_choose_jsons", {})
        trajectories = step_1_data.get("trajectories", [])
        existing_megajson = current_state.get(self.step_id, {}).get("megajson") or step_1_data.get("megajson")

        report = []
        junction_previews = []
        if existing_megajson:
            report = existing_megajson.get("stitching_alignment_report", [])
            junction_previews = _ensure_junction_previews(existing_megajson)
        elif len(trajectories) >= 2:
            try:
                existing_megajson = stitch_multistep_trajectories(trajectories, n_bridge_frames=8)
                report = existing_megajson.get("stitching_alignment_report", [])
                junction_previews = _ensure_junction_previews(existing_megajson)
            except Exception:
                existing_megajson = None

        saved_n_bridge = current_state.get(self.step_id, {}).get("n_bridge_frames")
        if saved_n_bridge is None and junction_previews:
            saved_n_bridge = junction_previews[0].get("n_bridge_frames", 8)
        if saved_n_bridge is None:
            saved_n_bridge = 8

        saved_trims = current_state.get(self.step_id, {}).get("junction_trims", {})

        return {
            "type": "stitch_align_ui",
            "n_steps": len(trajectories) or (existing_megajson.get("n_steps", 1) if existing_megajson else 0),
            "sequence": step_1_data.get("sequence", []),
            "stitching_report": report,
            "junction_previews": junction_previews,
            "megajson": existing_megajson,
            "has_megajson": bool(existing_megajson),
            "n_bridge_frames": saved_n_bridge,
            "junction_trims": saved_trims,
            "message": "Inspect boundary 3D alignment, choose Hero intermediate molecules, and verify zero RMSD continuity."
        }

    def _execute_stitch(
        self, global_state: Dict[str, Any], junction_choices: Dict[str, Any], n_bridge: int,
        junction_trims: Dict[str, Any] = None, auto_trim: bool = False
    ) -> Dict[str, Any]:
        step_1_data = global_state.get("step_1_choose_jsons", {})
        trajectories = step_1_data.get("trajectories", [])

        source_trajs = trajectories if len(trajectories) >= 2 else []
        existing = step_1_data.get("megajson") or global_state.get(self.step_id, {}).get("megajson")

        if not source_trajs and existing:
            source_trajs = _extract_elementary_trajectories_from_megajson(existing)

        trims = junction_trims or {}

        if len(source_trajs) >= 2:
            try:
                megajson = stitch_multistep_trajectories(
                    source_trajs, junction_choices=junction_choices, n_bridge_frames=n_bridge,
                    junction_trims=trims, auto_trim=auto_trim
                )
                if existing:
                    for key in ("solvent", "molecule_bundles", "stoichiometry_summary", "orbital_timeline", "svg_wiring", "mechanism_svg"):
                        if key in existing and key not in megajson:
                            megajson[key] = existing[key]

                frame_energies = [
                    round(float(fr.get("rel_energy_kcal", fr.get("energy_kcal", 0.0))), 2)
                    for fr in megajson.get("frames", [])
                ]
                megajson["energy_profile"] = {
                    "frame_energies_kcal": frame_energies,
                    "overall_barrier_kcal": megajson.get("overall_barrier_kcal", 0.0),
                    "overall_reaction_energy_kcal": megajson.get("overall_reaction_energy_kcal", 0.0)
                }

                if self.step_id not in global_state:
                    global_state[self.step_id] = {}
                global_state[self.step_id]["megajson"] = megajson
                global_state[self.step_id]["junction_trims"] = trims
                step_1_data["megajson"] = megajson

                return {
                    "success": True,
                    "megajson": megajson,
                    "report": megajson.get("stitching_alignment_report", []),
                    "junction_previews": _ensure_junction_previews(megajson),
                    "junction_choices": junction_choices,
                    "junction_trims": trims,
                    "n_bridge_frames": n_bridge
                }
            except Exception as e:
                return {"error": f"Stitching failed: {e}"}

        if existing:
            return {
                "success": True,
                "megajson": existing,
                "report": existing.get("stitching_alignment_report", []),
                "junction_previews": _ensure_junction_previews(existing),
                "junction_choices": junction_choices,
                "junction_trims": trims,
                "n_bridge_frames": n_bridge
            }

        return {"error": "Need at least 2 trajectories or a loaded cascade trajectory from Step 1."}

    def _handle_preview_stitch(self, global_state: Dict[str, Any], user_input: Dict[str, Any] = None) -> Dict[str, Any]:
        user_input = user_input or {}
        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))
        junction_trims = user_input.get("junction_trims", {})
        auto_trim = bool(user_input.get("auto_trim", False))
        return self._execute_stitch(global_state, junction_choices, n_bridge, junction_trims=junction_trims, auto_trim=auto_trim)

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")
        if action_type in ("preview_stitch", "preview_alignment", "recalculate", "apply_trim", "auto_trim"):
            return self._handle_preview_stitch(global_state, user_input)
        return {"error": f"Unknown action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        user_input = user_input or {}
        junction_choices = user_input.get("junction_choices", {})
        n_bridge = int(user_input.get("n_bridge_frames", 8))
        junction_trims = user_input.get("junction_trims", {})
        auto_trim = bool(user_input.get("auto_trim", False))
        res = self._execute_stitch(global_state, junction_choices, n_bridge, junction_trims=junction_trims, auto_trim=auto_trim)
        if not res.get("success"):
            return StepResult(success=False, error_message=res.get("error", "Stitching failed."))
        return StepResult(
            success=True,
            data={
                "megajson": res["megajson"],
                "report": res.get("report", []),
                "junction_previews": res.get("junction_previews", []),
                "junction_choices": junction_choices,
                "junction_trims": junction_trims,
                "n_bridge_frames": n_bridge
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("megajson"))
