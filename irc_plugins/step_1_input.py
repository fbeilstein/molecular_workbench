import os
import json
import shutil
from typing import Any, Dict
from .base import PipelineStep, StepResult
from .project_manager import (
    decompose_multistep_reaction,
    init_project,
    save_project_meta,
    _step_state_path,
    _step_traj_path,
    STEPS_DIR
)


class Step1Input(PipelineStep):
    @property
    def step_id(self) -> str:
        return "step_1_input"

    @property
    def step_name(self) -> str:
        return "1. Reaction Input"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s1 = current_state.get(self.step_id, {})
        return {
            "type": "ketcher",
            "message": "Draw the reactants on the left and products on the right, separated by a reaction arrow. Multi-arrow reactions (A >> B >> C) are automatically decomposed into steps!",
            "raw_smiles": s1.get("raw_smiles", ""),
            "raw_rxnfile": s1.get("raw_rxnfile", "")
        }

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        smiles = user_input.get("smiles", "")
        rxnfile = user_input.get("rxnfile", "")

        if not smiles or ">>" not in smiles:
            return StepResult(success=False, error_message="Please draw a valid reaction (with an arrow).")

        subreactions = decompose_multistep_reaction(smiles, rxnfile)
        if len(subreactions) > 1:
            proj = init_project()
            solvent = global_state.get("solvent", proj.get("solvent", "water"))

            # Remove any excess step directories from older runs
            existing_dirs = [d for d in os.listdir(STEPS_DIR) if d.startswith("step_")]
            for d in existing_dirs:
                try:
                    idx = int(d.replace("step_", ""))
                    if idx > len(subreactions):
                        shutil.rmtree(os.path.join(STEPS_DIR, d), ignore_errors=True)
                except ValueError:
                    pass

            new_steps = []
            for i, sub in enumerate(subreactions, start=1):
                sid = f"step_{i}"
                s_dir = os.path.join(STEPS_DIR, sid)
                os.makedirs(s_dir, exist_ok=True)
                s_state_file = _step_state_path(sid)
                
                cur_st = {}
                if os.path.exists(s_state_file):
                    try:
                        with open(s_state_file, 'r') as f:
                            cur_st = json.load(f)
                    except Exception:
                        pass
                
                cur_st["solvent"] = solvent
                cur_st["step_1_input"] = {
                    "raw_smiles": sub["smiles"],
                    "raw_rxnfile": sub["rxnfile"]
                }
                # Reset downstream steps for this subreaction
                for downstream in ["step_2_mapping", "step_3_ts", "step_4_irc"]:
                    cur_st.pop(downstream, None)

                # Remove old trajectory if present
                t_path = _step_traj_path(sid)
                if os.path.exists(t_path):
                    try:
                        os.remove(t_path)
                    except Exception:
                        pass

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
            save_project_meta(proj)

            return StepResult(success=True, data={
                "raw_smiles": subreactions[0]["smiles"],
                "raw_rxnfile": subreactions[0]["rxnfile"]
            })

        return StepResult(success=True, data={
            "raw_smiles": smiles,
            "raw_rxnfile": rxnfile
        })

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "raw_smiles" in step_data
