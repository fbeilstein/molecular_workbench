import json
import os
from typing import Any, Dict
import numpy as np
from ase import Atoms
from .base import PipelineStep, StepResult
from .ts_builder import (
    ReactionModel,
    molblocks_to_positions,
    run_irc_path
)

LOG_PATH = 'optimization.log'
TRAJECTORY_FILE = 'irc_trajectory.json'


class Step4IRC(PipelineStep):
    @property
    def step_id(self) -> str:
        return "step_4_irc"

    @property
    def step_name(self) -> str:
        return "4. IRC Trajectory"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        existing_traj = current_state.get("step_4_irc", {}).get("trajectory")
        if not existing_traj and os.path.exists(TRAJECTORY_FILE):
            try:
                with open(TRAJECTORY_FILE, 'r') as f:
                    existing_traj = json.load(f)
            except Exception:
                pass

        ts_step = current_state.get("step_3_ts")
        ts_molblocks = ts_step.get("docked_reactants", []) if ts_step else []

        if not ts_molblocks and not existing_traj:
            return {
                "message": "Please complete Step 3 (Transition State Optimization) first.",
                "type": "irc_viewer_ui",
                "ts_molblocks": [],
                "trajectory": None
            }

        return {
            "type": "irc_viewer_ui",
            "ts_molblocks": ts_molblocks,
            "trajectory": existing_traj,
            "solvent": current_state.get("solvent", "water"),
            "message": "Step 4: Intrinsic Reaction Coordinate (IRC) path connecting Reactants -> TS -> Products."
        }

    def _calculate_irc(self, user_input: Dict[str, Any], global_state: Dict[str, Any]) -> Dict[str, Any]:
        ts_blocks = user_input.get("ts_molblocks") or global_state.get("step_3_ts", {}).get("docked_reactants", [])
        if not ts_blocks:
            return {"error": "No Transition State geometry found. Complete Step 3 first."}

        raw_smiles = global_state.get("step_1_input", {}).get("raw_smiles", "")
        if not raw_smiles or ">>" not in raw_smiles:
            return {"error": "Reaction SMILES missing from state."}

        max_steps = int(user_input.get("max_steps", 30))
        dx = float(user_input.get("step_size", 0.08))
        solvent = user_input.get("solvent") or global_state.get("solvent", "water")
        relax_fragments = bool(user_input.get("relax_fragments", True))

        log_fh = open(LOG_PATH, 'w', buffering=1)

        def log(s=""):
            log_fh.write(s + "\n")
            log_fh.flush()

        try:
            mapping = global_state.get("step_2_mapping", {}).get("final_mapping", {})
            model = ReactionModel(raw_smiles, mapping, solvent=solvent)
            mols, pos_list = molblocks_to_positions(ts_blocks)
            positions = np.vstack(pos_list)
            symbols = [a.GetSymbol() for m in mols for a in m.GetAtoms()]

            atoms = Atoms(symbols=symbols, positions=positions)
            log(f"=== Starting IRC Calculation for {raw_smiles} ===")
            log(f"  System: {len(symbols)} atoms, charge = {model.charge}, multiplicity = {model.multiplicity}")
            log(f"  Environment: {solvent.upper() if solvent else 'GAS PHASE'}")

            traj = run_irc_path(atoms, model, log, max_steps=max_steps, dx=dx, relax_fragments=relax_fragments)

            # Persist trajectory JSON to root
            try:
                with open(TRAJECTORY_FILE, 'w') as f:
                    json.dump(traj, f, indent=2)
            except Exception as err:
                log(f"  Note: Failed to write {TRAJECTORY_FILE}: {err}")

            return {
                "success": True,
                "trajectory": traj,
                "message": f"IRC path generated: {len(traj['frames'])} frames (Forward barrier: {traj['barriers']['forward_barrier_kcal']:.1f} kcal/mol, Reverse barrier: {traj['barriers']['reverse_barrier_kcal']:.1f} kcal/mol)"
            }
        except Exception as e:
            import traceback
            log("\nIRC ERROR: " + traceback.format_exc())
            return {"error": f"IRC calculation failed: {e}"}
        finally:
            log_fh.close()

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")
        if action_type == "calculate_irc":
            return self._calculate_irc(user_input, global_state)
        return {"error": f"Unknown custom action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        traj = user_input.get("trajectory")
        return StepResult(success=True, data={"trajectory": traj})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("trajectory"))
