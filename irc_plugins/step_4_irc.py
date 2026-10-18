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

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
os.makedirs(OUT_DIR, exist_ok=True)
LOG_PATH = os.path.join(OUT_DIR, 'optimization.log')
TRAJECTORY_FILE = os.path.join(OUT_DIR, 'irc_trajectory.json')
LEGACY_TRAJECTORY_FILE = os.path.join(ROOT_DIR, 'irc_trajectory.json')


class Step4IRC(PipelineStep):
    @property
    def step_id(self) -> str:
        return "step_4_irc"

    @property
    def step_name(self) -> str:
        return "4. IRC Trajectory"

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        existing_traj = current_state.get("step_4_irc", {}).get("trajectory")
        if not existing_traj:
            for tf in (TRAJECTORY_FILE, LEGACY_TRAJECTORY_FILE):
                if os.path.exists(tf):
                    try:
                        with open(tf, 'r') as f:
                            existing_traj = json.load(f)
                            break
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
        engine = user_input.get("engine") or "xtb"
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
            log(f"  Engine: {engine.upper()} | Environment: {solvent.upper() if solvent else 'GAS PHASE'}")

            traj = run_irc_path(atoms, model, log, max_steps=max_steps, dx=dx, relax_fragments=relax_fragments, engine=engine)

            # Persist trajectory JSON with previous backup
            try:
                import shutil
                if os.path.exists(TRAJECTORY_FILE):
                    shutil.copy2(TRAJECTORY_FILE, os.path.join(OUT_DIR, "irc_trajectory_prev.json"))
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

    def _handle_prolong_trajectory(self, user_input: Dict[str, Any], global_state: Dict[str, Any]) -> Dict[str, Any]:
        import glob
        from .ts_engine.trajectory_stitcher import stitch_multistep_trajectories
        current_traj = user_input.get("trajectory") or global_state.get("step_4_irc", {}).get("trajectory")
        if not current_traj and os.path.exists(TRAJECTORY_FILE):
            with open(TRAJECTORY_FILE, "r") as f:
                current_traj = json.load(f)
        if not current_traj:
            return {"error": "Current IRC trajectory missing."}

        cand_files = [p for p in [os.path.join(OUT_DIR, "irc_trajectory_prev.json")] if os.path.exists(p)]
        t_dir = os.path.join(os.path.dirname(OUT_DIR), "trajectories")
        if os.path.exists(t_dir):
            cand_files.extend(sorted(glob.glob(os.path.join(t_dir, "*.json"))))

        prev_traj = None
        cur_smiles = current_traj.get("reaction_smiles")
        for c_path in cand_files:
            try:
                with open(c_path, "r") as f:
                    tj = json.load(f)
                if tj.get("reaction_smiles") != cur_smiles and "frames" in tj:
                    prev_traj = tj
                    break
            except Exception:
                continue
        if not prev_traj:
            return {"error": "No previous IRC trajectory found to prolong onto."}
        try:
            mega = stitch_multistep_trajectories([prev_traj, current_traj], n_bridge_frames=8)
            with open(TRAJECTORY_FILE, "w") as f:
                json.dump(mega, f, indent=2)
            prolong_p = os.path.join(os.path.dirname(OUT_DIR), "trajectories", "prolonged_pathway.json")
            try:
                with open(prolong_p, "w") as f:
                    json.dump(mega, f, indent=2)
            except Exception:
                pass
            return {
                "success": True,
                "trajectory": mega,
                "message": f"Pathway prolonged! Master trajectory now has {mega.get('total_frames', len(mega['frames']))} frames across {mega.get('n_steps', 2)} steps."
            }
        except Exception as e:
            return {"error": f"Failed to prolong trajectory: {e}"}

    def _handle_prolong_reaction(self, user_input: Dict[str, Any], global_state: Dict[str, Any]) -> Dict[str, Any]:
        """Take product of this reaction and seed it as reactant for next elementary reaction."""
        import glob
        from .ts_engine.trajectory_stitcher import extract_fragment_molblock
        current_traj = user_input.get("trajectory") or global_state.get("step_4_irc", {}).get("trajectory")
        if not current_traj and os.path.exists(TRAJECTORY_FILE):
            try:
                with open(TRAJECTORY_FILE, "r") as f:
                    current_traj = json.load(f)
            except Exception:
                pass
        if not current_traj:
            return {"error": "Current IRC trajectory missing."}

        t_dir = os.path.join(os.path.dirname(OUT_DIR), "trajectories")
        os.makedirs(t_dir, exist_ok=True)
        existing_steps = [p for p in sorted(glob.glob(os.path.join(t_dir, "*.json"))) if not p.endswith("prolonged_pathway.json")]
        step_num = len(existing_steps) + 1
        step_path = os.path.join(t_dir, f"irc_trajectory{step_num}.json")
        try:
            with open(step_path, "w") as f:
                json.dump(current_traj, f, indent=2)
        except Exception:
            pass

        prods = current_traj.get("fragments", {}).get("products", [])
        if not prods:
            return {"error": "No products identified in current trajectory."}

        prod_smiles = [p.get("smiles") for p in prods if p.get("smiles")]
        seeded_blocks = [extract_fragment_molblock(current_traj, p) for p in prods]
        next_smiles = ".".join(prod_smiles) + ">>"

        global_state["step_1_input"] = {"raw_smiles": next_smiles, "raw_rxnfile": ""}
        global_state.pop("step_2_mapping", None)
        global_state["step_3_ts"] = {"docked_reactants": seeded_blocks}
        global_state.pop("step_4_irc", None)

        state_file = os.path.join(OUT_DIR, "reaction_state.json")
        try:
            with open(state_file, "w") as f:
                json.dump(global_state, f, indent=2)
        except Exception:
            pass

        return {
            "success": True,
            "next_reactant_smiles": ".".join(prod_smiles),
            "message": f"Step {step_num} archived! Seeded product ({'.'.join(prod_smiles)}) for next reaction."
        }

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")
        if action_type == "calculate_irc":
            return self._calculate_irc(user_input, global_state)
        if action_type in ("prolong_to_next", "prolong_reaction"):
            return self._handle_prolong_reaction(user_input, global_state)
        if action_type == "prolong_trajectory":
            return self._handle_prolong_trajectory(user_input, global_state)
        return {"error": f"Unknown custom action: {action_type}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        traj = user_input.get("trajectory")
        return StepResult(success=True, data={"trajectory": traj})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("trajectory"))
