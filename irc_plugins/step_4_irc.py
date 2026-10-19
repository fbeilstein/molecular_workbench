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
        from .project_manager import init_project, _step_traj_path
        proj = init_project()
        active_id = proj.get("active_step_id", "step_1")
        step_traj_file = _step_traj_path(active_id)

        existing_traj = current_state.get("step_4_irc", {}).get("trajectory")
        if existing_traj and (existing_traj.get("n_steps", 1) > 1 or existing_traj.get("title") == "Master Multistep IRC Trajectory"):
            existing_traj = None

        if not existing_traj:
            cands = [step_traj_file]
            if len(proj.get("steps", [])) <= 1:
                cands.extend([TRAJECTORY_FILE, LEGACY_TRAJECTORY_FILE])
            for tf in cands:
                if os.path.exists(tf):
                    try:
                        with open(tf, 'r') as f:
                            tj = json.load(f)
                        if tj and tj.get("frames") and tj.get("n_steps", 1) == 1 and tj.get("title") != "Master Multistep IRC Trajectory":
                            existing_traj = tj
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
        engine = str(user_input.get("engine") or "xtb").lower()
        orca_mode = str(user_input.get("orca_mode") or "sp").lower()
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
            log(f"  Engine: {engine.upper()} | Mode: {orca_mode.upper() if engine == 'orca' else 'Standard'} | Solvent: {solvent.upper() if solvent else 'GAS PHASE'}")

            if engine == "orca":
                from .ts_engine.orca_runner import (
                    run_orca_thermo_profile,
                    run_orca_native_irc,
                    run_orca_neb_ci
                )
                nprocs = min(4, os.cpu_count() or 4)
                if orca_mode == "thermo":
                    traj = run_orca_thermo_profile(atoms, model, solvent, log, max_steps=max_steps, dx=dx, nprocs=nprocs)
                elif orca_mode == "irc":
                    traj = run_orca_native_irc(atoms, model, solvent, log, max_points=max_steps, nprocs=nprocs)
                elif orca_mode == "neb":
                    traj = run_orca_neb_ci(atoms, model, solvent, log, n_images=16, nprocs=nprocs)
                else: # "sp"
                    traj = run_irc_path(atoms, model, log, max_steps=max_steps, dx=dx, relax_fragments=relax_fragments, engine="orca")
            else:
                traj = run_irc_path(atoms, model, log, max_steps=max_steps, dx=dx, relax_fragments=relax_fragments, engine="xtb")

            if not traj or not isinstance(traj, dict) or not traj.get("frames"):
                return {"error": "Failed to calculate IRC trajectory. Check log for details."}

            # Persist trajectory JSON with previous backup
            try:
                import shutil
                if os.path.exists(TRAJECTORY_FILE):
                    shutil.copy2(TRAJECTORY_FILE, os.path.join(OUT_DIR, "irc_trajectory_prev.json"))
                with open(TRAJECTORY_FILE, 'w') as f:
                    json.dump(traj, f, indent=2)

                # Persist to active project step directory
                from .project_manager import init_project, _step_traj_path
                proj = init_project()
                active_id = proj.get("active_step_id", "step_1")
                step_traj_file = _step_traj_path(active_id)
                with open(step_traj_file, 'w') as f:
                    json.dump(traj, f, indent=2)

                # Also backup to trajectories directory with step index
                step_idx = 1
                for s in proj.get("steps", []):
                    if s["id"] == active_id:
                        step_idx = s.get("index", 1)
                        break
                traj_dir = os.path.join(ROOT_DIR, "trajectories")
                os.makedirs(traj_dir, exist_ok=True)
                with open(os.path.join(traj_dir, f"irc_trajectory{step_idx}.json"), 'w') as f:
                    json.dump(traj, f, indent=2)
            except Exception as err:
                log(f"  Note: Failed to write {TRAJECTORY_FILE}: {err}")

            is_gibbs = traj.get("barriers", {}).get("energy_type") == "gibbs"
            prefix = "Thermochemical ΔG" if is_gibbs else "IRC"
            b_unit = "ΔG‡" if is_gibbs else "ΔE‡"
            fwd_b = traj.get('barriers', {}).get('forward_barrier_kcal', 0.0)
            rev_b = traj.get('barriers', {}).get('reverse_barrier_kcal', 0.0)

            return {
                "success": True,
                "trajectory": traj,
                "message": f"{prefix} path generated: {len(traj['frames'])} frames ({b_unit} fwd: {fwd_b:.1f} kcal/mol, {b_unit} rev: {rev_b:.1f} kcal/mol)"
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
                if tj.get("reaction_smiles") and tj.get("reaction_smiles") != cur_smiles and "frames" in tj and tj.get("fragments"):
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
