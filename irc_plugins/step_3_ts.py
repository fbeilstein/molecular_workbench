import os
from typing import Any, Dict
from .base import PipelineStep, StepResult
from rdkit import Chem
from rdkit.Chem import AllChem

OUT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'out'))
os.makedirs(OUT_DIR, exist_ok=True)
LOG_PATH = os.path.join(OUT_DIR, 'optimization.log')


def _evaluate_frequency_mode(atoms, model, log):
    try:
        from .ts_builder import analyse_mode
        return analyse_mode(atoms, model, log)
    except Exception as e:
        log(f"  frequency check note: {e}")
        return None


def _format_frequency_verdict(ts_report):
    if not ts_report:
        return ""
    if ts_report.get("n_imag") == 1 and ts_report.get("mode_consistent", True):
        return f" | 1 imag freq ({ts_report['imag_freqs'][0]:.0f}i cm-1) - TS confirmed"
    if ts_report.get("n_imag") == 0:
        return " | ⚠️ 0 imag freqs (Local Minimum - Not a TS)"
    return f" | ⚠️ {ts_report.get('n_imag', 0)} imag freq(s) - check log"


class Step3TS(PipelineStep):
    _dock_cache: Dict[Any, Any] = {}

    @property
    def step_id(self) -> str:
        return "step_3_ts"

    @property
    def step_name(self) -> str:
        return "3. 3D Editor & TS Guess"

    def _generate_3d_molblock(self, smiles: str) -> str:
        from .ts_engine.docking import generate_lowest_energy_molblock
        return generate_lowest_energy_molblock(smiles)

    @staticmethod
    def _model(state: Dict[str, Any]):
        from .ts_builder import ReactionModel
        raw = state.get("step_1_input", {}).get("raw_smiles", "")
        if not raw or ">>" not in raw:
            raise ValueError(f"No valid reaction SMILES in state: '{raw}'")
        mapping = state.get("step_2_mapping", {}).get("final_mapping", {})
        if not mapping:
            from .step_2_mapping import Step2Mapping
            mapping = Step2Mapping().get_ui_schema(state).get("initial_mapping", {})
        return ReactionModel(raw, mapping, solvent=state.get("solvent", "water"))

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        mapping_data = current_state.get("step_2_mapping")
        if not mapping_data:
            return {"message": "Please complete Step 2 (Mapping) first."}

        raw_smiles = current_state.get("step_1_input", {}).get("raw_smiles", "")
        if ">>" not in raw_smiles:
            return {"message": "Invalid reaction smiles."}

        r_part, p_part = raw_smiles.split(">>")
        r_smiles_list = r_part.split(".") if r_part else []
        p_smiles_list = p_part.split(".") if p_part else []
        p_blocks = [self._generate_3d_molblock(s) for s in p_smiles_list]
        r_clean_blocks = self._side_by_side([self._generate_3d_molblock(s) for s in r_smiles_list])

        # Automatic docking from the mapping (cached: /api/plugins is polled often)
        mapping = mapping_data.get("final_mapping", {})
        if not mapping:
            from .step_2_mapping import Step2Mapping
            mapping = Step2Mapping().get_ui_schema(current_state).get("initial_mapping", {})
        key = (raw_smiles, tuple(sorted(mapping.items())))
        try:
            if key not in self._dock_cache:
                from .ts_builder import auto_dock_molblocks
                model = self._model(current_state)
                blocks, n_var, _ = auto_dock_molblocks(model, variant=0)
                self._dock_cache.clear()
                self._dock_cache[key] = (blocks, n_var, model.reaction_info())
            r_blocks, n_var, info = self._dock_cache[key]
            msg = ("Reactants docked from atom mapping. "
                   "Choose an optimization strategy or refine fragments if needed.")
        except Exception as e:
            r_blocks = r_clean_blocks
            n_var, info = 1, {"forming": [], "breaking": [], "order_changed": [], "summary": [f"Auto-dock note: {e}"]}
            msg = ("Auto-dock could not align fragments automatically. "
                   "Position them in 3D and select an optimization strategy.")

        return {
            "type": "ts_editor_ui",
            "reactants_3d": r_blocks,
            "clean_reactants_3d": r_clean_blocks,
            "products_3d": p_blocks,
            "reaction_info": info,
            "n_dock_variants": n_var,
            "solvent": current_state.get("solvent", "water"),
            "message": msg,
        }

    @staticmethod
    def _side_by_side(blocks, gap=2.5):
        """Fallback layout: separate fragments along x instead of overlapping."""
        import numpy as np
        from .ts_builder import set_positions
        out, x_cursor = [], 0.0
        for mb in blocks:
            m = Chem.MolFromMolBlock(mb, removeHs=False) if mb else None
            if m is None:
                out.append(mb)
                continue
            p = m.GetConformer().GetPositions()
            p = p - p.mean(axis=0)
            half = float(np.max(np.abs(p[:, 0]))) if len(p) else 0.0
            p[:, 0] += x_cursor + half
            x_cursor += 2 * half + gap
            out.append(Chem.MolToMolBlock(set_positions(m, p)))
        return out

    # ------------------------------------------------------------------
    def _handle_auto_dock(self, global_state, user_input):
        try:
            from .ts_builder import auto_dock_molblocks
            model = self._model(global_state)
            model.apply_override(user_input.get("reaction_override"))
            variant = int(user_input.get("variant", 0))
            blocks, n_var, score = auto_dock_molblocks(model, variant=variant)
            return {"success": True, "molblocks": blocks, "n_variants": n_var,
                    "message": f"Docking pose {variant % n_var + 1}/{n_var} (score {score:.2f})"}
        except Exception as e:
            return {"error": f"Auto-dock failed: {e}"}

    def _handle_refine_fragment(self, global_state, user_input):
        try:
            from .ts_builder import refine_fragment_shape
            frag_idx = int(user_input.get("fragment_index", 0))
            mb = user_input.get("molblock")
            method = user_input.get("method", "mmff")
            model = self._model(global_state)
            chg, mult = 0, 1
            if frag_idx < len(model.r_mols):
                m = model.r_mols[frag_idx]
                chg = int(sum(a.GetFormalCharge() for a in m.GetAtoms()))
                mult = int(sum(a.GetNumRadicalElectrons() for a in m.GetAtoms())) + 1
            new_mb = refine_fragment_shape(mb, method=method, charge=chg, mult=mult)
            return {"success": True, "molblock": new_mb, "fragment_index": frag_idx,
                    "message": f"Reactant {frag_idx + 1} refined via {method.upper()} in-place."}
        except Exception as e:
            return {"error": f"Fragment refinement failed: {e}"}

    def _handle_reset_clean(self, global_state, user_input):
        try:
            from .ts_builder import auto_dock_molblocks
            model = self._model(global_state)
            blocks, _, _ = auto_dock_molblocks(model, variant=0)
            return {"success": True, "molblocks": blocks,
                    "message": "Reset scene to clean reactant conformers."}
        except Exception:
            raw_smiles = global_state.get("step_1_input", {}).get("raw_smiles", "")
            r_part = raw_smiles.split(">>")[0]
            r_smiles_list = r_part.split(".") if r_part else []
            clean_blocks = self._side_by_side([self._generate_3d_molblock(s) for s in r_smiles_list])
            return {"success": True, "molblocks": clean_blocks,
                    "message": "Reset scene to clean reactant conformers."}

    def _handle_seed_from_previous(self, global_state, user_input):
        import glob, json
        from .ts_engine.trajectory_stitcher import extract_fragment_molblock
        cand_files = [p for p in [os.path.join(OUT_DIR, "irc_trajectory_prev.json"), os.path.join(OUT_DIR, "irc_trajectory.json")] if os.path.exists(p)]
        t_dir = os.path.join(os.path.dirname(OUT_DIR), "trajectories")
        if os.path.exists(t_dir):
            cand_files.extend(sorted(glob.glob(os.path.join(t_dir, "*.json")), reverse=True))
        if not cand_files:
            return {"error": "No previous IRC trajectory found in out/ or trajectories/."}
        try:
            raw_smiles = global_state.get("step_1_input", {}).get("raw_smiles", "")
            r_smiles = raw_smiles.split(">>")[0].split(".") if ">>" in raw_smiles else []
            current_blocks = list(user_input.get("molblocks", []))
            if not current_blocks:
                from .ts_builder import auto_dock_molblocks
                current_blocks, _, _ = auto_dock_molblocks(self._model(global_state), variant=0)
            for c_path in cand_files:
                try:
                    with open(c_path, "r") as f:
                        traj = json.load(f)
                    prods = traj.get("fragments", {}).get("products", [])
                    seeded = 0
                    for r_idx, r_smi in enumerate(r_smiles):
                        can_r = Chem.CanonSmiles(r_smi) if r_smi else ""
                        for fp in prods:
                            can_p = Chem.CanonSmiles(fp.get("smiles", "")) if fp.get("smiles") else ""
                            if can_r and can_r == can_p:
                                mb = extract_fragment_molblock(traj, fp)
                                if r_idx < len(current_blocks):
                                    current_blocks[r_idx] = mb
                                else:
                                    current_blocks.append(mb)
                                seeded += 1
                                break
                    if seeded > 0:
                        src = os.path.basename(c_path)
                        return {"success": True, "molblocks": current_blocks, "message": f"Seeded {seeded} reactant(s) from {src}!"}
                except Exception:
                    continue
            return {"error": "None of the current reactants match products in previous trajectories."}
        except Exception as e:
            return {"error": f"Seeding failed: {e}"}

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")
        handlers = {
            "auto_dock": self._handle_auto_dock,
            "refine_fragment": self._handle_refine_fragment,
            "reset_clean": self._handle_reset_clean,
            "seed_from_previous": self._handle_seed_from_previous,
            "optimize_ts": lambda gs, ui: self._optimize_ts(ui, gs),
        }
        handler = handlers.get(action_type)
        if not handler:
            return {"error": f"Unknown custom action: {action_type}"}
        return handler(global_state, user_input)

    @staticmethod
    def _run_strategy(strategy, atoms, model, log):
        from .ts_builder import (
            run_dimer_ts, run_drive_sella_ts, run_seeded_sella_ts, run_neb_ts,
            run_rits_ts, run_lego_ts, run_guided_manual_ts, run_genetic_swarm_ts, run_orca_dft_ts
        )
        dispatch = {
            "rits": lambda: (run_rits_ts(atoms, model, log), None),
            "lego": lambda: (run_lego_ts(atoms, model, log), None),
            "genetic_swarm": lambda: (run_genetic_swarm_ts(atoms, model, log), None),
            "manual_valley": lambda: (run_guided_manual_ts(atoms, model, log), None),
            "auto_dock": lambda: (run_drive_sella_ts(atoms, model, log), None),
            "orca_dft": lambda: (run_orca_dft_ts(atoms, model, log), None),
            "dimer": lambda: ([], run_dimer_ts(atoms, model, log)),
            "drive": lambda: (run_drive_sella_ts(atoms, model, log), None),
            "sella": lambda: (run_seeded_sella_ts(atoms, model, log), None),
            "neb": lambda: (run_neb_ts(atoms, model, log), None),
        }
        runner = dispatch.get(strategy, dispatch["rits"])
        profile, _ = runner()
        return profile or []

    def _optimize_ts(self, user_input, global_state):
        import numpy as np
        from ase import Atoms
        from .ts_builder import (molblocks_to_positions, positions_to_molblocks,
                                 check_geometry_sanity, analyse_mode)

        docked = user_input.get("docked_molblocks", [])
        mode = user_input.get("mode", "auto")
        strategy = user_input.get("strategy") or "rits"

        if not docked:
            return {"error": "No molblocks provided"}

        raw = global_state.get("step_1_input", {}).get("raw_smiles", "")
        if not raw or ">>" not in raw:
            return {"error": "Reaction SMILES missing. Please complete Step 1 (Input) and Step 2 (Mapping) first."}

        log_fh = open(LOG_PATH, 'w', buffering=1)

        def log(s=""):
            log_fh.write(s + "\n")
            log_fh.flush()

        try:
            if user_input.get("solvent"):
                global_state["solvent"] = user_input["solvent"]
            model = self._model(global_state)
            overridden = model.apply_override(user_input.get("reaction_override"))
            mols, pos_list = molblocks_to_positions(docked)
            positions = np.vstack(pos_list)
            symbols = [a.GetSymbol() for m in mols for a in m.GetAtoms()]
            if len(symbols) != model.n_atoms or symbols != model.symbols:
                return {"error": "3D atoms do not match the mapped reaction (re-run Step 3)."}
            atoms = Atoms(symbols=symbols, positions=positions)

            src = "MANUAL selection in 3D" if overridden else "Step 2 mapping"
            log(f"=== Reaction analysis (reacting bonds from {src}) ===")
            for line in model.summary_lines():
                log(line)

            log(f"\nOptimization Strategy: {strategy.upper()}")
            profile = self._run_strategy(strategy, atoms, model, log)

            # Sanity verification: protect against bond tearing / unphysical distortion
            is_sane = check_geometry_sanity(atoms, model, log)
            if not is_sane:
                log("\n❌ SANITY GATE TRIGGERED: Spectator bonds were broken or severely stretched.")
                log("   Retaining original input coordinates to protect your setup.")
                return {
                    "success": False,
                    "diverged": True,
                    "error": "Optimization distorted spectator bonds beyond sanity threshold. Previous coordinates preserved.",
                    "profile": profile,
                }

            if atoms.calc is None:
                from .ts_engine.calculators import make_xtb
                atoms.calc = make_xtb(model.charge, model.multiplicity, solvent=getattr(model, "solvent", "water"))
            try:
                e_ts = atoms.get_potential_energy()
            except Exception:
                e_ts = 0.0
            log(f"\nFinal geometry energy: E = {e_ts:.4f} eV")
            for i, j in model.formed + model.broken:
                d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
                log(f"  {model.label(i)}-{model.label(j)}: {d:.3f} A")
            for i, j, ro, po in model.order_changed:
                d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
                log(f"  Order changed {model.label(i)}-{model.label(j)} ({ro:.1f}->{po:.1f}): {d:.3f} A")

            log("\n=== Frequency check ===")
            ts_report = _evaluate_frequency_mode(atoms, model, log)
            opt_blocks = positions_to_molblocks(mols, atoms.get_positions())
            verdict = _format_frequency_verdict(ts_report)

            return {
                "success": True,
                "optimized_molblocks": opt_blocks,
                "profile": profile,
                "ts_report": ts_report,
                "message": f"TS optimization finished ({strategy}). E = {e_ts:.3f} eV{verdict}",
            }
        except Exception as e:
            import traceback
            log("\nERROR: " + traceback.format_exc())
            return {"error": f"Optimization failed: {e}"}
        finally:
            log_fh.close()

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        return StepResult(success=True, data={"docked_reactants": user_input.get("docked_molblocks", [])})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "docked_reactants" in step_data

