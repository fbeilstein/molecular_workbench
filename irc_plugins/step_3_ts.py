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
    _mb_cache: Dict[str, str] = {}

    @property
    def step_id(self) -> str:
        return "step_3_ts"

    @property
    def step_name(self) -> str:
        return "3. 3D Editor & TS Guess"

    def _generate_3d_molblock(self, smiles: str) -> str:
        if smiles not in self._mb_cache:
            from .ts_engine.docking import generate_lowest_energy_molblock
            self._mb_cache[smiles] = generate_lowest_energy_molblock(smiles)
        return self._mb_cache[smiles]

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
            return {"type": "empty", "message": "Please complete Step 2 (Mapping) first."}

        raw_smiles = current_state.get("step_1_input", {}).get("raw_smiles", "")
        if ">>" not in raw_smiles:
            return {"type": "empty", "message": "Invalid reaction smiles. Please complete Step 1 first."}

        r_part, p_part = raw_smiles.split(">>")
        r_smiles_list = r_part.split(".") if r_part else []
        p_smiles_list = p_part.split(".") if p_part else []
        p_blocks = [self._generate_3d_molblock(s) for s in p_smiles_list]
        r_clean_blocks = self._side_by_side([self._generate_3d_molblock(s) for s in r_smiles_list])

        existing_docked = current_state.get("step_3_ts", {}).get("docked_reactants", [])
        if existing_docked:
            r_blocks, n_var, info = existing_docked, 1, {"forming": [], "breaking": [], "order_changed": [], "summary": []}
            msg = "Using existing docked reactants."
        else:
            mapping = mapping_data.get("final_mapping") or Step2Mapping().get_ui_schema(current_state).get("initial_mapping", {})
            key = (raw_smiles, tuple(sorted(mapping.items())))
            try:
                if key not in self._dock_cache:
                    from .ts_builder import auto_dock_molblocks
                    model = self._model(current_state)
                    blocks, n_var, _ = auto_dock_molblocks(model, variant=0)
                    self._dock_cache.clear()
                    self._dock_cache[key] = (blocks, n_var, model.reaction_info())
                r_blocks, n_var, info = self._dock_cache[key]
                msg = "Reactants docked from atom mapping. Choose an optimization strategy."
            except Exception as e:
                r_blocks, n_var = r_clean_blocks, 1
                info = {"forming": [], "breaking": [], "order_changed": [], "summary": [f"Auto-dock note: {e}"]}
                msg = "Position fragments in 3D and select an optimization strategy."

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
            if m is None and mb:
                m = Chem.MolFromMolBlock(mb, sanitize=False, removeHs=False)
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
            blocks, _, _ = auto_dock_molblocks(self._model(global_state), variant=0)
            return {"success": True, "molblocks": blocks, "message": "Reset scene to clean reactant conformers."}
        except Exception:
            r_part = global_state.get("step_1_input", {}).get("raw_smiles", "").split(">>")[0]
            clean = self._side_by_side([self._generate_3d_molblock(s) for s in (r_part.split(".") if r_part else [])])
            return {"success": True, "molblocks": clean, "message": "Reset scene to clean reactant conformers."}

    def _handle_seed_from_previous(self, global_state, user_input):
        import glob, json
        from .ts_engine.trajectory_stitcher import extract_fragment_molblock
        from .project_manager import init_project, _step_traj_path

        proj = init_project()
        active_id = proj.get("active_step_id", "step_1")
        steps = proj.get("steps", [])

        active_idx = 1
        for s in steps:
            if s["id"] == active_id:
                active_idx = s.get("index", 1)
                break

        prev_steps = [s for s in steps if s.get("index", 1) < active_idx]
        prev_steps.sort(key=lambda s: s.get("index", 1), reverse=True)
        cand_files = [_step_traj_path(s["id"]) for s in prev_steps if os.path.exists(_step_traj_path(s["id"]))]

        for p in [os.path.join(OUT_DIR, "irc_trajectory_prev.json"), os.path.join(OUT_DIR, "irc_trajectory.json")]:
            if os.path.exists(p) and p not in cand_files:
                cand_files.append(p)
        t_dir = os.path.join(os.path.dirname(OUT_DIR), "trajectories")
        if os.path.exists(t_dir):
            for p in sorted(glob.glob(os.path.join(t_dir, "*.json")), reverse=True):
                if p not in cand_files and not p.endswith("prolonged_pathway.json"):
                    cand_files.append(p)

        if not cand_files:
            return {"error": "No previous IRC trajectory found in project steps or trajectories/."}
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
                                target_m = None
                                try: target_m = self._model(global_state).r_mols[r_idx]
                                except Exception: target_m = Chem.AddHs(Chem.MolFromSmiles(r_smi)) if r_smi else None
                                mb = extract_fragment_molblock(traj, fp, target_mol=target_m)
                                if r_idx < len(current_blocks):
                                    current_blocks[r_idx] = mb
                                else:
                                    current_blocks.append(mb)
                                seeded += 1
                                break
                    if seeded > 0:
                        src = os.path.basename(c_path)
                        global_state.setdefault("step_3_ts", {})["docked_reactants"] = current_blocks
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
        d = {
            "rits": run_rits_ts, "lego": run_lego_ts, "genetic_swarm": run_genetic_swarm_ts,
            "manual_valley": run_guided_manual_ts, "auto_dock": run_drive_sella_ts,
            "orca_dft": run_orca_dft_ts, "drive": run_drive_sella_ts,
            "sella": run_seeded_sella_ts, "neb": run_neb_ts, "dimer": run_dimer_ts
        }
        res = d.get(strategy, run_rits_ts)(atoms, model, log)
        return (res[0] if isinstance(res, tuple) else res) or []

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
            if symbols != model.symbols and sorted(symbols) == sorted(model.symbols) and len(mols) == len(model.r_mols):
                try:
                    reordered = []
                    for r_idx, ref_m in enumerate(model.r_mols):
                        mt = ref_m.GetSubstructMatch(mols[r_idx])
                        if not (mt and len(mt) == ref_m.GetNumAtoms()):
                            from rdkit.Chem import AllChem
                            try:
                                assigned = AllChem.AssignBondOrdersFromTemplate(ref_m, mols[r_idx])
                                mt = ref_m.GetSubstructMatch(assigned)
                            except Exception:
                                pass
                        if mt and len(mt) == ref_m.GetNumAtoms():
                            inv = {ref_i: tr_i for tr_i, ref_i in enumerate(mt)}
                            reordered.append(pos_list[r_idx][[inv[k] for k in range(len(ref_m.GetAtoms()))]])
                    if len(reordered) == len(model.r_mols):
                        positions, symbols, mols = np.vstack(reordered), list(model.symbols), list(model.r_mols)
                except Exception:
                    pass
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
            try: e_ts = atoms.get_potential_energy()
            except Exception: e_ts = 0.0
            log(f"\nFinal geometry energy: E = {e_ts:.4f} eV")
            for i, j in model.formed + model.broken:
                log(f"  {model.label(i)}-{model.label(j)}: {np.linalg.norm(atoms.positions[i] - atoms.positions[j]):.3f} A")
            for i, j, ro, po in model.order_changed:
                log(f"  Order changed {model.label(i)}-{model.label(j)} ({ro:.1f}->{po:.1f}): {np.linalg.norm(atoms.positions[i] - atoms.positions[j]):.3f} A")

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

