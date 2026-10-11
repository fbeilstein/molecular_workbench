from typing import Any, Dict
from .base import PipelineStep, StepResult
from rdkit import Chem
from rdkit.Chem import AllChem

LOG_PATH = 'optimization.log'


class Step3TS(PipelineStep):
    _dock_cache: Dict[Any, Any] = {}

    @property
    def step_id(self) -> str:
        return "step_3_ts"

    @property
    def step_name(self) -> str:
        return "3. 3D Editor & TS Guess"

    def _generate_3d_molblock(self, smiles: str) -> str:
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            return ""
        mol = Chem.AddHs(mol)
        params = AllChem.ETKDGv3()
        params.randomSeed = 42
        if AllChem.EmbedMolecule(mol, params) != 0:
            AllChem.Compute2DCoords(mol)
        else:
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except Exception:
                pass
        return Chem.MolToMolBlock(mol)

    @staticmethod
    def _model(state: Dict[str, Any]):
        from .ts_builder import ReactionModel
        raw = state.get("step_1_input", {}).get("raw_smiles", "")
        mapping = state.get("step_2_mapping", {}).get("final_mapping", {})
        return ReactionModel(raw, mapping)

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        mapping_data = current_state.get("step_2_mapping")
        if not mapping_data:
            return {"message": "Please complete Step 2 (Mapping) first."}

        raw_smiles = current_state.get("step_1_input", {}).get("raw_smiles", "")
        if ">>" not in raw_smiles:
            return {"message": "Invalid reaction smiles."}

        r_part, p_part = raw_smiles.split(">>")
        p_smiles_list = p_part.split(".") if p_part else []
        p_blocks = [self._generate_3d_molblock(s) for s in p_smiles_list]

        # Automatic docking from the mapping (cached: /api/plugins is polled often)
        key = (raw_smiles, tuple(sorted(mapping_data.get("final_mapping", {}).items())))
        try:
            if key not in self._dock_cache:
                from .ts_builder import auto_dock_molblocks
                model = self._model(current_state)
                blocks, n_var, _ = auto_dock_molblocks(model, variant=0)
                self._dock_cache.clear()
                self._dock_cache[key] = (blocks, n_var, model.reaction_info())
            r_blocks, n_var, info = self._dock_cache[key]
            msg = ("Reactants were docked automatically from the atom mapping "
                   "(green dashes = bonds that will form, red = bonds that will break). "
                   "Just press 'Find TS' - manual adjustment is optional.")
        except Exception as e:
            r_smiles_list = r_part.split(".") if r_part else []
            r_blocks = self._side_by_side([self._generate_3d_molblock(s) for s in r_smiles_list])
            n_var, info = 1, {"forming": [], "breaking": [], "summary": [f"Auto-dock failed: {e}"]}
            msg = ("Automatic docking failed - arrange the reactants manually, "
                   "pick the reacting bonds in 3D and use 'Refine my guess'.")

        return {
            "type": "ts_editor_ui",
            "reactants_3d": r_blocks,
            "products_3d": p_blocks,
            "reaction_info": info,
            "n_dock_variants": n_var,
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
    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        action_type = user_input.get("action_type")

        if action_type == "auto_dock":
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

        if action_type == "optimize_ts":
            return self._optimize_ts(user_input, global_state)

        return {"error": "Unknown custom action"}

    def _optimize_ts(self, user_input, global_state):
        import numpy as np
        from ase import Atoms
        from .ts_builder import (molblocks_to_positions, positions_to_molblocks,
                                 bond_drive, refine_guess, refine_ts, analyse_mode,
                                 EV_TO_KCAL, DOCK_EXTRA)

        docked = user_input.get("docked_molblocks", [])
        mode = user_input.get("mode", "auto")
        if not docked:
            return {"error": "No molblocks provided"}

        log_fh = open(LOG_PATH, 'w', buffering=1)

        def log(s=""):
            log_fh.write(s + "\n")
            log_fh.flush()

        try:
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

            ts_report = None
            profile = []
            if mode == "refine":
                log("\n=== Manual mode: your geometry is the TS guess ===")
                for i, j in model.formed:
                    d = np.linalg.norm(positions[i] - positions[j])
                    if d > model.d0(i, j) + DOCK_EXTRA + 0.2:
                        log(f"  WARNING: forming bond {model.label(i)}...{model.label(j)} is "
                            f"{d:.2f} A - this looks like a reactant complex, not a TS. "
                            f"Sella may fall back to a minimum; consider 'Find TS (scan)' "
                            f"or move the fragments closer (~{model.d0(i, j) + 0.45:.1f} A).")
                refine_guess(atoms, model, log)
            elif model.formed or model.broken:
                log("\n=== Bond drive: restrained GFN2-xTB scan R -> P ===")
                frames, energies, lams = bond_drive(atoms, model, log)
                profile = [[l, (e - energies[0]) * EV_TO_KCAL] for l, e in zip(lams, energies)]
                i_max = int(np.argmax(energies))
                if i_max == 0:
                    log("  WARNING: energy only decreases along the drive (barrierless?). "
                        "Using the second point.")
                    i_max = min(1, len(frames) - 1)
                elif i_max == len(frames) - 1:
                    log("  WARNING: maximum at the end of the scan - guess may be product-like.")
                atoms.set_positions(frames[i_max])
                log(f"  -> TS guess at lambda={lams[i_max]:.2f}, "
                    f"barrier estimate {(energies[i_max] - energies[0]) * EV_TO_KCAL:.1f} kcal/mol")
            else:
                log("\nWARNING: no bond changes found; running Sella directly from the input geometry.")

            log("\n=== Sella saddle-point refinement (order=1) ===")
            converged = refine_ts(atoms, model, log_fh, log)
            e_ts = atoms.get_potential_energy()
            log(f"  Sella {'converged' if converged else 'did NOT converge'}; E = {e_ts:.4f} eV")
            for i, j in model.formed + model.broken:
                d = np.linalg.norm(atoms.positions[i] - atoms.positions[j])
                log(f"  {model.label(i)}-{model.label(j)}: {d:.3f} A")

            log("\n=== Frequency check ===")
            try:
                ts_report = analyse_mode(atoms, model, log)
            except Exception as e:
                log(f"  frequency check failed: {e}")

            opt_blocks = positions_to_molblocks(mols, atoms.get_positions())
            verdict = ""
            if ts_report:
                if ts_report.get("n_imag") == 1 and ts_report.get("mode_consistent", True):
                    verdict = f" | 1 imaginary freq ({ts_report['imag_freqs'][0]:.0f}i cm-1) - looks like the right TS"
                else:
                    verdict = f" | {ts_report.get('n_imag')} imaginary freq(s) - check the log"
            return {
                "success": True,
                "optimized_molblocks": opt_blocks,
                "profile": profile,
                "ts_report": ts_report,
                "message": f"TS search finished. E = {e_ts:.3f} eV{verdict}",
            }
        except Exception as e:
            import traceback
            log("\nERROR: " + traceback.format_exc())
            return {"error": f"Optimization failed: {e}"}
        finally:
            log_fh.close()

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        # Save the final coordinates for the next step (IRC)
        docked_molblocks = user_input.get("docked_molblocks", [])
        return StepResult(success=True, data={
            "docked_reactants": docked_molblocks
        })

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "docked_reactants" in step_data
