from typing import Any, Dict
from .base import PipelineStep, StepResult
from rdkit import Chem
from rdkit.Chem import AllChem

class Step3TS(PipelineStep):
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
        # Attempt 3D embedding
        if AllChem.EmbedMolecule(mol, AllChem.ETKDGv3()) != 0:
            # Fallback if it fails
            AllChem.Compute2DCoords(mol)
        else:
            # Quick forcefield optimization to clean up geometry
            try:
                AllChem.MMFFOptimizeMolecule(mol)
            except:
                pass
        return Chem.MolToMolBlock(mol)

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        mapping_data = current_state.get("step_2_mapping")
        if not mapping_data:
            return {"message": "Please complete Step 2 (Mapping) first."}
            
        step_1 = current_state.get("step_1_input", {})
        raw_smiles = step_1.get("raw_smiles", "")
        
        if ">>" not in raw_smiles:
            return {"message": "Invalid reaction smiles."}
            
        r_part, p_part = raw_smiles.split(">>")
        r_smiles_list = r_part.split(".") if r_part else []
        p_smiles_list = p_part.split(".") if p_part else []
        
        r_blocks = [self._generate_3d_molblock(s) for s in r_smiles_list]
        p_blocks = [self._generate_3d_molblock(s) for s in p_smiles_list]
        
        return {
            "type": "ts_editor_ui",
            "reactants_3d": r_blocks,
            "products_3d": p_blocks,
            "message": "Manual 3D Docking Editor. Use the controls to arrange the reactant molecules into the pre-reaction complex."
        }

    def custom_action(self, user_input: Dict[str, Any]) -> Dict[str, Any]:
        action_type = user_input.get("action_type")
        if action_type == "optimize_ts":
            docked_molblocks = user_input.get("docked_molblocks", [])
            # TODO: Wire up Sella + ASE + PySCF/xTB here.
            # For now, we simulate an optimization by returning the exact same coordinates
            # but we could add a tiny bit of noise or run a fast MMFF94 pass.
            # Since MMFF94 doesn't support bonds breaking, returning the raw coords is safest for a mockup.
            import time
            time.sleep(1) # simulate compute time
            return {
                "success": True,
                "optimized_molblocks": docked_molblocks,
                "message": "Optimization complete (Mock: No physics engines installed yet)."
            }
        return {"error": "Unknown custom action"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        # Save the final coordinates for the next step (IRC)
        docked_molblocks = user_input.get("docked_molblocks", [])
        return StepResult(success=True, data={
            "docked_reactants": docked_molblocks
        })
        
    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "docked_reactants" in step_data
