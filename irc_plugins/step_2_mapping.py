from typing import Any, Dict
from .base import PipelineStep, StepResult
from .rxn_mapping import parse_reaction_smiles, auto_map_atoms, _validate_pt_mapping, _get_rdkit

class Step2Mapping(PipelineStep):
    @property
    def step_id(self) -> str:
        return "step_2_mapping"
        
    @property
    def step_name(self) -> str:
        return "2. Atom Mapping"

    def _generate_svg_and_coords(self, mols):
        Chem, _, _, _ = _get_rdkit()
        from rdkit.Chem.Draw import rdMolDraw2D
        from rdkit.Chem import rdDepictor
        
        svgs = []
        coords = []
        
        for mol in mols:
            if not mol:
                svgs.append("")
                coords.append({})
                continue
                
            # Compute 2D coordinates for the drawing
            rdDepictor.Compute2DCoords(mol)
            
            # Use fixed dimensions
            width = 400
            height = 300
            d2d = rdMolDraw2D.MolDraw2DSVG(width, height)
            
            opts = d2d.drawOptions()
            opts.clearBackground = False # transparent
            
            d2d.DrawMolecule(mol)
            d2d.FinishDrawing()
            
            frag_coords = {}
            for i in range(mol.GetNumAtoms()):
                pt = d2d.GetDrawCoords(i)
                frag_coords[i] = {"x": pt.x, "y": pt.y, "symbol": mol.GetAtomWithIdx(i).GetSymbol()}
                
            # Add a viewBox wrapper div or just return raw SVG
            svgs.append(d2d.GetDrawingText())
            coords.append(frag_coords)
            
        return svgs, coords

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        if "step_1_input" not in current_state:
            return {"message": "Please complete Step 1 first."}
            
        smiles = current_state["step_1_input"]["raw_smiles"]
        
        try:
            r_smi, p_smi = parse_reaction_smiles(smiles)
        except Exception as e:
            return {"message": f"Error parsing reaction: {e}"}
            
        Chem, _, _, _ = _get_rdkit()
        r_mols = [Chem.MolFromSmiles(s) for s in r_smi]
        p_mols = [Chem.MolFromSmiles(s) for s in p_smi]
        
        # Auto-map without explicit Hs to avoid symmetrical explosion
        mapping, _ = auto_map_atoms(r_smi, p_smi)
        mapping = _validate_pt_mapping(mapping, r_smi, p_smi)
        
        # Now add explicit Hs. RDKit preserves existing atom indices.
        r_mols_h = [Chem.AddHs(m) if m else None for m in r_mols]
        p_mols_h = [Chem.AddHs(m) if m else None for m in p_mols]
        
        # Extend mapping to the newly added attached hydrogens
        # Convert mapping to list of keys so we can iterate safely
        heavy_keys = list(mapping.keys())
        for r_key in heavy_keys:
            p_key = mapping[r_key]
            r_fi, r_ai = r_key[1], r_key[2]
            p_fi, p_ai = p_key[1], p_key[2]
            
            if not r_mols_h[r_fi] or not p_mols_h[p_fi]: continue
            
            r_atom = r_mols_h[r_fi].GetAtomWithIdx(r_ai)
            r_hs = [n.GetIdx() for n in r_atom.GetNeighbors() if n.GetAtomicNum() == 1]
            
            p_atom = p_mols_h[p_fi].GetAtomWithIdx(p_ai)
            p_hs = [n.GetIdx() for n in p_atom.GetNeighbors() if n.GetAtomicNum() == 1]
            
            # Map up to the minimum number of Hs on both sides
            for rh, ph in zip(r_hs, p_hs):
                mapping[('R', r_fi, rh)] = ('P', p_fi, ph)
        
        r_svgs, r_coords = self._generate_svg_and_coords(r_mols_h)
        p_svgs, p_coords = self._generate_svg_and_coords(p_mols_h)
        
        # Convert tuple keys to strings for JSON
        json_mapping = {}
        for k, v in mapping.items():
            k_str = f"{k[0]}_{k[1]}_{k[2]}" # e.g. "R_0_4"
            v_str = f"{v[0]}_{v[1]}_{v[2]}" # e.g. "P_0_5"
            json_mapping[k_str] = v_str
            
        return {
            "type": "mapping_ui",
            "reactants": {
                "smiles": r_smi,
                "svgs": r_svgs,
                "coords": r_coords
            },
            "products": {
                "smiles": p_smi,
                "svgs": p_svgs,
                "coords": p_coords
            },
            "initial_mapping": json_mapping,
            "message": "Verify and correct atom mappings. Click a wire to delete it. Click two atoms to create a new wire."
        }

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        mapping = user_input.get("mapping", {})
        
        # Keep keys as strings (e.g. "R_0_1") because JSON requires string keys
        return StepResult(success=True, data={
            "final_mapping": mapping
        })
        
    def validate(self, step_data: Dict[str, Any]) -> bool:
        return "final_mapping" in step_data
