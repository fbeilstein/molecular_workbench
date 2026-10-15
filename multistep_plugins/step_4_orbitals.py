from typing import Any, Dict, List, Optional
from rdkit import Chem
from .base import PipelineStep, StepResult
from irc_plugins.ts_engine import orbital_engine

_SCF_CACHE: Dict[int, Dict[str, Any]] = {}
_GRID_CACHE: Dict[str, Any] = {}
_CUBE_CACHE: Dict[str, str] = {}


def _get_system_charge(megajson: Dict[str, Any]) -> int:
    boundaries = megajson.get("step_boundaries", [])
    if boundaries and ">>" in boundaries[0].get("reaction_smiles", ""):
        reactants = boundaries[0]["reaction_smiles"].split(">>")[0].split(".")
        return sum(sum(a.GetFormalCharge() for a in Chem.MolFromSmiles(r).GetAtoms()) for r in reactants if Chem.MolFromSmiles(r))
    return -1


def _resolve_orbital_key(d: Dict[str, Any], orb_name: str, chem_key: str = None) -> Optional[str]:
    if not d or "vecs" not in d:
        return None
    if orb_name in d.get("vecs", {}):
        return orb_name
    all_metas = {**d.get("ibos_meta", {}), **d.get("canonical_pi_meta", {}), **d.get("antibonds_meta", {}), **d.get("meta", {})}
    for target in [chem_key, orb_name]:
        if not target:
            continue
        for k, meta in all_metas.items():
            if (meta.get("chem_key") == target or meta.get("name") == target) and k in d.get("vecs", {}):
                return k
    return None


def _extract_landmarks(megajson: Dict[str, Any]) -> List[Dict[str, Any]]:
    return []


def _build_dynamic_default_timeline(megajson: Dict[str, Any], scf_0: Dict[str, Any]) -> Dict[str, Any]:
    return {"default_isovalue": 0.03, "tracks": []}


class Step4Orbitals(PipelineStep):
    """Step 4: Interactive Molecular Orbital Studio and Schedule configuration."""

    @property
    def step_id(self) -> str:
        return "step_4_orbitals"

    @property
    def step_name(self) -> str:
        return "4. Orbitals"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s4 = current_state.get("step_4_orbitals", {})
        if s4.get("megajson"):
            return s4["megajson"]
        s3 = current_state.get("step_3_energy_stitch", {})
        return s3.get("megajson") or current_state.get("step_2_stitch_align", {}).get("megajson") or {}

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        megajson = self._get_megajson(current_state)
        if not megajson or not megajson.get("frames"):
            return {"type": "orbitals_ui", "has_megajson": False, "message": "No pathway found. Complete earlier steps first."}
        scf_0 = self._ensure_frame_computed(0, megajson)
        ep = megajson.get("energy_profile", {}).get("frame_energies_kcal", []) or [f.get("rel_energy_kcal", 0.0) for f in megajson["frames"]]
        return {
            "type": "orbitals_ui", "has_megajson": True, "total_frames": len(megajson["frames"]),
            "step_boundaries": megajson.get("step_boundaries", []), "symbols": megajson.get("symbols", []),
            "trajectory": [f["coords"] for f in megajson["frames"]],
            "landmarks": [], "energy_profile": ep,
            "timeline": megajson.get("orbital_timeline") or _build_dynamic_default_timeline(megajson, scf_0),
            "default_isovalue": 0.03, "message": "Reaction Orbital Studio: dynamic QM discovery & timeline."
        }

    def _get_grid(self, megajson: Dict[str, Any]):
        if "lo" not in _GRID_CACHE:
            all_pos = [f["coords"] for f in megajson["frames"]]
            _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"] = orbital_engine.build_global_grid(all_pos, grid_points=32, margin_angstrom=2.5)
        return _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"]

    def _ensure_frame_computed(self, frame_idx: int, megajson: Dict[str, Any]):
        if frame_idx in _SCF_CACHE:
            return _SCF_CACHE[frame_idx]
        symbols, charge = megajson["symbols"], _get_system_charge(megajson)
        if 0 not in _SCF_CACHE:
            mol0, mf0, mull0 = orbital_engine.compute_frame_scf(symbols, megajson["frames"][0]["coords"], charge=charge, spin=0)
            ibos_m0, star_m0, ibos_v0, pops0, ibo_c0 = orbital_engine.compute_frame_ibos(mol0, mf0, symbols)
            f_meta0, c_meta0, f_vecs0 = orbital_engine.extract_frontier_orbitals(mol0, mf0, None, None, ibos_m0, star_m0, ibos_v0)
            n_core0 = len([k for k, v in ibos_m0.items() if v.get("type") == "core"])
            _SCF_CACHE[0] = {"mol": mol0, "mf": mf0, "meta": f_meta0, "ibos_meta": ibos_m0, "antibonds_meta": star_m0, "canonical_pi_meta": c_meta0, "ibo_pops": pops0, "ibo_coeff": ibo_c0, "vecs": {**ibos_v0, **f_vecs0}, "charges": mull0, "n_core": n_core0}
            if frame_idx == 0:
                return _SCF_CACHE[0]

        ref_idx = min(_SCF_CACHE.keys(), key=lambda k: abs(k - frame_idx))
        ref_data = _SCF_CACHE[ref_idx]
        mol, mf, mull = orbital_engine.compute_frame_scf(symbols, megajson["frames"][frame_idx]["coords"], charge=charge, spin=0)
        n_c = ref_data.get("n_core", 11)
        ibos_m, star_m, ibos_v, pops, ibo_c = orbital_engine.compute_frame_ibos(mol, mf, symbols, {"mol": ref_data["mol"], "ibo_coeff": ref_data["ibo_coeff"], "n_core": n_c})
        f_meta, c_meta, f_vecs = orbital_engine.extract_frontier_orbitals(mol, mf, ref_data["mol"], ref_data["vecs"], ibos_m, star_m, ibos_v)
        _SCF_CACHE[frame_idx] = {"mol": mol, "mf": mf, "meta": f_meta, "ibos_meta": ibos_m, "antibonds_meta": star_m, "canonical_pi_meta": c_meta, "ibo_pops": pops, "ibo_coeff": ibo_c, "vecs": {**ibos_v, **f_vecs}, "charges": mull, "n_core": n_c}
        return _SCF_CACHE[frame_idx]

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        megajson = self._get_megajson(global_state or {})
        if not megajson:
            return {"error": "Pathway data is missing"}
        action = user_input.get("action_type")
        frame_idx = int(user_input.get("frame_index", 0))
        n_frames = len(megajson.get("frames", []))
        if frame_idx < 0 or frame_idx >= n_frames:
            return {"error": f"Frame {frame_idx} out of range"}

        if action == "get_frame_orbitals":
            data = self._ensure_frame_computed(frame_idx, megajson)
            return {"success": True, "frame_index": frame_idx, "charges": data["charges"], "frontier_orbitals": data["meta"], "all_ibos": list(data["ibos_meta"].values()), "antibonds": list(data.get("antibonds_meta", {}).values()), "canonical_pi": list(data.get("canonical_pi_meta", {}).values()), "coords": megajson["frames"][frame_idx]["coords"]}

        if action == "get_density_flux":
            sf, ef = max(0, int(user_input.get("start_frame", 0))), min(n_frames - 1, int(user_input.get("end_frame", n_frames - 1)))
            ds, de = self._ensure_frame_computed(sf, megajson), self._ensure_frame_computed(ef, megajson)
            return {"success": True, "start_frame": sf, "end_frame": ef, "fluxes": orbital_engine.compute_density_flux(ds["ibo_pops"], de["ibo_pops"], megajson["symbols"])}

        if action == "save_timeline":
            if user_input.get("timeline"):
                megajson["orbital_timeline"] = user_input["timeline"]
                return {"success": True, "message": "Timeline saved"}
            return {"error": "No timeline data"}

        if action == "batch_precompute_cubes":
            sf, sz, sel = int(user_input.get("start_frame", 0)), int(user_input.get("batch_size", 2)), user_input.get("selected_tracks", [])
            ef = min(n_frames - 1, sf + sz - 1)
            lo, step, npts = self._get_grid(megajson)
            batch_cubes = {}
            for fi in range(sf, ef + 1):
                d = self._ensure_frame_computed(fi, megajson)
                for trk in sel:
                    ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                    ck = trk.get("chem_key") if isinstance(trk, dict) else None
                    rk = _resolve_orbital_key(d, ok, ck)
                    if rk and rk in d.get("vecs", {}):
                        kr = f"{fi}_{rk}"
                        if kr not in _CUBE_CACHE:
                            v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
                            _CUBE_CACHE[kr] = orbital_engine.generate_cube_text(f"F{fi}_{rk}", megajson["symbols"], megajson["frames"][fi]["coords"], v, lo, step, npts)
                        _CUBE_CACHE[f"{fi}_{ok}"] = _CUBE_CACHE[kr]
                        batch_cubes[f"{fi}_{ok}"] = _CUBE_CACHE[kr]
                        if ck:
                            _CUBE_CACHE[f"{fi}_{ck}"] = _CUBE_CACHE[kr]
                            batch_cubes[f"{fi}_{ck}"] = _CUBE_CACHE[kr]
            return {"success": True, "start_frame": sf, "end_frame": ef, "total_frames": n_frames, "percent": int(round((ef + 1) / n_frames * 100)), "is_complete": bool(ef >= n_frames - 1), "cubes": batch_cubes}

        if action in ["get_orbital_mesh", "get_cube"]:
            ok, ck, iso = str(user_input.get("orbital", "HOMO")), user_input.get("chem_key"), float(user_input.get("isovalue", 0.03))
            for k in [f"{frame_idx}_{ok}", f"{frame_idx}_{ck}" if ck else None]:
                if k and k in _CUBE_CACHE:
                    return {"success": True, "frame_index": frame_idx, "orbital": ok, "isovalue": iso, "cube": _CUBE_CACHE[k], "filename": f"{ok}_f{frame_idx}.cube", "coords": megajson["frames"][frame_idx]["coords"]}
            d = self._ensure_frame_computed(frame_idx, megajson)
            rk = _resolve_orbital_key(d, ok, ck)
            if not rk or rk not in d.get("vecs", {}):
                return {"success": False, "available": False, "error": f"Orbital '{ok}' unavailable at frame {frame_idx}"}
            lo, step, npts = self._get_grid(megajson)
            kr = f"{frame_idx}_{rk}"
            if kr not in _CUBE_CACHE:
                v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
                _CUBE_CACHE[kr] = orbital_engine.generate_cube_text(f"F{frame_idx}_{rk}", megajson["symbols"], megajson["frames"][frame_idx]["coords"], v, lo, step, npts)
            _CUBE_CACHE[f"{frame_idx}_{ok}"] = _CUBE_CACHE[kr]
            if ck:
                _CUBE_CACHE[f"{frame_idx}_{ck}"] = _CUBE_CACHE[kr]
            return {"success": True, "frame_index": frame_idx, "orbital": ok, "isovalue": iso, "cube": _CUBE_CACHE[kr], "filename": f"{ok}_f{frame_idx}.cube", "coords": megajson["frames"][frame_idx]["coords"]}

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson:
            return StepResult(success=False, error_message="Pathway data missing.")
        tl = user_input.get("timeline") or megajson.get("orbital_timeline") or {}
        megajson["orbital_timeline"] = tl
        megajson.pop("landmark_cubes", None)
        return StepResult(success=True, data={"completed": True, "megajson": megajson, "orbital_timeline": tl, "landmarks": _extract_landmarks(megajson)})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed") and step_data.get("megajson"))
