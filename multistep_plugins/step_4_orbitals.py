import os
from typing import Any, Dict, List, Optional
import numpy as np
from rdkit import Chem
from .base import PipelineStep, StepResult
from irc_plugins.ts_engine import orbital_engine

# In-memory calculation cache across frames
_SCF_CACHE: Dict[int, Dict[str, Any]] = {}
_GRID_CACHE: Dict[str, Any] = {}
_CUBE_CACHE: Dict[str, str] = {}


def _get_system_charge(megajson: Dict[str, Any]) -> int:
    boundaries = megajson.get("step_boundaries", [])
    if boundaries:
        smi = boundaries[0].get("reaction_smiles", "")
        if ">>" in smi:
            reactants = smi.split(">>")[0].split(".")
            total_q = 0
            for r in reactants:
                m = Chem.MolFromSmiles(r)
                if m:
                    total_q += sum(a.GetFormalCharge() for a in m.GetAtoms())
            return total_q
    return -1


def _extract_landmarks(megajson: Dict[str, Any]) -> List[Dict[str, Any]]:
    frames = megajson.get("frames", [])
    n_frames = len(frames)
    boundaries = megajson.get("step_boundaries", [])
    landmarks = [{"name": "Reactants", "frame_index": 0, "desc": "Initial Reactants"}]

    for idx, b in enumerate(boundaries):
        ts_idx = b.get("ts_frame")
        if ts_idx is not None and 0 <= ts_idx < n_frames:
            landmarks.append({
                "name": f"TS {idx + 1}",
                "frame_index": ts_idx,
                "desc": f"Transition State {idx + 1}"
            })
        end_idx = b.get("end_frame")
        if end_idx is not None and 0 <= end_idx < n_frames - 1:
            landmarks.append({
                "name": f"Intermediate {idx + 1}",
                "frame_index": end_idx,
                "desc": f"Reaction Intermediate {idx + 1}"
            })

    if n_frames > 1:
        landmarks.append({"name": "Products", "frame_index": n_frames - 1, "desc": "Final Products"})
    return landmarks


class Step4Orbitals(PipelineStep):
    """Step 4: Interactive Molecular Orbital Studio and Schedule configuration."""

    @property
    def step_id(self) -> str:
        return "step_4_orbitals"

    @property
    def step_name(self) -> str:
        return "4. Orbitals"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s3 = current_state.get("step_3_energy_stitch", {})
        mj = s3.get("megajson")
        if not mj:
            mj = current_state.get("step_2_stitch_align", {}).get("megajson")
        return mj or {}

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        megajson = self._get_megajson(current_state)
        if not megajson or not megajson.get("frames"):
            return {
                "type": "orbitals_ui",
                "has_megajson": False,
                "message": "No stitched reaction pathway found. Complete earlier steps first."
            }

        landmarks = _extract_landmarks(megajson)
        existing_schedule = megajson.get("orbital_schedule") or {
            "default_isovalue": 0.03,
            "selected_orbital": "IBO_30",
            "stage_assignments": {
                "0": {"orbitals": ["IBO_30", "IBO_19", "LUMO"], "isovalue": 0.03},
                "1": {"orbitals": ["IBO_30", "IBO_20", "LUMO"], "isovalue": 0.03}
            }
        }

        return {
            "type": "orbitals_ui",
            "has_megajson": True,
            "total_frames": len(megajson.get("frames", [])),
            "step_boundaries": megajson.get("step_boundaries", []),
            "symbols": megajson.get("symbols", []),
            "landmarks": landmarks,
            "schedule": existing_schedule,
            "default_isovalue": 0.03,
            "message": "Intrinsic Reaction Orbitals (IBO Movie): visualize electron probability flow during bond formation."
        }

    def _get_grid(self, megajson: Dict[str, Any]):
        if "lo" not in _GRID_CACHE:
            all_pos = [f["coords"] for f in megajson["frames"]]
            lo, step, npts = orbital_engine.build_global_grid(all_pos, grid_points=32, margin_angstrom=2.5)
            _GRID_CACHE["lo"] = lo
            _GRID_CACHE["step"] = step
            _GRID_CACHE["npts"] = npts
        return _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"]

    def _ensure_frame_computed(self, frame_idx: int, megajson: Dict[str, Any]):
        if frame_idx in _SCF_CACHE:
            return _SCF_CACHE[frame_idx]

        symbols = megajson["symbols"]
        charge = _get_system_charge(megajson)
        coords = megajson["frames"][frame_idx]["coords"]

        prev_mol, prev_vecs = None, None
        prev_ibo_data = None
        sorted_keys = sorted(_SCF_CACHE.keys())
        if sorted_keys:
            closest = min(sorted_keys, key=lambda k: abs(k - frame_idx))
            prev_mol = _SCF_CACHE[closest]["mol"]
            prev_vecs = _SCF_CACHE[closest]["vecs"]
            prev_ibo_data = {
                "pops": _SCF_CACHE[closest]["ibo_pops"],
                "vecs": _SCF_CACHE[closest]["ibo_vecs"],
                "vec_list": list(_SCF_CACHE[closest]["ibo_vecs"].values())
            }

        mol, mf = orbital_engine.compute_frame_scf(symbols, coords, charge=charge, spin=0)
        frontier_meta, frontier_vecs = orbital_engine.extract_frontier_orbitals(mol, mf, prev_mol, prev_vecs)
        ibos_meta, ibos_vecs, ibo_pops = orbital_engine.compute_frame_ibos(mol, mf, symbols, prev_ibo_data)

        # Merge vectors for unified lookup
        all_vecs = {**ibos_vecs, **frontier_vecs}

        _SCF_CACHE[frame_idx] = {
            "mol": mol,
            "mf": mf,
            "meta": frontier_meta,
            "ibos_meta": ibos_meta,
            "ibo_pops": ibo_pops,
            "ibo_vecs": ibos_vecs,
            "vecs": all_vecs
        }
        return _SCF_CACHE[frame_idx]

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        megajson = self._get_megajson(global_state)
        if not megajson:
            return {"error": "Pathway data is missing"}

        action = user_input.get("action_type")
        frame_idx = int(user_input.get("frame_index", 0))
        n_frames = len(megajson.get("frames", []))
        if frame_idx < 0 or frame_idx >= n_frames:
            return {"error": f"Frame index {frame_idx} out of range (0..{n_frames - 1})"}

        if action == "get_frame_orbitals":
            data = self._ensure_frame_computed(frame_idx, megajson)
            reaction_flow = [
                {"key": "IBO_30", "role": "forming", "label": "⚡ Forming C-O Bond: LP(O9) → C1-O9", "badge": "Nucleophile"},
                {"key": "IBO_19", "role": "opening", "label": "⚡ Carbonyl Opening: C1=O3 (π) → LP(O3⁻)", "badge": "Electrophile"},
                {"key": "IBO_20", "role": "breaking", "label": "⚡ Leaving Group: C1-Cl2 → LP(Cl⁻)", "badge": "Leaving Group"},
                {"key": "LUMO", "role": "virtual", "label": "Acceptor Virtual Space (LUMO)", "badge": "Virtual"},
                {"key": "HOMO", "role": "frontier", "label": "Canonical HOMO", "badge": "Canonical"},
            ]
            return {
                "success": True,
                "frame_index": frame_idx,
                "reaction_flow": reaction_flow,
                "frontier_orbitals": data["meta"],
                "all_ibos": list(data["ibos_meta"].values()),
                "coords": megajson["frames"][frame_idx]["coords"]
            }

        if action == "get_orbital_mesh":
            orb_name = str(user_input.get("orbital", "IBO_30"))
            isovalue = float(user_input.get("isovalue", 0.03))
            data = self._ensure_frame_computed(frame_idx, megajson)
            if orb_name not in data["vecs"]:
                return {"error": f"Orbital '{orb_name}' not available"}

            lo, step, npts = self._get_grid(megajson)
            vol = orbital_engine.evaluate_mo_grid(data["mol"], data["vecs"][orb_name], lo, step, npts)
            mesh = orbital_engine.extract_isosurface_mesh(vol, lo, step, isovalue=isovalue)
            cube_text = orbital_engine.generate_cube_text(
                f"Frame {frame_idx} {orb_name}",
                megajson["symbols"],
                megajson["frames"][frame_idx]["coords"],
                vol, lo, step, npts
            )
            return {
                "success": True,
                "frame_index": frame_idx,
                "orbital": orb_name,
                "isovalue": isovalue,
                "mesh": mesh,
                "cube": cube_text,
                "coords": megajson["frames"][frame_idx]["coords"]
            }

        if action == "pregenerate_all_cubes":
            landmarks = _extract_landmarks(megajson)
            lo, step, npts = self._get_grid(megajson)
            precomputed = 0
            for lm in landmarks:
                f = lm["frame_index"]
                d = self._ensure_frame_computed(f, megajson)
                for o_name in ["IBO_30", "IBO_19", "IBO_20", "LUMO", "HOMO"]:
                    if o_name in d["vecs"]:
                        k = f"{f}_{o_name}"
                        if k not in _CUBE_CACHE:
                            vol = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][o_name], lo, step, npts)
                            c_txt = orbital_engine.generate_cube_text(
                                f"Frame {f} {o_name}", megajson["symbols"], megajson["frames"][f]["coords"],
                                vol, lo, step, npts
                            )
                            _CUBE_CACHE[k] = c_txt
                            precomputed += 1
            return {"success": True, "precomputed_count": precomputed, "cached_total": len(_CUBE_CACHE)}

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson:
            return StepResult(success=False, error_message="Pathway data missing.")

        schedule = user_input.get("schedule") or megajson.get("orbital_schedule") or {}
        megajson["orbital_schedule"] = schedule

        landmarks = _extract_landmarks(megajson)
        landmark_cubes = {}
        for lm in landmarks:
            f = lm["frame_index"]
            d = self._ensure_frame_computed(f, megajson)
            lo, step, npts = self._get_grid(megajson)
            for o_name in ["IBO_30", "IBO_19", "IBO_20", "LUMO", "HOMO"]:
                if o_name in d["vecs"]:
                    vol = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][o_name], lo, step, npts)
                    c_txt = orbital_engine.generate_cube_text(
                        f"Frame {f} {o_name}", megajson["symbols"], megajson["frames"][f]["coords"],
                        vol, lo, step, npts
                    )
                    landmark_cubes[f"{f}_{o_name}"] = c_txt

        megajson["landmark_cubes"] = landmark_cubes

        return StepResult(
            success=True,
            data={
                "completed": True,
                "megajson": megajson,
                "orbital_schedule": schedule,
                "landmarks": landmarks
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed") and step_data.get("megajson"))
