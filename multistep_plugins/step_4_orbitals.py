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
_TRACK_MAP: Dict[str, int] = {}


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
        existing_timeline = megajson.get("orbital_timeline") or {
            "default_isovalue": 0.03,
            "tracks": [
                {
                    "id": "trk_ibo23",
                    "orbital_key": "IBO_23",
                    "title": "C1-O9 σ (Forming)",
                    "color": "#06b6d4",
                    "default_on": True,
                    "clips": [{"id": "c1", "start": 20, "end": 68, "title": "Nu Addition"}]
                },
                {
                    "id": "trk_ibo21",
                    "orbital_key": "IBO_21",
                    "title": "C1=O3 π (Opening)",
                    "color": "#a855f7",
                    "default_on": True,
                    "clips": [{"id": "c2", "start": 0, "end": 68, "title": "π → Oxyanion"}]
                },
                {
                    "id": "trk_lumo",
                    "orbital_key": "LUMO",
                    "title": "LUMO (Carbonyl π*)",
                    "color": "#f59e0b",
                    "default_on": False,
                    "clips": [{"id": "c3", "start": 0, "end": 45, "title": "Acceptor Lobe"}]
                },
                {
                    "id": "trk_ibo20",
                    "orbital_key": "IBO_20",
                    "title": "C1-Cl2 σ (Cleaving)",
                    "color": "#10b981",
                    "default_on": True,
                    "clips": [{"id": "c4", "start": 68, "end": 146, "title": "Cl⁻ Elimination"}]
                }
            ]
        }

        ep_data = megajson.get("energy_profile", {}).get("frame_energies_kcal", [])
        if not ep_data and megajson.get("frames"):
            ep_data = [f.get("rel_energy_kcal", 0.0) for f in megajson["frames"]]

        return {
            "type": "orbitals_ui",
            "has_megajson": True,
            "total_frames": len(megajson.get("frames", [])),
            "step_boundaries": megajson.get("step_boundaries", []),
            "symbols": megajson.get("symbols", []),
            "landmarks": landmarks,
            "energy_profile": ep_data,
            "timeline": existing_timeline,
            "default_isovalue": 0.03,
            "message": "Reaction Orbital Studio (Kdenlive Timeline): audition, clip, and direct electron flow across frames."
        }

    def _get_grid(self, megajson: Dict[str, Any]):
        if "lo" not in _GRID_CACHE:
            all_pos = [f["coords"] for f in megajson["frames"]]
            lo, step, npts = orbital_engine.build_global_grid(all_pos, grid_points=50, margin_angstrom=2.5)
            _GRID_CACHE["lo"] = lo
            _GRID_CACHE["step"] = step
            _GRID_CACHE["npts"] = npts
        return _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"]

    def _ensure_frame_computed(self, frame_idx: int, megajson: Dict[str, Any]):
        if frame_idx in _SCF_CACHE:
            return _SCF_CACHE[frame_idx]

        symbols = megajson["symbols"]
        charge = _get_system_charge(megajson)
        if 0 not in _SCF_CACHE:
            coords0 = megajson["frames"][0]["coords"]
            mol0, mf0, mull0 = orbital_engine.compute_frame_scf(symbols, coords0, charge=charge, spin=0)
            ibos_m0, star_m0, ibos_v0, pops0, ibo_c0 = orbital_engine.compute_frame_ibos(mol0, mf0, symbols)
            f_meta0, c_meta0, f_vecs0 = orbital_engine.extract_frontier_orbitals(
                mol0, mf0, None, None, ibos_m0, star_m0, ibos_v0
            )
            _SCF_CACHE[0] = {
                "mol": mol0, "mf": mf0, "meta": f_meta0, "ibos_meta": ibos_m0,
                "antibonds_meta": star_m0, "canonical_pi_meta": c_meta0,
                "ibo_pops": pops0, "ibo_coeff": ibo_c0, "vecs": {**ibos_v0, **f_vecs0},
                "charges": mull0
            }
            if frame_idx == 0:
                return _SCF_CACHE[0]

        # Use nearest computed frame as reference for phase / orbital continuity
        cur = max([k for k in sorted(_SCF_CACHE.keys()) if k <= frame_idx], default=0)
        coords = megajson["frames"][frame_idx]["coords"]
        mol, mf, mull = orbital_engine.compute_frame_scf(symbols, coords, charge=charge, spin=0)
        prev_d = {"mol": _SCF_CACHE[cur]["mol"], "ibo_coeff": _SCF_CACHE[cur]["ibo_coeff"]}
        ibos_m, star_m, ibos_v, pops, ibo_c = orbital_engine.compute_frame_ibos(mol, mf, symbols, prev_d)
        f_meta, c_meta, f_vecs = orbital_engine.extract_frontier_orbitals(
            mol, mf, _SCF_CACHE[cur]["mol"], _SCF_CACHE[cur]["vecs"], ibos_m, star_m, ibos_v
        )
        _SCF_CACHE[frame_idx] = {
            "mol": mol, "mf": mf, "meta": f_meta, "ibos_meta": ibos_m,
            "antibonds_meta": star_m, "canonical_pi_meta": c_meta,
            "ibo_pops": pops, "ibo_coeff": ibo_c, "vecs": {**ibos_v, **f_vecs},
            "charges": mull
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
            return {
                "success": True,
                "frame_index": frame_idx,
                "charges": data["charges"],
                "frontier_orbitals": data["meta"],
                "all_ibos": list(data["ibos_meta"].values()),
                "antibonds": list(data.get("antibonds_meta", {}).values()),
                "canonical_pi": list(data.get("canonical_pi_meta", {}).values()),
                "coords": megajson["frames"][frame_idx]["coords"]
            }

        if action == "get_density_flux":
            start_f = max(0, min(n_frames - 1, int(user_input.get("start_frame", 0))))
            end_f = max(0, min(n_frames - 1, int(user_input.get("end_frame", n_frames - 1))))
            d_start = self._ensure_frame_computed(start_f, megajson)
            d_end = self._ensure_frame_computed(end_f, megajson)
            fluxes = orbital_engine.compute_density_flux(d_start["ibo_pops"], d_end["ibo_pops"], megajson["symbols"])
            return {
                "success": True,
                "start_frame": start_f,
                "end_frame": end_f,
                "fluxes": fluxes
            }

        if action == "save_timeline":
            timeline = user_input.get("timeline")
            if timeline:
                megajson["orbital_timeline"] = timeline
                return {"success": True, "message": "Timeline persisted successfully"}
            return {"error": "No timeline data provided"}

        if action == "batch_precompute_cubes":
            start_f = int(user_input.get("start_frame", 0))
            batch_sz = int(user_input.get("batch_size", 6))
            selected_tracks = user_input.get("selected_tracks", ["IBO_23", "IBO_21", "IBO_20", "LUMO"])
            end_f = min(n_frames - 1, start_f + batch_sz - 1)
            lo, step, npts = self._get_grid(megajson)
            batch_cubes = {}
            for fi in range(start_f, end_f + 1):
                d = self._ensure_frame_computed(fi, megajson)
                for trk in selected_tracks:
                    if trk in d["vecs"]:
                        k = f"{fi}_{trk}"
                        if k not in _CUBE_CACHE:
                            vol = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][trk], lo, step, npts)
                            _CUBE_CACHE[k] = orbital_engine.generate_cube_text(
                                f"Frame {fi} {trk}", megajson["symbols"], megajson["frames"][fi]["coords"], vol, lo, step, npts
                            )
                        batch_cubes[k] = _CUBE_CACHE[k]
            return {
                "success": True, "start_frame": start_f, "end_frame": end_f, "total_frames": n_frames,
                "percent": int(round((end_f + 1) / n_frames * 100)), "cubes": batch_cubes,
                "is_complete": bool(end_f >= n_frames - 1)
            }

        if action in ["get_orbital_mesh", "get_cube"]:
            orb_name = str(user_input.get("orbital", "IBO_23"))
            chem_key = user_input.get("chem_key")
            isovalue = float(user_input.get("isovalue", 0.03))
            d = self._ensure_frame_computed(frame_idx, megajson)

            resolved_key = None
            if orb_name in d["vecs"]:
                resolved_key = orb_name
            else:
                search_targets = [chem_key, orb_name]
                all_metas = {
                    **d.get("ibos_meta", {}),
                    **d.get("canonical_pi_meta", {}),
                    **d.get("antibonds_meta", {}),
                    **d.get("meta", {})
                }
                for target in search_targets:
                    if not target:
                        continue
                    for k, meta in all_metas.items():
                        if meta.get("chem_key") == target or meta.get("name") == target:
                            resolved_key = k
                            break
                    if resolved_key:
                        break

            if not resolved_key or resolved_key not in d["vecs"]:
                return {"error": f"Orbital '{orb_name}' not available at frame {frame_idx}"}

            orb_name = resolved_key

            lo, step, npts = self._get_grid(megajson)
            k = f"{frame_idx}_{orb_name}"
            if k not in _CUBE_CACHE:
                vol = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][orb_name], lo, step, npts)
                _CUBE_CACHE[k] = orbital_engine.generate_cube_text(
                    f"Frame {frame_idx} {orb_name}", megajson["symbols"], megajson["frames"][frame_idx]["coords"], vol, lo, step, npts
                )
            cube_text = _CUBE_CACHE[k]
            return {
                "success": True, "frame_index": frame_idx, "orbital": orb_name,
                "isovalue": isovalue, "cube": cube_text, "filename": f"{orb_name}_frame{frame_idx}.cube",
                "coords": megajson["frames"][frame_idx]["coords"]
            }

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson:
            return StepResult(success=False, error_message="Pathway data missing.")

        timeline = user_input.get("timeline") or megajson.get("orbital_timeline") or {}
        megajson["orbital_timeline"] = timeline

        landmarks = _extract_landmarks(megajson)
        landmark_cubes = {}
        track_keys = [t.get("orbital_key") for t in timeline.get("tracks", []) if t.get("orbital_key")]
        if not track_keys:
            track_keys = ["IBO_23", "IBO_21", "IBO_20", "LUMO", "HOMO"]

        for lm in landmarks:
            f = lm["frame_index"]
            d = self._ensure_frame_computed(f, megajson)
            lo, step, npts = self._get_grid(megajson)
            for o_name in track_keys:
                if o_name in d["vecs"]:
                    vol = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][o_name], lo, step, npts)
                    c_txt = orbital_engine.generate_cube_text(
                        f"Frame {f} {o_name}", megajson["symbols"], megajson["frames"][f]["coords"],
                        vol, lo, step, npts
                    )
                    landmark_cubes[f"{f}_{o_name}"] = c_txt

        megajson["landmark_cubes"] = {**_CUBE_CACHE, **landmark_cubes}

        return StepResult(
            success=True,
            data={
                "completed": True,
                "megajson": megajson,
                "orbital_timeline": timeline,
                "landmarks": landmarks
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed") and step_data.get("megajson"))
