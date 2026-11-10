"""
Step 6: Interactive Molecular Orbital Studio
Calculates and localizes molecular orbitals (frontier MOs and IBOs) and schedules orbital tracks.
Grid meshing and chemical selection delegated to orbital_cube_mesher.py and orbital_selection.py.
"""

import json
import os
import threading
from typing import Any, Dict, List, Optional

from .base import PipelineStep, StepResult
from .ts_engine import orbital_engine
from .orbital_selection import (
    get_system_charge,
    enrich_timeline_tracks,
    resolve_orbital_key as _resolve_orbital_key,
)
from .orbital_cube_mesher import (
    clear_orbital_caches,
    get_or_compute_cube,
    batch_precompute_cubes,
    get_cached_cubes,
)
from .orbital_tracker import (
    clear_orbital_store,
    get_existing_orbital_store,
    get_or_build_orbital_store,
    get_frame_pool_payload,
)

from .constants import ROOT_DIR, CASCADE_DIR, PROLONGED_PATH
TRAJECTORIES_DIR = CASCADE_DIR

_SCF_CACHE: Dict[int, Dict[str, Any]] = {}


class Step6Orbitals(PipelineStep):
    """Step 6: Interactive Molecular Orbital Studio with distinct bonding/antibonding palettes."""

    @property
    def step_id(self) -> str: return "step_6_orbitals"

    @property
    def step_name(self) -> str: return "6. Orbitals"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s6, s5 = current_state.get(self.step_id, {}), current_state.get("step_5_cascade", {})
        mj = s6.get("megajson") or s5.get("megajson")
        if not mj and os.path.exists(PROLONGED_PATH):
            try: mj = json.load(open(PROLONGED_PATH, "r"))
            except Exception: mj = {}
        if mj and os.path.exists(PROLONGED_PATH):
            try:
                disk = json.load(open(PROLONGED_PATH, "r"))
                if disk.get("global_electron_tracks"):
                    mj["global_electron_tracks"] = disk["global_electron_tracks"]
                if disk.get("orbital_timeline"):
                    mj["orbital_timeline"] = disk["orbital_timeline"]
                if disk.get("trajectory_atom_labels"):
                    mj["trajectory_atom_labels"] = disk["trajectory_atom_labels"]
                if disk.get("step_atom_labels"):
                    mj["step_atom_labels"] = disk["step_atom_labels"]
                if disk.get("frames"):
                    mj["frames"] = disk["frames"]
                if disk.get("solvent") and "solvent" not in mj:
                    mj["solvent"] = disk["solvent"]
            except Exception: pass
        if mj and current_state.get("solvent"):
            mj["solvent"] = current_state["solvent"]
        return mj or {}

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        megajson = self._get_megajson(current_state)
        if not megajson or not megajson.get("frames"):
            return {"type": "orbitals_ui", "has_megajson": False, "message": "No reaction pathway available. Please complete Step 5 first."}

        # Check if orbital store is available to attach ESP charges and electron tracks
        try:
            store = get_existing_orbital_store(megajson)
            if store and "npz" in store:
                npz = store["npz"]
                n_f = len(megajson["frames"])
                charges_traj = [npz[f"f{f}_charges"].tolist() for f in range(n_f) if f"f{f}_charges" in npz]
                if len(charges_traj) == n_f:
                    megajson["trajectory_charges"] = charges_traj
        except Exception:
            pass

        ep = megajson.get("energy_profile", {}).get("frame_energies_kcal", []) or [f.get("rel_energy_kcal", 0.0) for f in megajson["frames"]]
        timeline = megajson.get("orbital_timeline")
        if not isinstance(timeline, dict):
            timeline = {"default_isovalue": 0.035, "tracks": []}
        else:
            tracks = timeline.get("tracks", [])
            timeline["tracks"] = enrich_timeline_tracks(megajson, tracks if isinstance(tracks, list) else [])
            timeline.setdefault("default_isovalue", 0.035)

        return {
            "type": "orbitals_ui", "has_megajson": True, "total_frames": len(megajson["frames"]),
            "step_boundaries": megajson.get("step_boundaries", []), "symbols": megajson.get("symbols", []),
            "trajectory": [f.get("coords", f.get("positions", [])) for f in megajson["frames"]],
            "trajectory_symbols": [f.get("symbols") or megajson.get("symbols", []) for f in megajson["frames"]],
            "trajectory_atom_labels": [f.get("atom_labels") for f in megajson["frames"]],
            "trajectory_bonds": [f.get("bonds", []) for f in megajson["frames"]],
            "trajectory_charges": megajson.get("trajectory_charges", []), "energy_profile": ep,
            "timeline": timeline, "default_isovalue": 0.035,
            "color_palettes": {
                "bonding": {"pos": "#38bdf8", "neg": "#1d4ed8", "label": "Bonding (Cyan/Blue)"},
                "antibonding": {"pos": "#fb923c", "neg": "#db2777", "label": "Anti-Bonding (Orange/Magenta)"},
                "core": {"pos": "#94a3b8", "neg": "#475569", "label": "Core (Gray)"}
            },
            "message": "Reaction Orbital Studio: explore frontier MOs and localized IBOs across the reaction coordinate."
        }

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        megajson = self._get_megajson(global_state or {})
        if not megajson:
            return {"error": "Pathway data is missing. Please complete Step 5."}

        action = user_input.get("action_type") or user_input.get("action")
        if action in ["clear_cache", "purge_cache"]:
            clear_orbital_store()
            clear_orbital_caches()
            return {"success": True, "message": "Orbital store and cube caches cleared."}

        if action in ["calculate_orbitals", "build_orbital_store"]:
            force = bool(user_input.get("force", False))
            store = get_or_build_orbital_store(megajson, force_rebuild=force)
            return {"success": True, "message": "Orbital tracking computed successfully.", "total_tracks": len(store.get("registry", {}).get("tracks", {}))}

        frame_idx = int(user_input.get("frame_index", user_input.get("frame_idx", 0)))
        n_frames = len(megajson.get("frames", []))
        if frame_idx < 0 or frame_idx >= n_frames:
            return {"error": f"Frame {frame_idx} out of range (0-{n_frames-1})"}

        if action == "get_frame_orbitals":
            return get_frame_pool_payload(megajson, frame_idx)

        if action == "save_timeline":
            if user_input.get("timeline"):
                tl = user_input["timeline"]
                if isinstance(tl, dict) and "tracks" in tl and isinstance(tl["tracks"], list):
                    tl["tracks"] = enrich_timeline_tracks(megajson, tl["tracks"])
                megajson["orbital_timeline"] = tl
                try:
                    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
                    with open(PROLONGED_PATH, "w") as f:
                        json.dump(megajson, f, indent=2)
                except Exception:
                    pass
                return {"success": True, "message": "Orbital timeline saved.", "timeline": tl}
            return {"error": "No timeline data provided."}

        if action == "batch_precompute_cubes":
            sf, sz, sel = int(user_input.get("start_frame", 0)), int(user_input.get("batch_size", 2)), user_input.get("selected_tracks", [])
            force = bool(user_input.get("force", False) or user_input.get("overwrite", False))
            return batch_precompute_cubes(megajson, sf, sz, sel, force)

        if action == "get_cached_cubes":
            frames, tracks = [int(x) for x in user_input.get("frames", [])], user_input.get("tracks", [])
            return {"success": True, "cubes": get_cached_cubes(frames, tracks, megajson=megajson)}

        if action in ["get_orbital_mesh", "get_cube", "get_orbital_cube"]:
            tid = user_input.get("track_id") or user_input.get("orbital", "HOMO")
            ok = str(user_input.get("orbital", tid))
            ck_raw = user_input.get("chem_key")
            ck = str(ck_raw) if (ck_raw is not None and not isinstance(ck_raw, (int, float))) else None
            try:
                iso = float(user_input.get("isovalue", 0.035) or 0.035)
            except (ValueError, TypeError):
                iso = 0.035
            return get_or_compute_cube(frame_idx, ok, ck, iso, None, megajson, track_id=tid)

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson: return StepResult(success=False, error_message="Pathway data missing.")
        tl = user_input.get("timeline") or megajson.get("orbital_timeline") or {}
        megajson["orbital_timeline"] = tl
        try:
            os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
            with open(PROLONGED_PATH, "w") as f:
                json.dump(megajson, f, indent=2)
        except Exception:
            pass
        return StepResult(success=True, data={"completed": True, "megajson": megajson, "orbital_timeline": tl})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed") and step_data.get("megajson"))
