"""
Step 7: Reaction Cascade Bundler & Pure Data Export
Bundles trajectory frames, SVG mechanism diagrams, and orbitals for the Molecular Workbench Viewer.
SVG processing and ZIP packaging delegated to svg_processor.py and bundle_assembler.py.
"""

import json
import os
import re
from typing import Any, Dict

from .base import PipelineStep, StepResult
from .step_6_orbitals import Step6Orbitals, enrich_timeline_tracks
from .svg_processor import normalize_svg_text
from .bundle_assembler import (
    handle_molecule_bundle_upload,
    handle_molecule_bundle_query,
    handle_bundle_export,
)
from .bundle_utils import (
    recompute_highres_batch,
    get_cached_highres_cube,
    get_or_extract_surface_mesh,
    purge_cube_files,
    save_state_to_file,
)

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, 'trajectories')
PROLONGED_PATH = os.path.join(TRAJECTORIES_DIR, 'prolonged_pathway.json')


class Step7Bundle(PipelineStep):
    """Step 7: Reaction Cascade Bundler & Pure Data Export for Molecular Workbench Viewer."""

    @property
    def step_id(self) -> str:
        return "step_7_bundle"

    @property
    def step_name(self) -> str:
        return "7. Bundle & Export"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s7 = current_state.get(self.step_id, {})
        if s7.get("megajson"):
            return s7["megajson"]
        for sk in ("step_6_orbitals", "step_5_cascade"):
            st = current_state.get(sk, {})
            if st.get("megajson"):
                return st["megajson"]
        if os.path.exists(PROLONGED_PATH):
            try:
                with open(PROLONGED_PATH, "r") as f:
                    return json.load(f)
            except Exception:
                pass
        return {}

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        megajson = self._get_megajson(current_state)
        if not megajson or not megajson.get("frames"):
            return {"type": "bundle_viewer_ui", "has_megajson": False, "message": "No reaction cascade found. Complete earlier steps first."}

        frames = megajson.get("frames", [])
        timeline = megajson.get("orbital_timeline") or current_state.get("step_6_orbitals", {}).get("orbital_timeline", {})
        if timeline and isinstance(timeline, dict) and "tracks" in timeline:
            timeline["tracks"] = enrich_timeline_tracks(megajson, timeline["tracks"])

        molecule_bundles = megajson.get("molecule_bundles", {})
        molecule_summary = {
            m_id: {"filename": data.get("filename", "bundle.zip"), "orbitals": data.get("orbitals", []), "count": len(data.get("orbitals", []))}
            for m_id, data in molecule_bundles.items() if isinstance(data, dict)
        }

        ui_megajson = dict(megajson)
        if "highres_cubes" in ui_megajson:
            ui_megajson["highres_cubes"] = {k: (v if (isinstance(v, str) and "\n" not in v) else True) for k, v in ui_megajson["highres_cubes"].items()}
        ui_megajson["molecule_bundles"] = molecule_summary

        bundle_path = os.path.join(ROOT_DIR, "out", "bundles", "reaction_cascade_bundle.zip")
        has_bundle = os.path.exists(bundle_path)

        return {
            "type": "bundle_viewer_ui", "has_megajson": True, "megajson": ui_megajson,
            "total_frames": len(frames), "step_boundaries": megajson.get("step_boundaries", []),
            "overall_barrier_kcal": megajson.get("overall_barrier_kcal", 0.0),
            "overall_reaction_energy_kcal": megajson.get("overall_reaction_energy_kcal", 0.0),
            "energy_profile": megajson.get("energy_profile", {}), "symbols": megajson.get("symbols", []),
            "trajectory_symbols": [f.get("symbols") or megajson.get("symbols", []) for f in frames],
            "trajectory_bonds": [f.get("bonds", []) for f in frames],
            "trajectory_charges": megajson.get("trajectory_charges", []), "orbital_timeline": timeline,
            "orbital_schedule": megajson.get("orbital_schedule", {}),
            "landmarks": current_state.get("step_6_orbitals", {}).get("landmarks", []),
            "mechanism_svg": normalize_svg_text(megajson.get("mechanism_svg", "")),
            "svg_wiring": megajson.get("svg_wiring", {}), "orbital_resolution": megajson.get("orbital_resolution", 32),
            "highres_count": len(megajson.get("highres_cubes", {})), "molecule_bundles": molecule_summary,
            "existing_bundle": {
                "exists": has_bundle, "filename": "reaction_cascade_bundle.zip",
                "size_bytes": os.path.getsize(bundle_path) if has_bundle else 0,
                "download_url": "/out/bundles/reaction_cascade_bundle.zip",
                "viewer_url": "/viewer.html?bundle=/out/bundles/reaction_cascade_bundle.zip"
            },
            "message": "Export lightweight, max-compression pure data archive for Molecular Workbench Viewer."
        }

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        global_state = global_state or {}
        megajson = self._get_megajson(global_state)
        if not megajson:
            return {"error": "Reaction pathway data is missing"}

        action = user_input.get("action_type")

        if action == "get_frame_charges":
            frame_idx = int(user_input.get("frame_index", 0))
            tc = megajson.get("trajectory_charges")
            if tc and 0 <= frame_idx < len(tc):
                return {"success": True, "frame_index": frame_idx, "charges": tc[frame_idx]}
            fr = megajson.get("frames", [])
            if 0 <= frame_idx < len(fr) and fr[frame_idx].get("charges"):
                return {"success": True, "frame_index": frame_idx, "charges": fr[frame_idx]["charges"]}
            return Step6Orbitals().custom_action(user_input, global_state)

        if action == "recompute_highres_cubes":
            res = recompute_highres_batch(
                megajson, user_input.get("selected_tracks", []), int(user_input.get("start_frame", 0)),
                int(user_input.get("end_frame", 15)), int(user_input.get("resolution", 50)), float(user_input.get("isovalue", 0.03))
            )
            save_state_to_file(global_state)
            return res

        if action in ["get_cube", "get_orbital_mesh"]:
            frame_idx, ok = int(user_input.get("frame_index", 0)), str(user_input.get("orbital", "HOMO"))
            iso = float(user_input.get("isovalue", 0.03))
            pos, neg = get_or_extract_surface_mesh(frame_idx, ok, iso, megajson, global_state)
            hrc = get_cached_highres_cube(frame_idx, ok, megajson)
            if pos or neg:
                resp = {"success": True, "frame_index": frame_idx, "orbital": ok, "isovalue": iso, "is_mesh": True, "pos": pos, "neg": neg, "coords": megajson["frames"][frame_idx].get("coords", megajson["frames"][frame_idx].get("positions", []))}
                if hrc: resp["cube"] = hrc
                return resp
            if hrc:
                return {"success": True, "frame_index": frame_idx, "orbital": ok, "isovalue": iso, "is_highres": True, "cube": hrc, "filename": f"{ok}_f{frame_idx}_highres.cube", "coords": megajson["frames"][frame_idx].get("coords", megajson["frames"][frame_idx].get("positions", []))}
            return Step6Orbitals().custom_action(user_input, global_state)

        if action == "set_isovalue":
            iso = float(user_input.get("isovalue", 0.03))
            megajson["orbital_isovalue"] = megajson["isovalue"] = iso
            for sk in ("step_7_bundle", "step_6_orbitals", "step_5_cascade"):
                if sk in global_state and isinstance(global_state[sk], dict) and isinstance(global_state[sk].get("megajson"), dict):
                    global_state[sk]["megajson"]["orbital_isovalue"] = global_state[sk]["megajson"]["isovalue"] = iso
            save_state_to_file(global_state)
            return {"success": True, "isovalue": iso}

        if action in ["get_cached_cubes", "get_cached_meshes"]:
            frames, tracks, iso = [int(x) for x in user_input.get("frames", [])], user_input.get("tracks", []), float(user_input.get("isovalue", 0.03))
            megajson["orbital_isovalue"] = megajson["isovalue"] = iso
            meshes = {}
            for fi in frames:
                for trk in tracks:
                    ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                    if ok != "track_esp_charge":
                        p, n = get_or_extract_surface_mesh(fi, ok, iso, megajson, global_state)
                        if p or n: meshes[f"{fi}_{ok}"] = meshes[f"{fi}_{ok}_{round(iso * 1000)}"] = {"pos": p, "neg": n}
            if meshes: return {"success": True, "meshes": meshes, "is_mesh": True}
            res = Step6Orbitals().custom_action(user_input, global_state)
            if isinstance(res, dict) and res.get("success"):
                cubes = res.setdefault("cubes", {})
                for fi in frames:
                    for trk in tracks:
                        ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                        if f"{fi}_{ok}" not in cubes and (hrc := get_cached_highres_cube(fi, ok, megajson)): cubes[f"{fi}_{ok}"] = hrc
            return res

        if action == "save_mechanism_svg":
            svg = normalize_svg_text(user_input.get("svg_text", ""))
            megajson["mechanism_svg"] = svg
            for sk in ("step_7_bundle", "step_6_orbitals", "step_5_cascade"):
                if sk in global_state and isinstance(global_state[sk], dict) and isinstance(global_state[sk].get("megajson"), dict):
                    global_state[sk]["megajson"]["mechanism_svg"] = svg
            save_state_to_file(global_state)
            if os.path.exists(PROLONGED_PATH):
                try:
                    with open(PROLONGED_PATH, "r") as pf: pj = json.load(pf)
                    pj["mechanism_svg"] = svg
                    with open(PROLONGED_PATH, "w") as pf: json.dump(pj, pf, indent=2)
                except Exception: pass
            return {"success": True, "svg_length": len(svg), "svg_text": svg}

        if action == "save_svg_wiring":
            wiring = user_input.get("wiring", {})
            megajson["svg_wiring"] = wiring
            mb = megajson.get("molecule_bundles", {})
            for grp_id, cfg in wiring.items():
                if isinstance(cfg, dict) and cfg.get("moleculeName"):
                    m_name = cfg["moleculeName"]
                    norm_name = re.sub(r'[\s\-_]+', '', m_name).lower()
                    for k, v in list(mb.items()):
                        if k == grp_id or re.sub(r'[\s\-_]+', '', str(k)).lower() == norm_name or re.sub(r'[\s\-_]+', '', str(v.get("molecule_key", ""))).lower() == norm_name:
                            v["molecule_name"] = m_name
                            mb[m_name] = mb[grp_id] = v
                            break
            for sk in ("step_7_bundle", "step_6_orbitals", "step_5_cascade"):
                if sk in global_state and isinstance(global_state[sk], dict) and isinstance(global_state[sk].get("megajson"), dict):
                    global_state[sk]["megajson"]["svg_wiring"] = wiring
            save_state_to_file(global_state)
            if os.path.exists(PROLONGED_PATH):
                try:
                    with open(PROLONGED_PATH, "r") as pf: pj = json.load(pf)
                    pj["svg_wiring"] = wiring
                    with open(PROLONGED_PATH, "w") as pf: json.dump(pj, pf, indent=2)
                except Exception: pass
            return {"success": True, "wired_count": len(wiring)}

        if action == "upload_molecule_bundle":
            return handle_molecule_bundle_upload(megajson, user_input, global_state)

        if action == "get_molecule_bundle":
            return handle_molecule_bundle_query(megajson, user_input)

        if action == "delete_molecule_bundle":
            megajson.get("molecule_bundles", {}).pop(str(user_input.get("molecule_id", "")), None)
            save_state_to_file(global_state)
            return {"success": True, "molecule_id": str(user_input.get("molecule_id", ""))}

        if action in ["repack_bundle", "export_bundle"]:
            return handle_bundle_export(megajson, user_input, global_state, ROOT_DIR)

        if action == "purge_cubes":
            purged = purge_cube_files()
            megajson.pop("highres_cubes", None)
            save_state_to_file(global_state)
            return {"success": True, "purged": purged}

        if action == "get_megajson":
            return {"success": True, "megajson": megajson, "filename": "multistep_pathway.json"}

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson:
            return StepResult(success=False, error_message="No reaction cascade available.")
        return StepResult(success=True, data={"completed": True, "megajson": megajson})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed"))
