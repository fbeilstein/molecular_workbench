import json
import os
import re
import sys
import shutil
import tempfile
import subprocess
import zipfile
from typing import Any, Dict, Optional, Tuple
from .base import PipelineStep, StepResult
from .step_6_orbitals import Step6Orbitals
from .bundle_utils import (
    recompute_highres_batch,
    get_cached_highres_cube,
    create_maxcomp_bundle_zip,
    unpack_molecule_bundle,
    purge_cube_files,
    get_or_extract_surface_mesh,
    save_state_to_file,
)

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, 'trajectories')
PROLONGED_PATH = os.path.join(TRAJECTORIES_DIR, 'prolonged_pathway.json')


def normalize_svg_text(svg_text: str) -> str:
    """Normalize SVG mechanism text to match Inkscape rendering while preserving arrow strokes and markers."""
    if not svg_text or not svg_text.strip():
        return svg_text

    if "<text" in svg_text:
        inkscape_bin = shutil.which("inkscape")
        if inkscape_bin:
            in_path = None
            out_path = None
            try:
                with tempfile.NamedTemporaryFile("w", suffix=".svg", delete=False) as in_f:
                    in_f.write(svg_text)
                    in_path = in_f.name
                out_path = in_path + ".out.svg"
                cmd = [inkscape_bin, '--actions=select-by-element:text;object-to-path', in_path, '-o', out_path]
                res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=15)
                if res.returncode == 0 and os.path.exists(out_path):
                    with open(out_path, "r") as f:
                        converted = f.read()
                    if "<svg" in converted:
                        return converted
            except Exception as e:
                print(f"[Step7] Inkscape text-to-path error: {e}")
            finally:
                for p in (in_path, out_path):
                    if p and os.path.exists(p):
                        try:
                            os.remove(p)
                        except Exception:
                            pass

        # Fallback regex sanitization if inkscape is not available: only modify <text> elements!
        def fix_style(m):
            style = m.group(1)
            style = re.sub(r'stroke:\s*#[0-9a-fA-F]+;?', 'stroke:none;', style)
            if 'font-size' not in style:
                style += ';font-size:10.5px;'
            if 'font-family' not in style:
                style += ';font-family:-apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif;'
            return f'style="{style}"'

        svg_text = re.sub(r'<text\b([^>]*?)style="([^"]*?)"', lambda m: f'<text{m.group(1)}{fix_style(m)}', svg_text)

    return svg_text


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
        s6 = current_state.get("step_6_orbitals", {})
        if s6.get("megajson"):
            return s6["megajson"]
        s5 = current_state.get("step_5_cascade", {})
        if s5.get("megajson"):
            return s5["megajson"]
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
            return {
                "type": "bundle_viewer_ui",
                "has_megajson": False,
                "message": "No reaction cascade found. Complete earlier steps first."
            }

        frames = megajson.get("frames", [])
        timeline = megajson.get("orbital_timeline") or current_state.get("step_6_orbitals", {}).get("orbital_timeline", {})
        from .step_6_orbitals import enrich_timeline_tracks
        if timeline and isinstance(timeline, dict) and "tracks" in timeline:
            timeline["tracks"] = enrich_timeline_tracks(megajson, timeline["tracks"])
        molecule_bundles = megajson.get("molecule_bundles", {})
        molecule_summary = {
            m_id: {
                "filename": data.get("filename", "bundle.zip"),
                "orbitals": data.get("orbitals", []),
                "count": len(data.get("orbitals", []))
            }
            for m_id, data in molecule_bundles.items() if isinstance(data, dict)
        }

        ui_megajson = dict(megajson)
        if "highres_cubes" in ui_megajson:
            ui_megajson["highres_cubes"] = {
                k: (v if (isinstance(v, str) and "\n" not in v) else True)
                for k, v in ui_megajson["highres_cubes"].items()
            }
        ui_megajson["molecule_bundles"] = molecule_summary

        bundle_path = os.path.join(ROOT_DIR, "out", "bundles", "reaction_cascade_bundle.zip")
        has_bundle = os.path.exists(bundle_path)
        bundle_size = os.path.getsize(bundle_path) if has_bundle else 0

        return {
            "type": "bundle_viewer_ui",
            "has_megajson": True,
            "megajson": ui_megajson,
            "total_frames": len(frames),
            "step_boundaries": megajson.get("step_boundaries", []),
            "overall_barrier_kcal": megajson.get("overall_barrier_kcal", 0.0),
            "overall_reaction_energy_kcal": megajson.get("overall_reaction_energy_kcal", 0.0),
            "energy_profile": megajson.get("energy_profile", {}),
            "symbols": megajson.get("symbols", []),
            "trajectory_charges": megajson.get("trajectory_charges", []),
            "orbital_timeline": timeline,
            "orbital_schedule": megajson.get("orbital_schedule", {}),
            "landmarks": current_state.get("step_6_orbitals", {}).get("landmarks", []),
            "mechanism_svg": normalize_svg_text(megajson.get("mechanism_svg", "")),
            "svg_wiring": megajson.get("svg_wiring", {}),
            "orbital_resolution": megajson.get("orbital_resolution", 32),
            "highres_count": len(megajson.get("highres_cubes", {})),
            "molecule_bundles": molecule_summary,
            "existing_bundle": {
                "exists": has_bundle,
                "filename": "reaction_cascade_bundle.zip",
                "size_bytes": bundle_size,
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

        # 1. High-resolution orbital recomputation
        if action == "recompute_highres_cubes":
            resolution = int(user_input.get("resolution", 50))
            tracks = user_input.get("selected_tracks", [])
            sf = int(user_input.get("start_frame", 0))
            ef = int(user_input.get("end_frame", sf + 15))
            iso = float(user_input.get("isovalue", 0.03))

            res = recompute_highres_batch(megajson, tracks, sf, ef, resolution, iso)
            save_state_to_file(global_state)
            return res

        # 2. Orbital Surface Mesh / Cube Retrieval
        if action in ["get_cube", "get_orbital_mesh"]:
            frame_idx = int(user_input.get("frame_index", 0))
            ok = str(user_input.get("orbital", "HOMO"))
            iso = float(user_input.get("isovalue", 0.03))

            pos, neg = get_or_extract_surface_mesh(frame_idx, ok, iso, megajson, global_state)
            highres_cube = get_cached_highres_cube(frame_idx, ok, megajson)
            if not highres_cube:
                try:
                    from .step_6_orbitals import _CUBE_CACHE
                    highres_cube = _CUBE_CACHE.get(f"{frame_idx}_{ok}")
                except Exception:
                    pass

            if pos or neg:
                resp = {
                    "success": True, "frame_index": frame_idx, "orbital": ok, "isovalue": iso,
                    "is_mesh": True, "pos": pos, "neg": neg,
                    "coords": megajson["frames"][frame_idx].get("coords") or megajson["frames"][frame_idx].get("positions", [])
                }
                if highres_cube:
                    resp["cube"] = highres_cube
                return resp

            if highres_cube:
                return {
                    "success": True, "frame_index": frame_idx, "orbital": ok, "isovalue": iso,
                    "is_highres": True, "cube": highres_cube,
                    "filename": f"{ok}_f{frame_idx}_highres.cube",
                    "coords": megajson["frames"][frame_idx].get("coords") or megajson["frames"][frame_idx].get("positions", [])
                }

            return Step6Orbitals().custom_action(user_input, global_state)

        if action in ["get_cached_cubes", "get_cached_meshes"]:
            frames = [int(x) for x in user_input.get("frames", [])]
            tracks = user_input.get("tracks", [])
            iso = float(user_input.get("isovalue", 0.03))
            meshes = {}
            for fi in frames:
                for trk in tracks:
                    ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                    if ok == "track_esp_charge":
                        continue
                    pos, neg = get_or_extract_surface_mesh(fi, ok, iso, megajson, global_state)
                    if pos or neg:
                        meshes[f"{fi}_{ok}"] = {"pos": pos, "neg": neg}

            if meshes:
                return {"success": True, "meshes": meshes, "is_mesh": True}

            res = Step6Orbitals().custom_action(user_input, global_state)
            if isinstance(res, dict) and res.get("success"):
                cubes = res.setdefault("cubes", {})
                for fi in frames:
                    for trk in tracks:
                        ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                        if f"{fi}_{ok}" not in cubes:
                            hrc = get_cached_highres_cube(fi, ok, megajson)
                            if hrc:
                                cubes[f"{fi}_{ok}"] = hrc
            return res

        # 3. Mechanism SVG & Wiring persistence
        if action == "save_mechanism_svg":
            svg_text = user_input.get("svg_text", "")
            svg_text = normalize_svg_text(svg_text)
            megajson["mechanism_svg"] = svg_text
            for s_key in ("step_7_bundle", "step_6_orbitals", "step_5_cascade"):
                if s_key in global_state and isinstance(global_state[s_key], dict):
                    s_mj = global_state[s_key].get("megajson")
                    if isinstance(s_mj, dict):
                        s_mj["mechanism_svg"] = svg_text
            save_state_to_file(global_state)
            if os.path.exists(PROLONGED_PATH):
                try:
                    with open(PROLONGED_PATH, "r") as pf:
                        pj = json.load(pf)
                    pj["mechanism_svg"] = svg_text
                    with open(PROLONGED_PATH, "w") as pf:
                        json.dump(pj, pf, indent=2)
                except Exception as ex:
                    print(f"[Step7] Warning persisting svg to prolonged_pathway: {ex}")
            return {"success": True, "svg_length": len(svg_text), "svg_text": svg_text}

        if action == "save_svg_wiring":
            wiring = user_input.get("wiring", {})
            megajson["svg_wiring"] = wiring
            mb = megajson.get("molecule_bundles", {})
            for grp_id, cfg in wiring.items():
                if isinstance(cfg, dict) and cfg.get("moleculeName"):
                    m_name = cfg["moleculeName"]
                    norm_name = re.sub(r'[\s\-_]+', '', m_name).lower()
                    for k, v in list(mb.items()):
                        k_norm = re.sub(r'[\s\-_]+', '', str(k)).lower()
                        vk_norm = re.sub(r'[\s\-_]+', '', str(v.get("molecule_key", ""))).lower()
                        if k == grp_id or k_norm == norm_name or vk_norm == norm_name:
                            v["molecule_name"] = m_name
                            mb[m_name] = mb[grp_id] = v
                            break
            for s_key in ("step_7_bundle", "step_6_orbitals", "step_5_cascade"):
                if s_key in global_state and isinstance(global_state[s_key], dict):
                    s_mj = global_state[s_key].get("megajson")
                    if isinstance(s_mj, dict):
                        s_mj["svg_wiring"] = wiring
            save_state_to_file(global_state)
            if os.path.exists(PROLONGED_PATH):
                try:
                    with open(PROLONGED_PATH, "r") as pf:
                        pj = json.load(pf)
                    pj["svg_wiring"] = wiring
                    with open(PROLONGED_PATH, "w") as pf:
                        json.dump(pj, pf, indent=2)
                except Exception as ex:
                    print(f"[Step7] Warning persisting wiring to prolonged_pathway: {ex}")
            return {"success": True, "wired_count": len(wiring)}

        # 4. Molecular Workbench ZIP bundle upload & orbital catalog
        if action == "upload_molecule_bundle":
            m_id = str(user_input.get("molecule_id", "molecule"))
            m_name = str(user_input.get("molecule_name", ""))
            wiring = megajson.get("svg_wiring", {})
            if not m_name and m_id in wiring:
                m_name = wiring[m_id].get("moleculeName") or wiring[m_id].get("molName", "")
            filename = str(user_input.get("filename", "bundle.zip"))
            data_b64 = user_input.get("data_b64", "")

            try:
                b_res = unpack_molecule_bundle(data_b64)
            except Exception as e:
                return {"success": False, "error": f"Failed to parse ZIP bundle: {e}"}

            bundles = megajson.setdefault("molecule_bundles", {})
            entry = {
                "filename": filename,
                "molecule_id": m_id,
                "molecule_key": b_res.get("molecule_key", m_id),
                "molecule_name": m_name or b_res.get("molecule_name", m_id),
                "orbitals": b_res.get("orbitals", []),
                "orbitals_tree": b_res.get("orbitals_tree", []),
                "xyz": b_res.get("xyz", ""),
                "esp_surface": b_res.get("esp_surface"),
                "files": b_res.get("files", {}),
                "data_b64": data_b64
            }
            bundles[m_id] = entry
            if m_name:
                bundles[m_name] = entry
            if b_res.get("molecule_key") and b_res["molecule_key"] != m_id:
                bundles[b_res["molecule_key"]] = entry

            save_state_to_file(global_state)
            return {
                "success": True,
                "molecule_id": m_id,
                "molecule_key": b_res.get("molecule_key", m_id),
                "molecule_name": m_name or b_res.get("molecule_name", m_id),
                "orbitals": b_res.get("orbitals", []),
                "orbitals_tree": b_res.get("orbitals_tree", []),
                "xyz": b_res.get("xyz", ""),
                "esp_surface": b_res.get("esp_surface"),
                "files": b_res.get("files", {}),
                "count": len(b_res.get("orbitals", []))
            }

        if action == "get_molecule_bundle":
            m_id = str(user_input.get("molecule_id", ""))
            m_name = str(user_input.get("molecule_name", ""))
            wiring = megajson.get("svg_wiring", {})
            if not m_name and m_id in wiring:
                m_name = wiring[m_id].get("moleculeName") or wiring[m_id].get("molName", "")

            def norm_k(s):
                return re.sub(r'[\s\-_]+', '', str(s or '')).lower()
            target_keys = {norm_k(m_id), norm_k(m_name)} - {""}

            bundles = megajson.get("molecule_bundles", {})
            b = bundles.get(m_id) or bundles.get(m_name)
            if not b:
                for k, v in bundles.items():
                    k_n = norm_k(k)
                    vk_n = norm_k(v.get("molecule_key", ""))
                    vn_n = norm_k(v.get("molecule_name", ""))
                    if (k_n in target_keys) or (vk_n in target_keys) or (vn_n in target_keys):
                        b = v
                        break
            if not b:
                return {"success": False, "error": f"No bundle found for molecule '{m_name or m_id}'"}
            return {
                "success": True,
                "molecule_id": m_id,
                "molecule_key": b.get("molecule_key", m_id),
                "molecule_name": m_name or b.get("molecule_name", m_id),
                "filename": b.get("filename"),
                "orbitals": b.get("orbitals", []),
                "orbitals_tree": b.get("orbitals_tree", []),
                "xyz": b.get("xyz", ""),
                "esp_surface": b.get("esp_surface"),
                "files": b.get("files", {})
            }

        if action == "delete_molecule_bundle":
            m_id = str(user_input.get("molecule_id", ""))
            bundles = megajson.get("molecule_bundles", {})
            bundles.pop(m_id, None)
            save_state_to_file(global_state)
            return {"success": True, "molecule_id": m_id}

        # 5. Pure Data Archive Packager (NO standalone HTML, NO cubes, NO xyz)
        if action in ["repack_bundle", "export_bundle"]:
            esp_surfaces = user_input.get("esp_surfaces")
            if esp_surfaces and isinstance(esp_surfaces, dict):
                megajson.setdefault("esp_surfaces", {}).update(esp_surfaces)
                save_state_to_file(global_state)

            bundle_dir = os.path.join(ROOT_DIR, "out", "bundles")
            os.makedirs(bundle_dir, exist_ok=True)
            bundle_filename = "multistep_reaction_bundle_maxcomp.zip"
            bundle_path = os.path.join(bundle_dir, bundle_filename)
            alt_path = os.path.join(bundle_dir, "reaction_cascade_bundle.zip")

            zip_bytes, uncompressed_bytes = create_maxcomp_bundle_zip(megajson, global_state)
            with open(bundle_path, "wb") as f:
                f.write(zip_bytes)
            with open(alt_path, "wb") as f:
                f.write(zip_bytes)

            ratio = (1.0 - (len(zip_bytes) / max(1, uncompressed_bytes))) * 100.0
            return {
                "success": True,
                "download_url": f"/out/bundles/{bundle_filename}",
                "viewer_url": f"/viewer.html?bundle=/out/bundles/{bundle_filename}",
                "filename": bundle_filename,
                "size_bytes": len(zip_bytes),
                "uncompressed_bytes": uncompressed_bytes,
                "compression_ratio_pct": round(ratio, 1)
            }

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

        return StepResult(
            success=True,
            data={
                "completed": True,
                "megajson": megajson
            }
        )

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed"))
