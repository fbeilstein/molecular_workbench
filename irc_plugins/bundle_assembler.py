"""
Bundle Assembler Module
Coordinates molecule ZIP bundle uploads, lookups, and pure data archive packaging.
"""

import os
import re
from typing import Any, Dict
from .bundle_utils import (
    unpack_molecule_bundle,
    create_maxcomp_bundle_zip,
    save_state_to_file,
)


def handle_molecule_bundle_upload(
    megajson: Dict[str, Any],
    user_input: Dict[str, Any],
    global_state: Dict[str, Any]
) -> Dict[str, Any]:
    """Handle user upload and unpacking of a single-molecule bundle ZIP."""
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


def handle_molecule_bundle_query(
    megajson: Dict[str, Any],
    user_input: Dict[str, Any]
) -> Dict[str, Any]:
    """Retrieve cataloged bundle details for a requested molecule ID or name."""
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


def handle_bundle_export(
    megajson: Dict[str, Any],
    user_input: Dict[str, Any],
    global_state: Dict[str, Any],
    root_dir: str
) -> Dict[str, Any]:
    """Package the reaction cascade into maximum-compression ZIP bundle."""
    iso_input = user_input.get("isovalue")
    if iso_input is not None:
        iso_val = float(iso_input)
        megajson["orbital_isovalue"] = iso_val
        megajson["isovalue"] = iso_val
        for s_key in ("step_7_bundle", "step_6_orbitals", "step_5_cascade"):
            if s_key in global_state and isinstance(global_state[s_key], dict):
                s_mj = global_state[s_key].get("megajson")
                if isinstance(s_mj, dict):
                    s_mj["orbital_isovalue"] = iso_val
                    s_mj["isovalue"] = iso_val
        save_state_to_file(global_state)

    esp_surfaces = user_input.get("esp_surfaces")
    if esp_surfaces and isinstance(esp_surfaces, dict):
        megajson.setdefault("esp_surfaces", {}).update(esp_surfaces)
        save_state_to_file(global_state)

    bundle_dir = os.path.join(root_dir, "out", "bundles")
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
