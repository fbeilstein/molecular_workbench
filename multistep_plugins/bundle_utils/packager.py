import base64
from concurrent.futures import ThreadPoolExecutor
import io
import json
import os
import re
import sys
import zipfile
from typing import Any, Dict, List, Tuple
from .highres_engine import _HIGHRES_CUBE_CACHE, CUBE_DIR

TOOLS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "tools"))
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)
from cube_converter import CubeConverter


def purge_cube_files() -> int:
    """Purge temporary .cube files from CUBE_DIR and memory caches to save space."""
    count = 0
    legacy_d = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "trajectories", "cubes"))
    for d in (CUBE_DIR, legacy_d):
        if os.path.exists(d):
            for f in os.listdir(d):
                if f.endswith(".cube"):
                    try:
                        os.remove(os.path.join(d, f))
                        count += 1
                    except OSError:
                        pass
    _HIGHRES_CUBE_CACHE.clear()
    try:
        from ..step_4_orbitals import _CUBE_CACHE
        _CUBE_CACHE.clear()
    except Exception:
        pass
    return count


def build_bundle_manifest(megajson: Dict[str, Any], global_state: Dict[str, Any], surface_keys: List[str] = None, esp_keys: List[str] = None) -> Dict[str, Any]:
    """Build standardized bundle.json manifest with surface mesh listings."""
    boundaries = megajson.get("step_boundaries", [])
    title = " → ".join(b.get("step_name", f"Step {i+1}") for i, b in enumerate(boundaries)) or "Multistep Reaction"
    rxn_smi = " -> ".join(b.get("reaction_smiles", "") for b in boundaries)
    surfaces, esps = surface_keys or [], esp_keys or []
    seen_dirs = set()
    mol_list = []
    for m in megajson.get("molecule_bundles", {}).values():
        if isinstance(m, dict) and m.get("directory") and m["directory"] not in seen_dirs:
            seen_dirs.add(m["directory"])
            mol_list.append({"id": m.get("molecule_id"), "key": m.get("molecule_key"), "name": m.get("molecule_name") or m.get("molecule_key"), "directory": m.get("directory")})
    return {
        "title": title,
        "reaction_smiles": rxn_smi,
        "total_frames": len(megajson.get("frames", [])),
        "overall_barrier_kcal": megajson.get("overall_barrier_kcal", 0.0),
        "overall_reaction_energy_kcal": megajson.get("overall_reaction_energy_kcal", 0.0),
        "solvent": megajson.get("solvent", "water"),
        "orbital_resolution": megajson.get("orbital_resolution", 32),
        "step_boundaries": boundaries,
        "symbols": megajson.get("symbols", []),
        "energy_profile": megajson.get("energy_profile", {}),
        "orbital_timeline": megajson.get("orbital_timeline", {}),
        "has_mechanism_svg": bool(megajson.get("mechanism_svg")),
        "svg_wiring": megajson.get("svg_wiring", {}),
        "has_orbitals": len(surfaces) > 0,
        "orbitals_type": "surface_mesh_json",
        "orbital_surfaces_count": len(surfaces),
        "has_esp_surfaces": len(esps) > 0,
        "esp_surfaces_count": len(esps),
        "molecules": mol_list,
        "files": {
            "trajectory_json": "multistep_pathway.json",
            "mechanism_svg": "mechanism.svg" if megajson.get("mechanism_svg") else None,
            "svg_wiring": "svg_wiring.json" if megajson.get("svg_wiring") else None,
            "orbitals_directory": "orbitals/",
            "esp_directory": "esp/" if len(esps) > 0 else None,
            "molecules_directory": "molecules/"
        }
    }


class MoleculeBundleResult(dict):
    """Dictionary supporting both dict access and tuple unpacking for backward compatibility."""
    def __iter__(self):
        return iter([self.get("orbitals", []), self.get("files", {})])


def unpack_molecule_bundle(data_b64: str) -> MoleculeBundleResult:
    """Unpack a base64-encoded Molecular Workbench zip bundle and extract orbital meshes/esp."""
    raw_zip = base64.b64decode(data_b64)
    orbitals_found: List[str] = []
    files_dict: Dict[str, Any] = {}
    meta = {"xyz": "", "orbitals_tree": [], "esp_surface": None, "molecule_key": "", "smiles": ""}

    with zipfile.ZipFile(io.BytesIO(raw_zip), "r") as zf:
        for name in zf.namelist():
            try:
                content = zf.read(name).decode("utf-8", errors="ignore")
                files_dict[name] = files_dict[os.path.basename(name)] = content
                if name.endswith(".xyz") and not meta["xyz"]:
                    meta["xyz"] = content
                elif name.endswith(".cube") or name.endswith(("_pos.json", "_neg.json")):
                    base = os.path.basename(name).replace(".cube", "").replace("_pos.json", "").replace("_neg.json", "")
                    if base and base not in orbitals_found:
                        orbitals_found.append(base)
            except Exception:
                continue

        mol_data = None
        for name, content in files_dict.items():
            if (name.startswith("molecules/") and name.endswith(".json") and not name.endswith(("_esp.json", "_pos.json", "_neg.json"))) or name == "manifest.json":
                try:
                    d = json.loads(content)
                    mol_data = d.get("molecules", [d])[0] if (isinstance(d, dict) and "molecules" in d) else d
                    break
                except Exception:
                    pass

        if isinstance(mol_data, dict):
            meta["molecule_key"] = mol_data.get("key", mol_data.get("name", ""))
            meta["smiles"] = mol_data.get("smiles", "")
            meta["xyz"] = meta["xyz"] or mol_data.get("xyz", "")
            meta["esp_surface"] = mol_data.get("esp_surface")
            if mol_data.get("orbitals"):
                meta["orbitals_tree"] = mol_data["orbitals"]
                for it in [it for cat in meta["orbitals_tree"] for it in cat.get("items", [])]:
                    b = os.path.basename(it.get("file", "")).replace(".json", "").replace(".cube", "")
                    if b and b not in orbitals_found:
                        orbitals_found.append(b)

        if not meta["esp_surface"]:
            for fname in files_dict:
                if fname.endswith(("_esp.json", "/esp.json")) or fname == "esp.json":
                    meta["esp_surface"] = {"file": fname}
                    break

        mol_key = meta["molecule_key"]
        if mol_key:
            ws_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
            local_full_zip = os.path.join(ws_root, "output", mol_key, f"{mol_key}.rxnbundle.zip")
            if os.path.exists(local_full_zip):
                try:
                    with zipfile.ZipFile(local_full_zip, "r") as lzf:
                        for lname in lzf.namelist():
                            if lname not in files_dict and lname.endswith((".json", ".xyz")):
                                lc = lzf.read(lname).decode("utf-8", errors="ignore")
                                files_dict[lname] = files_dict[os.path.basename(lname)] = lc
                                if lname.endswith(("_pos.json", "_neg.json")):
                                    orb = os.path.basename(lname).replace("_pos.json", "").replace("_neg.json", "")
                                    if orb not in orbitals_found:
                                        orbitals_found.append(orb)
                except Exception:
                    pass

    return MoleculeBundleResult({
        "orbitals": orbitals_found, "files": files_dict, "xyz": meta["xyz"],
        "orbitals_tree": meta["orbitals_tree"], "esp_surface": meta["esp_surface"],
        "molecule_key": meta["molecule_key"], "smiles": meta["smiles"]
    })



from .orbital_packager import convert_and_pack_trajectory_orbitals


def _convert_and_pack_orbitals(megajson: Dict[str, Any], zf: zipfile.ZipFile, isovalue: float = 0.035, global_state: Dict[str, Any] = None) -> Tuple[List[str], int]:
    """Calculate lightweight 3D surface meshes across all frames, pack them into zip, and never pack .cube files."""
    return convert_and_pack_trajectory_orbitals(megajson, zf, isovalue, global_state)


def create_maxcomp_bundle_zip(megajson: Dict[str, Any], global_state: Dict[str, Any]) -> Tuple[bytes, int]:
    """Calculate surfaces and pack into zip with max compression (level 9). Never pack cubes, xyz, or index.html."""
    zip_buf = io.BytesIO()
    isovalue = float(megajson.get("orbital_isovalue") or megajson.get("isovalue") or 0.035)

    with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        # 1. Convert cubes to surfaces and pack into orbitals/
        surface_keys, orb_bytes = _convert_and_pack_orbitals(megajson, zf, isovalue, global_state)


        # 2. Pack ESP surfaces into esp/ with 4-decimal vertex precision
        esp_surfaces = megajson.get("esp_surfaces", {})
        esp_keys, esp_bytes = [], 0
        for f_key, mesh in esp_surfaces.items():
            if not mesh or not isinstance(mesh, dict):
                continue
            clean_f = str(f_key).replace("frame_", "")
            clean_mesh = dict(mesh)
            if "vertexArr" in clean_mesh and isinstance(clean_mesh["vertexArr"], list):
                clean_mesh["vertexArr"] = [
                    {'x': round(float(v.get('x', 0)), 4), 'y': round(float(v.get('y', 0)), 4), 'z': round(float(v.get('z', 0)), 4)}
                    if isinstance(v, dict) else [round(float(c), 4) for c in v]
                    for v in clean_mesh["vertexArr"]
                ]
            if "color" in clean_mesh and isinstance(clean_mesh["color"], list):
                clean_mesh["color"] = [
                    {'r': round(float(c.get('r', 0)), 2), 'g': round(float(c.get('g', 0)), 2), 'b': round(float(c.get('b', 0)), 2)}
                    if isinstance(c, dict) else round(float(c), 2)
                    for c in clean_mesh["color"]
                ]
            mesh_json = json.dumps(clean_mesh, separators=(',', ':')).encode('utf-8')
            zf.writestr(f"esp/frame_{clean_f}.json", mesh_json)
            esp_bytes += len(mesh_json)
            esp_keys.append(clean_f)

        # 3. Store uploaded molecule bundles as UNPACKED data (no nested zip archives, no cubes, no xyz)
        molecule_bundles = megajson.get("molecule_bundles") or (global_state.get("step_4_orbitals", {}).get("megajson", {}).get("molecule_bundles") if global_state else None) or {}
        clean_mol_bundles = {}
        mol_bytes = 0
        written_files = set()
        wiring = megajson.get("svg_wiring", {})

        def norm_name(s):
            return re.sub(r'[\s\-_]+', '', str(s or '')).lower()

        def resolve_molecule_name(m_id, mol_key, b_data):
            if isinstance(b_data, dict) and b_data.get("molecule_name"):
                return b_data["molecule_name"]
            target_ids = {m_id, mol_key, norm_name(m_id), norm_name(mol_key)}
            for grp_id, w_cfg in wiring.items():
                if isinstance(w_cfg, dict) and w_cfg.get("moleculeName"):
                    wn = w_cfg["moleculeName"]
                    if grp_id in target_ids or norm_name(wn) in target_ids or norm_name(grp_id) in target_ids:
                        return wn
            return (b_data.get("molecule_name") or b_data.get("title") or mol_key.replace('_', ' ').title()) if isinstance(b_data, dict) else mol_key

        for m_id, b_data in molecule_bundles.items():
            mol_key = b_data.get("molecule_key") or m_id
            mol_dir = f"molecules/{mol_key}"

            if "data_b64" in b_data:
                try:
                    with zipfile.ZipFile(io.BytesIO(base64.b64decode(b_data["data_b64"])), "r") as zin:
                        for it in zin.infolist():
                            if it.is_dir() or it.filename.endswith((".cube", ".xyz")):
                                continue
                            sub = it.filename
                            if sub.startswith("molecules/"):
                                sub = sub[10:]
                            if sub.startswith(f"{mol_key}/"):
                                sub = sub[len(mol_key) + 1:]
                            dest = f"{mol_dir}/{sub}"
                            if dest not in written_files:
                                content = zin.read(it.filename)
                                zf.writestr(dest, content)
                                mol_bytes += len(content)
                                written_files.add(dest)
                except Exception as e:
                    print(f"[Packager] Warning unpacking molecule {m_id}: {e}")

            for fname, fcontent in b_data.get("files", {}).items():
                if fname.endswith((".cube", ".xyz")):
                    continue
                dest = f"{mol_dir}/{os.path.basename(fname)}"
                if dest not in written_files:
                    c_bytes = fcontent.encode("utf-8") if isinstance(fcontent, str) else fcontent
                    zf.writestr(dest, c_bytes)
                    mol_bytes += len(c_bytes)
                    written_files.add(dest)

            mol_name = resolve_molecule_name(m_id, mol_key, b_data)
            entry = {
                "molecule_id": m_id, "molecule_key": mol_key, "molecule_name": mol_name, "directory": f"{mol_dir}/",
                "orbitals": b_data.get("orbitals", []), "orbitals_tree": b_data.get("orbitals_tree", []),
                "xyz": b_data.get("xyz", ""), "esp_surface": b_data.get("esp_surface"), "file_count": len(written_files)
            }
            clean_mol_bundles[m_id] = entry
            if mol_name:
                clean_mol_bundles[mol_name] = entry
                clean_mol_bundles[norm_name(mol_name)] = entry

        # Fallback to copy pre-packed molecules if molecule_bundles had no unpacked files
        if not written_files:
            ws_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
            for ref_zip in (os.path.join(ws_root, "bundles", "tetrahedral_intermediate_CO.zip"), os.path.join(ws_root, "out", "tetrahedral_intermediate_CO.zip")):
                if os.path.exists(ref_zip):
                    try:
                        with zipfile.ZipFile(ref_zip, "r") as rzf:
                            for r_name in rzf.namelist():
                                if r_name.startswith("molecules/") and not rzf.getinfo(r_name).is_dir():
                                    if r_name not in written_files:
                                        r_bytes = rzf.read(r_name)
                                        zf.writestr(r_name, r_bytes)
                                        mol_bytes += len(r_bytes)
                                        written_files.add(r_name)
                            if "multistep_pathway.json" in rzf.namelist():
                                ref_pw = json.loads(rzf.read("multistep_pathway.json").decode("utf-8"))
                                if ref_pw.get("molecule_bundles"):
                                    for rk, rv in ref_pw["molecule_bundles"].items():
                                        if isinstance(rv, dict):
                                            r_entry = dict(rv)
                                            r_name = resolve_molecule_name(rk, r_entry.get("molecule_key", rk), r_entry)
                                            r_entry["molecule_name"] = r_name
                                            clean_mol_bundles[rk] = r_entry
                                            if r_name:
                                                clean_mol_bundles[r_name] = r_entry
                                                clean_mol_bundles[norm_name(r_name)] = r_entry
                    except Exception:
                        pass
                if written_files:
                    break

        # 4. Canonicalize & clean megajson (strip dead keys, remove positions, round precision)
        clean_megajson = dict(megajson)
        for dead_k in ("combined_frames", "steps", "landmark_orbital_meshes", "stitching_alignment_report", "junction_previews", "highres_cubes", "molecule_bundles", "esp_surfaces"):
            clean_megajson.pop(dead_k, None)

        clean_megajson["orbital_surfaces"] = sorted(surface_keys)
        clean_megajson["esp_frame_keys"] = sorted(esp_keys)
        clean_megajson["molecule_bundles"] = clean_mol_bundles
        megajson["orbital_surfaces"] = sorted(surface_keys)

        clean_frames = []
        for f in clean_megajson.get("frames", []):
            rc = f.get("coords", f.get("positions", []))
            e = f.get("rel_energy_kcal") if f.get("rel_energy_kcal") is not None else f.get("energy_kcal", 0.0)
            cf = {
                "step_index": f.get("step_index", 1),
                "coords": [[round(float(c), 4) for c in pt] for pt in rc],
                "is_ts": bool(f.get("is_ts", False)),
                "rel_energy_kcal": round(float(e or 0.0), 2)
            }
            if "energy_kcal" in f and f["energy_kcal"] is not None:
                cf["energy_kcal"] = round(float(f["energy_kcal"]), 2)
            if "step" in f and f["step"] is not None:
                cf["step"] = f["step"]
            if "s_arc" in f and f["s_arc"] is not None:
                cf["s_arc"] = round(float(f["s_arc"]), 2)
            if f.get("charges"):
                cf["charges"] = [round(float(q), 3) for q in f["charges"]]
            clean_frames.append(cf)
        clean_megajson["frames"] = clean_frames

        if "trajectory_charges" in clean_megajson:
            clean_megajson["trajectory_charges"] = [[round(float(q), 3) for q in fq] for fq in clean_megajson["trajectory_charges"]]

        if isinstance(clean_megajson.get("energy_profile"), dict) and "frame_energies_kcal" in clean_megajson["energy_profile"]:
            clean_megajson["energy_profile"]["frame_energies_kcal"] = [round(float(e), 2) for e in clean_megajson["energy_profile"]["frame_energies_kcal"]]

        for k in ("overall_barrier_kcal", "overall_reaction_energy_kcal"):
            if clean_megajson.get(k) is not None:
                clean_megajson[k] = round(float(clean_megajson[k]), 2)

        # 5. Build manifest and metadata (no XYZ, no index.html)
        manifest = build_bundle_manifest(clean_megajson, global_state, surface_keys, esp_keys)
        json_content = json.dumps(clean_megajson, indent=2)

        uncompressed_bytes = orb_bytes + esp_bytes + mol_bytes + len(json_content.encode("utf-8"))

        zf.writestr("bundle.json", json.dumps(manifest, indent=2))
        zf.writestr("multistep_pathway.json", json_content)

        svg_content = clean_megajson.get("mechanism_svg")
        if svg_content:
            zf.writestr("mechanism.svg", svg_content)
            uncompressed_bytes += len(svg_content.encode("utf-8"))

        wiring = clean_megajson.get("svg_wiring")
        if wiring:
            wiring_str = json.dumps(wiring, indent=2)
            zf.writestr("svg_wiring.json", wiring_str)
            uncompressed_bytes += len(wiring_str.encode("utf-8"))

    # 6. Purge temporary cubes to save space!
    purged_count = purge_cube_files()
    print(f"[Packager] Purged {purged_count} temporary cube files after packaging.")

    zip_buf.seek(0)
    return zip_buf.getvalue(), uncompressed_bytes
