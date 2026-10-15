import base64
from concurrent.futures import ThreadPoolExecutor
import io
import json
import os
import sys
import zipfile
from typing import Any, Dict, List, Tuple
from .standalone_template import generate_standalone_html
from .highres_engine import _HIGHRES_CUBE_CACHE, CUBE_DIR

TOOLS_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "tools"))
if TOOLS_DIR not in sys.path:
    sys.path.insert(0, TOOLS_DIR)
from cube_converter import CubeConverter


def generate_xyz(megajson: Dict[str, Any]) -> str:
    """Generate multi-model XYZ text format from pathway frames."""
    symbols = megajson.get("symbols", [])
    n_atoms = len(symbols)
    lines = []
    for f_idx, frame in enumerate(megajson.get("frames", [])):
        coords = frame.get("coords", frame.get("positions", []))
        e = frame.get("rel_energy_kcal") or frame.get("energy_kcal") or 0.0
        lines.append(f"{n_atoms}\nFrame {f_idx} | Step {frame.get('step_index', 0)} | Relative Free Energy: {e:.2f} kcal/mol")
        for sym, pos in zip(symbols, coords):
            lines.append(f"{sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}")
    return "\n".join(lines) + "\n"


def purge_cube_files() -> int:
    """Purge temporary .cube files from CUBE_DIR and memory caches to save space."""
    count = 0
    if os.path.exists(CUBE_DIR):
        for f in os.listdir(CUBE_DIR):
            if f.endswith(".cube"):
                try:
                    os.remove(os.path.join(CUBE_DIR, f))
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
    surfaces = surface_keys or []
    esps = esp_keys or []
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
        "files": {
            "trajectory_xyz": "multistep_pathway.xyz",
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
    """Unpack a base64-encoded Molecular Workbench zip bundle and extract orbital meshes/xyz/esp."""
    raw_zip = base64.b64decode(data_b64)
    orbitals_found: List[str] = []
    files_dict: Dict[str, Any] = {}
    meta: Dict[str, Any] = {
        "xyz": "",
        "orbitals_tree": [],
        "esp_surface": None,
        "molecule_key": "",
        "smiles": ""
    }
    with zipfile.ZipFile(io.BytesIO(raw_zip), "r") as zf:
        for name in zf.namelist():
            try:
                content = zf.read(name).decode("utf-8", errors="ignore")
            except Exception:
                continue
            files_dict[name] = content
            files_dict[os.path.basename(name)] = content
            if name.endswith(".xyz") and not meta["xyz"]:
                meta["xyz"] = content
            elif name.endswith(".cube"):
                orb = os.path.basename(name).replace(".cube", "")
                if orb not in orbitals_found:
                    orbitals_found.append(orb)
            elif name.endswith("_pos.json") or name.endswith("_neg.json"):
                orb = os.path.basename(name).replace("_pos.json", "").replace("_neg.json", "")
                if orb not in orbitals_found:
                    orbitals_found.append(orb)

        mol_data = None
        for name in zf.namelist():
            if name.startswith("molecules/") and name.endswith(".json") and not any(name.endswith(x) for x in ["_esp.json", "_pos.json", "_neg.json"]):
                try:
                    mol_data = json.loads(files_dict.get(name, "{}"))
                    break
                except Exception:
                    pass

        if not mol_data and "manifest.json" in files_dict:
            try:
                man = json.loads(files_dict["manifest.json"])
                mols = man.get("molecules", [])
                mol_data = mols[0] if (mols and isinstance(mols, list)) else man
            except Exception:
                pass

        if mol_data:
            meta["molecule_key"] = mol_data.get("key", mol_data.get("name", ""))
            meta["smiles"] = mol_data.get("smiles", "")
            if not meta["xyz"] and mol_data.get("xyz"):
                meta["xyz"] = mol_data.get("xyz")
            if mol_data.get("orbitals"):
                meta["orbitals_tree"] = mol_data.get("orbitals")
                for cat in meta["orbitals_tree"]:
                    for it in cat.get("items", []):
                        f = it.get("file", "")
                        b = os.path.basename(f).replace(".json", "").replace(".cube", "")
                        if b and b not in orbitals_found:
                            orbitals_found.append(b)
            if mol_data.get("esp_surface"):
                meta["esp_surface"] = mol_data.get("esp_surface")

        if not meta["esp_surface"]:
            for fname in files_dict:
                if fname.endswith("_esp.json") or fname == "esp.json" or fname.endswith("/esp.json"):
                    meta["esp_surface"] = {"file": fname}
                    break

        mol_key = meta["molecule_key"]
        if mol_key:
            workspace_root = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
            local_full_zip = os.path.join(workspace_root, "output", mol_key, f"{mol_key}.rxnbundle.zip")
            if os.path.exists(local_full_zip):
                try:
                    with zipfile.ZipFile(local_full_zip, "r") as lzf:
                        for lname in lzf.namelist():
                            if lname not in files_dict and (lname.endswith(".json") or lname.endswith(".xyz")):
                                lcontent = lzf.read(lname).decode("utf-8", errors="ignore")
                                files_dict[lname] = lcontent
                                files_dict[os.path.basename(lname)] = lcontent
                                if lname.endswith("_pos.json") or lname.endswith("_neg.json"):
                                    orb = os.path.basename(lname).replace("_pos.json", "").replace("_neg.json", "")
                                    if orb not in orbitals_found:
                                        orbitals_found.append(orb)
                except Exception:
                    pass

    return MoleculeBundleResult({
        "orbitals": orbitals_found,
        "files": files_dict,
        "xyz": meta["xyz"],
        "orbitals_tree": meta["orbitals_tree"],
        "esp_surface": meta["esp_surface"],
        "molecule_key": meta["molecule_key"],
        "smiles": meta["smiles"]
    })



def _convert_and_pack_orbitals(megajson: Dict[str, Any], zf: zipfile.ZipFile, isovalue: float = 0.035) -> Tuple[List[str], int]:
    """Calculate lightweight 3D surface meshes from cubes, pack them into zip, and never pack .cube files."""
    sources: Dict[str, Any] = {}
    if os.path.exists(CUBE_DIR):
        sources.update({f[:-5]: os.path.join(CUBE_DIR, f) for f in os.listdir(CUBE_DIR) if f.endswith(".cube")})

    for k, text in _HIGHRES_CUBE_CACHE.items():
        if "_" in k and not k.startswith("res"):
            parts = k.split("_")
            clean_k = f"{parts[1]}_{'_'.join(parts[2:])}" if (len(parts) > 2 and parts[0].isdigit() and parts[1].isdigit()) else k
            sources.setdefault(clean_k, text)

    for k, val in megajson.get("highres_cubes", {}).items():
        if k not in sources:
            if isinstance(val, str) and "\n" in val:
                sources[k] = val
            elif isinstance(val, str) and val.endswith(".cube"):
                cp = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", val))
                if os.path.exists(cp):
                    sources[k] = cp

    # 4. From Step4 _CUBE_CACHE
    try:
        from ..step_4_orbitals import _CUBE_CACHE
        for k, text in _CUBE_CACHE.items():
            if k not in sources:
                sources[k] = text
    except Exception:
        pass

    if not sources:
        return [], 0

    print(f"[Packager] Converting {len(sources)} orbital cubes to 3D surface meshes (isovalue={isovalue})...")

    def _extract_mesh(item):
        key, source = item
        try:
            vol, origin, step = CubeConverter.read_cube(source)
            pos = CubeConverter.extract_mesh(vol, origin, step, isovalue)
            neg = CubeConverter.extract_mesh(vol, origin, step, -isovalue)
            return key, pos, neg
        except Exception as e:
            print(f"[Packager] Warning: failed to extract mesh for {key}: {e}")
            return key, None, None

    packed_keys = []
    uncompressed_bytes = 0

    with ThreadPoolExecutor(max_workers=min(8, os.cpu_count() or 4)) as executor:
        for key, pos, neg in executor.map(_extract_mesh, sources.items()):
            has_mesh = False
            if pos:
                p_bytes = json.dumps(pos, separators=(',', ':')).encode('utf-8')
                zf.writestr(f"orbitals/{key}_pos.json", p_bytes)
                uncompressed_bytes += len(p_bytes)
                has_mesh = True
            if neg:
                n_bytes = json.dumps(neg, separators=(',', ':')).encode('utf-8')
                zf.writestr(f"orbitals/{key}_neg.json", n_bytes)
                uncompressed_bytes += len(n_bytes)
                has_mesh = True
            if has_mesh:
                packed_keys.append(key)

    print(f"[Packager] Successfully packed {len(packed_keys)} orbital surface pairs into archive.")
    return packed_keys, uncompressed_bytes


def create_maxcomp_bundle_zip(megajson: Dict[str, Any], global_state: Dict[str, Any]) -> Tuple[bytes, int]:
    """Calculate surfaces and pack into zip with max compression (level 9). Never pack .cube files. Purge temporary cubes."""
    zip_buf = io.BytesIO()
    isovalue = float(megajson.get("orbital_isovalue") or megajson.get("isovalue") or 0.035)

    with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        # 1. Convert cubes to surfaces and pack into orbitals/
        surface_keys, orb_bytes = _convert_and_pack_orbitals(megajson, zf, isovalue)

        # 2. Pack ESP surfaces into esp/
        esp_surfaces = megajson.get("esp_surfaces", {})
        esp_keys, esp_bytes = [], 0
        for f_key, mesh in esp_surfaces.items():
            if not mesh or not isinstance(mesh, dict):
                continue
            clean_f = str(f_key).replace("frame_", "")
            mesh_json = json.dumps(mesh, separators=(',', ':')).encode('utf-8')
            zf.writestr(f"esp/frame_{clean_f}.json", mesh_json)
            esp_bytes += len(mesh_json)
            esp_keys.append(clean_f)

        # 3. Store uploaded molecule bundles as UNPACKED data (no nested zip archives!)
        molecule_bundles = megajson.get("molecule_bundles", {})
        clean_mol_bundles = {}
        mol_bytes = 0
        for m_id, b_data in molecule_bundles.items():
            mol_key = b_data.get("molecule_key") or m_id
            mol_dir = f"molecules/{mol_key}"
            written_files = set()

            if "data_b64" in b_data:
                try:
                    raw_b = base64.b64decode(b_data["data_b64"])
                    with zipfile.ZipFile(io.BytesIO(raw_b), "r") as zin:
                        for item in zin.infolist():
                            if item.is_dir() or item.filename.endswith(".cube"):
                                continue
                            sub = item.filename
                            if sub.startswith("molecules/"):
                                sub = sub[10:]
                            if sub.startswith(f"{mol_key}/"):
                                sub = sub[len(mol_key) + 1:]
                            dest = f"{mol_dir}/{sub}"
                            content = zin.read(item.filename)
                            zf.writestr(dest, content)
                            mol_bytes += len(content)
                            written_files.add(dest)
                except Exception as e:
                    print(f"[Packager] Warning unpacking molecule {m_id}: {e}")

            for fname, fcontent in b_data.get("files", {}).items():
                if fname.endswith(".cube"):
                    continue
                dest = f"{mol_dir}/{os.path.basename(fname)}"
                if dest not in written_files:
                    c_bytes = fcontent.encode("utf-8") if isinstance(fcontent, str) else fcontent
                    zf.writestr(dest, c_bytes)
                    mol_bytes += len(c_bytes)
                    written_files.add(dest)

            clean_mol_bundles[m_id] = {
                "molecule_id": m_id,
                "molecule_key": mol_key,
                "directory": f"{mol_dir}/",
                "orbitals": b_data.get("orbitals", []),
                "orbitals_tree": b_data.get("orbitals_tree", []),
                "xyz": b_data.get("xyz", ""),
                "esp_surface": b_data.get("esp_surface"),
                "file_count": len(written_files)
            }

        # 4. Clean megajson (remove cubes & raw b64 strings)
        clean_megajson = dict(megajson)
        clean_megajson.pop("highres_cubes", None)
        clean_megajson["orbital_surfaces"] = sorted(surface_keys)
        clean_megajson["esp_frame_keys"] = sorted(esp_keys)
        clean_megajson["molecule_bundles"] = clean_mol_bundles
        megajson.pop("highres_cubes", None)
        megajson["orbital_surfaces"] = sorted(surface_keys)

        # 5. Build manifest and metadata
        manifest = build_bundle_manifest(clean_megajson, global_state, surface_keys, esp_keys)
        xyz_content = generate_xyz(clean_megajson)
        json_content = json.dumps(clean_megajson, indent=2)
        standalone_html = generate_standalone_html(clean_megajson, manifest)

        uncompressed_bytes = (
            orb_bytes + esp_bytes + mol_bytes +
            len(xyz_content.encode("utf-8")) +
            len(json_content.encode("utf-8")) +
            len(standalone_html.encode("utf-8"))
        )

        zf.writestr("bundle.json", json.dumps(manifest, indent=2))
        zf.writestr("multistep_pathway.xyz", xyz_content)
        zf.writestr("multistep_pathway.json", json_content)
        zf.writestr("index.html", standalone_html)

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
    print(f"[Packager] Purged {purged_count} temporary cube files after surface packaging.")

    zip_buf.seek(0)
    return zip_buf.getvalue(), uncompressed_bytes
