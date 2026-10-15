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
    frames = megajson.get("frames", [])
    lines = []
    for f_idx, frame in enumerate(frames):
        coords = frame.get("coords", frame.get("positions", []))
        step_idx = frame.get("step_index", 0)
        energy = frame.get("rel_energy_kcal")
        if energy is None:
            energy = frame.get("energy_kcal")
        if energy is None:
            energy = 0.0
        lines.append(f"{n_atoms}")
        lines.append(f"Frame {f_idx} | Step {step_idx} | Relative Free Energy: {energy:.2f} kcal/mol")
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


def build_bundle_manifest(megajson: Dict[str, Any], global_state: Dict[str, Any], surface_keys: List[str] = None) -> Dict[str, Any]:
    """Build standardized bundle.json manifest with surface mesh listings."""
    boundaries = megajson.get("step_boundaries", [])
    title = " → ".join(b.get("step_name", f"Step {i+1}") for i, b in enumerate(boundaries)) or "Multistep Reaction"
    rxn_smi = " -> ".join(b.get("reaction_smiles", "") for b in boundaries)
    surfaces = surface_keys or []
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
        "files": {
            "trajectory_xyz": "multistep_pathway.xyz",
            "trajectory_json": "multistep_pathway.json",
            "mechanism_svg": "mechanism.svg" if megajson.get("mechanism_svg") else None,
            "svg_wiring": "svg_wiring.json" if megajson.get("svg_wiring") else None,
            "orbitals_directory": "orbitals/",
            "molecule_bundles": "molecule_bundles/"
        }
    }


def unpack_molecule_bundle(data_b64: str) -> Tuple[List[str], Dict[str, Any]]:
    """Unpack a base64-encoded Molecular Workbench zip bundle and extract orbital meshes/xyz."""
    raw_zip = base64.b64decode(data_b64)
    orbitals_found = []
    files_dict = {}
    with zipfile.ZipFile(io.BytesIO(raw_zip), "r") as zf:
        for name in zf.namelist():
            if name.endswith(".cube"):
                orb_name = os.path.basename(name).replace(".cube", "")
                orbitals_found.append(orb_name)
                files_dict[orb_name] = zf.read(name).decode("utf-8", errors="ignore")
            elif name.endswith("_pos.json") or name.endswith("_neg.json"):
                base_orb = os.path.basename(name).replace("_pos.json", "").replace("_neg.json", "")
                if base_orb not in orbitals_found:
                    orbitals_found.append(base_orb)
                files_dict[os.path.basename(name)] = zf.read(name).decode("utf-8", errors="ignore")
            elif name.endswith(".xyz"):
                files_dict["xyz"] = zf.read(name).decode("utf-8", errors="ignore")
            elif name.endswith(".json"):
                files_dict["json"] = zf.read(name).decode("utf-8", errors="ignore")
    return orbitals_found, files_dict


def _convert_and_pack_orbitals(megajson: Dict[str, Any], zf: zipfile.ZipFile, isovalue: float = 0.035) -> Tuple[List[str], int]:
    """Calculate lightweight 3D surface meshes from cubes, pack them into zip, and never pack .cube files."""
    sources: Dict[str, Any] = {}

    # 1. From CUBE_DIR
    if os.path.exists(CUBE_DIR):
        for f in os.listdir(CUBE_DIR):
            if f.endswith(".cube"):
                k = f[:-5]
                sources[k] = os.path.join(CUBE_DIR, f)

    # 2. From _HIGHRES_CUBE_CACHE
    for k, text in _HIGHRES_CUBE_CACHE.items():
        if "_" in k and not k.startswith("res"):
            parts = k.split("_")
            if len(parts) > 2 and parts[0].isdigit() and parts[1].isdigit():
                clean_k = f"{parts[1]}_{'_'.join(parts[2:])}"
                if clean_k not in sources:
                    sources[clean_k] = text
            else:
                sources[k] = text

    # 3. From megajson["highres_cubes"]
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

        # 2. Update megajson to reference orbital surfaces instead of cubes
        clean_megajson = dict(megajson)
        clean_megajson.pop("highres_cubes", None)
        clean_megajson["orbital_surfaces"] = sorted(surface_keys)
        megajson.pop("highres_cubes", None)
        megajson["orbital_surfaces"] = sorted(surface_keys)

        # 3. Build manifest and metadata
        manifest = build_bundle_manifest(clean_megajson, global_state, surface_keys)
        xyz_content = generate_xyz(clean_megajson)
        json_content = json.dumps(clean_megajson, indent=2)
        standalone_html = generate_standalone_html(clean_megajson, manifest)

        uncompressed_bytes = (
            orb_bytes +
            len(xyz_content.encode("utf-8")) +
            len(json_content.encode("utf-8")) +
            len(standalone_html.encode("utf-8"))
        )

        zf.writestr("bundle.json", json.dumps(manifest, indent=2))
        zf.writestr("multistep_pathway.xyz", xyz_content)
        zf.writestr("multistep_pathway.json", json_content)
        zf.writestr("index.html", standalone_html)

        # Store Mechanism SVG if available
        svg_content = clean_megajson.get("mechanism_svg")
        if svg_content:
            zf.writestr("mechanism.svg", svg_content)
            uncompressed_bytes += len(svg_content.encode("utf-8"))

        # Store SVG Wiring if available
        wiring = clean_megajson.get("svg_wiring")
        if wiring:
            wiring_str = json.dumps(wiring, indent=2)
            zf.writestr("svg_wiring.json", wiring_str)
            uncompressed_bytes += len(wiring_str.encode("utf-8"))

        # Store uploaded molecule bundles without any .cube files
        molecule_bundles = clean_megajson.get("molecule_bundles", {})
        for m_id, b_data in molecule_bundles.items():
            if "data_b64" in b_data:
                raw_b = base64.b64decode(b_data["data_b64"])
                try:
                    cbuf = io.BytesIO()
                    with zipfile.ZipFile(io.BytesIO(raw_b), "r") as zin, zipfile.ZipFile(cbuf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zout:
                        for item in zin.infolist():
                            if not item.filename.endswith(".cube"):
                                zout.writestr(item, zin.read(item.filename))
                    cbytes = cbuf.getvalue()
                except Exception:
                    cbytes = raw_b
                fn = f"molecule_bundles/{m_id}_{b_data.get('filename', 'bundle.zip')}"
                zf.writestr(fn, cbytes)
                uncompressed_bytes += len(cbytes)

    # 4. Purge temporary cubes to save space!
    purged_count = purge_cube_files()
    print(f"[Packager] Purged {purged_count} temporary cube files after surface packaging.")

    zip_buf.seek(0)
    return zip_buf.getvalue(), uncompressed_bytes
