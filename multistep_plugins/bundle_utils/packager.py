import base64
import io
import json
import os
import zipfile
from typing import Any, Dict, List, Tuple
from .standalone_template import generate_standalone_html
from .highres_engine import _HIGHRES_CUBE_CACHE


def generate_xyz(megajson: Dict[str, Any]) -> str:
    """Generate multi-model XYZ text format from pathway frames."""
    symbols = megajson.get("symbols", [])
    n_atoms = len(symbols)
    frames = megajson.get("frames", [])
    lines = []
    for f_idx, frame in enumerate(frames):
        coords = frame.get("coords", frame.get("positions", []))
        energy = frame.get("rel_energy_kcal", frame.get("energy_kcal", 0.0))
        step_idx = frame.get("step_index", 0)
        lines.append(f"{n_atoms}")
        lines.append(f"Frame {f_idx} | Step {step_idx} | Relative Free Energy: {energy:.2f} kcal/mol")
        for sym, pos in zip(symbols, coords):
            lines.append(f"{sym:<2} {pos[0]:12.6f} {pos[1]:12.6f} {pos[2]:12.6f}")
    return "\n".join(lines) + "\n"


def build_bundle_manifest(megajson: Dict[str, Any], global_state: Dict[str, Any]) -> Dict[str, Any]:
    """Build standardized bundle.json manifest."""
    boundaries = megajson.get("step_boundaries", [])
    title = " → ".join(b.get("step_name", f"Step {i+1}") for i, b in enumerate(boundaries)) or "Multistep Reaction"
    rxn_smi = " -> ".join(b.get("reaction_smiles", "") for b in boundaries)
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
    """Unpack a base64-encoded Molecular Workbench zip bundle and extract orbital cubes/xyz."""
    raw_zip = base64.b64decode(data_b64)
    orbitals_found = []
    files_dict = {}
    with zipfile.ZipFile(io.BytesIO(raw_zip), "r") as zf:
        for name in zf.namelist():
            if name.endswith(".cube"):
                orb_name = os.path.basename(name).replace(".cube", "")
                orbitals_found.append(orb_name)
                files_dict[orb_name] = zf.read(name).decode("utf-8", errors="ignore")
            elif name.endswith(".xyz"):
                files_dict["xyz"] = zf.read(name).decode("utf-8", errors="ignore")
            elif name.endswith(".json"):
                files_dict["json"] = zf.read(name).decode("utf-8", errors="ignore")
    return orbitals_found, files_dict


def create_maxcomp_bundle_zip(megajson: Dict[str, Any], global_state: Dict[str, Any]) -> Tuple[bytes, int]:
    """Repack all data, SVG, wiring, high-res cubes, molecule bundles, and standalone viewer with maximal compression (level 9)."""
    zip_buf = io.BytesIO()
    manifest = build_bundle_manifest(megajson, global_state)
    xyz_content = generate_xyz(megajson)
    json_content = json.dumps(megajson, indent=2)
    standalone_html = generate_standalone_html(megajson, manifest)

    uncompressed_bytes = (
        len(xyz_content.encode("utf-8")) +
        len(json_content.encode("utf-8")) +
        len(standalone_html.encode("utf-8"))
    )

    with zipfile.ZipFile(zip_buf, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        zf.writestr("bundle.json", json.dumps(manifest, indent=2))
        zf.writestr("multistep_pathway.xyz", xyz_content)
        zf.writestr("multistep_pathway.json", json_content)
        zf.writestr("index.html", standalone_html)

        # Store Mechanism SVG if available
        svg_content = megajson.get("mechanism_svg")
        if svg_content:
            zf.writestr("mechanism.svg", svg_content)
            uncompressed_bytes += len(svg_content.encode("utf-8"))

        # Store SVG Wiring if available
        wiring = megajson.get("svg_wiring")
        if wiring:
            wiring_str = json.dumps(wiring, indent=2)
            zf.writestr("svg_wiring.json", wiring_str)
            uncompressed_bytes += len(wiring_str.encode("utf-8"))

        # Store all high-resolution cubes
        highres_cubes = megajson.get("highres_cubes", {})
        for key, cube_text in highres_cubes.items():
            if cube_text:
                fn = f"orbitals/{key}.cube"
                zf.writestr(fn, cube_text)
                uncompressed_bytes += len(cube_text.encode("utf-8"))

        # Also store any cached cubes from _HIGHRES_CUBE_CACHE
        for k, cube_text in _HIGHRES_CUBE_CACHE.items():
            fn = f"orbitals/{k}.cube"
            if fn not in zf.namelist():
                zf.writestr(fn, cube_text)
                uncompressed_bytes += len(cube_text.encode("utf-8"))

        # Store uploaded molecule bundles from Molecular Workbench
        molecule_bundles = megajson.get("molecule_bundles", {})
        for m_id, b_data in molecule_bundles.items():
            if "data_b64" in b_data:
                raw_b = base64.b64decode(b_data["data_b64"])
                zf.writestr(f"molecule_bundles/{m_id}_{b_data.get('filename', 'bundle.zip')}", raw_b)
                uncompressed_bytes += len(raw_b)

    zip_buf.seek(0)
    return zip_buf.getvalue(), uncompressed_bytes
