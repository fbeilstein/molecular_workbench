"""
Project Archival Module
Handles packaging active reaction projects into portable .zip archives
and restoring archived projects into the active workspace.
"""

import io
import json
import os
import shutil
import time
import zipfile
from typing import Any, Dict, Optional, Tuple

from .constants import (
    ROOT_DIR,
    OUT_DIR,
    STEPS_DIR,
    PROJECT_FILE,
    CASCADE_DIR,
    PROLONGED_PATH,
)


def export_project_archive(out_dir: str = OUT_DIR) -> Tuple[bytes, str]:
    """
    Export the current multi-step project into a zip archive.
    Returns (zip_bytes, suggested_filename).
    """
    proj_path = os.path.join(out_dir, 'project.json')
    steps_path = os.path.join(out_dir, 'steps')

    proj_data = {}
    if os.path.exists(proj_path):
        try:
            with open(proj_path, 'r') as f:
                proj_data = json.load(f)
        except Exception:
            pass

    # Determine filename from project title / first step
    steps = proj_data.get("steps", [])
    slug = "reaction_project"
    if steps and steps[0].get("smiles"):
        smi = steps[0]["smiles"].split(">>")[0]
        safe_smi = "".join(c if c.isalnum() else "_" for c in smi)[:24].strip("_")
        if safe_smi:
            slug = f"project_{safe_smi}"

    timestamp = time.strftime("%Y%m%d_%H%M%S")
    filename = f"{slug}_{timestamp}.zip"

    buf = io.BytesIO()
    with zipfile.ZipFile(buf, 'w', compression=zipfile.ZIP_DEFLATED) as zf:
        # Write project.json
        if os.path.exists(proj_path):
            zf.write(proj_path, arcname="project.json")
        elif proj_data:
            zf.writestr("project.json", json.dumps(proj_data, indent=2))

        # Write each step directory
        if os.path.exists(steps_path):
            for root, dirs, files in os.walk(steps_path):
                for file in files:
                    full_p = os.path.join(root, file)
                    rel_p = os.path.relpath(full_p, out_dir)
                    zf.write(full_p, arcname=rel_p)

        # Write root master files (stitched cascade pathway, global state, orbital tracks)
        for extra_file in ["reaction_state.json", "master_pathway.json", "curated_orbitals_backup.json"]:
            full_p = os.path.join(out_dir, extra_file)
            if os.path.exists(full_p):
                zf.write(full_p, arcname=extra_file)

        # Write generated bundles directory if present
        bundles_path = os.path.join(out_dir, 'bundles')
        if os.path.exists(bundles_path):
            for root, dirs, files in os.walk(bundles_path):
                for file in files:
                    full_p = os.path.join(root, file)
                    rel_p = os.path.relpath(full_p, out_dir)
                    zf.write(full_p, arcname=rel_p)

        # Write precalculated cubes directory if present
        cubes_path = os.path.join(out_dir, 'cubes')
        if os.path.exists(cubes_path):
            for root, dirs, files in os.walk(cubes_path):
                for file in files:
                    full_p = os.path.join(root, file)
                    rel_p = os.path.relpath(full_p, out_dir)
                    zf.write(full_p, arcname=rel_p)

        # Write comprehensive metadata manifest
        manifest = {
            "version": "2.0",
            "exported_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "reaction_steps_count": len(steps),
            "reaction_steps": [s.get("id") for s in steps],
            "solvent": proj_data.get("solvent", "water"),
            "pipeline_phases_included": [
                "Phase 1: Reaction Input",
                "Phase 2: Atom Mapping",
                "Phase 3: Transition State (TS)",
                "Phase 4: IRC Trajectories",
                "Phase 5: Master Cascade Stitched Pathway",
                "Phase 6: Molecular Orbitals Timeline",
                "Phase 7: Packaging & Bundle Viewer"
            ],
            "has_master_pathway": os.path.exists(os.path.join(out_dir, "master_pathway.json")),
            "has_curated_orbitals": os.path.exists(os.path.join(out_dir, "curated_orbitals_backup.json")),
            "has_bundles": os.path.exists(bundles_path) and len(os.listdir(bundles_path)) > 0
        }
        zf.writestr("ARCHIVE_MANIFEST.json", json.dumps(manifest, indent=2))

    buf.seek(0)
    return buf.getvalue(), filename


def import_project_archive(zip_bytes: bytes, out_dir: str = OUT_DIR) -> Tuple[bool, Dict[str, Any], Dict[str, Any], str]:
    """
    Restore a project from a zip archive into out_dir.
    Returns (success, project_summary, active_state, message).
    """
    from . import project_manager

    if not zip_bytes:
        return False, {}, {}, "Empty archive payload."

    buf = io.BytesIO(zip_bytes)
    if not zipfile.is_zipfile(buf):
        return False, {}, {}, "Uploaded file is not a valid zip archive."

    buf.seek(0)
    try:
        with zipfile.ZipFile(buf, 'r') as zf:
            namelist = zf.namelist()
            if "project.json" not in namelist and not any(n.startswith("steps/") for n in namelist):
                return False, {}, {}, "Archive missing project.json or steps/ directory."

            # Clean existing steps directory to prevent old step leakage
            steps_path = os.path.join(out_dir, 'steps')
            if os.path.exists(steps_path):
                shutil.rmtree(steps_path, ignore_errors=True)
            os.makedirs(steps_path, exist_ok=True)

            allowed_prefixes = ("steps/", "bundles/", "cubes/")
            allowed_files = {"project.json", "reaction_state.json", "master_pathway.json", "curated_orbitals_backup.json", "irc_trajectory.json"}

            # Extract recognized project files
            for member in namelist:
                # Prevent directory traversal attacks
                norm = os.path.normpath(member)
                if norm.startswith("..") or os.path.isabs(norm):
                    continue
                if norm in allowed_files or any(norm.startswith(p) for p in allowed_prefixes):
                    dest_p = os.path.join(out_dir, norm)
                    os.makedirs(os.path.dirname(dest_p), exist_ok=True)
                    with zf.open(member) as src, open(dest_p, "wb") as dst:
                        shutil.copyfileobj(src, dst)

        # Sync master_pathway.json to CASCADE_DIR/prolonged_pathway.json and out/irc_trajectory.json
        master_p = os.path.join(out_dir, "master_pathway.json")
        if os.path.exists(master_p):
            os.makedirs(CASCADE_DIR, exist_ok=True)
            shutil.copyfile(master_p, PROLONGED_PATH)
            shutil.copyfile(master_p, os.path.join(out_dir, "irc_trajectory.json"))

        # Re-initialize project and load active state
        proj = project_manager.init_project()
        active_state = project_manager.load_active_state()

        # Ensure cascade phases (Step 5, Step 6, Step 7) from top-level reaction_state.json
        # are preserved in active_state
        main_state_p = os.path.join(out_dir, "reaction_state.json")
        if os.path.exists(main_state_p):
            try:
                with open(main_state_p, 'r') as f:
                    main_st = json.load(f)
                for k in ["step_5_cascade", "step_6_orbitals", "step_7_bundle"]:
                    if k in main_st:
                        active_state[k] = main_st[k]
            except Exception:
                pass

        summary = project_manager.get_project_summary()

        n_steps = len(proj.get("steps", []))
        msg = f"Project successfully loaded with {n_steps} reaction step{'s' if n_steps != 1 else ''}!"
        return True, summary, active_state, msg
    except Exception as e:
        return False, {}, {}, f"Failed to import project archive: {e}"
