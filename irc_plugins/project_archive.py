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

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
STEPS_DIR = os.path.join(OUT_DIR, 'steps')
PROJECT_FILE = os.path.join(OUT_DIR, 'project.json')


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

        # Write metadata manifest
        manifest = {
            "version": "1.0",
            "exported_at": time.strftime("%Y-%m-%d %H:%M:%S"),
            "steps_count": len(steps),
            "solvent": proj_data.get("solvent", "water")
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

            # Extract recognized project files
            for member in namelist:
                # Prevent directory traversal attacks
                norm = os.path.normpath(member)
                if norm.startswith("..") or os.path.isabs(norm):
                    continue
                if norm == "project.json" or norm.startswith("steps" + os.sep) or norm.startswith("steps/"):
                    dest_p = os.path.join(out_dir, norm)
                    os.makedirs(os.path.dirname(dest_p), exist_ok=True)
                    with zf.open(member) as src, open(dest_p, "wb") as dst:
                        shutil.copyfileobj(src, dst)

        # Re-initialize project and load active state
        proj = project_manager.init_project()
        active_state = project_manager.load_active_state()
        summary = project_manager.get_project_summary()

        n_steps = len(proj.get("steps", []))
        msg = f"Project successfully loaded with {n_steps} reaction step{'s' if n_steps != 1 else ''}!"
        return True, summary, active_state, msg
    except Exception as e:
        return False, {}, {}, f"Failed to import project archive: {e}"
