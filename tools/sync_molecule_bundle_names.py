#!/usr/bin/env python3
"""
sync_molecule_bundle_names.py

Checks and synchronizes internal naming within single-molecule bundles
in bundles/molecules/ so that internal manifest fields, molecule JSON files,
orbital item file references, and mesh file prefixes match the bundle's key.

Usage:
    python tools/sync_molecule_bundle_names.py --check
    python tools/sync_molecule_bundle_names.py --fix [--dry-run]
    python tools/sync_molecule_bundle_names.py --fix-zip-names [--dry-run]
"""

import argparse
import glob
import io
import json
import os
import shutil
import sys
import zipfile

# Recommended corrections for zip filenames (typos and naming inconsistencies)
RECOMMENDED_ZIP_RENAMES = {
    "cloride_ion": "chloride",          # Typo fix (cloride -> chloride) and match systematic naming (no _ion)
    "methanoate_ion": "methanoate",    # Match ethanoate, butanoate (no _ion)
    "1,1,1,3,3,3-hexafluoro-2-propanol": "1,1,1,3,3,3-hexafluoropropan-2-ol", # IUPAC PIN (matches conjugate base)
    "methyl_bromide": "bromomethane",  # IUPAC PIN substitutive nomenclature
}


def audit_bundle(zip_path, target_key=None):
    """
    Inspect a single-molecule bundle and return a report of discrepancies.
    """
    base_name = os.path.splitext(os.path.basename(zip_path))[0]
    key = target_key or base_name

    issues = []
    namelist = []
    manifest = {}

    try:
        with zipfile.ZipFile(zip_path, "r") as zf:
            namelist = zf.namelist()
            if "manifest.json" not in namelist:
                issues.append("Missing manifest.json")
                return {
                    "zip_path": zip_path,
                    "base_name": base_name,
                    "target_key": key,
                    "old_key": None,
                    "issues": issues,
                    "is_synced": False,
                }

            manifest = json.loads(zf.read("manifest.json").decode("utf-8"))
    except Exception as e:
        return {
            "zip_path": zip_path,
            "base_name": base_name,
            "target_key": key,
            "old_key": None,
            "issues": [f"Failed to read archive: {e}"],
            "is_synced": False,
        }

    manifest_name = manifest.get("name")
    manifest_title = manifest.get("title")

    mols = manifest.get("molecules", [])
    m0_key = mols[0].get("key") if mols else None

    # Detect old key
    old_key = manifest_name or m0_key

    # Check manifest fields
    if manifest_name != key:
        issues.append(f"manifest.name is '{manifest_name}' (expected '{key}')")
    if manifest_title != key:
        issues.append(f"manifest.title is '{manifest_title}' (expected '{key}')")
    if m0_key != key:
        issues.append(f"manifest.molecules[0].key is '{m0_key}' (expected '{key}')")

    # Check main molecule json
    expected_mol_json = f"molecules/{key}.json"
    if expected_mol_json not in namelist:
        found_mol_jsons = [
            n for n in namelist
            if n.startswith("molecules/") and n.count("/") == 1
            and not (n.endswith("_pos.json") or n.endswith("_neg.json") or n.endswith("_esp.json"))
        ]
        issues.append(f"Main mol JSON '{expected_mol_json}' not found (present: {found_mol_jsons})")

    # Check internal file references in manifest
    if mols:
        m0 = mols[0]
        esp = m0.get("esp_surface")
        if isinstance(esp, dict):
            esp_file = esp.get("file", "")
            expected_esp = f"molecules/{key}_esp.json"
            if esp_file != expected_esp:
                issues.append(f"ESP surface reference '{esp_file}' != '{expected_esp}'")

        bad_orb_refs = 0
        for cat in m0.get("orbitals", []):
            for it in cat.get("items", []):
                f = it.get("file", "")
                if f and not f.startswith(f"molecules/{key}_"):
                    bad_orb_refs += 1
        if bad_orb_refs > 0:
            issues.append(f"{bad_orb_refs} orbital item file references do not match prefix 'molecules/{key}_'")

    # Check mesh files prefix
    mesh_files = [n for n in namelist if n.endswith("_pos.json") or n.endswith("_neg.json") or n.endswith("_esp.json")]
    bad_mesh_files = [n for n in mesh_files if not n.startswith(f"molecules/{key}_")]
    if bad_mesh_files:
        issues.append(f"{len(bad_mesh_files)} mesh files do not start with 'molecules/{key}_' (sample: '{bad_mesh_files[0]}')")

    return {
        "zip_path": zip_path,
        "base_name": base_name,
        "target_key": key,
        "old_key": old_key,
        "issues": issues,
        "is_synced": len(issues) == 0,
        "total_files": len(namelist),
    }


def _remap_mol_json_content(content: bytes, old_key: str, new_key: str) -> bytes:
    """Update internal molecule key and orbital file paths inside molecule JSON."""
    try:
        mol_data = json.loads(content.decode("utf-8"))
        mol_data["key"] = new_key
        old_pfx = f"molecules/{old_key}_"
        new_pfx = f"molecules/{new_key}_"

        esp = mol_data.get("esp_surface")
        if isinstance(esp, dict):
            f = esp.get("file", "")
            if f.startswith(old_pfx):
                esp["file"] = f.replace(old_pfx, new_pfx, 1)

        for cat in mol_data.get("orbitals", []):
            for it in cat.get("items", []):
                f = it.get("file", "")
                if f.startswith(old_pfx):
                    it["file"] = f.replace(old_pfx, new_pfx, 1)

        return json.dumps(mol_data, indent=2).encode("utf-8")
    except Exception as e:
        print(f"Warning: could not parse molecule JSON: {e}")
        return content


def fix_bundle(zip_path, target_key=None, dry_run=False, backup=False):
    """
    Rewrite bundle archive so that manifest and all internal files/references
    consistently use target_key.
    """
    base_name = os.path.splitext(os.path.basename(zip_path))[0]
    new_key = target_key or base_name

    with zipfile.ZipFile(zip_path, "r") as zf_in:
        namelist = zf_in.namelist()
        if "manifest.json" not in namelist:
            raise ValueError(f"No manifest.json in {zip_path}")

        manifest = json.loads(zf_in.read("manifest.json").decode("utf-8"))
        old_key = manifest.get("name")
        if not old_key and manifest.get("molecules"):
            old_key = manifest["molecules"][0].get("key")

        if not old_key:
            # Fallback to base name of any non-mesh mol json
            for n in namelist:
                if n.startswith("molecules/") and n.count("/") == 1 and not (
                    n.endswith("_pos.json") or n.endswith("_neg.json") or n.endswith("_esp.json")
                ):
                    old_key = os.path.splitext(os.path.basename(n))[0]
                    break

        if not old_key:
            old_key = base_name

        # Update manifest dict
        manifest["name"] = new_key
        manifest["title"] = new_key

        if manifest.get("molecules"):
            m0 = manifest["molecules"][0]
            m0["key"] = new_key
            if isinstance(m0.get("esp_surface"), dict):
                f = m0["esp_surface"].get("file", "")
                if f.startswith(f"molecules/{old_key}_"):
                    m0["esp_surface"]["file"] = f.replace(f"molecules/{old_key}_", f"molecules/{new_key}_", 1)
                elif f == f"molecules/{old_key}.json":
                    m0["esp_surface"]["file"] = f"molecules/{new_key}_esp.json"
            for cat in m0.get("orbitals", []):
                for it in cat.get("items", []):
                    f = it.get("file", "")
                    if f.startswith(f"molecules/{old_key}_"):
                        it["file"] = f.replace(f"molecules/{old_key}_", f"molecules/{new_key}_", 1)

        if dry_run:
            return True

        out_buf = io.BytesIO()
        with zipfile.ZipFile(out_buf, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf_out:
            # 1. Write updated manifest
            zf_out.writestr("manifest.json", json.dumps(manifest, indent=2))

            # 2. Process all other files
            for fname in namelist:
                if fname == "manifest.json":
                    continue

                content = zf_in.read(fname)
                new_fname = fname

                # Main molecule json
                if fname == f"molecules/{old_key}.json" or (
                    fname.startswith("molecules/") and fname.endswith(".json")
                    and not (fname.endswith("_pos.json") or fname.endswith("_neg.json") or fname.endswith("_esp.json"))
                ):
                    new_fname = f"molecules/{new_key}.json"
                    content = _remap_mol_json_content(content, old_key, new_key)

                # Mesh files: molecules/<old_key>_...
                elif fname.startswith(f"molecules/{old_key}_"):
                    new_fname = fname.replace(f"molecules/{old_key}_", f"molecules/{new_key}_", 1)

                zf_out.writestr(new_fname, content)

    if dry_run:
        return True

    if backup:
        shutil.copy2(zip_path, zip_path + ".bak")

    out_buf.seek(0)
    with open(zip_path, "wb") as f_out:
        f_out.write(out_buf.getvalue())

    return True


def rename_bundle_files(directory, rename_map, dry_run=False, update_readme=True):
    """
    Renames zip and svg files in directory based on rename_map: {old_base: new_base}.
    Also updates references in README.md if requested.
    """
    renamed = []
    for old_base, new_base in rename_map.items():
        old_zip = os.path.join(directory, f"{old_base}.zip")
        new_zip = os.path.join(directory, f"{new_base}.zip")
        old_svg = os.path.join(directory, f"{old_base}.svg")
        new_svg = os.path.join(directory, f"{new_base}.svg")

        if os.path.exists(old_zip):
            print(f"  Rename: {old_base}.zip -> {new_base}.zip")
            if not dry_run:
                os.rename(old_zip, new_zip)
            renamed.append((old_base, new_base))

        if os.path.exists(old_svg):
            print(f"  Rename: {old_base}.svg -> {new_base}.svg")
            if not dry_run:
                os.rename(old_svg, new_svg)

    if update_readme and renamed and not dry_run:
        readme_path = "README.md"
        if os.path.exists(readme_path):
            with open(readme_path, "r", encoding="utf-8") as f:
                content = f.read()
            for old_base, new_base in renamed:
                content = content.replace(f"bundles/molecules/{old_base}.zip", f"bundles/molecules/{new_base}.zip")
                content = content.replace(f"bundles/molecules/{old_base}.svg", f"bundles/molecules/{new_base}.svg")
            
            # Specific cleanups in README names
            content = content.replace("**chloride ion**", "**chloride** or **chloride ion**")
            content = content.replace("**methanoate ion**", "**methanoate** or **formate**")
            content = content.replace(
                "**1,1,1,3,3,3-hexafluoro-2-propanol** or **hexafluoroisopropanol**",
                "**1,1,1,3,3,3-hexafluoropropan-2-ol**, or **1,1,1,3,3,3-hexafluoro-2-propanol**, or **hexafluoroisopropanol**"
            )
            content = content.replace(
                "**methyl bromide** or **bromomethane**",
                "**bromomethane**, or **methyl bromide**"
            )
            content = content.replace("SMILES: CBr||", "SMILES: CBr|")
            with open(readme_path, "w", encoding="utf-8") as f:
                f.write(content)
            print("  Updated references and titles in README.md")

    return renamed


def main():
    parser = argparse.ArgumentParser(description="Check and synchronize molecule bundle naming.")
    parser.add_argument("--dir", default="bundles/molecules", help="Directory containing molecule bundles")
    parser.add_argument("--check", action="store_true", help="Only check and report status (default)")
    parser.add_argument("--fix", action="store_true", help="Fix and synchronize internal bundle contents")
    parser.add_argument("--fix-zip-names", action="store_true", help="Rename zip and svg files for known typos/IUPAC fixes")
    parser.add_argument("--dry-run", action="store_true", help="Simulate actions without writing files")
    parser.add_argument("--backup", action="store_true", help="Create .bak backups before modifying zip files")

    args = parser.parse_args()

    bundles = sorted(glob.glob(os.path.join(args.dir, "*.zip")))
    if not bundles:
        print(f"No zip bundles found in {args.dir}")
        sys.exit(1)

    print(f"Found {len(bundles)} molecule bundles in {args.dir}\n")

    if args.fix_zip_names:
        print("=== Step 1: Renaming zip and svg files for known typos/naming fixes ===")
        renamed = rename_bundle_files(args.dir, RECOMMENDED_ZIP_RENAMES, dry_run=args.dry_run)
        print(f"Renamed {len(renamed)} bundle pairs.\n")
        # Refresh bundle list
        bundles = sorted(glob.glob(os.path.join(args.dir, "*.zip")))

    print("=== Checking bundle synchronization ===")
    synced_count = 0
    desynced_count = 0
    desynced_reports = []

    for b in bundles:
        rep = audit_bundle(b)
        if rep["is_synced"]:
            synced_count += 1
        else:
            desynced_count += 1
            desynced_reports.append(rep)

    print(f"Synchronized bundles: {synced_count} / {len(bundles)}")
    print(f"Desynchronized bundles: {desynced_count} / {len(bundles)}\n")

    if desynced_reports:
        print("Desynchronized bundles details:")
        for r in desynced_reports:
            print(f"  • {r['base_name']}: old_key='{r['old_key']}' ({len(r['issues'])} issues)")
            for iss in r["issues"][:3]:
                print(f"      - {iss}")
            if len(r["issues"]) > 3:
                print(f"      ... and {len(r['issues']) - 3} more issues")

    if args.fix:
        print(f"\n=== Fixing {len(desynced_reports)} desynchronized bundles{' (DRY RUN)' if args.dry_run else ''} ===")
        for r in desynced_reports:
            print(f"  Fixing: {r['base_name']} (migrating '{r['old_key']}' -> '{r['target_key']}')...")
            fix_bundle(r["zip_path"], target_key=r["target_key"], dry_run=args.dry_run, backup=args.backup)

        print("\nFix complete. Re-verifying...")
        recheck_synced = 0
        for b in bundles:
            rep = audit_bundle(b)
            if rep["is_synced"]:
                recheck_synced += 1
        print(f"Post-fix synchronized: {recheck_synced} / {len(bundles)}")

    elif not args.fix_zip_names:
        print("\nTip: Run with --fix to update internal naming, or --fix --dry-run to preview.")


if __name__ == "__main__":
    main()
