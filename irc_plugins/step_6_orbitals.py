"""
Step 6: Interactive Molecular Orbital Studio
Calculates and localizes molecular orbitals (HOMO/LUMO frontier MOs and localized IBOs),
providing distinct visual palettes for bonding vs. anti-bonding orbitals,
and scheduling orbital tracks across reaction coordinates.
"""

import json
import os
import re
import threading
from typing import Any, Dict, List, Optional
from rdkit import Chem

from .base import PipelineStep, StepResult
from .ts_engine import orbital_engine

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
OUT_DIR = os.path.join(ROOT_DIR, 'out')
CUBES_DIR = os.path.join(OUT_DIR, 'cubes')
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, 'trajectories')
PROLONGED_PATH = os.path.join(TRAJECTORIES_DIR, 'prolonged_pathway.json')

_SCF_CACHE: Dict[int, Dict[str, Any]] = {}
_GRID_CACHE: Dict[str, Any] = {}
_CUBE_CACHE: Dict[str, str] = {}
_CACHE_LOCK = threading.Lock()


def _get_system_charge(megajson: Dict[str, Any]) -> int:
    """Infer molecular system net formal charge from reaction SMILES or megajson."""
    for ck in ("charge", "net_charge", "formal_charge"):
        if ck in megajson:
            try:
                return int(megajson[ck])
            except (ValueError, TypeError):
                pass
    boundaries = megajson.get("step_boundaries", [])
    if boundaries and ">>" in boundaries[0].get("reaction_smiles", ""):
        reactants = boundaries[0]["reaction_smiles"].split(">>")[0].split(".")
        try:
            return sum(
                sum(a.GetFormalCharge() for a in Chem.MolFromSmiles(r).GetAtoms())
                for r in reactants if Chem.MolFromSmiles(r)
            )
        except Exception:
            return 0
    return 0


def _resolve_orbital_key(d: Dict[str, Any], orb_name: str, chem_key: str = None) -> Optional[str]:
    if not d or "vecs" not in d:
        return None
    if orb_name in d.get("vecs", {}):
        return orb_name
    all_metas = {
        **d.get("ibos_meta", {}),
        **d.get("canonical_pi_meta", {}),
        **d.get("antibonds_meta", {}),
        **d.get("meta", {})
    }
    for target in [orb_name, chem_key]:
        if not target:
            continue
        for k, meta in all_metas.items():
            if (meta.get("key") == target or meta.get("name") == target or meta.get("chem_key") == target or meta.get("label") == target) and k in d.get("vecs", {}):
                return k
def _clean_role_label(meta: Optional[Dict[str, Any]]) -> Optional[str]:
    if not meta:
        return None
    typ = meta.get("type", "")
    atoms = meta.get("atoms", [])
    atoms_str = "-".join(atoms) if atoms else ""
    if typ == "lone_pair":
        return f"LP({atoms[0]})" if atoms else "LP"
    elif "sigma_star" in typ or "sigmastar" in typ:
        return f"σ*({atoms_str})" if atoms_str else "σ*"
    elif "pistar" in typ or "pi_star" in typ or typ == "delocalized_pistar":
        return f"π*({atoms_str})" if atoms_str else "π*"
    elif "pi" in typ or typ == "delocalized_pi":
        return f"π({atoms_str})" if atoms_str else "π"
    elif "sigma" in typ or typ == "delocalized_sigma":
        return f"σ({atoms_str})" if atoms_str else "σ"
    return meta.get("label") or typ


def enrich_timeline_tracks(megajson: Dict[str, Any], tracks: list, s6_instance=None) -> list:
    """Enrich orbital timeline tracks with normal chemical nomenclature and role transitions."""
    if not tracks:
        return tracks

    n_frames = len(megajson.get("frames", []))
    f0 = _SCF_CACHE.get(0) if n_frames > 0 else None
    flast = _SCF_CACHE.get(n_frames - 1) if n_frames > 0 else None

    for trk in tracks:
        if not isinstance(trk, dict):
            continue
        ok = trk.get("orbital_key", "")
        if ok == "track_esp_charge" or trk.get("id") == "trk_esp_charge":
            trk["label"] = "ESP Surface"
            trk["title"] = "ESP Surface"
            trk["short_label"] = "ESP"
            trk["role_transition"] = "ESP"
            trk["category"] = "ESP"
            continue

        start_role = trk.get("start_role")
        end_role = trk.get("end_role")

        if (not start_role or not end_role) and f0 is not None and flast is not None:
            r0 = _resolve_orbital_key(f0, ok, trk.get("chem_key"))
            m0 = {**f0.get("ibos_meta", {}), **f0.get("canonical_pi_meta", {}), **f0.get("antibonds_meta", {})}.get(r0, {}) if r0 else {}
            rlast = _resolve_orbital_key(flast, ok, trk.get("chem_key"))
            mlast = {**flast.get("ibos_meta", {}), **flast.get("canonical_pi_meta", {}), **flast.get("antibonds_meta", {})}.get(rlast, {}) if rlast else {}
            start_role = _clean_role_label(m0)
            end_role = _clean_role_label(mlast)

        ck = trk.get("chem_key", "")
        if not start_role and ck:
            if ck.startswith("lp_"):
                start_role = f"LP({ck.replace('lp_', '')})"
            elif ck.startswith("sigma_"):
                start_role = f"σ({ck.replace('sigma_', '')})"
            elif ck.startswith("pi_"):
                start_role = f"π({ck.replace('pi_', '')})"
            elif "pistar" in ck:
                start_role = f"π*({ck.split('_')[1] if '_' in ck else ''})"

        if start_role and end_role and start_role != end_role:
            label = f"{start_role} → {end_role}"
            role_trans = f"{start_role.split('(')[0]} → {end_role.split('(')[0]}"
        else:
            label = start_role or end_role or trk.get("title") or ok
            label = re.sub(r"^Non-Bonding\s+", "", label)
            label = re.sub(r"\[Antibond\]", "", label).strip()
            role_trans = label.split("(")[0] if "(" in label else label

        trk["label"] = label
        trk["title"] = label
        trk["role_transition"] = role_trans
        trk["short_label"] = role_trans
        if start_role:
            trk["start_role"] = start_role
        if end_role:
            trk["end_role"] = end_role
        trk["category"] = role_trans

    return tracks


class Step6Orbitals(PipelineStep):
    """Step 6: Interactive Molecular Orbital Studio with distinct bonding/antibonding palettes."""

    @property
    def step_id(self) -> str:
        return "step_6_orbitals"

    @property
    def step_name(self) -> str:
        return "6. Orbitals"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        """Fetch pathway megajson, which is guaranteed to be generated by Step 5."""
        s6 = current_state.get(self.step_id, {})
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
                "type": "orbitals_ui",
                "has_megajson": False,
                "message": "No reaction pathway available. Please complete Step 5 first."
            }

        ep = megajson.get("energy_profile", {}).get("frame_energies_kcal", []) or [
            f.get("rel_energy_kcal", 0.0) for f in megajson["frames"]
        ]
        timeline = megajson.get("orbital_timeline")
        if not isinstance(timeline, dict):
            timeline = {
                "default_isovalue": 0.03,
                "tracks": []
            }
        else:
            if "tracks" not in timeline or not isinstance(timeline["tracks"], list):
                timeline["tracks"] = []
            else:
                timeline["tracks"] = enrich_timeline_tracks(megajson, timeline["tracks"], self)
            if "default_isovalue" not in timeline:
                timeline["default_isovalue"] = 0.03

        return {
            "type": "orbitals_ui",
            "has_megajson": True,
            "total_frames": len(megajson["frames"]),
            "step_boundaries": megajson.get("step_boundaries", []),
            "symbols": megajson.get("symbols", []),
            "trajectory": [f.get("coords", f.get("positions", [])) for f in megajson["frames"]],
            "trajectory_charges": megajson.get("trajectory_charges", []),
            "energy_profile": ep,
            "timeline": timeline,
            "default_isovalue": 0.03,
            "color_palettes": {
                "bonding": {"pos": "#38bdf8", "neg": "#1d4ed8", "label": "Bonding (Cyan/Blue)"},
                "antibonding": {"pos": "#fb923c", "neg": "#db2777", "label": "Anti-Bonding (Orange/Magenta)"},
                "core": {"pos": "#94a3b8", "neg": "#475569", "label": "Core (Gray)"}
            },
            "message": "Reaction Orbital Studio: explore frontier MOs and localized IBOs across the reaction coordinate."
        }

    def _get_grid(self, megajson: Dict[str, Any]):
        if "lo" not in _GRID_CACHE:
            all_pos = [f.get("coords", f.get("positions", [])) for f in megajson["frames"]]
            _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"] = orbital_engine.build_global_grid(
                all_pos, grid_points=32, margin_angstrom=2.5
            )
        return _GRID_CACHE["lo"], _GRID_CACHE["step"], _GRID_CACHE["npts"]

    def _tag_orbital_classes(self, f_meta: Dict[str, Any], ibos_m: Dict[str, Any], star_m: Dict[str, Any], c_meta: Dict[str, Any]):
        """Tag every orbital with its chemical nature and distinct color palette."""
        # 1. Frontier MOs
        for k, m in f_meta.items():
            name = m.get("name", "").upper()
            if "HOMO" in name:
                m["orbital_class"] = "bonding"
                m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}
            elif "LUMO" in name:
                m["orbital_class"] = "antibonding"
                m["palette"] = {"pos": "#fb923c", "neg": "#db2777"}
            else:
                m["orbital_class"] = "bonding"
                m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}

        # 2. Localized IBOs
        for k, m in ibos_m.items():
            t = m.get("type", "")
            if t == "core":
                m["orbital_class"] = "core"
                m["palette"] = {"pos": "#94a3b8", "neg": "#475569"}
            else:
                m["orbital_class"] = "bonding"
                m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}

        # 3. Antibonds
        for k, m in star_m.items():
            m["orbital_class"] = "antibonding"
            m["palette"] = {"pos": "#fb923c", "neg": "#db2777"}

        # 4. Canonical pi
        for k, m in c_meta.items():
            m["orbital_class"] = "bonding"
            m["palette"] = {"pos": "#38bdf8", "neg": "#1d4ed8"}

    def _compute_single_frame(self, fi: int, prev_data: Optional[Dict[str, Any]], megajson: Dict[str, Any]) -> Dict[str, Any]:
        symbols = megajson["symbols"]
        charge = _get_system_charge(megajson)
        coords = megajson["frames"][fi].get("coords", megajson["frames"][fi].get("positions", []))

        mol, mf, mull = orbital_engine.compute_frame_scf(symbols, coords, charge=charge, spin=0)
        if prev_data is None:
            ibos_m, star_m, ibos_v, pops, ibo_c = orbital_engine.compute_frame_ibos(mol, mf, symbols)
            f_meta, c_meta, f_vecs = orbital_engine.extract_frontier_orbitals(mol, mf, None, None, ibos_m, star_m, ibos_v)
            n_c = len([k for k, v in ibos_m.items() if v.get("type") == "core"])
        else:
            n_c = prev_data.get("n_core")
            if n_c is None:
                n_c = len([k for k, v in prev_data.get("ibos_meta", {}).items() if v.get("type") == "core"])
            ref_dict = {
                "mol": prev_data["mol"],
                "ibo_coeff": prev_data["ibo_coeff"],
                "n_core": n_c,
                "val_pops": prev_data.get("val_pops") if prev_data.get("val_pops") is not None else prev_data["ibo_pops"][n_c:]
            }
            ibos_m, star_m, ibos_v, pops, ibo_c = orbital_engine.compute_frame_ibos(mol, mf, symbols, ref_dict)
            f_meta, c_meta, f_vecs = orbital_engine.extract_frontier_orbitals(mol, mf, prev_data["mol"], prev_data["vecs"], ibos_m, star_m, ibos_v)

        self._tag_orbital_classes(f_meta, ibos_m, star_m, c_meta)

        frame_dict = {
            "mol": mol, "mf": mf, "meta": f_meta, "ibos_meta": ibos_m,
            "antibonds_meta": star_m, "canonical_pi_meta": c_meta,
            "ibo_pops": pops, "val_pops": pops[n_c:], "ibo_coeff": ibo_c, "vecs": {**ibos_v, **f_vecs},
            "charges": mull, "n_core": n_c
        }
        _SCF_CACHE[fi] = frame_dict
        return frame_dict

    def _ensure_frame_computed(self, frame_idx: int, megajson: Dict[str, Any]) -> Dict[str, Any]:
        with _CACHE_LOCK:
            if frame_idx in _SCF_CACHE:
                return _SCF_CACHE[frame_idx]
            if 0 not in _SCF_CACHE:
                self._compute_single_frame(0, None, megajson)
            if frame_idx == 0:
                return _SCF_CACHE[0]

            prev_fi = max([k for k in _SCF_CACHE.keys() if k < frame_idx], default=0)
            for fi in range(prev_fi + 1, frame_idx + 1):
                if fi not in _SCF_CACHE:
                    self._compute_single_frame(fi, _SCF_CACHE[fi - 1], megajson)
            return _SCF_CACHE[frame_idx]

    def custom_action(self, user_input: Dict[str, Any], global_state: Dict[str, Any] = None) -> Dict[str, Any]:
        megajson = self._get_megajson(global_state or {})
        if not megajson:
            return {"error": "Pathway data is missing. Please complete Step 5."}

        action = user_input.get("action_type")
        frame_idx = int(user_input.get("frame_index", 0))
        n_frames = len(megajson.get("frames", []))
        if frame_idx < 0 or frame_idx >= n_frames:
            return {"error": f"Frame {frame_idx} out of range (0-{n_frames-1})"}

        if action == "clear_cache":
            with _CACHE_LOCK:
                _SCF_CACHE.clear()
                _CUBE_CACHE.clear()
                _GRID_CACHE.clear()
                if os.path.exists(CUBES_DIR):
                    for cf in os.listdir(CUBES_DIR):
                        if cf.endswith(".cube"):
                            try:
                                os.remove(os.path.join(CUBES_DIR, cf))
                            except Exception:
                                pass
            return {"success": True, "message": "Orbital caches and disk cubes cleared."}

        if action == "get_frame_orbitals":
            data = self._ensure_frame_computed(frame_idx, megajson)
            coords = megajson["frames"][frame_idx].get("coords", megajson["frames"][frame_idx].get("positions", []))
            return {
                "success": True,
                "frame_index": frame_idx,
                "charges": data["charges"],
                "frontier_orbitals": data["meta"],
                "all_ibos": list(data["ibos_meta"].values()),
                "antibonds": list(data.get("antibonds_meta", {}).values()),
                "canonical_pi": list(data.get("canonical_pi_meta", {}).values()),
                "coords": coords
            }

        if action == "get_density_flux":
            sf = max(0, int(user_input.get("start_frame", 0)))
            ef = min(n_frames - 1, int(user_input.get("end_frame", n_frames - 1)))
            ds = self._ensure_frame_computed(sf, megajson)
            de = self._ensure_frame_computed(ef, megajson)
            return {
                "success": True,
                "start_frame": sf,
                "end_frame": ef,
                "fluxes": orbital_engine.compute_density_flux(ds["ibo_pops"], de["ibo_pops"], megajson["symbols"])
            }

        if action == "save_timeline":
            if user_input.get("timeline"):
                tl = user_input["timeline"]
                if isinstance(tl, dict) and "tracks" in tl and isinstance(tl["tracks"], list):
                    tl["tracks"] = enrich_timeline_tracks(megajson, tl["tracks"], self)
                megajson["orbital_timeline"] = tl
                # Persist to prolonged_pathway.json
                try:
                    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
                    with open(PROLONGED_PATH, "w") as f:
                        json.dump(megajson, f, indent=2)
                except Exception:
                    pass
                return {"success": True, "message": "Orbital timeline saved.", "timeline": tl}
            return {"error": "No timeline data provided."}

        if action == "batch_precompute_cubes":
            sf = int(user_input.get("start_frame", 0))
            sz = int(user_input.get("batch_size", 2))
            sel = user_input.get("selected_tracks", [])
            force = bool(user_input.get("force", False) or user_input.get("overwrite", False))
            ef = min(n_frames - 1, sf + sz - 1)
            lo, step, npts = self._get_grid(megajson)
            batch_cubes = {}

            os.makedirs(CUBES_DIR, exist_ok=True)
            for fi in range(sf, ef + 1):
                d = self._ensure_frame_computed(fi, megajson)
                coords = megajson["frames"][fi].get("coords", megajson["frames"][fi].get("positions", []))
                for trk in sel:
                    ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                    ck = trk.get("chem_key") if isinstance(trk, dict) else None
                    if not ok or ok == "track_esp_charge":
                        continue
                    cache_key = f"{fi}_{ok}"

                    rk = _resolve_orbital_key(d, ok, ck)
                    if not rk or rk not in d.get("vecs", {}):
                        # Orbital / bond is not formed at this frame
                        continue

                    disk_cube_path = os.path.join(CUBES_DIR, f"{fi}_{rk}.cube")

                    if not force:
                        if cache_key in _CUBE_CACHE:
                            cached_val = _CUBE_CACHE[cache_key]
                            if cached_val.startswith(f"F{fi}_{rk}"):
                                batch_cubes[cache_key] = cached_val
                                continue

                        if os.path.exists(disk_cube_path):
                            try:
                                with open(disk_cube_path, "r") as cf:
                                    first_line = cf.readline().strip()
                                    if first_line == f"F{fi}_{rk}":
                                        cf.seek(0)
                                        cube_text = cf.read()
                                        _CUBE_CACHE[cache_key] = cube_text
                                        batch_cubes[cache_key] = cube_text
                                        continue
                            except Exception:
                                pass

                    v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
                    cube_text = orbital_engine.generate_cube_text(
                        f"F{fi}_{rk}", megajson["symbols"], coords, v, lo, step, npts
                    )
                    _CUBE_CACHE[cache_key] = cube_text
                    batch_cubes[cache_key] = cube_text
                    try:
                        with open(disk_cube_path, "w") as cf:
                            cf.write(cube_text)
                    except Exception:
                        pass

            return {
                "success": True,
                "start_frame": sf,
                "end_frame": ef,
                "total_frames": n_frames,
                "percent": int(round((ef + 1) / n_frames * 100)),
                "is_complete": bool(ef >= n_frames - 1),
                "cubes": batch_cubes
            }

        if action == "get_cached_cubes":
            frames = [int(x) for x in user_input.get("frames", [])]
            tracks = user_input.get("tracks", [])
            found = {}
            for fi in frames:
                d = _SCF_CACHE.get(fi)
                for trk in tracks:
                    ok = trk.get("orbital_key") if isinstance(trk, dict) else str(trk)
                    ck = trk.get("chem_key") if isinstance(trk, dict) else None
                    if not ok or ok == "track_esp_charge":
                        continue
                    cache_key = f"{fi}_{ok}"
                    if cache_key in _CUBE_CACHE:
                        found[cache_key] = _CUBE_CACHE[cache_key]
                        continue
                    rk = _resolve_orbital_key(d, ok, ck) if d else ok
                    if rk:
                        disk_cube_path = os.path.join(CUBES_DIR, f"{fi}_{rk}.cube")
                        if os.path.exists(disk_cube_path):
                            try:
                                with open(disk_cube_path, "r") as cf:
                                    first_line = cf.readline().strip()
                                    if first_line == f"F{fi}_{rk}":
                                        cf.seek(0)
                                        cube_text = cf.read()
                                        _CUBE_CACHE[cache_key] = cube_text
                                        found[cache_key] = cube_text
                            except Exception:
                                pass
            return {"success": True, "cubes": found}

        if action in ["get_orbital_mesh", "get_cube"]:
            ok = str(user_input.get("orbital", "HOMO"))
            ck_raw = user_input.get("chem_key")
            ck = str(ck_raw) if (ck_raw is not None and not isinstance(ck_raw, (int, float))) else None
            iso_raw = user_input.get("isovalue")
            try:
                iso = float(iso_raw) if (iso_raw is not None and iso_raw != "") else 0.03
            except (ValueError, TypeError):
                iso = 0.03
            coords = megajson["frames"][frame_idx].get("coords", megajson["frames"][frame_idx].get("positions", []))

            cache_key = f"{frame_idx}_{ok}"

            d = self._ensure_frame_computed(frame_idx, megajson)
            rk = _resolve_orbital_key(d, ok, ck)
            if not rk or rk not in d.get("vecs", {}):
                return {"success": False, "available": False, "error": f"Orbital '{ck or ok}' not formed at frame {frame_idx}"}

            disk_cube_path = os.path.join(CUBES_DIR, f"{frame_idx}_{rk}.cube")

            # Check memory cache
            if cache_key in _CUBE_CACHE:
                cube_text = _CUBE_CACHE[cache_key]
                if cube_text.startswith(f"F{frame_idx}_{rk}"):
                    return {
                        "success": True,
                        "cached": True,
                        "frame_index": frame_idx,
                        "orbital": ok,
                        "resolved_orbital": rk,
                        "chem_key": ck,
                        "isovalue": iso,
                        "cube": cube_text,
                        "filename": f"{ok}_f{frame_idx}.cube",
                        "coords": coords
                    }

            # Check disk cache with header validation
            if os.path.exists(disk_cube_path):
                try:
                    with open(disk_cube_path, "r") as cf:
                        first_line = cf.readline().strip()
                        if first_line == f"F{frame_idx}_{rk}":
                            cf.seek(0)
                            cube_text = cf.read()
                            _CUBE_CACHE[cache_key] = cube_text
                            return {
                                "success": True,
                                "cached": True,
                                "frame_index": frame_idx,
                                "orbital": ok,
                                "resolved_orbital": rk,
                                "chem_key": ck,
                                "isovalue": iso,
                                "cube": cube_text,
                                "filename": f"{ok}_f{frame_idx}.cube",
                                "coords": coords
                            }
                except Exception:
                    pass

            lo, step, npts = self._get_grid(megajson)
            v = orbital_engine.evaluate_mo_grid(d["mol"], d["vecs"][rk], lo, step, npts)
            cube_text = orbital_engine.generate_cube_text(
                f"F{frame_idx}_{rk}", megajson["symbols"], coords, v, lo, step, npts
            )
            _CUBE_CACHE[cache_key] = cube_text

            # Persist to disk
            try:
                os.makedirs(CUBES_DIR, exist_ok=True)
                with open(disk_cube_path, "w") as cf:
                    cf.write(cube_text)
            except Exception:
                pass

            return {
                "success": True,
                "cached": False,
                "frame_index": frame_idx,
                "orbital": ok,
                "resolved_orbital": rk,
                "chem_key": ck,
                "isovalue": iso,
                "cube": cube_text,
                "filename": f"{ok}_f{frame_idx}.cube",
                "coords": coords
            }

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson:
            return StepResult(success=False, error_message="Pathway data missing.")
        tl = user_input.get("timeline") or megajson.get("orbital_timeline") or {}
        megajson["orbital_timeline"] = tl
        return StepResult(success=True, data={"completed": True, "megajson": megajson, "orbital_timeline": tl})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed") and step_data.get("megajson"))
