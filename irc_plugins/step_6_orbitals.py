"""
Step 6: Interactive Molecular Orbital Studio
Calculates and localizes molecular orbitals (frontier MOs and IBOs) and schedules orbital tracks.
Grid meshing and chemical selection delegated to orbital_cube_mesher.py and orbital_selection.py.
"""

import json
import os
import threading
from typing import Any, Dict, List, Optional

from .base import PipelineStep, StepResult
from .ts_engine import orbital_engine
from .orbital_selection import (
    get_system_charge,
    resolve_orbital_key,
    tag_orbital_classes,
    enrich_timeline_tracks,
)
from .orbital_indexer import get_or_build_catalog
from .orbital_cube_mesher import (
    clear_orbital_caches,
    get_or_compute_cube,
    batch_precompute_cubes,
    get_cached_cubes,
)

ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
TRAJECTORIES_DIR = os.path.join(ROOT_DIR, 'trajectories')
PROLONGED_PATH = os.path.join(TRAJECTORIES_DIR, 'prolonged_pathway.json')

_SCF_CACHE: Dict[int, Dict[str, Any]] = {}
_CACHE_LOCK = threading.Lock()
_resolve_orbital_key = resolve_orbital_key


class Step6Orbitals(PipelineStep):
    """Step 6: Interactive Molecular Orbital Studio with distinct bonding/antibonding palettes."""

    @property
    def step_id(self) -> str: return "step_6_orbitals"

    @property
    def step_name(self) -> str: return "6. Orbitals"

    def _get_megajson(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        s6, s5 = current_state.get(self.step_id, {}), current_state.get("step_5_cascade", {})
        mj = s6.get("megajson") or s5.get("megajson")
        if not mj and os.path.exists(PROLONGED_PATH):
            try: mj = json.load(open(PROLONGED_PATH, "r"))
            except Exception: mj = {}
        if mj and os.path.exists(PROLONGED_PATH) and not mj.get("global_electron_tracks"):
            try:
                disk = json.load(open(PROLONGED_PATH, "r"))
                mj.setdefault("global_electron_tracks", disk.get("global_electron_tracks", {}))
                mj.setdefault("orbital_timeline", disk.get("orbital_timeline", {}))
            except Exception: pass
        return mj or {}

    def get_ui_schema(self, current_state: Dict[str, Any]) -> Dict[str, Any]:
        megajson = self._get_megajson(current_state)
        if not megajson or not megajson.get("frames"):
            return {"type": "orbitals_ui", "has_megajson": False, "message": "No reaction pathway available. Please complete Step 5 first."}

        ep = megajson.get("energy_profile", {}).get("frame_energies_kcal", []) or [f.get("rel_energy_kcal", 0.0) for f in megajson["frames"]]
        timeline = megajson.get("orbital_timeline")
        if not isinstance(timeline, dict):
            timeline = {"default_isovalue": 0.03, "tracks": []}
        else:
            tracks = timeline.get("tracks", [])
            timeline["tracks"] = enrich_timeline_tracks(megajson, tracks if isinstance(tracks, list) else [], _SCF_CACHE, self._ensure_frame_computed)
            timeline.setdefault("default_isovalue", 0.03)

        return {
            "type": "orbitals_ui", "has_megajson": True, "total_frames": len(megajson["frames"]),
            "step_boundaries": megajson.get("step_boundaries", []), "symbols": megajson.get("symbols", []),
            "trajectory": [f.get("coords", f.get("positions", [])) for f in megajson["frames"]],
            "trajectory_symbols": [f.get("symbols") or megajson.get("symbols", []) for f in megajson["frames"]],
            "trajectory_bonds": [f.get("bonds", []) for f in megajson["frames"]],
            "trajectory_charges": megajson.get("trajectory_charges", []), "energy_profile": ep,
            "timeline": timeline, "default_isovalue": 0.03,
            "color_palettes": {
                "bonding": {"pos": "#38bdf8", "neg": "#1d4ed8", "label": "Bonding (Cyan/Blue)"},
                "antibonding": {"pos": "#fb923c", "neg": "#db2777", "label": "Anti-Bonding (Orange/Magenta)"},
                "core": {"pos": "#94a3b8", "neg": "#475569", "label": "Core (Gray)"}
            },
            "message": "Reaction Orbital Studio: explore frontier MOs and localized IBOs across the reaction coordinate."
        }

    def _compute_single_frame(self, fi: int, prev_data: Optional[Dict[str, Any]], megajson: Dict[str, Any]) -> Dict[str, Any]:
        fr = megajson["frames"][fi]
        symbols = fr.get("symbols") or megajson.get("symbols", [])
        charge = get_system_charge(megajson, fi)
        coords = fr.get("coords", fr.get("positions", []))

        mol, mf, mull = orbital_engine.compute_frame_scf(symbols, coords, charge=charge, spin=0)
        same_sys = prev_data is not None and prev_data.get("symbols") == symbols
        if not same_sys:
            ibos_m, star_m, ibos_v, pops, ibo_c = orbital_engine.compute_frame_ibos(mol, mf, symbols)
            f_meta, c_meta, f_vecs = orbital_engine.extract_frontier_orbitals(mol, mf, None, None, ibos_m, star_m, ibos_v)
            n_c = len([k for k, v in ibos_m.items() if v.get("type") == "core"])
        else:
            n_c = prev_data.get("n_core") or len([k for k, v in prev_data.get("ibos_meta", {}).items() if v.get("type") == "core"])
            ref_dict = {
                "mol": prev_data["mol"], "ibo_coeff": prev_data["ibo_coeff"], "n_core": n_c,
                "val_pops": prev_data.get("val_pops") if prev_data.get("val_pops") is not None else prev_data["ibo_pops"][n_c:]
            }
            ibos_m, star_m, ibos_v, pops, ibo_c = orbital_engine.compute_frame_ibos(mol, mf, symbols, ref_dict)
            f_meta, c_meta, f_vecs = orbital_engine.extract_frontier_orbitals(mol, mf, prev_data["mol"], prev_data["vecs"], ibos_m, star_m, ibos_v)

        st_idx = fr.get("step_index")
        for m in (ibos_m, star_m, c_meta):
            for itm in m.values():
                itm["step_index"] = st_idx
        tag_orbital_classes(f_meta, ibos_m, star_m, c_meta)

        frame_dict = {
            "mol": mol, "mf": mf, "meta": f_meta, "ibos_meta": ibos_m,
            "antibonds_meta": star_m, "canonical_pi_meta": c_meta,
            "ibo_pops": pops, "val_pops": pops[n_c:], "ibo_coeff": ibo_c, "vecs": {**ibos_v, **f_vecs},
            "charges": mull, "n_core": n_c, "symbols": symbols, "step_index": st_idx
        }
        _SCF_CACHE[fi] = frame_dict
        return frame_dict

    def _ensure_frame_computed(self, frame_idx: int, megajson: Dict[str, Any]) -> Dict[str, Any]:
        with _CACHE_LOCK:
            if frame_idx in _SCF_CACHE: return _SCF_CACHE[frame_idx]
            from .orbital_indexer import get_frame_groups
            groups = get_frame_groups(megajson)
            g_idx = next((i for i, (gs, ge) in enumerate(groups) if gs <= frame_idx <= ge), 0)
            g_start, g_end = groups[g_idx]
            if frame_idx == g_start: return self._compute_single_frame(g_start, None, megajson)
            grp_cached = [k for k in _SCF_CACHE.keys() if g_start <= k < frame_idx]
            if not grp_cached:
                self._compute_single_frame(g_start, None, megajson)
                grp_cached = [g_start]
            prev_fi = max(grp_cached)
            for fi in range(prev_fi + 1, frame_idx + 1):
                if fi not in _SCF_CACHE: self._compute_single_frame(fi, _SCF_CACHE[fi - 1], megajson)
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
            clear_orbital_caches()
            return {"success": True, "message": "Orbital caches and disk cubes cleared."}

        if action == "get_frame_orbitals":
            data = self._ensure_frame_computed(frame_idx, megajson)
            fr = megajson["frames"][frame_idx]
            catalog = get_or_build_catalog(megajson, _SCF_CACHE, self._ensure_frame_computed)

            def _enrich_pool(metas):
                res = []
                for m in metas:
                    k = m.get("key") or m.get("name")
                    trk = next((t for t in catalog.values() if t.get("frame_map", [None]*n_frames)[frame_idx] == k), None)
                    itm = dict(m)
                    itm["track_id"] = trk["track_id"] if trk else k
                    if trk:
                        itm["label"] = trk["label"]
                        itm["decoupled_label"] = trk["label"]
                    itm["local_key"] = k
                    res.append(itm)
                return res

            f_mos = {}
            for fk, fm in data["meta"].items():
                itm = dict(fm)
                itm["track_id"] = "trk_homo" if "HOMO" in fk else ("trk_lumo" if "LUMO" in fk else fk)
                itm["local_key"] = fk
                f_mos[fk] = itm

            return {
                "success": True, "frame_index": frame_idx, "charges": data["charges"],
                "frontier_orbitals": f_mos, "all_ibos": _enrich_pool(list(data["ibos_meta"].values())),
                "antibonds": _enrich_pool(list(data.get("antibonds_meta", {}).values())),
                "canonical_pi": _enrich_pool(list(data.get("canonical_pi_meta", {}).values())),
                "coords": fr.get("coords", fr.get("positions", [])),
                "symbols": fr.get("symbols") or megajson.get("symbols", []),
                "bonds": fr.get("bonds", []), "step_index": fr.get("step_index")
            }

        if action == "get_density_flux":
            sf, ef = max(0, int(user_input.get("start_frame", 0))), min(n_frames - 1, int(user_input.get("end_frame", n_frames - 1)))
            ds, de = self._ensure_frame_computed(sf, megajson), self._ensure_frame_computed(ef, megajson)
            return {"success": True, "start_frame": sf, "end_frame": ef, "fluxes": orbital_engine.compute_density_flux(ds["ibo_pops"], de["ibo_pops"], megajson["symbols"])}

        if action == "save_timeline":
            if user_input.get("timeline"):
                tl = user_input["timeline"]
                if isinstance(tl, dict) and "tracks" in tl and isinstance(tl["tracks"], list):
                    tl["tracks"] = enrich_timeline_tracks(megajson, tl["tracks"], _SCF_CACHE, self._ensure_frame_computed)
                megajson["orbital_timeline"] = tl
                try:
                    os.makedirs(TRAJECTORIES_DIR, exist_ok=True)
                    with open(PROLONGED_PATH, "w") as f:
                        json.dump(megajson, f, indent=2)
                except Exception:
                    pass
                return {"success": True, "message": "Orbital timeline saved.", "timeline": tl}
            return {"error": "No timeline data provided."}

        if action == "batch_precompute_cubes":
            sf, sz, sel = int(user_input.get("start_frame", 0)), int(user_input.get("batch_size", 2)), user_input.get("selected_tracks", [])
            force = bool(user_input.get("force", False) or user_input.get("overwrite", False))
            return batch_precompute_cubes(megajson, sf, sz, sel, force, self._ensure_frame_computed)

        if action == "get_cached_cubes":
            frames, tracks = [int(x) for x in user_input.get("frames", [])], user_input.get("tracks", [])
            return {"success": True, "cubes": get_cached_cubes(frames, tracks, _SCF_CACHE, megajson)}

        if action in ["get_orbital_mesh", "get_cube"]:
            tid = user_input.get("track_id") or user_input.get("orbital", "HOMO")
            ok = str(user_input.get("orbital", tid))
            ck_raw = user_input.get("chem_key")
            ck = str(ck_raw) if (ck_raw is not None and not isinstance(ck_raw, (int, float))) else None
            try:
                iso = float(user_input.get("isovalue", 0.03) or 0.03)
            except (ValueError, TypeError):
                iso = 0.03
            d = self._ensure_frame_computed(frame_idx, megajson)
            return get_or_compute_cube(frame_idx, ok, ck, iso, d, megajson, track_id=tid)

        return {"error": f"Unknown action: {action}"}

    def run(self, global_state: Dict[str, Any], user_input: Dict[str, Any]) -> StepResult:
        megajson = self._get_megajson(global_state)
        if not megajson: return StepResult(success=False, error_message="Pathway data missing.")
        tl = user_input.get("timeline") or megajson.get("orbital_timeline") or {}
        megajson["orbital_timeline"] = tl
        return StepResult(success=True, data={"completed": True, "megajson": megajson, "orbital_timeline": tl})

    def validate(self, step_data: Dict[str, Any]) -> bool:
        return bool(step_data and step_data.get("completed") and step_data.get("megajson"))
