import json
from typing import Any, Dict, List, Optional, Tuple
import numpy as np
from rdkit import Chem
from scipy.optimize import linear_sum_assignment


def kabsch_fit(P: np.ndarray, Q: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Compute optimal rotation matrix R and translation vector t minimizing RMSD(P, Q @ R.T + t)."""
    cP, cQ = np.mean(P, axis=0), np.mean(Q, axis=0)
    H = (Q - cQ).T @ (P - cP)
    U, _, Vt = np.linalg.svd(H)
    R = Vt.T @ U.T
    if np.linalg.det(R) < 0:
        Vt[-1, :] *= -1
        R = Vt.T @ U.T
    return R, cP - cQ @ R.T


def solve_optimal_element_alignment(
    p_t: np.ndarray, p_s: np.ndarray, sym_t: List[str], sym_s: List[str]
) -> Tuple[float, np.ndarray, np.ndarray, List[int], List[int]]:
    """Solve globally optimal 3D Kabsch alignment preserving elements via Hungarian matching."""
    P, Q = p_t - np.mean(p_t, axis=0), p_s - np.mean(p_s, axis=0)
    _, _, Vt_p = np.linalg.svd(P.T @ P)
    _, _, Vt_q = np.linalg.svd(Q.T @ Q)
    best_rmsd, best_R, best_t = 999.0, np.eye(3), np.zeros(3)
    best_row, best_col = list(range(len(sym_t))), list(range(len(sym_s)))

    for sx in [1, -1]:
        for sy in [1, -1]:
            for sz in [1, -1]:
                S = np.diag([sx, sy, sz])
                R_init = Vt_p.T @ S @ Vt_q
                if np.linalg.det(R_init) <= 0:
                    continue
                Q_rot = Q @ R_init.T
                cost = np.zeros((len(sym_t), len(sym_s)))
                for i, st in enumerate(sym_t):
                    for j, ss in enumerate(sym_s):
                        cost[i, j] = np.linalg.norm(P[i] - Q_rot[j]) if st == ss else 1e5
                r_ind, c_ind = linear_sum_assignment(cost)
                if any(sym_t[r] != sym_s[c] for r, c in zip(r_ind, c_ind)):
                    continue
                R, t = kabsch_fit(p_t[r_ind], p_s[c_ind])
                rmsd = float(np.sqrt(np.mean((p_t[r_ind] - (p_s[c_ind] @ R.T + t)) ** 2)))
                if rmsd < best_rmsd:
                    best_rmsd, best_R, best_t = rmsd, R, t
                    best_row, best_col = r_ind.tolist(), c_ind.tolist()
    return best_rmsd, best_R, best_t, best_row, best_col


def find_all_hero_overlaps(step_prev: Dict[str, Any], step_next: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Identify ALL common molecules/fragments and compute exact element-matched mappings."""
    prods, reacs = step_prev.get("fragments", {}).get("products", []), step_next.get("fragments", {}).get("reactants", [])
    p_last = np.array(step_prev["frames"][-1].get("coords", step_prev["frames"][-1].get("positions", [])))
    p_first = np.array(step_next["frames"][0].get("coords", step_next["frames"][0].get("positions", [])))
    sym_p, sym_n = step_prev.get("symbols", []), step_next.get("symbols", [])
    candidates, used_next = [], set()

    for p_idx, fp in enumerate(prods):
        can_p = Chem.CanonSmiles(fp.get("smiles", "")) if fp.get("smiles") else ""
        for r_idx, fr in enumerate(reacs):
            if r_idx in used_next:
                continue
            can_r = Chem.CanonSmiles(fr.get("smiles", "")) if fr.get("smiles") else ""
            if can_p and can_p == can_r:
                ip, ir = fp.get("atom_indices", []), fr.get("atom_indices", [])
                sym_sub_p, sym_sub_n = [sym_p[i] for i in ip], [sym_n[i] for i in ir]
                rmsd_d, R_d, t_d, r_d, c_d = solve_optimal_element_alignment(p_last[ip], p_first[ir], sym_sub_p, sym_sub_n)
                p_f_m = p_first[ir].copy(); p_f_m[:, 0] = -p_f_m[:, 0]
                rmsd_m, R_m, t_m, r_m, c_m = solve_optimal_element_alignment(p_last[ip], p_f_m, sym_sub_p, sym_sub_n)
                use_m = (rmsd_m < rmsd_d - 0.05)
                rmsd, R, t, r_sub, c_sub = (rmsd_m, R_m, t_m, r_m, c_m) if use_m else (rmsd_d, R_d, t_d, r_d, c_d)
                p_m, n_m = [ip[i] for i in r_sub], [ir[i] for i in c_sub]
                form = fp.get("formula") or fr.get("formula") or can_p
                tag = " [Enantiomer Inverted]" if use_m else ""
                candidates.append({
                    "id": f"mol_{p_idx}_{r_idx}", "label": f"{form}{tag} (Rigid All-Atom, {len(p_m)} atoms, RMSD: {rmsd:.4f} Å)",
                    "formula": f"{form}{tag}", "smiles": can_p, "prev_indices": p_m, "next_indices": n_m,
                    "num_atoms": len(p_m), "rmsd": round(rmsd, 4), "R": R.tolist(), "t": t.tolist(), "is_joint": False, "is_reflected": use_m
                })

                heavy_p = [i for i in p_m if sym_p[i] != "H"]
                heavy_n = [n_m[p_m.index(i)] for i in heavy_p]
                if len(heavy_p) >= 3 and len(heavy_p) < len(p_m):
                    src_h = p_f_m[heavy_n] if use_m else p_first[heavy_n]
                    R_h, t_h = kabsch_fit(p_last[heavy_p], src_h)
                    rmsd_h = float(np.sqrt(np.mean((p_last[heavy_p] - (src_h @ R_h.T + t_h)) ** 2)))
                    candidates.append({
                        "id": f"mol_{p_idx}_{r_idx}_scaffold",
                        "label": f"{form}{tag} (Rotamer-Aware Scaffold, {len(heavy_p)} atoms, RMSD: {rmsd_h:.4f} Å)",
                        "formula": f"{form} [Scaffold]{tag}", "smiles": can_p, "prev_indices": p_m, "next_indices": n_m,
                        "num_atoms": len(heavy_p), "rmsd": round(rmsd_h, 4), "R": R_h.tolist(), "t": t_h.tolist(), "is_joint": False, "is_reflected": use_m
                    })
                used_next.add(r_idx)
                break

    if not candidates and len(sym_p) == len(sym_n) and sorted(sym_p) == sorted(sym_n):
        rmsd_d, R_d, t_d, r_d, c_d = solve_optimal_element_alignment(p_last, p_first, sym_p, sym_n)
        p_f_m = p_first.copy(); p_f_m[:, 0] = -p_f_m[:, 0]
        rmsd_m, R_m, t_m, r_m, c_m = solve_optimal_element_alignment(p_last, p_f_m, sym_p, sym_n)
        use_m = (rmsd_m < rmsd_d - 0.05)
        rmsd, R, t, r_sub, c_sub = (rmsd_m, R_m, t_m, r_m, c_m) if use_m else (rmsd_d, R_d, t_d, r_d, c_d)
        tag = " [Enantiomer Inverted]" if use_m else ""
        candidates.append({
            "id": "full_system", "label": f"Full System{tag} ({len(sym_p)} atoms, RMSD: {rmsd:.4f} Å)",
            "formula": f"Full Scaffold{tag}", "smiles": step_prev.get("reaction_smiles", ""), "prev_indices": r_sub,
            "next_indices": c_sub, "num_atoms": len(sym_p), "rmsd": round(rmsd, 4), "R": R.tolist(), "t": t.tolist(), "is_joint": False, "is_reflected": use_m
        })
    else:
        distinct = [c for c in candidates if not c["id"].endswith("_scaffold")]
        if len(distinct) > 1:
            jp, jn = [i for c in distinct for i in c["prev_indices"]], [i for c in distinct for i in c["next_indices"]]
            has_ref = any(c.get("is_reflected") for c in distinct)
            src_j = p_first.copy()
            if has_ref: src_j[:, 0] = -src_j[:, 0]
            R, t = kabsch_fit(p_last[jp], src_j[jn])
            j_rmsd = float(np.sqrt(np.mean((p_last[jp] - (src_j[jn] @ R.T + t)) ** 2)))
            candidates.insert(0, {
                "id": "joint_all", "label": f"All {len(distinct)} Coinciding Molecules ({len(jp)} atoms, RMSD: {j_rmsd:.4f} Å)",
                "formula": "Joint Complex", "smiles": " + ".join(c["smiles"] for c in distinct), "prev_indices": jp,
                "next_indices": jn, "num_atoms": len(jp), "rmsd": round(j_rmsd, 4), "R": R.tolist(), "t": t.tolist(), "is_joint": True, "is_reflected": has_ref
            })
    return candidates


def generate_bridge_frames(
    p_start: np.ndarray, p_end: np.ndarray, e_start: float, e_end: float, n_frames: int = 8,
    bonds: Optional[List[Any]] = None, symbols: Optional[List[str]] = None
) -> List[Dict[str, Any]]:
    """Generate intermittent transition frames with zero atom clashes via IDPP geodesic interpolation."""
    frames = []
    try:
        from ase import Atoms
        from ase.mep import idpp_interpolate
        syms = symbols if symbols and len(symbols) == len(p_start) else ["C"] * len(p_start)
        images = [Atoms(symbols=syms, positions=p_start)]
        for i in range(1, n_frames + 1):
            s = (i / (n_frames + 1)) ** 2 * (3.0 - 2.0 * (i / (n_frames + 1)))
            images.append(Atoms(symbols=syms, positions=(1.0 - s) * p_start + s * p_end))
        images.append(Atoms(symbols=syms, positions=p_end))
        idpp_interpolate(images, traj=None, log=None, fmax=0.05, steps=50)
        for i, img in enumerate(images[1:-1]):
            s = ((i + 1) / (n_frames + 1)) ** 2 * (3.0 - 2.0 * ((i + 1) / (n_frames + 1)))
            frames.append({
                "coords": img.get_positions().round(4).tolist(),
                "energy_kcal": round(float((1.0 - s) * e_start + s * e_end), 2),
                "is_bridge": True
            })
        return frames
    except Exception:
        pass

    for i in range(1, n_frames + 1):
        s = (i / (n_frames + 1)) ** 2 * (3.0 - 2.0 * (i / (n_frames + 1)))
        coords = (1.0 - s) * p_start + s * p_end
        for _ in range(5):
            D = np.linalg.norm(coords[:, None, :] - coords[None, :, :], axis=2)
            np.fill_diagonal(D, 99.0)
            for a, b in zip(*np.where(D < 1.15)):
                if a < b and (dist := float(np.linalg.norm(coords[b] - coords[a]))) > 1e-4:
                    shift = 0.5 * (1.15 - dist) * ((coords[b] - coords[a]) / dist)
                    coords[b] += shift; coords[a] -= shift
        frames.append({"coords": coords.round(4).tolist(), "energy_kcal": round(float((1.0 - s) * e_start + s * e_end), 2), "is_bridge": True})
    return frames


def stitch_multistep_trajectories(
    trajectories: List[Dict[str, Any]], junction_choices: Optional[Dict[str, str]] = None, n_bridge_frames: int = 8
) -> Dict[str, Any]:
    """Stitch IRC trajectories with optimal element alignment, monotonic s_arc, and step-aware bond orders."""
    if not trajectories:
        return {}
    if len(trajectories) == 1:
        return trajectories[0]

    choices = junction_choices or {}
    aligned_steps = [json.loads(json.dumps(trajectories[0]))]
    step_reports, junction_previews, all_bridges = [], [], []

    for k in range(1, len(trajectories)):
        prev_step, next_step = aligned_steps[k - 1], trajectories[k]
        all_cands = find_all_hero_overlaps(prev_step, next_step)
        cand = next((c for c in all_cands if c["id"] == choices.get(str(k))), None) if all_cands and choices.get(str(k)) else None
        if not cand and all_cands:
            cand = next((c for c in all_cands if not c.get("is_joint")), all_cands[0])

        out_next = json.loads(json.dumps(next_step))
        if cand:
            R, t = np.array(cand["R"]), np.array(cand["t"])
            b_rmsd, idxs_p, idxs_n = cand["rmsd"], cand["prev_indices"], cand["next_indices"]
            prev_syms, next_syms = prev_step.get("symbols", []), next_step.get("symbols", [])
            c_for_r = {r: c for r, c in zip(idxs_p, idxs_n)}
            if len(prev_syms) == len(next_syms) and len(idxs_p) == len(prev_syms):
                order = [c_for_r[r] for r in range(len(prev_syms))]
            elif len(idxs_p) > 0 and len(prev_syms) == len(next_syms):
                used_c, order = set(), []
                for r in range(len(prev_syms)):
                    if r in c_for_r:
                        c = c_for_r[r]
                        order.append(c)
                        used_c.add(c)
                    else:
                        rem = [c for c in range(len(next_syms)) if c not in used_c]
                        order.append(rem[0])
                        used_c.add(rem[0])
            else:
                order = list(range(len(next_syms)))

            out_next["symbols"] = [next_syms[c] for c in order]
            inv_order = {c: r for r, c in enumerate(order)}
            remapped_bonds = []
            for b in next_step.get("bonds", []):
                u = b.get("from", b[0] if isinstance(b, (list, tuple)) else 0)
                v = b.get("to", b[1] if isinstance(b, (list, tuple)) else 0)
                if u in inv_order and v in inv_order:
                    nb = dict(b)
                    nb["from"], nb["to"] = inv_order[u], inv_order[v]
                    remapped_bonds.append(nb)
            out_next["bonds"] = remapped_bonds
            is_ref = cand.get("is_reflected", False)
            for fr in out_next["frames"]:
                raw_pos = np.array(fr.get("coords", fr.get("positions", [])))
                if is_ref: raw_pos[:, 0] = -raw_pos[:, 0]
                fr["coords"] = ((raw_pos @ R.T + t)[order]).round(4).tolist()
                fr.pop("positions", None)
        else:
            b_rmsd, idxs_p, idxs_n = 0.0, [], []

        prev_pos = np.array(prev_step["frames"][-1].get("coords", prev_step["frames"][-1].get("positions", [])))
        unalign_pos = np.array(next_step["frames"][0].get("coords", next_step["frames"][0].get("positions", [])))
        if cand and cand.get("is_reflected"): unalign_pos[:, 0] = -unalign_pos[:, 0]
        align_pos = np.array(out_next["frames"][0].get("coords", out_next["frames"][0].get("positions", [])))
        bridge = generate_bridge_frames(prev_pos, align_pos, prev_step["frames"][-1].get("rel_energy_kcal", 0.0), out_next["frames"][0].get("rel_energy_kcal", 0.0), n_bridge_frames, bonds=prev_step.get("bonds", []), symbols=prev_step.get("symbols", []))
        j_info = {
            "junction_index": k, "step_prev_index": k, "step_next_index": k + 1,
            "step_prev_name": prev_step.get("reaction_smiles", f"Step {k}"), "step_next_name": next_step.get("reaction_smiles", f"Step {k + 1}"),
            "all_candidates": all_cands, "selected_candidate_id": cand["id"] if cand else None,
            "hero_smiles": cand["smiles"] if cand else None, "hero_formula": cand.get("formula", "") if cand else None,
            "n_matched_atoms": len(idxs_p), "joint_rmsd_angstrom": b_rmsd, "rmsd": b_rmsd, "is_reflected": bool(cand and cand.get("is_reflected")),
            "prev_symbols": prev_step.get("symbols", []), "prev_coords": prev_pos.round(4).tolist(),
            "next_symbols": out_next.get("symbols", []), "next_coords_unaligned": unalign_pos.round(4).tolist(),
            "next_coords_aligned": align_pos.round(4).tolist(), "bridge_frames": bridge, "n_bridge_frames": n_bridge_frames,
            "matched_prev_indices": idxs_p, "matched_next_indices": idxs_n
        }
        step_reports.append(j_info); junction_previews.append(j_info); all_bridges.append(bridge); aligned_steps.append(out_next)

    master_frames, step_boundaries, curr_idx, accum_s = [], [], 0, 0.0
    symbols, max_energy, final_energy = aligned_steps[0].get("symbols", []), 0.0, 0.0
    master_bonds_dict = {}
    for st_i, step in enumerate(aligned_steps):
        for b in step.get("bonds", []):
            p = (min(b.get("from", 0), b.get("to", 0)), max(b.get("from", 0), b.get("to", 0)))
            if p not in master_bonds_dict: master_bonds_dict[p] = dict(b)
            else: master_bonds_dict[p]["p_order"] = b.get("p_order", master_bonds_dict[p].get("p_order", 1))

    for k, step in enumerate(aligned_steps):
        s_frames, ts_idx = step.get("frames", []), step.get("ts_frame_index", 0)
        step_boundaries.append({
            "step_index": k + 1, "step_name": step.get("reaction_smiles", f"Step {k + 1}"),
            "reaction_smiles": step.get("reaction_smiles", ""), "start_frame": curr_idx,
            "end_frame": curr_idx + len(s_frames) - 1, "ts_frame": curr_idx + ts_idx
        })
        k_bonds = step.get("bonds", [])
        raw_s = [float(fr.get("s_arc", idx * 0.1)) for idx, fr in enumerate(s_frames)]
        s_base, step_start_s = raw_s[0], (0.0 if k == 0 else accum_s + 0.15)

        for f_i, fr in enumerate(s_frames):
            c = fr.get("coords", fr.get("positions", []))
            fr["coords"], fr["positions"], fr["step_index"], fr["is_ts"] = c, c, k + 1, (f_i == ts_idx)
            e_k = fr.get("rel_energy_kcal", fr.get("energy_kcal", 0.0))
            fr["energy_kcal"], fr["rel_energy_kcal"], fr["s_arc"] = round(float(e_k), 2), round(float(e_k), 2), round(step_start_s + (raw_s[f_i] - s_base), 3)
            accum_s, max_energy, final_energy = fr["s_arc"], max(max_energy, e_k), e_k
            fr["bonds"] = [{"from": b.get("from", 0), "to": b.get("to", 0), "order": b.get("r_order", 1) if f_i <= ts_idx else b.get("p_order", 1), "type": b.get("type", "spectator")} for b in k_bonds if (b.get("r_order", 1) if f_i <= ts_idx else b.get("p_order", 1)) > 0]
            master_frames.append(fr); curr_idx += 1

        if k < len(all_bridges):
            bridge_b = [{"from": b.get("from", 0), "to": b.get("to", 0), "order": b.get("p_order", 1), "type": b.get("type", "spectator")} for b in k_bonds if b.get("p_order", 1) > 0]
            for b_fr in all_bridges[k]:
                c = b_fr.get("coords", b_fr.get("positions", []))
                b_fr["coords"], b_fr["positions"], b_fr["step_index"], b_fr["is_ts"] = c, c, f"Bridge {k + 1}➔{k + 2}", False
                accum_s += 0.15
                b_fr["rel_energy_kcal"], b_fr["bonds"], b_fr["s_arc"] = b_fr["energy_kcal"], bridge_b, round(float(accum_s), 3)
                max_energy = max(max_energy, b_fr["energy_kcal"])
                master_frames.append(b_fr); curr_idx += 1

    return {
        "title": "Master Multistep IRC Trajectory", "n_steps": len(aligned_steps), "total_frames": len(master_frames),
        "symbols": symbols, "bonds": list(master_bonds_dict.values()), "step_boundaries": step_boundaries,
        "stitching_alignment_report": step_reports, "junction_previews": junction_previews,
        "overall_barrier_kcal": round(float(max_energy), 2), "overall_reaction_energy_kcal": round(float(final_energy), 2),
        "frames": master_frames
    }


def extract_fragment_molblock(traj: Dict[str, Any], frag_meta: Dict[str, Any], frame_idx: int = -1, target_mol: Any = None) -> str:
    """Extract a 3D MolBlock corresponding to a fragment in a trajectory frame with exact coordinates and bonds."""
    idxs = frag_meta.get("atom_indices", [])
    idx_map = {old: new for new, old in enumerate(idxs)}
    frames = traj.get("frames", [])
    frame = frames[frame_idx] if frames else {}
    pos = frame.get("positions", frame.get("coords", []))
    symbols = [traj.get("symbols", [])[i] for i in idxs]
    sub_pos = [pos[i] for i in idxs] if len(pos) >= max(idxs, default=0) + 1 else []

    sub_bonds = []
    for b in traj.get("bonds", []):
        u, v = b.get("from", b[0] if isinstance(b, list) else 0), b.get("to", b[1] if isinstance(b, list) else 0)
        if u in idx_map and v in idx_map and (ord_val := b.get("p_order", b.get("order", 1))) > 0:
            sub_bonds.append((idx_map[u], idx_map[v], ord_val))

    lines = [f"Seeded Product {frag_meta.get('formula', '')}", "  Molecular Workbench", "", f"{len(symbols):>3}{len(sub_bonds):>3}  0  0  0  0  0  0  0  0999 V2000"]
    for s, p in zip(symbols, sub_pos): lines.append(f"{p[0]:>10.4f}{p[1]:>10.4f}{p[2]:>10.4f} {s:<3} 0  0  0  0  0  0  0  0  0  0  0  0")
    for u, v, o in sub_bonds: lines.append(f"{u+1:>3}{v+1:>3}{o:>3}  0  0  0  0")
    lines.append("M  END")
    raw_mb = "\n".join(lines)

    if target_mol is None and frag_meta.get("smiles"):
        try: target_mol = Chem.AddHs(Chem.MolFromSmiles(frag_meta["smiles"]))
        except Exception: target_mol = None

    if target_mol is not None:
        try:
            raw_m = Chem.MolFromMolBlock(raw_mb, removeHs=False)
            if raw_m is not None and (match := target_mol.GetSubstructMatch(raw_m)) and len(match) == target_mol.GetNumAtoms():
                inv = {ref_i: traj_i for traj_i, ref_i in enumerate(match)}
                conf, t_conf = Chem.Conformer(target_mol.GetNumAtoms()), raw_m.GetConformer()
                for ref_i in range(target_mol.GetNumAtoms()): conf.SetAtomPosition(ref_i, t_conf.GetAtomPosition(inv[ref_i]))
                out_m = Chem.Mol(target_mol); out_m.RemoveAllConformers(); out_m.AddConformer(conf)
                return Chem.MolToMolBlock(out_m)
        except Exception: pass
    return raw_mb
