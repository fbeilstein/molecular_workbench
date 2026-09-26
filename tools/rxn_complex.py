import os
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.transform import Rotation
from rxn_mapping import _get_rdkit
from rxn_mechanism import detect_driving_coordinates, map_driving_coords_to_complex


def _read_xyz(path):
    """Read xyz file → list of (symbol, x, y, z)."""
    with open(path) as f:
        lines = f.readlines()
    n = int(lines[0].strip())
    atoms = []
    for line in lines[2:2 + n]:
        parts = line.split()
        atoms.append((parts[0], float(parts[1]), float(parts[2]), float(parts[3])))
    return atoms


def _write_xyz_atoms(atoms, path, comment='Complex'):
    """Write list of (sym, x, y, z) tuples to xyz file."""
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for sym, x, y, z in atoms:
            f.write(f"{sym:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def _build_neb_complexes(r_xyz_paths, p_xyz_paths, mapping, reactants_smi,
                         products_smi, out_dir, name, spacing=2.5):
    """
    Build matched reactant/product complex xyz files for NEB.
    
    NEB requires: same atom count, same atom types, same order.
    To avoid tearing product molecules apart, it rigidly translates whole
    product fragments to positions minimizing distance to mapped R atoms.
    """
    Chem, _, _, _ = _get_rdkit()

    # Detect driving coordinates (formed/broken bonds)
    driving_coords = detect_driving_coordinates(reactants_smi, products_smi, mapping)
    formed_bonds = driving_coords['formed_bonds'] + driving_coords.get('proton_transfers', [])
    broken_bonds = driving_coords['broken_bonds']

    # For proton-transfer-only reactions, use wider spacing so the scan
    # trajectory shows the full approach of the proton acceptor
    is_pt_only = (not driving_coords['formed_bonds'] and
                  not driving_coords['broken_bonds'] and
                  bool(driving_coords.get('proton_transfers')))
    if is_pt_only:
        spacing = 4.0

    # Read all fragment geometries
    r_frags = [_read_xyz(p) for p in r_xyz_paths]
    p_frags = [_read_xyz(p) for p in p_xyz_paths]

    # Build reactant complex intelligently
    r_complex = []
    r_atom_origins = []  # Track (frag_idx, atom_idx) for each complex atom
    placed_frags = set()
    
    def get_coord(fi, ai):
        for ci, (f, a) in enumerate(r_atom_origins):
            if f == fi and a == ai:
                return r_complex[ci][1:4]
        return None

    if len(r_frags) > 0:
        fi = 0
        xs = [a[1] for a in r_frags[0]]
        cx = (min(xs) + max(xs)) / 2.0
        for ai, (sym, x, y, z) in enumerate(r_frags[0]):
            r_complex.append((sym, x - cx, y, z))
            r_atom_origins.append((fi, ai))
        placed_frags.add(0)

    for fi in range(1, len(r_frags)):
        frag = r_frags[fi]
        if not frag: continue
        
        target_f, target_a, source_a = None, None, None
        for (r1, r2) in formed_bonds:
            if r1[0] == fi and r2[0] in placed_frags:
                source_a, target_f, target_a = r1[1], r2[0], r2[1]
                break
            elif r2[0] == fi and r1[0] in placed_frags:
                source_a, target_f, target_a = r2[1], r1[0], r1[1]
                break
                
        if target_f is not None:
             leave_a = None
             for (b1, b2) in broken_bonds:
                 if b1 == (target_f, target_a): leave_a = b2[1]
                 elif b2 == (target_f, target_a): leave_a = b1[1]
                 
             target_coord = get_coord(target_f, target_a)
             if leave_a is not None and target_coord is not None:
                 leave_coord = get_coord(target_f, leave_a)
                 v_dist = sum((a-b)**2 for a,b in zip(target_coord, leave_coord))**0.5
                 vec = [(t-l)/v_dist for t, l in zip(target_coord, leave_coord)] if v_dist > 0 else [1.0, 0.0, 0.0]
             else:
                 import numpy as np
                 t_frag = r_frags[target_f]
                 t_coord = np.array(t_frag[target_a][1:4])
                 neighbors = []
                 for i, (sym, x, y, z) in enumerate(t_frag):
                     if i != target_a:
                         pos = np.array([x, y, z])
                         if np.linalg.norm(pos - t_coord) < 1.8:
                             neighbors.append(pos)
                 
                 vec = [1.0, 0.0, 0.0]
                 if len(neighbors) >= 2:
                     v1 = neighbors[0] - t_coord
                     v2 = neighbors[1] - t_coord
                     n = np.cross(v1, v2)
                     norm = np.linalg.norm(n)
                     if norm > 1e-3:
                         vec = (n / norm).tolist()
                 elif len(neighbors) == 1:
                     v = t_coord - neighbors[0]
                     norm = np.linalg.norm(v)
                     if norm > 1e-3:
                         vec = (v / norm).tolist()
                 
             import numpy as np
             desired_source_pos = np.array([t + v*(spacing) for t, v in zip(target_coord, vec)])
             curr_source_pos = np.array(frag[source_a][1:4])
             
             # Calculate rotation to point source_a toward target
             # 1. Get the centroid of the fragment EXCLUDING source_a
             if len(frag) > 1:
                 other_pts = np.array([frag[i][1:4] for i in range(len(frag)) if i != source_a])
                 centroid = np.mean(other_pts, axis=0)
                 local_vec = centroid - curr_source_pos
                 if np.linalg.norm(local_vec) > 0:
                     local_vec = local_vec / np.linalg.norm(local_vec)
                     # We want the fragment's bulk to point OUTWARD (aligned with vec)
                     target_vec = np.array(vec)
                     target_vec = target_vec / np.linalg.norm(target_vec)
                     axis = np.cross(local_vec, target_vec)
                     axis_len = np.linalg.norm(axis)
                     dot = np.dot(local_vec, target_vec)
                     if axis_len > 1e-5:
                         axis = axis / axis_len
                         angle = np.arccos(np.clip(dot, -1.0, 1.0))
                         K = np.array([[0, -axis[2], axis[1]],
                                       [axis[2], 0, -axis[0]],
                                       [-axis[1], axis[0], 0]])
                         R = np.eye(3) + np.sin(angle)*K + (1-np.cos(angle))*(K@K)
                     elif dot < 0:
                         R = -np.eye(3)
                     else:
                         R = np.eye(3)
                 else:
                     R = np.eye(3)
             else:
                 R = np.eye(3)

             # Apply rotation and translation
             for ai, (sym, x, y, z) in enumerate(frag):
                 pt = np.array([x, y, z]) - curr_source_pos
                 pt_rot = R @ pt
                 pt_final = pt_rot + desired_source_pos
                 r_complex.append((sym, pt_final[0], pt_final[1], pt_final[2]))
                 r_atom_origins.append((fi, ai))
        else:
             offset = 5.0 * fi
             for ai, (sym, x, y, z) in enumerate(frag):
                 r_complex.append((sym, x + offset, y, z))
                 r_atom_origins.append((fi, ai))
        placed_frags.add(fi)

    # Align product fragments using Kabsch algorithm
    import numpy as np

    def kabsch_align(P, Q, allow_reflection=False):
        # P and Q are Nx3 matrices. Aligns P onto Q.
        P_centroid = np.mean(P, axis=0)
        Q_centroid = np.mean(Q, axis=0)
        P_centered = P - P_centroid
        Q_centered = Q - Q_centroid
        
        H = P_centered.T @ Q_centered
        U, S, Vt = np.linalg.svd(H)
        R = Vt.T @ U.T
        if not allow_reflection and np.linalg.det(R) < 0:
            Vt[2, :] *= -1
            R = Vt.T @ U.T
        t = Q_centroid - P_centroid @ R.T
        return R, t

    translated_p_frags = []
    for p_fi, frag in enumerate(p_frags):
        p_sources = []
        r_targets = []
        
        is_reactive = False
        
        # Collect all mapped atoms for this product fragment
        for p_ai in range(len(frag)):
            for rk, pv in mapping.items():
                if pv[1] == p_fi and pv[2] == p_ai:
                    r_fi, r_ai = rk[1], rk[2]
                    
                    # Check if atom is involved in heavy-atom formed/broken bonds (not PT)
                    for (b1, b2) in (driving_coords['formed_bonds'] + driving_coords['broken_bonds']):
                        if (r_fi, r_ai) == b1 or (r_fi, r_ai) == b2:
                            is_reactive = True
                            
                    # Find Rcoord in r_complex
                    for ci, (r_f, r_a) in enumerate(r_atom_origins):
                        if r_f == r_fi and r_a == r_ai:
                            p_sources.append(frag[p_ai][1:4])
                            r_targets.append(r_complex[ci][1:4])
                            break
                            
        P_pts = np.array(p_sources)
        Q_pts = np.array(r_targets)
        
        if len(P_pts) >= 2:
            # Full 3D alignment (or 2-point vector alignment)
            # Allow reflection if fragment is a spectator or only does PT
            R, t = kabsch_align(P_pts, Q_pts, allow_reflection=(not is_reactive))
            t_frag = []
            for sym, x, y, z in frag:
                p_vec = np.array([x, y, z])
                p_rot = R @ p_vec + t
                t_frag.append((sym, p_rot[0], p_rot[1], p_rot[2]))
        elif len(P_pts) == 1:
            # Translation only (1 atom)
            centroid_P = np.mean(P_pts, axis=0)
            centroid_Q = np.mean(Q_pts, axis=0)
            shift = centroid_Q - centroid_P
            t_frag = [(sym, x+shift[0], y+shift[1], z+shift[2]) for sym, x, y, z in frag]
        else:
            # Fragment has no mapped atoms (e.g. leaving [H+] or [Br-])
            # Translate it away from the rest of the placed product fragments
            if len(translated_p_frags) > 0:
                all_pts = []
                for prev_frag in translated_p_frags:
                    all_pts.extend([a[1:4] for a in prev_frag])
                if len(all_pts) > 0:
                    centroid = np.mean(all_pts, axis=0)
                    frag_centroid = np.mean([a[1:4] for a in frag], axis=0)
                    # Push it 8.0 A away from the centroid
                    push_dir = frag_centroid - centroid
                    norm = np.linalg.norm(push_dir)
                    if norm > 1e-3:
                         push_dir = push_dir / norm
                    else:
                         push_dir = np.array([0.0, 0.0, 1.0])
                    shift = (centroid + push_dir * 8.0) - frag_centroid
                    t_frag = [(sym, x+shift[0], y+shift[1], z+shift[2]) for sym, x, y, z in frag]
                else:
                    t_frag = frag.copy()
            else:
                t_frag = frag.copy()
            
        translated_p_frags.append(t_frag)

    # Now build matched p_complex taking atoms from translated_p_frags
    p_complex = [None] * len(r_complex)
    used_p_atoms = set()

    # Step 1: Place mapped heavy atoms
    for rk, pv in mapping.items():
        r_fi, r_ai = rk[1], rk[2]
        p_fi, p_ai = pv[1], pv[2]
        for ci, (fi, ai) in enumerate(r_atom_origins):
            if fi == r_fi and ai == r_ai:
                if p_fi < len(translated_p_frags) and p_ai < len(translated_p_frags[p_fi]):
                    p_atom = translated_p_frags[p_fi][p_ai]
                    p_complex[ci] = p_atom
                    used_p_atoms.add((p_fi, p_ai))
                break

    # Step 2: Match remaining atoms (hydrogens) globally minimizing distance
    import numpy as np
    from scipy.optimize import linear_sum_assignment
    
    unmapped_r = []
    for ci, (sym, rx, ry, rz) in enumerate(r_complex):
        if p_complex[ci] is None:
            unmapped_r.append((ci, sym, np.array([rx, ry, rz])))
            
    unmapped_p = []
    for p_fi in range(len(translated_p_frags)):
        for p_ai, p_atom in enumerate(translated_p_frags[p_fi]):
            if (p_fi, p_ai) not in used_p_atoms:
                unmapped_p.append(((p_fi, p_ai), p_atom[0], np.array(p_atom[1:4])))
                
    for symbol in set(u[1] for u in unmapped_r):
        r_group = [u for u in unmapped_r if u[1] == symbol]
        p_group = [u for u in unmapped_p if u[1] == symbol]
        
        if len(r_group) != len(p_group):
            print(f"  Warning: Unmapped count mismatch for {symbol}: R={len(r_group)}, P={len(p_group)}")
            continue
            
        cost_matrix = np.zeros((len(r_group), len(p_group)))
        for i, r in enumerate(r_group):
            for j, p in enumerate(p_group):
                cost_matrix[i, j] = np.linalg.norm(r[2] - p[2])
                
        row_ind, col_ind = linear_sum_assignment(cost_matrix)
        for r_idx, p_idx in zip(row_ind, col_ind):
            ci = r_group[r_idx][0]
            p_fi, p_ai = p_group[p_idx][0]
            p_complex[ci] = translated_p_frags[p_fi][p_ai]
            used_p_atoms.add((p_fi, p_ai))
            
    # Fallback greedy for any remaining (due to count mismatch)
    for ci, (sym, rx, ry, rz) in enumerate(r_complex):
        if p_complex[ci] is not None:
            continue
            
        r_pt = np.array([rx, ry, rz])
        best_p_fi, best_p_ai = None, None
        min_dist = float('inf')
        
        for p_fi in range(len(translated_p_frags)):
            for p_ai, p_atom in enumerate(translated_p_frags[p_fi]):
                if (p_fi, p_ai) in used_p_atoms:
                    continue
                if p_atom[0] == sym:
                    p_pt = np.array(p_atom[1:4])
                    dist = np.linalg.norm(r_pt - p_pt)
                    if dist < min_dist:
                        min_dist = dist
                        best_p_fi = p_fi
                        best_p_ai = p_ai
                        
        if best_p_fi is not None:
            p_complex[ci] = translated_p_frags[best_p_fi][best_p_ai]
            used_p_atoms.add((best_p_fi, best_p_ai))
        else:
            p_complex[ci] = (sym, rx, ry, rz)
            print(f"  Warning: no P match for R atom {ci} ({sym}), using R position")

    # Post-process: for dissociation products (broken bonds), push the leaving
    # group atom away from its former partner so the optimizer can't recombine them.
    # This is critical for E1, SN1, etc. where ion pairs form.
    import math as _math
    MIN_LEAVING_DIST = 8.0  # Å — far beyond electrostatic recombination range
    
    separations_to_push = list(broken_bonds)
    if 'proton_transfers' in driving_coords:
        separations_to_push.extend(driving_coords['proton_transfers'])
        
    for (b1, b2) in separations_to_push:
        # b1 and b2 are (frag_idx, atom_idx) in reactant frame
        # Find their complex indices
        ci1, ci2 = None, None
        for ci, (fi, ai) in enumerate(r_atom_origins):
            if fi == b1[0] and ai == b1[1]:
                ci1 = ci
            if fi == b2[0] and ai == b2[1]:
                ci2 = ci
        if ci1 is None or ci2 is None:
            continue
        # Check if they're in different product fragments
        p1_frag, p2_frag = None, None
        for rk, pv in mapping.items():
            if rk[1] == b1[0] and rk[2] == b1[1]:
                p1_frag = pv[1]
            if rk[1] == b2[0] and rk[2] == b2[1]:
                p2_frag = pv[1]
        if p1_frag is not None and p2_frag is not None and p1_frag != p2_frag:
            # These atoms are in different product fragments — push apart
            pa1 = p_complex[ci1]
            pa2 = p_complex[ci2]
            dx = pa2[1] - pa1[1]
            dy = pa2[2] - pa1[2]
            dz = pa2[3] - pa1[3]
            dist = _math.sqrt(dx*dx + dy*dy + dz*dz)
            if dist < MIN_LEAVING_DIST and dist > 0.01:
                # Push the smaller fragment (fewer atoms) away
                # Determine which complex index is the "leaving" atom
                frag1_size = sum(1 for _, pv in mapping.items() if pv[1] == p1_frag)
                frag2_size = sum(1 for _, pv in mapping.items() if pv[1] == p2_frag)
                if frag2_size <= frag1_size:
                    # Push ci2 (and its fragment) away from ci1
                    push_atom, anchor_atom = ci2, ci1
                    push_frag = p2_frag
                else:
                    push_atom, anchor_atom = ci1, ci2
                    push_frag = p1_frag
                # Direction: from anchor toward push atom
                vec = [dx/dist, dy/dist, dz/dist] if push_atom == ci2 else [-dx/dist, -dy/dist, -dz/dist]
                extra = MIN_LEAVING_DIST - dist
                # Shift ALL atoms belonging to this product fragment
                push_indices = set()
                for rk, pv in mapping.items():
                    if pv[1] == push_frag:
                        for ci_s, (fi_s, ai_s) in enumerate(r_atom_origins):
                            if fi_s == rk[1] and ai_s == rk[2]:
                                push_indices.add(ci_s)
                # Also include unmapped atoms assigned to this fragment
                for ci_s in range(len(p_complex)):
                    if ci_s in push_indices:
                        continue
                    # Check if this unmapped atom is near any pushed atom
                    if p_complex[ci_s] is not None and ci_s not in push_indices:
                        for pi in list(push_indices):
                            pa_p = p_complex[pi]
                            pa_c = p_complex[ci_s]
                            d = _math.sqrt((pa_p[1]-pa_c[1])**2 + (pa_p[2]-pa_c[2])**2 + (pa_p[3]-pa_c[3])**2)
                            if d < 2.0:  # Close enough to be same fragment
                                push_indices.add(ci_s)
                                break
                if not push_indices:
                    push_indices = {push_atom}
                for ci_s in push_indices:
                    pa = p_complex[ci_s]
                    p_complex[ci_s] = (pa[0],
                                       pa[1] + vec[0]*extra,
                                       pa[2] + vec[1]*extra,
                                       pa[3] + vec[2]*extra)
                sym1 = p_complex[anchor_atom][0]
                sym2 = p_complex[push_atom][0]
                new_dist = dist + extra
                print(f"  ⚡ Pushed leaving group {sym2}({push_atom}) away from "
                      f"{sym1}({anchor_atom}): {dist:.2f} → {new_dist:.2f} Å")

    # Verify atom types
    mismatches = []
    for i, (ra, pa) in enumerate(zip(r_complex, p_complex)):
        if ra[0] != pa[0]:
            mismatches.append(f"  Atom {i}: R={ra[0]} vs P={pa[0]}")
    if mismatches:
        print(f"  WARNING: Atom type mismatches in complex!")
        for m in mismatches:
            print(m)

    r_path = os.path.join(out_dir, f"{name}_reactant_complex.xyz")
    p_path = os.path.join(out_dir, f"{name}_product_complex.xyz")
    _write_xyz_atoms(r_complex, r_path, f"Reactant complex ({len(r_complex)} atoms)")
    _write_xyz_atoms(p_complex, p_path, f"Product complex ({len(p_complex)} atoms)")

    # Map driving coordinates to complex atom indices
    complex_driving = map_driving_coords_to_complex(driving_coords, r_atom_origins)

    return r_path, p_path, complex_driving