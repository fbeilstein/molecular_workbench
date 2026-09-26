import os
import numpy as np
from pyscf import lo
from pyscf.tools import cubegen
from .geometry import make_atom_ids
from .classification import _CORE_ELECTRONS, _atom_populations, _classify_orbital, _is_core_orbital

_COVALENT_R_BOHR = {
    'H': 0.6 * 1.88973,
    'C': 1.4 * 1.88973,
    'N': 1.3 * 1.88973,
    'O': 1.25 * 1.88973,
    'F': 1.15 * 1.88973,
}

def _find_rings(mol, max_size=8):
    """Find small rings (3-8 atoms) in the molecule using DFS.
    Returns list of frozensets of atom indices, excluding H atoms.
    """
    coords = mol.atom_coords()  # Bohr
    natm = mol.natm
    syms = [mol.atom_symbol(i) for i in range(natm)]

    # Build adjacency from covalent radii
    adj = [[] for _ in range(natm)]
    for i in range(natm):
        for j in range(i + 1, natm):
            ri = _COVALENT_R_BOHR.get(syms[i], 1.45)
            rj = _COVALENT_R_BOHR.get(syms[j], 1.45)
            if np.linalg.norm(coords[i] - coords[j]) < (ri + rj) * 1.3:
                adj[i].append(j)
                adj[j].append(i)

    rings = set()
    for start in range(natm):
        if syms[start] == 'H':
            continue
        _dfs_rings(adj, start, [start], {start}, rings, syms, max_size)

    # Deduplicate: keep only minimal (non-superset) rings
    unique = []
    for ring in sorted(rings, key=len):
        if not any(prev.issubset(ring) and len(prev) < len(ring)
                   for prev in unique):
            unique.append(ring)
    return unique


def _dfs_rings(adj, current, path, visited, rings, syms, max_size):
    if len(path) > max_size:
        return
    for nb in adj[current]:
        if syms[nb] == 'H':
            continue
        if nb == path[0] and len(path) >= 3:
            rings.add(frozenset(path))
            continue
        if nb in visited:
            continue
        visited.add(nb)
        path.append(nb)
        _dfs_rings(adj, nb, path, visited, rings, syms, max_size)
        path.pop()
        visited.discard(nb)


def _ring_plane_normals(mol, rings):
    """Compute plane normal for each ring via PCA.
    Returns list of (frozenset_of_atom_indices, unit_normal_vector).
    """
    coords = mol.atom_coords()
    result = []
    for ring in rings:
        ring_coords = coords[list(ring)]
        centered = ring_coords - ring_coords.mean(axis=0)
        _, _, Vt = np.linalg.svd(centered)
        normal = Vt[-1]
        normal /= np.linalg.norm(normal)
        result.append((ring, normal))
    return result


def _find_pi_canonical(mol, mf, atom_labels):
    """Identify canonical MOs with π character using PCA-based ring planes.

    For each ring system, finds the plane normal via SVD.  Then for each
    valence canonical MO, checks whether its AO coefficients on ring atoms
    are dominated by p-orbitals aligned with the plane normal (= π character).

    Returns list of (mo_index, energy_ev) for each π MO, sorted by energy.
    """
    rings = _find_rings(mol)
    if not rings:
        return []

    ring_planes = _ring_plane_normals(mol, rings)

    # Union of all ring atom indices
    all_ring_atoms = set()
    for ring, _ in ring_planes:
        all_ring_atoms |= ring

    # Collect all unique plane normals (group coplanar rings)
    normals = []
    for _, n in ring_planes:
        # Check if this normal is already represented (parallel or anti-parallel)
        duplicate = False
        for existing in normals:
            if abs(abs(np.dot(n, existing)) - 1.0) < 0.1:
                duplicate = True
                break
        if not duplicate:
            normals.append(n)

    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy
    if isinstance(mo_coeff, (tuple, list)) or (isinstance(mo_coeff, np.ndarray) and mo_coeff.ndim == 3):
        mo_coeff = mo_coeff[0]
        mo_energy = mo_energy[0]
        n_occ = mol.nelec[0]
    else:
        n_occ = mol.nelectron // 2

    # Core count
    core_electrons = sum(_CORE_ELECTRONS.get(s, 0) for s in atom_labels)
    n_core = core_electrons // 2

    ao_labels = mol.ao_labels(fmt=False)
    p_dirs = {'x': np.array([1, 0, 0]),
              'y': np.array([0, 1, 0]),
              'z': np.array([0, 0, 1])}

    pi_mos = []
    for i in range(n_core, n_occ):
        mo = mo_coeff[:, i]

        # For each unique ring plane normal, accumulate perpendicular vs parallel
        # p-character on ring atoms
        is_pi = False
        for normal in normals:
            p_perp, p_along, total_ring = 0.0, 0.0, 0.0
            for mu, (aidx, *rest) in enumerate(ao_labels):
                if aidx not in all_ring_atoms:
                    continue
                c2 = mo[mu] ** 2
                total_ring += c2
                sublabel = rest[-1] if rest else ''
                if sublabel in p_dirs:
                    proj = abs(np.dot(p_dirs[sublabel], normal))
                    p_perp += c2 * proj
                    p_along += c2 * (1 - proj)

            # Significant ring density + dominated by perpendicular p
            if total_ring > 0.1 and p_perp > p_along * 2:
                is_pi = True
                break

        if is_pi:
            energy_ev = float(mo_energy[i]) * 27.2114
            pi_mos.append((i, energy_ev))

    pi_mos.sort(key=lambda x: x[1])
    return pi_mos


def compute_canonical(mol, mf, atom_labels, name, out_dir, grid_points=50):
    """Generate HOMO, LUMO, and canonical π system cube files.

    Instead of blindly exporting HOMO-N, detects which canonical MOs are
    π-type using PCA-based ring plane analysis and exports only those.
    """
    from pyscf import tools
    import numpy as np
    from pyscf.tools import cubegen

    # Handle open-shell (UHF/UKS) where mo_coeff is a tuple
    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy
    if isinstance(mo_coeff, (tuple, list)) or (isinstance(mo_coeff, np.ndarray) and mo_coeff.ndim == 3):
        mo_coeff = mo_coeff[0]
        mo_energy = mo_energy[0]
        n_alpha = mol.nelec[0]
        homo_idx = n_alpha - 1
    else:
        homo_idx = mol.nelectron // 2 - 1
    lumo_idx = homo_idx + 1
    results = {}

    # Always export HOMO and LUMO
    for label, idx in [('homo', homo_idx), ('lumo', lumo_idx)]:
        if idx < 0 or idx >= mo_coeff.shape[1]:
            continue
        energy_ev = mo_energy[idx] * 27.2114
        cube_name = f"{name}_{label}.cube"
        cube_path = os.path.join(out_dir, cube_name)
        cubegen.orbital(mol, cube_path, mo_coeff[:, idx],
                              nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
        results[label] = {
            'file': cube_name,
            'energy_ev': round(energy_ev, 2),
            'index': idx,
        }
        print(f"  ✓ {label.upper()} (E={energy_ev:+.2f} eV) → {cube_name}")

    # Detect and export canonical π MOs
    pi_mos = _find_pi_canonical(mol, mf, atom_labels)
    pi_results = []
    for mo_idx, energy_ev in pi_mos:
        depth = homo_idx - mo_idx
        if depth == 0:
            # Already exported as HOMO
            pi_results.append({
                'file': results['homo']['file'],
                'energy_ev': round(energy_ev, 2),
                'index': mo_idx,
                'homo_label': 'HOMO',
            })
            continue
        label = f"pi_{mo_idx}"
        cube_name = f"{name}_{label}.cube"
        cube_path = os.path.join(out_dir, cube_name)
        cubegen.orbital(mol, cube_path, mo_coeff[:, mo_idx],
                              nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
        homo_label = f"HOMO-{depth}" if depth > 0 else "HOMO"
        pi_results.append({
            'file': cube_name,
            'energy_ev': round(energy_ev, 2),
            'index': mo_idx,
            'homo_label': homo_label,
        })
        print(f"  ✓ π ({homo_label}, E={energy_ev:+.2f} eV) → {cube_name}")

    if pi_results:
        results['pi_system'] = pi_results
        print(f"  Canonical π system: {len(pi_results)} MOs detected")

    return results


def compute_localized(mol, mf, atom_labels, name, out_dir, grid_points=50):
    """Localize MOs using hybrid Boys/IBO and export to cube/json."""
    import os
    from pyscf.tools import cubegen
    from pyscf import lo
    
    n_core = 0
    for sym in atom_labels:
        n_core += _CORE_ELECTRONS.get(sym, 0)
    n_core //= 2
    n_occ = mol.nelectron // 2
    
    mo_coeff = mf.mo_coeff
    ovlp = mol.intor_symmetric('int1e_ovlp')
    atom_ids = make_atom_ids(atom_labels)
    lp_counts = {}

    sigma, pi, lone_pairs = [], [], []

    # -- Localize Occupieds (IBO is perfect for occupieds) --
    # Construct full IAOs from all occupied orbitals to form proper minimal basis
    full_occ = mo_coeff[:, :n_occ]
    iaos = lo.iao.iao(mol, full_occ)

    # Localize ONLY valence orbitals using the full IAOs to prevent core mixing
    val_occ = mo_coeff[:, n_core:n_occ]
    print(f"  Localizing {val_occ.shape[1]} valence MOs (IBO)...")
    try:
        occ_loc = lo.ibo.ibo(mol, val_occ, iaos=iaos)
        for i in range(occ_loc.shape[1]):
            mo = occ_loc[:, i]
            pop = _atom_populations(mol, mo, ovlp)
            info = _classify_orbital(mol, mo, pop, atom_labels, atom_ids, ovlp)
            
            if info['type'] == 'lone_pair':
                a_id = info['atom']
                cube_name = f"{name}_lp_{a_id}_{lp_counts.get(a_id,0)+1}_{i}.cube"
                lp_counts[a_id] = lp_counts.get(a_id, 0) + 1
                cube_path = os.path.join(out_dir, cube_name)
                cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                lone_pairs.append({'atom': a_id, 'index': lp_counts[a_id], 'file': cube_name})
                print(f"  ✓ LP({a_id} #{lp_counts[a_id]}) → {cube_name}")
            elif info['type'] == 'sigma':
                a1, a2 = info['atoms']
                cube_name = f"{name}_sigma_{a1}_{a2}_{i}.cube"
                cube_path = os.path.join(out_dir, cube_name)
                cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                sigma.append({'atoms': [a1, a2], 'file': cube_name})
                print(f"  ✓ σ({a1}–{a2}) → {cube_name}")
            elif info['type'] == 'pi':
                a1, a2 = info['atoms']
                cube_name = f"{name}_pi_{a1}_{a2}_{i}.cube"
                cube_path = os.path.join(out_dir, cube_name)
                cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                pi.append({'atoms': [a1, a2], 'file': cube_name})
                print(f"  ✓ π({a1}={a2}) → {cube_name}")

    except Exception as e:
        print(f"  ⚠ Localization failed: {e}")

    # -- Localize Virtuals (IBO with Occupied IAOs) --
    sigma_star, pi_star = [], []
    n_virt_to_loc = len(sigma) + len(pi)
    if n_virt_to_loc > 0 and n_occ + n_virt_to_loc <= mo_coeff.shape[1]:
        virt_coeff = mo_coeff[:, n_occ : n_occ + n_virt_to_loc]
        print(f"  Localizing {virt_coeff.shape[1]} lowest virtual MOs (IBO)...")
        try:
            # Re-use the exact same IAOs from the occupied space to project virtuals!
            # This is the secret to avoiding leakage and preventing unphysical banana bonds.
            virt_loc = lo.ibo.ibo(mol, virt_coeff, iaos=iaos)
            
            sigma_star_idx, pi_star_idx = 0, 0
            for i in range(virt_loc.shape[1]):
                mo = virt_loc[:, i]
                pop = _atom_populations(mol, mo, ovlp)
                info = _classify_orbital(mol, mo, pop, atom_labels, atom_ids, ovlp)
                
                if info['type'] == 'sigma':
                    a1, a2 = info['atoms']
                    cube_name = f"{name}_sigmastar_{a1}_{a2}_{sigma_star_idx}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    sigma_star.append({'atoms': [a1, a2], 'file': cube_name})
                    print(f"  ✓ σ*({a1}–{a2}) → {cube_name} [IBO]")
                    sigma_star_idx += 1
                
                elif info['type'] == 'pi':
                    a1, a2 = info['atoms']
                    cube_name = f"{name}_pistar_{a1}_{a2}_{pi_star_idx}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    pi_star.append({'atoms': [a1, a2], 'file': cube_name})
                    print(f"  ✓ π*({a1}–{a2}) → {cube_name} [IBO]")
                    pi_star_idx += 1

        except Exception as e:
            print(f"  ⚠ Virtual localization failed: {e}")

    return {
        'sigma': sigma,
        'pi': pi,
        'lone_pairs': lone_pairs,
        'sigma_star': sigma_star,
        'pi_star': pi_star
    }


