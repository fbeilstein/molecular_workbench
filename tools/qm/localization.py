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


def _find_pi_canonical(mol, mf, atom_labels, n_virt_scan=None):
    """Identify canonical MOs with π character using PCA-based ring planes.

    For each ring system, finds the plane normal via SVD.  Then for each
    valence canonical MO (occupied AND virtual), checks whether its AO
    coefficients on ring atoms are dominated by p-orbitals aligned with
    the plane normal (= π character).

    Args:
        n_virt_scan: How many virtual MOs to scan for π*. If None, scans
                     as many virtuals as there are valence occupieds.

    Returns:
        (pi_occ, pi_virt) where each is a list of (mo_index, energy_ev),
        sorted by energy. pi_occ = bonding π, pi_virt = antibonding π*.
    """
    rings = _find_rings(mol)
    if not rings:
        return [], []

    ring_planes = _ring_plane_normals(mol, rings)

    # Union of all ring atom indices
    all_ring_atoms = set()
    for ring, _ in ring_planes:
        all_ring_atoms |= ring

    # Collect all unique plane normals (group coplanar rings)
    normals = []
    for _, n in ring_planes:
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

    core_electrons = sum(_CORE_ELECTRONS.get(s, 0) for s in atom_labels)
    n_core = core_electrons // 2
    n_valence = n_occ - n_core

    ao_labels = mol.ao_labels(fmt=False)
    p_dirs = {'x': np.array([1, 0, 0]),
              'y': np.array([0, 1, 0]),
              'z': np.array([0, 0, 1])}

    def _is_pi_mo(mo_idx):
        """Check if canonical MO at mo_idx has π character."""
        mo = mo_coeff[:, mo_idx]
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

            if total_ring > 0.1 and p_perp > p_along * 2:
                return True
        return False

    # Scan occupied valence MOs for bonding π
    pi_occ = []
    for i in range(n_core, n_occ):
        if _is_pi_mo(i):
            energy_ev = float(mo_energy[i]) * 27.2114
            pi_occ.append((i, energy_ev))
    pi_occ.sort(key=lambda x: x[1])

    # Scan virtual MOs for antibonding π*
    if n_virt_scan is None:
        n_virt_scan = n_valence
    n_total = mo_coeff.shape[1]
    virt_end = min(n_occ + n_virt_scan, n_total)

    pi_virt = []
    for i in range(n_occ, virt_end):
        if _is_pi_mo(i):
            energy_ev = float(mo_energy[i]) * 27.2114
            pi_virt.append((i, energy_ev))
    pi_virt.sort(key=lambda x: x[1])

    return pi_occ, pi_virt


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

    # Detect and export canonical π MOs (occupied bonding)
    pi_occ, pi_virt = _find_pi_canonical(mol, mf, atom_labels)
    pi_results = []
    for mo_idx, energy_ev in pi_occ:
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
        print(f"  Canonical π system: {len(pi_results)} bonding MOs detected")

    # Detect and export canonical π* MOs (virtual antibonding)
    pi_star_results = []
    for mo_idx, energy_ev in pi_virt:
        offset = mo_idx - lumo_idx
        if offset == 0:
            # Already exported as LUMO
            pi_star_results.append({
                'file': results['lumo']['file'],
                'energy_ev': round(energy_ev, 2),
                'index': mo_idx,
                'lumo_label': 'LUMO',
            })
            continue
        label = f"pistar_{mo_idx}"
        cube_name = f"{name}_{label}.cube"
        cube_path = os.path.join(out_dir, cube_name)
        cubegen.orbital(mol, cube_path, mo_coeff[:, mo_idx],
                              nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
        lumo_label = f"LUMO+{offset}" if offset > 0 else "LUMO"
        pi_star_results.append({
            'file': cube_name,
            'energy_ev': round(energy_ev, 2),
            'index': mo_idx,
            'lumo_label': lumo_label,
        })
        print(f"  ✓ π* ({lumo_label}, E={energy_ev:+.2f} eV) → {cube_name}")

    if pi_star_results:
        results['pi_star_system'] = pi_star_results
        print(f"  Canonical π* system: {len(pi_star_results)} antibonding MOs detected")

    return results


def compute_localized(mol, mf, atom_labels, name, out_dir, grid_points=50):
    """Localize MOs and export to cube files.

    Strategy:
    - Detect canonical π/π* MOs (ring plane analysis) and export them directly
      (canonical MOs look correct for π systems — textbook shapes).
    - Exclude π MOs from the localization subspace.
    - Localize only the σ/LP occupied MOs (IBO) and σ* virtual MOs (PM).
    """
    import os
    from pyscf.tools import cubegen
    from pyscf import lo

    n_core = 0
    for sym in atom_labels:
        n_core += _CORE_ELECTRONS.get(sym, 0)
    n_core //= 2
    n_occ = mol.nelectron // 2

    mo_coeff = mf.mo_coeff
    mo_energy = mf.mo_energy
    ovlp = mol.intor_symmetric('int1e_ovlp')
    atom_ids = make_atom_ids(atom_labels)
    lp_counts = {}

    sigma, pi, lone_pairs = [], [], []
    sigma_star, pi_star = [], []

    # ── Step 1: Detect canonical π/π* MOs ──
    n_valence = n_occ - n_core
    pi_occ_mos, pi_virt_mos = _find_pi_canonical(mol, mf, atom_labels,
                                                   n_virt_scan=n_valence)
    pi_occ_indices = set(idx for idx, _ in pi_occ_mos)
    pi_virt_indices = set(idx for idx, _ in pi_virt_mos)
    n_pi_occ = len(pi_occ_indices)
    n_pi_virt = len(pi_virt_indices)

    if n_pi_occ > 0:
        print(f"  Detected {n_pi_occ} canonical π MOs — exporting directly")
    if n_pi_virt > 0:
        print(f"  Detected {n_pi_virt} canonical π* MOs — exporting directly")

    # ── Step 2: Export canonical π bonds directly ──
    for mo_idx, energy_ev in pi_occ_mos:
        depth = (n_occ - 1) - mo_idx
        homo_label = f"HOMO-{depth}" if depth > 0 else "HOMO"
        cube_name = f"{name}_pi_canonical_{mo_idx}.cube"
        cube_path = os.path.join(out_dir, cube_name)
        cubegen.orbital(mol, cube_path, mo_coeff[:, mo_idx],
                        nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
        pi.append({
            'atoms': ['ring'],  # canonical — spread over ring
            'file': cube_name,
            'energy_ev': round(energy_ev, 2),
            'canonical_label': f"π ({homo_label})",
        })
        print(f"  ✓ π ({homo_label}, E={energy_ev:+.2f} eV) → {cube_name}")

    # ── Step 3: Export canonical π* antibonds directly ──
    lumo_idx = n_occ
    for mo_idx, energy_ev in pi_virt_mos:
        offset = mo_idx - lumo_idx
        lumo_label = f"LUMO+{offset}" if offset > 0 else "LUMO"
        cube_name = f"{name}_pistar_canonical_{mo_idx}.cube"
        cube_path = os.path.join(out_dir, cube_name)
        cubegen.orbital(mol, cube_path, mo_coeff[:, mo_idx],
                        nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
        pi_star.append({
            'atoms': ['ring'],
            'file': cube_name,
            'energy_ev': round(energy_ev, 2),
            'canonical_label': f"π* ({lumo_label})",
        })
        print(f"  ✓ π* ({lumo_label}, E={energy_ev:+.2f} eV) → {cube_name}")

    # ── Step 4: Localize σ/LP occupied MOs (excluding π) ──
    # Build the non-π valence occupied subspace
    non_pi_cols = [i for i in range(n_core, n_occ) if i not in pi_occ_indices]
    if non_pi_cols:
        sigma_occ = mo_coeff[:, non_pi_cols]
        print(f"  Localizing {sigma_occ.shape[1]} σ/LP valence MOs (IBO)...")
        try:
            full_occ = mo_coeff[:, :n_occ]
            iaos = lo.iao.iao(mol, full_occ)
            occ_loc = lo.ibo.ibo(mol, sigma_occ, iaos=iaos, max_iter=500)

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
                elif info['type'] in ('sigma', 'delocalized_sigma'):
                    atoms = info['atoms']
                    a1, a2 = atoms[0], atoms[1]
                    cube_name = f"{name}_sigma_{a1}_{a2}_{i}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    sigma.append({'atoms': [a1, a2], 'file': cube_name})
                    print(f"  ✓ σ({a1}–{a2}) → {cube_name}")
                elif info['type'] in ('pi', 'delocalized_pi'):
                    # Residual π that leaked through — shouldn't happen often
                    # after removing canonical π, but handle gracefully
                    atoms = info['atoms']
                    atoms_str = "_".join(atoms)
                    cube_name = f"{name}_pi_residual_{atoms_str}_{i}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    pi.append({'atoms': atoms, 'file': cube_name})
                    print(f"  ⚠ residual π({','.join(atoms)}) → {cube_name}")

        except Exception as e:
            print(f"  ⚠ σ/LP localization failed: {e}")

    # ── Step 5: Localize σ* virtual MOs (excluding π*) ──
    n_sigma_bonds = len(sigma)
    if n_sigma_bonds > 0:
        # Determine number of valence virtual orbitals
        n_min_basis = sum(5 if mol.atom_symbol(i) != 'H' else 1 for i in range(mol.natm))
        n_val_virt = max(n_sigma_bonds, n_min_basis - n_occ)
        
        # Take only the lowest valence virtual MOs to avoid diffuse/Rydberg mixing
        non_pi_virt_cols = [i for i in range(n_occ, mo_coeff.shape[1]) if i not in pi_virt_indices][:n_val_virt]

        if non_pi_virt_cols:
            virt_coeff = mo_coeff[:, non_pi_virt_cols]
            print(f"  Localizing {virt_coeff.shape[1]} σ* virtual MOs (Pipek-Mezey)...")
            try:
                loc_virt = lo.PipekMezey(mol, virt_coeff)
                loc_virt.init_guess = 'random'
                virt_loc = loc_virt.kernel()

                # We will find many virtuals; keep the most symmetric one for each bond
                best_sigmastar = {}
                for i in range(virt_loc.shape[1]):
                    mo = virt_loc[:, i]
                    pop = _atom_populations(mol, mo, ovlp)
                    info = _classify_orbital(mol, mo, pop, atom_labels, atom_ids, ovlp)

                    atoms = info.get('atoms', [])
                    if len(atoms) == 2:
                        a1, a2 = atoms[0], atoms[1]
                        pair = tuple(sorted([a1, a2]))
                        
                        # Calculate symmetry score (how close to 50/50 population)
                        abs_pop = np.abs(pop)
                        frac = abs_pop / abs_pop.sum()
                        idx1 = atom_ids.index(a1)
                        idx2 = atom_ids.index(a2)
                        sym_score = abs(frac[idx1] - frac[idx2])

                        # Keep the one that is most symmetrically shared between the 2 atoms
                        if pair not in best_sigmastar or sym_score < best_sigmastar[pair]['sym_score']:
                            best_sigmastar[pair] = {
                                'mo': mo,
                                'sym_score': sym_score,
                                'atoms': [a1, a2],
                                'idx': i
                            }

                # Now export only the best one for each bonded pair
                for pair, data in best_sigmastar.items():
                    a1, a2 = data['atoms']
                    cube_name = f"{name}_sigmastar_{a1}_{a2}_{data['idx']}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    cubegen.orbital(mol, cube_path, data['mo'], nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    sigma_star.append({'atoms': [a1, a2], 'file': cube_name})
                    print(f"  ✓ σ*({a1}–{a2}) → {cube_name} (sym: {data['sym_score']:.2f})")

            except Exception as e:
                print(f"  ⚠ σ* localization failed: {e}")

    return {
        'sigma': sigma,
        'pi': pi,
        'lone_pairs': lone_pairs,
        'sigma_star': sigma_star,
        'pi_star': pi_star,
    }
