#!/usr/bin/env python3
"""
mol_orbitals — Compute molecular orbitals and generate cube files

Takes an optimized .xyz file and computes:
  - HOMO/LUMO canonical orbitals
  - Localized σ bonds, π bonds, and lone pairs (Pipek-Mezey)
  - orbitals.json manifest describing all generated cubes

Uses PySCF for quantum chemistry calculations.

Usage:
    python mol_orbitals.py molecule.xyz --name aspirin -o output/
    python mol_orbitals.py molecule.xyz --charge -1 --basis cc-pvdz
"""

import argparse
import json
import os
import sys

import numpy as np


# ── XYZ parsing ──────────────────────────────────────────────────────────────

def parse_xyz(filepath):
    """Parse XYZ file → PySCF atom string + atom labels list."""
    with open(filepath, 'r') as f:
        lines = f.readlines()
    n_atoms = int(lines[0].strip())
    atoms = []
    labels = []
    for line in lines[2:2 + n_atoms]:
        parts = line.split()
        if len(parts) >= 4:
            sym = parts[0]
            x, y, z = float(parts[1]), float(parts[2]), float(parts[3])
            atoms.append(f'{sym} {x} {y} {z}')
            labels.append(sym)
    return '; '.join(atoms), labels


def make_atom_ids(atom_labels):
    """Create unique atom IDs: ['C','C','O','H'] → ['C1','C2','O1','H1']."""
    ids = []
    count = {}
    for sym in atom_labels:
        count[sym] = count.get(sym, 0) + 1
        ids.append(f"{sym}{count[sym]}")
    return ids


# ── SCF calculation ──────────────────────────────────────────────────────────

def run_scf(xyz_file, charge=0, spin=0, basis='6-31g*', method='b3lyp',
            chkfile=None, read_chk=False):
    """Run DFT/HF calculation, return (mol, mf, atom_labels)."""
    from pyscf import gto, scf, dft

    atom_str, labels = parse_xyz(xyz_file)
    n_atoms = len(labels)
    print(f"  SCF: {n_atoms} atoms, basis={basis}")

    # Auto-detect odd electron count before building molecule
    _Z = {'H':1,'He':2,'Li':3,'Be':4,'B':5,'C':6,'N':7,'O':8,'F':9,'Ne':10,
          'Na':11,'Mg':12,'Al':13,'Si':14,'P':15,'S':16,'Cl':17,'Ar':18,
          'K':19,'Ca':20,'Br':35,'I':53}
    total_z = sum(_Z.get(s, 0) for s in labels)
    n_elec_check = total_z - charge
    if n_elec_check % 2 != 0 and spin == 0:
        spin = 1  # doublet for odd-electron systems
        print(f"  ⚠ Odd electron count ({n_elec_check}), auto-setting spin=1 (doublet)")

    mol = gto.M(atom=atom_str, basis=basis, charge=charge, spin=spin, verbose=0)
    n_elec = mol.nelectron
    print(f"  {n_elec} electrons (charge={charge}, spin={spin})")

    if spin > 0:
        # Open-shell: use unrestricted methods
        if method.lower() in ('hf', 'rhf', 'uhf'):
            mf = scf.UHF(mol)
        else:
            mf = dft.UKS(mol)
            mf.xc = method
        print(f"  Using unrestricted (open-shell) method")
    else:
        if method.lower() in ('hf', 'rhf'):
            mf = scf.RHF(mol)
        else:
            mf = dft.RKS(mol)
            mf.xc = method

    if chkfile:
        mf.chkfile = chkfile
        if read_chk and os.path.exists(chkfile):
            mf.init_guess = 'chkfile'

    mf.verbose = 0
    print(f"  Running {method.upper()} calculation...")
    energy = mf.kernel()

    if not mf.converged:
        print(f"  WARNING: SCF did not converge!", file=sys.stderr)

    print(f"  Total energy: {energy:.6f} Hartree")
    return mol, mf, labels


# ── Ring detection helpers ────────────────────────────────────────────────────

_COVALENT_R_BOHR = {  # covalent radii in Bohr
    'H': 0.586, 'C': 1.436, 'N': 1.341, 'O': 1.247, 'F': 1.077,
    'S': 1.984, 'P': 2.022, 'Cl': 1.928, 'Br': 2.268, 'I': 2.627,
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


# ── Canonical orbitals (HOMO/LUMO + π system) ────────────────────────────────

def compute_canonical(mol, mf, atom_labels, name, out_dir, grid_points=50):
    """Generate HOMO, LUMO, and canonical π system cube files.

    Instead of blindly exporting HOMO-N, detects which canonical MOs are
    π-type using PCA-based ring plane analysis and exports only those.
    """
    from pyscf import tools
    import numpy as np

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
        tools.cubegen.orbital(mol, cube_path, mo_coeff[:, idx],
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
        tools.cubegen.orbital(mol, cube_path, mo_coeff[:, mo_idx],
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


# ── Localized orbitals (σ/π/LP) ──────────────────────────────────────────────

def compute_localized(mol, mf, atom_labels, name, out_dir, grid_points=50):
    """IBO (Intrinsic Bond Orbital) localization → classify σ/π/LP, write cubes.

    IBO localization maximizes atomic orbital populations, naturally preserving
    the σ/π separation critical for aromatic systems.  Unlike Boys localization
    (which produces "banana bonds" by mixing σ and π), IBO keeps pi orbitals
    as pure p_z contributions, making classification straightforward.
    """
    from pyscf import lo, tools
    import numpy as np

    # Handle open-shell (UHF/UKS)
    mo_coeff = mf.mo_coeff
    if isinstance(mo_coeff, (tuple, list)) or (isinstance(mo_coeff, np.ndarray) and mo_coeff.ndim == 3):
        mo_coeff = mo_coeff[0]  # Use alpha MOs
        n_occ = mol.nelec[0]
    else:
        n_occ = mol.nelectron // 2

    # Freeze core orbitals by stripping lowest energy MOs
    core_electrons = 0
    for sym in atom_labels:
        core_electrons += _CORE_ELECTRONS.get(sym, 0)
    n_core = core_electrons // 2

    # Localize ONLY valence OCCUPIED orbitals
    occ_coeff = mo_coeff[:, n_core:n_occ]

    print(f"  Localizing {occ_coeff.shape[1]} valence MOs (IBO)...")
    loc_coeff = lo.ibo.ibo(mol, occ_coeff, verbose=0)
    print(f"  Localization converged")

    sigma, pi, lone_pairs = [], [], []
    ovlp = mol.intor_symmetric('int1e_ovlp')
    atom_ids = make_atom_ids(atom_labels)
    lp_count = {}

    for i in range(loc_coeff.shape[1]):
        mo = loc_coeff[:, i]
        pop = _atom_populations(mol, mo, ovlp)
        info = _classify_orbital(mol, mo, pop, atom_labels, atom_ids, ovlp)

        if info['type'] == 'core':
            continue

        elif info['type'] == 'lone_pair':
            atom = info['atom']
            lp_count[atom] = lp_count.get(atom, 0) + 1
            idx = lp_count[atom]
            cube_name = f"{name}_lp_{atom}_{idx}_{i}.cube"
            cube_path = os.path.join(out_dir, cube_name)
            tools.cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
            lone_pairs.append({'atom': atom, 'index': idx, 'file': cube_name})
            print(f"  ✓ LP({atom} #{idx}) → {cube_name}")

        elif info['type'] == 'sigma':
            a1, a2 = info['atoms']
            cube_name = f"{name}_sigma_{a1}_{a2}_{i}.cube"
            cube_path = os.path.join(out_dir, cube_name)
            tools.cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
            sigma.append({'atoms': [a1, a2], 'file': cube_name})
            print(f"  ✓ σ({a1}–{a2}) → {cube_name}")

        elif info['type'] == 'pi':
            a1, a2 = info['atoms']
            cube_name = f"{name}_pi_{a1}_{a2}_{i}.cube"
            cube_path = os.path.join(out_dir, cube_name)
            tools.cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
            pi.append({'atoms': [a1, a2], 'file': cube_name})
            print(f"  ✓ π({a1}={a2}) → {cube_name}")

    print(f"  Localized: {len(sigma)} σ, {len(pi)} π, {len(lone_pairs)} LP")

    # -- Localize Virtuals (antibonding) --
    sigma_star, pi_star = [], []
    n_virt_to_loc = len(sigma) + len(pi)
    if n_virt_to_loc > 0 and n_occ + n_virt_to_loc <= mo_coeff.shape[1]:
        virt_coeff = mo_coeff[:, n_occ : n_occ + n_virt_to_loc]
        print(f"  Localizing {virt_coeff.shape[1]} lowest virtual MOs (Boys)...")
        try:
            virt_loc = lo.Boys(mol, virt_coeff).kernel()
            for i in range(virt_loc.shape[1]):
                mo = virt_loc[:, i]
                pop = _atom_populations(mol, mo, ovlp)
                info = _classify_orbital(mol, mo, pop, atom_labels, atom_ids, ovlp)
                
                if info['type'] == 'sigma':
                    a1, a2 = info['atoms']
                    cube_name = f"{name}_sigmastar_{a1}_{a2}_{i}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    tools.cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    sigma_star.append({'atoms': [a1, a2], 'file': cube_name})
                    print(f"  ✓ σ*({a1}–{a2}) → {cube_name}")
                elif info['type'] == 'pi':
                    a1, a2 = info['atoms']
                    cube_name = f"{name}_pistar_{a1}_{a2}_{i}.cube"
                    cube_path = os.path.join(out_dir, cube_name)
                    tools.cubegen.orbital(mol, cube_path, mo, nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)
                    pi_star.append({'atoms': [a1, a2], 'file': cube_name})
                    print(f"  ✓ π*({a1}–{a2}) → {cube_name}")
        except Exception as e:
            print(f"  ⚠ Virtual localization failed: {e}")

    return {
        'sigma': sigma,
        'pi': pi,
        'lone_pairs': lone_pairs,
        'sigma_star': sigma_star,
        'pi_star': pi_star
    }



# ── Orbital classification helpers ───────────────────────────────────────────

def _atom_populations(mol, mo_coeff, ovlp):
    """Mulliken population of a single MO on each atom."""
    dm_mo = np.outer(mo_coeff, mo_coeff)
    ps = dm_mo * ovlp
    ao_labels = mol.ao_labels(fmt=False)
    pop = np.zeros(mol.natm)
    for mu, (atom_idx, *_) in enumerate(ao_labels):
        pop[atom_idx] += ps[mu, :].sum()
    return pop


def _classify_orbital(mol, mo_coeff, pop, atom_labels, atom_ids, ovlp):
    """Classify a localized MO as lone_pair, sigma, or pi."""
    total = pop.sum()
    if total < 1e-6:
        return {'type': 'core'}

    frac = pop / total
    major = [(i, frac[i]) for i in range(mol.natm) if frac[i] > 0.04]
    major.sort(key=lambda x: -x[1])

    # Core is handled pre-localization, so if >95% on one heavy atom, it's a Valence Lone Pair
    if len(major) == 1 and major[0][1] > 0.95:
        aidx = major[0][0]
        sym = atom_labels[aidx]
        if sym != 'H':
            return {'type': 'lone_pair', 'atom': atom_ids[aidx]}

    # Bond: two+ atoms with >10% population
    if len(major) >= 2:
        a1_idx, a2_idx = major[0][0], major[1][0]
        bond_type = _classify_sigma_pi(mol, mo_coeff, a1_idx, a2_idx)
        return {'type': bond_type, 'atoms': [atom_ids[a1_idx], atom_ids[a2_idx]]}

    # Lone pair fallback
    if len(major) == 1:
        idx = major[0][0]
        if atom_labels[idx] == 'H':
            import numpy as np
            coords = mol.atom_coords()
            dists = np.linalg.norm(coords - coords[idx], axis=1)
            dists[idx] = 999.9
            closest_idx = np.argmin(dists)
            return {'type': 'sigma', 'atoms': [atom_ids[closest_idx], atom_ids[idx]]}
        return {'type': 'lone_pair', 'atom': atom_ids[idx]}

    return {'type': 'sigma', 'atoms': ['?', '?']}


# Core electron counts per element (number of core ELECTRONS, not MOs)
# Core = all electrons below the valence shell
_CORE_ELECTRONS = {
    'H': 0, 'He': 0,
    'Li': 2, 'Be': 2, 'B': 2, 'C': 2, 'N': 2, 'O': 2, 'F': 2, 'Ne': 2,
    'Na': 10, 'Mg': 10, 'Al': 10, 'Si': 10, 'P': 10, 'S': 10, 'Cl': 10, 'Ar': 10,
    'K': 18, 'Ca': 18,
    'Sc': 18, 'Ti': 18, 'V': 18, 'Cr': 18, 'Mn': 18, 'Fe': 18, 'Co': 18, 'Ni': 18,
    'Cu': 18, 'Zn': 18,
    'Ga': 28, 'Ge': 28, 'As': 28, 'Se': 28, 'Br': 28, 'Kr': 28,
    'Rb': 36, 'Sr': 36,
    'I': 46, 'Xe': 46,
}

def _is_core_orbital(mol, mo_coeff, atom_idx, element):
    """Check if a localized MO on this atom is a core (inner-shell) orbital
    by computing its orbital energy proxy (kinetic energy contribution)."""
    core_e = _CORE_ELECTRONS.get(element, 0)
    if core_e == 0:
        return False

    # Check if the orbital is dominated by low-n AOs (core-like)
    ao_labels = mol.ao_labels(fmt=False)
    core_ao_weight = 0.0
    valence_ao_weight = 0.0

    # Determine valence shell n for this element
    # Row 1: n=1, Row 2: n=2, Row 3: n=3, Row 4: n=4 etc
    _VALENCE_N = {
        'H': 1, 'He': 1,
        'Li': 2, 'Be': 2, 'B': 2, 'C': 2, 'N': 2, 'O': 2, 'F': 2, 'Ne': 2,
        'Na': 3, 'Mg': 3, 'Al': 3, 'Si': 3, 'P': 3, 'S': 3, 'Cl': 3, 'Ar': 3,
        'K': 4, 'Ca': 4, 'Sc': 4, 'Ti': 4, 'V': 4, 'Cr': 4, 'Mn': 4, 'Fe': 4,
        'Co': 4, 'Ni': 4, 'Cu': 4, 'Zn': 4,
        'Ga': 4, 'Ge': 4, 'As': 4, 'Se': 4, 'Br': 4, 'Kr': 4,
        'Rb': 5, 'Sr': 5, 'I': 5, 'Xe': 5,
    }
    val_n = _VALENCE_N.get(element, 1)

    for mu, (ai, *rest) in enumerate(ao_labels):
        if ai != atom_idx:
            continue
        c2 = mo_coeff[mu] ** 2
        # AO label format: (atom_idx, element, shell_label, sublabel)
        # shell_label is '1s', '2s', '2p', '3s', '3p', '3d', '4s', etc.
        shell_label = rest[1] if len(rest) > 1 else ''
        n = 0
        for ch in shell_label:
            if ch.isdigit():
                n = int(ch)
                break
        if n > 0 and n < val_n:
            core_ao_weight += c2
        else:
            valence_ao_weight += c2

    total_w = core_ao_weight + valence_ao_weight
    if total_w < 1e-10:
        return False

    return core_ao_weight / total_w > 0.5


def _classify_sigma_pi(mol, mo_coeff, atom1_idx, atom2_idx):
    """Distinguish σ vs π by AO angular momentum analysis."""
    coords = mol.atom_coords()
    
    # Hydrogens cannot form pi bonds. Prevent polarization functions from falsely triggering pi character.
    sym1 = mol.atom_symbol(atom1_idx)
    sym2 = mol.atom_symbol(atom2_idx)
    if sym1 == 'H' or sym2 == 'H':
        return 'sigma'
        
    bond_vec = coords[atom2_idx] - coords[atom1_idx]
    bond_len = np.linalg.norm(bond_vec)
    if bond_len < 1e-6:
        return 'sigma'
    bond_hat = bond_vec / bond_len

    ao_labels = mol.ao_labels(fmt=False)
    p_dirs = {'x': np.array([1, 0, 0]),
              'y': np.array([0, 1, 0]),
              'z': np.array([0, 0, 1])}

    p_along, p_perp = 0.0, 0.0
    for mu, (aidx, *rest) in enumerate(ao_labels):
        if aidx not in (atom1_idx, atom2_idx):
            continue
        sublabel = rest[-1] if rest else ''
        coeff = abs(mo_coeff[mu])
        if coeff < 1e-4:
            continue
        if sublabel in p_dirs:
            proj = abs(np.dot(p_dirs[sublabel], bond_hat))
            p_along += coeff * proj
            p_perp += coeff * (1 - proj)

    return 'pi' if p_perp > p_along * 1.5 else 'sigma'


# ── Main pipeline ────────────────────────────────────────────────────────────

def compute_all_orbitals(xyz_file, name=None, out_dir=None, charge=0, spin=0,
                         basis='6-31g*', method='b3lyp', grid_points=50,
                         chkfile=None, read_chk=False):
    """Full pipeline: SCF → canonical + localized orbitals → cube files + JSON."""
    if name is None:
        name = os.path.splitext(os.path.basename(xyz_file))[0]
    if out_dir is None:
        out_dir = os.path.dirname(xyz_file) or '.'
    os.makedirs(out_dir, exist_ok=True)

    # Run SCF
    mol, mf, labels = run_scf(xyz_file, charge, spin, basis, method,
                               chkfile, read_chk)

    # Canonical orbitals
    canonical = compute_canonical(mol, mf, labels, name, out_dir, grid_points)

    # Localized orbitals
    localized = compute_localized(mol, mf, labels, name, out_dir, grid_points)

    # ESP surface
    esp_data = compute_esp_surface(mol, mf, name, out_dir, grid_points)

    # Write manifest
    manifest = {
        'name': name,
        'method': method,
        'basis': basis,
        'charge': charge,
        'spin': spin,
        'n_electrons': mol.nelectron,
        'energy_hartree': float(mf.e_tot),
        'canonical': canonical,
        'localized': localized,
    }
    if esp_data:
        manifest['esp_surface'] = esp_data

    manifest_path = os.path.join(out_dir, f"{name}_orbitals.json")
    with open(manifest_path, 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"  Manifest: {manifest_path}")

    return manifest


# ── Electrostatic Potential Surface ──────────────────────────────────────────

ANG2BOHR = 1.8897259886

def compute_esp_surface(mol, mf, name, out_dir, grid_points=50, iso_density=0.01):
    """Compute molecular electrostatic potential mapped onto electron density isosurface.

    Steps:
      1. Compute total electron density on a 3D grid
      2. Extract isosurface at iso_density using marching cubes
      3. Evaluate the electrostatic potential at each surface vertex
      4. Save mesh + per-vertex ESP values as compact JSON

    The frontend can then color the surface using a red↔blue diverging map.
    """
    from pyscf import gto, tools
    from skimage import measure
    import trimesh

    print(f"  Computing ESP surface (isodensity={iso_density})...")

    # Build the density matrix
    dm = mf.make_rdm1()
    if isinstance(dm, (tuple, list)):
        dm = dm[0] + dm[1]  # alpha + beta for open-shell

    # Generate electron density cube
    dens_cube = os.path.join(out_dir, f"{name}_density.cube")
    tools.cubegen.density(mol, dens_cube, dm,
                          nx=grid_points, ny=grid_points, nz=grid_points, margin=5.0)

    # Read the cube file
    with open(dens_cube, 'r') as f:
        lines = f.readlines()

    natoms_line = lines[2].split()
    natoms = int(natoms_line[0])
    origin = np.array([float(x) for x in natoms_line[1:4]])

    nx_line = lines[3].split()
    ny_line = lines[4].split()
    nz_line = lines[5].split()
    npts = (int(nx_line[0]), int(ny_line[0]), int(nz_line[0]))
    step = np.array([float(nx_line[1]), float(ny_line[2]), float(nz_line[3])])

    vals = []
    for line in lines[6 + natoms:]:
        vals.extend([float(x) for x in line.split()])
    vol = np.array(vals).reshape(npts)

    # Clean up cube file
    os.remove(dens_cube)

    # Extract isosurface
    try:
        verts, faces, normals, _ = measure.marching_cubes(vol, iso_density)
    except Exception as e:
        print(f"  ✗ ESP surface extraction failed: {e}")
        return None

    # Clean with trimesh
    mesh = trimesh.Trimesh(vertices=verts, faces=faces, vertex_normals=normals, process=True)
    verts = mesh.vertices
    faces = mesh.faces
    normals = mesh.vertex_normals

    # Convert from grid indices to Bohr, then to Angstrom
    verts_bohr = verts * step + origin
    verts_ang = verts_bohr / ANG2BOHR

    # Evaluate ESP at each surface vertex (in Bohr coordinates)
    print(f"  Evaluating ESP at {len(verts)} surface vertices...")

    # Nuclear contribution (vectorized)
    nuc_charges = mol.atom_charges()
    nuc_coords = mol.atom_coords()  # Bohr

    esp_values = np.zeros(len(verts_bohr))
    for a in range(mol.natm):
        dists = np.linalg.norm(verts_bohr - nuc_coords[a], axis=1)
        dists = np.maximum(dists, 1e-6)  # avoid division by zero
        esp_values += nuc_charges[a] / dists

    # Electronic contribution: V_elec(r) = -Tr[D * (1/|r-R|)]
    # Using PySCF's with_rinv_origin + int1e_rinv integral
    print(f"  Evaluating electronic ESP...")
    for i, coord in enumerate(verts_bohr):
        with mol.with_rinv_origin(coord):
            ints = mol.intor('int1e_rinv')
        esp_values[i] -= np.einsum('ij,ij->', dm, ints)

    # Convert ESP from Hartree/e to kcal/mol for chemical intuition
    HARTREE2KCAL = 627.509
    esp_kcal = esp_values * HARTREE2KCAL

    # Build compact JSON mesh
    v_list = [{'x': round(float(v[0]), 4), 'y': round(float(v[1]), 4), 'z': round(float(v[2]), 4)}
              for v in verts_ang]
    f_list = [int(idx) for face in faces for idx in face]
    n_list = [{'x': -round(float(n[0]), 4), 'y': -round(float(n[1]), 4), 'z': -round(float(n[2]), 4)}
              for n in normals]
    esp_list = [round(float(v), 2) for v in esp_kcal]

    esp_data = {
        'vertices': v_list,
        'faces': f_list,
        'normals': n_list,
        'esp_values': esp_list,  # kcal/mol per vertex
        'esp_min': round(float(esp_kcal.min()), 2),
        'esp_max': round(float(esp_kcal.max()), 2),
    }

    esp_path = os.path.join(out_dir, f"{name}_esp.json")
    with open(esp_path, 'w') as f:
        json.dump(esp_data, f, separators=(',', ':'))

    sz = os.path.getsize(esp_path) / 1024
    print(f"  ✓ ESP surface: {len(verts)} vertices, ESP range [{esp_kcal.min():.1f}, {esp_kcal.max():.1f}] kcal/mol → {name}_esp.json ({sz:.0f} KB)")

    return {'file': f"{name}_esp.json", 'esp_min': esp_data['esp_min'], 'esp_max': esp_data['esp_max']}


# ── CLI ──────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description='Compute molecular orbitals (HOMO/LUMO/σ/π/LP) and generate cube files',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s water.xyz --name water
  %(prog)s aspirin.xyz --charge 0 --basis 6-31g* -o output/
  %(prog)s chloride.xyz --charge -1 --name Cl
        """
    )
    parser.add_argument('xyz_file', help='Input .xyz file')
    parser.add_argument('--name', '-n', default=None, help='Base name (default: from filename)')
    parser.add_argument('-o', '--output-dir', default=None, help='Output directory')
    parser.add_argument('--charge', type=int, default=0, help='Molecular charge')
    parser.add_argument('--spin', type=int, default=0, help='Spin (2S, e.g. 0=singlet)')
    parser.add_argument('--basis', default='6-31g*', help='Basis set')
    parser.add_argument('--method', default='b3lyp', help='DFT functional or hf')
    parser.add_argument('--grid', type=int, default=50, help='Cube grid points per dimension')
    parser.add_argument('--chkfile', default=None, help='PySCF checkpoint file')
    parser.add_argument('--read-chk', action='store_true', help='Read guess from chkfile')

    args = parser.parse_args()
    if not os.path.exists(args.xyz_file):
        print(f"ERROR: {args.xyz_file} not found", file=sys.stderr)
        sys.exit(1)

    compute_all_orbitals(args.xyz_file, args.name, args.output_dir,
                         args.charge, args.spin, args.basis, args.method,
                         args.grid, args.chkfile, args.read_chk)


if __name__ == '__main__':
    main()
