import numpy as np
from qm.constants import (
    BOHR_TO_ANG,
    CORE_ELECTRONS,
    HARTREE_TO_EV,
    LP_DOMINANT_POPULATION_THRESHOLD,
    LP_HYPERCONJUGATION_PRIMARY_MIN,
    LP_HYPERCONJUGATION_SECONDARY_MAX,
    NODAL_AMPLITUDE_FLOOR,
    NODAL_MIN_SIGNIFICANT_POINTS,
    NODAL_SCAN_RADII_BOHR,
    NODAL_SCAN_SAMPLE_POINTS,
    NODAL_SIGNIFICANT_AMP_RATIO,
    SP3_MAX_NEIGHBOR_DIST_ANG,
    TWO_CENTER_BOND_MIN_POPULATION,
)

_CORE_ELECTRONS = CORE_ELECTRONS

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
    abs_pop = np.abs(pop)
    total = abs_pop.sum()
    if total < 1e-6:
        return {'type': 'core'}

    frac = abs_pop / total
    major = [(i, frac[i]) for i in range(mol.natm) if frac[i] > 0.04]
    major.sort(key=lambda x: -x[1])

    # Valence Lone Pair:
    # Lone pairs on electronegative heteroatoms (O, N, F, Cl, Br, I, S, P) frequently participate
    # in negative hyperconjugation (e.g. oxyanion n_O -> sigma*(C-Cl) or sigma*(C-O)), yielding
    # 75-90% on the heteroatom and 10-25% on adjacent atoms.
    # If >=75% is on a single non-H atom, or if the primary atom holds >=68% with >=3:1 ratio
    # over the secondary atom, this is a Valence Lone Pair (non-bonding), NOT a covalent 2-center bond.
    if len(major) >= 1:
        top_idx, top_frac = major[0]
        sym = atom_labels[top_idx]
        if sym != 'H':
            is_lp = False
            second_frac = major[1][1] if len(major) >= 2 else 0.0
            if top_frac >= LP_DOMINANT_POPULATION_THRESHOLD:
                is_lp = True
            elif top_frac >= LP_HYPERCONJUGATION_PRIMARY_MIN and second_frac < LP_HYPERCONJUGATION_SECONDARY_MAX:
                is_lp = True
            elif len(major) == 1:
                is_lp = True

            if is_lp:
                return {'type': 'lone_pair', 'atom': atom_ids[top_idx]}

    # Bond: two+ atoms with >10% population
    if len(major) >= 2:
        sig_atoms = [m for m in major if m[1] > TWO_CENTER_BOND_MIN_POPULATION]
        if len(sig_atoms) > 2:
            # Find a genuinely bonded pair among the significant atoms to determine sigma/pi symmetry
            coords = mol.atom_coords()
            a1_idx, a2_idx = sig_atoms[0][0], sig_atoms[1][0]
            for i in range(len(sig_atoms)):
                for j in range(i + 1, len(sig_atoms)):
                    idx_i, idx_j = sig_atoms[i][0], sig_atoms[j][0]
                    dist = np.linalg.norm(coords[idx_i] - coords[idx_j])
                    # Ensure atoms are bonded (ignoring H-H since they can't form pi anyway)
                    # Typical C-C bond is ~1.2-1.6 Angstroms. We use 1.8 as a safe upper bound.
                    if 0.5 < dist < 1.8:
                        a1_idx, a2_idx = idx_i, idx_j
                        break
                else:
                    continue
                break

            bond_type = _classify_sigma_pi(mol, mo_coeff, a1_idx, a2_idx)
            if bond_type == 'pi':
                return {'type': 'delocalized_pi', 'atoms': [atom_ids[m[0]] for m in sig_atoms]}
            else:
                return {'type': 'delocalized_sigma', 'atoms': [atom_ids[m[0]] for m in sig_atoms]}
        
        if len(sig_atoms) == 2:
            a1_idx, a2_idx = sig_atoms[0][0], sig_atoms[1][0]
            bond_type = _classify_sigma_pi(mol, mo_coeff, a1_idx, a2_idx)
            return {'type': bond_type, 'atoms': [atom_ids[a1_idx], atom_ids[a2_idx]]}

    # Fallback (if failed bond checks)
    idx = major[0][0]
    if atom_labels[idx] == 'H':
        coords = mol.atom_coords()
        dists = np.linalg.norm(coords - coords[idx], axis=1)
        dists[idx] = 999.9
        closest_idx = np.argmin(dists)
        return {'type': 'sigma', 'atoms': [atom_ids[closest_idx], atom_ids[idx]]}
    return {'type': 'lone_pair', 'atom': atom_ids[idx]}





def count_core_orbitals(mol, mo_energy: np.ndarray, nocc: int, cutoff_ev: float = -50.0) -> int:
    """Determine number of closed core orbitals combining energy threshold and noble-gas core table."""
    if nocc <= 0:
        return 0
    n_core_energy = int(np.sum((mo_energy[:nocc] * HARTREE_TO_EV) < cutoff_ev))
    n_core_table = sum(_CORE_ELECTRONS.get(mol.atom_symbol(i), 0) for i in range(mol.natm)) // 2
    return min(int(nocc), max(n_core_energy, n_core_table))



def _classify_sigma_pi(mol, mo_coeff, atom1_idx, atom2_idx):
    """Distinguish σ vs π using a nodal-plane test.

    Samples the MO on concentric circles perpendicular to the bond axis at
    the bond midpoint. A σ orbital is cylindrically symmetric (same sign
    all around), while a π orbital has a nodal plane through the bond
    axis (opposite signs on diametrically opposed points).
    """
    coords_bohr = mol.atom_coords()
    coords_ang = coords_bohr * BOHR_TO_ANG

    # Hydrogens cannot form pi bonds.
    if mol.atom_symbol(atom1_idx) == 'H' or mol.atom_symbol(atom2_idx) == 'H':
        return 'sigma'

    # Tetrahedral / 4-coordinate atoms (like sp3 carbon) cannot form pi bonds without violating octet rules.
    for a_idx in (atom1_idx, atom2_idx):
        if mol.atom_symbol(a_idx) in ('C', 'Si'):
            dists_ang = np.linalg.norm(coords_ang - coords_ang[a_idx], axis=1)
            # count bonded neighbors within 1.95 A (excluding self)
            n_neighbors = int(np.sum((dists_ang > 0.2) & (dists_ang < SP3_MAX_NEIGHBOR_DIST_ANG)))
            if n_neighbors >= 4:
                return 'sigma'

    bond_vec = coords_bohr[atom2_idx] - coords_bohr[atom1_idx]
    bond_len = np.linalg.norm(bond_vec)
    if bond_len < 1e-6:
        return 'sigma'
    bond_hat = bond_vec / bond_len

    # Build an orthogonal frame perpendicular to the bond
    ref = np.array([1., 0., 0.]) if abs(bond_hat[0]) < 0.9 else np.array([0., 1., 0.])
    perp1 = np.cross(bond_hat, ref)
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(bond_hat, perp1)

    midpoint = (coords_bohr[atom1_idx] + coords_bohr[atom2_idx]) / 2.0
    angles = 2 * np.pi * np.arange(NODAL_SCAN_SAMPLE_POINTS) / NODAL_SCAN_SAMPLE_POINTS
    radii = NODAL_SCAN_RADII_BOHR
    is_pi = False

    for radius in radii:
        sample_coords = midpoint + radius * (
            np.outer(np.cos(angles), perp1) + np.outer(np.sin(angles), perp2)
        )
        ao_vals = mol.eval_gto("GTOval", sample_coords)
        mo_vals = ao_vals @ mo_coeff
        max_abs = np.max(np.abs(mo_vals))
        if max_abs < NODAL_AMPLITUDE_FLOOR:
            continue

        sig_mask = np.abs(mo_vals) > max_abs * NODAL_SIGNIFICANT_AMP_RATIO
        n_sig = int(np.sum(sig_mask))
        if n_sig < NODAL_MIN_SIGNIFICANT_POINTS:
            continue

        n_pos = int(np.sum(mo_vals[sig_mask] > 0))
        n_neg = n_sig - n_pos
        minority_frac = min(n_pos, n_neg) / n_sig

        opposite_inversions = 0
        opposite_tested = 0
        for j in range(6):
            if sig_mask[j] and sig_mask[j + 6]:
                opposite_tested += 1
                if mo_vals[j] * mo_vals[j + 6] < 0:
                    opposite_inversions += 1

        if minority_frac > 0.25 and (opposite_tested >= 3 and opposite_inversions >= 3):
            is_pi = True
            break
        elif minority_frac > 0.35 and opposite_inversions >= 2:
            is_pi = True
            break

    return 'pi' if is_pi else 'sigma'


