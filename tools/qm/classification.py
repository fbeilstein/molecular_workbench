import numpy as np

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

    # Core is handled pre-localization, so if >95% on one heavy atom, it's a Valence Lone Pair
    if len(major) == 1 and major[0][1] > 0.95:
        aidx = major[0][0]
        sym = atom_labels[aidx]
        if sym != 'H':
            return {'type': 'lone_pair', 'atom': atom_ids[aidx]}

    # Bond: two+ atoms with >10% population
    if len(major) >= 2:
        sig_atoms = [m for m in major if m[1] > 0.10]
        if len(sig_atoms) > 2:
            bond_type = _classify_sigma_pi(mol, mo_coeff, sig_atoms[0][0], sig_atoms[1][0])
            if bond_type == 'pi':
                return {'type': 'delocalized_pi', 'atoms': [atom_ids[m[0]] for m in sig_atoms]}
            else:
                return {'type': 'delocalized_sigma', 'atoms': [atom_ids[m[0]] for m in sig_atoms]}
        
        if len(sig_atoms) == 2:
            a1_idx, a2_idx = sig_atoms[0][0], sig_atoms[1][0]
            bond_type = _classify_sigma_pi(mol, mo_coeff, a1_idx, a2_idx)
            return {'type': bond_type, 'atoms': [atom_ids[a1_idx], atom_ids[a2_idx]]}

    # Lone pair fallback (if only 1 atom has >10% population, or it failed bond checks)
    idx = major[0][0]
    if atom_labels[idx] == 'H':
        coords = mol.atom_coords()
        dists = np.linalg.norm(coords - coords[idx], axis=1)
        dists[idx] = 999.9
        closest_idx = np.argmin(dists)
        return {'type': 'sigma', 'atoms': [atom_ids[closest_idx], atom_ids[idx]]}
    return {'type': 'lone_pair', 'atom': atom_ids[idx]}


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
    """Distinguish σ vs π using a nodal-plane test.

    Samples the MO on a small circle perpendicular to the bond axis at
    the bond midpoint.  A σ orbital is cylindrically symmetric (same sign
    all around), while a π orbital has a nodal plane through the bond
    axis (half positive, half negative).

    This approach is basis-set agnostic — works equally well for
    s, p, d, and f orbital contributions.
    """
    coords = mol.atom_coords()

    # Hydrogens cannot form pi bonds.
    if mol.atom_symbol(atom1_idx) == 'H' or mol.atom_symbol(atom2_idx) == 'H':
        return 'sigma'

    bond_vec = coords[atom2_idx] - coords[atom1_idx]
    bond_len = np.linalg.norm(bond_vec)
    if bond_len < 1e-6:
        return 'sigma'
    bond_hat = bond_vec / bond_len

    # Build an orthogonal frame perpendicular to the bond
    ref = np.array([1., 0., 0.]) if abs(bond_hat[0]) < 0.9 else np.array([0., 1., 0.])
    perp1 = np.cross(bond_hat, ref)
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(bond_hat, perp1)

    # Sample 12 points on a circle at the bond midpoint
    midpoint = (coords[atom1_idx] + coords[atom2_idx]) / 2
    n_samples = 12
    radius = 0.5  # Bohr (~0.26 Å)
    angles = 2 * np.pi * np.arange(n_samples) / n_samples
    sample_coords = midpoint + radius * (
        np.outer(np.cos(angles), perp1) + np.outer(np.sin(angles), perp2)
    )

    ao_vals = mol.eval_gto("GTOval", sample_coords)
    mo_vals = ao_vals @ mo_coeff

    max_abs = np.max(np.abs(mo_vals))
    if max_abs < 1e-6:
        return 'sigma'

    # Among samples with significant amplitude, count positive vs negative.
    # σ → all same sign (minority ≈ 0), π → half-and-half (minority ≈ 0.5).
    sig_mask = np.abs(mo_vals) > max_abs * 0.1
    n_sig = int(np.sum(sig_mask))
    if n_sig < 4:
        return 'sigma'

    n_pos = int(np.sum(mo_vals[sig_mask] > 0))
    n_neg = n_sig - n_pos
    minority_frac = min(n_pos, n_neg) / n_sig

    return 'pi' if minority_frac > 0.25 else 'sigma'


