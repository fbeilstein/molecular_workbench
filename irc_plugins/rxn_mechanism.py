"""
Reaction mechanism analysis: driving coordinates, push-arrow parser,
and bond change detection.
"""

import json
import math
import os
from .rxn_mapping import _get_rdkit


def _dist(p1, p2):
    return math.sqrt((p1['x'] - p2[0]) ** 2 + (p1['y'] - p2[1]) ** 2)


def _try_kekulize(mol, Chem):
    try:
        Chem.Kekulize(mol, clearAromaticFlags=True)
    except Exception:
        pass


def _match_arrow_to_features(arrow, mols):
    """Match arrow head and tail to nearest atom or bond midpoints."""
    start = arrow['data']['pos'][0]
    end = arrow['data']['pos'][1]
    start_map = None
    start_is_bond = False
    min_ds = 999.0
    end_map = None
    min_de = 999.0

    for mol in mols.values():
        if mol.get('type') != 'molecule':
            continue

        for atom in mol.get('atoms', []):
            loc = atom['location']
            ds, de = _dist(start, loc), _dist(end, loc)
            if ds < min_ds:
                min_ds, start_map, start_is_bond = ds, atom.get("mapping"), False
            if de < min_de:
                min_de, end_map = de, atom.get("mapping")

        for bond in mol.get('bonds', []):
            atoms = mol['atoms']
            a1, a2 = atoms[bond['atoms'][0]]['location'], atoms[bond['atoms'][1]]['location']
            mid = [(a1[0] + a2[0]) / 2, (a1[1] + a2[1]) / 2]
            ds = _dist(start, mid)
            if ds < min_ds:
                min_ds = ds
                start_map = (atoms[bond['atoms'][0]].get('mapping'), atoms[bond['atoms'][1]].get('mapping'))
                start_is_bond = True

    return start_map, start_is_bond, end_map


def parse_ket_mechanism(ket_path):
    """Parse curved push-arrows from Ketcher KET file to identify formed and broken bonds."""
    with open(ket_path) as f:
        ket = json.load(f)

    arrows = [
        node for node in ket.get('root', {}).get('nodes', [])
        if node.get('type') == 'arrow' and 'elliptical-arc-arrow' in node.get('data', {}).get('mode', '')
    ]
    if not arrows:
        return None

    mols = {k: v for k, v in ket.items() if k.startswith('mol')}
    formed_bonds, broken_bonds = [], []

    for a in arrows:
        start_map, start_is_bond, end_map = _match_arrow_to_features(a, mols)
        if start_is_bond and isinstance(start_map, tuple):
            if start_map[0] is not None and start_map[1] is not None:
                broken_bonds.append(sorted([start_map[0], start_map[1]]))
            continue
        if not start_is_bond and start_map is not None and end_map is not None:
            formed_bonds.append(sorted([start_map, end_map]))

    return {
        'formed_bonds': formed_bonds,
        'broken_bonds': broken_bonds
    }


def _extract_ket_coords(ket_mech, map_num_to_r, r_mols):
    """Convert KET map numbers to fragment atom indices."""
    formed_bonds, formed_labels = [], []
    for b1, b2 in ket_mech['formed_bonds']:
        r1, r2 = map_num_to_r.get(b1), map_num_to_r.get(b2)
        if not r1 or not r2:
            continue
        formed_bonds.append((r1, r2))
        sym1 = r_mols[r1[0]].GetAtomWithIdx(r1[1]).GetSymbol()
        sym2 = r_mols[r2[0]].GetAtomWithIdx(r2[1]).GetSymbol()
        formed_labels.append(f"{sym1}(R{r1[0]+1}:{r1[1]})–{sym2}(R{r2[0]+1}:{r2[1]})")

    broken_bonds, broken_labels = [], []
    for b1, b2 in ket_mech['broken_bonds']:
        r1, r2 = map_num_to_r.get(b1), map_num_to_r.get(b2)
        if not r1 or not r2:
            continue
        broken_bonds.append((r1, r2))
        sym1 = r_mols[r1[0]].GetAtomWithIdx(r1[1]).GetSymbol()
        sym2 = r_mols[r2[0]].GetAtomWithIdx(r2[1]).GetSymbol()
        broken_labels.append(f"{sym1}(R{r1[0]+1}:{r1[1]})–{sym2}(R{r2[0]+1}:{r2[1]})")

    pt_labels, proton_transfers = [], []
    for fb in formed_bonds:
        for bb in broken_bonds:
            shared = set(fb).intersection(bb)
            if not shared:
                continue
            shared_atom = list(shared)[0]
            donor = bb[0] if bb[1] == shared_atom else bb[1]
            acceptor = fb[0] if fb[1] == shared_atom else fb[1]
            sym1 = r_mols[donor[0]].GetAtomWithIdx(donor[1]).GetSymbol()
            sym2 = r_mols[acceptor[0]].GetAtomWithIdx(acceptor[1]).GetSymbol()
            pt_labels.append(f"Proton Hop via {sym1}–{sym2}")
            proton_transfers.append((donor, acceptor))

    return {
        'formed_bonds': formed_bonds,
        'broken_bonds': broken_bonds,
        'proton_transfers': proton_transfers,
        'formed_labels': formed_labels,
        'broken_labels': broken_labels,
        'pt_labels': pt_labels
    }


def _detect_formed_bonds(p_mols, p_to_r, r_mols, Chem):
    formed_bonds, formed_labels = [], []
    for p_fi, p_mol in enumerate(p_mols):
        if not p_mol:
            continue
        _try_kekulize(p_mol, Chem)
        for bond in p_mol.GetBonds():
            r1 = p_to_r.get((p_fi, bond.GetBeginAtomIdx()))
            r2 = p_to_r.get((p_fi, bond.GetEndAtomIdx()))
            if not r1 or not r2 or r1[0] == r2[0]:
                continue
            formed_bonds.append((r1, r2))
            sym1 = r_mols[r1[0]].GetAtomWithIdx(r1[1]).GetSymbol() if r_mols[r1[0]] else '?'
            sym2 = r_mols[r2[0]].GetAtomWithIdx(r2[1]).GetSymbol() if r_mols[r2[0]] else '?'
            formed_labels.append(f"{sym1}(R{r1[0]+1}:{r1[1]})–{sym2}(R{r2[0]+1}:{r2[1]})")
    return formed_bonds, formed_labels


def _detect_broken_bonds(r_mols, r_to_p, Chem):
    broken_bonds, broken_labels = [], []
    for r_fi, r_mol in enumerate(r_mols):
        if not r_mol:
            continue
        _try_kekulize(r_mol, Chem)
        for bond in r_mol.GetBonds():
            b_idx = bond.GetBeginAtomIdx()
            e_idx = bond.GetEndAtomIdx()
            p1 = r_to_p.get((r_fi, b_idx))
            p2 = r_to_p.get((r_fi, e_idx))
            if p1 is None or p2 is None or p1[0] == p2[0]:
                continue
            broken_bonds.append(((r_fi, b_idx), (r_fi, e_idx)))
            sym1 = r_mol.GetAtomWithIdx(b_idx).GetSymbol()
            sym2 = r_mol.GetAtomWithIdx(e_idx).GetSymbol()
            broken_labels.append(f"{sym1}(R{r_fi+1}:{b_idx})–{sym2}(R{r_fi+1}:{e_idx})")
    return broken_bonds, broken_labels


def _count_atom_hs(atom):
    return atom.GetTotalNumHs() + sum(1 for n in atom.GetNeighbors() if n.GetAtomicNum() == 1)


def _detect_proton_transfers(r_mols, p_mols, r_to_p):
    donors, acceptors = [], []
    for (r_fi, r_ai), (p_fi, p_ai) in r_to_p.items():
        r_mol = r_mols[r_fi] if r_fi < len(r_mols) else None
        p_mol = p_mols[p_fi] if p_fi < len(p_mols) else None
        if not r_mol or not p_mol:
            continue
        r_atom = r_mol.GetAtomWithIdx(r_ai)
        p_atom = p_mol.GetAtomWithIdx(p_ai)
        r_h = _count_atom_hs(r_atom)
        p_h = _count_atom_hs(p_atom)
        if r_h > p_h:
            donors.append((r_fi, r_ai))
        elif r_h < p_h:
            acceptors.append((r_fi, r_ai))

    if len(donors) != 1 or len(acceptors) != 1 or donors[0][0] == acceptors[0][0]:
        return [], []
    d, a = donors[0], acceptors[0]
    sym_d = r_mols[d[0]].GetAtomWithIdx(d[1]).GetSymbol()
    sym_a = r_mols[a[0]].GetAtomWithIdx(a[1]).GetSymbol()
    return [(d, a)], [f"PT: {sym_d}(R{d[0]+1}:{d[1]})→{sym_a}(R{a[0]+1}:{a[1]})"]


def detect_driving_coordinates(reactants_smi, products_smi, mapping, ket_path=None):
    """
    Detect bonds that form and break during a reaction.

    If a KET file with push-arrows is provided, extracts driving coords exactly.
    Otherwise, uses atom mapping to compare connectivity in reactants vs products.
    """
    Chem, _, _, _ = _get_rdkit()

    p_to_r = {(pv[1], pv[2]): (rk[1], rk[2]) for rk, pv in mapping.items()}
    r_to_p = {(rk[1], rk[2]): (pv[1], pv[2]) for rk, pv in mapping.items()}

    params = Chem.SmilesParserParams()
    params.removeHs = False
    r_mols = [Chem.MolFromSmiles(s, params) for s in reactants_smi]
    p_mols = [Chem.MolFromSmiles(s, params) for s in products_smi]

    map_num_to_r = {}
    for r_fi, r_mol in enumerate(r_mols):
        if not r_mol:
            continue
        for a in r_mol.GetAtoms():
            m = a.GetAtomMapNum()
            if m > 0:
                map_num_to_r[m] = (r_fi, a.GetIdx())

    if ket_path and os.path.exists(ket_path):
        ket_mech = parse_ket_mechanism(ket_path)
        if ket_mech:
            print("  ⚡ Detected explicit push-arrows in KET file.")
            return _extract_ket_coords(ket_mech, map_num_to_r, r_mols)

    formed_bonds, formed_labels = _detect_formed_bonds(p_mols, p_to_r, r_mols, Chem)
    broken_bonds, broken_labels = _detect_broken_bonds(r_mols, r_to_p, Chem)
    pt_bonds, pt_labels = _detect_proton_transfers(r_mols, p_mols, r_to_p)

    return {
        'formed_bonds': formed_bonds,
        'broken_bonds': broken_bonds,
        'proton_transfers': pt_bonds,
        'formed_labels': formed_labels,
        'broken_labels': broken_labels,
        'pt_labels': pt_labels,
    }


def _map_pairs_to_complex(pairs, lookup):
    c_pairs, c_1based = [], []
    for r1, r2 in pairs:
        c1 = lookup.get((r1[0], r1[1]))
        c2 = lookup.get((r2[0], r2[1]))
        if c1 is not None and c2 is not None:
            c_pairs.append((c1, c2))
            c_1based.append((c1 + 1, c2 + 1))
    return c_pairs, c_1based


def map_driving_coords_to_complex(driving_coords, r_atom_origins):
    """Convert fragment-local driving coordinates to complex-global atom indices."""
    lookup = {(fi, ai): ci for ci, (fi, ai) in enumerate(r_atom_origins)}

    formed, formed_1based = _map_pairs_to_complex(driving_coords['formed_bonds'], lookup)
    broken, broken_1based = _map_pairs_to_complex(driving_coords['broken_bonds'], lookup)
    pt, pt_1based = _map_pairs_to_complex(driving_coords.get('proton_transfers', []), lookup)

    return {
        'formed': formed,
        'broken': broken,
        'pt': pt,
        'formed_1based': formed_1based,
        'broken_1based': broken_1based,
        'pt_1based': pt_1based,
    }