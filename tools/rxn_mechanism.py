import json
import os
import math
from rxn_mapping import _get_rdkit


def parse_ket_mechanism(ket_path):
    import json
    import math
    with open(ket_path) as f:
        ket = json.load(f)

    arrows = []
    for node in ket.get('root', {}).get('nodes', []):
        if node.get('type') == 'arrow' and 'elliptical-arc-arrow' in node.get('data', {}).get('mode', ''):
            arrows.append(node)

    if not arrows:
        return None

    mols = {}
    for k, v in ket.items():
        if k.startswith('mol'):
            mols[k] = v

    def dist(p1, p2):
        return math.sqrt((p1['x']-p2[0])**2 + (p1['y']-p2[1])**2)

    formed_bonds = []
    broken_bonds = []

    for a in arrows:
        start = a['data']['pos'][0]
        end = a['data']['pos'][1]
        
        start_map = None
        start_is_bond = False
        min_ds = 999
        
        end_map = None
        min_de = 999
        
        for mk, mol in mols.items():
            if mol.get('type') != 'molecule': continue
            
            # Check atoms
            for atom in mol.get('atoms', []):
                loc = atom['location']
                ds = dist(start, loc)
                de = dist(end, loc)
                if ds < min_ds:
                    min_ds = ds
                    start_map = atom.get("mapping")
                    start_is_bond = False
                if de < min_de:
                    min_de = de
                    end_map = atom.get("mapping")
                    
            # Check bonds
            for bond in mol.get('bonds', []):
                a1 = mol['atoms'][bond['atoms'][0]]['location']
                a2 = mol['atoms'][bond['atoms'][1]]['location']
                mid = [(a1[0]+a2[0])/2, (a1[1]+a2[1])/2]
                ds = dist(start, mid)
                if ds < min_ds:
                    min_ds = ds
                    map1 = mol['atoms'][bond['atoms'][0]].get('mapping')
                    map2 = mol['atoms'][bond['atoms'][1]].get('mapping')
                    start_map = (map1, map2)
                    start_is_bond = True

        if start_is_bond and isinstance(start_map, tuple):
            if start_map[0] is not None and start_map[1] is not None:
                broken_bonds.append(sorted([start_map[0], start_map[1]]))
        elif not start_is_bond and start_map is not None and end_map is not None:
            formed_bonds.append(sorted([start_map, end_map]))
            
    return {
        'formed_bonds': formed_bonds,
        'broken_bonds': broken_bonds
    }


def detect_driving_coordinates(reactants_smi, products_smi, mapping, ket_path=None):
    """
    Detect bonds that form and break during a reaction.

    If a KET file with push-arrows is provided, extracts driving coords exactly.
    Otherwise, uses atom mapping to compare connectivity in reactants vs products.

    Args:
        reactants_smi: list of reactant SMILES
        products_smi: list of product SMILES
        mapping: dict {('R', frag_idx, atom_idx): ('P', frag_idx, atom_idx)}
        ket_path: optional path to reaction.ket file

    Returns:
        dict with keys:
          - formed_bonds: [(r_atom1, r_atom2), ...] where r_atom = ('R', frag_idx, atom_idx)
          - broken_bonds: [(r_atom1, r_atom2), ...]
          - formed_labels: ["Cl(0)–C(1)", ...] human-readable
          - broken_labels: ["C(1)–Br(2)", ...]
    """
    Chem, _, _, _ = _get_rdkit()

    # Build reverse mapping: product → reactant coordinates
    p_to_r = {(pv[1], pv[2]): (rk[1], rk[2])
              for rk, pv in mapping.items()}
    r_to_p = {(rk[1], rk[2]): (pv[1], pv[2])
              for rk, pv in mapping.items()}

    params = Chem.SmilesParserParams()
    params.removeHs = False
    r_mols = [Chem.MolFromSmiles(s, params) for s in reactants_smi]
    p_mols = [Chem.MolFromSmiles(s, params) for s in products_smi]

    # Map number -> (frag_idx, atom_idx)
    map_num_to_r_atom = {}
    for r_fi, r_mol in enumerate(r_mols):
        if not r_mol: continue
        for a in r_mol.GetAtoms():
            m = a.GetAtomMapNum()
            if m > 0:
                map_num_to_r_atom[m] = (r_fi, a.GetIdx())

    # Check for KET explicit mechanism
    if ket_path and os.path.exists(ket_path):
        ket_mech = parse_ket_mechanism(ket_path)
        if ket_mech:
            print("  ⚡ Detected explicit push-arrows in KET file.")
            formed_bonds = []
            formed_labels = []
            for b1, b2 in ket_mech['formed_bonds']:
                r1 = map_num_to_r_atom.get(b1)
                r2 = map_num_to_r_atom.get(b2)
                if r1 and r2:
                    formed_bonds.append((r1, r2))
                    sym1 = r_mols[r1[0]].GetAtomWithIdx(r1[1]).GetSymbol()
                    sym2 = r_mols[r2[0]].GetAtomWithIdx(r2[1]).GetSymbol()
                    formed_labels.append(f"{sym1}(R{r1[0]+1}:{r1[1]})–{sym2}(R{r2[0]+1}:{r2[1]})")
            
            broken_bonds = []
            broken_labels = []
            for b1, b2 in ket_mech['broken_bonds']:
                r1 = map_num_to_r_atom.get(b1)
                r2 = map_num_to_r_atom.get(b2)
                if r1 and r2:
                    broken_bonds.append((r1, r2))
                    sym1 = r_mols[r1[0]].GetAtomWithIdx(r1[1]).GetSymbol()
                    sym2 = r_mols[r2[0]].GetAtomWithIdx(r2[1]).GetSymbol()
                    broken_labels.append(f"{sym1}(R{r1[0]+1}:{r1[1]})–{sym2}(R{r2[0]+1}:{r2[1]})")

            # Check if it's a proton transfer (broken and formed share an atom)
            pt_labels = []
            proton_transfers = []
            for fb in formed_bonds:
                for bb in broken_bonds:
                    shared = set(fb).intersection(bb)
                    if shared:
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

    # Heuristic fallback
    formed_bonds = []
    formed_labels = []
    for p_fi, p_mol in enumerate(p_mols):
        if not p_mol:
            continue
        try:
            Chem.Kekulize(p_mol, clearAromaticFlags=True)
        except Exception:
            pass
        for bond in p_mol.GetBonds():
            r1 = p_to_r.get((p_fi, bond.GetBeginAtomIdx()))
            r2 = p_to_r.get((p_fi, bond.GetEndAtomIdx()))
            if r1 and r2 and r1[0] != r2[0]:
                formed_bonds.append((r1, r2))
                # Human-readable label
                sym1 = r_mols[r1[0]].GetAtomWithIdx(r1[1]).GetSymbol() if r_mols[r1[0]] else '?'
                sym2 = r_mols[r2[0]].GetAtomWithIdx(r2[1]).GetSymbol() if r_mols[r2[0]] else '?'
                formed_labels.append(f"{sym1}(R{r1[0]+1}:{r1[1]})–{sym2}(R{r2[0]+1}:{r2[1]})")

    # Detect broken bonds: bonded in reactant, cross-fragment in product
    broken_bonds = []
    broken_labels = []
    for r_fi, r_mol in enumerate(r_mols):
        if not r_mol:
            continue
        try:
            Chem.Kekulize(r_mol, clearAromaticFlags=True)
        except Exception:
            pass
        for bond in r_mol.GetBonds():
            p1 = r_to_p.get((r_fi, bond.GetBeginAtomIdx()))
            p2 = r_to_p.get((r_fi, bond.GetEndAtomIdx()))
            if (p1 is not None) and (p2 is not None) and (p1[0] != p2[0]):
                broken_bonds.append((
                    (r_fi, bond.GetBeginAtomIdx()),
                    (r_fi, bond.GetEndAtomIdx())
                ))
                sym1 = r_mol.GetAtomWithIdx(bond.GetBeginAtomIdx()).GetSymbol()
                sym2 = r_mol.GetAtomWithIdx(bond.GetEndAtomIdx()).GetSymbol()
                broken_labels.append(f"{sym1}(R{r_fi+1}:{bond.GetBeginAtomIdx()})–{sym2}(R{r_fi+1}:{bond.GetEndAtomIdx()})")

    # Detect proton transfers (implicit hydrogen changes on mapped heavy atoms)
    proton_donors = []
    proton_acceptors = []
    
    for r_fi, r_mol in enumerate(r_mols):
        if not r_mol: continue
        for r_atom in r_mol.GetAtoms():
            r_ai = r_atom.GetIdx()
            p_match = r_to_p.get((r_fi, r_ai))
            if p_match:
                p_fi, p_ai = p_match
                p_mol = p_mols[p_fi]
                if p_mol:
                    p_atom = p_mol.GetAtomWithIdx(p_ai)
                    # When removeHs=False, explicit hydrogens are full atoms in the graph, so we must count them
                    r_h_implicit = r_atom.GetTotalNumHs()
                    r_h_explicit = sum(1 for n in r_atom.GetNeighbors() if n.GetAtomicNum() == 1)
                    r_h = r_h_implicit + r_h_explicit
                    
                    p_h_implicit = p_atom.GetTotalNumHs()
                    p_h_explicit = sum(1 for n in p_atom.GetNeighbors() if n.GetAtomicNum() == 1)
                    p_h = p_h_implicit + p_h_explicit
                    
                    if r_h > p_h:
                        proton_donors.append((r_fi, r_ai))
                    elif r_h < p_h:
                        proton_acceptors.append((r_fi, r_ai))

    proton_transfers = []
    pt_labels = []
    if len(proton_donors) == 1 and len(proton_acceptors) == 1:
        d_fi, d_ai = proton_donors[0]
        a_fi, a_ai = proton_acceptors[0]
        if d_fi != a_fi:
            proton_transfers.append((proton_donors[0], proton_acceptors[0]))
            sym_d = r_mols[d_fi].GetAtomWithIdx(d_ai).GetSymbol()
            sym_a = r_mols[a_fi].GetAtomWithIdx(a_ai).GetSymbol()
            pt_labels.append(f"PT: {sym_d}(R{d_fi+1}:{d_ai})→{sym_a}(R{a_fi+1}:{a_ai})")

    return {
        'formed_bonds': formed_bonds,
        'broken_bonds': broken_bonds,
        'proton_transfers': proton_transfers,
        'formed_labels': formed_labels,
        'broken_labels': broken_labels,
        'pt_labels': pt_labels,
    }


def map_driving_coords_to_complex(driving_coords, r_atom_origins):
    """
    Convert fragment-local driving coordinates to complex-global atom indices.

    Args:
        driving_coords: dict from detect_driving_coordinates()
        r_atom_origins: list of (frag_idx, atom_idx) for each complex atom

    Returns:
        dict with:
          - formed: [(complex_idx_1, complex_idx_2), ...] 0-based
          - broken: [(complex_idx_1, complex_idx_2), ...] 0-based
          - formed_1based: [(i, j), ...] for xTB input
          - broken_1based: [(i, j), ...] for xTB input
    """
    def frag_to_complex(frag_idx, atom_idx):
        for ci, (fi, ai) in enumerate(r_atom_origins):
            if fi == frag_idx and ai == atom_idx:
                return ci
        return None

    formed = []
    formed_1based = []
    for (r1, r2) in driving_coords['formed_bonds']:
        ci1 = frag_to_complex(r1[0], r1[1])
        ci2 = frag_to_complex(r2[0], r2[1])
        if ci1 is not None and ci2 is not None:
            formed.append((ci1, ci2))
            formed_1based.append((ci1 + 1, ci2 + 1))

    broken = []
    broken_1based = []
    for (r1, r2) in driving_coords['broken_bonds']:
        ci1 = frag_to_complex(r1[0], r1[1])
        ci2 = frag_to_complex(r2[0], r2[1])
        if ci1 is not None and ci2 is not None:
            broken.append((ci1, ci2))
            broken_1based.append((ci1 + 1, ci2 + 1))

    pt = []
    pt_1based = []
    for (r1, r2) in driving_coords.get('proton_transfers', []):
        ci1 = frag_to_complex(r1[0], r1[1])
        ci2 = frag_to_complex(r2[0], r2[1])
        if ci1 is not None and ci2 is not None:
            pt.append((ci1, ci2))
            pt_1based.append((ci1 + 1, ci2 + 1))

    return {
        'formed': formed,
        'broken': broken,
        'pt': pt,
        'formed_1based': formed_1based,
        'broken_1based': broken_1based,
        'pt_1based': pt_1based,
    }