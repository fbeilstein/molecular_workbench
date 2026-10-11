"""
Reaction atom mapping utilities: SMILES parsing, MCS-based automatic mapping,
and proton transfer mapping correction.
"""

import json
import os


def _get_rdkit():
    from rdkit import Chem
    from rdkit.Chem import AllChem, rdFMCS, rdmolops
    return Chem, AllChem, rdFMCS, rdmolops


def parse_reaction_smiles(rxn_smiles):
    """Parse 'R1.R2>>P1.P2' into reactant and product SMILES lists."""
    if '>>' not in rxn_smiles:
        raise ValueError(f"Not a reaction SMILES (no '>>'): {rxn_smiles}")

    parts = rxn_smiles.split('>>')
    if len(parts) != 2:
        raise ValueError(f"Invalid reaction SMILES: {rxn_smiles}")

    reactant_str, product_str = parts
    reactants = [s.strip() for s in reactant_str.split('.') if s.strip()]
    products = [s.strip() for s in product_str.split('.') if s.strip()]

    if not reactants or not products:
        raise ValueError("Reaction must have both reactants and products")

    return reactants, products


def detect_charge(smiles):
    """Detect formal charge from SMILES using RDKit."""
    Chem, _, _, _ = _get_rdkit()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return 0
    return sum(a.GetFormalCharge() for a in mol.GetAtoms())


def _has_atom_maps(mols):
    for mol in mols:
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                return True
    return False


def _map_from_atom_numbers(r_mols, p_mols):
    r_map = {}
    for fi, mol in enumerate(r_mols):
        if not mol:
            continue
        for atom in mol.GetAtoms():
            mn = atom.GetAtomMapNum()
            if mn > 0:
                r_map[mn] = (fi, atom.GetIdx(), atom.GetSymbol())

    mapping = {}
    summary_lines = []
    for fi, mol in enumerate(p_mols):
        if not mol:
            continue
        for atom in mol.GetAtoms():
            mn = atom.GetAtomMapNum()
            if mn <= 0 or mn not in r_map:
                continue
            rfi, rai, rsym = r_map[mn]
            mapping[('R', rfi, rai)] = ('P', fi, atom.GetIdx())
            summary_lines.append(
                f"  R{rfi+1}:{rsym}({rai}) → P{fi+1}:{atom.GetSymbol()}({atom.GetIdx()})  [map #{mn}]"
            )
    return mapping, summary_lines


def _bond_change_cost(rb, pb):
    if rb is not None and pb is not None:
        return abs(rb.GetBondTypeAsDouble() - pb.GetBondTypeAsDouble())
    if (rb is None) != (pb is None):
        return 2.0
    return 0.0


def _score_alignment_pairs(pairs, rmol, pmol):
    cost = 0.0
    for ra, pa in pairs:
        cost += 2.0 * abs(rmol.GetAtomWithIdx(ra).GetTotalNumHs() - pmol.GetAtomWithIdx(pa).GetTotalNumHs())

    n_pairs = len(pairs)
    for x in range(n_pairs):
        for y in range(x + 1, n_pairs):
            rb = rmol.GetBondBetweenAtoms(pairs[x][0], pairs[y][0])
            pb = pmol.GetBondBetweenAtoms(pairs[x][1], pairs[y][1])
            cost += _bond_change_cost(rb, pb)
    return cost


def _find_best_mcs_alignment(rmol, pmol, r_match, p_matches, ri, pi, mapping, claimed_p):
    best_pairs, best_cost = None, None
    for p_match in p_matches:
        pairs = [(ra, pa) for ra, pa in zip(r_match, p_match) if ('R', ri, ra) not in mapping]
        if not pairs or any((pi, pa) in claimed_p for _, pa in pairs):
            continue
        cost = _score_alignment_pairs(pairs, rmol, pmol)
        if best_cost is None or cost < best_cost:
            best_pairs, best_cost = pairs, cost
    return best_pairs


def _align_mol_pair_mcs(rmol, pmol, ri, pi, mapping, claimed_p, rdFMCS, Chem):
    try:
        mcs = rdFMCS.FindMCS([rmol, pmol],
                             atomCompare=rdFMCS.AtomCompare.CompareElements,
                             bondCompare=rdFMCS.BondCompare.CompareAny,
                             timeout=10)
    except Exception as e:
        return [f"  Warning: MCS failed for R{ri+1}↔P{pi+1}: {e}"]

    if mcs.numAtoms == 0:
        return []

    patt = Chem.MolFromSmarts(mcs.smartsString)
    qp = Chem.AdjustQueryParameters.NoAdjustments()
    qp.makeBondsGeneric = True
    patt = Chem.AdjustQueryProperties(patt, qp)
    r_match = rmol.GetSubstructMatch(patt)
    p_matches = pmol.GetSubstructMatches(patt, uniquify=False, maxMatches=5000)

    best_pairs = _find_best_mcs_alignment(rmol, pmol, r_match, p_matches, ri, pi, mapping, claimed_p)
    if not r_match or not best_pairs:
        return []

    new_mappings = []
    for r_ai, p_ai in best_pairs:
        key = ('R', ri, r_ai)
        p_key = (pi, p_ai)
        if key in mapping or p_key in claimed_p:
            continue
        mapping[key] = ('P', pi, p_ai)
        claimed_p.add(p_key)
        rsym = rmol.GetAtomWithIdx(r_ai).GetSymbol()
        psym = pmol.GetAtomWithIdx(p_ai).GetSymbol()
        new_mappings.append(f"  R{ri+1}:{rsym}({r_ai}) → P{pi+1}:{psym}({p_ai})")
    return new_mappings


def auto_map_atoms(reactants_smi, products_smi):
    """
    Attempt automatic atom mapping using RDKit atom map numbers or MCS-based matching.
    """
    Chem, _, rdFMCS, _ = _get_rdkit()

    params = Chem.SmilesParserParams()
    params.removeHs = False
    r_mols = [Chem.MolFromSmiles(s, params) for s in reactants_smi]
    p_mols = [Chem.MolFromSmiles(s, params) for s in products_smi]

    if _has_atom_maps(r_mols + p_mols):
        mapping, summary_lines = _map_from_atom_numbers(r_mols, p_mols)
        return mapping, '\n'.join(summary_lines) if summary_lines else "  No mapping found"

    summary_lines = ["  (Auto-mapped via MCS — please verify)"]
    r_mols_sorted = sorted([(i, m) for i, m in enumerate(r_mols) if m is not None],
                           key=lambda x: x[1].GetNumHeavyAtoms(), reverse=True)
    p_mols_sorted = sorted([(i, m) for i, m in enumerate(p_mols) if m is not None],
                           key=lambda x: x[1].GetNumHeavyAtoms(), reverse=True)

    claimed_p = set()
    mapping = {}
    for ri, rmol in r_mols_sorted:
        for pi, pmol in p_mols_sorted:
            msgs = _align_mol_pair_mcs(rmol, pmol, ri, pi, mapping, claimed_p, rdFMCS, Chem)
            summary_lines.extend(msgs)

    summary = '\n'.join(summary_lines) if summary_lines else "  No mapping found"
    return mapping, summary


def _collect_pt_candidates(r_to_p, r_mols, p_mols):
    donors, acceptors = [], []
    for (r_fi, r_ai), (p_fi, p_ai) in r_to_p.items():
        r_mol = r_mols[r_fi] if r_fi < len(r_mols) else None
        p_mol = p_mols[p_fi] if p_fi < len(p_mols) else None
        if r_mol is None or p_mol is None:
            continue
        r_atom = r_mol.GetAtomWithIdx(r_ai)
        p_atom = p_mol.GetAtomWithIdx(p_ai)
        r_h, p_h = r_atom.GetTotalNumHs(), p_atom.GetTotalNumHs()
        if r_h > p_h:
            donors.append((r_fi, r_ai, p_fi, p_ai, r_h - p_h))
        elif r_h < p_h:
            acceptors.append((r_fi, r_ai, p_fi, p_ai, p_h - r_h))
    return donors, acceptors


def _is_pt_mapping_consistent(d_donor, a_acceptor, r_mols, p_mols):
    d_r_fi, d_r_ai, d_p_fi, d_p_ai, _ = d_donor
    a_r_fi, a_r_ai, a_p_fi, a_p_ai, _ = a_acceptor

    d_r_atom = r_mols[d_r_fi].GetAtomWithIdx(d_r_ai)
    d_p_atom = p_mols[d_p_fi].GetAtomWithIdx(d_p_ai)
    if d_p_atom.GetTotalNumHs() < d_r_atom.GetTotalNumHs():
        return True

    if d_p_atom.GetTotalNumHs() == d_r_atom.GetTotalNumHs():
        a_r_atom = r_mols[a_r_fi].GetAtomWithIdx(a_r_ai)
        a_p_atom = p_mols[a_p_fi].GetAtomWithIdx(a_p_ai)
        if a_p_atom.GetTotalNumHs() > a_r_atom.GetTotalNumHs():
            return True

    return False


def _validate_pt_mapping(mapping, reactants_smi, products_smi):
    """Fix atom mapping for proton transfer reactions."""
    Chem, _, _, _ = _get_rdkit()
    params = Chem.SmilesParserParams()
    params.removeHs = False
    r_mols = [Chem.MolFromSmiles(s, params) for s in reactants_smi]
    p_mols = [Chem.MolFromSmiles(s, params) for s in products_smi]

    r_to_p = {(rk[1], rk[2]): (pv[1], pv[2]) for rk, pv in mapping.items()}
    donors, acceptors = _collect_pt_candidates(r_to_p, r_mols, p_mols)

    if len(donors) != 1 or len(acceptors) != 1:
        return mapping

    d_donor, a_acceptor = donors[0], acceptors[0]
    if d_donor[0] == a_acceptor[0]:
        return mapping

    if _is_pt_mapping_consistent(d_donor, a_acceptor, r_mols, p_mols):
        return mapping

    # Swap product assignments
    d_rk = ('R', d_donor[0], d_donor[1])
    a_rk = ('R', a_acceptor[0], a_acceptor[1])
    d_pv, a_pv = mapping[d_rk], mapping[a_rk]

    new_mapping = dict(mapping)
    new_mapping[d_rk] = a_pv
    new_mapping[a_rk] = d_pv

    print(f"  ⚡ PT mapping fix: swapped R{d_donor[0]+1}:{d_donor[1]}↔R{a_acceptor[0]+1}:{a_acceptor[1]} product assignments")
    return new_mapping