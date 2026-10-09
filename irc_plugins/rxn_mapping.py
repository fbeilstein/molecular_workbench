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


def auto_map_atoms(reactants_smi, products_smi):
    """
    Attempt automatic atom mapping using RDKit atom map numbers
    or MCS-based matching.

    Returns a mapping dict and human-readable summary.
    """
    Chem, _, rdFMCS, _ = _get_rdkit()

    params = Chem.SmilesParserParams()
    params.removeHs = False
    r_mols = [Chem.MolFromSmiles(s, params) for s in reactants_smi]
    p_mols = [Chem.MolFromSmiles(s, params) for s in products_smi]

    # Check if atom map numbers are already present
    has_maps = False
    for mol in r_mols + p_mols:
        if mol is None:
            continue
        for atom in mol.GetAtoms():
            if atom.GetAtomMapNum() > 0:
                has_maps = True
                break

    mapping = {}
    summary_lines = []

    if has_maps:
        r_map = {}
        for fi, mol in enumerate(r_mols):
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                mn = atom.GetAtomMapNum()
                if mn > 0:
                    r_map[mn] = (fi, atom.GetIdx(), atom.GetSymbol())

        for fi, mol in enumerate(p_mols):
            if mol is None:
                continue
            for atom in mol.GetAtoms():
                mn = atom.GetAtomMapNum()
                if mn > 0 and mn in r_map:
                    rfi, rai, rsym = r_map[mn]
                    mapping[('R', rfi, rai)] = ('P', fi, atom.GetIdx())
                    summary_lines.append(
                        f"  R{rfi+1}:{rsym}({rai}) → P{fi+1}:{atom.GetSymbol()}({atom.GetIdx()})  [map #{mn}]"
                    )
    else:
        summary_lines.append("  (Auto-mapped via MCS — please verify)")
        # Sort by number of heavy atoms (descending) to match largest backbones first
        r_mols_sorted = sorted([(i, m) for i, m in enumerate(r_mols) if m is not None],
                               key=lambda x: x[1].GetNumHeavyAtoms(), reverse=True)
        p_mols_sorted = sorted([(i, m) for i, m in enumerate(p_mols) if m is not None],
                               key=lambda x: x[1].GetNumHeavyAtoms(), reverse=True)

        claimed_p = set()
        for ri, rmol in r_mols_sorted:
            for pi, pmol in p_mols_sorted:
                try:
                    mcs = rdFMCS.FindMCS([rmol, pmol],
                                         atomCompare=rdFMCS.AtomCompare.CompareElements,
                                         bondCompare=rdFMCS.BondCompare.CompareAny,
                                         timeout=10)
                    if mcs.numAtoms == 0:
                        continue

                    patt = Chem.MolFromSmarts(mcs.smartsString)
                    r_match = rmol.GetSubstructMatch(patt)
                    p_match = pmol.GetSubstructMatch(patt)

                    if r_match and p_match:
                        for r_ai, p_ai in zip(r_match, p_match):
                            key = ('R', ri, r_ai)
                            p_key = (pi, p_ai)
                            if key not in mapping and p_key not in claimed_p:
                                mapping[key] = ('P', pi, p_ai)
                                claimed_p.add(p_key)
                                rsym = rmol.GetAtomWithIdx(r_ai).GetSymbol()
                                psym = pmol.GetAtomWithIdx(p_ai).GetSymbol()
                                summary_lines.append(
                                    f"  R{ri+1}:{rsym}({r_ai}) → P{pi+1}:{psym}({p_ai})"
                                )
                except Exception as e:
                    summary_lines.append(f"  Warning: MCS failed for R{ri+1}↔P{pi+1}: {e}")

    summary = '\n'.join(summary_lines) if summary_lines else "  No mapping found"
    return mapping, summary


def _validate_pt_mapping(mapping, reactants_smi, products_smi):
    """Fix atom mapping for proton transfer reactions.

    MCS can swap donor/acceptor when both are the same element (e.g. two
    oxygens).  Detect this by checking hydrogen-count consistency: the
    reactant atom that *loses* an H must map to the product atom that
    has *fewer* H, not more.

    Returns corrected mapping dict.
    """
    Chem, _, _, _ = _get_rdkit()
    params = Chem.SmilesParserParams()
    params.removeHs = False
    r_mols = [Chem.MolFromSmiles(s, params) for s in reactants_smi]
    p_mols = [Chem.MolFromSmiles(s, params) for s in products_smi]

    r_to_p = {(rk[1], rk[2]): (pv[1], pv[2]) for rk, pv in mapping.items()}

    # Collect atoms whose H-count changes
    donors = []   # (r_fi, r_ai, p_fi, p_ai, delta_h)
    acceptors = []
    for (r_fi, r_ai), (p_fi, p_ai) in r_to_p.items():
        r_mol = r_mols[r_fi]
        p_mol = p_mols[p_fi]
        if r_mol is None or p_mol is None:
            continue
        r_atom = r_mol.GetAtomWithIdx(r_ai)
        p_atom = p_mol.GetAtomWithIdx(p_ai)
        r_h = r_atom.GetTotalNumHs()
        p_h = p_atom.GetTotalNumHs()
        if r_h > p_h:
            donors.append((r_fi, r_ai, p_fi, p_ai, r_h - p_h))
        elif r_h < p_h:
            acceptors.append((r_fi, r_ai, p_fi, p_ai, p_h - r_h))

    if len(donors) != 1 or len(acceptors) != 1:
        return mapping  # Not a simple 1-to-1 proton transfer

    d_r_fi, d_r_ai, d_p_fi, d_p_ai, _ = donors[0]
    a_r_fi, a_r_ai, a_p_fi, a_p_ai, _ = acceptors[0]

    # The donor (loses H) and acceptor (gains H) must be in DIFFERENT
    # reactant fragments for this to be an inter-fragment proton transfer
    if d_r_fi == a_r_fi:
        return mapping  # Intra-fragment rearrangement, mapping is fine

    # Check: is the mapping consistent?
    # The donor atom loses H, so its product partner should also have
    # fewer H than before.  If the product partner of the donor actually
    # gained H, the mapping is swapped.
    d_r_atom = r_mols[d_r_fi].GetAtomWithIdx(d_r_ai)
    d_p_atom = p_mols[d_p_fi].GetAtomWithIdx(d_p_ai)
    if d_p_atom.GetTotalNumHs() > d_r_atom.GetTotalNumHs():
        # Already correct (product has more H — wait, donor should lose H!)
        # If product partner has MORE H, that means the donor was mapped
        # to the product that gained H — WRONG!
        pass  # Fall through to swap
    elif d_p_atom.GetTotalNumHs() < d_r_atom.GetTotalNumHs():
        # Product partner has fewer H than reactant — consistent!
        return mapping
    else:
        # Same H count — check acceptor side
        a_r_atom = r_mols[a_r_fi].GetAtomWithIdx(a_r_ai)
        a_p_atom = p_mols[a_p_fi].GetAtomWithIdx(a_p_ai)
        if a_p_atom.GetTotalNumHs() > a_r_atom.GetTotalNumHs():
            return mapping  # Acceptor gained H — consistent
        # else fall through to swap

    # Mapping is swapped — the donor was matched to the acceptor's product.
    # Swap the product assignments of the two atoms.
    d_rk = ('R', d_r_fi, d_r_ai)
    a_rk = ('R', a_r_fi, a_r_ai)
    d_pv = mapping[d_rk]  # currently wrong
    a_pv = mapping[a_rk]  # currently wrong

    new_mapping = dict(mapping)
    new_mapping[d_rk] = a_pv  # donor → acceptor's old product
    new_mapping[a_rk] = d_pv  # acceptor → donor's old product

    print(f"  ⚡ PT mapping fix: swapped "
          f"R{d_r_fi+1}:{d_r_ai}↔R{a_r_fi+1}:{a_r_ai} product assignments")
    return new_mapping