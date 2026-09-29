import os
from rdkit import Chem

class TopologyDetector:
    def __init__(self, smiles: str):
        self.smiles = smiles
        self.mol = Chem.MolFromSmiles(smiles)
        if self.mol:
            # We need to map the SMILES to the exact same hydrogen-added graph used in the XYZ.
            # mol_prep uses AddHs, so we do the same.
            self.mol = Chem.AddHs(self.mol)
            Chem.Kekulize(self.mol, clearAromaticFlags=False)
            
        # Also keep a raw unsanitized version. This preserves explicit double bonds
        # (e.g. in hypervalent atoms like [O-][Cl]=O) that RDKit's sanitizer might downgrade.
        self.mol_raw = Chem.MolFromSmiles(smiles, sanitize=False)

    def get_aromatic_atoms(self):
        """Return a list of atom indices that are aromatic."""
        if not self.mol:
            return []
        
        aromatic_indices = []
        for atom in self.mol.GetAtoms():
            if atom.GetIsAromatic():
                aromatic_indices.append(atom.GetIdx())
        
        return aromatic_indices

    def get_conjugated_fragments(self, min_size=3):
        """Return a list of conjugated fragments.
        Each fragment is a list of atom indices that form a continuous conjugated system.
        Only fragments with at least `min_size` atoms are returned.
        """
        if not self.mol or not self.mol_raw:
            return []
            
        adj = {}
        
        # 1. Aromatic and standard conjugated bonds (from sanitized mol)
        for bond in self.mol.GetBonds():
            if bond.GetIsConjugated() or bond.GetIsAromatic():
                u = bond.GetBeginAtomIdx()
                v = bond.GetEndAtomIdx()
                adj.setdefault(u, set()).add(v)
                adj.setdefault(v, set()).add(u)
                
        # 2. Hypervalent Resonance Systems
        # RDKit misses resonance in hypervalent atoms (e.g. [O-][Cl]=O).
        # We manually group them ONLY if there is a true resonance structure 
        # (a double bond mixing with a charged single bond).
        # This prevents falsely grouping isolated double bonds on tetrahedral centers (e.g. HClO3).
        for atom in self.mol_raw.GetAtoms():
            if atom.GetSymbol() in ['S', 'P', 'Cl', 'Br', 'I', 'N']:
                double_neighbors = []
                charged_single_neighbors = []
                for bond in atom.GetBonds():
                    neighbor = bond.GetOtherAtom(atom)
                    if bond.GetBondType() == Chem.BondType.DOUBLE:
                        double_neighbors.append(neighbor.GetIdx())
                    elif bond.GetBondType() == Chem.BondType.SINGLE and neighbor.GetFormalCharge() != 0:
                        charged_single_neighbors.append(neighbor.GetIdx())
                
                # If resonance is possible between the double bond(s) and charged single bond(s), group them
                if double_neighbors and charged_single_neighbors:
                    u = atom.GetIdx()
                    for v in double_neighbors + charged_single_neighbors:
                        adj.setdefault(u, set()).add(v)
                        adj.setdefault(v, set()).add(u)
                
        visited = set()
        fragments = []
        
        for node in adj:
            if node not in visited:
                comp = []
                queue = [node]
                visited.add(node)
                while queue:
                    curr = queue.pop(0)
                    comp.append(curr)
                    for neighbor in adj.get(curr, []):
                        if neighbor not in visited:
                            visited.add(neighbor)
                            queue.append(neighbor)
                if len(comp) >= min_size:
                    fragments.append(sorted(comp))
                    
        return fragments
