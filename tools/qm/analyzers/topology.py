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

    def get_aromatic_atoms(self):
        """Return a list of atom indices that are aromatic."""
        if not self.mol:
            return []
        
        aromatic_indices = []
        for atom in self.mol.GetAtoms():
            if atom.GetIsAromatic():
                aromatic_indices.append(atom.GetIdx())
        
        return aromatic_indices

    def get_conjugated_atoms(self):
        """Return a list of atom indices that are part of conjugated pi systems."""
        if not self.mol:
            return []
            
        conjugated_indices = []
        for atom in self.mol.GetAtoms():
            # Aromatic atoms are conjugated
            if atom.GetIsAromatic():
                conjugated_indices.append(atom.GetIdx())
                continue
                
            # Non-aromatic atoms with pi bonds are conjugated if adjacent to another pi bond
            # In RDKit, an atom might be flagged as conjugated if it's involved in a conjugated system
            if hasattr(atom, 'GetIsConjugated') and atom.GetIsConjugated():
                conjugated_indices.append(atom.GetIdx())
                
        return conjugated_indices
