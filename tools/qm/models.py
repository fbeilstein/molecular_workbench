from dataclasses import dataclass, field
from typing import List, Optional
import numpy as np

@dataclass
class Orbital:
    """Represents a classified molecular orbital before 3D rendering."""
    type: str  # e.g., 'sigma', 'pi', 'sigma_star', 'pi_star', 'lone_pair', 'core'
    atoms: List[str]  # e.g., ['C1', 'C2'] or ['C1'] or ['ring']
    mo_coeff: np.ndarray  # The 1D array of orbital coefficients
    energy_ev: Optional[float] = None
    canonical_label: Optional[str] = None
    sym_score: Optional[float] = None  # For anti-bonds to track symmetry

    def to_dict(self, cube_file_path: str) -> dict:
        """Export to the JSON manifest format."""
        d = {
            'atoms': self.atoms,
            'file': cube_file_path
        }
        if self.energy_ev is not None:
            d['energy_ev'] = round(self.energy_ev, 2)
        if self.canonical_label is not None:
            d['canonical_label'] = self.canonical_label
        return d
