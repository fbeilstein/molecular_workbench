import numpy as np
from rdkit import Chem
from .common import rcov, _unit


def _parse_mapping(mapping_json, offsets):
    """Parse raw UI mapping json into reactant global index -> (pfi, pai)."""
    r2p = {}
    for k, v in (mapping_json or {}).items():
        try:
            _, rfi, rai = k.split("_")
            _, pfi, pai = v.split("_")
            ci = offsets[int(rfi)] + int(rai)
            r2p[ci] = (int(pfi), int(pai))
        except (ValueError, IndexError):
            continue

    seen, clean = set(), {}
    for ci in sorted(r2p):
        target = r2p[ci]
        if target in seen:
            continue
        seen.add(target)
        clean[ci] = target
    return clean


def _extract_product_bonds(p_mols, p2r):
    p_bonds = set()
    p_bond_orders = {}
    for pfi, m in enumerate(p_mols):
        for b in m.GetBonds():
            i = p2r.get((pfi, b.GetBeginAtomIdx()))
            j = p2r.get((pfi, b.GetEndAtomIdx()))
            if i is None or j is None:
                continue
            pair = (min(i, j), max(i, j))
            p_bonds.add(pair)
            p_bond_orders[pair] = float(b.GetBondTypeAsDouble())
    return p_bonds, p_bond_orders


def _match_free_atoms_by_symbol(symbols, free, r2p, n_atoms):
    for ci in range(n_atoms):
        if ci in r2p:
            continue
        for k, (pfi, pai, sym) in enumerate(free):
            if sym != symbols[ci]:
                continue
            r2p[ci] = (pfi, pai)
            free.pop(k)
            break


class ReactionModel:
    def __init__(self, raw_smiles, mapping_json, solvent: str = "water"):
        if not raw_smiles or ">>" not in raw_smiles:
            raise ValueError(f"Invalid reaction SMILES: '{raw_smiles}' (must contain '>>')")
        self.solvent = solvent
        r_part, p_part = raw_smiles.split(">>")
        self.r_smiles = [s.strip() for s in r_part.split(".") if s.strip()]
        self.p_smiles = [s.strip() for s in p_part.split(".") if s.strip()]

        self.r_mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in self.r_smiles]
        self.p_mols = [Chem.AddHs(Chem.MolFromSmiles(s)) for s in self.p_smiles]
        if any(m is None for m in self.r_mols + self.p_mols):
            raise ValueError("Could not parse reaction SMILES")

        for m in self.r_mols + self.p_mols:
            try:
                Chem.Kekulize(m, clearAromaticFlags=True)
            except Exception:
                pass

        self.offsets, off = [], 0
        for m in self.r_mols:
            self.offsets.append(off)
            off += m.GetNumAtoms()
        self.n_atoms = off

        self.frag_of, self.symbols = [], []
        for fi, m in enumerate(self.r_mols):
            for a in m.GetAtoms():
                self.frag_of.append(fi)
                self.symbols.append(a.GetSymbol())

        self._build_reactant_graph()

        clean_mapping = _parse_mapping(mapping_json, self.offsets)
        self.n_explicit = len(clean_mapping)
        self.r2p = self._complete_mapping(clean_mapping)
        self.p2r = {v: k for k, v in self.r2p.items()}

        self._build_product_graph()
        self._classify_bonds()

        self.charge = int(sum(a.GetFormalCharge() for m in self.r_mols for a in m.GetAtoms()))
        n_rad = int(sum(a.GetNumRadicalElectrons() for m in self.r_mols for a in m.GetAtoms()))
        self.multiplicity = n_rad + 1
        self.overridden = False

    def _build_reactant_graph(self):
        self.r_bonds = set()
        self.r_bond_orders = {}
        self.neighbors = [[] for _ in range(self.n_atoms)]
        for fi, m in enumerate(self.r_mols):
            o = self.offsets[fi]
            for b in m.GetBonds():
                i, j = o + b.GetBeginAtomIdx(), o + b.GetEndAtomIdx()
                pair = (min(i, j), max(i, j))
                self.r_bonds.add(pair)
                self.r_bond_orders[pair] = float(b.GetBondTypeAsDouble())
                self.neighbors[i].append(j)
                self.neighbors[j].append(i)
    def _build_product_graph(self):
        self.p_bonds, self.p_bond_orders = _extract_product_bonds(self.p_mols, self.p2r)

    def _classify_bonds(self):
        mapped = set(self.r2p)
        self.formed = sorted(b for b in self.p_bonds - self.r_bonds)
        self.broken = sorted(b for b in self.r_bonds - self.p_bonds if b[0] in mapped and b[1] in mapped)

        self.order_changed = []
        for pair in self.r_bonds.intersection(self.p_bonds):
            r_ord = self.r_bond_orders.get(pair, 1.0)
            p_ord = self.p_bond_orders.get(pair, 1.0)
            if abs(r_ord - p_ord) > 0.4:
                self.order_changed.append((pair[0], pair[1], r_ord, p_ord))

        changing = set(self.formed).union(set(self.broken)).union({(i, j) for i, j, _, _ in self.order_changed})
        self.spectator_bonds = sorted([b for b in self.r_bonds if b not in changing])

    def _find_neighbor_match(self, nb, ci, r2p, used):
        if nb not in r2p:
            return None
        pfi, pai = r2p[nb]
        patom = self.p_mols[pfi].GetAtomWithIdx(pai)
        for pn in patom.GetNeighbors():
            key = (pfi, pn.GetIdx())
            if pn.GetSymbol() == self.symbols[ci] and key not in used:
                return key
        return None

    def _find_mapped_neighbor(self, ci, r2p, used):
        for nb in self.neighbors[ci]:
            match = self._find_neighbor_match(nb, ci, r2p, used)
            if match:
                return match
        return None

    def _complete_mapping(self, r2p):
        r2p = dict(r2p)
        used = set(r2p.values())

        # Pass 1: Topological expansion from already-mapped neighbors
        for ci in range(self.n_atoms):
            if ci in r2p:
                continue
            match = self._find_mapped_neighbor(ci, r2p, used)
            if match:
                r2p[ci] = match
                used.add(match)

        # Pass 2: Remaining free atoms matched by symbol
        free = [(pfi, a.GetIdx(), a.GetSymbol())
                for pfi, m in enumerate(self.p_mols) for a in m.GetAtoms()
                if (pfi, a.GetIdx()) not in used]
        _match_free_atoms_by_symbol(self.symbols, free, r2p, self.n_atoms)
        return r2p

    def local(self, ci: int):
        fi = self.frag_of[ci]
        return fi, ci - self.offsets[fi]

    def label(self, ci: int) -> str:
        fi, li = self.local(ci)
        return f"{self.symbols[ci]}(R{fi + 1}:{li})"

    def d0(self, i: int, j: int) -> float:
        return rcov(self.symbols[i]) + rcov(self.symbols[j])

    def summary_lines(self):
        lines = [
            f"Charge {self.charge:+d}, multiplicity {self.multiplicity}",
            f"Mapped atoms: {len(self.r2p)}/{self.n_atoms} ({self.n_explicit} from Step 2, rest completed automatically)"
        ]
        for i, j in self.formed:
            lines.append(f"  FORMING  {self.label(i)} ... {self.label(j)}")
        for i, j in self.broken:
            lines.append(f"  BREAKING {self.label(i)} -- {self.label(j)}")
        for i, j, ro, po in self.order_changed:
            lines.append(f"  PI/ORDER {self.label(i)}--{self.label(j)}: order {ro:.1f} -> {po:.1f}")
        if not self.formed and not self.broken and not self.order_changed:
            lines.append("  (no bond changes detected - check the mapping)")
        return lines

    def reaction_info(self):
        return {
            "forming": [[list(self.local(i)), list(self.local(j))] for i, j in self.formed],
            "breaking": [[list(self.local(i)), list(self.local(j))] for i, j in self.broken],
            "order_changed": [[[list(self.local(i)), list(self.local(j))], ro, po]
                              for i, j, ro, po in self.order_changed],
            "summary": self.summary_lines(),
            "charge": self.charge,
            "multiplicity": self.multiplicity,
        }

    def _convert_pairs(self, pairs):
        out = []
        for (fi, ai), (fj, aj) in pairs:
            i = self.offsets[int(fi)] + int(ai)
            j = self.offsets[int(fj)] + int(aj)
            if i == j or not (0 <= i < self.n_atoms and 0 <= j < self.n_atoms):
                continue
            out.append((min(i, j), max(i, j)))
        return sorted(set(out))

    def apply_override(self, override):
        if not override:
            return False
        forming = self._convert_pairs(override.get("forming", []))
        breaking = self._convert_pairs(override.get("breaking", []))
        if not forming and not breaking:
            return False

        self.formed, self.broken = forming, breaking
        changing = set(self.formed).union(set(self.broken)).union({(i, j) for i, j, _, _ in self.order_changed})
        self.spectator_bonds = sorted([b for b in self.r_bonds if b not in changing])
        self.overridden = True
        return True

    def attack_mode(self, ci: int, coords_of):
        x = coords_of(ci)
        nbrs = self.neighbors[ci]
        leaving = [j for j in nbrs if (min(ci, j), max(ci, j)) in self.broken]
        if leaving:
            v = -sum(_unit(coords_of(j) - x) for j in leaving)
            if np.linalg.norm(v) > 1e-3:
                return ("vec", _unit(v))

        if not nbrs:
            return None

        units = [_unit(coords_of(j) - x) for j in nbrs]
        s = np.sum(units, axis=0)
        if len(nbrs) == 1 or np.linalg.norm(s) > 0.35:
            return ("vec", _unit(-s))
        if len(nbrs) == 2:
            return ("perp", units[0])
        if len(nbrs) == 3:
            n = np.cross(units[1] - units[0], units[2] - units[0])
            if np.linalg.norm(n) > 1e-6:
                return ("axis", _unit(n))
        return None

    def embed_fragments(self, seed=42):
        from .docking import embed_fragments
        return embed_fragments(self, seed=seed)

    def dock(self, frag_coords, variant=0, n_starts=32, seed=7):
        from .docking import dock_fragments
        return dock_fragments(self, frag_coords, variant=variant, n_starts=n_starts, seed=seed)
