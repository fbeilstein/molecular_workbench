# Agent Knowledge & Technical Gotchas
**WARNING: Read this file immediately if you are an AI agent returning to this project after losing context.**

## 1. PySCF Symmetry & xTB Noise
- **DO NOT ENABLE SYMMETRY.** The `gto.M(..., symmetry=True)` argument must remain `False`. 
- **Why?** The pipeline optimizes geometries with xTB first. xTB produces atomic coordinates that often break strict point-group symmetry by negligible amounts (e.g., $10^{-6}$ Å). PySCF's symmetric SCF builder enforces rigid symmetry groups; these microscopic numerical deviations will cause PySCF to crash with matrix dimension mismatch errors when computing overlap matrices (e.g., `ValueError: operands could not be broadcast together with shapes (10,10) (11,11)`). Keep `symmetry=False` in `scf.py` and `pipeline.py`.

## 2. Core Electron Handling (`chemcore`)
- **NEVER localize the entire occupied space.** 
- **Why?** Core electrons (e.g., the $1s$ electrons of Carbon and Oxygen) are tightly bound and $>99\%$ localized on a single atom. If they are passed into the `ValenceLocalizer`, the `classification.py` logic will falsely flag them as "Valence Lone Pairs" because they satisfy the population checks (`>95% on one heavy atom`). 
- **The Fix:** The pipeline uses `from pyscf.data.elements import chemcore` to explicitly calculate the number of core orbitals `n_core = chemcore(mol)`. The occupied array passed to the localizer must be strictly sliced to exclude them: `occ_coeff = mf.mo_coeff[:, n_core:n_occ]`.

## 3. UI Label Redundancy (`app.js` and `bundle_exporter.py`)
- **DO NOT append `(Energy eV)` to labels in the Python backend.**
- **Why?** The Javascript frontend (`static/app.js`) is hardcoded to automatically append the energy string `(${orb.energy_ev} eV)` to the label of any UI item that possesses an `energy_ev` property in its JSON block. 
- **Gotcha:** If you write backend code that sets a manifest label to `f"π (HOMO-4), {energy} eV"`, the UI will render it as `π (HOMO-4), -9.88 eV (-9.88 eV)`. Keep the python-generated labels clean and let the frontend handle the energy append.

## 4. Manifest Keys and Renderer Synchronization
- **The Dict Chain:** The pipeline creates `Orbital` objects with a specific `.type` (e.g., `pi_star`). `renderer.py` checks this `.type` and assigns it to a specific array in the manifest dictionary (e.g., `manifest['localized']['pi_star']`). `bundle_exporter.py` then reads that exact key to format the UI grouping.
- **Gotcha:** If you invent a new `.type` string in `valence.py` (e.g., `rydberg`), you MUST explicitly add an `elif orb.type == 'rydberg':` branch in `renderer.py` and a corresponding rendering block in `bundle_exporter.py`. If you don't map it through all three files, the orbital will be silently dropped and will never appear in the web UI.

## 5. Lone Pair vs Bond Thresholds
- **Gotcha:** The classification logic in `classification.py` requires a secondary atom to possess $>10\%$ of the orbital's electron density to qualify as a bond partner. Do not lower this threshold to 4-5%. 
- **Why?** Polar lone pairs (like the in-plane $p_y$ lone pair on Oxygen in $C=O$) possess a tiny $\sim4.5\%$ delocalization tail onto the adjacent Carbon atom. If the bond threshold is too low, the algorithm will detect 2 atoms, check the angular momentum, see it's perpendicular to the axis, and falsely categorize the Oxygen lone pair as a second $\pi$ bond (resulting in two $\pi$ bonds and two $\pi^*$ antibonds for Formaldehyde!). The 10% threshold forces it to safely fall back to the `Lone Pair` logic.

## 6. Stereochemistry in SMILES (The `\C` JSON Bug)
- **Gotcha:** Never manually construct JSON string values that contain SMILES strings using f-strings (e.g., `f' "smiles": "{smiles}" '`).
- **Why?** SMILES strings containing stereochemistry use backslashes (e.g., `C(/C)(\C)=C(\C)/C` for tetramethylethylene). If injected directly into a JSON template, the backslash creates an invalid escape sequence (like `\C`) causing a `JSONDecodeError`. 
- **The Fix:** Always use `json.dumps(smiles)` to properly escape backslashes, and if writing a script via bash heredoc, always quote the EOF delimiter (`cat <<'EOF'`) to prevent bash from prematurely expanding the backslashes.

## 7. The Hypervalent "Triple Bond" Demotion
- **Gotcha:** Do not remove the `seen_pi` interceptor loop in `valence.py`!
- **Why?** In hypervalent molecules (like $ClO_4^-$ or $SO_4^{2-}$), terminal oxygens have two degenerate perpendicular lone pairs. Both of these donate heavily ($>10\%$) to the central atom's $d$-orbitals. Without the interceptor, the localizer flags *both* as $\pi$ bonds, resulting in a physically accurate but visually ridiculous "triple bond" for every terminal oxygen. The `seen_pi` logic elegantly forces the second degenerate pair back down to a localized oxygen lone pair, enforcing textbook Lewis structure aesthetics.

## 8. Frontier Orbitals vs. Chemical Groups (`pipeline.py`)
- **DO NOT** use `continue` statements to skip `homo_idx` or `lumo_idx` when iterating over canonical pi systems (`pi_occ`/`pi_virt`) in `pipeline.py` or when localizing orbitals.
- **Why?** The UI groups orbitals chemically (`π Bonds`, `Lone Pairs`, etc.) and also explicitly highlights the `Frontier` (HOMO/LUMO). A previous developer attempted to avoid "duplicating" the HOMO in the UI by explicitly excluding it from the `π Bonds` list (via `if idx == homo_idx: continue`). 
- **Gotcha:** This caused the HOMO to be entirely dropped from the underlying chemical group list. Consequently, molecules where the HOMO is a $\pi$ orbital (e.g. cyclopentadienyl anion, annulene) incorrectly rendered one fewer $\pi$ bond in the frontend checkboxes. The correct behavior is to let the HOMO/LUMO be processed normally into their chemical systems; they will intentionally be duplicated in the JSON manifest (once as `_homo.json` for the Frontier UI section, and once as `_pi_...json` for the chemical group UI section) so both UI checklists are mathematically complete.

## 9. Strict File Length Limit (< 400 Lines)
- **Constraint:** All non-vendor source code files (Python, JS, CSS) must strictly remain **under 400 lines of code** (target < 350 lines).
- **Why?** Modular, readable architecture prevents monolithic sprawl.
- **Gotcha:** When editing or extending features (e.g., in UI managers, packagers, or plugins), always check line counts (`wc -l <file>`) before and after edits. Decompose large handlers into dedicated helper modules.

## 10. Bundle Structure & Direct URL Resolution
- **Bundle Directories:**
  - `bundles/reactions/`: Multistep reaction pathway archives (`.zip`) containing `multistep_pathway.json`, `mechanism.svg`, `svg_wiring.json`, and precomputed surface meshes.
  - `bundles/molecules/`: Single-molecule orbital archives (`.zip`) containing molecule definitions and surface meshes.
- **Explicit URL Fetching:**
  - `viewer_loader.js` performs direct, transparent `fetch(url)` requests without proxy hacks or implicit path rewriting.
  - Links in `README.md` and documentation must explicitly point to committed paths (e.g. `bundles/reactions/...`, `bundles/molecules/...`).
