# Molecular Workbench

A comprehensive computational chemistry suite for exploring molecular structures, quantum orbitals, transition states, reaction paths (IRC), and generating lightweight, zero-backend 3D visualizations (`.rxnbundle.zip`).

---

## Prerequisites

- **Python**: Version **3.9 to 3.12** is required (to ensure compatibility with PySCF, tblite, and ASE).
- **ORCA** (Optional): Version 5.x / 6.x if using ORCA-driven transition state and IRC calculations.

---

## Installation

This repository respects your preferred Python environment (e.g. `venv`, `conda`).

1. **Clone the repository**:
   ```bash
   git clone git@github.com:fbeilstein/molecular_workbench.git
   cd molecular_workbench
   ```

2. **Create and activate your environment**:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. **Install Dependencies & Tooling**:
   ```bash
   pip install -r requirements.txt
   bash setup.sh --with-ml
   ```
   *This installs RDKit, ASE, Sella, PySCF, tblite (xTB) and downloads the required `xtb` engine binary into the `tools/` folder.*

---

## Running the Workbenches

The repository provides three dedicated workbenches and a zero-backend standalone viewer:

### 1. Single-Molecule Quantum Workbench
Compute 3D geometries, Intrinsic Bond Orbitals (IBO), and Electrostatic Potential (ESP) maps for individual molecules:
```bash
python server.py --port 4321
```
Open [**http://localhost:4321**](http://localhost:4321) in your browser.


### 2. Transition State & IRC Workbench
Find reaction transition states, analyze vibrational normal modes, integrate full IRC reaction pathways, calculate orbitals, and prepare bundle:
```bash
python server_irc.py --port 8082
```
Open [**http://localhost:8082**](http://localhost:8082) in your browser.


### 3. Standalone 3D Web Viewer (`viewer.html`)
A standalone client-side viewer that unpacks `.rxnbundle.zip` archives directly in the browser's memory using `JSZip` and renders them with `3Dmol.js`.

- You can simply [OPEN MOLECULAR VIEWER](https://fbeilstein.github.io/molecular_workbench/viewer.html) in browser, choose zip bundle, and view.
- You can embed it into your page, the bundle can be passed as `https://fbeilstein.github.io/molecular_workbench/viewer.html?bundle=<path-or-url>`. Be aware of COARS.
- You can simply [DOWNLOAD IT](https://github.com/fbeilstein/molecular_workbench/releases/download/irc_v1/viewer.zip) in archive, unpack, and open viewer.html in any browser.

---

## Documentation

- [Transition State Engines & Search Strategies](docs/ts_engines.md): Details on Genetic Swarm, Dual-Endpoint IDPP, and Sella DLC optimization.
- [Machine Learning TS Models (RitS & TS-Diff)](docs/ml_models.md): Setup instructions and weight checkpoint downloads for generative flow-matching TS models.
- [Architecture Overview](docs/architecture.md): High-level system structure, components, and data pipelines.
- [Scientific Choices & Quantum Methodology](docs/scientific_choices.md): Level of theory, Intrinsic Bond Orbitals (IBO), and exact virtual subspace antisymmetrization.
- [Agent Knowledge & Technical Gotchas](docs/agent_knowledge.md): Developer notes, PySCF symmetry gotchas, and bundle structures.

---

## Gallery

<details>
<summary><strong>Reactions</strong></summary>

|Reaction|Bundle|
|---|---|
|![](https://github.com/fbeilstein/molecular_workbench/blob/master/bundles/reactions/tetrahedral_intermediate_CO.svg)|[View Reaction](https://fbeilstein.github.io/molecular_workbench/viewer.html?bundle=bundles/reactions/tetrahedral_intermediate_CO.zip&download)

</details>

[Molecules](https://fbeilstein.github.io/molecular_workbench/gallery.html)

