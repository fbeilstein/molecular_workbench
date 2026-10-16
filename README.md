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

### 2. Multistep Reaction Pathway Studio
Build, align, and render complex multistep reaction mechanisms with synchronized 3D orbitals, Potential Energy Surfaces (PES), and interactive 2D SVG wiring:
```bash
python server_multistep.py --port 8083
```
Open [**http://localhost:8083**](http://localhost:8083) in your browser.
- **Stage 1**: Ingest reaction step JSONs or IRC trajectories.
- **Stage 2**: Automatic coordinate frame alignment across reaction steps.
- **Stage 3**: Harmonized energy profiling and activation barrier calculation.
- **Stage 4**: Synchronized orbital schedules and continuous Electrostatic Potential (ESP) surfaces.
- **Stage 5**: Interactive SVG mechanism wiring and maximal-compression `.rxnbundle.zip` exporter.

### 3. Transition State & IRC Workbench
Find reaction transition states, analyze vibrational normal modes, and integrate full IRC reaction pathways:
```bash
python server_irc.py --port 8082
```
Open [**http://localhost:8082**](http://localhost:8082) in your browser.
- **Step 1**: Draw or input reactant and product SMILES.
- **Step 2**: Interactive Ketcher atom mapping with automated graph validation.
- **Step 3**: 3D TS search (Genetic Swarm, Universal IDPP, ORCA, or RitS), imaginary frequency validation, and IRC integration.

### 4. Standalone 3D Web Viewer (`viewer.html`)
A standalone client-side viewer that unpacks `.rxnbundle.zip` archives directly in the browser's memory using `JSZip` and renders them with `3Dmol.js`:
- Open `viewer.html` directly or serve via GitHub Pages.
- Drag & drop any `.zip` bundle, or pass it via URL query parameter: `viewer.html?bundle=<path-or-url>`.

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
|![](https://github.com/fbeilstein/molecular_workbench/blob/master/bundles/reactions/tetrahedral_intermediate_CO.svg)|[View Reaction](https://fbeilstein.github.io/molecular_workbench/viewer.html?bundle=bundles/reactions/tetrahedral_intermediate_CO.zip)

</details>

<details>
<summary><strong>Molecules</strong></summary>

|Molecule|Name|Bundle|Notes|
|--|--|--|--|
|![](https://github.com/fbeilstein/molecular_workbench/blob/master/bundles/molecules/%5B18%5Dannulene.svg)|**Cyclooctadecanonaene** or **[18]annulene**|[View Molecule](https://fbeilstein.github.io/molecular_workbench/viewer.html?bundle=bundles/molecules/%5B18%5Dannulene.zip)|Aromatic: $4n+2$ electrons + planar. SMILES (NOT WORKING IN Ketcher!): C1([H])C([H])=C([H])C([H])=C([H])C([H])=C([H])C([H])=C([H])C([H])=C([H])C([H])=C([H])C([H])=C([H])C([H])=C([H])C=1[H] |c:7,11,19,23,31,34,t:3,15,27||
|![](https://github.com/fbeilstein/molecular_workbench/blob/master/bundles/molecules/all-cis-%5B10%5Dannulene.svg)|**Cyclodecapentaene** or **[10]annulene**, or **all-cis-[10]annulene**|[View Molecule](https://fbeilstein.github.io/molecular_workbench/viewer.html?bundle=bundles/molecules/all-cis-%5B10%5Dannulene.zip)|NOT Aromatic: $4n+2$ electrons BUT non-planar. SMILES: C1C=CC=CC=CC=CC=1 |c:1,3,5,7,9||

</details>
