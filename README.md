# Molecular Workbench

A comprehensive computational chemistry suite for exploring molecular structures, orbitals, reaction paths (IRC), and generating embeddable 3D visualizations (`.rxnbundle`).

---

## Prerequisites

- **Python**: Version **3.9 to 3.12** is required (to ensure compatibility with PySCF, tblite, and ASE).
- **ORCA** (Optional): Version 5.x / 6.x if using ORCA-driven transition state and IRC calculations.

---

## Installation

This repository respects your preferred Python environment (e.g. `venv`, `conda`).

1. **Clone the repository**:
   ```bash
   git clone <your-repo-url>
   cd molecular_workbench
   ```

2. **Create and activate your environment**:
   ```bash
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. **Install Core Dependencies**:
   ```bash
   pip install -r requirements.txt
   bash setup.sh --with-ml
   ```
   *This installs RDKit, ASE, Sella, PySCF, tblite (xTB) and downloads the required `xtb` engine binary into the `tools/` folder.*

---

## Transition State Engines & Strategies

The Transition State Workbench (`server_irc.py`) supports three robust, transparent TS search strategies:

1. **🧬 Genetic Swarm (Default Intermolecular)**:
   - A multi-generation genetic algorithm exploring $SE(3) \times SO(3)$ (both 3D translation and 3D rigid-body quaternions for both reactants).
   - Followed by Stage 2 Quantum Valley relaxation and Stage 3 Sella Delocalized Internal Coordinate (DLC) saddle polish.

2. **🧱 Dual-Endpoint IDPP (Universal All Reaction Classes)**:
   - Image-Dependent Pair Potential (IDPP) geometric interpolation between Reactant ($R$) and Product ($P$) endpoints.
   - **Universal**: Operates on any number of forming and breaking bonds simultaneously (substitutions, Diels-Alder cycloadditions, additions, eliminations, and rearrangements) without reaction-specific heuristics.

3. **⚡ Machine Learning Flow-Matching (RitS / TS-Diff)**:
   - **Prerequisites Required**: Flow-matching neural network models require PyTorch, PyTorch Geometric, and external model weight checkpoints (~1.5 GB).
   - **Zero Silent Failures**: If selected without installed packages or checkpoint files, the engine **never silently fakes or falls back to ground states**. It immediately raises an informative error pointing to setup instructions.
   - For full setup and weights download instructions, see [docs/ML_MODELS.md](docs/ML_MODELS.md).

---

## Usage

### 1. IRC & Transition State Workbench
To run the full 3-step Reaction & Transition State Workbench:
```bash
source .venv/bin/activate
python server_irc.py --port 8082
```
Open **http://localhost:8082** in your browser.
- **Step 1**: Draw or input reactant and product SMILES.
- **Step 2**: Interactive Ketcher atom mapping with automated validation.
- **Step 3**: 3D Transition State Search (Genetic Swarm, IDPP, ORCA, or RitS) with real-time vibrational mode analysis and IRC path generation.

### 2. Orbital & Surface Workbench
To run the orbital and electrostatic potential (ESP) workbench:
```bash
source .venv/bin/activate
python server.py --port 4321
```
Open **http://localhost:4321** in your browser.

---

## Project Structure

- `server_irc.py`: IRC & Transition State HTTP server.
- `irc_plugins/`: Reaction model, mapping engine, TS optimization engine, and mode analysis.
  - `irc_plugins/ts_engine/`: Genetic swarm, Universal IDPP interpolator, RitS wrapper, calculators.
- `static_irc/`: Three.js 3D viewer, docking controls, real-time optimization log, and HUD.
- `docs/ML_MODELS.md`: Step-by-step installation instructions for optional ML models and weight checkpoints.
- `tools/`: Computational backend binaries and helper scripts.
