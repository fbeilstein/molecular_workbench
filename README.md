# Molecular Workbench

A comprehensive suite for exploring molecular structures, orbitals, reaction paths, and generating embeddable 3D visualizations (`.rxnbundle`).

## Prerequisites

- **Python**: Version **3.9 to 3.12** is required (to ensure compatibility with PySCF and other numerical libraries).
- **ORCA** (Optional but recommended for reaction path fallback): ORCA is a quantum chemistry program suite.
  - Download: [ORCA Forum](https://orcaforum.kofo.mpg.de/) (Registration required).
  - Setup: Ensure the `orca` binary is available in your system `$PATH` (e.g. `/data/orca/orca`).

## Installation

This repository uses a clean installation philosophy that respects your preferred Python environment (e.g., conda, venv, poetry).

1. **Clone the repository**:
   ```bash
   git clone <your-repo-url>
   cd molecular_workbench
   ```

2. **Activate your preferred Python environment**:
   ```bash
   # Example using venv:
   python3 -m venv .venv
   source .venv/bin/activate
   ```

3. **Run the Setup Script**:
   ```bash
   bash setup.sh
   ```
   *This script simply runs `pip install -r requirements.txt` into your active environment and downloads the required `xtb` engine binary into the `tools/` folder. It does not bloat your system with unnecessary defaults.*

## Usage

### Starting the Server

The application runs a local web server to provide an interactive 3D workbench UI. To start it:

```bash
# Start the server on the default port (4321)
python3 server.py --port 4321
```

Once running, open your web browser and navigate to:
**http://localhost:4321**

### Features

- **3D Viewer**: Uses 3Dmol.js to render molecular structures and volumetric data.
- **Orbital Computation**: Generates MOs and ESP surfaces using PySCF.
- **Reaction Paths**: Explores transition states and reaction coordinates using XTB/GSM.
- **Packaging**: Once you analyze a structure, use the "Download Bundle" button in the UI to generate a lightweight `.rxnbundle.zip` containing your selected ESPs and Orbitals. This bundle can be directly embedded into slide presentations.

## Project Structure

- `server.py`: Main HTTP server and API endpoints.
- `static/`: Frontend JS, CSS, and 3Dmol UI logic.
- `tools/`: Computational backend scripts (PySCF wrappers, XTB pathing, format conversion).
- `output/`: Generated `.xyz`, `.cube`, and `.rxnbundle` data. (Ignored by git)
