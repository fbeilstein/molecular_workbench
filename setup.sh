#!/usr/bin/env bash
# ============================================================================
# Molecular Workbench: Automated Setup & Fine-Tuning
#
# Configures the local .venv, installs core & optional computational chemistry
# dependencies, downloads xTB binaries, tunes thread/runtime parameters,
# and performs automated verification.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
VENV_DIR="$SCRIPT_DIR/.venv"
TOOLS_DIR="$SCRIPT_DIR/tools"
XTB_DIR="$TOOLS_DIR/xtb"
OUTPUT_DIR="$SCRIPT_DIR/output"
ML_DIR="$TOOLS_DIR/ml_models/rits"

# UI Colors & Formatting
GREEN='\033[0;32m'
BLUE='\033[0;34m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
BOLD='\033[1m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
step()  { echo -e "\n${BLUE}${BOLD}==>${NC} ${BOLD}$1${NC}"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# Parse Command-line Options
WITH_ML=false
SKIP_TESTS=false

for arg in "$@"; do
    case "$arg" in
        --with-ml|--ml)
            WITH_ML=true
            ;;
        --skip-tests)
            SKIP_TESTS=true
            ;;
        -h|--help)
            echo "Molecular Workbench Setup"
            echo "Usage: ./setup.sh [OPTIONS]"
            echo ""
            echo "Options:"
            echo "  --with-ml     Also install PyTorch, PyG, and RitS package for ML transition states"
            echo "  --skip-tests  Skip post-installation sanity tests"
            echo "  -h, --help    Show this help message"
            exit 0
            ;;
        *)
            warn "Unknown option '$arg' (ignoring). Run with -h for help."
            ;;
    esac
done

echo -e "${BOLD}════════════════════════════════════════════════════════════════════${NC}"
echo -e "${BOLD}       Molecular Workbench: Environment Setup & Fine-Tuning       ${NC}"
echo -e "${BOLD}════════════════════════════════════════════════════════════════════${NC}"

# ── 1. Initialize or Detect Local .venv ──────────────────────────────────────
step "Step 1: Setting up local virtual environment (.venv)..."

if [ ! -d "$VENV_DIR" ]; then
    info "No .venv found. Creating fresh virtual environment at $VENV_DIR..."
    if ! command -v python3 >/dev/null 2>&1; then
        error "python3 is required to initialize .venv, but was not found on PATH."
    fi
    python3 -m venv "$VENV_DIR"
    info "✓ Created virtual environment."
else
    info "✓ Found existing virtual environment at $VENV_DIR."
fi

VENV_PY="$VENV_DIR/bin/python3"
VENV_PIP="$VENV_DIR/bin/pip"

if [ ! -x "$VENV_PY" ] || [ ! -x "$VENV_PIP" ]; then
    error "Virtual environment binaries missing at $VENV_DIR/bin."
fi

PY_VER=$("$VENV_PY" -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}")')
info "Using Python $PY_VER ($VENV_PY)"

info "Upgrading pip, setuptools, and wheel..."
"$VENV_PIP" install --quiet --upgrade pip setuptools wheel

# ── 2. Install Core Dependencies ─────────────────────────────────────────────
step "Step 2: Installing core dependencies from requirements.txt..."

if [ -f "$SCRIPT_DIR/requirements.txt" ]; then
    "$VENV_PIP" install -r "$SCRIPT_DIR/requirements.txt"
    info "✓ Core requirements installed."
else
    error "requirements.txt not found in $SCRIPT_DIR!"
fi

# ── 3. Semi-empirical Engine Setup (xTB) ─────────────────────────────────────
step "Step 3: Setting up xTB semi-empirical quantum engine..."

mkdir -p "$TOOLS_DIR"
if [ -x "$XTB_DIR/bin/xtb" ]; then
    info "✓ xTB standalone binary already present at $XTB_DIR/bin/xtb"
else
    info "Downloading precompiled Grimme xTB 6.7.1 binary for Linux x86_64..."
    mkdir -p "$XTB_DIR"
    TMP_ARCHIVE="/tmp/xtb-6.7.1.tar.xz"
    wget -qO "$TMP_ARCHIVE" "https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz"
    tar xf "$TMP_ARCHIVE" -C "$XTB_DIR" --strip-components=1
    rm -f "$TMP_ARCHIVE"
    chmod +x "$XTB_DIR/bin/xtb"
    info "✓ xTB installed successfully at $XTB_DIR/bin/xtb"
fi

# ── 4. Setup Directories & Frontend Libraries ────────────────────────────────
step "Step 4: Preparing directories & frontend libraries..."
mkdir -p "$OUTPUT_DIR"
mkdir -p "$ML_DIR"
mkdir -p "$SCRIPT_DIR/static_irc/libs"

if [ ! -f "$SCRIPT_DIR/static_irc/libs/3Dmol-min.js" ]; then
    info "Downloading 3Dmol.js (3Dmol-min.js) for offline molecular viewer..."
    curl -sSL "https://3Dmol.org/build/3Dmol-min.js" -o "$SCRIPT_DIR/static_irc/libs/3Dmol-min.js" || \
    wget -qO "$SCRIPT_DIR/static_irc/libs/3Dmol-min.js" "https://3Dmol.org/build/3Dmol-min.js" || \
    warn "Could not download 3Dmol-min.js offline; browser will use CDN fallback."
else
    info "✓ 3Dmol.js library present at static_irc/libs/3Dmol-min.js"
fi
info "✓ Verified directories and frontend assets."

# ── 5. Optional Machine Learning Stack (RitS) ────────────────────────────────
if [ "$WITH_ML" = true ]; then
    step "Step 5: Installing Machine Learning TS stack (PyTorch, PyG, RitS)..."
    
    # Check GPU availability
    if command -v nvidia-smi >/dev/null 2>&1; then
        info "NVIDIA GPU detected. Installing PyTorch..."
    else
        info "No NVIDIA GPU detected. Installing PyTorch CPU wheel..."
    fi
    "$VENV_PIP" install torch --index-url https://download.pytorch.org/whl/cpu
    "$VENV_PIP" install torch-geometric

    RITS_REPO_DIR="$TOOLS_DIR/rits_repo"
    if [ -d "$RITS_REPO_DIR" ]; then
        info "✓ Found local RitS repository at $RITS_REPO_DIR"
        info "Installing RitS (megalodon) package..."
        "$VENV_PIP" install --no-deps -e "$RITS_REPO_DIR"
        info "✓ Machine Learning dependencies installed."
    else
        info "Cloning isayevlab/RitS into $RITS_REPO_DIR..."
        if git clone https://github.com/isayevlab/RitS.git "$RITS_REPO_DIR"; then
            info "Installing RitS (megalodon) package..."
            "$VENV_PIP" install --no-deps -e "$RITS_REPO_DIR"
            info "✓ Machine Learning dependencies installed."
        else
            warn "Could not clone RitS from GitHub (network or firewall issue)."
            warn "Core TS engines (🧬 Genetic Swarm & 🧱 Dual-Endpoint IDPP) remain fully functional."
            warn "You can clone RitS manually when connected:"
            warn "  git clone https://github.com/isayevlab/RitS.git tools/rits_repo"
            warn "  pip install --no-deps -e tools/rits_repo"
        fi
    fi

    # Automated weights download via gdown
    "$VENV_PIP" install --quiet gdown
    RITS_CKPT_FILE="$ML_DIR/rits.ckpt"
    if [ -f "$RITS_CKPT_FILE" ] && [ -s "$RITS_CKPT_FILE" ]; then
        info "✓ RitS weights checkpoint already present at $RITS_CKPT_FILE"
    else
        info "Downloading RitS model weights checkpoint (~364 MB) via gdown..."
        if "$VENV_DIR/bin/gdown" "https://drive.google.com/uc?id=1tmCmIuzgDsZYg3DeMeA0fcKDqvbCyqye" -O "$RITS_CKPT_FILE"; then
            info "✓ Downloaded RitS weights checkpoint to $RITS_CKPT_FILE"
        else
            warn "Automated checkpoint download failed. You can download manually from:"
            warn "  https://drive.google.com/drive/folders/1xGX0uZgaGHU_e-a8hDdflV4wWApczNVU"
            warn "and place into $RITS_CKPT_FILE"
        fi
    fi
else
    info "Note: ML flow-matching dependencies (RitS) skipped (run with --with-ml to install)."
fi

# ── 6. Environment Fine-Tuning ───────────────────────────────────────────────
step "Step 6: Fine-tuning computational runtime environment..."

DOTENV_FILE="$SCRIPT_DIR/.env"
cat <<EOF > "$DOTENV_FILE"
# Molecular Workbench Environment Configuration
OMP_NUM_THREADS=2
MKL_NUM_THREADS=2
OPENBLAS_NUM_THREADS=2
XTBPATH=$XTB_DIR/share/xtb
PATH=$XTB_DIR/bin:\$PATH
PYTHONPATH=$SCRIPT_DIR:\$PYTHONPATH
RITS_CHECKPOINT=$ML_DIR/rits_checkpoint.pt
EOF
info "✓ Wrote tuned runtime parameters to .env"

# Inject into .venv/bin/activate if not already present
ACTIVATE_SCRIPT="$VENV_DIR/bin/activate"
if [ -f "$ACTIVATE_SCRIPT" ]; then
    if ! grep -q "Molecular Workbench Fine-Tuning" "$ACTIVATE_SCRIPT"; then
        cat <<'EOF' >> "$ACTIVATE_SCRIPT"

# --- Molecular Workbench Fine-Tuning ---
export OMP_NUM_THREADS=2
export MKL_NUM_THREADS=2
export OPENBLAS_NUM_THREADS=2
export XTBPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../tools/xtb/share/xtb" 2>/dev/null && pwd || true)"
export PATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../tools/xtb/bin" 2>/dev/null && pwd):$PATH"
export PYTHONPATH="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." 2>/dev/null && pwd):$PYTHONPATH"
export RITS_CHECKPOINT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../tools/ml_models/rits" 2>/dev/null && pwd)/rits_checkpoint.pt"
# ----------------------------------------
EOF
        info "✓ Applied environment fine-tuning hooks to .venv/bin/activate"
    fi
fi

# ── 7. Verification & Sanity Checks ─────────────────────────────────────────
if [ "$SKIP_TESTS" = false ]; then
    step "Step 7: Verifying fine-tuned environment and running sanity tests..."

    PYTHONPATH="$SCRIPT_DIR" "$VENV_PY" - <<'EOF'
import sys

# 1. RDKit Conformer Generation
from rdkit import Chem
from rdkit.Chem import AllChem
m = Chem.MolFromSmiles("CCO")
m = Chem.AddHs(m)
AllChem.EmbedMolecule(m, randomSeed=42)
assert m.GetNumConformers() > 0, "RDKit conformer embedding failed"
print("  ✓ RDKit 3D conformer embedding: OK")

# 2. PySCF Engine
import pyscf
print("  ✓ PySCF quantum chemistry engine: OK")

# 3. ASE & tblite (xTB)
from ase import Atoms
from tblite.ase import TBLite
h2o = Atoms("H2O", positions=[[0, 0, 0], [0, 0.757, 0.586], [0, -0.757, 0.586]])
h2o.calc = TBLite(method="GFN2-xTB")
e = h2o.get_potential_energy()
assert e < 0, "tblite energy evaluation failed"
print(f"  ✓ tblite (GFN2-xTB) calculation: OK ({e:.2f} eV)")

# 4. Sella Delocalized Internal Coordinate Optimizer
from sella import Sella
opt = Sella(h2o, internal=True, order=0, logfile=None)
print("  ✓ Sella DLC optimizer: OK")

# 5. Universal Transition State IDPP Interpolation Engine
from irc_plugins.ts_engine.lego_builder import build_lego_ts_guess
print("  ✓ Universal Dual-Endpoint IDPP TS engine: OK")
EOF
    info "✓ All sanity tests passed successfully!"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo -e "${GREEN}${BOLD}════════════════════════════════════════════════════════════════════${NC}"
echo -e "${GREEN}${BOLD}             Molecular Workbench Setup & Tuning Complete!           ${NC}"
echo -e "${GREEN}${BOLD}════════════════════════════════════════════════════════════════════${NC}"
echo ""
echo "To activate your fine-tuned environment, run:"
echo "  source .venv/bin/activate"
echo ""
echo "To launch the Reaction & Transition State Workbench (Step 1-3 + IRC):"
echo "  .venv/bin/python3 server_irc.py --port 8082"
echo "  (Open http://localhost:8082)"
echo ""
echo "To launch the Molecular Orbitals & ESP Surface Workbench:"
echo "  .venv/bin/python3 server.py --port 4321"
echo "  (Open http://localhost:4321)"
echo ""
if [ "$WITH_ML" = false ]; then
    echo "Note: To install optional ML flow-matching models (RitS), run:"
    echo "  ./setup.sh --with-ml"
    echo "See docs/ML_MODELS.md for model checkpoint weights download."
fi
echo ""
