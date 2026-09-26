#!/usr/bin/env bash
# ============================================================================
# Biochemistry Workbench Setup
# Creates a virtualenv and installs all required Python packages for:
#   - 3D molecular geometry generation (RDKit, xTB)
#   - Molecular orbital calculation (PySCF)
#   - Reaction trajectory generation (xTB path, GSM, ORCA fallback)
#   - ChemDraw export (via RDKit)
#
# Usage:
#   bash workbench/setup.sh            # First install
#   bash workbench/setup.sh --force    # Recreate venv from scratch
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
VENV_DIR="$PROJECT_DIR/.venv"
TOOLS_DIR="$SCRIPT_DIR/tools"
XTB_DIR="$TOOLS_DIR/xtb"
GSM_DIR="$TOOLS_DIR/gsm"
OUTPUT_DIR="$SCRIPT_DIR/output"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# ── Check prerequisites ─────────────────────────────────────────────────────
command -v python3 >/dev/null 2>&1 || error "python3 not found. Install Python 3.9+ first."

PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
info "Found Python $PYTHON_VERSION"

# ── Create virtual environment ───────────────────────────────────────────────
FORCE=false
[[ "${1:-}" == "--force" ]] && FORCE=true

if [ -d "$VENV_DIR" ]; then
    if [ "$FORCE" = true ]; then
        rm -rf "$VENV_DIR"
        python3 -m venv "$VENV_DIR"
        info "Recreated virtual environment"
    else
        info "Using existing virtual environment at $VENV_DIR"
    fi
else
    python3 -m venv "$VENV_DIR"
    info "Created virtual environment at $VENV_DIR"
fi

# Activate
source "$VENV_DIR/bin/activate"
pip install --upgrade pip setuptools wheel -q

# ── Create output directory ──────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"

# ── Install packages ─────────────────────────────────────────────────────────

info "═══ Stage 1: Core cheminformatics (RDKit) ═══"
pip install rdkit -q
info "  ✓ RDKit installed"

info "═══ Stage 2: Quantum chemistry (PySCF) ═══"
pip install pyscf -q
info "  ✓ PySCF installed"

info "═══ Stage 3: xTB semi-empirical engine ═══"
if [ -f "$XTB_DIR/bin/xtb" ]; then
    info "  ✓ xTB already present at $XTB_DIR/bin/xtb"
else
    info "  → Downloading xTB binary..."
    mkdir -p "$XTB_DIR"
    wget -qO /tmp/xtb.tar.xz "https://github.com/grimme-lab/xtb/releases/download/v6.7.1/xtb-6.7.1-linux-x86_64.tar.xz"
    tar xf /tmp/xtb.tar.xz -C "$XTB_DIR" --strip-components=1
    rm /tmp/xtb.tar.xz
    info "  ✓ xTB installed at $XTB_DIR/bin/xtb"
fi

info "═══ Stage 4: GSM reaction path engine ═══"
if [ -f "$GSM_DIR/gsm.orca" ]; then
    info "  ✓ GSM already present at $GSM_DIR/gsm.orca"
else
    info "  → Downloading GSM binary from Grimme fork..."
    mkdir -p "$GSM_DIR"
    wget -qO /tmp/gsm.tar.xz "https://github.com/grimme-lab/molecularGSM/releases/download/rev1/gsm.tar.xz"
    tar xf /tmp/gsm.tar.xz -C "$GSM_DIR" --strip-components=1
    rm /tmp/gsm.tar.xz
    chmod +x "$GSM_DIR/gsm.orca"
    info "  ✓ GSM installed at $GSM_DIR/gsm.orca"
fi

info "═══ Stage 5: Transition state engine (autodE) ═══"
pip install git+https://github.com/duartegroup/autodE.git -q 2>/dev/null || warn "  ⚠ autodE install failed. TS search will be unavailable."
info "  ✓ autodE installed (or skipped)"

# NWChem — open-source high-level QM fallback (optional)
if command -v apt-get >/dev/null 2>&1; then
    if ! command -v nwchem >/dev/null 2>&1; then
        info "  → Installing NWChem (open-source DFT engine, requires sudo)"
        sudo apt-get install -y nwchem -q 2>/dev/null || warn "  ⚠ NWChem install failed."
    fi
fi

info "═══ Stage 6: Utilities ═══"
pip install numpy matplotlib -q
info "  ✓ NumPy, Matplotlib installed"

# ── Verify ───────────────────────────────────────────────────────────────────
echo ""
info "═══ Verification ═══"

python3 -c "from rdkit import Chem; print('  ✓ RDKit', Chem.rdBase.rdkitVersion)" 2>/dev/null || warn "  ✗ RDKit import failed"
python3 -c "import pyscf; print('  ✓ PySCF', pyscf.__version__)" 2>/dev/null || warn "  ✗ PySCF import failed"

if [ -x "$XTB_DIR/bin/xtb" ]; then
    XTB_VER=$("$XTB_DIR/bin/xtb" --version 2>&1 | head -1 || echo "unknown")
    info "  ✓ xTB: $XTB_VER"
else
    warn "  ✗ xTB binary not found or not executable"
fi

if [ -x "$GSM_DIR/gsm.orca" ]; then
    info "  ✓ GSM (gsm.orca) found"
else
    warn "  ✗ GSM binary not found or not executable"
fi

if [ -f "$GSM_DIR/tm2orca.py" ]; then
    info "  ✓ tm2orca.py found"
else
    warn "  ✗ tm2orca.py not found — GSM won't work"
fi

if command -v nwchem >/dev/null 2>&1; then
    info "  ✓ NWChem found"
fi

if [ -x "/data/orca/orca" ]; then
    info "  ✓ ORCA found at /data/orca/orca"
else
    warn "  ⚠ ORCA not found — overnight IRC scripts will not be runnable"
fi

# ── Summary ──────────────────────────────────────────────────────────────────
echo ""
echo "═══════════════════════════════════════════════════════════════"
info "Setup complete!"
echo ""
echo "  Reaction path cascade:"
echo "    Tier 1: xTB --path     (fast, ~1 min)"
echo "    Tier 2: Relaxed scan   (driving coordinates, ~2 min)"
echo "    Tier 3: DE-GSM + xTB   (robust fallback, ~5 min)"
echo ""
echo "  Activate the environment:"
echo "    source $VENV_DIR/bin/activate"
echo ""
echo "  Launch the workbench:"
echo "    python $SCRIPT_DIR/server.py"
echo "    → http://localhost:8080"
echo ""
echo "═══════════════════════════════════════════════════════════════"
