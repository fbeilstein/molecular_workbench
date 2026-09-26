#!/usr/bin/env bash
# ============================================================================
# Molecular Workbench Setup
# Installs required Python packages and downloads the xTB binary.
#
# Note: Please activate your preferred Python environment (venv, conda, etc.) 
# before running this script.
# ============================================================================

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
TOOLS_DIR="$SCRIPT_DIR/tools"
XTB_DIR="$TOOLS_DIR/xtb"
OUTPUT_DIR="$SCRIPT_DIR/output"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

info()  { echo -e "${GREEN}[INFO]${NC} $1"; }
warn()  { echo -e "${YELLOW}[WARN]${NC} $1"; }
error() { echo -e "${RED}[ERROR]${NC} $1"; exit 1; }

# ── Check Python environment ────────────────────────────────────────────────
if ! command -v python3 >/dev/null 2>&1; then
    error "python3 not found."
fi

PYTHON_VERSION=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
info "Using Python $PYTHON_VERSION from $(which python3)"

# ── Install dependencies ────────────────────────────────────────────────────
if [ -f "$SCRIPT_DIR/requirements.txt" ]; then
    info "Installing Python dependencies from requirements.txt..."
    python3 -m pip install -r "$SCRIPT_DIR/requirements.txt"
else
    warn "requirements.txt not found! Skipping Python dependency installation."
fi

# ── Create output directory ──────────────────────────────────────────────────
mkdir -p "$OUTPUT_DIR"

# ── Download xTB ─────────────────────────────────────────────────────────────
info "═══ Setting up xTB semi-empirical engine ═══"
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

# ── Verify ───────────────────────────────────────────────────────────────────
echo ""
info "═══ Verification ═══"

python3 -c "import rdkit; print('  ✓ RDKit installed')" 2>/dev/null || warn "  ✗ RDKit import failed"
python3 -c "import pyscf; print('  ✓ PySCF installed')" 2>/dev/null || warn "  ✗ PySCF import failed"

if [ -x "$XTB_DIR/bin/xtb" ]; then
    XTB_VER=$("$XTB_DIR/bin/xtb" --version 2>&1 | head -1 || echo "unknown")
    info "  ✓ xTB: $XTB_VER"
else
    warn "  ✗ xTB binary not found or not executable"
fi

echo ""
info "Setup complete! You can now start the server:"
echo "  python3 server.py --port 4321"
