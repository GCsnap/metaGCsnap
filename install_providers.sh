#!/usr/bin/env bash
# install_providers.sh – create the gcsnap conda/mamba environment and optionally
# add provider-specific dependencies.
#
# Usage:
#   bash install_providers.sh --base
#   bash install_providers.sh --ncbi
#   bash install_providers.sh --mgnify
#   bash install_providers.sh --local
#   bash install_providers.sh --complete
#
# A flag is always required.  The gcsnap base environment is always created
# first; provider flags add their extra packages on top via env update.
# The gcsnap package itself is installed in editable mode (pip install -e .)
# at the end regardless of which flag is used.
# mamba is used automatically if available, otherwise falls back to conda.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ENVS_DIR="${SCRIPT_DIR}/envs"

# ── usage ─────────────────────────────────────────────────────────────────────

usage() {
    echo ""
    echo "Usage: bash install_providers.sh <flag>"
    echo ""
    echo "  --base       Install base environment only"
    echo "  --ncbi       Install base + NCBI provider dependencies"
    echo "  --mgnify     Install base + MGnify provider dependencies"
    echo "  --local      Install base + local provider dependencies (mpi4py)"
    echo "  --complete   Install base + all provider dependencies"
    echo ""
    exit 1
}

# ── require exactly one flag ──────────────────────────────────────────────────

if [[ $# -ne 1 ]]; then
    echo "Error: exactly one flag is required."
    usage
fi

FLAG="$1"

case "$FLAG" in
    --base|--ncbi|--mgnify|--local|--complete) ;;
    *) echo "Error: unknown flag '$FLAG'."; usage ;;
esac

# ── helpers ───────────────────────────────────────────────────────────────────

# Auto-detect mamba or conda.
if command -v mamba &>/dev/null; then
    PKG_MANAGER="mamba"
elif command -v conda &>/dev/null; then
    PKG_MANAGER="conda"
else
    echo "Error: neither mamba nor conda found. Please install one and try again."
    exit 1
fi

echo "==> Using ${PKG_MANAGER}"

# Initialise the package manager so that 'activate' works inside the script.
PKG_BASE="$(${PKG_MANAGER} info --base 2>/dev/null)"
source "${PKG_BASE}/etc/profile.d/conda.sh"
if [[ "${PKG_MANAGER}" == "mamba" ]]; then
    source "${PKG_BASE}/etc/profile.d/mamba.sh" 2>/dev/null || true
fi

update_provider() {
    local yml="$1"
    echo "  -> updating with ${yml} ..."
    ${PKG_MANAGER} env update -n gcsnap -f "${yml}"
}

# ── step 1: base environment ──────────────────────────────────────────────────

echo ""
echo "==> Creating base environment from envs/GCsnap_base.yml ..."
${PKG_MANAGER} env create -f "${ENVS_DIR}/GCsnap_base.yml" || {
    echo "  (environment already exists – updating instead)"
    ${PKG_MANAGER} env update -n gcsnap -f "${ENVS_DIR}/GCsnap_base.yml" --prune
}

# ── step 2: provider extras ───────────────────────────────────────────────────

echo ""
echo "==> Adding provider dependencies (flag: ${FLAG}) ..."

case "$FLAG" in
    --base)
        echo "  -> base only, no provider extras."
        ;;
    --ncbi)
        update_provider "${ENVS_DIR}/provider_ncbi.yml"
        ;;
    --mgnify)
        update_provider "${ENVS_DIR}/provider_mgnify.yml"
        ;;
    --local)
        update_provider "${ENVS_DIR}/provider_local.yml"
        ;;
    --complete)
        update_provider "${ENVS_DIR}/provider_ncbi.yml"
        update_provider "${ENVS_DIR}/provider_mgnify.yml"
        update_provider "${ENVS_DIR}/provider_local.yml"
        ;;
esac

# ── step 3: install gcsnap in editable mode ───────────────────────────────────

echo ""
echo "==> Installing gcsnap in editable mode ..."
${PKG_MANAGER} run -n gcsnap pip install -e "${SCRIPT_DIR}"

echo ""
echo "==> Done.  Activate your environment with:"
echo "      ${PKG_MANAGER} activate gcsnap"
echo ""
