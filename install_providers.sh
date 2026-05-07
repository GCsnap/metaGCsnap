#!/usr/bin/env bash
# install_providers.sh – create the gcsnap conda environment and optionally
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
# first; provider flags add their extra packages on top via conda env update.
# The gcsnap package itself is installed in editable mode (pip install -e .)
# at the end regardless of which flag is used.

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

# Initialise conda so that 'conda activate' works inside the script.
CONDA_BASE="$(conda info --base 2>/dev/null)" || {
    echo "Error: conda not found. Please install conda/mamba and try again."
    exit 1
}
source "${CONDA_BASE}/etc/profile.d/conda.sh"

update_provider() {
    local yml="$1"
    echo "  -> updating with ${yml} ..."
    conda env update -n gcsnap -f "${yml}"
}

# ── step 1: base environment ──────────────────────────────────────────────────

echo ""
echo "==> Creating base environment from envs/GCsnap_base.yml ..."
conda env create -f "${ENVS_DIR}/GCsnap_base.yml" || {
    echo "  (environment already exists – updating instead)"
    conda env update -n gcsnap -f "${ENVS_DIR}/GCsnap_base.yml" --prune
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
conda run -n gcsnap pip install -e "${SCRIPT_DIR}"

echo ""
echo "==> Done.  Activate your environment with:"
echo "      conda activate gcsnap"
echo ""
