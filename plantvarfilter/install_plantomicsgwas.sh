#!/usr/bin/env bash
# ==============================================================================
# PlantOmicsGWAS - One-shot HPC installer
#
# What this does, in order:
#   1. Finds conda/mamba, or installs Miniforge if neither exists
#   2. Creates a dedicated environment with every external bioinformatics
#      tool the package needs (samtools, bcftools, bowtie2, minimap2,
#      plink, plink2, htslib) plus a build toolchain (cmake, compilers)
#   3. Installs the plantomicsgwas Python package + optional extras
#      (fastlmm, pysnptools, xgboost, geneview, cyvcf2)
#   4. Verifies every tool and the CLI actually work
#
# Usage:
#   bash install_plantomicsgwas.sh [environment_name]
#
# Default environment name: plantomicsgwas
# ==============================================================================

set -euo pipefail

ENV_NAME="${1:-plantomicsgwas}"
PYTHON_VERSION="3.11"   # chosen for best current compatibility with fastlmm

echo "=============================================="
echo " PlantOmicsGWAS - One-shot HPC installer"
echo " Environment name : $ENV_NAME"
echo " Python version   : $PYTHON_VERSION"
echo "=============================================="
echo ""

# ------------------------------------------------------------------
# 1. Find or install conda/mamba
# ------------------------------------------------------------------
if command -v mamba >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    CONDA_CMD=mamba
elif command -v conda >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base)"
    CONDA_CMD=conda
else
    echo "[INFO] No conda/mamba found on this system."
    echo "[INFO] Installing Miniforge into \$HOME/miniforge3 (no admin rights needed) ..."

    OS_NAME="$(uname)"
    ARCH_NAME="$(uname -m)"
    INSTALLER_URL="https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-${OS_NAME}-${ARCH_NAME}.sh"

    curl -L -o "$HOME/miniforge_installer.sh" "$INSTALLER_URL"
    bash "$HOME/miniforge_installer.sh" -b -p "$HOME/miniforge3"
    rm -f "$HOME/miniforge_installer.sh"

    CONDA_BASE="$HOME/miniforge3"
    CONDA_CMD=mamba
fi

# shellcheck disable=SC1090
source "${CONDA_BASE}/etc/profile.d/conda.sh"
command -v mamba >/dev/null 2>&1 || CONDA_CMD=conda

echo "[INFO] Using package manager : $CONDA_CMD"
echo "[INFO] Conda base            : $CONDA_BASE"
echo ""

# ------------------------------------------------------------------
# 2. Create (or reuse) the environment with all external tools
# ------------------------------------------------------------------
# NOTE (Issue 1.2 fix): cmake + a C/C++ compiler toolchain are included here
# so that any package needing a source build later (e.g. xgboost, fastlmm,
# pysnptools on platforms without a prebuilt wheel) has a working toolchain
# from the start, instead of failing mid-way through step 3 with something
# like "PermissionError: [Errno 13] Permission denied: 'cmake'".
if conda env list | grep -qE "^${ENV_NAME}[[:space:]]"; then
    echo "[INFO] Environment '$ENV_NAME' already exists — reusing it."
else
    echo "[INFO] Creating environment '$ENV_NAME' with bioinformatics tools..."
    "$CONDA_CMD" create -y -n "$ENV_NAME" \
        -c conda-forge -c bioconda \
        python="$PYTHON_VERSION" \
        cmake compilers \
        samtools bcftools bowtie2 minimap2 plink plink2 htslib
fi

conda activate "$ENV_NAME"
echo "[INFO] Environment activated: $(python3 --version)"
echo ""

# ------------------------------------------------------------------
# 3. Install the Python package + optional extras
# ------------------------------------------------------------------
echo "[INFO] Upgrading pip..."
pip install --upgrade pip --quiet

echo "[INFO] Installing plantomicsgwas (this may take a few minutes,"
echo "       fastlmm and other extras build from source)..."
pip install "plantomicsgwas[all]"

echo ""

# ------------------------------------------------------------------
# 4. Verify everything actually works
# ------------------------------------------------------------------
echo "=============================================="
echo " Verifying installation"
echo "=============================================="

python3 - <<'PYEOF'
try:
    from plantvarfilter.linux import resolve_tool
except ImportError:
    from plantvarfilter.windows import resolve_tool  # type: ignore

tools = ["bcftools", "samtools", "tabix", "bgzip", "plink", "plink2", "bowtie2", "minimap2"]
print()
for tool in tools:
    try:
        path = resolve_tool(tool)
    except Exception as exc:
        path = None
    status = "OK     " if path else "MISSING"
    print(f"  [{status}] {tool:10} -> {path}")
print()
PYEOF

echo "[INFO] Checking plantomicsgwas-compute CLI..."
if plantomicsgwas-compute list-steps >/tmp/plantomicsgwas_check.log 2>&1; then
    echo "  [OK     ] plantomicsgwas-compute list-steps"
else
    echo "  [FAILED ] plantomicsgwas-compute list-steps -- see /tmp/plantomicsgwas_check.log"
fi

echo ""
echo "=============================================="
echo " Done! Environment is ready."
echo "=============================================="
echo ""
echo "IMPORTANT for multi-node HPC clusters:"
echo "  This environment must be visible from every compute node that will"
echo "  run a job -- not just the login/head node. This is normally true"
echo "  automatically if \$HOME (or wherever you installed conda) is on a"
echo "  shared filesystem (NFS / Lustre / GPFS), which is the standard"
echo "  setup on most HPC clusters. If your cluster instead gives every"
echo "  node a separate local disk, re-run this script once per node, or"
echo "  ask your cluster admin how shared software installs are normally"
echo "  handled on your system."
echo ""
echo "For SLURM/PBS/LSF job scripts, add this to the config's hpc: section"
echo "so every submitted job activates the environment automatically:"
echo ""
echo "    hpc:"
echo "      conda_env: ${CONDA_BASE}/envs/${ENV_NAME}"
echo ""
echo "=============================================="
echo " Dropping you into a ready-to-use shell now..."
echo " (type 'exit' to leave this environment)"
echo "=============================================="
echo ""

# Land the person directly in an activated, ready-to-run shell instead of
# just printing instructions they'd have to copy/paste themselves.
exec bash --rcfile <(echo "source \"${CONDA_BASE}/etc/profile.d/conda.sh\"; conda activate ${ENV_NAME}; PS1='(${ENV_NAME}) \\u@\\h:\\w\\$ '")