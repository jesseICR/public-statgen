#!/usr/bin/env bash
#
# setup_admixture.sh
# ------------------
# Installs the ADMIXTURE 1.3.0 software for ancestry estimation.
#
# ADMIXTURE (Alexander et al. 2009) is a model-based ancestry estimation tool
# that uses maximum likelihood to estimate individual ancestry proportions from
# SNP genotype data.
#
# Outputs:
#   tools/bin/admixture — ADMIXTURE binary
#
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOLS_DIR="${PROJECT_DIR}/tools"
ADMIXTURE_DIR="${TOOLS_DIR}/admixture"
BIN_DIR="${TOOLS_DIR}/bin"

# ---------------------------------------------------------------------------
# Skip if already installed
# ---------------------------------------------------------------------------
if [[ -x "${BIN_DIR}/admixture" ]]; then
    echo "  [skip] ADMIXTURE already installed at ${BIN_DIR}/admixture"
    "${BIN_DIR}/admixture" --version 2>&1 | head -1 || true
    exit 0
fi

# ---------------------------------------------------------------------------
# Download ADMIXTURE 1.3.0
# ---------------------------------------------------------------------------
mkdir -p "${ADMIXTURE_DIR}" "${BIN_DIR}"

echo "  [download] Downloading ADMIXTURE 1.3.0 ..."
TARBALL="${ADMIXTURE_DIR}/admixture_linux-1.3.0.tar.gz"

curl -fSL -o "${TARBALL}" \
    "https://dalexander.github.io/admixture/binaries/admixture_linux-1.3.0.tar.gz"

tar -xzf "${TARBALL}" -C "${ADMIXTURE_DIR}" --strip-components=2
rm -f "${TARBALL}"

# Symlink into tools/bin
ln -sf "${ADMIXTURE_DIR}/admixture" "${BIN_DIR}/admixture"
chmod +x "${BIN_DIR}/admixture"

# Verify
if "${BIN_DIR}/admixture" --version 2>&1 | head -1; then
    echo "  [ok] ADMIXTURE installed successfully."
else
    echo "  Error: ADMIXTURE binary failed to run." >&2
    exit 1
fi
