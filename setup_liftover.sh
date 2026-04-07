#!/usr/bin/env bash
#
# setup_liftover.sh
# -----------------
# Downloads and installs the UCSC liftOver binary into ./tools/bin/.
# Works on macOS (Intel & Apple Silicon) and Linux (x86_64).
#
# Outputs:
#   tools/bin/liftOver   — UCSC liftOver binary
#
set -euo pipefail

TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)/tools"
BIN_DIR="${TOOLS_DIR}/bin"

mkdir -p "${BIN_DIR}"

if [[ -x "${BIN_DIR}/liftOver" ]]; then
    echo "  [skip] liftOver already installed at ${BIN_DIR}/liftOver"
    exit 0
fi

OS="$(uname -s)"
ARCH="$(uname -m)"

case "${OS}-${ARCH}" in
    Darwin-arm64)   PLATFORM="macOSX.arm64" ;;
    Darwin-x86_64)  PLATFORM="macOSX.x86_64" ;;
    Linux-x86_64)   PLATFORM="linux.x86_64" ;;
    *)
        echo "Error: unsupported OS/architecture '${OS}-${ARCH}' for liftOver" >&2
        exit 1
        ;;
esac

URL="https://hgdownload.soe.ucsc.edu/admin/exe/${PLATFORM}/liftOver"

echo "  [download] liftOver ..."
curl -fSL -o "${BIN_DIR}/liftOver" "${URL}"
chmod +x "${BIN_DIR}/liftOver"

echo "  [installed] ${BIN_DIR}/liftOver"
