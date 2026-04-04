#!/usr/bin/env bash
#
# setup_plink.sh
# --------------
# Downloads and installs PLINK 1.9 and PLINK 2.0 into ./tools/bin/.
# Works on macOS (Intel & Apple Silicon) and Linux (x86_64).
#
# Outputs:
#   tools/bin/plink1   — PLINK 1.9 binary
#   tools/bin/plink2   — PLINK 2.0 binary
#
set -euo pipefail

PLINK1_DATE="20250819"
PLINK2_DATE="20260311"

TOOLS_DIR="$(cd "$(dirname "$0")" && pwd)/tools"
BIN_DIR="${TOOLS_DIR}/bin"
TMP_DIR="${TOOLS_DIR}/tmp"

mkdir -p "${BIN_DIR}" "${TMP_DIR}"

OS="$(uname -s)"
ARCH="$(uname -m)"

# ---------------------------------------------------------------------------
# Determine download URLs
# ---------------------------------------------------------------------------
case "${OS}" in
    Darwin) PLINK1_URL="https://s3.amazonaws.com/plink1-assets/plink_mac_${PLINK1_DATE}.zip" ;;
    Linux)  PLINK1_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_${PLINK1_DATE}.zip" ;;
    *)      echo "Error: unsupported OS '${OS}' for PLINK 1.9" >&2; exit 1 ;;
esac

case "${OS}-${ARCH}" in
    Darwin-arm64)   PLINK2_URL="https://s3.amazonaws.com/plink2-assets/plink2_mac_arm64_${PLINK2_DATE}.zip" ;;
    Darwin-x86_64)  PLINK2_URL="https://s3.amazonaws.com/plink2-assets/plink2_mac_avx2_${PLINK2_DATE}.zip" ;;
    Linux-x86_64)   PLINK2_URL="https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_${PLINK2_DATE}.zip" ;;
    Linux-aarch64)
        echo "Error: PLINK 2.0 does not provide a Linux aarch64 binary." >&2
        echo "See https://www.cog-genomics.org/plink/2.0/ for available builds." >&2
        exit 1
        ;;
    *) echo "Error: unsupported OS/architecture '${OS}-${ARCH}' for PLINK 2.0" >&2; exit 1 ;;
esac

# ---------------------------------------------------------------------------
# Install helper
# ---------------------------------------------------------------------------
install_plink() {
    local name="$1" url="$2" src_bin="$3" dest_bin="$4"

    if [[ -x "${dest_bin}" ]]; then
        echo "  [skip] ${name} already installed at ${dest_bin}"
        return 0
    fi

    echo "  [download] ${name} ..."
    curl -fSL -o "${TMP_DIR}/${name}.zip" "${url}"
    python3 -c "
import zipfile, os, sys
with zipfile.ZipFile(sys.argv[1]) as z:
    os.makedirs(sys.argv[2], exist_ok=True)
    for m in z.infolist():
        if not m.is_dir():
            m.filename = os.path.basename(m.filename)
            z.extract(m, sys.argv[2])
" "${TMP_DIR}/${name}.zip" "${TMP_DIR}/${name}"
    cp "${TMP_DIR}/${name}/${src_bin}" "${dest_bin}"
    chmod +x "${dest_bin}"
    echo "  [installed] ${dest_bin}"
}

# ---------------------------------------------------------------------------
# Download and install
# ---------------------------------------------------------------------------
echo "==> Installing PLINK binaries ..."

install_plink "plink1" "${PLINK1_URL}" "plink"  "${BIN_DIR}/plink1"
install_plink "plink2" "${PLINK2_URL}" "plink2" "${BIN_DIR}/plink2"

# ---------------------------------------------------------------------------
# Clean up
# ---------------------------------------------------------------------------
rm -rf "${TMP_DIR}"

echo "==> PLINK setup complete. Binaries in ${BIN_DIR}/"
