#!/usr/bin/env bash
#
# setup_python.sh
# ----------------
# Creates a Python virtual environment in ./tools/venv/ and installs
# pipeline dependencies from requirements.txt.
#
# Outputs:
#   tools/venv/   — Python virtual environment with dependencies installed
#
set -euo pipefail

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
VENV_DIR="${PROJECT_DIR}/tools/venv"
REQUIREMENTS="${PROJECT_DIR}/requirements.txt"

# ---------------------------------------------------------------------------
# Create venv (skip if it already exists)
# ---------------------------------------------------------------------------
if [[ -d "${VENV_DIR}" ]]; then
    echo "  [skip] Python venv already exists at ${VENV_DIR}"
else
    echo "  [create] Python venv at ${VENV_DIR} ..."
    python3 -m venv "${VENV_DIR}"
fi

# ---------------------------------------------------------------------------
# Install / update dependencies
# ---------------------------------------------------------------------------
echo "  [install] Python dependencies from requirements.txt ..."
"${VENV_DIR}/bin/pip" install --quiet --upgrade pip
"${VENV_DIR}/bin/pip" install --quiet -r "${REQUIREMENTS}"

echo "==> Python setup complete."
