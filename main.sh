#!/usr/bin/env bash
#
# main.sh — Public Statistical Genetics Pipeline
# ====================================
# Orchestrates the full public statistical genetics pipeline.
#
# Usage:
#   bash main.sh
#
# Steps (current):
#   1. Install PLINK 1.9 & 2.0
#   2. Download KG, HGDP, SGDP, and Neural ADMIXTURE data
#   3. QC KG and HGDP data
#   4. Set up Python virtual environment
#   5. QC SGDP data (liftover hg19→hg38, match to KG, assign rsIDs)
#   6. Merge KG + HGDP + SGDP into a single fileset
#   7. Build merged metadata CSV with Neural ADMIXTURE ancestry labels
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
TOOLS_BIN="${PROJECT_DIR}/tools/bin"

export DOWNLOADS_DIR="${PROJECT_DIR}/downloads"
export QC_DIR="${PROJECT_DIR}/qc"
export MERGE_DIR="${PROJECT_DIR}/merge"
export PLINK1="${TOOLS_BIN}/plink1"
export PLINK2="${TOOLS_BIN}/plink2"
export SNPS_FILE="${PROJECT_DIR}/rsids_dense_chr1_22.txt"
export PLINK_MEMORY=14000
export PLINK_THREADS=6
export GENO_THRESHOLD=0.03
export PYTHON="${PROJECT_DIR}/tools/venv/bin/python"

# ---------------------------------------------------------------------------
# Step 1 — Install PLINK
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 1: Install PLINK"
echo "============================================"
bash "${PROJECT_DIR}/setup_plink.sh"

# Verify the binaries are functional
"${PLINK1}" --version > /dev/null 2>&1 || { echo "Error: plink1 failed to run" >&2; exit 1; }
"${PLINK2}" --version > /dev/null 2>&1 || { echo "Error: plink2 failed to run" >&2; exit 1; }
echo ""

# ---------------------------------------------------------------------------
# Step 2 — Download reference panel data
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 2: Download KG, HGDP, SGDP, and Neural ADMIXTURE data"
echo "============================================"
mkdir -p "${DOWNLOADS_DIR}"
bash "${PROJECT_DIR}/download_files.sh"
echo ""

# ---------------------------------------------------------------------------
# Step 3 — QC KG and HGDP
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 3: QC KG and HGDP data"
echo "============================================"

if [[ -f "${QC_DIR}/kg_qc.bed" && -f "${QC_DIR}/hgdp_qc.bed" ]]; then
    echo "  [skip] QC output already exists in ${QC_DIR}/"
else
    mkdir -p "${QC_DIR}"
    bash "${PROJECT_DIR}/qc_kg_hgdp.sh"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 4 — Set up Python virtual environment
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 4: Set up Python environment"
echo "============================================"
bash "${PROJECT_DIR}/setup_python.sh"
echo ""

# ---------------------------------------------------------------------------
# Step 5 — QC SGDP (liftover + rsID matching)
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 5: QC SGDP data"
echo "============================================"

if [[ -f "${QC_DIR}/sgdp_qc.bed" ]]; then
    echo "  [skip] SGDP QC output already exists in ${QC_DIR}/"
else
    mkdir -p "${QC_DIR}"

    # Unzip bim and copy bed/fam into QC directory
    unzip -p "${DOWNLOADS_DIR}/sgdp_all.bim.zip" > "${QC_DIR}/sgdp_all.bim"
    cp -n "${DOWNLOADS_DIR}/sgdp_all.bed" "${QC_DIR}/"
    cp -n "${DOWNLOADS_DIR}/sgdp_all.fam" "${QC_DIR}/"

    "${PYTHON}" "${PROJECT_DIR}/qc_sgdp.py"

    # Clean up raw SGDP files from QC directory
    rm -f "${QC_DIR}"/sgdp_all.*
fi
echo ""

# ---------------------------------------------------------------------------
# Step 6 — Merge KG + HGDP + SGDP
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 6: Merge KG + HGDP + SGDP"
echo "============================================"

if [[ -f "${MERGE_DIR}/merged_kg_hgdp_sgdp.bed" ]]; then
    echo "  [skip] Merged fileset already exists in ${MERGE_DIR}/"
else
    mkdir -p "${MERGE_DIR}"
    bash "${PROJECT_DIR}/merge_kg_hgdp_sgdp.sh"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 7 — Build metadata CSV
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 7: Build merged metadata CSV"
echo "============================================"

if [[ -f "${PROJECT_DIR}/summary/metadata.csv" ]]; then
    echo "  [skip] Metadata CSV already exists in ${PROJECT_DIR}/summary/"
else
    "${PYTHON}" "${PROJECT_DIR}/build_metadata.py"
fi
echo ""

echo "============================================"
echo "Pipeline steps 1-7 complete."
echo "============================================"
