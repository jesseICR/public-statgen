#!/usr/bin/env bash
#
# main.sh — Public Statistical Genetics Pipeline
# ====================================
# Orchestrates the full public statistical genetics pipeline.
#
# Usage:
#   bash main.sh
#
# Steps:
#   1. Install PLINK 1.9 & 2.0
#   2. Download KG, HGDP, SGDP, Neural ADMIXTURE, and GIAB data
#   3. QC KG and HGDP data
#   4. Set up Python virtual environment
#   5. QC SGDP data (liftover hg19->hg38, match to KG, assign rsIDs)
#   6. Merge KG + HGDP + SGDP into a single fileset
#   7. Prepare and merge GIAB Ashkenazi parents (HG003, HG004)
#   8. Build merged metadata CSV with Neural ADMIXTURE ancestry labels
#   9. Build supervised reference population assignments (K=6)
#  10. Install ADMIXTURE software
#  11. QC merged panel for ADMIXTURE (geno/MAF0.03/HWE/LD/mind/kinship)
#  12. Run ADMIXTURE supervised ancestry (3-fold CV + final projection)
#  13. Analyze ADMIXTURE results (structure plots, metadata, allele freqs)
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
export SUPERVISED_ADMIXTURE="${PROJECT_DIR}/supervised_admixture"
export ADMIXTURE="${TOOLS_BIN}/admixture"

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
# Step 7 — Prepare and merge GIAB Ashkenazi parents
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 7: Prepare and merge GIAB Ashkenazi parents"
echo "============================================"

# Check if GIAB samples are already in the merged fam
if grep -qw "HG003" "${MERGE_DIR}/merged_kg_hgdp_sgdp.fam" 2>/dev/null; then
    echo "  [skip] GIAB samples already in merged panel"
else
    "${PYTHON}" "${PROJECT_DIR}/prepare_giab.py"
    bash "${PROJECT_DIR}/merge_giab.sh"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 8 — Build metadata CSV
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 8: Build merged metadata CSV"
echo "============================================"

if [[ -f "${PROJECT_DIR}/summary/metadata.csv" ]]; then
    echo "  [skip] Metadata CSV already exists in ${PROJECT_DIR}/summary/"
else
    "${PYTHON}" "${PROJECT_DIR}/build_metadata.py"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 9 — Build supervised reference populations (K=6)
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 9: Build supervised reference populations (K=6)"
echo "============================================"

if [[ -f "${PROJECT_DIR}/summary/supervised.csv" ]]; then
    echo "  [skip] Supervised CSV already exists in ${PROJECT_DIR}/summary/"
else
    "${PYTHON}" "${PROJECT_DIR}/build_supervised.py"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 10 — Install ADMIXTURE
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 10: Install ADMIXTURE"
echo "============================================"
bash "${PROJECT_DIR}/setup_admixture.sh"

"${ADMIXTURE}" --version 2>&1 | head -1 || { echo "Error: admixture failed to run" >&2; exit 1; }
echo ""

# ---------------------------------------------------------------------------
# Step 11 — QC merged panel for ADMIXTURE (MAF 0.03)
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 11: QC merged panel for ADMIXTURE (MAF 0.03)"
echo "============================================"

if [[ -f "${SUPERVISED_ADMIXTURE}/ancestry_qc.bed" ]]; then
    echo "  [skip] ADMIXTURE QC already complete in ${SUPERVISED_ADMIXTURE}/"
else
    mkdir -p "${SUPERVISED_ADMIXTURE}/scrap"
    bash "${PROJECT_DIR}/qc_admixture.sh"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 12 — Run ADMIXTURE supervised ancestry (3-fold CV + final)
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 12: Run ADMIXTURE supervised ancestry estimation"
echo "============================================"

# Detect K from supervised.csv
K=$("${PYTHON}" -c "
import pandas as pd
print(pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')['reference_population'].nunique())
")

if [[ -f "${SUPERVISED_ADMIXTURE}/admixture_final.${K}.Q" ]]; then
    echo "  [skip] ADMIXTURE output already exists in ${SUPERVISED_ADMIXTURE}/"
else
    "${PYTHON}" "${PROJECT_DIR}/run_admixture_supervised.py"
fi
echo ""

# ---------------------------------------------------------------------------
# Step 13 — Analyze ADMIXTURE results
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 13: Analyze ADMIXTURE results"
echo "============================================"

if [[ -f "${PROJECT_DIR}/summary/admixture-global-${K}/metadata_ancestry.csv" ]]; then
    echo "  [skip] ADMIXTURE analysis output already exists"
else
    "${PYTHON}" "${PROJECT_DIR}/analyze_admixture_results.py"
fi
echo ""

echo "============================================"
echo "Pipeline steps 1-13 complete."
echo "============================================"
