#!/usr/bin/env bash
#
# main.sh — Public Statistical Genetics Pipeline
# ====================================
# Orchestrates the full public statistical genetics pipeline.
#
# Usage:
#   bash main.sh
#
# The pipeline will prompt for:
#   - ADMIXTURE K model (K=3, K=5, or K=6)
#   - Minor allele frequency (MAF) threshold
#
# Steps:
#   1. Install PLINK 1.9 & 2.0 and UCSC liftOver
#   2. Download KG, HGDP, SGDP, Neural ADMIXTURE, GIAB, and population data
#   3. QC KG and HGDP data
#   4. Set up Python virtual environment
#   5. QC SGDP data (UCSC liftOver hg19->hg38, match to KG, assign rsIDs)
#   6. Merge KG + HGDP + SGDP into a single fileset
#   7. Prepare and merge GIAB Ashkenazi parents (HG003, HG004)
#   8. Build merged metadata CSV with Neural ADMIXTURE ancestry labels
#   9. Build supervised reference population assignments
#  10. Install ADMIXTURE software
#  11. QC merged panel for ADMIXTURE (geno/MAF/LD/mind/kinship/HWE)
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
export PYTHON="${PROJECT_DIR}/tools/venv/bin/python"
export SUPERVISED_ADMIXTURE="${PROJECT_DIR}/supervised_admixture"
export ADMIXTURE="${TOOLS_BIN}/admixture"
export LIFTOVER="${TOOLS_BIN}/liftOver"
export CHAIN_FILE="${DOWNLOADS_DIR}/hg19ToHg38.over.chain.gz"

# ---------------------------------------------------------------------------
# Pipeline constants
# ---------------------------------------------------------------------------
# Post-merge genotype missingness
export GENO_THRESHOLD=0.03

# ADMIXTURE QC — genotype missingness
export GENO_ADMIXTURE=0.01

# ADMIXTURE QC — Hardy-Weinberg p-value threshold
export HWE_PVALUE="1e-50"

# ADMIXTURE QC — LD pruning parameters
export LD_WINDOW=50
export LD_STEP=10
export LD_R2=0.1

# ADMIXTURE QC — individual missingness
export MIND_ADMIXTURE=0.01

# ADMIXTURE QC — kinship (KING cutoffs)
export KING_CUTOFF_AMR=0.088          # AMR: ~3rd-degree relatives
export KING_CUTOFF_NONAMR=0.05        # Non-AMR: ~2nd-degree relatives

# ADMIXTURE run parameters
export N_FOLDS=3
export ADMIXTURE_SEED=42
export FLAG_THRESHOLD=0.95

# ---------------------------------------------------------------------------
# Configuration — interactive or via environment variables
# ---------------------------------------------------------------------------
# Non-interactive mode: set K_MODEL and MAF_ADMIXTURE as environment variables
# to skip the interactive prompts (useful for Docker / CI).
#   Example: K_MODEL=6 MAF_ADMIXTURE=0.0100 bash main.sh
# ---------------------------------------------------------------------------
if [[ -n "${K_MODEL:-}" && -n "${MAF_ADMIXTURE:-}" ]]; then
    echo "============================================"
    echo "Pipeline Configuration (from environment)"
    echo "============================================"
    if [[ "$K_MODEL" != "3" && "$K_MODEL" != "5" && "$K_MODEL" != "6" ]]; then
        echo "Error: K_MODEL must be 3, 5, or 6." >&2
        exit 1
    fi
    echo "  K_MODEL=${K_MODEL}"
    echo "  MAF_ADMIXTURE=${MAF_ADMIXTURE}"
    echo ""
else
    echo "============================================"
    echo "Pipeline Configuration"
    echo "============================================"
    echo ""
    echo "Select ADMIXTURE model:"
    echo ""
    echo "  K=3  African, American, European"
    echo "       Basic 3-way continental ancestry. Use for populations of"
    echo "       primarily African, American, and European descent."
    echo ""
    echo "  K=5  African, American, East Asian, European, South Asian"
    echo "       Adds East Asian and South Asian ancestry components."
    echo "       Use when samples may include or be admixed with"
    echo "       these ancestries."
    echo ""
    echo "  K=6  African, American, East Asian, European, Oceanian, South Asian"
    echo "       Full global ancestry including Oceanian / Pacific Islander."
    echo ""
    read -r -p "Enter K (3, 5, or 6) [default: 6]: " K_INPUT
    K_INPUT="${K_INPUT:-6}"
    if [[ "$K_INPUT" != "3" && "$K_INPUT" != "5" && "$K_INPUT" != "6" ]]; then
        echo "Error: K must be 3, 5, or 6." >&2
        exit 1
    fi
    export K_MODEL="${K_INPUT}"
    echo "  -> K=${K_MODEL}"
    echo ""

    echo "Select minor allele frequency (MAF) threshold for ADMIXTURE QC."
    echo "  Recommended range: 0.5% to 5%."
    echo "  Higher MAF = fewer, more common SNPs (faster, less noise)."
    echo "  Lower MAF  = more SNPs retained (more resolution, noisier)."
    echo ""
    read -r -p "Enter MAF as a percentage, e.g. 1 for 1% [default: 1]: " MAF_INPUT
    MAF_INPUT="${MAF_INPUT:-1}"
    if ! [[ "$MAF_INPUT" =~ ^[0-9]*\.?[0-9]+$ ]]; then
        echo "Error: MAF must be a number." >&2
        exit 1
    fi
    MAF_ADMIXTURE=$(awk "BEGIN {printf \"%.4f\", ${MAF_INPUT} / 100}")
    export MAF_ADMIXTURE
    echo "  -> MAF=${MAF_ADMIXTURE} (${MAF_INPUT}%)"
    echo ""
fi

# ---------------------------------------------------------------------------
# Volume mount support (Docker)
# ---------------------------------------------------------------------------
# When a volume is mounted at pipeline-output/ (the documented Docker usage),
# redirect all generated data directories there via symlinks so that results
# persist after the container exits.
#
#   docker run --rm -v $(pwd)/data:/app/pipeline-output public-statgen
#
PIPELINE_OUTPUT="${PROJECT_DIR}/pipeline-output"
if [[ -d "${PIPELINE_OUTPUT}" ]]; then
    echo "==> Volume detected at ${PIPELINE_OUTPUT} — redirecting output directories"
    for subdir in downloads qc merge supervised_admixture summary logs; do
        target="${PIPELINE_OUTPUT}/${subdir}"
        link="${PROJECT_DIR}/${subdir}"
        mkdir -p "${target}"
        if [[ ! -L "${link}" ]]; then
            rm -rf "${link}"
            ln -s "${target}" "${link}"
        fi
    done
    echo ""
fi

# ---------------------------------------------------------------------------
# Logging — all subsequent output goes to both terminal and log file
# ---------------------------------------------------------------------------
mkdir -p "${PROJECT_DIR}/logs"
LOG_FILE="${PROJECT_DIR}/logs/main_run_$(date +%Y%m%d_%H%M%S).log"
exec > >(tee -a "${LOG_FILE}") 2>&1

echo "============================================"
echo "Pipeline Run — $(date)"
echo "============================================"
echo "  K_MODEL=${K_MODEL}"
echo "  MAF_ADMIXTURE=${MAF_ADMIXTURE}"
echo "  GENO_ADMIXTURE=${GENO_ADMIXTURE}"
echo "  HWE_PVALUE=${HWE_PVALUE}"
echo "  LD_WINDOW=${LD_WINDOW}, LD_STEP=${LD_STEP}, LD_R2=${LD_R2}"
echo "  MIND_ADMIXTURE=${MIND_ADMIXTURE}"
echo "  KING_CUTOFF_AMR=${KING_CUTOFF_AMR}"
echo "  KING_CUTOFF_NONAMR=${KING_CUTOFF_NONAMR}"
echo "  N_FOLDS=${N_FOLDS}, ADMIXTURE_SEED=${ADMIXTURE_SEED}"
echo "  LOG_FILE=${LOG_FILE}"
echo ""

K="${K_MODEL}"

# ---------------------------------------------------------------------------
# Step 1 — Install PLINK and liftOver
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 1: Install PLINK and liftOver"
echo "============================================"
bash "${PROJECT_DIR}/setup_plink.sh"
bash "${PROJECT_DIR}/setup_liftover.sh"

# Verify the binaries are functional
"${PLINK1}" --version > /dev/null 2>&1 || { echo "Error: plink1 failed to run" >&2; exit 1; }
"${PLINK2}" --version > /dev/null 2>&1 || { echo "Error: plink2 failed to run" >&2; exit 1; }
[[ -x "${LIFTOVER}" ]] || { echo "Error: liftOver not found at ${LIFTOVER}" >&2; exit 1; }
echo ""

# ---------------------------------------------------------------------------
# Step 2 — Download reference panel data
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 2: Download KG, HGDP, SGDP, Neural ADMIXTURE, GIAB, and population data"
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
# Step 5 — QC SGDP (UCSC liftOver + rsID matching)
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 5: QC SGDP data"
echo "============================================"

if [[ -f "${QC_DIR}/sgdp_qc.bed" ]]; then
    echo "  [skip] SGDP QC output already exists in ${QC_DIR}/"
else
    mkdir -p "${QC_DIR}"

    # Unzip bim and copy bed/fam into QC directory
    python3 -c "
import zipfile, sys
with zipfile.ZipFile(sys.argv[1]) as z:
    names = [n for n in z.namelist() if not n.endswith('/')]
    sys.stdout.buffer.write(z.read(names[0]))
" "${DOWNLOADS_DIR}/sgdp_all.bim.zip" > "${QC_DIR}/sgdp_all.bim"
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
# Step 9 — Build supervised reference populations
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 9: Build supervised reference populations"
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
# Step 11 — QC merged panel for ADMIXTURE
# ---------------------------------------------------------------------------
echo "============================================"
echo "Step 11: QC merged panel for ADMIXTURE (MAF ${MAF_ADMIXTURE})"
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
echo "Step 12: Run ADMIXTURE supervised ancestry estimation (K=${K})"
echo "============================================"

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
echo "Step 13: Analyze ADMIXTURE results (K=${K})"
echo "============================================"

if [[ -f "${PROJECT_DIR}/summary/admixture-global-${K}/metadata_ancestry.csv" ]]; then
    echo "  [skip] ADMIXTURE analysis output already exists"
else
    "${PYTHON}" "${PROJECT_DIR}/analyze_admixture_results.py"
fi
echo ""

echo "============================================"
echo "Pipeline steps 1-13 complete."
echo "Log saved to ${LOG_FILE}"
echo "============================================"
