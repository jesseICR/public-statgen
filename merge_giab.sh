#!/usr/bin/env bash
#
# merge_giab.sh
# -------------
# Merges GIAB Ashkenazi parent samples (HG003, HG004) into the reference panel.
#
# Steps:
#   1. Merge GIAB bed with reference panel using plink1
#   2. Handle any allele-mismatch SNPs (exclude + retry)
#   3. Normalize fam file (FID=0, PAT=0, MAT=0) for all samples
#   4. Report SNP and sample counts
#   5. Clean up intermediates
#
# Expected environment:
#   MERGE_DIR      — directory containing merged_kg_hgdp_sgdp.{bed,bim,fam}
#   QC_DIR         — directory containing giab_qc.{bed,bim,fam}
#   PLINK1         — path to plink1 binary
#   PLINK_MEMORY   — memory limit in MB
#   PLINK_THREADS  — number of threads
#
# Outputs (in MERGE_DIR, replacing existing merged files):
#   merged_kg_hgdp_sgdp.{bed,bim,fam}   — now includes GIAB samples
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
for var in MERGE_DIR QC_DIR PLINK1 PLINK_MEMORY PLINK_THREADS; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: ${var} is not set." >&2
        exit 1
    fi
done

PLINK_FLAGS=(--memory "${PLINK_MEMORY}" --threads "${PLINK_THREADS}")

REF="${MERGE_DIR}/merged_kg_hgdp_sgdp"
GIAB="${QC_DIR}/giab_qc"
OUT="${MERGE_DIR}/merged_with_giab"

SNPS_BEFORE=$(wc -l < "${REF}.bim")
SAMPLES_BEFORE=$(wc -l < "${REF}.fam")
echo "  Reference panel: ${SNPS_BEFORE} SNPs, ${SAMPLES_BEFORE} samples"

GIAB_SNPS=$(wc -l < "${GIAB}.bim")
GIAB_SAMPLES=$(wc -l < "${GIAB}.fam")
echo "  GIAB data: ${GIAB_SNPS} SNPs, ${GIAB_SAMPLES} samples"

# ---------------------------------------------------------------------------
# 1. Merge GIAB into reference panel
# ---------------------------------------------------------------------------
echo "  Merging GIAB into reference panel ..."

"${PLINK1}" --bfile "${REF}" \
    --bmerge "${GIAB}.bed" "${GIAB}.bim" "${GIAB}.fam" \
    --make-bed --allow-no-sex \
    --out "${OUT}" \
    "${PLINK_FLAGS[@]}" || true

# Handle merge errors (allele mismatches)
if [[ -f "${OUT}-merge.missnp" ]]; then
    N_MISSNP=$(wc -l < "${OUT}-merge.missnp")
    echo "  Excluding ${N_MISSNP} problem SNPs and retrying ..."

    "${PLINK1}" --bfile "${REF}" \
        --exclude "${OUT}-merge.missnp" \
        --make-bed \
        --out "${MERGE_DIR}/ref_clean" \
        "${PLINK_FLAGS[@]}"

    "${PLINK1}" --bfile "${GIAB}" \
        --exclude "${OUT}-merge.missnp" \
        --make-bed \
        --out "${MERGE_DIR}/giab_clean" \
        "${PLINK_FLAGS[@]}"

    "${PLINK1}" --bfile "${MERGE_DIR}/ref_clean" \
        --bmerge "${MERGE_DIR}/giab_clean.bed" "${MERGE_DIR}/giab_clean.bim" "${MERGE_DIR}/giab_clean.fam" \
        --make-bed --allow-no-sex \
        --out "${OUT}" \
        "${PLINK_FLAGS[@]}"

    rm -f "${MERGE_DIR}"/ref_clean.* "${MERGE_DIR}"/giab_clean.*
fi

# ---------------------------------------------------------------------------
# 2. Replace original with merged
# ---------------------------------------------------------------------------
rm -f "${REF}".{bed,bim,fam}
for ext in bed bim fam; do
    mv "${OUT}.${ext}" "${REF}.${ext}"
done

# ---------------------------------------------------------------------------
# 3. Normalize fam file (FID=0, PAT=0, MAT=0 for all samples)
# ---------------------------------------------------------------------------
echo "  Normalizing fam file ..."
awk 'BEGIN{OFS="\t"} {$1=0; $3=0; $4=0; $6=-9; print}' "${REF}.fam" > "${REF}.fam.tmp"
mv "${REF}.fam.tmp" "${REF}.fam"

# ---------------------------------------------------------------------------
# 4. Report
# ---------------------------------------------------------------------------
SNPS_AFTER=$(wc -l < "${REF}.bim")
SAMPLES_AFTER=$(wc -l < "${REF}.fam")
SNPS_LOST=$((SNPS_BEFORE - SNPS_AFTER))

echo "  Merged panel: ${SNPS_AFTER} SNPs, ${SAMPLES_AFTER} samples"
echo "  SNPs lost from reference panel: ${SNPS_LOST}"
echo "  Samples added: $((SAMPLES_AFTER - SAMPLES_BEFORE))"

# ---------------------------------------------------------------------------
# 5. Clean up
# ---------------------------------------------------------------------------
rm -f "${OUT}"* "${MERGE_DIR}"/*-merge.* "${MERGE_DIR}"/*.nosex "${MERGE_DIR}"/*.log
rm -f "${QC_DIR}"/giab_qc.*

echo "  GIAB merge complete."
