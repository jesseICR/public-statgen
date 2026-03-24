#!/usr/bin/env bash
#
# merge_kg_hgdp_sgdp.sh
# ----------------------
# Merges QC'd KG, HGDP, and SGDP bed/bim/fam filesets into a single fileset.
#
# Steps:
#   1. Align SGDP alleles to KG (fix A1/A2 flips)
#   2. Merge KG + HGDP (trial merge → exclude problem SNPs → clean merge)
#   3. Remove overlapping samples from SGDP (prefer KG/HGDP versions)
#   4. Three-way merge: (KG+HGDP) + unique SGDP
#   5. Post-merge QC: remove ambiguous SNPs (A/T, C/G) and apply geno filter
#   6. Clean up intermediates
#
# Expected environment:
#   QC_DIR         — directory containing kg_qc, hgdp_qc, sgdp_qc filesets
#   MERGE_DIR      — output directory for merged fileset (must exist)
#   PLINK1         — path to plink1 binary
#   PLINK2         — path to plink2 binary
#   PLINK_MEMORY   — memory limit in MB for PLINK
#   PLINK_THREADS  — number of threads for PLINK
#   GENO_THRESHOLD — max per-variant missingness (e.g. 0.03)
#
# Outputs (in MERGE_DIR):
#   merged_kg_hgdp_sgdp.{bed,bim,fam}
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
for var in QC_DIR MERGE_DIR PLINK1 PLINK2 PLINK_MEMORY PLINK_THREADS GENO_THRESHOLD; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: ${var} is not set." >&2
        exit 1
    fi
done

PLINK_FLAGS=(--memory "${PLINK_MEMORY}" --threads "${PLINK_THREADS}")

# ---------------------------------------------------------------------------
# 1. Align SGDP alleles to KG reference
# ---------------------------------------------------------------------------
echo "  1/5 Aligning SGDP alleles to KG ..."

# Build rsID → A1 reference from KG
awk '{print $2, $5}' "${QC_DIR}/kg_qc.bim" > "${MERGE_DIR}/kg_a1_ref.txt"

"${PLINK1}" --bfile "${QC_DIR}/sgdp_qc" \
    --a1-allele "${MERGE_DIR}/kg_a1_ref.txt" 2 1 \
    --make-bed \
    --out "${MERGE_DIR}/sgdp_aligned" \
    "${PLINK_FLAGS[@]}"

# ---------------------------------------------------------------------------
# 2. Merge KG + HGDP
# ---------------------------------------------------------------------------
echo "  2/5 Merging KG + HGDP ..."

# Find common SNPs
awk '{print $2}' "${QC_DIR}/kg_qc.bim"   | sort > "${MERGE_DIR}/kg_snps.txt"
awk '{print $2}' "${QC_DIR}/hgdp_qc.bim" | sort > "${MERGE_DIR}/hgdp_snps.txt"
comm -12 "${MERGE_DIR}/kg_snps.txt" "${MERGE_DIR}/hgdp_snps.txt" > "${MERGE_DIR}/common_kg_hgdp.txt"

# Trial merge — expected to fail on a few problem SNPs
"${PLINK1}" --bfile "${QC_DIR}/kg_qc" \
    --extract "${MERGE_DIR}/common_kg_hgdp.txt" \
    --bmerge "${QC_DIR}/hgdp_qc.bed" "${QC_DIR}/hgdp_qc.bim" "${QC_DIR}/hgdp_qc.fam" \
    --make-bed \
    --out "${MERGE_DIR}/temp_merge1" \
    "${PLINK_FLAGS[@]}" || true

# Exclude problem SNPs and re-extract from both datasets
"${PLINK1}" --bfile "${QC_DIR}/kg_qc" \
    --extract "${MERGE_DIR}/common_kg_hgdp.txt" \
    --exclude "${MERGE_DIR}/temp_merge1-merge.missnp" \
    --make-bed \
    --out "${MERGE_DIR}/kg_qc_clean" \
    "${PLINK_FLAGS[@]}"

"${PLINK1}" --bfile "${QC_DIR}/hgdp_qc" \
    --extract "${MERGE_DIR}/common_kg_hgdp.txt" \
    --exclude "${MERGE_DIR}/temp_merge1-merge.missnp" \
    --make-bed \
    --out "${MERGE_DIR}/hgdp_qc_clean" \
    "${PLINK_FLAGS[@]}"

# Clean merge
"${PLINK1}" --bfile "${MERGE_DIR}/kg_qc_clean" \
    --bmerge "${MERGE_DIR}/hgdp_qc_clean.bed" "${MERGE_DIR}/hgdp_qc_clean.bim" "${MERGE_DIR}/hgdp_qc_clean.fam" \
    --make-bed \
    --out "${MERGE_DIR}/merged_kg_hgdp" \
    "${PLINK_FLAGS[@]}"

# ---------------------------------------------------------------------------
# 3. Remove overlapping samples from SGDP
# ---------------------------------------------------------------------------
echo "  3/5 Removing overlapping SGDP samples ..."

awk '{print $2}' "${MERGE_DIR}/merged_kg_hgdp.fam" | sort -u > "${MERGE_DIR}/kg_hgdp_iids.txt"

# Build remove list: SGDP samples whose IID already exists in KG/HGDP
awk 'NR==FNR {ids[$1]=1; next} $2 in ids {print $1, $2}' \
    "${MERGE_DIR}/kg_hgdp_iids.txt" "${MERGE_DIR}/sgdp_aligned.fam" \
    > "${MERGE_DIR}/sgdp_remove.txt"

echo "    Overlapping samples to remove: $(wc -l < "${MERGE_DIR}/sgdp_remove.txt")"

"${PLINK1}" --bfile "${MERGE_DIR}/sgdp_aligned" \
    --remove "${MERGE_DIR}/sgdp_remove.txt" \
    --make-bed \
    --out "${MERGE_DIR}/sgdp_unique" \
    "${PLINK_FLAGS[@]}"

# ---------------------------------------------------------------------------
# 4. Final three-way merge
# ---------------------------------------------------------------------------
echo "  4/5 Three-way merge: (KG+HGDP) + SGDP ..."

"${PLINK1}" --bfile "${MERGE_DIR}/merged_kg_hgdp" \
    --bmerge "${MERGE_DIR}/sgdp_unique.bed" "${MERGE_DIR}/sgdp_unique.bim" "${MERGE_DIR}/sgdp_unique.fam" \
    --make-bed \
    --allow-no-sex \
    --out "${MERGE_DIR}/merged_kg_hgdp_sgdp" \
    "${PLINK_FLAGS[@]}"

# ---------------------------------------------------------------------------
# 5. Post-merge QC: remove ambiguous SNPs and apply geno filter
# ---------------------------------------------------------------------------
echo "  5/5 Post-merge QC (ambiguous SNPs + geno filter) ..."

# Identify ambiguous (A/T, C/G) SNPs
awk '($5 == "A" && $6 == "T") || ($5 == "T" && $6 == "A") || \
     ($5 == "C" && $6 == "G") || ($5 == "G" && $6 == "C") \
     {print $2}' "${MERGE_DIR}/merged_kg_hgdp_sgdp.bim" \
     > "${MERGE_DIR}/ambiguous_snps.txt"

echo "    Ambiguous SNPs to remove: $(wc -l < "${MERGE_DIR}/ambiguous_snps.txt")"

# Apply both filters in a single PLINK2 pass
"${PLINK2}" --bfile "${MERGE_DIR}/merged_kg_hgdp_sgdp" \
    --exclude "${MERGE_DIR}/ambiguous_snps.txt" \
    --geno "${GENO_THRESHOLD}" \
    --no-input-missing-phenotype \
    --make-bed \
    --out "${MERGE_DIR}/merged_kg_hgdp_sgdp_qc" \
    "${PLINK_FLAGS[@]}"

# Replace raw merge with QC'd version
rm -f "${MERGE_DIR}"/merged_kg_hgdp_sgdp.{bed,bim,fam}
for ext in bed bim fam; do
    mv "${MERGE_DIR}/merged_kg_hgdp_sgdp_qc.${ext}" "${MERGE_DIR}/merged_kg_hgdp_sgdp.${ext}"
done

echo "    Samples: $(wc -l < "${MERGE_DIR}/merged_kg_hgdp_sgdp.fam")"
echo "    Variants: $(wc -l < "${MERGE_DIR}/merged_kg_hgdp_sgdp.bim")"

# ---------------------------------------------------------------------------
# 6. Clean up intermediates
# ---------------------------------------------------------------------------
echo "  Cleaning up intermediate files ..."
rm -f "${MERGE_DIR}"/temp_merge1* "${MERGE_DIR}"/*_clean* "${MERGE_DIR}"/*-merge.*
rm -f "${MERGE_DIR}"/sgdp_aligned* "${MERGE_DIR}"/sgdp_unique*
rm -f "${MERGE_DIR}"/merged_kg_hgdp.{bed,bim,fam,log,nosex}
rm -f "${MERGE_DIR}"/merged_kg_hgdp_sgdp_qc.*
rm -f "${MERGE_DIR}"/{kg_a1_ref,kg_snps,hgdp_snps,common_kg_hgdp,kg_hgdp_iids,sgdp_remove,ambiguous_snps}.txt
rm -f "${MERGE_DIR}"/*.nosex "${MERGE_DIR}"/*.log

echo "  Merge complete."
