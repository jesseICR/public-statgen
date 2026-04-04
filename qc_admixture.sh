#!/usr/bin/env bash
#
# qc_admixture.sh
# ----------------
# Quality-control pipeline preparing the merged reference panel for
# supervised ancestry estimation with ADMIXTURE.
#
# SNP-based filters (geno, MAF, LD) are computed ONLY on supervised
# samples — the individuals that will actually train the model. The resulting
# SNP list is then applied to all samples so non-supervised individuals
# remain in the panel for projection.
#
# Hardy-Weinberg filtering is performed AFTER kinship pruning so that it
# runs only on unrelated supervised individuals, avoiding HWE inflation
# from cryptic relatedness.
#
# QC order:
#   0. Extract supervised samples for QC
#   1. Genotype missingness filter (--geno)                  [supervised only]
#   2. Minor allele frequency filter (--maf)                 [supervised only]
#   3. Exclude long-range LD regions                          [supervised only]
#   4. LD pruning (--indep-pairwise)                          [supervised only]
#   5. Individual missingness filter (--mind)                 [supervised only]
#   6. Kinship: remove related individuals                    [supervised only]
#        AMR samples:     KING cutoff (env KING_CUTOFF_AMR)
#        Non-AMR samples: KING cutoff (env KING_CUTOFF_NONAMR)
#   7. Hardy-Weinberg equilibrium exact test (--hwe)          [unrelated supervised]
#   8. Apply final supervised SNP+sample list to all samples for projection
#
# Expected environment:
#   MERGE_DIR              — directory containing merged_kg_hgdp_sgdp.{bed,bim,fam}
#   PLINK1                 — path to plink1 binary
#   PLINK2                 — path to plink2 binary
#   PYTHON                 — path to venv python
#   PLINK_MEMORY           — memory limit in MB
#   PLINK_THREADS          — number of threads
#   SUPERVISED_ADMIXTURE   — output directory for all ADMIXTURE work
#   GENO_ADMIXTURE         — genotype missingness threshold (e.g. 0.01)
#   MAF_ADMIXTURE          — minor allele frequency threshold (e.g. 0.02)
#   HWE_PVALUE             — Hardy-Weinberg p-value threshold (e.g. 1e-50)
#   LD_WINDOW              — LD pruning window size in SNPs (e.g. 50)
#   LD_STEP                — LD pruning step size (e.g. 10)
#   LD_R2                  — LD pruning r-squared threshold (e.g. 0.1)
#   MIND_ADMIXTURE         — individual missingness threshold (e.g. 0.01)
#   KING_CUTOFF_AMR        — kinship cutoff for AMR samples (e.g. 0.088)
#   KING_CUTOFF_NONAMR     — kinship cutoff for non-AMR samples (e.g. 0.05)
#
# Input files (read-only):
#   summary/supervised.csv  — sample_id,population_id,reference_population
#   summary/metadata.csv    — sample_id,population_id,superpopulation,...
#
# Outputs (in SUPERVISED_ADMIXTURE):
#   ancestry_qc.{bed,bim,fam}    — final QC'd fileset (all samples, supervised SNPs)
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
for var in MERGE_DIR PLINK1 PLINK2 PYTHON PLINK_MEMORY PLINK_THREADS SUPERVISED_ADMIXTURE \
           GENO_ADMIXTURE MAF_ADMIXTURE HWE_PVALUE LD_WINDOW LD_STEP LD_R2 \
           MIND_ADMIXTURE KING_CUTOFF_AMR KING_CUTOFF_NONAMR; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: ${var} is not set." >&2
        exit 1
    fi
done

PROJECT_DIR="$(cd "$(dirname "$0")" && pwd)"
SCRAP="${SUPERVISED_ADMIXTURE}/scrap"
PLINK_FLAGS=(--memory "${PLINK_MEMORY}" --threads "${PLINK_THREADS}" --no-input-missing-phenotype)

mkdir -p "${SUPERVISED_ADMIXTURE}" "${SCRAP}"

# ---------------------------------------------------------------------------
# Helper: format a number with comma thousands separators
# ---------------------------------------------------------------------------
fmt() { printf "%'d" "$1"; }

# ---------------------------------------------------------------------------
# Helper: count SNPs and samples in a PLINK bed fileset
# ---------------------------------------------------------------------------
count_snps()   { wc -l < "$1.bim"; }
count_samples(){ wc -l < "$1.fam"; }

# ---------------------------------------------------------------------------
# Starting point
# ---------------------------------------------------------------------------
INPUT="${MERGE_DIR}/merged_kg_hgdp_sgdp"
SNPS_BEFORE=$(count_snps "${INPUT}")
SAMPLES_BEFORE=$(count_samples "${INPUT}")
echo "  Starting data: $(fmt ${SNPS_BEFORE}) SNPs, $(fmt ${SAMPLES_BEFORE}) samples"
echo ""

# ===================================================================
# STEP 0 — Extract supervised samples for SNP filtering
# ===================================================================
echo "  Step 0: Extract supervised samples for SNP-based QC"

"${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
fam = pd.read_csv('${INPUT}.fam', sep='\s+', header=None, names=['FID','IID','PAT','MAT','SEX','PHENO'])
sup_ids = set(sup['sample_id'])
keep = fam[fam['IID'].isin(sup_ids)][['FID','IID']]
keep.to_csv('${SCRAP}/supervised_keep.txt', sep='\t', header=False, index=False)
print(f'    Supervised samples in merged panel: {len(keep):,} of {len(sup):,}')
"

"${PLINK2}" --bfile "${INPUT}" \
    --keep "${SCRAP}/supervised_keep.txt" \
    --make-bed \
    --out "${SCRAP}/supervised_only" \
    "${PLINK_FLAGS[@]}"

SUP_SAMPLES=$(count_samples "${SCRAP}/supervised_only")
echo "    Extracted $(fmt ${SUP_SAMPLES}) supervised samples for SNP filtering"
echo ""

# ===================================================================
# STEP 1 — Genotype missingness (--geno) [supervised only]
# ===================================================================
echo "  Step 1/8: Genotype missingness filter (--geno ${GENO_ADMIXTURE}) [supervised only]"

"${PLINK2}" --bfile "${SCRAP}/supervised_only" \
    --geno "${GENO_ADMIXTURE}" \
    --make-bed \
    --out "${SCRAP}/geno_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/geno_filtered")
SNPS_REMOVED=$((SNPS_BEFORE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs with >${GENO_ADMIXTURE} missingness in supervised samples"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 2 — Minor allele frequency (--maf) [supervised only]
# ===================================================================
echo "  Step 2/8: MAF filter (--maf ${MAF_ADMIXTURE}) [supervised only]"

SNPS_BEFORE_MAF=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/geno_filtered" \
    --maf "${MAF_ADMIXTURE}" \
    --make-bed \
    --out "${SCRAP}/maf_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/maf_filtered")
SNPS_REMOVED=$((SNPS_BEFORE_MAF - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs with MAF < ${MAF_ADMIXTURE} in supervised samples"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 3 — Exclude long-range LD regions [supervised only]
# ===================================================================
echo "  Step 3/8: Exclude long-range LD regions (Price et al. 2008, hg38) [supervised only]"

HIGH_LD_BED="${DOWNLOADS_DIR}/high_ld_regions_hg38.bed"
if [[ ! -f "${HIGH_LD_BED}" ]]; then
    echo "Error: high-LD regions file not found: ${HIGH_LD_BED}" >&2
    echo "Run download_files.sh first." >&2
    exit 1
fi

# Convert BED (0-based) to PLINK range format (1-based)
awk 'BEGIN{OFS="\t"} {
    chr = $1; sub(/^chr/, "", chr)
    print chr, $2+1, $3, $4
}' "${HIGH_LD_BED}" > "${SCRAP}/high_ld_ranges.txt"

SNPS_BEFORE_LD=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/maf_filtered" \
    --exclude range "${SCRAP}/high_ld_ranges.txt" \
    --make-bed \
    --out "${SCRAP}/ldregion_excluded" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/ldregion_excluded")
SNPS_REMOVED=$((SNPS_BEFORE_LD - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs in long-range LD regions"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 4 — LD pruning (--indep-pairwise) [supervised only]
# ===================================================================
echo "  Step 4/8: LD pruning (window=${LD_WINDOW}, step=${LD_STEP}, r2=${LD_R2}) [supervised only]"

"${PLINK2}" --bfile "${SCRAP}/ldregion_excluded" \
    --indep-pairwise "${LD_WINDOW}" "${LD_STEP}" "${LD_R2}" \
    --out "${SCRAP}/ld_prune" \
    "${PLINK_FLAGS[@]}"

SNPS_BEFORE_PRUNE=${SNPS_AFTER}

"${PLINK2}" --bfile "${SCRAP}/ldregion_excluded" \
    --extract "${SCRAP}/ld_prune.prune.in" \
    --make-bed \
    --out "${SCRAP}/ld_pruned" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/ld_pruned")
SNPS_REMOVED=$((SNPS_BEFORE_PRUNE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs by LD pruning"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 5 — Individual missingness (--mind) [supervised only]
# ===================================================================
echo "  Step 5/8: Individual missingness filter (--mind ${MIND_ADMIXTURE}) [supervised only]"

SAMPLES_BEFORE_MIND=$(count_samples "${SCRAP}/ld_pruned")

"${PLINK2}" --bfile "${SCRAP}/ld_pruned" \
    --mind "${MIND_ADMIXTURE}" \
    --make-bed \
    --out "${SCRAP}/mind_filtered" \
    "${PLINK_FLAGS[@]}"

SAMPLES_AFTER_MIND=$(count_samples "${SCRAP}/mind_filtered")
SAMPLES_REMOVED_MIND=$((SAMPLES_BEFORE_MIND - SAMPLES_AFTER_MIND))

if [[ ${SAMPLES_REMOVED_MIND} -gt 0 ]]; then
    echo "    Removed $(fmt ${SAMPLES_REMOVED_MIND}) supervised individuals with >${MIND_ADMIXTURE} missingness"

    awk '{print $2}' "${SCRAP}/ld_pruned.fam" | sort > "${SCRAP}/before_mind_iids.txt"
    awk '{print $2}' "${SCRAP}/mind_filtered.fam" | sort > "${SCRAP}/after_mind_iids.txt"
    comm -23 "${SCRAP}/before_mind_iids.txt" "${SCRAP}/after_mind_iids.txt" > "${SCRAP}/removed_mind_iids.txt"

    "${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
removed = set(open('${SCRAP}/removed_mind_iids.txt').read().split())
sup_removed = sup[sup['sample_id'].isin(removed)]
print('    Removed supervised individuals:')
for _, row in sup_removed.iterrows():
    print(f'      {row[\"sample_id\"]} — {row[\"population_id\"]} ({row[\"reference_population\"]})')
"
else
    echo "    No supervised individuals removed"
fi
echo "    Remaining: $(fmt ${SAMPLES_AFTER_MIND}) supervised individuals"
echo ""

# ===================================================================
# STEP 6 — Kinship (separate AMR and non-AMR) [supervised only]
# ===================================================================
echo "  Step 6/8: Kinship filtering [supervised only] (AMR=${KING_CUTOFF_AMR}, non-AMR=${KING_CUTOFF_NONAMR})"

"${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
fam = pd.read_csv('${SCRAP}/mind_filtered.fam', sep='\s+', header=None, names=['FID','IID','PAT','MAT','SEX','PHENO'])
sup_ref = dict(zip(sup['sample_id'], sup['reference_population']))
amr_ids = set(sup.loc[sup['reference_population'] == 'American', 'sample_id'])
amr_in_fam = fam[fam['IID'].isin(amr_ids)]
nonamr_in_fam = fam[~fam['IID'].isin(amr_ids)]
amr_in_fam[['FID', 'IID']].to_csv('${SCRAP}/amr_samples.txt', sep='\t', header=False, index=False)
nonamr_in_fam[['FID', 'IID']].to_csv('${SCRAP}/nonamr_samples.txt', sep='\t', header=False, index=False)
print(f'    Supervised AMR: {len(amr_in_fam):,}')
print(f'    Supervised non-AMR: {len(nonamr_in_fam):,}')
"

# --- 6a. Kinship on AMR ---
echo ""
echo "    6a. AMR kinship (KING cutoff = ${KING_CUTOFF_AMR})"

AMR_COUNT_BEFORE=$(wc -l < "${SCRAP}/amr_samples.txt")

"${PLINK2}" --bfile "${SCRAP}/mind_filtered" \
    --keep "${SCRAP}/amr_samples.txt" \
    --king-cutoff "${KING_CUTOFF_AMR}" \
    --out "${SCRAP}/amr_king" \
    "${PLINK_FLAGS[@]}"

if [[ -f "${SCRAP}/amr_king.king.cutoff.out.id" ]]; then
    AMR_REMOVED=$(wc -l < "${SCRAP}/amr_king.king.cutoff.out.id")
    if head -1 "${SCRAP}/amr_king.king.cutoff.out.id" | grep -q "^#"; then
        AMR_REMOVED=$((AMR_REMOVED - 1))
    fi
else
    AMR_REMOVED=0
fi

echo "        AMR removed (related): $(fmt ${AMR_REMOVED}) of $(fmt ${AMR_COUNT_BEFORE})"

if [[ ${AMR_REMOVED} -gt 0 ]]; then
    "${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
removed = pd.read_csv('${SCRAP}/amr_king.king.cutoff.out.id', sep='\t')
removed_ids = set(removed.iloc[:, 1] if removed.columns[0] == '#FID' else removed.iloc[:, 0])
sup_removed = sup[sup['sample_id'].isin(removed_ids)]
print('        Removed supervised AMR individuals:')
for _, row in sup_removed.iterrows():
    print(f'          {row[\"sample_id\"]} — {row[\"population_id\"]} ({row[\"reference_population\"]})')
"
fi

# --- 6b. Kinship on non-AMR ---
echo ""
echo "    6b. Non-AMR kinship (KING cutoff = ${KING_CUTOFF_NONAMR})"

NONAMR_COUNT_BEFORE=$(wc -l < "${SCRAP}/nonamr_samples.txt")

"${PLINK2}" --bfile "${SCRAP}/mind_filtered" \
    --keep "${SCRAP}/nonamr_samples.txt" \
    --king-cutoff "${KING_CUTOFF_NONAMR}" \
    --out "${SCRAP}/nonamr_king" \
    "${PLINK_FLAGS[@]}"

if [[ -f "${SCRAP}/nonamr_king.king.cutoff.out.id" ]]; then
    NONAMR_REMOVED=$(wc -l < "${SCRAP}/nonamr_king.king.cutoff.out.id")
    if head -1 "${SCRAP}/nonamr_king.king.cutoff.out.id" | grep -q "^#"; then
        NONAMR_REMOVED=$((NONAMR_REMOVED - 1))
    fi
else
    NONAMR_REMOVED=0
fi

echo "        Non-AMR removed (related): $(fmt ${NONAMR_REMOVED}) of $(fmt ${NONAMR_COUNT_BEFORE})"

if [[ ${NONAMR_REMOVED} -gt 0 ]]; then
    "${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
removed = pd.read_csv('${SCRAP}/nonamr_king.king.cutoff.out.id', sep='\t')
removed_ids = set(removed.iloc[:, 1] if removed.columns[0] == '#FID' else removed.iloc[:, 0])
sup_removed = sup[sup['sample_id'].isin(removed_ids)]
if len(sup_removed) > 0:
    print('        Removed supervised non-AMR individuals:')
    for _, row in sup_removed.iterrows():
        print(f'          {row[\"sample_id\"]} — {row[\"population_id\"]} ({row[\"reference_population\"]})')
"
fi

# --- Save supervised QC'd sample list ---
echo ""
tail -n +2 "${SCRAP}/amr_king.king.cutoff.in.id" > "${SCRAP}/supervised_qc_keep.txt"
tail -n +2 "${SCRAP}/nonamr_king.king.cutoff.in.id" >> "${SCRAP}/supervised_qc_keep.txt"

TOTAL_REMOVED=$((AMR_REMOVED + NONAMR_REMOVED))
SUP_FINAL=$(wc -l < "${SCRAP}/supervised_qc_keep.txt")
echo "    Total supervised removed by kinship: $(fmt ${TOTAL_REMOVED})"
echo "    Supervised samples after kinship: $(fmt ${SUP_FINAL})"
echo ""

# ===================================================================
# STEP 7 — HWE on unrelated supervised samples
# ===================================================================
echo "  Step 7/8: HWE exact test (--hwe ${HWE_PVALUE}) [unrelated supervised only]"

# Extract unrelated supervised from LD-pruned fileset
"${PLINK2}" --bfile "${SCRAP}/ld_pruned" \
    --keep "${SCRAP}/supervised_qc_keep.txt" \
    --make-bed \
    --out "${SCRAP}/unrelated_supervised" \
    "${PLINK_FLAGS[@]}"

SNPS_BEFORE_HWE=$(count_snps "${SCRAP}/unrelated_supervised")

"${PLINK2}" --bfile "${SCRAP}/unrelated_supervised" \
    --hwe "${HWE_PVALUE}" \
    --make-bed \
    --out "${SCRAP}/hwe_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER_HWE=$(count_snps "${SCRAP}/hwe_filtered")
SNPS_REMOVED_HWE=$((SNPS_BEFORE_HWE - SNPS_AFTER_HWE))
echo "    Removed $(fmt ${SNPS_REMOVED_HWE}) SNPs failing HWE (p < ${HWE_PVALUE}) in unrelated supervised"
echo "    Remaining: $(fmt ${SNPS_AFTER_HWE}) SNPs"
echo ""

# ===================================================================
# STEP 8 — Apply final SNPs + add non-supervised back for projection
# ===================================================================
echo "  Step 8/8: Build final panel (supervised QC'd SNPs, all samples)"

# Get the SNP list from the HWE-filtered unrelated supervised data
awk '{print $2}' "${SCRAP}/hwe_filtered.bim" > "${SCRAP}/final_snps.txt"
SNPS_FINAL=$(wc -l < "${SCRAP}/final_snps.txt")

# Apply SNP list to the full merged panel (brings non-supervised back)
# then remove only the supervised individuals who failed QC
"${PYTHON}" -c "
import pandas as pd
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
sup_ids = set(sup['sample_id'])
keep = set(open('${SCRAP}/supervised_qc_keep.txt').read().split())
# Supervised samples that were removed by QC
sup_removed = sup_ids - keep
fam = pd.read_csv('${INPUT}.fam', sep='\s+', header=None, names=['FID','IID','PAT','MAT','SEX','PHENO'])
# Keep: all non-supervised + supervised that passed QC
keep_fam = fam[~fam['IID'].isin(sup_removed)]
keep_fam[['FID','IID']].to_csv('${SCRAP}/final_keep.txt', sep='\t', header=False, index=False)
n_nonsup = (~keep_fam['IID'].isin(sup_ids)).sum()
n_sup = keep_fam['IID'].isin(sup_ids).sum()
print(f'    Supervised (QC passed): {n_sup:,}')
print(f'    Non-supervised (all): {n_nonsup:,}')
print(f'    Total: {len(keep_fam):,}')
"

"${PLINK2}" --bfile "${INPUT}" \
    --extract "${SCRAP}/final_snps.txt" \
    --keep "${SCRAP}/final_keep.txt" \
    --make-bed \
    --out "${SUPERVISED_ADMIXTURE}/ancestry_qc" \
    "${PLINK_FLAGS[@]}"

SAMPLES_FINAL=$(count_samples "${SUPERVISED_ADMIXTURE}/ancestry_qc")
SNPS_FINAL_OUT=$(count_snps "${SUPERVISED_ADMIXTURE}/ancestry_qc")
echo "    Output: $(fmt ${SAMPLES_FINAL}) samples, $(fmt ${SNPS_FINAL_OUT}) SNPs"
echo ""

# ===================================================================
# Summary
# ===================================================================
echo "  ================================================"
echo "  QC Summary"
echo "  ================================================"
echo "  All filters computed on supervised samples only (n=$(fmt ${SUP_SAMPLES}))"
echo "  Parameters:"
echo "    Geno=${GENO_ADMIXTURE}, MAF=${MAF_ADMIXTURE}, HWE=${HWE_PVALUE}"
echo "    LD window=${LD_WINDOW} step=${LD_STEP} r2=${LD_R2}"
echo "    Mind=${MIND_ADMIXTURE}, KING AMR=${KING_CUTOFF_AMR}, KING non-AMR=${KING_CUTOFF_NONAMR}"
echo "  SNPs:    $(fmt ${SNPS_BEFORE}) -> $(fmt ${SNPS_FINAL_OUT})"
echo "  Supervised samples: $(fmt ${SUP_SAMPLES}) -> $(fmt ${SUP_FINAL})"
echo "  Final panel: $(fmt ${SAMPLES_FINAL}) samples ($(fmt ${SUP_FINAL}) supervised + non-supervised)"
echo "  Output:  ${SUPERVISED_ADMIXTURE}/ancestry_qc.{bed,bim,fam}"
echo ""

# Clean up logs
rm -f "${SCRAP}"/*.log "${SCRAP}"/*.nosex
rm -f "${SCRAP}"/before_mind_iids.txt "${SCRAP}"/after_mind_iids.txt

echo "  ADMIXTURE QC complete."
