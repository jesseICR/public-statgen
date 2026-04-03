#!/usr/bin/env bash
#
# qc_admixture.sh
# ----------------
# Quality-control pipeline preparing the merged reference panel for
# supervised ancestry estimation with ADMIXTURE.
#
# SNP-based filters (geno, MAF, HWE, LD) are computed ONLY on supervised
# samples — the individuals that will actually train the model. The resulting
# SNP list is then applied to all samples so non-supervised individuals
# remain in the panel for projection.
#
# QC order:
#   0. Extract supervised samples for SNP filtering
#   1. Genotype missingness filter (--geno 0.01)          [supervised only]
#   2. Minor allele frequency filter (--maf 0.03)         [supervised only]
#   3. Hardy-Weinberg equilibrium exact test (--hwe 1e-50) [supervised only]
#   4. Exclude long-range LD regions                       [supervised only]
#   5. LD pruning (--indep-pairwise 50 10 0.1)            [supervised only]
#   6. Apply SNP list to ALL samples
#   7. Individual missingness filter (--mind 0.01)         [all samples]
#   8. Kinship: remove related individuals                 [all samples]
#        AMR samples:     KING cutoff 0.088  (approx 3rd-degree)
#        Non-AMR samples: KING cutoff 0.05   (approx 2nd-degree)
#
# Expected environment:
#   MERGE_DIR              — directory containing merged_kg_hgdp_sgdp.{bed,bim,fam}
#   PLINK1                 — path to plink1 binary
#   PLINK2                 — path to plink2 binary
#   PYTHON                 — path to venv python
#   PLINK_MEMORY           — memory limit in MB
#   PLINK_THREADS          — number of threads
#   SUPERVISED_ADMIXTURE   — output directory for all ADMIXTURE work
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
for var in MERGE_DIR PLINK1 PLINK2 PYTHON PLINK_MEMORY PLINK_THREADS SUPERVISED_ADMIXTURE; do
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
# STEP 1 — Genotype missingness (--geno 0.01) [supervised only]
# ===================================================================
echo "  Step 1/8: Genotype missingness filter (--geno 0.01) [supervised only]"

"${PLINK2}" --bfile "${SCRAP}/supervised_only" \
    --geno 0.01 \
    --make-bed \
    --out "${SCRAP}/geno_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/geno_filtered")
SNPS_REMOVED=$((SNPS_BEFORE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs with >1% missingness in supervised samples"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 2 — Minor allele frequency (--maf 0.03) [supervised only]
# ===================================================================
echo "  Step 2/8: MAF filter (--maf 0.03) [supervised only]"

SNPS_BEFORE_MAF=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/geno_filtered" \
    --maf 0.03 \
    --make-bed \
    --out "${SCRAP}/maf_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/maf_filtered")
SNPS_REMOVED=$((SNPS_BEFORE_MAF - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs with MAF < 0.03 in supervised samples"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 3 — Hardy-Weinberg equilibrium (--hwe 1e-50) [supervised only]
# ===================================================================
echo "  Step 3/8: HWE exact test (--hwe 1e-50) [supervised only]"

SNPS_BEFORE_HWE=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/maf_filtered" \
    --hwe 1e-50 \
    --make-bed \
    --out "${SCRAP}/hwe_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_AFTER=$(count_snps "${SCRAP}/hwe_filtered")
SNPS_REMOVED=$((SNPS_BEFORE_HWE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs failing HWE (p < 1e-50)"
echo "    Remaining: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 4 — Exclude long-range LD regions [supervised only]
# ===================================================================
echo "  Step 4/8: Exclude long-range LD regions (Price et al. 2008, hg38) [supervised only]"

# Download the hg38 high-LD regions BED file
HIGH_LD_BED="${SCRAP}/high_ld_regions_hg38.bed"
if [[ ! -f "${HIGH_LD_BED}" ]]; then
    echo "    Downloading high-LD regions ..."
    curl -fSL -o "${HIGH_LD_BED}" \
        "https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed"
fi

# Convert BED (0-based) to PLINK range format (1-based)
awk 'BEGIN{OFS="\t"} {
    chr = $1; sub(/^chr/, "", chr)
    print chr, $2+1, $3, $4
}' "${HIGH_LD_BED}" > "${SCRAP}/high_ld_ranges.txt"

SNPS_BEFORE_LD=${SNPS_AFTER}
"${PLINK2}" --bfile "${SCRAP}/hwe_filtered" \
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
# STEP 5 — LD pruning (--indep-pairwise 50 10 0.1) [supervised only]
# ===================================================================
echo "  Step 5/8: LD pruning (window=50, step=10, r2=0.1) [supervised only]"

"${PLINK2}" --bfile "${SCRAP}/ldregion_excluded" \
    --indep-pairwise 50 10 0.1 \
    --out "${SCRAP}/ld_prune" \
    "${PLINK_FLAGS[@]}"

SNPS_BEFORE_PRUNE=${SNPS_AFTER}

# Save the final SNP list — this is what we'll apply to all samples
cp "${SCRAP}/ld_prune.prune.in" "${SCRAP}/supervised_snps.txt"
SNPS_AFTER=$(wc -l < "${SCRAP}/supervised_snps.txt")
SNPS_REMOVED=$((SNPS_BEFORE_PRUNE - SNPS_AFTER))
echo "    Removed $(fmt ${SNPS_REMOVED}) SNPs by LD pruning"
echo "    Final supervised SNP list: $(fmt ${SNPS_AFTER}) SNPs"
echo ""

# ===================================================================
# STEP 6 — Apply supervised SNP list to ALL samples
# ===================================================================
echo "  Step 6/8: Apply supervised SNP list to all samples"

"${PLINK2}" --bfile "${INPUT}" \
    --extract "${SCRAP}/supervised_snps.txt" \
    --make-bed \
    --out "${SCRAP}/all_snp_filtered" \
    "${PLINK_FLAGS[@]}"

SNPS_FINAL=$(count_snps "${SCRAP}/all_snp_filtered")
SAMPLES_ALL=$(count_samples "${SCRAP}/all_snp_filtered")
echo "    All samples with supervised SNPs: $(fmt ${SAMPLES_ALL}) samples, $(fmt ${SNPS_FINAL}) SNPs"
echo ""

# ===================================================================
# STEP 7 — Individual missingness (--mind 0.01) [all samples]
# ===================================================================
echo "  Step 7/8: Individual missingness filter (--mind 0.01) [all samples]"

SAMPLES_BEFORE_MIND=$(count_samples "${SCRAP}/all_snp_filtered")

"${PLINK2}" --bfile "${SCRAP}/all_snp_filtered" \
    --mind 0.01 \
    --make-bed \
    --out "${SCRAP}/mind_filtered" \
    "${PLINK_FLAGS[@]}"

SAMPLES_AFTER_MIND=$(count_samples "${SCRAP}/mind_filtered")
SAMPLES_REMOVED_MIND=$((SAMPLES_BEFORE_MIND - SAMPLES_AFTER_MIND))

if [[ ${SAMPLES_REMOVED_MIND} -gt 0 ]]; then
    echo "    Removed $(fmt ${SAMPLES_REMOVED_MIND}) individuals with >1% missingness"

    awk '{print $2}' "${SCRAP}/all_snp_filtered.fam" | sort > "${SCRAP}/before_mind_iids.txt"
    awk '{print $2}' "${SCRAP}/mind_filtered.fam" | sort > "${SCRAP}/after_mind_iids.txt"
    comm -23 "${SCRAP}/before_mind_iids.txt" "${SCRAP}/after_mind_iids.txt" > "${SCRAP}/removed_mind_iids.txt"

    "${PYTHON}" -c "
import pandas as pd, sys
meta = pd.read_csv('${PROJECT_DIR}/summary/metadata.csv')
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
removed = set(open('${SCRAP}/removed_mind_iids.txt').read().split())
meta_removed = meta[meta['sample_id'].isin(removed)]
print('    Removed individuals by population:')
for pop, group in meta_removed.groupby('population_id'):
    sup_match = sup[sup['population_id'] == pop]
    ref_pop = sup_match['reference_population'].iloc[0] if len(sup_match) > 0 else 'non-reference'
    ids = group['sample_id'].tolist()
    if len(ids) <= 5:
        print(f'      {pop} ({ref_pop}): {\", \".join(ids)}')
    else:
        print(f'      {pop} ({ref_pop}): {len(ids)} samples')
"
else
    echo "    No individuals removed (all pass missingness threshold)"
fi
echo "    Remaining: $(fmt ${SAMPLES_AFTER_MIND}) individuals"
echo ""

# ===================================================================
# STEP 8 — Kinship (separate AMR and non-AMR) [all samples]
# ===================================================================
echo "  Step 8/8: Kinship filtering (AMR threshold=0.088, non-AMR threshold=0.05)"

"${PYTHON}" -c "
import pandas as pd
meta = pd.read_csv('${PROJECT_DIR}/summary/metadata.csv')
fam = pd.read_csv('${SCRAP}/mind_filtered.fam', sep='\s+', header=None, names=['FID','IID','PAT','MAT','SEX','PHENO'])
amr_ids = set(meta.loc[meta['superpopulation'] == 'American', 'sample_id'])
amr_in_fam = fam[fam['IID'].isin(amr_ids)]
nonamr_in_fam = fam[~fam['IID'].isin(amr_ids)]
amr_in_fam[['FID', 'IID']].to_csv('${SCRAP}/amr_samples.txt', sep='\t', header=False, index=False)
nonamr_in_fam[['FID', 'IID']].to_csv('${SCRAP}/nonamr_samples.txt', sep='\t', header=False, index=False)
print(f'    AMR samples: {len(amr_in_fam):,}')
print(f'    Non-AMR samples: {len(nonamr_in_fam):,}')
"

# --- 8a. Kinship on AMR (threshold 0.088) ---
echo ""
echo "    8a. AMR kinship (KING cutoff = 0.088)"

AMR_COUNT_BEFORE=$(wc -l < "${SCRAP}/amr_samples.txt")

"${PLINK2}" --bfile "${SCRAP}/mind_filtered" \
    --keep "${SCRAP}/amr_samples.txt" \
    --king-cutoff 0.088 \
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
meta = pd.read_csv('${PROJECT_DIR}/summary/metadata.csv')
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
removed = pd.read_csv('${SCRAP}/amr_king.king.cutoff.out.id', sep='\t')
removed_ids = set(removed.iloc[:, 1] if removed.columns[0] == '#FID' else removed.iloc[:, 0])
meta_removed = meta[meta['sample_id'].isin(removed_ids)]
print('        Removed AMR individuals by population:')
for pop, group in meta_removed.groupby('population_id'):
    sup_match = sup[sup['population_id'] == pop]
    ref_pop = sup_match['reference_population'].iloc[0] if len(sup_match) > 0 else 'non-reference'
    ids = group['sample_id'].tolist()
    if len(ids) <= 5:
        print(f'          {pop} ({ref_pop}): {\", \".join(ids)}')
    else:
        print(f'          {pop} ({ref_pop}): {len(ids)} samples')
"
fi

# --- 8b. Kinship on non-AMR (threshold 0.05) ---
echo ""
echo "    8b. Non-AMR kinship (KING cutoff = 0.05)"

NONAMR_COUNT_BEFORE=$(wc -l < "${SCRAP}/nonamr_samples.txt")

"${PLINK2}" --bfile "${SCRAP}/mind_filtered" \
    --keep "${SCRAP}/nonamr_samples.txt" \
    --king-cutoff 0.05 \
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
meta = pd.read_csv('${PROJECT_DIR}/summary/metadata.csv')
sup = pd.read_csv('${PROJECT_DIR}/summary/supervised.csv')
removed = pd.read_csv('${SCRAP}/nonamr_king.king.cutoff.out.id', sep='\t')
removed_ids = set(removed.iloc[:, 1] if removed.columns[0] == '#FID' else removed.iloc[:, 0])
meta_removed = meta[meta['sample_id'].isin(removed_ids)]
print('        Removed non-AMR individuals by population:')
for pop, group in meta_removed.groupby('population_id'):
    sup_match = sup[sup['population_id'] == pop]
    ref_pop = sup_match['reference_population'].iloc[0] if len(sup_match) > 0 else 'non-reference'
    ids = group['sample_id'].tolist()
    if len(ids) <= 5:
        print(f'          {pop} ({ref_pop}): {\", \".join(ids)}')
    else:
        print(f'          {pop} ({ref_pop}): {len(ids)} samples')
"
fi

# --- Combine kept samples from both kinship runs ---
echo ""
echo "    Combining kept samples from AMR and non-AMR kinship ..."

tail -n +2 "${SCRAP}/amr_king.king.cutoff.in.id" > "${SCRAP}/kinship_keep.txt"
tail -n +2 "${SCRAP}/nonamr_king.king.cutoff.in.id" >> "${SCRAP}/kinship_keep.txt"

TOTAL_REMOVED=$((AMR_REMOVED + NONAMR_REMOVED))
echo "    Total removed by kinship: $(fmt ${TOTAL_REMOVED})"

"${PLINK2}" --bfile "${SCRAP}/mind_filtered" \
    --keep "${SCRAP}/kinship_keep.txt" \
    --make-bed \
    --out "${SUPERVISED_ADMIXTURE}/ancestry_qc" \
    "${PLINK_FLAGS[@]}"

SAMPLES_FINAL=$(count_samples "${SUPERVISED_ADMIXTURE}/ancestry_qc")
SNPS_FINAL=$(count_snps "${SUPERVISED_ADMIXTURE}/ancestry_qc")
echo "    After kinship: $(fmt ${SAMPLES_FINAL}) individuals, $(fmt ${SNPS_FINAL}) SNPs"
echo ""

# ===================================================================
# Summary
# ===================================================================
TOTAL_SNPS_REMOVED=$((SNPS_BEFORE - SNPS_FINAL))
TOTAL_SAMPLES_REMOVED=$((SAMPLES_BEFORE - SAMPLES_FINAL))

echo "  ================================================"
echo "  QC Summary"
echo "  ================================================"
echo "  SNP filters computed on supervised samples only (n=$(fmt ${SUP_SAMPLES}))"
echo "  SNPs:    $(fmt ${SNPS_BEFORE}) -> $(fmt ${SNPS_FINAL})  (removed $(fmt ${TOTAL_SNPS_REMOVED}))"
echo "  Samples: $(fmt ${SAMPLES_BEFORE}) -> $(fmt ${SAMPLES_FINAL})  (removed $(fmt ${TOTAL_SAMPLES_REMOVED}))"
echo "  Output:  ${SUPERVISED_ADMIXTURE}/ancestry_qc.{bed,bim,fam}"
echo ""

# Clean up logs
rm -f "${SCRAP}"/*.log "${SCRAP}"/*.nosex
rm -f "${SCRAP}"/before_mind_iids.txt "${SCRAP}"/after_mind_iids.txt

echo "  ADMIXTURE QC complete."
