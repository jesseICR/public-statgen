#!/usr/bin/env bash
#
# qc_kg_hgdp.sh
# --------------
# QC and format conversion for 1000 Genomes (KG) and HGDP pfiles.
#
# Steps:
#   1. Decompress .pgen.zst / .pvar.zst from DOWNLOADS_DIR into QC_DIR
#   2. Filter (autosomal, biallelic, SNP extract, remove KG relatives)
#      and convert directly to bed/bim/fam
#   3. Clean up intermediate decompressed files
#
# Expected environment:
#   DOWNLOADS_DIR  — directory containing the raw downloaded pfiles
#   QC_DIR         — output directory for QC'd bed/bim/fam files
#   PLINK2         — path to plink2 binary
#   SNPS_FILE      — path to SNP rsID list for --extract
#   PLINK_MEMORY   — memory limit in MB for PLINK
#   PLINK_THREADS  — number of threads for PLINK
#
# Outputs (in QC_DIR):
#   kg_qc.{bed,bim,fam}
#   hgdp_qc.{bed,bim,fam}
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
for var in DOWNLOADS_DIR QC_DIR PLINK2 SNPS_FILE PLINK_MEMORY PLINK_THREADS; do
    if [[ -z "${!var:-}" ]]; then
        echo "Error: ${var} is not set." >&2
        exit 1
    fi
done

# ---------------------------------------------------------------------------
# Shared PLINK2 flags
# ---------------------------------------------------------------------------
PLINK2_COMMON=(
    --chr 1-22
    --max-alleles 2
    --extract "${SNPS_FILE}"
    --set-missing-var-ids '@_#_$1_$2'
    --allow-extra-chr
    --memory "${PLINK_MEMORY}"
    --threads "${PLINK_THREADS}"
)

# ---------------------------------------------------------------------------
# 1. Decompress zst files into QC_DIR
# ---------------------------------------------------------------------------
echo "  Decompressing KG pfiles ..."
"${PLINK2}" --zst-decompress "${DOWNLOADS_DIR}/kg_all.pvar.zst"   > "${QC_DIR}/kg_all.pvar"
"${PLINK2}" --zst-decompress "${DOWNLOADS_DIR}/kg_all.pgen.zst"   > "${QC_DIR}/kg_all.pgen"
cp "${DOWNLOADS_DIR}/kg_all.psam" "${QC_DIR}/kg_all.psam"

echo "  Decompressing HGDP pfiles ..."
"${PLINK2}" --zst-decompress "${DOWNLOADS_DIR}/hgdp_all.pvar.zst" > "${QC_DIR}/hgdp_all.pvar"
"${PLINK2}" --zst-decompress "${DOWNLOADS_DIR}/hgdp_all.pgen.zst" > "${QC_DIR}/hgdp_all.pgen"
cp "${DOWNLOADS_DIR}/hgdp_all.psam" "${QC_DIR}/hgdp_all.psam"

# ---------------------------------------------------------------------------
# 2. QC filter → bed/bim/fam (single pass)
# ---------------------------------------------------------------------------
echo "  QC filtering KG (removing related individuals) → bed ..."
"${PLINK2}" --pfile "${QC_DIR}/kg_all" \
    --remove "${DOWNLOADS_DIR}/deg2_hg38.king.cutoff.out.id" \
    "${PLINK2_COMMON[@]}" \
    --make-bed \
    --out "${QC_DIR}/kg_qc"

echo "  QC filtering HGDP → bed ..."
"${PLINK2}" --pfile "${QC_DIR}/hgdp_all" \
    "${PLINK2_COMMON[@]}" \
    --make-bed \
    --out "${QC_DIR}/hgdp_qc"

# ---------------------------------------------------------------------------
# 3. Clean up intermediate decompressed files
# ---------------------------------------------------------------------------
echo "  Cleaning up intermediate files ..."
rm -f "${QC_DIR}"/kg_all.{pgen,pvar,psam}
rm -f "${QC_DIR}"/hgdp_all.{pgen,pvar,psam}
rm -f "${QC_DIR}"/kg_qc.log "${QC_DIR}"/hgdp_qc.log

echo "  QC complete."
