#!/usr/bin/env bash
#
# download_files.sh
# ------------------
# Downloads 1000 Genomes (KG), HGDP, SGDP, and Neural ADMIXTURE reference files.
#
# Expected environment:
#   DOWNLOADS_DIR  — destination directory for all downloaded files (must exist)
#
# Outputs (written to DOWNLOADS_DIR):
#   kg_all.{pgen.zst,pvar.zst,psam}     — 1000 Genomes pfiles (hg38)
#   deg2_hg38.king.cutoff.out.id         — KG related-sample exclusion list
#   hgdp_all.{pgen.zst,pvar.zst,psam}   — HGDP pfiles (hg38)
#   sgdp_all.{bed,bim.zip,fam}           — SGDP bed/bim/fam (hg19)
#   sgdp_metadata.txt                    — SGDP sample metadata
#   neural/data/                         — Neural ADMIXTURE pretrained data
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Validate environment
# ---------------------------------------------------------------------------
if [[ -z "${DOWNLOADS_DIR:-}" ]]; then
    echo "Error: DOWNLOADS_DIR is not set." >&2
    exit 1
fi

if [[ ! -d "${DOWNLOADS_DIR}" ]]; then
    echo "Error: DOWNLOADS_DIR does not exist: ${DOWNLOADS_DIR}" >&2
    exit 1
fi

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
download() {
    local dest="$1" url="$2"
    if [[ -f "${dest}" ]]; then
        echo "  [skip] $(basename "${dest}") already exists"
        return 0
    fi
    echo "  [download] $(basename "${dest}")"
    curl -sSf -o "${dest}" "${url}"
}

# ---------------------------------------------------------------------------
# 1000 Genomes (KG) — hg38 pfiles
# ---------------------------------------------------------------------------
echo "==> Downloading 1000 Genomes (KG) data ..."

download "${DOWNLOADS_DIR}/kg_all.pgen.zst" \
    "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst"

download "${DOWNLOADS_DIR}/kg_all.pvar.zst" \
    "https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz"

download "${DOWNLOADS_DIR}/kg_all.psam" \
    "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x"

download "${DOWNLOADS_DIR}/deg2_hg38.king.cutoff.out.id" \
    "https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id"

# ---------------------------------------------------------------------------
# HGDP — hg38 pfiles
# ---------------------------------------------------------------------------
echo "==> Downloading HGDP data ..."

download "${DOWNLOADS_DIR}/hgdp_all.pgen.zst" \
    "https://www.dropbox.com/s/1hssuhgqg4r345f/hgdp_statphase.pgen.zst"

download "${DOWNLOADS_DIR}/hgdp_all.pvar.zst" \
    "https://www.dropbox.com/s/y7iza47zf75fhy4/hgdp_statphase.pvar.zst"

download "${DOWNLOADS_DIR}/hgdp_all.psam" \
    "https://www.dropbox.com/s/0zg57558fqpj3w1/hgdp.psam"

# ---------------------------------------------------------------------------
# SGDP — hg19 bed/bim/fam
# ---------------------------------------------------------------------------
echo "==> Downloading SGDP data ..."

download "${DOWNLOADS_DIR}/sgdp_all.bim.zip" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bim.zip"

download "${DOWNLOADS_DIR}/sgdp_all.bed" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bed"

download "${DOWNLOADS_DIR}/sgdp_all.fam" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.fam"

download "${DOWNLOADS_DIR}/sgdp_metadata.txt" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt"

# ---------------------------------------------------------------------------
# Neural ADMIXTURE — pretrained data from Figshare
# ---------------------------------------------------------------------------
echo "==> Downloading Neural ADMIXTURE data ..."

NEURAL_DIR="${DOWNLOADS_DIR}/neural"
if [[ -d "${NEURAL_DIR}/allchms" ]]; then
    echo "  [skip] Neural ADMIXTURE data already exists"
else
    mkdir -p "${NEURAL_DIR}"
    echo "  [download] Neural ADMIXTURE data.tar.gz"
    curl -LJO --output-dir "${NEURAL_DIR}" https://ndownloader.figshare.com/files/34438760
    # tar may warn about hardlinks on macOS — ignore non-fatal errors
    tar -xzf "${NEURAL_DIR}/data.tar.gz" -C "${NEURAL_DIR}" || true
    rm -f "${NEURAL_DIR}/data.tar.gz"
fi

# ---------------------------------------------------------------------------
# GIAB Ashkenazi Jewish parents — hg38 benchmark VCFs + high-confidence BEDs
# ---------------------------------------------------------------------------
echo "==> Downloading GIAB Ashkenazi parent data ..."

GIAB_DIR="${DOWNLOADS_DIR}/giab"
mkdir -p "${GIAB_DIR}"

GIAB_BASE="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio"

download "${GIAB_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "${GIAB_BASE}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

download "${GIAB_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "${GIAB_BASE}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

download "${GIAB_DIR}/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "${GIAB_BASE}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"

download "${GIAB_DIR}/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "${GIAB_BASE}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"

# ---------------------------------------------------------------------------
# 1000 Genomes population descriptions (for plot labels)
# ---------------------------------------------------------------------------
echo "==> Downloading 1000 Genomes population descriptions ..."

download "${DOWNLOADS_DIR}/kg_population_names.tsv" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv"

echo "==> Downloads complete."
