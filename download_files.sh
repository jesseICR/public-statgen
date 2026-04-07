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
#   high_ld_regions_hg38.bed             — Long-range LD regions (hg38)
#
set -euo pipefail

# ---------------------------------------------------------------------------
# Generate checksums mode (run after a known-good download to populate blanks)
# Usage: DOWNLOADS_DIR=./downloads bash download_files.sh --generate-checksums
# ---------------------------------------------------------------------------
if [[ "${1:-}" == "--generate-checksums" ]]; then
    if [[ -z "${DOWNLOADS_DIR:-}" ]]; then
        echo "Error: DOWNLOADS_DIR is not set." >&2; exit 1
    fi
    echo "SHA-256 checksums for files in ${DOWNLOADS_DIR}:"
    for f in \
        kg_all.pgen.zst kg_all.pvar.zst kg_all.psam \
        deg2_hg38.king.cutoff.out.id \
        hgdp_all.pgen.zst hgdp_all.pvar.zst hgdp_all.psam \
        sgdp_all.bim.zip sgdp_all.bed sgdp_all.fam sgdp_metadata.txt \
        kg_population_names.tsv high_ld_regions_hg38.bed \
        giab/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        giab/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
        giab/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz \
        giab/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed; do
        path="${DOWNLOADS_DIR}/${f}"
        if [[ -f "${path}" ]]; then
            if command -v sha256sum &>/dev/null; then
                hash="$(sha256sum "${path}" | awk '{print $1}')"
            elif command -v shasum &>/dev/null; then
                hash="$(shasum -a 256 "${path}" | awk '{print $1}')"
            fi
            printf "  %-55s %s\n" "${f}" "${hash}"
        else
            printf "  %-55s [missing]\n" "${f}"
        fi
    done
    exit 0
fi

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

# Portable SHA-256: macOS has shasum, Linux has sha256sum
if command -v sha256sum &>/dev/null; then
    sha256() { sha256sum "$1" | awk '{print $1}'; }
elif command -v shasum &>/dev/null; then
    sha256() { shasum -a 256 "$1" | awk '{print $1}'; }
else
    sha256() { python3 -c "import hashlib,sys;print(hashlib.sha256(open(sys.argv[1],'rb').read()).hexdigest())" "$1"; }
fi

download() {
    local dest="$1" url="$2" expected_sha256="${3:-}"
    if [[ -f "${dest}" ]]; then
        echo "  [skip] $(basename "${dest}") already exists"
        return 0
    fi
    echo "  [download] $(basename "${dest}")"
    curl -sSfL -o "${dest}" "${url}"
    if [[ -n "${expected_sha256}" ]]; then
        local actual
        actual="$(sha256 "${dest}")"
        if [[ "${actual}" != "${expected_sha256}" ]]; then
            echo "  [CHECKSUM FAILED] $(basename "${dest}")" >&2
            echo "    Expected: ${expected_sha256}" >&2
            echo "    Actual:   ${actual}" >&2
            rm -f "${dest}"
            return 1
        fi
        echo "  [verified] $(basename "${dest}")"
    fi
}

# ---------------------------------------------------------------------------
# 1000 Genomes (KG) — hg38 pfiles
# ---------------------------------------------------------------------------
echo "==> Downloading 1000 Genomes (KG) data ..."

download "${DOWNLOADS_DIR}/kg_all.pgen.zst" \
    "https://www.dropbox.com/s/j72j6uciq5zuzii/all_hg38.pgen.zst" \
    ""  # ~5 GB — run 'bash download_files.sh --generate-checksums' to populate

download "${DOWNLOADS_DIR}/kg_all.pvar.zst" \
    "https://www.dropbox.com/scl/fi/fn0bcm5oseyuawxfvkcpb/all_hg38_rs.pvar.zst?rlkey=przncwb78rhz4g4ukovocdxaz" \
    ""  # ~0.7 GB

download "${DOWNLOADS_DIR}/kg_all.psam" \
    "https://www.dropbox.com/scl/fi/u5udzzaibgyvxzfnjcvjc/hg38_corrected.psam?rlkey=oecjnk4vmbhc8b1p202l0ih4x" \
    "727a01b637d5e1ef02516a53069dc552b1fa46f442004ad50ea3d069c7741cdd"

download "${DOWNLOADS_DIR}/deg2_hg38.king.cutoff.out.id" \
    "https://www.dropbox.com/s/4zhmxpk5oclfplp/deg2_hg38.king.cutoff.out.id" \
    "22094d03c593276ba88ea5c564e106bab2a755f0d2a3e6c77aaff2f8bcd83b66"

# ---------------------------------------------------------------------------
# HGDP — hg38 pfiles
# ---------------------------------------------------------------------------
echo "==> Downloading HGDP data ..."

download "${DOWNLOADS_DIR}/hgdp_all.pgen.zst" \
    "https://www.dropbox.com/s/1hssuhgqg4r345f/hgdp_statphase.pgen.zst" \
    ""  # ~2 GB

download "${DOWNLOADS_DIR}/hgdp_all.pvar.zst" \
    "https://www.dropbox.com/s/y7iza47zf75fhy4/hgdp_statphase.pvar.zst" \
    ""  # ~0.5 GB

download "${DOWNLOADS_DIR}/hgdp_all.psam" \
    "https://www.dropbox.com/s/0zg57558fqpj3w1/hgdp.psam" \
    "d4ffa5e1307303858ad0a420b2a6682e71bb00cd057791d8a5a872143f738cc5"

# ---------------------------------------------------------------------------
# SGDP — hg19 bed/bim/fam
# ---------------------------------------------------------------------------
echo "==> Downloading SGDP data ..."

download "${DOWNLOADS_DIR}/sgdp_all.bim.zip" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bim.zip" \
    "949dc6013e6c102a0c00555fa30b36ab2b4257e4a1bacd349ae7cc626986923e"

download "${DOWNLOADS_DIR}/sgdp_all.bed" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.bed" \
    ""  # ~1.3 GB

download "${DOWNLOADS_DIR}/sgdp_all.fam" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/variant_set/cteam_extended.v4.maf0.1perc.fam" \
    "6f4e33ae2e2f29d59b503dcd3a7cf43afc75f42d0d8ec9c52df9d7ae5cc0535c"

download "${DOWNLOADS_DIR}/sgdp_metadata.txt" \
    "https://sharehost.hms.harvard.edu/genetics/reich_lab/sgdp/SGDP_metadata.279public.21signedLetter.44Fan.samples.txt" \
    "1e4f0076ddadc57ce9feb252c1f9a782ea9d5ba4e6606627fcd60e9c5bf2ef33"

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
    "${GIAB_BASE}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    ""  # ~0.4 GB

download "${GIAB_DIR}/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "${GIAB_BASE}/HG003_NA24149_father/latest/GRCh38/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "652afd3046705af3200f9c87c255fef11bb212dd76c75a19999c9b2df8a3180c"

download "${GIAB_DIR}/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    "${GIAB_BASE}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz" \
    ""  # ~0.4 GB

download "${GIAB_DIR}/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "${GIAB_BASE}/HG004_NA24143_mother/latest/GRCh38/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed" \
    "88d1c926fdc8abd9c39b3e9fe2af3fc43dd67d284690a3dc6127b369d96bccad"

# ---------------------------------------------------------------------------
# 1000 Genomes population descriptions (for plot labels)
# ---------------------------------------------------------------------------
echo "==> Downloading 1000 Genomes population descriptions ..."

download "${DOWNLOADS_DIR}/kg_population_names.tsv" \
    "https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/20131219.populations.tsv" \
    "420932693a5997e1d47f634b93eab861b87d5e5d23ac141e75163426f8076aaa"

# ---------------------------------------------------------------------------
# High-LD regions BED (for ADMIXTURE QC LD exclusion)
# ---------------------------------------------------------------------------
echo "==> Downloading high-LD regions ..."

download "${DOWNLOADS_DIR}/high_ld_regions_hg38.bed" \
    "https://raw.githubusercontent.com/meyer-lab-cshl/plinkQC/master/inst/extdata/high-LD-regions-hg38-GRCh38.bed" \
    "6e2ba6abf0151209ed71f632c736f841a2f84e25aff17d1d169ffe52a62ee177"

echo "==> Downloads complete."
