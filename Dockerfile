# =============================================================================
# Public Statistical Genetics Pipeline — Docker Image
# =============================================================================
# Builds a reproducible Linux x86_64 image with all bioinformatics tools
# (PLINK 1.9, PLINK 2.0, UCSC liftOver, ADMIXTURE) and Python dependencies
# pre-installed. Data downloads happen at runtime.
#
# Build:
#   docker build -t public-statgen .
#
# Run (defaults: K=6, MAF=1%):
#   docker run --rm -v $(pwd)/data:/app/pipeline-output public-statgen
#
# Run with custom parameters:
#   docker run --rm -e K_MODEL=3 -e MAF_ADMIXTURE=0.0200 \
#     -v $(pwd)/data:/app/pipeline-output public-statgen
#
# =============================================================================

FROM python:3.11-slim-bookworm

LABEL org.opencontainers.image.source="https://github.com/jesseICR/public-statgen"
LABEL org.opencontainers.image.description="Reproducible statistical genetics pipeline for ancestry estimation"
LABEL org.opencontainers.image.licenses="MIT"

# ---- System dependencies ----------------------------------------------------
RUN apt-get update && apt-get install -y --no-install-recommends \
        curl \
        unzip \
        gawk \
        bc \
        zstd \
        gzip \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# ---- Python dependencies (cached layer) -------------------------------------
COPY requirements.txt .
RUN python3 -m venv /app/tools/venv \
    && /app/tools/venv/bin/pip install --no-cache-dir --upgrade pip \
    && /app/tools/venv/bin/pip install --no-cache-dir -r requirements.txt

# ---- Pipeline code -----------------------------------------------------------
COPY . .

# ---- Pre-install bioinformatics tools ----------------------------------------
# These scripts detect Linux x86_64 and download the correct binaries.
RUN bash setup_plink.sh \
    && bash setup_liftover.sh \
    && bash setup_admixture.sh

# ---- Verify tools work -------------------------------------------------------
RUN tools/bin/plink1 --version > /dev/null 2>&1 \
    && tools/bin/plink2 --version > /dev/null 2>&1 \
    && test -x tools/bin/liftOver \
    && tools/bin/admixture --version 2>&1 | head -1

# ---- Default configuration (override with -e at runtime) ---------------------
ENV K_MODEL=6
ENV MAF_ADMIXTURE=0.0100

# ---- Runtime -----------------------------------------------------------------
# The pipeline generates ~15 GB of intermediate data and needs ~91 GB peak.
# Mount a volume at /app/pipeline-output to persist results, or run in /app
# and accept that data lives inside the container.
ENTRYPOINT ["bash", "main.sh"]
