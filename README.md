# Public Statistical Genetics Pipeline

A fully reproducible pipeline for harmonizing, aligning, and performing rigorous quality control on the 1000 Genomes (KG), Human Genome Diversity Project (HGDP), and Simons Genome Diversity Project (SGDP) reference panels.

## Quick Start

```bash
bash main.sh
```

That's it. The pipeline is idempotent — re-running skips any steps that have already completed.

## Requirements

- **Python 3** (for SGDP liftover / QC; a venv is created automatically)
- **wget** and **curl** (for downloading data and tools)
- **unzip** (for extracting PLINK archives and SGDP bim)
- macOS (Intel or Apple Silicon) or Linux (x86_64)

PLINK 1.9 and 2.0 are installed automatically by the pipeline.
Python dependencies (`pandas`, `liftover`) are installed into a local venv.

## Pipeline Overview

`main.sh` is the top-level orchestrator. It runs each stage in order and exports shared paths (`DOWNLOADS_DIR`, `QC_DIR`, `PLINK1`, `PLINK2`, `PYTHON`, `SNPS_FILE`, `PLINK_MEMORY`, `PLINK_THREADS`) so that all downstream scripts use consistent locations.

### Current stages

| Stage | Script | Description |
|-------|--------|-------------|
| 1 | `setup_plink.sh` | Downloads and installs PLINK 1.9 & 2.0 into `tools/bin/` |
| 2 | `download_files.sh` | Downloads KG, HGDP, SGDP, and Neural ADMIXTURE files into `downloads/` |
| 3 | `qc_kg_hgdp.sh` | Decompress, filter (autosomal, biallelic, SNP extract, remove KG relatives), and convert to bed/bim/fam in `qc/` |
| 4 | `setup_python.sh` | Creates a Python venv in `tools/venv/` and installs dependencies from `requirements.txt` |
| 5 | `qc_sgdp.py` | Lifts SGDP from hg19→hg38, matches to KG by (chrom, pos, alleles), assigns rsIDs, outputs `qc/sgdp_qc.{bed,bim,fam}` |
| 6 | `merge_kg_hgdp_sgdp.sh` | Aligns alleles, merges KG+HGDP, deduplicates SGDP samples, three-way merge into `merge/merged_kg_hgdp_sgdp.{bed,bim,fam}` |
| 7 | `build_metadata.py` | Merges KG, HGDP, SGDP metadata with Neural ADMIXTURE ancestry labels into `summary/metadata.csv` |
| 8 | `build_supervised.py` | Assigns samples to K=8 supervised ADMIXTURE reference populations, writes `summary/supervised.csv` |

Additional stages (projection, ADMIXTURE) will be added to `main.sh` as the pipeline grows.

## Directory Structure

```
public-statgen/
├── main.sh                      # Run this
├── setup_plink.sh               # Stage 1: install PLINK
├── download_files.sh             # Stage 2: download reference panels + neural data
├── qc_kg_hgdp.sh               # Stage 3: QC KG and HGDP
├── setup_python.sh              # Stage 4: set up Python venv
├── qc_sgdp.py                   # Stage 5: QC SGDP (liftover + rsID match)
├── merge_kg_hgdp_sgdp.sh       # Stage 6: three-way merge
├── build_metadata.py            # Stage 7: merge metadata + neural ancestry
├── build_supervised.py          # Stage 8: supervised ADMIXTURE reference populations
├── requirements.txt             # Python dependencies
├── rsids_dense_chr1_22.txt      # SNP list for filtering
├── tools/
│   ├── bin/                     # PLINK binaries (created by setup_plink.sh)
│   │   ├── plink1
│   │   └── plink2
│   └── venv/                    # Python venv (created by setup_python.sh)
├── downloads/                   # Raw downloaded data (created by main.sh)
│   ├── kg_all.{pgen.zst,pvar.zst,psam}
│   ├── deg2_hg38.king.cutoff.out.id
│   ├── hgdp_all.{pgen.zst,pvar.zst,psam}
│   ├── sgdp_all.{bed,bim.zip,fam}
│   ├── sgdp_metadata.txt
│   └── neural/                   # Neural ADMIXTURE pretrained data
│       ├── allchms/data/
│       ├── chm22/data/
│       └── chm22sim/data/
├── qc/                          # QC'd bed/bim/fam
│   ├── kg_qc.{bed,bim,fam}
│   ├── hgdp_qc.{bed,bim,fam}
│   └── sgdp_qc.{bed,bim,fam}
├── merge/                       # Merged fileset
│   └── merged_kg_hgdp_sgdp.{bed,bim,fam}
└── literature_reference/        # Sample-level info extracted from publications
    ├── sharma_all_of_us_2025.csv
    ├── marino_creatinine_2022.csv
    ├── koenig_harmonized_outliers_2024.txt
    ├── ancestry_martin_outliers_2017.csv
    ├── other_spanish_outliers.txt
    └── american_admixed_outliers.txt
```

## Data Sources

- **1000 Genomes (KG)** — hg38 pfiles from the [PLINK 2.0 resources page](https://www.cog-genomics.org/plink/2.0/resources)
- **HGDP** — hg38 pfiles (statistically phased) from the same source
- **SGDP** — hg19 bed/bim/fam from the [Reich Lab](https://reichdata.hms.harvard.edu/pub/datasets/sgdp/)

## Literature Reference

The `literature_reference/` directory contains sample-level information extracted from publications that characterize the KG, HGDP, and SGDP reference panels. These files document ancestry outliers, population compositions, and other metadata used in the pipeline's quality control and interpretation.

| File | Description |
|------|-------------|
| `sharma_all_of_us_2025.csv` | Reference populations and sample counts used in the All of Us ancestry analysis |
| `marino_creatinine_2022.csv` | Reference populations used in creatinine ancestry analysis |
| `koenig_harmonized_outliers_2024.txt` | Sample IDs of outliers identified during harmonization of diverse human genomes |
| `ancestry_martin_outliers_2017.csv` | Samples with considerable admixture identified in genetic risk prediction analysis |
| `other_spanish_outliers.txt` | Additional Spanish-ancestry outlier samples |
| `american_admixed_outliers.txt` | American reference samples excluded due to admixture |

### Publications

- **Koenig et al.** — *A harmonized public resource of deeply sequenced diverse human genomes.* Genome Research (2024). [doi:10.1101/gr.278378.123](https://doi.org/10.1101/gr.278378.123)

- **Martin et al.** — *Human Demographic History Impacts Genetic Risk Prediction across Diverse Populations.* American Journal of Human Genetics (2017). [doi:10.1016/j.ajhg.2017.03.004](https://doi.org/10.1016/j.ajhg.2017.03.004)

- **Sharma et al.** — *Genetic ancestry and population structure in the All of Us Research Program cohort.* bioRxiv (2024). [doi:10.1101/2024.12.21.629909](https://doi.org/10.1101/2024.12.21.629909)

- **Marino-Ramirez et al.** — *Effects of genetic ancestry and socioeconomic deprivation on ethnic differences in serum creatinine.* Gene (2022). [doi:10.1016/j.gene.2022.146709](https://doi.org/10.1016/j.gene.2022.146709)

- **Dominguez Mantes et al.** — *Neural ADMIXTURE for rapid genomic clustering.* Nature Computational Science (2023). [doi:10.1038/s43588-023-00482-7](https://doi.org/10.1038/s43588-023-00482-7) — Used in the pipeline for ancestry label assignment but not directly represented in the literature reference files.
