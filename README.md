# Public Statistical Genetics Pipeline

A fully reproducible bioinformatics pipeline for harmonizing, aligning, and performing rigorous quality control on the 1000 Genomes (KG), Human Genome Diversity Project (HGDP), Simons Genome Diversity Project (SGDP), and Genome in a Bottle (GIAB) Ashkenazi Jewish reference panels, followed by supervised ADMIXTURE ancestry estimation. The pipeline merges four major human genetic diversity datasets into a single reference panel of 4,324 samples, runs supervised ADMIXTURE with cross-validation, and produces continental ancestry fractions, structure plots, and estimated allele frequencies — all from a single `bash main.sh` command.

## Quick Start

```bash
bash main.sh
```

The pipeline will prompt for two configuration choices:

1. **ADMIXTURE K model** — K=3, K=5, or K=6 (default: K=6)
2. **MAF threshold** — minor allele frequency percentage for ADMIXTURE QC (default: 1%)

After that it runs unattended. The pipeline is idempotent — re-running skips any steps that have already completed.

## Requirements

- **Python 3** (for SGDP liftover / QC; a venv is created automatically)
- **curl** (for downloading data and tools)
- macOS (Intel or Apple Silicon), Linux (x86_64), or Windows via WSL

PLINK 1.9, PLINK 2.0, and ADMIXTURE are installed automatically by the pipeline.
Python dependencies (`pandas`, `liftover`, `numpy`, `matplotlib`) are installed into a local venv.

## Pipeline Overview

`main.sh` is the top-level orchestrator. It runs each stage in order and exports shared paths (`DOWNLOADS_DIR`, `QC_DIR`, `MERGE_DIR`, `PLINK1`, `PLINK2`, `PYTHON`, `SNPS_FILE`, `PLINK_MEMORY`, `PLINK_THREADS`, `SUPERVISED_ADMIXTURE`, `ADMIXTURE`) so that all downstream scripts use consistent locations.

### Stages

| Stage | Script | Description |
|-------|--------|-------------|
| 1 | `setup_plink.sh` | Downloads and installs PLINK 1.9 & 2.0 into `tools/bin/` |
| 2 | `download_files.sh` | Downloads KG, HGDP, SGDP, Neural ADMIXTURE, and GIAB files into `downloads/` |
| 3 | `qc_kg_hgdp.sh` | Decompress, filter (autosomal, biallelic, SNP extract, remove KG relatives), and convert to bed/bim/fam in `qc/` |
| 4 | `setup_python.sh` | Creates a Python venv in `tools/venv/` and installs dependencies from `requirements.txt` |
| 5 | `qc_sgdp.py` | Lifts SGDP from hg19 to hg38, matches to KG by (chrom, pos, alleles), assigns rsIDs, outputs `qc/sgdp_qc.{bed,bim,fam}` |
| 6 | `merge_kg_hgdp_sgdp.sh` | Aligns alleles, merges KG+HGDP, deduplicates SGDP samples, three-way merge, removes ambiguous SNPs, applies geno filter |
| 7 | `prepare_giab.py` + `merge_giab.sh` | Converts GIAB Ashkenazi parent VCFs to PLINK (filling hom-ref from high-confidence BED), merges into reference panel, normalizes fam |
| 8 | `build_metadata.py` | Merges KG, HGDP, SGDP, GIAB metadata with Neural ADMIXTURE ancestry labels into `summary/metadata.csv` |
| 9 | `build_supervised.py` | Assigns samples to K=6 supervised ADMIXTURE reference populations, writes `summary/supervised.csv` |
| 10 | `setup_admixture.sh` | Downloads and installs ADMIXTURE into `tools/bin/` |
| 11 | `qc_admixture.sh` | QC on supervised samples only: geno, MAF, long-range LD exclusion, LD pruning, mind, kinship, HWE (on unrelated). Applies resulting SNP list to all samples for projection. |
| 12 | `run_admixture_supervised.py` | 3-fold stratified cross-validation + final supervised ADMIXTURE run (K=6) |
| 13 | `analyze_admixture_results.py` | Structure plots, cross-validation accuracy, augmented metadata with ancestry fractions, formatted allele frequency file |

## Directory Structure

```
public-statgen/
├── main.sh                      # Run this
├── setup_plink.sh               # Stage 1: install PLINK
├── download_files.sh            # Stage 2: download reference panels + GIAB + neural data
├── qc_kg_hgdp.sh               # Stage 3: QC KG and HGDP
├── setup_python.sh              # Stage 4: set up Python venv
├── qc_sgdp.py                   # Stage 5: QC SGDP (liftover + rsID match)
├── merge_kg_hgdp_sgdp.sh       # Stage 6: three-way merge
├── prepare_giab.py              # Stage 7a: convert GIAB VCFs to PLINK
├── merge_giab.sh                # Stage 7b: merge GIAB into reference panel
├── build_metadata.py            # Stage 8: merge metadata + neural ancestry
├── build_supervised.py          # Stage 9: supervised ADMIXTURE reference populations
├── setup_admixture.sh           # Stage 10: install ADMIXTURE
├── qc_admixture.sh              # Stage 11: QC for ADMIXTURE
├── run_admixture_supervised.py  # Stage 12: run ADMIXTURE supervised
├── analyze_admixture_results.py # Stage 13: analyze ADMIXTURE output
├── requirements.txt             # Python dependencies
├── rsids_dense_chr1_22.txt      # SNP list for filtering
├── tools/
│   ├── bin/                     # Binaries (created by setup scripts)
│   │   ├── plink1
│   │   ├── plink2
│   │   └── admixture
│   └── venv/                    # Python venv (created by setup_python.sh)
├── downloads/                   # Raw downloaded data (created by main.sh)
│   ├── kg_all.{pgen.zst,pvar.zst,psam}
│   ├── deg2_hg38.king.cutoff.out.id
│   ├── hgdp_all.{pgen.zst,pvar.zst,psam}
│   ├── sgdp_all.{bed,bim.zip,fam}
│   ├── sgdp_metadata.txt
│   ├── neural/                  # Neural ADMIXTURE pretrained data
│   └── giab/                    # GIAB Ashkenazi parents (VCFs + BEDs)
│       ├── HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
│       ├── HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
│       ├── HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
│       └── HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
├── qc/                          # QC'd bed/bim/fam
│   ├── kg_qc.{bed,bim,fam}
│   ├── hgdp_qc.{bed,bim,fam}
│   └── sgdp_qc.{bed,bim,fam}
├── merge/                       # Merged fileset (KG + HGDP + SGDP + GIAB)
│   └── merged_kg_hgdp_sgdp.{bed,bim,fam}
├── supervised_admixture/        # ADMIXTURE working directory
│   ├── ancestry_qc.{bed,bim,fam}
│   ├── fold_assignments.csv
│   ├── admixture_fold{1,2,3}.K.Q
│   ├── admixture_final.K.{Q,P}
│   └── scrap/                   # Intermediate QC files
├── summary/                     # Final outputs
│   ├── metadata.csv
│   ├── supervised.csv
│   ├── koenig_harmonized_outliers_2024.txt
│   └── admixture-global-K/
│       ├── metadata_ancestry.csv
│       ├── structure_holdout.png
│       ├── structure_projected.png
│       └── admixture_allele_freqs.tsv
├── outputs/                     # Git-tracked deterministic outputs
│   └── admixture-global-6/
│       ├── metadata_ancestry.csv
│       ├── structure_holdout.png
│       ├── structure_projected.png
│       └── admixture_allele_freqs.tsv
└── literature_reference/        # Sample-level info extracted from publications
```

## Runtime and Storage

**Benchmarked on:** 128 cores, 503 GiB RAM | PLINK threads: 6, PLINK memory: 14 GB

### Per-Step Runtime

| Step | Description | Runtime |
|------|-------------|---------|
| 1 | Install PLINK (1.9 + 2.0) | 1 s |
| 2 | Download KG, HGDP, SGDP, Neural ADMIXTURE, GIAB | ~7 m 30 s |
| 3 | QC KG and HGDP (decompress zst, filter, convert) | 9 m 49 s |
| 4 | Set up Python venv + dependencies | 11 s |
| 5 | QC SGDP (liftover hg19 → hg38, filter) | 1 m 38 s |
| 6 | Merge KG + HGDP + SGDP | 16 s |
| 7 | Prepare and merge GIAB Ashkenazi parents | 16 s |
| 8 | Build metadata CSV | < 1 s |
| 9 | Build supervised reference populations | < 1 s |
| 10 | Install ADMIXTURE | < 1 s |
| 11 | QC for ADMIXTURE (geno, MAF, HWE, LD prune, kinship) | 4 s |
| 12 | Run ADMIXTURE supervised (3-fold CV + final, K=6) | 70 m 15 s |
| 13 | Analyze ADMIXTURE results (CV stats, plots, metadata) | 16 s |

**Total pipeline runtime: ~90 minutes.**

The two dominant steps are Step 3 (QC KG + HGDP, ~10 min) and Step 12 (ADMIXTURE runs, ~70 min), which together account for ~89% of total runtime. Step 12 runs ADMIXTURE 4 times (3 CV folds + 1 final), each taking ~17–19 min wall time using 6 threads.

### Storage

| Directory | Size | Contents |
|-----------|------|----------|
| `downloads/` | 13 GB | Raw downloaded files (KG, HGDP pgen.zst, SGDP bed, Neural ADMIXTURE, GIAB VCFs) |
| `tools/` | 296 MB | PLINK 1.9 + 2.0, ADMIXTURE, Python venv |
| `qc/` | 496 MB | QC'd BED/BIM/FAM for KG, HGDP, SGDP |
| `merge/` | 389 MB | Merged three-dataset BED/BIM/FAM |
| `supervised_admixture/` | 858 MB | ADMIXTURE QC panel, fold .Q/.P files, final .Q/.P |
| `summary/` | 7.5 MB | metadata, supervised CSV, structure plots, allele freqs |
| **Total** | **~15 GB** | |

**Peak transient storage:** Step 3 decompresses KG and HGDP pgen.zst files (~5 GB each → ~8.9 GB uncompressed). Peak project size during step 3 is approximately **91 GB** before intermediates are removed.

The `downloads/` directory can be deleted after QC (steps 3 + 5) to reclaim ~13 GB, reducing the final footprint to ~2 GB.

## Data Sources

- **1000 Genomes (KG)** — hg38 pfiles from the [PLINK 2.0 resources page](https://www.cog-genomics.org/plink/2.0/resources)
- **HGDP** — hg38 pfiles (statistically phased) from the same source
- **SGDP** — hg19 bed/bim/fam from the [Reich Lab](https://reichdata.hms.harvard.edu/pub/datasets/sgdp/)
- **GIAB Ashkenazi Jewish trio** — hg38 benchmark VCFs from [NIST GIAB](https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/) (parents HG003 and HG004 only)
- **Neural ADMIXTURE** — pretrained ancestry model from [Figshare](https://doi.org/10.6084/m9.figshare.19387538.v1)

## GIAB Integration

The GIAB Ashkenazi Jewish parents (HG003 = father, HG004 = mother) are included as non-supervised samples for ancestry projection. They are not used in ADMIXTURE training — their ancestry fractions are estimated by the trained model.

Because GIAB benchmark VCFs only contain variant calls (not homozygous-reference sites), `prepare_giab.py` uses the high-confidence BED files to fill in reference-homozygous genotypes at panel SNP positions that fall within callable regions. This preserves the full SNP set during the merge.

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
| `oceanian_admixed_outliers.txt` | Oceanian reference samples excluded due to admixture |

### Publications

- **Koenig et al.** — *A harmonized public resource of deeply sequenced diverse human genomes.* Genome Research (2024). [doi:10.1101/gr.278378.123](https://doi.org/10.1101/gr.278378.123)

- **Martin et al.** — *Human Demographic History Impacts Genetic Risk Prediction across Diverse Populations.* American Journal of Human Genetics (2017). [doi:10.1016/j.ajhg.2017.03.004](https://doi.org/10.1016/j.ajhg.2017.03.004)

- **Sharma et al.** — *Genetic ancestry and population structure in the All of Us Research Program cohort.* bioRxiv (2024). [doi:10.1101/2024.12.21.629909](https://doi.org/10.1101/2024.12.21.629909)

- **Marino-Ramirez et al.** — *Effects of genetic ancestry and socioeconomic deprivation on ethnic differences in serum creatinine.* Gene (2022). [doi:10.1016/j.gene.2022.146709](https://doi.org/10.1016/j.gene.2022.146709)

- **Dominguez Mantes et al.** — *Neural ADMIXTURE for rapid genomic clustering.* Nature Computational Science (2023). [doi:10.1038/s43588-023-00482-7](https://doi.org/10.1038/s43588-023-00482-7)

- **Zook et al.** — *An open resource for accurately benchmarking small variant and reference calls.* Nature Biotechnology (2019). [doi:10.1038/s41587-019-0074-6](https://doi.org/10.1038/s41587-019-0074-6)

## ADMIXTURE Results (K=6)

The results below were obtained using the pipeline defaults (K=6, MAF=0.01). The pipeline produces a supervised ADMIXTURE analysis with six continental ancestry components: **African**, **American**, **East Asian**, **European**, **Oceanian**, and **South Asian**. A copy of the output files is checked into [`outputs/admixture-global-6/`](outputs/admixture-global-6/) so results are visible directly from the repository.

### Structure Plots

#### Cross-Validation Holdout

![Structure plot — supervised holdout](outputs/admixture-global-6/structure_holdout.png)

The holdout plot shows out-of-sample ancestry estimates for supervised reference samples. The pipeline uses 3-fold stratified cross-validation: in each fold, one-third of supervised samples are held out from training and their ancestry is predicted by the model trained on the remaining two-thirds. The three folds are combined into a single plot. Clean, single-color bars indicate the model accurately recovers known ancestry labels; mixed bars would indicate misclassification or genuine admixture within the reference panel. Samples are grouped by superpopulation and sorted geographically (roughly out-of-Africa order: African, European, East Asian, Oceanian, South Asian).

#### Projected Ancestry (Full Model)

![Structure plot — projected ancestry](outputs/admixture-global-6/structure_projected.png)

The projected plot shows ancestry estimates from the final ADMIXTURE model (trained on all supervised samples) applied to the full panel of 3,695 individuals, including unsupervised samples such as the two GIAB Ashkenazi Jewish parents (HG003 and HG004). Each vertical bar represents one individual; bar height shows the proportion assigned to each of the six ancestry components. Populations are grouped by superpopulation and dataset (1000 Genomes, HGDP, SGDP, GIAB).

The GIAB Ashkenazi parents appear in the European block and show a mixed profile (~78% European, ~12% South Asian, ~6% African). This does not indicate literal sub-Saharan African or South Asian ancestry. Rather, it reflects the smooth allele frequency cline across Western Eurasia, the Middle East, and North Africa. With only six components, the model has no dedicated Middle Eastern or North African cluster, so the Levantine and Near Eastern ancestry signal in Ashkenazi Jews is partitioned across the nearest available components: allele frequencies shared with populations along the Mediterranean and North African coast are absorbed by the African component, while those shared with populations along the Central and Western Asian gradient are captured by the South Asian component. This is a well-understood limitation of supervised clustering methods applied to populations that fall along continuous geographic clines rather than between discrete genetic clusters.

### Limitations of the South Asian Reference Panel

The South Asian supervised component is trained on three 1000 Genomes populations: **GIH** (Gujarati Indian, Houston), **ITU** (Indian Telugu, UK), and **STU** (Sri Lankan Tamil, UK). These 324 samples are all drawn from the Indian subcontinent's southern and western Dravidian-speaking and Indo-European-speaking groups, and two of the three are diaspora samples collected in the US and UK. This creates several limitations:

- **Geographic bias.** The reference panel has no representation from the northern tier of South Asia (e.g., Punjabi, Pashtun, Balochi, Sindhi, Bengali populations), nor from Central Asian groups that form the eastern end of the West Eurasian cline. The "South Asian" component is therefore anchored to a geographically narrow slice of the subcontinent's genetic diversity.
- **Cline truncation.** South Asian genetic variation is structured along a north–south and west–east Ancestral North Indian (ANI) to Ancestral South Indian (ASI) cline. By sampling only the middle-to-southern portion of this cline, the supervised labels train the model on a restricted allele frequency range. Populations at the ANI-heavy end of the cline (e.g., northwestern South Asians) will have their ANI-associated alleles partially absorbed by the European component, inflating European fractions and deflating South Asian fractions for these groups.
- **Diaspora sampling effects.** GIH, ITU, and STU were recruited from immigrant communities, which may not be representative of the source populations due to founder effects, selective migration, and community endogamy in the diaspora.

### Output Data Files

- [**metadata_ancestry.csv**](outputs/admixture-global-6/metadata_ancestry.csv) — Per-sample metadata augmented with the six ancestry fraction columns, max ancestry, and assigned group for all 4,324 samples in the panel.
- [**admixture_allele_freqs.tsv**](outputs/admixture-global-6/admixture_allele_freqs.tsv) — Estimated allele frequencies at each SNP for each of the six ancestry clusters (the ADMIXTURE P matrix), formatted with rsID and allele columns.

### MAF Threshold Selection

We evaluated four minor allele frequency thresholds (0.005, 0.01, 0.02, 0.05) during the ADMIXTURE QC step and observed the following:

- **Higher MAF thresholds (0.02, 0.05)** inflate European ancestry fractions in South Asian, admixed American, Middle Eastern, and African groups. At MAF 0.05, low-frequency variants that distinguish these populations from Europeans are discarded, collapsing real structure into the European component.
- **Very low MAF (0.005)** preserves more rare variation but introduces degeneracies — the model assigns excess Oceanian ancestry to European samples, likely due to noise from very rare variants or convergence artifacts.
- **MAF 0.01** provides the best balance: it retains enough low-frequency variants to separate closely related continental groups while avoiding the noise that degrades the model at MAF 0.005. This is the threshold used in the results checked into this repository.

The final QC-pruned SNP set at MAF 0.01 contains **135,020 SNPs** (after genotype missingness, MAF, long-range LD exclusion, LD pruning, sample missingness, kinship, and HWE filters).

### SNP Density Sensitivity

We also tested expanding the input SNP list from the default ~500K rsID set (`rsids_dense_chr1_22.txt`) to the SBayesRC array of ~7 million SNPs. After the same QC pipeline, the denser set yielded approximately 245K post-QC SNPs — roughly 1.8x the 135K from the default list. The resulting ancestry fractions were nearly identical: differences were negligible across all populations and samples. The population structure captured by the K=6 model is fully saturated by the LD-pruned SNPs derived from the original ~500K starter list, and increasing marker density provides no meaningful improvement in ancestry resolution. The denser set did produce slightly cleaner results in specific cases — it eliminated spurious Oceanian ancestry fractions for Mbuti and reduced South Asian noise for ACB and ASW — but these are minor refinements rather than substantive changes to the overall ancestry estimates.

## Keywords

population genetics, statistical genetics, genetic ancestry, global ancestry estimation, ancestry inference, population structure analysis, ADMIXTURE, supervised ADMIXTURE, ancestry fractions, structure plot, reference panel, 1000 Genomes, HGDP, Human Genome Diversity Project, SGDP, Simons Genome Diversity Project, GIAB, Genome in a Bottle, Ashkenazi Jewish genetics, PLINK, PLINK2, bioinformatics pipeline, reproducible genomics, SNP quality control, genotype QC, allele frequency estimation, LD pruning, linkage disequilibrium, Hardy-Weinberg equilibrium, kinship filtering, liftover, hg19 to hg38, genome build conversion, GRCh38, continental ancestry, population stratification, cross-validation, K=6, minor allele frequency, MAF filtering, human genetic diversity, genomic data harmonization, merge reference panels, ancestry estimation pipeline, open source genetics
