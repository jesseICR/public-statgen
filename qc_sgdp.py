#!/usr/bin/env python3
"""
qc_sgdp.py
----------
QC pipeline for SGDP data:
  1. Lift SGDP variant positions from hg19 → hg38
  2. Match to KG variants on (chrom, pos, alleles) to assign rsIDs
  3. Extract matched variants and rewrite bim with hg38 positions + rsIDs

Expected environment variables:
  QC_DIR         — directory containing sgdp_all.{bed,bim,fam} and kg_qc.bim
  PLINK2         — path to plink2 binary
  PLINK_MEMORY   — memory limit in MB for PLINK
  PLINK_THREADS  — number of threads for PLINK
  LIFTOVER       — path to UCSC liftOver binary
  CHAIN_FILE     — path to hg19ToHg38.over.chain.gz

Inputs (in QC_DIR):
  sgdp_all.{bed,bim,fam}   — SGDP data (hg19 coordinates)
  kg_qc.bim                — KG bim (hg38, used for rsID matching)

Outputs (in QC_DIR):
  sgdp_qc.{bed,bim,fam}    — QC'd SGDP data (hg38 coordinates, KG rsIDs)
"""

import os
import subprocess

import pandas as pd

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
QC_DIR = os.environ["QC_DIR"]
PLINK2 = os.environ["PLINK2"]
PLINK_MEMORY = os.environ["PLINK_MEMORY"]
PLINK_THREADS = os.environ["PLINK_THREADS"]
LIFTOVER = os.environ["LIFTOVER"]
CHAIN_FILE = os.environ["CHAIN_FILE"]


# ---------------------------------------------------------------------------
# 1. Read SGDP bim and lift positions to hg38
# ---------------------------------------------------------------------------
print("  Reading SGDP bim and lifting hg19 → hg38 ...")
sgdp_bim = pd.read_csv(
    os.path.join(QC_DIR, "sgdp_all.bim"),
    sep="\t", header=None,
    names=["chrom", "rsid", "cm", "pos_hg19", "a1", "a2"],
)

# Keep autosomes only
sgdp_bim = sgdp_bim[sgdp_bim["chrom"].isin(range(1, 23))].copy()

# Write temporary BED for UCSC liftOver (BIM is 1-based; BED is 0-based half-open)
bed_in = os.path.join(QC_DIR, "_liftover_input.bed")
bed_out = os.path.join(QC_DIR, "_liftover_mapped.bed")
bed_unmapped = os.path.join(QC_DIR, "_liftover_unmapped.bed")

sgdp_bim = sgdp_bim.reset_index(drop=True)
sgdp_bim["idx"] = sgdp_bim.index
bed_df = pd.DataFrame({
    "chrom": "chr" + sgdp_bim["chrom"].astype(str),
    "start": (sgdp_bim["pos_hg19"] - 1).astype(int),
    "end": sgdp_bim["pos_hg19"].astype(int),
    "idx": sgdp_bim["idx"].astype(int),
})
bed_df.to_csv(bed_in, sep="\t", header=False, index=False)

subprocess.run([LIFTOVER, bed_in, CHAIN_FILE, bed_out, bed_unmapped], check=True)

# Parse mapped output and convert 0-based start back to 1-based position
mapped = pd.read_csv(
    bed_out, sep="\t", header=None,
    names=["chrom_hg38", "start", "end", "idx"],
)
mapped["pos_hg38"] = mapped["start"] + 1
# Drop variants that mapped to alt contigs (e.g. chr7_KI270803v1_alt)
mapped = mapped[mapped["chrom_hg38"].str.match(r"^chr\d+$")].copy()
mapped["chrom_hg38_num"] = mapped["chrom_hg38"].str.replace("chr", "").astype(int)

# Merge back on index and drop variants that changed chromosome
sgdp_bim = sgdp_bim.merge(mapped[["idx", "pos_hg38", "chrom_hg38_num"]], on="idx", how="inner")
sgdp_bim = sgdp_bim[sgdp_bim["chrom"] == sgdp_bim["chrom_hg38_num"]].copy()
sgdp_bim.drop(columns=["idx", "chrom_hg38_num"], inplace=True)
sgdp_bim["pos_hg38"] = sgdp_bim["pos_hg38"].astype(int)

# ---------------------------------------------------------------------------
# 2. Match to KG on (chrom, pos, alleles) to get rsIDs
# ---------------------------------------------------------------------------
print("  Matching SGDP variants to KG by (chrom, pos, alleles) ...")
kg_bim = pd.read_csv(
    os.path.join(QC_DIR, "kg_qc.bim"),
    sep="\t", header=None,
    names=["chrom", "rsid", "cm", "pos", "a1", "a2"],
)

matched = sgdp_bim.merge(
    kg_bim,
    left_on=["chrom", "pos_hg38"], right_on=["chrom", "pos"],
    suffixes=("_sgdp", "_kg"), how="inner",
)

# Keep only allele-concordant variants (allowing strand flips)
alleles_match = (
    ((matched["a1_sgdp"] == matched["a1_kg"]) & (matched["a2_sgdp"] == matched["a2_kg"]))
    | ((matched["a1_sgdp"] == matched["a2_kg"]) & (matched["a2_sgdp"] == matched["a1_kg"]))
)
matched = matched[alleles_match].copy()
print(f"  Matched {len(matched):,} variants")

# ---------------------------------------------------------------------------
# 3. Extract matched variants with PLINK2
# ---------------------------------------------------------------------------
print("  Extracting matched variants with PLINK2 ...")
extract_ids_path = os.path.join(QC_DIR, "sgdp_select_ids.txt")
matched["rsid_sgdp"].to_csv(extract_ids_path, index=False, header=False)

subprocess.run(
    [
        PLINK2,
        "--bfile", os.path.join(QC_DIR, "sgdp_all"),
        "--extract", extract_ids_path,
        "--make-bed",
        "--out", os.path.join(QC_DIR, "sgdp_qc"),
        "--memory", PLINK_MEMORY,
        "--threads", PLINK_THREADS,
    ],
    check=True,
)

# ---------------------------------------------------------------------------
# 4. Update bim with hg38 positions and KG rsIDs
# ---------------------------------------------------------------------------
print("  Updating sgdp_qc.bim with hg38 positions and rsIDs ...")
rsid_map = matched.set_index("rsid_sgdp")["rsid_kg"].to_dict()
pos_map = matched.set_index("rsid_sgdp")["pos"].to_dict()

qc_bim = pd.read_csv(
    os.path.join(QC_DIR, "sgdp_qc.bim"),
    sep="\t", header=None,
    names=["chrom", "rsid", "cm", "pos", "a1", "a2"],
)
qc_bim["pos"] = qc_bim["rsid"].map(pos_map)
qc_bim["rsid"] = qc_bim["rsid"].map(rsid_map)
qc_bim.to_csv(
    os.path.join(QC_DIR, "sgdp_qc.bim"),
    sep="\t", header=False, index=False,
)

# ---------------------------------------------------------------------------
# 5. Clean up intermediates
# ---------------------------------------------------------------------------
for name in [
    "sgdp_select_ids.txt", "sgdp_qc.log",
    "_liftover_input.bed", "_liftover_mapped.bed", "_liftover_unmapped.bed",
]:
    path = os.path.join(QC_DIR, name)
    if os.path.exists(path):
        os.remove(path)

print("  SGDP QC complete.")
