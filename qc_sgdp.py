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

Inputs (in QC_DIR):
  sgdp_all.{bed,bim,fam}   — SGDP data (hg19 coordinates)
  kg_qc.bim                — KG bim (hg38, used for rsID matching)

Outputs (in QC_DIR):
  sgdp_qc.{bed,bim,fam}    — QC'd SGDP data (hg38 coordinates, KG rsIDs)
"""

import os
import subprocess

import pandas as pd
from liftover import get_lifter

# ---------------------------------------------------------------------------
# Environment
# ---------------------------------------------------------------------------
QC_DIR = os.environ["QC_DIR"]
PLINK2 = os.environ["PLINK2"]
PLINK_MEMORY = os.environ["PLINK_MEMORY"]
PLINK_THREADS = os.environ["PLINK_THREADS"]

# ---------------------------------------------------------------------------
# Liftover helper
# ---------------------------------------------------------------------------
converter = get_lifter("hg19", "hg38", one_based=True)


def liftover_pos(chrom: int, pos_hg19: int):
    """Lift a single (chrom, pos) from hg19 → hg38. Returns None on failure."""
    result = converter.convert_coordinate(str(chrom), pos_hg19)
    if len(result) != 1:
        return None
    lifted_chrom, lifted_pos = result[0][0], result[0][1]
    if lifted_chrom != f"chr{chrom}":
        return None
    return lifted_pos


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

# Liftover
sgdp_bim["pos_hg38"] = [
    liftover_pos(row.chrom, row.pos_hg19) for row in sgdp_bim.itertuples()
]
sgdp_bim = sgdp_bim.dropna(subset=["pos_hg38"])
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
for name in ["sgdp_select_ids.txt", "sgdp_qc.log"]:
    path = os.path.join(QC_DIR, name)
    if os.path.exists(path):
        os.remove(path)

print("  SGDP QC complete.")
