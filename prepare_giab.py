"""
prepare_giab.py — Convert GIAB Ashkenazi parent VCFs to PLINK format.

Reads the GIAB benchmark VCFs for HG003 (father) and HG004 (mother),
fills in homozygous-reference genotypes at reference panel positions within
high-confidence callable regions, and converts to PLINK bed/bim/fam.

The benchmark VCFs only contain variant calls.  Positions in the panel that
fall within GIAB high-confidence regions but are NOT variant are homozygous
reference — this script fills those in so the merged panel retains all SNPs.

Expected environment:
  PLINK2        — path to plink2 binary
  MERGE_DIR     — directory containing merged_kg_hgdp_sgdp.bim (panel SNPs)
  QC_DIR        — output directory for giab_qc.{bed,bim,fam}
  DOWNLOADS_DIR — directory with GIAB downloads (VCFs + BEDs)
  PLINK_MEMORY  — memory limit in MB
  PLINK_THREADS — number of threads

Outputs (in QC_DIR):
  giab_qc.{bed,bim,fam}
"""

import gzip
import os
import subprocess
import sys
from bisect import bisect_right


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PLINK2 = os.environ["PLINK2"]
MERGE_DIR = os.environ["MERGE_DIR"]
QC_DIR = os.environ["QC_DIR"]
DOWNLOADS_DIR = os.environ["DOWNLOADS_DIR"]
PLINK_MEMORY = os.environ["PLINK_MEMORY"]
PLINK_THREADS = os.environ["PLINK_THREADS"]

GIAB_DIR = os.path.join(DOWNLOADS_DIR, "giab")

SAMPLES = ["HG003", "HG004"]
SAMPLE_FILES = {
    "HG003": {
        "vcf": os.path.join(GIAB_DIR, "HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"),
        "bed": os.path.join(GIAB_DIR, "HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"),
    },
    "HG004": {
        "vcf": os.path.join(GIAB_DIR, "HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz"),
        "bed": os.path.join(GIAB_DIR, "HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed"),
    },
}


# ---------------------------------------------------------------------------
# Load reference panel positions from bim
# ---------------------------------------------------------------------------
def load_panel_snps(bim_path):
    """Return (dict of (chr,pos)->(rsid,a1,a2), ordered list of keys)."""
    snps = {}
    ordered = []
    with open(bim_path) as f:
        for line in f:
            parts = line.split()
            chr_num = int(parts[0])
            rsid = parts[1]
            pos = int(parts[3])
            a1 = parts[4]   # ALT (plink2 convention)
            a2 = parts[5]   # REF (plink2 convention)
            key = (chr_num, pos)
            if key not in snps:
                snps[key] = (rsid, a1, a2)
                ordered.append(key)
    return snps, ordered


# ---------------------------------------------------------------------------
# Load GIAB high-confidence BED intervals
# ---------------------------------------------------------------------------
def load_callable_regions(bed_path):
    """Return dict of chr_num -> sorted list of (start, end) tuples (0-based)."""
    regions = {}
    with open(bed_path) as f:
        for line in f:
            parts = line.split("\t")
            chrom = parts[0]
            if not chrom.startswith("chr"):
                continue
            chr_str = chrom[3:]
            if not chr_str.isdigit():
                continue
            chr_num = int(chr_str)
            start = int(parts[1])
            end = int(parts[2])
            regions.setdefault(chr_num, []).append((start, end))
    for v in regions.values():
        v.sort()
    return regions


def is_callable(intervals, pos):
    """Check if a 1-based VCF position falls within a callable BED interval."""
    if not intervals:
        return False
    pos_0 = pos - 1
    idx = bisect_right(intervals, (pos_0, float("inf"))) - 1
    if idx < 0:
        return False
    start, end = intervals[idx]
    return start <= pos_0 < end


# ---------------------------------------------------------------------------
# Parse GIAB VCF for biallelic SNPs at panel positions
# ---------------------------------------------------------------------------
def load_giab_variants(vcf_path, panel_snps):
    """Return dict of (chr,pos) -> (ref, alt, gt_string)."""
    variants = {}
    with gzip.open(vcf_path, "rt") as f:
        for line in f:
            if line[0] == "#":
                continue
            parts = line.split("\t", 10)
            chrom = parts[0]
            if not chrom.startswith("chr"):
                continue
            chr_str = chrom[3:]
            if not chr_str.isdigit():
                continue
            chr_num = int(chr_str)
            pos = int(parts[1])
            key = (chr_num, pos)
            if key not in panel_snps:
                continue
            ref = parts[3]
            alt = parts[4]
            if len(ref) != 1 or len(alt) != 1 or "," in alt:
                continue
            fmt_fields = parts[8].split(":")
            gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else 0
            gt = parts[9].split(":")[gt_idx].strip().replace("|", "/")
            if gt in ("0/0", "0/1", "1/0", "1/1"):
                variants[key] = (ref, alt, gt)
    return variants


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    bim_path = os.path.join(MERGE_DIR, "merged_kg_hgdp_sgdp.bim")
    print("Loading reference panel SNPs ...")
    panel_snps, ordered_keys = load_panel_snps(bim_path)
    print(f"  Panel SNPs: {len(panel_snps):,}")

    # Load per-sample data
    sample_data = {}
    for sid in SAMPLES:
        paths = SAMPLE_FILES[sid]
        print(f"\nProcessing {sid} ...")
        regions = load_callable_regions(paths["bed"])
        n_intervals = sum(len(v) for v in regions.values())
        print(f"  Callable intervals: {n_intervals:,}")
        variants = load_giab_variants(paths["vcf"], panel_snps)
        print(f"  Panel-matching variants: {len(variants):,}")
        sample_data[sid] = {"regions": regions, "variants": variants}

    # Write combined VCF at all panel positions
    vcf_path = os.path.join(QC_DIR, "giab_combined.vcf")
    print(f"\nWriting combined VCF → {vcf_path}")

    stats = {"variant": 0, "hom_ref": 0, "missing": 0, "allele_mismatch": 0}

    with open(vcf_path, "w") as f:
        f.write("##fileformat=VCFv4.2\n")
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
        header = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        header += SAMPLES
        f.write("\t".join(header) + "\n")

        for key in ordered_keys:
            chr_num, pos = key
            rsid, a1, a2 = panel_snps[key]  # a1=ALT, a2=REF in plink2

            # Collect variant calls at this position
            var_info = {}
            for sid in SAMPLES:
                if key in sample_data[sid]["variants"]:
                    var_info[sid] = sample_data[sid]["variants"][key]

            # Determine allele orientation
            flip_needed = False
            skip = False
            if var_info:
                vcf_ref, vcf_alt, _ = next(iter(var_info.values()))
                if {vcf_ref, vcf_alt} != {a1, a2}:
                    stats["allele_mismatch"] += 1
                    skip = True
                elif vcf_ref != a2:
                    flip_needed = True

            if skip:
                continue

            # Build genotypes
            gts = []
            for sid in SAMPLES:
                if sid in var_info:
                    _, _, gt = var_info[sid]
                    if flip_needed:
                        flip_map = {"0": "1", "1": "0", ".": "."}
                        gt = "/".join(flip_map[x] for x in gt.split("/"))
                    gts.append(gt)
                    stats["variant"] += 1
                elif is_callable(sample_data[sid]["regions"].get(chr_num, []), pos):
                    gts.append("0/0")
                    stats["hom_ref"] += 1
                else:
                    gts.append("./.")
                    stats["missing"] += 1

            # REF=a2 (panel REF), ALT=a1 (panel ALT)
            record = [f"chr{chr_num}", str(pos), rsid, a2, a1, ".", ".", ".", "GT"]
            record += gts
            f.write("\t".join(record) + "\n")

    total_gt = stats["variant"] + stats["hom_ref"] + stats["missing"]
    n_snps = total_gt // len(SAMPLES) if SAMPLES else 0
    print(f"\n  SNPs written: {n_snps:,}")
    print(f"  Per-sample genotype counts:")
    print(f"    Variant:          {stats['variant']:,}")
    print(f"    Hom-ref (filled): {stats['hom_ref']:,}")
    print(f"    Missing:          {stats['missing']:,}")
    print(f"  Allele mismatches (skipped): {stats['allele_mismatch']:,}")
    if total_gt > 0:
        print(f"  Missing rate: {stats['missing'] / total_gt * 100:.1f}%")

    # Convert VCF to PLINK bed/bim/fam
    print("\nConverting to PLINK format with plink2 ...")
    out_prefix = os.path.join(QC_DIR, "giab_qc")
    cmd = [
        PLINK2,
        "--vcf", vcf_path,
        "--chr", "1-22",
        "--snps-only", "just-acgt",
        "--max-alleles", "2",
        "--make-bed",
        "--out", out_prefix,
        "--memory", PLINK_MEMORY,
        "--threads", PLINK_THREADS,
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"plink2 stderr:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"plink2 failed with return code {result.returncode}")

    bim_count = sum(1 for _ in open(out_prefix + ".bim"))
    fam_count = sum(1 for _ in open(out_prefix + ".fam"))
    print(f"  Output: {bim_count:,} SNPs, {fam_count} samples")

    os.remove(vcf_path)
    print("\nGIAB preparation complete.")


if __name__ == "__main__":
    main()
