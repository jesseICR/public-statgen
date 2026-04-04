"""
run_admixture_supervised.py — Supervised ancestry estimation with ADMIXTURE.

Performs 3-fold stratified cross-validation on supervised reference samples,
then a final run using all supervised samples as reference with projection
onto every sample in the dataset.

K model selection (env K_MODEL):
  K=3: African, American, European
  K=5: African, American, East Asian, European, South Asian
  K=6: African, American, East Asian, European, Oceanian, South Asian

The supervised.csv always contains all populations, but only populations
in the active K model are labeled in the .pop file. Samples from inactive
populations are projected (unlabeled).

Cross-validation:
  - Splits active supervised samples into N stratified folds (balanced by
    reference_population, seed from env).
  - For each fold: creates a .pop file labeling training samples with
    their reference_population, holdout and non-supervised with "-".
  - Runs ADMIXTURE in --supervised mode.

Final run:
  - All active supervised samples labeled in .pop.
  - ADMIXTURE estimates ancestry for every sample (supervised + non-supervised).

Expected environment:
  SUPERVISED_ADMIXTURE  — directory with ancestry_qc.{bed,bim,fam}
  K_MODEL               — ADMIXTURE model (3, 5, or 6)
  N_FOLDS               — number of cross-validation folds (default 3)
  ADMIXTURE_SEED        — random seed (default 42)
  PLINK_THREADS         — thread count for ADMIXTURE
  ADMIXTURE             — path to admixture binary (optional)

Input files:
  summary/supervised.csv     — sample_id, population_id, reference_population
  {SUPERVISED_ADMIXTURE}/ancestry_qc.{bed,bim,fam}

Output files (in SUPERVISED_ADMIXTURE):
  admixture_fold{1..N}.K.Q   — fold ancestry proportions
  admixture_final.K.Q         — final ancestry proportions (all samples)
  admixture_final.K.P         — allele frequency matrix (for projection)
  fold_assignments.csv        — which fold each supervised sample belongs to
"""

import os
import subprocess
import sys

import pandas as pd
import numpy as np

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPERVISED_ADMIXTURE = os.environ["SUPERVISED_ADMIXTURE"]
SCRAP = os.path.join(SUPERVISED_ADMIXTURE, "scrap")
THREADS = os.environ.get("PLINK_THREADS", "6")
ADMIXTURE = os.environ.get("ADMIXTURE", os.path.join(PROJECT_DIR, "tools", "bin", "admixture"))
N_FOLDS = int(os.environ.get("N_FOLDS", "3"))
SEED = int(os.environ.get("ADMIXTURE_SEED", "42"))

K_MODEL = int(os.environ.get("K_MODEL", "6"))
K_MODEL_POPS = {
    3: ["African", "American", "European"],
    5: ["African", "American", "East Asian", "European", "South Asian"],
    6: ["African", "American", "East Asian", "European", "Oceanian", "South Asian"],
}
active_pops = set(K_MODEL_POPS[K_MODEL])

os.makedirs(SCRAP, exist_ok=True)


def fmt(n):
    return f"{n:,}"


# ---------------------------------------------------------------------------
# 1. Load supervised labels and QC'd fam
# ---------------------------------------------------------------------------
print("=" * 60)
print("ADMIXTURE supervised ancestry estimation")
print("=" * 60)

supervised = pd.read_csv(os.path.join(PROJECT_DIR, "summary", "supervised.csv"))

# Filter to active populations for this K model
supervised_active = supervised[supervised["reference_population"].isin(active_pops)].copy()
ref_pops = sorted(active_pops)
K = len(ref_pops)

print(f"\nK model: K={K_MODEL}")
print(f"Active reference populations (K={K}):")
for rp in ref_pops:
    n = (supervised_active["reference_population"] == rp).sum()
    print(f"  {rp}: {fmt(n)}")
print(f"  Total active supervised: {fmt(len(supervised_active))}")

inactive = supervised[~supervised["reference_population"].isin(active_pops)]
if len(inactive) > 0:
    print(f"\nInactive populations (projected, not labeled):")
    for rp in sorted(inactive["reference_population"].unique()):
        n = (inactive["reference_population"] == rp).sum()
        print(f"  {rp}: {fmt(n)}")

# Map sample_id → reference_population (active only)
sid_to_ref = dict(zip(supervised_active["sample_id"], supervised_active["reference_population"]))

# Read the QC'd fam to get sample order
fam_path = os.path.join(SUPERVISED_ADMIXTURE, "ancestry_qc.fam")
fam = pd.read_csv(fam_path, sep=r"\s+", header=None,
                   names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"])
sample_order = fam["IID"].tolist()
print(f"\nQC'd samples: {fmt(len(sample_order))}")

# Which active supervised samples survived QC?
supervised_in_qc = supervised_active[supervised_active["sample_id"].isin(sample_order)].copy()
n_survived = len(supervised_in_qc)
n_lost = len(supervised_active) - n_survived
print(f"Active supervised samples in QC'd data: {fmt(n_survived)} ({fmt(n_lost)} removed by QC)")

if n_lost > 0:
    lost = supervised_active[~supervised_active["sample_id"].isin(sample_order)]
    print("  Lost by reference population:")
    for rp, grp in lost.groupby("reference_population"):
        pops = grp.groupby("population_id").size()
        for pop, cnt in pops.items():
            if cnt <= 5:
                ids = grp[grp["population_id"] == pop]["sample_id"].tolist()
                print(f"    {pop} ({rp}): {', '.join(ids)}")
            else:
                print(f"    {pop} ({rp}): {cnt} samples")

# ---------------------------------------------------------------------------
# 2. Stratified N-fold split
# ---------------------------------------------------------------------------
print(f"\nSplitting {fmt(n_survived)} active supervised samples into {N_FOLDS} stratified folds ...")

rng = np.random.RandomState(SEED)
supervised_in_qc["fold"] = -1

for rp in ref_pops:
    mask = supervised_in_qc["reference_population"] == rp
    indices = supervised_in_qc.index[mask].tolist()
    rng.shuffle(indices)
    for i, idx in enumerate(indices):
        supervised_in_qc.loc[idx, "fold"] = (i % N_FOLDS) + 1

for fold in range(1, N_FOLDS + 1):
    fold_samples = supervised_in_qc[supervised_in_qc["fold"] == fold]
    counts = fold_samples["reference_population"].value_counts().sort_index()
    print(f"  Fold {fold}: {fmt(len(fold_samples))} samples")
    for rp, n in counts.items():
        print(f"    {rp}: {fmt(n)}")

fold_path = os.path.join(SUPERVISED_ADMIXTURE, "fold_assignments.csv")
supervised_in_qc.to_csv(fold_path, index=False)
print(f"\nWrote {fold_path}")


# ---------------------------------------------------------------------------
# 3. Helper to write .pop file and run ADMIXTURE
# ---------------------------------------------------------------------------
def write_pop_file(labeled_ids, pop_path):
    """Write .pop file for ADMIXTURE supervised mode.

    labeled_ids: dict of sample_id → reference_population for reference samples.
    Each line in .pop corresponds to the same line in the .fam file.
    Reference samples get their group label; others get "-".

    Spaces in population names are replaced with underscores to avoid parsing
    issues in ADMIXTURE output.
    """
    with open(pop_path, "w") as f:
        for sid in sample_order:
            if sid in labeled_ids:
                label = labeled_ids[sid].replace(" ", "_")
                f.write(label + "\n")
            else:
                f.write("-\n")


def run_admixture(prefix, run_name):
    """Run ADMIXTURE --supervised on prefix.bed with K populations.

    ADMIXTURE outputs {basename}.{K}.Q and {basename}.{K}.P in the CWD.
    We run from SUPERVISED_ADMIXTURE so output lands there.
    """
    bed_path = os.path.relpath(f"{prefix}.bed", SUPERVISED_ADMIXTURE)
    cmd = [
        ADMIXTURE, bed_path, str(K),
        "--supervised",
        f"-j{THREADS}",
        f"--seed={SEED}",
    ]
    print(f"    Running: admixture {bed_path} {K} --supervised -j{THREADS} --seed={SEED}")
    result = subprocess.run(cmd, capture_output=True, text=True,
                            cwd=SUPERVISED_ADMIXTURE)
    if result.returncode != 0:
        print(f"    ADMIXTURE STDERR:\n{result.stderr}", file=sys.stderr)
        raise RuntimeError(f"ADMIXTURE failed with return code {result.returncode}")

    # Print select output lines
    for line in result.stdout.strip().split("\n"):
        if any(kw in line.lower() for kw in ["loglikelihood", "converge", "iteration"]):
            print(f"    {line}")

    # ADMIXTURE outputs basename.K.Q in the CWD
    basename = os.path.basename(prefix)
    q_file = os.path.join(SUPERVISED_ADMIXTURE, f"{basename}.{K}.Q")
    if os.path.exists(q_file):
        # Rename to the desired output name
        target = os.path.join(SUPERVISED_ADMIXTURE, f"{run_name}.{K}.Q")
        os.rename(q_file, target)
        # Also rename .P file
        p_file = os.path.join(SUPERVISED_ADMIXTURE, f"{basename}.{K}.P")
        if os.path.exists(p_file):
            os.rename(p_file, os.path.join(SUPERVISED_ADMIXTURE, f"{run_name}.{K}.P"))
        return target
    else:
        raise FileNotFoundError(f"Expected output {q_file} not found")


def setup_fold_files(fold_name, labeled_ids):
    """Create symlinked bed/bim/fam and .pop file in scrap/ for a fold."""
    prefix = os.path.join(SCRAP, fold_name)
    qc_prefix = os.path.join(SUPERVISED_ADMIXTURE, "ancestry_qc")

    for ext in [".bed", ".bim", ".fam"]:
        target = prefix + ext
        if os.path.exists(target):
            os.remove(target)
        os.symlink(os.path.abspath(qc_prefix + ext), target)

    write_pop_file(labeled_ids, prefix + ".pop")
    return prefix


# ---------------------------------------------------------------------------
# 4. Cross-validation runs
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print(f"{N_FOLDS}-fold cross-validation")
print("=" * 60)

labeled_all = {row["sample_id"]: row["reference_population"]
               for _, row in supervised_in_qc.iterrows()}

for fold in range(1, N_FOLDS + 1):
    print(f"\n--- Fold {fold}: holdout = fold {fold}, training = folds {[f for f in range(1, N_FOLDS+1) if f != fold]} ---")

    holdout_ids = set(
        supervised_in_qc.loc[supervised_in_qc["fold"] == fold, "sample_id"]
    )
    training_ids = {
        sid: ref for sid, ref in labeled_all.items() if sid not in holdout_ids
    }

    print(f"    Training: {fmt(len(training_ids))} samples")
    print(f"    Holdout:  {fmt(len(holdout_ids))} samples")

    prefix = setup_fold_files(f"fold{fold}", training_ids)
    q_file = run_admixture(prefix, f"admixture_fold{fold}")
    print(f"    Output: {q_file}")

# ---------------------------------------------------------------------------
# 5. Final run: all supervised labeled, project onto everyone
# ---------------------------------------------------------------------------
print("\n" + "=" * 60)
print("Final run: all active supervised as reference")
print("=" * 60)

prefix = setup_fold_files("final", labeled_all)
q_file = run_admixture(prefix, "admixture_final")
print(f"Output: {q_file}")

print("\n" + "=" * 60)
print("ADMIXTURE runs complete.")
print("=" * 60)
