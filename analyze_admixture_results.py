"""
analyze_admixture_results.py — Post-hoc analysis of ADMIXTURE supervised results.

Reads the cross-validation and final ADMIXTURE .Q/.P files, then produces:

  1. Structure plot of all supervised holdout samples (OOS predictions combined)
  2. Structure plot of non-supervised samples (full-model predictions)
  3. Per-fold cross-validation accuracy summary
  4. Augmented metadata CSV with ancestry fractions and supervised flag
  5. Formatted allele frequency file (.P) for downstream projection

Ancestry prediction logic:
  - Supervised samples (used to train the full model): their ancestry comes
    from the fold where they were HELD OUT. Each supervised sample appears in
    exactly one of the 3 folds as holdout, so we get a true out-of-sample
    estimate for every supervised sample.
  - Non-supervised samples (never used in any training): their ancestry comes
    from the FINAL run where all supervised samples are labeled. The final
    .P file is the inference artifact for projecting new data.

Expected environment:
  SUPERVISED_ADMIXTURE — directory with ADMIXTURE output files

Output files (in summary/admixture-global-K/):
  structure_holdout.png       — supervised samples, OOS ancestry estimates
  structure_projected.png     — non-supervised samples, full-model estimates
  crossval_summary.csv        — per-fold, per-population accuracy stats
  metadata_ancestry.csv       — full metadata + supervised flag + ancestry fractions
  admixture_allele_freqs.tsv  — final .P file with SNP IDs and population headers
"""

import os

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
SUPERVISED_ADMIXTURE = os.environ["SUPERVISED_ADMIXTURE"]

FLAG_THRESHOLD = 0.95
N_FOLDS = 3

PALETTE = {
    "African":      "#E69F00",
    "American":     "#D55E00",
    "East Asian":   "#009E73",
    "European":     "#0072B2",
    "Oceanian":     "#CC79A7",
    "South Asian":  "#56B4E9",
}


def fmt(n):
    return f"{n:,}"


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
supervised = pd.read_csv(os.path.join(PROJECT_DIR, "summary", "supervised.csv"))
metadata = pd.read_csv(os.path.join(PROJECT_DIR, "summary", "metadata.csv"))
folds = pd.read_csv(os.path.join(SUPERVISED_ADMIXTURE, "fold_assignments.csv"))

ref_pops = sorted(supervised["reference_population"].unique())
K = len(ref_pops)

OUTPUT_DIR = os.path.join(PROJECT_DIR, "summary", f"admixture-global-{K}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

sid_to_ref = dict(zip(supervised["sample_id"], supervised["reference_population"]))
sid_to_pop = dict(zip(metadata["sample_id"], metadata["population_id"]))

# Read fam for sample order
fam = pd.read_csv(
    os.path.join(SUPERVISED_ADMIXTURE, "ancestry_qc.fam"),
    sep=r"\s+", header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
sample_ids = fam["IID"].tolist()
sample_id_set = set(sample_ids)


# ---------------------------------------------------------------------------
# Detect ADMIXTURE column order
# ---------------------------------------------------------------------------
def detect_column_order(q_path):
    """Detect which .Q column corresponds to which reference population.

    Uses the final run where all supervised samples are labeled.  For each
    reference population, find the column with ~1.0 mean among its reference
    samples.  Returns a list of population names in column order.
    """
    q = pd.read_csv(q_path, sep=r"\s+", header=None)
    q.index = sample_ids

    col_map = [None] * K
    used_cols = set()
    for rp in ref_pops:
        rp_ids = [sid for sid in sample_ids if sid_to_ref.get(sid) == rp]
        if not rp_ids:
            continue
        means = q.loc[rp_ids].mean()
        for col_idx in means.sort_values(ascending=False).index:
            if col_idx not in used_cols:
                col_map[col_idx] = rp
                used_cols.add(col_idx)
                break
    return col_map


final_q_path = os.path.join(SUPERVISED_ADMIXTURE, f"admixture_final.{K}.Q")
col_names = detect_column_order(final_q_path)
print(f"Detected ADMIXTURE column order: {col_names}")


def load_q_file(q_path):
    """Load an ADMIXTURE .Q file with detected column names."""
    q = pd.read_csv(q_path, sep=r"\s+", header=None, names=col_names)
    q.index = sample_ids
    q.index.name = "sample_id"
    return q


# ---------------------------------------------------------------------------
# Structure plot with population_id labels
# ---------------------------------------------------------------------------
def structure_plot(q_df, title, output_path, group_col, pop_col,
                   figsize_per_sample=0.02, min_width=14, max_width=50):
    """Stacked-bar structure plot sorted by group then population.

    group_col: Series mapping sample_id -> group (for thick dividers + top labels)
    pop_col:   Series mapping sample_id -> population_id (for thin dividers + bottom labels)
    """
    ancestry_cols = [c for c in ref_pops if c in q_df.columns]
    df = q_df[ancestry_cols].copy()
    df["_group"] = group_col.values if hasattr(group_col, 'values') else df.index.map(group_col)
    df["_pop"] = pop_col.values if hasattr(pop_col, 'values') else df.index.map(pop_col)
    df = df.sort_values(["_group", "_pop"])

    groups = df["_group"].values
    pops = df["_pop"].values
    df = df.drop(columns=["_group", "_pop"])

    n = len(df)
    width = max(min_width, min(max_width, n * figsize_per_sample))
    fig, ax = plt.subplots(figsize=(width, 4.5))

    x = np.arange(n)
    bottom = np.zeros(n)
    for col in ancestry_cols:
        color = PALETTE.get(col, "#999999")
        ax.bar(x, df[col].values, bottom=bottom, width=1.0, color=color,
               edgecolor="none", label=col, align="edge")
        bottom += df[col].values

    ax.set_xlim(0, n)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Ancestry proportion")
    ax.set_title(title, fontsize=11, pad=30)
    ax.set_xticks([])

    # --- Group dividers (thick) and top labels ---
    group_starts = [0]
    for i in range(1, n):
        if groups[i] != groups[i - 1]:
            ax.axvline(i, color="black", linewidth=2, zorder=10)
            group_starts.append(i)
    group_starts.append(n)

    for j in range(len(group_starts) - 1):
        mid = (group_starts[j] + group_starts[j + 1]) / 2
        grp = groups[group_starts[j]]
        ax.text(mid, 1.02, grp, ha="center", va="bottom", fontsize=8,
                fontweight="bold", transform=ax.get_xaxis_transform())

    # --- Population dividers (thin) and bottom labels ---
    pop_starts = [0]
    for i in range(1, n):
        if pops[i] != pops[i - 1]:
            if groups[i] == groups[i - 1]:
                ax.axvline(i, color="black", linewidth=0.5, zorder=5)
            pop_starts.append(i)
    pop_starts.append(n)

    for j in range(len(pop_starts) - 1):
        mid = (pop_starts[j] + pop_starts[j + 1]) / 2
        pop_label = pops[pop_starts[j]]
        pop_n = pop_starts[j + 1] - pop_starts[j]
        if pop_n >= 2:
            ax.text(mid, -0.01, pop_label, ha="right", va="top", fontsize=6,
                    rotation=45, transform=ax.get_xaxis_transform())

    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {output_path}")


# ===================================================================
# 1. Build out-of-sample predictions for supervised samples
# ===================================================================
print("=" * 60)
print("Building out-of-sample predictions")
print("=" * 60)

# For each supervised sample, get their ancestry from the fold where
# they were held out (not used in training). This gives a true OOS estimate.
oos_frames = []
for fold in range(1, N_FOLDS + 1):
    q = load_q_file(os.path.join(SUPERVISED_ADMIXTURE, f"admixture_fold{fold}.{K}.Q"))
    holdout_sids = set(folds.loc[folds["fold"] == fold, "sample_id"])
    holdout_sup = [s for s in holdout_sids if s in sid_to_ref]
    oos_frames.append(q.loc[holdout_sup])
    print(f"  Fold {fold}: {fmt(len(holdout_sup))} holdout supervised samples")

oos = pd.concat(oos_frames)
print(f"  Total OOS predictions: {fmt(len(oos))}")

# ===================================================================
# 2. Cross-validation summary
# ===================================================================
print("\n" + "=" * 60)
print("Cross-validation accuracy")
print("=" * 60)

oos["reference_population"] = oos.index.map(sid_to_ref)
oos["expected_ancestry"] = oos.apply(
    lambda row: row.get(row["reference_population"], 0.0), axis=1
)
oos["fold"] = oos.index.map(dict(zip(folds["sample_id"], folds["fold"])))

cv_rows = []
for fold in range(1, N_FOLDS + 1):
    for rp in ref_pops:
        sub = oos[(oos["fold"] == fold) & (oos["reference_population"] == rp)]
        if len(sub) == 0:
            continue
        n_pass = (sub["expected_ancestry"] >= FLAG_THRESHOLD).sum()
        cv_rows.append({
            "fold": fold,
            "reference_population": rp,
            "n_holdout": len(sub),
            "mean_expected_ancestry": round(sub["expected_ancestry"].mean(), 4),
            "median_expected_ancestry": round(sub["expected_ancestry"].median(), 4),
            "min_expected_ancestry": round(sub["expected_ancestry"].min(), 4),
            "n_pass_95pct": n_pass,
            "pct_pass_95pct": round(n_pass / len(sub) * 100, 1),
        })

cv_summary = pd.DataFrame(cv_rows)
cv_path = os.path.join(OUTPUT_DIR, "crossval_summary.csv")
cv_summary.to_csv(cv_path, index=False)
print(f"Wrote {cv_path}\n")

for rp in ref_pops:
    rp_rows = cv_summary[cv_summary["reference_population"] == rp]
    mean_acc = rp_rows["mean_expected_ancestry"].mean()
    total_n = rp_rows["n_holdout"].sum()
    total_pass = rp_rows["n_pass_95pct"].sum()
    pct = total_pass / total_n * 100 if total_n > 0 else 0
    print(f"  {rp}: mean OOS={mean_acc:.4f}, pass>=95%: {fmt(total_pass)}/{fmt(total_n)} ({pct:.1f}%)")

# ===================================================================
# 3. Structure plot — supervised holdout (OOS)
# ===================================================================
print("\n" + "=" * 60)
print("Structure plots")
print("=" * 60)

ancestry_cols = [c for c in ref_pops if c in oos.columns]
oos_group = oos.index.map(sid_to_ref)
oos_pop = oos.index.map(sid_to_pop)

structure_plot(
    oos[ancestry_cols],
    f"ADMIXTURE K={K} — Supervised Holdout, Out-of-Sample (n={fmt(len(oos))})",
    os.path.join(OUTPUT_DIR, "structure_holdout.png"),
    group_col=oos_group,
    pop_col=oos_pop,
)

# ===================================================================
# 4. Structure plot — non-supervised (full model)
# ===================================================================
q_final = load_q_file(final_q_path)
sup_ids = set(supervised["sample_id"])
nonsup_mask = ~q_final.index.isin(sup_ids)
q_nonsup = q_final.loc[nonsup_mask]

nonsup_assigned = q_nonsup[ancestry_cols].idxmax(axis=1)
nonsup_pop = q_nonsup.index.map(sid_to_pop)

structure_plot(
    q_nonsup[ancestry_cols],
    f"ADMIXTURE K={K} — Non-Supervised Samples, Full Model (n={fmt(len(q_nonsup))})",
    os.path.join(OUTPUT_DIR, "structure_projected.png"),
    group_col=nonsup_assigned,
    pop_col=nonsup_pop,
)

# ===================================================================
# 5. Augmented metadata CSV
# ===================================================================
print("\n" + "=" * 60)
print("Building metadata with ancestry fractions")
print("=" * 60)

# OOS ancestry for supervised samples
oos_anc = oos[ancestry_cols].copy()
oos_anc.index.name = "sample_id"

# Full-model ancestry for non-supervised samples
final_nonsup = q_final.loc[nonsup_mask, ancestry_cols].copy()
final_nonsup.index.name = "sample_id"

# Combine: OOS for supervised, final for non-supervised
ancestry_combined = pd.concat([oos_anc, final_nonsup])

# Build output
meta_out = metadata.copy()
meta_out["supervised"] = meta_out["sample_id"].isin(sup_ids).map({True: "yes", False: "no"})
meta_out["in_qc_panel"] = meta_out["sample_id"].isin(sample_id_set).map({True: "yes", False: "no"})

ancestry_reset = ancestry_combined.reset_index()
meta_out = meta_out.merge(ancestry_reset, on="sample_id", how="left")

meta_out["max_ancestry"] = meta_out[ancestry_cols].max(axis=1)
meta_out["assigned_group"] = meta_out[ancestry_cols].idxmax(axis=1)

id_cols = ["sample_id", "population_id", "superpopulation", "dataset", "country",
           "neural_ancestry", "supervised", "in_qc_panel"]
tail_cols = ["max_ancestry", "assigned_group"]
col_order = [c for c in id_cols if c in meta_out.columns] + ancestry_cols + tail_cols
meta_out = meta_out[col_order]

meta_path = os.path.join(OUTPUT_DIR, "metadata_ancestry.csv")
meta_out.to_csv(meta_path, index=False)
print(f"Wrote {meta_path}")
print(f"  Total samples: {fmt(len(meta_out))}")
print(f"  Supervised (OOS predictions): {(meta_out['supervised'] == 'yes').sum()}")
print(f"  Non-supervised (full-model predictions): "
      f"{((meta_out['supervised'] == 'no') & (meta_out['in_qc_panel'] == 'yes')).sum()}")
print(f"  Not in QC panel (no predictions): "
      f"{(meta_out['in_qc_panel'] == 'no').sum()}")

# ===================================================================
# 6. Format and copy the final .P file
# ===================================================================
print("\n" + "=" * 60)
print("Formatting allele frequency file (.P)")
print("=" * 60)

# The .P file has one row per SNP, K columns of allele frequencies.
# Column order matches col_names (detected from .Q file).
# We add SNP metadata and population headers so users can project new data.
p_path = os.path.join(SUPERVISED_ADMIXTURE, f"admixture_final.{K}.P")
bim = pd.read_csv(
    os.path.join(SUPERVISED_ADMIXTURE, "ancestry_qc.bim"),
    sep="\t", header=None,
    names=["chr", "snp_id", "cm", "pos", "a1", "a2"],
)

p_raw = pd.read_csv(p_path, sep=r"\s+", header=None, names=col_names)
p_out = pd.concat([bim[["chr", "snp_id", "pos", "a1", "a2"]], p_raw], axis=1)

p_out_path = os.path.join(OUTPUT_DIR, "admixture_allele_freqs.tsv")
p_out.to_csv(p_out_path, sep="\t", index=False)
print(f"Wrote {p_out_path}")
print(f"  {fmt(len(p_out))} SNPs x {K} populations")
print(f"  Population column order: {col_names}")

print("\n" + "=" * 60)
print("Analysis complete.")
print(f"All outputs in {OUTPUT_DIR}/")
print("=" * 60)
