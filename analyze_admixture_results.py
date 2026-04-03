"""
analyze_admixture_results.py — Post-hoc analysis of ADMIXTURE supervised results.

Reads the cross-validation and final ADMIXTURE .Q/.P files, then produces:

  1. Structure plot of all supervised holdout samples (OOS predictions combined)
  2. Structure plot of projected samples (full-model predictions), organized
     by metadata superpopulation with within-population ancestry sorting
  3. Augmented metadata CSV with ancestry fractions and supervised flag
  4. Formatted allele frequency file (.P) for downstream projection

Ancestry prediction logic:
  - Active supervised samples (used to train the model): their ancestry comes
    from the fold where they were HELD OUT. Each sample appears in exactly one
    fold as holdout, yielding a true out-of-sample estimate.
  - All other samples (non-supervised + inactive supervised): their ancestry
    comes from the FINAL run where all active supervised samples are labeled.

Expected environment:
  SUPERVISED_ADMIXTURE — directory with ADMIXTURE output files
  K_MODEL              — ADMIXTURE model (3, 5, or 6)
  FLAG_THRESHOLD       — ancestry fraction pass threshold (default 0.95)
  N_FOLDS              — cross-validation folds (default 3)

Output files (in summary/admixture-global-K/):
  structure_holdout.png       — supervised samples, OOS ancestry estimates
  structure_projected.png     — projected samples, full-model estimates
  metadata_ancestry.csv       — full metadata + supervised flag + ancestry fractions
  admixture_allele_freqs.tsv  — final .P file with SNP IDs and population headers
"""

import os
import re

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

FLAG_THRESHOLD = float(os.environ.get("FLAG_THRESHOLD", "0.95"))
N_FOLDS = int(os.environ.get("N_FOLDS", "3"))

K_MODEL = int(os.environ.get("K_MODEL", "6"))
K_MODEL_POPS = {
    3: ["African", "American", "European"],
    5: ["African", "American", "East Asian", "European", "South Asian"],
    6: ["African", "American", "East Asian", "European", "Oceanian", "South Asian"],
}
active_pops = set(K_MODEL_POPS[K_MODEL])

PALETTE = {
    "African":      "#E69F00",
    "American":     "#D55E00",
    "East Asian":   "#009E73",
    "European":     "#0072B2",
    "Oceanian":     "#CC79A7",
    "South Asian":  "#56B4E9",
}

# Geographic display order for superpopulation grouping (roughly out-of-Africa)
SUPERPOP_DISPLAY_ORDER = [
    "African",
    "Middle Eastern",
    "European",
    "West Eurasian",
    "South Asian",
    "Central South Asian",
    "Central Asian Siberian",
    "East Asian",
    "Oceanian",
    "American",
]

DATASET_ABBREV = {
    "kg": "1000G",
    "hgdp": "HGDP",
    "sgdp": "SGDP",
    "giab": "GIAB",
}
DATASET_FULLNAME = {
    "1000G": "1000 Genomes Project",
    "HGDP": "Human Genome Diversity Project",
    "SGDP": "Simons Genome Diversity Project",
    "GIAB": "Genome in a Bottle",
}
DATASET_ORDER = ["kg", "hgdp", "sgdp", "giab"]


def fmt(n):
    return f"{n:,}"


def superpop_sort_key(sp):
    """Return sort key for superpopulation display order."""
    try:
        return SUPERPOP_DISPLAY_ORDER.index(sp)
    except ValueError:
        return len(SUPERPOP_DISPLAY_ORDER)


# ---------------------------------------------------------------------------
# Load data
# ---------------------------------------------------------------------------
supervised = pd.read_csv(os.path.join(PROJECT_DIR, "summary", "supervised.csv"))
metadata = pd.read_csv(os.path.join(PROJECT_DIR, "summary", "metadata.csv"))
folds = pd.read_csv(os.path.join(SUPERVISED_ADMIXTURE, "fold_assignments.csv"))

ref_pops = sorted(active_pops)
K = len(ref_pops)

OUTPUT_DIR = os.path.join(PROJECT_DIR, "summary", f"admixture-global-{K}")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Active supervised: only populations used for labeling in this K model
supervised_active = supervised[supervised["reference_population"].isin(active_pops)]
active_sup_ids = set(supervised_active["sample_id"])
all_sup_ids = set(supervised["sample_id"])

sid_to_ref = dict(zip(supervised_active["sample_id"], supervised_active["reference_population"]))
sid_to_pop = dict(zip(metadata["sample_id"], metadata["population_id"]))
sid_to_superpop = dict(zip(metadata["sample_id"], metadata["superpopulation"]))
sid_to_dataset = dict(zip(metadata["sample_id"], metadata["dataset"]))

# Read fam for sample order
fam = pd.read_csv(
    os.path.join(SUPERVISED_ADMIXTURE, "ancestry_qc.fam"),
    sep=r"\s+", header=None,
    names=["FID", "IID", "PAT", "MAT", "SEX", "PHENO"],
)
sample_ids = fam["IID"].tolist()
sample_id_set = set(sample_ids)

# ---------------------------------------------------------------------------
# Load 1KG population descriptions for formatted labels
# ---------------------------------------------------------------------------
pop_desc_path = os.path.join(PROJECT_DIR, "downloads", "kg_population_names.tsv")
pop_code_to_name = {}
if os.path.exists(pop_desc_path):
    try:
        pop_desc_df = pd.read_csv(pop_desc_path, sep="\t")
        if "Population Code" in pop_desc_df.columns and "Population Description" in pop_desc_df.columns:
            pop_code_to_name = dict(zip(
                pop_desc_df["Population Code"],
                pop_desc_df["Population Description"],
            ))
            print(f"Loaded {len(pop_code_to_name)} population descriptions from {pop_desc_path}")
    except Exception as e:
        print(f"Warning: could not load population descriptions: {e}")


def humanize_label(text):
    """Make a population label more readable.

    Replaces underscores with spaces and inserts spaces before CamelCase
    transitions (e.g. 'BergamoItalian' -> 'Bergamo Italian',
    'PapuanHighlands' -> 'Papuan Highlands').
    """
    text = text.replace("_", " ")
    text = re.sub(r'([a-z])([A-Z])', r'\1 \2', text)
    return text


def format_pop_label(pop_id):
    """Format a population label with description if available."""
    if pop_id in pop_code_to_name:
        return f"{pop_code_to_name[pop_id]} ({pop_id})"
    return humanize_label(pop_id)


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
# Structure plot with population_id labels and within-population sorting
# ---------------------------------------------------------------------------
def structure_plot(q_df, title, output_path, group_col, pop_col,
                   group_order=None, dataset_col=None, format_labels=False,
                   figsize_per_sample=0.02, min_width=14, max_width=50):
    """Stacked-bar structure plot with hierarchical sorting and variable-width bars.

    group_col:    Series mapping sample_id -> group (thick dividers + top labels)
    pop_col:      Series mapping sample_id -> population_id (thin dividers + bottom labels)
    group_order:  optional list of group names in display order
    dataset_col:  optional Series mapping sample_id -> dataset (medium dividers + mid labels)
    format_labels: if True, use population descriptions for bottom labels

    Small populations get wider bars for visibility: n<3 → 3x width, n<8 → 2x.
    All dividers use a white-gap underlay so colored bars don't bleed through.
    All population labels are rendered using two stagger levels with tick marks.
    """
    ancestry_cols = [c for c in ref_pops if c in q_df.columns]
    df = q_df[ancestry_cols].copy()
    df["_group"] = group_col.values if hasattr(group_col, "values") else df.index.map(group_col)
    df["_pop"] = pop_col.values if hasattr(pop_col, "values") else df.index.map(pop_col)
    has_dataset = dataset_col is not None
    if has_dataset:
        df["_dataset"] = dataset_col.values if hasattr(dataset_col, "values") else df.index.map(dataset_col)

    if group_order is None:
        group_order = sorted(df["_group"].dropna().unique())

    def _sort_pop_block(block):
        means = block[ancestry_cols].mean().sort_values(ascending=False)
        return block.sort_values(by=means.index.tolist(),
                                 ascending=[False] * len(ancestry_cols))

    # Build sorted dataframe: group → dataset → population → ancestry sort
    sorted_parts = []
    for grp in group_order:
        grp_df = df[df["_group"] == grp]
        if grp_df.empty:
            continue
        if has_dataset:
            ds_order = [d for d in DATASET_ORDER if d in grp_df["_dataset"].values]
            ds_order += sorted(set(grp_df["_dataset"].unique()) - set(ds_order))
            for ds in ds_order:
                ds_df = grp_df[grp_df["_dataset"] == ds]
                for pop in sorted(ds_df["_pop"].unique()):
                    block = ds_df[ds_df["_pop"] == pop].copy()
                    sorted_parts.append(_sort_pop_block(block))
        else:
            for pop in sorted(grp_df["_pop"].unique()):
                block = grp_df[grp_df["_pop"] == pop].copy()
                sorted_parts.append(_sort_pop_block(block))

    missing = df[df["_group"].isna()]
    if not missing.empty:
        sorted_parts.append(missing)

    df = pd.concat(sorted_parts)
    groups = df["_group"].values
    pops = df["_pop"].values
    datasets = df["_dataset"].values if has_dataset else None
    drop_cols = ["_group", "_pop"] + (["_dataset"] if has_dataset else [])
    df = df.drop(columns=drop_cols)

    n = len(df)

    # --- Compute variable widths per individual ---
    # Count population sizes and assign wider bars to small populations
    pop_sizes = pd.Series(pops).groupby(pops).transform("size").values
    ind_widths = np.where(pop_sizes < 3, 3.0, np.where(pop_sizes < 8, 2.0, 1.0))
    x_pos = np.concatenate([[0.0], np.cumsum(ind_widths[:-1])])
    total_x = x_pos[-1] + ind_widths[-1] if n > 0 else 0

    fig_width = max(min_width, min(max_width, total_x * figsize_per_sample))
    fig, ax = plt.subplots(figsize=(fig_width, 5))

    bottom = np.zeros(n)
    for col in ancestry_cols:
        color = PALETTE.get(col, "#999999")
        ax.bar(x_pos, df[col].values, bottom=bottom, width=ind_widths,
               color=color, edgecolor="none", label=col, align="edge")
        bottom += df[col].values

    ax.set_xlim(0, total_x)
    ax.set_ylim(0, 1)
    ax.set_ylabel("Ancestry proportion")
    ax.set_title(title, fontsize=11, pad=45 if has_dataset else 30)
    ax.set_xticks([])

    def _block_mid(start_idx, end_idx):
        """Midpoint in x-space for a block of bars [start_idx, end_idx)."""
        left = x_pos[start_idx]
        right = x_pos[end_idx - 1] + ind_widths[end_idx - 1]
        return (left + right) / 2

    def _divider_x(idx):
        """X coordinate for a divider at the left edge of bar idx."""
        return x_pos[idx]

    # --- Group dividers (white gap + black line) and top labels ---
    group_starts = [0]
    for i in range(1, n):
        if groups[i] != groups[i - 1]:
            dx = _divider_x(i)
            ax.axvline(dx, color="white", linewidth=4, zorder=9)
            ax.axvline(dx, color="black", linewidth=1.5, zorder=10)
            group_starts.append(i)
    group_starts.append(n)

    top_y = 1.07 if has_dataset else 1.02
    for j in range(len(group_starts) - 1):
        mid = _block_mid(group_starts[j], group_starts[j + 1])
        grp = groups[group_starts[j]]
        ax.text(mid, top_y, humanize_label(str(grp)), ha="center", va="bottom",
                fontsize=7, fontweight="bold",
                transform=ax.get_xaxis_transform())

    # --- Dataset dividers (white gap + dashed gray) and labels ---
    if datasets is not None:
        ds_starts = [0]
        for i in range(1, n):
            if datasets[i] != datasets[i - 1]:
                if groups[i] == groups[i - 1]:
                    dx = _divider_x(i)
                    ax.axvline(dx, color="white", linewidth=2.5, zorder=6)
                    ax.axvline(dx, color="#444444", linewidth=0.8,
                               linestyle=(0, (4, 2)), zorder=7)
                ds_starts.append(i)
        ds_starts.append(n)

        for j in range(len(ds_starts) - 1):
            mid = _block_mid(ds_starts[j], ds_starts[j + 1])
            ds = datasets[ds_starts[j]]
            ds_label = DATASET_ABBREV.get(ds, ds)
            ax.text(mid, 1.015, ds_label, ha="center", va="bottom",
                    fontsize=5, fontstyle="italic", color="#333333",
                    transform=ax.get_xaxis_transform())

    # --- Population dividers (white gap + thin line) and bottom labels ---
    pop_starts = [0]
    for i in range(1, n):
        if pops[i] != pops[i - 1]:
            same_group = groups[i] == groups[i - 1]
            same_dataset = datasets is None or datasets[i] == datasets[i - 1]
            if same_group and same_dataset:
                dx = _divider_x(i)
                ax.axvline(dx, color="white", linewidth=1.5, zorder=4)
                ax.axvline(dx, color="black", linewidth=0.3, zorder=5)
            pop_starts.append(i)
    pop_starts.append(n)

    # Stagger labels at two y-levels with tick marks to avoid overlap
    trans = ax.get_xaxis_transform()
    stagger_y = [-0.03, -0.14]

    for j in range(len(pop_starts) - 1):
        mid = _block_mid(pop_starts[j], pop_starts[j + 1])
        pop_label = pops[pop_starts[j]]
        display_label = format_pop_label(pop_label) if format_labels else humanize_label(pop_label)

        level = j % 2
        label_y = stagger_y[level]

        # Tick mark from bar bottom to label
        ax.plot([mid, mid], [-0.005, label_y + 0.008], transform=trans,
                color="#aaaaaa", linewidth=0.3, clip_on=False, zorder=1)
        # Label text
        ax.text(mid, label_y, display_label, ha="right", va="top", fontsize=4,
                rotation=90, transform=trans, clip_on=False)

    # --- Legends ---
    ax.legend(loc="upper left", bbox_to_anchor=(1.01, 1), fontsize=8, frameon=False)

    # Dataset key (if datasets are shown)
    if has_dataset:
        key_lines = [f"{abbr} = {DATASET_FULLNAME[abbr]}"
                     for abbr in DATASET_ABBREV.values() if abbr in DATASET_FULLNAME]
        key_text = "\n".join(key_lines)
        ax.text(1.01, 0.38, key_text, fontsize=5.5, ha="left", va="top",
                transform=ax.transAxes, linespacing=1.5,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#f8f8f8",
                          edgecolor="#cccccc", linewidth=0.5))

    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved {output_path}")


# ===================================================================
# 1. Build out-of-sample predictions for supervised samples
# ===================================================================
print("=" * 60)
print("Building out-of-sample predictions")
print("=" * 60)

# For each active supervised sample, get their ancestry from the fold where
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
# 2. Structure plot — supervised holdout (OOS)
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
    group_order=ref_pops,
    format_labels=True,
)

# ===================================================================
# 3. Structure plot — projected samples (full model)
# ===================================================================
q_final = load_q_file(final_q_path)

# Projected = everything NOT in the active supervised holdout set
# (includes non-supervised + inactive supervised populations)
oos_sids = set(oos.index)
projected_mask = ~q_final.index.isin(oos_sids)
q_projected = q_final.loc[projected_mask]

# Group by metadata superpopulation (not majority ancestry)
projected_group = q_projected.index.map(sid_to_superpop)
projected_pop = q_projected.index.map(sid_to_pop)
projected_dataset = q_projected.index.map(sid_to_dataset)

# Determine display order based on superpopulations present
present_superpops = sorted(projected_group.dropna().unique(), key=superpop_sort_key)

structure_plot(
    q_projected[ancestry_cols],
    f"ADMIXTURE K={K} — Projected Ancestry, Full Model (n={fmt(len(q_projected))})",
    os.path.join(OUTPUT_DIR, "structure_projected.png"),
    group_col=projected_group,
    pop_col=projected_pop,
    group_order=present_superpops,
    dataset_col=projected_dataset,
    format_labels=True,
)

# ===================================================================
# 4. Augmented metadata CSV
# ===================================================================
print("\n" + "=" * 60)
print("Building metadata with ancestry fractions")
print("=" * 60)

# OOS ancestry for active supervised samples
oos_anc = oos[ancestry_cols].copy()
oos_anc.index.name = "sample_id"

# Full-model ancestry for all other samples in the QC panel
final_rest = q_final.loc[projected_mask, ancestry_cols].copy()
final_rest.index.name = "sample_id"

# Combine: OOS for active supervised, final for everyone else
ancestry_combined = pd.concat([oos_anc, final_rest])

# Build output
meta_out = metadata.copy()
meta_out["supervised"] = meta_out["sample_id"].isin(all_sup_ids).map({True: "yes", False: "no"})
meta_out["in_qc_panel"] = meta_out["sample_id"].isin(sample_id_set).map({True: "yes", False: "no"})

ancestry_reset = ancestry_combined.reset_index()
meta_out = meta_out.merge(ancestry_reset, on="sample_id", how="left")

meta_out["max_ancestry"] = meta_out[ancestry_cols].max(axis=1)
has_ancestry = meta_out[ancestry_cols].notna().any(axis=1)
meta_out.loc[has_ancestry, "assigned_group"] = meta_out.loc[has_ancestry, ancestry_cols].idxmax(axis=1)

id_cols = ["sample_id", "population_id", "superpopulation", "dataset", "country",
           "neural_ancestry", "supervised", "in_qc_panel"]
tail_cols = ["max_ancestry", "assigned_group"]
col_order = [c for c in id_cols if c in meta_out.columns] + ancestry_cols + tail_cols
meta_out = meta_out[col_order]

meta_path = os.path.join(OUTPUT_DIR, "metadata_ancestry.csv")
meta_out.to_csv(meta_path, index=False)
print(f"Wrote {meta_path}")
print(f"  Total samples: {fmt(len(meta_out))}")
print(f"  Supervised (OOS predictions): {len(oos)}")
print(f"  Projected (full-model predictions): "
      f"{(meta_out['in_qc_panel'] == 'yes').sum() - len(oos)}")
print(f"  Not in QC panel (no predictions): "
      f"{(meta_out['in_qc_panel'] == 'no').sum()}")

# ===================================================================
# 5. Format and copy the final .P file
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
