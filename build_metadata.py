"""
build_metadata.py — Merge KG, HGDP, and SGDP metadata into a single CSV.

Produces summary/metadata.csv with standardized superpopulation labels and
Neural ADMIXTURE single-ancestry classifications.
"""

import os

import pandas as pd


# ---------------------------------------------------------------------------
# Superpopulation mappings (dataset-specific labels → standardized names)
# ---------------------------------------------------------------------------

KG_SUPERPOPULATION_MAP = {
    "AFR": "African",
    "EUR": "European",
    "SAS": "South Asian",
    "EAS": "East Asian",
    "AMR": "American",
}

HGDP_SUPERPOPULATION_MAP = {
    "AFRICA": "African",
    "EUROPE": "European",
    "EAST_ASIA": "East Asian",
    "CENTRAL_SOUTH_ASIA": "Central South Asian",
    "MIDDLE_EAST": "Middle Eastern",
    "AMERICA": "American",
    "OCEANIA": "Oceanian",
}

SGDP_SUPERPOPULATION_MAP = {
    "Africa": "African",
    "WestEurasia": "West Eurasian",
    "SouthAsia": "South Asian",
    "CentralAsiaSiberia": "Central Asian Siberian",
    "America": "American",
    "EastAsia": "East Asian",
    "Oceania": "Oceanian",
}

# Neural ADMIXTURE pop column → ancestry label (empirically determined by
# cross-referencing Q matrix dominant components with known KG/HGDP labels).
NEURAL_ANCESTRY_LABELS = {
    "pop1": "African",
    "pop2": "American",
    "pop3": "East Asian",
    "pop4": "European",
    "pop5": "Oceanian",
    "pop6": "South Asian",
    "pop7": "West Asian",
}

FINAL_COLUMNS = ["sample_id", "population_id", "superpopulation", "dataset", "country", "neural_ancestry"]


# ---------------------------------------------------------------------------
# Metadata loaders
# ---------------------------------------------------------------------------

def load_kg_metadata():
    """Load 1000 Genomes metadata and standardize superpopulation labels."""
    df = pd.read_csv("downloads/kg_all.psam", sep="\t")
    df.rename(columns={"#IID": "sample_id", "Population": "population_id", "SuperPop": "superpop_id"}, inplace=True)
    df["superpopulation"] = df["superpop_id"].map(KG_SUPERPOPULATION_MAP)
    df["dataset"] = "kg"
    df["country"] = pd.NA
    return df[["sample_id", "population_id", "superpopulation", "dataset", "country"]]


def load_hgdp_metadata():
    """Load HGDP metadata and standardize superpopulation labels."""
    df = pd.read_csv("downloads/hgdp_all.psam", sep="\t")
    df.rename(columns={"#IID": "sample_id", "population": "population_id", "region": "region_id"}, inplace=True)
    df["superpopulation"] = df["region_id"].map(HGDP_SUPERPOPULATION_MAP)
    df["dataset"] = "hgdp"
    df["country"] = pd.NA
    return df[["sample_id", "population_id", "superpopulation", "dataset", "country"]]


def load_sgdp_metadata(kg_sample_ids, hgdp_sample_ids):
    """Load SGDP metadata, resolve ID overlaps, and standardize labels.

    The SGDP has complex ID handling: samples may appear under Illumina IDs,
    Sample IDs, or Sample ID aliases. This function also deduplicates against
    KG and HGDP samples that appear in the SGDP fam file.

    Returns:
        tuple: (metadata DataFrame, illumina_to_sample dict, sample_to_alias dict)
            The two dicts are needed downstream for Neural ADMIXTURE ID resolution.
    """
    # Load SGDP metadata (tab-separated, latin1 encoding for special characters)
    sgdp_metadata = pd.read_csv("downloads/sgdp_metadata.txt", sep="\t", encoding="latin1")
    sgdp_metadata.rename(columns={
        "Illumina_ID": "illumina_id",
        "Sample_ID": "sample_id",
        "Sample_ID(Aliases)": "sample_id_aliases",
        "SGDP_ID": "sgdp_id",
        "Population_ID": "population_id",
        "Region": "region_id",
        "Country": "country",
    }, inplace=True)

    # Load SGDP fam file (space-separated, no header)
    sgdp_fam = pd.read_csv("downloads/sgdp_all.fam", sep=" ", header=None)
    sgdp_fam.columns = ["population_id", "sample_id", "paternal_id", "maternal_id", "sex", "phenotype"]

    # Build ID-mapping dicts BEFORE any filtering (needed for Neural ADMIXTURE)
    illumina_to_sample = sgdp_metadata.set_index("illumina_id")["sample_id"].to_dict()
    sample_to_alias = sgdp_metadata.set_index("sample_id")["sample_id_aliases"].to_dict()

    # Remove samples that overlap with KG or HGDP
    kg_hgdp_ids = set(kg_sample_ids) | set(hgdp_sample_ids)
    sgdp_overlap_ids = kg_hgdp_ids.intersection(sgdp_fam["sample_id"])

    sgdp_metadata = sgdp_metadata[~sgdp_metadata["sample_id_aliases"].isin(sgdp_overlap_ids)]
    sgdp_metadata.reset_index(drop=True, inplace=True)
    sgdp_fam = sgdp_fam[~sgdp_fam["sample_id"].isin(sgdp_overlap_ids)]
    sgdp_fam.reset_index(drop=True, inplace=True)

    # Drop constant-value fam columns (paternal_id, maternal_id, phenotype)
    for col in ["paternal_id", "maternal_id", "phenotype"]:
        assert sgdp_fam[col].nunique(dropna=False) == 1
        sgdp_fam.drop(col, axis=1, inplace=True)

    # Swap sample_id and sample_id_aliases so aliases become the canonical IDs
    sgdp_metadata.rename(columns={
        "sample_id": "sample_id_old",
        "sample_id_aliases": "sample_id",
    }, inplace=True)

    # Add fam-only entries (samples in the fam file but not in the metadata)
    sgdp_fam_only_entries = sgdp_fam[~sgdp_fam["sample_id"].isin(sgdp_metadata["sample_id"])]
    sgdp_metadata = pd.concat([sgdp_metadata, sgdp_fam_only_entries])

    # Remove any remaining overlaps with HGDP
    sgdp_metadata = sgdp_metadata[~sgdp_metadata["sample_id"].isin(hgdp_sample_ids)]
    sgdp_metadata.reset_index(drop=True, inplace=True)

    # Standardize superpopulation labels
    sgdp_metadata["superpopulation"] = sgdp_metadata["region_id"].map(SGDP_SUPERPOPULATION_MAP)
    sgdp_metadata["dataset"] = "sgdp"

    result = sgdp_metadata[["sample_id", "population_id", "superpopulation", "dataset", "country"]].copy()
    return result, illumina_to_sample, sample_to_alias


def load_neural_ancestry(illumina_to_sample, sample_to_alias, merged_or_related_ids):
    """Load Neural ADMIXTURE Q matrix and assign single-ancestry labels.

    Each sample's dominant ancestry component (highest Q value) is mapped to
    a named population label using NEURAL_ANCESTRY_LABELS.

    The ID resolution follows a two-step process:
      1. Illumina IDs → Sample IDs
      2. Sample IDs → aliases (only for IDs not already in merged_or_related_ids)
    """
    # Load the Q matrix (7 ancestry components per sample)
    q_matrix = pd.read_csv("downloads/neural/allchms/data/GT_train.Q", sep=" ", header=None)
    pop_columns = [f"pop{i}" for i in range(1, q_matrix.shape[1] + 1)]
    q_matrix.columns = pop_columns

    # Load the neural fam file to get sample IDs
    neural_fam = pd.read_csv("downloads/neural/allchms/data/whole_filt_ld_single.fam", sep="\t", header=None)
    original_neural_ids = neural_fam.iloc[:, 1].tolist()

    # Two-step ID resolution
    resolved_ids = [illumina_to_sample.get(sid, sid) for sid in original_neural_ids]
    resolved_ids = [
        sample_to_alias.get(sid, sid) if sid not in merged_or_related_ids else sid
        for sid in resolved_ids
    ]
    assert len(original_neural_ids) == len(resolved_ids)

    # Assign dominant ancestry label via argmax across pop columns
    dominant_pop = q_matrix[pop_columns].idxmax(axis=1)
    neural_ancestry = dominant_pop.map(NEURAL_ANCESTRY_LABELS)

    result = pd.DataFrame({"sample_id": resolved_ids, "neural_ancestry": neural_ancestry})
    result.drop_duplicates(subset="sample_id", inplace=True)
    return result


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    # Load the three metadata sources
    kg_metadata = load_kg_metadata()
    hgdp_metadata = load_hgdp_metadata()
    sgdp_metadata, illumina_to_sample, sample_to_alias = load_sgdp_metadata(
        kg_sample_ids=kg_metadata["sample_id"].tolist(),
        hgdp_sample_ids=hgdp_metadata["sample_id"].tolist(),
    )

    # Build merged + related ID set (needed for Neural ADMIXTURE ID resolution)
    merged_fam = pd.read_csv("merge/merged_kg_hgdp_sgdp.fam", sep="\t", header=None)
    merged_ids = merged_fam.iloc[:, 1].tolist()
    kg_related = pd.read_csv("downloads/deg2_hg38.king.cutoff.out.id", sep="\t")
    kg_related_ids = kg_related.iloc[:, 0].tolist()
    merged_or_related_ids = set(merged_ids + kg_related_ids)

    # Load Neural ADMIXTURE ancestry classifications
    neural_ancestry = load_neural_ancestry(illumina_to_sample, sample_to_alias, merged_or_related_ids)

    # Merge all three datasets vertically, then join neural ancestry
    metadata = pd.concat([kg_metadata, hgdp_metadata, sgdp_metadata], ignore_index=True)
    metadata = metadata.merge(neural_ancestry, on="sample_id", how="left")

    # Impute missing superpopulations from other samples with the same population_id
    pop_to_superpop = (
        metadata.dropna(subset="superpopulation")
        .groupby("population_id")["superpopulation"]
        .unique()
    )
    for pop_id, superpops in pop_to_superpop.items():
        if len(superpops) == 1:
            mask = (metadata["population_id"] == pop_id) & metadata["superpopulation"].isna()
            metadata.loc[mask, "superpopulation"] = superpops[0]

    assert list(metadata.columns) == FINAL_COLUMNS
    assert set(metadata["sample_id"]) == merged_or_related_ids
    assert metadata["superpopulation"].notna().all(), "Some samples still missing superpopulation"

    # Write output
    os.makedirs("summary", exist_ok=True)
    metadata.to_csv("summary/metadata.csv", index=False)
    print(f"Wrote {len(metadata)} samples to summary/metadata.csv")
    print(f"  Datasets: {metadata['dataset'].value_counts().to_dict()}")
    print(f"  Neural ancestry coverage: {metadata['neural_ancestry'].notna().sum()} / {len(metadata)}")

    # Print breakdown of each neural ancestry group
    neural_samples = metadata.dropna(subset="neural_ancestry")
    for ancestry, group in sorted(neural_samples.groupby("neural_ancestry")):
        print(f"\n{'=' * 60}")
        print(f"Neural Ancestry: {ancestry} (n={len(group)})")
        print(f"{'=' * 60}")

        print(f"\n  By dataset:")
        for val, count in group["dataset"].value_counts().items():
            print(f"    {val}: {count}")

        print(f"\n  By superpopulation:")
        for val, count in group["superpopulation"].value_counts().items():
            print(f"    {val}: {count}")

        print(f"\n  By population_id:")
        for val, count in group["population_id"].value_counts().items():
            print(f"    {val}: {count}")

        countries = group["country"].dropna()
        if len(countries) > 0:
            print(f"\n  By country:")
            for val, count in countries.value_counts().items():
                print(f"    {val}: {count}")

        mismatches = group[group["neural_ancestry"] != group["superpopulation"]]
        if len(mismatches) > 0:
            print(f"\n  Superpopulation mismatches ({len(mismatches)} samples):")
            for (pop_id, superpop), count in (
                mismatches.groupby(["population_id", "superpopulation"]).size()
                .sort_values(ascending=False).items()
            ):
                print(f"    {pop_id} ({superpop}): {count}")


if __name__ == "__main__":
    main()
