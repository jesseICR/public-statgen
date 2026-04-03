"""
build_supervised.py — Assign samples to K=6 supervised ADMIXTURE reference populations.

Rules derived from union of Sharma (All of Us) and Marino (creatinine) reference panels:
  African:     Union of both, minus LWK
  American:    Admix3way subset of PEL, Colombian, Maya, Pima, Karitiana, Surui + SGDP: Mixe, Quechua, Zapotec, Chane, Nahua
  East Asian:  Union of both, minus Cambodian and KHV, plus CDX and CHS from 1KGP
  European:    Union of both
  Oceanian:    All Sharma
  South Asian: GIH, ITU, STU

Outlier exclusions:
  - Koenig harmonized outliers (summary/koenig_harmonized_outliers_2024.txt)
  - Other Spanish outliers (literature_reference/other_spanish_outliers.txt)
  - American admixed outliers (literature_reference/american_admixed_outliers.txt)
  - Oceanian admixed outliers (literature_reference/oceanian_admixed_outliers.txt)
"""

import pandas as pd

REFERENCE_POPULATIONS = {
    "African": ["ESN", "MSL", "YRI", "Mbuti", "Biaka"],
    "American": ["PEL", "Colombian", "Maya", "Pima", "Karitiana", "Surui", "Mixe", "Quechua", "Zapotec", "Chane", "Nahua"],
    "East Asian": ["CHB", "CHS", "CDX", "Dai", "JPT", "She", "Tujia"],
    "European": ["FIN", "GBR", "IBS", "TSI", "Tuscan", "French", "Basque", "BergamoItalian"],
    "Oceanian": ["Bougainville", "PapuanHighlands", "PapuanSepik"],
    "South Asian": ["GIH", "ITU", "STU"],
}

OUTLIER_FILES = {
    "Koenig harmonized outliers": "summary/koenig_harmonized_outliers_2024.txt",
    "Other Spanish outliers": "literature_reference/other_spanish_outliers.txt",
    "American admixed outliers": "literature_reference/american_admixed_outliers.txt",
    "Oceanian admixed outliers": "literature_reference/oceanian_admixed_outliers.txt",
}

# Invert to population_id -> reference_population
POP_TO_REF = {}
for ref_pop, pop_ids in REFERENCE_POPULATIONS.items():
    for pid in pop_ids:
        POP_TO_REF[pid] = ref_pop


def load_outlier_ids():
    """Load all outlier sample IDs from exclusion lists."""
    all_outliers = set()
    for label, path in OUTLIER_FILES.items():
        with open(path, "r") as f:
            ids = {line.strip() for line in f if line.strip()}
        all_outliers |= ids
    return all_outliers


def main():
    metadata = pd.read_csv("summary/metadata.csv")
    outlier_ids = load_outlier_ids()

    # Filter to reference populations only
    mask = metadata["population_id"].isin(POP_TO_REF)
    supervised = metadata.loc[mask, ["sample_id", "population_id"]].copy()
    supervised["reference_population"] = supervised["population_id"].map(POP_TO_REF)

    # Exclude outliers and report each one
    outlier_mask = supervised["sample_id"].isin(outlier_ids)
    if outlier_mask.any():
        excluded = supervised[outlier_mask]
        print("Outlier exclusions:")
        for label, path in OUTLIER_FILES.items():
            with open(path, "r") as f:
                ids = {line.strip() for line in f if line.strip()}
            hits = excluded[excluded["sample_id"].isin(ids)]
            if len(hits) > 0:
                pop_counts = hits.groupby(["population_id", "reference_population"]).size()
                print(f"  {label} ({path}): {len(hits)} samples")
                for (pop, ref), count in pop_counts.items():
                    if count <= 5:
                        sids = hits[(hits["population_id"] == pop)]["sample_id"].tolist()
                        print(f"    {', '.join(sids)} — {pop} ({ref})")
                    else:
                        print(f"    {count} samples — {pop} ({ref})")
        print()
        supervised = supervised[~outlier_mask].copy()

    # Verify no unmapped samples
    assert supervised["reference_population"].notna().all()

    supervised.to_csv("summary/supervised.csv", index=False)
    print(f"Wrote {len(supervised)} samples to summary/supervised.csv\n")

    # Print population_ids in each reference population
    for ref_pop, pop_ids in sorted(REFERENCE_POPULATIONS.items()):
        counts = supervised[supervised["reference_population"] == ref_pop]["population_id"].value_counts()
        total = counts.sum()
        print(f"{ref_pop} (n={total}):")
        for pid in sorted(pop_ids):
            n = counts.get(pid, 0)
            print(f"  {pid}: {n}")
        print()


if __name__ == "__main__":
    main()
