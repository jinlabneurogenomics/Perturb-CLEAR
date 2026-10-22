"""Summarize soft morphology composition within mapped RNA subpopulations."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


RNA_SUBPOPULATIONS = ["T0", "T1", "T2", "T3"]
M_WEIGHTS = [f"morphology_weight_M{index}" for index in range(1, 6)]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize morphology composition within mapped RNA subpopulations.")
    parser.add_argument("output_dir", type=Path)
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("--mapping-workbook", type=Path)
    source.add_argument("--rna-scores", type=Path)
    parser.add_argument("--morphology-scores", type=Path)
    parser.add_argument("--datasets", default="Sorensen2026")
    args = parser.parse_args()
    if args.rna_scores is not None and args.morphology_scores is None:
        parser.error("--morphology-scores is required with --rna-scores")
    if args.mapping_workbook is not None and args.morphology_scores is not None:
        parser.error("Do not combine --mapping-workbook with --morphology-scores")
    return args


def canonical_dataset(value: object) -> str:
    text = str(value).lower()
    if "sorensen" in text:
        return "Sorensen2026"
    if "scala" in text:
        return "Scala2021"
    return str(value)


def prepare_rna_scores(rna: pd.DataFrame) -> pd.DataFrame:
    """Return one row per neuron with T0 to T3 global scores."""
    if {"RNA_subpopulation", "probability", "source_cell_id"}.issubset(rna.columns):
        rna = rna.copy()
        rna["RNA_subpopulation"] = (
            rna["RNA_subpopulation"].astype(str).str.replace("STref C", "T", regex=False)
        )
        rna["dataset"] = rna["dataset"].map(canonical_dataset)
        rna["source_cell_id"] = rna["source_cell_id"].astype(str)
        wide = rna.pivot_table(
            index=["dataset", "source_cell_id"],
            columns="RNA_subpopulation",
            values="probability",
            aggfunc="first",
        ).reset_index()
    elif "neuron_ID" in rna.columns:
        required = ["dataset", "neuron_ID"] + [
            f"probability_{subpopulation}_all4" for subpopulation in RNA_SUBPOPULATIONS
        ]
        missing = [column for column in required if column not in rna.columns]
        if missing:
            raise ValueError(f"Wide RNA mapping table is missing columns: {missing}")
        wide = rna[required].copy()
        wide["dataset"] = wide["dataset"].map(canonical_dataset)
        wide["source_cell_id"] = wide.pop("neuron_ID").astype(str)
        wide = wide.rename(
            columns={
                f"probability_{subpopulation}_all4": subpopulation
                for subpopulation in RNA_SUBPOPULATIONS
            }
        )
    else:
        raise ValueError(
            "RNA scores must be either a long probability table or the released wide mapping table"
        )
    for subpopulation in RNA_SUBPOPULATIONS:
        if subpopulation not in wide.columns:
            wide[subpopulation] = 0.0
    return wide


def prepare_morphology_scores(morphology: pd.DataFrame) -> pd.DataFrame:
    """Standardize column names from released or newly generated mapping tables."""
    morphology = morphology.copy()
    aliases = {
        "neuron_ID": "source_cell_id",
        "cell_id": "source_cell_id",
        "animal_ID": "animal_id",
        "published_broad_class": "broad_cell_type",
        "outside_reference_range": "outside_reference_global_0_05",
    }
    for old, new in aliases.items():
        if old in morphology.columns and new not in morphology.columns:
            morphology = morphology.rename(columns={old: new})
    for index in range(1, 6):
        old = f"weight_M{index}"
        new = f"morphology_weight_M{index}"
        if old in morphology.columns and new not in morphology.columns:
            morphology = morphology.rename(columns={old: new})
    morphology["dataset"] = morphology["dataset"].map(canonical_dataset)
    morphology["source_cell_id"] = morphology["source_cell_id"].astype(str)
    return morphology


def load_score_tables(args: argparse.Namespace) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Load Table 11 or fresh RNA and morphology mapping outputs."""
    if args.mapping_workbook is None:
        return (
            pd.read_csv(args.rna_scores, low_memory=False),
            pd.read_csv(args.morphology_scores, low_memory=False),
        )

    metadata = pd.read_excel(
        args.mapping_workbook, sheet_name="Published neurons", header=0
    )
    rna = pd.read_excel(args.mapping_workbook, sheet_name="RNA mapping", header=0)
    morphology = pd.read_excel(
        args.mapping_workbook, sheet_name="Morphology mapping", header=0
    )
    required_metadata = ["neuron_ID", "source"]
    missing_metadata = [column for column in required_metadata if column not in metadata.columns]
    if missing_metadata:
        raise ValueError(f"Published-neuron metadata is missing columns: {missing_metadata}")
    lookup = metadata[required_metadata].drop_duplicates("neuron_ID").copy()
    lookup["dataset"] = lookup["source"].map(canonical_dataset)
    lookup = lookup.drop(columns="source")
    rna = rna.merge(lookup, on="neuron_ID", how="left", validate="one_to_one")
    morphology = morphology.merge(
        lookup, on="neuron_ID", how="left", validate="one_to_one"
    )
    if rna["dataset"].isna().any() or morphology["dataset"].isna().any():
        raise ValueError("Every mapping row must match the Published neurons worksheet")
    if "OOD_p_global" in morphology.columns:
        morphology["outside_reference_global_0_05"] = (
            pd.to_numeric(morphology["OOD_p_global"], errors="coerce") < 0.05
        )
    return rna, morphology


def main() -> None:
    args = parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    datasets = [item.strip() for item in args.datasets.split(",") if item.strip()]

    rna_raw, morphology_raw = load_score_tables(args)
    wide = prepare_rna_scores(rna_raw)
    morphology = prepare_morphology_scores(morphology_raw)
    required = ["dataset", "source_cell_id", "animal_id", "broad_cell_type"] + M_WEIGHTS
    missing = [column for column in required if column not in morphology.columns]
    if missing:
        raise ValueError(f"Morphology mapping is missing columns: {missing}")

    merged = morphology.merge(wide, on=["dataset", "source_cell_id"], how="inner", validate="one_to_one")
    merged = merged.loc[merged["dataset"].isin(datasets)].copy()
    if "outside_reference_global_0_05" in merged.columns:
        merged = merged.loc[~merged["outside_reference_global_0_05"].fillna(True)].copy()

    layer_45 = merged["broad_cell_type"].astype(str).eq("L4/5 IT")
    layer_23 = merged["broad_cell_type"].astype(str).eq("L2/3 IT")
    merged["global_top_subpopulation"] = merged[RNA_SUBPOPULATIONS].idxmax(axis=1)
    merged["rna_layer_concordant"] = (
        (layer_45 & merged["global_top_subpopulation"].isin(["T0", "T1"]))
        | (layer_23 & merged["global_top_subpopulation"].isin(["T2", "T3"]))
    )
    merged = merged.loc[merged["rna_layer_concordant"]].copy()
    layer_45 = merged["broad_cell_type"].astype(str).eq("L4/5 IT")
    layer_23 = merged["broad_cell_type"].astype(str).eq("L2/3 IT")
    merged["within_layer_mass"] = np.where(
        layer_45, merged["T0"] + merged["T1"],
        np.where(layer_23, merged["T2"] + merged["T3"], np.nan),
    )
    for subpopulation in RNA_SUBPOPULATIONS:
        merged[f"conditional_{subpopulation}"] = merged[subpopulation] / merged["within_layer_mass"].clip(lower=np.finfo(float).eps)
    merged["RNA_subpopulation"] = np.where(
        layer_45,
        np.where(merged["conditional_T0"] >= merged["conditional_T1"], "T0", "T1"),
        np.where(
            layer_23,
            np.where(merged["conditional_T2"] >= merged["conditional_T3"], "T2", "T3"),
            "unresolved",
        ),
    )
    merged = merged.loc[merged["RNA_subpopulation"].isin(RNA_SUBPOPULATIONS)].copy()

    cell_columns = [
        "dataset", "source_cell_id", "animal_id", "broad_cell_type", "RNA_subpopulation",
        "within_layer_mass", "global_top_subpopulation", "rna_layer_concordant",
    ] + [f"conditional_{subpopulation}" for subpopulation in RNA_SUBPOPULATIONS] + M_WEIGHTS
    merged[cell_columns].to_csv(args.output_dir / "cell_level_rna_morphology_mapping.csv.gz", index=False)

    cell_summary = (
        merged.groupby(["dataset", "broad_cell_type", "RNA_subpopulation"], observed=True)
        .agg(
            n_cells=("source_cell_id", "nunique"),
            n_animals=("animal_id", "nunique"),
            **{f"mean_{column}": (column, "mean") for column in M_WEIGHTS},
        )
        .reset_index()
    )
    cell_summary.to_csv(args.output_dir / "rna_subpopulation_morphology_composition.csv", index=False)

    animal = (
        merged.groupby(["dataset", "animal_id", "broad_cell_type", "RNA_subpopulation"], observed=True)[M_WEIGHTS]
        .mean()
        .reset_index()
    )
    animal.to_csv(args.output_dir / "animal_level_morphology_composition.csv", index=False)

    pd.DataFrame(
        {
            "parameter": [
                "datasets",
                "rna_subpopulations",
                "morphology_weights",
                "require_layer_concordance",
                "exclude_outside_reference_global_0_05",
                "retained_cells",
                "retained_animals",
            ],
            "value": [
                ";".join(datasets),
                ";".join(RNA_SUBPOPULATIONS),
                ";".join(M_WEIGHTS),
                True,
                "outside_reference_global_0_05" in morphology.columns,
                merged["source_cell_id"].nunique(),
                merged["animal_id"].nunique(),
            ],
        }
    ).to_csv(args.output_dir / "run_parameters.csv", index=False)
    print(f"Composition summaries written to {args.output_dir}")


if __name__ == "__main__":
    main()
