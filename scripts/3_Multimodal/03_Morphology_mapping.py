"""Map published neurons to fixed Perturb-CLEAR morphology profiles.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import minimize_scalar
from scipy.spatial.distance import cdist


M_ORDER = ["M1", "M2", "M3", "M4", "M5"]
REFERENCE_DATASET = "Perturb_CLEAR"
SEED = 20260831
META_COLUMNS = ["cell_id", "dataset", "broad_cell_type", "animal_id", "official_M"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Map published neurons to the fixed Perturb-CLEAR morphology reference."
    )
    parser.add_argument("reference_features", type=Path, help="Perturb-CLEAR reference-neuron feature CSV")
    parser.add_argument("published_features", type=Path, help="Published-neuron feature CSV")
    parser.add_argument("feature_manifest", type=Path, help="CSV listing the 36 morphology features")
    parser.add_argument("combat_script", type=Path, help="Layer-protected ComBat R script")
    parser.add_argument("output_dir", type=Path, help="Output directory")
    parser.add_argument("--rscript", default="Rscript", help="Rscript executable")
    return parser.parse_args()


def normalize_m_label(value: object) -> str:
    text = str(value).strip().replace("M_", "M")
    return text if text in M_ORDER else ""


def row_mean_fill(frame: pd.DataFrame) -> pd.DataFrame:
    numeric = frame.astype(float).copy()
    row_means = numeric.mean(axis=1, skipna=True).fillna(0.0)
    return numeric.T.fillna(row_means).T


def reference_z(data: pd.DataFrame, features: list[str]) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    raw = data[features].apply(pd.to_numeric, errors="coerce")
    reference = data["dataset"].eq(REFERENCE_DATASET)
    means = raw.loc[reference].mean(axis=0, skipna=True)
    standard_deviations = raw.loc[reference].std(axis=0, ddof=1, skipna=True)
    invalid = (~np.isfinite(standard_deviations)) | (standard_deviations <= 0)
    if invalid.any():
        raise ValueError(f"Non-usable reference features: {standard_deviations.index[invalid].tolist()}")
    return row_mean_fill((raw - means) / standard_deviations), means, standard_deviations


def rms_distances(matrix: np.ndarray, profiles: np.ndarray) -> np.ndarray:
    return cdist(matrix, profiles, metric="euclidean") / math.sqrt(matrix.shape[1])


def soft_weights(distances: np.ndarray, temperature: float) -> np.ndarray:
    logits = -(distances**2) / max(float(temperature), 1e-12)
    logits -= logits.max(axis=1, keepdims=True)
    values = np.exp(logits)
    return values / values.sum(axis=1, keepdims=True)


def multiclass_log_loss(truth: np.ndarray, probabilities: np.ndarray) -> float:
    chosen = probabilities[np.arange(len(truth)), truth]
    return float(-np.mean(np.log(np.clip(chosen, 1e-15, 1.0))))


def calibrate_reference(
    reference_raw: pd.DataFrame, features: list[str]
) -> tuple[float, np.ndarray, dict[str, np.ndarray], pd.DataFrame]:
    distance_blocks: list[np.ndarray] = []
    truths: list[str] = []
    calibration_rows: list[dict[str, object]] = []
    layer_calibration: dict[str, list[float]] = {}
    for animal in sorted(reference_raw["animal_id"].astype(str).unique()):
        train = reference_raw.loc[~reference_raw["animal_id"].astype(str).eq(animal)].copy()
        test = reference_raw.loc[reference_raw["animal_id"].astype(str).eq(animal)].copy()
        train_raw = train[features].apply(pd.to_numeric, errors="coerce")
        test_raw = test[features].apply(pd.to_numeric, errors="coerce")
        means = train_raw.mean(axis=0)
        standard_deviations = train_raw.std(axis=0, ddof=1)
        if (~np.isfinite(standard_deviations) | (standard_deviations <= 0)).any():
            raise ValueError(f"Invalid fold scaling while holding out animal {animal}")
        train_z = row_mean_fill((train_raw - means) / standard_deviations)
        test_z = row_mean_fill((test_raw - means) / standard_deviations)
        profiles = (
            train_z.assign(official_M=train["official_M"].to_numpy())
            .groupby("official_M", observed=True)
            .mean()
            .reindex(M_ORDER)
        )
        if profiles.isna().any().any():
            raise ValueError(f"A morphology type is absent while holding out animal {animal}")
        distances = rms_distances(test_z.to_numpy(), profiles.to_numpy())
        distance_blocks.append(distances)
        truths.extend(test["official_M"].tolist())
        nearest = distances.min(axis=1)
        for cell_id, layer, distance in zip(test["cell_id"], test["broad_cell_type"], nearest):
            layer_calibration.setdefault(str(layer), []).append(float(distance))
            calibration_rows.append(
                {
                    "held_out_animal": animal,
                    "cell_id": cell_id,
                    "broad_cell_type": layer,
                    "nearest_distance": float(distance),
                }
            )
    stacked = np.vstack(distance_blocks)
    truth_index = np.asarray([M_ORDER.index(value) for value in truths], dtype=int)
    fit = minimize_scalar(
        lambda log_temperature: multiclass_log_loss(
            truth_index, soft_weights(stacked, math.exp(log_temperature))
        ),
        bounds=(-8, 6),
        method="bounded",
    )
    return (
        math.exp(float(fit.x)),
        stacked.min(axis=1),
        {key: np.asarray(value) for key, value in layer_calibration.items()},
        pd.DataFrame(calibration_rows),
    )


def empirical_upper_tail(values: np.ndarray, calibration: np.ndarray) -> np.ndarray:
    return np.asarray(
        [(1 + np.sum(calibration >= value)) / (len(calibration) + 1) for value in values],
        dtype=float,
    )


def main() -> None:
    args = parse_args()
    np.random.seed(SEED)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    table_dir = args.output_dir / "tables"
    intermediate_dir = args.output_dir / "intermediate"
    table_dir.mkdir(exist_ok=True)
    intermediate_dir.mkdir(exist_ok=True)

    reference_data = pd.read_csv(args.reference_features, low_memory=False)
    published_data = pd.read_csv(args.published_features, low_memory=False)
    if not reference_data.columns.equals(published_data.columns):
        missing_from_published = sorted(set(reference_data.columns) - set(published_data.columns))
        missing_from_reference = sorted(set(published_data.columns) - set(reference_data.columns))
        if missing_from_published or missing_from_reference:
            raise ValueError(
                "Reference and published-neuron tables must contain the same columns. "
                f"Missing from published: {missing_from_published}; "
                f"missing from reference: {missing_from_reference}"
            )
        published_data = published_data.loc[:, reference_data.columns]
    if "dataset" not in reference_data.columns or "dataset" not in published_data.columns:
        raise ValueError("Both input tables need a dataset column")
    if not reference_data["dataset"].eq(REFERENCE_DATASET).all():
        raise ValueError("The reference table may contain only Perturb_CLEAR neurons")
    if published_data["dataset"].eq(REFERENCE_DATASET).any():
        raise ValueError("The published-neuron table must not contain Perturb_CLEAR neurons")
    data = pd.concat([reference_data, published_data], ignore_index=True, sort=False)
    missing_metadata = [column for column in META_COLUMNS if column not in data.columns]
    if missing_metadata:
        raise ValueError(f"Missing metadata columns: {missing_metadata}")
    data["cell_id"] = data["cell_id"].astype(str)
    data["animal_id"] = data["animal_id"].astype(str)
    data["official_M"] = data["official_M"].map(normalize_m_label)
    if data["cell_id"].duplicated().any():
        raise ValueError("cell_id must be unique")

    feature_manifest = pd.read_csv(args.feature_manifest)
    feature_column = next(
        (column for column in ["feature", "metric", "feature_name", "Metrics"] if column in feature_manifest.columns),
        None,
    )
    if feature_column is None:
        raise ValueError("The feature manifest needs a feature, metric, feature_name, or Metrics column")
    features = feature_manifest[feature_column].dropna().astype(str).tolist()
    if len(features) != 36 or len(set(features)) != 36:
        raise ValueError("Exactly 36 unique morphology features are required")
    missing_features = [feature for feature in features if feature not in data.columns]
    if missing_features:
        raise ValueError(f"Missing morphology features: {missing_features}")

    reference = data["dataset"].eq(REFERENCE_DATASET)
    if reference.sum() == 0:
        raise ValueError("The Perturb_CLEAR reference batch is absent")
    if (data.loc[reference, "official_M"] == "").any():
        raise ValueError("Every Perturb_CLEAR reference neuron needs an M1-M5 label")
    if set(data.loc[reference, "official_M"]) != set(M_ORDER):
        raise ValueError("The reference must contain all five morphology types")
    if "condition" in data.columns:
        condition = data.loc[reference, "condition"].astype(str).str.lower()
        accepted = condition.isin(["ctrl", "control", "st2", "st_2"])
        if not accepted.all():
            raise ValueError("Only control Perturb_CLEAR neurons may define the morphology reference")

    standardized, means, standard_deviations = reference_z(data, features)
    combat_input = pd.concat(
        [data[META_COLUMNS].reset_index(drop=True), standardized.reset_index(drop=True)], axis=1
    )
    combat_input_path = intermediate_dir / "reference_standardized_features.csv"
    combat_output_path = intermediate_dir / "layer_protected_combat_features.csv"
    combat_input.to_csv(combat_input_path, index=False)
    subprocess.run(
        [args.rscript, str(args.combat_script), str(combat_input_path), str(combat_output_path)],
        check=True,
    )
    corrected = pd.read_csv(combat_output_path, low_memory=False)
    corrected_matrix = corrected[features].astype(float)

    reference_raw = data.loc[reference, META_COLUMNS + features].copy()
    temperature, global_calibration, layer_calibration, calibration_detail = calibrate_reference(
        reference_raw, features
    )
    profiles = (
        corrected_matrix.loc[reference]
        .assign(official_M=data.loc[reference, "official_M"].to_numpy())
        .groupby("official_M", observed=True)
        .mean()
        .reindex(M_ORDER)
    )
    distances = rms_distances(corrected_matrix.to_numpy(), profiles.to_numpy())
    weights = soft_weights(distances, temperature)
    order = np.argsort(distances, axis=1)
    nearest = distances[np.arange(len(data)), order[:, 0]]
    second = distances[np.arange(len(data)), order[:, 1]]
    global_p = empirical_upper_tail(nearest, global_calibration)
    layer_p = np.asarray(
        [
            empirical_upper_tail(np.asarray([distance]), layer_calibration[str(layer)])[0]
            if str(layer) in layer_calibration
            else np.nan
            for distance, layer in zip(nearest, data["broad_cell_type"])
        ]
    )

    result = data[META_COLUMNS].copy()
    result["mapping_role"] = np.where(reference, "fixed_reference", "published_query")
    result["nearest_M_profile_diagnostic"] = [M_ORDER[index] for index in order[:, 0]]
    result["nearest_distance"] = nearest
    result["second_distance"] = second
    result["distance_margin"] = second - nearest
    result["maximum_relative_similarity"] = weights.max(axis=1)
    result["normalized_entropy"] = -(
        weights * np.log(np.clip(weights, 1e-15, 1.0))
    ).sum(axis=1) / math.log(len(M_ORDER))
    result["ood_p_global"] = global_p
    result["ood_p_layer"] = layer_p
    result["outside_reference_global_0_05"] = global_p < 0.05
    for index, label in enumerate(M_ORDER):
        result[f"distance_{label}"] = distances[:, index]
        result[f"morphology_weight_{label}"] = weights[:, index]
    result.to_csv(table_dir / "morphology_mapping_scores.csv.gz", index=False)

    profile_table = profiles.reset_index().rename(columns={"official_M": "morphology_type"})
    profile_table.to_csv(table_dir / "fixed_morphology_profiles.csv", index=False)
    calibration_detail["temperature"] = temperature
    calibration_detail.to_csv(table_dir / "leave_one_animal_out_calibration.csv", index=False)
    pd.DataFrame(
        {
            "feature": features,
            "reference_mean": means.reindex(features).to_numpy(),
            "reference_standard_deviation": standard_deviations.reindex(features).to_numpy(),
        }
    ).to_csv(table_dir / "reference_feature_scaling.csv", index=False)

    (table_dir / "run_parameters.json").write_text(
        json.dumps(
            {
                "seed": SEED,
                "reference_dataset": REFERENCE_DATASET,
                "morphology_types": M_ORDER,
                "n_features": len(features),
                "temperature": temperature,
                "outside_reference_p_threshold": 0.05,
                "reference_cells": int(reference.sum()),
                "query_cells": int((~reference).sum()),
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    print(f"Morphology mapping complete: {table_dir}")


if __name__ == "__main__":
    main()
