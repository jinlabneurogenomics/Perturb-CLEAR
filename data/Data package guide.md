# Data package guide

The `data` directory contains project-owned morphology data and compact inputs for the multimodal workflows.

## Morphology

- `Perturb-CLEAR screen raw metrics.xlsx`: pooled-screen morphometric measurements.
- `Adnp validation raw metrics.xlsx`: independent Adnp validation-cohort measurements.
- `Perturb-CLEAR screen and Adnp validation sample summary.xlsx`: sample metadata and counts.
- `MetricsSet1.csv`: the prespecified 36-feature morphometric panel.

## Multimodal

- `Published multimodally profiled neuron mapping.xlsx`: Similar to supplementary Table 11, containing metadata and the submitted RNA-subpopulation and morphology-type mapping scores for Sorensen and Scala neurons.
- `Published neuron metadata.csv`: identifiers, source datasets, published broad classes, mapping availability, and analytical roles for 205 published neurons.
- `Published neuron morphometrics.csv`: the 36 morphology features for 91 Sorensen, 42 Scala, and 72 Gouwens neurons. Gouwens neurons contribute only to normalization.
- `Perturb-CLEAR reference neurons morphometrics.csv`: the same 36 features for all 162 Perturb-CLEAR control neurons used to define the fixed morphology reference.
- `Published neuron morphology separation features.csv`: the 82 Sorensen neurons and standardized 36-feature representation used for the RMS-shift and held-out-animal analyses.

`03_Morphology_mapping.py` requires five metadata fields in both morphometric input tables: `cell_id`, `dataset`, `broad_cell_type`, `animal_id`, and `official_M`. The published-neuron `official_M` values may be empty. The script reads the reference and published-neuron tables separately and combines them in memory. Fitted profiles, calibration parameters, and per-neuron mapping scores are generated as outputs.

`04_Summarize_RNA_and_Mtype_labels.py` can read the submitted mapping scores directly from Table 11. It can also accept fresh RNA and morphology mapping outputs when the full expression objects are available.

`05_Morphology_separation_analysis.R` uses `Published neuron morphology separation features.csv` together with `MetricsSet1.csv`.