# Perturb-CLEAR analysis code and processed data

This repository contains the analysis scripts and compact data tables used for the Perturb-CLEAR morphology, Perturb-seq, and multimodal analyses.

## Contents

- `scripts/1_Morphology`: morphology extraction, clustering, screen analysis, classifier validation, and Adnp follow-up analyses.
- `scripts/2_Perturb_seq`: Perturb-seq object construction, perturbation relabeling, cell-composition analysis, and differential expression.
- `scripts/3_Multimodal`: fixed-reference mapping, RNA-subpopulation and morphology summaries, morphology separation, marker analysis, enrichment, and Adnp module analyses.
- `data/Morphology`: project-owned morphometric measurements, sample summaries, and the 36-feature panel.
- `data/Multimodal`: processed published-neuron morphometrics, the Perturb-CLEAR morphology reference, mapping scores, metadata, and the morphology-separation input.
- `Environment and package requirements.md`, `data/Published datasets.md`, and `data/Data package guide.md`: software requirements, source-data links, and the released data inventory.

## Getting started

Install the packages listed in `Environment and package requirements.md`. Run scripts for each relevant workflow. Each script lists its command-line arguments at the beginning.

## Data provenance

Sorensen et al. neurons are the primary paired RNA and morphology cohort. Scala et al. neurons provide a sensitivity cohort. Gouwens et al. morphology data contribute only to normalization. Source links are provided in `data/Published datasets.md`.