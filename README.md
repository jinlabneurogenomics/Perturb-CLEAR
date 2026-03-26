# Perturb-CLEAR

We developed Perturb-CLEAR, which integrates pooled CRISPR screening and whole-mount imaging to quantify brain-wide cytoarchitecture, and paired it with Perturb-seq to link structural phenotypes to transcriptomic changes. Applying Perturb-CLEAR to the developing mouse cortex revealed morphogenesis trajectories accompanied by transcriptomic dynamics. Moreover, systematic perturbation of NDD risk genes uncovered gene-specific multimodal phenotypes. Combined morphology and transcriptome analyses link NDD risk genes to concordant multimodal cellular phenotypes in the developing brain, highlighting diverse paths of perturbation effect propagation across modalities.

## Image segmentation
**generate_output_windows.py**: Image segmentation to create crops of individual neurons using soma coordinates

## Morphology data QC and quantification
**swc_qc.py**: script-based quality control of .swc files to remove tracing artifacts  
**morphometric_feature_extraction.py**: Extract standard morphometrics from swc files

## Perturb-seq clustering and QC
See PerturbSeqAnalysis and the associated documentation for the code used in upstream analysis of the Perturb-seq data (starting with CellRanger output through clustering and cell type identification).

## Perturb-seq downstream analysis
**propeller.R**: Method to detect cell type proportion changes  
**run_DEG.R**: Identify DEGs using glmGamPoi

## Joint analysis
**TWI_calculation.R**: Method to calculate transcriptome wide impact
**eDist_calculation.R**: Method to calculate eDist

