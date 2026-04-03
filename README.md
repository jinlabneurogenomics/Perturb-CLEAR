# Perturb-CLEAR

We developed Perturb-CLEAR, which integrates pooled CRISPR screening and whole-mount imaging to quantify brain-wide cytoarchitecture, and paired it with Perturb-seq to link structural phenotypes to transcriptomic changes. Applying Perturb-CLEAR to the developing mouse cortex revealed morphogenesis trajectories accompanied by transcriptomic dynamics. Moreover, systematic perturbation of NDD risk genes uncovered gene-specific multimodal phenotypes. Combined morphology and transcriptome analyses link NDD risk genes to concordant multimodal cellular phenotypes in the developing brain, highlighting diverse paths of perturbation effect propagation across modalities.

## System requirements and installation guide
For softwares or packages, the system requiresments and installation guide can be found on their source websites: 
- **Imaris 10.2.0**: https://imaris.oxinst.com/support/installation-instructions
- **neuTube v1.0z**: https://neutracing.com/download/
- **Imagej v1.54r**: https://imagej.net/downloads
- **R v4.0.3**: https://cran.r-project.org/doc/manuals/r-release/R-admin.html
- **Python 3.12** https://www.python.org/downloads
The scripts will run on any system that with R and Python installed.

## Scripts
### Image segmentation
- **generate_output_windows.py**: Image processing to create crops of individual neurons using soma coordinates.

### Morphology data QC and quantification
- **swc_qc.py**: script-based quality control of .swc files to remove tracing artifacts  
- **morphometric_feature_extraction.py**: Extract standard morphometrics from swc files

### Perturb-CLEAR downstream analysis
- **Analysis.Data.R**: morphology-based clustering, cell type labeling, perturbation testing, and NMF analysis  
- **Cluster.py**: Contains helper functions used in the Analysis.Data.R

### Perturb-seq downstream analysis
- **DETest.R**: Identify DEGs using glmGamPoi

### Developmental Omics
- **ComplexHeatmap_L2-3_L4-5_IT.R**: Visualize gene expression patterns across neuronal subtypes
- **DESeq2_Testing_Mouse_DevVIS_Age.R**: Perform differential expression analysis across developmental ages
- **DEGs_Grouping_by_Expression_Trend.R**: Cluster DEGs based on shared expression dynamics
- **gProfileR_Enrich_by_Expression_Trend.R**: Run pathway and gene set enrichment analysis for genes grouped by expression trends
- **NDD_RiskGenes_FisherEnrichTest.R**: Test enrichment of neurodevelopmental disorder risk genes using Fisher’s exact test
