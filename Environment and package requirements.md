# Environment and package requirements

## Environment requirements

| Environment | Runtime | Main use |
|---|---|---|
| Morphology and mapping | R 4.4.0 with Bioconductor 3.20 | Morphology statistics and ComBat normalization |
| RNA analysis | R 4.4.3 with Bioconductor 3.20 | Perturb-seq preprocessing, mapping, differential expression, and enrichment |
| Coexpression analysis | R 4.3.3 | hdWGCNA network construction and module summaries |
| Python analysis | Python 3.12 | Image processing, morphometric extraction, morphology mapping, and table summaries |

We reported the environments used for original analyses. Their differences reflect analysis history and package compatibility, not different analysis methods.

## R package requirements

| Package | Version | Released workflow |
|---|---:|---|
| Seurat | 5.4.0 | Perturb-seq object construction and RNA reference mapping |
| SeuratObject | 5.x | RNA assay, count, and metadata access |
| Azimuth | Record installed version | Perturb-seq cell annotation |
| harmony | 1.2.4 | Perturb-seq integration |
| glmGamPoi | 1.18.0 | Differential expression and subpopulation-resolved perturbation models |
| fgsea | 1.32.4 | Directional Gene Ontology enrichment |
| org.Mm.eg.db | 3.20.0 | Mouse gene annotation |
| GO.db | 3.20.0 | Gene Ontology annotation |
| sva | 3.54.0 | ComBat normalization with broad cell class protected |
| speckle | Record installed version | Perturbation-associated cell-composition testing |
| limma | 3.x | Cell-composition testing |
| edgeR | 3.x or later | Count processing in the cell-composition workflow |
| scater | Record installed version | Single-cell data handling for cell-composition analysis |
| SingleCellExperiment | 1.28.x | Cell-composition analysis and h5ad conversion |
| patchwork | 1.x | Figure assembly in the cell-composition workflow |
| e1071 | 1.7-4 | Morphology SVM classification |
| lmerTest | 3.x | Morphology mixed models |
| nlme | 3.1.x | Animal-aware morphology models |
| RcppML | Record installed version | Adnp morphology program analysis |
| umap | Record installed version | Morphology embedding |
| hdWGCNA | 0.4.11 | Coexpression network construction |
| WGCNA | 1.74 | Network and module calculations |
| zellkonverter | 1.16.x | h5ad import |
| data.table | 1.16.x | Tabular processing |
| Matrix | 1.7.x | Sparse matrix calculations |
| qs | 0.27.x | R object serialization |
| future | 1.34.x | Sequential execution control during Seurat mapping |
| DBI | 1.2.x | Annotation database access |
| RSQLite | 2.3.x | Annotation database queries |
| dplyr | 1.1.x | Tabular processing |
| tidyr | 1.3.x | Table reshaping |
| readr | 2.1.x | Delimited-file input |
| readxl | 1.4.x | Workbook input |
| openxlsx | 4.x | Workbook output |
| tibble | 3.2.x | Tabular output |
| purrr | 1.x | Iteration and list processing |
| ggplot2 | 3.x | Plotting |
| cowplot | 1.x | Figure assembly |

We recommend installing the Bioconductor release that matches the selected R version before installing Bioconductor packages.

## Python package requirements

| Package | Version | Released workflow |
|---|---:|---|
| numpy | 2.1.3 | Morphometric extraction and morphology mapping |
| pandas | 2.2.3 | Tabular processing and Table 11 input |
| scipy | 1.14.1 | Morphometric calculations and morphology mapping |
| openpyxl | NA | Reading Table 11 `.xlsx` worksheets through pandas |
| tifffile | NA | Image-window generation |
| tqdm | NA | Image-processing progress reporting |