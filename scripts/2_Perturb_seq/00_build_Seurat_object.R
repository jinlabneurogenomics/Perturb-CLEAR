# Portable implementation of the Perturb-seq Seurat object construction.
# Laboratory-specific paths are replaced with command-line inputs.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Azimuth)
  library(harmony)
  library(qs)
  library(Seurat)
})

# Command-line arguments, in order:
# 1. Channel 1 Cell Ranger filtered-feature-matrix directory.
# 2. Channel 2 Cell Ranger filtered-feature-matrix directory.
# 3. Output directory.
# 4. Optional Azimuth reference name; defaults to mousecortexref.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Expected 3 required arguments and at most 1 optional argument in the order documented above.")
}

channel1_dir <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
channel2_dir <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
azimuth_reference <- if (length(args) == 4L) args[[4]] else "mousecortexref"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_gene_expression <- function(matrix_dir) {
  matrices <- Read10X(data.dir = matrix_dir)
  if (is.list(matrices)) {
    if (!"Gene Expression" %in% names(matrices)) {
      stop("Gene Expression matrix is absent from: ", matrix_dir)
    }
    return(matrices[["Gene Expression"]])
  }
  matrices
}

channel1 <- CreateSeuratObject(
  counts = read_gene_expression(channel1_dir),
  project = "Ch1"
)
channel2 <- CreateSeuratObject(
  counts = read_gene_expression(channel2_dir),
  project = "Ch2"
)

combined <- merge(channel1, y = channel2, project = "combined")
combined[["percent.mt"]] <- PercentageFeatureSet(combined, pattern = "^mt-")

qc_summary <- as.data.table(combined@meta.data, keep.rownames = "cell_id")[, .(
  neurons_before_filter = .N,
  median_detected_genes = median(nFeature_RNA),
  median_UMI_count = median(nCount_RNA),
  median_mitochondrial_percent = median(percent.mt)
), by = orig.ident]

combined <- subset(combined, subset = nFeature_RNA > 500 & percent.mt < 10)
combined <- NormalizeData(combined)
combined <- FindVariableFeatures(combined, assay = "RNA")
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 40, verbose = FALSE)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:30)
combined <- FindNeighbors(combined, reduction = "pca", dims = 1:30)
combined <- FindClusters(combined, resolution = 0.3)

combined <- IntegrateLayers(
  object = combined,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE
)
combined <- JoinLayers(combined)
combined <- RunAzimuth(combined, reference = azimuth_reference)

qc_after <- as.data.table(combined@meta.data, keep.rownames = "cell_id")[, .(
  neurons_after_filter = .N
), by = orig.ident]
qc_summary <- merge(qc_summary, qc_after, by = "orig.ident", all = TRUE)

qsave(combined, file.path(output_dir, "perturbseq_seurat_object.qs"))
fwrite(qc_summary, file.path(output_dir, "cell_filtering_summary.csv"))
fwrite(data.table(
  parameter = c(
    "minimum_detected_genes", "maximum_mitochondrial_percent",
    "PCA_components", "UMAP_and_neighbor_dimensions", "cluster_resolution",
    "Azimuth_reference"
  ),
  value = c("500", "10", "40", "1:30", "0.3", azimuth_reference)
), file.path(output_dir, "processing_parameters.csv"))
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))

message("Perturb-seq object construction complete: ", output_dir)
