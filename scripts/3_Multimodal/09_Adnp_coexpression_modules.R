rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(zellkonverter)
  library(SingleCellExperiment)
  library(Seurat)
  library(dplyr)
  library(WGCNA)
  library(hdWGCNA)
})
set.seed(12345)

# Command-line arguments, in order:
# 1. L4/5 IT Adnp and ST2 expression object in h5ad format.
# 2. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Expected 2 command-line arguments in the order documented above.")
}
input_h5ad <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_root <- args[[2]]
module_dir <- file.path(output_root, "modules")
dir.create(module_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(module_dir, "TOM"), recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 16)
required_conditions <- c("ST2", "Adnp")

message("Reading ", input_h5ad)
sce <- readH5AD(input_h5ad)
stopifnot(all(c("X", "counts") %in% assayNames(sce)))
seurat_obj <- as.Seurat(sce, counts = "counts", data = "X")
rm(sce)
gc()

required_metadata <- c("condition", "animal_id", "broad_cell_type")
missing_metadata <- setdiff(required_metadata, colnames(seurat_obj@meta.data))
if (length(missing_metadata) > 0) {
  stop("Missing required metadata: ", paste(missing_metadata, collapse = ", "))
}
if (!all(as.character(seurat_obj$broad_cell_type) == "L4/5 IT")) {
  stop("Input contains cells outside L4/5 IT.")
}

seurat_obj$condition <- factor(seurat_obj$condition, levels = required_conditions)
seurat_obj$animal_id <- factor(seurat_obj$animal_id)
condition_by_animal <- table(seurat_obj$animal_id, seurat_obj$condition)
paired_animals <- rownames(condition_by_animal)[apply(condition_by_animal > 0, 1, all)]
if (length(paired_animals) < 3L) {
  stop("Fewer than three animals contain both ST2 and Adnp neurons.")
}
seurat_obj <- subset(seurat_obj, subset = animal_id %in% paired_animals)
seurat_obj$animal_id <- droplevels(seurat_obj$animal_id)

cell_counts <- seurat_obj@meta.data %>%
  count(animal_id, condition, name = "n_cells") %>%
  arrange(animal_id, condition)
write.csv(
  cell_counts, file.path(module_dir, "cell_counts_by_animal_condition.csv"),
  row.names = FALSE
)

DefaultAssay(seurat_obj) <- "originalexp"
available_reductions <- Reductions(seurat_obj)
reference_pca <- grep(
  "st2.*reference.*pca|reference.*pca",
  available_reductions, value = TRUE, ignore.case = TRUE
)
if (length(reference_pca) != 1L) {
  stop(
    "Expected one transferred ST2-reference PCA reduction. Found: ",
    paste(available_reductions, collapse = ", ")
  )
}
metacell_reduction <- reference_pca[[1]]
message("Using transferred ST2-reference reduction: ", metacell_reduction)

seurat_obj <- SetupForWGCNA(
  seurat_obj,
  gene_select = "fraction",
  fraction = 0.05,
  wgcna_name = "Adnp_vs_ST2_paired"
)
seurat_obj <- MetacellsByGroups(
  seurat_obj = seurat_obj,
  group.by = c("condition", "animal_id"),
  reduction = metacell_reduction,
  min_cells = 20,
  k = 5,
  ident.group = "condition"
)
seurat_obj <- NormalizeMetacells(seurat_obj)
seurat_obj <- SetDatExpr(
  seurat_obj,
  group_name = required_conditions,
  group.by = "condition",
  assay = "originalexp",
  layer = "data"
)

seurat_obj <- TestSoftPowers(seurat_obj, networkType = "signed")
power_table <- GetPowerTable(seurat_obj)
write.csv(
  power_table, file.path(module_dir, "soft_power_table.csv"), row.names = FALSE
)
selected_power <- power_table$Power[which(power_table$SFT.R.sq >= 0.8)[1]]
if (!identical(as.integer(selected_power), 7L)) {
  stop("The smallest soft power with scale-free fit >= 0.8 is no longer 7.")
}

old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(module_dir)
seurat_obj <- ConstructNetwork(
  seurat_obj,
  soft_power = 7,
  corType = "pearson",
  networkType = "signed",
  TOMType = "signed",
  minModuleSize = 50,
  mergeCutHeight = 0.2,
  randomSeed = 12345,
  tom_name = "Adnp_ST2_paired_TOM",
  tom_outdir = file.path(module_dir, "TOM"),
  overwrite_tom = TRUE
)

seurat_obj <- ModuleEigengenes(seurat_obj, assay = "originalexp")
seurat_obj <- ModuleConnectivity(
  seurat_obj,
  group.by = "condition",
  group_name = required_conditions,
  harmonized = FALSE,
  assay = "originalexp"
)
seurat_obj <- ResetModuleNames(seurat_obj, new_name = "Adnp-M")
seurat_obj <- ModuleEigengenes(seurat_obj, assay = "originalexp")

module_eigengenes <- GetMEs(seurat_obj, harmonized = FALSE)
old_me_columns <- intersect(colnames(module_eigengenes), colnames(seurat_obj@meta.data))
if (length(old_me_columns) > 0) {
  seurat_obj@meta.data <- seurat_obj@meta.data[
    , !colnames(seurat_obj@meta.data) %in% old_me_columns, drop = FALSE
  ]
}
seurat_obj <- AddMetaData(seurat_obj, metadata = module_eigengenes)

modules <- GetModules(seurat_obj) %>% filter(module != "grey")
hubs <- GetHubGenes(seurat_obj, n_hubs = 100) %>% filter(module != "grey")
expected_modules <- paste0("Adnp-M", 1:6)
observed_modules <- sort(unique(as.character(modules$module)))
if (!identical(observed_modules, expected_modules)) {
  stop(
    "Expected exactly six renamed modules (Adnp-M1 to Adnp-M6); found: ",
    paste(observed_modules, collapse = ", ")
  )
}
if (!all(expected_modules %in% unique(as.character(hubs$module)))) {
  stop("At least one expected module is missing from the hub-gene table.")
}
write.csv(
  modules, file.path(module_dir, "hdWGCNA_module_assignments.csv"),
  row.names = FALSE
)
write.csv(
  hubs, file.path(module_dir, "hdWGCNA_hub_genes_top100.csv"),
  row.names = FALSE
)

saveRDS(
  seurat_obj,
  file.path(module_dir, "hdWGCNA_final_Adnp_vs_ST2_paired.rds")
)
writeLines(
  capture.output(sessionInfo()), file.path(module_dir, "sessionInfo.txt")
)
message("Module construction complete: ", module_dir)
