rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(SeuratObject)
  library(dplyr)
  library(tidyr)
  library(tibble)
})
set.seed(12345)

# Command-line arguments, in order:
# 1. hdWGCNA object in RDS format.
# 2. T0/T1 subpopulation-assignment CSV.
# 3. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Expected 3 command-line arguments in the order documented above.")
}
input_rds <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
subpopulation_assignment_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_root <- args[[3]]
subpopulation_dir <- file.path(output_root, "subpopulation_effects")
dir.create(subpopulation_dir, recursive = TRUE, showWarnings = FALSE)
seurat_obj <- readRDS(input_rds)
meta <- seurat_obj@meta.data
required_metadata <- c("condition", "animal_id")
missing_metadata <- setdiff(required_metadata, colnames(meta))
if (length(missing_metadata) > 0) {
  stop("Missing required metadata: ", paste(missing_metadata, collapse = ", "))
}

subpopulations <- read.csv(subpopulation_assignment_file, stringsAsFactors = FALSE, check.names = FALSE)
if (!all(c("cell_id", "RNA_subpopulation") %in% colnames(subpopulations))) {
  stop("Subpopulation-assignment table must contain cell_id and RNA_subpopulation columns.")
}
if (anyDuplicated(subpopulations$cell_id)) stop("Subpopulation-assignment table has duplicated cell_id values.")
if (!all(subpopulations$RNA_subpopulation %in% c("T0", "T1"))) stop("RNA_subpopulation must contain only T0 or T1.")

module_cols <- grep("^Adnp-M[0-9]+$", colnames(meta), value = TRUE)
module_cols <- module_cols[order(as.integer(sub("Adnp-M", "", module_cols)))]
if (length(module_cols) != 6L) {
  stop("Expected six Adnp-M module eigengenes; found ", length(module_cols))
}

cell_me <- meta %>%
  rownames_to_column("cell_id") %>%
  select(cell_id, condition, animal_id, all_of(module_cols)) %>%
  left_join(subpopulations[, c("cell_id", "RNA_subpopulation")], by = "cell_id") %>%
  mutate(
    condition = factor(as.character(condition), levels = c("ST2", "Adnp")),
    animal_id = factor(animal_id),
    subpopulation = factor(as.character(RNA_subpopulation), levels = c("T0", "T1"))
  ) %>%
  filter(!is.na(condition), !is.na(subpopulation)) %>%
  pivot_longer(
    cols = all_of(module_cols), names_to = "module", values_to = "ME"
  ) %>%
  mutate(module = factor(module, levels = module_cols))

subpopulation_pseudobulk <- cell_me %>%
  group_by(animal_id, condition, subpopulation, module) %>%
  summarise(ME_mean = mean(ME, na.rm = TRUE), n_cells = n(), .groups = "drop")
write.csv(
  subpopulation_pseudobulk,
  file.path(subpopulation_dir, "subpopulation_pseudobulk_module_activity.csv"),
  row.names = FALSE
)

# This manuscript panel is descriptive. Calculate the Adnp-minus-ST2 change
# within each animal before averaging animals with equal weight. No P values or
# multiple-testing adjustments are calculated.
animal_deltas <- subpopulation_pseudobulk %>%
  select(animal_id, condition, subpopulation, module, ME_mean, n_cells) %>%
  pivot_wider(
    names_from = condition,
    values_from = c(ME_mean, n_cells),
    names_sep = "_"
  ) %>%
  filter(!is.na(ME_mean_ST2), !is.na(ME_mean_Adnp)) %>%
  mutate(delta_Adnp_minus_ST2 = ME_mean_Adnp - ME_mean_ST2)

write.csv(
  animal_deltas,
  file.path(subpopulation_dir, "subpopulation_animal_module_deltas_Adnp_vs_ST2.csv"),
  row.names = FALSE
)

effects <- animal_deltas %>%
  group_by(subpopulation, module) %>%
  summarise(
    estimate = mean(delta_Adnp_minus_ST2, na.rm = TRUE),
    n_paired_animals = n_distinct(animal_id),
    .groups = "drop"
  ) %>%
  mutate(comparison = paste0(subpopulation, ": Adnp - ST2"))

write.csv(
  effects,
  file.path(subpopulation_dir, "subpopulation_module_deltas_Adnp_vs_ST2.csv"),
  row.names = FALSE
)
writeLines(
  capture.output(sessionInfo()), file.path(subpopulation_dir, "sessionInfo.txt")
)
message("Descriptive subpopulation-dependent module summaries complete: ", subpopulation_dir)

