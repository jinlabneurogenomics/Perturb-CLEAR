# Subpopulation-resolved Adnp response model.
# The fixed four-subpopulation mapping is treated as an input. ST2 reference neurons
# receive one-hot scores; Adnp neurons retain their complete soft score vector.

options(stringsAsFactors = FALSE)
compiler::enableJIT(0)
suppressPackageStartupMessages({
  library(Matrix)
  library(SeuratObject)
  library(data.table)
  library(qs)
  library(glmGamPoi)
})

# Command-line arguments, in order:
# 1. Perturb-seq object in qs format.
# 2. All-cell reference-mapping CSV, optionally gzip-compressed.
# 3. Four-subpopulation mapping-score CSV, optionally gzip-compressed.
# 4. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Expected 4 command-line arguments in the order documented above.")
}

object_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
mapping_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
score_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[4]], winslash = "/", mustWork = FALSE)
table_dir <- file.path(output_dir, "tables")
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 260909L
control_assignment <- "ST_2"
adnp_assignment <- "Adnp_1"
fdr_cutoff <- 0.05
finite_effect_boundary <- 9.9
set.seed(seed)

object <- qread(object_path)
counts <- LayerData(object[["RNA"]], layer = "counts")
object_meta <- as.data.table(object@meta.data, keep.rownames = "source_cell_id")
object_meta[, `:=`(
  source_cell_id = as.character(source_cell_id),
  Assign = as.character(Assign),
  animal_id = as.character(orig.ident),
  CT = as.character(CT)
)]
rm(object)
invisible(gc())

mapping <- fread(mapping_path)
mapping <- mapping[dataset %chin% c("Perturbseq_ST2_reference", "Perturbseq_query")]
mapping[, source_cell_id := as.character(source_cell_id)]
mapping[, RNA_subpopulation := sub("STref C", "T", as.character(RNA_subpopulation))]

score_long <- fread(score_path)[dataset == "Perturbseq_query"]
score_long[, source_cell_id := as.character(source_cell_id)]
score_long[, RNA_subpopulation := sub("STref C", "T", as.character(RNA_subpopulation))]
score_wide <- dcast(score_long, source_cell_id ~ RNA_subpopulation, value.var = "probability")
for (subpopulation in paste0("T", 0:3)) if (!subpopulation %in% names(score_wide)) score_wide[, (subpopulation) := 0]
setnames(score_wide, paste0("T", 0:3), paste0("prob_T", 0:3))

cell_info <- merge(
  mapping[, .(
    source_cell_id, dataset, mapping_role, Assign, animal_id,
    source_broad_cell_type, RNA_subpopulation
  )],
  score_wide,
  by = "source_cell_id", all.x = TRUE, sort = FALSE
)
for (index in 0:3) {
  column <- paste0("prob_T", index)
  cell_info[Assign == control_assignment, (column) := as.numeric(RNA_subpopulation == paste0("T", index))]
}
probability_columns <- paste0("prob_T", 0:3)
if (anyNA(cell_info[, ..probability_columns])) stop("Incomplete four-subpopulation score vectors.")
if (max(abs(rowSums(cell_info[, ..probability_columns]) - 1)) > 1e-8) {
  stop("Four-subpopulation score vectors do not sum to one.")
}
cell_info[, probability_source := fifelse(
  Assign == control_assignment,
  "one-hot fixed-reference membership",
  "frozen RPCA transfer score"
)]
cell_info[, normalized_entropy :=
  -rowSums(as.matrix(.SD) * log(pmax(as.matrix(.SD), .Machine$double.eps))) / log(4),
  .SDcols = probability_columns]
fwrite(cell_info, file.path(table_dir, "cell_subpopulation_scores.csv.gz"), na = "NA")

analysis_meta <- cell_info[Assign %chin% c(control_assignment, adnp_assignment)]
analysis_meta <- merge(
  analysis_meta,
  object_meta[, .(source_cell_id, object_animal_id = animal_id, CT)],
  by = "source_cell_id", all.x = TRUE, sort = FALSE
)
analysis_meta[, animal_id := fifelse(is.na(animal_id) | !nzchar(animal_id), object_animal_id, animal_id)]
analysis_meta <- analysis_meta[match(colnames(counts), source_cell_id), nomatch = 0]
analysis_meta <- analysis_meta[Assign %chin% c(control_assignment, adnp_assignment)]
analysis_meta[, condition := relevel(factor(Assign), control_assignment)]
analysis_meta[, animal_id := factor(animal_id)]
counts_fit <- counts[, analysis_meta$source_cell_id, drop = FALSE]
keep_gene <- Matrix::rowMeans(counts_fit > 0) > 0.01
counts_fit <- counts_fit[keep_gene, , drop = FALSE]

base_design <- model.matrix(
  ~ animal_id + prob_T0 + prob_T1 + prob_T2 + condition,
  data = analysis_meta
)
condition_indicator <- as.numeric(analysis_meta$Assign == adnp_assignment)
interaction_design <- cbind(
  INT_Adnp_T0 = condition_indicator * analysis_meta$prob_T0,
  INT_Adnp_T1 = condition_indicator * analysis_meta$prob_T1,
  INT_Adnp_T2 = condition_indicator * analysis_meta$prob_T2
)
design <- cbind(base_design, interaction_design)
storage.mode(design) <- "double"
if (qr(design)$rank != ncol(design)) stop("The subpopulation-resolved Adnp design is rank deficient.")

condition_coefficient <- paste0("condition", adnp_assignment)
if (!condition_coefficient %in% colnames(design)) stop("Adnp condition coefficient is absent.")
fit_path <- file.path(checkpoint_dir, "adnp_four_subpopulation_glmGamPoi.qs")
fit <- if (file.exists(fit_path)) {
  qread(fit_path)
} else {
  fitted <- glm_gp(
    counts_fit, design = design, verbose = TRUE, on_disk = FALSE,
    subsample = TRUE, overdispersion_shrinkage = FALSE,
    do_cox_reid_adjustment = FALSE
  )
  qsave(fitted, fit_path, preset = "balanced")
  fitted
}

subpopulation_contrasts <- c(
  T0 = paste(condition_coefficient, "+ INT_Adnp_T0"),
  T1 = paste(condition_coefficient, "+ INT_Adnp_T1"),
  T2 = paste(condition_coefficient, "+ INT_Adnp_T2"),
  T3 = condition_coefficient
)
interaction_contrasts <- c(
  T0_vs_T1 = "INT_Adnp_T0 - INT_Adnp_T1",
  T2_vs_T3 = "INT_Adnp_T2"
)

run_test <- function(contrast, test_type, label) {
  result <- as.data.table(test_de(
    fit, contrast = contrast, pval_adjust_method = "BH",
    compute_lfc_se = TRUE, verbose = FALSE
  ))
  setnames(result, "name", "Gene")
  result[, `:=`(
    perturbation = "Adnp", test_type = test_type,
    subpopulation_or_pair = label, contrast = contrast,
    p_value = as.numeric(pval), FDR = as.numeric(adj_pval),
    log2FC = as.numeric(lfc),
    stable_call = is.finite(adj_pval) & adj_pval < fdr_cutoff &
      is.finite(lfc) & abs(lfc) < finite_effect_boundary
  )]
  result
}

result_rows <- list()
for (subpopulation in names(subpopulation_contrasts)) {
  result_rows[[paste0("subpopulation_", subpopulation)]] <- run_test(
    subpopulation_contrasts[[subpopulation]], "subpopulation_specific_Adnp_effect", subpopulation
  )
}
for (pair in names(interaction_contrasts)) {
  result_rows[[paste0("interaction_", pair)]] <- run_test(
    interaction_contrasts[[pair]], "Adnp_by_subpopulation_interaction", pair
  )
}
results <- rbindlist(result_rows, fill = TRUE)
fwrite(results, file.path(table_dir, "adnp_subpopulation_resolved_gene_results.csv.gz"), na = "NA")
fwrite(results[test_type == "subpopulation_specific_Adnp_effect", .(
  genes_tested = .N,
  stable_DEGs_FDR_0_05 = sum(stable_call),
  median_absolute_log2FC = median(abs(log2FC), na.rm = TRUE)
), by = subpopulation_or_pair], file.path(table_dir, "adnp_subpopulation_response_summary.csv"), na = "NA")
fwrite(results[stable_call == TRUE, .(
  Gene, test_type, subpopulation_or_pair, log2FC, p_value, FDR, stable_call
)], file.path(table_dir, "adnp_stable_DEGs.csv"), na = "NA")

fwrite(data.table(
  parameter = c(
    "seed", "control_assignment", "perturbation_assignment",
    "gene_detection_fraction", "FDR_cutoff", "absolute_log2FC_boundary",
    "reference_subpopulation", "multiple_testing"
  ),
  value = c(
    seed, control_assignment, adnp_assignment, 0.01, fdr_cutoff,
    finite_effect_boundary, "T3", "BH_within_contrast"
  )
), file.path(table_dir, "run_parameters.csv"))
fwrite(data.table(
  n_cells = nrow(analysis_meta), n_genes = nrow(counts_fit),
  n_control_cells = sum(analysis_meta$Assign == control_assignment),
  n_Adnp_cells = sum(analysis_meta$Assign == adnp_assignment),
  n_animals = uniqueN(analysis_meta$animal_id),
  n_design_columns = ncol(design), design_rank = qr(design)$rank,
  seed = seed
), file.path(table_dir, "model_audit.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "R_sessionInfo.txt"))
message("Subpopulation-resolved Adnp analysis complete: ", table_dir)
