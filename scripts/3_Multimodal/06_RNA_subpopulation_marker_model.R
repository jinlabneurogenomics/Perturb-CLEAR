# Primary RNA-subpopulation marker model used for the manuscript.

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(qs)
  library(SeuratObject)
  library(glmGamPoi)
})

# Command-line arguments, in order:
# 1. Fixed four-subpopulation RNA reference in qs format.
# 2. Mapped query-object list in qs format.
# 3. Mapping-score CSV, optionally gzip-compressed.
# 4. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Expected 4 command-line arguments in the order documented above.")
}

reference_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
mapped_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
score_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[4]], winslash = "/", mustWork = FALSE)
table_dir <- file.path(output_dir, "tables")
checkpoint_dir <- file.path(output_dir, "checkpoints")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(checkpoint_dir, recursive = TRUE, showWarnings = FALSE)

fdr_cutoff <- 0.05
finite_effect_boundary <- 9.9

reference <- qread(reference_path)
mapped_objects <- qread(mapped_path)
query <- mapped_objects[["Sorensen2026"]]
if (is.null(query)) stop("Sorensen2026 is absent from the mapped-query object list.")

reference_counts <- LayerData(reference[["RNA"]], layer = "counts")
query_counts <- LayerData(query[["RNA"]], layer = "counts")
shared_genes <- intersect(rownames(reference_counts), rownames(query_counts))
if (length(shared_genes) <= 10000L) stop("Unexpectedly small shared-gene universe.")
reference_counts <- reference_counts[shared_genes, , drop = FALSE]
query_counts <- query_counts[shared_genes, , drop = FALSE]

reference_meta <- as.data.table(reference@meta.data, keep.rownames = "source_cell_id")
reference_meta[, `:=`(
  source_cell_id = as.character(source_cell_id),
  dataset = "Perturbseq_ST2",
  animal_id = as.character(orig.ident),
  broad_layer = fifelse(grepl("^L2/3", as.character(CT)), "L2/3 IT",
                        fifelse(grepl("^L4/5", as.character(CT)), "L4/5 IT", NA_character_)),
  subpopulation = as.character(RNA_subpopulation),
  subpopulation_probability = 1,
  label_source = "fixed ST2 training label"
)]
reference_meta[, valid_layer_subpopulation :=
  (broad_layer == "L4/5 IT" & subpopulation %chin% c("T0", "T1")) |
  (broad_layer == "L2/3 IT" & subpopulation %chin% c("T2", "T3"))]
reference_meta[, cell_key := paste0("PSEQ::", source_cell_id)]
colnames(reference_counts) <- reference_meta$cell_key[match(colnames(reference_counts), reference_meta$source_cell_id)]

query_meta <- as.data.table(query@meta.data, keep.rownames = "object_cell_id")
query_meta[, `:=`(
  source_cell_id = as.character(source_cell_id),
  dataset = "Sorensen2026",
  animal_id = as.character(animal_id),
  broad_layer = as.character(broad_cell_type)
)]

score_long <- fread(score_path)[dataset == "Sorensen2026"]
score_long[, source_cell_id := as.character(source_cell_id)]
score_wide <- dcast(score_long, source_cell_id ~ RNA_subpopulation, value.var = "probability")
internal_subpopulations <- paste0("T", 0:3)
if (!all(internal_subpopulations %in% names(score_wide))) stop("The complete four-subpopulation score vector is required.")
setnames(score_wide, internal_subpopulations, paste0("probability_T", 0:3))
query_meta <- merge(query_meta, score_wide, by = "source_cell_id", all.x = TRUE, sort = FALSE)
if (anyNA(query_meta$probability_T0)) stop("Mapping-score merge produced missing values.")

query_meta[, layer_pair_mass := fifelse(
  broad_layer == "L4/5 IT", probability_T0 + probability_T1,
  fifelse(broad_layer == "L2/3 IT", probability_T2 + probability_T3, NA_real_)
)]
query_meta[, `:=`(
  conditional_T0 = fifelse(broad_layer == "L4/5 IT", probability_T0 / pmax(layer_pair_mass, .Machine$double.eps), NA_real_),
  conditional_T1 = fifelse(broad_layer == "L4/5 IT", probability_T1 / pmax(layer_pair_mass, .Machine$double.eps), NA_real_),
  conditional_T2 = fifelse(broad_layer == "L2/3 IT", probability_T2 / pmax(layer_pair_mass, .Machine$double.eps), NA_real_),
  conditional_T3 = fifelse(broad_layer == "L2/3 IT", probability_T3 / pmax(layer_pair_mass, .Machine$double.eps), NA_real_)
)]
query_meta[, subpopulation := fifelse(
  broad_layer == "L4/5 IT", fifelse(conditional_T0 >= conditional_T1, "T0", "T1"),
  fifelse(broad_layer == "L2/3 IT", fifelse(conditional_T2 >= conditional_T3, "T2", "T3"), NA_character_)
)]
query_meta[, subpopulation_probability := fcase(
  subpopulation == "T0", conditional_T0,
  subpopulation == "T1", conditional_T1,
  subpopulation == "T2", conditional_T2,
  subpopulation == "T3", conditional_T3,
  default = NA_real_
)]
query_meta[, `:=`(
  label_source = "four-subpopulation scores renormalized within published broad layer",
  valid_layer_subpopulation = is.finite(layer_pair_mass) & layer_pair_mass > 0 & !is.na(subpopulation),
  cell_key = paste0("SORENSEN::", source_cell_id)
)]
colnames(query_counts) <- query_meta$cell_key[match(colnames(query_counts), query_meta$source_cell_id)]

manifest_columns <- c(
  "cell_key", "source_cell_id", "dataset", "animal_id", "broad_layer",
  "subpopulation", "label_source", "subpopulation_probability", "valid_layer_subpopulation"
)
cell_manifest <- rbindlist(list(
  reference_meta[, ..manifest_columns],
  query_meta[, ..manifest_columns]
), fill = TRUE)
all_counts <- cbind(reference_counts, query_counts)
cell_manifest <- cell_manifest[match(colnames(all_counts), cell_key)]
stopifnot(identical(cell_manifest$cell_key, colnames(all_counts)), !anyDuplicated(cell_manifest$cell_key))

fwrite(cell_manifest, file.path(table_dir, "cell_manifest.csv"))
fwrite(
  data.table(gene = shared_genes, gene_universe = "shared symbols in ST2 and Sorensen count assays"),
  file.path(table_dir, "shared_gene_universe.csv")
)

contrasts <- data.table(
  contrast_id = c("T0_vs_T1", "T2_vs_T3"),
  broad_layer = c("L4/5 IT", "L2/3 IT"),
  positive_subpopulation = c("T0", "T2"),
  negative_subpopulation = c("T1", "T3")
)

weighted_metrics <- function(counts, subpopulation_score) {
  library_size <- Matrix::colSums(counts)
  cpm <- counts %*% Diagonal(x = 1e6 / pmax(library_size, 1))
  negative_score <- 1 - subpopulation_score
  data.table(
    Gene = rownames(counts),
    weighted_detected_positive = as.numeric((counts > 0) %*% subpopulation_score) / sum(subpopulation_score),
    weighted_detected_negative = as.numeric((counts > 0) %*% negative_score) / sum(negative_score),
    weighted_mean_CPM_positive = as.numeric(cpm %*% subpopulation_score) / sum(subpopulation_score),
    weighted_mean_CPM_negative = as.numeric(cpm %*% negative_score) / sum(negative_score)
  )
}

profile_rows <- list()
result_rows <- list()
audit_rows <- list()

for (index in seq_len(nrow(contrasts))) {
  contrast <- contrasts[index]
  keep <- cell_manifest$valid_layer_subpopulation & cell_manifest$broad_layer == contrast$broad_layer
  meta <- copy(cell_manifest[keep])
  counts <- all_counts[, meta$cell_key, drop = FALSE]
  meta[, subpopulation_score := fifelse(subpopulation == contrast$positive_subpopulation, subpopulation_probability, 1 - subpopulation_probability)]
  stopifnot(all(is.finite(meta$subpopulation_score)), all(meta$subpopulation_score >= 0), all(meta$subpopulation_score <= 1))

  library_size <- Matrix::colSums(counts)
  cpm <- counts %*% Diagonal(x = 1e6 / pmax(library_size, 1))
  profile_rows[[index]] <- rbindlist(list(
    data.table(
      Gene = rownames(counts), contrast_id = contrast$contrast_id,
      subpopulation = contrast$positive_subpopulation,
      subpopulation_probability_mass = sum(meta$subpopulation_score),
      weighted_mean_CPM = as.numeric(cpm %*% meta$subpopulation_score) / sum(meta$subpopulation_score)
    ),
    data.table(
      Gene = rownames(counts), contrast_id = contrast$contrast_id,
      subpopulation = contrast$negative_subpopulation,
      subpopulation_probability_mass = sum(1 - meta$subpopulation_score),
      weighted_mean_CPM = as.numeric(cpm %*% (1 - meta$subpopulation_score)) / sum(1 - meta$subpopulation_score)
    )
  ))
  profile_rows[[index]][, `:=`(
    broad_layer = contrast$broad_layer,
    n_neurons = nrow(meta), n_animals = uniqueN(meta$animal_id),
    cohort = "ST2 plus Sorensen",
    normalization = "cell CPM weighted by within-layer continuous subpopulation score"
  )]

  keep_gene <- Matrix::rowMeans(counts > 0) > 0.01
  counts_fit <- counts[keep_gene, , drop = FALSE]
  meta[, technical_batch := factor(fifelse(dataset == "Sorensen2026", "Sorensen2026", paste0("Perturbseq_", animal_id)))]
  design <- model.matrix(~ technical_batch + subpopulation_score, data = meta)
  if (qr(design)$rank < ncol(design)) stop("Rank-deficient design for ", contrast$contrast_id)
  checkpoint <- file.path(checkpoint_dir, paste0(contrast$contrast_id, "_glmGamPoi.qs"))
  fit <- if (file.exists(checkpoint)) qread(checkpoint) else {
    fitted <- glm_gp(
      counts_fit, design = design, verbose = FALSE, on_disk = FALSE,
      subsample = TRUE, overdispersion_shrinkage = FALSE,
      do_cox_reid_adjustment = FALSE
    )
    qsave(fitted, checkpoint, preset = "balanced")
    fitted
  }
  result <- as.data.table(as.data.frame(test_de(
    fit, contrast = "subpopulation_score", pval_adjust_method = "BH", compute_lfc_se = FALSE
  )))
  result[, Gene := as.character(name)]
  result <- merge(result, weighted_metrics(counts_fit, meta$subpopulation_score), by = "Gene", all.x = TRUE, sort = FALSE)
  result[, `:=`(
    cohort = "ST2 plus Sorensen", contrast_id = contrast$contrast_id,
    broad_layer = contrast$broad_layer,
    positive_subpopulation = contrast$positive_subpopulation,
    negative_subpopulation = contrast$negative_subpopulation,
    p_value = as.numeric(pval), FDR = as.numeric(adj_pval), lfc = as.numeric(lfc),
    n_neurons = nrow(meta), n_animals = uniqueN(meta$animal_id),
    positive_probability_mass = sum(meta$subpopulation_score),
    negative_probability_mass = sum(1 - meta$subpopulation_score),
    model_formula = "~ technical_batch + subpopulation_score",
    gene_filter = "detected in more than 1 percent of fitted neurons",
    inference_role = "exploratory marker model"
  )]
  result_rows[[index]] <- result
  audit_rows[[index]] <- data.table(
    contrast_id = contrast$contrast_id, n_neurons = nrow(meta),
    n_ST2 = sum(meta$dataset == "Perturbseq_ST2"),
    n_Sorensen = sum(meta$dataset == "Sorensen2026"),
    n_genes_tested = nrow(result), design_columns = ncol(design),
    design_rank = qr(design)$rank,
    subpopulation_score_min = min(meta$subpopulation_score), subpopulation_score_max = max(meta$subpopulation_score)
  )
}

profiles <- rbindlist(profile_rows)
results <- rbindlist(result_rows, fill = TRUE)
results[, stable_marker := is.finite(FDR) & FDR < fdr_cutoff & is.finite(lfc) & abs(lfc) < finite_effect_boundary]
results[, marker_subpopulation := fifelse(
  stable_marker & lfc > 0, positive_subpopulation,
  fifelse(stable_marker & lfc < 0, negative_subpopulation, NA_character_)
)]

fwrite(profiles, file.path(table_dir, "subpopulation_expression_profiles.csv.gz"), compress = "gzip", na = "NA")
fwrite(results, file.path(table_dir, "marker_results.csv.gz"), compress = "gzip", na = "NA")
fwrite(results[stable_marker == TRUE, .(
  Gene, contrast_id, broad_layer, marker_subpopulation, log2FC = lfc, p_value, FDR,
  n_neurons, n_animals,
  selection_rule = "BH FDR < 0.05, finite estimate, absolute log2FC < 9.9"
)], file.path(table_dir, "stable_marker_sets.csv"), na = "NA")
fwrite(rbindlist(audit_rows), file.path(table_dir, "model_audit.csv"))
fwrite(data.table(
  parameter = c(
    "cohort", "reference_score_type", "query_score_type", "gene_detection_fraction",
    "model_formula", "FDR_cutoff", "absolute_log2FC_boundary"
  ),
  value = c(
    "ST2;Sorensen2026", "one_hot", "within_layer_conditional",
    0.01, "~ technical_batch + subpopulation_score", fdr_cutoff, finite_effect_boundary
  )
), file.path(table_dir, "run_parameters.csv"))
fwrite(data.table(
  software = c("R", "data.table", "Matrix", "qs", "SeuratObject", "glmGamPoi"),
  version = c(
    as.character(getRversion()), as.character(packageVersion("data.table")),
    as.character(packageVersion("Matrix")), as.character(packageVersion("qs")),
    as.character(packageVersion("SeuratObject")), as.character(packageVersion("glmGamPoi"))
  )
), file.path(table_dir, "software_versions.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "R_sessionInfo.txt"))
message("Marker analysis complete: ", table_dir)
