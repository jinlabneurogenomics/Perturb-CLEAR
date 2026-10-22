options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(qs)
  library(Matrix)
  library(glmGamPoi)
  library(dplyr)
  library(readr)
})

# Command-line arguments, in order:
# 1. Perturb-seq Seurat object in qs format.
# 2. Output directory.
# 3. Optional blocking-covariate name: Sample or orig.ident; defaults to orig.ident.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L || length(args) > 3L) {
  stop("Expected 2 required arguments and at most 1 optional argument in the order documented above.")
}
input_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
blocking_covariate <- if (length(args) == 3L) args[[3]] else "orig.ident"
if (!blocking_covariate %in% c("Sample", "orig.ident")) {
  stop("Blocking covariate must be 'Sample' or 'orig.ident'.")
}
blocking_tag <- if (blocking_covariate == "Sample") "Sample" else "origIdent"
design_label <- paste0("~Assign+", blocking_covariate)
checkpoint_dir <- file.path(output_dir, "checkpoints")
nominal_root <- file.path(output_dir, "nominal_csv")
for (d in c(output_dir, checkpoint_dir, nominal_root)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}
ct_include <- c("L2/3 IT CTX-1", "L4/5 IT CTX")
references <- c("ST2" = "ST_2", "NT2" = "NT_2")
min_num_assign <- 10L
filter_detection_fraction <- 0.01
log_file <- file.path(output_dir, "shared_fit_progress.log")

log_msg <- function(...) {
  msg <- paste0(
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ",
    paste(..., collapse = "")
  )
  cat(msg, "\n")
  cat(msg, "\n", file = log_file, append = TRUE)
}

sanitize_file <- function(x) gsub("[^A-Za-z0-9]+", "_", x)
sheet_name_for <- function(ct, assign) {
  x <- paste0(gsub("/", "_", ct, fixed = TRUE), "_Assign", assign)
  substr(x, 1, 31)
}

normalize_result_columns <- function(res, gene_names) {
  out <- as.data.frame(res, stringsAsFactors = FALSE)
  required_numeric <- c("pval", "adj_pval", "f_statistic", "df1", "df2", "lfc")
  for (nm in required_numeric) {
    if (!nm %in% names(out)) out[[nm]] <- NA_real_
    out[[nm]] <- as.numeric(out[[nm]])
  }
  if (!"name" %in% names(out)) out$name <- gene_names
  out$Gene <- gene_names
  out <- out[, c("name", "pval", "adj_pval", "f_statistic", "df1", "df2", "lfc", "Gene")]
  out[order(out$pval, na.last = TRUE), , drop = FALSE]
}

contrast_for <- function(assign, reference, coefficient_names) {
  coefficient <- paste0("Assign", assign)
  nt_coefficient <- "AssignNT_2"
  if (reference == "ST_2") {
    if (!coefficient %in% coefficient_names) {
      stop("Missing coefficient for ST2 contrast: ", coefficient)
    }
    return(coefficient)
  }
  if (reference != "NT_2") stop("Unsupported reference: ", reference)
  if (!nt_coefficient %in% coefficient_names) stop("Missing AssignNT_2 coefficient.")
  if (assign == "ST_2") return("-1 * AssignNT_2")
  if (!coefficient %in% coefficient_names) {
    stop("Missing coefficient for NT2 contrast: ", coefficient)
  }
  paste(coefficient, "-", nt_coefficient)
}

write_nominal_csv <- function(tab, reference_tag, ct, assign) {
  outdir <- file.path(nominal_root, reference_tag)
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  path <- file.path(outdir, paste0(sheet_name_for(ct, assign), ".csv"))
  write_csv(tab %>% filter(is.finite(pval), pval < 0.05), path, na = "NA")
  path
}

log_msg("Loading current subset object: ", input_file)
obj <- qread(input_file)
if (!blocking_covariate %in% colnames(obj@meta.data)) {
  stop("Current object lacks blocking covariate: ", blocking_covariate)
}
if ("BioSamp" %in% colnames(obj@meta.data)) {
  stop("Unexpected BioSamp found; this script is intentionally blocked by ", blocking_covariate, ".")
}
counts <- obj@assays[["RNA"]]@counts
meta <- obj@meta.data
meta$cell_id <- rownames(meta)
meta$CT <- as.character(meta$CT)
meta$Assign <- as.character(meta$Assign)
meta[[blocking_covariate]] <- as.character(meta[[blocking_covariate]])
meta <- meta %>% filter(
  CT %in% ct_include,
  !is.na(Assign),
  !is.na(.data[[blocking_covariate]])
)
counts <- counts[, meta$cell_id, drop = FALSE]
rm(obj)
invisible(gc())

cell_counts <- meta %>% count(CT, Assign, name = "n_cells") %>% arrange(CT, Assign)
sample_balance <- meta %>%
  count(CT, .data[[blocking_covariate]], Assign, name = "n_cells") %>%
  arrange(CT, .data[[blocking_covariate]], Assign)
write_csv(cell_counts, file.path(output_dir, "cell_counts_by_CT_Assign.csv"))
write_csv(
  sample_balance,
  file.path(output_dir, paste0("sample_balance_by_CT_", blocking_tag, "_Assign.csv"))
)

all_results <- list(ST2 = list(), NT2 = list())
all_summaries <- list()

for (ct_index in seq_along(ct_include)) {
  ct <- ct_include[[ct_index]]
  ct_tag <- sanitize_file(ct)
  checkpoint_file <- file.path(checkpoint_dir, paste0(ct_tag, "__shared_results.qs"))
  if (file.exists(checkpoint_file)) {
    saved <- qread(checkpoint_file)
    all_results$ST2 <- c(all_results$ST2, saved$results$ST2)
    all_results$NT2 <- c(all_results$NT2, saved$results$NT2)
    all_summaries[[ct]] <- saved$summary
    log_msg("Reused completed CT checkpoint: ", ct)
    next
  }

  meta_ct <- meta %>% filter(CT == ct)
  assignment_counts <- table(meta_ct$Assign)
  included_assignments <- names(assignment_counts)[assignment_counts > min_num_assign]
  if (!all(unname(references) %in% included_assignments)) {
    stop("Both ST_2 and NT_2 must have >", min_num_assign, " cells in ", ct)
  }
  meta_ct <- meta_ct %>% filter(Assign %in% included_assignments) %>% droplevels()
  counts_ct <- counts[, meta_ct$cell_id, drop = FALSE]
  keep_genes <- Matrix::rowMeans(counts_ct > 0) > filter_detection_fraction
  counts_ct <- counts_ct[keep_genes, , drop = FALSE]

  meta_ct$Assign <- relevel(factor(meta_ct$Assign), "ST_2")
  meta_ct[[blocking_covariate]] <- factor(meta_ct[[blocking_covariate]])
  design_formula <- reformulate(c("Assign", blocking_covariate))
  design <- model.matrix(design_formula, meta_ct)
  if (qr(design)$rank < ncol(design)) {
    stop("Shared ", design_label, " design is rank-deficient for ", ct,
         "; no automatic fallback is allowed for this controlled run.")
  }
  coefficient_names <- grep("^Assign", colnames(design), value = TRUE)
  fit_checkpoint <- file.path(checkpoint_dir, paste0(ct_tag, "__shared_fit.qs"))
  fit_metadata_file <- file.path(checkpoint_dir, paste0(ct_tag, "__shared_fit_metadata.rds"))
  if (!file.exists(fit_checkpoint)) {
    set.seed(260716L + ct_index)
    log_msg(
      "Shared fit start: CT=", ct,
      " cells=", ncol(counts_ct),
      " genes=", nrow(counts_ct),
      " assignments=", length(levels(meta_ct$Assign)),
      " blocking_covariate=", blocking_covariate,
      " block_levels=", length(levels(meta_ct[[blocking_covariate]])),
      " design_columns=", ncol(design)
    )
    fit_time <- system.time({
      fit <- glm_gp(
        counts_ct,
        design = design,
        verbose = FALSE,
        on_disk = FALSE,
        subsample = TRUE,
        overdispersion_shrinkage = FALSE,
        do_cox_reid_adjustment = FALSE
      )
    })
    fit_elapsed_sec <- unname(fit_time[["elapsed"]])
    log_msg("Shared fit complete: CT=", ct, " elapsed_sec=", round(fit_elapsed_sec, 1))
    log_msg("Saving shared fit checkpoint: ", fit_checkpoint)
    qsave(fit, fit_checkpoint, preset = "fast")
    saveRDS(list(fit_elapsed_sec = fit_elapsed_sec), fit_metadata_file)
    rm(fit)
    invisible(gc())
  } else {
    fit_elapsed_sec <- if (file.exists(fit_metadata_file)) {
      readRDS(fit_metadata_file)$fit_elapsed_sec
    } else {
      NA_real_
    }
    log_msg("Reusing shared fit checkpoint: ", fit_checkpoint)
  }

  contrast_checkpoint_dir <- file.path(checkpoint_dir, ct_tag, "contrast_results")
  dir.create(contrast_checkpoint_dir, recursive = TRUE, showWarnings = FALSE)
  tasks <- list()
  task_index <- 0L
  for (reference_tag in names(references)) {
    reference <- references[[reference_tag]]
    tested_assignments <- sort(setdiff(levels(meta_ct$Assign), reference))
    for (assign in tested_assignments) {
      task_index <- task_index + 1L
      sheet_name <- sheet_name_for(ct, assign)
      nominal_dir <- file.path(nominal_root, reference_tag)
      dir.create(nominal_dir, recursive = TRUE, showWarnings = FALSE)
      tasks[[task_index]] <- list(
        reference_tag = reference_tag,
        reference = reference,
        CT = ct,
        Assign = assign,
        sheet_name = sheet_name,
        contrast = contrast_for(assign, reference, coefficient_names),
        n_cells_CT_fit = ncol(counts_ct),
        n_cells_assign = sum(meta_ct$Assign == assign),
        n_cells_reference = sum(meta_ct$Assign == reference),
        n_assignments_in_fit = length(levels(meta_ct$Assign)),
        n_samples_in_fit = length(levels(meta_ct[[blocking_covariate]])),
        n_block_levels_in_fit = length(levels(meta_ct[[blocking_covariate]])),
        blocking_covariate = blocking_covariate,
        fit_elapsed_sec = fit_elapsed_sec,
        nominal_file = file.path(nominal_dir, paste0(sheet_name, ".csv")),
        result_checkpoint = file.path(
          contrast_checkpoint_dir,
          paste0(reference_tag, "__", sanitize_file(assign), ".qs")
        )
      )
    }
  }

  run_contrast_task <- function(task) {
    if (file.exists(task$result_checkpoint)) return(task$result_checkpoint)
    test_time <- system.time({
      res <- glmGamPoi::test_de(
        fit_worker,
        contrast = task$contrast,
        pval_adjust_method = "BH",
        compute_lfc_se = FALSE
      )
    })
    full <- normalize_result_columns(res, gene_names_worker)
    nominal <- full[is.finite(full$pval) & full$pval < 0.05, , drop = FALSE]
    utils::write.csv(nominal, task$nominal_file, row.names = FALSE, na = "NA")
    summary <- data.frame(
      reference = task$reference,
      CT = task$CT,
      Assign = task$Assign,
      sheet_name = task$sheet_name,
      contrast = task$contrast,
      n_cells_CT_fit = task$n_cells_CT_fit,
      n_cells_assign = task$n_cells_assign,
      n_cells_reference = task$n_cells_reference,
      n_assignments_in_fit = task$n_assignments_in_fit,
      n_samples_in_fit = task$n_samples_in_fit,
      n_block_levels_in_fit = task$n_block_levels_in_fit,
      blocking_covariate = task$blocking_covariate,
      n_genes_tested = nrow(full),
      n_nominal_p_lt_0.05 = sum(full$pval < 0.05, na.rm = TRUE),
      n_FDR_lt_0.05 = sum(full$adj_pval < 0.05, na.rm = TRUE),
      n_FDR_lt_0.10 = sum(full$adj_pval < 0.10, na.rm = TRUE),
      n_FDR_lt_0.20 = sum(full$adj_pval < 0.20, na.rm = TRUE),
      design = design_label,
      fit_elapsed_sec = task$fit_elapsed_sec,
      test_elapsed_sec = unname(test_time[["elapsed"]]),
      nominal_csv = task$nominal_file,
      stringsAsFactors = FALSE
    )
    qs::qsave(
      list(
        reference_tag = task$reference_tag,
        sheet_name = task$sheet_name,
        result = full,
        summary = summary
      ),
      task$result_checkpoint,
      preset = "fast"
    )
    task$result_checkpoint
  }

  pending_tasks <- tasks[!vapply(tasks, function(x) file.exists(x$result_checkpoint), logical(1))]
  if (length(pending_tasks)) {
    worker_count <- min(8L, length(pending_tasks))
    log_msg(
      "Parallel contrast extraction start: CT=", ct,
      " pending=", length(pending_tasks),
      " workers=", worker_count
    )
    worker_log <- file.path(output_dir, paste0("parallel_workers_", ct_tag, ".log"))
    cl <- parallel::makeCluster(worker_count, outfile = worker_log)
    on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
    parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())
    gene_names_worker <- rownames(counts_ct)
    parallel::clusterExport(
      cl,
      c(
        "fit_checkpoint", "gene_names_worker", "normalize_result_columns",
        "run_contrast_task", "design_label"
      ),
      envir = environment()
    )
    parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(qs))
      suppressPackageStartupMessages(library(glmGamPoi))
      fit_worker <- qs::qread(fit_checkpoint)
      NULL
    })
    completed_files <- parallel::parLapplyLB(cl, pending_tasks, run_contrast_task)
    parallel::stopCluster(cl)
    on.exit(NULL, add = FALSE)
    log_msg("Parallel contrast extraction complete: CT=", ct)
  } else {
    log_msg("All contrast checkpoints already exist: CT=", ct)
  }

  expected_files <- vapply(tasks, `[[`, character(1), "result_checkpoint")
  if (!all(file.exists(expected_files))) {
    stop("Missing contrast checkpoints after extraction for ", ct)
  }
  completed <- lapply(expected_files, qread)
  ct_results <- list(ST2 = list(), NT2 = list())
  ct_summary <- list()
  for (i in seq_along(completed)) {
    item <- completed[[i]]
    ct_results[[item$reference_tag]][[item$sheet_name]] <- item$result
    ct_summary[[i]] <- item$summary
  }
  ct_summary <- bind_rows(ct_summary)
  qsave(
    list(results = ct_results, summary = ct_summary),
    checkpoint_file,
    preset = "fast"
  )
  all_results$ST2 <- c(all_results$ST2, ct_results$ST2)
  all_results$NT2 <- c(all_results$NT2, ct_results$NT2)
  all_summaries[[ct]] <- ct_summary
  rm(counts_ct, meta_ct, design, ct_results, completed)
  invisible(gc())
}

summary_table <- bind_rows(all_summaries) %>% arrange(reference, CT, Assign)
write_csv(summary_table, file.path(output_dir, paste0("DEG_statistics_sharedFit_", blocking_tag, "_NT2_ST2_L23_L45.csv")))
write_csv(
  summary_table %>%
    group_by(reference, CT) %>%
    summarise(
      n_comparisons = n(),
      total_nominal = sum(n_nominal_p_lt_0.05),
      total_FDR_0.05 = sum(n_FDR_lt_0.05),
      total_FDR_0.10 = sum(n_FDR_lt_0.10),
      total_FDR_0.20 = sum(n_FDR_lt_0.20),
      median_nominal = median(n_nominal_p_lt_0.05),
      median_FDR_0.05 = median(n_FDR_lt_0.05),
      median_FDR_0.20 = median(n_FDR_lt_0.20),
      .groups = "drop"
    ),
  file.path(output_dir, paste0("DEG_statistics_aggregate_sharedFit_", blocking_tag, "_NT2_ST2_L23_L45.csv"))
)

qsave(all_results$ST2, file.path(output_dir, paste0("DE.sharedFit.", blocking_tag, ".ST2_L23_L45.qs")), preset = "high")
qsave(all_results$NT2, file.path(output_dir, paste0("DE.sharedFit.", blocking_tag, ".NT2_L23_L45.qs")), preset = "high")
qsave(
  list(
    run_parameters = list(
      design = design_label,
      blocking_covariate = blocking_covariate,
      reference_assignments = unname(references),
      gene_detection_fraction = filter_detection_fraction,
      minimum_assignment_cell_count_exclusive = min_num_assign,
      subsample = TRUE,
      overdispersion_shrinkage = FALSE,
      do_cox_reid_adjustment = FALSE,
      compute_lfc_se = FALSE
    ),
    results = all_results,
    summary = summary_table,
    cell_counts = cell_counts,
    blocking_balance = sample_balance
  ),
  file.path(output_dir, paste0("DE.sharedFit.", blocking_tag, ".NT2_ST2.combined_L23_L45.qs")),
  preset = "high"
)

run_parameters <- tibble::tibble(
  parameter = c(
    "input_file", "cell_types", "blocking_covariate", "design",
    "reference_assignments", "gene_detection_fraction",
    "minimum_assignment_cell_count_exclusive", "subsample",
    "overdispersion_shrinkage", "do_cox_reid_adjustment", "compute_lfc_se"
  ),
  value = c(
    normalizePath(input_file, winslash = "/", mustWork = TRUE),
    paste(ct_include, collapse = ";"), blocking_covariate, design_label,
    paste(unname(references), collapse = ";"), filter_detection_fraction,
    min_num_assign, TRUE, FALSE, FALSE, FALSE
  )
)
write_csv(
  run_parameters,
  file.path(output_dir, paste0("run_parameters_sharedFit_", blocking_tag, "_L23_L45.csv"))
)

log_msg("All shared-fit outputs completed: ", output_dir)
