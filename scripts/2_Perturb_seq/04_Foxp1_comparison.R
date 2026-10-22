suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
})

# Command-line arguments, in order:
# 1. Perturb-seq differential-expression CSV.
# 2. Ortiz et al. Foxp1 differential-expression workbook.
# 3. Output directory.
# 4. Optional published cell-type label; defaults to the L4/5 IT label used in the analysis.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3L || length(args) > 4L) {
  stop("Expected 3 required arguments and at most 1 optional argument in the order documented above.")
}
pseq_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
ortiz_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
ortiz_celltype <- if (length(args) == 4L) args[[4]] else "006 L4_5 IT CTX Glut.xlsx"
fdr_cutoff <- 0.05
example_labels <- c("Pdzrn3", "Ppfibp1", "Grid2", "Ncam2")

pseq <- read.csv(pseq_file, check.names = FALSE, stringsAsFactors = FALSE) %>%
  transmute(
    Gene = as.character(Gene),
    pseq_lfc = as.numeric(lfc),
    pseq_pval = as.numeric(pval),
    pseq_fdr = as.numeric(adj_pval)
  )
ortiz_all <- read_excel(ortiz_file, sheet = "combined", skip = 1)
ortiz <- ortiz_all %>%
  filter(celltype == ortiz_celltype) %>%
  transmute(
    Gene = as.character(Gene),
    ortiz_lfc = as.numeric(logFC),
    ortiz_pval = as.numeric(PValue),
    ortiz_fdr = as.numeric(padj)
  )

stopifnot(!anyDuplicated(pseq$Gene), !anyDuplicated(ortiz$Gene))
stopifnot(all(pseq$pseq_pval < 0.05, na.rm = TRUE))
stopifnot(all(ortiz$ortiz_pval < 0.05, na.rm = TRUE))

joined_all <- inner_join(pseq, ortiz, by = "Gene") %>%
  mutate(
    capped_or_unstable = abs(pseq_lfc) >= 9.5 | abs(ortiz_lfc) >= 9.5,
    significant_this_study = is.finite(pseq_fdr) & pseq_fdr < fdr_cutoff,
    significant_ortiz = is.finite(ortiz_fdr) & ortiz_fdr < fdr_cutoff
  )

plot_data <- joined_all %>%
  filter(!capped_or_unstable, is.finite(pseq_lfc), is.finite(ortiz_lfc)) %>%
  mutate(
    significance_group = case_when(
      significant_this_study & significant_ortiz ~ "Significant in both",
      significant_ortiz ~ "Significant in Ortiz et al.",
      significant_this_study ~ "Significant in this study",
      TRUE ~ "Not significant"
    ),
    significance_group = factor(
      significance_group,
      levels = c("Not significant", "Significant in Ortiz et al.", "Significant in this study", "Significant in both")
    ),
    inherited_example_label = Gene %in% example_labels
  )

analysis_set <- plot_data %>%
  filter(significant_this_study, significant_ortiz)

pearson_estimable <- nrow(analysis_set) >= 3
sign_concordance <- if (nrow(analysis_set)) mean(sign(analysis_set$pseq_lfc) == sign(analysis_set$ortiz_lfc)) else NA_real_
if (pearson_estimable) {
  cor_result <- cor.test(
    analysis_set$pseq_lfc,
    analysis_set$ortiz_lfc,
    method = "pearson",
    alternative = "two.sided"
  )
  pearson_r <- unname(cor_result$estimate)
  pearson_p <- cor_result$p.value
  fit <- lm(ortiz_lfc ~ pseq_lfc, data = analysis_set)
  prediction_grid <- data.frame(
    pseq_lfc = seq(min(analysis_set$pseq_lfc), max(analysis_set$pseq_lfc), length.out = 200)
  )
  prediction <- predict(fit, newdata = prediction_grid, interval = "confidence", level = 0.95)
  prediction_grid$fit <- prediction[, "fit"]
  prediction_grid$lwr <- prediction[, "lwr"]
  prediction_grid$upr <- prediction[, "upr"]
} else {
  pearson_r <- NA_real_
  pearson_p <- NA_real_
  prediction_grid <- data.frame(pseq_lfc = numeric(), fit = numeric(), lwr = numeric(), upr = numeric())
}

primary_summary <- data.frame(
  model = "glmGamPoi with animal blocking",
  comparison = "L4/5 IT Foxp1_1 versus ST_2 compared with Ortiz et al. Foxp1 transgenic L4/5 IT",
  displayed_gene_rule = "nominal P < 0.05 in both datasets; capped or unstable |log2FC| >= 9.5 excluded",
  analysis_gene_rule = "BH FDR < 0.05 in both datasets",
  shared_nominal_before_stability_filter = nrow(joined_all),
  displayed_shared_nominal = nrow(plot_data),
  excluded_capped_or_unstable = sum(joined_all$capped_or_unstable),
  shared_fdr_supported = nrow(analysis_set),
  pearson_r = pearson_r,
  pearson_p_two_sided = pearson_p,
  sign_concordance = sign_concordance,
  stringsAsFactors = FALSE
)
write.csv(primary_summary, file.path(out_dir, "Foxp1_L45_consistency_summary.csv"), row.names = FALSE)
write.csv(plot_data %>% arrange(Gene), file.path(out_dir, "Foxp1_L45_consistency_plotdata.csv"), row.names = FALSE, na = "")
write.csv(analysis_set %>% arrange(Gene), file.path(out_dir, "Foxp1_L45_shared_FDR005.csv"), row.names = FALSE, na = "")
