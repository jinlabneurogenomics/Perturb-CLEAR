# Compare Adnp effect estimates between the pooled screen and validation cohort.
# The script exports joined estimates, Pearson correlations, and linear-model summaries.
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

# Command-line arguments, in order:
# 1. Screen-versus-validation effect-size CSV.
# 2. Output-file prefix.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Expected 2 command-line arguments in the order documented above.")
}
infile <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
out_prefix <- args[[2]]
# --- parameters you may want to edit ---
padj_cutoff <- 0.05
pval_floor  <- 1e-300  # prevents -Inf when p=0

# --- read data ---
df <- readr::read_csv(infile, show_col_types = FALSE)

required_cols <- c(
  "Metric",
  "Estimate_valid", "Estimate_pool",
  "t.value_valid",  "t.value_pool",
  "Pr...t.._valid",
  "padj_byTest_valid", "padj_byTest_pool"
)
missing_cols <- setdiff(required_cols, colnames(df))
if (length(missing_cols) > 0) {
  stop(paste0(
    "Missing required columns: ", paste(missing_cols, collapse = ", "), "\n",
    "Found columns: ", paste(colnames(df), collapse = ", ")
  ))
}

# derive aesthetics
df2 <- df %>%
  mutate(
    sig_validation = padj_byTest_valid < padj_cutoff,
    sig_pool       = padj_byTest_pool  < padj_cutoff,
    sig_group = case_when(
      sig_validation & sig_pool      ~ "Significant in both",
      sig_validation != sig_pool     ~ "Significant in one",
      TRUE                           ~ "Not significant"
    ),
    neglog10_pval_validation = -log10(pmax(Pr...t.._valid, pval_floor))
  )

# long format for plotting (one row per metric per stat)
# NOTE: "valid" becomes the y-axis column "validation"
plot_df <- df2 %>%
  select(
    Metric,
    sig_group,
    neglog10_pval_validation,
    Pr...t.._valid,
    padj_byTest_valid,
    padj_byTest_pool,
    Estimate_valid, Estimate_pool,
    t.value_valid,  t.value_pool
  ) %>%
  pivot_longer(
    cols = c(Estimate_valid, Estimate_pool, t.value_valid, t.value_pool),
    names_to = c("stat", "sample"),
    names_pattern = "^(Estimate|t\\.value)_(valid|pool)$",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = sample, values_from = value) %>%
  rename(validation = valid) %>%
  mutate(
    stat = recode(
      stat,
      "Estimate" = "Estimate (Cohen's d)",
      "t.value"  = "t-statistic (t.value)"
    )
  )

# Pearson correlation summary
cor_summary <- plot_df %>%
  group_by(stat) %>%
  summarise(
    n = sum(complete.cases(pool, validation)),
    pearson_r = cor(pool, validation, use = "complete.obs", method = "pearson"),
    pearson_p = cor.test(pool, validation, method = "pearson")$p.value,
    .groups = "drop"
  )

# linear regression parameters: validation ~ pool (exact coefficients + 95% CI)
lm_models <- plot_df %>%
  filter(complete.cases(pool, validation)) %>%
  group_by(stat) %>%
  summarise(
    model = list(lm(validation ~ pool, data = cur_data())),
    .groups = "drop"
  )

lm_params <- lm_models %>%
  mutate(
    intercept = sapply(model, function(m) unname(coef(m)[1])),
    slope     = sapply(model, function(m) unname(coef(m)[2])),
    r_squared = sapply(model, function(m) summary(m)$r.squared),
    
    intercept_se = sapply(model, function(m) summary(m)$coefficients[1, 2]),
    slope_se     = sapply(model, function(m) summary(m)$coefficients[2, 2]),
    
    intercept_p  = sapply(model, function(m) summary(m)$coefficients[1, 4]),
    slope_p      = sapply(model, function(m) summary(m)$coefficients[2, 4]),
    
    intercept_ci_low  = sapply(model, function(m) confint(m)[1, 1]),
    intercept_ci_high = sapply(model, function(m) confint(m)[1, 2]),
    slope_ci_low      = sapply(model, function(m) confint(m)[2, 1]),
    slope_ci_high     = sapply(model, function(m) confint(m)[2, 2])
  ) %>%
  select(-model)

# Combine Pearson + LM outputs into one summary table
summary_df <- cor_summary %>%
  left_join(lm_params, by = "stat")

# Write outputs (data + correlation/regression summary)
readr::write_csv(plot_df, paste0(out_prefix, "_joined_effects.csv"))
readr::write_csv(summary_df, paste0(out_prefix, "_correlation_summary.csv"))
