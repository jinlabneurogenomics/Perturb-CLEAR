# Adnp morphology comparison between L2/3 IT and L4/5 IT neurons.


if (.Platform$OS.type == "windows") {
  r_minor <- strsplit(R.version[["minor"]], ".", fixed = TRUE)[[1]][1]
  user_library <- file.path(
    Sys.getenv("LOCALAPPDATA"), "R", "win-library",
    paste0(R.version[["major"]], ".", r_minor)
  )
  if (dir.exists(user_library)) .libPaths(c(user_library, .libPaths()))
}

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
})

if (!requireNamespace("nlme", quietly = TRUE)) stop("Package nlme is required")

set.seed(20260811)

# Inputs
# Command-line arguments, in order:
# 1. Morphometric measurements CSV.
# 2. Metric-label dictionary CSV.
# 3. Cell-metadata CSV.
# 4. Prespecified basal-feature list CSV.
# 5. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5L) {
  stop("Expected 5 command-line arguments in the order documented above.")
}
full_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
dictionary_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
metadata_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
metric_file <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
out_dir <- args[[5]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
for (path in c(full_file, dictionary_file, metadata_file, metric_file)) {
  if (!file.exists(path)) stop("Missing input: ", path)
}

fdr_symbol <- function(x) {
  case_when(
    is.na(x) ~ "n.s.",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

# Inputs and exact metric nomenclature
metrics <- read_csv(metric_file, show_col_types = FALSE)$Metrics
if (length(metrics) != 19) stop("Expected exactly 19 metrics, found ", length(metrics))

dictionary <- read_csv(dictionary_file, show_col_types = FALSE) %>%
  distinct(Old, .keep_all = TRUE)
label_lookup <- setNames(dictionary$New, dictionary$Old)
if (any(!metrics %in% names(label_lookup))) {
  stop("Dictionary is missing: ", paste(metrics[!metrics %in% names(label_lookup)], collapse = ", "))
}

block_lookup <- c(
  Basal_dendrite_total_length = "Size",
  Basal_branch_len_mean = "Size",
  Basal_extent_height = "Size",
  Basal_extent_depth = "Size",
  Basal_extent_width = "Size",
  Basal_alpha_3D_volume = "Size",
  Basal_convex_3D_volume = "Size",
  Basal_n_branches = "Complexity",
  Basal_n_terminals = "Complexity",
  Basal_n_bifurcations = "Complexity",
  Basal_sholl_peak = "Complexity",
  Basal_max_branch_order = "Complexity",
  Basal_sholl_auc = "Complexity",
  Basal_sholl_entropy = "Complexity",
  Basal_sholl_radius_at_peak = "Complexity",
  Basal_sholl_radius_com = "Complexity",
  Basal_bif_plane_torque_mean = "Angle",
  Basal_bif_parent_child_tilt_min_mean = "Angle",
  Basal_bif_childpair_angle_min_mean = "Angle"
)
if (any(!metrics %in% names(block_lookup))) stop("Block lookup is incomplete")

full_data <- read_csv(full_file, show_col_types = FALSE, name_repair = "unique")
metadata <- read_csv(metadata_file, show_col_types = FALSE) %>%
  select(Neuron, Label) %>%
  distinct(Neuron, .keep_all = TRUE)

dat <- full_data %>%
  left_join(metadata, by = c("file" = "Neuron")) %>%
  filter(Label %in% c("L2/3 IT", "L4/5 IT"), Perturbation %in% c("ctrl", "Adnp")) %>%
  transmute(
    Neuron = file,
    Brain = factor(Brain),
    CellType = factor(Label, levels = c("L2/3 IT", "L4/5 IT")),
    Group = factor(Perturbation, levels = c("ctrl", "Adnp")),
    across(all_of(metrics), as.numeric)
  )

if (nrow(dat) != 243) stop("Expected 243 target neurons, found ", nrow(dat))
if (any(!is.finite(as.matrix(dat[metrics])))) stop("Non-finite values remain in the 19-metric matrix")

sample_counts <- dat %>%
  count(CellType, Group, Brain, name = "n_neurons") %>%
  arrange(CellType, Group, Brain)
write_csv(sample_counts, file.path(out_dir, "sample_counts_by_brain.csv"))

# Standardize each metric to the matching cell-type control distribution ------
metric_long <- dat %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Raw_value") %>%
  group_by(CellType, Metric) %>%
  mutate(
    control_mean = mean(Raw_value[Group == "ctrl"]),
    control_sd = sd(Raw_value[Group == "ctrl"]),
    Z = (Raw_value - control_mean) / control_sd
  ) %>%
  ungroup() %>%
  mutate(
    Metric_label = unname(label_lookup[Metric]),
    Block = factor(unname(block_lookup[Metric]), levels = c("Size", "Complexity", "Angle"))
  )

write_csv(metric_long, file.path(out_dir, "metric_neuron_data.csv"))

# Mixed models
fit_term <- function(d, formula, term_pattern) {
  fit <- tryCatch(
    nlme::lme(
      formula, random = ~1 | Brain, data = d, method = "REML",
      na.action = na.omit,
      control = nlme::lmeControl(returnObject = TRUE, msMaxIter = 200)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    lm_fit <- lm(formula, data = d)
    tab <- summary(lm_fit)$coefficients
    term <- grep(term_pattern, rownames(tab), value = TRUE)[1]
    return(tibble(
      estimate = tab[term, "Estimate"], SE = tab[term, "Std. Error"],
      DF = df.residual(lm_fit), p_value = tab[term, "Pr(>|t|)"], model = "lm fallback"
    ))
  }
  tab <- summary(fit)$tTable
  term <- grep(term_pattern, rownames(tab), value = TRUE)[1]
  tibble(
    estimate = tab[term, "Value"], SE = tab[term, "Std.Error"],
    DF = tab[term, "DF"], p_value = tab[term, "p-value"], model = "lme: random intercept for Brain"
  )
}

metric_effects <- metric_long %>%
  group_by(CellType, Metric, Metric_label, Block) %>%
  group_modify(~ fit_term(.x, Z ~ Group, "^GroupAdnp$")) %>%
  ungroup() %>%
  group_by(CellType) %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    lower = estimate - qt(0.975, pmax(DF, 1)) * SE,
    upper = estimate + qt(0.975, pmax(DF, 1)) * SE,
    Significance = fdr_symbol(FDR)
  ) %>%
  ungroup()

metric_interactions <- metric_long %>%
  group_by(Metric, Metric_label, Block) %>%
  group_modify(~ fit_term(
    .x, Z ~ CellType * Group,
    "CellType.*:GroupAdnp|GroupAdnp:CellType"
  )) %>%
  ungroup() %>%
  mutate(
    FDR = p.adjust(p_value, method = "BH"),
    lower = estimate - qt(0.975, pmax(DF, 1)) * SE,
    upper = estimate + qt(0.975, pmax(DF, 1)) * SE,
    Significance = fdr_symbol(FDR),
    Contrast = "Adnp effect in L4/5 IT minus Adnp effect in L2/3 IT"
  )

write_csv(metric_effects, file.path(out_dir, "metric_celltype_effects.csv"))
write_csv(metric_interactions, file.path(out_dir, "metric_celltype_interactions.csv"))
