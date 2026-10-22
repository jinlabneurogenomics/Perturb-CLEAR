# Adnp morphology programs from consensus non-negative matrix factorization.

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
  library(nlme)
})

analysis_seed <- 20261009L
set.seed(analysis_seed)

# Inputs
# Command-line arguments, in order:
# 1. Morphometric measurements CSV.
# 2. Metric-label dictionary CSV.
# 3. Cell-metadata CSV.
# 4. Prespecified basal-feature list CSV.
# 5. Prior program-loading CSV used only to align program names.
# 6. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 6L) {
  stop("Expected 6 command-line arguments in the order documented above.")
}
full_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
dictionary_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
metadata_file <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
metric_file <- normalizePath(args[[4]], winslash = "/", mustWork = TRUE)
prior_loading_file <- normalizePath(args[[5]], winslash = "/", mustWork = TRUE)
output_dir <- args[[6]]
data_dir <- file.path(output_dir, "data")
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

input_files <- c(full_file, dictionary_file, metadata_file, metric_file, prior_loading_file)
if (any(!file.exists(input_files))) {
  stop("Missing inputs: ", paste(input_files[!file.exists(input_files)], collapse = ", "))
}
# Fixed settings and user-facing labels
candidate_ranks <- 2:5
selected_rank <- 3L
n_resamples <- 20L
n_per_cell_class <- 27L
n_random_starts <- 12L
nmf_max_iter <- 5000L
projection_max_iter <- 15000L
convergence_tolerance <- 1e-6

program_biological_names <- c(
  "Program 1" = "Branch complexity",
  "Program 2" = "Branch order and geometry",
  "Program 3" = "Arbor scale and radial extent"
)
fdr_symbol <- function(x) {
  case_when(
    is.na(x) ~ "n.s.",
    x < 0.001 ~ "***",
    x < 0.01 ~ "**",
    x < 0.05 ~ "*",
    TRUE ~ "n.s."
  )
}

# Input validation
metrics <- read_csv(metric_file, show_col_types = FALSE)$Metrics
if (length(metrics) != 19L || anyDuplicated(metrics)) {
  stop("Expected 19 unique prespecified metrics")
}

dictionary <- read_csv(dictionary_file, show_col_types = FALSE) %>%
  distinct(Old, .keep_all = TRUE)
label_lookup <- setNames(dictionary$New, dictionary$Old)
if (any(!metrics %in% names(label_lookup))) {
  stop("Metric dictionary is missing: ", paste(setdiff(metrics, names(label_lookup)), collapse = ", "))
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
  Basal_bif_plane_torque_mean = "Geometry",
  Basal_bif_parent_child_tilt_min_mean = "Geometry",
  Basal_bif_childpair_angle_min_mean = "Geometry"
)
if (any(!metrics %in% names(block_lookup))) stop("Metric block mapping is incomplete")

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

if (nrow(dat) != 243L) stop("Expected 243 target neurons, found ", nrow(dat))
if (anyDuplicated(dat$Neuron)) stop("Neuron identifiers are not unique")
if (any(!is.finite(as.matrix(dat[metrics])))) stop("Non-finite metric values remain")
if (any(as.matrix(dat[metrics]) < 0)) stop("NMF inputs must be non-negative")
if (any(is.na(dat$Brain)) || any(is.na(dat$CellType)) || any(is.na(dat$Group))) {
  stop("Missing experimental metadata remain")
}

sample_counts <- dat %>%
  count(CellType, Group, Brain, name = "N_neurons") %>%
  arrange(CellType, Group, Brain)
write_csv(sample_counts, file.path(data_dir, "sample_counts_by_brain.csv"))

sample_summary <- dat %>%
  group_by(CellType, Group) %>%
  summarise(N_neurons = n(), N_brains = n_distinct(Brain), .groups = "drop")
write_csv(sample_summary, file.path(data_dir, "sample_summary.csv"))

# Matching-control scaling
raw_matrix <- as.matrix(dat[metrics])
log_matrix <- log1p(raw_matrix)
nmf_matrix <- log_matrix
nmf_control_means <- list()
for (ct in levels(dat$CellType)) {
  rows <- dat$CellType == ct
  ctrl_rows <- rows & dat$Group == "ctrl"
  ref <- colMeans(log_matrix[ctrl_rows, , drop = FALSE])
  if (any(!is.finite(ref)) || any(ref <= 0)) stop("Invalid NMF control reference for ", ct)
  nmf_matrix[rows, ] <- sweep(log_matrix[rows, , drop = FALSE], 2, ref, FUN = "/")
  nmf_control_means[[ct]] <- ref
}
if (any(!is.finite(nmf_matrix)) || any(nmf_matrix < 0)) stop("Invalid scaled NMF matrix")

nmf_input_export <- bind_cols(
  dat %>% select(Neuron, Brain, CellType, Group),
  as.data.frame(nmf_matrix)
)
write_csv(nmf_input_export, file.path(data_dir, "nmf_preprocessed_matrix.csv"))

metric_long <- dat %>%
  pivot_longer(all_of(metrics), names_to = "Metric", values_to = "Raw_value") %>%
  group_by(CellType, Metric) %>%
  mutate(
    control_mean = mean(Raw_value[Group == "ctrl"]),
    control_sd = sd(Raw_value[Group == "ctrl"]),
    Metric_z = (Raw_value - control_mean) / control_sd
  ) %>%
  ungroup() %>%
  mutate(
    Metric_label = unname(label_lookup[Metric]),
    Block = unname(block_lookup[Metric])
  )
if (any(!is.finite(metric_long$Metric_z))) stop("Invalid matching-control metric z scores")
write_csv(metric_long, file.path(data_dir, "metric_neuron_data.csv"))

# Corrected, asserted control resampling
sample_one_value <- function(x) {
  if (!length(x)) stop("Cannot sample from an empty vector")
  x[sample.int(length(x), size = 1L)]
}

sample_l23_controls <- function(seed, n_target = n_per_cell_class) {
  set.seed(seed)
  candidates <- which(dat$CellType == "L2/3 IT" & dat$Group == "ctrl")
  by_brain <- split(candidates, droplevels(dat$Brain[candidates]))
  first <- vapply(by_brain, sample_one_value, integer(1))
  if (length(first) > n_target) stop("More L2/3 control brains than target neurons")
  remaining <- setdiff(candidates, first)
  n_fill <- n_target - length(first)
  fill <- if (n_fill > 0L) remaining[sample.int(length(remaining), n_fill)] else integer(0)
  rows <- c(unname(first), fill)
  if (length(rows) != n_target || anyDuplicated(rows)) stop("Invalid L2/3 control sample")
  if (any(dat$CellType[rows] != "L2/3 IT") || any(dat$Group[rows] != "ctrl")) {
    stop("L2/3 sampling contamination detected")
  }
  rows
}

l45_ctrl <- which(dat$CellType == "L4/5 IT" & dat$Group == "ctrl")
if (length(l45_ctrl) != n_per_cell_class || anyDuplicated(l45_ctrl)) {
  stop("Expected 27 unique L4/5 control neurons")
}

training_rows <- vector("list", n_resamples)
training_membership <- vector("list", n_resamples)
for (resample_id in seq_len(n_resamples)) {
  l23_rows <- sample_l23_controls(analysis_seed + resample_id)
  rows <- c(l23_rows, l45_ctrl)
  if (length(rows) != 2L * n_per_cell_class || anyDuplicated(rows)) {
    stop("Training resample contains duplicate or missing rows: ", resample_id)
  }
  if (any(dat$Group[rows] != "ctrl")) stop("Adnp neuron entered training resample ", resample_id)
  tab <- table(dat$CellType[rows])
  if (any(tab != n_per_cell_class)) stop("Unbalanced cell classes in resample ", resample_id)
  training_rows[[resample_id]] <- rows
  training_membership[[resample_id]] <- tibble(
    Resample = resample_id,
    Neuron = dat$Neuron[rows],
    Brain = as.character(dat$Brain[rows]),
    CellType = as.character(dat$CellType[rows]),
    Group = as.character(dat$Group[rows])
  )
}
training_membership <- bind_rows(training_membership)

membership_audit <- training_membership %>%
  group_by(Resample) %>%
  summarise(
    N_rows = n(),
    N_unique_neurons = n_distinct(Neuron),
    N_L23 = sum(CellType == "L2/3 IT"),
    N_L45 = sum(CellType == "L4/5 IT"),
    N_ctrl = sum(Group == "ctrl"),
    N_Adnp = sum(Group == "Adnp"),
    N_brains = n_distinct(Brain),
    .groups = "drop"
  )
if (any(membership_audit$N_rows != 54L) ||
    any(membership_audit$N_unique_neurons != 54L) ||
    any(membership_audit$N_L23 != 27L) ||
    any(membership_audit$N_L45 != 27L) ||
    any(membership_audit$N_ctrl != 54L) ||
    any(membership_audit$N_Adnp != 0L)) {
  stop("Training membership audit failed")
}
write_csv(training_membership, file.path(data_dir, "nmf_training_membership.csv"))
write_csv(membership_audit, file.path(data_dir, "nmf_training_membership_audit.csv"))

# NMF implementation
# Euclidean multiplicative updates follow Lee and Seung (NIPS 2000).
nmf_fit_mu <- function(v, rank, n_run, max_iter, tolerance, seed) {
  eps <- 1e-10
  best <- NULL
  best_error <- Inf
  for (run_id in seq_len(n_run)) {
    set.seed(seed + run_id + rank * 10000L)
    w <- matrix(runif(nrow(v) * rank, 0.2, 1), nrow(v), rank)
    h <- matrix(runif(rank * ncol(v), 0.2, 1), rank, ncol(v))
    previous_error <- Inf
    converged <- FALSE
    final_iter <- max_iter
    for (iter in seq_len(max_iter)) {
      h <- h * (t(w) %*% v) / (t(w) %*% w %*% h + eps)
      w <- w * (v %*% t(h)) / (w %*% h %*% t(h) + eps)
      if (iter %% 25L == 0L || iter == max_iter) {
        current_error <- sqrt(sum((v - w %*% h)^2) / sum(v^2))
        relative_change <- abs(previous_error - current_error) / max(previous_error, eps)
        if (is.finite(relative_change) && relative_change < tolerance) {
          converged <- TRUE
          final_iter <- iter
          break
        }
        previous_error <- current_error
      }
    }
    err <- sqrt(sum((v - w %*% h)^2) / sum(v^2))
    if (err < best_error) {
      h_sums <- rowSums(h)
      best <- list(
        w = sweep(w, 2, h_sums, FUN = "*"),
        h = sweep(h, 1, h_sums, FUN = "/"),
        error = err,
        converged = converged,
        iterations = final_iter,
        selected_start = run_id
      )
      best_error <- err
    }
  }
  best
}

nmf_project <- function(v, h, max_iter, tolerance) {
  eps <- 1e-10
  w <- matrix(1, nrow(v), nrow(h))
  previous_error <- Inf
  converged <- FALSE
  final_iter <- max_iter
  for (iter in seq_len(max_iter)) {
    w <- w * (v %*% t(h)) / (w %*% h %*% t(h) + eps)
    if (iter %% 25L == 0L || iter == max_iter) {
      current_error <- sqrt(sum((v - w %*% h)^2) / sum(v^2))
      relative_change <- abs(previous_error - current_error) / max(previous_error, eps)
      if (is.finite(relative_change) && relative_change < tolerance) {
        converged <- TRUE
        final_iter <- iter
        break
      }
      previous_error <- current_error
    }
  }
  list(w = w, converged = converged, iterations = final_iter)
}

all_permutations <- function(values) {
  if (length(values) == 1L) return(matrix(values, nrow = 1L))
  do.call(rbind, lapply(seq_along(values), function(i) {
    cbind(values[i], all_permutations(values[-i]))
  }))
}

best_alignment <- function(reference_h, target_h) {
  rank <- nrow(reference_h)
  permutations <- all_permutations(seq_len(rank))
  total_scores <- apply(permutations, 1, function(perm) {
    sum(vapply(seq_len(rank), function(k) {
      cor(reference_h[k, ], target_h[perm[k], ], method = "pearson")
    }, numeric(1)))
  })
  best_perm <- permutations[which.max(total_scores), ]
  aligned <- target_h[best_perm, , drop = FALSE]
  correlations <- vapply(seq_len(rank), function(k) {
    cor(reference_h[k, ], aligned[k, ], method = "pearson")
  }, numeric(1))
  list(permutation = best_perm, aligned_h = aligned, correlations = correlations)
}

align_fit_collection <- function(fits, rank) {
  n <- length(fits)
  pairwise_similarity <- matrix(1, n, n)
  for (i in seq_len(n - 1L)) {
    for (j in (i + 1L):n) {
      match <- best_alignment(fits[[i]]$h, fits[[j]]$h)
      pairwise_similarity[i, j] <- mean(match$correlations)
      pairwise_similarity[j, i] <- pairwise_similarity[i, j]
    }
  }
  reference_index <- which.max(rowMeans(pairwise_similarity))
  reference_h <- fits[[reference_index]]$h
  aligned_array <- array(NA_real_, dim = c(n, rank, ncol(reference_h)))
  stability <- vector("list", n)
  for (i in seq_len(n)) {
    match <- best_alignment(reference_h, fits[[i]]$h)
    aligned_array[i, , ] <- match$aligned_h
    stability[[i]] <- tibble(
      Resample = i,
      Rank = rank,
      Reference_resample = reference_index,
      Reference_program = seq_len(rank),
      Aligned_program = match$permutation,
      Loading_correlation = match$correlations,
      Reconstruction_error = fits[[i]]$error,
      Converged = fits[[i]]$converged,
      Iterations = fits[[i]]$iterations,
      Selected_start = fits[[i]]$selected_start
    )
  }
  consensus_h <- apply(aligned_array, c(2, 3), mean)
  consensus_h <- sweep(consensus_h, 1, rowSums(consensus_h), FUN = "/")
  list(
    consensus_h = consensus_h,
    stability = bind_rows(stability),
    reference_resample = reference_index,
    pairwise_similarity = pairwise_similarity
  )
}

# Fit ranks 2-5 using the same corrected balanced resamples
fits_by_rank <- setNames(vector("list", length(candidate_ranks)), candidate_ranks)
aligned_by_rank <- setNames(vector("list", length(candidate_ranks)), candidate_ranks)

for (rank_value in candidate_ranks) {
  message("Fitting rank ", rank_value, " across ", n_resamples, " corrected resamples")
  fits <- vector("list", n_resamples)
  for (resample_id in seq_len(n_resamples)) {
    rows <- training_rows[[resample_id]]
    fits[[resample_id]] <- nmf_fit_mu(
      nmf_matrix[rows, , drop = FALSE],
      rank = rank_value,
      n_run = n_random_starts,
      max_iter = nmf_max_iter,
      tolerance = convergence_tolerance,
      seed = analysis_seed + rank_value * 100000L + resample_id * 100L
    )
  }
  fits_by_rank[[as.character(rank_value)]] <- fits
  aligned_by_rank[[as.character(rank_value)]] <- align_fit_collection(fits, rank_value)
}

rank_fit_data <- bind_rows(lapply(candidate_ranks, function(rank_value) {
  aligned_by_rank[[as.character(rank_value)]]$stability
}))
write_csv(rank_fit_data, file.path(data_dir, "rank_resample_fit_and_stability.csv"))

rank_summary <- rank_fit_data %>%
  group_by(Rank) %>%
  summarise(
    Median_error = median(Reconstruction_error),
    Error_Q1 = quantile(Reconstruction_error, 0.25),
    Error_Q3 = quantile(Reconstruction_error, 0.75),
    Median_component_stability = median(Loading_correlation),
    Minimum_component_median_stability = min(
      tapply(Loading_correlation, Reference_program, median)
    ),
    Fraction_converged = mean(Converged),
    .groups = "drop"
  ) %>%
  mutate(
    Prespecified_primary = Rank == selected_rank,
    Error_improvement_from_previous = lag(Median_error) - Median_error,
    Relative_error_improvement_from_previous = Error_improvement_from_previous / lag(Median_error)
  )
write_csv(rank_summary, file.path(data_dir, "rank_summary.csv"))

# Final rank-3 consensus and continuity mapping
consensus_h_unordered <- aligned_by_rank[[as.character(selected_rank)]]$consensus_h
colnames(consensus_h_unordered) <- metrics

prior_loadings <- read_csv(prior_loading_file, show_col_types = FALSE) %>%
  select(Program, Metric, Loading) %>%
  filter(Program %in% names(program_biological_names), Metric %in% metrics) %>%
  pivot_wider(names_from = Metric, values_from = Loading)
prior_programs <- prior_loadings$Program
prior_h <- as.matrix(prior_loadings[metrics])
rownames(prior_h) <- prior_programs

continuity_match <- best_alignment(prior_h, consensus_h_unordered)
consensus_h <- continuity_match$aligned_h
rownames(consensus_h) <- names(program_biological_names)
colnames(consensus_h) <- metrics

label_alignment <- tibble(
  Program = names(program_biological_names),
  Proposed_name = unname(program_biological_names),
  Corrected_unordered_component = continuity_match$permutation,
  Correlation_with_superseded_loading_profile = continuity_match$correlations,
  Label_source = "Inherited from prior manuscript labels after loading-profile matching; post hoc and not independently validated"
)
write_csv(label_alignment, file.path(data_dir, "program_label_alignment.csv"))

projection <- nmf_project(
  nmf_matrix, consensus_h,
  max_iter = projection_max_iter,
  tolerance = convergence_tolerance
)
if (!projection$converged) warning("Projection reached the iteration limit")
nmf_scores <- projection$w
colnames(nmf_scores) <- names(program_biological_names)

# Loading tables and row order
loading_data <- as.data.frame(consensus_h) %>%
  tibble::rownames_to_column("Program") %>%
  pivot_longer(-Program, names_to = "Metric", values_to = "Loading") %>%
  mutate(
    Metric_label = unname(label_lookup[Metric]),
    Block = unname(block_lookup[Metric]),
    Program = factor(Program, levels = names(program_biological_names))
  ) %>%
  group_by(Program) %>%
  mutate(Relative_loading = Loading / max(Loading)) %>%
  ungroup()

profile_matrix <- loading_data %>%
  select(Metric, Program, Relative_loading) %>%
  mutate(Program = as.character(Program)) %>%
  pivot_wider(names_from = Program, values_from = Relative_loading) %>%
  arrange(Metric)
profile_values <- as.matrix(profile_matrix[names(program_biological_names)])
rownames(profile_values) <- profile_matrix$Metric
profile_cor <- cor(t(profile_values), method = "pearson")
profile_cor[!is.finite(profile_cor)] <- 0
diag(profile_cor) <- 1
metric_hclust <- hclust(as.dist(1 - profile_cor), method = "average")
nmf_metric_order <- rownames(profile_values)[metric_hclust$order]

loading_data <- loading_data %>%
  mutate(
    Metric_order = match(Metric, nmf_metric_order),
    Proposed_name = unname(program_biological_names[as.character(Program)])
  ) %>%
  arrange(Program, Metric_order)
write_csv(loading_data, file.path(data_dir, "program_loadings.csv"))

loading_matrix <- loading_data %>%
  select(Metric, Metric_label, Metric_order, Program, Loading) %>%
  pivot_wider(names_from = Program, values_from = Loading) %>%
  arrange(Metric_order) %>%
  select(-Metric_order)
write_csv(loading_matrix, file.path(data_dir, "program_loading_matrix.csv"))

metric_order_table <- tibble(
  Metric = nmf_metric_order,
  Metric_label = unname(label_lookup[nmf_metric_order]),
  Plot_order_top_to_bottom = seq_along(nmf_metric_order)
)
write_csv(metric_order_table, file.path(data_dir, "metric_plot_order.csv"))

program_definitions <- loading_data %>%
  arrange(Program, desc(Loading)) %>%
  group_by(Program) %>%
  slice_head(n = 5) %>%
  summarise(
    Proposed_name = first(Proposed_name),
    Top_metrics = paste(Metric_label, collapse = "; "),
    .groups = "drop"
  ) %>%
  left_join(label_alignment, by = c("Program", "Proposed_name"))
write_csv(program_definitions, file.path(data_dir, "program_definitions.csv"))

# Scores and mixed models
score_data <- bind_cols(
  dat %>% select(Neuron, Brain, CellType, Group),
  as.data.frame(nmf_scores)
) %>%
  pivot_longer(
    all_of(names(program_biological_names)),
    names_to = "Program", values_to = "Score"
  ) %>%
  mutate(Program = factor(Program, levels = names(program_biological_names))) %>%
  group_by(CellType, Program) %>%
  mutate(
    control_mean = mean(Score[Group == "ctrl"]),
    control_sd = sd(Score[Group == "ctrl"]),
    Score_z = (Score - control_mean) / control_sd
  ) %>%
  ungroup()
if (any(!is.finite(score_data$Score_z))) stop("Invalid standardized NMF scores")
write_csv(score_data, file.path(data_dir, "neuron_program_scores.csv"))

fit_lme_term <- function(d, formula, term_pattern) {
  fit <- nlme::lme(
    fixed = formula,
    random = ~1 | Brain,
    data = d,
    method = "REML",
    na.action = na.fail,
    control = nlme::lmeControl(returnObject = FALSE, msMaxIter = 300)
  )
  tab <- summary(fit)$tTable
  terms <- grep(term_pattern, rownames(tab), value = TRUE)
  if (length(terms) != 1L) stop("Expected one model term matching ", term_pattern)
  term <- terms[[1]]
  tibble(
    Term = term,
    Estimate = tab[term, "Value"],
    SE = tab[term, "Std.Error"],
    DF = tab[term, "DF"],
    P_value = tab[term, "p-value"],
    Model = "nlme::lme; random intercept for Brain; REML"
  )
}

score_effects <- score_data %>%
  group_by(CellType, Program) %>%
  group_modify(~ fit_lme_term(.x, Score_z ~ Group, "^GroupAdnp$")) %>%
  ungroup() %>%
  group_by(CellType) %>%
  mutate(FDR_within_celltype_3 = p.adjust(P_value, method = "BH")) %>%
  ungroup() %>%
  mutate(
    FDR_figure_family_6 = p.adjust(P_value, method = "BH"),
    CI_lower = Estimate - qt(0.975, pmax(DF, 1)) * SE,
    CI_upper = Estimate + qt(0.975, pmax(DF, 1)) * SE,
    Significance = fdr_symbol(FDR_figure_family_6),
    Primary_FDR_family = "All six within-cell-class program comparisons in the score figure"
  )
write_csv(score_effects, file.path(data_dir, "celltype_program_effects.csv"))

score_interactions <- score_data %>%
  group_by(Program) %>%
  group_modify(~ fit_lme_term(
    .x,
    Score_z ~ CellType * Group,
    "CellTypeL4/5 IT:GroupAdnp|GroupAdnp:CellTypeL4/5 IT"
  )) %>%
  ungroup() %>%
  mutate(
    FDR_interaction_family_3 = p.adjust(P_value, method = "BH"),
    CI_lower = Estimate - qt(0.975, pmax(DF, 1)) * SE,
    CI_upper = Estimate + qt(0.975, pmax(DF, 1)) * SE,
    Significance = fdr_symbol(FDR_interaction_family_3),
    Interpretation = "Adnp effect in L4/5 IT minus Adnp effect in L2/3 IT"
  )
write_csv(score_interactions, file.path(data_dir, "celltype_interactions.csv"))

score_group_summary <- score_data %>%
  group_by(CellType, Group, Program) %>%
  summarise(
    Mean_score_z = mean(Score_z),
    Median_score_z = median(Score_z),
    SD = sd(Score_z),
    SEM = SD / sqrt(n()),
    N_neurons = n(),
    N_brains = n_distinct(Brain),
    .groups = "drop"
  )
write_csv(score_group_summary, file.path(data_dir, "program_score_group_summary.csv"))

# Pearson-R interpretation matrices
score_metric_data <- inner_join(
  score_data %>% select(Neuron, Brain, CellType, Group, Program, Score_z),
  metric_long %>% select(Neuron, Brain, CellType, Group, Metric, Metric_label, Metric_z),
  by = c("Neuron", "Brain", "CellType", "Group"),
  relationship = "many-to-many"
)

safe_cor <- function(x, y) {
  out <- suppressWarnings(cor(x, y, method = "pearson", use = "complete.obs"))
  ifelse(is.finite(out), out, NA_real_)
}

balanced_correlation_resamples <- training_membership %>%
  inner_join(
    score_metric_data,
    by = c("Neuron", "Brain", "CellType", "Group"),
    relationship = "many-to-many"
  ) %>%
  group_by(Resample, Program, Metric, Metric_label) %>%
  summarise(
    Pearson_R = safe_cor(Score_z, Metric_z),
    N_neurons = n_distinct(Neuron),
    N_brains = n_distinct(Brain),
    .groups = "drop"
  )
if (any(balanced_correlation_resamples$N_neurons != 54L) ||
    any(!is.finite(balanced_correlation_resamples$Pearson_R))) {
  stop("Balanced control correlation audit failed")
}
write_csv(
  balanced_correlation_resamples,
  file.path(data_dir, "Pearson_R_balanced_control_resamples.csv")
)

balanced_correlation <- balanced_correlation_resamples %>%
  group_by(Program, Metric, Metric_label) %>%
  summarise(
    Pearson_R = tanh(mean(atanh(pmax(pmin(Pearson_R, 0.999999), -0.999999)))),
    Resample_SD = sd(Pearson_R),
    Resample_min = min(Pearson_R),
    Resample_max = max(Pearson_R),
    N_resamples = n(),
    N_neurons_per_resample = first(N_neurons),
    .groups = "drop"
  ) %>%
  mutate(
    Cohort = "Balanced controls",
    Definition = "Fisher-z mean of Pearson correlations across control resamples"
  )

all_neuron_correlation <- score_metric_data %>%
  group_by(Program, Metric, Metric_label) %>%
  summarise(
    Pearson_R = safe_cor(Score_z, Metric_z),
    N_neurons = n_distinct(Neuron),
    N_brains = n_distinct(Brain),
    .groups = "drop"
  ) %>%
  mutate(
    Cohort = "All neurons",
    Definition = "Pearson correlations across control and Adnp neurons",
    Resample_SD = NA_real_,
    Resample_min = NA_real_,
    Resample_max = NA_real_,
    N_resamples = NA_integer_,
    N_neurons_per_resample = NA_integer_
  )

correlation_data <- bind_rows(
  balanced_correlation %>% mutate(N_neurons = N_neurons_per_resample, N_brains = NA_integer_),
  all_neuron_correlation
) %>%
  mutate(
    Program = factor(Program, levels = names(program_biological_names)),
    Proposed_name = unname(program_biological_names[as.character(Program)]),
    Metric_order = match(Metric, nmf_metric_order)
  ) %>%
  arrange(Cohort, Program, Metric_order)
write_csv(correlation_data, file.path(data_dir, "program_metric_Pearson_R.csv"))

write_correlation_matrix <- function(data, filename) {
  out <- data %>%
    select(Metric, Metric_label, Metric_order, Program, Pearson_R) %>%
    pivot_wider(names_from = Program, values_from = Pearson_R) %>%
    arrange(Metric_order) %>%
    select(-Metric_order)
  write_csv(out, file.path(data_dir, filename))
}
write_correlation_matrix(
  correlation_data %>% filter(Cohort == "Balanced controls"),
  "Pearson_R_matrix.csv"
)
write_correlation_matrix(
  correlation_data %>% filter(Cohort == "All neurons"),
  "Pearson_R_all_neurons_sensitivity_matrix.csv"
)