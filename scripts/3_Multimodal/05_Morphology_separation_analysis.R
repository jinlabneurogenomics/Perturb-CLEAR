# Exploratory morphology separation of RNA subpopulations in published neurons.

options(stringsAsFactors = FALSE)

# Command-line arguments, in order:
# 1. Published-neuron morphology-separation table, optionally gzip-compressed.
# 2. Feature-manifest CSV.
# 3. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Expected 3 command-line arguments in the order documented above.")
}

feature_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
manifest_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[3]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20261001L
n_bootstrap <- 20000L
n_randomization <- 100000L
max_exact <- 50000L
set.seed(seed)

features <- read.csv(feature_path, check.names = FALSE, stringsAsFactors = FALSE)
manifest <- read.csv(manifest_path, check.names = FALSE, stringsAsFactors = FALSE)
feature_name_column <- intersect(c("feature", "metric", "feature_name", "Metrics"), names(manifest))[1]
if (is.na(feature_name_column)) stop("Feature manifest needs a feature, metric, feature_name, or Metrics column.")
feature_columns <- as.character(manifest[[feature_name_column]])
if (length(feature_columns) != 36L || length(unique(feature_columns)) != 36L) {
  stop("Exactly 36 unique morphology features are required.")
}

weight_columns <- paste0("morphology_weight_M", 1:5)
if ("neuron_ID" %in% names(features) && !"source_cell_id" %in% names(features)) {
  names(features)[names(features) == "neuron_ID"] <- "source_cell_id"
}
if ("published_broad_class" %in% names(features) && !"broad_cell_type" %in% names(features)) {
  names(features)[names(features) == "published_broad_class"] <- "broad_cell_type"
}
for (index in 1:5) {
  old_name <- paste0("weight_M", index)
  new_name <- paste0("morphology_weight_M", index)
  if (old_name %in% names(features) && !new_name %in% names(features)) {
    names(features)[names(features) == old_name] <- new_name
  }
}
required_columns <- c(
  "source_cell_id", "dataset", "animal_id", "broad_cell_type",
  "RNA_subpopulation", weight_columns, feature_columns
)
if (!all(required_columns %in% names(features))) {
  stop("The published-neuron morphology-separation table is incomplete.")
}
cells <- features[, required_columns, drop = FALSE]
cells$source_cell_id <- as.character(cells$source_cell_id)
cells <- cells[cells$dataset == "Sorensen2026", , drop = FALSE]
if (any(!is.finite(as.matrix(cells[, feature_columns, drop = FALSE])))) {
  stop("Morphology features must be finite after preprocessing.")
}

pair_definition <- data.frame(
  pair_id = c("T0_vs_T1", "T2_vs_T3"),
  broad_layer = c("L4/5 IT", "L2/3 IT"),
  subpopulation_a = c("T0", "T2"),
  subpopulation_b = c("T1", "T3"),
  stringsAsFactors = FALSE
)
cells$pair_id <- ifelse(
  cells$broad_cell_type == "L4/5 IT", "T0_vs_T1",
  ifelse(cells$broad_cell_type == "L2/3 IT", "T2_vs_T3", NA_character_)
)
cells <- cells[!is.na(cells$pair_id), , drop = FALSE]

profile_columns <- c(feature_columns, weight_columns)
animal_subpopulation <- aggregate(
  cells[, profile_columns, drop = FALSE],
  by = cells[, c("animal_id", "broad_cell_type", "RNA_subpopulation", "pair_id"), drop = FALSE],
  FUN = mean
)
cell_counts <- aggregate(
  source_cell_id ~ animal_id + broad_cell_type + RNA_subpopulation + pair_id,
  data = cells, FUN = length
)
names(cell_counts)[names(cell_counts) == "source_cell_id"] <- "n_cells"
animal_subpopulation <- merge(
  animal_subpopulation, cell_counts,
  by = c("animal_id", "broad_cell_type", "RNA_subpopulation", "pair_id"),
  all.x = TRUE, sort = FALSE
)
write.csv(animal_subpopulation, file.path(output_dir, "animal_subpopulation_morphology_profiles.csv"), row.names = FALSE, na = "")

auc_rank <- function(score, positive) {
  positive <- as.integer(positive)
  n_positive <- sum(positive == 1L)
  n_negative <- sum(positive == 0L)
  if (n_positive == 0L || n_negative == 0L) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[positive == 1L]) - n_positive * (n_positive + 1) / 2) / (n_positive * n_negative)
}

rms_shift <- function(x, y) sqrt(mean((colMeans(x) - colMeans(y))^2))

js_distance <- function(p, q) {
  p <- as.numeric(p) / sum(p)
  q <- as.numeric(q) / sum(q)
  midpoint <- (p + q) / 2
  kl <- function(x, y) {
    keep <- x > 0
    sum(x[keep] * log2(x[keep] / y[keep]))
  }
  sqrt(0.5 * kl(p, midpoint) + 0.5 * kl(q, midpoint))
}

loao_predictions <- function(data, subpopulation_a, subpopulation_b) {
  data <- data[data$RNA_subpopulation %in% c(subpopulation_a, subpopulation_b), , drop = FALSE]
  output <- list()
  for (held_out_animal in unique(data$animal_id)) {
    train <- data[data$animal_id != held_out_animal, , drop = FALSE]
    test <- data[data$animal_id == held_out_animal, , drop = FALSE]
    if (sum(train$RNA_subpopulation == subpopulation_a) < 2L || sum(train$RNA_subpopulation == subpopulation_b) < 2L) next
    training_matrix <- as.matrix(train[, feature_columns, drop = FALSE])
    test_matrix <- as.matrix(test[, feature_columns, drop = FALSE])
    mean_a <- colMeans(training_matrix[train$RNA_subpopulation == subpopulation_a, , drop = FALSE])
    mean_b <- colMeans(training_matrix[train$RNA_subpopulation == subpopulation_b, , drop = FALSE])
    direction <- mean_a - mean_b
    midpoint <- (mean_a + mean_b) / 2
    score <- as.numeric(sweep(test_matrix, 2, midpoint, "-") %*% direction)
    output[[length(output) + 1L]] <- data.frame(
      animal_id = test$animal_id,
      broad_layer = test$broad_cell_type,
      pair_id = test$pair_id,
      RNA_subpopulation = test$RNA_subpopulation,
      score = score,
      heldout_positive = test$RNA_subpopulation == subpopulation_a,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, output)
}

prediction_rows <- list()
summary_rows <- list()
feature_effect_rows <- list()
for (index in seq_len(nrow(pair_definition))) {
  pair <- pair_definition[index, ]
  profile <- animal_subpopulation[animal_subpopulation$pair_id == pair$pair_id, , drop = FALSE]
  group_a <- as.matrix(profile[profile$RNA_subpopulation == pair$subpopulation_a, feature_columns, drop = FALSE])
  group_b <- as.matrix(profile[profile$RNA_subpopulation == pair$subpopulation_b, feature_columns, drop = FALSE])
  predictions <- loao_predictions(profile, pair$subpopulation_a, pair$subpopulation_b)
  predictions$pair_id <- pair$pair_id
  prediction_rows[[index]] <- predictions
  summary_rows[[index]] <- data.frame(
    pair_id = pair$pair_id, broad_layer = pair$broad_layer,
    subpopulation_a = pair$subpopulation_a, subpopulation_b = pair$subpopulation_b,
    RMS_shift = rms_shift(group_a, group_b),
    held_out_AUC = auc_rank(predictions$score, predictions$heldout_positive),
    n_subpopulation_a_profiles = nrow(group_a), n_subpopulation_b_profiles = nrow(group_b),
    n_subpopulation_a_cells = sum(cells$RNA_subpopulation == pair$subpopulation_a),
    n_subpopulation_b_cells = sum(cells$RNA_subpopulation == pair$subpopulation_b),
    stringsAsFactors = FALSE
  )
  mean_a <- colMeans(group_a)
  mean_b <- colMeans(group_b)
  feature_effect_rows[[index]] <- data.frame(
    pair_id = pair$pair_id, feature = feature_columns,
    subpopulation_a_mean = mean_a, subpopulation_b_mean = mean_b,
    subpopulation_a_minus_subpopulation_b = mean_a - mean_b,
    stringsAsFactors = FALSE
  )
}
predictions <- do.call(rbind, prediction_rows)
summary_table <- do.call(rbind, summary_rows)
feature_effects <- do.call(rbind, feature_effect_rows)

bootstrap_pair <- function(data, prediction_data, subpopulation_a, subpopulation_b, iterations) {
  animals <- unique(data$animal_id)
  data_animal_index <- match(data$animal_id, animals)
  prediction_animal_index <- match(prediction_data$animal_id, animals)
  if (anyNA(prediction_animal_index)) stop("Profiles and held-out predictions do not align.")
  matrix <- as.matrix(data[, feature_columns, drop = FALSE])
  subpopulation_a_flag <- data$RNA_subpopulation == subpopulation_a
  subpopulation_b_flag <- data$RNA_subpopulation == subpopulation_b
  prediction_a <- prediction_data$heldout_positive
  prediction_b <- !prediction_data$heldout_positive
  score_a <- prediction_data$score[prediction_a]
  score_b <- prediction_data$score[prediction_b]
  auc_comparison <- outer(
    score_a, score_b,
    FUN = function(a, b) as.numeric(a > b) + 0.5 * as.numeric(a == b)
  )
  output <- matrix(NA_real_, nrow = iterations, ncol = 2L)
  colnames(output) <- c("RMS", "AUC")
  for (iteration in seq_len(iterations)) {
    repeat {
      animal_count <- tabulate(
        sample.int(length(animals), length(animals), replace = TRUE),
        nbins = length(animals)
      )
      row_weight <- animal_count[data_animal_index]
      weight_a <- row_weight * subpopulation_a_flag
      weight_b <- row_weight * subpopulation_b_flag
      if (sum(weight_a) > 0 && sum(weight_b) > 0) break
    }
    mean_a <- colSums(matrix * weight_a) / sum(weight_a)
    mean_b <- colSums(matrix * weight_b) / sum(weight_b)
    output[iteration, "RMS"] <- sqrt(mean((mean_a - mean_b)^2))
    prediction_weight <- animal_count[prediction_animal_index]
    weight_prediction_a <- prediction_weight[prediction_a]
    weight_prediction_b <- prediction_weight[prediction_b]
    output[iteration, "AUC"] <- sum(
      auc_comparison * outer(weight_prediction_a, weight_prediction_b)
    ) / (sum(weight_prediction_a) * sum(weight_prediction_b))
  }
  output
}

set.seed(seed)
bootstrap <- data.frame(iteration = seq_len(n_bootstrap))
for (index in seq_len(nrow(pair_definition))) {
  pair <- pair_definition[index, ]
  profile <- animal_subpopulation[animal_subpopulation$pair_id == pair$pair_id, , drop = FALSE]
  pair_predictions <- predictions[predictions$pair_id == pair$pair_id, , drop = FALSE]
  draws <- bootstrap_pair(profile, pair_predictions, pair$subpopulation_a, pair$subpopulation_b, n_bootstrap)
  bootstrap[[paste0(pair$pair_id, "__RMS")]] <- draws[, "RMS"]
  bootstrap[[paste0(pair$pair_id, "__AUC")]] <- draws[, "AUC"]
}

interval <- function(values) unname(quantile(values, c(0.025, 0.975), na.rm = TRUE))
summary_table$RMS_lower_95 <- NA_real_
summary_table$RMS_upper_95 <- NA_real_
summary_table$AUC_lower_95 <- NA_real_
summary_table$AUC_upper_95 <- NA_real_
for (index in seq_len(nrow(summary_table))) {
  pair_id <- summary_table$pair_id[index]
  summary_table[index, c("RMS_lower_95", "RMS_upper_95")] <- interval(bootstrap[[paste0(pair_id, "__RMS")]])
  summary_table[index, c("AUC_lower_95", "AUC_upper_95")] <- interval(bootstrap[[paste0(pair_id, "__AUC")]])
}

enumerate_or_sample <- function(n, n_positive, iterations, max_enumeration) {
  total <- choose(n, n_positive)
  if (is.finite(total) && total <= max_enumeration) {
    combinations <- combn(n, n_positive)
    list(indices = combinations, mode = "exact", draws = ncol(combinations), total = total)
  } else {
    sampled <- replicate(iterations, sample.int(n, n_positive, replace = FALSE))
    list(indices = sampled, mode = "Monte Carlo", draws = ncol(sampled), total = total)
  }
}

set.seed(seed + 1L)
permutation_objects <- list()
permutation_rows <- list()
for (index in seq_len(nrow(pair_definition))) {
  pair <- pair_definition[index, ]
  profile <- animal_subpopulation[animal_subpopulation$pair_id == pair$pair_id, , drop = FALSE]
  animal_frequency <- table(profile$animal_id)
  shared_animals <- names(animal_frequency)[animal_frequency > 1L]
  unique_only <- profile[!profile$animal_id %in% shared_animals, , drop = FALSE]
  unique_only <- unique_only[order(unique_only$animal_id), , drop = FALSE]
  observed_positive <- unique_only$RNA_subpopulation == pair$subpopulation_a
  assignments <- enumerate_or_sample(
    nrow(unique_only), sum(observed_positive), n_randomization, max_exact
  )
  probability_matrix <- as.matrix(unique_only[, weight_columns, drop = FALSE])
  feature_matrix <- as.matrix(unique_only[, feature_columns, drop = FALSE])
  metrics_for_assignment <- function(positive) {
    c(
      JS = js_distance(colMeans(probability_matrix[positive, , drop = FALSE]),
                       colMeans(probability_matrix[!positive, , drop = FALSE])),
      RMS = rms_shift(feature_matrix[positive, , drop = FALSE],
                      feature_matrix[!positive, , drop = FALSE])
    )
  }
  observed <- metrics_for_assignment(observed_positive)
  null <- matrix(NA_real_, nrow = assignments$draws, ncol = 2L, dimnames = list(NULL, c("JS", "RMS")))
  for (draw in seq_len(assignments$draws)) {
    positive <- seq_len(nrow(unique_only)) %in% assignments$indices[, draw]
    null[draw, ] <- metrics_for_assignment(positive)
  }
  permutation_objects[[pair$pair_id]] <- list(observed = observed, null = null)
  for (metric in c("JS", "RMS")) {
    exceed <- sum(null[, metric] >= observed[metric] - 1e-12)
    p_value <- if (assignments$mode == "exact") {
      exceed / assignments$draws
    } else {
      (exceed + 1) / (assignments$draws + 1)
    }
    permutation_rows[[length(permutation_rows) + 1L]] <- data.frame(
      test_scope = "within_pair_subpopulation_association", pair_id = pair$pair_id,
      metric = metric, observed = observed[metric], p_value = p_value,
      mode = assignments$mode, draws = assignments$draws,
      n_profiles = nrow(unique_only), excluded_shared_subpopulation_animals = length(shared_animals),
      stringsAsFactors = FALSE
    )
  }
}

set.seed(seed + 2L)
for (metric in c("JS", "RMS")) {
  null_a <- permutation_objects[["T0_vs_T1"]]$null[, metric]
  null_b <- permutation_objects[["T2_vs_T3"]]$null[, metric]
  null_contrast <- sample(null_a, n_randomization, replace = TRUE) -
    sample(null_b, n_randomization, replace = TRUE)
  observed_contrast <- permutation_objects[["T0_vs_T1"]]$observed[metric] -
    permutation_objects[["T2_vs_T3"]]$observed[metric]
  p_value <- (sum(null_contrast >= observed_contrast - 1e-12) + 1) /
    (n_randomization + 1)
  permutation_rows[[length(permutation_rows) + 1L]] <- data.frame(
    test_scope = "T0T1_greater_than_T2T3", pair_id = "T0_vs_T1_minus_T2_vs_T3",
    metric = metric, observed = observed_contrast, p_value = p_value,
    mode = "Monte Carlo from independent within-layer assignment nulls",
    draws = n_randomization, n_profiles = NA_integer_,
    excluded_shared_subpopulation_animals = NA_integer_, stringsAsFactors = FALSE
  )
}
permutation_summary <- do.call(rbind, permutation_rows)
permutation_summary$p_BH_six_tests <- p.adjust(permutation_summary$p_value, method = "BH")
permutation_summary$manuscript_role <- ifelse(
  permutation_summary$test_scope == "T0T1_greater_than_T2T3" & permutation_summary$metric == "RMS",
  "reported post hoc contrast", "multiplicity-family audit only"
)

roc_rows <- list()
for (index in seq_len(nrow(pair_definition))) {
  pair <- pair_definition[index, ]
  prediction <- predictions[predictions$pair_id == pair$pair_id, , drop = FALSE]
  prediction <- prediction[order(prediction$score, decreasing = TRUE), , drop = FALSE]
  positive <- prediction$heldout_positive
  true_positive <- c(0, cumsum(positive))
  false_positive <- c(0, cumsum(!positive))
  roc_rows[[index]] <- data.frame(
    pair_id = pair$pair_id, threshold_rank = 0:nrow(prediction),
    false_positive_rate = false_positive / sum(!positive),
    true_positive_rate = true_positive / sum(positive),
    threshold = c(Inf, prediction$score),
    AUC = auc_rank(prediction$score, positive),
    stringsAsFactors = FALSE
  )
}

write.csv(summary_table, file.path(output_dir, "morphology_separation_summary.csv"), row.names = FALSE, na = "")
write.csv(feature_effects, file.path(output_dir, "feature_level_subpopulation_differences.csv"), row.names = FALSE, na = "")
write.csv(predictions, file.path(output_dir, "leave_one_animal_out_predictions.csv"), row.names = FALSE, na = "")
write.csv(do.call(rbind, roc_rows), file.path(output_dir, "roc_coordinates.csv"), row.names = FALSE, na = "")
write.csv(permutation_summary, file.path(output_dir, "animal_label_randomization.csv"), row.names = FALSE, na = "")
bootstrap_connection <- gzfile(file.path(output_dir, "whole_animal_bootstrap.csv.gz"), "wt")
write.csv(bootstrap, bootstrap_connection, row.names = FALSE, na = "")
close(bootstrap_connection)
write.csv(data.frame(
  parameter = c(
    "seed", "dataset", "feature_count", "bootstrap_draws",
    "bootstrap_interval", "randomization_draws", "maximum_exact_permutations",
    "subpopulation_pairs"
  ),
  value = c(
    seed, "Sorensen2026", length(feature_columns), n_bootstrap,
    "0.025;0.975", n_randomization, max_exact, "T0_vs_T1;T2_vs_T3"
  ),
  stringsAsFactors = FALSE
), file.path(output_dir, "run_parameters.csv"), row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
message("Morphology separation analysis complete: ", output_dir)
