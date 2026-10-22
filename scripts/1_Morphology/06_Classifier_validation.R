# Command-line arguments, in order:
# 1. Classifier-label CSV containing Neuron, Label, and predSVM.
# 2. Expert-label CSV containing Neuron and CellType.
# 3. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Expected 3 command-line arguments in the order documented above.")
}
classifier_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
expert_file <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
output_dir <- args[[3]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

old <- read.csv(
  classifier_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)
manual <- read.csv(
  expert_file,
  check.names = FALSE,
  stringsAsFactors = FALSE
)

if (!all(c("Neuron", "Label", "predSVM") %in% colnames(old))) {
  stop("old_labels.csv must contain Neuron, Label, and predSVM")
}
if (!all(c("Neuron", "CellType") %in% colnames(manual))) {
  stop("cellType_manual_annotation.csv must contain Neuron and CellType")
}
if (anyDuplicated(old$Neuron)) stop("old_labels.csv has duplicated Neuron identifiers")
if (anyDuplicated(manual$Neuron)) stop("Manual annotations have duplicated Neuron identifiers")

match_index <- match(manual$Neuron, old$Neuron)
if (anyNA(match_index)) {
  stop("Manual neurons absent from old_labels.csv: ", paste(manual$Neuron[is.na(match_index)], collapse = ", "))
}
joined <- cbind(
  manual,
  old[match_index, setdiff(colnames(old), c("", "Neuron")), drop = FALSE]
)

valid_classes <- c("L2_3_IT", "L4_5_IT")
if (!all(joined$CellType %in% valid_classes)) stop("Unexpected manual CellType value")
if (!all(joined$Label %in% valid_classes)) stop("Unexpected old Label value")
if (any(!is.finite(joined$predSVM))) stop("predSVM contains missing/non-finite values")
if (any(joined$predSVM < 0 | joined$predSVM > 1)) stop("predSVM is outside [0, 1]")

joined$ManualL4_5 <- as.integer(joined$CellType == "L4_5_IT")
joined$PredictedL4_5 <- as.integer(joined$Label == "L4_5_IT")
joined$ThresholdLabel <- ifelse(joined$predSVM > 0.5, "L4_5_IT", "L2_3_IT")
joined$Correct <- joined$CellType == joined$Label
joined$Confidence <- pmax(joined$predSVM, 1 - joined$predSVM)
joined$AbsoluteCalibrationError <- abs(joined$ManualL4_5 - joined$predSVM)
if (!all(joined$ThresholdLabel == joined$Label)) {
  stop("Stored Label is not fully consistent with predSVM > 0.5")
}

safe_divide <- function(numerator, denominator) {
  if (!is.finite(denominator) || denominator == 0) return(NA_real_)
  numerator / denominator
}

roc_auc <- function(y, score) {
  positive <- sum(y == 1)
  negative <- sum(y == 0)
  if (positive == 0 || negative == 0) return(NA_real_)
  ranks <- rank(score, ties.method = "average")
  (sum(ranks[y == 1]) - positive * (positive + 1) / 2) / (positive * negative)
}

average_precision <- function(y, score) {
  if (sum(y == 1) == 0) return(NA_real_)
  order_index <- order(score, decreasing = TRUE)
  y_ordered <- y[order_index]
  precision <- cumsum(y_ordered) / seq_along(y_ordered)
  mean(precision[y_ordered == 1])
}

calculate_metrics <- function(data) {
  y <- data$ManualL4_5
  pred <- data$PredictedL4_5
  score <- pmin(pmax(data$predSVM, 1e-8), 1 - 1e-8)
  tp <- sum(y == 1 & pred == 1)
  tn <- sum(y == 0 & pred == 0)
  fp <- sum(y == 0 & pred == 1)
  fn <- sum(y == 1 & pred == 0)
  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  precision <- safe_divide(tp, tp + fp)
  npv <- safe_divide(tn, tn + fn)
  accuracy <- safe_divide(tp + tn, length(y))
  balanced_accuracy <- mean(c(sensitivity, specificity), na.rm = TRUE)
  f1 <- safe_divide(2 * tp, 2 * tp + fp + fn)
  mcc_denominator <- sqrt(
    as.numeric(tp + fp) * as.numeric(tp + fn) *
      as.numeric(tn + fp) * as.numeric(tn + fn)
  )
  mcc <- safe_divide(tp * tn - fp * fn, mcc_denominator)
  observed <- accuracy
  expected <- safe_divide(
    (tp + fp) * (tp + fn) + (fn + tn) * (fp + tn),
    length(y)^2
  )
  kappa <- safe_divide(observed - expected, 1 - expected)

  logit_score <- qlogis(score)
  calibration_fit <- suppressWarnings(glm(y ~ logit_score, family = binomial()))
  calibration_intercept_fit <- suppressWarnings(glm(
    y ~ 1,
    offset = logit_score,
    family = binomial()
  ))

  c(
    N = length(y),
    Manual_L4_5_prevalence = mean(y),
    Accuracy = accuracy,
    Balanced_accuracy = balanced_accuracy,
    Sensitivity_L4_5 = sensitivity,
    Specificity_L2_3 = specificity,
    Precision_L4_5 = precision,
    Negative_predictive_value_L2_3 = npv,
    F1_L4_5 = f1,
    Matthews_correlation = mcc,
    Cohen_kappa = kappa,
    ROC_AUC = roc_auc(y, score),
    Average_precision_PR = average_precision(y, score),
    Brier_score = mean((score - y)^2),
    Log_loss = -mean(y * log(score) + (1 - y) * log(1 - score)),
    Calibration_intercept = unname(coef(calibration_intercept_fit)[1]),
    Calibration_slope = unname(coef(calibration_fit)[2])
  )
}

observed_metrics <- calculate_metrics(joined)

# Cluster bootstrap: sample brains rather than individual neurons.
set.seed(20260716)
bootstrap_replicates <- 2000L
brains <- unique(joined$Brain)
bootstrap_matrix <- matrix(
  NA_real_,
  nrow = bootstrap_replicates,
  ncol = length(observed_metrics),
  dimnames = list(NULL, names(observed_metrics))
)
for (b in seq_len(bootstrap_replicates)) {
  sampled_brains <- sample(brains, length(brains), replace = TRUE)
  bootstrap_data <- do.call(rbind, lapply(sampled_brains, function(brain) {
    joined[joined$Brain == brain, , drop = FALSE]
  }))
  bootstrap_matrix[b, ] <- calculate_metrics(bootstrap_data)
}

ci_lower <- apply(bootstrap_matrix, 2, quantile, probs = 0.025, na.rm = TRUE)
ci_upper <- apply(bootstrap_matrix, 2, quantile, probs = 0.975, na.rm = TRUE)
ci_lower["N"] <- NA_real_
ci_upper["N"] <- NA_real_
summary_metrics <- data.frame(
  Metric = names(observed_metrics),
  Estimate = as.numeric(observed_metrics),
  Cluster_bootstrap_CI_lower = as.numeric(ci_lower[names(observed_metrics)]),
  Cluster_bootstrap_CI_upper = as.numeric(ci_upper[names(observed_metrics)]),
  Direction = ifelse(names(observed_metrics) %in% c("N", "Manual_L4_5_prevalence"), "Descriptive",
    ifelse(names(observed_metrics) %in% c("Brier_score", "Log_loss"), "Lower is better",
      ifelse(names(observed_metrics) %in% c("Calibration_intercept"), "Ideal = 0",
        ifelse(names(observed_metrics) %in% c("Calibration_slope"), "Ideal = 1", "Higher is better")
      )
    )
  ),
  stringsAsFactors = FALSE
)
write.csv(
  summary_metrics,
  file.path(output_dir, "validation_summary_metrics.csv"),
  row.names = FALSE
)

confusion <- as.data.frame.matrix(table(
  Manual = joined$CellType,
  Predicted = joined$Label
))
confusion <- cbind(Manual = rownames(confusion), confusion)
rownames(confusion) <- NULL
write.csv(confusion, file.path(output_dir, "confusion_matrix.csv"), row.names = FALSE)

thresholds <- sort(unique(c(seq(0, 1, by = 0.01), joined$predSVM)))
threshold_performance <- do.call(rbind, lapply(thresholds, function(threshold) {
  y <- joined$ManualL4_5
  pred <- as.integer(joined$predSVM > threshold)
  tp <- sum(y == 1 & pred == 1)
  tn <- sum(y == 0 & pred == 0)
  fp <- sum(y == 0 & pred == 1)
  fn <- sum(y == 1 & pred == 0)
  sensitivity <- safe_divide(tp, tp + fn)
  specificity <- safe_divide(tn, tn + fp)
  precision <- safe_divide(tp, tp + fp)
  data.frame(
    Threshold = threshold,
    Accuracy = safe_divide(tp + tn, length(y)),
    Balanced_accuracy = mean(c(sensitivity, specificity), na.rm = TRUE),
    Sensitivity_L4_5 = sensitivity,
    Specificity_L2_3 = specificity,
    Precision_L4_5 = precision,
    F1_L4_5 = safe_divide(2 * tp, 2 * tp + fp + fn),
    Youden_J = sensitivity + specificity - 1,
    Predicted_L4_5 = tp + fp,
    stringsAsFactors = FALSE
  )
}))
write.csv(
  threshold_performance,
  file.path(output_dir, "threshold_performance.csv"),
  row.names = FALSE
)

confidence_thresholds <- seq(0.5, 0.95, by = 0.05)
confidence_coverage <- do.call(rbind, lapply(confidence_thresholds, function(threshold) {
  retained <- joined$Confidence >= threshold
  subset <- joined[retained, , drop = FALSE]
  metrics <- if (nrow(subset) > 0) calculate_metrics(subset) else rep(NA_real_, length(observed_metrics))
  data.frame(
    Confidence_threshold = threshold,
    Retained_cells = sum(retained),
    Coverage = mean(retained),
    Accuracy = unname(metrics["Accuracy"]),
    Balanced_accuracy = unname(metrics["Balanced_accuracy"]),
    Sensitivity_L4_5 = unname(metrics["Sensitivity_L4_5"]),
    Specificity_L2_3 = unname(metrics["Specificity_L2_3"]),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  confidence_coverage,
  file.path(output_dir, "confidence_coverage.csv"),
  row.names = FALSE
)

calibration_breaks <- seq(0, 1, by = 0.1)
joined$Calibration_bin <- cut(
  joined$predSVM,
  breaks = calibration_breaks,
  include.lowest = TRUE,
  right = TRUE
)
calibration_bins <- do.call(rbind, lapply(levels(joined$Calibration_bin), function(bin) {
  subset <- joined[joined$Calibration_bin == bin, , drop = FALSE]
  if (nrow(subset) == 0) return(NULL)
  data.frame(
    Score_bin = bin,
    N = nrow(subset),
    Mean_predicted_L4_5_score = mean(subset$predSVM),
    Observed_manual_L4_5_fraction = mean(subset$ManualL4_5),
    Bin_accuracy = mean(subset$Correct),
    stringsAsFactors = FALSE
  )
}))
write.csv(
  calibration_bins,
  file.path(output_dir, "calibration_bins.csv"),
  row.names = FALSE
)

group_summary <- function(group_column) {
  groups <- unique(joined[[group_column]])
  do.call(rbind, lapply(groups, function(group) {
    subset <- joined[joined[[group_column]] == group, , drop = FALSE]
    data.frame(
      Stratum = group_column,
      Group = as.character(group),
      N = nrow(subset),
      Brains = length(unique(subset$Brain)),
      Manual_L4_5 = sum(subset$ManualL4_5),
      Predicted_L4_5 = sum(subset$PredictedL4_5),
      Accuracy = mean(subset$Correct),
      stringsAsFactors = FALSE
    )
  }))
}
stratified_summary <- do.call(rbind, lapply(
  c("Brain", "Perturbation", "Type", "Age"),
  group_summary
))
write.csv(
  stratified_summary,
  file.path(output_dir, "stratified_agreement.csv"),
  row.names = FALSE
)

joined_export <- joined[, c(
  "Neuron", "Brain", "Perturbation", "Age", "Type", "CellType", "Label",
  "predSVM", "Confidence", "Correct", "AbsoluteCalibrationError"
)]
colnames(joined_export)[colnames(joined_export) == "CellType"] <- "ManualLabel"
colnames(joined_export)[colnames(joined_export) == "Label"] <- "OldSVMLabel"
colnames(joined_export)[colnames(joined_export) == "predSVM"] <- "OldSVM_L4_5_score"
write.csv(
  joined_export,
  file.path(output_dir, "manual_validation_joined.csv"),
  row.names = FALSE
)

best_threshold_row <- threshold_performance[
  which.max(threshold_performance$Youden_J),
  ,
  drop = FALSE
]
run_diagnostics <- data.frame(
  Check = c(
    "Old_label_rows",
    "Manual_annotation_rows",
    "Manual_unique_neurons",
    "Manual_neurons_matched",
    "Manual_brains",
    "Manual_L2_3_IT",
    "Manual_L4_5_IT",
    "Stored_labels_matching_score_threshold",
    "Score_minimum",
    "Score_maximum",
    "Exploratory_Youden_threshold",
    "Exploratory_Youden_balanced_accuracy",
    "Brain_cluster_bootstrap_replicates"
  ),
  Value = c(
    nrow(old),
    nrow(manual),
    length(unique(manual$Neuron)),
    nrow(joined),
    length(unique(joined$Brain)),
    sum(joined$CellType == "L2_3_IT"),
    sum(joined$CellType == "L4_5_IT"),
    sum(joined$ThresholdLabel == joined$Label),
    min(joined$predSVM),
    max(joined$predSVM),
    best_threshold_row$Threshold,
    best_threshold_row$Balanced_accuracy,
    bootstrap_replicates
  ),
  stringsAsFactors = FALSE
)
write.csv(
  run_diagnostics,
  file.path(output_dir, "run_diagnostics.csv"),
  row.names = FALSE
)

