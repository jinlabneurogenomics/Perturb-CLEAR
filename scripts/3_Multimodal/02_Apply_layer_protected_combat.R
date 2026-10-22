# Apply the manuscript morphology normalization only: parametric ComBat with
# Perturb-CLEAR controls as the reference batch and broad layer protected.

suppressPackageStartupMessages({
  library(sva)
})

# Command-line arguments, in order:
# 1. Standardized reference and query morphometrics CSV.
# 2. Output CSV.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Expected 2 command-line arguments in the order documented above.")
}

input_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_path <- args[[2]]
data <- read.csv(input_path, check.names = FALSE, stringsAsFactors = FALSE)
metadata_columns <- c("cell_id", "dataset", "broad_cell_type", "animal_id", "official_M")
if (!all(metadata_columns %in% names(data))) stop("Required metadata columns are missing.")
feature_columns <- setdiff(names(data), metadata_columns)
matrix <- as.matrix(data[, feature_columns, drop = FALSE])
storage.mode(matrix) <- "double"
if (any(!is.finite(matrix))) stop("The standardized feature matrix must be finite.")

batch <- factor(data$dataset)
if (!"Perturb_CLEAR" %in% levels(batch)) stop("Perturb_CLEAR reference batch is absent.")
layer_model <- model.matrix(~ broad_cell_type, data = data)
corrected <- ComBat(
  dat = t(matrix), batch = batch, mod = layer_model,
  par.prior = TRUE, prior.plots = FALSE, mean.only = FALSE,
  ref.batch = "Perturb_CLEAR"
)

output <- cbind(data[, metadata_columns, drop = FALSE], as.data.frame(t(corrected), check.names = FALSE))
write.csv(output, output_path, row.names = FALSE, na = "")
writeLines(capture.output(sessionInfo()), paste0(output_path, ".sessionInfo.txt"))
