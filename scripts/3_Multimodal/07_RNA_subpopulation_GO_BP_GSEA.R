# Build the Gene Ontology Biological Process mapping and run directional GSEA
# for the T0 versus T1 RNA-subpopulation marker model.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(fgsea)
  library(DBI)
  library(RSQLite)
})

# Command-line arguments, in order:
# 1. RNA-subpopulation marker-result CSV, optionally gzip-compressed.
# 2. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Expected 2 command-line arguments in the order documented above.")
}
marker_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[2]], winslash = "/", mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 20260918L
min_size <- 10L
max_size <- 500L
n_perm_simple <- 10000L
effect_boundary <- 9.9

marker <- fread(marker_file)
accepted_contrast_ids <- "T0_vs_T1"
selected_contrast_id <- "not supplied"
if ("contrast_id" %in% names(marker)) {
  available_contrasts <- unique(as.character(marker$contrast_id))
  selected_contrasts <- intersect(accepted_contrast_ids, available_contrasts)
  if (length(selected_contrasts) != 1L) {
    stop(
      "Expected exactly one T0-versus-T1 contrast identifier. Available values: ",
      paste(available_contrasts, collapse = ", ")
    )
  }
  selected_contrast_id <- selected_contrasts[[1]]
  marker <- marker[contrast_id == selected_contrast_id]
}
required_columns <- c("Gene", "lfc", "f_statistic", "pval", "df1")
if (!all(required_columns %in% names(marker))) {
  stop("Marker table is missing required columns: ", paste(setdiff(required_columns, names(marker)), collapse = ", "))
}
if (!nrow(marker) || anyDuplicated(marker$Gene)) {
  stop("The T0-versus-T1 marker table must contain one row per gene.")
}

gene_symbols <- unique(as.character(marker$Gene))
gene_symbols <- gene_symbols[nzchar(gene_symbols)]
org_package <- find.package("org.Mm.eg.db", quiet = TRUE)
go_package <- find.package("GO.db", quiet = TRUE)
if (!nzchar(org_package) || !nzchar(go_package)) {
  stop("Packages org.Mm.eg.db and GO.db are required.")
}
org_database <- file.path(org_package, "extdata", "org.Mm.eg.sqlite")
go_database <- file.path(go_package, "extdata", "GO.sqlite")

org_connection <- dbConnect(SQLite(), org_database)
on.exit(if (dbIsValid(org_connection)) dbDisconnect(org_connection), add = TRUE)
dbWriteTable(
  org_connection, "selected_symbols", data.frame(symbol = gene_symbols),
  temporary = TRUE, overwrite = TRUE
)
mapping <- as.data.table(dbGetQuery(
  org_connection,
  paste(
    "SELECT DISTINCT bp.go_id AS GO_ID, gi.symbol AS Gene",
    "FROM go_bp_all AS bp",
    "INNER JOIN gene_info AS gi ON bp._id = gi._id",
    "INNER JOIN selected_symbols AS selected ON gi.symbol = selected.symbol",
    "WHERE bp.go_id IS NOT NULL AND bp.go_id != ''",
    "ORDER BY bp.go_id, gi.symbol"
  )
))
dbDisconnect(org_connection)

ontology_connection <- dbConnect(SQLite(), go_database)
on.exit(if (dbIsValid(ontology_connection)) dbDisconnect(ontology_connection), add = TRUE)
go_terms <- as.data.table(dbGetQuery(
  ontology_connection,
  paste(
    "SELECT go_id AS GO_ID, term AS Description",
    "FROM go_term WHERE ontology = 'BP' ORDER BY go_id"
  )
))
dbDisconnect(ontology_connection)
mapping <- merge(mapping, go_terms, by = "GO_ID", all = FALSE, sort = FALSE)
mapping <- unique(mapping[, .(GO_ID, Gene, Description)])
setorder(mapping, GO_ID, Gene)

stable <- marker[
  is.finite(lfc) & is.finite(f_statistic) & is.finite(pval) &
    df1 == 1 & f_statistic > 0 & abs(lfc) < effect_boundary
]
stable[, signed_sqrt_F := sign(lfc) * sqrt(f_statistic)]
rank_values <- stable$signed_sqrt_F
names(rank_values) <- stable$Gene
rank_values <- sort(rank_values, decreasing = TRUE)
if (anyDuplicated(names(rank_values))) stop("Duplicated genes in the ranked vector.")

pathways <- lapply(split(mapping$Gene, mapping$GO_ID), unique)
go_terms <- unique(mapping[, .(GO_ID, Description)])
set.seed(seed)
result <- suppressWarnings(fgseaMultilevel(
  pathways = pathways,
  stats = rank_values,
  minSize = min_size,
  maxSize = max_size,
  eps = 0,
  nPermSimple = n_perm_simple,
  scoreType = "std",
  nproc = 1
))
result <- as.data.table(result)
result[, leadingEdge := NULL]
setnames(result, "pathway", "GO_ID")
result <- merge(result, go_terms, by = "GO_ID", all.x = TRUE, sort = FALSE)
result[, `:=`(
  direction = fifelse(NES > 0, "T0", "T1"),
  n_ranked_genes = length(rank_values),
  rank_method = "sign(log2FC) * sqrt(F)",
  seed = seed,
  minSize = min_size,
  maxSize = max_size
)]
result[, abs_NES_sort := abs(NES)]
setorderv(result, c("padj", "abs_NES_sort", "GO_ID"), c(1L, -1L, 1L))
result[, abs_NES_sort := NULL]

fwrite(mapping, file.path(output_dir, "GO_BP_term_to_gene.csv.gz"), compress = "gzip")
fwrite(data.table(
  software = c("R", "RSQLite", "org.Mm.eg.db", "GO.db", "fgsea"),
  version = c(
    as.character(getRversion()), as.character(packageVersion("RSQLite")),
    as.character(packageVersion("org.Mm.eg.db")), as.character(packageVersion("GO.db")),
    as.character(packageVersion("fgsea"))
  ),
  n_input_genes = length(gene_symbols),
  n_mapped_genes = uniqueN(mapping$Gene),
  n_BP_terms = uniqueN(mapping$GO_ID)
), file.path(output_dir, "GO_BP_mapping_manifest.csv"))
fwrite(result, file.path(output_dir, "go_bp_gsea_all.csv.gz"), compress = "gzip", na = "")
fwrite(data.table(
  n_tested_genes = nrow(marker),
  n_ranked_genes = length(rank_values),
  n_excluded_genes = nrow(marker) - length(rank_values)
), file.path(output_dir, "go_bp_gsea_rank_audit.csv"), na = "")
fwrite(data.table(
  parameter = c(
    "seed", "minSize", "maxSize", "nPermSimple", "eps", "scoreType",
    "nproc", "absolute_log2FC_boundary", "input_contrast_id"
  ),
  value = c(seed, min_size, max_size, n_perm_simple, 0, "std", 1L, effect_boundary, selected_contrast_id)
), file.path(output_dir, "run_parameters.csv"), na = "")
writeLines(capture.output(sessionInfo()), file.path(output_dir, "R_sessionInfo.txt"))
