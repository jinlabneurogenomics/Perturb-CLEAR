# Build the fixed four-subpopulation RNA reference from ST2 Perturb-seq neurons and
# map other Perturb-seq cells plus published paired-modality neurons by RPCA.

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(qs)
  library(Seurat)
})

# Command-line arguments, in order:
# 1. Perturb-seq object in qs format.
# 2. Published query-object list in qs format.
# 3. Shared-gene manifest CSV.
# 4. Output directory.
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Expected 4 command-line arguments in the order documented above.")
}

perturbseq_path <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
published_path <- normalizePath(args[[2]], winslash = "/", mustWork = TRUE)
gene_manifest_path <- normalizePath(args[[3]], winslash = "/", mustWork = TRUE)
output_dir <- normalizePath(args[[4]], winslash = "/", mustWork = FALSE)
table_dir <- file.path(output_dir, "tables")
data_dir <- file.path(output_dir, "data")
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)

seed <- 260909L
dims <- 30L
resolution_candidates <- c(0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
min_cluster_cells <- 30L
min_animals <- 3L
max_animal_fraction <- 0.70
min_clusters <- 3L
max_clusters <- 6L
set.seed(seed)
options(future.globals.maxSize = 100 * 1024^3)
future::plan(future::sequential)

gene_manifest <- fread(gene_manifest_path)
gene_column <- intersect(c("gene_name", "gene", "Gene"), names(gene_manifest))[1]
if (is.na(gene_column)) stop("The shared-gene manifest needs a gene_name, gene, or Gene column.")
shared_genes <- unique(as.character(gene_manifest[[gene_column]]))
shared_genes <- shared_genes[nzchar(shared_genes)]

make_object <- function(counts, metadata) {
  metadata <- metadata[match(colnames(counts), metadata$source_cell_id), , drop = FALSE]
  rownames(metadata) <- metadata$source_cell_id
  object <- CreateSeuratObject(
    counts = counts, assay = "RNA", meta.data = metadata,
    min.cells = 0, min.features = 0
  )
  NormalizeData(
    object, assay = "RNA", normalization.method = "LogNormalize",
    scale.factor = 10000, verbose = FALSE
  )
}

message("Loading Perturb-seq cells")
perturbseq <- qread(perturbseq_path)
perturbseq_counts <- LayerData(perturbseq[["RNA"]], layer = "counts")
shared_genes <- intersect(shared_genes, rownames(perturbseq_counts))
if (length(shared_genes) <= 10000L) stop("Unexpectedly small shared-gene universe.")
perturbseq_counts <- perturbseq_counts[shared_genes, , drop = FALSE]
perturbseq_meta <- as.data.frame(perturbseq@meta.data)
perturbseq_meta$source_cell_id <- rownames(perturbseq_meta)
required_metadata <- c("Assign", "orig.ident", "CT")
if (!all(required_metadata %in% names(perturbseq_meta))) {
  stop("Perturb-seq metadata needs Assign, orig.ident, and CT.")
}
perturbseq_meta$Assign <- as.character(perturbseq_meta$Assign)
perturbseq_meta$orig.ident <- as.character(perturbseq_meta$orig.ident)
perturbseq_meta$broad_cell_type <- as.character(perturbseq_meta$CT)
rm(perturbseq)
invisible(gc())

reference_meta <- perturbseq_meta[perturbseq_meta$Assign == "ST_2", , drop = FALSE]
perturbseq_query_meta <- perturbseq_meta[perturbseq_meta$Assign != "ST_2", , drop = FALSE]
reference <- make_object(
  perturbseq_counts[, reference_meta$source_cell_id, drop = FALSE],
  reference_meta
)

message("Fitting the ST2-only shared-gene reference")
reference <- ScaleData(
  reference, assay = "RNA", features = shared_genes,
  do.center = TRUE, do.scale = TRUE, block.size = 1000,
  min.cells.to.block = 0, verbose = FALSE
)
reference <- RunPCA(
  reference, assay = "RNA", features = shared_genes, npcs = dims,
  reduction.name = "pca.st2", reduction.key = "ST2PC_",
  approx = TRUE, seed.use = seed, verbose = FALSE
)
reference <- FindNeighbors(
  reference, reduction = "pca.st2", dims = seq_len(dims),
  k.param = 30, nn.method = "annoy", n.trees = 50,
  graph.name = c("ST2_nn", "ST2_snn"), verbose = FALSE
)

resolution_rows <- list()
for (resolution in resolution_candidates) {
  column <- paste0("ST2_res_", gsub("\\.", "_", format(resolution, nsmall = 1)))
  reference <- FindClusters(
    reference, graph.name = "ST2_snn", resolution = resolution,
    algorithm = 1, random.seed = seed, cluster.name = column,
    verbose = FALSE
  )
  balance <- data.frame(
    cluster = as.character(reference@meta.data[[column]]),
    animal = as.character(reference@meta.data$orig.ident),
    stringsAsFactors = FALSE
  )
  cluster_size <- table(balance$cluster)
  animals_per_cluster <- tapply(balance$animal, balance$cluster, function(x) length(unique(x)))
  animal_dominance <- tapply(balance$animal, balance$cluster, function(x) max(table(x)) / length(x))
  n_clusters <- length(cluster_size)
  passes <- min(cluster_size) >= min_cluster_cells &&
    min(animals_per_cluster) >= min_animals &&
    max(animal_dominance) <= max_animal_fraction &&
    n_clusters >= min_clusters && n_clusters <= max_clusters
  resolution_rows[[as.character(resolution)]] <- data.frame(
    resolution = resolution, cluster_column = column, n_clusters = n_clusters,
    min_cluster_cells = min(cluster_size),
    min_animals_per_cluster = min(animals_per_cluster),
    max_single_animal_fraction = max(animal_dominance),
    pass_prespecified_rule = passes,
    stringsAsFactors = FALSE
  )
}

resolution_table <- rbindlist(resolution_rows)
fwrite(resolution_table, file.path(table_dir, "reference_resolution_sweep.csv"))
eligible <- resolution_table[pass_prespecified_rule == TRUE]
if (!nrow(eligible)) stop("No clustering resolution passed the prespecified balance rule.")
selected_column <- eligible[which.max(resolution)]$cluster_column
raw_cluster <- as.character(reference@meta.data[[selected_column]])
reference$RNA_subpopulation <- factor(
  paste0("T", raw_cluster),
  levels = paste0("T", sort(unique(as.integer(raw_cluster))))
)
if (nlevels(reference$RNA_subpopulation) != 4L) {
  stop("The selected reference is not a four-subpopulation solution.")
}
reference$broad_cell_type <- as.character(reference$CT)
reference <- RunUMAP(
  reference, reduction = "pca.st2", dims = seq_len(dims),
  reduction.name = "umap.st2", reduction.key = "ST2UMAP_",
  n.neighbors = 30, min.dist = 0.3, metric = "cosine",
  return.model = TRUE, seed.use = seed, verbose = FALSE
)

selected_balance <- as.data.table(reference@meta.data)[, .N, by = .(RNA_subpopulation, orig.ident, broad_cell_type)]
selected_balance[, cluster_total := sum(N), by = RNA_subpopulation]
selected_balance[, fraction_of_cluster := N / cluster_total]
fwrite(selected_balance, file.path(table_dir, "reference_subpopulation_balance.csv"))

map_rpca <- function(reference_object, query_object, genes, dataset_name) {
  effective_dims <- min(dims, ncol(query_object) - 2L)
  dims_use <- seq_len(effective_dims)
  rpca_reduction <- paste0("pca.rpca.", effective_dims)
  reference_object[[rpca_reduction]] <- CreateDimReducObject(
    embeddings = Embeddings(reference_object, "pca.st2")[, dims_use, drop = FALSE],
    loadings = Loadings(reference_object, "pca.st2")[, dims_use, drop = FALSE],
    stdev = Stdev(reference_object[["pca.st2"]])[dims_use],
    key = "rpcaRefPC_", assay = "RNA"
  )
  anchors <- FindTransferAnchors(
    reference = reference_object, query = query_object,
    normalization.method = "LogNormalize",
    reference.assay = "RNA", query.assay = "RNA", reduction = "rpca",
    reference.reduction = rpca_reduction, features = genes, scale = TRUE,
    npcs = effective_dims, dims = dims_use, k.anchor = 5, k.filter = NA,
    k.score = 30, max.features = length(genes), nn.method = "annoy",
    n.trees = 50, approx.pca = TRUE, verbose = FALSE
  )
  reference_pca <- t(Embeddings(reference_object, "pca.st2")[, dims_use, drop = FALSE])
  rownames(reference_pca) <- paste0("PC_", dims_use)
  mapped <- MapQuery(
    anchorset = anchors, query = query_object, reference = reference_object,
    refdata = list(
      RNA_subpopulation = "RNA_subpopulation",
      broad_cell_type = "broad_cell_type",
      reference_PCA = reference_pca
    ),
    new.reduction.name = "ref.pca", reference.reduction = rpca_reduction,
    reference.dims = dims_use, store.weights = TRUE,
    transferdata.args = list(k.weight = 20),
    integrateembeddings.args = list(k.weight = 20),
    reduction.model = NULL, verbose = FALSE
  )
  mapped <- RunUMAP(
    mapped, reduction = "ref.pca", dims = dims_use,
    reduction.model = reference_object[["umap.st2"]],
    reduction.name = "ref.umap", reduction.key = "refUMAP_",
    verbose = FALSE
  )
  anchor_table <- as.data.frame(slot(anchors, "anchors"))
  list(
    mapped = mapped,
    anchor_summary = data.frame(
      dataset = dataset_name, query_cells = ncol(query_object),
      reference_cells = ncol(reference_object), shared_genes = length(genes),
      dimensions = effective_dims, anchors = nrow(anchor_table),
      mean_anchor_score = mean(anchor_table$score),
      median_anchor_score = median(anchor_table$score),
      stringsAsFactors = FALSE
    )
  )
}

extract_query <- function(mapped, metadata, dataset_name) {
  md <- as.data.frame(mapped@meta.data)
  md$source_cell_id <- rownames(md)
  pcs <- as.data.frame(Embeddings(mapped, "ref.pca"))
  pcs$source_cell_id <- rownames(pcs)
  colnames(pcs)[seq_len(ncol(pcs) - 1L)] <- paste0("PC_", seq_len(ncol(pcs) - 1L))
  umap <- as.data.frame(Embeddings(mapped, "ref.umap"))
  umap$source_cell_id <- rownames(umap)
  colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")
  result <- data.frame(
    source_cell_id = md$source_cell_id,
    dataset = dataset_name,
    mapping_role = "query",
    Assign = metadata$Assign[match(md$source_cell_id, metadata$source_cell_id)],
    animal_id = metadata$orig.ident[match(md$source_cell_id, metadata$source_cell_id)],
    source_broad_cell_type = metadata$broad_cell_type[match(md$source_cell_id, metadata$source_cell_id)],
    predicted_broad_cell_type = as.character(md$predicted.broad_cell_type),
    RNA_subpopulation = as.character(md$predicted.RNA_subpopulation),
    RNA_subpopulation_top_score = as.numeric(md$predicted.RNA_subpopulation.score),
    stringsAsFactors = FALSE
  )
  result <- merge(result, pcs, by = "source_cell_id", sort = FALSE)
  merge(result, umap, by = "source_cell_id", sort = FALSE)
}

query_objects <- list()
query_metadata <- list()
query_objects[["Perturbseq_query"]] <- make_object(
  perturbseq_counts[, perturbseq_query_meta$source_cell_id, drop = FALSE],
  perturbseq_query_meta
)
query_metadata[["Perturbseq_query"]] <- perturbseq_query_meta

message("Loading published raw-count query objects")
published_saved <- qread(published_path)
for (dataset in c("Sorensen2026", "Scala2021")) {
  key_candidates <- c(dataset, paste0(dataset, "__Seurat_RPCA_anchors"))
  key <- key_candidates[key_candidates %in% names(published_saved)][1]
  if (is.na(key)) next
  old_object <- published_saved[[key]]
  published_counts <- LayerData(old_object[["RNA"]], layer = "counts")
  published_counts <- published_counts[shared_genes, , drop = FALSE]
  published_meta <- as.data.frame(old_object@meta.data)
  published_meta$source_cell_id <- rownames(published_meta)
  published_meta$Assign <- if ("Assign" %in% names(published_meta)) as.character(published_meta$Assign) else "WT"
  published_meta$orig.ident <- if ("animal_id" %in% names(published_meta)) as.character(published_meta$animal_id) else dataset
  published_meta$broad_cell_type <- if ("broad_cell_type" %in% names(published_meta)) as.character(published_meta$broad_cell_type) else NA_character_
  query_objects[[dataset]] <- make_object(published_counts, published_meta)
  query_metadata[[dataset]] <- published_meta
}
if (!"Sorensen2026" %in% names(query_objects)) stop("Sorensen2026 query object is required.")
rm(published_saved, perturbseq_counts)
invisible(gc())

query_frames <- list()
anchor_rows <- list()
mapped_objects <- list()
for (dataset in names(query_objects)) {
  message("RPCA mapping: ", dataset)
  mapped_result <- map_rpca(reference, query_objects[[dataset]], shared_genes, dataset)
  mapped_objects[[dataset]] <- mapped_result$mapped
  query_frames[[dataset]] <- extract_query(mapped_result$mapped, query_metadata[[dataset]], dataset)
  anchor_rows[[dataset]] <- mapped_result$anchor_summary
}

reference_pc <- as.data.frame(Embeddings(reference, "pca.st2"))
reference_pc$source_cell_id <- rownames(reference_pc)
colnames(reference_pc)[seq_len(ncol(reference_pc) - 1L)] <- paste0("PC_", seq_len(ncol(reference_pc) - 1L))
reference_umap <- as.data.frame(Embeddings(reference, "umap.st2"))
reference_umap$source_cell_id <- rownames(reference_umap)
colnames(reference_umap)[1:2] <- c("UMAP_1", "UMAP_2")
reference_metadata <- as.data.frame(reference@meta.data)
reference_metadata$source_cell_id <- rownames(reference_metadata)
reference_frame <- data.frame(
  source_cell_id = reference_metadata$source_cell_id,
  dataset = "Perturbseq_ST2_reference",
  mapping_role = "reference",
  Assign = "ST_2",
  animal_id = as.character(reference_metadata$orig.ident),
  source_broad_cell_type = as.character(reference_metadata$broad_cell_type),
  predicted_broad_cell_type = as.character(reference_metadata$broad_cell_type),
  RNA_subpopulation = as.character(reference_metadata$RNA_subpopulation),
  RNA_subpopulation_top_score = 1,
  stringsAsFactors = FALSE
)
reference_frame <- merge(reference_frame, reference_pc, by = "source_cell_id", sort = FALSE)
reference_frame <- merge(reference_frame, reference_umap, by = "source_cell_id", sort = FALSE)
all_mapping <- rbindlist(c(list(reference_frame), query_frames), fill = TRUE)
fwrite(all_mapping, file.path(table_dir, "all_cell_reference_mapping.csv.gz"), na = "NA")
fwrite(rbindlist(anchor_rows), file.path(table_dir, "anchor_summary.csv"), na = "NA")

score_rows <- list()
for (dataset in names(mapped_objects)) {
  md <- as.data.frame(mapped_objects[[dataset]]@meta.data)
  md$source_cell_id <- rownames(md)
  score_columns <- grep("^prediction.score.RNA_subpopulation\\.", names(md), value = TRUE)
  if (!length(score_columns)) next
  score_matrix <- as.matrix(md[, score_columns, drop = FALSE])
  colnames(score_matrix) <- sub("^prediction.score.RNA_subpopulation\\.", "", score_columns)
  score_matrix <- score_matrix / pmax(rowSums(score_matrix), .Machine$double.eps)
  entropy <- -rowSums(ifelse(score_matrix > 0, score_matrix * log(score_matrix), 0)) / log(ncol(score_matrix))
  long <- as.data.table(score_matrix, keep.rownames = "source_cell_id")
  long <- melt(long, id.vars = "source_cell_id", variable.name = "RNA_subpopulation", value.name = "probability")
  long[, `:=`(
    dataset = dataset,
    normalized_entropy = entropy[match(source_cell_id, rownames(score_matrix))]
  )]
  score_rows[[dataset]] <- long
}
fwrite(rbindlist(score_rows), file.path(table_dir, "four_subpopulation_mapping_scores.csv.gz"), na = "NA")

qsave(reference, file.path(data_dir, "four_subpopulation_reference.qs"), preset = "balanced")
qsave(mapped_objects, file.path(data_dir, "mapped_queries.qs"), preset = "high")

fwrite(data.table(
  parameter = c(
    "seed", "reference_assignment", "shared_gene_count", "PCA_dimensions",
    "neighbor_k", "resolution_candidates", "selected_resolution",
    "minimum_cluster_cells", "minimum_animals_per_cluster",
    "maximum_single_animal_fraction", "minimum_clusters", "maximum_clusters",
    "anchor_k", "score_k", "weight_k", "annoy_trees"
  ),
  value = c(
    seed, "ST_2", length(shared_genes), dims, 30L,
    paste(resolution_candidates, collapse = ";"),
    eligible[which.max(resolution)]$resolution,
    min_cluster_cells, min_animals, max_animal_fraction, min_clusters, max_clusters,
    5L, 30L, 20L, 50L
  )
), file.path(table_dir, "run_parameters.csv"))
writeLines(capture.output(sessionInfo()), file.path(table_dir, "R_sessionInfo.txt"))
message("Four-subpopulation RNA mapping complete: ", output_dir)
