source("src/pbmc_10k/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1] %||% "config/pbmc_10k_multi.yaml"
cfg <- read_config(config_path)
ensure_dirs(cfg)

object <- readRDS(file.path(cfg$output_dir, "objects", "01_filtered_seurat.rds"))

message("Running standard scRNA workflow")
object <- NormalizeData(object, verbose = FALSE)
object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = cfg$clustering$nfeatures, verbose = FALSE)
vars_to_regress <- if (isTRUE(cfg$clustering$regress_percent_mt)) "percent.mt" else NULL
object <- ScaleData(object, vars.to.regress = vars_to_regress, verbose = FALSE)
object <- RunPCA(object, npcs = cfg$clustering$npcs, verbose = FALSE)
object <- FindNeighbors(object, dims = seq_len(cfg$clustering$npcs), verbose = FALSE)
object <- FindClusters(object, resolution = cfg$clustering$resolution, verbose = FALSE)
object <- RunUMAP(object, dims = seq_len(cfg$clustering$npcs), verbose = FALSE)

cluster_counts <- as.data.table(table(object$seurat_clusters))
setnames(cluster_counts, c("cluster", "n_cells"))
fwrite(cluster_counts, file.path(cfg$output_dir, "02_clustering", "cluster_counts.csv"))

cluster_umap <- DimPlot(object, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("PBMC 10k clusters")
pca_elbow <- ElbowPlot(object, ndims = cfg$clustering$npcs)
save_plot(cluster_umap, file.path(cfg$output_dir, "02_clustering", "umap_clusters.png"), width = 8, height = 6)
save_plot(pca_elbow, file.path(cfg$output_dir, "02_clustering", "elbow_plot.png"), width = 7, height = 5)
plot_marker_panel(
  object,
  features = c("IL7R", "LTB", "NKG7", "MS4A1", "LST1", "FCER1A", "PPBP"),
  path = file.path(cfg$output_dir, "02_clustering", "canonical_marker_featureplots.png")
)

message("Summarizing cluster-average marker signals")
avg_expr <- AverageExpression(object, assays = "RNA", slot = "data", verbose = FALSE)$RNA
variable_genes <- intersect(VariableFeatures(object), rownames(avg_expr))
avg_dt <- as.data.table(avg_expr[variable_genes, , drop = FALSE], keep.rownames = "gene")
avg_long <- melt(avg_dt, id.vars = "gene", variable.name = "cluster", value.name = "avg_expression")
avg_long[, cluster := sub("^g", "", cluster)]
top_avg <- avg_long[order(cluster, -avg_expression)][, head(.SD, 10L), by = cluster]
fwrite(avg_long, file.path(cfg$output_dir, "02_clustering", "cluster_average_expression_variable_features.csv"))
fwrite(top_avg, file.path(cfg$output_dir, "02_clustering", "top10_average_features_per_cluster.csv"))

saveRDS(object, file.path(cfg$output_dir, "objects", "02_clustered_seurat.rds"))
message("Saved clustered object")
