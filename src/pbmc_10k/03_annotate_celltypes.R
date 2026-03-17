source("src/pbmc_10k/00_utils.R")

assign_label_from_markers <- function(top_markers, score_row) {
  genes <- unique(trimws(unlist(strsplit(top_markers %||% "", ","))))

  if (any(genes %in% c("PPBP", "PF4", "SDPR", "NRGN", "TUBB1"))) {
    return("Platelet")
  }
  if (any(genes %in% c("MKI67", "TOP2A", "TYMS", "STMN1", "PCNA"))) {
    return("Proliferating")
  }
  if (any(genes %in% c("FCER1A", "CLEC10A", "CST3", "IFI30", "HLA-DPA1"))) {
    return("DC")
  }
  if (any(genes %in% c("FCGR3A", "IFITM3", "AIF1", "TYMP")) ||
      (("LST1" %in% genes) && ("FCGR3A" %in% genes))) {
    return("FCGR3A_Monocyte")
  }
  if (any(genes %in% c("S100A8", "S100A9", "FCN1", "LYZ", "CTSS"))) {
    return("CD14_Monocyte")
  }
  if (any(genes %in% c("MS4A1", "CD79A", "BANK1")) ||
      (("CD74" %in% genes || "HLA-DRA" %in% genes) && !any(genes %in% c("LYZ", "S100A8", "S100A9")))) {
    return("B")
  }
  if (any(genes %in% c("GNLY", "NKG7", "PRF1", "CTSW", "CCL5", "GZMK"))) {
    if (("FCGR3A" %in% genes) && !(any(genes %in% c("IL7R", "LTB", "TRAC", "CD3D", "IL32")))) {
      return("NK")
    }
    return("Cytotoxic_T")
  }
  if (any(genes %in% c("IL7R", "LTB", "TCF7", "TRAC", "CD3D", "IL32", "LINC00861", "ETS1"))) {
    return("CD4_T")
  }

  score_labels <- setdiff(names(score_row), c("cluster", "top_markers", "suggested_label", "refined_label"))
  best_score <- score_labels[which.max(as.numeric(score_row[, ..score_labels]))]
  best_score
}

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1] %||% "config/pbmc_10k_multi.yaml"
cfg <- read_config(config_path)
ensure_dirs(cfg)

object <- readRDS(file.path(cfg$output_dir, "objects", "02_clustered_seurat.rds"))
markers_dt <- fread(file.path(cfg$output_dir, "02_clustering", "top10_average_features_per_cluster.csv"))
avg_expr <- AverageExpression(object, assays = "RNA", slot = "data", verbose = FALSE)$RNA
score_dt <- score_clusters(avg_expr, cluster_marker_sets())
score_dt[, cluster := as.character(cluster)]

top_marker_text <- markers_dt[
  order(cluster, -avg_expression)
][
  ,
  .(top_markers = paste(gene, collapse = ", ")),
  by = cluster
]
top_marker_text[, cluster := as.character(cluster)]
annotation_dt <- merge(score_dt, top_marker_text, by = "cluster", all.x = TRUE)
score_columns <- setdiff(names(score_dt), "cluster")
annotation_dt[, suggested_label := score_columns[max.col(.SD, ties.method = "first")], .SDcols = score_columns]
annotation_dt[, refined_label := vapply(
  seq_len(.N),
  function(i) assign_label_from_markers(top_markers[i], annotation_dt[i]),
  character(1)
)]
fwrite(annotation_dt, file.path(cfg$output_dir, "03_annotation", "cluster_annotation_scores.csv"))

cluster_to_label <- annotation_dt[, .(cluster, refined_label)]
label_map <- setNames(cluster_to_label$refined_label, cluster_to_label$cluster)
object$celltype_broad <- unname(label_map[as.character(object$seurat_clusters)])
object$celltype_broad[is.na(object$celltype_broad)] <- "Unassigned"

annotation_table <- as.data.table(object[[]])[, .N, by = .(cluster = seurat_clusters, celltype_broad)]
setnames(annotation_table, "N", "n_cells")
fwrite(annotation_table, file.path(cfg$output_dir, "03_annotation", "cluster_celltype_table.csv"))

celltype_umap <- DimPlot(object, reduction = "umap", group.by = "celltype_broad", label = TRUE, repel = TRUE) +
  ggtitle("PBMC 10k broad cell types")
save_plot(celltype_umap, file.path(cfg$output_dir, "03_annotation", "umap_celltypes.png"), width = 9, height = 7)

dot_features <- intersect(
  c("IL7R", "LTB", "NKG7", "GNLY", "MS4A1", "CD79A", "LST1", "FCER1A", "PPBP", "MKI67"),
  rownames(object)
)
if (length(dot_features) > 0L) {
  dot_plot <- DotPlot(object, features = dot_features, group.by = "celltype_broad") +
    RotatedAxis()
  save_plot(dot_plot, file.path(cfg$output_dir, "03_annotation", "celltype_marker_dotplot.png"), width = 11, height = 6)
}

saveRDS(object, file.path(cfg$output_dir, "objects", "03_annotated_seurat.rds"))
message("Saved annotated object")
