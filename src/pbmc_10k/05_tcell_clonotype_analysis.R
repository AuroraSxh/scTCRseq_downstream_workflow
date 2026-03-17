source("src/pbmc_10k/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1] %||% "config/pbmc_10k_multi.yaml"
cfg <- read_config(config_path)
ensure_dirs(cfg)

object <- readRDS(file.path(cfg$output_dir, "objects", "04_tcr_integrated_seurat.rds"))
tcell_labels <- c("CD4_T", "Cytotoxic_T")
t_cells <- subset(object, subset = celltype_broad %in% tcell_labels)

legacy_outputs <- file.path(
  cfg$output_dir,
  "05_tcell",
  c(
    "tcell_clone_state_table.csv",
    "umap_tcell_clone_size.png",
    "clone_size_fraction_by_tcell_cluster.png"
  )
)
legacy_outputs <- legacy_outputs[file.exists(legacy_outputs)]
if (length(legacy_outputs) > 0L) {
  invisible(file.remove(legacy_outputs))
}

if (ncol(t_cells) < 50L) {
  stop("Too few T cells detected for T-cell focused analysis.")
}

message("Re-embedding T-cell compartment")
t_cells <- NormalizeData(t_cells, verbose = FALSE)
t_cells <- FindVariableFeatures(t_cells, selection.method = "vst", nfeatures = min(cfg$clustering$nfeatures, 2500), verbose = FALSE)
vars_to_regress <- if (isTRUE(cfg$clustering$regress_percent_mt)) "percent.mt" else NULL
t_cells <- ScaleData(t_cells, vars.to.regress = vars_to_regress, verbose = FALSE)
t_cells <- RunPCA(t_cells, npcs = min(20, cfg$clustering$npcs), verbose = FALSE)
t_cells <- FindNeighbors(t_cells, dims = 1:20, verbose = FALSE)
t_cells <- FindClusters(t_cells, resolution = 0.4, verbose = FALSE)
t_cells <- RunUMAP(t_cells, dims = 1:20, verbose = FALSE)

meta_dt <- as.data.table(t_cells[[]], keep.rownames = "barcode")
meta_dt[, seurat_clusters := as.character(seurat_clusters)]
cluster_levels <- sort(unique(meta_dt$seurat_clusters))

tcell_cluster_counts <- as.data.table(table(meta_dt$seurat_clusters))
setnames(tcell_cluster_counts, c("cluster", "n_cells"))
tcell_cluster_counts[, cluster := as.character(cluster)]
fwrite(tcell_cluster_counts, file.path(cfg$output_dir, "05_tcell", "tcell_cluster_counts.csv"))

clone_meta <- meta_dt[
  !is.na(clonotype_id) & !clonotype_multiplet,
  .(barcode, celltype_broad, seurat_clusters, clonotype_id, clonotype_cdr3s_aa)
]

top_clones <- data.table(
  top_rank = integer(),
  clonotype_id = character(),
  clonotype_cdr3s_aa = character(),
  n_cells = integer(),
  n_tcell_clusters = integer(),
  tcell_labels = character(),
  dominant_tcell_cluster = character(),
  dominant_cluster_cells = integer()
)
diversity_dt <- data.table(
  cluster = character(),
  n_cells = integer(),
  n_tcr_cells = integer(),
  tcr_detection_rate = numeric(),
  unique_clonotypes = integer(),
  shannon_diversity = numeric(),
  evenness = numeric(),
  expanded_clonotypes = integer(),
  expanded_cells = integer(),
  expanded_fraction = numeric(),
  dominant_clone_cells = integer(),
  dominant_clone_fraction = numeric()
)
top_overlap <- data.table(
  clonotype_id = character(),
  top_rank = integer(),
  seurat_clusters = character(),
  N = integer(),
  n_tcr_cells = integer(),
  cluster_fraction = numeric()
)

if (nrow(clone_meta) > 0L) {
  clone_cluster_counts <- clone_meta[, .N, by = .(clonotype_id, clonotype_cdr3s_aa, seurat_clusters)]
  dominant_cluster_dt <- clone_cluster_counts[
    order(clonotype_id, -N, seurat_clusters),
    .SD[1],
    by = clonotype_id
  ][, .(
    clonotype_id,
    dominant_tcell_cluster = seurat_clusters,
    dominant_cluster_cells = N
  )]

  top_clones <- clone_meta[
    ,
    .(
      n_cells = .N,
      n_tcell_clusters = uniqueN(seurat_clusters),
      tcell_labels = collapse_unique(celltype_broad)
    ),
    by = .(clonotype_id, clonotype_cdr3s_aa)
  ][order(-n_cells, clonotype_id)]
  if (any(top_clones$n_cells > 1L)) {
    top_clones <- top_clones[n_cells > 1L]
  }
  top_clones <- merge(top_clones, dominant_cluster_dt, by = "clonotype_id", all.x = TRUE)
  top_clones <- top_clones[1:min(cfg$clonotype$top_n, .N)]
  top_clones[, top_rank := seq_len(.N)]
  setcolorder(
    top_clones,
    c(
      "top_rank",
      "clonotype_id",
      "clonotype_cdr3s_aa",
      "n_cells",
      "n_tcell_clusters",
      "tcell_labels",
      "dominant_tcell_cluster",
      "dominant_cluster_cells"
    )
  )

  cluster_clone_sizes <- clone_meta[, .N, by = .(seurat_clusters, clonotype_id)]
  diversity_dt <- cluster_clone_sizes[
    ,
    {
      counts <- N
      total <- sum(counts)
      p <- counts / total
      richness <- .N
      shannon <- -sum(p * log(p))
      evenness <- if (richness > 1L) shannon / log(richness) else 1
      expanded_cells <- sum(counts[counts > 1L])
      dominant_clone_cells <- max(counts)
      list(
        n_tcr_cells = total,
        unique_clonotypes = richness,
        shannon_diversity = round(shannon, 3),
        evenness = round(evenness, 3),
        expanded_clonotypes = sum(counts > 1L),
        expanded_cells = expanded_cells,
        expanded_fraction = round(expanded_cells / total, 3),
        dominant_clone_cells = dominant_clone_cells,
        dominant_clone_fraction = round(dominant_clone_cells / total, 3)
      )
    },
    by = seurat_clusters
  ]
  diversity_dt <- merge(
    tcell_cluster_counts,
    diversity_dt,
    by.x = "cluster",
    by.y = "seurat_clusters",
    all.x = TRUE
  )
  diversity_dt[is.na(n_tcr_cells), c(
    "n_tcr_cells",
    "unique_clonotypes",
    "shannon_diversity",
    "evenness",
    "expanded_clonotypes",
    "expanded_cells",
    "expanded_fraction",
    "dominant_clone_cells",
    "dominant_clone_fraction"
  ) := list(0L, 0L, 0, 0, 0L, 0L, 0, 0L, 0)]
  diversity_dt[, tcr_detection_rate := round(ifelse(n_cells > 0, n_tcr_cells / n_cells, 0), 3)]
  setcolorder(
    diversity_dt,
    c(
      "cluster",
      "n_cells",
      "n_tcr_cells",
      "tcr_detection_rate",
      "unique_clonotypes",
      "shannon_diversity",
      "evenness",
      "expanded_clonotypes",
      "expanded_cells",
      "expanded_fraction",
      "dominant_clone_cells",
      "dominant_clone_fraction"
    )
  )

  if (nrow(top_clones) > 0L) {
    cluster_totals <- clone_meta[, .(n_tcr_cells = .N), by = seurat_clusters]
    top_ids <- top_clones$clonotype_id
    top_overlap <- merge(
      CJ(seurat_clusters = cluster_levels, clonotype_id = top_ids, unique = TRUE),
      clone_cluster_counts[clonotype_id %in% top_ids, .(seurat_clusters, clonotype_id, N)],
      by = c("seurat_clusters", "clonotype_id"),
      all.x = TRUE
    )
    top_overlap[is.na(N), N := 0L]
    top_overlap <- merge(top_overlap, top_clones[, .(clonotype_id, top_rank)], by = "clonotype_id", all.x = TRUE)
    top_overlap <- merge(top_overlap, cluster_totals, by = "seurat_clusters", all.x = TRUE)
    top_overlap[is.na(n_tcr_cells), n_tcr_cells := 0L]
    top_overlap[, cluster_fraction := round(ifelse(n_tcr_cells > 0, N / n_tcr_cells, 0), 4)]
    setcolorder(top_overlap, c("clonotype_id", "top_rank", "seurat_clusters", "N", "n_tcr_cells", "cluster_fraction"))
    setorder(top_overlap, top_rank, seurat_clusters)
  }
}

fwrite(top_clones, file.path(cfg$output_dir, "05_tcell", "tcell_top_clonotypes.csv"))
fwrite(diversity_dt, file.path(cfg$output_dir, "05_tcell", "tcell_clonotype_diversity.csv"))
fwrite(top_overlap, file.path(cfg$output_dir, "05_tcell", "tcell_top_clonotype_overlap.csv"))

top_ids <- top_clones$clonotype_id
t_cells$top_clonotype <- ifelse(!is.na(t_cells$clonotype_id) & t_cells$clonotype_id %in% top_ids, t_cells$clonotype_id, "Other_T_cells")
t_cells$top_clonotype <- factor(t_cells$top_clonotype, levels = c(top_ids, "Other_T_cells"))

tcell_umap <- DimPlot(t_cells, reduction = "umap", group.by = "seurat_clusters", label = TRUE) +
  ggtitle("T-cell subset clusters")
top_clone_plot <- DimPlot(t_cells, reduction = "umap", group.by = "top_clonotype") +
  ggtitle("Top clonotypes in T-cell subset")
save_plot(tcell_umap, file.path(cfg$output_dir, "05_tcell", "umap_tcell_clusters.png"), width = 8, height = 6)
save_plot(top_clone_plot, file.path(cfg$output_dir, "05_tcell", "umap_tcell_top_clonotypes.png"), width = 10, height = 7)

if (nrow(top_clones) > 0L) {
  expansion_plot <- ggplot(top_clones, aes(x = reorder(clonotype_id, n_cells), y = n_cells)) +
    geom_col(fill = "#2F6B8A") +
    geom_text(aes(label = n_cells), hjust = -0.1, size = 3) +
    coord_flip() +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    labs(x = NULL, y = "Cells", title = "Top clonotype expansion in T-cell subset") +
    theme_bw(base_size = 11)
  save_plot(expansion_plot, file.path(cfg$output_dir, "05_tcell", "top_clonotype_expansion_barplot.png"), width = 8, height = 5)

  overlap_plot <- ggplot(
    top_overlap,
    aes(
      x = factor(seurat_clusters, levels = cluster_levels),
      y = factor(clonotype_id, levels = rev(top_ids)),
      fill = N
    )
  ) +
    geom_tile(color = "white", linewidth = 0.2) +
    scale_fill_gradient(low = "#F4F1E8", high = "#2F6B8A") +
    labs(x = "T-cell cluster", y = "Top clonotype", fill = "Cells", title = "Top clonotype overlap across T-cell clusters") +
    theme_bw(base_size = 11)
  save_plot(overlap_plot, file.path(cfg$output_dir, "05_tcell", "top_clonotype_overlap_by_tcell_cluster.png"), width = 8, height = 5.5)
}

if (nrow(diversity_dt) > 0L) {
  diversity_plot <- ggplot(
    diversity_dt,
    aes(x = factor(cluster, levels = cluster_levels), y = shannon_diversity, fill = expanded_fraction)
  ) +
    geom_col() +
    geom_text(aes(label = paste0("n=", unique_clonotypes)), vjust = -0.35, size = 3) +
    scale_fill_gradient(low = "#DCE8D0", high = "#B24A3A") +
    scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    labs(
      x = "T-cell cluster",
      y = "Shannon diversity",
      fill = "Expanded fraction",
      title = "Clonotype diversity across T-cell clusters"
    ) +
    theme_bw(base_size = 11)
  save_plot(diversity_plot, file.path(cfg$output_dir, "05_tcell", "clonotype_diversity_by_tcell_cluster.png"), width = 8, height = 5)
}

saveRDS(t_cells, file.path(cfg$output_dir, "objects", "05_tcell_seurat.rds"))
message("Saved T-cell focused object")
