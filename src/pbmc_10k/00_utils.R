suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(yaml)
})

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0L) {
    y
  } else {
    x
  }
}

read_config <- function(path) {
  cfg <- yaml::read_yaml(path)
  stopifnot(file.exists(cfg$input$gex_h5))
  stopifnot(file.exists(cfg$input$tcr_contigs))
  stopifnot(file.exists(cfg$input$tcr_clonotypes))
  cfg
}

ensure_dirs <- function(cfg) {
  dirs <- c(
    cfg$output_dir,
    cfg$log_dir,
    file.path(cfg$output_dir, "objects"),
    file.path(cfg$output_dir, "01_qc"),
    file.path(cfg$output_dir, "02_clustering"),
    file.path(cfg$output_dir, "03_annotation"),
    file.path(cfg$output_dir, "04_tcr"),
    file.path(cfg$output_dir, "05_tcell"),
    file.path(cfg$output_dir, "06_report")
  )
  invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))
}

save_plot <- function(plot_obj, path, width = 8, height = 6) {
  ggplot2::ggsave(
    filename = path,
    plot = plot_obj,
    width = width,
    height = height,
    dpi = 300
  )
}

write_metadata <- function(object, path) {
  meta <- as.data.table(object[[]], keep.rownames = "barcode")
  fwrite(meta, path)
}

collapse_unique <- function(values) {
  values <- unique(stats::na.omit(values))
  if (length(values) == 0L) {
    return(NA_character_)
  }
  paste(sort(values), collapse = ";")
}

clone_size_bin <- function(freq) {
  out <- rep("No_TCR", length(freq))
  out[!is.na(freq) & freq == 1] <- "Singleton"
  out[!is.na(freq) & freq >= 2 & freq <= 4] <- "Small"
  out[!is.na(freq) & freq >= 5 & freq <= 9] <- "Medium"
  out[!is.na(freq) & freq >= 10] <- "Large"
  factor(out, levels = c("No_TCR", "Singleton", "Small", "Medium", "Large"))
}

cluster_marker_sets <- function() {
  list(
    CD4_T = c("IL7R", "LTB", "IL32", "LTB", "MALAT1", "TRAC"),
    Cytotoxic_T = c("NKG7", "CCL5", "CTSW", "GZMK", "TRAC", "CD3D"),
    NK = c("NKG7", "GNLY", "PRF1", "TYROBP", "FCGR3A"),
    B = c("MS4A1", "CD79A", "CD74", "HLA-DRA", "BANK1"),
    CD14_Monocyte = c("LST1", "S100A8", "S100A9", "FCN1", "CTSS"),
    FCGR3A_Monocyte = c("FCGR3A", "LST1", "IFITM3", "AIF1", "TYMP"),
    DC = c("FCER1A", "CST3", "HLA-DPA1", "CLEC10A", "SAT1"),
    Platelet = c("PPBP", "PF4", "SDPR", "NRGN", "GNG11"),
    Proliferating = c("MKI67", "TOP2A", "TYMS", "STMN1", "PCNA")
  )
}

score_clusters <- function(avg_expr, marker_sets) {
  cluster_names <- sub("^g", "", colnames(avg_expr))
  scores <- lapply(names(marker_sets), function(set_name) {
    genes <- intersect(marker_sets[[set_name]], rownames(avg_expr))
    if (length(genes) == 0L) {
      return(rep(NA_real_, length(cluster_names)))
    }
    colMeans(avg_expr[genes, , drop = FALSE])
  })
  score_dt <- data.table(cluster = cluster_names)
  for (index in seq_along(marker_sets)) {
    score_dt[[names(marker_sets)[index]]] <- scores[[index]]
  }
  score_dt
}

load_gene_expression_matrix <- function(gex_h5) {
  matrix_or_list <- Read10X_h5(gex_h5)
  if (is.list(matrix_or_list)) {
    if (!"Gene Expression" %in% names(matrix_or_list)) {
      stop("Gene Expression assay not found in H5 input.")
    }
    return(matrix_or_list[["Gene Expression"]])
  }
  matrix_or_list
}

plot_marker_panel <- function(object, features, path, ncol = 3) {
  keep <- intersect(features, rownames(object))
  if (length(keep) == 0L) {
    return(invisible(NULL))
  }
  plot_obj <- FeaturePlot(object, features = keep, ncol = ncol, order = TRUE)
  save_plot(plot_obj, path, width = 12, height = 8)
}

top_markers_by_cluster <- function(markers_dt, top_n = 10L) {
  markers_dt <- as.data.table(markers_dt)
  markers_dt <- markers_dt[p_val_adj < 0.05]
  markers_dt[order(cluster, -avg_log2FC)][, head(.SD, top_n), by = cluster]
}
