source("src/pbmc_10k/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1] %||% "config/pbmc_10k_multi.yaml"
cfg <- read_config(config_path)
ensure_dirs(cfg)

annotated <- readRDS(file.path(cfg$output_dir, "objects", "04_tcr_integrated_seurat.rds"))
t_cells <- readRDS(file.path(cfg$output_dir, "objects", "05_tcell_seurat.rds"))
qc_summary <- fread(file.path(cfg$output_dir, "01_qc", "qc_summary.csv"))
cluster_counts <- fread(file.path(cfg$output_dir, "02_clustering", "cluster_counts.csv"))
tcr_detection <- fread(file.path(cfg$output_dir, "04_tcr", "tcr_detection_by_celltype.csv"))
top_clones <- fread(file.path(cfg$output_dir, "04_tcr", "top_clonotypes.csv"))
tcell_top_clones <- fread(file.path(cfg$output_dir, "05_tcell", "tcell_top_clonotypes.csv"))
tcell_diversity <- fread(file.path(cfg$output_dir, "05_tcell", "tcell_clonotype_diversity.csv"))

celltype_counts <- as.data.table(annotated[[]])[, .N, by = celltype_broad][order(-N)]
tcell_counts <- as.data.table(t_cells[[]])[, .N, by = .(celltype_broad, seurat_clusters)][order(celltype_broad, seurat_clusters)]

report_lines <- c(
  paste0("# PBMC 10k scRNA + scTCR downstream analysis report"),
  "",
  paste0("Generated: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  "",
  "## Input",
  paste0("- GEX matrix: `", cfg$input$gex_h5, "`"),
  paste0("- TCR contigs: `", cfg$input$tcr_contigs, "`"),
  paste0("- TCR clonotypes: `", cfg$input$tcr_clonotypes, "`"),
  "",
  "## QC summary",
  paste0("- Raw cells after CreateSeuratObject: ", qc_summary[metric == "raw_cells", value]),
  paste0("- Cells retained after QC: ", qc_summary[metric == "filtered_cells", value]),
  paste0("- Retained fraction: ", qc_summary[metric == "retained_fraction", value]),
  "",
  "## Clustering summary",
  paste0("- Broad RNA clusters: ", nrow(cluster_counts)),
  paste0("- Broad annotated cell types: ", length(unique(annotated$celltype_broad))),
  "",
  "### Broad cell-type counts"
)
report_lines <- c(report_lines, paste0("- ", celltype_counts$celltype_broad, ": ", celltype_counts$N))
report_lines <- c(
  report_lines,
  "",
  "## TCR integration summary",
  paste0("- Cells with detected productive TCR: ", sum(annotated$has_tcr, na.rm = TRUE)),
  paste0("- Cells with paired TRA/TRB: ", sum(annotated$paired_tcr, na.rm = TRUE)),
  paste0("- Unique clonotypes observed: ", length(unique(stats::na.omit(annotated$clonotype_id)))),
  "",
  "### TCR detection by broad cell type"
)
report_lines <- c(
  report_lines,
  paste0("- ", tcr_detection$celltype_broad, ": ", round(100 * tcr_detection$tcr_rate, 2), "%")
)
report_lines <- c(report_lines, "", "### Global top clonotypes")
if (nrow(top_clones) > 0L) {
  report_lines <- c(
    report_lines,
    paste0("- ", top_clones$clonotype_id, ": ", top_clones$clonotype_frequency, " cells; ", top_clones$clonotype_cdr3s_aa)
  )
} else {
  report_lines <- c(report_lines, "- No clonotypes were detected.")
}

report_lines <- c(
  report_lines,
  "",
  "## T-cell focused clonotype view",
  paste0("- T-cell subset cells: ", ncol(t_cells)),
  paste0("- T-cell subset clusters: ", length(unique(t_cells$seurat_clusters))),
  "",
  "### T-cell subset composition"
)
report_lines <- c(
  report_lines,
  paste0("- ", tcell_counts$celltype_broad, " / cluster ", tcell_counts$seurat_clusters, ": ", tcell_counts$N)
)
report_lines <- c(report_lines, "", "### Top clonotypes in T-cell subset")
if (nrow(tcell_top_clones) > 0L) {
  report_lines <- c(
    report_lines,
    paste0(
      "- ", tcell_top_clones$clonotype_id,
      ": ", tcell_top_clones$n_cells,
      " cells across ", tcell_top_clones$n_tcell_clusters,
      " T-cell clusters; dominant cluster ", tcell_top_clones$dominant_tcell_cluster,
      " (", tcell_top_clones$dominant_cluster_cells, " cells); ",
      tcell_top_clones$clonotype_cdr3s_aa
    )
  )
} else {
  report_lines <- c(report_lines, "- No T-cell clonotypes passed the reporting filters.")
}
report_lines <- c(report_lines, "", "### Top clonotype overlap across T-cell clusters")
if (nrow(tcell_top_clones[n_tcell_clusters > 1]) > 0L) {
  overlap_summary <- tcell_top_clones[n_tcell_clusters > 1]
  report_lines <- c(
    report_lines,
    paste0(
      "- ", overlap_summary$clonotype_id,
      " spans ", overlap_summary$n_tcell_clusters,
      " T-cell clusters; dominant cluster ", overlap_summary$dominant_tcell_cluster,
      " contains ", overlap_summary$dominant_cluster_cells,
      " of ", overlap_summary$n_cells, " cells."
    )
  )
} else {
  report_lines <- c(report_lines, "- Top reported T-cell clonotypes were confined to single T-cell clusters.")
}
report_lines <- c(report_lines, "", "### T-cell clonotype diversity by cluster")
if (nrow(tcell_diversity) > 0L) {
  report_lines <- c(
    report_lines,
    paste0(
      "- Cluster ", tcell_diversity$cluster,
      ": ", tcell_diversity$n_tcr_cells, "/", tcell_diversity$n_cells, " TCR+ cells, ",
      tcell_diversity$unique_clonotypes, " clonotypes, Shannon ", sprintf("%.3f", tcell_diversity$shannon_diversity),
      ", expanded fraction ", round(100 * tcell_diversity$expanded_fraction, 1), "%"
    )
  )
} else {
  report_lines <- c(report_lines, "- No cluster-level diversity metrics were available.")
}

report_lines <- c(
  report_lines,
  "",
  "## Output locations",
  paste0("- QC outputs: `", file.path(cfg$output_dir, "01_qc"), "`"),
  paste0("- Clustering outputs: `", file.path(cfg$output_dir, "02_clustering"), "`"),
  paste0("- Annotation outputs: `", file.path(cfg$output_dir, "03_annotation"), "`"),
  paste0("- TCR outputs: `", file.path(cfg$output_dir, "04_tcr"), "`"),
  paste0("- T-cell outputs: `", file.path(cfg$output_dir, "05_tcell"), "`"),
  paste0("- Saved objects: `", file.path(cfg$output_dir, "objects"), "`"),
  paste0("- Step logs: `", cfg$log_dir, "`")
)

writeLines(report_lines, con = file.path(cfg$output_dir, "06_report", "analysis_summary.md"))
fwrite(celltype_counts, file.path(cfg$output_dir, "06_report", "celltype_counts.csv"))
message("Wrote summary report")
