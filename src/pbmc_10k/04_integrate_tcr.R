source("src/pbmc_10k/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1] %||% "config/pbmc_10k_multi.yaml"
cfg <- read_config(config_path)
ensure_dirs(cfg)

object <- readRDS(file.path(cfg$output_dir, "objects", "03_annotated_seurat.rds"))
contigs <- fread(cfg$input$tcr_contigs)
clonotypes <- fread(cfg$input$tcr_clonotypes)

message("Summarizing productive high-confidence TCR calls")
contigs_use <- contigs[
  is_cell == TRUE &
    high_confidence == TRUE &
    productive == TRUE &
    chain %in% c("TRA", "TRB")
]

tcr_meta <- contigs_use[
  ,
  .(
    has_tcr = TRUE,
    paired_tcr = all(c("TRA", "TRB") %in% unique(chain)),
    clonotype_id = {
      ids <- unique(na.omit(raw_clonotype_id))
      if (length(ids) == 1L) ids else NA_character_
    },
    clonotype_multiplet = uniqueN(na.omit(raw_clonotype_id)) > 1L,
    tcr_chain_count = uniqueN(chain),
    tcr_chain_types = collapse_unique(chain),
    cdr3s = collapse_unique(paste(chain, cdr3, sep = ":")),
    total_tcr_reads = sum(reads, na.rm = TRUE),
    total_tcr_umis = sum(umis, na.rm = TRUE)
  ),
  by = barcode
]

clonotypes_use <- clonotypes[, .(
  clonotype_id,
  clonotype_frequency = frequency,
  clonotype_proportion = proportion,
  clonotype_cdr3s_aa = cdr3s_aa
)]
tcr_meta <- merge(tcr_meta, clonotypes_use, by = "clonotype_id", all.x = TRUE)
tcr_meta$clone_size_bin <- clone_size_bin(tcr_meta$clonotype_frequency)
tcr_meta_df <- as.data.frame(tcr_meta)
rownames(tcr_meta_df) <- tcr_meta_df$barcode

common <- intersect(colnames(object), rownames(tcr_meta_df))
object <- AddMetaData(object, metadata = tcr_meta_df[common, setdiff(colnames(tcr_meta_df), "barcode"), drop = FALSE])
object$has_tcr[is.na(object$has_tcr)] <- FALSE
object$paired_tcr[is.na(object$paired_tcr)] <- FALSE
object$clonotype_multiplet[is.na(object$clonotype_multiplet)] <- FALSE
object$clone_size_bin <- clone_size_bin(object$clonotype_frequency)
object$has_tcr_label <- ifelse(object$has_tcr, "TCR_detected", "No_TCR")

meta_dt <- as.data.table(object[[]], keep.rownames = "barcode")
fwrite(meta_dt, file.path(cfg$output_dir, "04_tcr", "seurat_metadata_with_tcr.csv"))

tcr_detection <- meta_dt[, .(
  n_cells = .N,
  n_tcr = sum(has_tcr, na.rm = TRUE),
  tcr_rate = round(mean(has_tcr, na.rm = TRUE), 4)
), by = celltype_broad][order(-tcr_rate)]
fwrite(tcr_detection, file.path(cfg$output_dir, "04_tcr", "tcr_detection_by_celltype.csv"))

top_clones <- unique(meta_dt[!is.na(clonotype_id) & !clonotype_multiplet, .(
  clonotype_id,
  clonotype_frequency,
  clonotype_proportion,
  clonotype_cdr3s_aa
)])[order(-clonotype_frequency)][1:min(cfg$clonotype$top_n, .N)]
fwrite(top_clones, file.path(cfg$output_dir, "04_tcr", "top_clonotypes.csv"))

clone_by_celltype <- meta_dt[
  !is.na(clonotype_id) & !clonotype_multiplet,
  .N,
  by = .(celltype_broad, clone_size_bin)
]
fwrite(clone_by_celltype, file.path(cfg$output_dir, "04_tcr", "clone_size_by_celltype.csv"))

has_tcr_plot <- DimPlot(object, reduction = "umap", group.by = "has_tcr_label") +
  ggtitle("TCR detection on RNA UMAP")
clone_size_plot <- DimPlot(object, reduction = "umap", group.by = "clone_size_bin") +
  ggtitle("Clone size on RNA UMAP")
save_plot(has_tcr_plot, file.path(cfg$output_dir, "04_tcr", "umap_tcr_detected.png"), width = 8, height = 6)
save_plot(clone_size_plot, file.path(cfg$output_dir, "04_tcr", "umap_clone_size.png"), width = 8, height = 6)

detection_plot <- ggplot(tcr_detection, aes(x = reorder(celltype_broad, -tcr_rate), y = tcr_rate, fill = celltype_broad)) +
  geom_col(show.legend = FALSE) +
  coord_flip() +
  labs(x = NULL, y = "TCR detection rate", title = "TCR detection by broad cell type") +
  theme_bw(base_size = 11)
save_plot(detection_plot, file.path(cfg$output_dir, "04_tcr", "tcr_detection_by_celltype.png"), width = 8, height = 5)

if (nrow(top_clones) > 0L) {
  top_clone_plot <- ggplot(top_clones, aes(x = reorder(clonotype_id, clonotype_frequency), y = clonotype_frequency)) +
    geom_col(fill = "#2F6B8A") +
    coord_flip() +
    labs(x = NULL, y = "Cells", title = "Top expanded clonotypes") +
    theme_bw(base_size = 11)
  save_plot(top_clone_plot, file.path(cfg$output_dir, "04_tcr", "top_clonotypes_barplot.png"), width = 8, height = 5)
}

saveRDS(object, file.path(cfg$output_dir, "objects", "04_tcr_integrated_seurat.rds"))
message("Saved TCR-integrated object")
