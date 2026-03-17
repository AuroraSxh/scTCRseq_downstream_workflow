source("src/pbmc_10k/00_utils.R")

args <- commandArgs(trailingOnly = TRUE)
config_path <- args[1] %||% "config/pbmc_10k_multi.yaml"
cfg <- read_config(config_path)
ensure_dirs(cfg)

message("Loading gene expression matrix")
gex <- load_gene_expression_matrix(cfg$input$gex_h5)
object <- CreateSeuratObject(
  counts = gex,
  project = cfg$sample_name,
  min.cells = 3,
  min.features = cfg$qc$min_features
)
object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
object$sample_id <- cfg$sample_name

raw_meta <- as.data.table(object[[]], keep.rownames = "barcode")
fwrite(raw_meta, file.path(cfg$output_dir, "01_qc", "raw_metadata.csv"))

qc_violin <- VlnPlot(
  object,
  features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
  ncol = 3,
  pt.size = 0.05
)
qc_scatter <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt")
save_plot(qc_violin, file.path(cfg$output_dir, "01_qc", "qc_violin_raw.png"), width = 12, height = 4.5)
save_plot(qc_scatter, file.path(cfg$output_dir, "01_qc", "qc_scatter_raw.png"), width = 10, height = 5)

message("Applying QC filters")
filtered <- subset(
  object,
  subset =
    nFeature_RNA >= cfg$qc$min_features &
    nFeature_RNA <= cfg$qc$max_features &
    nCount_RNA <= cfg$qc$max_counts &
    percent.mt <= cfg$qc$max_percent_mt
)

qc_summary <- data.table(
  metric = c("raw_cells", "filtered_cells", "retained_fraction"),
  value = c(ncol(object), ncol(filtered), round(ncol(filtered) / ncol(object), 4))
)
fwrite(qc_summary, file.path(cfg$output_dir, "01_qc", "qc_summary.csv"))
write_metadata(filtered, file.path(cfg$output_dir, "01_qc", "filtered_metadata.csv"))
saveRDS(filtered, file.path(cfg$output_dir, "objects", "01_filtered_seurat.rds"))
message("Saved filtered object with ", ncol(filtered), " cells")
