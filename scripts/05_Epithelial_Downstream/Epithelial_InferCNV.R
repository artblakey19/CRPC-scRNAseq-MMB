# Epithelial_InferCNV.R
# Run inferCNV on filtered/annotated epithelial cells.
#
# Observations:
#   - Results/05_Epithelial_Downstream/epi_annotated.rds
# References:
#   - Non-epithelial cells from Results/01_Integrated/combined_CRPC.rds
#
# Writes:
#   - Results/05_Epithelial_Downstream/InferCNV/inputs/
#   - Results/05_Epithelial_Downstream/InferCNV/infercnv_run/
#   - Results/05_Epithelial_Downstream/InferCNV/epithelial_cnv_scores.tsv
#   - Results/05_Epithelial_Downstream/InferCNV/UMAP_infercnv_cnv_score.png
#   - Results/05_Epithelial_Downstream/InferCNV/VlnPlot_infercnv_cnv_score.png

suppressPackageStartupMessages({
    library(Seurat)
    library(infercnv)
    library(Matrix)
    library(dplyr)
    library(ggplot2)
})

COMBINED_RDS <- "Results/01_Integrated/combined_CRPC.rds"
EPI_RDS      <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
GENE_ORDER   <- "Resources/hg38_gencode_v27.txt"
OUT_DIR      <- "Results/05_Epithelial_Downstream/InferCNV"
INPUT_DIR    <- file.path(OUT_DIR, "inputs")
RUN_DIR      <- file.path(OUT_DIR, "infercnv_run")

dir.create(INPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(RUN_DIR, showWarnings = FALSE, recursive = TRUE)

threads <- as.integer(Sys.getenv("INFERCNV_THREADS", "16"))
if (is.na(threads) || threads < 1) threads <- 16

run_hmm <- tolower(Sys.getenv("INFERCNV_HMM", "false")) %in%
    c("1", "true", "t", "yes", "y")

clean_group <- function(x) {
    x <- gsub("[^A-Za-z0-9]+", "_", x)
    x <- gsub("^_+|_+$", "", x)
    x[x == "" | is.na(x)] <- "NA"
    x
}

message("Loading input objects")
combined <- readRDS(COMBINED_RDS)
epi <- readRDS(EPI_RDS)

required_combined <- c("orig.ident", "celltype")
missing_combined <- setdiff(required_combined, colnames(combined@meta.data))
if (length(missing_combined) > 0) {
    stop("combined object is missing metadata columns: ",
         paste(missing_combined, collapse = ", "), call. = FALSE)
}
if (!"annotation" %in% colnames(epi@meta.data)) {
    stop("epithelial object is missing metadata column: annotation",
         call. = FALSE)
}

epi_cells <- intersect(colnames(epi), colnames(combined))
if (length(epi_cells) == 0) {
    stop("No overlap between epithelial and combined object cell names",
         call. = FALSE)
}

ref_cells <- colnames(combined)[
    !is.na(combined$celltype) & combined$celltype != "Epithelial"
]
if (length(ref_cells) == 0) {
    stop("No non-epithelial reference cells found in combined object",
         call. = FALSE)
}

analysis_cells <- c(epi_cells, ref_cells)

message("Building inferCNV cell annotations")
epi_meta <- data.frame(
    cell = epi_cells,
    class = "observation",
    sample = as.character(combined$orig.ident[epi_cells]),
    annotation = as.character(epi$annotation[epi_cells]),
    stringsAsFactors = FALSE
)
epi_meta$group <- paste0(
    "obs_", clean_group(epi_meta$annotation), "__", clean_group(epi_meta$sample)
)

ref_meta <- data.frame(
    cell = ref_cells,
    class = "reference",
    sample = as.character(combined$orig.ident[ref_cells]),
    annotation = as.character(combined$celltype[ref_cells]),
    stringsAsFactors = FALSE
)
ref_meta$group <- paste0("ref_", clean_group(ref_meta$annotation))

cell_meta <- rbind(epi_meta, ref_meta)
cell_meta <- cell_meta[match(analysis_cells, cell_meta$cell), ]
stopifnot(identical(cell_meta$cell, analysis_cells))

ref_group_names <- sort(unique(ref_meta$group))

annotations_file <- file.path(INPUT_DIR, "cell_annotations.tsv")
cell_meta_file <- file.path(INPUT_DIR, "cell_metadata.tsv")
group_counts_file <- file.path(INPUT_DIR, "cell_group_counts.tsv")

write.table(
    cell_meta[, c("cell", "group")],
    file = annotations_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
)
write.table(
    cell_meta,
    file = cell_meta_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
)
write.table(
    as.data.frame(table(class = cell_meta$class, group = cell_meta$group)) |>
        dplyr::filter(Freq > 0),
    file = group_counts_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
)

message("Preparing gene order")
gene_order <- read.delim(
    GENE_ORDER, header = FALSE, stringsAsFactors = FALSE,
    col.names = c("gene", "chr", "start", "end")
)
canonical_chr <- paste0("chr", c(1:22, "X", "Y"))
gene_order <- gene_order[
    gene_order$chr %in% canonical_chr & !duplicated(gene_order$gene),
]

DefaultAssay(combined) <- "RNA"
counts <- GetAssayData(combined, assay = "RNA", layer = "counts")
counts <- counts[, analysis_cells, drop = FALSE]

gene_order <- gene_order[gene_order$gene %in% rownames(counts), ]
gene_order$chr <- factor(gene_order$chr, levels = canonical_chr)
gene_order <- gene_order[order(gene_order$chr, gene_order$start,
                               gene_order$end), ]

gene_order_file <- file.path(INPUT_DIR, "gene_order_hg38_gencode_v27.tsv")
write.table(
    gene_order[, c("gene", "chr", "start", "end")],
    file = gene_order_file,
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = FALSE
)

counts <- counts[gene_order$gene, , drop = FALSE]

summary_file <- file.path(INPUT_DIR, "run_summary.txt")
writeLines(c(
    paste("combined_rds:", COMBINED_RDS),
    paste("epithelial_rds:", EPI_RDS),
    paste("gene_order_source:", GENE_ORDER),
    paste("observation_cells:", nrow(epi_meta)),
    paste("reference_cells:", nrow(ref_meta)),
    paste("total_cells:", ncol(counts)),
    paste("genes_in_count_matrix:", nrow(counts)),
    paste("reference_groups:", paste(ref_group_names, collapse = ",")),
    paste("cutoff:", 0.1),
    paste("denoise:", TRUE),
    paste("HMM:", run_hmm),
    paste("threads:", threads)
), con = summary_file)

message("Observation cells: ", nrow(epi_meta))
message("Reference cells: ", nrow(ref_meta))
message("Genes used before inferCNV chr exclusions: ", nrow(counts))
message("Reference groups: ", paste(ref_group_names, collapse = ", "))
message("Running inferCNV; output directory: ", RUN_DIR)

set.seed(42)
infercnv_obj <- CreateInfercnvObject(
    raw_counts_matrix = counts,
    annotations_file = annotations_file,
    gene_order_file = gene_order_file,
    ref_group_names = ref_group_names
)

infercnv_obj <- infercnv::run(
    infercnv_obj,
    cutoff = 0.1,
    out_dir = RUN_DIR,
    cluster_by_groups = TRUE,
    cluster_references = TRUE,
    analysis_mode = "samples",
    denoise = TRUE,
    HMM = run_hmm,
    num_threads = threads,
    output_format = "png",
    png_res = 300,
    plot_steps = FALSE,
    save_rds = TRUE,
    save_final_rds = TRUE,
    useRaster = TRUE
)

message("Computing per-epithelial inferCNV burden score")
obs_cells <- colnames(infercnv_obj@expr.data)[
    colnames(infercnv_obj@expr.data) %in% epi_cells
]
obs_expr <- infercnv_obj@expr.data[, obs_cells, drop = FALSE]
cnv_score <- Matrix::colMeans(abs(obs_expr - 1))

score_df <- data.frame(
    cell = names(cnv_score),
    infercnv_cnv_score = as.numeric(cnv_score),
    sample = as.character(combined$orig.ident[names(cnv_score)]),
    annotation = as.character(epi$annotation[names(cnv_score)]),
    stringsAsFactors = FALSE
)
write.table(
    score_df,
    file = file.path(OUT_DIR, "epithelial_cnv_scores.tsv"),
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
)

score_summary <- score_df |>
    dplyr::group_by(annotation) |>
    dplyr::summarise(
        n = dplyr::n(),
        median = stats::median(infercnv_cnv_score),
        mean = mean(infercnv_cnv_score),
        q75 = stats::quantile(infercnv_cnv_score, 0.75),
        max = max(infercnv_cnv_score),
        .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(median))
write.table(
    score_summary,
    file = file.path(OUT_DIR, "epithelial_cnv_score_by_annotation.tsv"),
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
)

sample_summary <- score_df |>
    dplyr::group_by(sample) |>
    dplyr::summarise(
        n = dplyr::n(),
        median = stats::median(infercnv_cnv_score),
        mean = mean(infercnv_cnv_score),
        q75 = stats::quantile(infercnv_cnv_score, 0.75),
        max = max(infercnv_cnv_score),
        .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(median))
write.table(
    sample_summary,
    file = file.path(OUT_DIR, "epithelial_cnv_score_by_sample.tsv"),
    sep = "\t", quote = FALSE,
    row.names = FALSE, col.names = TRUE
)

epi$infercnv_cnv_score <- NA_real_
epi$infercnv_cnv_score[names(cnv_score)] <- as.numeric(cnv_score)

p_umap <- FeaturePlot(
    epi, features = "infercnv_cnv_score", pt.size = 0.25
) +
    scale_color_viridis_c(option = "magma", na.value = "grey85") +
    ggtitle("Epithelial inferCNV burden score")
ggsave(
    file.path(OUT_DIR, "UMAP_infercnv_cnv_score.png"),
    plot = p_umap, width = 10, height = 8, bg = "white"
)

p_vln <- VlnPlot(
    epi, features = "infercnv_cnv_score", group.by = "annotation",
    pt.size = 0
) +
    RotatedAxis() +
    ggtitle("Epithelial inferCNV burden score by annotation")
ggsave(
    file.path(OUT_DIR, "VlnPlot_infercnv_cnv_score.png"),
    plot = p_vln, width = 13, height = 7, bg = "white"
)

saveRDS(
    epi,
    file.path(OUT_DIR, "epi_annotated_with_infercnv_score.rds")
)

message("inferCNV complete")
