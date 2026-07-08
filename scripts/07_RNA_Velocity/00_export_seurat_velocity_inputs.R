#!/usr/bin/env Rscript
# Export epithelial Seurat metadata and UMAP coordinates for scVelo.

suppressPackageStartupMessages({
    library(Seurat)
})

IN_RDS <- Sys.getenv(
    "SEURAT_RDS",
    "Results/05_Epithelial_Downstream/epi_annotated.rds"
)
OUT_DIR <- Sys.getenv("VELOCITY_OUT_DIR", "Results/07_RNA_Velocity")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

if (!file.exists(IN_RDS)) {
    stop("Missing Seurat RDS: ", IN_RDS, call. = FALSE)
}

message("Loading ", IN_RDS)
epi <- readRDS(IN_RDS)

required_meta <- c("orig.ident", "seurat_clusters", "annotation")
missing_meta <- setdiff(required_meta, colnames(epi@meta.data))
if (length(missing_meta) > 0) {
    stop(
        "Missing metadata column(s): ",
        paste(missing_meta, collapse = ", "),
        call. = FALSE
    )
}

if (!"umap" %in% Reductions(epi)) {
    stop("The Seurat object does not contain a 'umap' reduction.", call. = FALSE)
}

cells <- colnames(epi)
umap <- Embeddings(epi, "umap")
umap <- umap[cells, seq_len(min(2, ncol(umap))), drop = FALSE]
if (ncol(umap) < 2) {
    stop("UMAP embedding has fewer than two dimensions.", call. = FALSE)
}
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")

sample_prefix <- sub("_.*$", "", cells)
sample_map <- c(P1 = "CRPC1", P2 = "CRPC2", P3 = "CRPC3")
sample <- unname(sample_map[sample_prefix])
sample[is.na(sample)] <- as.character(epi$orig.ident)[is.na(sample)]

meta <- data.frame(
    cell = cells,
    sample_prefix = sample_prefix,
    sample = sample,
    orig.ident = as.character(epi$orig.ident),
    seurat_cluster = as.character(epi$seurat_clusters),
    annotation = as.character(epi$annotation),
    UMAP_1 = umap[, 1],
    UMAP_2 = umap[, 2],
    stringsAsFactors = FALSE,
    check.names = FALSE
)

optional_cols <- intersect(
    c(
        "Phase", "S.Score", "G2M.Score",
        "nFeature_RNA", "nCount_RNA", "percent.mt",
        "scDblFinder.score", "scDblFinder.class"
    ),
    colnames(epi@meta.data)
)
for (col in optional_cols) {
    meta[[col]] <- epi@meta.data[cells, col, drop = TRUE]
}

metadata_csv <- file.path(OUT_DIR, "seurat_velocity_metadata.csv")
umap_csv <- file.path(OUT_DIR, "seurat_umap.csv")
cells_txt <- file.path(OUT_DIR, "epithelial_cells.txt")
summary_csv <- file.path(OUT_DIR, "seurat_annotation_summary.csv")

write.csv(meta, metadata_csv, row.names = FALSE, quote = TRUE)
write.csv(meta[, c("cell", "UMAP_1", "UMAP_2")], umap_csv, row.names = FALSE)
writeLines(cells, cells_txt)

summary_df <- as.data.frame(table(
    sample = meta$sample,
    annotation = meta$annotation
))
summary_df <- summary_df[summary_df$Freq > 0, ]
write.csv(summary_df, summary_csv, row.names = FALSE)

message("Wrote: ", metadata_csv)
message("Wrote: ", umap_csv)
message("Wrote: ", cells_txt)
message("Wrote: ", summary_csv)
