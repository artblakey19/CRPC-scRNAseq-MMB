# scVelo step 0 — export Seurat metadata + UMAP coords to CSV for Python side.
#
# scVelo (Python) 가 epi_filtered_clustered.rds 를 직접 못 읽기 때문에, cell name
# 별 metadata + UMAP coordinate 를 CSV 로 내보낸다. cell name 은 Seurat 의
# "P{idx}_{barcode}-1" 형식 (Integrated_Analysis.R add.cell.ids 와 동일).
#
# Output: Results/05_Epithelial_Downstream/scVelo/epi_metadata_for_scvelo.csv

library(Seurat)
library(dplyr)

IN_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/scVelo"
OUT_CSV <- file.path(OUT_DIR, "epi_metadata_for_scvelo.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

epi <- readRDS(IN_RDS)

umap <- as.data.frame(Embeddings(epi, "umap"))
colnames(umap) <- c("UMAP_1", "UMAP_2")

meta <- epi@meta.data %>%
    dplyr::select(any_of(c(
        "orig.ident", "seurat_clusters", "Phase",
        "nFeature_RNA", "nCount_RNA", "percent.mt",
        "S.Score", "G2M.Score"
    )))

# Attach barcode-keyed annotation (annotation_per_cell.csv, written by
# Epithelial_Annotation.R) so the Python side never re-derives labels from a
# hardcoded cluster->label map — renumber-safe across reseeded reruns.
ANNOT_CSV <- "Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv"
if (!file.exists(ANNOT_CSV)) {
    stop("annotation_per_cell.csv not found — run Epithelial_Annotation.R first: ",
         ANNOT_CSV)
}
annot <- read.csv(ANNOT_CSV, stringsAsFactors = FALSE)

out <- cbind(cell = rownames(meta), meta, umap)
out$annotation <- annot$annotation[match(out$cell, annot$cell)]
if (any(is.na(out$annotation))) {
    stop(sum(is.na(out$annotation)), " cells have no annotation match — ",
         "annotation_per_cell.csv is out of sync with the Stage 4 RDS.")
}
write.csv(out, OUT_CSV, row.names = FALSE)

message("Exported ", nrow(out), " cells to ", OUT_CSV)
