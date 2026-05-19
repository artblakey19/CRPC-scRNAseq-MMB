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
        "S.Score", "G2M.Score", "copykat_prediction"   # legacy carry-over
    )))

out <- cbind(cell = rownames(meta), meta, umap)
write.csv(out, OUT_CSV, row.names = FALSE)

message("Exported ", nrow(out), " cells to ", OUT_CSV)
