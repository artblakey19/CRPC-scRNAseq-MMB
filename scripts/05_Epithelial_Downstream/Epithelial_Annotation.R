# Epithelial_Annotation_decontX.R
# decontX filtered epi cluster -> cell-type label mapping.
# 원본 Epithelial_Annotation.R 와 동일 라벨 체계(12 labels). 매핑은 decontX
# res-0.3 canonical 클러스터를 원본 annotated 객체와 marker-profile 상관
# (Spearman 0.94–0.99, 1:1)으로 도출하고 사용자가 확정(2026-07-11).
#
# Reads:  Results/04_Epithelial_Filtered/epi_filtered_clustered.rds
# Writes: Results/05_Epithelial_Downstream/epi_annotated.rds
#         Results/05_Epithelial_Downstream/Annotation/*

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})
source("scripts/00_utils/scRNA_utils.R")

IN_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
OUT_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# decontX res-0.3 (12 clusters) 전용 매핑. cluster 정수 ID 는 rerun 간 불안정하므로
# 이 맵은 현재 Results/04_..._decontX RDS 에만 유효.
#   BE 5 = cl3 (AP-1/heat-stress), OE = cl6 (CFTR+), LE(ARPC) = cl9 (KLK2/FOLH1),
#   Ionocyte = cl11 (FOXI1/ATP6V), Club = cl2 (MMP7/PIGR).
cluster_to_label <- c(
    "0"  = "BE 1",
    "1"  = "BE 4",
    "2"  = "Club",
    "3"  = "BE 5",
    "4"  = "Hillock 1",
    "5"  = "BE 2",
    "6"  = "OE",
    "7"  = "Hillock 2",
    "8"  = "BE 6",
    "9"  = "LE(ARPC)",
    "10" = "BE 3",
    "11" = "Ionocyte"
)

label_levels <- c(
    "LE(ARPC)", "Club", "Hillock 1", "Hillock 2",
    "BE 1", "BE 2", "BE 3",
    "BE 4", "BE 5", "BE 6",
    "OE",
    "Ionocyte"
)

epi <- readRDS(IN_RDS)
clu <- as.character(epi$seurat_clusters)

missing <- setdiff(unique(clu), names(cluster_to_label))
if (length(missing) > 0) {
    stop("Unmapped clusters: ", paste(missing, collapse = ", "))
}

epi$annotation <- factor(unname(cluster_to_label[clu]), levels = label_levels)
Idents(epi) <- "annotation"

n_lab <- nlevels(epi$annotation)
pal <- utils_cb_palette(n_lab)

p1 <- DimPlot(epi, group.by = "annotation", label = TRUE, repel = TRUE,
              pt.size = 0.3, cols = pal) +
    ggtitle("Filtered epithelial — annotation (decontX)")
ggsave(file.path(OUT_DIR, "UMAP_annotation.png"),
       plot = p1, width = 10, height = 8, bg = "white")

p2 <- DimPlot(epi, group.by = "annotation", split.by = "orig.ident",
              label = TRUE, repel = TRUE, pt.size = 0.3, cols = pal)
ggsave(file.path(OUT_DIR, "UMAP_annotation_split.png"),
       plot = p2, width = 24, height = 8, bg = "white")

tab <- as.data.frame(table(cluster = epi$seurat_clusters,
                           annotation = epi$annotation))
tab <- tab[tab$Freq > 0, ]
write.csv(tab, file.path(OUT_DIR, "cluster_to_label.csv"), row.names = FALSE)

write.csv(
    data.frame(
        cell = colnames(epi),
        seurat_cluster = as.character(epi$seurat_clusters),
        annotation = as.character(epi$annotation)
    ),
    file.path(OUT_DIR, "annotation_per_cell.csv"),
    row.names = FALSE
)

writeLines(label_levels, file.path(OUT_DIR, "label_levels.txt"))

saveRDS(epi, OUT_RDS)
message("Annotated object (decontX): ", OUT_RDS)
message("Figures + mapping: ", OUT_DIR)
