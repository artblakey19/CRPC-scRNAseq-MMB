# Epithelial_Annotation.R
# Filtered epi cluster -> cell-type label mapping (user-curated).
# Reads:  Results/04_Epithelial_Filtered/epi_filtered_clustered.rds
# Writes: Results/05_Epithelial_Downstream/epi_annotated.rds
#         Results/05_Epithelial_Downstream/Annotation/UMAP_annotation.png
#         Results/05_Epithelial_Downstream/Annotation/UMAP_annotation_split.png
#         Results/05_Epithelial_Downstream/Annotation/cluster_to_label.csv
#
# Source filtered RDS is left untouched.

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

cluster_to_label <- c(
    "0"  = "BE 1",
    "1"  = "OE 1",
    "2"  = "OE 2",
    "3"  = "BE 2",
    "4"  = "Club-like",
    "5"  = "Hillock-like 1",
    "6"  = "OE 3",
    "7"  = "Hillock-like 2",
    "8"  = "OE 4",
    "9"  = "ARPC",
    "10" = "Ionocyte-like"
)

# Display order: lineage groups, OE 1..4 numeric, then Ionocyte-like as its own
# rare-cell category (FOXI1+ ATP6V+ — see Epithelial_Ionocyte_FeaturePlot.R).
label_levels <- c(
    "ARPC", "Club-like", "Hillock-like 1", "Hillock-like 2",
    "BE 1", "BE 2",
    "OE 1", "OE 2", "OE 3", "OE 4",
    "Ionocyte-like"
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
    ggtitle("Filtered epithelial — annotation")
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

# Per-cell mapping table used by downstream plotting scripts to avoid
# re-loading the full annotated RDS just to recover labels.
write.csv(
    data.frame(
        cell = colnames(epi),
        seurat_cluster = as.character(epi$seurat_clusters),
        annotation = as.character(epi$annotation)
    ),
    file.path(OUT_DIR, "annotation_per_cell.csv"),
    row.names = FALSE
)

# CSV loses factor level order; persist the display ordering separately so
# downstream plotting scripts can restore it without hardcoding.
writeLines(label_levels, file.path(OUT_DIR, "label_levels.txt"))

saveRDS(epi, OUT_RDS)
message("Annotated object: ", OUT_RDS)
message("Figures + mapping: ", OUT_DIR)
