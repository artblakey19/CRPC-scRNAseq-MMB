# Epithelial_ClusterMarkers.R
# Per-cell-type marker genes for the filtered epithelial object.
#
# Stage 04's Results/04_Epithelial_Filtered/all_markers.csv is keyed by the raw
# integer seurat_clusters (0-11). This script re-derives markers keyed by the
# CURRENT annotation labels (BE 1-6, LE(ARPC), Club, Hillock 1/2, OE, Ionocyte)
# so downstream figures/tables read in cell-type terms. Clustering is unchanged
# — annotation is a label-only renaming — so results are 1:1 with the integer
# markers, just relabelled and re-ordered for display.
#
# Markers are computed on the SCT assay (PrepSCTFindMarkers), matching the
# convention in utils_save_all_markers() used at Stage 04.
#
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds
# Writes: Results/05_Epithelial_Downstream/ClusterMarkers/all_markers_by_celltype.csv
#         Results/05_Epithelial_Downstream/ClusterMarkers/top_markers_by_celltype.csv
#         Results/05_Epithelial_Downstream/ClusterMarkers/top_markers_DotPlot.png
#         Results/05_Epithelial_Downstream/ClusterMarkers/top_markers_Heatmap.png

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})
source("scripts/00_utils/scRNA_utils.R")

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/ClusterMarkers"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

TOP_N_TABLE <- 25   # top markers per cell type written to the top-markers CSV
TOP_N_PLOT  <- 5    # top markers per cell type shown in the dotplot/heatmap

epi <- readRDS(IN_RDS)

# Display order defined by the annotation script; fall back to sorted levels.
lvl_file <- "Results/05_Epithelial_Downstream/Annotation/label_levels.txt"
if (file.exists(lvl_file)) {
    label_levels <- readLines(lvl_file)
    epi$annotation <- factor(as.character(epi$annotation), levels = label_levels)
}
Idents(epi) <- "annotation"

# SCT assay + PrepSCTFindMarkers, matching utils_save_all_markers() convention.
DefaultAssay(epi) <- "SCT"
epi <- PrepSCTFindMarkers(epi)

all_markers <- FindAllMarkers(
    epi,
    only.pos        = TRUE,
    min.pct         = 0.25,
    logfc.threshold = 0.5
)

# Keep FindAllMarkers' cluster ordering (= annotation factor levels).
write.csv(all_markers, file.path(OUT_DIR, "all_markers_by_celltype.csv"),
          row.names = FALSE)

top_tbl <- all_markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = TOP_N_TABLE, with_ties = FALSE) %>%
    ungroup()
write.csv(top_tbl, file.path(OUT_DIR, "top_markers_by_celltype.csv"),
          row.names = FALSE)

# ---- Figures: top-N markers per cell type ----
top_plot <- all_markers %>%
    filter(p_val_adj < 0.05) %>%
    group_by(cluster) %>%
    slice_max(order_by = avg_log2FC, n = TOP_N_PLOT, with_ties = FALSE) %>%
    ungroup()
plot_features <- unique(top_plot$gene)

p_dot <- DotPlot(epi, features = plot_features) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = paste0("Top ", TOP_N_PLOT, " markers per epithelial cell type"),
         x = NULL, y = NULL)
ggsave(file.path(OUT_DIR, "top_markers_DotPlot.png"),
       plot = p_dot, width = max(16, length(plot_features) * 0.28), height = 8,
       bg = "white", limitsize = FALSE)

# Heatmap on scaled SCT data; downsample cells per group so the raster stays legible.
epi_hm <- subset(epi, downsample = 300)
epi_hm <- ScaleData(epi_hm, features = plot_features, verbose = FALSE)
n_grp  <- nlevels(droplevels(Idents(epi_hm)))
p_hm <- DoHeatmap(epi_hm, features = plot_features,
                  group.colors = utils_cb_palette(n_grp)) +
    scale_fill_gradient2(low = "#1F77B4", mid = "grey95", high = "#D7261E") +
    theme(axis.text.y = element_text(size = 6))
ggsave(file.path(OUT_DIR, "top_markers_Heatmap.png"),
       plot = p_hm, width = 14, height = max(10, length(plot_features) * 0.14),
       bg = "white", limitsize = FALSE)

message("Markers by cell type → ", OUT_DIR)
message("  all_markers_by_celltype.csv  (", nrow(all_markers), " rows)")
message("  top_markers_by_celltype.csv  (top ", TOP_N_TABLE, "/cell type)")
