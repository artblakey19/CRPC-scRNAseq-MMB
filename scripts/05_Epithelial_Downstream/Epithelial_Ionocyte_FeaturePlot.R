# Epithelial_Ionocyte_FeaturePlot.R
# FeaturePlot for candidate ionocyte markers on filtered epithelial UMAP
# to confirm OE 5 (cluster 10) identity.
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds
# Writes: Results/05_Epithelial_Downstream/Annotation/FeaturePlot_ionocyte_markers.png
#         Results/05_Epithelial_Downstream/Annotation/VlnPlot_ionocyte_markers.png

suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

markers <- c("FOXI1", "ATP6V1G3", "ATP6V1B1", "KRT7", "KIT")

epi <- readRDS(IN_RDS)
# SCT assay has normalized 'data' layer; RNA in this object holds counts only.
DefaultAssay(epi) <- "SCT"

present <- markers %in% rownames(epi)
if (any(!present)) {
    message("Genes not in object: ", paste(markers[!present], collapse = ", "))
}
markers <- markers[present]

fp_list <- lapply(markers, function(g) {
    FeaturePlot(epi, features = g, order = TRUE, pt.size = 0.25) +
        ggtitle(g) +
        theme(plot.title = element_text(face = "bold"))
})
fp <- wrap_plots(fp_list, ncol = 3)
ggsave(file.path(OUT_DIR, "FeaturePlot_ionocyte_markers.png"),
       plot = fp, width = 15, height = 10, bg = "white", dpi = 200)

vp <- VlnPlot(epi, features = markers, group.by = "annotation",
              pt.size = 0, stack = TRUE, flip = TRUE) +
    NoLegend() +
    ggtitle("Ionocyte candidate markers by annotation")
ggsave(file.path(OUT_DIR, "VlnPlot_ionocyte_markers.png"),
       plot = vp, width = 9, height = 7, bg = "white", dpi = 200)

message("Wrote: ", file.path(OUT_DIR, "FeaturePlot_ionocyte_markers.png"))
message("Wrote: ", file.path(OUT_DIR, "VlnPlot_ionocyte_markers.png"))
