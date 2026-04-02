#### Pseudotime Trajectory ####

#### Add necessary tools to library ####

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)
setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK")

#### Combined_ARKO_FibSM ####
Idents(object = Combined_ARKO_FibSM) <- "seurat_clusters"
DimPlot(Combined_ARKO_FibSM, reduction = "umap", pt.size = 0.3)

tiff(file = "Combined_ARKO_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined_ARKO_FibSM, reduction = "umap", pt.size = 0.3)
dev.off()

DimPlot(Combined_ARKO_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim")

#Cell count (split)
Idents(object = Combined_ARKO_FibSM) <- "stim"
Combined_ARKO_FibSM$stim.seurat_clusters <- paste(Idents(Combined_ARKO_FibSM), Combined_ARKO_FibSM$seurat_clusters, sep = "_")
Idents(object = Combined_ARKO_FibSM) <- "stim.seurat_clusters"
table(Idents(Combined_ARKO_FibSM))

FeaturePlot(Combined_ARKO_FibSM, reduction = "umap", features = c("EGFP"), max.cutoff = "q90", cols = c("light grey", "purple"), pt.size = 0.3, split.by = "stim")
FeaturePlot(Combined_ARKO_FibSM, reduction = "umap", features = c("Ar"), max.cutoff = "q90", cols = c("light grey", "purple"), pt.size = 0.3, split.by = "stim")
FeaturePlot(Combined_ARKO_FibSM, reduction = "umap", features = c("Esr1"), max.cutoff = "q90", cols = c("light grey", "purple"), pt.size = 0.3, split.by = "stim")
FeaturePlot(Combined_ARKO_FibSM, reduction = "umap", features = c("Inhbb"), max.cutoff = "q90", cols = c("light grey", "purple"), pt.size = 0.3, split.by = "stim")
FeaturePlot(Combined_ARKO_FibSM, reduction = "umap", features = c("Wnt7b"), max.cutoff = "q90", cols = c("light grey", "purple"), pt.size = 0.3, split.by = "stim")
FeaturePlot(Combined_ARKO_FibSM, reduction = "umap", features = c("Igfbp3"), max.cutoff = "q90", cols = c("light grey", "purple"), pt.size = 0.3, split.by = "stim")


#DEGs
DefaultAssay(Combined_ARKO_FibSM) <- "RNA"
Idents(object = Combined_ARKO_FibSM) <- "seurat_clusters"
Combined_ARKO_FibSM <- ScaleData(Combined_ARKO_FibSM, features = rownames(Combined_ARKO_FibSM))
Combined_ARKO_FibSM.markers <- FindAllMarkers(Combined_ARKO_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Combined_ARKO_FibSM.markers, "Combined_ARKO_FibSM.markers.csv")

#Pseudotime plots
plot_cell_trajectory(ARKO_FibSM_Pseudo, color_by = "seurat_clusters", show_branch_points = FALSE)
plot_cell_trajectory(ARKO_FibSM_Pseudo, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
plot_cell_trajectory(ARKO_FibSM_Pseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
plot_cell_trajectory(ARKO_FibSM_Pseudo, markers = "Esr1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
plot_cell_trajectory(ARKO_FibSM_Pseudo, markers = "Inhbb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
plot_cell_trajectory(ARKO_FibSM_Pseudo, markers = "Igfbp3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))