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

#### E18_CtrlvARKO.combined Fib/SM####
Idents(object = E18_CtrlvARKO.combined) <- "CellTypes"
DimPlot(E18_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3)

#subclustering Epi/Stro
Combined_E18_Epi <- subset(E18_CtrlvARKO.combined, idents = c("UGE", "BE", "WD"))
Combined_E18_Stro <- subset(E18_CtrlvARKO.combined, idents = c("Fibroblast", "Prolif Stro", "SM", "Glia", "Neur", "Pericyte", "Endothelial", "Myoblast", "Leukocyte"))

#Clustering
Idents(object = Combined_E18_Epi) <- "seurat_clusters"
DefaultAssay(Combined_E18_Epi) <- "integrated"
Combined_E18_Epi <- ScaleData(Combined_E18_Epi, verbose = FALSE)
Combined_E18_Epi <- RunPCA(Combined_E18_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_E18_Epi, ndims = 50)
Combined_E18_Epi <- FindNeighbors(Combined_E18_Epi, reduction = "pca", dims = 1:20)
Combined_E18_Epi <- FindClusters(Combined_E18_Epi, resolution = 0.5)
Combined_E18_Epi <- RunUMAP(Combined_E18_Epi, reduction = "pca", dims = 1:20)
DimPlot(Combined_E18_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_E18_Stro) <- "seurat_clusters"
DefaultAssay(Combined_E18_Stro) <- "integrated"
Combined_E18_Stro <- ScaleData(Combined_E18_Stro, verbose = FALSE)
Combined_E18_Stro <- RunPCA(Combined_E18_Stro, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_E18_Epi, ndims = 50)
Combined_E18_Stro <- FindNeighbors(Combined_E18_Stro, reduction = "pca", dims = 1:20)
Combined_E18_Stro <- FindClusters(Combined_E18_Stro, resolution = 0.5)
Combined_E18_Stro <- RunUMAP(Combined_E18_Stro, reduction = "pca", dims = 1:20)
DimPlot(Combined_E18_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(Combined_E18_Stro) <- "RNA"
FeaturePlot(Combined_E18_Stro, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
Idents(object = Combined_E18_Stro) <- "CellTypes"
DimPlot(Combined_E18_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

Combined_E18_FibSM <- subset(Combined_E18_Stro, idents = c("Fibroblast", "SM"))

#clustering
Idents(object = Combined_E18_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_E18_FibSM) <- "integrated"
Combined_E18_FibSM <- ScaleData(Combined_E18_FibSM, verbose = FALSE)
Combined_E18_FibSM <- RunPCA(Combined_E18_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_E18_FibSM, ndims = 50)
Combined_E18_FibSM <- FindNeighbors(Combined_E18_FibSM, reduction = "pca", dims = 1:20)
Combined_E18_FibSM <- FindClusters(Combined_E18_FibSM, resolution = 0.5)
Combined_E18_FibSM <- RunUMAP(Combined_E18_FibSM, reduction = "pca", dims = 1:20)
DimPlot(Combined_E18_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_E18_FibSM) <- "CellTypes"
DimPlot(Combined_E18_FibSM, reduction = "umap", pt.size = 0.3)
Idents(object = Combined_E18_FibSM) <- "seurat_clusters"
DimPlot(Combined_E18_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE)

#Cell count (split)
Idents(object = Combined_E18_FibSM) <- "stim"
Combined_E18_FibSM$stim.seurat_clusters <- paste(Idents(Combined_E18_FibSM), Combined_E18_FibSM$seurat_clusters, sep = "_")
Idents(object = Combined_E18_FibSM) <- "stim.seurat_clusters"
table(Idents(Combined_E18_FibSM))

#Featureplot
DefaultAssay(Combined_E18_FibSM) <- "RNA"
FeaturePlot(Combined_E18_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

tiff(file = "Combined_E18_FibSM Rspo3 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_E18_FibSM, reduction = "umap", features = c("Rspo3"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_E18_FibSM Ptn Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_E18_FibSM, reduction = "umap", features = c("Ptn"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_E18_FibSM Wnt5a Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_E18_FibSM, reduction = "umap", features = c("Wnt5a"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_E18_FibSM Inhbb Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_E18_FibSM, reduction = "umap", features = c("Inhbb"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = Combined_E18_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_E18_FibSM) <- "RNA"
Combined_E18_FibSM <- ScaleData(Combined_E18_FibSM, features = rownames(Combined_E18_FibSM))
Combined_E18_FibSM.markers <- FindAllMarkers(Combined_E18_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined_E18_FibSMTop10 <- Combined_E18_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Combined_E18_FibSM, features = c(Combined_E18_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
Combined_E18_FibSM.0.1markers <- FindAllMarkers(Combined_E18_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(Combined_E18_FibSM.0.1markers, "Combined_E18_FibSM.0.1markers-1.csv")

#### Combined_P11_FibSM ####

Idents(object = P11_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P11_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3)

#subclustering Epi/Stro
Combined_P11_Epi <- subset(P11_CtrlvARKO.combined, idents = c("UGE", "BE", "WD"))
Combined_P11_Stro <- subset(P11_CtrlvARKO.combined, idents = c("Fibroblast", "Prolif Stro", "SM", "Glia", "Neur", "Pericyte", "Endothelial", "Myoblast", "Leukocyte"))

#Clustering
Idents(object = Combined_E18_Epi) <- "seurat_clusters"
DefaultAssay(Combined_E18_Epi) <- "integrated"
Combined_E18_Epi <- ScaleData(Combined_E18_Epi, verbose = FALSE)
Combined_E18_Epi <- RunPCA(Combined_E18_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_E18_Epi, ndims = 50)
Combined_E18_Epi <- FindNeighbors(Combined_E18_Epi, reduction = "pca", dims = 1:20)
Combined_E18_Epi <- FindClusters(Combined_E18_Epi, resolution = 0.5)
Combined_E18_Epi <- RunUMAP(Combined_E18_Epi, reduction = "pca", dims = 1:20)
DimPlot(Combined_E18_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_E18_Stro) <- "seurat_clusters"
DefaultAssay(Combined_E18_Stro) <- "integrated"
Combined_E18_Stro <- ScaleData(Combined_E18_Stro, verbose = FALSE)
Combined_E18_Stro <- RunPCA(Combined_E18_Stro, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_E18_Epi, ndims = 50)
Combined_E18_Stro <- FindNeighbors(Combined_E18_Stro, reduction = "pca", dims = 1:20)
Combined_E18_Stro <- FindClusters(Combined_E18_Stro, resolution = 0.5)
Combined_E18_Stro <- RunUMAP(Combined_E18_Stro, reduction = "pca", dims = 1:20)
DimPlot(Combined_E18_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(Combined_E18_Stro) <- "RNA"
FeaturePlot(Combined_E18_Stro, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
Idents(object = Combined_E18_Stro) <- "CellTypes"
DimPlot(Combined_E18_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

Combined_E18_FibSM <- subset(Combined_E18_Stro, idents = c("Fibroblast", "SM"))

#clustering
Idents(object = Combined_E18_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_E18_FibSM) <- "integrated"
Combined_E18_FibSM <- ScaleData(Combined_E18_FibSM, verbose = FALSE)
Combined_E18_FibSM <- RunPCA(Combined_E18_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_E18_FibSM, ndims = 50)
Combined_E18_FibSM <- FindNeighbors(Combined_E18_FibSM, reduction = "pca", dims = 1:20)
Combined_E18_FibSM <- FindClusters(Combined_E18_FibSM, resolution = 0.5)
Combined_E18_FibSM <- RunUMAP(Combined_E18_FibSM, reduction = "pca", dims = 1:20)
DimPlot(Combined_E18_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_E18_FibSM) <- "CellTypes"
DimPlot(Combined_E18_FibSM, reduction = "umap", pt.size = 0.3)
Idents(object = Combined_E18_FibSM) <- "seurat_clusters"
DimPlot(Combined_E18_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE)

#Cell count (split)
Idents(object = Combined_E18_FibSM) <- "stim"
Combined_E18_FibSM$stim.seurat_clusters <- paste(Idents(Combined_E18_FibSM), Combined_E18_FibSM$seurat_clusters, sep = "_")
Idents(object = Combined_E18_FibSM) <- "stim.seurat_clusters"
table(Idents(Combined_E18_FibSM))

#Featureplot
DefaultAssay(Combined_E18_FibSM) <- "RNA"
FeaturePlot(Combined_E18_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Heatmap
Idents(object = Combined_E18_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_E18_FibSM) <- "RNA"
Combined_E18_FibSM <- ScaleData(Combined_E18_FibSM, features = rownames(Combined_E18_FibSM))
Combined_E18_FibSM.markers <- FindAllMarkers(Combined_E18_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined_E18_FibSMTop10 <- Combined_E18_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Combined_E18_FibSM, features = c(Combined_E18_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
Combined_E18_FibSM.0.1markers <- FindAllMarkers(Combined_E18_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(Combined_E18_FibSM.0.1markers, "Combined_E18_FibSM.0.1markers-1.csv")