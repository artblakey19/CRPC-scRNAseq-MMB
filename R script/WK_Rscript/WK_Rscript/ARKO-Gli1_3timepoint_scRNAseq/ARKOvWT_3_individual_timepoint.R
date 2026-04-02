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
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/E18.5_only_Analysis/E18only_Rworkspace.RData")

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


#### P11_CtrlvARKO.combined Fib/SM####
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P11_only_Analysis/P11only_Rworkspace.RData")

Idents(object = P11_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P11_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3)

#subclustering Epi/Stro
Combined_P11_Epi <- subset(P11_CtrlvARKO.combined, idents = c("Basal", "Luminal", "Intermediate", "SV", "Ductus Deferens", "Neuroendocrine"))
Combined_P11_Stro <- subset(P11_CtrlvARKO.combined, idents = c("Fibroblast", "Prolif Stro", "SM", "Glia", "Neur", "Pericytes", "Endothelial", "Leukocytes", "Adipocytes"))

#Clustering
Idents(object = Combined_P11_Epi) <- "seurat_clusters"
DefaultAssay(Combined_P11_Epi) <- "integrated"
Combined_P11_Epi <- ScaleData(Combined_P11_Epi, verbose = FALSE)
Combined_P11_Epi <- RunPCA(Combined_P11_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P11_Epi, ndims = 50)
Combined_P11_Epi <- FindNeighbors(Combined_P11_Epi, reduction = "pca", dims = 1:20)
Combined_P11_Epi <- FindClusters(Combined_P11_Epi, resolution = 0.5)
Combined_P11_Epi <- RunUMAP(Combined_P11_Epi, reduction = "pca", dims = 1:20)
DimPlot(Combined_P11_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_P11_Stro) <- "seurat_clusters"
DefaultAssay(Combined_P11_Stro) <- "integrated"
Combined_P11_Stro <- ScaleData(Combined_P11_Stro, verbose = FALSE)
Combined_P11_Stro <- RunPCA(Combined_P11_Stro, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P11_Epi, ndims = 50)
Combined_P11_Stro <- FindNeighbors(Combined_P11_Stro, reduction = "pca", dims = 1:20)
Combined_P11_Stro <- FindClusters(Combined_P11_Stro, resolution = 0.5)
Combined_P11_Stro <- RunUMAP(Combined_P11_Stro, reduction = "pca", dims = 1:20)
DimPlot(Combined_P11_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

#subclustering Fib/SM
Idents(object = Combined_P11_Stro) <- "CellTypes"
DimPlot(Combined_P11_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

Combined_P11_FibSM <- subset(Combined_P11_Stro, idents = c("Fibroblast", "SM"))

#clustering
Idents(object = Combined_P11_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_P11_FibSM) <- "integrated"
Combined_P11_FibSM <- ScaleData(Combined_P11_FibSM, verbose = FALSE)
Combined_P11_FibSM <- RunPCA(Combined_P11_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P11_FibSM, ndims = 50)
Combined_P11_FibSM <- FindNeighbors(Combined_P11_FibSM, reduction = "pca", dims = 1:20)
Combined_P11_FibSM <- FindClusters(Combined_P11_FibSM, resolution = 0.5)
Combined_P11_FibSM <- RunUMAP(Combined_P11_FibSM, reduction = "pca", dims = 1:20)
DimPlot(Combined_P11_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_P11_FibSM) <- "CellTypes"
DimPlot(Combined_P11_FibSM, reduction = "umap", pt.size = 0.3)
Idents(object = Combined_P11_FibSM) <- "seurat_clusters"
DimPlot(Combined_P11_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE)

#Cell count (split)
Idents(object = Combined_P11_FibSM) <- "stim"
Combined_P11_FibSM$stim.seurat_clusters <- paste(Idents(Combined_P11_FibSM), Combined_P11_FibSM$seurat_clusters, sep = "_")
Idents(object = Combined_P11_FibSM) <- "stim.seurat_clusters"
table(Idents(Combined_P11_FibSM))

#Featureplot
DefaultAssay(Combined_P11_FibSM) <- "RNA"
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Heatmap
Idents(object = Combined_P11_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_P11_FibSM) <- "RNA"
Combined_P11_FibSM <- ScaleData(Combined_P11_FibSM, features = rownames(Combined_P11_FibSM))
Combined_P11_FibSM.markers <- FindAllMarkers(Combined_P11_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined_P11_FibSMTop10 <- Combined_P11_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Combined_P11_FibSM, features = c(Combined_P11_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
Combined_P11_FibSM.0.1markers <- FindAllMarkers(Combined_P11_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(Combined_P11_FibSM.0.1markers, "Combined_P11_FibSM.0.1markers-1.csv")

#FeaturePlots
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Igfbp3"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Sfrp2"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Inhbb"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Dcn"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")


FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Cxcl1"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Ctgf"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Igf1"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Ptn"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Ccl2"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Cxcl2"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Csf1"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P11_FibSM, reduction = "umap", features = c("Rspo3"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")



#### P42_CtrlvARKO.combined Fib/SM####
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P42_only_Analysis/P42only_Rworkspace.RData")

Idents(object = P42_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P42_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3)

#subclustering Epi/Stro
Combined_P42_Epi <- subset(P42_CtrlvARKO.combined, idents = c("Basal", "Luminal", "SV", "Prolif Epi"))
Combined_P42_Stro <- subset(P42_CtrlvARKO.combined, idents = c("Fibroblast", "SM", "Lymphocytes", "Leukocytes", "Endothelial", "Pericyte", "Glia"))

#Clustering
Idents(object = Combined_P42_Epi) <- "seurat_clusters"
DefaultAssay(Combined_P42_Epi) <- "integrated"
Combined_P42_Epi <- ScaleData(Combined_P42_Epi, verbose = FALSE)
Combined_P42_Epi <- RunPCA(Combined_P42_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P42_Epi, ndims = 50)
Combined_P42_Epi <- FindNeighbors(Combined_P42_Epi, reduction = "pca", dims = 1:15)
Combined_P42_Epi <- FindClusters(Combined_P42_Epi, resolution = 0.5)
Combined_P42_Epi <- RunUMAP(Combined_P42_Epi, reduction = "pca", dims = 1:15)
DimPlot(Combined_P42_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_P42_Stro) <- "seurat_clusters"
DefaultAssay(Combined_P42_Stro) <- "integrated"
Combined_P42_Stro <- ScaleData(Combined_P42_Stro, verbose = FALSE)
Combined_P42_Stro <- RunPCA(Combined_P42_Stro, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P42_Epi, ndims = 50)
Combined_P42_Stro <- FindNeighbors(Combined_P42_Stro, reduction = "pca", dims = 1:24)
Combined_P42_Stro <- FindClusters(Combined_P42_Stro, resolution = 0.5)
Combined_P42_Stro <- RunUMAP(Combined_P42_Stro, reduction = "pca", dims = 1:24)
DimPlot(Combined_P42_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

#subclustering Fib/SM
Idents(object = Combined_P42_Stro) <- "CellTypes"
DimPlot(Combined_P42_Stro, reduction = "umap", pt.size = 0.3, label = TRUE)

Combined_P42_FibSM <- subset(Combined_P42_Stro, idents = c("Fibroblast", "SM"))

#clustering
Idents(object = Combined_P42_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_P42_FibSM) <- "integrated"
Combined_P42_FibSM <- ScaleData(Combined_P42_FibSM, verbose = FALSE)
Combined_P42_FibSM <- RunPCA(Combined_P42_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P42_FibSM, ndims = 50)
Combined_P42_FibSM <- FindNeighbors(Combined_P42_FibSM, reduction = "pca", dims = 1:20)
Combined_P42_FibSM <- FindClusters(Combined_P42_FibSM, resolution = 0.5)
Combined_P42_FibSM <- RunUMAP(Combined_P42_FibSM, reduction = "pca", dims = 1:20)
DimPlot(Combined_P42_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

#Remove cluster8 (Pbsn-positive)
Idents(object = Combined_P42_FibSM) <- "seurat_clusters"
Combined_P42_FibSM <- subset(Combined_P42_FibSM, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "9", "10"))

Idents(object = Combined_P42_FibSM ) <- "seurat_clusters"
DefaultAssay(Combined_P42_FibSM ) <- "integrated"
Combined_P42_FibSM  <- ScaleData(Combined_P42_FibSM , verbose = FALSE)
Combined_P42_FibSM  <- RunPCA(Combined_P42_FibSM , npcs = 50, verbose = FALSE)
ElbowPlot(Combined_P42_FibSM , ndims = 50)
Combined_P42_FibSM  <- FindNeighbors(Combined_P42_FibSM , reduction = "pca", dims = 1:20)
Combined_P42_FibSM  <- FindClusters(Combined_P42_FibSM , resolution = 0.5)
Combined_P42_FibSM  <- RunUMAP(Combined_P42_FibSM , reduction = "pca", dims = 1:20)
DimPlot(Combined_P42_FibSM , reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Combined_P42_FibSM) <- "CellTypes"
DimPlot(Combined_P42_FibSM, reduction = "umap", pt.size = 0.3)
Idents(object = Combined_P42_FibSM) <- "seurat_clusters"
DimPlot(Combined_P42_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE)

#Cell count (split)
Idents(object = Combined_P42_FibSM) <- "stim"
Combined_P42_FibSM$stim.seurat_clusters <- paste(Idents(Combined_P42_FibSM), Combined_P42_FibSM$seurat_clusters, sep = "_")
Idents(object = Combined_P42_FibSM) <- "stim.seurat_clusters"
table(Idents(Combined_P42_FibSM))

#Featureplot
DefaultAssay(Combined_P42_FibSM) <- "RNA"
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Heatmap
Idents(object = Combined_P42_FibSM) <- "seurat_clusters"
DefaultAssay(Combined_P42_FibSM) <- "RNA"
Combined_P42_FibSM <- ScaleData(Combined_P42_FibSM, features = rownames(Combined_P42_FibSM))
Combined_P42_FibSM.markers <- FindAllMarkers(Combined_P42_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Combined_P42_FibSMTop10 <- Combined_P42_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(Combined_P42_FibSM, features = c(Combined_P42_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
Combined_P42_FibSM.0.1markers <- FindAllMarkers(Combined_P42_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(Combined_P42_FibSM.0.1markers, "Combined_P42_FibSM.0.1markers.csv")

#FeaturePlots
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Igfbp3"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Sfrp2"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Inhbb"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Dcn"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")


FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Cxcl1"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Ctgf"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Igf1"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Ptn"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Ccl2"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Cxcl2"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Csf1"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_P42_FibSM, reduction = "umap", features = c("Rspo3"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
