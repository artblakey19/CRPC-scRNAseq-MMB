#### 6 samples Fib/SM ####

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
setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 FibSM")

load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/E18.5_only_Analysis/E18only_Rworkspace.RData")
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P42_only_Analysis/P42only_Rworkspace.RData")
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P11_only_Analysis/P11only_Rworkspace.RData")

#### E18.5 Ctrl####
Idents(object = E18_Ctrl) <- "seurat_clusters"
DefaultAssay(E18_Ctrl) <- "RNA"
E18_Ctrl <- ScaleData(E18_Ctrl, verbose = FALSE)
E18_Ctrl <- RunPCA(E18_Ctrl, npcs = 50, verbose = FALSE)
ElbowPlot(E18_Ctrl)

E18_Ctrl <- FindNeighbors(E18_Ctrl, reduction = "pca", dims = 1:20)
E18_Ctrl <- FindClusters(E18_Ctrl, resolution = 0.5)

E18_Ctrl <- RunUMAP(E18_Ctrl, reduction = "pca", dims = 1:20)
DimPlot(E18_Ctrl, reduction = "umap", pt.size = 0.3)

DefaultAssay(E18_Ctrl) <- "RNA"
FeaturePlot(E18_Ctrl, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_Ctrl, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_Ctrl, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
E18_Ctrl_FibSM <- subset(E18_Ctrl, idents = c("0", "1", "3", "4", "7", "6"))
Idents(object = E18_Ctrl_FibSM) <- "seurat_clusters"
E18_Ctrl_FibSM <- RenameIdents(object = E18_Ctrl_FibSM, '0' = "Fibroblast", '1' = "Fibroblast", '3' = "Fibroblast", '4' = "Fibroblast", '7' = "Fibroblast", '6' = "SM")
DimPlot(E18_Ctrl_FibSM, reduction = "umap", pt.size = 1)
E18_Ctrl_FibSM[["CellTypes"]] <- Idents(object = E18_Ctrl_FibSM)
Idents(object = E18_Ctrl_FibSM) <- "CellTypes"
DimPlot(E18_Ctrl_FibSM, reduction = "umap", pt.size = 1)

#### E18.5 ARKO####
Idents(object = E18_ARKO) <- "seurat_clusters"
E18_ARKO <- ScaleData(E18_ARKO, verbose = FALSE)
E18_ARKO <- RunPCA(E18_ARKO, npcs = 50, verbose = FALSE)
ElbowPlot(E18_ARKO)

E18_ARKO <- FindNeighbors(E18_ARKO, reduction = "pca", dims = 1:20)
E18_ARKO <- FindClusters(E18_ARKO, resolution = 0.5)
E18_ARKO <- RunUMAP(E18_ARKO, reduction = "pca", dims = 1:20)
DimPlot(E18_ARKO, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(E18_ARKO) <- "RNA"
FeaturePlot(E18_ARKO, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_ARKO, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_ARKO, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
E18_ARKO_FibSM <- subset(E18_ARKO, idents = c("0", "1", "4", "5", "6", "7", "8", "9"))
Idents(object = E18_ARKO_FibSM) <- "seurat_clusters"
E18_ARKO_FibSM <- RenameIdents(object = E18_ARKO_FibSM, '0' = "Fibroblast", '1' = "Fibroblast", '4' = "Fibroblast", '5' = "Fibroblast", '6' = "Fibroblast", '7' = "Fibroblast", '8' = "Fibroblast", '9' = "SM")
DimPlot(E18_ARKO_FibSM, reduction = "umap", pt.size = 1)
E18_ARKO_FibSM[["CellTypes"]] <- Idents(object = E18_ARKO_FibSM)
Idents(object = E18_ARKO_FibSM) <- "CellTypes"
DimPlot(E18_ARKO_FibSM, reduction = "umap", pt.size = 1)

#### P11 Ctrl####
Idents(object = P11_Ctrl) <- "seurat_clusters"
P11_Ctrl <- ScaleData(P11_Ctrl, verbose = FALSE)
P11_Ctrl <- RunPCA(P11_Ctrl, npcs = 50, verbose = FALSE)
ElbowPlot(P11_Ctrl)

P11_Ctrl <- FindNeighbors(P11_Ctrl, reduction = "pca", dims = 1:15)
P11_Ctrl <- FindClusters(P11_Ctrl, resolution = 0.5)
P11_Ctrl <- RunUMAP(P11_Ctrl, reduction = "pca", dims = 1:15)
DimPlot(P11_Ctrl, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(P11_Ctrl) <- "RNA"
FeaturePlot(P11_Ctrl, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_Ctrl, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_Ctrl, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
P11_Ctrl_FibSM <- subset(P11_Ctrl, idents = c("4", "3", "0", "11", "13", "2", "5"))
Idents(object = P11_Ctrl_FibSM) <- "seurat_clusters"
P11_Ctrl_FibSM <- RenameIdents(object = P11_Ctrl_FibSM, '4' = "Fibroblast", '3' = "Fibroblast", '0' = "Fibroblast", '11' = "Fibroblast", '13' = "Fibroblast", '2' = "SM", '5' = "SM")
P11_Ctrl_FibSM[["CellTypes"]] <- Idents(object = P11_Ctrl_FibSM)
Idents(object = P11_Ctrl_FibSM) <- "CellTypes"
DimPlot(P11_Ctrl_FibSM, reduction = "umap", pt.size = 1)

#### P11 ARKO####
Idents(object = P11_ARKO) <- "seurat_clusters"
P11_ARKO <- ScaleData(P11_ARKO, verbose = FALSE)
P11_ARKO <- RunPCA(P11_ARKO, npcs = 50, verbose = FALSE)
ElbowPlot(P11_ARKO)

P11_ARKO <- FindNeighbors(P11_ARKO, reduction = "pca", dims = 1:20)
P11_ARKO <- FindClusters(P11_ARKO, resolution = 0.5)
P11_ARKO <- RunUMAP(P11_ARKO, reduction = "pca", dims = 1:20)
DimPlot(P11_ARKO, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(P11_ARKO) <- "RNA"
FeaturePlot(P11_ARKO, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_ARKO, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_ARKO, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
P11_ARKO_FibSM <- subset(P11_ARKO, idents = c("1", "3", "5", "7", "0", "2", "6"))
Idents(object = P11_ARKO_FibSM) <- "seurat_clusters"
P11_ARKO_FibSM <- RenameIdents(object = P11_ARKO_FibSM, '1' = "Fibroblast", '3' = "Fibroblast", '5' = "Fibroblast", '7' = "Fibroblast", '0' = "SM", '2' = "SM", '6' = "SM")
P11_ARKO_FibSM[["CellTypes"]] <- Idents(object = P11_ARKO_FibSM)
Idents(object = P11_ARKO_FibSM) <- "CellTypes"
DimPlot(P11_ARKO_FibSM, reduction = "umap", pt.size = 1)

#### P42 Ctrl####
Idents(object = P42_Ctrl) <- "seurat_clusters"
P42_Ctrl <- ScaleData(P42_Ctrl, verbose = FALSE)
P42_Ctrl <- RunPCA(P42_Ctrl, npcs = 50, verbose = FALSE)
ElbowPlot(P42_Ctrl)

P42_Ctrl <- FindNeighbors(P42_Ctrl, reduction = "pca", dims = 1:20)
P42_Ctrl <- FindClusters(P42_Ctrl, resolution = 0.5)
P42_Ctrl <- RunUMAP(P42_Ctrl, reduction = "pca", dims = 1:20)
DimPlot(P42_Ctrl, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(P42_Ctrl) <- "RNA"
FeaturePlot(P42_Ctrl, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_Ctrl, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_Ctrl, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
P42_Ctrl_FibSM <- subset(P42_Ctrl, idents = c("1", "4", "16", "12", "13"))
Idents(object = P42_Ctrl_FibSM) <- "seurat_clusters"
P42_Ctrl_FibSM <- RenameIdents(object = P42_Ctrl_FibSM, '1' = "Fibroblast", '4' = "Fibroblast", '16' = "Fibroblast", '12' = "SM", '13' = "SM")
P42_Ctrl_FibSM[["CellTypes"]] <- Idents(object = P42_Ctrl_FibSM)
Idents(object = P42_Ctrl_FibSM) <- "CellTypes"
DimPlot(P42_Ctrl_FibSM, reduction = "umap", pt.size = 1)

#### P42 ARKO####
Idents(object = P42_ARKO) <- "seurat_clusters"
P42_ARKO <- ScaleData(P42_ARKO, verbose = FALSE)
P42_ARKO <- RunPCA(P42_ARKO, npcs = 50, verbose = FALSE)
ElbowPlot(P42_ARKO)

P42_ARKO <- FindNeighbors(P42_ARKO, reduction = "pca", dims = 1:15)
P42_ARKO <- FindClusters(P42_ARKO, resolution = 0.5)
P42_ARKO <- RunUMAP(P42_ARKO, reduction = "pca", dims = 1:15)
DimPlot(P42_ARKO, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(P42_ARKO) <- "RNA"
FeaturePlot(P42_ARKO, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_ARKO, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_ARKO, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#subclustering Fib/SM
P42_ARKO_FibSM <- subset(P42_ARKO, idents = c("0", "2", "5", "1", "9"))
Idents(object = P42_ARKO_FibSM) <- "seurat_clusters"
P42_ARKO_FibSM <- RenameIdents(object = P42_ARKO_FibSM, '0' = "Fibroblast", '2' = "Fibroblast", '5' = "Fibroblast", '1' = "SM", '9' = "SM")
P42_ARKO_FibSM[["CellTypes"]] <- Idents(object = P42_ARKO_FibSM)
Idents(object = P42_ARKO_FibSM) <- "CellTypes"
DimPlot(P42_ARKO_FibSM, reduction = "umap", pt.size = 1)

#### Merging 6 FibSM ####

#Set Current idents
Idents(object = E18_Ctrl_FibSM) <- "seurat_clusters"
Idents(object = E18_ARKO_FibSM) <- "seurat_clusters"
Idents(object = P11_Ctrl_FibSM) <- "seurat_clusters"
Idents(object = P11_ARKO_FibSM) <- "seurat_clusters"
Idents(object = P42_Ctrl_FibSM) <- "seurat_clusters"
Idents(object = P42_ARKO_FibSM) <- "seurat_clusters"
E18_Ctrl_FibSM$stim <- "E18_Ctrl"
E18_ARKO_FibSM$stim <- "E18_ARKO"
P11_Ctrl_FibSM$stim <- "P11_Ctrl"
P11_ARKO_FibSM$stim <- "P11_ARKO"
P42_Ctrl_FibSM$stim <- "P42_Ctrl"
P42_ARKO_FibSM$stim <- "P42_ARKO"
ARKOvCtrl_FibSM.anchors <- FindIntegrationAnchors(object.list = list(E18_Ctrl_FibSM, E18_ARKO_FibSM, P11_Ctrl_FibSM, P11_ARKO_FibSM, P42_Ctrl_FibSM, P42_ARKO_FibSM), dims = 1:20)
ARKOvCtrl_FibSM<- IntegrateData(anchorset = ARKOvCtrl_FibSM.anchors, dims = 1:20)
DefaultAssay(ARKOvCtrl_FibSM) <- "integrated"

#Run the standard workflow for visualization and clustering
ARKOvCtrl_FibSM <- ScaleData(ARKOvCtrl_FibSM, verbose = FALSE)
ARKOvCtrl_FibSM <- RunPCA(ARKOvCtrl_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_FibSM)

# umap and Clustering
ARKOvCtrl_FibSM <- FindNeighbors(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:20)
ARKOvCtrl_FibSM <- FindClusters(ARKOvCtrl_FibSM, resolution = 0.5)
ARKOvCtrl_FibSM <- RunUMAP(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:20)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 

Idents(object = ARKOvCtrl_FibSM) <- "CellTypes"
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim")

#Cell count (split)
Idents(object = ARKOvCtrl_FibSM) <- "stim"
ARKOvCtrl_FibSM$stim.seurat_clusters <- paste(Idents(ARKOvCtrl_FibSM), ARKOvCtrl_FibSM$seurat_clusters, sep = "_")
Idents(object = ARKOvCtrl_FibSM) <- "stim.seurat_clusters"
table(Idents(ARKOvCtrl_FibSM))

#Featureplot
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Heatmap
Idents(object = ARKOvCtrl_FibSM) <- "seurat_clusters"
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
ARKOvCtrl_FibSM <- ScaleData(ARKOvCtrl_FibSM, features = rownames(ARKOvCtrl_FibSM))
ARKOvCtrl_FibSM.markers <- FindAllMarkers(ARKOvCtrl_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ARKOvCtrl_FibSMTop10 <- ARKOvCtrl_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ARKOvCtrl_FibSM, features = c(ARKOvCtrl_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
ARKOvCtrl_FibSM.0.1markers <- FindAllMarkers(ARKOvCtrl_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(ARKOvCtrl_FibSM.0.1markers, "ARKOvCtrl_FibSM.0.1markers.csv")

#### eliminate Epi cells ####
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Add ARQ info
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
ARKOvCtrl_FibSMCdh1Pos <- subset(x=ARKOvCtrl_FibSM, subset = Cdh1 > 0)
ARKOvCtrl_FibSMCdh1Neg <- subset(x=ARKOvCtrl_FibSM, subset = Cdh1 == 0)
Idents(object = ARKOvCtrl_FibSMCdh1Pos) <- "Cdh1Pos"
Idents(object = ARKOvCtrl_FibSMCdh1Neg) <- "Cdh1Neg"
ARKOvCtrl_FibSMCdh1Pos[["Cdh1Exp"]] <- Idents(object = ARKOvCtrl_FibSMCdh1Pos)
ARKOvCtrl_FibSMCdh1Neg[["Cdh1Exp"]] <- Idents(object = ARKOvCtrl_FibSMCdh1Neg)
ARKOvCtrl_FibSMCdh1  <- merge(x = ARKOvCtrl_FibSMCdh1Pos, y = ARKOvCtrl_FibSMCdh1Neg)
Idents(object = ARKOvCtrl_FibSMCdh1) <- "Cdh1Exp"
ARKOvCtrl_FibSM$Cdh1Exp <- Idents(object = ARKOvCtrl_FibSMCdh1)

Idents(object = ARKOvCtrl_FibSM) <- "Cdh1Exp"
ARKOvCtrl_FibSM <- subset(ARKOvCtrl_FibSM, idents = c("Cdh1Neg"))
DimPlot(ARKOvCtrl_FibSM, reduction = "tsne", pt.size = 0.3)

#Run the standard workflow for visualization and clustering
DefaultAssay(ARKOvCtrl_FibSM) <- "integrated"

ARKOvCtrl_FibSM <- ScaleData(ARKOvCtrl_FibSM, verbose = FALSE)
ARKOvCtrl_FibSM <- RunPCA(ARKOvCtrl_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_FibSM)

# umap and Clustering
ARKOvCtrl_FibSM <- FindNeighbors(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:15)
ARKOvCtrl_FibSM <- FindClusters(ARKOvCtrl_FibSM, resolution = 0.5)
ARKOvCtrl_FibSM <- RunUMAP(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:15)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = ARKOvCtrl_FibSM) <- "seurat_clusters"
tiff(file = "ARKOvCtrl_FibSM_umap.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_FibSM_split_umap.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

Idents(object = ARKOvCtrl_FibSM) <- "CellTypes"
tiff(file = "ARKOvCtrl_FibSM_cellytype_umap.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3) 
dev.off()
tiff(file = "ARKOvCtrl_FibSM_cellytype_split_umap.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell count (split)
Idents(object = ARKOvCtrl_FibSM) <- "stim"
ARKOvCtrl_FibSM$stim.seurat_clusters <- paste(Idents(ARKOvCtrl_FibSM), ARKOvCtrl_FibSM$seurat_clusters, sep = "_")
Idents(object = ARKOvCtrl_FibSM) <- "stim.seurat_clusters"
table(Idents(ARKOvCtrl_FibSM))

#Featureplot
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
tiff(file = "ARKOvCtrl_FibSM_split_EGFP.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM_split_Ar.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM_split_Gli1.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "ARKOvCtrl_FibSM_EGFP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM_Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM_Gli1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = ARKOvCtrl_FibSM) <- "seurat_clusters"
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
ARKOvCtrl_FibSM <- ScaleData(ARKOvCtrl_FibSM, features = rownames(ARKOvCtrl_FibSM))
ARKOvCtrl_FibSM.markers <- FindAllMarkers(ARKOvCtrl_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ARKOvCtrl_FibSMTop10 <- ARKOvCtrl_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(ARKOvCtrl_FibSM, features = c(ARKOvCtrl_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
ARKOvCtrl_FibSM.0.1markers <- FindAllMarkers(ARKOvCtrl_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(ARKOvCtrl_FibSM.0.1markers, "ARKOvCtrl_FibSM.0.1markers.csv")


#### Pseudotime 6 FibSM ####

DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
FibSMPseudo <- as.CellDataSet(ARKOvCtrl_FibSM)
FibSMPseudo <- detectGenes(FibSMPseudo, min_expr = 0.1)
print(head(fData(FibSMPseudo)))

expressed_genes <- row.names(subset(fData(FibSMPseudo),
                                    num_cells_expressed >= 10))

pData(FibSMPseudo)$Total_mRNAs <- Matrix::colSums(exprs(FibSMPseudo))
FibSMPseudo <- FibSMPseudo[,pData(FibSMPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(FibSMPseudo), geom =
        "density")

FibSMPseudo <- estimateSizeFactors(FibSMPseudo)
FibSMPseudo <- estimateDispersions(FibSMPseudo)

disp_table <- dispersionTable(FibSMPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
FibSMPseudo <- setOrderingFilter(FibSMPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(FibSMPseudo)

#FibSMPseudo@auxClusteringData[["umap"]]$variance_explained <- NULL
plot_pc_variance_explained(FibSMPseudo, return_all = F) # norm_method='log'

FibSMPseudo <- reduceDimension(FibSMPseudo, max_components = 2, num_dim = 20,
                             reduction_method = 'umap', verbose = T)
FibSMPseudo <- clusterCells(FibSMPseudo, num_clusters = 2)

plot_cell_clusters(FibSMPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(FibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
FibSMPseudo <- setOrderingFilter(FibSMPseudo, ordering_genes)
plot_ordering_genes(FibSMPseudo)

FibSMPseudo <- reduceDimension(FibSMPseudo, max_components = 2,
                             method = 'DDRTree')

FibSMPseudo <- orderCells(FibSMPseudo)

GM_state <- function(FibSMPseudo){
  if (length(unique(pData(FibSMPseudo)$State)) > 1){
    T0_counts <- table(pData(FibSMPseudo)$State, pData(FibSMPseudo)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

FibSMPseudo <- orderCells(FibSMPseudo, root_state = GM_state(FibSMPseudo))

#Visualization
plot_cell_trajectory(FibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

plot_cell_trajectory(FibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))

