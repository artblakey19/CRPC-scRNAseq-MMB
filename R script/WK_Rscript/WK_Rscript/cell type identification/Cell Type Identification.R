#Cell type identification

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(RColorBrewer)
library(SeuratData)
library(BiocManager)
library(SeuratWrappers)
library(monocle3)
library(clustree)

####WT_E18.5####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Cell Type Identification/E18.5")

DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1 <- ScaleData(E18_Ctrl1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_Ctrl1))
E18_Ctrl1 <- RunPCA(E18_Ctrl1, features = VariableFeatures(E18_Ctrl1))
ElbowPlot(E18_Ctrl1, ndims = 50)

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 0.3)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res0.3.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 0.5)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 0.8)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 1)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 1.2)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res1.2.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 1.5)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res1.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 1.8)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res1.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 2.0)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res2.0.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 2.2)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res2.2.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 2.5)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res2.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 2.8)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res2.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 3.0)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res3.0.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#clustree
clustree(E18_Ctrl1, prefix = "RNA_snn_res.")
tiff(file = "E18_Ctrl1 clustree.tiff", width = 12, height = 10, units = "in", compression = "lzw", res = 200)
clustree(E18_Ctrl1, prefix = "RNA_snn_res.")
dev.off()

#DEGs
Idents(object = E18_Ctrl1) <- "seurat_clusters"
DefaultAssay(E18_Ctrl1) <- "RNA"
all.genes <- rownames(E18_Ctrl1)
E18_Ctrl1 <- ScaleData(E18_Ctrl1, features = all.genes)
E18_Ctrl1.seurat_cluster.markers <- FindAllMarkers(E18_Ctrl1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(E18_Ctrl1.seurat_cluster.markers, file = "E18_Ctrl1.seurat_cluster.markers.csv")

#Heatmap
DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1.seurat_clusterTop30 <- E18_Ctrl1.seurat_cluster.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
tiff(file = "E18_Ctrl1 seurat res3 Heatmap Top30.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 200)
DoHeatmap(E18_Ctrl1, features = c(E18_Ctrl1.seurat_clusterTop30$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####Optimized resolution####
DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1 <- ScaleData(E18_Ctrl1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_Ctrl1))
E18_Ctrl1 <- RunPCA(E18_Ctrl1, features = VariableFeatures(E18_Ctrl1))
ElbowPlot(E18_Ctrl1, ndims = 50)

tiff(file = "E18_Ctrl1 ElbowPlot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
ElbowPlot(E18_Ctrl1, ndims = 50)
dev.off()

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 1.6)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:32)
E18_Ctrl1 <- RunTSNE(E18_Ctrl1, reduction = "pca", dims = 1:32)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims32 res1.6.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#DEGs
Idents(object = E18_Ctrl1) <- "seurat_clusters"
DefaultAssay(E18_Ctrl1) <- "RNA"
all.genes <- rownames(E18_Ctrl1)
E18_Ctrl1 <- ScaleData(E18_Ctrl1, features = all.genes)
E18_Ctrl1.seurat_cluster.res1.6.markers <- FindAllMarkers(E18_Ctrl1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(E18_Ctrl1.seurat_cluster.res1.6.markers, file = "E18_Ctrl1.seurat_cluster.res1.6.markers.csv")

#Heatmap
DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1.seurat_clusterTop30 <- E18_Ctrl1.seurat_cluster.res1.6.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
tiff(file = "E18_Ctrl1 seurat res1.6 Heatmap Top30.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 200)
DoHeatmap(E18_Ctrl1, features = c(E18_Ctrl1.seurat_clusterTop30$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#FeaturePlot
DefaultAssay(E18_Ctrl1)<-"RNA"
tiff(file = "E18_Ctrl1 celltype marker expression plots.tiff", width = 12, height = 40, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Epcam", "Cdh1", 
                                                                          "Krt14", "Krt15", "Krt19", "Anxa1", "Psca", 
                                                                          "Upk3a", "Upk1a", "Foxq1", "Krt7", "Krt20",
                                                                          "Plac8", "Cldn3", "Cldn4", "Hoxb9", "Pax8", "Wfdc2", "Gata3", "Pax2",
                                                                          "Snai1", "Crabp1", "Wnt9b", "Foxf1", "Acta2", "Scmh1",
                                                                          "Vim",  "Fbln1", "Myh11", "Tagln", "Actg2", 
                                                                          "Syp", "Chgb", "Scg2", 
                                                                          "Rgs5", "Apold1", "Esam", 
                                                                          "Plp1", "Mpz", "Cdh19",
                                                                          "Cdh5", "Cldn5", "Pecam1",
                                                                          "Myog", "Actc1", "Mylpf",
                                                                          "Tyrobp", "C1qa","Rgs1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

##mesonephric derivatives markers
#Featureplots
DefaultAssay(E18_Ctrl1)<-"RNA"
tiff(file = "E18_Ctrl1 Hoxb7 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Hoxb7"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Wfdc2 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Wfdc2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Gsta3 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Gsta3"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Sox17 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Sox17"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Pax2 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Pax2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Pax8 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Pax8"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Lhx1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Lhx1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Cdh16 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Cdh16"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Hoxb9 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Hoxb9"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

##Bladder Urothelium markers
#Featureplots
tiff(file = "E18_Ctrl1 Upk3a expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Upk3a"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Foxq1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Foxq1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt7 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt7"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt20 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt20"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Sptssb expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Sptssb"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Upk2 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Upk2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Snx31 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Snx31"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

##UGE markers
#Featureplots
tiff(file = "E18_Ctrl1 Krt15 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt15"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Cxcl15 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Cxcl15"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Tspan1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Tspan1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt14 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt14"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Hoxb13 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Hoxb13"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Hoxd13 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Hoxd13"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

##Distinguish Superficial, Intermediate, Basal, Budding
#Featureplots
tiff(file = "E18_Ctrl1 Trp63 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Trp63"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt5 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt13 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt13"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt14 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt14"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Krt20 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Krt20"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Foxa2 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Foxa2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Clca1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Clca1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Col17a1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Col17a1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Wif1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Wif1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Etv4 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Etv4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Nkx3-1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Nkx3-1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q80")
dev.off()
tiff(file = "E18_Ctrl1 Upk1b expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Upk1b"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
#Violin plot
DefaultAssay(E18_Ctrl1)<-"RNA"
tiff(file = "E18_Ctrl1 Upk3b Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Upk3b"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Upk3a Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Upk3a"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Trp63 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Trp63"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Krt5 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Krt5"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Krt13 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Krt13"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Krt14 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Krt14"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Krt20 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Krt20"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Foxa2 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Foxa2"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Nkx3-1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Nkx3-1"), pt.size = 0)
dev.off()

##SVM, Lamina propria, Muscularis muscosa, Submucosa, Humcularis propria, Mesenchymal pad, Urethra
#Featureplots
tiff(file = "E18_Ctrl1 Crabp1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Crabp1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Snai1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Snai1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Foxf1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Foxf1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Acta2 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Fgf10 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Fgf10"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
#Vlnplots
tiff(file = "E18_Ctrl1 Crabp1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Crabp1"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Snai1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Snai1"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Foxf1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Foxf1"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Acta2 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Acta2"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Myh11 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Myh11"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Fgf10 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Fgf10"), pt.size = 0)
dev.off()

##Telocyte markers
#Featureplots
tiff(file = "E18_Ctrl1 Pdgfra expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Cd34 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Kit expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Kit"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Foxl1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 GFP expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = 2.0, max.cutoff = "q90")
dev.off()
#Vlnplots
tiff(file = "E18_Ctrl1 Pdgfra Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Pdgfra"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Cd34 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Cd34"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Kit Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Kit"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1 Foxl1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1, features = c("Foxl1"), pt.size = 0)
dev.off()

#Rename
Idents(object = E18_Ctrl1) <- "seurat_clusters"
E18_Ctrl1 <- RenameIdents(object = E18_Ctrl1, 
                                      '29' = "UGE_Budding", '14' ="UGE_Basal", '18' ="UGE_Superficial",
                          '11' ="Bladder_Basal",'26' ="Bladder_Superficial",
                                      '7' ="Wolffian Duct",
                          '4' ="Telocyte",'3' ="UGM_Lamina propria1",'17' ="UGM_Lamina propria2",
                          '28' ="UGM_Muscularis mucosa",'5' ="UGM_Submucosa1",
                                      '8' ="UGM_Submucosa2",'1' ="UGM_Muscularis propria1",'0' ="UGM_Muscularis propria2",
                          '12' ="UGM_MP",'21' ="Urethra_Mesenchyme1",'22' ="Urethra_Mesenchyme2",
                          '13' ="Urethra_Muscularis", '6'="WD_Mesenchyme1", '2'="WD_Mesenchyme2",
                          '16'="WD_Muscularis1", '10'="WD_Muscularis2",
                                      '20' = "Pericyte", '9' = "Glia1", '27' = "Glia2",'19' = "Endothelia",
                                      '24' = "Myoblast",'15' = "Neurocrest1",
                                      '23' = "Neurocrest2",'25' = "Immune")
E18_Ctrl1[["CellTypes"]] <- Idents(object = E18_Ctrl1)

#Umap
Idents(object = E18_Ctrl1) <- "CellTypes"
tiff(file = "E18_Ctrl1 UMAP CellTypes.tiff", width = 10, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "E18_Ctrl1 UMAP CellTypes label.tiff", width = 12, height = 8, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

        

####Subset FBSM####
Idents(object = E18_Ctrl1) <- "CellTypes"
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE) 

E18_Ctrl1_FBSM <- subset(E18_Ctrl1, idents = c("Telocyte",
                                           "UGM_Lamina propria1","UGM_Lamina propria2",
                                           "UGM_Muscularis mucosa",
                                           "UGM_Submucosa1","UGM_Submucosa2",
                                           "UGM_Muscularis propria1","UGM_Muscularis propria2",
                                           "UGM_MP","UGM_Lamina propria1"))

DefaultAssay(E18_Ctrl1_FBSM) <- "RNA"
#Run the standard workflow for visualization and clustering
E18_Ctrl1_FBSM <- ScaleData(E18_Ctrl1_FBSM, verbose = FALSE)
E18_Ctrl1_FBSM <- RunPCA(E18_Ctrl1_FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(E18_Ctrl1_FBSM, ndims = 50)
# UMAP and Clustering
E18_Ctrl1_FBSM <- FindNeighbors(E18_Ctrl1_FBSM, reduction = "pca", dims = 1:23)
E18_Ctrl1_FBSM <- FindClusters(E18_Ctrl1_FBSM, resolution = 0.8)
E18_Ctrl1_FBSM <- RunUMAP(E18_Ctrl1_FBSM, reduction = "pca", dims = 1:23)

Idents(object = E18_Ctrl1_FBSM) <- "seurat_clusters"
DimPlot(E18_Ctrl1_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell Cycle scoring
DefaultAssay(E18_Ctrl1_FBSM) <- "RNA"
all.genes <- rownames(E18_Ctrl1_FBSM)
E18_Ctrl1_FBSM <- ScaleData(E18_Ctrl1_FBSM, features = all.genes)
E18_Ctrl1_FBSM <- CellCycleScoring(E18_Ctrl1_FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_Ctrl1_FBSM) <- "Phase"
E18_Ctrl1_FBSM <- RenameIdents(object = E18_Ctrl1_FBSM, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
E18_Ctrl1_FBSM[["Phase"]] <- Idents(object = E18_Ctrl1_FBSM)
DimPlot(E18_Ctrl1_FBSM, reduction = "umap")

tiff(file = "E18_Ctrl1_FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM, reduction = "umap", pt.size = 0.5, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#Cell Cycle regression
E18_Ctrl1_FBSM1 <- E18_Ctrl1_FBSM
DefaultAssay(E18_Ctrl1_FBSM1) <- "RNA"
E18_Ctrl1_FBSM1 <- ScaleData(E18_Ctrl1_FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_Ctrl1_FBSM1))
E18_Ctrl1_FBSM1 <- RunPCA(E18_Ctrl1_FBSM1, features = VariableFeatures(E18_Ctrl1_FBSM1))
ElbowPlot(E18_Ctrl1_FBSM1, ndims = 50)

E18_Ctrl1_FBSM1 <- FindNeighbors(E18_Ctrl1_FBSM1, reduction = "pca", dims = 1:26)
E18_Ctrl1_FBSM1 <- FindClusters(E18_Ctrl1_FBSM1, resolution = 0.8)
E18_Ctrl1_FBSM1 <- RunUMAP(E18_Ctrl1_FBSM1, reduction = "pca", dims = 1:26)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_Ctrl1_FBSM1) <- "seurat_clusters"
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E18_Ctrl1_FBSM1 seurat_clusters UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()


Idents(object = E18_Ctrl1_FBSM1) <- "Phase"
E18_Ctrl1_FBSM1 <- RenameIdents(object = E18_Ctrl1_FBSM1, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
E18_Ctrl1_FBSM1[["Phase"]] <- Idents(object = E18_Ctrl1_FBSM1)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap")

Idents(object = E18_Ctrl1_FBSM1) <- "Phase"
tiff(file = "E18_Ctrl1_FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.5, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#Dotplot
Idents(object = E18_Ctrl1_FBSM1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1_FBSM1 CellTypes stromal markers DotPlot.tiff", width =10 , height = 4, units = "in", compression = "lzw", res = 200)
DotPlot(E18_Ctrl1_FBSM1, features = c("Epcam", "Cdh1", 
                                      "Snai1", "Crabp1", "Foxf1", "Acta2", "Scmh1",
                                      "Pdgfra", "Cd34", "Kit", "Foxl1", 
                                      "Vim",  "Fbln1", "Myh11", "Tagln", "Actg2", "Fgf10"),
        cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#FeaturePlot
DefaultAssay(E18_Ctrl1_FBSM1)<-"RNA"
tiff(file = "E18_Ctrl1_FBSM1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Epcam", "Cdh1", 
                                                        "Snai1", "Crabp1", "Foxf1", "Acta2", "Scmh1",
                                                        "Pdgfra", "Cd34", "Kit", "Foxl1",
                                                        "Vim",  "Fbln1", "Myh11", "Tagln", "Fgf10"
                                                        ) 
                                                        , cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = E18_Ctrl1_FBSM1) <- "seurat_clusters"
E18_Ctrl1_FBSM1 <- RenameIdents(object = E18_Ctrl1_FBSM1, 
                             '1' = "Telocyte1", '9' ="Lamina propria1", 
                             '7' ="Lamina propria2",'2' ="Lamina propria3",
                             '12'="Muscularis Mucosa",
                             '4' ="Submucosa1",'10' ="Submucosa2",
                             '3'= "Muscularis propria1",'0'= "Muscularis propria2",
                             '5'= "Muscularis propria3",'6'= "Muscularis propria4",
                             '11'= "Muscularis propria5", '8'="OS")
E18_Ctrl1_FBSM1[["FBSMCellTypes"]] <- Idents(object = E18_Ctrl1_FBSM1)

#Umap
Idents(object = E18_Ctrl1_FBSM1) <- "FBSMCellTypes"
tiff(file = "E18_Ctrl1_FBSM1 UMAP FBSMCellTypes.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 UMAP FBSMCellTypes label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#DEGs
Idents(object = E18_Ctrl1_FBSM1) <- "FBSMCellTypes"
DefaultAssay(E18_Ctrl1_FBSM1) <- "RNA"
all.genes <- rownames(E18_Ctrl1_FBSM1)
E18_Ctrl1_FBSM1 <- ScaleData(E18_Ctrl1_FBSM1, features = all.genes)
E18_Ctrl1_FBSM1.FBSMCellTypes.markers <- FindAllMarkers(E18_Ctrl1_FBSM1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(E18_Ctrl1_FBSM1.FBSMCellTypes.markers, file = "E18_Ctrl1_FBSM1.FBSMCellTypes.markers.csv")

#Featureplot
tiff(file = "E18_Ctrl1_FBSM1 Pdgfra expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q80")
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Cd34 expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Foxl1 expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Kit expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Kit"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

#Vlnplot
tiff(file = "E18_Ctrl1_FBSM1 Pdgfra Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1_FBSM1, features = c("Pdgfra"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Cd34 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1_FBSM1, features = c("Cd34"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Foxl1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1_FBSM1, features = c("Foxl1"), pt.size = 0)
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Kit Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl1_FBSM1, features = c("Kit"), pt.size = 0)
dev.off()

####WT_P60####
mouse_cell_cycle_genes <- readRDS(file = "/Volumes/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Cell Type Identification/P60")

#loading data
WT_P60.unfiltered.data <- Read10X("/Volumes/user_data/zjsun/seq/200108_Tgen_scRNA/33580_count_EGFPmm10_v3/outs/filtered_feature_bc_matrix")
WT_P60.unfiltered <- CreateSeuratObject(counts = WT_P60.unfiltered.data,  min.cells = 3, min.features = 250, project = "WT_P60.unfiltered")
WT_P60.unfiltered <- NormalizeData(WT_P60.unfiltered)

WT_P60.unfiltered[["percent.mt"]] <- PercentageFeatureSet(WT_P60.unfiltered, pattern = "^mt-")

table(Idents(WT_P60.unfiltered))

tiff(file = "WT_P60 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()

tiff(file = "WT_P60 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(WT_P60.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "WT_P60 pre filteration")
dev.off()

tiff(file = "WT_P60 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(WT_P60.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "WT_P60 pre filteration")
dev.off()

#Filtering
WT_P60 <- subset(WT_P60.unfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7200 & percent.mt < 15)
Idents(object = WT_P60) <- 'WT_P60'

table(Idents(WT_P60))
tiff(file = "WT_P60 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3,)
dev.off()
tiff(file = "WT_P60 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(WT_P60@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "WT_P60 post filteration")
dev.off()
tiff(file = "WT_P60 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(WT_P60@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "WT_P60 post filteration")
dev.off()

WT_P60 <- FindVariableFeatures(WT_P60, selection.method = "vst", nfeatures = 2500)
tiff(file = "WT_P60 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VariableFeaturePlot(WT_P60)
dev.off()

all.genes <- rownames(WT_P60)
WT_P60 <- ScaleData(WT_P60, features = all.genes)
WT_P60 <- RunPCA(WT_P60, features = VariableFeatures(object = WT_P60))
DimPlot(WT_P60, reduction = "pca")
ElbowPlot(WT_P60, ndims = 50)

WT_P60 <- FindNeighbors(WT_P60, dims = 1:20)
WT_P60 <- FindClusters(WT_P60, resolution = 0.5)
WT_P60 <- RunUMAP(WT_P60, dims = 1:20)
DimPlot(WT_P60, reduction = "umap", pt.size = 0.5, label = TRUE)

#Cell cycle assignment
DefaultAssay(WT_P60) <- "RNA"
all.genes <- rownames(WT_P60)
WT_P60 <- ScaleData(WT_P60, features = all.genes)
WT_P60 <- CellCycleScoring(WT_P60, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = WT_P60) <- "Phase"
DimPlot(WT_P60, reduction = "umap")
tiff(file = "WT_P60 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT_P60, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell cycle regression
WT_P60_1 <- WT_P60
DefaultAssay(WT_P60_1) <- "RNA"
WT_P60_1 <- ScaleData(WT_P60_1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(WT_P60_1))
WT_P60_1 <- RunPCA(WT_P60_1, features = VariableFeatures(WT_P60_1))
ElbowPlot(WT_P60_1, ndims = 50)

tiff(file = "WT_P60_1 ElbowPlot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
ElbowPlot(WT_P60_1, ndims = 50)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 0.5)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 0.8)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 1.0)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res1.0.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 1.2)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res1.2.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 1.5)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res1.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 1.8)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res1.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 2.0)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res2.0.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 2.2)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res2.2.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 2.5)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res2.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 2.8)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res2.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 3.0)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res3.0.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#clustree
clustree(WT_P60_1, prefix = "RNA_snn_res.")
tiff(file = "WT_P60_1 clustree.tiff", width = 12, height = 10, units = "in", compression = "lzw", res = 200)
clustree(WT_P60_1, prefix = "RNA_snn_res.")
dev.off()

#DEGs
Idents(object = WT_P60_1) <- "seurat_clusters"
DefaultAssay(WT_P60_1) <- "RNA"
all.genes <- rownames(WT_P60_1)
WT_P60_1 <- ScaleData(WT_P60_1, features = all.genes)
WT_P60_1.seurat_cluster.markers <- FindAllMarkers(WT_P60_1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(WT_P60_1.seurat_cluster.markers, file = "WT_P60_1.seurat_cluster.markers.csv")

#Heatmap
DefaultAssay(WT_P60_1) <- "RNA"
WT_P60_1.seurat_clusterTop30 <- WT_P60_1.seurat_cluster.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
tiff(file = "WT_P60_1 seurat res3 Heatmap Top30.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 200)
DoHeatmap(WT_P60_1, features = c(WT_P60_1.seurat_clusterTop30$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####Optimized resolution####
WT_P60_1 <- FindNeighbors(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- FindClusters(WT_P60_1, resolution = 1.8)
WT_P60_1 <- RunUMAP(WT_P60_1, reduction = "pca", dims = 1:32)
WT_P60_1 <- RunTSNE(WT_P60_1, reduction = "pca", dims = 1:32)

Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 UMAP dims32 res1.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = WT_P60_1) <- "Phase"
DimPlot(WT_P60_1, reduction = "umap")
tiff(file = "WT_P60_1 res1.8 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#DEGs
Idents(object = WT_P60_1) <- "seurat_clusters"
DefaultAssay(WT_P60_1) <- "RNA"
all.genes <- rownames(WT_P60_1)
WT_P60_1 <- ScaleData(WT_P60_1, features = all.genes)
WT_P60_1.seurat_cluster.res1.8.markers <- FindAllMarkers(WT_P60_1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(WT_P60_1.seurat_cluster.res1.8.markers, file = "WT_P60_1.seurat_cluster.res1.8.markers.csv")

#Heatmap
DefaultAssay(WT_P60_1) <- "RNA"
WT_P60_1.seurat_clusterTop30 <- WT_P60_1.seurat_cluster.res1.8.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
tiff(file = "WT_P60_1 seurat res1.8 Heatmap Top30.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 200)
DoHeatmap(WT_P60_1, features = c(WT_P60_1.seurat_clusterTop30$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Dotplot
Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 CellTypes markers DotPlot.tiff", width =20 , height = 6, units = "in", compression = "lzw", res = 200)
DotPlot(WT_P60_1, features = c("Gsdma",	"Tgm4",	"Fut9",	"C1s2",	"Mt3",	"Ren1",	
                               "Glb1l2",	"Tcaf3", "Nlrp10",	"2210407C18Rik",	"Ms4a5",	
                               "A630095E13Rik",	"Mgll",	"Defb10",	"Sp5",	"H2-Q10",
                               "Pcp4",	"Upb1",	"Msmb",	"Abo", "Trpv6",	"Scgb2b20",
                               "Scgb2b7",	"Nupr1",	"Slc26a4",	"Csta1",	"Dcaf12l1",	
                               "Sbp",	"Sbpl",	"Crabp1",	"Spink1",	"Ppp1r1b",	"Spns2",
                               "Scnn1a",	"Pglyrp1",	"Gprc5a",	"Psca",	"Krt4",	
                               "Tacstd2",	"Ly6a",	"Lgals7",	"Aqp3",	"Grp",	"Dapl1",
                               "Anxa8",	"Fxyd4",	"Dlk2",	"Col18a1",	"Fmod",	"Krt14",	"Scpep1",	"Trp63",
                               "Upk3a", "Upk1a", "Foxq1",
                               "Plac8", "Pax8", "Svs6", "Svs2", "Svs4"),
        cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Dotplot
Idents(object = WT_P60_1) <- "seurat_clusters"
tiff(file = "WT_P60_1 CellTypes stromal markers DotPlot.tiff", width =20 , height = 6, units = "in", compression = "lzw", res = 200)
DotPlot(WT_P60_1, features = c("Rorb", "F2r", "Adm", "Rdh10", "Fzd1",
                               "Wnt2", "Wif1", "Ifitm1", "Srd5a2",
                               "Sult1e1", "Lum", "Ogn", "Dpt", "C3",
                               "Ctgf", "Clec3b", "Dcn", "Ebf1", "Gpx3", "Igf1",
                               "Lgr5", "Apoe", "Osr1", "Sfrp2", "Mfap4",
                               "Foxf1", "Foxl1", "Pdgfra", "Cd34", "Kit", "Gli1", "Ar",
                               "C1qa", "Cd68",
                               "Mcemp1", "Ly75",
                               "Cd3g", "Cd2",
                               "Cdh5", "Pecam1",
                               "Rgs4", "Pcp4l1",
                               "Plp1", "Cdh19"),
        cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#FeaturePlot
DefaultAssay(WT_P60_1)<-"RNA"
tiff(file = "WT_P60_1 celltype marker expression plots.tiff", width = 16, height = 40, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_1, reduction = "umap", features = c("Epcam", "Cdh1",
                                                       "Krt14", "Lgals7",
                                                       "C1s2", "Ren1",
                                                       "Glb1l2", "Nlrp10",
                                                       "Defb10", "Sp5",
                                                       "Trpv6", "Slc26a4",
                                                       "Ppp1r1b", "Psca",
                                                       "Fbln1", "Fbn1",
                                                       "Acta2", "Myh11",
                                                       "Pdgfra", "Foxl1", "Spon1",
                                                       "Rorb", "Wif1",
                                                       "Sult1e1", "Lum",
                                                       "Lgr5", "Osr1",
                                                      "C1qa", "Cd68",
                                                      "Mcemp1", "Ly75",
                                                      "Cd3g", "Cd2",
                                                      "Cdh5", "Pecam1",
                                                      "Rgs5", "Gja4",
                                                      "Plp1", "Cdh19"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = WT_P60_1) <- "seurat_clusters"
WT_P60_1 <- RenameIdents(object = WT_P60_1, 
                          '7' = "BE1", '3' ="BE2", '2' ="BE3",'5' ="BE4",
                          '1' ="AP1_Hoxb9+",'17' ="AP2_Defb50+", '28'= "AP3_Krt20+",
                          '9' ="DP1", '21'="DP2_Fos", 
                          '13' ="LP",'4' ="VP1_Trpv6+",'10' ="VP2_Rpl10+",'11' ="VP3_Cxcl17+", '24'="Prog/UrLE",
                          '29' ="OE",'6' ="Telocyte",
                          '0' ="SubEpithelial FB1", '19' ="SubEpithelial FB2",'8' ="Intestitial FB1",
                          '12' ="Intestitial FB2",'23' ="Urethral FB",'20' ="SM",
                          '18' ="Glia", '25'="Pericyte", '27'="Pericyte",
                          '16'="non-lymphatic VE", '30'="lymphatic VE",
                          '14' = "Monocyte", '22' = "Monocyte-derived DC", '26' = "Mast Cell",'15' = "T cell")
WT_P60_1[["CellTypes"]] <- Idents(object = WT_P60_1)

#Umap
Idents(object = WT_P60_1) <- "CellTypes"
tiff(file = "WT_P60_1 UMAP CellTypes.tiff", width = 10, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "WT_P60_1 UMAP CellTypes label.tiff", width = 12, height = 8, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#Rename_CellTypes1
Idents(object = WT_P60_1) <- "seurat_clusters"
WT_P60_1 <- RenameIdents(object = WT_P60_1, 
                         '7' = "BE", '3' ="BE", '2' ="BE",'5' ="BE",
                         '1' ="AP",'17' ="AP", '28'= "AP",
                         '9' ="DP", '21'="DP", 
                         '13' ="LP",'4' ="VP",'10' ="VP",'11' ="VP", '24'="Prog/UrLE",
                         '29' ="OE",'6' ="Telocyte",
                         '0' ="SubEpithelial FB", '19' ="SubEpithelial FB",'8' ="Intestitial FB",
                         '12' ="Intestitial FB",'23' ="Urethral FB",'20' ="SM",
                         '18' ="Glia", '25'="Pericyte", '27'="Pericyte",
                         '16'="VE", '30'="VE",
                         '14' = "Monocyte", '22' = "Monocyte", '26' = "Mast Cell",'15' = "T cell")
WT_P60_1[["CellTypes1"]] <- Idents(object = WT_P60_1)

#Umap
Idents(object = WT_P60_1) <- "CellTypes1"
tiff(file = "WT_P60_1 UMAP CellTypes1.tiff", width = 10, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "WT_P60_1 UMAP CellTypes1 label.tiff", width = 12, height = 8, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#DEGs_CellTypes1
Idents(object = WT_P60_1) <- "CellTypes1"
DefaultAssay(WT_P60_1) <- "RNA"
all.genes <- rownames(WT_P60_1)
WT_P60_1 <- ScaleData(WT_P60_1, features = all.genes)
WT_P60_1.CellTypes1.markers <- FindAllMarkers(WT_P60_1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(WT_P60_1.CellTypes1.markers, file = "WT_P60_1.CellTypes1.markers.csv")

#Dotplot
Idents(object = WT_P60_1) <- "CellTypes1"
tiff(file = "WT_P60_1 CellTypes1 markers DotPlot.tiff", width =20 , height = 6, units = "in", compression = "lzw", res = 200)
DotPlot(WT_P60_1, features = c("Krt14", "Krt17", "Lgals7", "Lcn2", "Plet1",
                               "Ren1", "Glb1l3", "Fgl1", "C1s2", "Mt3", "Fut9", "Atp12a",
                               "Tcaf3", "Nlrp10", "Glb1l2", "Rimklb", "Slc26a3", "Greb1",
                               "Pcp4", "Aqp8", "Upb1", "Ceacam2", "Defb10", "Spink11", "Spc25",
                               "5430419D17Rik", "Abo", "Lrrc26", "Fam25c", "Timp4", "Csta1", "Kcnk3", "Slc26a4", "Pdzk1", "Gucy2g", "Atp2c2", "Trpv6",
                               "Pglyrp1", "Gsdmc2", "Ppp1r1b", "Gsdmc3", "Acsm3", "Krt4", "Psca", "Wfdc2",
                               "Cdca8", "Hmmr", "Kif4", "Nuf2", "Racgap1", "Tpx2",
                               "Mfap4", "Pcolce", "Hsd11b1", "Itgbl1", "Fxyd6", "Rcn3", "Cpxm2", "Edn3", 
                               "Pdgfrl", "Ogn", "Fbln7", "Pdgfra", "Lamc3", "Podn",
                               "Col6a5", "Mdk", "Tmem204", "Scarf2", "Dner", "Ctsk", "Dkk2", "Mfap2", "Calml4", "Atp1a2", "Twist1",
                               "Wif1", "Rorb", "Gpc1", "Wnt2", "Nptx2", "Sfrp2", "Osr2", 
                               "Wfdc1", "Wt1", 
                               "Sult1e1", "C3", "Dpt", "Thbd", "Clec3b", "Lum", "Cilp",
                               "Lgr5", "Pi15", "Aqp1", "Osr1", "C7", "Pdlim3", "Sod3", "Spock2", "Ltbp1", 
                               "Cnn1", "Myh11", "Actg2", "Acta1", "Tnnc2", "Lmod1",
                               "Plp1", "Fabp7", "Kcna1", "Cdh19", "S100b",
                               "Rgs5", "Gja4", "Ndufa4l2", "Rgs4", "Vtn", "Cox4i2",
                               "Flt1", "Pecam1", "Egfl7", "Cdh5", "Cd93", "Plvap", "Ptprb", "Apold1", "Adgrl4", "Sox18",
                               "Ccl4", "H2-Aa", "Lyz2", "C1qb", "C1qa",
                               "Ccr7", "Ccl22", "Timd4", "Pkib", "Plbd1", "Napsa", "Myo1g",
                               "Rgs1", "Cd3g", "Nkg7", "Laptm5", "Sh2d2a", "Cd53", "Ptprc", "Cytip"),
        cols = c("light grey", "red")) + RotatedAxis()
dev.off()

Idents(object = WT_P60_1) <- "CellTypes1"
tiff(file = "WT_P60_1 CellTypes1 markers DotPlot.tiff", width =20 , height = 6, units = "in", compression = "lzw", res = 200)
DotPlot(WT_P60_1, features = c("Krt14", "Krt17", "Lgals7", "Lcn2", "Plet1",
                               "Ren1", "Glb1l3", "C1s2", "Mt3", "Fut9",
                               "Tcaf3", "Nlrp10", "Glb1l2", "Rimklb", "Slc26a3", 
                               "Aqp8", "Upb1", "Ceacam2", "Defb10", "Spink11", 
                               "Fam25c", "Csta1", "Slc26a4", "Gucy2g", "Trpv6",
                               "Ppp1r1b", "Acsm3", "Krt4", "Psca", "Wfdc2",
                               "Cdca8", "Kif4", "Nuf2", "Racgap1", "Tpx2",
                               "Hsd11b1", "Itgbl1", "Fxyd6", "Edn3", "Podn",
                               "Rorb", "Gpc1", "Wnt2", "Wt1", "F2r", 
                               "Sult1e1", "C3", "Cilp", "Ifi205", "Dpep1", 
                               "Lgr5", "Pi15", "Osr1", "Sod3", "Ltbp1", 
                               "Cnn1", "Actg2", "Acta1", "Tnnc2", "Lmod1",
                               "Plp1", "Fabp7", "Kcna1", "Cdh19", "S100b",
                               "Rgs5", "Gja4", "Ndufa4l2", "Rgs4", "Cox4i2",
                               "Pecam1", "Egfl7", "Cdh5", "Plvap", "Ptprb", 
                               "Ccl4", "C1qb", "C1qa", "Ccl3", "C1qc",
                               "Ccr7", "Ccl22", "Timd4", "Pkib", "Napsa",
                               "Cd3g", "Nkg7", "Sh2d2a", "Il2rb", "Ptpn22"),
        cols = c("light grey", "red")) + RotatedAxis()
dev.off()

##telocyte markers
#Featureplots
DefaultAssay(WT_P60_1)<-"RNA"
tiff(file = "WT_P60_1 Pdgfra expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "WT_P60_1 Cd34 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "WT_P60_1 Foxl1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "WT_P60_1 Kit expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_1, reduction = "umap", features = c("Kit"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "WT_P60_1 Pdgfra Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_1, features = c("Pdgfra"), pt.size = 0)
dev.off()
tiff(file = "WT_P60_1 Cd34 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_1, features = c("Cd34"), pt.size = 0)
dev.off()
tiff(file = "WT_P60_1 Kit Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_1, features = c("Kit"), pt.size = 0)
dev.off()
tiff(file = "WT_P60_1 Foxl1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_1, features = c("Foxl1"), pt.size = 0)
dev.off()

####Subset FBSM####
Idents(object = WT_P60_1) <- "CellTypes"
WT_P60_FBSM <- subset(WT_P60_1, idents = c("SubEpithelial FB1",
                                           "SubEpithelial FB2",
                                           "Telocyte", "Intestitial FB1",
                                           "Intestitial FB2", "Urethral FB",
                                           "SM"))

DefaultAssay(WT_P60_FBSM) <- "RNA"
#Run the standard workflow for visualization and clustering
WT_P60_FBSM <- ScaleData(WT_P60_FBSM, verbose = FALSE)
WT_P60_FBSM <- RunPCA(WT_P60_FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(WT_P60_FBSM, ndims = 50)
# UMAP and Clustering
WT_P60_FBSM <- FindNeighbors(WT_P60_FBSM, reduction = "pca", dims = 1:23)
WT_P60_FBSM <- FindClusters(WT_P60_FBSM, resolution = 0.8)
WT_P60_FBSM <- RunUMAP(WT_P60_FBSM, reduction = "pca", dims = 1:23)

Idents(object = WT_P60_FBSM) <- "seurat_clusters"
DimPlot(WT_P60_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
DimPlot(WT_P60_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim") 

#Cell Cycle scoring
DefaultAssay(WT_P60_FBSM) <- "RNA"
all.genes <- rownames(WT_P60_FBSM)
WT_P60_FBSM <- ScaleData(WT_P60_FBSM, features = all.genes)
WT_P60_FBSM <- CellCycleScoring(WT_P60_FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = WT_P60_FBSM) <- "Phase"
WT_P60_FBSM <- RenameIdents(object = WT_P60_FBSM, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
WT_P60_FBSM[["Phase"]] <- Idents(object = WT_P60_FBSM)
DimPlot(WT_P60_FBSM, reduction = "umap")

tiff(file = "WT_P60_FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_FBSM, reduction = "umap", pt.size = 0.5, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#Cell Cycle regression
WT_P60_FBSM1 <- WT_P60_FBSM
DefaultAssay(WT_P60_FBSM1) <- "RNA"
WT_P60_FBSM1 <- ScaleData(WT_P60_FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(WT_P60_FBSM1))
WT_P60_FBSM1 <- RunPCA(WT_P60_FBSM1, features = VariableFeatures(WT_P60_FBSM1))
ElbowPlot(WT_P60_FBSM1, ndims = 50)

WT_P60_FBSM1 <- FindNeighbors(WT_P60_FBSM1, reduction = "pca", dims = 1:24)
WT_P60_FBSM1 <- FindClusters(WT_P60_FBSM1, resolution = 0.8)
WT_P60_FBSM1 <- RunUMAP(WT_P60_FBSM1, reduction = "pca", dims = 1:24)
DimPlot(WT_P60_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = WT_P60_FBSM1) <- "seurat_clusters"
DimPlot(WT_P60_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "WT_P60_FBSM1 seurat_clusters UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_FBSM1, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()


Idents(object = WT_P60_FBSM1) <- "Phase"
WT_P60_FBSM1 <- RenameIdents(object = WT_P60_FBSM1, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
WT_P60_FBSM1[["Phase"]] <- Idents(object = WT_P60_FBSM1)
DimPlot(WT_P60_FBSM1, reduction = "umap")

Idents(object = WT_P60_FBSM1) <- "Phase"
tiff(file = "WT_P60_FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_FBSM1, reduction = "umap", pt.size = 0.5, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#DEGs
Idents(object = WT_P60_FBSM1) <- "seurat_clusters"
DefaultAssay(WT_P60_FBSM1) <- "RNA"
all.genes <- rownames(WT_P60_FBSM1)
WT_P60_FBSM1 <- ScaleData(WT_P60_FBSM1, features = all.genes)
WT_P60_FBSM1.seurat_cluster.res0.8.markers <- FindAllMarkers(WT_P60_FBSM1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(WT_P60_FBSM1.seurat_cluster.res0.8.markers, file = "WT_P60_FBSM1.seurat_cluster.res0.8.markers.csv")

#Dotplot
Idents(object = WT_P60_FBSM1) <- "seurat_clusters"
tiff(file = "WT_P60_FBSM1 CellTypes stromal markers DotPlot.tiff", width =15 , height = 4, units = "in", compression = "lzw", res = 200)
DotPlot(WT_P60_FBSM1, features = c("Rorb", "F2r", "Adm", "Rdh10", "Fzd1",
                               "Wnt2", "Wif1", "Ifitm1", "Srd5a2",
                               "Sult1e1", "Lum", "Ogn", "Dpt", "C3",
                               "Ctgf", "Clec3b", "Dcn", "Ebf1", "Gpx3", "Igf1",
                               "Lgr5", "Apoe", "Osr1", "Sfrp2", "Mfap4",
                               "Foxf1", "Foxl1", "Pdgfra", "Cd34", "Kit", "Gli1", 
                               "Ar", "Fbln1", "Acta2", "Tagln", "Myh11"),
        cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#FeaturePlot
DefaultAssay(WT_P60_FBSM1)<-"RNA"
tiff(file = "WT_P60_FBSM1 celltype marker expression plots.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_FBSM1, reduction = "umap", features = c(
                                                       "Fbln1", "Fbn1",
                                                       "Acta2", "Myh11",
                                                       "Pdgfra", "Foxl1", "Cd34", "Kit",
                                                       "Rorb", "Wif1",
                                                       "Sult1e1", "Lum",
                                                       "Lgr5", "Osr1"
                                                       ), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = WT_P60_FBSM1) <- "seurat_clusters"
WT_P60_FBSM1 <- RenameIdents(object = WT_P60_FBSM1, 
                         '2' = "SubEpithelial FB1", '0' ="SubEpithelial FB2", 
                         '7' ="SubEpithelial FB3",'5' ="Telocyte",
                         '3' ="MyoFB",'4' ="Intestitial FB1", '1'= "Intestitial FB2",
                         '9' ="VP/Proximal FB", '6'="SM", '8'="OS")
WT_P60_FBSM1[["FBSMCellTypes"]] <- Idents(object = WT_P60_FBSM1)

#Umap
Idents(object = WT_P60_FBSM1) <- "FBSMCellTypes"
tiff(file = "WT_P60_FBSM1 UMAP FBSMCellTypes.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_FBSM1, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "WT_P60_FBSM1 UMAP FBSMCellTypes label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(WT_P60_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#DEGs
Idents(object = WT_P60_FBSM1) <- "FBSMCellTypes"
DefaultAssay(WT_P60_FBSM1) <- "RNA"
all.genes <- rownames(WT_P60_FBSM1)
WT_P60_FBSM1 <- ScaleData(WT_P60_FBSM1, features = all.genes)
WT_P60_FBSM1.FBSMCellTypes.res0.8.markers <- FindAllMarkers(WT_P60_FBSM1, only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.1)
write.csv(WT_P60_FBSM1.FBSMCellTypes.res0.8.markers, file = "WT_P60_FBSM1.FBSMCellTypes.res0.8.markers.csv")

#Featureplot
tiff(file = "WT_P60_FBSM1 Pdgfra expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_FBSM1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q80")
dev.off()
tiff(file = "WT_P60_FBSM1 Cd34 expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_FBSM1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "WT_P60_FBSM1 Foxl1 expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_FBSM1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q80")
dev.off()
tiff(file = "WT_P60_FBSM1 Kit expression plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
FeaturePlot(WT_P60_FBSM1, reduction = "umap", features = c("Kit"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

#Vlnplot
tiff(file = "WT_P60_FBSM1 Pdgfra Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_FBSM1, features = c("Pdgfra"), pt.size = 0)
dev.off()
tiff(file = "WT_P60_FBSM1 Cd34 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_FBSM1, features = c("Cd34"), pt.size = 0)
dev.off()
tiff(file = "WT_P60_FBSM1 Foxl1 Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_FBSM1, features = c("Foxl1"), pt.size = 0)
dev.off()
tiff(file = "WT_P60_FBSM1 Kit Vln.tiff", width = 10, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(WT_P60_FBSM1, features = c("Kit"), pt.size = 0)
dev.off()