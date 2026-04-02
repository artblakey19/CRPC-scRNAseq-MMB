####Slingshot####

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.6/pseudotime/")

list.files()
ARQPosEpi <- readRDS("ARQPosEpi.rds")
ARQNegEpi <- readRDS("ARQNegEpi.rds")

ARQPosEpi1 <- as.SingleCellExperiment(ARQPosEpi)
ARQNegEpi1 <- as.SingleCellExperiment(ARQNegEpi)

####ARQPosEpi####
DefaultAssay(ARQPosEpi) <- "RNA"
Idents(object = ARQPosEpi) <- "seurat_clusters"
DimPlot(ARQPosEpi, reduction = "umap")

ARQPosEpi1 <- ARQPosEpi

DefaultAssay(ARQPosEpi1) <- "integrated"

#Run the standard workflow for visualization and clustering
ARQPosEpi1 <- ScaleData(ARQPosEpi1, verbose = FALSE)
ARQPosEpi1 <- RunPCA(ARQPosEpi1, npcs = 30, verbose = FALSE)
ElbowPlot(ARQPosEpi1)
#Umap and Clustering
ARQPosEpi1 <- RunTSNE(ARQPosEpi1, reduction = "pca", dims = 1:20)
ARQPosEpi1 <- RunUMAP(ARQPosEpi1, reduction = "pca", dims = 1:20)
DimPlot(ARQPosEpi1, reduction = "umap")

Idents(object = ARQPosEpi1) <- "EpiCellType"
DimPlot(ARQPosEpi1, reduction = "umap")

ARQPosEpi1_sds <- slingshot(Embeddings(ARQPosEpi1, "umap"), clusterLabels = ARQPosEpi1$EpiCellType, 
                           start.clus = "BE1", end.clus = "LE1", stretch = 2)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_clust <- cell_pal(ARQPosEpi$EpiCellType, hue_pal()) 

plot(reducedDim(ARQPosEpi1_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQPosEpi1_sds, lwd = 2, col = 'black')

tiff(file = "ARQPosEpi slingshot Pseudotime-final.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot(reducedDim(ARQPosEpi1_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQPosEpi1_sds, lwd = 2, col = 'black')
dev.off()

#
nc <- 3
pt <- slingPseudotime(ARQPosEpi1_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)

tiff(file = "ARQPosEpi slingshot Pseudotime split-final.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQPosEpi1_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQPosEpi1_sds, lwd = 2, col = 'black')
  
}
dev.off()

tiff(file = "ARQPosEpi1 slingshot UMAP label.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQPosEpi1, reduction = "umap", pt.size = 1, label = TRUE)
dev.off()

####ARQNegEpi####

ARQNegEpi1 <- ARQNegEpi

DefaultAssay(ARQNegEpi1) <- "integrated"

#Run the standard workflow for visualization and clustering
ARQNegEpi1 <- ScaleData(ARQNegEpi1, verbose = FALSE)
ARQNegEpi1 <- RunPCA(ARQNegEpi1, npcs = 30, verbose = FALSE)
ElbowPlot(ARQNegEpi1)
#Umap and Clustering
ARQNegEpi1 <- FindNeighbors(ARQNegEpi1, reduction = "pca", dims = 1:5)
ARQNegEpi1 <- FindClusters(ARQNegEpi1, resolution = 0.3)
ARQNegEpi1 <- RunTSNE(ARQNegEpi1, reduction = "pca", dims = 1:20)
ARQNegEpi1 <- RunUMAP(ARQNegEpi1, reduction = "pca", dims = 1:20)
DimPlot(ARQNegEpi1, reduction = "umap")

Idents(object = ARQNegEpi1) <- "EpiCellType"
DimPlot(ARQNegEpi1, reduction = "umap")

Idents(object = ARQNegEpi1) <- "seurat_clusters"
DimPlot(ARQNegEpi1, reduction = "umap", label = TRUE)

ARQNegEpi1_sds <- slingshot(Embeddings(ARQNegEpi1, "umap"), clusterLabels = ARQNegEpi1$seurat_clusters, 
                            start.clus = "6")

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_clust <- cell_pal(ARQNegEpi1$EpiCellType, hue_pal()) 

plot(reducedDim(ARQNegEpi1_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQNegEpi1_sds, lwd = 2, col = 'black')

tiff(file = "ARQNegEpi slingshot Pseudotime-final.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot(reducedDim(ARQNegEpi1_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQNegEpi1_sds, lwd = 2, col = 'black')
dev.off()

#
nc <- 3
pt <- slingPseudotime(ARQNegEpi1_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)

tiff(file = "ARQNegEpi slingshot Pseudotime split-final.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQNegEpi1_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQNegEpi1_sds, lwd = 2, col = 'black')
  
}
dev.off()

tiff(file = "ARQNegEpi1 slingshot UMAP label.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQNegEpi1, reduction = "umap", pt.size = 1, label = TRUE)
dev.off()