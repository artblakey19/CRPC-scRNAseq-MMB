#Human scRNAseq

#Add necessary tools to library
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(scales)
library(viridis)
library(scater)
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)
library(irlba)

####Load data GSE157703 PCa1####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE157703/GSE157703_RAW")

GSE157703_PCa1_raw <- read.table(file = "GSM4773521_PCa1_gene_counts_matrix.csv", sep = ",", header = T, row.names=1, as.is=T)
GSE157703_PCa1 <- CreateSeuratObject(counts = GSE157703_PCa1_raw, min.cells = 3, min.features = 500)
GSE157703_PCa1 <- NormalizeData(GSE157703_PCa1)
GSE157703_PCa1 <- FindVariableFeatures(GSE157703_PCa1, selection.method = "vst", nfeatures = 2000)

#Run the standard workflow for visualization and clustering Ctrl2
GSE157703_PCa1 <- ScaleData(GSE157703_PCa1, verbose = FALSE)
GSE157703_PCa1 <- RunPCA(GSE157703_PCa1, npcs = 30, verbose = FALSE)
ElbowPlot(GSE157703_PCa1, ndims = 30)
# UMAP and Clustering
GSE157703_PCa1 <- FindNeighbors(GSE157703_PCa1, reduction = "pca", dims = 1:15)
GSE157703_PCa1 <- FindClusters(GSE157703_PCa1, resolution = 0.7)
GSE157703_PCa1 <- RunUMAP(GSE157703_PCa1, reduction = "pca", dims = 1:15)
DimPlot(GSE157703_PCa1, reduction = "umap", pt.size = 0.5, label = TRUE) 

#DEGs
DefaultAssay(GSE157703_PCa1)<-"RNA"
Idents(object = GSE157703_PCa1) <- "seurat_clusters"
all.genes <- rownames(GSE157703_PCa1)
GSE157703_PCa1 <- ScaleData(GSE157703_PCa1, features = all.genes)
GSE157703_PCa1.seurat.markers <- FindAllMarkers(GSE157703_PCa1, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(GSE157703_PCa1.seurat.markers, file = "GSE157703_PCa1.seurat.markers.csv")

tiff(file = "GSE157703_PCa1 celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa1, reduction = "umap", features = c("AR", "EPCAM", "VIM",
                                                             "KRT5", "KRT8", "KRT18",
                                                             "FBLN1", "ACTA2", "DES",
                                                             "PECAM1", "TYROBP", "CCL5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSE157703_PCa1) <- "seurat_clusters"
DimPlot(GSE157703_PCa1, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE157703_PCa1 <- RenameIdents(object = GSE157703_PCa1, 
                                '11' = "BE", '5' = "LE", '3' = "LE",
                               '10' = "FB", '4' = "MyoFB", '7' = "SM",
                               '1' = "Endo", '6' = "Endo", 
                               '9' = "Immune", '8' = "Immune", '0' = "Immune", '2' = "Immune")  
GSE157703_PCa1[["CellTypes"]] <- Idents(object = GSE157703_PCa1)

#UMAP_CellTypes
Idents(object = GSE157703_PCa1) <- "CellTypes"
tiff(file = "GSE157703_PCa1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE157703_PCa1, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

#Cell count
Idents(object = GSE157703_PCa1) <- "stim"
table(Idents(GSE157703_PCa1))

#FeaturePlots
DefaultAssay(GSE157703_PCa1) <- "RNA"
tiff(file = "GSE157703_PCa1 GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa1, reduction = "umap", features = c("GLI1"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "GSE157703_PCa1 SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa1, reduction = "umap", features = c("SHH"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "GSE157703_PCa1 AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa1, reduction = "umap", features = c("AR"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

####Load data GSE157703 PCa2####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE157703/GSE157703_RAW")

GSE157703_PCa2_raw <- read.table(file = "GSM4773522_PCa2_gene_counts_matrix.txt", sep = "", header = T, row.names=1, as.is=T)
GSE157703_PCa2 <- CreateSeuratObject(counts = GSE157703_PCa2_raw, min.cells = 3, min.features = 500)
GSE157703_PCa2 <- NormalizeData(GSE157703_PCa2)
GSE157703_PCa2 <- FindVariableFeatures(GSE157703_PCa2, selection.method = "vst", nfeatures = 2000)

#Run the standard workflow for visualization and clustering Ctrl2
GSE157703_PCa2 <- ScaleData(GSE157703_PCa2, verbose = FALSE)
GSE157703_PCa2 <- RunPCA(GSE157703_PCa2, npcs = 30, verbose = FALSE)
ElbowPlot(GSE157703_PCa2, ndims = 30)
# UMAP and Clustering
GSE157703_PCa2 <- FindNeighbors(GSE157703_PCa2, reduction = "pca", dims = 1:15)
GSE157703_PCa2 <- FindClusters(GSE157703_PCa2, resolution = 0.7)
GSE157703_PCa2 <- RunUMAP(GSE157703_PCa2, reduction = "pca", dims = 1:15)
DimPlot(GSE157703_PCa2, reduction = "umap", pt.size = 0.5, label = TRUE) 

#DEGs
DefaultAssay(GSE157703_PCa2)<-"RNA"
Idents(object = GSE157703_PCa2) <- "seurat_clusters"
all.genes <- rownames(GSE157703_PCa2)
GSE157703_PCa2 <- ScaleData(GSE157703_PCa2, features = all.genes)
GSE157703_PCa2.seurat.markers <- FindAllMarkers(GSE157703_PCa2, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(GSE157703_PCa2.seurat.markers, file = "GSE157703_PCa2.seurat.markers.csv")

tiff(file = "GSE157703_PCa2 celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa2, reduction = "umap", features = c("AR", "EPCAM", "VIM",
                                                             "KRT5", "KRT8", "PLP1",
                                                             "FBLN1", "ACTA2", "DES",
                                                             "PECAM1", "TYROBP", "CCL5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSE157703_PCa2) <- "seurat_clusters"
DimPlot(GSE157703_PCa2, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE157703_PCa2 <- RenameIdents(object = GSE157703_PCa2, 
                               '11' = "BE", '2' = "LE", '3' = "LE", '6' = "LE", '10' = "LE",
                               '15' = "FB", '4' = "MyoFB", '9' = "SM",
                               '5' = "Endo", '7' = "Endo", '8' = "Endo",
                               '12' = "Immune", '13' = "Immune", '1' = "Immune", '0' = "Immune",
                               '14' = "Immune", '16' = "Glia")  
GSE157703_PCa2[["CellTypes"]] <- Idents(object = GSE157703_PCa2)

#UMAP_CellTypes
Idents(object = GSE157703_PCa2) <- "CellTypes"
tiff(file = "GSE157703_PCa2 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE157703_PCa2, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

#Cell count
Idents(object = GSE157703_PCa2) <- "stim"
table(Idents(GSE157703_PCa2))

#FeaturePlots
DefaultAssay(GSE157703_PCa2) <- "RNA"
tiff(file = "GSE157703_PCa2 GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa2, reduction = "umap", features = c("GLI1"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "GSE157703_PCa2 SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa2, reduction = "umap", features = c("SHH"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "GSE157703_PCa2 AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE157703_PCa2, reduction = "umap", features = c("AR"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

####Load data GSK141445####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE141445/GSE141445_RAW")

GSE141445_PCa_raw <- read.table(file = "data.raw.matrix.txt", sep = "", header = T, row.names=1, as.is=T)
GSE141445_PCa <- CreateSeuratObject(counts = GSE141445_PCa_raw,  min.cells = 3, min.features = 500)
GSE141445_PCa <- NormalizeData(GSE141445_PCa)
GSE141445_PCa <- FindVariableFeatures(GSE141445_PCa, selection.method = "vst", nfeatures = 5000)
#Run the standard workflow for visualization and clustering Ctrl2
GSE141445_PCa <- ScaleData(GSE141445_PCa, verbose = FALSE)
GSE141445_PCa <- RunPCA(GSE141445_PCa, npcs = 30, verbose = FALSE)
ElbowPlot(GSE141445_PCa, ndims = 50)
# UMAP and Clustering
GSE141445_PCa <- FindNeighbors(GSE141445_PCa, reduction = "pca", dims = 1:20)
GSE141445_PCa <- FindClusters(GSE141445_PCa, resolution = 0.5)
GSE141445_PCa <- RunUMAP(GSE141445_PCa, reduction = "pca", dims = 1:20)
DimPlot(GSE141445_PCa, reduction = "umap", pt.size = 0.5, label = TRUE) 

#DEGs
DefaultAssay(GSE141445_PCa)<-"RNA"
Idents(object = GSE141445_PCa) <- "seurat_clusters"
all.genes <- rownames(GSE141445_PCa)
GSE141445_PCa <- ScaleData(GSE141445_PCa, features = all.genes)
GSE141445_PCa.seurat.markers <- FindAllMarkers(GSE141445_PCa, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(GSE141445_PCa.seurat.markers, file = "GSE141445_PCa.seurat.markers.csv")

tiff(file = "GSE141445_PCa celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE141445_PCa, reduction = "umap", features = c("AR", "EPCAM", "VIM",
                                                             "KRT5", "KRT8", "CHGB",
                                                             "FBLN1", "ACTA2", "DES",
                                                             "PECAM1", "TYROBP", "CCL5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSE141445_PCa) <- "seurat_clusters"
DimPlot(GSE141445_PCa, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE141445_PCa <- RenameIdents(object = GSE141445_PCa, 
                               '16' = "BE", '10' = "LE", '5' = "LE", '2' = "LE", '27' = "LE",
                              '3' = "LE", '0' = "LE", '17' = "LE", '22' = "LE",
                              '15' = "LE", '21' = "LE", '19' = "LE", '12' = "LE",
                              '4' = "LE", '8' = "LE", '7' = "LE", '13' = "LE",'18' = "LE",
                               '25' = "FB", '8' = "SM",
                               '9' = "Endo", '14' = "Endo", '23' = "Endo", '24' = "Endo",
                               '11' = "Immune", '6' = "Immune", '1' = "Immune", '26' = "Immune",
                               '28' = "Immune", '20' = "Immune")  
GSE141445_PCa[["CellTypes"]] <- Idents(object = GSE141445_PCa)

#UMAP_CellTypes
Idents(object = GSE141445_PCa) <- "CellTypes"
tiff(file = "GSE141445_PCa CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE141445_PCa, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

#Cell count
Idents(object = GSE141445_PCa) <- "stim"
table(Idents(GSE141445_PCa))

#FeaturePlots
DefaultAssay(GSE141445_PCa) <- "RNA"
tiff(file = "GSE141445_PCa GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE141445_PCa, reduction = "umap", features = c("GLI1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE141445_PCa SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE141445_PCa, reduction = "umap", features = c("SHH"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE141445_PCa AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE141445_PCa, reduction = "umap", features = c("AR"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()



####Load data GSE117403_D17####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE120716_Normal/GSE117403_D17")
GSE117403_D17.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE120716_Normal/GSE117403_D17")
GSE117403_D17 <- CreateSeuratObject(counts = GSE117403_D17.data, min.cells = 3, min.features = 500)
GSE117403_D17 <- NormalizeData(GSE117403_D17)
GSE117403_D17 <- FindVariableFeatures(GSE117403_D17, selection.method = "vst", nfeatures = 2000)

#Run the standard workflow for visualization and clustering Ctrl2
GSE117403_D17 <- ScaleData(GSE117403_D17, verbose = FALSE)
GSE117403_D17 <- RunPCA(GSE117403_D17, npcs = 30, verbose = FALSE)
ElbowPlot(GSE117403_D17, ndims = 30)
# UMAP and Clustering
GSE117403_D17 <- FindNeighbors(GSE117403_D17, reduction = "pca", dims = 1:15)
GSE117403_D17 <- FindClusters(GSE117403_D17, resolution = 0.7)
GSE117403_D17 <- RunUMAP(GSE117403_D17, reduction = "pca", dims = 1:15)
DimPlot(GSE117403_D17, reduction = "umap", pt.size = 0.5, label = TRUE) 

DefaultAssay(GSE117403_D17) <- "RNA"
tiff(file = "GSE117403_D17 celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("AR", "EPCAM", "VIM",
                                                             "KRT5", "KRT8", "CHGA",
                                                             "FBLN1", "ACTA2", "DES",
                                                             "CLDN5", "RGS1", "CCL5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSE117403_D17) <- "seurat_clusters"
DimPlot(GSE117403_D17, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE117403_D17 <- RenameIdents(object = GSE117403_D17, 
                                '0' = "BE",'1' = "BE",'11' = "BE",
                              '5' = "OE", '2' = "OE", '6' = "OE",'15' = "OE", '3'="LE", '10'="LE",
                               '8' = "FB", '4' = "FB",  '9' = "SM", '12' = "SM",
                               '7' = "Endo", '13' = "Endo", '14'= "Endo"
                               )  
GSE117403_D17[["CellTypes"]] <- Idents(object = GSE117403_D17)

#UMAP_CellTypes
Idents(object = GSE117403_D17) <- "CellTypes"
tiff(file = "GSE117403_D17 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D17, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

Idents(object = GSE117403_D17) <- "CellTypes"
tiff(file = "GSE117403_D17 CellTypes UMAP without label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D17, reduction = "umap", pt.size = 0.5)
dev.off()

#Cell count
Idents(object = GSE117403_D17) <- "stim"
table(Idents(GSE117403_D17))

#FeaturePlots
DefaultAssay(GSE117403_D17) <- "RNA"
tiff(file = "GSE117403_D17 GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("GLI1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17 SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("SHH"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17 AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("AR"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17 VIM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("VIM"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17 CDH1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("CDH1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17 FBLN1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17, reduction = "umap", features = c("FBLN1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Stroma
Idents(object = GSE117403_D17) <- "CellTypes"
tiff(file = "GSE117403_D17 FBSM Highlighted UMAP without label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D17, reduction = "umap", pt.size = 0.5, cols = c("light grey", "light grey","light grey",
                                                                   "darkblue", "darkblue", "light grey"))
dev.off()

GSE117403_D17_FBSM <- subset(GSE117403_D17, idents = c("FB", "SM"))

#Run the standard workflow for visualization and clustering
GSE117403_D17_FBSM <- ScaleData(GSE117403_D17_FBSM, verbose = FALSE)
GSE117403_D17_FBSM <- RunPCA(GSE117403_D17_FBSM, npcs = 50, verbose = FALSE)
ElbowPlot(GSE117403_D17_FBSM, ndims = 30)

#Umap and Clustering
GSE117403_D17_FBSM <- FindNeighbors(GSE117403_D17_FBSM, reduction = "pca", dims = 1:17)
GSE117403_D17_FBSM <- FindClusters(GSE117403_D17_FBSM, resolution = 0.5)
GSE117403_D17_FBSM <- RunUMAP(GSE117403_D17_FBSM, reduction = "pca", dims = 1:17)
DimPlot(GSE117403_D17_FBSM, reduction = "umap", pt.size = 0.5, label = TRUE)

Idents(object = GSE117403_D17_FBSM) <- "seurat_clusters"
GSE117403_D17_FBSM <- subset(GSE117403_D17_FBSM, idents = c("0", "1", "2", "3", "4",
                                                            "5", "6", "7", "8", "10"))

GSE117403_D17_FBSM <- FindNeighbors(GSE117403_D17_FBSM, reduction = "pca", dims = 1:17)
GSE117403_D17_FBSM <- FindClusters(GSE117403_D17_FBSM, resolution = 0.5)
GSE117403_D17_FBSM <- RunUMAP(GSE117403_D17_FBSM, reduction = "pca", dims = 1:17)
DimPlot(GSE117403_D17_FBSM, reduction = "umap", pt.size = 0.5, label = TRUE)

Idents(object = GSE117403_D17_FBSM) <- "seurat_clusters"
GSE117403_D17_FBSM <- subset(GSE117403_D17_FBSM, idents = c("0", "1", "2", "3", "4",
                                                            "5", "6", "7"))
GSE117403_D17_FBSM <- FindNeighbors(GSE117403_D17_FBSM, reduction = "pca", dims = 1:15)
GSE117403_D17_FBSM <- FindClusters(GSE117403_D17_FBSM, resolution = 0.5)
GSE117403_D17_FBSM <- RunUMAP(GSE117403_D17_FBSM, reduction = "pca", dims = 1:15)
DimPlot(GSE117403_D17_FBSM, reduction = "umap", pt.size = 0.5, label = TRUE)

#Rename CellTypes
Idents(object = GSE117403_D17_FBSM) <- "seurat_clusters"
DimPlot(GSE117403_D17_FBSM, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE117403_D17_FBSM <- RenameIdents(object = GSE117403_D17_FBSM, 
                              '1' = "FB1",'3' = "FB2",'5' = "FB3",
                              '0' = "FB4", '7' = "FB5", '8' = "FB6", '2' = "SM1",'4' = "SM2", '6'="SM3"
                            )  
GSE117403_D17_FBSM[["FBSMCellTypes"]] <- Idents(object = GSE117403_D17_FBSM)


#UMAP_CellTypes
Idents(object = GSE117403_D17_FBSM) <- "FBSMCellTypes"
tiff(file = "GSE117403_D17_FBSM FBSMCellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D17_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#UMAP_CellTypes
Idents(object = GSE117403_D17_FBSM) <- "FBSMCellTypes"
tiff(file = "GSE117403_D17_FBSM FBSMCellTypes UMAP without label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D17_FBSM, reduction = "umap", pt.size = 0.3)
dev.off()

#FeaturePlots
DefaultAssay(GSE117403_D17_FBSM) <- "RNA"
tiff(file = "GSE117403_D17_FBSM GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17_FBSM, reduction = "umap", features = c("GLI1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17_FBSM SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17_FBSM, reduction = "umap", features = c("SHH"), cols = c("light grey", "light grey"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17_FBSM AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17_FBSM, reduction = "umap", features = c("AR"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D17_FBSM VIM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D17_FBSM, reduction = "umap", features = c("VIM"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####Load data GSE117403_D27####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE120716_Normal/GSE117403_D27")
GSE117403_D27.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE120716_Normal/GSE117403_D27")
GSE117403_D27 <- CreateSeuratObject(counts = GSE117403_D27.data, min.cells = 3, min.features = 500)
GSE117403_D27 <- NormalizeData(GSE117403_D27)
GSE117403_D27 <- FindVariableFeatures(GSE117403_D27, selection.method = "vst", nfeatures = 2000)

#Run the standard workflow for visualization and clustering Ctrl2
GSE117403_D27 <- ScaleData(GSE117403_D27, verbose = FALSE)
GSE117403_D27 <- RunPCA(GSE117403_D27, npcs = 30, verbose = FALSE)
ElbowPlot(GSE117403_D27, ndims = 30)
# UMAP and Clustering
GSE117403_D27 <- FindNeighbors(GSE117403_D27, reduction = "pca", dims = 1:15)
GSE117403_D27 <- FindClusters(GSE117403_D27, resolution = 0.7)
GSE117403_D27 <- RunUMAP(GSE117403_D27, reduction = "pca", dims = 1:15)
DimPlot(GSE117403_D27, reduction = "umap", pt.size = 0.5, label = TRUE) 

DefaultAssay(GSE117403_D27) <- "RNA"
tiff(file = "GSE117403_D27 celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D27, reduction = "umap", features = c("AR", "EPCAM", "VIM",
                                                            "KRT5", "KLK3", "PIGR", "KRT13",
                                                            "CALCA",
                                                            "FBLN1", "ACTA2", 
                                                            "CLDN5", "RGS1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSE117403_D27) <- "seurat_clusters"
DimPlot(GSE117403_D27, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE117403_D27 <- RenameIdents(object = GSE117403_D27, 
                              '0' = "BE",'2' = "BE",'6' = "BE", '13' = "BE",
                              '11' = "OE", '8' = "OE", '14' = "OE", '3' = "OE",
                              '7'="LE", '10'="LE",
                              '15' = "FB", '4' = "FB", '9' = "FB", "5" = "FB",
                              '1' = "SM", '12' = "SM",
                              '16' = "Endo"
)  
GSE117403_D27[["CellTypes"]] <- Idents(object = GSE117403_D27)

#UMAP_CellTypes
Idents(object = GSE117403_D27) <- "CellTypes"
tiff(file = "GSE117403_D27 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D27, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

Idents(object = GSE117403_D27) <- "CellTypes"
tiff(file = "GSE117403_D27 CellTypes UMAP without label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_D27, reduction = "umap", pt.size = 0.5)
dev.off()

#Cell count
Idents(object = GSE117403_D27) <- "stim"
table(Idents(GSE117403_D27))

#FeaturePlots
DefaultAssay(GSE117403_D27) <- "RNA"
tiff(file = "GSE117403_D27 GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D27, reduction = "umap", features = c("GLI1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D27 SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D27, reduction = "umap", features = c("SHH"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D27 AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D27, reduction = "umap", features = c("AR"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_D27 VIM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_D27, reduction = "umap", features = c("VIM"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
####Load data GSE117403_D27####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE120716_Normal/GSE117403_Pd")
GSE117403_Pd.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE120716_Normal/GSE117403_Pd")
GSE117403_Pd <- CreateSeuratObject(counts = GSE117403_Pd.data, min.cells = 3, min.features = 500)
GSE117403_Pd <- NormalizeData(GSE117403_Pd)
GSE117403_Pd <- FindVariableFeatures(GSE117403_Pd, selection.method = "vst", nfeatures = 2000)

#Run the standard workflow for visualization and clustering Ctrl2
GSE117403_Pd <- ScaleData(GSE117403_Pd, verbose = FALSE)
GSE117403_Pd <- RunPCA(GSE117403_Pd, npcs = 30, verbose = FALSE)
ElbowPlot(GSE117403_Pd, ndims = 30)
# UMAP and Clustering
GSE117403_Pd <- FindNeighbors(GSE117403_Pd, reduction = "pca", dims = 1:15)
GSE117403_Pd <- FindClusters(GSE117403_Pd, resolution = 0.7)
GSE117403_Pd <- RunUMAP(GSE117403_Pd, reduction = "pca", dims = 1:15)
DimPlot(GSE117403_Pd, reduction = "umap", pt.size = 0.5, label = TRUE) 

DefaultAssay(GSE117403_Pd) <- "RNA"
tiff(file = "GSE117403_Pd celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_Pd, reduction = "umap", features = c("AR", "EPCAM", "VIM",
                                                            "KRT5", "KLK3", "PIGR", "KRT13",
                                                            "CALCA",
                                                            "FBLN1", "ACTA2", 
                                                            "CLDN5", "RGS1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSE117403_Pd) <- "seurat_clusters"
DimPlot(GSE117403_Pd, reduction = "umap", pt.size = 0.3, label=TRUE)
GSE117403_Pd <- RenameIdents(object = GSE117403_Pd, 
                              '1' = "BE",'5' = "BE",'4' = "BE", '0' = "BE", '4' = "BE", '2' = "BE", '3' = "BE",
                              '9' = "OE", '11' = "OE", '15' = "OE", '8' = "OE", '17' = "OE", '10' = "OE",
                              '6'="LE", '13'="LE",
                              '7' = "FB", '20' = "SM", '12' = "SM", 
                              '18' = "Immune", '19' = "Immune",
                              '21' = "Endo", '14' = "Endo", '16' = "Endo"
)  
GSE117403_Pd[["CellTypes"]] <- Idents(object = GSE117403_Pd)

#UMAP_CellTypes
Idents(object = GSE117403_Pd) <- "CellTypes"
tiff(file = "GSE117403_Pd CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_Pd, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

Idents(object = GSE117403_Pd) <- "CellTypes"
tiff(file = "GSE117403_Pd CellTypes UMAP without label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSE117403_Pd, reduction = "umap", pt.size = 0.5)
dev.off()

#Cell count
Idents(object = GSE117403_Pd) <- "stim"
table(Idents(GSE117403_Pd))

#FeaturePlots
DefaultAssay(GSE117403_Pd) <- "RNA"
tiff(file = "GSE117403_Pd GLI1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_Pd, reduction = "umap", features = c("GLI1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_Pd SHH UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_Pd, reduction = "umap", features = c("SHH"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_Pd AR UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_Pd, reduction = "umap", features = c("AR"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "GSE117403_Pd VIM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSE117403_Pd, reduction = "umap", features = c("VIM"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####SRR13644609####

#Bam file to count matrix
BiocManager::install(c('GenomicFeatures', 'GenomicAlignments', 'DBI'), force = TRUE)
library(GenomicFeatures)
library(GenomicAlignments)
library(DBI)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/SRR13644609")

txdb = makeTxDbFromGFF(file="gencode.v27lift37.annotation.gtf.gz", format="gtf",  organism="Homo sapiens")
saveDb(txdb, file="Homo_sapiens_gencode_v27.sqlite")

txdb  = loadDb("Homo_sapiens_gencode_v27.sqlite")
genes = exonsBy(txdb, by="gene")

filenames = list.files(pattern="possorted_genome_bam10.bam.1")
bamfiles  = BamFileList(filenames, yieldSize=1000000)

se = summarizeOverlaps(features=genes, reads=bamfiles, mode="Union",
                      ignore.strand=FALSE, inter.feature=TRUE, param=ScanBamParam(), preprocess.reads=NULL, singleEnd=TRUE, fragments=FALSE)

as1= assay(se)

write.table(as1, file = "SRR13644609_gene_level_counts.txt", append = FALSE,
            quote = FALSE, sep = "\t", eol = "\n", row.names = TRUE,
            col.names = TRUE)


####RNAseq####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/GSE157703/GSE157703_RAW")

GSE157703_PCa1_raw <- read.table(file = "GSM4773521_PCa1_gene_counts_matrix.csv", sep = ",", header = T, row.names=1, as.is=T)
