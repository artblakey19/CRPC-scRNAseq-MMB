####Load required packages####
library(reshape2)
library(vegan)
library(rgl)
library(gplots)
library(grid)
library(gridExtra)
library(GenomicFeatures)
library(ggplot2)
library(statmod)
library(edgeR)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(tidyverse)
####SU2C_PCF_mCRPC####
setwd("//Users/won/Desktop/Sun Lab/Data/HGF-MET-Bcat/Human/2019_JCI_PNelson")
countdata <- read.table(file = "COUNTS.csv", sep = ",", header = T, row.names=1, as.is=T)

a <- as.matrix(countdata)
#Log Transform
b=log2(a/10+1)
c=log2(a+1)
CRPC2019 =  CreateSeuratObject(counts = c, project = "2019PNelson")

ARNE_meta1 <- CRPC2019@meta.data
ARNE_meta1<- rownames_to_column(ARNE_meta1)
write_csv(ARNE_meta1, "ARNE_meta1.csv")
ARNE_meta <- column_to_rownames(ARNE_meta, var = "rowname")
CRPC2019 <- AddMetaData(CRPC2019,ARNE_meta)

#Clustering
CRPC2019 <- FindVariableFeatures(CRPC2019, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(CRPC2019)
CRPC2019 <- ScaleData(CRPC2019, features = all.genes)
CRPC2019 <- RunPCA(CRPC2019, features = VariableFeatures(object = CRPC2019))

ElbowPlot(CRPC2019, ndims = 50)

CRPC2019 <- FindNeighbors(CRPC2019, reduction = "pca", dims = 1:20)
CRPC2019 <- FindClusters(CRPC2019, resolution = 0.5)
CRPC2019 <- RunUMAP(CRPC2019, reduction = "pca", dims = 1:20)
DimPlot(CRPC2019, reduction = "umap", pt.size = 0.3) 

Idents(object = CRPC2019) <- "ARNE"
DimPlot(CRPC2019, reduction = "umap", pt.size = 0.3) 

#DEGs_
DefaultAssay(CRPC2019) <- "RNA"
all.genes <- rownames(CRPC2019)
CRPC2019 <- ScaleData(CRPC2019, features = all.genes)
DNNEPCvARPC.Markers <- FindMarkers(CRPC2019, ident.1 = c("DNPC1", "DNPC2", "NEPC"), 
                                                              ident.2 = c("ARPC"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(DNNEPCvARPC.Markers, "DNNEPCvARPC.Markers.csv")

DefaultAssay(CRPC2019) <- "RNA"
all.genes <- rownames(CRPC2019)
CRPC2019 <- ScaleData(CRPC2019, features = all.genes)
DNPCvARPC.Markers <- FindMarkers(CRPC2019, ident.1 = c("DNPC1", "DNPC2"), 
                                   ident.2 = c("ARPC"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(DNPCvARPC.Markers, "DNPCvARPC.Markers.csv")

DefaultAssay(CRPC2019) <- "RNA"
all.genes <- rownames(CRPC2019)
CRPC2019 <- ScaleData(CRPC2019, features = all.genes)
DNPC2vARPC.Markers <- FindMarkers(CRPC2019, ident.1 = c("DNPC2"), 
                                 ident.2 = c("ARPC"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(DNPC2vARPC.Markers, "DNPC2vARPC.Markers.csv")

#Rename
Idents(object = CRPC2019) <- "ARNE"
CRPC2019 <- RenameIdents(object = CRPC2019, 
                      'ARPC'="ARPC", 'NEPC'="NEPC",'DNPC1'="DNPC",
                      'DNPC2'="DNPC")  
CRPC2019[["ARNE1"]] <- Idents(object = CRPC2019)

Idents(object = CRPC2019) <- "ARNE1"
DefaultAssay(CRPC2019) <- "RNA"
CRPC2019 <- ScaleData(CRPC2019, features = rownames(CRPC2019))
tiff(file = "CRPC2019 ARNE1 XPO RP family selected.tiff", width = 4.5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(CRPC2019, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                              "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3",
                              "TCF7", "TCF4", "PLAUR", "MMP7","LGR6", "FST", "BMP2", "SNAI1",
                              "XPO1", "XPO5", "XPO7",
                              "MRPL2", "MRPL52", "MRPL34", 
                              "RPS6KA2", "RPS6KL1", "MRPS6", 
                              "EIF5A2", "EIF1B", "EIF1"
), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()
