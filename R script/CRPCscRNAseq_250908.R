#download packages
install.packages('Seurat')
install.packages('devtools')
install.packages('GGally')
install.packages('devtools')
install.packages('devtools')
install.packages('devtools')
install.packages('devtools')
install.packages('devtools')
install.packages('devtools')

#Add necessary tools to library
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(Seurat)
library(devtools)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 

####load or save data####
setwd("D:/CRPC scRNAseq/Results")

####Loading data####
##CRPC1
CRPC1unfiltered.data <- Read10X("D:/CRPC scRNAseq/Raw data/CRPC1")
CRPC1unfiltered <- CreateSeuratObject(counts = CRPC1unfiltered.data,  min.cells = 3, min.features = 500, project = "CRPC1unfiltered")
CRPC1unfiltered <- NormalizeData(CRPC1unfiltered)
####Initial processing, Filtering and Clustering####
#CRPC1unfiltered
CRPC1unfiltered[["percent.mt"]] <- PercentageFeatureSet(CRPC1unfiltered, pattern = "^MT.")

tiff(file = "CRPC1 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC1unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

tiff(file = "CRPC1 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC1unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "CRPC1 Pre-filteration")
dev.off()
tiff(file = "CRPC1 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC1unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "CRPC1 Pre-filteration")
dev.off()

#CRPC1unfiltered
CRPC1 <- subset(CRPC1unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 10)
tiff(file = "CRPC1 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "CRPC1 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC1@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "CRPC1 Post-filteration")
dev.off()
tiff(file = "CRPC1 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC1@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "CRPC1 Post-filteration")
dev.off()

#Genes and UMI counts per cell
mean(CRPC1$nCount_RNA)
mean(CRPC1$nFeature_RNA)

#Clustering
CRPC1 <- FindVariableFeatures(CRPC1, selection.method = "vst", nfeatures = 5000)
tiff(file = "CRPC1 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(CRPC1)
dev.off()

CRPC1 <- ScaleData(CRPC1, verbose = FALSE)
CRPC1 <- RunPCA(CRPC1, npcs = 50, verbose = FALSE)
tiff(file = "CRPC1 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC1, ndims = 50)
dev.off()

CRPC1 <- FindNeighbors(CRPC1, reduction = "pca", dims = 1:20)
CRPC1 <- FindClusters(CRPC1, resolution = 0.5)
CRPC1 <- RunUMAP(CRPC1, reduction = "pca", dims = 1:20)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3) 

#Cell cycle assignment
exp.mat <- read.table(file = "D:/CRPC scRNAseq/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(CRPC1) <- "RNA"
all.genes <- rownames(CRPC1)
CRPC1 <- ScaleData(CRPC1, features = all.genes)
CRPC1 <- CellCycleScoring(CRPC1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = CRPC1) <- "Phase"
DimPlot(CRPC1, reduction = "umap")
tiff(file = "CRPC1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#FeaturePlots for cell markers
DefaultAssay(CRPC1)<-"RNA"
FeaturePlot(CRPC1, reduction = "umap", features = c("AR", "PSA",
                                                       "KRT5", "TRP63", "KRT8", "CD24A", 
                                                       "FBLN1", "MYH11", "PLP1", "PECAM1",
                                                       "RGS5", "TYROBP", "CCL5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#DEGs
DefaultAssay(CRPC1) <- "RNA"
Idents(object = CRPC1) <- "seurat_clusters"
DimPlot(CRPC1, reduction = "umap", label = TRUE)
CRPC1 <- ScaleData(CRPC1, features = rownames(CRPC1))
CRPC1.allMarkers <- FindAllMarkers(CRPC1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(CRPC1.allMarkers, "CRPC1.allMarkers.csv")

#Cell Type Identification
Idents(object = CRPC1) <- "seurat_clusters"
CRPC1 <- RenameIdents(object = CRPC1, 
                         '23'="BE", '24' = "BE", '0'= "BE", '4'="BE",
                         '21'="LE", '10'="LE",'5'="LE",'1'="LE",
                         '3'="LE",'2'="LE",'20'="LE", '7'="LE", '16'="LE", '6'="LE", '11'="LE", '13'="LE",'19'="SV",
                         '15'="FB", '8'="FB",'17'="SM", '22'="Pericyte",
                         '25'="Glia",'14'="VE", '12'="Immune", '18'="Immune",
                         '9'="Immune")  
CRPC1[["CellTypes"]] <- Idents(object = CRPC1)

####Loading data####
##CRPC2
CRPC2unfiltered.data <- Read10X("D:/CRPC scRNAseq/Raw data/CRPC1")
CRPC2unfiltered <- CreateSeuratObject(counts = CRPC2unfiltered,  min.cells = 3, min.features = 500, project = "CRPC2unfiltered")
CRPC2unfiltered <- NormalizeData(CRPC2unfiltered)
####Initial processing, Filtering and Clustering####
#CRPC2unfiltered
CRPC2unfiltered[["percent.mt"]] <- PercentageFeatureSet(CRPC2unfiltered, pattern = "^MT.")

tiff(file = "CRPC2 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC2unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

tiff(file = "CRPC2 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC2unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "CRPC2 Pre-filteration")
dev.off()
tiff(file = "CRPC1 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC1unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "CRPC1 Pre-filteration")
dev.off()

#CRPC1unfiltered
CRPC1 <- subset(CRPC1unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 10)
tiff(file = "CRPC2 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "CRPC2 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "CRPC1 Post-filteration")
dev.off()
tiff(file = "CRPC1 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC1@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "CRPC1 Post-filteration")
dev.off()

#Genes and UMI counts per cell
mean(CRPC1$nCount_RNA)
mean(CRPC1$nFeature_RNA)

#Clustering
CRPC1 <- FindVariableFeatures(CRPC1, selection.method = "vst", nfeatures = 5000)
tiff(file = "CRPC1 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(CRPC1)
dev.off()

CRPC1 <- ScaleData(CRPC1, verbose = FALSE)
CRPC1 <- RunPCA(CRPC1, npcs = 50, verbose = FALSE)
tiff(file = "CRPC1 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC1, ndims = 50)
dev.off()

CRPC1 <- FindNeighbors(CRPC1, reduction = "pca", dims = 1:20)
CRPC1 <- FindClusters(CRPC1, resolution = 0.5)
CRPC1 <- RunUMAP(CRPC1, reduction = "pca", dims = 1:20)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3) 

#Cell cycle assignment
exp.mat <- read.table(file = "D:/CRPC scRNAseq/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(CRPC1) <- "RNA"
all.genes <- rownames(CRPC1)
CRPC1 <- ScaleData(CRPC1, features = all.genes)
CRPC1 <- CellCycleScoring(CRPC1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = CRPC1) <- "Phase"
DimPlot(CRPC1, reduction = "umap")
tiff(file = "CRPC1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#FeaturePlots for cell markers
DefaultAssay(CRPC1)<-"RNA"
FeaturePlot(CRPC1, reduction = "umap", features = c("AR", "PSA",
                                                    "KRT5", "TRP63", "KRT8", "CD24A", 
                                                    "FBLN1", "MYH11", "PLP1", "PECAM1",
                                                    "RGS5", "TYROBP", "CCL5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#DEGs
DefaultAssay(CRPC1) <- "RNA"
Idents(object = CRPC1) <- "seurat_clusters"
DimPlot(CRPC1, reduction = "umap")
CRPC1 <- ScaleData(CRPC1, features = rownames(CRPC1))
CRPC1.allMarkers <- FindAllMarkers(CRPC1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(CRPC1.clusters.allMarkers, "CRPC1.Markers.csv")

#Cell Type Identification
Idents(object = CRPC1) <- "seurat_clusters"
CRPC1 <- RenameIdents(object = CRPC1, 
                      '23'="BE", '24' = "BE", '0'= "BE", '4'="BE",
                      '21'="LE", '10'="LE",'5'="LE",'1'="LE",
                      '3'="LE",'2'="LE",'20'="LE", '7'="LE", '16'="LE", '6'="LE", '11'="LE", '13'="LE",'19'="SV",
                      '15'="FB", '8'="FB",'17'="SM", '22'="Pericyte",
                      '25'="Glia",'14'="VE", '12'="Immune", '18'="Immune",
                      '9'="Immune")  
CRPC1[["CellTypes"]] <- Idents(object = CRPC1)