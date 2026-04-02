#download packages
install.packages('Seurat')
install.packages('SeuratObject')
install.packages('devtools')
install.packages('GGally')
install.packages('patchwork')
install.packages('SeuratWrappers')
install.packages('Azimuth')
install.packages('future')

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
library(devtools)
library(R.utils)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 
library(patchwork)
library(SeuratWrappers)
library(future)

####CRPC1####
setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC1")

####Loading data####
##CRPC1
CRPC1unfiltered.data <- Read10X("D:/CRPC scRNAseq/Raw data/CRPC1/filtered_feature_bc_matrix")
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

#Cell Type Identification
Idents(object = CRPC1) <- "seurat_clusters"
CRPC1 <- RenameIdents(object = CRPC1, 
                      '23'="stim", '24' = "stim", '0'= "stim", '4'="stim",
                      '21'="stim", '10'="stim",'5'="stim",'1'="stim",
                      '3'="stim",'2'="stim",'20'="stim", '7'="stim", '16'="stim", '6'="stim", '11'="stim", '13'="stim",'19'="stim",
                      '15'="stim", '8'="stim",'17'="stim", '22'="stim",
                      '25'="stim",'14'="stim", '12'="stim", '18'="stim",
                      '9'="stim")  
CRPC1[["Original"]] <- Idents(object = CRPC1)

Idents(object = CRPC1) <- "Original"
tiff(file = "CRPC1 lightgrey UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3, cols = c("lightgrey")) + NoLegend()
dev.off()

####CRPC2####
setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC2")

####Loading data####
##CRPC2
CRPC2unfiltered.data <- Read10X("D:/CRPC scRNAseq/Raw data/CRPC2/filtered_feature_bc_matrix")
CRPC2unfiltered <- CreateSeuratObject(counts = CRPC2unfiltered.data,  min.cells = 3, min.features = 500, project = "CRPC2unfiltered")
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
tiff(file = "CRPC2 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC2unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "CRPC2 Pre-filteration")
dev.off()

#CRPC2unfiltered
CRPC2 <- subset(CRPC2unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 10)
tiff(file = "CRPC2 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "CRPC2 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "CRPC2 Post-filteration")
dev.off()
tiff(file = "CRPC2 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "CRPC2 Post-filteration")
dev.off()

#Genes and UMI counts per cell
mean(CRPC2$nCount_RNA)
mean(CRPC2$nFeature_RNA)

#Clustering
CRPC2 <- FindVariableFeatures(CRPC2, selection.method = "vst", nfeatures = 5000)
tiff(file = "CRPC2 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(CRPC2)
dev.off()

CRPC2 <- ScaleData(CRPC2, verbose = FALSE)
CRPC2 <- RunPCA(CRPC2, npcs = 50, verbose = FALSE)
tiff(file = "CRPC2 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC2, ndims = 50)
dev.off()

CRPC2 <- FindNeighbors(CRPC2, reduction = "pca", dims = 1:20)
CRPC2 <- FindClusters(CRPC2, resolution = 0.5)
CRPC2 <- RunUMAP(CRPC2, reduction = "pca", dims = 1:20)
DimPlot(CRPC2, reduction = "umap", pt.size = 0.3) 

#Cell cycle assignment
exp.mat <- read.table(file = "D:/CRPC scRNAseq/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(CRPC2) <- "RNA"
all.genes <- rownames(CRPC2)
CRPC2 <- ScaleData(CRPC2, features = all.genes)
CRPC2 <- CellCycleScoring(CRPC2, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = CRPC2) <- "Phase"
DimPlot(CRPC2, reduction = "umap")
tiff(file = "CRPC2 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC2, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#DEGs
DefaultAssay(CRPC2) <- "RNA"
Idents(object = CRPC2) <- "seurat_clusters"
DimPlot(CRPC2, reduction = "umap")
CRPC2 <- ScaleData(CRPC2, features = rownames(CRPC2))
CRPC2.allMarkers <- FindAllMarkers(CRPC2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(CRPC2.allMarkers, "CRPC2.allMarkers.csv")

#Cell Type Identification
Idents(object = CRPC2) <- "seurat_clusters"
CRPC2 <- RenameIdents(object = CRPC2, 
                      '23'="BE", '24' = "BE", '0'= "BE", '4'="BE",
                      '21'="LE", '10'="LE",'5'="LE",'1'="LE",
                      '3'="LE",'2'="LE",'20'="LE", '7'="LE", '16'="LE", '6'="LE", '11'="LE", '13'="LE",'19'="SV",
                      '15'="FB", '8'="FB",'17'="SM", '22'="Pericyte",
                      '25'="Glia",'14'="VE", '12'="Immune", '18'="Immune",
                      '9'="Immune")  
CRPC2[["CellTypes"]] <- Idents(object = CRPC2)

Idents(object = CRPC2) <- "seurat_clusters"
CRPC2 <- RenameIdents(object = CRPC2, 
                      '23'="Original", '24' = "Original", '0'= "Original", '4'="Original",
                      '21'="Original", '10'="Original",'5'="Original",'1'="Original",
                      '3'="Original",'2'="Original",'20'="Original", '7'="Original", '16'="Original", '6'="Original", 
                      '11'="Original", '13'="Original",'19'="Original",
                      '15'="Original", '8'="Original",'17'="Original", '22'="Original",
                      '25'="Original",'14'="Original", '12'="Original", '18'="Original",
                      '9'="Original")  
CRPC2[["Original"]] <- Idents(object = CRPC2)

Idents(object = CRPC2) <- "Original"
tiff(file = "CRPC2 darkblue UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC2, reduction = "umap", pt.size = 0.3, cols = c("darkblue")) + NoLegend()
dev.off()


####CRPC3####
setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC3")

####Loading data####
##CRPC3
CRPC3unfiltered.data <- Read10X("D:/CRPC scRNAseq/Raw data/CRPC3/filtered_feature_bc_matrix")
CRPC3unfiltered <- CreateSeuratObject(counts = CRPC3unfiltered.data,  min.cells = 3, min.features = 500, project = "CRPC3unfiltered")
CRPC3unfiltered <- NormalizeData(CRPC3unfiltered)

####Initial processing, Filtering and Clustering####
#CRPC3unfiltered
CRPC3unfiltered[["percent.mt"]] <- PercentageFeatureSet(CRPC3unfiltered, pattern = "^MT.")

tiff(file = "CRPC3 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC3unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()

tiff(file = "CRPC3 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC3unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "CRPC3 Pre-filteration")
dev.off()
tiff(file = "CRPC3 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC3unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "CRPC3 Pre-filteration")
dev.off()

#CRPC3unfiltered
CRPC3 <- subset(CRPC3unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 10)
tiff(file = "CRPC3 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CRPC3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "CRPC3 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC3@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "CRPC3 Post-filteration")
dev.off()
tiff(file = "CRPC3 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(CRPC3@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "CRPC3 Post-filteration")
dev.off()

#Genes and UMI counts per cell
mean(CRPC3$nCount_RNA)
mean(CRPC3$nFeature_RNA)

#Clustering
CRPC3 <- FindVariableFeatures(CRPC3, selection.method = "vst", nfeatures = 5000)
tiff(file = "CRPC3 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(CRPC3)
dev.off()

CRPC3 <- ScaleData(CRPC3, verbose = FALSE)
CRPC3 <- RunPCA(CRPC3, npcs = 50, verbose = FALSE)
tiff(file = "CRPC3 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC3, ndims = 50)
dev.off()

CRPC3 <- FindNeighbors(CRPC3, reduction = "pca", dims = 1:20)
CRPC3 <- FindClusters(CRPC3, resolution = 0.5)
CRPC3 <- RunUMAP(CRPC3, reduction = "pca", dims = 1:20)
DimPlot(CRPC3, reduction = "umap", pt.size = 0.3) 

#Cell cycle assignment
exp.mat <- read.table(file = "D:/CRPC scRNAseq/cell_cycle_vignette_files/nestorawa_forcellcycle_expressionMatrix.txt",
                      header = TRUE, as.is = TRUE, row.names = 1)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(CRPC3) <- "RNA"
all.genes <- rownames(CRPC3)
CRPC3 <- ScaleData(CRPC3, features = all.genes)
CRPC3 <- CellCycleScoring(CRPC3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = CRPC3) <- "Phase"
DimPlot(CRPC3, reduction = "umap")
tiff(file = "CRPC3 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC3, reduction = "umap", pt.size = 0.3, cols = c("goldenrod2", "grey75", "deepskyblue2"))
dev.off()

#DEGs
DefaultAssay(CRPC3) <- "RNA"
Idents(object = CRPC3) <- "seurat_clusters"
DimPlot(CRPC3, reduction = "umap")
CRPC3 <- ScaleData(CRPC3, features = rownames(CRPC3))
CRPC3.allMarkers <- FindAllMarkers(CRPC3, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(CRPC3.allMarkers, "CRPC3.allMarkers.csv")

#Cell Type Identification
Idents(object = CRPC3) <- "seurat_clusters"
CRPC3 <- RenameIdents(object = CRPC3, 
                      '23'="BE", '24' = "BE", '0'= "BE", '4'="BE",
                      '21'="LE", '10'="LE",'5'="LE",'1'="LE",
                      '3'="LE",'2'="LE",'20'="LE", '7'="LE", '16'="LE", '6'="LE", '11'="LE", '13'="LE",'19'="SV",
                      '15'="FB", '8'="FB",'17'="SM", '22'="Pericyte",
                      '25'="Glia",'14'="VE", '12'="Immune", '18'="Immune",
                      '9'="Immune")  
CRPC3[["CellTypes"]] <- Idents(object = CRPC3)

Idents(object = CRPC3) <- "seurat_clusters"
CRPC3 <- RenameIdents(object = CRPC3, 
                      '23'="Original", '24' = "Original", '0'= "Original", '4'="Original",
                      '21'="Original", '10'="Original",'5'="Original",'1'="Original",
                      '3'="Original",'2'="Original",'20'="Original", '7'="Original", 
                      '16'="Original", '6'="Original", '11'="Original", '13'="Original",'19'="Original",
                      '15'="Original", '8'="Original",'17'="Original", '22'="Original",
                      '25'="Original",'14'="Original", '12'="Original", '18'="Original",
                      '9'="Original")
CRPC3[["Original"]] <- Idents(object = CRPC3)

Idents(object = CRPC3) <- "Original"
tiff(file = "CRPC3 salmon UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC3, reduction = "umap", pt.size = 0.3, cols = c("salmon")) + NoLegend()
dev.off()
####Merging Datasets####

setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC.combined")

CRPC1[["orig.clusters"]] <- Idents(object = CRPC1)
CRPC2[["orig.clusters"]] <- Idents(object = CRPC2)
CRPC3[["orig.clusters"]] <- Idents(object = CRPC3)

CRPC1$orig.ident <- "CRPC1"
CRPC2$orig.ident <- "CRPC2"
CRPC3$orig.ident <- "CRPC3"

Idents(object = CRPC1) <- "seurat_clusters"
Idents(object = CRPC2) <- "seurat_clusters"
Idents(object = CRPC3) <- "seurat_clusters"

CRPC1$stim <- "CRPC1"
CRPC2$stim <- "CRPC2"
CRPC3$stim <- "CRPC3"

CRPC.combined <- merge(
  x = CRPC1,
  y = list(CRPC2, CRPC3),
  add.cell.ids = c("CRPC1", "CRPC2", "CRPC3"),
  project = "CRPC"
)

DefaultAssay(CRPC.combined) <- "RNA"
CRPC.combined <- ScaleData(CRPC.combined, verbose = FALSE)
CRPC.combined <- RunPCA(CRPC.combined, npcs = 50, verbose = FALSE)

CRPC.combined <- JoinLayers(CRPC.combined, assay = "RNA", layers = c("counts", "data"))
CRPC.combined[["RNA"]] <- split(
  CRPC.combined[["RNA"]],
  f = CRPC.combined$stim,
  layers = c("counts", "data")
)

## 7) IntegrateLayers (RPCA)
CRPC.combined <- IntegrateLayers(
  object = CRPC.combined,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  dims = 1:30,
  verbose = FALSE
)

CRPC.combined <- FindNeighbors(CRPC.combined, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)
CRPC.combined <- FindClusters(CRPC.combined, resolution = 0.5, verbose = FALSE)
CRPC.combined <- RunUMAP(CRPC.combined, reduction = "integrated.rpca", dims = 1:30, verbose = FALSE)

## 9) 확인 플롯
DimPlot(CRPC.combined, group.by = "stim", reduction = "umap")

Idents(object = CRPC.combined) <- "stim"
tiff(file = "CRPC.combined stim UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined, reduction = "umap", pt.size = 0.3, cols = c("lightgrey", "darkblue", "salmon")) + NoLegend()
dev.off()
tiff(file = "CRPC.combined split stim UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("lightgrey", "darkblue", "salmon")) + NoLegend()
dev.off()

#FeaturePlots for cell markers
DefaultAssay(CRPC.combined)<-"RNA"
tiff(file = "CRPC.combined Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("KLK2", "MSMB", "TMPRSS2",
                                                            "EPCAM", "KRT8", "KRT5",
                                                            "ACPP", "NKX3-1", "AZGP1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined Endo expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("VWF", "SELE", "IFI27",
                                                            "FLT1", "SPARCL1", "PTPRB",
                                                            "SDPR", "DARC", "PLVAP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined FB expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("DCN", "FBLN1", "COL1A2",
                                                            "C7", "IGF1", "IGFBP5",
                                                            "CCDC80", "CFD", "LTBP4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined MYELOID expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("IL1B", "HLA-DRA", "HLA-DPA1",
                                                            "HLA-DPB1", "HLA-DRB1", "CD74",
                                                            "IL8", "IFI30", "LYZ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined PLASMA expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("IGKC", "IGHA1", "IGJ",
                                                            "IGHA2", "MZB1", "IGHG3",
                                                            "SLAMF7", "IGHG4", "IGHG1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined SM expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("MYH11", "RGS5", "ACTA2",
                                                            "TAGLN", "MYL9", "MYLK",
                                                            "MCAM", "CALD1", "LMOD1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined T CELL expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("IL7R", "TRBC2", "CCL5",
                                                            "IFNG", "CD8A", "ETS1",
                                                            "CCL4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined MAST CELL expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("CPA3", "TPSAB1", "KIT",
                                                            "VWA5A", "IL1RL1", "CTSG",
                                                            "ACSL4", "MS4A2", "GATA2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined B CELL expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined, reduction = "umap", features = c("MS4A1", "CXCR5", "CD79A",
                                                            "BANK1", "LY9", "CCR7",
                                                            "IRF8"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#New Labeling
Idents(object = CRPC.combined) <- "seurat_clusters"

CRPC.combined <- RenameIdents(object = CRPC.combined, 
                              '0' = "Basal Epithelium",'11' = "Basal Epithelium",
                              '8' = "Basal Epithelium",'4' = "Basal Epithelium",
                              '10' = "Basal Epithelium",'1' = "Basal Epithelium",
                              '6' = "Basal Epithelium",'15' = "Basal Epithelium",
                              '22' = "Basal Epithelium",
                              '13' = "Luminal Epithelium",'5' = "Luminal Epithelium",
                              '14' = "Luminal Epithelium",'16' = "Luminal Epithelium",
                              '20' = "Other Epithelium",
                              '9' = "Fibroblast",'3' = "Fibroblast",'18' = "Fibroblast",
                              '12' = "Smooth Muscle",
                              '2' = "Endothelium",
                              '7' = "Lymphocyte",'19' = "Lymphocyte",
                              '17' = "Myeloid",'21' = "Mast Cells", '23' = "Other stroma"
                              )
CRPC.combined[["celltype"]] <- Idents(object = CRPC.combined)

Idents(object = CRPC.combined) <- "celltype"
DimPlot(CRPC.combined, reduction = "umap", pt.size = 0.3) 
Idents(object = CRPC.combined) <- "celltype"
DimPlot(CRPC.combined, reduction = "umap", pt.size = 0.3,  split.by = "stim") 



saveRDS(CRPC.combined, file = "C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC.combined/CRPC.combined.rds")

#Cell cycle assignment
DefaultAssay(CRPC.combined) <- "RNA"
all.genes <- rownames(CRPC.combined)
CRPC.combined <- ScaleData(CRPC.combined, features = all.genes)
CRPC.combined <- CellCycleScoring(CRPC.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = CRPC.combined) <- "Phase"
DimPlot(CRPC.combined, reduction = "umap")
tiff(file = "CRPC.combined Cell Cycle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell Cycle regression
CRPC.combined1 <- CRPC.combined
DefaultAssay(CRPC.combined1) <- "integrated"
CRPC.combined1 <- ScaleData(CRPC.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CRPC.combined1))
CRPC.combined1 <- RunPCA(CRPC.combined1, features = VariableFeatures(CRPC.combined1))
ElbowPlot(CRPC.combined1, ndims = 50)

CRPC.combined1 <- FindNeighbors(CRPC.combined1, reduction = "pca", dims = 1:30)
CRPC.combined1 <- FindClusters(CRPC.combined1, resolution = 0.8)
CRPC.combined1 <- RunUMAP(CRPC.combined1, reduction = "pca", dims = 1:30)

Idents(object = CRPC.combined1) <- "seurat_clusters"
DimPlot(CRPC.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(CRPC.combined1, reduction = "umap", pt.size = 0.3, split.by='stim',label = TRUE)

Idents(object = CRPC.combined1) <- "seurat_clusters"
tiff(file = "CRPC.combined1 UMAP dim30.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = CRPC.combined1) <- "Phase"
tiff(file = "CRPC.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

#FeaturePlots for cell markers
DefaultAssay(CRPC.combined1)<-"RNA"
tiff(file = "CRPC.combined1 Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("KLK2", "MSMB", "TMPRSS2",
                                                            "EPCAM", "KRT8", "KRT5",
                                                            "ACPP", "NKX3-1", "AZGP1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 Endo expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("VWF", "SELE", "IFI27",
                                                            "FLT1", "SPARCL1", "PTPRB",
                                                            "SDPR", "DARC", "PLVAP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 FB expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("DCN", "FBLN1", "COL1A2",
                                                            "C7", "IGF1", "IGFBP5",
                                                            "CCDC80", "CFD", "LTBP4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 MYELOID expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("IL1B", "HLA-DRA", "HLA-DPA1",
                                                            "HLA-DPB1", "HLA-DRB1", "CD74",
                                                            "IL8", "IFI30", "LYZ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 PLASMA expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("IGKC", "IGHA1", "IGJ",
                                                            "IGHA2", "MZB1", "IGHG3",
                                                            "SLAMF7", "IGHG4", "IGHG1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 SM expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("MYH11", "RGS5", "ACTA2",
                                                            "TAGLN", "MYL9", "MYLK",
                                                            "MCAM", "CALD1", "LMOD1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 T CELL expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("IL7R", "TRBC2", "CCL5",
                                                            "IFNG", "CD8A", "ETS1",
                                                            "CCL4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 MAST CELL expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("CPA3", "TPSAB1", "KIT",
                                                            "VWA5A", "IL1RL1", "CTSG",
                                                            "ACSL4", "MS4A2", "GATA2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "CRPC.combined1 B CELL expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC.combined1, reduction = "umap", features = c("MS4A1", "CXCR5", "CD79A",
                                                            "BANK1", "LY9", "CCR7",
                                                            "IRF8"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

FeaturePlot(CRPC.combined1, reduction = "umap", features = c("CHGA"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#New Labeling
Idents(object = CRPC.combined1) <- "seurat_clusters"
CRPC.combined1 <- RenameIdents(object = CRPC.combined1, 
                              '17' = "Epithelium",'5' = "Epithelium",'11' = "Epithelium",
                              '4' = "Epithelium",'7' = "Epithelium",'0' = "Epithelium",
                              '2' = "Epithelium",'20' = "Epithelium",'22' = "Epithelium",
                              '13' = "Epithelium",'12' = "Epithelium",'1' = "Epithelium",
                              '15' = "Epithelium",
                              '3' = "Fibroblast",'8' = "Fibroblast",'16' = "Fibroblast",
                              '10' = "Smooth Muscle",
                              '9' = "Endothelium",'14' = "Endothelium",
                              '24' = "Plasma",
                              '6' = "Lymphocyte",'19' = "Lymphocyte",
                              '18' = "Myeloids",'21' = "Mast Cells", '23' = "Other stroma"
)
CRPC.combined1[["celltype"]] <- Idents(object = CRPC.combined1)

Idents(object = CRPC.combined1) <- "celltype"
DimPlot(CRPC.combined1, reduction = "umap", pt.size = 0.3) 
tiff(file = "CRPC.combined1 celltype UMAP.tiff", width = 7, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined1, reduction = "umap", pt.size = 0.3) 
dev.off()
tiff(file = "CRPC.combined1 celltype split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC.combined1, reduction = "umap",split.by = "stim", pt.size = 0.3) 
dev.off()

#Cell counts
Idents(object = CRPC.combined1) <- "stim"
table(Idents(CRPC.combined1))

#####Subset Epithelial cells####

setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC.combined.epi")

#subset primary
Idents(object = CRPC.combined1) <- "celltype"
DimPlot(CRPC.combined1, reduction = "umap")
CRPC.combined1.Epi <- subset(CRPC.combined1, idents = c("Epithelium"))

CRPC.combined1.Epi[["RNA"]] <- split(CRPC.combined1.Epi[["RNA"]], f = CRPC.combined1.Epi$stim)
CRPC.combined1.Epi

#Clustering
CRPC.combined1.Epi <- NormalizeData(CRPC.combined1.Epi)
CRPC.combined1.Epi <- FindVariableFeatures(CRPC.combined1.Epi)
CRPC.combined1.Epi <- ScaleData(CRPC.combined1.Epi)
CRPC.combined1.Epi <- RunPCA(CRPC.combined1.Epi)
ElbowPlot(CRPC.combined1.Epi, ndims = 50)

CRPC.combined1.Epi <- FindNeighbors(CRPC.combined1.Epi, reduction = "pca", dims = 1:20)
CRPC.combined1.Epi <- FindClusters(CRPC.combined1.Epi, resolution = 1, cluster.name = "unintegrated_clusters")
CRPC.combined1.Epi <- RunUMAP(CRPC.combined1.Epi, reduction = "pca", dims = 1:20, reduction.name = "umap.unintegrated")
DimPlot(CRPC.combined1.Epi, reduction = "umap.unintegrated", group.by = c("stim", "seurat_clusters"))

#CCAIntegration
CRPC.combined1.Epi <- IntegrateLayers(object = CRPC.combined1.Epi, 
                                      method = CCAIntegration, orig.reduction = "pca", 
                                      new.reduction = "integrated.cca",
                        verbose = FALSE)

#RPCAIntegration
CRPC.combined1.Epi <- IntegrateLayers(
  object = CRPC.combined1.Epi,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE
)

CRPC.combined1.Epi <- FindNeighbors(CRPC.combined1.Epi, reduction = "integrated.cca", dims = 1:30)
CRPC.combined1.Epi <- FindClusters(CRPC.combined1.Epi, resolution = 1, cluster.name = "cca_clusters")
CRPC.combined1.Epi <- RunUMAP(CRPC.combined1.Epi, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
p1 <- DimPlot(
  CRPC.combined1.Epi,
  reduction = "umap.cca",
  group.by = c("stim", "seurat_clusters", "cca_clusters"),
  combine = FALSE, label.size = 2
)

CRPC.combined1.Epi <- FindNeighbors(CRPC.combined1.Epi, reduction = "integrated.rpca", dims = 1:30)
CRPC.combined1.Epi <- FindClusters(CRPC.combined1.Epi, resolution = 2, cluster.name = "rpca_clusters")
CRPC.combined1.Epi <- RunUMAP(CRPC.combined1.Epi, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
p2 <- DimPlot(
  CRPC.combined1.Epi,
  reduction = "umap.rpca",
  group.by = c("stim", "seurat_clusters", "rpca_clusters"),
  combine = FALSE, label.size = 2
)

wrap_plots(c(p1, p2), ncol = 2, byrow = F)

# re-join layers after integration
CRPC.combined1.Epi[["RNA"]] <- JoinLayers(CRPC.combined1.Epi[["RNA"]])

CRPC.combined1.Epi <- FindNeighbors(CRPC.combined1.Epi, reduction = "integrated.cca", dims = 1:20)
CRPC.combined1.Epi <- FindClusters(CRPC.combined1.Epi, resolution = 1)
DimPlot(CRPC.combined1.Epi, reduction = "umap", split.by = "stim", label = TRUE)


