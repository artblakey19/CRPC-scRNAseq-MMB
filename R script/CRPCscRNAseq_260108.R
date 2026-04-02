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
CRPC1unfiltered.data <- Read10X("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC1/filtered_feature_bc_matrix")
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
CRPC1 <- subset(CRPC1unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 15)
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

all.genes <- rownames(CRPC1)
CRPC1 <- ScaleData(CRPC1, features = all.genes)
CRPC1 <- RunPCA(CRPC1, features = VariableFeatures(object = CRPC1))

tiff(file = "CRPC1 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC1, ndims = 50)
dev.off()

CRPC1 <- FindNeighbors(CRPC1, reduction = "pca", dims = 1:30)
CRPC1 <- FindClusters(CRPC1, resolution = 0.3)
CRPC1 <- RunUMAP(CRPC1, reduction = "pca", dims = 1:30)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3, label = TRUE) 

#FeaturePlots for cell markers
DefaultAssay(CRPC1)<-"RNA"
tiff(file = "CRPC1 Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC1, reduction = "umap", features = c("AR", "KLK2", "MSMB", "TMPRSS2",
                                                            "EPCAM", "KRT8", "KRT5",
                                                            "CDH1", "NKX3-1", "AZGP1", "VIM"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell Type Identification
Idents(object = CRPC1) <- "seurat_clusters"
CRPC1 <- RenameIdents(object = CRPC1, 
                      '0'="CRPC1", '1' = "CRPC1", '2'= "CRPC1", '3'="CRPC1",
                      '4'="CRPC1", '5'="CRPC1",'6'="CRPC1",'7'="CRPC1",
                      '8'="CRPC1",'9'="CRPC1",'10'="CRPC1", '11'="CRPC1", '12'="CRPC1", 
                      '13'="CRPC1", '14'="CRPC1", '13'="CRPC1",'19'="CRPC1",
                      '15'="CRPC1", '16'="CRPC1")  
CRPC1[["Original"]] <- Idents(object = CRPC1)

Idents(object = CRPC1) <- "Original"
tiff(file = "CRPC1 lightgrey UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC1, reduction = "umap", pt.size = 0.3, cols = c("lightgrey")) + NoLegend()
dev.off()

#subset epi
Idents(object = CRPC1) <- "seurat_clusters"
CRPC1.epi <- subset(CRPC1, idents = c("1", "5", "0", "13", "6", "10", "14", "2", "9"))

####CRPC2####
setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC2")

##CRPC2
CRPC2unfiltered.data <- Read10X("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC2/filtered_feature_bc_matrix")
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
CRPC2 <- subset(CRPC2unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 8000 & percent.mt < 15)
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

all.genes <- rownames(CRPC2)
CRPC2 <- ScaleData(CRPC2, features = all.genes)
CRPC2 <- RunPCA(CRPC2, features = VariableFeatures(object = CRPC2))

tiff(file = "CRPC2 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC2, ndims = 50)
dev.off()

CRPC2 <- FindNeighbors(CRPC2, reduction = "pca", dims = 1:30)
CRPC2 <- FindClusters(CRPC2, resolution = 0.3)
CRPC2 <- RunUMAP(CRPC2, reduction = "pca", dims = 1:30)
DimPlot(CRPC2, reduction = "umap", pt.size = 0.3, label = TRUE) 

#FeaturePlots for cell markers
DefaultAssay(CRPC2)<-"RNA"
tiff(file = "CRPC2 Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC2, reduction = "umap", features = c("AR", "KLK2", "MSMB", "TMPRSS2",
                                                    "EPCAM", "KRT8", "KRT5",
                                                    "CDH1", "NKX3-1", "AZGP1", "VIM"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell Type Identification
Idents(object = CRPC2) <- "seurat_clusters"
CRPC2 <- RenameIdents(object = CRPC2, 
                      '0'="CRPC2", '1' = "CRPC2", '2'= "CRPC2", '3'="CRPC2",
                      '4'="CRPC2", '5'="CRPC2",'6'="CRPC2",'7'="CRPC2",
                      '8'="CRPC2",'9'="CRPC2",'10'="CRPC2", '11'="CRPC2", '12'="CRPC2", 
                      '13'="CRPC2", '14'="CRPC2")  
CRPC2[["Original"]] <- Idents(object = CRPC2)

Idents(object = CRPC2) <- "Original"
tiff(file = "CRPC2 lightseagreen UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen")) + NoLegend()
dev.off()

#subset epi
Idents(object = CRPC2) <- "seurat_clusters"
CRPC2.epi <- subset(CRPC2, idents = c("11", "5", "6", "3", "0", "2", "10"))

####CRPC3####
setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC3")

##CRPC3
CRPC3unfiltered.data <- Read10X("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC3/filtered_feature_bc_matrix")
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
CRPC3 <- subset(CRPC3unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 15)
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

all.genes <- rownames(CRPC3)
CRPC3 <- ScaleData(CRPC3, features = all.genes)
CRPC3 <- RunPCA(CRPC3, features = VariableFeatures(object = CRPC3))

tiff(file = "CRPC3 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(CRPC3, ndims = 50)
dev.off()

CRPC3 <- FindNeighbors(CRPC3, reduction = "pca", dims = 1:30)
CRPC3 <- FindClusters(CRPC3, resolution = 0.3)
CRPC3 <- RunUMAP(CRPC3, reduction = "pca", dims = 1:30)
DimPlot(CRPC3, reduction = "umap", pt.size = 0.3, label = TRUE) 

#FeaturePlots for cell markers
DefaultAssay(CRPC3)<-"RNA"
tiff(file = "CRPC3 Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(CRPC3, reduction = "umap", features = c("AR", "KLK2", "MSMB", "TMPRSS2",
                                                    "EPCAM", "KRT8", "KRT5",
                                                    "CDH1", "NKX3-1", "AZGP1", "VIM"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell Type Identification
Idents(object = CRPC3) <- "seurat_clusters"
CRPC3 <- RenameIdents(object = CRPC3, 
                      '0'="CRPC3", '1' = "CRPC3", '2'= "CRPC3", '3'="CRPC3",
                      '4'="CRPC3", '5'="CRPC3",'6'="CRPC3",'7'="CRPC3",
                      '8'="CRPC3",'9'="CRPC3",'10'="CRPC3", '11'="CRPC3", '12'="CRPC3", 
                      '13'="CRPC3", '14'="CRPC3", '15'="CRPC3")  
CRPC3[["Original"]] <- Idents(object = CRPC3)

Idents(object = CRPC3) <- "Original"
tiff(file = "CRPC3 orange UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC3, reduction = "umap", pt.size = 0.3, cols = c("orange")) + NoLegend()
dev.off()

#subset epi
Idents(object = CRPC3) <- "seurat_clusters"
CRPC3.epi <- subset(CRPC3, idents = c("4", "0", "10", "7", "1", "3", "13", "2"))

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

####Merge Epithelial cells####

#CRPC1.epi
setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC1")

Idents(object = CRPC1) <- "seurat_clusters"
DimPlot(CRPC1, reduction = "umap", label = TRUE)

DefaultAssay(CRPC1)<-"RNA"
FeaturePlot(CRPC1, reduction = "umap", features = c("KLK2", "MSMB", "TMPRSS2",
                                                             "EPCAM", "KRT8", "KRT5",
                                                             "ACPP", "NKX3-1", "AZGP1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")


CRPC1 <- subset(CRPC1, idents = c("Epithelium"))


library(Seurat)

setwd("C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC.combined.epi")

# (A) 메타데이터 세팅(질문에서 작성한 그대로)
CRPC1.epi[["orig.clusters"]] <- Idents(CRPC1.epi)
CRPC2.epi[["orig.clusters"]] <- Idents(CRPC2.epi)
CRPC3.epi[["orig.clusters"]] <- Idents(CRPC3.epi)

CRPC1.epi$orig.ident <- "CRPC1.epi"
CRPC2.epi$orig.ident <- "CRPC2.epi"
CRPC3.epi$orig.ident <- "CRPC3.epi"

Idents(CRPC1.epi) <- "seurat_clusters"
Idents(CRPC2.epi) <- "seurat_clusters"
Idents(CRPC3.epi) <- "seurat_clusters"

CRPC1.epi$stim <- "CRPC1.epi"
CRPC2.epi$stim <- "CRPC2.epi"
CRPC3.epi$stim <- "CRPC3.epi"

# RNA assay 존재 확인 및 기본 assay 통일
DefaultAssay(CRPC1.epi) <- "RNA"
DefaultAssay(CRPC2.epi) <- "RNA"
DefaultAssay(CRPC3.epi) <- "RNA"

# (권장) counts만 남기고 data/scale 등은 제거하여 merge 부담 최소화
CRPC1.epi <- DietSeurat(CRPC1.epi, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
CRPC2.epi <- DietSeurat(CRPC2.epi, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)
CRPC3.epi <- DietSeurat(CRPC3.epi, assays = "RNA", layers = "counts", dimreducs = NULL, graphs = NULL)

gc()

# (B) Merge
epi.merged <- merge(
  x = CRPC1.epi,
  y = list(CRPC2.epi, CRPC3.epi),
  add.cell.ids = c("CRPC1", "CRPC2", "CRPC3"),
  project = "CRPC_epi"
)
DefaultAssay(epi.merged) <- "RNA"

# (C) Split
epi.list <- SplitObject(epi.merged, split.by = "orig.ident")

# (D) LogNormalize + HVG (샘플당 2000 HVG)
epi.list <- lapply(epi.list, function(x) {
  x <- NormalizeData(x, normalization.method = "LogNormalize", scale.factor = 1e4, verbose = FALSE)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  x
})

# (E) 통합 feature 선택: 3000까지 뽑되, 실제 통합은 2000으로 제한(아래에서)
features <- SelectIntegrationFeatures(object.list = epi.list, nfeatures = 3000)

# (F) RPCA용 PCA 선행(동일 features로 scale + PCA)
epi.list <- lapply(epi.list, function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, npcs = 50, verbose = FALSE)
  x
})

# 메모리 회수(중간 객체 정리)
gc()

# (G) RPCA anchors
anchors <- FindIntegrationAnchors(
  object.list = epi.list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:20,     # 1:30 -> 1:20
  k.anchor = 5
)

# 통합 직전: 가장 큰 메모리 먹는 list 제거(anchors만 남기고)
rm(epi.list)
gc()

# (H) IntegrateData: 핵심은 features.to.integrate로 통합 feature 제한
features.to.integrate <- features[1:1500]  # 1500까지 줄여도 됨

epi.integrated <- IntegrateData(
  anchorset = anchors,
  dims = 1:20,
  features.to.integrate = features.to.integrate,
  k.weight = 50,                 # 기본 100 -> 50 (필요시 30)
  new.assay.name = "integrated"
)

# anchors도 정리
rm(anchors)
gc()

# (I) Downstream
DefaultAssay(epi.integrated) <- "integrated"
epi.integrated <- ScaleData(epi.integrated, verbose = FALSE)
epi.integrated <- RunPCA(epi.integrated, npcs = 50, verbose = FALSE)
epi.integrated <- RunUMAP(epi.integrated, dims = 1:20)
epi.integrated <- FindNeighbors(epi.integrated, dims = 1:20)
epi.integrated <- FindClusters(epi.integrated, resolution = 0.5)

DimPlot(epi.integrated, group.by = "orig.ident")

Idents(object = epi.integrated) <- "seurat_clusters"
tiff(file = "epi.integrated UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(epi.integrated, reduction = "umap", pt.size = 0.3, label = T)
dev.off()

##Subset Epi
Idents(object = epi.integrated) <- "seurat_clusters"
epi.integrated1 <- subset(epi.integrated, idents = c("0", "1", "2", "3", "4", "5", 
                                                     "6", "7", "8", "9", "10", "11"
                                                     , "13", "14", "15", "16"))

#Downstream
DefaultAssay(epi.integrated1) <- "integrated"
epi.integrated1 <- ScaleData(epi.integrated1, verbose = FALSE)
epi.integrated1 <- RunPCA(epi.integrated1, npcs = 50, verbose = FALSE)
epi.integrated1 <- RunUMAP(epi.integrated1, dims = 1:20)
epi.integrated1 <- FindNeighbors(epi.integrated1, dims = 1:20)
epi.integrated1 <- FindClusters(epi.integrated1, resolution = 0.5)

Idents(object = epi.integrated1) <- "seurat_clusters"
DimPlot(epi.integrated1, reduction = "umap", pt.size = 0.3, split.by = "stim", label = T)

#FeaturePlots for cell markers
DefaultAssay(epi.integrated1)<-"RNA"
tiff(file = "epi.integrated1 Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated1, reduction = "umap", features = c("AR", "EPCAM", "CDH1", "VIM",
                                                              "KRT5", "KRT14", "KRT15", "TP63", "KRT17",
                                                              "KRT8", "KRT19", "NKX3-1", "MSMB", "TMPRSS2", "KLK3", "KLK2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell Type Identification
Idents(object = epi.integrated1) <- "seurat_clusters"
epi.integrated1 <- RenameIdents(object = epi.integrated1, 
                      '0'="CRPC3", '1' = "CRPC3", '2'= "CRPC3", '3'="CRPC3",
                      '4'="CRPC3", '5'="CRPC3",'6'="CRPC3",'7'="CRPC3",
                      '8'="CRPC3",'9'="CRPC3",'10'="CRPC3", '11'="CRPC3", '12'="CRPC3", 
                      '13'="CRPC3", '14'="CRPC3", '15'="CRPC3")  
CRPC3[["Original"]] <- Idents(object = CRPC3)

Idents(object = CRPC3) <- "Original"
tiff(file = "CRPC3 orange UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CRPC3, reduction = "umap", pt.size = 0.3, cols = c("orange")) + NoLegend()
dev.off()

##Subset Epi2
Idents(object = epi.integrated1) <- "seurat_clusters"
epi.integrated2 <- subset(epi.integrated1, idents = c("0", "1", "2", "3", "4", "5", 
                                                     "6", "7", "8", "9", "10", "11"
                                                     , "12", "14", "15"))

#Downstream
DefaultAssay(epi.integrated2) <- "integrated"
epi.integrated2 <- ScaleData(epi.integrated2, verbose = FALSE)
epi.integrated2 <- RunPCA(epi.integrated2, npcs = 50, verbose = FALSE)
epi.integrated2 <- RunUMAP(epi.integrated2, dims = 1:20)
epi.integrated2 <- FindNeighbors(epi.integrated2, dims = 1:20)
epi.integrated2 <- FindClusters(epi.integrated2, resolution = 0.5)

Idents(object = epi.integrated2) <- "seurat_clusters"
DimPlot(epi.integrated2, reduction = "umap", pt.size = 0.3, split.by = "stim", label = T)

#FeaturePlots for cell markers
DefaultAssay(epi.integrated2)<-"RNA"
tiff(file = "epi.integrated2 Epi expression plots.tiff", width = 12.5, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = c("AR", "TMPRSS2", "KLK3", "NKX3-1",
                                                              "PSCA", "KRT4",  "TACSTD2", "PIGR",
                                                              "KRT15", "KRT14", "KRT6A", "KRT5", "KRT16", "KRT17", 
                                                              "CD24", 
                                                              "SERPINB1", "LCN2", "KRT13", "KRT19",
                                                              "KRT8", "KRT18"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap for AR signature
DefaultAssay(epi.integrated2)<-"RNA"
ARlist <- list(c("NAP1L2", "TMPRSS2", "CHRNA2", "TARP", "STEAP4", "ALDH1A3", "PART1", "PLPP1", "KLK3", "KLK2", 'NKX3-1', "AR", "PMEPA1", "SLC45A3", "FKBP5"))
epi.integrated2 <- AddModuleScore(object = epi.integrated2, features = ARlist, name = "ARDownstream") 
epi.integrated2[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = epi.integrated2, vars = 'ARDownstream1')))
all.genes <- rownames(epi.integrated2)
epi.integrated2 <- ScaleData(epi.integrated2, features = all.genes)

tiff(file = "epi.integrated2 ARscore expression plots.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'ARDownstream1', cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q20", max.cutoff = "q99")
dev.off()
tiff(file = "epi.integrated2 ARscore split expression plots.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'ARDownstream1', cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q20", max.cutoff = "q99", split.by = "stim")
dev.off()

#Heatmap for NE signature
DefaultAssay(epi.integrated2)<-"RNA"
NElist <- list(c('NKX2-1', "INSM1", "SOX2", "LMO3", "POU3F2", "ASCL1", "CHGA", "CHGB", "SCN3A", "CELF3", "ELAVL4", "PCSK1", "SRRM4", "SNAP25", "SYP", "CHRNB2", "ACTL6B", "SCG3", "ENO2"))
epi.integrated2 <- AddModuleScore(object = epi.integrated2, features = NElist, name = "NEmodule") 
epi.integrated2[['NEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = epi.integrated2, vars = 'NEmodule1')))
all.genes <- rownames(epi.integrated2)
epi.integrated2 <- ScaleData(epi.integrated2, features = all.genes)

tiff(file = "epi.integrated2 NEscore expression plots.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'NEmodule1', cols = c("light grey", "red"), min.cutoff = "q30", pt.size = 0.3)
dev.off()
tiff(file = "epi.integrated2 NEscore split expression plots.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'NEmodule1', cols = c("light grey", "red"), min.cutoff = "q30", pt.size = 0.3, split.by = "stim")
dev.off()

tiff(file = "epi.integrated2 CHGB expression plots.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = "CHGB", cols = c("light grey", "red"), max.cutoff = "q90", pt.size = 0.3)
dev.off()
tiff(file = "epi.integrated2 CHGB split expression plots.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = "CHGB", cols = c("light grey", "red"), max.cutoff = "q90", pt.size = 0.3, split.by = "stim")
dev.off()

#Heatmap for BE signature
BElist <- list(c('KRT5', "KRT15", "KRT13", "SLC14A1", "DST", "LAMB3", "NTN4", "FLRT3", "TIMP3", "TP63", "CLU"))
epi.integrated2 <- AddModuleScore(object = epi.integrated2, features = BElist, name = "BEmodule") 
epi.integrated2[['BEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = epi.integrated2, vars = 'BEmodule1')))
all.genes <- rownames(epi.integrated2)
epi.integrated2 <- ScaleData(epi.integrated2, features = all.genes)

tiff(file = "epi.integrated2 BEscore expression plots.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'BEmodule1', cols = c("light grey", "red"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3)
dev.off()
tiff(file = "epi.integrated2 BEscore split expression plots.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'BEmodule1', cols = c("light grey", "red"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3, split.by = "stim")
dev.off()

#Heatmap for Club signature
Clublist <- list(c('MMP7', "PIGR", "CP", "RARRES1", "IGFBP3", "KRT7", "LCN2", "KRT4", "LXN", "NCOA7", "CEACAM6"))
epi.integrated2 <- AddModuleScore(object = epi.integrated2, features = Clublist, name = "Clubmodule") 
epi.integrated2[['Clubmodule']] <- CreateAssayObject(data = t(x = FetchData(object = epi.integrated2, vars = 'Clubmodule1')))
all.genes <- rownames(epi.integrated2)
epi.integrated2 <- ScaleData(epi.integrated2, features = all.genes)

tiff(file = "epi.integrated2 Clubmodule expression plots.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'Clubmodule1', cols = c("light grey", "red"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3)
dev.off()
tiff(file = "epi.integrated2 Clubmodule split expression plots.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'Clubmodule1', cols = c("light grey", "red"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3, split.by = "stim")
dev.off()


#Heatmap for LE signature
LElist <- list(c('MSMB', "ACPP", "NEFH", "ANPEP", "LSAMP", "AXGP1", "MME", "KRT8", "KRT18", "KLK3", "CD38"))
epi.integrated2 <- AddModuleScore(object = epi.integrated2, features = LElist, name = "LEmodule") 
epi.integrated2[['LEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = epi.integrated2, vars = 'LEmodule1')))
all.genes <- rownames(epi.integrated2)
epi.integrated2 <- ScaleData(epi.integrated2, features = all.genes)

tiff(file = "epi.integrated2 LEmodule expression plots.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'LEmodule1', cols = c("light grey", "red"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3)
dev.off()
tiff(file = "epi.integrated2 LEmodule split expression plots.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(epi.integrated2, reduction = "umap", features = 'LEmodule1', cols = c("light grey", "red"), min.cutoff = "q10", max.cutoff = "q90", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell Type Identification
Idents(object = epi.integrated2) <- "seurat_clusters"
DimPlot(epi.integrated2, reduction = "umap", pt.size = 0.3, split.by = "stim", label = T)

epi.integrated2 <- RenameIdents(object = epi.integrated2, 
                                '4'="BE1", '7' = "BE2", 
                                '5'= "LE1", '13'= "LE2", '12'= "LE3",
                                '3'= "CRPC1", '9'= "CRPC2", '10'= "CRPC3",'11'= "CRPC4",
                                '1'= "CRPC5",'0'= "CRPC6",'6'= "CRPC7",'8'= "CRPC8",
                                '2'= "Club",'14'= "OE")  
epi.integrated2[["Celltype"]] <- Idents(object = epi.integrated2)

Idents(object = epi.integrated2) <- "Celltype"
tiff(file = "epi.integrated2 Celltype UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(epi.integrated2, reduction = "umap", pt.size = 0.3)
dev.off()
Idents(object = epi.integrated2) <- "Celltype"
tiff(file = "epi.integrated2 Celltype split UMAP.tiff", width = 12.5, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(epi.integrated2, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

Idents(object = epi.integrated2) <- "Celltype"
epi.integrated2$stim.Celltype <- paste(Idents(epi.integrated2), epi.integrated2$stim, sep = "_")
Idents(object = epi.integrated2) <- "stim.Celltype"
table(Idents(epi.integrated2))


saveRDS(epi.integrated2, file = "C:/Users/김원경/Desktop/WK Kim Lab/CRPC Project/CRPC.combined.epi/epi.integrated2.rds")


#Heatmap
#Heatmap
DefaultAssay(epi.integrated2) <- "RNA"
Idents(object = epi.integrated2) <- "Celltype"
epi.integrated2 <- ScaleData(epi.integrated2, features = rownames(epi.integrated2))
epi.integrated2.all.markers <- FindAllMarkers(epi.integrated2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
epi.integrated2.all.markers.Top50 <- epi.integrated2.all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
tiff(file = "epi.integrated2 Heatmap Top50 purple.tiff", width = 6, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(epi.integrated2, features = c(epi.integrated2.all.markers.Top50$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

