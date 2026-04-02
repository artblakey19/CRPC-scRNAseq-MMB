####ARKO_P35####

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
library(scDblFinder)

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_ARKO_P35")

ARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/ARKO-Gli1/Single_Cell_Seq/190219_scRNA/27348_ARKO1_count_EGFPmm10/outs/filtered_feature_bc_matrix")
ARKOunfiltered <- CreateSeuratObject(counts = ARKOunfiltered.data,  min.cells = 3, min.features = 500, project = "ARKOunfiltered")
ARKOunfiltered <- NormalizeData(ARKOunfiltered)

ARKOunfiltered[["percent.mt"]] <- PercentageFeatureSet(ARKOunfiltered, pattern = "^mt-")

#new filtering paramaters

table(Idents(ARKOunfiltered))
VlnPlot(ARKOunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
hist(ARKOunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARKO pre filteration")
hist(ARKOunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARKO pre filteration")

ARKO <- subset(ARKOunfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)

table(Idents(ARKO))
VlnPlot(ARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
hist(ARKO@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARKO post filteration")
hist(ARKO@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARKO post filteration")

ARKO <- FindVariableFeatures(ARKO, selection.method = "vst", nfeatures = 5000)

all.genes <- rownames(ARKO)
ARKO <- ScaleData(ARKO, features = all.genes)
ARKO <- RunPCA(ARKO, features = VariableFeatures(object = ARKO))
print(ARKO[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ARKO, dims = 1:2, reduction = "pca")
DimPlot(ARKO, reduction = "pca")
ElbowPlot(ARKO, ndims = 50)
ARKO <- FindNeighbors(ARKO, dims = 1:25)
ARKO <- FindClusters(ARKO, resolution = 0.5)
head(Idents(ARKO), 5)
ARKO <- RunTSNE(ARKO, dims = 1:25)
ARKO <- RunUMAP(ARKO, dims = 1:25)
DimPlot(ARKO, reduction = "umap", pt.size = 1, label = TRUE)

Idents(object = TwoM) <- "seurat_clusters"
Idents(object = WT) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"
tiff(file = "ARKO UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKO, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "WT UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "TwoM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3)
dev.off()

#### Merging Datasets TwoMvWT-P35vWT-P60 ####

#Stash old idents
TwoM[["orig.clusters"]] <- Idents(object = TwoM)
WT[["orig.clusters"]] <- Idents(object = WT)
ARKO[["orig.clusters"]] <- Idents(object = ARKO)

#Set Current idents
Idents(object = TwoM) <- "seurat_clusters"
Idents(object = WT) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"

TwoM$stim <- "ARQ9"
WT$stim <- "WT"
ARKO$stim <- "ARKO"

TwoMvWTvARKO.anchors <- FindIntegrationAnchors(object.list = list(TwoM, WT, ARKO), dims = 1:20)
TwoMvWTvARKO.combined <- IntegrateData(anchorset = TwoMvWTvARKO.anchors, dims = 1:20)

DefaultAssay(TwoMvWTvARKO.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWTvARKO.combined <- ScaleData(TwoMvWTvARKO.combined, verbose = FALSE)
TwoMvWTvARKO.combined <- RunPCA(TwoMvWTvARKO.combined, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWTvARKO.combined, ndims = 50)

#Umap and Clustering
TwoMvWTvARKO.combined <- FindNeighbors(TwoMvWTvARKO.combined, reduction = "pca", dims = 1:16)
TwoMvWTvARKO.combined <- FindClusters(TwoMvWTvARKO.combined, resolution = 0.5)
TwoMvWTvARKO.combined <- RunTSNE(TwoMvWTvARKO.combined, reduction = "pca", dims = 1:16)
TwoMvWTvARKO.combined <- RunUMAP(TwoMvWTvARKO.combined, reduction = "pca", dims = 1:16)
DimPlot(TwoMvWTvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "TwoMvWTvARKO.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWTvARKO.combined) <- "RNA"
all.genes <- rownames(TwoMvWTvARKO.combined)
TwoMvWTvARKO.combined <- ScaleData(TwoMvWTvARKO.combined, features = all.genes)
TwoMvWTvARKO.combined <- CellCycleScoring(TwoMvWTvARKO.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWTvARKO.combined) <- "Phase"
DimPlot(TwoMvWTvARKO.combined, reduction = "umap")
tiff(file = "TwoMvWTvARKO.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWTvARKO.combined1 <- TwoMvWTvARKO.combined
DefaultAssay(TwoMvWTvARKO.combined1) <- "integrated"
TwoMvWTvARKO.combined1 <- ScaleData(TwoMvWTvARKO.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWTvARKO.combined1))
TwoMvWTvARKO.combined1 <- RunPCA(TwoMvWTvARKO.combined1, features = VariableFeatures(TwoMvWTvARKO.combined1))
ElbowPlot(TwoMvWTvARKO.combined1, ndims = 50)

TwoMvWTvARKO.combined1 <- FindNeighbors(TwoMvWTvARKO.combined1, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined1 <- FindClusters(TwoMvWTvARKO.combined1, resolution = 0.5)
TwoMvWTvARKO.combined1 <- RunUMAP(TwoMvWTvARKO.combined1, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined1 <- RunTSNE(TwoMvWTvARKO.combined1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWTvARKO.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWTvARKO.combined1) <- "Phase"
tiff(file = "TwoMvWTvARKO.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
Idents(object = TwoMvWTvARKO.combined1) <- "seurat_clusters"
tiff(file = "TwoMvWTvARKO.combined1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWTvARKO.combined1) <- "stim"
TwoMvWTvARKO.combined1 <- RenameIdents(object = TwoMvWTvARKO.combined1, 'ARQ9' = "ARQ9", 'WT' = "WT", 'ARKO' = "ARKO")
TwoMvWTvARKO.combined1[["stim"]] <- Idents(object = TwoMvWTvARKO.combined1)
Idents(object = TwoMvWTvARKO.combined1) <- "seurat_clusters"
tiff(file = "TwoMvWTvARKO.combined1 seurat split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

DimPlot(TwoMvWTvARKO.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell type identification
DefaultAssay(TwoMvWTvARKO.combined1) <- "RNA"
tiff(file = "TwoMvWTvARKO.combined1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWTvARKO.combined1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2", "Epcam",
                                                                "Vim", "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                                "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####Re-clustering Stro####

Idents(object = TwoMvWTvARKO.combined1) <- "seurat_clusters"
TwoMvWTvARKO.combined.Stro <- subset(TwoMvWTvARKO.combined1, idents = c("13", "16", "4", "17", "3", "1", "11", "9"))
Idents(object = TwoMvWTvARKO.combined.Stro) <- "seurat_clusters"
DefaultAssay(TwoMvWTvARKO.combined.Stro) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWTvARKO.combined.Stro <- ScaleData(TwoMvWTvARKO.combined.Stro, verbose = FALSE)
TwoMvWTvARKO.combined.Stro <- RunPCA(TwoMvWTvARKO.combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWTvARKO.combined.Stro, ndims = 50)

#Umap and Clustering
TwoMvWTvARKO.combined.Stro <- FindNeighbors(TwoMvWTvARKO.combined.Stro, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined.Stro <- FindClusters(TwoMvWTvARKO.combined.Stro, resolution = 0.5)
TwoMvWTvARKO.combined.Stro <- RunTSNE(TwoMvWTvARKO.combined.Stro, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined.Stro <- RunUMAP(TwoMvWTvARKO.combined.Stro, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWTvARKO.combined.Stro, reduction = "umap", pt.size = 0.3)

tiff(file = "TwoMvWTvARKO.combined.Stro UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.Stro, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWTvARKO.combined.Stro) <- "RNA"
all.genes <- rownames(TwoMvWTvARKO.combined.Stro)
TwoMvWTvARKO.combined.Stro <- ScaleData(TwoMvWTvARKO.combined.Stro, features = all.genes)
TwoMvWTvARKO.combined.Stro <- CellCycleScoring(TwoMvWTvARKO.combined.Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWTvARKO.combined.Stro) <- "Phase"
DimPlot(TwoMvWTvARKO.combined.Stro, reduction = "umap")
tiff(file = "TwoMvWTvARKO.combined.Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.Stro, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWTvARKO.combined.Stro1 <- TwoMvWTvARKO.combined.Stro
DefaultAssay(TwoMvWTvARKO.combined.Stro1) <- "integrated"
TwoMvWTvARKO.combined.Stro1 <- ScaleData(TwoMvWTvARKO.combined.Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWTvARKO.combined.Stro1))
TwoMvWTvARKO.combined.Stro1 <- RunPCA(TwoMvWTvARKO.combined.Stro1, features = VariableFeatures(TwoMvWTvARKO.combined.Stro1))
ElbowPlot(TwoMvWTvARKO.combined.Stro1, ndims = 50)

TwoMvWTvARKO.combined.Stro1 <- FindNeighbors(TwoMvWTvARKO.combined.Stro1, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined.Stro1 <- FindClusters(TwoMvWTvARKO.combined.Stro1, resolution = 0.5)
TwoMvWTvARKO.combined.Stro1 <- RunUMAP(TwoMvWTvARKO.combined.Stro1, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined.Stro1 <- RunTSNE(TwoMvWTvARKO.combined.Stro1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWTvARKO.combined.Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWTvARKO.combined.Stro1) <- "Phase"
tiff(file = "TwoMvWTvARKO.combined.Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.Stro1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
Idents(object = TwoMvWTvARKO.combined.Stro1) <- "seurat_clusters"
tiff(file = "TwoMvWTvARKO.combined.Stro1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

#Rename
DefaultAssay(TwoMvWTvARKO.combined.Stro1) <- "RNA"
tiff(file = "TwoMvWTvARKO.combined.Stro1 expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWTvARKO.combined.Stro1, reduction = "umap", features = c("Fbln1", "Myh11", "Pecam1", "Plp1", "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 1)
dev.off()

Idents(object = TwoMvWTvARKO.combined.Stro1) <- "seurat_clusters"
TwoMvWTvARKO.combined.Stro1 <- RenameIdents(object = TwoMvWTvARKO.combined.Stro1, '5' = "FB1", '0' = "FB2", '4' = "FB3", '2' = "FB4", '11' = "FB5",
                                       '1' = "SM1", '3' = "SM2", '6' = "VE", '15' = "Pericyte", '9' = "Glia", '13' = "Glia", '10' = "Leu",
                                       '8' = "Leu", '12' = "Lym", '7' = "Lym", '14' = "OF")
TwoMvWTvARKO.combined.Stro1[["StroCellType"]] <- Idents(object = TwoMvWTvARKO.combined.Stro1)

#Umap
Idents(object = TwoMvWTvARKO.combined.Stro1) <- "StroCellType"
tiff(file = "TwoMvWTvARKO.combined.Stro1 StroCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

#FB only
TwoMvWTvARKO.combined.FB <- subset(TwoMvWTvARKO.combined.Stro1, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
TwoMvWTvARKO.combined.FB <- RunTSNE(TwoMvWTvARKO.combined.FB, reduction = "pca", dims = 1:20)
TwoMvWTvARKO.combined.FB <- RunUMAP(TwoMvWTvARKO.combined.FB, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWTvARKO.combined.FB, reduction = "umap", pt.size = 0.3)


tiff(file = "TwoMvWTvARKO.combined.FB UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.FB, reduction = "umap", pt.size = 1, label = TRUE)
dev.off()

tiff(file = "TwoMvWTvARKO.combined.FB split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWTvARKO.combined.FB, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWTvARKO.combined.FB) <- "stim"
TwoMvWTvARKO.combined.FB$stim.StroCellType <- paste(Idents(TwoMvWTvARKO.combined.FB), TwoMvWTvARKO.combined.FB$StroCellType, sep = "_")
Idents(object = TwoMvWTvARKO.combined.FB) <- "stim.StroCellType"
table(Idents(TwoMvWTvARKO.combined.FB))

#Add EGFP info
DefaultAssay(TwoMvWTvARKO.combined.FB) <- "RNA"
TwoMvWTvARKO.combined.FBEGFPNeg <- subset(x=TwoMvWTvARKO.combined.FB, subset = EGFP < 1)
TwoMvWTvARKO.combined.FBEGFPPos <- subset(x=TwoMvWTvARKO.combined.FB, subset = EGFP >= 1)
Idents(object = TwoMvWTvARKO.combined.FBEGFPNeg) <- "EGFPNeg"
Idents(object = TwoMvWTvARKO.combined.FBEGFPPos) <- "EGFPPos"
TwoMvWTvARKO.combined.FBEGFPNeg[["EGFPExp"]] <- Idents(object = TwoMvWTvARKO.combined.FBEGFPNeg)
TwoMvWTvARKO.combined.FBEGFPPos[["EGFPExp"]] <- Idents(object = TwoMvWTvARKO.combined.FBEGFPPos)
TwoMvWTvARKO.combined.FBEGFP <- merge(x = TwoMvWTvARKO.combined.FBEGFPNeg, y = TwoMvWTvARKO.combined.FBEGFPPos)
Idents(object = TwoMvWTvARKO.combined.FBEGFP) <- "EGFPExp"
TwoMvWTvARKO.combined.FB$EGFPExp <- Idents(object = TwoMvWTvARKO.combined.FBEGFP)

#Add ARQ info
DefaultAssay(TwoMvWTvARKO.combined.FB) <- "RNA"
TwoMvWTvARKO.combined.FBARQNeg <- subset(x=TwoMvWTvARKO.combined.FB, subset = ARQ == 0)
TwoMvWTvARKO.combined.FBARQPos <- subset(x=TwoMvWTvARKO.combined.FB, subset = ARQ > 0)
Idents(object = TwoMvWTvARKO.combined.FBARQNeg) <- "ARQNeg"
Idents(object = TwoMvWTvARKO.combined.FBARQPos) <- "ARQPos"
TwoMvWTvARKO.combined.FBARQNeg[["ARQExp"]] <- Idents(object = TwoMvWTvARKO.combined.FBARQNeg)
TwoMvWTvARKO.combined.FBARQPos[["ARQExp"]] <- Idents(object = TwoMvWTvARKO.combined.FBARQPos)
TwoMvWTvARKO.combined.FBARQ <- merge(x = TwoMvWTvARKO.combined.FBARQNeg, y = TwoMvWTvARKO.combined.FBARQPos)
Idents(object = TwoMvWTvARKO.combined.FBARQ) <- "ARQExp"
TwoMvWTvARKO.combined.FB$ARQExp <- Idents(object = TwoMvWTvARKO.combined.FBARQ)

#Cellcounts
Idents(object = TwoMvWTvARKO.combined.FB) <- "EGFPExp."
TwoMvWTvARKO.combined.FB$EGFPExp.StroCellType <- paste(Idents(TwoMvWTvARKO.combined.FB), TwoMvWTvARKO.combined.FB$StroCellType, sep = "_")
Idents(object = TwoMvWTvARKO.combined.FB) <- "EGFPExp"
TwoMvWTvARKO.combined.FB$EGFPExp.ARQExp <- paste(Idents(TwoMvWTvARKO.combined.FB), TwoMvWTvARKO.combined.FB$ARQExp, sep = "_")
Idents(object = TwoMvWTvARKO.combined.FB) <- "EGFPExp.ARQExp"
TwoMvWTvARKO.combined.FB$EGFPExp.ARQExp.StroCellType <- paste(Idents(TwoMvWTvARKO.combined.FB), TwoMvWTvARKO.combined.FB$StroCellType, sep = "_")
Idents(object = TwoMvWTvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType"
TwoMvWTvARKO.combined.FB$EGFPExp.ARQExp.StroCellType.stim <- paste(Idents(TwoMvWTvARKO.combined.FB), TwoMvWTvARKO.combined.FB$stim, sep = "_")
Idents(object = TwoMvWTvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"

Idents(object = TwoMvWTvARKO.combined.FB) <- "stim"
TwoMvARKO.combined.FB <- subset(TwoMvWTvARKO.combined.FB, idents = c("ARQ9", "ARKO"))

#FB1vARKO
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.EGFPFB <- subset(TwoMvARKO.combined.FB, idents = c("EGFPPos_ARQPos_FB1_ARQ9", "EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"))
TwoMvARKO.combined.EGFPFB <- RenameIdents(object = TwoMvARKO.combined.EGFPFB, 'EGFPPos_ARQPos_FB1_ARQ9' = "EGFPPos_ARQPos_FB1_ARQ9", 'EGFPPos_ARQNeg_FB1_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB2_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB3_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB4_ARKO' = "EGFPPos_ARKO"
                                          , 'EGFPPos_ARQNeg_FB5_ARKO' = "EGFPPos_ARKO")
TwoMvARKO.combined.EGFPFB[["FB1vARKO"]] <- Idents(object = TwoMvARKO.combined.EGFPFB)

DefaultAssay(TwoMvARKO.combined.EGFPFB) <- "RNA"
Idents(object = TwoMvARKO.combined.EGFPFB) <- "FB1vARKO"
TwoMvARKO.combined.EGFPFB <- ScaleData(TwoMvARKO.combined.EGFPFB, features = rownames(TwoMvARKO.combined.EGFPFB))
TwoMvARKO.combined.EGFPFB.markers <- FindAllMarkers(TwoMvARKO.combined.EGFPFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvARKO.combined.EGFPFB.markers.Top50 <- TwoMvARKO.combined.EGFPFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB1vsARKO_mGFP+_All Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvARKO.combined.EGFPFB, features = c(TwoMvARKO.combined.EGFPFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#FB1vARKO
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.EGFPFB <- subset(TwoMvARKO.combined.FB, idents = c("EGFPPos_ARQPos_FB1_ARQ9", "EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"))
TwoMvARKO.combined.EGFPFB <- RenameIdents(object = TwoMvARKO.combined.EGFPFB, 'EGFPPos_ARQPos_FB1_ARQ9' = "EGFPPos_ARQPos_FB1_ARQ9", 'EGFPPos_ARQNeg_FB1_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB2_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB3_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB4_ARKO' = "EGFPPos_ARKO"
                                          , 'EGFPPos_ARQNeg_FB5_ARKO' = "EGFPPos_ARKO")
TwoMvARKO.combined.EGFPFB[["FB1vARKO"]] <- Idents(object = TwoMvARKO.combined.EGFPFB)

DefaultAssay(TwoMvARKO.combined.EGFPFB) <- "RNA"
Idents(object = TwoMvARKO.combined.EGFPFB) <- "FB1vARKO"
TwoMvARKO.combined.EGFPFB <- ScaleData(TwoMvARKO.combined.EGFPFB, features = rownames(TwoMvARKO.combined.EGFPFB))
TwoMvARKO.combined.EGFPFB.markers <- FindAllMarkers(TwoMvARKO.combined.EGFPFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvARKO.combined.EGFPFB.markers.Top50 <- TwoMvARKO.combined.EGFPFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB1vsARKO_mGFP+_All Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvARKO.combined.EGFPFB, features = c(TwoMvARKO.combined.EGFPFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#FB2vARKO
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.EGFPFB <- subset(TwoMvARKO.combined.FB, idents = c("EGFPPos_ARQPos_FB2_ARQ9", "EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"))
TwoMvARKO.combined.EGFPFB <- RenameIdents(object = TwoMvARKO.combined.EGFPFB, 'EGFPPos_ARQPos_FB2_ARQ9' = "EGFPPos_ARQPos_FB2_ARQ9", 'EGFPPos_ARQNeg_FB1_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB2_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB3_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB4_ARKO' = "EGFPPos_ARKO"
                                          , 'EGFPPos_ARQNeg_FB5_ARKO' = "EGFPPos_ARKO")
TwoMvARKO.combined.EGFPFB[["FB2vARKO"]] <- Idents(object = TwoMvARKO.combined.EGFPFB)

DefaultAssay(TwoMvARKO.combined.EGFPFB) <- "RNA"
Idents(object = TwoMvARKO.combined.EGFPFB) <- "FB2vARKO"
TwoMvARKO.combined.EGFPFB <- ScaleData(TwoMvARKO.combined.EGFPFB, features = rownames(TwoMvARKO.combined.EGFPFB))
TwoMvARKO.combined.EGFPFB.markers <- FindAllMarkers(TwoMvARKO.combined.EGFPFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvARKO.combined.EGFPFB.markers.Top50 <- TwoMvARKO.combined.EGFPFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB2vsARKO_mGFP+_All Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvARKO.combined.EGFPFB, features = c(TwoMvARKO.combined.EGFPFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#FB3vARKO
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.EGFPFB <- subset(TwoMvARKO.combined.FB, idents = c("EGFPPos_ARQPos_FB3_ARQ9", "EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"))
TwoMvARKO.combined.EGFPFB <- RenameIdents(object = TwoMvARKO.combined.EGFPFB, 'EGFPPos_ARQPos_FB3_ARQ9' = "EGFPPos_ARQPos_FB3_ARQ9", 'EGFPPos_ARQNeg_FB1_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB2_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB3_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB4_ARKO' = "EGFPPos_ARKO"
                                          , 'EGFPPos_ARQNeg_FB5_ARKO' = "EGFPPos_ARKO")
TwoMvARKO.combined.EGFPFB[["FB3vARKO"]] <- Idents(object = TwoMvARKO.combined.EGFPFB)

DefaultAssay(TwoMvARKO.combined.EGFPFB) <- "RNA"
Idents(object = TwoMvARKO.combined.EGFPFB) <- "FB3vARKO"
TwoMvARKO.combined.EGFPFB <- ScaleData(TwoMvARKO.combined.EGFPFB, features = rownames(TwoMvARKO.combined.EGFPFB))
TwoMvARKO.combined.EGFPFB.markers <- FindAllMarkers(TwoMvARKO.combined.EGFPFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvARKO.combined.EGFPFB.markers.Top50 <- TwoMvARKO.combined.EGFPFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB3vsARKO_mGFP+_All Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvARKO.combined.EGFPFB, features = c(TwoMvARKO.combined.EGFPFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#FB4vARKO
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.EGFPFB <- subset(TwoMvARKO.combined.FB, idents = c("EGFPPos_ARQPos_FB4_ARQ9", "EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"))
TwoMvARKO.combined.EGFPFB <- RenameIdents(object = TwoMvARKO.combined.EGFPFB, 'EGFPPos_ARQPos_FB4_ARQ9' = "EGFPPos_ARQPos_FB4_ARQ9", 'EGFPPos_ARQNeg_FB1_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB2_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB3_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB4_ARKO' = "EGFPPos_ARKO"
                                          , 'EGFPPos_ARQNeg_FB5_ARKO' = "EGFPPos_ARKO")
TwoMvARKO.combined.EGFPFB[["FB4vARKO"]] <- Idents(object = TwoMvARKO.combined.EGFPFB)

DefaultAssay(TwoMvARKO.combined.EGFPFB) <- "RNA"
Idents(object = TwoMvARKO.combined.EGFPFB) <- "FB4vARKO"
TwoMvARKO.combined.EGFPFB <- ScaleData(TwoMvARKO.combined.EGFPFB, features = rownames(TwoMvARKO.combined.EGFPFB))
TwoMvARKO.combined.EGFPFB.markers <- FindAllMarkers(TwoMvARKO.combined.EGFPFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvARKO.combined.EGFPFB.markers.Top50 <- TwoMvARKO.combined.EGFPFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB4vsARKO_mGFP+_All Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvARKO.combined.EGFPFB, features = c(TwoMvARKO.combined.EGFPFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#FB5vARKO
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.EGFPFB <- subset(TwoMvARKO.combined.FB, idents = c("EGFPPos_ARQPos_FB5_ARQ9", "EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"))
TwoMvARKO.combined.EGFPFB <- RenameIdents(object = TwoMvARKO.combined.EGFPFB, 'EGFPPos_ARQPos_FB5_ARQ9' = "EGFPPos_ARQPos_FB5_ARQ9", 'EGFPPos_ARQNeg_FB1_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB2_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB3_ARKO' = "EGFPPos_ARKO", 'EGFPPos_ARQNeg_FB4_ARKO' = "EGFPPos_ARKO"
                                          , 'EGFPPos_ARQNeg_FB5_ARKO' = "EGFPPos_ARKO")
TwoMvARKO.combined.EGFPFB[["FB5vARKO"]] <- Idents(object = TwoMvARKO.combined.EGFPFB)

DefaultAssay(TwoMvARKO.combined.EGFPFB) <- "RNA"
Idents(object = TwoMvARKO.combined.EGFPFB) <- "FB5vARKO"
TwoMvARKO.combined.EGFPFB <- ScaleData(TwoMvARKO.combined.EGFPFB, features = rownames(TwoMvARKO.combined.EGFPFB))
TwoMvARKO.combined.EGFPFB.markers <- FindAllMarkers(TwoMvARKO.combined.EGFPFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvARKO.combined.EGFPFB.markers.Top50 <- TwoMvARKO.combined.EGFPFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB5vsARKO_mGFP+_All Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvARKO.combined.EGFPFB, features = c(TwoMvARKO.combined.EGFPFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#DEGs_TwoM_mGFP+ARQ+FB2vsmGFP+ARQ+FB1-5

DefaultAssay(TwoMvARKO.combined.FB) <- "RNA"
Idents(object = TwoMvARKO.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
TwoMvARKO.combined.FB <- ScaleData(TwoMvARKO.combined.FB, features = rownames(TwoMvARKO.combined.FB))
TwoM_EGFPPosARQPos_FBvARKO_EGFPPos_FB.Markers <- FindMarkers(TwoMvARKO.combined.FB, ident.1 = c("EGFPPos_ARQPos_FB1_ARQ9", "EGFPPos_ARQPos_FB2_ARQ9", "EGFPPos_ARQPos_FB3_ARQ9", "EGFPPos_ARQPos_FB4_ARQ9", "EGFPPos_ARQPos_FB5_ARQ9"), 
                                                      ident.2 = c("EGFPPos_ARQNeg_FB1_ARKO", "EGFPPos_ARQNeg_FB2_ARKO", "EGFPPos_ARQNeg_FB3_ARKO", "EGFPPos_ARQNeg_FB4_ARKO", "EGFPPos_ARQNeg_FB5_ARKO"), min.pct = 0.1, logfc.threshold = 0.1)
write.csv(TwoM_EGFPPosARQPos_FBvARKO_EGFPPos_FB.Markers, "TwoM_EGFPPosARQPos_FBvARKO_EGFPPos_FB.Markers.csv")
