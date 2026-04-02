#ARQ-Osr1 SCseq Work Flow

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
library(ggpubr)
library(GGally)

#### Initial Filtering and Clustering ####

#Setup workspace to make file calling & saving easy
setwd("W:/ARQ9_scRNAseq_WK/AllNew_UMAP")

###PIN### 

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

PIN.data <- Read10X("//isi-dcnl/user_data/zjsun/BIC/AR Transgene Tumorigenesis/ARQ9-Osr1 SingleCell Seq/PCa_AR_transgene_analysis/PIN4_outs/filtered_feature_bc_matrix")
PIN <- CreateSeuratObject(counts = PIN.data,  min.cells = 3, min.features = 500, project = "PIN")
PIN <- NormalizeData(PIN)

#Initial processing & filtering

PIN[["percent.mt"]] <- PercentageFeatureSet(PIN, pattern = "^mt-")

tiff(file = "PIN Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(PIN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "PIN Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Pre-filteration")
dev.off()
tiff(file = "PIN Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Pre-filteration")
dev.off()

PIN2 <- subset(PIN, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
tiff(file = "PIN Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(PIN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "PIN Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN2@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Post-filteration")
dev.off()
tiff(file = "PIN Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN2@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Post-filteration")
dev.off()

PIN2 <- FindVariableFeatures(PIN2, selection.method = "vst", nfeatures = 5000)
tiff(file = "PIN Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(PIN2)
dev.off()

#Clustering
PIN2 <- ScaleData(PIN2, verbose = FALSE)
PIN2 <- RunPCA(PIN2, npcs = 50, verbose = FALSE)
tiff(file = "PIN ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(PIN2, ndims = 50)
dev.off()

PIN2 <- FindNeighbors(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- FindClusters(PIN2, resolution = 0.5)
PIN2 <- JackStraw(PIN2, num.replicate = 100)
PIN2 <- ScoreJackStraw(PIN2, dims = 1:20)
JackStrawPlot(PIN2, dims = 1:20)
PIN2 <- RunTSNE(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- RunUMAP(PIN2, reduction = "pca", dims = 1:20)
tiff(file = "PIN2 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PIN2, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("#B71800")) 
dev.off()

###Tumor###

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

Tumor.data <- Read10X("//isi-dcnl/user_data/zjsun/BIC/AR Transgene Tumorigenesis/ARQ9-Osr1 SingleCell Seq/PCa_AR_transgene_analysis/Tumor4_outs/filtered_feature_bc_matrix")
Tumor <- CreateSeuratObject(counts = Tumor.data,  min.cells = 3, min.features = 500, project = "Tumor")
Tumor <- NormalizeData(Tumor)

#Initial processing & filtering

Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^mt-")

tiff(file = "Tumor Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "Tumor Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Pre-filteration")
dev.off()
tiff(file = "Tumor Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Pre-filteration")
dev.off()

Tumor2 <- subset(Tumor, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
tiff(file = "Tumor Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Tumor2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "Tumor Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Post-filteration")
dev.off()
tiff(file = "Tumor Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Post-filteration")
dev.off()

Tumor2 <- FindVariableFeatures(Tumor2, selection.method = "vst", nfeatures = 5000)
tiff(file = "Tumor Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(Tumor2)
dev.off()

Tumor2 <- ScaleData(Tumor2, verbose = FALSE)
Tumor2 <- RunPCA(Tumor2, npcs = 50, verbose = FALSE)
tiff(file = "Tumor ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(Tumor2, ndims = 50)
dev.off()

Tumor2 <- FindNeighbors(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- FindClusters(Tumor2, resolution = 0.5)
Tumor2 <- JackStraw(Tumor2, num.replicate = 100)
Tumor2 <- ScoreJackStraw(Tumor2, dims = 1:20)
JackStrawPlot(Tumor2, dims = 1:20)
Tumor2 <- RunTSNE(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- RunUMAP(Tumor2, reduction = "pca", dims = 1:20)
tiff(file = "Tumor2 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Tumor2, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("#009FB7")) 
dev.off()

#### Merging Datasets ####

#Stash old idents
Tumor2[["orig.clusters"]] <- Idents(object = Tumor2)
PIN2[["orig.clusters"]] <- Idents(object = PIN2)

#Set Current idents
Idents(object = PIN2) <- "seurat_clusters"
Idents(object = Tumor2) <- "seurat_clusters"
PIN2$stim <- "PIN"
Tumor2$stim <- "Tumor"
PINvTumor.anchors <- FindIntegrationAnchors(object.list = list(PIN2, Tumor2), dims = 1:20)
PINvTumor.combined <- IntegrateData(anchorset = PINvTumor.anchors, dims = 1:20)
DefaultAssay(PINvTumor.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
PINvTumor.combined <- ScaleData(PINvTumor.combined, verbose = FALSE)
PINvTumor.combined <- RunPCA(PINvTumor.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
PINvTumor.combined <- FindNeighbors(PINvTumor.combined, reduction = "pca", dims = 1:20)
PINvTumor.combined <- FindClusters(PINvTumor.combined, resolution = 0.5)
PINvTumor.combined <- RunTSNE(PINvTumor.combined, reduction = "pca", dims = 1:20)
PINvTumor.combined <- RunUMAP(PINvTumor.combined, reduction = "pca", dims = 1:20)
tiff(file = "PINvTumor.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("#B71800", "#009FB7")) 
dev.off()

#New labeling
new.cluster.ids <- c("LE", "BE", "LE", "LE", "LE", "LE", "LE", "Lym", "LE", "Fib", "LE", "LE", "Leu", "Lym", "Endo", "SM", "BE", "BE")
names(new.cluster.ids) <- levels(PINvTumor.combined)
PINvTumor.combined <- RenameIdents(PINvTumor.combined, new.cluster.ids)
PINvTumor.combined[["CellType"]] <- Idents(object = PINvTumor.combined)
tiff(file = "PINvTumor.combined CellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined, reduction = "umap", pt.size = 0.3) 
dev.off()

Idents(object = PINvTumor.combined) <- "stim"
PINonly <- subset(PINvTumor.combined, idents = c("PIN"))
Idents(object = PINonly) <- "CellType"
tiff(file = "PINonly UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINonly, reduction = "umap", pt.size = 0.3) 
dev.off()

Idents(object = PINvTumor.combined) <- "stim"
Tumoronly <- subset(PINvTumor.combined, idents = c("Tumor"))
Idents(object = PINonly) <- "CellType"
tiff(file = "Tumoronly UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Tumoronly, reduction = "umap", pt.size = 0.3) 
dev.off()

#FeaturePlot
DefaultAssay(PINvTumor.combined) <- "RNA"
tiff(file = "ARQ combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "EGFP combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#CoexpressionPlot
DefaultAssay(PINonly) <- "RNA"
tiff(file = "PINonly EGFP hAR Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "ARQ"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "PINonly EGFP Epcam Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Epcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "PINonly EGFP Vim Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Vim"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "PINonly EGFP Myh11 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Myh11"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "PINonly EGFP Krt5 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Krt5"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "PINonly EGFP Krt8 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Krt8"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

DefaultAssay(Tumoronly) <- "RNA"
tiff(file = "Tumoronly EGFP hAR Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Tumoronly, reduction = "umap", features = c("EGFP", "ARQ"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "Tumoronly EGFP Epcam Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Tumoronly, reduction = "umap", features = c("EGFP", "Epcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "Tumoronly EGFP Vim Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Tumoronly, reduction = "umap", features = c("EGFP", "Vim"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "Tumoronly EGFP Myh11 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Tumoronly, reduction = "umap", features = c("EGFP", "Myh11"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "Tumoronly EGFP Krt5 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Tumoronly, reduction = "umap", features = c("EGFP", "Krt5"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "Tumoronly EGFP Krt8 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Tumoronly, reduction = "umap", features = c("EGFP", "Krt8"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

#Dotplot
Idents(object = PINvTumor.combined) <- "CellType"
PINvTumor.combined <- RenameIdents(object = PINvTumor.combined, 'BE' = "BE", 'LE' = "LE", 'Fib' = "Fib", 'SM' = "SM", 'Leu' = "Leu", 'Endo' = "Endo", 'Lym' = "Lym")
DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
tiff(file = "Combined DotPlot.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Cell counts
Idents(object = PINvTumor.combined) <- "stim"
PINvTumor.combined$stim.seurat_clusters <- paste(Idents(PINvTumor.combined), PINvTumor.combined$seurat_clusters, sep = "_")
Idents(object = PINvTumor.combined) <- "stim.seurat_clusters"
table(Idents(PINvTumor.combined))

#### Reclustering Epithelial cells ####

Idents(object = PINvTumor.combined) <- "CellType"
PINvTumor.combined.Epi <- subset(PINvTumor.combined, idents = c("BE", "LE"))
Idents(object = PINvTumor.combined.Epi) <- "seurat_clusters"
DefaultAssay(PINvTumor.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
PINvTumor.combined.Epi <- ScaleData(PINvTumor.combined.Epi, verbose = FALSE)
PINvTumor.combined.Epi <- RunPCA(PINvTumor.combined.Epi, npcs = 30, verbose = FALSE)

#Umap and Clustering
PINvTumor.combined.Epi <- FindNeighbors(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi <- FindClusters(PINvTumor.combined.Epi, resolution = 0.5)
PINvTumor.combined.Epi <- RunTSNE(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi <- RunUMAP(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)

tiff(file = "PINvTumor.combined.Epi UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
dev.off()

Idents(object = PINvTumor.combined.Epi) <- "stim"
PINonlyEpi <- subset(PINvTumor.combined.Epi, idents = c("PIN"))
Idents(object = PINonlyEpi) <- "CellType"
tiff(file = "PINonlyEpi UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINonlyEpi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
dev.off()

TumoronlyEpi <- subset(PINvTumor.combined.Epi, idents = c("Tumor"))
Idents(object = TumoronlyEpi) <- "CellType"
tiff(file = "TumoronlyEpi UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TumoronlyEpi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
dev.off()

#DEGs
DefaultAssay(PINvTumor.combined.Epi) <- "RNA"
Idents(object = PINvTumor.combined.Epi) <- "CellType"
PINvTumor.combined.Epi <- ScaleData(PINvTumor.combined.Epi, features = rownames(PINvTumor.combined.Epi))
PINvTumor.combined.Epi.markers <- FindAllMarkers(PINvTumor.combined.Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(PINvTumor.combined.Epi.Markers, "PINvTumor.combined.Epi.Markers.csv")
PINvTumor.combined.Epi.ARQ.EGFP.Ar.Pbsn.markers <- FindAllMarkers(PINvTumor.combined.Epi, features = c("Pbsn", "Ar", "EGFP", "ARQ"), logfc.threshold = 0)
write.csv(PINvTumor.combined.Epi.ARQ.EGFP.Ar.Pbsn.markers, "PINvTumor.combined.Epi.ARQ.EGFP.Ar.Pbsn.marker.csv")

#DotPlot
Idents(object = PINvTumor.combined.Epi) <- "seurat_clusters"
PINvTumor.combined.Epi <- RenameIdents(object = PINvTumor.combined.Epi, '4' = "BE1", '7' = "BE2", '0' = "LE1", '2' = "LE2", '6' = "LE3", '1' = "LE4", '8' = "LE5", '3' = "LE6", '5' = "LE7", '10' = "LE8", '9' = "UrLE", '11' = "OE")
PINvTumor.combined.Epi[["CellType"]] <- Idents(object = PINvTumor.combined.Epi)
tiff(file = "Epi CellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined.Epi, features = c("Cd53", "Sh2d2a", "Cytip", "Srgn", "Rgs1", "Myof", "Mgat4a", "Bace2", "Cyba", "Wfdc2", "Gjb2", "Timp4", "Cxcl17", "Oit1", "Ppp1r1b", "Pnliprp1", "C1s2", "C1rb", "Pbsn", "Tgm4", "Pmaip1", "Crip1", "Cd55", "Kctd14", "Sbspon", "Aldh1a1", "Hmgcs2", "Ceacam2", "Mme", "Apof", "Klk1b24", "Ptgds","Gpx3", "Svs3a", "Wfdc15b", "Cdca3", "Birc5", "Ube2c", "Stmn1", "Mki67", "Apoc4", "Nkain4", "Gulo", "Npl", "Mt3", "Nxf7", "Syngr1", "Msx2", "Defa20", "Wif1", "Lamb3", "Htra1", "Lbp", "Ltbp4", "Sult5a1", "Col17a1", "Cldn1", "Tpm2", "Krt17", "Krt14", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Cell counts
Idents(object = PINvTumor.combined.Epi) <- "stim"
PINvTumor.combined.Epi$stim.seurat_clusters <- paste(Idents(PINvTumor.combined.Epi), PINvTumor.combined.Epi$seurat_clusters, sep = "_")
Idents(object = PINvTumor.combined.Epi) <- "stim.seurat_clusters"
table(Idents(PINvTumor.combined.Epi))

#FeaturePlots
DefaultAssay(PINonlyEpi) <- "RNA"
tiff(file = "PINonlyEpi hAR Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi mGFP Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi Krt5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi Krt8 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
DefaultAssay(TumoronlyEpi) <- "RNA"
tiff(file = "TumoronlyEpi hAR Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi mGFP Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi Krt5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi Krt8 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####Cell Cycle Regression####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/PINvTumor.combined.Epi")

DefaultAssay(PINvTumor.combined.Epi) <- "RNA"
all.genes <- rownames(PINvTumor.combined.Epi)
PINvTumor.combined.Epi <- ScaleData(PINvTumor.combined.Epi, features = all.genes)
PINvTumor.combined.Epi <- CellCycleScoring(PINvTumor.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = PINvTumor.combined.Epi) <- "Phase"
DimPlot(PINvTumor.combined.Epi, reduction = "umap")
tiff(file = "PINvTumor.combined.Epi Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
tiff(file = "PINvTumor.combined.Epi Cell Cyle stim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
PINvTumor.combined.Epi1 <- PINvTumor.combined.Epi
DefaultAssay(PINvTumor.combined.Epi1) <- "integrated"
PINvTumor.combined.Epi1 <- ScaleData(PINvTumor.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(PINvTumor.combined.Epi1))
PINvTumor.combined.Epi1 <- RunPCA(PINvTumor.combined.Epi1, features = VariableFeatures(PINvTumor.combined.Epi1))
ElbowPlot(PINvTumor.combined.Epi1, ndims = 30)
PINvTumor.combined.Epi1 <- FindNeighbors(PINvTumor.combined.Epi1, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi1 <- FindClusters(PINvTumor.combined.Epi1, resolution = 0.5)
PINvTumor.combined.Epi1 <- RunUMAP(PINvTumor.combined.Epi1, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi1 <- RunTSNE(PINvTumor.combined.Epi1, reduction = "pca", dims = 1:20)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = PINvTumor.combined.Epi1) <- "Phase"
tiff(file = "PINvTumor.combined.Epi1 Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
tiff(file = "PINvTumor.combined.Epi1 Cell Cyle stim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = PINvTumor.combined.Epi1) <- "seurat_clusters"
tiff(file = "PINvTumor.combined.Epi1 seurat UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "PINvTumor.combined.Epi1 ProE highlight after cellcycle regression UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "grey", "grey", "purple", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

tiff(file = "PINvTumor.combined.Epi ProE highlight before cellcycle regression UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "grey", "grey", "purple", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

####Cell Type Identification####

#DEGs
DefaultAssay(PINvTumor.combined.Epi1) <- "RNA"
Idents(object = PINvTumor.combined.Epi1) <- "seurat_clusters"
PINvTumor.combined.Epi1 <- ScaleData(PINvTumor.combined.Epi1, features = rownames(PINvTumor.combined.Epi1))
PINvTumor.combined.Epi1.markers <- FindAllMarkers(PINvTumor.combined.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(PINvTumor.combined.Epi1.markers, "PINvTumor.combined.Epi1.markers.csv")

#Rename
Idents(object = PINvTumor.combined.Epi1) <- "seurat_clusters"
PINvTumor.combined.Epi1 <- RenameIdents(object = PINvTumor.combined.Epi1, '4' = "BE1", '6' = "BE2", '0' = "LE1", '2' = "LE2", '1' = "LE3", '9' = "LE4", '3' = "LE5", '5' = "LE6", '8' = "LE7", '7' = "UrLE", '10' = "OE")
PINvTumor.combined.Epi1[["EpiCellType"]] <- Idents(object = PINvTumor.combined.Epi1)

#Umap
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
tiff(file = "PINvTumor.combined.Epi1 UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
dev.off()
tiff(file = "PINvTumor.combined.Epi1 stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
dev.off()

#Dotplots
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
tiff(file = "Epi1 EpiCellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined.Epi1, features = c("Ar", "ARQ", "EGFP", "Krt14", "Krt17", "Tpm2", "Krt16", "Col17a1", "Sult5a1", "Ltbp4", "Lbp", "Htra1", "Lamb3", "Wif1", "Defa20", "Defa22", "Msx2", "Syngr1", "Mt3", "Npl", "Gulo", "Mt4", "Nkain4", "Svs3a", "Wfdc15b", "Gpx3", "Ptgds", "Klk1b24", "Pcp4", "Smgc", "Snhg11", "Itln1", "Ceacam2", "Crip1", "Cd55", "Pmaip1", "Sbspon", "Kctd14", "Tgm4", "C1rb", "Gm5615", "Pnliprp1", "C1s2", "Oit1", "Cxcl17", "Timp4", "Aqp5", "Gjb2", "Wfdc2", "Msln", "Gsdmc2", "Atp10b", "Mgat4a", "Rgs1", "Srgn", "Tnfrsf9", "Cytip", "Coro1a"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

tiff(file = "Epi1 EpiCellType hAR GFP Ar Pbsn DotPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined.Epi1, features = c("ARQ", "EGFP", "Ar", "Pbsn"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Cell counts
Idents(object = PINvTumor.combined.Epi1) <- "stim"
PINvTumor.combined.Epi1$stim.EpiCellType <- paste(Idents(PINvTumor.combined.Epi1), PINvTumor.combined.Epi1$EpiCellType, sep = "_")
Idents(object = PINvTumor.combined.Epi1) <- "stim.EpiCellType"
table(Idents(PINvTumor.combined.Epi1))

#FeaturePlots
Idents(object = PINvTumor.combined.Epi1) <- "stim"
PINonlyEpi <- subset(PINvTumor.combined.Epi1, idents = c("PIN"))
TumoronlyEpi <- subset(PINvTumor.combined.Epi1, idents = c("Tumor"))

DefaultAssay(PINonlyEpi) <- "RNA"
tiff(file = "PINonlyEpi hAR Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi mGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi Krt5 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi Krt8 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi Ki67 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "PINonlyEpi Pcna Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Pcna"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

DefaultAssay(TumoronlyEpi) <- "RNA"
tiff(file = "TumoronlyEpi hAR Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi mGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi Krt5 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi Krt8 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi Ki67 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TumoronlyEpi Pcna Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Pcna"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Add Ki67 info
DefaultAssay(PINvTumor.combined.Epi1) <- "RNA"
PINvTumor.combined.Epi1Mki67Pos <- subset(x=PINvTumor.combined.Epi1, subset = Mki67 > 0)
PINvTumor.combined.Epi1Mki67Neg <- subset(x=PINvTumor.combined.Epi1, subset = Mki67 == 0)
Idents(object = PINvTumor.combined.Epi1Mki67Pos) <- "Mki67Pos"
Idents(object = PINvTumor.combined.Epi1Mki67Neg) <- "Mki67Neg"
PINvTumor.combined.Epi1Mki67Pos[["Mki67Exp"]] <- Idents(object = PINvTumor.combined.Epi1Mki67Pos)
PINvTumor.combined.Epi1Mki67Neg[["Mki67Exp"]] <- Idents(object = PINvTumor.combined.Epi1Mki67Neg)
PINvTumor.combined.Epi1Mki67 <- merge(x = PINvTumor.combined.Epi1Mki67Pos, y = PINvTumor.combined.Epi1Mki67Neg)
Idents(object = PINvTumor.combined.Epi1Mki67) <- "Mki67Exp"
PINvTumor.combined.Epi1$Mki67Exp <- Idents(object = PINvTumor.combined.Epi1Mki67)
Idents(object = PINvTumor.combined.Epi1) <- "Mki67Exp"
PINvTumor.combined.Epi1$Mki67Exp.EpiCellType <- paste(Idents(PINvTumor.combined.Epi1), PINvTumor.combined.Epi1$EpiCellType, sep = "_")
Idents(object = PINvTumor.combined.Epi1) <- "Mki67Exp.EpiCellType"
PINvTumor.combined.Epi1$stim.Mki67Exp.EpiCellType <- paste(Idents(PINvTumor.combined.Epi1), PINvTumor.combined.Epi1$stim, sep = "_")
Idents(object = PINvTumor.combined.Epi1) <- "stim.Mki67Exp.EpiCellType"
table(Idents(PINvTumor.combined.Epi1))

#Add Pcna info
DefaultAssay(PINvTumor.combined.Epi1) <- "RNA"
PINvTumor.combined.Epi1PcnaPos <- subset(x=PINvTumor.combined.Epi1, subset = Pcna > 0.5)
PINvTumor.combined.Epi1PcnaNeg <- subset(x=PINvTumor.combined.Epi1, subset = Pcna < 0.5)
Idents(object = PINvTumor.combined.Epi1PcnaPos) <- "PcnaPos"
Idents(object = PINvTumor.combined.Epi1PcnaNeg) <- "PcnaNeg"
PINvTumor.combined.Epi1PcnaPos[["PcnaExp"]] <- Idents(object = PINvTumor.combined.Epi1PcnaPos)
PINvTumor.combined.Epi1PcnaNeg[["PcnaExp"]] <- Idents(object = PINvTumor.combined.Epi1PcnaNeg)
PINvTumor.combined.Epi1Pcna <- merge(x = PINvTumor.combined.Epi1PcnaPos, y = PINvTumor.combined.Epi1PcnaNeg)
Idents(object = PINvTumor.combined.Epi1Pcna) <- "PcnaExp"
PINvTumor.combined.Epi1$PcnaExp <- Idents(object = PINvTumor.combined.Epi1Pcna)
Idents(object = PINvTumor.combined.Epi1) <- "PcnaExp"
PINvTumor.combined.Epi1$PcnaExp.EpiCellType <- paste(Idents(PINvTumor.combined.Epi1), PINvTumor.combined.Epi1$EpiCellType, sep = "_")
Idents(object = PINvTumor.combined.Epi1) <- "PcnaExp.EpiCellType"
PINvTumor.combined.Epi1$stim.PcnaExp.EpiCellType <- paste(Idents(PINvTumor.combined.Epi1), PINvTumor.combined.Epi1$stim, sep = "_")
Idents(object = PINvTumor.combined.Epi1) <- "stim.PcnaExp.EpiCellType"
table(Idents(PINvTumor.combined.Epi1))


#### ARQ+vARQ- BE in Integrated ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.3")

Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
tiff(file = "onlyBE Highlighted TSNE.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("skyblue", "blue", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

#Add ARQ info
DefaultAssay(PINvTumor.combined.Epi1) <- "RNA"
PINvTumor.combined.Epi1ARQPos <- subset(x=PINvTumor.combined.Epi1, subset = ARQ > 0)
PINvTumor.combined.Epi1ARQNeg <- subset(x=PINvTumor.combined.Epi1, subset = ARQ == 0)
Idents(object = PINvTumor.combined.Epi1ARQPos) <- "ARQPos"
Idents(object = PINvTumor.combined.Epi1ARQNeg) <- "ARQNeg"
PINvTumor.combined.Epi1ARQPos[["ARQExp"]] <- Idents(object = PINvTumor.combined.Epi1ARQPos)
PINvTumor.combined.Epi1ARQNeg[["ARQExp"]] <- Idents(object = PINvTumor.combined.Epi1ARQNeg)
PINvTumor.combined.Epi1ARQ <- merge(x = PINvTumor.combined.Epi1ARQPos, y = PINvTumor.combined.Epi1ARQNeg)
Idents(object = PINvTumor.combined.Epi1ARQ) <- "ARQExp"
PINvTumor.combined.Epi1$ARQExp <- Idents(object = PINvTumor.combined.Epi1ARQ)
Idents(object = PINvTumor.combined.Epi1) <- "ARQExp"
PINvTumor.combined.Epi1$ARQExp.EpiCellType <- paste(Idents(PINvTumor.combined.Epi1), PINvTumor.combined.Epi1$EpiCellType, sep = "_")
Idents(object = PINvTumor.combined.Epi1) <- "ARQExp.EpiCellType"


Idents(object = PINvTumor.combined.Epi1) <- "ARQExp"
PINvTumor.combined.Epi2 <- subset(PINvTumor.combined.Epi1, idents = c("ARQPos", "ARQNeg"))

Idents(object = PINvTumor.combined.Epi2) <- "ARQExp.EpiCellType"
tiff(file = "onlyBE ARQ+vARQ- UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi2, reduction = "umap", pt.size = 0.3, cols = c("#3399FF", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "#E06666", "grey", "grey", "grey", "grey", "grey", "#E06666", "#3399FF", "grey", "grey", "grey", "grey", "grey")) + NoLegend()
dev.off

Idents(object = PINvTumor.combined.Epi1) <- "CellType"
onlyBE2 <- subset(PINvTumor.combined.Epi1, idents = c("BE1", "BE2"))

Idents(object = onlyBE2) <- "ARQExp"
onlyBE2 <- subset(onlyBE2, idents = c("ARQPos", "ARQNeg"))
onlyBE2 <- RenameIdents(object = onlyBE2, 'ARQNeg' = "ARQNeg", 'ARQPos' = "ARQPos")
DimPlot(onlyBE2, reduction = "umap", pt.size = 0.3)

#DEGs
onlyBE2.0.1.Markers <- FindMarkers(onlyBE2, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(onlyBE2.0.1.Markers, "onlyBE2.0.1.Markers.csv")
onlyBE2.0.Markers <- FindMarkers(onlyBE2, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(onlyBE2.0.Markers, "onlyBE2.0.Markers.csv")

#p.adjust
DEG_onlyBE2 <- read.csv("onlyBE2.0.Markers.csv") 
DEG_onlyBE2_pvalue <- DEG_onlyBE2$p_val
DEG_onlyBE2_pvalue=as.numeric(DEG_onlyBE2_pvalue)
DEG_onlyBE2_BH = p.adjust(DEG_onlyBE2_pvalue, "BH")
write.csv(DEG_onlyBE2_BH, "DEG_onlyBE2_BH.csv")

#Heatmap
DefaultAssay(onlyBE2) <- "RNA"
onlyBE2 <- ScaleData(onlyBE2, features = rownames(onlyBE2))
onlyBE2.markers <- FindAllMarkers(onlyBE2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
onlyBE2Top50 <- onlyBE2.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "onlyBE2 Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(onlyBE2, features = c(onlyBE2Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs-1
Idents(object = PINvTumor.combined.Epi1) <- "ARQExp.EpiCellType"
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3)
        
onlyBE3 <- subset(PINvTumor.combined.Epi1, idents = c("ARQNeg_BE1", "ARQPos_BE2"))

onlyBE3.0.1.Markers <- FindMarkers(onlyBE3, ident.1 = "ARQPos_BE2", ident.2 = "ARQNeg_BE1", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(onlyBE3.0.1.Markers, "onlyBE3.0.1.Markers.csv")
onlyBE3.0.Markers <- FindMarkers(onlyBE3, ident.1 = "ARQPos_BE2", ident.2 = "ARQNeg_BE1", min.pct = 0, logfc.threshold = 0)
write.csv(onlyBE3.0.Markers, "onlyBE3.0.Markers.csv")

#p.adjust
DEG_onlyBE3 <- read.csv("onlyBE3.0.Markers.csv") 
DEG_onlyBE3_pvalue <- DEG_onlyBE3$p_val
DEG_onlyBE3_pvalue=as.numeric(DEG_onlyBE3_pvalue)
DEG_onlyBE3_BH = p.adjust(DEG_onlyBE3_pvalue, "BH")
write.csv(DEG_onlyBE3_BH, "DEG_onlyBE3_BH.csv")

#Heatmap
Idents(object = onlyBE3) <- "ARQExp.EpiCellType"
DefaultAssay(onlyBE3) <- "RNA"
onlyBE3 <- ScaleData(onlyBE3, features = rownames(onlyBE3))
onlyBE3.markers <- FindAllMarkers(onlyBE3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
onlyBE3Top50 <- onlyBE3.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "onlyBE3 Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(onlyBE3, features = c(onlyBE3Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Boxplot Generation

#Add ARQ info
DefaultAssay(onlyBE2) <- "RNA"
onlyBE2ARQNeg <- subset(x=onlyBE2, subset = ARQ == 0)
onlyBE2ARQPos <- subset(x=onlyBE2, subset = ARQ > 0)
Idents(object = onlyBE2ARQNeg) <- "ARQNeg"
Idents(object = onlyBE2ARQPos) <- "ARQPos"
onlyBE2ARQNeg[["ARQExp"]] <- Idents(object = onlyBE2ARQNeg)
onlyBE2ARQPos[["ARQExp"]] <- Idents(object = onlyBE2ARQPos)
onlyBE2ARQ <- merge(x = onlyBE2ARQNeg, y = onlyBE2ARQPos)
Idents(object = onlyBE2ARQ) <- "ARQExp"
onlyBE2$ARQExp <- Idents(object = onlyBE2ARQ)

boxdata = FetchData(onlyBE2, c("ARQExp", "ARQ", "Igf1r", "Fos", "Jak2", "Mapk13"))
tail(boxdata,6)

tiff(file = "onlyBE2 hARtg Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=ARQ, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "onlyBE2 Igf1r Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Igf1r, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "onlyBE2 Fos Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Fos, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
dev.off()
tiff(file = "onlyBE2 Jak2 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Jak2, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
dev.off()
tiff(file = "onlyBE2 Mapk13 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Mapk13, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2.25))
dev.off()

#Boxplot Generation

#Add ARQ info
DefaultAssay(onlyBE3) <- "RNA"
onlyBE3ARQNeg <- subset(x=onlyBE3, subset = ARQ == 0)
onlyBE3ARQPos <- subset(x=onlyBE3, subset = ARQ > 0)
Idents(object = onlyBE3ARQNeg) <- "ARQNeg"
Idents(object = onlyBE3ARQPos) <- "ARQPos"
onlyBE3ARQNeg[["ARQExp"]] <- Idents(object = onlyBE3ARQNeg)
onlyBE3ARQPos[["ARQExp"]] <- Idents(object = onlyBE3ARQPos)
onlyBE3ARQ <- merge(x = onlyBE3ARQNeg, y = onlyBE3ARQPos)
Idents(object = onlyBE3ARQ) <- "ARQExp"
onlyBE3$ARQExp <- Idents(object = onlyBE3ARQ)

boxdata = FetchData(onlyBE3, c("ARQExp", "ARQ", "Igf1r", "Fos", "Jak2", "Mapk13", "Jun"))
tail(boxdata,6)

#w/o dots
tiff(file = "onlyBE3 hARtg Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=ARQ, fill = ARQExp)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "onlyBE3 Igf1r Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Igf1r, fill = ARQExp)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "onlyBE3 Fos Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Fos, fill = ARQExp)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
dev.off()
tiff(file = "onlyBE3 Jak2 Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Jak2, fill = ARQExp)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
dev.off() 
tiff(file = "onlyBE3 Mapk13 Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Mapk13, fill = ARQExp)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2.25))
dev.off() 

#w/ dots
tiff(file = "onlyBE3 hARtg Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=ARQ, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "onlyBE3 Igf1r Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Igf1r, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "onlyBE3 Fos Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Fos, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
dev.off()
tiff(file = "onlyBE3 Jak2 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Jak2, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
dev.off()
tiff(file = "onlyBE3 Mapk13 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=ARQExp, y=Mapk13, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2.25))
dev.off()

####Gene-Gene Spearman correlation####
install.packages("GGally")
library(GGally)

Idents(object = PINvTumor.combined.Epi1) <- "CellType"
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3)
onlyBE <- subset(PINvTumor.combined.Epi1, ident = c("BE1", "BE2"))

#Seurat object.
GOI <- c('ARQ','Igf1r','Fos','Jak2', 'Mapk13')  
GOI_index <- is.element(rownames(onlyBE),GOI)
Cell_index <- is.element(Idents(onlyBE),c('BE1','BE2'))

expr_GOI <- onlyBE@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- onlyBE@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI))

tiff(file = "spearman correlation BE.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "spearman correlation BE NoLabel.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))
dev.off()

head(ggcorr(t(expr_GOI))$data, 10)
head(ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))$data, 10)

Idents(object = PINvTumor.combined.Epi1) <- "ARQExp.EpiCellType"
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3)
onlyBE3 <- subset(PINvTumor.combined.Epi1, idents = c("ARQNeg_BE1", "ARQPos_BE2"))
Idents(object = onlyBE3) <- "EpiCellType"
DimPlot(onlyBE3, reduction = "umap", pt.size = 0.3)


#Seurat object.
GOI <- c('ARQ','Igf1r','Fos','Jak2', 'Mapk13')  
GOI_index <- is.element(rownames(onlyBE3),GOI)
Cell_index <- is.element(Idents(onlyBE3),c('BE1','BE2'))

expr_GOI <- onlyBE3@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- onlyBE3@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI))

tiff(file = "spearman correlation BE.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "spearman correlation BE NoLabel.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))
dev.off()

head(ggcorr(t(expr_GOI))$data, 10)
head(ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))$data, 10)

#Seurat object.
GOI1 <- c('ARQ','Igf1r', 'Jak2', 'Mapk13')  
GOI_index1 <- is.element(rownames(onlyBE3),GOI1)
Cell_index1 <- is.element(Idents(onlyBE3),c('BE1','BE2'))

expr_GOI1 <- onlyBE3@assays$RNA@data[GOI_index1,Cell_index1] 
expr_GOI1 <- onlyBE3@assays$RNA@counts[GOI_index1,Cell_index1]
ggcorr(t(expr_GOI1))

tiff(file = "spearman correlation BE except Fos.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI1), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "spearman correlation BE NoLabel except Fos.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI1), method = c("pairwise", "spearman"))
dev.off()

head(ggcorr(t(expr_GOI1), method = c("pairwise", "spearman"))$data, 10)

#Seurat object.
GOI2 <- c('ARQ','Igf1r', 'Jak2', 'Mapk13', 'Akt1', 'Foxo1', 'Hras', 'Socs2', 'Socs3', 'Stat3', 'Elk1')  
GOI_index2 <- is.element(rownames(onlyBE3),GOI2)
Cell_index2 <- is.element(Idents(onlyBE3),c('BE1','BE2'))

expr_GOI2 <- onlyBE3@assays$RNA@data[GOI_index2,Cell_index2] 
expr_GOI2 <- onlyBE3@assays$RNA@counts[GOI_index2,Cell_index2]
ggcorr(t(expr_GOI2))

tiff(file = "spearman correlation BE except Fos.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI2), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "spearman correlation BE NoLabel except Fos.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI2), method = c("pairwise", "spearman"))
dev.off()

head(ggcorr(t(expr_GOI2), method = c("pairwise", "spearman"))$data, 50)


####Merge with WT and PIN and Tumor####

WTunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191214_32848_33580_test/counts/33580_WT/outs/filtered_feature_bc_matrix")
WTunfiltered <- CreateSeuratObject(counts = WTunfiltered.data,  min.cells = 3, min.features = 500, project = "WTunfiltered")
WTunfiltered <- NormalizeData(WTunfiltered)

#Initial processing & filtering
WTunfiltered[["percent.mt"]] <- PercentageFeatureSet(WTunfiltered, pattern = "^mt-")

tiff(file = "WT Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WTunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "darkgreen", xlab = "nFeature_RNA", main = "WT Pre-filteration")
dev.off()
tiff(file = "WT Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$percent.mt, breaks = 100, col = "darkgreen", xlab = "percent.mt", main = "WT Pre-filteration")
dev.off()

WT <- subset(WTunfiltered, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
tiff(file = "WT Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$nFeature_RNA, breaks = 100, col = "darkgreen", xlab = "nFeature_RNA", main = "WT Post-filteration")
dev.off()
tiff(file = "WT Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$percent.mt, breaks = 100, col = "darkgreen", xlab = "percent.mt", main = "WT Post-filteration")
dev.off()

WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 5000)
tiff(file = "WT Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT)
dev.off()

table(Idents(WTunfiltered))
table(Idents(WT))

#Run the standard workflow for visualization and clustering Ctrl2
WT <- ScaleData(WT, verbose = FALSE)
WT <- RunPCA(WT, npcs = 30, verbose = FALSE)
ElbowPlot(WT, ndims = 50)
# t-SNE and Clustering
WT <- FindNeighbors(WT, reduction = "pca", dims = 1:20)
WT <- FindClusters(WT, resolution = 0.5)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:20)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:20)
DimPlot(WT, reduction = "umap", pt.size = 0.3, label = TRUE) 

WT.Epi <- subset(WT, idents = c("0", "4", "2", "7", "11", "6", "3", "9"))
DimPlot(WT.Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 

DefaultAssay(WT.Epi) <- "RNA"
FeaturePlot(WT.Epi, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
FeaturePlot(WT.Epi, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")

#Merge
Idents(object = PINvTumor.combined.Epi1) <- "stim"
PIN.Epi <- subset(PINvTumor.combined.Epi1, idents = c("PIN"))
Tumor.Epi <- subset(PINvTumor.combined.Epi1, idents = c("Tumor"))

Idents(object = PIN.Epi) <- "seurat_clusters"
Idents(object = Tumor.Epi) <- "seurat_clusters"
Idents(object = WT.Epi) <- "seurat_clusters"

PIN.Epi[["orig.clusters"]] <- Idents(object = PIN.Epi)
Tumor.Epi[["orig.clusters"]] <- Idents(object = Tumor.Epi)
WT.Epi[["orig.clusters"]] <- Idents(object = WT.Epi)

Idents(object = PIN.Epi) <- "seurat_clusters"
Idents(object = Tumor.Epi) <- "seurat_clusters"
Idents(object = WT.Epi) <- "seurat_clusters"

PIN.Epi$stim <- "PIN"
Tumor.Epi$stim <- "Tumor"
WT.Epi$stim <- "WT"

ARQ9vWT.Epi.anchors <- FindIntegrationAnchors(object.list = list(PIN.Epi, Tumor.Epi, WT.Epi), dims = 1:20)
ARQ9vWT.Epi.combined <- IntegrateData(anchorset = ARQ9vWT.Epi.anchors, dims = 1:20)

ARQ9vWT.Epi.combined1 <- ARQ9vWT.Epi.combined

#without re-clustering
DefaultAssay(ARQ9vWT.Epi.combined) <- "integrated"
ARQ9vWT.Epi.combined <- ScaleData(ARQ9vWT.Epi.combined, verbose = FALSE)
ARQ9vWT.Epi.combined <- RunPCA(ARQ9vWT.Epi.combined, npcs = 30, verbose = FALSE)
ARQ9vWT.Epi.combined <- RunUMAP(ARQ9vWT.Epi.combined, reduction = "pca", dims = 1:20)
DimPlot(ARQ9vWT.Epi.combined, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = ARQ9vWT.Epi.combined) <- "EpiCellType"
ARQ9vWT.Epi.combined <- RenameIdents(object = ARQ9vWT.Epi.combined, 'BE1' = "BE1", 'BE2' = "BE2",
                                     'LE1' = "LE1", 'LE2' = "LE2", 
                                      'LE3' = "LE3", 'LE4' = "LE4", 'LE5' = "LE5", 'LE6' = "LE6", 'LE7' = "LE7",
                                      'ProE' = "ProE", 'OE' = "OE")
ARQ9vWT.Epi.combined[["WTEpiCellType"]] <- Idents(object = ARQ9vWT.Epi.combined)

tiff(file = "ARQ9vWT.Epi.combined UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50")) 
dev.off()
tiff(file = "ARQ9vWT.Epi.combined split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50")) 
dev.off()

#with re-clustering
DefaultAssay(ARQ9vWT.Epi.combined1) <- "integrated"
ARQ9vWT.Epi.combined1 <- ScaleData(ARQ9vWT.Epi.combined1, verbose = FALSE)
ARQ9vWT.Epi.combined1 <- RunPCA(ARQ9vWT.Epi.combined1, npcs = 30, verbose = FALSE)
ElbowPlot(ARQ9vWT.Epi.combined1, ndims = 50)
ARQ9vWT.Epi.combined1 <- FindNeighbors(ARQ9vWT.Epi.combined1, reduction = "pca", dims = 1:20)
ARQ9vWT.Epi.combined1 <- FindClusters(ARQ9vWT.Epi.combined1, resolution = 0.7)
ARQ9vWT.Epi.combined1 <- RunUMAP(ARQ9vWT.Epi.combined1, reduction = "pca", dims = 1:20)
DimPlot(ARQ9vWT.Epi.combined1, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = ARQ9vWT.Epi.combined1) <- "seurat_clusters"
ARQ9vWT.Epi.combined1 <- RenameIdents(object = ARQ9vWT.Epi.combined1, '0' = "BE1", '11' = "BE2",
                                     '1' = "LE1", '8' = "LE1", '3' = "LE1", 
                                     '2' = "LE2", '9' = "LE3", '13' = "LE4", '4' = "LE5", '5' = "LE5", '6' = "LE6",
                                     '7' = "LE7", '10' = "LE7", '12' = "UrLE", '14' = "OE", '15' = "OE")
ARQ9vWT.Epi.combined1[["WTEpiCellType"]] <- Idents(object = ARQ9vWT.Epi.combined1)

Idents(object = ARQ9vWT.Epi.combined1) <- "WTEpiCellType"
tiff(file = "ARQ9vWT.Epi.combined1 UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50")) 
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1 split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50")) 
dev.off()

Idents(object = ARQ9vWT.Epi.combined1) <- "stim"
ARQ9vWT.Epi.combined1_PIN <- subset(ARQ9vWT.Epi.combined1, idents = c("PIN"))
ARQ9vWT.Epi.combined1_Tumor <- subset(ARQ9vWT.Epi.combined1, idents = c("Tumor"))
ARQ9vWT.Epi.combined1_WT <- subset(ARQ9vWT.Epi.combined1, idents = c("WT"))

DefaultAssay(ARQ9vWT.Epi.combined1_PIN) <- "RNA"
tiff(file = "ARQ9vWT.Epi.combined1_PIN hARtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_PIN, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_PIN Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_PIN, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_PIN Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_PIN, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_PIN Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_PIN, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()

DefaultAssay(ARQ9vWT.Epi.combined1_Tumor) <- "RNA"
tiff(file = "ARQ9vWT.Epi.combined1_Tumor hARtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_Tumor, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_Tumor Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_Tumor, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_Tumor Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_Tumor, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_Tumor Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_Tumor, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()

DefaultAssay(ARQ9vWT.Epi.combined1_WT) <- "RNA"
tiff(file = "ARQ9vWT.Epi.combined1_WT hARtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_WT, reduction = "umap", features = c("ARQ"), cols = c("light grey", "light grey"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_WT Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_WT, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_WT Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_WT, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined1_WT Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined1_WT, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = "q90")
dev.off()


####Merge with WT and ARQ########
PINvTumor.combined.Epi1[["orig.clusters"]] <- Idents(object = PINvTumor.combined.Epi1)
WT.Epi[["orig.clusters"]] <- Idents(object = WT.Epi)

Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
Idents(object = WT.Epi) <- "seurat_clusters"

PINvTumor.combined.Epi1$stim.EpiCellType <- "ARQ"
WT.Epi$stim.EpiCellType <- "WT"
WT.Epi$stim <- "WT"

ARQ9vWT.Epi.anchors1 <- FindIntegrationAnchors(object.list = list(PINvTumor.combined.Epi1, WT.Epi), dims = 1:20)
ARQ9vWT.Epi.combined3 <- IntegrateData(anchorset = ARQ9vWT.Epi.anchors1, dims = 1:20)

mouse_cell_cycle_genes <- readRDS("mouse_cell_cycle_genes.Rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(ARQ9vWT.Epi.combined3) <- "RNA"
all.genes <- rownames(ARQ9vWT.Epi.combined3)
ARQ9vWT.Epi.combined3 <- ScaleData(ARQ9vWT.Epi.combined3, features = all.genes)
ARQ9vWT.Epi.combined3 <- CellCycleScoring(ARQ9vWT.Epi.combined3, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Take Cell cycle out 
DefaultAssay(ARQ9vWT.Epi.combined3) <- "integrated"
ARQ9vWT.Epi.combined3 <- ScaleData(ARQ9vWT.Epi.combined3, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ARQ9vWT.Epi.combined3))
ARQ9vWT.Epi.combined3 <- RunPCA(ARQ9vWT.Epi.combined3, features = VariableFeatures(ARQ9vWT.Epi.combined3))
ElbowPlot(ARQ9vWT.Epi.combined3, ndims = 30)
ARQ9vWT.Epi.combined3 <- RunUMAP(ARQ9vWT.Epi.combined3, reduction = "pca", dims = 1:20)
DimPlot(ARQ9vWT.Epi.combined3, reduction = "umap", pt.size = 0.3, label = TRUE)

ARQ9vWT.Epi.combined4 <- ARQ9vWT.Epi.combined3

#without re-clustering
Idents(object = ARQ9vWT.Epi.combined3) <- "EpiCellType"
ARQ9vWT.Epi.combined3 <- RenameIdents(object = ARQ9vWT.Epi.combined3, 'BE1' = "BE1", 'BE2' = "BE2",
                                     'LE1' = "LE1", 'LE2' = "LE2", 
                                     'LE3' = "LE3", 'LE4' = "LE4", 'LE5' = "LE5", 'LE6' = "LE6", 'LE7' = "LE7", 'UrLE' = "UrLE",
                                     'OE' = "OE")
ARQ9vWT.Epi.combined3[["WTEpiCellType"]] <- Idents(object = ARQ9vWT.Epi.combined3)
DimPlot(ARQ9vWT.Epi.combined3, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "dark grey", "light grey")) 

tiff(file = "ARQ9vWT.Epi.combined3 UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined3, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "dark grey", "light grey")) 
dev.off()
tiff(file = "ARQ9vWT.Epi.combined3 split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined3, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "dark grey", "light grey")) 
dev.off()

#with re-clustering
DefaultAssay(ARQ9vWT.Epi.combined4) <- "RNA"
all.genes <- rownames(ARQ9vWT.Epi.combined4)
ARQ9vWT.Epi.combined4 <- ScaleData(ARQ9vWT.Epi.combined4, features = all.genes)
ARQ9vWT.Epi.combined4 <- CellCycleScoring(ARQ9vWT.Epi.combined4, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#Take Cell cycle out 
DefaultAssay(ARQ9vWT.Epi.combined4) <- "integrated"
ARQ9vWT.Epi.combined4 <- ScaleData(ARQ9vWT.Epi.combined4, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ARQ9vWT.Epi.combined4))
ARQ9vWT.Epi.combined4 <- RunPCA(ARQ9vWT.Epi.combined4, features = VariableFeatures(ARQ9vWT.Epi.combined4))
ElbowPlot(ARQ9vWT.Epi.combined4, ndims = 30)

ARQ9vWT.Epi.combined4 <- FindNeighbors(ARQ9vWT.Epi.combined4, reduction = "pca", dims = 1:20)
ARQ9vWT.Epi.combined4 <- FindClusters(ARQ9vWT.Epi.combined4, resolution = 0.8)
ARQ9vWT.Epi.combined4 <- RunUMAP(ARQ9vWT.Epi.combined4, reduction = "pca", dims = 1:20)
DimPlot(ARQ9vWT.Epi.combined4, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = ARQ9vWT.Epi.combined4) <- "seurat_clusters"
ARQ9vWT.Epi.combined4 <- subset(ARQ9vWT.Epi.combined4, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18"))
ARQ9vWT.Epi.combined4 <- RenameIdents(object = ARQ9vWT.Epi.combined4, '1' = "BE1", '5' = "BE1", '9' = "BE2",
                                      '0' = "LE1", '4' = "LE1", 
                                      '8' = "LE2", '12' = "LE3", '2' = "LE3", '10' = "LE3", '11' = "LE4", '3' = "LE5", '13' = "LE6", '17' = "LE6", '15' = "LE6", '6' = "LE6",
                                      '7' = "LE7", '14' = "LE7", '16' = "UrLE", '18' = "OE")
ARQ9vWT.Epi.combined4[["WTEpiCellType"]] <- Idents(object = ARQ9vWT.Epi.combined4)

Idents(object = ARQ9vWT.Epi.combined4) <- "WTEpiCellType"
tiff(file = "ARQ9vWT.Epi.combined4 UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined4, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "dark grey", "light grey")) 
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4 split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vWT.Epi.combined4, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "dark grey", "light grey")) 
dev.off()

Idents(object = ARQ9vWT.Epi.combined4) <- "stim"
ARQ9vWT.Epi.combined4_PIN <- subset(ARQ9vWT.Epi.combined4, idents = c("PIN"))
ARQ9vWT.Epi.combined4_Tumor <- subset(ARQ9vWT.Epi.combined4, idents = c("Tumor"))
ARQ9vWT.Epi.combined4_WT <- subset(ARQ9vWT.Epi.combined4, idents = c("WT"))

DefaultAssay(ARQ9vWT.Epi.combined4_PIN) <- "RNA"
tiff(file = "ARQ9vWT.Epi.combined4_PIN hARtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_PIN, reduction = "umap", features = c("ARQ"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_PIN mGFP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_PIN, reduction = "umap", features = c("EGFP"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_PIN Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_PIN, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_PIN Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_PIN, reduction = "umap", features = c("Krt19"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_PIN Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_PIN, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

DefaultAssay(ARQ9vWT.Epi.combined4_Tumor) <- "RNA"
tiff(file = "ARQ9vWT.Epi.combined4_Tumor hARtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_Tumor, reduction = "umap", features = c("ARQ"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_Tumor mGFP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_Tumor, reduction = "umap", features = c("EGFP"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_Tumor Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_Tumor, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_Tumor Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_Tumor, reduction = "umap", features = c("Krt19"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_Tumor Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_Tumor, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

DefaultAssay(ARQ9vWT.Epi.combined4_WT) <- "RNA"
tiff(file = "ARQ9vWT.Epi.combined4_WT hARtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_WT, reduction = "umap", features = c("ARQ"), split.by = "stim", cols = c("light grey", "light grey"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_WT mGFP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_WT, reduction = "umap", features = c("EGFP"), split.by = "stim", cols = c("light grey", "light grey"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_WT Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_WT, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_WT Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_WT, reduction = "umap", features = c("Krt19"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()
tiff(file = "ARQ9vWT.Epi.combined4_WT Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vWT.Epi.combined4_WT, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

#DEGs
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.3")

Idents(object = ARQ9vWT.Epi.combined4) <- "WTEpiCellType"
ARQ9vWT.Epi.combined4$WTEpiCellType.stim.EpiCellType <- paste(Idents(ARQ9vWT.Epi.combined4), ARQ9vWT.Epi.combined4$stim.EpiCellType, sep = "_")
Idents(object = ARQ9vWT.Epi.combined4) <- "WTEpiCellType.stim.EpiCellType"
DimPlot(ARQ9vWT.Epi.combined4, reduction = "umap", pt.size = 0.3)

BE2_ARQvBE_WT <- subset(ARQ9vWT.Epi.combined4, idents = c("BE2_ARQ", "BE1_WT", "BE2_WT"))
BE_ARQvBE_WT <- subset(ARQ9vWT.Epi.combined4, idents = c("BE2_ARQ", "BE1_ARQ", "BE1_WT", "BE2_WT"))

BE2_ARQvBE_WT  <- RenameIdents(object = BE2_ARQvBE_WT, 'BE2_ARQ' = "BE2_ARQ", 'BE1_WT' = "BE_WT", 'BE2_WT' = "BE_WT")
BE2_ARQvBE_WT[["BE2_ARQvBE_WT"]] <- Idents(object = BE2_ARQvBE_WT)

BE_ARQvBE_WT  <- RenameIdents(object = BE_ARQvBE_WT, 'BE2_ARQ' = "BE_ARQ", 'BE1_ARQ' = "BE_ARQ", 'BE1_WT' = "BE_WT", 'BE2_WT' = "BE_WT")
BE_ARQvBE_WT[["BE_ARQvBE_WT"]] <- Idents(object = BE_ARQvBE_WT)

#BE2_ARQvBE_WT
Idents(object = BE2_ARQvBE_WT) <- "BE2_ARQvBE_WT"

BE2_ARQvBE_WT.0.1.Markers <- FindMarkers(BE2_ARQvBE_WT, ident.1 = "BE2_ARQ", ident.2 = "BE_WT", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(BE2_ARQvBE_WT.0.1.Markers, "BE2_ARQvBE_WT.0.1.Markers.csv")
BE2_ARQvBE_WT.0.Markers <- FindMarkers(BE2_ARQvBE_WT, ident.1 = "BE2_ARQ", ident.2 = "BE_WT", min.pct = 0, logfc.threshold = 0)
write.csv(BE2_ARQvBE_WT.0.Markers, "BE2_ARQvBE_WT.0.Markers.csv")

#p.adjust
DEG_BE2_ARQvBE_WT <- read.csv("BE2_ARQvBE_WT.0.Markers.csv") 
DEG_BE2_ARQvBE_WT_pvalue <- DEG_BE2_ARQvBE_WT$p_val
DEG_BE2_ARQvBE_WT_pvalue=as.numeric(DEG_BE2_ARQvBE_WT_pvalue)
DEG_BE2_ARQvBE_WT_BH = p.adjust(DEG_BE2_ARQvBE_WT_pvalue, "BH")
write.csv(DEG_BE2_ARQvBE_WT_BH, "DEG_BE2_ARQvBE_WT_BH.csv")

#BE_ARQvBE_WT
Idents(object = BE_ARQvBE_WT) <- "BE_ARQvBE_WT"

BE_ARQvBE_WT.0.1.Markers <- FindMarkers(BE_ARQvBE_WT, ident.1 = "BE_ARQ", ident.2 = "BE_WT", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(BE_ARQvBE_WT.0.1.Markers, "BE_ARQvBE_WT.0.1.Markers.csv")
BE_ARQvBE_WT.0.Markers <- FindMarkers(BE_ARQvBE_WT, ident.1 = "BE_ARQ", ident.2 = "BE_WT", min.pct = 0, logfc.threshold = 0)
write.csv(BE_ARQvBE_WT.0.Markers, "BE_ARQvBE_WT.0.Markers.csv")

#p.adjust
DEG_BE_ARQvBE_WT <- read.csv("BE_ARQvBE_WT.0.Markers.csv") 
DEG_BE_ARQvBE_WT_pvalue <- DEG_BE_ARQvBE_WT$p_val
DEG_BE_ARQvBE_WT_pvalue=as.numeric(DEG_BE_ARQvBE_WT_pvalue)
DEG_BE_ARQvBE_WT_BH = p.adjust(DEG_BE_ARQvBE_WT_pvalue, "BH")
write.csv(DEG_BE_ARQvBE_WT_BH, "DEG_BE_ARQvBE_WT_BH.csv")

####LE1-2vLE5-7 in PINvTumor.combined.Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.4")

Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
tiff(file = "LE1-2 LE5-7 Highlighted TSNE.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "#E06666", "#E06666", "grey", "grey", "#3399FF", "#3399FF", "#3399FF", "grey", "grey"))
dev.off()

TumorvNormal.Epi1 <- subset(PINvTumor.combined.Epi1, idents = c("LE1", "LE2", "LE5", "LE6", "LE7"))
TumorvNormal.Epi1  <- RenameIdents(object = TumorvNormal.Epi1, 'LE5' = "Normal", 'LE6' = "Normal", 'LE7' = "Normal", 'LE1' = "Tumor", 'LE2' = "Tumor")
TumorvNormal.Epi1[["NormalTumor"]] <- Idents(object = TumorvNormal.Epi1)

#DEGs
TumorvNormal.Epi1.0.01.Markers <- FindMarkers(TumorvNormal.Epi1, ident.1 = "Tumor", ident.2 = "Normal", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(TumorvNormal.Epi1.0.01.Markers, "TumorvNormal.Epi1.0.01.Markers.csv")

TumorvNormal.Epi1.0.Markers <- FindMarkers(TumorvNormal.Epi1, ident.1 = "Tumor", ident.2 = "Normal", min.pct = 0., logfc.threshold = 0.)
write.csv(TumorvNormal.Epi1.0.Markers, "TumorvNormal.Epi1.0.Markers.csv")

#p.adjust
DEG_Tumor <- read.csv("TumorvNormal.Epi1.0.Markers.csv") 
DEG_Tumor_pvalue <- DEG_Tumor$p_val
DEG_Tumor_pvalue=as.numeric(DEG_Tumor_pvalue)
DEG_Tumor_BH = p.adjust(DEG_Tumor_pvalue, "BH")
write.csv(DEG_Tumor_BH, "DEG_Tumor_BH.csv")

#Heatmap
DefaultAssay(TumorvNormal.Epi1) <- "RNA"
TumorvNormal.Epi1 <- ScaleData(TumorvNormal.Epi1, features = rownames(TumorvNormal.Epi1))
TumorvNormal.Epi1.markers <- FindAllMarkers(TumorvNormal.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TumorvNormal.Epi1Top50 <- TumorvNormal.Epi1.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "LE1-2vLE5-7 Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TumorvNormal.Epi1, features = c(TumorvNormal.Epi1Top50$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Boxplots
boxdata = FetchData(TumorvNormal.Epi1, c("NormalTumor", "ARQ", "Tcf4", "Myc", "Ccnd1", "Axin2", "Lgr5"))
tail(boxdata,6)

#w/o dots
tiff(file = "TumorvNormal.Epi1 hARtg Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=ARQ, fill = NormalTumor)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Tcf4 Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Tcf4, fill = NormalTumor)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Myc Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Myc, fill = NormalTumor)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) #coord_cartesian(ylim=c(0, 2))
dev.off()
tiff(file = "TumorvNormal.Epi1 Ccnd1 Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Ccnd1, fill = NormalTumor)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) #coord_cartesian(ylim=c(0, 2))
dev.off() 
tiff(file = "TumorvNormal.Epi1 Axin2 Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Axin2, fill = NormalTumor)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) #coord_cartesian(ylim=c(0, 2.25))
dev.off() 
tiff(file = "TumorvNormal.Epi1 Lgr5 Boxplot with Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Lgr5, fill = NormalTumor)) + geom_boxplot(outlier.color = NA) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) #coord_cartesian(ylim=c(0, 2.25))
dev.off() 

#w/ dots
tiff(file = "TumorvNormal.Epi1 hARtg Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=ARQ, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Tcf4 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Tcf4, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Myc Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Myc, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Ccnd1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Ccnd1, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Axin2 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Axin2, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "TumorvNormal.Epi1 Lgr5 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=NormalTumor, y=Lgr5, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()

#expression Plots
DefaultAssay(PINvTumor.combined.Epi1) <- "RNA"
tiff(file = "PINvTumor.combined.Epi1 Tcf4.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Tcf4"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
dev.off()
tiff(file = "PINvTumor.combined.Epi1 Ccnd1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Ccnd1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
dev.off()
tiff(file = "PINvTumor.combined.Epi1 Axin2.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Axin2"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
dev.off()
tiff(file = "PINvTumor.combined.Epi1 Lgr5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Lgr5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
dev.off()

#Seurat object.
DefaultAssay(TumorvNormal.Epi1) <- "RNA"
Idents(object = TumorvNormal.Epi1) <- "NormalTumor"
DimPlot(TumorvNormal.Epi1, reduction = "umap", pt.size = 0.3)
GOI3 <- c('ARQ','Tcf4', 'Myc', 'Ccnd1', 'Axin2', 'Lgr5')  
GOI_index3 <- is.element(rownames(TumorvNormal.Epi1),GOI3)
Cell_index3 <- is.element(Idents(TumorvNormal.Epi1),c('Normal','Tumor'))

expr_GOI3 <- TumorvNormal.Epi1@assays$RNA@data[GOI_index3,Cell_index3] 
expr_GOI3 <- TumorvNormal.Epi1@assays$RNA@counts[GOI_index3,Cell_index3]
ggcorr(t(expr_GOI3))

tiff(file = "spearman correlation TumorvNormal.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI3), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "spearman correlation TumorvNormal.tiff NoLabel.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI3), method = c("pairwise", "spearman"))
dev.off()

head(ggcorr(t(expr_GOI3), method = c("pairwise", "spearman"))$data, 15)

#Except Myc
GOI4 <- c('ARQ','Tcf4', 'Ccnd1', 'Axin2', 'Lgr5')  
GOI_index4 <- is.element(rownames(TumorvNormal.Epi1),GOI4)
Cell_index4 <- is.element(Idents(TumorvNormal.Epi1),c('Normal','Tumor'))

expr_GOI4 <- TumorvNormal.Epi1@assays$RNA@data[GOI_index4,Cell_index4] 
expr_GOI4 <- TumorvNormal.Epi1@assays$RNA@counts[GOI_index4,Cell_index4]
ggcorr(t(expr_GOI4))

tiff(file = "spearman correlation TumorvNormal except Myc.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI4), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "spearman correlation TumorvNormal NoLabel except Myc.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI4), method = c("pairwise", "spearman"))
dev.off()

head(ggcorr(t(expr_GOI4), method = c("pairwise", "spearman"))$data, 15)

####Pseudotime using Monocle2 with hARtg+ vs hARtg-####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.6")

#Rename
Idents(object = PINvTumor.combined.Epi1) <- "seurat_clusters"
PINvTumor.combined.Epi1 <- RenameIdents(object = PINvTumor.combined.Epi1, '4' = "BE1", '6' = "BE2", '0' = "LE1", '2' = "LE2", '1' = "LE3", '9' = "LE4", '3' = "LE5", '5' = "LE6", '8' = "LE7", '7' = "UrLE", '10' = "OE")
PINvTumor.combined.Epi1[["EpiCellType"]] <- Idents(object = PINvTumor.combined.Epi1)


Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
PINvTumor.combined.Epi2 <- subset(PINvTumor.combined.Epi1, idents = c("BE1", "BE2", "LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7"))

Idents(object = PINvTumor.combined.Epi2) <- "seurat_clusters"
DimPlot(PINvTumor.combined.Epi2, reduction = "umap")

Idents(object = PINvTumor.combined.Epi2) <- "ARQExp"
ARQPosEpi <- subset(PINvTumor.combined.Epi2, idents = c("ARQPos"))
ARQNegEpi <- subset(PINvTumor.combined.Epi2, idents = c("ARQNeg"))

#ARQPosEpi
Idents(object = ARQPosEpi) <- "seurat_clusters"
DimPlot(ARQPosEpi, reduction = "tsne")


DefaultAssay(ARQPosEpi) <- "RNA"
EpiPseudo <- as.CellDataSet(ARQPosEpi)
EpiPseudo <- detectGenes(EpiPseudo, min_expr = 0.1)
print(head(fData(EpiPseudo)))

expressed_genes <- row.names(subset(fData(EpiPseudo),
                                    num_cells_expressed >= 10))

pData(EpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo))
EpiPseudo <- EpiPseudo[,pData(EpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo), geom =
        "density")

EpiPseudo <- estimateSizeFactors(EpiPseudo)
EpiPseudo <- estimateDispersions(EpiPseudo)

disp_table <- dispersionTable(EpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo <- setOrderingFilter(EpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo)

#EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo, return_all = F) # norm_method='log'

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
EpiPseudo <- clusterCells(EpiPseudo, num_clusters = 2)

plot_cell_clusters(EpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(EpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.1))
EpiPseudo <- setOrderingFilter(EpiPseudo, ordering_genes)
plot_ordering_genes(EpiPseudo)

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2,
                             method = 'DDRTree')

EpiPseudo <- orderCells(EpiPseudo)

GM_state <- function(EpiPseudo){
  if (length(unique(pData(EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo)$State, pData(EpiPseudo)$seurat_clusters)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo <- orderCells(EpiPseudo, root_state = GM_state(EpiPseudo))

#Visualization
tiff(file = "EpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "EpiPseudo Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo EpiCellType.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) 
dev.off()
tiff(file = "EpiPseudo stim EpiCellType stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo CellMarkers Wnt targets purple.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt5", "Krt14", "Krt8", "Krt18", "Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo CellMarkers Wnt targets.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt5", "Krt14", "Krt8", "Krt18", "Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "EpiPseudo Wnt targets.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo Wnt targets redcolor.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "EpiPseudo IGF1R Wnt stem cell targets.tiff", width = 16, height = 21, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Igf1r", "Mapk13", "Jak2", "Tcf4", "Ccnd1", "Axin2", "Lgr5", "Trp63", "Ly6a", "Psca", "Tacstd2", "Itga6", "Osr1", "Pbsn"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "EpiPseudo celltype IGF1R Wnt stem cell targets.tiff", width = 16, height = 21, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Krt5", "Krt8", "Igf1r", "Mapk13", "Jak2", "Tcf4", "Ccnd1", "Axin2", "Lgr5", "Trp63", "Ly6a", "Psca", "Tacstd2", "Itga6", "Osr1", "Pbsn"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

#ARQNegEpi
Idents(object = ARQNegEpi) <- "seurat_clusters"
DefaultAssay(ARQNegEpi) <- "RNA"
EpiPseudo1 <- as.CellDataSet(ARQNegEpi)
EpiPseudo1 <- detectGenes(EpiPseudo1, min_expr = 0.1)
print(head(fData(EpiPseudo1)))

expressed_genes <- row.names(subset(fData(EpiPseudo1),
                                    num_cells_expressed >= 10))

pData(EpiPseudo1)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo1))
EpiPseudo1 <- EpiPseudo1[,pData(EpiPseudo1)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo1), geom =
        "density")

EpiPseudo1 <- estimateSizeFactors(EpiPseudo1)
EpiPseudo1 <- estimateDispersions(EpiPseudo1)

disp_table <- dispersionTable(EpiPseudo1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo1 <- setOrderingFilter(EpiPseudo1, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo1)

#EpiPseudo1@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo1, return_all = F) # norm_method='log'

EpiPseudo1 <- reduceDimension(EpiPseudo1, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
EpiPseudo1 <- clusterCells(EpiPseudo1, num_clusters = 2)

plot_cell_clusters(EpiPseudo1, color_by = "seurat_clusters")

diff_test_res1 <- differentialGeneTest(EpiPseudo1[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes1 <- row.names (subset(diff_test_res1, qval < 0.01))
EpiPseudo1 <- setOrderingFilter(EpiPseudo1, ordering_genes1)
plot_ordering_genes(EpiPseudo1)

EpiPseudo1 <- reduceDimension(EpiPseudo1, max_components = 2,
                             method = 'DDRTree')

EpiPseudo1 <- orderCells(EpiPseudo1)

GM_state <- function(EpiPseudo1){
  if (length(unique(pData(EpiPseudo1)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo1)$State, pData(EpiPseudo1)$seurat_clusters)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo1 <- orderCells(EpiPseudo1, root_state = GM_state(EpiPseudo1))

#Visualization
tiff(file = "EpiPseudo1 Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "EpiPseudo1 Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo1 EpiCellType.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) 
dev.off()
tiff(file = "EpiPseudo1 stim EpiCellType stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo1 CellMarkers Wnt targets purple.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt5", "Krt14", "Krt18", "Mki67", "Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo1 CellMarkers Wnt targets.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt5", "Krt14", "Krt8", "Krt18", "Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "EpiPseudo1 Wnt targets.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo1 IGF1R targets.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, markers = c("Igf1r", "Mapk13", "Jak2", "Fos"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "EpiPseudo1 celltype IGF1R Wnt stem cell targets.tiff", width = 16, height = 21, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, markers = c("ARQ", "EGFP", "Ar", "Krt5", "Krt8", "Igf1r", "Mapk13", "Jak2", "Tcf4", "Ccnd1", "Axin2", "Lgr5", "Trp63", "Ly6a", "Psca", "Tacstd2", "Itga6", "Osr1", "Pbsn"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

#Expression level of Wnt downstreams in plot_genes_branched_pseudotime
EpiPseudo_genes1 <- row.names(subset(fData(EpiPseudo1), gene_short_name %in% c("Tcf4", "Axin2")))
tiff(file = "ARQPosEpi plot_genes_branched_pseudotime Tcf4 Axin2.tiff", width = 4, height = 6, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(EpiPseudo1[EpiPseudo_genes1,],
                               branch_point = 1,
                               color_by = "EpiCellType",
                               ncol = 1)+ scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
dev.off()

EpiPseudo_genes1 <- row.names(subset(fData(EpiPseudo1), gene_short_name %in% c("Ccnd1", "Lgr5")))
tiff(file = "ARQPosEpi plot_genes_branched_pseudotime Ccnd1 Lgr5.tiff", width = 4, height = 6, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(EpiPseudo1[EpiPseudo_genes1,],
                               branch_point = 1,
                               color_by = "EpiCellType",
                               ncol = 1)+ scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
dev.off()

#WT.Epi
Idents(object = WT.Epi) <- "seurat_clusters"
DimPlot(WT.Epi, reduction = "tsne")

DefaultAssay(WT.Epi) <- "RNA"
EpiPseudo1 <- as.CellDataSet(WT.Epi)
EpiPseudo1 <- detectGenes(EpiPseudo1, min_expr = 0.1)
print(head(fData(EpiPseudo1)))

expressed_genes <- row.names(subset(fData(EpiPseudo1),
                                    num_cells_expressed >= 10))

pData(EpiPseudo1)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo1))
EpiPseudo1 <- EpiPseudo1[,pData(EpiPseudo1)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo1), geom =
        "density")

EpiPseudo1 <- estimateSizeFactors(EpiPseudo1)
EpiPseudo1 <- estimateDispersions(EpiPseudo1)

disp_table <- dispersionTable(EpiPseudo1)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo1 <- setOrderingFilter(EpiPseudo1, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo1)

#EpiPseudo1@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo1, return_all = F) # norm_method='log'

EpiPseudo1 <- reduceDimension(EpiPseudo1, max_components = 2, num_dim = 20,
                              reduction_method = 'tSNE', verbose = T)
EpiPseudo1 <- clusterCells(EpiPseudo1, num_clusters = 2)

plot_cell_clusters(EpiPseudo1, color_by = "seurat_clusters")

diff_test_res1 <- differentialGeneTest(EpiPseudo1[expressed_genes,],
                                       fullModelFormulaStr = "~seurat_clusters")
ordering_genes1 <- row.names (subset(diff_test_res1, qval < 0.001))
EpiPseudo1 <- setOrderingFilter(EpiPseudo1, ordering_genes1)
plot_ordering_genes(EpiPseudo1)

EpiPseudo1 <- reduceDimension(EpiPseudo1, max_components = 2,
                              method = 'DDRTree')

EpiPseudo1 <- orderCells(EpiPseudo1)

GM_state <- function(EpiPseudo1){
  if (length(unique(pData(EpiPseudo1)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo1)$State, pData(EpiPseudo1)$seurat_clusters)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo1 <- orderCells(EpiPseudo1, root_state = GM_state(EpiPseudo1))

#Visualization
tiff(file = "EpiPseudo1 Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "EpiPseudo1 seurat_clusters.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "EpiPseudo1 celltype IGF1R Wnt stem cell targets.tiff", width = 16, height = 21, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo1, markers = c("ARQ", "EGFP", "Ar", "Krt5", "Krt8", "Igf1r", "Mapk13", "Jak2", "Tcf4", "Ccnd1", "Axin2", "Lgr5", "Trp63", "Ly6a", "Psca", "Tacstd2", "Itga6", "Osr1", "Pbsn"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

#ARQNegEpi
Idents(object = ARQNegEpi) <- "seurat_clusters"
DefaultAssay(ARQNegEpi) <- "RNA"
EpiPseudo1 <- as.CellDataSet(ARQNegEpi)
EpiPseudo1 <- detectGenes(EpiPseudo1, min_expr = 0.1)
print(head(fData(EpiPseudo1)))


expressed_genes <- row.names(subset(fData(EpiPseudo2),
                                    num_cells_expressed >= 10))

pData(EpiPseudo2)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo2))
EpiPseudo2 <- EpiPseudo2[,pData(EpiPseudo2)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo2), geom =
        "density")

EpiPseudo2 <- estimateSizeFactors(EpiPseudo2)
EpiPseudo2 <- estimateDispersions(EpiPseudo2)

disp_table <- dispersionTable(EpiPseudo2)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo2 <- setOrderingFilter(EpiPseudo2, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo2)

#EpiPseudo2@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo2, return_all = F) # norm_method='log'

EpiPseudo2 <- reduceDimension(EpiPseudo2, max_components = 2, num_dim = 20,
                              reduction_method = 'tSNE', verbose = T)
EpiPseudo2 <- clusterCells(EpiPseudo2, num_clusters = 2)

plot_cell_clusters(EpiPseudo2, color_by = "seurat_clusters")

diff_test_res2 <- differentialGeneTest(EpiPseudo2[expressed_genes,],
                                       fullModelFormulaStr = "~seurat_clusters")
ordering_genes2 <- row.names (subset(diff_test_res2, qval < 0.1))
EpiPseudo2 <- setOrderingFilter(EpiPseudo2, ordering_genes2)
plot_ordering_genes(EpiPseudo2)

EpiPseudo2 <- reduceDimension(EpiPseudo2, max_components = 2,
                              method = 'DDRTree')

EpiPseudo2 <- orderCells(EpiPseudo2)

GM_state <- function(EpiPseudo2){
  if (length(unique(pData(EpiPseudo2)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo2)$State, pData(EpiPseudo2)$seurat_clusters)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo2 <- orderCells(EpiPseudo2, root_state = GM_state(EpiPseudo1))

#Visualization
tiff(file = "EpiPseudo2 Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo2, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()


####Pseudotime using Slingslot####

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.6/pseudotime")

#ARQPosEpi
DefaultAssay(ARQPosEpi) <- "RNA"
Idents(object = ARQPosEpi) <- "seurat_clusters"
DimPlot(ARQPosEpi, reduction = "umap")

ARQPosEpi_sds <- slingshot(Embeddings(ARQPosEpi, "umap"), clusterLabels = ARQPosEpi$seurat_clusters, 
                    start.clus = 4, end.clus = 0, stretch = 2)

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

cell_colors_clust <- cell_pal(ARQPosEpi$seurat_clusters, hue_pal())

plot(reducedDim(ARQPosEpi_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQPosEpi_sds, lwd = 2, type = 'lineages', col = 'black')

tiff(file = "ARQPosEpi slingshot Pseudotime.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot(reducedDim(ARQPosEpi_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQPosEpi_sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()

count.data <- GetAssayData(object = ARQPosEpi[["RNA"]], slot = "counts")

tiff(file = "ARQPosEpi slingshot Ccnd1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plotGenePseudotime(ARQPosEpi_sds, "Ccnd1", count.data)
dev.off()

#
nc <- 3
pt <- slingPseudotime(ARQPosEpi_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQPosEpi_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQPosEpi_sds, lwd = 2, col = 'black', type = 'lineages')
}

tiff(file = "ARQPosEpi slingshot curve1-3.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQPosEpi_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQPosEpi_sds, lwd = 2, col = 'black', type = 'lineages')
}
dev.off()

#
nc <- 3
pt <- slingPseudotime(ARQPosEpi_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQPosEpi_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQPosEpi_sds, lwd = 2, col = 'black', type = 'lineages')
}

tiff(file = "ARQPosEpi slingshot curve1-3.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQPosEpi_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQPosEpi_sds, lwd = 2, col = 'black', type = 'lineages')
}
dev.off()

#ARQNegEpi
DefaultAssay(ARQNegEpi) <- "RNA"
Idents(object = ARQNegEpi) <- "seurat_clusters"
DimPlot(ARQNegEpi, reduction = "umap")

ARQNegEpi_sds <- slingshot(Embeddings(ARQNegEpi, "umap"), clusterLabels = ARQNegEpi$seurat_clusters, 
                           start.clus = 4, stretch = 2)

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

cell_colors_clust <- cell_pal(ARQNegEpi$seurat_clusters, hue_pal())

plot(reducedDim(ARQNegEpi_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQNegEpi_sds, lwd = 2, type = 'lineages', col = 'black')

tiff(file = "ARQNegEpi slingshot Pseudotime.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot(reducedDim(ARQNegEpi_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(ARQNegEpi_sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()

#
nc <- 3
pt <- slingPseudotime(ARQNegEpi_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQNegEpi_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQNegEpi_sds, lwd = 2, col = 'black', type = 'lineages')
}

#
nc <- 1
pt <- slingPseudotime(ARQNegEpi_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(ARQNegEpi_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(ARQNegEpi_sds, lwd = 2, col = 'black', type = 'lineages')
}

saveRDS(ARQPosEpi, file = "ARQPosEpi.rds")
saveRDS(ARQNegEpi, file = "ARQNegEpi.rds")

####Fetchdata####
PINvTumor.combined.Epi2.FetchData  <- FetchData(PINvTumor.combined.Epi2, vars = c("orig.ident","seurat_clusters"))
