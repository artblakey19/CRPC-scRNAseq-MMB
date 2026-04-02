library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HGF-hMET-Pb")

HGF_hMET.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/201120_TGen/new_counts/count_for_163bp_hMETtg2_mm10/count_39553_HGF_hMET_Pb_8M/outs/filtered_feature_bc_matrix")
HGF_hMET <- CreateSeuratObject(counts = HGF_hMET.data,  min.cells = 3, min.features = 200, project = "HGF_hMET")
HGF_hMET <- NormalizeData(HGF_hMET)

ARQ_PCa.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/201120_TGen/new_counts/count_for_163bp_hMETtg2_mm10/28801_s7181_Tumor_count_ref5/outs/filtered_feature_bc_matrix")
ARQ_PCa <- CreateSeuratObject(counts = ARQ_PCa.data,  min.cells = 3, min.features = 200, project = "ARQ_PCa")
ARQ_PCa <- NormalizeData(ARQ_PCa)

//isi-dcnl/user_data/zjsun/seq/201120_TGen/new_counts/count_for_163bp_hMETtg2_mm10/28801_s7181_Tumor_count_ref5/outs/filtered_feature_bc_matrix 

HGF_hMET[["percent.mt"]] <- PercentageFeatureSet(HGF_hMET, pattern = "^mt-")
VlnPlot(HGF_hMET, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(HGF_hMET@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "HGF_hMET Pre-filteration")

hist(HGF_hMET@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "HGF_hMET Pre-filteration")

plot1 <- FeatureScatter(HGF_hMET, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HGF_hMET, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

HGF_hMET_1 <- subset(HGF_hMET, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 10)
VlnPlot(HGF_hMET_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(HGF_hMET_1@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "HGF_hMET_1 Post-filteration")

hist(HGF_hMET_1@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "HGF_hMET_1 Post-filteration")

plot1 <- FeatureScatter(HGF_hMET_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(HGF_hMET_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

HGF_hMET_1 <- FindVariableFeatures(HGF_hMET_1, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(HGF_hMET_1)

#Clustering
HGF_hMET_1 <- ScaleData(HGF_hMET_1, verbose = FALSE)
HGF_hMET_1 <- RunPCA(HGF_hMET_1, npcs = 50, verbose = FALSE)
ElbowPlot(HGF_hMET_1, ndims = 50)

HGF_hMET_1 <- FindNeighbors(HGF_hMET_1, reduction = "pca", dims = 1:19)
HGF_hMET_1 <- FindClusters(HGF_hMET_1, resolution = 0.5)
HGF_hMET_1 <- RunTSNE(HGF_hMET_1, reduction = "pca", dims = 1:19)
HGF_hMET_1 <- RunUMAP(HGF_hMET_1, reduction = "pca", dims = 1:19)
DimPlot(HGF_hMET_1, reduction = "umap", pt.size = 0.3) 

DefaultAssay(HGF_hMET_1) <- "RNA"
FeaturePlot(HGF_hMET_1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

tiff(file = "HGF_hMET hMETtg Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_hMET_1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "HGF_hMET Met Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_hMET_1, reduction = "umap", features = c("Met"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####ARQ_PCa####
ARQ_PCa[["percent.mt"]] <- PercentageFeatureSet(ARQ_PCa, pattern = "^mt-")
VlnPlot(ARQ_PCa, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(ARQ_PCa@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARQ_PCa Pre-filteration")

hist(ARQ_PCa@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARQ_PCa Pre-filteration")

plot1 <- FeatureScatter(ARQ_PCa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ARQ_PCa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ARQ_PCa_1 <- subset(ARQ_PCa, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
VlnPlot(ARQ_PCa_1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(ARQ_PCa_1@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARQ_PCa_1 Post-filteration")

hist(ARQ_PCa_1@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARQ_PCa_1 Post-filteration")

plot1 <- FeatureScatter(ARQ_PCa_1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ARQ_PCa_1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ARQ_PCa_1 <- FindVariableFeatures(ARQ_PCa_1, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(ARQ_PCa_1)

#Clustering
ARQ_PCa_1 <- ScaleData(ARQ_PCa_1, verbose = FALSE)
ARQ_PCa_1 <- RunPCA(ARQ_PCa_1, npcs = 50, verbose = FALSE)
ElbowPlot(ARQ_PCa_1, ndims = 50)

ARQ_PCa_1 <- FindNeighbors(ARQ_PCa_1, reduction = "pca", dims = 1:19)
ARQ_PCa_1 <- FindClusters(ARQ_PCa_1, resolution = 0.5)
ARQ_PCa_1 <- RunTSNE(ARQ_PCa_1, reduction = "pca", dims = 1:19)
ARQ_PCa_1 <- RunUMAP(ARQ_PCa_1, reduction = "pca", dims = 1:19)
DimPlot(ARQ_PCa_1, reduction = "umap", pt.size = 0.3) 

DefaultAssay(ARQ_PCa_1) <- "RNA"
FeaturePlot(ARQ_PCa_1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")


tiff(file = "ARQ_PCa hMETtg Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ_PCa_1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARQ_PCa Met Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ_PCa_1, reduction = "umap", features = c("Met"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()