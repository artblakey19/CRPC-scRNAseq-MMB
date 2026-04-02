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
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(scater)
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)
library(irlba)

####Loading data####
##E18.5
E181.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/counts/36982_6-E18-WT-Gli1/outs/filtered_feature_bc_matrix")
E181 <- CreateSeuratObject(counts = E181.data,  min.cells = 3, min.features = 200, project = "E181")
E181 <- NormalizeData(E181)
E182.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/20190417_Tgen/28126_scRNA/counts/28126_count/outs/filtered_feature_bc_matrix")
E182 <- CreateSeuratObject(counts = E182.data,  min.cells = 3, min.features = 200, project = "E182")
E182 <- NormalizeData(E182)
E183.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Adam Olson/AO Files/ARKO Gli1/Single Cell/WT/filtered_feature_bc_matrix")
E183 <- CreateSeuratObject(counts = E183.data,  min.cells = 3, min.features = 200, project = "E183")
E183 <- NormalizeData(E183)

##P14
P141.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191214_32848_33580_test/counts/32848_ARCreER_P14/outs/filtered_feature_bc_matrix")
P141 <- CreateSeuratObject(counts = P141.data,  min.cells = 3, min.features = 200, project = "P141")
P141 <- NormalizeData(P141)
P142.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/counts/36977_1-P14-WT-Gli1/outs/filtered_feature_bc_matrix")
P142 <- CreateSeuratObject(counts = P142.data,  min.cells = 3, min.features = 200, project = "P142")
P142 <- NormalizeData(P142)
P143.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/190927/31923_ARCrER/counts/31923_ARCrER/outs/filtered_feature_bc_matrix")
P143 <- CreateSeuratObject(counts = P143.data,  min.cells = 3, min.features = 200, project = "P143")
P143 <- NormalizeData(P143)

##P35
P351.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/190219_scRNA/27347_WT_count_EGFPmm10/outs/filtered_feature_bc_matrix")
P351 <- CreateSeuratObject(counts = P351.data,  min.cells = 3, min.features = 200, project = "P351")
P351 <- NormalizeData(P351)
P352.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191011_32414_test/32414_ARCreER_P35/32414_ARCreER_P35/outs/filtered_feature_bc_matrix")
P352 <- CreateSeuratObject(counts = P352.data,  min.cells = 3, min.features = 200, project = "P352")
P352 <- NormalizeData(P352)
P353.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/36980_4-P35-WT-Gli1/36980_count_ref5/outs/filtered_feature_bc_matrix")
P353 <- CreateSeuratObject(counts = P353.data,  min.cells = 3, min.features = 200, project = "P353")
P353 <- NormalizeData(P353)

####Initial processing, Filtering and Clustering####
##E181
E181[["percent.mt"]] <- PercentageFeatureSet(E181, pattern = "^mt-")
VlnPlot(E181, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(E181@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E181 Pre-filteration")
hist(E181@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E181 Pre-filteration")
E1811 <- subset(E181, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
VlnPlot(E1811, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(E1811@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E1811 Post-filteration")
hist(E1811@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E1811 Post-filteration")
E1811 <- FindVariableFeatures(E1811, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(E1811)
E1811 <- ScaleData(E1811, features = all.genes)
E1811 <- RunPCA(E1811, features = VariableFeatures(object = E1811))
print(E1811[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E1811, dims = 1:2, reduction = "pca")
DimPlot(E1811, reduction = "pca")
ElbowPlot(E1811, ndims = 50)
E1811 <- FindNeighbors(E1811, dims = 1:20)
E1811 <- FindClusters(E1811, resolution = 0.5)
head(Idents(E1811), 5)
E1811 <- RunTSNE(E1811, dims = 1:20)
E1811 <- RunUMAP(E1811, dims = 1:20)
DimPlot(E1811, reduction = "umap", pt.size = 1, label = TRUE)

##E182
E182[["percent.mt"]] <- PercentageFeatureSet(E182, pattern = "^mt-")
VlnPlot(E182, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(E182@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E182 Pre-filteration")
hist(E182@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E182 Pre-filteration")
plot1 <- FeatureScatter(E182, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E182, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
E1821 <- subset(E182, subset = nFeature_RNA> 450 & nFeature_RNA < 8000 & percent.mt < 10)
VlnPlot(E1821, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(E1821@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E1821 Post-filteration")
hist(E1821@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E1821 Post-filteration")
plot1 <- FeatureScatter(E1821, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E1821, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
E1821 <- FindVariableFeatures(E1821, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(E1821)
E1821 <- ScaleData(E1821, features = all.genes)
E1821 <- RunPCA(E1821, features = VariableFeatures(object = E1821))
print(E1821[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E1821, dims = 1:2, reduction = "pca")
DimPlot(E1821, reduction = "pca")
ElbowPlot(E1821, ndims = 50)
E1821 <- FindNeighbors(E1821, dims = 1:20)
E1821 <- FindClusters(E1821, resolution = 0.5)
head(Idents(E1821), 5)
E1821 <- RunTSNE(E1821, dims = 1:20)
E1821 <- RunUMAP(E1821, dims = 1:20)
DimPlot(E1821, reduction = "umap", pt.size = 1, label = TRUE)

##E183
E183[["percent.mt"]] <- PercentageFeatureSet(E183, pattern = "^mt-")
VlnPlot(E183, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(E183@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E183 Pre-filteration")
hist(E183@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E183 Pre-filteration")
plot1 <- FeatureScatter(E183, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E183, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
E1831 <- subset(E183, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 10)
VlnPlot(E1831, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(E1831@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E1831 Post-filteration")
hist(E1831@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E1831 Post-filteration")
plot1 <- FeatureScatter(E1831, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E1831, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
DefaultAssay(E1831) <- "RNA"
E1831 <- NormalizeData(E1831, verbose = FALSE)
E1831 <- FindVariableFeatures(E1831, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(E1831)
E1831 <- ScaleData(E1831, features = all.genes)
E1831 <- RunPCA(E1831, features = VariableFeatures(E1831))
print(E1831[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E1831, dims = 1:2, reduction = "pca")
DimPlot(E1831, reduction = "pca")
ElbowPlot(E1831, ndims = 50)
E1831 <- FindNeighbors(E1831, dims = 1:20)
E1831 <- FindClusters(E1831, resolution = 0.5)
head(Idents(E1831), 5)
E1831 <- RunTSNE(E1831, dims = 1:20)
E1831 <- RunUMAP(E1831, dims = 1:20)
DimPlot(E1831, reduction = "umap", pt.size = 1, label = TRUE)

##P141
P141[["percent.mt"]] <- PercentageFeatureSet(P141, pattern = "^mt-")
VlnPlot(P141, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P141@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P141 Pre-filteration")
hist(P141@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P141 Pre-filteration")
P1411 <- subset(P141, subset =  nFeature_RNA < 5000 & percent.mt < 10)
VlnPlot(P1411, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P1411@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P1411 Post-filteration")
hist(P1411@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P1411 Post-filteration")
plot1 <- FeatureScatter(P1411, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P1411, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P1411 <- FindVariableFeatures(P1411, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P1411)
P1411 <- ScaleData(P1411, features = all.genes)
P1411 <- RunPCA(P1411, features = VariableFeatures(object = E1811))
print(P1411[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P1411, dims = 1:2, reduction = "pca")
DimPlot(P1411, reduction = "pca")
ElbowPlot(P1411, ndims = 50)
P1411 <- FindNeighbors(P1411, dims = 1:20)
P1411 <- FindClusters(P1411, resolution = 0.5)
head(Idents(P1411), 5)
P1411 <- RunTSNE(P1411, dims = 1:20)
P1411 <- RunUMAP(P1411, dims = 1:20)
DimPlot(P1411, reduction = "umap", pt.size = 1, label = TRUE)

##P142
P142[["percent.mt"]] <- PercentageFeatureSet(P142, pattern = "^mt-")
VlnPlot(P142, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P142@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P142 Pre-filteration")
hist(P142@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P142 Pre-filteration")
P1421 <- subset(P142, subset =  nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(P1421, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P1421@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P1421 Post-filteration")
hist(P1421@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P1421 Post-filteration")
plot1 <- FeatureScatter(P1421, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P1421, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
DefaultAssay(P1421) <- "RNA"
P1421 <- NormalizeData(P1421, verbose = FALSE)
P1421 <- FindVariableFeatures(P1421, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P1421)
P1421 <- ScaleData(P1421, features = all.genes)
P1421 <- RunPCA(P1421, features = VariableFeatures(P1421))
print(P1421[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P1421, dims = 1:2, reduction = "pca")
DimPlot(P1421, reduction = "pca")
ElbowPlot(P1421, ndims = 50)
P1421 <- FindNeighbors(P1421, dims = 1:20)
P1421 <- FindClusters(P1421, resolution = 0.5)
head(Idents(P1421), 5)
P1421 <- RunTSNE(P1421, dims = 1:20)
P1421 <- RunUMAP(P1421, dims = 1:20)
DimPlot(P1421, reduction = "umap", pt.size = 1, label = TRUE)

##P143
P143[["percent.mt"]] <- PercentageFeatureSet(P143, pattern = "^mt-")
VlnPlot(P143, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P143@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P143 Pre-filteration")
hist(P143@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P143 Pre-filteration")
P1431 <- subset(P143, subset =  nFeature_RNA < 2000 & percent.mt < 10)
VlnPlot(P1431, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P1431@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P1431 Post-filteration")
hist(P1431@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P1431 Post-filteration")
plot1 <- FeatureScatter(P1431, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P1431, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P1431 <- FindVariableFeatures(P1431, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P1431)
P1431 <- ScaleData(P1431, features = all.genes)
P1431 <- RunPCA(P1431, features = VariableFeatures(object = E1811))
print(P1431[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P1431, dims = 1:2, reduction = "pca")
DimPlot(P1431, reduction = "pca")
ElbowPlot(P1431, ndims = 50)
P1431 <- FindNeighbors(P1431, dims = 1:20)
P1431 <- FindClusters(P1431, resolution = 0.5)
head(Idents(P1431), 5)
P1431 <- RunTSNE(P1431, dims = 1:20)
P1431 <- RunUMAP(P1431, dims = 1:20)
DimPlot(P1431, reduction = "umap", pt.size = 1, label = TRUE)

##P351
P351[["percent.mt"]] <- PercentageFeatureSet(P351, pattern = "^mt-")
VlnPlot(P351, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P351@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P351 Pre-filteration")
hist(P351@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P351 Pre-filteration")
plot1 <- FeatureScatter(P351, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P351, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P3511 <- subset(P351, subset = nFeature_RNA > 500 & nFeature_RNA < 8500 & percent.mt < 15)
VlnPlot(P3511, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P3511@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P3511 Post-filteration")
hist(P3511@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P3511 Post-filteration")
plot1 <- FeatureScatter(P3511, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3511, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P3511 <- FindVariableFeatures(P3511, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P3511)
P3511 <- ScaleData(P3511, features = all.genes)
P3511 <- RunPCA(P3511, features = VariableFeatures(object = E1811))
print(P3511[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P3511, dims = 1:2, reduction = "pca")
DimPlot(P3511, reduction = "pca")
ElbowPlot(P3511, ndims = 50)
P3511 <- FindNeighbors(P3511, dims = 1:20)
P3511 <- FindClusters(P3511, resolution = 0.5)
head(Idents(P3511), 5)
P3511 <- RunTSNE(P3511, dims = 1:20)
P3511 <- RunUMAP(P3511, dims = 1:20)
DimPlot(P3511, reduction = "umap", pt.size = 1, label = TRUE)

##P352
P352[["percent.mt"]] <- PercentageFeatureSet(P352, pattern = "^mt-")
VlnPlot(P352, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P352@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P352 Pre-filteration")
hist(P352@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P352 Pre-filteration")
plot1 <- FeatureScatter(P352, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P352, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P3521 <- subset(P352, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(P3521, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P3521@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P3521 Post-filteration")
hist(P3521@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P3521 Post-filteration")
plot1 <- FeatureScatter(P3521, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3521, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
DefaultAssay(P3521) <- "RNA"
P3521 <- NormalizeData(P3521, verbose = FALSE)
P3521 <- FindVariableFeatures(P3521, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P3521)
P3521 <- ScaleData(P3521, features = all.genes)
P3521 <- RunPCA(P3521, features = VariableFeatures(P3521))
print(P3521[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P3521, dims = 1:2, reduction = "pca")
DimPlot(P3521, reduction = "pca")
ElbowPlot(P3521, ndims = 50)
P3521 <- FindNeighbors(P3521, dims = 1:20)
P3521 <- FindClusters(P3521, resolution = 0.5)
head(Idents(P3521), 5)
P3521 <- RunTSNE(P3521, dims = 1:20)
P3521 <- RunUMAP(P3521, dims = 1:20)
DimPlot(P3521, reduction = "umap", pt.size = 1, label = TRUE)

#P353
P353[["percent.mt"]] <- PercentageFeatureSet(P353, pattern = "^mt-")
VlnPlot(P353, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P353@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P353 Pre-filteration")
hist(P353@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P353 Pre-filteration")
plot1 <- FeatureScatter(P353, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P353, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P3531 <- subset(P353, subset = nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(P3531, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(P3531@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P3531 Post-filteration")
hist(P3531@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P3531 Post-filteration")
plot1 <- FeatureScatter(P3531, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3531, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
P3531 <- FindVariableFeatures(P3531, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P3531)
P3531 <- ScaleData(P3531, features = all.genes)
P3531 <- RunPCA(P3531, features = VariableFeatures(object = E1811))
print(P3531[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P3531, dims = 1:2, reduction = "pca")
DimPlot(P3531, reduction = "pca")
ElbowPlot(P3531, ndims = 50)
P3531 <- FindNeighbors(P3531, dims = 1:20)
P3531 <- FindClusters(P3531, resolution = 0.5)
head(Idents(P3531), 5)
P3531 <- RunTSNE(P3531, dims = 1:20)
P3531 <- RunUMAP(P3531, dims = 1:20)
DimPlot(P3531, reduction = "umap", pt.size = 1, label = TRUE)

####Merging Datasets####
#Stash old idents
E1831[["orig.clusters"]] <- Idents(object = E1831)
P1121[["orig.clusters"]] <- Idents(object = P1121)
P3521[["orig.clusters"]] <- Idents(object = P3521)
#Set Current idents
Idents(object = E1831) <- "seurat_clusters"
Idents(object = P1121) <- "seurat_clusters"
Idents(object = P3521) <- "seurat_clusters"
E1831$stim <- "E18.5"
P1121$stim <- "P14"
P3521$stim <- "P35"
Osr1.combined.anchors <- FindIntegrationAnchors(object.list = list(E1831, P1121, P3521), dims = 1:20)
Osr1.combined <- IntegrateData(anchorset = Osr1.combined.anchors, dims = 1:20)
DefaultAssay(Osr1.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
Osr1.combined <- ScaleData(Osr1.combined, verbose = FALSE)
Osr1.combined <- RunPCA(Osr1.combined, npcs = 30, verbose = FALSE)
ElbowPlot(Osr1.combined, ndims = 50)
#tSNE, UMAP and Clustering
Osr1.combined <- FindNeighbors(Osr1.combined, reduction = "pca", dims = 1:16)
Osr1.combined <- FindClusters(Osr1.combined, resolution = 0.5)
Osr1.combined <- RunTSNE(Osr1.combined, reduction = "pca", dims = 1:16)
Osr1.combined <- RunUMAP(Osr1.combined, reduction = "pca", dims = 1:16)
DimPlot(Osr1.combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("darkred", "darkblue", "grey")) 
#Cell Cycle Regression
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes
DefaultAssay(Osr1.combined) <- "RNA"
all.genes <- rownames(Osr1.combined)
Osr1.combined <- ScaleData(Osr1.combined, features = all.genes)
Osr1.combined <- CellCycleScoring(Osr1.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(object = Osr1.combined) <- "Phase"
DimPlot(Osr1.combined, reduction = "umap", pt.size = 0.3, cols = c("purple", "#FFD966","#1D762E"))
Osr1.combined2 <- Osr1.combined
DefaultAssay(Osr1.combined2) <- "integrated"
Osr1.combined2 <- ScaleData(Osr1.combined2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Osr1.combined1))
Osr1.combined2 <- RunPCA(Osr1.combined2, features = VariableFeatures(Osr1.combined2))
ElbowPlot(Osr1.combined2, ndims = 30)
Osr1.combined2 <- FindNeighbors(Osr1.combined2, reduction = "pca", dims = 1:16)
Osr1.combined2 <- FindClusters(Osr1.combined2, resolution = 0.5)
Osr1.combined2 <- RunUMAP(Osr1.combined2, reduction = "pca", dims = 1:16)
Osr1.combined2 <- RunTSNE(Osr1.combined2, reduction = "pca", dims = 1:16)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3, label = TRUE)
Idents(object = Osr1.combined2) <- "Phase"
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3, cols = c("purple", "#FFD966","#1D762E"))
#New labeling
Idents(object = Osr1.combined2) <- "seurat_clusters"
Osr1.combined2 <- RenameIdents(object = Osr1.combined2, '2' = "BE", '5' = "BE", '4' = "LE", '1' = "FB", '3' = "FB", '8' = "FB", '0' = "SM", '6' = "SM", '11' = "Peri", '10' = "Leu", '9' = "VE", '7' = "Gli1", '12' = "SV", '13' = "NE")
Osr1.combined2[["CellType"]] <- Idents(object = Osr1.combined2)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3)
#UMAP
Idents(object = Osr1.combined2) <- "CellType"
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.5, cols = c("#1D762E", "#E06666",  "red", "#FF9933", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.5, split.by = "stim", cols = c("#1D762E", "#E06666",  "red", "#FF9933", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
#Dotplots
Idents(object = Osr1.combined2) <- "CellType"
tiff(file = "Osr1.combined2 CellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(Osr1.combined2, features = c("Krt14", "Krt15", "Krt17", "Krt5", "Krt6a", "Krt19", "Tspan1", "Krt8", "Clu", "Cldn3", "Igfbp3", "Fbln1", "Lum", "Bgn", "Col15a1", "Myl9", "Acta2", "Tagln", "Actg2", "Myh11", "Rgs5", "Ndufa4l2", "Apold1", "Mustn1",  "Lrrc32", "Cd74", "Ccl4", "C1qa", "C1qb", "Tyrobp", "Cldn5", "Cdh5", "Pecam1", "Plvap", "Aqp1", "Mpz", "Plp1", "Fabp7", "Cryab", "Gpm6b", "Pate4", "Svs2", "Svs4", "Wfdc15b", "Plac8", "Scg2", "Chgb", "Chga", "Syp", "Scg3"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()
#Featureplots
DefaultAssay(Osr1.combined2) <- "RNA"
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Myh11"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Osr1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
#Blended expression UMAP plots
Idents(object = Osr1.combined2) <- "stim"
E18only <- subset(Osr1.combined2, idents = c("E18.5"))
DefaultAssay(E18only) <- "RNA"
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Ly6a"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Tacstd2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Psca"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Itga6"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)

####% of indicated markers####
##E181
#Add Osr1 info
DefaultAssay(E1811) <- "RNA"
E1811Osr1Pos <- subset(x=E1811, subset = Osr1 > 0)
E1811Osr1Neg <- subset(x=E1811, subset = Osr1 == 0)
Idents(object = E1811Osr1Pos) <- "Osr1Pos"
Idents(object = E1811Osr1Neg) <- "Osr1Neg"
E1811Osr1Pos[["Osr1Exp"]] <- Idents(object = E1811Osr1Pos)
E1811Osr1Neg[["Osr1Exp"]] <- Idents(object = E1811Osr1Neg)
E1811Osr1 <- merge(x = E1811Osr1Pos, y = E1811Osr1Neg)
Idents(object = E1811Osr1) <- "Osr1Exp"
E1811$Osr1Exp <- Idents(object = E1811Osr1)
#Add Ly6a info
DefaultAssay(E1811) <- "RNA"
E1811Ly6aPos <- subset(x=E1811, subset = Ly6a > 0)
E1811Ly6aNeg <- subset(x=E1811, subset = Ly6a == 0)
Idents(object = E1811Ly6aPos) <- "Ly6aPos"
Idents(object = E1811Ly6aNeg) <- "Ly6aNeg"
E1811Ly6aPos[["Ly6aExp"]] <- Idents(object = E1811Ly6aPos)
E1811Ly6aNeg[["Ly6aExp"]] <- Idents(object = E1811Ly6aNeg)
E1811Ly6a <- merge(x = E1811Ly6aPos, y = E1811Ly6aNeg)
Idents(object = E1811Ly6a) <- "Ly6aExp"
E1811$Ly6aExp <- Idents(object = E1811Ly6a)
#Add Tacstd2 info
DefaultAssay(E1811) <- "RNA"
E1811Tacstd2Pos <- subset(x=E1811, subset = Tacstd2 > 0)
E1811Tacstd2Neg <- subset(x=E1811, subset = Tacstd2 == 0)
Idents(object = E1811Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = E1811Tacstd2Neg) <- "Tacstd2Neg"
E1811Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = E1811Tacstd2Pos)
E1811Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = E1811Tacstd2Neg)
E1811Tacstd2 <- merge(x = E1811Tacstd2Pos, y = E1811Tacstd2Neg)
Idents(object = E1811Tacstd2) <- "Tacstd2Exp"
E1811$Tacstd2Exp <- Idents(object = E1811Tacstd2)
#Add Psca info
DefaultAssay(E1811) <- "RNA"
E1811PscaPos <- subset(x=E1811, subset = Psca > 0)
E1811PscaNeg <- subset(x=E1811, subset = Psca == 0)
Idents(object = E1811PscaPos) <- "PscaPos"
Idents(object = E1811PscaNeg) <- "PscaNeg"
E1811PscaPos[["PscaExp"]] <- Idents(object = E1811PscaPos)
E1811PscaNeg[["PscaExp"]] <- Idents(object = E1811PscaNeg)
E1811Psca <- merge(x = E1811PscaPos, y = E1811PscaNeg)
Idents(object = E1811Psca) <- "PscaExp"
E1811$PscaExp <- Idents(object = E1811Psca)
#Total
Idents(object = E1811) <- "Osr1Exp"
table(Idents(E1811))
Idents(object = E1811) <- "Tacstd2Exp"
table(Idents(E1811))
Idents(object = E1811) <- "PscaExp"
table(Idents(E1811))
Idents(object = E1811) <- "Ly6aExp"
table(Idents(E1811))
#Epithelia
Idents(object = E1811) <- "seurat_clusters"
E1811.Epi <- subset(E1811, idents = c("4", "6", "7"))
Idents(object = E1811.Epi) <- "Osr1Exp"
table(Idents(E1811.Epi))
Idents(object = E1811.Epi) <- "Tacstd2Exp"
table(Idents(E1811.Epi))
Idents(object = E1811.Epi) <- "PscaExp"
table(Idents(E1811.Epi))
Idents(object = E1811.Epi) <- "Ly6aExp"
table(Idents(E1811.Epi))

##E182
#Add Osr1 info
DefaultAssay(E1821) <- "RNA"
E1821Osr1Pos <- subset(x=E1821, subset = Osr1 > 0)
E1821Osr1Neg <- subset(x=E1821, subset = Osr1 == 0)
Idents(object = E1821Osr1Pos) <- "Osr1Pos"
Idents(object = E1821Osr1Neg) <- "Osr1Neg"
E1821Osr1Pos[["Osr1Exp"]] <- Idents(object = E1821Osr1Pos)
E1821Osr1Neg[["Osr1Exp"]] <- Idents(object = E1821Osr1Neg)
E1821Osr1 <- merge(x = E1821Osr1Pos, y = E1821Osr1Neg)
Idents(object = E1821Osr1) <- "Osr1Exp"
E1821$Osr1Exp <- Idents(object = E1821Osr1)
#Add Ly6a info
DefaultAssay(E1821) <- "RNA"
E1821Ly6aPos <- subset(x=E1821, subset = Ly6a > 0)
E1821Ly6aNeg <- subset(x=E1821, subset = Ly6a == 0)
Idents(object = E1821Ly6aPos) <- "Ly6aPos"
Idents(object = E1821Ly6aNeg) <- "Ly6aNeg"
E1821Ly6aPos[["Ly6aExp"]] <- Idents(object = E1821Ly6aPos)
E1821Ly6aNeg[["Ly6aExp"]] <- Idents(object = E1821Ly6aNeg)
E1821Ly6a <- merge(x = E1821Ly6aPos, y = E1821Ly6aNeg)
Idents(object = E1821Ly6a) <- "Ly6aExp"
E1821$Ly6aExp <- Idents(object = E1821Ly6a)
#Add Tacstd2 info
DefaultAssay(E1821) <- "RNA"
E1821Tacstd2Pos <- subset(x=E1821, subset = Tacstd2 > 0)
E1821Tacstd2Neg <- subset(x=E1821, subset = Tacstd2 == 0)
Idents(object = E1821Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = E1821Tacstd2Neg) <- "Tacstd2Neg"
E1821Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = E1821Tacstd2Pos)
E1821Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = E1821Tacstd2Neg)
E1821Tacstd2 <- merge(x = E1821Tacstd2Pos, y = E1821Tacstd2Neg)
Idents(object = E1821Tacstd2) <- "Tacstd2Exp"
E1821$Tacstd2Exp <- Idents(object = E1821Tacstd2)
#Add Psca info
DefaultAssay(E1821) <- "RNA"
E1821PscaPos <- subset(x=E1821, subset = Psca > 0)
E1821PscaNeg <- subset(x=E1821, subset = Psca == 0)
Idents(object = E1821PscaPos) <- "PscaPos"
Idents(object = E1821PscaNeg) <- "PscaNeg"
E1821PscaPos[["PscaExp"]] <- Idents(object = E1821PscaPos)
E1821PscaNeg[["PscaExp"]] <- Idents(object = E1821PscaNeg)
E1821Psca <- merge(x = E1821PscaPos, y = E1821PscaNeg)
Idents(object = E1821Psca) <- "PscaExp"
E1821$PscaExp <- Idents(object = E1821Psca)
#Total
Idents(object = E1821) <- "Osr1Exp"
table(Idents(E1821))
Idents(object = E1821) <- "Tacstd2Exp"
table(Idents(E1821))
Idents(object = E1821) <- "PscaExp"
table(Idents(E1821))
Idents(object = E1821) <- "Ly6aExp"
table(Idents(E1821))
#Epithelia
Idents(object = E1821) <- "seurat_clusters"
E1821.Epi <- subset(E1821, idents = c("5", "9"))
Idents(object = E1821.Epi) <- "Osr1Exp"
table(Idents(E1821.Epi))
Idents(object = E1821.Epi) <- "Tacstd2Exp"
table(Idents(E1821.Epi))
Idents(object = E1821.Epi) <- "PscaExp"
table(Idents(E1821.Epi))
Idents(object = E1821.Epi) <- "Ly6aExp"
table(Idents(E1821.Epi))

##E183
#Add Osr1 info
DefaultAssay(E1831) <- "RNA"
E1831Osr1Pos <- subset(x=E1831, subset = Osr1 > 0)
E1831Osr1Neg <- subset(x=E1831, subset = Osr1 == 0)
Idents(object = E1831Osr1Pos) <- "Osr1Pos"
Idents(object = E1831Osr1Neg) <- "Osr1Neg"
E1831Osr1Pos[["Osr1Exp"]] <- Idents(object = E1831Osr1Pos)
E1831Osr1Neg[["Osr1Exp"]] <- Idents(object = E1831Osr1Neg)
E1831Osr1 <- merge(x = E1831Osr1Pos, y = E1831Osr1Neg)
Idents(object = E1831Osr1) <- "Osr1Exp"
E1831$Osr1Exp <- Idents(object = E1831Osr1)
#Add Ly6a info
DefaultAssay(E1831) <- "RNA"
E1831Ly6aPos <- subset(x=E1831, subset = Ly6a > 0)
E1831Ly6aNeg <- subset(x=E1831, subset = Ly6a == 0)
Idents(object = E1831Ly6aPos) <- "Ly6aPos"
Idents(object = E1831Ly6aNeg) <- "Ly6aNeg"
E1831Ly6aPos[["Ly6aExp"]] <- Idents(object = E1831Ly6aPos)
E1831Ly6aNeg[["Ly6aExp"]] <- Idents(object = E1831Ly6aNeg)
E1831Ly6a <- merge(x = E1831Ly6aPos, y = E1831Ly6aNeg)
Idents(object = E1831Ly6a) <- "Ly6aExp"
E1831$Ly6aExp <- Idents(object = E1831Ly6a)
#Add Tacstd2 info
DefaultAssay(E1831) <- "RNA"
E1831Tacstd2Pos <- subset(x=E1831, subset = Tacstd2 > 0)
E1831Tacstd2Neg <- subset(x=E1831, subset = Tacstd2 == 0)
Idents(object = E1831Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = E1831Tacstd2Neg) <- "Tacstd2Neg"
E1831Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = E1831Tacstd2Pos)
E1831Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = E1831Tacstd2Neg)
E1831Tacstd2 <- merge(x = E1831Tacstd2Pos, y = E1831Tacstd2Neg)
Idents(object = E1831Tacstd2) <- "Tacstd2Exp"
E1831$Tacstd2Exp <- Idents(object = E1831Tacstd2)
#Add Psca info
DefaultAssay(E1831) <- "RNA"
E1831PscaPos <- subset(x=E1831, subset = Psca > 0)
E1831PscaNeg <- subset(x=E1831, subset = Psca == 0)
Idents(object = E1831PscaPos) <- "PscaPos"
Idents(object = E1831PscaNeg) <- "PscaNeg"
E1831PscaPos[["PscaExp"]] <- Idents(object = E1831PscaPos)
E1831PscaNeg[["PscaExp"]] <- Idents(object = E1831PscaNeg)
E1831Psca <- merge(x = E1831PscaPos, y = E1831PscaNeg)
Idents(object = E1831Psca) <- "PscaExp"
E1831$PscaExp <- Idents(object = E1831Psca)
#Total
Idents(object = E1831) <- "Osr1Exp"
table(Idents(E1831))
Idents(object = E1831) <- "Tacstd2Exp"
table(Idents(E1831))
Idents(object = E1831) <- "PscaExp"
table(Idents(E1831))
Idents(object = E1831) <- "Ly6aExp"
table(Idents(E1831))
#Epithelia
Idents(object = E1831) <- "seurat_clusters"
E1831.Epi <- subset(E1831, idents = c("1"))
Idents(object = E1831.Epi) <- "Osr1Exp"
table(Idents(E1831.Epi))
Idents(object = E1831.Epi) <- "Tacstd2Exp"
table(Idents(E1831.Epi))
Idents(object = E1831.Epi) <- "PscaExp"
table(Idents(E1831.Epi))
Idents(object = E1831.Epi) <- "Ly6aExp"
table(Idents(E1831.Epi))

##P141
#Add Osr1 info
DefaultAssay(P1411) <- "RNA"
P1411Osr1Pos <- subset(x=P1411, subset = Osr1 > 0)
P1411Osr1Neg <- subset(x=P1411, subset = Osr1 == 0)
Idents(object = P1411Osr1Pos) <- "Osr1Pos"
Idents(object = P1411Osr1Neg) <- "Osr1Neg"
P1411Osr1Pos[["Osr1Exp"]] <- Idents(object = P1411Osr1Pos)
P1411Osr1Neg[["Osr1Exp"]] <- Idents(object = P1411Osr1Neg)
P1411Osr1 <- merge(x = P1411Osr1Pos, y = P1411Osr1Neg)
Idents(object = P1411Osr1) <- "Osr1Exp"
P1411$Osr1Exp <- Idents(object = P1411Osr1)
#Add Ly6a info
DefaultAssay(P1411) <- "RNA"
P1411Ly6aPos <- subset(x=P1411, subset = Ly6a > 0)
P1411Ly6aNeg <- subset(x=P1411, subset = Ly6a == 0)
Idents(object = P1411Ly6aPos) <- "Ly6aPos"
Idents(object = P1411Ly6aNeg) <- "Ly6aNeg"
P1411Ly6aPos[["Ly6aExp"]] <- Idents(object = P1411Ly6aPos)
P1411Ly6aNeg[["Ly6aExp"]] <- Idents(object = P1411Ly6aNeg)
P1411Ly6a <- merge(x = P1411Ly6aPos, y = P1411Ly6aNeg)
Idents(object = P1411Ly6a) <- "Ly6aExp"
P1411$Ly6aExp <- Idents(object = P1411Ly6a)
#Add Tacstd2 info
DefaultAssay(P1411) <- "RNA"
P1411Tacstd2Pos <- subset(x=P1411, subset = Tacstd2 > 0)
P1411Tacstd2Neg <- subset(x=P1411, subset = Tacstd2 == 0)
Idents(object = P1411Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = P1411Tacstd2Neg) <- "Tacstd2Neg"
P1411Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = P1411Tacstd2Pos)
P1411Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = P1411Tacstd2Neg)
P1411Tacstd2 <- merge(x = P1411Tacstd2Pos, y = P1411Tacstd2Neg)
Idents(object = P1411Tacstd2) <- "Tacstd2Exp"
P1411$Tacstd2Exp <- Idents(object = P1411Tacstd2)
#Add Psca info
DefaultAssay(P1411) <- "RNA"
P1411PscaPos <- subset(x=P1411, subset = Psca > 0)
P1411PscaNeg <- subset(x=P1411, subset = Psca == 0)
Idents(object = P1411PscaPos) <- "PscaPos"
Idents(object = P1411PscaNeg) <- "PscaNeg"
P1411PscaPos[["PscaExp"]] <- Idents(object = P1411PscaPos)
P1411PscaNeg[["PscaExp"]] <- Idents(object = P1411PscaNeg)
P1411Psca <- merge(x = P1411PscaPos, y = P1411PscaNeg)
Idents(object = P1411Psca) <- "PscaExp"
P1411$PscaExp <- Idents(object = P1411Psca)
#Total
Idents(object = P1411) <- "Osr1Exp"
table(Idents(P1411))
Idents(object = P1411) <- "Tacstd2Exp"
table(Idents(P1411))
Idents(object = P1411) <- "PscaExp"
table(Idents(P1411))
Idents(object = P1411) <- "Ly6aExp"
table(Idents(P1411))
#Epithelia
Idents(object = P1411) <- "seurat_clusters"
P1411.Epi <- subset(P1411, idents = c("1", "4", "13"))
Idents(object = P1411.Epi) <- "Osr1Exp"
table(Idents(P1411.Epi))
Idents(object = P1411.Epi) <- "Tacstd2Exp"
table(Idents(P1411.Epi))
Idents(object = P1411.Epi) <- "PscaExp"
table(Idents(P1411.Epi))
Idents(object = P1411.Epi) <- "Ly6aExp"
table(Idents(P1411.Epi))

##P142
#Add Osr1 info
DefaultAssay(P142) <- "RNA"
P142Osr1Pos <- subset(x=P142, subset = Osr1 > 0)
P1421Osr1Neg <- subset(x=P142, subset = Osr1 == 0)
Idents(object = P1421Osr1Pos) <- "Osr1Pos"
Idents(object = P1421Osr1Neg) <- "Osr1Neg"
P1421Osr1Pos[["Osr1Exp"]] <- Idents(object = P1421Osr1Pos)
P1421Osr1Neg[["Osr1Exp"]] <- Idents(object = P1421Osr1Neg)
P1421Osr1 <- merge(x = P1421Osr1Pos, y = P1421Osr1Neg)
Idents(object = P1421Osr1) <- "Osr1Exp"
P1421$Osr1Exp <- Idents(object = P1421Osr1)
#Add Ly6a info
DefaultAssay(P1421) <- "RNA"
P1421Ly6aPos <- subset(x=P1421, subset = Ly6a > 0)
P1421Ly6aNeg <- subset(x=P1421, subset = Ly6a == 0)
Idents(object = P1421Ly6aPos) <- "Ly6aPos"
Idents(object = P1421Ly6aNeg) <- "Ly6aNeg"
P1421Ly6aPos[["Ly6aExp"]] <- Idents(object = P1421Ly6aPos)
P1421Ly6aNeg[["Ly6aExp"]] <- Idents(object = P1421Ly6aNeg)
P1421Ly6a <- merge(x = P1421Ly6aPos, y = P1421Ly6aNeg)
Idents(object = P1421Ly6a) <- "Ly6aExp"
P1421$Ly6aExp <- Idents(object = P1421Ly6a)
#Add Tacstd2 info
DefaultAssay(P1421) <- "RNA"
P1421Tacstd2Pos <- subset(x=P1421, subset = Tacstd2 > 0)
P1421Tacstd2Neg <- subset(x=P1421, subset = Tacstd2 == 0)
Idents(object = P1421Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = P1421Tacstd2Neg) <- "Tacstd2Neg"
P1421Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = P1421Tacstd2Pos)
P1421Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = P1421Tacstd2Neg)
P1421Tacstd2 <- merge(x = P1421Tacstd2Pos, y = P1421Tacstd2Neg)
Idents(object = P1421Tacstd2) <- "Tacstd2Exp"
P1421$Tacstd2Exp <- Idents(object = P1421Tacstd2)
#Add Psca info
DefaultAssay(P1421) <- "RNA"
P1421PscaPos <- subset(x=P1421, subset = Psca > 0)
P1421PscaNeg <- subset(x=P1421, subset = Psca == 0)
Idents(object = P1421PscaPos) <- "PscaPos"
Idents(object = P1421PscaNeg) <- "PscaNeg"
P1421PscaPos[["PscaExp"]] <- Idents(object = P1421PscaPos)
P1421PscaNeg[["PscaExp"]] <- Idents(object = P1421PscaNeg)
P1421Psca <- merge(x = P1421PscaPos, y = P1421PscaNeg)
Idents(object = P1421Psca) <- "PscaExp"
P1421$PscaExp <- Idents(object = P1421Psca)
#Total
Idents(object = P1421) <- "Osr1Exp"
table(Idents(P1421))
Idents(object = P1421) <- "Tacstd2Exp"
table(Idents(P1421))
Idents(object = P1421) <- "PscaExp"
table(Idents(P1421))
Idents(object = P1421) <- "Ly6aExp"
table(Idents(P1421))
#Epithelia
Idents(object = P1421) <- "seurat_clusters"
P1421.Epi <- subset(P1421, idents = c("8", "6", "5", "4"))
Idents(object = P1421.Epi) <- "Osr1Exp"
table(Idents(P1421.Epi))
Idents(object = P1421.Epi) <- "Tacstd2Exp"
table(Idents(P1421.Epi))
Idents(object = P1421.Epi) <- "PscaExp"
table(Idents(P1421.Epi))
Idents(object = P1421.Epi) <- "Ly6aExp"
table(Idents(P1421.Epi))

##P143
DefaultAssay(P1431) <- "RNA"
P1431Osr1Pos <- subset(x=P1431, subset = Osr1 > 0)
P1431Osr1Neg <- subset(x=P1431, subset = Osr1 == 0)
Idents(object = P1431Osr1Pos) <- "Osr1Pos"
Idents(object = P1431Osr1Neg) <- "Osr1Neg"
P1431Osr1Pos[["Osr1Exp"]] <- Idents(object = P1431Osr1Pos)
P1431Osr1Neg[["Osr1Exp"]] <- Idents(object = P1431Osr1Neg)
P1431Osr1 <- merge(x = P1431Osr1Pos, y = P1431Osr1Neg)
Idents(object = P1431Osr1) <- "Osr1Exp"
P1431$Osr1Exp <- Idents(object = P1431Osr1)
#Add Ly6a info
DefaultAssay(P1431) <- "RNA"
P1431Ly6aPos <- subset(x=P1431, subset = Ly6a > 0)
P1431Ly6aNeg <- subset(x=P1431, subset = Ly6a == 0)
Idents(object = P1431Ly6aPos) <- "Ly6aPos"
Idents(object = P1431Ly6aNeg) <- "Ly6aNeg"
P1431Ly6aPos[["Ly6aExp"]] <- Idents(object = P1431Ly6aPos)
P1431Ly6aNeg[["Ly6aExp"]] <- Idents(object = P1431Ly6aNeg)
P1431Ly6a <- merge(x = P1421Ly6aPos, y = P1431Ly6aNeg)
Idents(object = P1431Ly6a) <- "Ly6aExp"
P1431$Ly6aExp <- Idents(object = P1431Ly6a)
#Add Tacstd2 info
DefaultAssay(P1431) <- "RNA"
P1431Tacstd2Pos <- subset(x=P1431, subset = Tacstd2 > 0)
P1431Tacstd2Neg <- subset(x=P1431, subset = Tacstd2 == 0)
Idents(object = P1431Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = P1431Tacstd2Neg) <- "Tacstd2Neg"
P1431Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = P1431Tacstd2Pos)
P1431Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = P1431Tacstd2Neg)
P1431Tacstd2 <- merge(x = P1431Tacstd2Pos, y = P1431Tacstd2Neg)
Idents(object = P1431Tacstd2) <- "Tacstd2Exp"
P1431$Tacstd2Exp <- Idents(object = P1431Tacstd2)
#Add Psca info
DefaultAssay(P1431) <- "RNA"
P1431PscaPos <- subset(x=P1431, subset = Psca > 0)
P1431PscaNeg <- subset(x=P1431, subset = Psca == 0)
Idents(object = P1431PscaPos) <- "PscaPos"
Idents(object = P1431PscaNeg) <- "PscaNeg"
P1431PscaPos[["PscaExp"]] <- Idents(object = P1431PscaPos)
P1431PscaNeg[["PscaExp"]] <- Idents(object = P1431PscaNeg)
P1431Psca <- merge(x = P1431PscaPos, y = P1431PscaNeg)
Idents(object = P1431Psca) <- "PscaExp"
P1431$PscaExp <- Idents(object = P1431Psca)
#Total
Idents(object = P1431) <- "Osr1Exp"
table(Idents(P1431))
Idents(object = P1431) <- "Tacstd2Exp"
table(Idents(P1431))
Idents(object = P1431) <- "PscaExp"
table(Idents(P1431))
Idents(object = P1431) <- "Ly6aExp"
table(Idents(P1431))
#Epithelia
Idents(object = P1431) <- "seurat_clusters"
P1431.Epi <- subset(P1431, idents = c("0", "2", "4"))
Idents(object = P1431.Epi) <- "Osr1Exp"
table(Idents(P1431.Epi))
Idents(object = P1431.Epi) <- "Tacstd2Exp"
table(Idents(P1431.Epi))
Idents(object = P1431.Epi) <- "PscaExp"
table(Idents(P1431.Epi))
Idents(object = P1431.Epi) <- "Ly6aExp"
table(Idents(P1431.Epi))

##P351
#Add Osr1 info
DefaultAssay(P3511) <- "RNA"
P3511Osr1Pos <- subset(x=P3511, subset = Osr1 > 0)
P3511Osr1Neg <- subset(x=P3511, subset = Osr1 == 0)
Idents(object = P3511Osr1Pos) <- "Osr1Pos"
Idents(object = P3511Osr1Neg) <- "Osr1Neg"
P3511Osr1Pos[["Osr1Exp"]] <- Idents(object = P3511Osr1Pos)
P3511Osr1Neg[["Osr1Exp"]] <- Idents(object = P3511Osr1Neg)
P3511Osr1 <- merge(x = P3511Osr1Pos, y = P3511Osr1Neg)
Idents(object = P3511Osr1) <- "Osr1Exp"
P3511$Osr1Exp <- Idents(object = P3511Osr1)
#Add Tacstd2 info
DefaultAssay(P3511) <- "RNA"
P3511Tacstd2Pos <- subset(x=P3511, subset = Tacstd2 > 0)
P3511Tacstd2Neg <- subset(x=P3511, subset = Tacstd2 == 0)
Idents(object = P3511Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = P3511Tacstd2Neg) <- "Tacstd2Neg"
P3511Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = P3511Tacstd2Pos)
P3511Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = P3511Tacstd2Neg)
P3511Tacstd2 <- merge(x = P3511Tacstd2Pos, y = P3511Tacstd2Neg)
Idents(object = P3511Tacstd2) <- "Tacstd2Exp"
P3511$Tacstd2Exp <- Idents(object = P3511Tacstd2)
#Add Psca info
DefaultAssay(P3511) <- "RNA"
P3511PscaPos <- subset(x=P3511, subset = Psca > 0)
P3511PscaNeg <- subset(x=P3511, subset = Psca == 0)
Idents(object = P3511PscaPos) <- "PscaPos"
Idents(object = P3511PscaNeg) <- "PscaNeg"
P3511PscaPos[["PscaExp"]] <- Idents(object = P3511PscaPos)
P3511PscaNeg[["PscaExp"]] <- Idents(object = P3511PscaNeg)
P3511Psca <- merge(x = P3511PscaPos, y = P3511PscaNeg)
Idents(object = P3511Psca) <- "PscaExp"
P3511$PscaExp <- Idents(object = P3511Psca)
#Add Ly6a info
DefaultAssay(P3511) <- "RNA"
P3511Ly6aPos <- subset(x=P3511, subset = Ly6a > 0)
P3511Ly6aNeg <- subset(x=P3511, subset = Ly6a == 0)
Idents(object = P3511Ly6aPos) <- "Ly6aPos"
Idents(object = P3511Ly6aNeg) <- "Ly6aNeg"
P3511Ly6aPos[["Ly6aExp"]] <- Idents(object = P3511Ly6aPos)
P3511Ly6aNeg[["Ly6aExp"]] <- Idents(object = P3511Ly6aNeg)
P3511Ly6a <- merge(x = P3511Ly6aPos, y = P3511Ly6aNeg)
Idents(object = P3511Ly6a) <- "Ly6aExp"
P3511$Ly6aExp <- Idents(object = P3511Ly6a)
#Total
Idents(object = P3511) <- "Osr1Exp"
table(Idents(P3511))
Idents(object = P3511) <- "Tacstd2Exp"
table(Idents(P3511))
Idents(object = P3511) <- "PscaExp"
table(Idents(P3511))
Idents(object = P3511) <- "Ly6aExp"
table(Idents(P3511))
#Epitheia
Idents(object = P3511) <- "seurat_clusters"
P3511.Epi <- subset(P3511, idents = c("0", "1", "4", "5", "8"))
Idents(object = P3511.Epi) <- "Osr1Exp"
table(Idents(P3511.Epi))
Idents(object = P3511.Epi) <- "Tacstd2Exp"
table(Idents(P3511.Epi))
Idents(object = P3511.Epi) <- "PscaExp"
table(Idents(P3511.Epi))
Idents(object = P3511.Epi) <- "Ly6aExp"
table(Idents(P3511.Epi))

#P352
#Add Osr1 info
DefaultAssay(P3521) <- "RNA"
P3521Osr1Pos <- subset(x=P3521, subset = Osr1 > 0)
P3521Osr1Neg <- subset(x=P3521, subset = Osr1 == 0)
Idents(object = P3521Osr1Pos) <- "Osr1Pos"
Idents(object = P3521Osr1Neg) <- "Osr1Neg"
P3521Osr1Pos[["Osr1Exp"]] <- Idents(object = P3521Osr1Pos)
P3521Osr1Neg[["Osr1Exp"]] <- Idents(object = P3521Osr1Neg)
P3521Osr1 <- merge(x = P3521Osr1Pos, y = P3521Osr1Neg)
Idents(object = P3521Osr1) <- "Osr1Exp"
P3521$Osr1Exp <- Idents(object = P3521Osr1)
#Add Tacstd2 info
DefaultAssay(P3521) <- "RNA"
P3521Tacstd2Pos <- subset(x=P3521, subset = Tacstd2 > 0)
P3521Tacstd2Neg <- subset(x=P3521, subset = Tacstd2 == 0)
Idents(object = P3521Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = P3521Tacstd2Neg) <- "Tacstd2Neg"
P3521Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = P3521Tacstd2Pos)
P3521Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = P3521Tacstd2Neg)
P3521Tacstd2 <- merge(x = P3521Tacstd2Pos, y = P3521Tacstd2Neg)
Idents(object = P3521Tacstd2) <- "Tacstd2Exp"
P3521$Tacstd2Exp <- Idents(object = P3521Tacstd2)
#Add Psca info
DefaultAssay(P3521) <- "RNA"
P3521PscaPos <- subset(x=P3521, subset = Psca > 0)
P3521PscaNeg <- subset(x=P3521, subset = Psca == 0)
Idents(object = P3521PscaPos) <- "PscaPos"
Idents(object = P3521PscaNeg) <- "PscaNeg"
P3521PscaPos[["PscaExp"]] <- Idents(object = P3521PscaPos)
P3521PscaNeg[["PscaExp"]] <- Idents(object = P3521PscaNeg)
P3521Psca <- merge(x = P3521PscaPos, y = P3521PscaNeg)
Idents(object = P3521Psca) <- "PscaExp"
P3521$PscaExp <- Idents(object = P3521Psca)
#Add Ly6a info
DefaultAssay(P3521) <- "RNA"
P3521Ly6aPos <- subset(x=P3521, subset = Ly6a > 0)
P3521Ly6aNeg <- subset(x=P3521, subset = Ly6a == 0)
Idents(object = P3521Ly6aPos) <- "Ly6aPos"
Idents(object = P3521Ly6aNeg) <- "Ly6aNeg"
P3521Ly6aPos[["Ly6aExp"]] <- Idents(object = P3521Ly6aPos)
P3521Ly6aNeg[["Ly6aExp"]] <- Idents(object = P3521Ly6aNeg)
P3521Ly6a <- merge(x = P3521Ly6aPos, y = P3521Ly6aNeg)
Idents(object = P3521Ly6a) <- "Ly6aExp"
P3521$Ly6aExp <- Idents(object = P3521Ly6a)
#Total
Idents(object = P3521) <- "Osr1Exp"
table(Idents(P3521))
Idents(object = P3521) <- "Tacstd2Exp"
table(Idents(P3521))
Idents(object = P3521) <- "PscaExp"
table(Idents(P3521))
Idents(object = P3521) <- "Ly6aExp"
table(Idents(P3521))
#Epithelia
Idents(object = P3521) <- "seurat_clusters"
P3521.Epi <- subset(P3521, idents = c("1", "7", "13", "3", "11"))
Idents(object = P3521.Epi) <- "Osr1Exp"
table(Idents(P3521.Epi))
Idents(object = P3521.Epi) <- "Tacstd2Exp"
table(Idents(P3521.Epi))
Idents(object = P3521.Epi) <- "PscaExp"
table(Idents(P3521.Epi))
Idents(object = P3521.Epi) <- "Ly6aExp"
table(Idents(P3521.Epi))

##P353
#Add Osr1 info
DefaultAssay(P3531) <- "RNA"
P3531Osr1Pos <- subset(x=P3531, subset = Osr1 > 0)
P3531Osr1Neg <- subset(x=P3531, subset = Osr1 == 0)
Idents(object = P3531Osr1Pos) <- "Osr1Pos"
Idents(object = P3531Osr1Neg) <- "Osr1Neg"
P3531Osr1Pos[["Osr1Exp"]] <- Idents(object = P3531Osr1Pos)
P3531Osr1Neg[["Osr1Exp"]] <- Idents(object = P3531Osr1Neg)
P3531Osr1 <- merge(x = P3531Osr1Pos, y = P3531Osr1Neg)
Idents(object = P3531Osr1) <- "Osr1Exp"
P3531$Osr1Exp <- Idents(object = P3531Osr1)
#Add Tacstd2 info
DefaultAssay(P3531) <- "RNA"
P3531Tacstd2Pos <- subset(x=P3531, subset = Tacstd2 > 0)
P3531Tacstd2Neg <- subset(x=P3531, subset = Tacstd2 == 0)
Idents(object = P3531Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = P3531Tacstd2Neg) <- "Tacstd2Neg"
P3531Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = P3531Tacstd2Pos)
P3531Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = P3531Tacstd2Neg)
P3531Tacstd2 <- merge(x = P3531Tacstd2Pos, y = P3531Tacstd2Neg)
Idents(object = P3531Tacstd2) <- "Tacstd2Exp"
P3531$Tacstd2Exp <- Idents(object = P3531Tacstd2)
#Add Psca info
DefaultAssay(P3531) <- "RNA"
P3531PscaPos <- subset(x=P3531, subset = Psca > 0)
P3531PscaNeg <- subset(x=P3531, subset = Psca == 0)
Idents(object = P3531PscaPos) <- "PscaPos"
Idents(object = P3531PscaNeg) <- "PscaNeg"
P3531PscaPos[["PscaExp"]] <- Idents(object = P3531PscaPos)
P3531PscaNeg[["PscaExp"]] <- Idents(object = P3531PscaNeg)
P3531Psca <- merge(x = P3531PscaPos, y = P3531PscaNeg)
Idents(object = P3531Psca) <- "PscaExp"
P3531$PscaExp <- Idents(object = P3531Psca)
#Add Ly6a info
DefaultAssay(P3531) <- "RNA"
P3531Ly6aPos <- subset(x=P3531, subset = Ly6a > 0)
P3531Ly6aNeg <- subset(x=P3531, subset = Ly6a == 0)
Idents(object = P3531Ly6aPos) <- "Ly6aPos"
Idents(object = P3531Ly6aNeg) <- "Ly6aNeg"
P3531Ly6aPos[["Ly6aExp"]] <- Idents(object = P3531Ly6aPos)
P3531Ly6aNeg[["Ly6aExp"]] <- Idents(object = P3531Ly6aNeg)
P3531Ly6a <- merge(x = P3531Ly6aPos, y = P3531Ly6aNeg)
Idents(object = P3531Ly6a) <- "Ly6aExp"
P3531$Ly6aExp <- Idents(object = P3531Ly6a)
#Total
Idents(object = P3531) <- "Osr1Exp"
table(Idents(P3531))
Idents(object = P3531) <- "Tacstd2Exp"
table(Idents(P3531))
Idents(object = P3531) <- "PscaExp"
table(Idents(P3531))
Idents(object = P3531) <- "Ly6aExp"
table(Idents(P3531))
#Epithelia
Idents(object = P3531) <- "seurat_clusters"
P3531.Epi <- subset(P3531, idents = c("0", "1", "2", "5", "8", "14"))
Idents(object = P3531.Epi) <- "Osr1Exp"
table(Idents(P3531.Epi))
Idents(object = P3531.Epi) <- "Tacstd2Exp"
table(Idents(P3531.Epi))
Idents(object = P3531.Epi) <- "PscaExp"
table(Idents(P3531.Epi))
Idents(object = P3531.Epi) <- "Ly6aExp"
table(Idents(P3531.Epi))

####Gene-Gene Spearman correlation####
#Epithelia
GOI <- c('Osr1','Psca','Ly6a','Tacstd2', 'Itga6')  
GOI_index <- is.element(rownames(E1831),GOI)
Cell_index <- is.element(Idents(E1831),c('1'))
expr_GOI <- E1831@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- E1831@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
#Stroma
GOI <- c('Osr1','Psca','Ly6a','Tacstd2', 'Itga6')  
GOI_index <- is.element(rownames(E1831),GOI)
Cell_index <- is.element(Idents(E1831),c("0", "5", "6", "2", "4", "7"))
expr_GOI <- E1831@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- E1831@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)

####Loading data####
##PIN
PIN.data <- Read10X("//isi-dcnl/user_data/zjsun/BIC/AR Transgene Tumorigenesis/ARQ9-Osr1 SingleCell Seq/PCa_AR_transgene_analysis/PIN4_outs/filtered_feature_bc_matrix")
PIN <- CreateSeuratObject(counts = PIN.data,  min.cells = 3, min.features = 500, project = "PIN")
PIN <- NormalizeData(PIN)
##Tumor
Tumor.data <- Read10X("//isi-dcnl/user_data/zjsun/BIC/AR Transgene Tumorigenesis/ARQ9-Osr1 SingleCell Seq/PCa_AR_transgene_analysis/Tumor4_outs/filtered_feature_bc_matrix")
Tumor <- CreateSeuratObject(counts = Tumor.data,  min.cells = 3, min.features = 500, project = "Tumor")
Tumor <- NormalizeData(Tumor)

####Intitial processing, Filtering and Clustering####
##PIN
PIN[["percent.mt"]] <- PercentageFeatureSet(PIN, pattern = "^mt-")
VlnPlot(PIN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(PIN@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Pre-filteration")
hist(PIN@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Pre-filteration")
PIN2 <- subset(PIN, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
VlnPlot(PIN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(PIN2@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Post-filteration")
hist(PIN2@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Post-filteration")
PIN2 <- FindVariableFeatures(PIN2, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(PIN2)
PIN2 <- ScaleData(PIN2, verbose = FALSE)
PIN2 <- RunPCA(PIN2, npcs = 50, verbose = FALSE)
ElbowPlot(PIN2, ndims = 50)
PIN2 <- FindNeighbors(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- FindClusters(PIN2, resolution = 0.5)
PIN2 <- JackStraw(PIN2, num.replicate = 100)
PIN2 <- ScoreJackStraw(PIN2, dims = 1:20)
JackStrawPlot(PIN2, dims = 1:20)
PIN2 <- RunTSNE(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- RunUMAP(PIN2, reduction = "pca", dims = 1:20)
DimPlot(PIN2, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("#B71800")) 

##Tumor
Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^mt-")
VlnPlot(Tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Tumor@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Pre-filteration")
hist(Tumor@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Pre-filteration")
Tumor2 <- subset(Tumor, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
VlnPlot(Tumor2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Tumor2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Post-filteration")
hist(Tumor2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Post-filteration")
Tumor2 <- FindVariableFeatures(Tumor2, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(Tumor2)
Tumor2 <- ScaleData(Tumor2, verbose = FALSE)
Tumor2 <- RunPCA(Tumor2, npcs = 50, verbose = FALSE)
ElbowPlot(Tumor2, ndims = 50)
Tumor2 <- FindNeighbors(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- FindClusters(Tumor2, resolution = 0.5)
Tumor2 <- JackStraw(Tumor2, num.replicate = 100)
Tumor2 <- ScoreJackStraw(Tumor2, dims = 1:20)
JackStrawPlot(Tumor2, dims = 1:20)
Tumor2 <- RunTSNE(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- RunUMAP(Tumor2, reduction = "pca", dims = 1:20)
DimPlot(Tumor2, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("#009FB7")) 

####Merging Datasets####
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
# tSNE, UMAP and Clustering
PINvTumor.combined <- FindNeighbors(PINvTumor.combined, reduction = "pca", dims = 1:20)
PINvTumor.combined <- FindClusters(PINvTumor.combined, resolution = 0.5)
PINvTumor.combined <- RunTSNE(PINvTumor.combined, reduction = "pca", dims = 1:20)
PINvTumor.combined <- RunUMAP(PINvTumor.combined, reduction = "pca", dims = 1:20)
DimPlot(PINvTumor.combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("#B71800", "#009FB7")) 
#New labeling
new.cluster.ids <- c("LE", "BE", "LE", "LE", "LE", "LE", "LE", "Lym", "LE", "Fib", "LE", "LE", "Leu", "Lym", "Endo", "SM", "BE", "BE")
names(new.cluster.ids) <- levels(PINvTumor.combined)
PINvTumor.combined <- RenameIdents(PINvTumor.combined, new.cluster.ids)
PINvTumor.combined[["CellType"]] <- Idents(object = PINvTumor.combined)
DimPlot(PINvTumor.combined, reduction = "umap", pt.size = 0.3) 
DimPlot(PINvTumor.combined, split.by = "stim", reduction = "umap", pt.size = 0.3) 
#Blended expression UMAP plots
Idents(object = PINvTumor.combined) <- "stim"
PINonly <- subset(PINvTumor.combined, idents = c("PIN"))
DefaultAssay(PINonly) <- "RNA"
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "ARQ"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Epcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Vim"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Myh11"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Krt5"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Krt8"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
#Dotplot
Idents(object = PINvTumor.combined) <- "CellType"
PINvTumor.combined <- RenameIdents(object = PINvTumor.combined, 'BE' = "BE", 'LE' = "LE", 'Fib' = "Fib", 'SM' = "SM", 'Leu' = "Leu", 'Endo' = "Endo", 'Lym' = "Lym")
DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
tiff(file = "Combined DotPlot.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#### Reclustering Epithelial cells ####
Idents(object = PINvTumor.combined) <- "CellType"
PINvTumor.combined.Epi <- subset(PINvTumor.combined, idents = c("BE", "LE"))
Idents(object = PINvTumor.combined.Epi) <- "seurat_clusters"
DefaultAssay(PINvTumor.combined.Epi) <- "integrated"
#Run the standard workflow for visualization and clustering
PINvTumor.combined.Epi <- ScaleData(PINvTumor.combined.Epi, verbose = FALSE)
PINvTumor.combined.Epi <- RunPCA(PINvTumor.combined.Epi, npcs = 30, verbose = FALSE)
#tSNE, Umap and Clustering
PINvTumor.combined.Epi <- FindNeighbors(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi <- FindClusters(PINvTumor.combined.Epi, resolution = 0.5)
PINvTumor.combined.Epi <- RunTSNE(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi <- RunUMAP(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
DimPlot(PINvTumor.combined.Epi, reduction = "umap", pt.size = 0.3)
#Cell Cycle Regression
DefaultAssay(PINvTumor.combined.Epi) <- "RNA"
all.genes <- rownames(PINvTumor.combined.Epi)
PINvTumor.combined.Epi <- ScaleData(PINvTumor.combined.Epi, features = all.genes)
PINvTumor.combined.Epi <- CellCycleScoring(PINvTumor.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Idents(object = PINvTumor.combined.Epi) <- "Phase"
DimPlot(PINvTumor.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
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
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
#Rename
Idents(object = PINvTumor.combined.Epi1) <- "seurat_clusters"
PINvTumor.combined.Epi1 <- RenameIdents(object = PINvTumor.combined.Epi1, '4' = "BE1", '6' = "BE2", '0' = "LE1", '2' = "LE2", '1' = "LE3", '9' = "LE4", '3' = "LE5", '5' = "LE6", '8' = "LE7", '7' = "UrLE", '10' = "OE")
PINvTumor.combined.Epi1[["EpiCellType"]] <- Idents(object = PINvTumor.combined.Epi1)
#Umap
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey50"))
#Dotplots
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
DotPlot(PINvTumor.combined.Epi1, features = c("Ar", "ARQ", "EGFP", "Krt14", "Krt17", "Tpm2", "Krt16", "Col17a1", "Sult5a1", "Ltbp4", "Lbp", "Htra1", "Lamb3", "Wif1", "Defa20", "Defa22", "Msx2", "Syngr1", "Mt3", "Npl", "Gulo", "Mt4", "Nkain4", "Svs3a", "Wfdc15b", "Gpx3", "Ptgds", "Klk1b24", "Pcp4", "Smgc", "Snhg11", "Itln1", "Ceacam2", "Crip1", "Cd55", "Pmaip1", "Sbspon", "Kctd14", "Tgm4", "C1rb", "Gm5615", "Pnliprp1", "C1s2", "Oit1", "Cxcl17", "Timp4", "Aqp5", "Gjb2", "Wfdc2", "Msln", "Gsdmc2", "Atp10b", "Mgat4a", "Rgs1", "Srgn", "Tnfrsf9", "Cytip", "Coro1a"), cols = c("light grey", "red")) + RotatedAxis()
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
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(PINonlyEpi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
DefaultAssay(TumoronlyEpi) <- "RNA"
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TumoronlyEpi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

#### ARQ+vARQ- BE in Integrated ####
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("skyblue", "blue", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
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
Idents(object = PINvTumor.combined.Epi1) <- "ARQExp"
PINvTumor.combined.Epi2 <- subset(PINvTumor.combined.Epi1, idents = c("ARQPos", "ARQNeg"))
Idents(object = PINvTumor.combined.Epi2) <- "ARQExp.EpiCellType"
DimPlot(PINvTumor.combined.Epi2, reduction = "umap", pt.size = 0.3, cols = c("#3399FF", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "#E06666", "grey", "grey", "grey", "grey", "grey", "#E06666", "#3399FF", "grey", "grey", "grey", "grey", "grey")) + NoLegend()
#DEGs from hARtg+ vs hARtg- BE
Idents(object = PINvTumor.combined.Epi2) <- "CellType"
onlyBE2 <- subset(PINvTumor.combined.Epi2, idents = c("BE1", "BE2"))
onlyBE2 <- RenameIdents(object = onlyBE2, 'ARQNeg' = "ARQNeg", 'ARQPos' = "ARQPos")
onlyBE2.0.1.Markers <- FindMarkers(onlyBE2, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(onlyBE2.0.1.Markers, "onlyBE2.0.1.Markers.csv")
onlyBE2.0.Markers <- FindMarkers(onlyBE2, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(onlyBE2.0.Markers, "onlyBE2.0.Markers.csv")
DEG_onlyBE2 <- read.csv("onlyBE2.0.Markers.csv") 
DEG_onlyBE2_pvalue <- DEG_onlyBE2$p_val
DEG_onlyBE2_pvalue=as.numeric(DEG_onlyBE2_pvalue)
DEG_onlyBE2_BH = p.adjust(DEG_onlyBE2_pvalue, "BH")
write.csv(DEG_onlyBE2_BH, "DEG_onlyBE2_BH.csv")
#Heatmap
Idents(object = PINvTumor.combined.Epi2) <- "CellType"
onlyBE2 <- subset(PINvTumor.combined.Epi2, idents = c("BE1", "BE2"))
onlyBE2 <- RenameIdents(object = onlyBE2, 'ARQNeg' = "ARQNeg", 'ARQPos' = "ARQPos")
DefaultAssay(onlyBE2) <- "RNA"
onlyBE2 <- ScaleData(onlyBE2, features = rownames(onlyBE2))
onlyBE2.markers <- FindAllMarkers(onlyBE2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
onlyBE2Top50 <- onlyBE2.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "onlyBE2 Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(onlyBE2, features = c(onlyBE2Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
#Boxplots
boxdata = FetchData(onlyBE2, c("ARQExp", "ARQ", "Igf1r", "Jak2", "Mapk13"))
tail(boxdata,5)
ggplot(boxdata, aes(x=ARQExp, y=ARQ, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=ARQExp, y=Igf1r, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=ARQExp, y=Fos, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
ggplot(boxdata, aes(x=ARQExp, y=Jak2, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2))
ggplot(boxdata, aes(x=ARQExp, y=Mapk13, fill = ARQExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666")) + coord_cartesian(ylim=c(0, 2.25))
#Gene-Gene Spearman correlation
GOI1 <- c('ARQ','Igf1r', 'Jak2', 'Mapk13')  
GOI_index1 <- is.element(rownames(onlyBE2),GOI1)
Cell_index1 <- is.element(Idents(onlyBE2),c('BE1','BE2'))
expr_GOI1 <- onlyBE2@assays$RNA@data[GOI_index1,Cell_index1] 
expr_GOI1 <- onlyBE2@assays$RNA@counts[GOI_index1,Cell_index1]
ggcorr(t(expr_GOI1), method = c("pairwise", "spearman"), label = TRUE)

####LE1-2vLE5-7 in PINvTumor.combined.Epi####
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
DimPlot(PINvTumor.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "#E06666", "#E06666", "grey", "grey", "#3399FF", "#3399FF", "#3399FF", "grey", "grey"))
#New labeling
TumorvNormal.Epi1 <- subset(PINvTumor.combined.Epi1, idents = c("LE1", "LE2", "LE5", "LE6", "LE7"))
TumorvNormal.Epi1  <- RenameIdents(object = TumorvNormal.Epi1, 'LE5' = "Normal", 'LE6' = "Normal", 'LE7' = "Normal", 'LE1' = "Tumor", 'LE2' = "Tumor")
TumorvNormal.Epi1[["NormalTumor"]] <- Idents(object = TumorvNormal.Epi1)
#Heatmap
DefaultAssay(TumorvNormal.Epi1) <- "RNA"
TumorvNormal.Epi1 <- ScaleData(TumorvNormal.Epi1, features = rownames(TumorvNormal.Epi1))
TumorvNormal.Epi1.markers <- FindAllMarkers(TumorvNormal.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TumorvNormal.Epi1Top50 <- TumorvNormal.Epi1.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
DoHeatmap(TumorvNormal.Epi1, features = c(TumorvNormal.Epi1Top50$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
#Boxplots
boxdata = FetchData(TumorvNormal.Epi1, c("NormalTumor", "ARQ", "Tcf4", "Ccnd1", "Axin2", "Lgr5"))
tail(boxdata,6)
ggplot(boxdata, aes(x=NormalTumor, y=ARQ, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=NormalTumor, y=Tcf4, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=NormalTumor, y=Myc, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=NormalTumor, y=Ccnd1, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=NormalTumor, y=Axin2, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
ggplot(boxdata, aes(x=NormalTumor, y=Lgr5, fill = NormalTumor)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
#Featureplots
DefaultAssay(PINvTumor.combined.Epi1) <- "RNA"
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Tcf4"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Ccnd1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Axin2"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
FeaturePlot(PINvTumor.combined.Epi1, reduction = "umap", features = c("Lgr5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 3)
#Gene-Gene Spearman correlation
DefaultAssay(TumorvNormal.Epi1) <- "RNA"
Idents(object = TumorvNormal.Epi1) <- "NormalTumor"
GOI <- c('ARQ','Tcf4', 'Ccnd1', 'Axin2', 'Lgr5')  
GOI_index <- is.element(rownames(TumorvNormal.Epi1),GOI)
Cell_index <- is.element(Idents(TumorvNormal.Epi1),c('Normal','Tumor'))
expr_GOI <- TumorvNormal.Epi1@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- TumorvNormal.Epi1@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
#DEGs LE1-2 vs LE5-7
DefaultAssay(TumorvNormal.Epi1) <- "RNA"
TumorvNormal.Epi1.0.Markers <- FindMarkers(TumorvNormal.Epi1, ident.1 = "Tumor", ident.2 = "Normal", min.pct = 0., logfc.threshold = 0.)
write.csv(TumorvNormal.Epi1.0.Markers, "TumorvNormal.Epi1.0.Markers.csv")
DEG_Tumor <- read.csv("TumorvNormal.Epi1.0.Markers.csv") 
DEG_Tumor_pvalue <- DEG_Tumor$p_val
DEG_Tumor_pvalue=as.numeric(DEG_Tumor_pvalue)
DEG_Tumor_BH = p.adjust(DEG_Tumor_pvalue, "BH")
write.csv(DEG_Tumor_BH, "DEG_Tumor_BH.csv")


####Pseudotime of hARtg+ vs hARtg- Epithelial cells using Monocle2 ####
Idents(object = PINvTumor.combined.Epi1) <- "EpiCellType"
PINvTumor.combined.Epi2 <- subset(PINvTumor.combined.Epi1, idents = c("BE1", "BE2", "LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7"))
Idents(object = PINvTumor.combined.Epi2) <- "ARQExp"
ARQPosEpi <- subset(PINvTumor.combined.Epi2, idents = c("ARQPos"))
ARQNegEpi <- subset(PINvTumor.combined.Epi2, idents = c("ARQNeg"))
##ARQPosEpi
Idents(object = ARQPosEpi) <- "seurat_clusters"
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
plot_pc_variance_explained(EpiPseudo, return_all = F) 
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
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
plot_cell_trajectory(EpiPseudo, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) 
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt14", "Krt18"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
plot_cell_trajectory(EpiPseudo, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
#Expression level of Wnt downstreams in plot_genes_branched_pseudotime
EpiPseudo_genes1 <- row.names(subset(fData(EpiPseudo1), gene_short_name %in% c("Tcf4", "Axin2")))
plot_genes_branched_pseudotime(EpiPseudo1[EpiPseudo_genes1,],
                               branch_point = 1,
                               color_by = "EpiCellType",
                               ncol = 1)+ scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
EpiPseudo_genes1 <- row.names(subset(fData(EpiPseudo1), gene_short_name %in% c("Ccnd1", "Lgr5")))
plot_genes_branched_pseudotime(EpiPseudo1[EpiPseudo_genes1,],
                               branch_point = 1,
                               color_by = "EpiCellType",
                               ncol = 1)+ scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))

##ARQNegEpi
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
plot_pc_variance_explained(EpiPseudo1, return_all = F)
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
EpiPseudo1 <- orderCells(EpiPseudo1, root_state = GM_state(EpiPseudo1))-
#Visualization
plot_cell_trajectory(EpiPseudo1, color_by = "Pseudotime", show_branch_points = FALSE)
plot_cell_trajectory(EpiPseudo1, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) 
plot_cell_trajectory(EpiPseudo1, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt14", "Krt18"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
plot_cell_trajectory(EpiPseudo1, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))

####Pseudotime of hARtg+ vs hARtg- Epithelial cells using Slingshot####
##ARQPosEpi
ARQPosEpi1 <- ARQPosEpi
DefaultAssay(ARQPosEpi1) <- "integrated"
ARQPosEpi1 <- ScaleData(ARQPosEpi1, verbose = FALSE)
ARQPosEpi1 <- RunPCA(ARQPosEpi1, npcs = 30, verbose = FALSE)
ElbowPlot(ARQPosEpi1)
ARQPosEpi1 <- FindNeighbors(ARQPosEpi1, reduction = "pca", dims = 1:5)
ARQPosEpi1 <- FindClusters(ARQPosEpi1, resolution = 0.3)
ARQPosEpi1 <- RunTSNE(ARQPosEpi1, reduction = "pca", dims = 1:20)
ARQPosEpi1 <- RunUMAP(ARQPosEpi1, reduction = "pca", dims = 1:20)
#Slingshot UMAP
Idents(object = ARQPosEpi1) <- "EpiCellType"
ARQPosEpi1_sds <- slingshot(Embeddings(ARQPosEpi1, "umap"), clusterLabels = ARQPosEpi1$EpiCellType)
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
#Slingshot Pseudotime
ARQPosEpi1.1 <- as.SingleCellExperiment(ARQPosEpi1)
sce.sling2 <- slingshot(ARQPosEpi1.1, clusterLabels = ARQPosEpi1$EpiCellType, reducedDim='UMAP')
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)
reducedDim(sce.sling2, "UMAP") <- reducedDim(ARQPosEpi1.1, "UMAP")
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)
gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo), pt.size = 1)
embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), size=1.2)
}
gg

##ARQNegEpi
ARQNegEpi1 <- ARQNegEpi
DefaultAssay(ARQNegEpi1) <- "integrated"
ARQNegEpi1 <- ScaleData(ARQNegEpi1, verbose = FALSE)
ARQNegEpi1 <- RunPCA(ARQNegEpi1, npcs = 30, verbose = FALSE)
ElbowPlot(ARQNegEpi1)
ARQNegEpi1 <- FindNeighbors(ARQNegEpi1, reduction = "pca", dims = 1:5)
ARQNegEpi1 <- FindClusters(ARQNegEpi1, resolution = 0.3)
ARQNegEpi1 <- RunTSNE(ARQNegEpi1, reduction = "pca", dims = 1:5)
ARQNegEpi1 <- RunUMAP(ARQNegEpi1, reduction = "pca", dims = 1:5)
#Slingshot UMAP
Idents(object = ARQNegEpi1) <- "EpiCellType"
ARQNegEpi1_sds <- slingshot(Embeddings(ARQNegEpi1, "umap"), clusterLabels = ARQNegEpi1$seurat_clusters)
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
#Slingshot Pseudotime
ARQNegEpi1.1 <- as.SingleCellExperiment(ARQNegEpi1)
sce.sling3 <- slingshot(ARQNegEpi1.1, reducedDim='UMAP', clusterLabels = ARQNegEpi1$seurat_clusters)
pseudo.paths1 <- slingPseudotime(sce.sling3)
head(pseudo.paths1)
shared.pseudo1 <- rowMeans(pseudo.paths1, na.rm=TRUE)
gg <- plotUMAP(sce.sling3, colour_by=I(shared.pseudo1))
embedded1 <- embedCurves(sce.sling3, "UMAP")
embedded1 <- slingCurves(embedded1)
for (path in embedded1) {
  embedded1 <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded1, aes(x=UMAP_1, y=UMAP_2), size=1.2)
}
gg

####PseudotimeDE####
#Perform pseudotime inference on the original dataset
rd <- irlba::prcomp_irlba(t(logcounts(ARQPosEpi1)), scale. = FALSE)$x[, 1:2]
reducedDims(ARQPosEpi1) <- SimpleList(UMAP = rd)
colData(ARQPosEpi1)$cl <- 1
ARQPosEpi_ori <- slingshot(ARQPosEpi1, reducedDim = 'UMAP', clusterLabels = "cl")
ARQPosEpi_ori_tbl <- tibble(cell = colnames(ARQPosEpi1), pseudotime = rescale(colData(ARQPosEpi_ori)$slingPseudotime_1))
head(ARQPosEpi_ori_tbl)
#Perform pseudotime inference on subsamples
set.seed(123)
options(mc.cores = 2)
n = 1000
#Ganerate random subsamples
ARQPosEpi_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(ARQPosEpi1)[2]), size = 0.8*dim(ARQPosEpi1)[2], replace = FALSE)
})
ARQPosEpi_sub_tbl <- mclapply(ARQPosEpi_index, function(x, sce) {
  sce <- sce[, x]
  rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(UMAP = rd)
  fit <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = "cl")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  merge.tbl <- left_join(tbl, ARQPosEpi_ori_tbl, by = "cell")
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = ARQPosEpi1)
PseudotimeDE::plotUncertainty(ARQPosEpi_ori_tbl, ARQPosEpi_sub_tbl)
#Perform DE test
system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5", "Krt5", "Krt19"),
                                                 ori.tbl = ARQPosEpi_ori_tbl,
                                                 sub.tbl = ARQPosEpi_sub_tbl[1:100], 
                                                 mat = ARQPosEpi1, 
                                                 model = "nb",
                                                 mc.cores = 2))

print(res)
#Visualization
ARQPosEpi2 <- GetAssayData(object = ARQPosEpi, slot = "counts")
ARQPosEpi2 <- as.matrix(ARQPosEpi2)
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi2,
                        model.fit = res$gam.fit)