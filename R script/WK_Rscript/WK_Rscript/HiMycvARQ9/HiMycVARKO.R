##HiMycvARKO-Gli1##

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(RColorBrewer)
library(SeuratData)
library(BiocManager)
library(SeuratWrappers)
library(monocle3)
library(clustree)

setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/HiMycvARQ9")

#Load ARKO
ARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/ARKO-Gli1/Single_Cell_Seq/190219_scRNA/27348_ARKO1_count_EGFPmm10/outs/filtered_feature_bc_matrix/")
ARKOunfiltered <- CreateSeuratObject(counts = ARKOunfiltered.data,  min.cells = 3, min.features = 500, project = "ARKOunfiltered")
ARKOunfiltered <- NormalizeData(ARKOunfiltered)
ARKOunfiltered[["percent.mt"]] <- PercentageFeatureSet(ARKOunfiltered, pattern = "^mt-")
ARKO <- subset(ARKOunfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)
VlnPlot(ARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
hist(ARKO@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARKO post filteration")
hist(ARKO@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARKO post filteration")
ARKO <- FindVariableFeatures(ARKO, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(ARKO)
ARKO <- ScaleData(ARKO, features = all.genes)
ARKO <- RunPCA(ARKO, features = VariableFeatures(object = ARKO))
ElbowPlot(ARKO, ndims = 50)
ARKO <- FindNeighbors(ARKO, reduction = "pca", dims = 1:25)
ARKO <- FindClusters(ARKO, resolution = 0.5)
ARKO <- RunUMAP(ARKO, reduction = "pca", dims = 1:25)
DimPlot(ARKO, reduction = "umap", pt.size = 1, label = TRUE)

#Load HiMyc
##HiMYC
Ctrlunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/2nd Run 6-4-2020/Sun Lab Analysis/Project_COHP_36595_1_XCTC1_count/outs/filtered_feature_bc_matrix")
Ctrlunfiltered <- CreateSeuratObject(counts = Ctrlunfiltered.data,  min.cells = 3, min.features = 500, project = "Ctrlunfiltered")
Ctrlunfiltered <- NormalizeData(Ctrlunfiltered)
Ctrlunfiltered[["percent.mt"]] <- PercentageFeatureSet(Ctrlunfiltered, pattern = "^mt-")
VlnPlot(Ctrlunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrlunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Pre-filteration")
hist(Ctrlunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Pre-filteration")
Ctrl <- subset(Ctrlunfiltered, subset = nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 15)
VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Post-filteration")
hist(Ctrl@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Post-filteration")
Ctrl <- FindVariableFeatures(Ctrl, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(Ctrl)
Ctrl <- ScaleData(Ctrl, verbose = FALSE)
Ctrl <- RunPCA(Ctrl, npcs = 30, verbose = FALSE)
ElbowPlot(Ctrl, ndims = 50)
Ctrl <- FindNeighbors(Ctrl, reduction = "pca", dims = 1:20)
Ctrl <- FindClusters(Ctrl, resolution = 0.5)
Ctrl <- RunUMAP(Ctrl, reduction = "pca", dims = 1:20)
DimPlot(Ctrl, reduction = "umap", pt.size = 0.3)

####Integration####
#Set Current idents
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
ARKO[["orig.clusters"]] <- Idents(object = ARKO)
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
ARKO$stim <- "ARKO"
CtrlvARKO.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, ARKO), dims = 1:20)
CtrlvARKO.combined <- IntegrateData(anchorset = CtrlvARKO.anchors, dims = 1:20)
DefaultAssay(CtrlvARKO.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, verbose = FALSE)
CtrlvARKO.combined <- RunPCA(CtrlvARKO.combined, npcs = 30, verbose = FALSE)
CtrlvARKO.combined <- FindNeighbors(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
CtrlvARKO.combined <- FindClusters(CtrlvARKO.combined, resolution = 0.5)
CtrlvARKO.combined <- RunUMAP(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
Idents(object = CtrlvARKO.combined) <- "stim"
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 

#Cell type identification
DefaultAssay(CtrlvARKO.combined)<-"RNA"
tiff(file = "CtrlvARKO.combined celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                 "Fbln1", "Myh11", "Pecam1",
                                                                 "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

Idents(object = CtrlvARKO.combined) <- "seurat_clusters"
tiff(file = "CtrlvARKO.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = T) 
dev.off()

####Subclustering FBSM####
Idents(object = CtrlvARKO.combined) <- "seurat_clusters"
CtrlvARKO.combined.FBSM <- subset(CtrlvARKO.combined, idents = c("1", "2"))
CtrlvARKO.combined.FBSM <- RunUMAP(CtrlvARKO.combined.FBSM, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(CtrlvARKO.combined.FBSM) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, verbose = FALSE)
CtrlvARKO.combined.FBSM <- RunPCA(CtrlvARKO.combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.FBSM, ndims = 50)
#Umap and Clustering
DefaultAssay(CtrlvARKO.combined.FBSM) <- "integrated"
CtrlvARKO.combined.FBSM <- FindNeighbors(CtrlvARKO.combined.FBSM, reduction = "pca", dims = 1:10)
CtrlvARKO.combined.FBSM <- FindClusters(CtrlvARKO.combined.FBSM, resolution = 2)
CtrlvARKO.combined.FBSM <- RunUMAP(CtrlvARKO.combined.FBSM, reduction = "pca", dims = 1:10)

Idents(object = CtrlvARKO.combined.FBSM) <- "seurat_clusters"
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(CtrlvARKO.combined.FBSM)<-"RNA"
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#CtrlvARKO.combined.FBSM1
Idents(object = CtrlvARKO.combined.FBSM) <- "seurat_clusters"
CtrlvARKO.combined.FBSM1 <- subset(CtrlvARKO.combined.FBSM, idents = c("0", "1", "2", "4", "5", "6", "7", 
                                                                       "8", "9", "10", "12", "13", "14",
                                                                       "15", "16", "18", "19", "20",
                                                                       "21", "22"))
CtrlvARKO.combined.FBSM1 <- RunUMAP(CtrlvARKO.combined.FBSM1, reduction = "pca", dims = 1:9)
DimPlot(CtrlvARKO.combined.FBSM1, reduction = "umap", label = TRUE, split.by = "stim")

DefaultAssay(CtrlvARKO.combined.FBSM1)<-"RNA"
FeaturePlot(CtrlvARKO.combined.FBSM1, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM1, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM1, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM1, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM1, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#higher resolution
CtrlvARKO.combined.FBSM2 <- CtrlvARKO.combined.FBSM1
DefaultAssay(CtrlvARKO.combined.FBSM2) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.FBSM2 <- ScaleData(CtrlvARKO.combined.FBSM2, verbose = FALSE)
CtrlvARKO.combined.FBSM2 <- RunPCA(CtrlvARKO.combined.FBSM2, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.FBSM2, ndims = 50)
#Umap and Clustering
DefaultAssay(CtrlvARKO.combined.FBSM2) <- "integrated"
CtrlvARKO.combined.FBSM2 <- FindNeighbors(CtrlvARKO.combined.FBSM2, reduction = "pca", dims = 1:5)
CtrlvARKO.combined.FBSM2 <- FindClusters(CtrlvARKO.combined.FBSM2, resolution = 3)
CtrlvARKO.combined.FBSM2 <- RunUMAP(CtrlvARKO.combined.FBSM2, reduction = "pca", dims = 1:5)

DimPlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(CtrlvARKO.combined.FBSM2)<-"RNA"
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q99")
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")


#CtrlvARKO.combined.FBSM3
Idents(object = CtrlvARKO.combined.FBSM2) <- "seurat_clusters"
CtrlvARKO.combined.FBSM3 <- subset(CtrlvARKO.combined.FBSM2, idents = c("0", "1", "2", "3", "4", "5", "6", "7", 
                                                                       "8", "9", "10", "11", "12", "13", "14",
                                                                       "15", "16", "17", "18", "19", 
                                                                       "21", "22", "24", "25", "26", "27", "28", "29"))
CtrlvARKO.combined.FBSM3 <- RunUMAP(CtrlvARKO.combined.FBSM3, reduction = "pca", dims = 1:5)
DimPlot(CtrlvARKO.combined.FBSM3, reduction = "umap", label = TRUE, split.by = "stim")

DefaultAssay(CtrlvARKO.combined.FBSM3)<-"RNA"
FeaturePlot(CtrlvARKO.combined.FBSM3, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM3, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

Idents(object = CtrlvARKO.combined.FBSM3) <- "seurat_clusters"
CtrlvARKO.combined.FBSM4 <- subset(CtrlvARKO.combined.FBSM3, idents = c("1", "2", "3", "4", "5", "6", "7", 
                                                                        "8", "9", "10", "11", "12", "13", "14",
                                                                        "15", "16", "17", "18", "19", 
                                                                        "21", "22", "24", "25", "26", "27", "28", "29"))
DimPlot(CtrlvARKO.combined.FBSM4, reduction = "umap", label = TRUE, split.by = "stim")

DefaultAssay(CtrlvARKO.combined.FBSM4)<-"RNA"
FeaturePlot(CtrlvARKO.combined.FBSM4, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.FBSM4, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Rename
Idents(object = CtrlvARKO.combined.FBSM4) <- "seurat_clusters"
CtrlvARKO.combined.FBSM4 <- RenameIdents(object = CtrlvARKO.combined.FBSM4, 
                                         '7' = "Telocyte", '12' = "Telocyte",
                                         '1' = "FB",'2' = "FB",'3' = "FB",'4' = "FB",
                                         '5' = "FB",'6' = "FB",'7' = "FB",'8' = "FB",'9' = "FB",
                                         '10' = "FB",'11' = "FB",'12' = "FB",'13' = "FB",'14' = "FB",'15' = "FB",
                                         '16' = "FB",'17' = "FB",'18' = "FB",'19' = "FB",'21' = "FB",'22' = "FB",
                                         '24' = "FB", '25' = "FB",'26' = "FB",'27' = "FB",'28' = "FB",
                                         '29' = "FB")
CtrlvARKO.combined.FBSM4[["FBSMCellTypes1"]] <- Idents(object = CtrlvARKO.combined.FBSM4)
DimPlot(CtrlvARKO.combined.FBSM4, split.by = "stim", reduction = "umap", cols = c("red","blue","darkorange","skyblue1", "bisque3","palevioletred3",
                                                               "green4"))
DefaultAssay(CtrlvARKO.combined.FBSM4) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM4) <- "FBSMCellTypes1"
tiff(file = "CtrlvARKO.combined.FBSM4 vlnplot EGFP split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM4, features = "EGFP", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM4 vlnplot Foxl1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM4, features = "Foxl1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM4 vlnplot Pdgfra split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM4, features = "Pdgfra", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM4 vlnplot Gli1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM4, features = "Gli1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM4 vlnplot CD34 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM4, features = "Cd34", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()

####SCT####

setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/HiMycvARQ9")

#Load ARKO
ARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/ARKO-Gli1/Single_Cell_Seq/190219_scRNA/27348_ARKO1_count_EGFPmm10/outs/filtered_feature_bc_matrix/")
ARKOunfiltered <- CreateSeuratObject(counts = ARKOunfiltered.data,  min.cells = 3, min.features = 500, project = "ARKOunfiltered")
ARKOunfiltered[["percent.mt"]] <- PercentageFeatureSet(ARKOunfiltered, pattern = "^mt-")
ARKO <- subset(ARKOunfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)
VlnPlot(ARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
hist(ARKO@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARKO post filteration")
hist(ARKO@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARKO post filteration")

#Apply sctransform normalization
ARKO <- SCTransform(ARKO, vars.to.regress = "percent.mt", verbose = FALSE)

#Perform dimensionality reduction by PCA and UMAP embedding
ARKO <- RunPCA(ARKO, verbose = FALSE)
ElbowPlot(ARKO, ndims = 50)
ARKO <- FindNeighbors(ARKO, dims = 1:25, verbose = FALSE)
ARKO <- FindClusters(ARKO, resolution = 0.5, verbose = FALSE)
ARKO <- RunUMAP(ARKO, dims = 1:25, verbose = FALSE)
DimPlot(ARKO, reduction = "umap", pt.size = 0.3, label = TRUE)

#Load HiMyc
##HiMYC
Ctrlunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/2nd Run 6-4-2020/Sun Lab Analysis/Project_COHP_36595_1_XCTC1_count/outs/filtered_feature_bc_matrix")
Ctrlunfiltered <- CreateSeuratObject(counts = Ctrlunfiltered.data,  min.cells = 3, min.features = 500, project = "Ctrlunfiltered")
Ctrlunfiltered[["percent.mt"]] <- PercentageFeatureSet(Ctrlunfiltered, pattern = "^mt-")
VlnPlot(Ctrlunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrlunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Pre-filteration")
hist(Ctrlunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Pre-filteration")
Ctrl <- subset(Ctrlunfiltered, subset = nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 15)
VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Post-filteration")
hist(Ctrl@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Post-filteration")
#Apply sctransform normalization
Ctrl <- SCTransform(Ctrl, vars.to.regress = "percent.mt", verbose = FALSE)

#Perform dimensionality reduction by PCA and UMAP embedding
Ctrl <- RunPCA(Ctrl, verbose = FALSE)
ElbowPlot(Ctrl, ndims = 50)
Ctrl <- FindNeighbors(Ctrl, dims = 1:25, verbose = FALSE)
Ctrl <- FindClusters(Ctrl, resolution = 0.5, verbose = FALSE)
Ctrl <- RunUMAP(Ctrl, dims = 1:25, verbose = FALSE)
DimPlot(Ctrl, reduction = "umap", pt.size = 0.3, label = TRUE)

#Stash old idents
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
ARKO[["orig.clusters"]] <- Idents(object = ARKO)

#Set Current idents
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"

Ctrl$stim <- "Ctrl"
ARKO$stim <- "ARKO"

#Perform integration using pearson residuals
ifnb.list <- list(Ctrl = Ctrl, ARKO = ARKO)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
CtrlvARKO.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                          anchor.features = features)
CtrlvARKO.combined.sct <- IntegrateData(anchorset = CtrlvARKO.anchors, normalization.method = "SCT")

CtrlvARKO.combined.sct <- RunPCA(CtrlvARKO.combined.sct, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.sct, ndims = 50)
CtrlvARKO.combined.sct <- RunUMAP(CtrlvARKO.combined.sct, reduction = "pca", dims = 1:25, verbose = FALSE)
CtrlvARKO.combined.sct <- FindNeighbors(CtrlvARKO.combined.sct, reduction = "pca", dims = 1:25)
CtrlvARKO.combined.sct <- FindClusters(CtrlvARKO.combined.sct, resolution = 0.5)
DimPlot(CtrlvARKO.combined.sct, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = CtrlvARKO.combined.sct) <- "stim"
CtrlvARKO.combined.sct <- RenameIdents(object = CtrlvARKO.combined.sct, 'Ctrl' = "Ctrl", 'ARKO' = "ARKO")
CtrlvARKO.combined.sct[["stim"]] <- Idents(object = CtrlvARKO.combined.sct)

tiff(file = "CtrlvARKO.combined.sct UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.sct, reduction = "umap", pt.size = 0.3, cols = c("darkblue", "grey"))
dev.off()

#Cell type identification
DefaultAssay(CtrlvARKO.combined.sct)<-"SCT"
tiff(file = "CtrlvARKO.combined.sct celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.sct, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                               "Fbln1", "Myh11", "Pecam1",
                                                               "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

Idents(object = CtrlvARKO.combined.sct) <- "seurat_clusters"
tiff(file = "CtrlvARKO.combined.sct UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.sct, reduction = "umap", pt.size = 0.3, label = T) 
dev.off()

####Subclustering FBSM####
Idents(object = CtrlvARKO.combined.sct) <- "seurat_clusters"
CtrlvARKO.combined.sct.FB <- subset(CtrlvARKO.combined.sct, idents = c("1", "2", "7"))
CtrlvARKO.combined.sct.FB <- RunUMAP(CtrlvARKO.combined.sct.FB, reduction = "pca", dims = 1:10)
DimPlot(CtrlvARKO.combined.sct.FB, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(CtrlvARKO.combined.sct.FB) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.sct.FB <- ScaleData(CtrlvARKO.combined.sct.FB, verbose = FALSE)
CtrlvARKO.combined.sct.FB <- RunPCA(CtrlvARKO.combined.sct.FB, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.sct.FB, ndims = 50)
#Umap and Clustering
DefaultAssay(CtrlvARKO.combined.sct.FB) <- "integrated"
CtrlvARKO.combined.sct.FB <- FindNeighbors(CtrlvARKO.combined.sct.FB, reduction = "pca", dims = 1:10)
CtrlvARKO.combined.sct.FB <- FindClusters(CtrlvARKO.combined.sct.FB, resolution = 2.5)
CtrlvARKO.combined.sct.FB <- RunUMAP(CtrlvARKO.combined.sct.FB, reduction = "pca", dims = 1:10)

Idents(object = CtrlvARKO.combined.sct.FB) <- "seurat_clusters"
DimPlot(CtrlvARKO.combined.sct.FB, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(CtrlvARKO.combined.sct.FB)<-"SCT"
FeaturePlot(CtrlvARKO.combined.sct.FB, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.sct.FB, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.sct.FB, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.sct.FB, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvARKO.combined.sct.FB, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#CtrlvARKO.combined.sct.FB1
Idents(object = CtrlvARKO.combined.sct.FB) <- "seurat_clusters"
CtrlvARKO.combined.sct.FB1 <- subset(CtrlvARKO.combined.sct.FB, idents = c("0", "1", "2", "3", "4", "5", "6", "7", 
                                                                        "8", "9", "10", "11", "12", "13", "14", "15",
                                                                        "16", 
                                                                          "18", "19", "20",
                                                                        "21", "22",  "23", "24", "25", "26"))
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap", label = TRUE, split.by = "stim")

#Rename
Idents(object = CtrlvARKO.combined.sct.FB1) <- "seurat_clusters"
CtrlvARKO.combined.sct.FB1 <- RenameIdents(object = CtrlvARKO.combined.sct.FB1, 
                                         '4' = "Telocyte", 
                                         '26' = "FB1",'6' = "FB1",'13' = "FB1",'3' = "FB1",'5' = "FB1",
                                         '24' = "FB2",'20' = "FB2",'23' = "FB2",'8' = "FB2",
                                         '0' = "FB3",'2' = "FB3",'7' = "FB3",'18' = "FB3",
                                         '10' = "FB4",'25' = "FB4",
                                         '21' = "FB4",'1' = "FB4",'11' = "FB4",
                                         '14' = "FB5",'19' = "FB5",'12' = "FB5",'9' = "FB5"
                                         ,'16' = "FB6",'15' = "FB6",'22' = "FB6")
CtrlvARKO.combined.sct.FB1[["FBSMCellTypes1"]] <- Idents(object = CtrlvARKO.combined.sct.FB1)
DimPlot(CtrlvARKO.combined.sct.FB1, split.by = "stim", reduction = "umap", cols = c("red","skyblue1","darkorange","green4", "bisque3","palevioletred3",
                                                                                    "steelblue3","grey"))

Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBSMCellTypes1"
tiff(file = "CtrlvARKO.combined.sct.FB1 FBSMCellTypes1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap", cols = c("red","skyblue1","darkorange","green4", "bisque3","palevioletred3",
                                                                 "steelblue3","grey"))
dev.off()
tiff(file = "CtrlvARKO.combined.sct.FB1 FBSMCellTypes1 split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap", split.by="stim", cols = c("red","skyblue1","darkorange","green4", "bisque3","palevioletred3",
                                                                                  "steelblue3","grey"))
dev.off()

Idents(object = CtrlvARKO.combined.sct.FB1) <- "stim"
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap")
Ctrl.combined.sct.FB1 <- subset(CtrlvARKO.combined.sct.FB1, idents = c("Ctrl"))
ARKO.combined.sct.FB1 <- subset(CtrlvARKO.combined.sct.FB1, idents = c("ARKO"))


DefaultAssay(Ctrl.combined.sct.FB1) <- "SCT"
tiff(file = "Ctrl.combined.sct.FB1 PDGFRa.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.sct.FB1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "Ctrl.combined.sct.FB1 Cd34.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.sct.FB1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "Ctrl.combined.sct.FB1 Foxl1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.sct.FB1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q50")
dev.off()
tiff(file = "Ctrl.combined.sct.FB1 Ar.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.sct.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "Ctrl.combined.sct.FB1 Gli1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.sct.FB1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q70")
dev.off()
tiff(file = "Ctrl.combined.sct.FB1 EGFP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.sct.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()

DefaultAssay(ARKO.combined.sct.FB1) <- "SCT"
tiff(file = "ARKO.combined.sct.FB1 PDGFRa.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Cd34.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Foxl1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Ar.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Gli1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 EGFP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99", min.cutoff = 2)
dev.off()

#split
Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBSMCellTypes1"
CtrlvARKO.combined.sct.FB1 <- RenameIdents(object = CtrlvARKO.combined.sct.FB1, 
                                         'Telocyte' = "Telocyte", 'FB1' = "FB",'FB2' = "FB",'FB3' = "FB",
                                         'FB4' = "FB",'FB5' = "FB", "FB6" = "FB")
CtrlvARKO.combined.sct.FB1[["FBvTC"]] <- Idents(object = CtrlvARKO.combined.sct.FB1)
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap")

Idents(object = CtrlvARKO.combined.sct.FB1) <- "stim"
CtrlvARKO.combined.sct.FB1 <- RenameIdents(object = CtrlvARKO.combined.sct.FB1, 
                                           'Ctrl' = "Ctrl", 'ARKO' = "ARKO")
CtrlvARKO.combined.sct.FB1[["stim"]] <- Idents(object = CtrlvARKO.combined.sct.FB1)
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap")

DefaultAssay(CtrlvARKO.combined.sct.FB1) <- "SCT"
Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBvTC"
tiff(file = "CtrlvARKO.combined.sct.FB1 vlnplot EGFP split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.sct.FB1, features = "EGFP", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.sct.FB1 vlnplot Foxl1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.sct.FB1, features = "Foxl1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.sct.FB1 vlnplot Pdgfra split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.sct.FB1, features = "Pdgfra", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.sct.FB1 vlnplot Gli1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.sct.FB1, features = "Gli1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.sct.FB1 vlnplot CD34 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.sct.FB1, features = "Cd34", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.sct.FB1 vlnplot Ar split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.sct.FB1, features = "Ar", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBvTC"
CtrlvARKO.combined.sct.FB1$FBvTC.stim <- paste(Idents(CtrlvARKO.combined.sct.FB1), CtrlvARKO.combined.sct.FB1$stim, sep = "_")
Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBvTC.stim"
table(Idents(CtrlvARKO.combined.sct.FB1))


Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBvTC.stim"
CtrlvARKO.combined.sct.FB1 <- RenameIdents(object = CtrlvARKO.combined.sct.FB1, 
                                           'Telocyte_Ctrl' = "Telocyte_Ctrl", 'FB_Ctrl' = "FB_Ctrl",
                                           'Telocyte_ARKO' = "Telocyte_ARKO",'FB_ARKO' = "FB_ARKO")
CtrlvARKO.combined.sct.FB1[["FBvTC.stim"]] <- Idents(object = CtrlvARKO.combined.sct.FB1)
DimPlot(CtrlvARKO.combined.sct.FB1, reduction = "umap")


#Dotplot
Idents(object = CtrlvARKO.combined.sct.FB1) <- "FBvTC.stim"
DefaultAssay(CtrlvARKO.combined.sct.FB1) <- "SCT"



