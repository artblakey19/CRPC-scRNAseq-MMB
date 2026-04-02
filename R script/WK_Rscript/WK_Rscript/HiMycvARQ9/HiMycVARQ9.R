##HiMycvARQ9-Gli1##

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

####Load Data ARQ9 2M####

setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/HiMycvARQ9")

#Setup the seurat object

TwoM.unfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/count_new_3ref/count_Project_COHP_46245_1_X3SC4/outs/filtered_feature_bc_matrix")
TwoM.unfiltered <- CreateSeuratObject(counts = TwoM.unfiltered.data,  min.cells = 3, min.features = 200, project = "ARQ9-2M")

#Initial processing & filtering

TwoM.unfiltered[["percent.mt"]] <- PercentageFeatureSet(TwoM.unfiltered, pattern = "^mt-")
tiff(file = "TwoM Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(TwoM.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "TwoM Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "TwoM Pre-filteration")
dev.off()
tiff(file = "TwoM Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "TwoM Pre-filteration")
dev.off()

plot1 <- FeatureScatter(TwoM.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TwoM.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

TwoM <- subset(TwoM.unfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & percent.mt < 15)
tiff(file = "TwoM Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(TwoM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "TwoM Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "TwoM Post-filteration")
dev.off()
tiff(file = "TwoM Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "TwoM Post-filteration")
dev.off()

plot1 <- FeatureScatter(TwoM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TwoM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
TwoM <- NormalizeData(TwoM)

#Identification of highly variable features
TwoM <- FindVariableFeatures(TwoM, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(TwoM)

#Scaling the data
all.genes <- rownames(TwoM)
TwoM <- ScaleData(TwoM, features = all.genes)

#Perform linear dimensional reduction
TwoM <- RunPCA(TwoM, features = VariableFeatures(object = TwoM))

#Determine the dimensionality of the dataset
ElbowPlot(TwoM, ndims = 50)

#Cluster the cells
TwoM <- FindNeighbors(TwoM, dims = 1:30)
TwoM <- FindClusters(TwoM, resolution = 0.8)

#Run non-linear dimensional reduction
TwoM <- RunUMAP(TwoM, dims = 1:30)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoM) <- "stim"
tiff(file = "TwoM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3, cols = c("darkblue")) 
dev.off()

table(Idents(TwoM))

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
TwoM[["orig.clusters"]] <- Idents(object = TwoM)
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = TwoM) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
TwoM$stim <- "ARQ9"
CtrlvTwoM.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, TwoM), dims = 1:20)
CtrlvTwoM.combined <- IntegrateData(anchorset = CtrlvTwoM.anchors, dims = 1:20)
DefaultAssay(CtrlvTwoM.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvTwoM.combined <- ScaleData(CtrlvTwoM.combined, verbose = FALSE)
CtrlvTwoM.combined <- RunPCA(CtrlvTwoM.combined, npcs = 30, verbose = FALSE)
CtrlvTwoM.combined <- FindNeighbors(CtrlvTwoM.combined, reduction = "pca", dims = 1:20)
CtrlvTwoM.combined <- FindClusters(CtrlvTwoM.combined, resolution = 0.5)
CtrlvTwoM.combined <- RunUMAP(CtrlvTwoM.combined, reduction = "pca", dims = 1:20)
Idents(object = CtrlvTwoM.combined) <- "stim"
DimPlot(CtrlvTwoM.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 

#Cell type identification
DefaultAssay(CtrlvTwoM.combined)<-"RNA"
tiff(file = "CtrlvTwoM.combined celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvTwoM.combined, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                      "Fbln1", "Myh11", "Pecam1",
                                                                      "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

Idents(object = CtrlvTwoM.combined) <- "seurat_clusters"
tiff(file = "CtrlvTwoM.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvTwoM.combined, reduction = "umap", pt.size = 0.3, label = T) 
dev.off()

#New Labeling
CtrlvTwoM.combined <- RenameIdents(object = CtrlvTwoM.combined, '0' = "BE",  '1' = "BE",  '17' = "BE",  '7' = "BE", 
                                   '4' = "LE",  '13' = "LE", '12' = "LE",'3' = "LE",'5' = "LE",
                                   '20' = "SV",'19' = "SV",'10' = "SV",'18' = "Glia",
                                   '16' = "SM",'14' = "Pericyte",'2' = "FB",'6' = "FB",'11' = "VE",
                                   '15' = "Immune",'9' = "Immune",'8' = "Immune")
CtrlvTwoM.combined[["celltype"]] <- Idents(object = CtrlvTwoM.combined)
DimPlot(CtrlvTwoM.combined, reduction = "umap", pt.size = 0.3, label = T) 

####Subset Stro####
Idents(object = CtrlvTwoM.combined) <- "celltype"
CtrlvTwoM.combined.Stro <- subset(CtrlvTwoM.combined, idents = c("FB", "SM", "Pericyte", "Immune", "VE", "Glia"))
DefaultAssay(CtrlvTwoM.combined.Stro) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvTwoM.combined.Stro <- ScaleData(CtrlvTwoM.combined.Stro, verbose = FALSE)
CtrlvTwoM.combined.Stro <- RunPCA(CtrlvTwoM.combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvTwoM.combined.Stro, ndims = 50)
#Umap and Clustering
CtrlvTwoM.combined.Stro <- FindNeighbors(CtrlvTwoM.combined.Stro, reduction = "pca", dims = 1:15)
CtrlvTwoM.combined.Stro <- FindClusters(CtrlvTwoM.combined.Stro, resolution = 0.5)
CtrlvTwoM.combined.Stro <- RunUMAP(CtrlvTwoM.combined.Stro, reduction = "pca", dims = 1:15)
DimPlot(CtrlvTwoM.combined.Stro, reduction = "umap", label = TRUE)

#Cell type identification
DefaultAssay(CtrlvTwoM.combined.Stro)<-"RNA"
tiff(file = "CtrlvTwoM.combined.Stro celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvTwoM.combined.Stro, reduction = "umap", features = c("Krt5",  "Krt8", "Cd24a", "Plp1",
                                                                 "Fbln1", "Myh11", "Pecam1", "Rgs5",
                                                                 "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####Subclustering FBSM####
Idents(object = CtrlvTwoM.combined.Stro) <- "seurat_clusters"
CtrlvTwoM.combined.FBSM <- subset(CtrlvTwoM.combined.Stro, idents = c("0", "5", "2", "1", "9"))
CtrlvTwoM.combined.FBSM <- RunUMAP(CtrlvTwoM.combined.FBSM, reduction = "pca", dims = 1:15)
DimPlot(CtrlvTwoM.combined.FBSM, reduction = "umap", label = TRUE)

#Cell type identification
DefaultAssay(CtrlvTwoM.combined.FBSM)<-"RNA"
tiff(file = "CtrlvTwoM.combined.FBSM celltype marker expression plots.tiff", width = 10, height = 10, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvTwoM.combined.FBSM, reduction = "umap", features = c("EGFP", "Pdgfra", "Cd34", "Foxl1", "Gli1", "Myh11"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#higher resolution
DefaultAssay(CtrlvTwoM.combined.FBSM) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvTwoM.combined.FBSM <- ScaleData(CtrlvTwoM.combined.FBSM, verbose = FALSE)
CtrlvTwoM.combined.FBSM <- RunPCA(CtrlvTwoM.combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvTwoM.combined.FBSM, ndims = 50)
#Umap and Clustering
DefaultAssay(CtrlvTwoM.combined.FBSM) <- "integrated"
CtrlvTwoM.combined.FBSM <- FindNeighbors(CtrlvTwoM.combined.FBSM, reduction = "pca", dims = 1:15)
CtrlvTwoM.combined.FBSM <- FindClusters(CtrlvTwoM.combined.FBSM, resolution = 1)
CtrlvTwoM.combined.FBSM <- RunUMAP(CtrlvTwoM.combined.FBSM, reduction = "pca", dims = 1:15)

Idents(object = CtrlvTwoM.combined.FBSM) <- "seurat_clusters"
DimPlot(CtrlvTwoM.combined.FBSM, reduction = "umap", label = TRUE)
DimPlot(CtrlvTwoM.combined.FBSM, reduction = "umap", split.by="stim", label = TRUE)

Idents(object = CtrlvTwoM.combined.FBSM) <- "seurat_clusters"
CtrlvTwoM.combined.FBSM1 <- subset(CtrlvTwoM.combined.FBSM, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12", "13"))
CtrlvTwoM.combined.FBSM1 <- RunUMAP(CtrlvTwoM.combined.FBSM1, reduction = "pca", dims = 1:15)
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", label = TRUE)
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(CtrlvTwoM.combined.FBSM1) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvTwoM.combined.FBSM1 <- ScaleData(CtrlvTwoM.combined.FBSM1, verbose = FALSE)
CtrlvTwoM.combined.FBSM1 <- RunPCA(CtrlvTwoM.combined.FBSM1, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvTwoM.combined.FBSM1, ndims = 50)
#Umap and Clustering
DefaultAssay(CtrlvTwoM.combined.FBSM1) <- "integrated"
CtrlvTwoM.combined.FBSM1 <- FindNeighbors(CtrlvTwoM.combined.FBSM1, reduction = "pca", dims = 1:20)
CtrlvTwoM.combined.FBSM1 <- FindClusters(CtrlvTwoM.combined.FBSM1, resolution = 2.5)
CtrlvTwoM.combined.FBSM1 <- RunUMAP(CtrlvTwoM.combined.FBSM1, reduction = "pca", dims = 1:20)

Idents(object = CtrlvTwoM.combined.FBSM1) <- "seurat_clusters"
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", label = TRUE)
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", split.by="stim", label = TRUE)


DefaultAssay(CtrlvTwoM.combined.FBSM1)<-"RNA"
FeaturePlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")


#Rename
Idents(object = CtrlvTwoM.combined.FBSM1) <- "seurat_clusters"
CtrlvTwoM.combined.FBSM1 <- RenameIdents(object = CtrlvTwoM.combined.FBSM1, 
                                         '5' = "Telocyte", 
                                         '0' = "FB",  '1' = "FB",'2' = "FB",'3' = "FB",'4' = "FB",
                                        '11' = "FB",'6' = "FB",'7' = "FB",'8' = "FB",'9' = "FB",
                                        '10' = "FB",'13' = "FB",'14' = "FB",'15' = "FB",'16' = "FB",
                                        '17' = "FB",'18' = "FB",'19' = "FB",'22' = "FB",'23' = "FB",'24' = "FB",
                                        '20'="SM", '12'="SM", '21'="SM")
CtrlvTwoM.combined.FBSM1[["FBSMCellTypes1"]] <- Idents(object = CtrlvTwoM.combined.FBSM1)
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", cols = c("red","blue","darkorange","skyblue1", "bisque3","palevioletred3",
                                                               "green4"))


Idents(object = CtrlvTwoM.combined.FBSM1) <- "FBSMCellTypes1"
tiff(file = "CtrlvTwoM.combined.FBSM1 FBSMCellTypes1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", cols = c("red","blue","darkorange","skyblue1", "bisque3","palevioletred3",
                                                               "green4"))
dev.off()
tiff(file = "CtrlvTwoM.combined.FBSM1 FBSMCellTypes1 split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvTwoM.combined.FBSM1, reduction = "umap", split.by="stim", cols = c("red","blue","darkorange","lightblue", "bisque3","palevioletred3",
                                                                                "green4"))
dev.off()

##Featurplots
#integrated

Idents(object = CtrlvTwoM.combined.FBSM1) <- "stim"
Ctrl.combined.FBSM1 <- subset(CtrlvTwoM.combined.FBSM1, idents = c("Ctrl"))
TwoM.combined.FBSM1 <- subset(CtrlvTwoM.combined.FBSM1, idents = c("ARQ9"))


DefaultAssay(Ctrl.combined.FBSM1) <- "RNA"
tiff(file = "Ctrl.combined.FBSM1 PDGFRa.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.FBSM1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "Ctrl.combined.FBSM1 Cd34.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.FBSM1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q60")
dev.off()
tiff(file = "Ctrl.combined.FBSM1 Foxl1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.FBSM1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q50")
dev.off()
tiff(file = "Ctrl.combined.FBSM1 Ar.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q80")
dev.off()
tiff(file = "Ctrl.combined.FBSM1 Gli1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Ctrl.combined.FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q70")
dev.off()

DefaultAssay(TwoM.combined.FBSM1) <- "RNA"
tiff(file = "TwoM.combined.FBSM1 PDGFRa.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM1 Cd34.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "TwoM.combined.FBSM1 Foxl1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "TwoM.combined.FBSM1 Ar.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM1 Gli1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q80")
dev.off()

Idents(object = CtrlvTwoM.combined.FBSM1) <- "FBSMCellTypes1"
CtrlvTwoM.combined.FBTC <- subset(CtrlvTwoM.combined.FBSM1, idents = c("FB", "Telocyte"))

DefaultAssay(CtrlvTwoM.combined.FBTC) <- "RNA"
Idents(object = CtrlvTwoM.combined.FBTC) <- "FBSMCellTypes1"
tiff(file = "CtrlvTwoM.combined.FBTC vlnplot EGFP split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvTwoM.combined.FBTC, features = "EGFP", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvTwoM.combined.FBTC vlnplot Foxl1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvTwoM.combined.FBTC, features = "Foxl1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvTwoM.combined.FBTC vlnplot Pdgfra split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvTwoM.combined.FBTC, features = "Pdgfra", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvTwoM.combined.FBTC vlnplot Gli1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvTwoM.combined.FBTC, features = "Gli1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "CtrlvTwoM.combined.FBTC vlnplot CD34 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvTwoM.combined.FBTC, features = "Cd34", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()


####TwoMvWT####
Idents(object = TwoMvWT.combined1.FBSM3) <- "FBSMCellTypes"
DimPlot(TwoMvWT.combined1.FBSM3, reduction = "umap", label = TRUE)
DimPlot(TwoMvWT.combined1.FBSM3, reduction = "umap", split.by="stim", label = TRUE)

#higher resolution
DefaultAssay(TwoMvWT.combined1.FBSM3) <- "integrated"
#Run the standard workflow for visualization and clustering
TwoMvWT.combined1.FBSM3 <- ScaleData(TwoMvWT.combined1.FBSM3, verbose = FALSE)
TwoMvWT.combined1.FBSM3 <- RunPCA(TwoMvWT.combined1.FBSM3, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined1.FBSM3, ndims = 50)
#Umap and Clustering
DefaultAssay(TwoMvWT.combined1.FBSM3) <- "integrated"
TwoMvWT.combined1.FBSM3 <- FindNeighbors(TwoMvWT.combined1.FBSM3, reduction = "pca", dims = 1:25)
TwoMvWT.combined1.FBSM3 <- FindClusters(TwoMvWT.combined1.FBSM3, resolution = 2)
TwoMvWT.combined1.FBSM3 <- RunUMAP(TwoMvWT.combined1.FBSM3, reduction = "pca", dims = 1:25)

Idents(object = TwoMvWT.combined1.FBSM3) <- "seurat_clusters"
DimPlot(TwoMvWT.combined1.FBSM3, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(TwoMvWT.combined1.FBSM3)<-"RNA"
FeaturePlot(TwoMvWT.combined1.FBSM3, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TwoMvWT.combined1.FBSM3, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TwoMvWT.combined1.FBSM3, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TwoMvWT.combined1.FBSM3, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TwoMvWT.combined1.FBSM3, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

