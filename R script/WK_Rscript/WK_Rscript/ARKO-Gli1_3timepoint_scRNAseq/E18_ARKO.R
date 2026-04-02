#Gli1-lineage cells_E18.5_ARKO

#### Add necessary tools to library ####

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

#### E18_Initial Filtering and Clustering ####
setwd("/Volumes/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/E18.5/ARKO")

E18_ARKO.unfiltered.data <- Read10X("/Volumes/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/CellRanger_Output/count_Project_COHP_36981_1_X3SC3/outs/filtered_feature_bc_matrix")
E18_ARKO.unfiltered <- CreateSeuratObject(counts = E18_ARKO.unfiltered.data,  min.cells = 3, min.features = 200, project = "E18_ARKO")
E18_ARKO.unfiltered <- NormalizeData(E18_ARKO.unfiltered)

E18_ARKO.unfiltered[["percent.mt"]] <- PercentageFeatureSet(E18_ARKO.unfiltered, pattern = "^mt-")
tiff(file = "E18_ARKO Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_ARKO.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "E18_ARKO Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_ARKO.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E18_ARKO Pre-filteration")
dev.off()
tiff(file = "E18_ARKO Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_ARKO.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E18_ARKO Pre-filteration")
dev.off()

plot1 <- FeatureScatter(E18_ARKO.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E18_ARKO.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E18_ARKO <- subset(E18_ARKO.unfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 10)
tiff(file = "E18_ARKO Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_ARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "E18_ARKO Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_ARKO@meta.data$nFeature_RNA, breaks = 100, col = "skyblue", xlab = "nFeature_RNA", main = "E18_ARKO Post-filteration")
dev.off()
tiff(file = "E18_ARKO Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_ARKO@meta.data$percent.mt, breaks = 100, col = "skyblue", xlab = "percent.mt", main = "E18_ARKO Post-filteration")
dev.off()

plot1 <- FeatureScatter(E18_ARKO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E18_ARKO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E18_ARKO <- FindVariableFeatures(E18_ARKO, selection.method = "vst", nfeatures = 2500)
VariableFeaturePlot(E18_ARKO)

#new filtering paramaters
DefaultAssay(E18_ARKO) <- "RNA"
all.genes <- rownames(E18_ARKO)
E18_ARKO <- ScaleData(E18_ARKO, features = all.genes)
E18_ARKO <- RunPCA(E18_ARKO, features = VariableFeatures(object = E18_ARKO))
ElbowPlot(E18_ARKO, ndims = 50)

E18_ARKO <- FindNeighbors(E18_ARKO, dims = 1:25)
E18_ARKO <- FindClusters(E18_ARKO, resolution = 0.5)
E18_ARKO <- RunUMAP(E18_ARKO, dims = 1:25)
DimPlot(E18_ARKO, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_ARKO) <- "stim"
tiff(file = "E18_ARKO UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO, reduction = "umap", pt.size = 0.3, cols = c("dark blue")) 
dev.off()

#Cell cycle scoring
DefaultAssay(E18_ARKO) <- "RNA"
all.genes <- rownames(E18_ARKO)
E18_ARKO <- ScaleData(E18_ARKO, features = all.genes)
E18_ARKO <- CellCycleScoring(E18_ARKO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_ARKO) <- "Phase"
DimPlot(E18_ARKO, reduction = "umap")
tiff(file = "E18_ARKO Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO, reduction = "umap", pt.size = 0.3, cols = c("orange", "magenta2","lightseagreen"))
dev.off()

#Cell Cycle regression
E18_ARKO1 <- E18_ARKO
DefaultAssay(E18_ARKO1) <- "RNA"
E18_ARKO1 <- ScaleData(E18_ARKO1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_ARKO1))
E18_ARKO1 <- RunPCA(E18_ARKO1, features = VariableFeatures(E18_ARKO1))
ElbowPlot(E18_ARKO1, ndims = 50)

E18_ARKO1 <- FindNeighbors(E18_ARKO1, reduction = "pca", dims = 1:23)
E18_ARKO1 <- FindClusters(E18_ARKO1, resolution = 0.8)
E18_ARKO1 <- RunUMAP(E18_ARKO1, reduction = "pca", dims = 1:23)

Idents(object = E18_ARKO1) <- "seurat_clusters"
DimPlot(E18_ARKO1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_ARKO1) <- "seurat_clusters"
tiff(file = "E18_ARKO1 UMAP dims23 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_ARKO1) <- "Phase"
tiff(file = "E18_ARKO1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1, reduction = "umap", pt.size = 0.3, cols = c( "orange", "magenta2", "lightseagreen"))
dev.off()

#UGS identification
DefaultAssay(E18_ARKO1)<-"RNA"
tiff(file = "E18_ARKO1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt14", 
                                                        "Upk1b", "Upk3b", "Pax8", 
                                                        "Vim", "Acta2", "Myh11", "Rgs5", "Plp1",
                                                        "Pecam1", "Tyrobp", "Rgs1", "Myog", "Syp"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#FeaturePlots
DefaultAssay(E18_ARKO1)<-"RNA"
tiff(file = "E18_ARKO1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_ARKO1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "E18_ARKO1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
Idents(object = E18_ARKO1) <- "seurat_clusters"
DefaultAssay(E18_ARKO1) <- "RNA"
all.genes <- rownames(E18_ARKO1)
E18_ARKO1 <- ScaleData(E18_ARKO1, features = all.genes)
E18_ARKO1.seurat_cluster.markers <- FindAllMarkers(E18_ARKO1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(E18_ARKO1.seurat_cluster.markers, file = "E18_ARKO1.seurat_cluster.markers.csv")

#Rename
Idents(object = E18_ARKO1) <- "seurat_clusters"
E18_ARKO1 <- RenameIdents(object = E18_ARKO1, 
                          '6' = "UGE", '12' ="UGE",  
                          '16' = "Urothelium", 
                          '19' = "WD",
                          '0' = "Mesenchyme",'1' = "Mesenchyme",'5' = "Mesenchyme",'4' = "Mesenchyme",
                          '11' = "Mesenchyme",'3' = "Mesenchyme",'2' = "Mesenchyme",'9' = "Mesenchyme",
                          '10' = "Mesenchyme",'8' = "Mesenchyme",'7' = "Mesenchyme",
                          '18' = "MyoBlast", 
                          '15' = "Pericyte",'14' = "Glia", '17' = "VE",  
                          '13' = "NE", 
                          '20'="Immune")
E18_ARKO1[["CellTypes"]] <- Idents(object = E18_ARKO1)

#CellTypes UMAP
Idents(object = E18_ARKO1) <- "CellTypes"
tiff(file = "E18_ARKO1 UMAP CellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red",  "green4", "bisque3", 
                                                                             "black"))
dev.off()

#EGFP expression
DefaultAssay(E18_ARKO1) <- "RNA"
E18_ARKO1EGFPPos <- subset(x=E18_ARKO1,  subset = `EGFP` > 1.3)
E18_ARKO1EGFPNeg <- subset(x=E18_ARKO1,  subset = `EGFP` < 1.3)
Idents(object = E18_ARKO1EGFPPos) <- "EGFPPos"
Idents(object = E18_ARKO1EGFPNeg) <- "EGFPNeg"
E18_ARKO1EGFPPos[["EGFPExp"]] <- Idents(object = E18_ARKO1EGFPPos)
E18_ARKO1EGFPNeg[["EGFPExp"]] <- Idents(object = E18_ARKO1EGFPNeg)
E18_ARKO1EGFP <- merge(x = E18_ARKO1EGFPPos, y = E18_ARKO1EGFPNeg)
Idents(object = E18_ARKO1EGFP) <- "EGFPExp"
E18_ARKO1$EGFPExp <- Idents(object = E18_ARKO1EGFP)
Idents(object = E18_ARKO1) <- "EGFPExp"
DimPlot(E18_ARKO1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Gli1 expression
DefaultAssay(E18_ARKO1) <- "RNA"
E18_ARKO1Gli1Pos <- subset(x=E18_ARKO1,  subset = `Gli1` > 0)
E18_ARKO1Gli1Neg <- subset(x=E18_ARKO1,  subset = `Gli1` == 0)
Idents(object = E18_ARKO1Gli1Pos) <- "Gli1Pos"
Idents(object = E18_ARKO1Gli1Neg) <- "Gli1Neg"
E18_ARKO1Gli1Pos[["Gli1Exp"]] <- Idents(object = E18_ARKO1Gli1Pos)
E18_ARKO1Gli1Neg[["Gli1Exp"]] <- Idents(object = E18_ARKO1Gli1Neg)
E18_ARKO1Gli1 <- merge(x = E18_ARKO1Gli1Pos, y = E18_ARKO1Gli1Neg)
Idents(object = E18_ARKO1Gli1) <- "Gli1Exp"
E18_ARKO1$Gli1Exp <- Idents(object =E18_ARKO1Gli1)
Idents(object = E18_ARKO1) <- "Gli1Exp"
DimPlot(E18_ARKO1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Cell counts
Idents(object = E18_ARKO1) <- "CellTypes"
E18_ARKO1$EGFPExp.CellTypes <- paste(Idents(E18_ARKO1), E18_ARKO1$EGFPExp, sep = "_")
Idents(object = E18_ARKO1) <- "EGFPExp.CellTypes"
table(Idents(E18_ARKO1))

#Cell counts
Idents(object = E18_ARKO1) <- "CellTypes"
E18_ARKO1$Gli1Exp.CellTypes <- paste(Idents(E18_ARKO1), E18_ARKO1$Gli1Exp, sep = "_")
Idents(object = E18_ARKO1) <- "Gli1Exp.CellTypes"
table(Idents(E18_ARKO1))

####Subclustering Stro####
Idents(object = E18_ARKO1) <- "CellTypes"
E18_ARKO1_Stro <- subset(E18_ARKO1, idents = c("Mesenchyme", "MyoBlast", "Pericyte", "Glia", "VE", "Immune"))

DefaultAssay(E18_ARKO1_Stro) <- "RNA"
#Run the standard workflow for visualization and clustering
E18_ARKO1_Stro <- ScaleData(E18_ARKO1_Stro, verbose = FALSE)
E18_ARKO1_Stro <- RunPCA(E18_ARKO1_Stro, npcs = 30, verbose = FALSE)
ElbowPlot(E18_ARKO1_Stro, ndims = 50)
# UMAP and Clustering
E18_ARKO1_Stro <- FindNeighbors(E18_ARKO1_Stro, reduction = "pca", dims = 1:23)
E18_ARKO1_Stro <- FindClusters(E18_ARKO1_Stro, resolution = 0.8)
E18_ARKO1_Stro <- RunUMAP(E18_ARKO1_Stro, reduction = "pca", dims = 1:23)
DimPlot(E18_ARKO1_Stro, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(E18_ARKO1_Stro) <- "RNA"
all.genes <- rownames(E18_ARKO1_Stro)
E18_ARKO1_Stro <- ScaleData(E18_ARKO1_Stro, features = all.genes)
E18_ARKO1_Stro <- CellCycleScoring(E18_ARKO1_Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_ARKO1_Stro) <- "Phase"
DimPlot(E18_ARKO1_Stro, reduction = "umap")

tiff(file = "E18_ARKO1_Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_Stro, reduction = "umap", pt.size = 0.3, cols = c("orange", "magenta2", "lightseagreen"))
dev.off()

#Cell Cycle regression
E18_ARKO1_Stro1 <- E18_ARKO1_Stro
DefaultAssay(E18_ARKO1_Stro1) <- "RNA"
E18_ARKO1_Stro1 <- ScaleData(E18_ARKO1_Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_ARKO1_Stro1))
E18_ARKO1_Stro1 <- RunPCA(E18_ARKO1_Stro1, features = VariableFeatures(E18_ARKO1_Stro1))
ElbowPlot(E18_ARKO1_Stro1, ndims = 50)

E18_ARKO1_Stro1 <- FindNeighbors(E18_ARKO1_Stro1, reduction = "pca", dims = 1:24)
E18_ARKO1_Stro1 <- FindClusters(E18_ARKO1_Stro1, resolution = 1.3)
E18_ARKO1_Stro1 <- RunUMAP(E18_ARKO1_Stro1, reduction = "pca", dims = 1:24)
DimPlot(E18_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_ARKO1_Stro1) <- "seurat_clusters"
tiff(file = "E18_ARKO1_Stro1 UMAP dims24 res1.3.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_ARKO1_Stro1) <- "Phase"
DimPlot(E18_ARKO1_Stro, reduction = "umap")

tiff(file = "E18_ARKO1_Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, cols = c("orange", "magenta2","lightseagreen"))
dev.off()

#E18_ARKO1_Stro1 celltype identification
DefaultAssay(E18_ARKO1_Stro1)<-"RNA"
tiff(file = "E18_ARKO1_Stro1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1_Stro1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt14", 
                                                              "Upk1b", "Upk3b", "Pax8", 
                                                              "Vim", "Acta2", "Myh11", "Rgs5", "Plp1",
                                                              "Pecam1", "Tyrobp", "Crabp1", "Myog", "Syp"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = E18_ARKO1_Stro1) <- "seurat_clusters"
E18_ARKO1_Stro1 <- RenameIdents(object = E18_ARKO1_Stro1, 
                                '7' = "FB", '8' = "FB",  '10' = "FB", '15' = "FB", '13' = "FB",
                                '17' = "FB", '5' = "FB",'1' = "FB",'2' = "FB",
                                '3' = "FB",'4' = "FB",'9' = "FB",'12' = "FB", '0'="WDM",
                                '18'="Myoblast",
                                '6' = "SM",'19' = "SM",
                                '14' = "Pericyte", '11' = "Glia",'16' = "VE", 
                                '20'="Immune",
                                '21'="OS")
E18_ARKO1_Stro1[["StroCellTypes"]] <- Idents(object = E18_ARKO1_Stro1)
DimPlot(E18_ARKO1_Stro1, reduction = "umap")

#CellTypes UMAP
Idents(object = E18_ARKO1_Stro1) <- "StroCellTypes"
tiff(file = "E18_ARKO1_Stro1 UMAP StroCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

####Subclustering FBSM####
Idents(object = E18_ARKO1_Stro1) <- "StroCellTypes"
E18_ARKO1_FBSM <- subset(E18_ARKO1_Stro1, idents = c("FB", "SM"))

DefaultAssay(E18_ARKO1_FBSM) <- "RNA"
#Run the standard workflow for visualization and clustering
E18_ARKO1_FBSM <- ScaleData(E18_ARKO1_FBSM, verbose = FALSE)
E18_ARKO1_FBSM <- RunPCA(E18_ARKO1_FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(E18_ARKO1_FBSM, ndims = 50)
# UMAP and Clustering
E18_ARKO1_FBSM <- FindNeighbors(E18_ARKO1_FBSM, reduction = "pca", dims = 1:22)
E18_ARKO1_FBSM <- FindClusters(E18_ARKO1_FBSM, resolution = 0.5)
E18_ARKO1_FBSM <- RunUMAP(E18_ARKO1_FBSM, reduction = "pca", dims = 1:22)
DimPlot(E18_ARKO1_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(E18_ARKO1_FBSM) <- "RNA"
all.genes <- rownames(E18_ARKO1_FBSM)
E18_ARKO1_FBSM <- ScaleData(E18_ARKO1_FBSM, features = all.genes)
E18_ARKO1_FBSM <- CellCycleScoring(E18_ARKO1_FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_ARKO1_FBSM) <- "Phase"
DimPlot(E18_ARKO1_FBSM, reduction = "umap")
tiff(file = "E18_ARKO1_FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_FBSM, reduction = "umap", pt.size = 0.3, cols = c("orange","magenta2","lightseagreen"))
dev.off()

#Cell Cycle regression
E18_ARKO1_FBSM1 <- E18_ARKO1_FBSM
DefaultAssay(E18_ARKO1_FBSM1) <- "RNA"
E18_ARKO1_FBSM1 <- ScaleData(E18_ARKO1_FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_ARKO1_FBSM1))
E18_ARKO1_FBSM1 <- RunPCA(E18_ARKO1_FBSM1, features = VariableFeatures(E18_ARKO1_FBSM1))
ElbowPlot(E18_ARKO1_FBSM1, ndims = 50)

E18_ARKO1_FBSM1 <- FindNeighbors(E18_ARKO1_FBSM1, reduction = "pca", dims = 1:25)
E18_ARKO1_FBSM1 <- FindClusters(E18_ARKO1_FBSM1, resolution = 0.8)
E18_ARKO1_FBSM1 <- RunUMAP(E18_ARKO1_FBSM1, reduction = "pca", dims = 1:25)

Idents(object = E18_ARKO1_FBSM1) <- "seurat_clusters"
DimPlot(E18_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_ARKO1_FBSM1) <- "seurat_clusters"
tiff(file = "E18_ARKO1_FBSM1 UMAP dims25 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_ARKO1_FBSM1) <- "Phase"
tiff(file = "E18_ARKO1_FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, cols = c("orange","magenta2","lightseagreen"))
dev.off()

#E18_ARKO1_FBSM1 celltype identification
DefaultAssay(E18_ARKO1_FBSM1)<-"RNA"
tiff(file = "E18_ARKO1_FBSM1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1_FBSM1, reduction = "umap", features = c(
  "Fbln1", "Dcn", "Vim", "Col1a1", 
  "Myh11", "Acta2", "Actg2", "Tagln", "Thy1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = E18_ARKO1_FBSM1) <- "seurat_clusters"
E18_ARKO1_FBSM1 <- RenameIdents(object = E18_ARKO1_FBSM1, 
                                '1' = "FB1",'2' = "FB1", '5' = "FB2", '4' = "FB3", 
                                '9' = "FB4", '0' = "FB5",'10' = "FB6",
                                '8' = "FB7", '3' = "FB8",'6' = "FB9",'11' = "FB10",
                                '12' = "SM1", '7' = "SM2" 
                                )
E18_ARKO1_FBSM1[["FBSMCellTypes"]] <- Idents(object = E18_ARKO1_FBSM1)

#CellTypes UMAP
Idents(object = E18_ARKO1_FBSM1) <- "FBSMCellTypes"
tiff(file = "E18_ARKO1_FBSM1 UMAP FBSMCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

#Cell counts
Idents(object = E18_ARKO1_FBSM1) <- "FBSMCellTypes"
E18_ARKO1_FBSM1$EGFPExp.FBSMCellTypes <- paste(Idents(E18_ARKO1_FBSM1), E18_ARKO1_FBSM1$EGFPExp, sep = "_")
Idents(object = E18_ARKO1_FBSM1) <- "EGFPExp.FBSMCellTypes"
table(Idents(E18_ARKO1_FBSM1))

#Cell counts
Idents(object = E18_ARKO1_FBSM1) <- "FBSMCellTypes"
E18_ARKO1_FBSM1$Gli1Exp.FBSMCellTypes <- paste(Idents(E18_ARKO1_FBSM1), E18_ARKO1_FBSM1$Gli1Exp, sep = "_")
Idents(object = E18_ARKO1_FBSM1) <- "Gli1Exp.FBSMCellTypes"
table(Idents(E18_ARKO1_FBSM1))

#FeaturePlots
DefaultAssay(E18_ARKO1_FBSM1)<-"RNA"
tiff(file = "E18_ARKO1_FBSM1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1_FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_ARKO1_FBSM1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1_FBSM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "E18_ARKO1_FBSM1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_ARKO1_FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

