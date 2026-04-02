#Gli1-lineage cells_P42_ARKO

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

#### P11_Initial Filtering and Clustering ####
setwd("/Volumes/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/P42/ARKO")

P42_ARKO.unfiltered.data <- Read10X("/Volumes/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/CellRanger_Output/count_Project_COHP_36979_1_X3SC3/outs/filtered_feature_bc_matrix")
P42_ARKO.unfiltered <- CreateSeuratObject(counts = P42_ARKO.unfiltered.data,  min.cells = 3, min.features = 200, project = "P42_ARKO")
P42_ARKO.unfiltered <- NormalizeData(P42_ARKO.unfiltered)

#P42_ARKO
#Initial processing & filtering
P42_ARKO.unfiltered[["percent.mt"]] <- PercentageFeatureSet(P42_ARKO.unfiltered, pattern = "^mt-")
tiff(file = "P42_ARKO Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(P42_ARKO.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "P42_ARKO Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_ARKO.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P42_ARKO Pre-filteration")
dev.off()
tiff(file = "P42_ARKO Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_ARKO.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P42_ARKO Pre-filteration")
dev.off()

plot1 <- FeatureScatter(P42_ARKO.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P42_ARKO.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P42_ARKO <- subset(P42_ARKO.unfiltered, subset = nFeature_RNA > 800 & nFeature_RNA < 7500 & percent.mt < 15)
tiff(file = "P42_ARKO Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(P42_ARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "P42_ARKO Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_ARKO@meta.data$nFeature_RNA, breaks = 100, col = "skyblue", xlab = "nFeature_RNA", main = "P42_ARKO Post-filteration")
dev.off()
tiff(file = "P42_ARKO Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_ARKO@meta.data$percent.mt, breaks = 100, col = "skyblue", xlab = "percent.mt", main = "P42_ARKO Post-filteration")
dev.off()

plot1 <- FeatureScatter(P42_ARKO, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P42_ARKO, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P42_ARKO <- FindVariableFeatures(P42_ARKO, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(P42_ARKO)

#new filtering paramaters
#Clustering
DefaultAssay(P42_ARKO) <- "RNA"
all.genes <- rownames(P42_ARKO)
P42_ARKO <- ScaleData(P42_ARKO, features = all.genes)
P42_ARKO <- RunPCA(P42_ARKO, features = VariableFeatures(object = P42_ARKO))
ElbowPlot(P42_ARKO, ndims = 50)

P42_ARKO <- FindNeighbors(P42_ARKO, dims = 1:26)
P42_ARKO <- FindClusters(P42_ARKO, resolution = 0.5)
P42_ARKO <- RunUMAP(P42_ARKO, dims = 1:26)
DimPlot(P42_ARKO, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle scoring
DefaultAssay(P42_ARKO) <- "RNA"
all.genes <- rownames(P42_ARKO)
P42_ARKO <- ScaleData(P42_ARKO, features = all.genes)
P42_ARKO <- CellCycleScoring(P42_ARKO, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P42_ARKO) <- "Phase"
DimPlot(P42_ARKO, reduction = "umap")
tiff(file = "P42_ARKO Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO, reduction = "umap", pt.size = 0.3, cols = c("orange","magenta2","lightseagreen"))
dev.off()

#Cell Cycle regression
P42_ARKO1 <- P42_ARKO
DefaultAssay(P42_ARKO1) <- "RNA"
P42_ARKO1 <- ScaleData(P42_ARKO1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P42_ARKO1))
P42_ARKO1 <- RunPCA(P42_ARKO1, features = VariableFeatures(P42_ARKO1))
ElbowPlot(P42_ARKO1, ndims = 50)

P42_ARKO1 <- FindNeighbors(P42_ARKO1, reduction = "pca", dims = 1:28)
P42_ARKO1 <- FindClusters(P42_ARKO1, resolution = 1)
P42_ARKO1 <- RunUMAP(P42_ARKO1, reduction = "pca", dims = 1:28)
DimPlot(P42_ARKO1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_ARKO1) <- "seurat_clusters"
tiff(file = "P42_ARKO1 UMAP dims28 res1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_ARKO1) <- "Phase"
tiff(file = "P42_ARKO1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1, reduction = "umap", pt.size = 0.3, cols = c( "orange", "magenta2","lightseagreen"))
dev.off()

#P11 celltype identification
DefaultAssay(P42_ARKO1)<-"RNA"
tiff(file = "P42_ARKO1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Cdo1", 
                                                        "Vim",  "Fbln1", "Myh11", "Acta2", "Tagln", "Actg2",
                                                        "Pecam1", "Tyrobp", "Plp1", "Syp", "Rgs5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#FeaturePlots
DefaultAssay(P42_ARKO1)<-"RNA"
tiff(file = "P42_ARKO1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_ARKO1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "P42_ARKO1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
Idents(object = P42_ARKO1) <- "seurat_clusters"
DefaultAssay(P42_ARKO1) <- "RNA"
all.genes <- rownames(P42_ARKO1)
P42_ARKO1 <- ScaleData(P42_ARKO1, features = all.genes)
P42_ARKO1.seurat_cluster.markers <- FindAllMarkers(P42_ARKO1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(P42_ARKO1.seurat_cluster.markers, file = "P42_ARKO1.seurat_cluster.markers.csv")

#Rename
Idents(object = P42_ARKO1) <- "seurat_clusters"
P42_ARKO1 <- RenameIdents(object = P42_ARKO1, 
                          '7' = "BE", '9' ="BE",  '15' ="BE",'10' ="BE",
                          '18' = "LE", 
                          '8' = "SV",'20' = "SV",'13' = "SV",
                          '3' = "FB",'2' = "FB",'0' = "FB",'6' = "FB",
                          '1' = "SM", '17' = "SM",
                          '16' = "Pericyte",'19' = "Pericyte", '25' = "Pericyte",'21' = "Glia",  '12' = "VE",  
                          '4'="Immune", '14'="Immune", '11'="Immune", '5'="Immune", '24'="Immune"
                          , '22'="Immune", '23'="Immune")
P42_ARKO1[["CellTypes"]] <- Idents(object = P42_ARKO1)

#CellTypes UMAP
Idents(object = P42_ARKO1) <- "CellTypes"
tiff(file = "P42_ARKO1 UMAP CellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red",  "green4", "bisque3", 
                                                                             "black"))
dev.off()

#EGFP expression
DefaultAssay(P42_ARKO1) <- "RNA"
P42_ARKO1EGFPPos <- subset(x=P42_ARKO1,  subset = `EGFP` > 1.3)
P42_ARKO1EGFPNeg <- subset(x=P42_ARKO1,  subset = `EGFP` < 1.3)
Idents(object = P42_ARKO1EGFPPos) <- "EGFPPos"
Idents(object = P42_ARKO1EGFPNeg) <- "EGFPNeg"
P42_ARKO1EGFPPos[["EGFPExp"]] <- Idents(object = P42_ARKO1EGFPPos)
P42_ARKO1EGFPNeg[["EGFPExp"]] <- Idents(object = P42_ARKO1EGFPNeg)
P42_ARKO1EGFP <- merge(x = P42_ARKO1EGFPPos, y = P42_ARKO1EGFPNeg)
Idents(object = P42_ARKO1EGFP) <- "EGFPExp"
P42_ARKO1$EGFPExp <- Idents(object = P42_ARKO1EGFP)
Idents(object = P42_ARKO1) <- "EGFPExp"
DimPlot(P42_ARKO1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Gli1 expression
DefaultAssay(P42_ARKO1) <- "RNA"
P42_ARKO1Gli1Pos <- subset(x=P42_ARKO1,  subset = `Gli1` > 0)
P42_ARKO1Gli1Neg <- subset(x=P42_ARKO1,  subset = `Gli1` == 0)
Idents(object = P42_ARKO1Gli1Pos) <- "Gli1Pos"
Idents(object = P42_ARKO1Gli1Neg) <- "Gli1Neg"
P42_ARKO1Gli1Pos[["Gli1Exp"]] <- Idents(object = P42_ARKO1Gli1Pos)
P42_ARKO1Gli1Neg[["Gli1Exp"]] <- Idents(object = P42_ARKO1Gli1Neg)
P42_ARKO1Gli1 <- merge(x = P42_ARKO1Gli1Pos, y = P42_ARKO1Gli1Neg)
Idents(object = P42_ARKO1Gli1) <- "Gli1Exp"
P42_ARKO1$Gli1Exp <- Idents(object = P42_ARKO1Gli1)
Idents(object = P42_ARKO1) <- "Gli1Exp"
DimPlot(P42_ARKO1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Cell counts
Idents(object = P42_ARKO1) <- "CellTypes"
P42_ARKO1$EGFPExp.CellTypes <- paste(Idents(P42_ARKO1), P42_ARKO1$EGFPExp, sep = "_")
Idents(object = P42_ARKO1) <- "EGFPExp.CellTypes"
table(Idents(P42_ARKO1))

#Cell counts
Idents(object = P42_ARKO1) <- "CellTypes"
P42_ARKO1$Gli1Exp.CellTypes <- paste(Idents(P42_ARKO1), P42_ARKO1$Gli1Exp, sep = "_")
Idents(object = P42_ARKO1) <- "Gli1Exp.CellTypes"
table(Idents(P42_ARKO1))

####Subclustering Stro####
Idents(object = P42_ARKO1) <- "CellTypes"
P42_ARKO1_Stro <- subset(P42_ARKO1, idents = c("FB", "SM", "Pericyte", "Glia", "VE", "Immune"))

DefaultAssay(P42_ARKO1_Stro) <- "RNA"
#Run the standard workflow for visualization and clustering
P42_ARKO1_Stro <- ScaleData(P42_ARKO1_Stro, verbose = FALSE)
P42_ARKO1_Stro <- RunPCA(P42_ARKO1_Stro, npcs = 30, verbose = FALSE)
ElbowPlot(P42_ARKO1_Stro, ndims = 50)
# UMAP and Clustering
P42_ARKO1_Stro <- FindNeighbors(P42_ARKO1_Stro, reduction = "pca", dims = 1:21)
P42_ARKO1_Stro <- FindClusters(P42_ARKO1_Stro, resolution = 0.8)
P42_ARKO1_Stro <- RunUMAP(P42_ARKO1_Stro, reduction = "pca", dims = 1:21)
DimPlot(P42_ARKO1_Stro, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(P42_ARKO1_Stro) <- "RNA"
all.genes <- rownames(P42_ARKO1_Stro)
P42_ARKO1_Stro <- ScaleData(P42_ARKO1_Stro, features = all.genes)
P42_ARKO1_Stro <- CellCycleScoring(P42_ARKO1_Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P42_ARKO1_Stro) <- "Phase"
DimPlot(P42_ARKO1_Stro, reduction = "umap")
tiff(file = "P42_ARKO1_Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_Stro, reduction = "umap", pt.size = 0.3, cols = c("orange","lightseagreen", "magenta2"))
dev.off()

#Cell Cycle regression
P42_ARKO1_Stro1 <- P42_ARKO1_Stro
DefaultAssay(P42_ARKO1_Stro1) <- "RNA"
P42_ARKO1_Stro1 <- ScaleData(P42_ARKO1_Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P42_ARKO1_Stro1))
P42_ARKO1_Stro1 <- RunPCA(P42_ARKO1_Stro1, features = VariableFeatures(P42_ARKO1_Stro1))
ElbowPlot(P42_ARKO1_Stro1, ndims = 50)

P42_ARKO1_Stro1 <- FindNeighbors(P42_ARKO1_Stro1, reduction = "pca", dims = 1:21)
P42_ARKO1_Stro1 <- FindClusters(P42_ARKO1_Stro1, resolution = 1.5)
P42_ARKO1_Stro1 <- RunUMAP(P42_ARKO1_Stro1, reduction = "pca", dims = 1:21)
DimPlot(P42_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_ARKO1_Stro1) <- "seurat_clusters"
tiff(file = "P42_ARKO1_Stro1 UMAP dims21 res1.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_ARKO1_Stro1) <- "Phase"
tiff(file = "P42_ARKO1_Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#P42_ARKO1_Stro1 celltype identification
DefaultAssay(P42_ARKO1_Stro1)<-"RNA"
tiff(file = "P42_ARKO1_Stro1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1_Stro1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Cdo1", 
                                                              "Vim",  "Fbln1", "Myh11", "Acta2", "Tagln", "Actg2",
                                                              "Pecam1", "Tyrobp", "Plp1", "Syp", "Rgs5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = P42_ARKO1_Stro1) <- "seurat_clusters"
P42_ARKO1_Stro1 <- RenameIdents(object = P42_ARKO1_Stro1, 
                                '16' = "FB", '15' = "FB",  '7' = "FB", '1' = "FB", '2' = "FB",
                                '0' = "FB", '6' = "FB", '11'="MyoFB",
                                '4' = "SM",'5' = "SM",
                                '14' = "Pericyte", '13' = "Pericyte",'17' = "Glia", '10' = "VE", 
                                '21'="Immune",'3'="Immune",'12'="Immune",'23'="Immune",
                                '9'="Immune",'18'="Immune",'8'="Immune",'19'="Immune",
                                '20'="OS", '22'="OS")
P42_ARKO1_Stro1[["StroCellTypes"]] <- Idents(object = P42_ARKO1_Stro1)

#CellTypes UMAP
Idents(object = P42_ARKO1_Stro1) <- "StroCellTypes"
tiff(file = "P42_ARKO1_Stro1 UMAP StroCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

####Subclustering FBSM####
Idents(object = P42_ARKO1_Stro1) <- "StroCellTypes"
P42_ARKO1_FBSM <- subset(P42_ARKO1_Stro1, idents = c("FB", "MyoFB", "SM"))

DefaultAssay(P42_ARKO1_FBSM) <- "RNA"
#Run the standard workflow for visualization and clustering
P42_ARKO1_FBSM <- ScaleData(P42_ARKO1_FBSM, verbose = FALSE)
P42_ARKO1_FBSM <- RunPCA(P42_ARKO1_FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(P42_ARKO1_FBSM, ndims = 50)
# UMAP and Clustering
P42_ARKO1_FBSM <- FindNeighbors(P42_ARKO1_FBSM, reduction = "pca", dims = 1:18)
P42_ARKO1_FBSM <- FindClusters(P42_ARKO1_FBSM, resolution = 0.8)
P42_ARKO1_FBSM <- RunUMAP(P42_ARKO1_FBSM, reduction = "pca", dims = 1:18)
DimPlot(P42_ARKO1_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(P42_ARKO1_FBSM) <- "RNA"
all.genes <- rownames(P42_ARKO1_FBSM)
P42_ARKO1_FBSM <- ScaleData(P42_ARKO1_FBSM, features = all.genes)
P42_ARKO1_FBSM <- CellCycleScoring(P42_ARKO1_FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P42_ARKO1_FBSM) <- "Phase"
DimPlot(P42_ARKO1_FBSM, reduction = "umap")
tiff(file = "P42_ARKO1_FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_FBSM, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "orange","magenta2"))
dev.off()

#Cell Cycle regression
P42_ARKO1_FBSM1 <- P42_ARKO1_FBSM
DefaultAssay(P42_ARKO1_FBSM1) <- "RNA"
P42_ARKO1_FBSM1 <- ScaleData(P42_ARKO1_FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P42_ARKO1_FBSM1))
P42_ARKO1_FBSM1 <- RunPCA(P42_ARKO1_FBSM1, features = VariableFeatures(P42_ARKO1_FBSM1))
ElbowPlot(P42_ARKO1_FBSM1, ndims = 50)

P42_ARKO1_FBSM1 <- FindNeighbors(P42_ARKO1_FBSM1, reduction = "pca", dims = 1:25)
P42_ARKO1_FBSM1 <- FindClusters(P42_ARKO1_FBSM1, resolution = 0.6)
P42_ARKO1_FBSM1 <- RunUMAP(P42_ARKO1_FBSM1, reduction = "pca", dims = 1:25)

Idents(object = P42_ARKO1_FBSM1) <- "seurat_clusters"
DimPlot(P42_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_ARKO1_FBSM1) <- "seurat_clusters"
tiff(file = "P42_ARKO1_FBSM1 UMAP dims25 res0.6.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_ARKO1_FBSM1) <- "Phase"
tiff(file = "P42_ARKO1_FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "orange","magenta2"))
dev.off()

#P42_ARKO1_FBSM1 celltype identification
DefaultAssay(P42_ARKO1_FBSM1)<-"RNA"
tiff(file = "P42_ARKO1_FBSM1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1_FBSM1, reduction = "umap", features = c(
  "Fbln1", "Dcn", "Vim", "Col1a1", 
  "Myh11", "Acta2", "Actg2", "Tagln"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = P42_ARKO1_FBSM1) <- "seurat_clusters"
P42_ARKO1_FBSM1 <- RenameIdents(object = P42_ARKO1_FBSM1, 
                                '4' = "FB1", '3' = "FB2", '8' = "FB3", 
                                '2' = "FB4", '0' = "FB5",'6' = "FB6",
                                '7' = "MyoFB1", 
                                '5' = "SM1", '1' = "SM2", 
                                '9' = "OS")
P42_ARKO1_FBSM1[["FBSMCellTypes"]] <- Idents(object = P42_ARKO1_FBSM1)

#CellTypes UMAP
Idents(object = P42_ARKO1_FBSM1) <- "FBSMCellTypes"
tiff(file = "P42_ARKO1_FBSM1 UMAP FBSMCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

#Cell counts
Idents(object = P42_ARKO1_FBSM1) <- "FBSMCellTypes"
P42_ARKO1_FBSM1$EGFPExp.FBSMCellTypes <- paste(Idents(P42_ARKO1_FBSM1), P42_ARKO1_FBSM1$EGFPExp, sep = "_")
Idents(object = P42_ARKO1_FBSM1) <- "EGFPExp.FBSMCellTypes"
table(Idents(P42_ARKO1_FBSM1))

#Cell counts
Idents(object = P42_ARKO1_FBSM1) <- "FBSMCellTypes"
P42_ARKO1_FBSM1$Gli1Exp.FBSMCellTypes <- paste(Idents(P42_ARKO1_FBSM1), P42_ARKO1_FBSM1$Gli1Exp, sep = "_")
Idents(object = P42_ARKO1_FBSM1) <- "Gli1Exp.FBSMCellTypes"
table(Idents(P42_ARKO1_FBSM1))

#FeaturePlots
DefaultAssay(P42_ARKO1_FBSM1)<-"RNA"
tiff(file = "P42_ARKO1_FBSM1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1_FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_ARKO1_FBSM1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1_FBSM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "P42_ARKO1_FBSM1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_ARKO1_FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

