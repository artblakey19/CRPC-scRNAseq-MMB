#Gli1-lineage cells_P42_WT

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
setwd("/Volumes/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/P42/Ctrl")

P42_Ctrl.unfiltered.data <- Read10X("/Volumes/user_data-1/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/CellRanger_Output/count_Project_COHP_36980_1_X3SC3/outs/filtered_feature_bc_matrix")
P42_Ctrl.unfiltered <- CreateSeuratObject(counts = P42_Ctrl.unfiltered.data,  min.cells = 3, min.features = 200, project = "P42_Ctrl")
P42_Ctrl.unfiltered <- NormalizeData(P42_Ctrl.unfiltered)

#P42_Ctrl
#Initial processing & filtering
P42_Ctrl.unfiltered[["percent.mt"]] <- PercentageFeatureSet(P42_Ctrl.unfiltered, pattern = "^mt-")
tiff(file = "P42_Ctrl Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(P42_Ctrl.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "P42_Ctrl Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_Ctrl.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P42_Ctrl Pre-filteration")
dev.off()
tiff(file = "P42_Ctrl Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_Ctrl.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P42_Ctrl Pre-filteration")
dev.off()

plot1 <- FeatureScatter(P42_Ctrl.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P42_Ctrl.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P42_Ctrl <- subset(P42_Ctrl.unfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 6100 & percent.mt < 15)
tiff(file = "P42_Ctrl Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(P42_Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "P42_Ctrl Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "skyblue", xlab = "nFeature_RNA", main = "P42_Ctrl Post-filteration")
dev.off()
tiff(file = "P42_Ctrl Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(P42_Ctrl@meta.data$percent.mt, breaks = 100, col = "skyblue", xlab = "percent.mt", main = "P42_Ctrl Post-filteration")
dev.off()

plot1 <- FeatureScatter(P42_Ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P42_Ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P42_Ctrl <- FindVariableFeatures(P42_Ctrl, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(P42_Ctrl)

#new filtering paramaters
#Clustering
DefaultAssay(P42_Ctrl) <- "RNA"
all.genes <- rownames(P42_Ctrl)
P42_Ctrl <- ScaleData(P42_Ctrl, features = all.genes)
P42_Ctrl <- RunPCA(P42_Ctrl, features = VariableFeatures(object = P42_Ctrl))
ElbowPlot(P42_Ctrl, ndims = 50)

P42_Ctrl <- FindNeighbors(P42_Ctrl, dims = 1:28)
P42_Ctrl <- FindClusters(P42_Ctrl, resolution = 0.5)
P42_Ctrl <- RunUMAP(P42_Ctrl, dims = 1:28)
DimPlot(P42_Ctrl, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle scoring
DefaultAssay(P42_Ctrl) <- "RNA"
all.genes <- rownames(P42_Ctrl)
P42_Ctrl <- ScaleData(P42_Ctrl, features = all.genes)
P42_Ctrl <- CellCycleScoring(P42_Ctrl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P42_Ctrl) <- "Phase"
DimPlot(P42_Ctrl, reduction = "umap")
tiff(file = "P42_Ctrl Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen","orange","magenta2"))
dev.off()

#Cell Cycle regression
P42_Ctrl1 <- P42_Ctrl
DefaultAssay(P42_Ctrl1) <- "RNA"
P42_Ctrl1 <- ScaleData(P42_Ctrl1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P42_Ctrl1))
P42_Ctrl1 <- RunPCA(P42_Ctrl1, features = VariableFeatures(P42_Ctrl1))
ElbowPlot(P42_Ctrl1, ndims = 50)

P42_Ctrl1 <- FindNeighbors(P42_Ctrl1, reduction = "pca", dims = 1:27)
P42_Ctrl1 <- FindClusters(P42_Ctrl1, resolution = 1)
P42_Ctrl1 <- RunUMAP(P42_Ctrl1, reduction = "pca", dims = 1:27)
DimPlot(P42_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_Ctrl1) <- "seurat_clusters"
tiff(file = "P42_Ctrl1 UMAP dims27 res2.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_Ctrl1) <- "Phase"
tiff(file = "P42_Ctrl1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#P11 celltype identification
DefaultAssay(P42_Ctrl1)<-"RNA"
tiff(file = "P42_Ctrl1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Cdo1", 
                                                        "Vim",  "Fbln1", "Myh11", "Acta2", "Tagln", "Actg2",
                                                        "Pecam1", "Tyrobp", "Plp1", "Syp", "Rgs5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#FeaturePlots
DefaultAssay(P42_Ctrl1)<-"RNA"
tiff(file = "P42_Ctrl1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_Ctrl1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "P42_Ctrl1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
Idents(object = P42_Ctrl1) <- "seurat_clusters"
DefaultAssay(P42_Ctrl1) <- "RNA"
all.genes <- rownames(P42_Ctrl1)
P42_Ctrl1 <- ScaleData(P42_Ctrl1, features = all.genes)
P42_Ctrl1.seurat_cluster.markers <- FindAllMarkers(P42_Ctrl1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(P42_Ctrl1.seurat_cluster.markers, file = "P42_Ctrl1.seurat_cluster.markers.csv")


#Rename
Idents(object = P42_Ctrl1) <- "seurat_clusters"
P42_Ctrl1 <- RenameIdents(object = P42_Ctrl1, 
                          '24' = "BE", '0' ="BE",  '4' ="BE",
                          '3' = "LE", '2' = "LE",'19' = "LE",'22' = "LE",
                          '6' = "LE",'7' = "LE",'17' = "LE",'9' = "LE",
                          '18' = "MyoE",
                          '27' = "SV",'14' = "SV",'25' = "SV",
                          '8' = "FB",'1' = "FB",'20' = "FB",'15' = "FB",'5' = "FB",
                          '12' = "SM", 
                          '16' = "Pericyte", '21' = "Glia",  '10' = "VE", '28' = "VE", 
                          '13'="Immune", '11'="Immune", '23'="Immune", '26'="Immune")
P42_Ctrl1[["CellTypes"]] <- Idents(object = P42_Ctrl1)

#CellTypes UMAP
Idents(object = P42_Ctrl1) <- "CellTypes"
tiff(file = "P42_Ctrl1 UMAP CellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red",  "green4", "bisque3", 
                                                                             "black"))
dev.off()

#EGFP expression
DefaultAssay(P42_Ctrl1) <- "RNA"
P42_Ctrl1EGFPPos <- subset(x=P42_Ctrl1,  subset = `EGFP` > 1.4)
P42_Ctrl1EGFPNeg <- subset(x=P42_Ctrl1,  subset = `EGFP` < 1.4)
Idents(object = P42_Ctrl1EGFPPos) <- "EGFPPos"
Idents(object = P42_Ctrl1EGFPNeg) <- "EGFPNeg"
P42_Ctrl1EGFPPos[["EGFPExp"]] <- Idents(object = P42_Ctrl1EGFPPos)
P42_Ctrl1EGFPNeg[["EGFPExp"]] <- Idents(object = P42_Ctrl1EGFPNeg)
P42_Ctrl1EGFP <- merge(x = P42_Ctrl1EGFPPos, y = P42_Ctrl1EGFPNeg)
Idents(object = P42_Ctrl1EGFP) <- "EGFPExp"
P42_Ctrl1$EGFPExp <- Idents(object = P42_Ctrl1EGFP)
Idents(object = P42_Ctrl1) <- "EGFPExp"
DimPlot(P42_Ctrl1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Gli1 expression
DefaultAssay(P42_Ctrl1) <- "RNA"
P42_Ctrl1Gli1Pos <- subset(x=P42_Ctrl1,  subset = `Gli1` > 0)
P42_Ctrl1Gli1Neg <- subset(x=P42_Ctrl1,  subset = `Gli1` == 0)
Idents(object = P42_Ctrl1Gli1Pos) <- "Gli1Pos"
Idents(object = P42_Ctrl1Gli1Neg) <- "Gli1Neg"
P42_Ctrl1Gli1Pos[["Gli1Exp"]] <- Idents(object = P42_Ctrl1Gli1Pos)
P42_Ctrl1Gli1Neg[["Gli1Exp"]] <- Idents(object = P42_Ctrl1Gli1Neg)
P42_Ctrl1Gli1 <- merge(x = P42_Ctrl1Gli1Pos, y = P42_Ctrl1Gli1Neg)
Idents(object = P42_Ctrl1Gli1) <- "Gli1Exp"
P42_Ctrl1$Gli1Exp <- Idents(object = P42_Ctrl1Gli1)
Idents(object = P42_Ctrl1) <- "Gli1Exp"
DimPlot(P42_Ctrl1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Cell counts
Idents(object = P42_Ctrl1) <- "CellTypes"
P42_Ctrl1$EGFPExp.CellTypes <- paste(Idents(P42_Ctrl1), P42_Ctrl1$EGFPExp, sep = "_")
Idents(object = P42_Ctrl1) <- "EGFPExp.CellTypes"
table(Idents(P42_Ctrl1))

#Cell counts
Idents(object = P42_Ctrl1) <- "CellTypes"
P42_Ctrl1$Gli1Exp.CellTypes <- paste(Idents(P42_Ctrl1), P42_Ctrl1$Gli1Exp, sep = "_")
Idents(object = P42_Ctrl1) <- "Gli1Exp.CellTypes"
table(Idents(P42_Ctrl1))

####Subclustering Stro####
Idents(object = P42_Ctrl1) <- "CellTypes"
P42_Ctrl1_Stro <- subset(P42_Ctrl1, idents = c("FB", "SM", "Pericyte", "Glia", "VE", "Immune"))

DefaultAssay(P42_Ctrl1_Stro) <- "RNA"
#Run the standard workflow for visualization and clustering
P42_Ctrl1_Stro <- ScaleData(P42_Ctrl1_Stro, verbose = FALSE)
P42_Ctrl1_Stro <- RunPCA(P42_Ctrl1_Stro, npcs = 30, verbose = FALSE)
ElbowPlot(P42_Ctrl1_Stro, ndims = 50)
# UMAP and Clustering
P42_Ctrl1_Stro <- FindNeighbors(P42_Ctrl1_Stro, reduction = "pca", dims = 1:21)
P42_Ctrl1_Stro <- FindClusters(P42_Ctrl1_Stro, resolution = 0.8)
P42_Ctrl1_Stro <- RunUMAP(P42_Ctrl1_Stro, reduction = "pca", dims = 1:21)
DimPlot(P42_Ctrl1_Stro, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(P42_Ctrl1_Stro) <- "RNA"
all.genes <- rownames(P42_Ctrl1_Stro)
P42_Ctrl1_Stro <- ScaleData(P42_Ctrl1_Stro, features = all.genes)
P42_Ctrl1_Stro <- CellCycleScoring(P42_Ctrl1_Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P42_Ctrl1_Stro) <- "Phase"
DimPlot(P42_Ctrl1_Stro, reduction = "umap")
tiff(file = "P42_Ctrl1_Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_Stro, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "magenta2","orange"))
dev.off()

#Cell Cycle regression
P42_Ctrl1_Stro1 <- P42_Ctrl1_Stro
DefaultAssay(P42_Ctrl1_Stro1) <- "RNA"
P42_Ctrl1_Stro1 <- ScaleData(P42_Ctrl1_Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P42_Ctrl1_Stro1))
P42_Ctrl1_Stro1 <- RunPCA(P42_Ctrl1_Stro1, features = VariableFeatures(P42_Ctrl1_Stro1))
ElbowPlot(P42_Ctrl1_Stro1, ndims = 50)

P42_Ctrl1_Stro1 <- FindNeighbors(P42_Ctrl1_Stro1, reduction = "pca", dims = 1:27)
P42_Ctrl1_Stro1 <- FindClusters(P42_Ctrl1_Stro1, resolution = 1.5)
P42_Ctrl1_Stro1 <- RunUMAP(P42_Ctrl1_Stro1, reduction = "pca", dims = 1:27)
DimPlot(P42_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_Ctrl1_Stro1) <- "seurat_clusters"
tiff(file = "P42_Ctrl1_Stro1 UMAP dims27 res1.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_Ctrl1_Stro1) <- "Phase"
tiff(file = "P42_Ctrl1_Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "magenta2", "orange"))
dev.off()

#P42_Ctrl1_Stro1 celltype identification
DefaultAssay(P42_Ctrl1_Stro1)<-"RNA"
tiff(file = "P42_Ctrl1_Stro1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1_Stro1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Cdo1", 
                                                              "Vim",  "Fbln1", "Myh11", "Acta2", "Tagln", "Actg2",
                                                              "Pecam1", "Tyrobp", "Plp1", "Syp", "Rgs5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = P42_Ctrl1_Stro1) <- "seurat_clusters"
P42_Ctrl1_Stro1 <- RenameIdents(object = P42_Ctrl1_Stro1, 
                                '0' = "FB", '5' = "FB",  '3' = "FB", '19' = "FB", '9' = "FB",
                                '1' = "FB", '2' = "FB", 
                                '10' = "SM",'15' = "SM",
                                '12' = "Pericyte", '17' = "Pericyte",'13' = "Glia", '7' = "VE", 
                                '22' = "VE",'16' = "VE",
                                '21'="Immune",'6'="Immune",'4'="Immune",'20'="Immune",
                                '8'="OS",'18'="OS",'14'="OS",'11'="OS")
P42_Ctrl1_Stro1[["StroCellTypes"]] <- Idents(object = P42_Ctrl1_Stro1)

#CellTypes UMAP
Idents(object = P42_Ctrl1_Stro1) <- "StroCellTypes"
tiff(file = "P42_Ctrl1_Stro1 UMAP StroCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

####Subclustering FBSM####
Idents(object = P42_Ctrl1_Stro1) <- "StroCellTypes"
P42_Ctrl1_FBSM <- subset(P42_Ctrl1_Stro1, idents = c("FB", "SM"))

DefaultAssay(P42_Ctrl1_FBSM) <- "RNA"
#Run the standard workflow for visualization and clustering
P42_Ctrl1_FBSM <- ScaleData(P42_Ctrl1_FBSM, verbose = FALSE)
P42_Ctrl1_FBSM <- RunPCA(P42_Ctrl1_FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(P42_Ctrl1_FBSM, ndims = 50)
# UMAP and Clustering
P42_Ctrl1_FBSM <- FindNeighbors(P42_Ctrl1_FBSM, reduction = "pca", dims = 1:20)
P42_Ctrl1_FBSM <- FindClusters(P42_Ctrl1_FBSM, resolution = 0.8)
P42_Ctrl1_FBSM <- RunUMAP(P42_Ctrl1_FBSM, reduction = "pca", dims = 1:20)
DimPlot(P42_Ctrl1_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(P42_Ctrl1_FBSM) <- "RNA"
all.genes <- rownames(P42_Ctrl1_FBSM)
P42_Ctrl1_FBSM <- ScaleData(P42_Ctrl1_FBSM, features = all.genes)
P42_Ctrl1_FBSM <- CellCycleScoring(P42_Ctrl1_FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P42_Ctrl1_FBSM) <- "Phase"
DimPlot(P42_Ctrl1_FBSM, reduction = "umap")
tiff(file = "P42_Ctrl1_FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_FBSM, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "magenta2","orange"))
dev.off()

#Cell Cycle regression
P42_Ctrl1_FBSM1 <- P42_Ctrl1_FBSM
DefaultAssay(P42_Ctrl1_FBSM1) <- "RNA"
P42_Ctrl1_FBSM1 <- ScaleData(P42_Ctrl1_FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P42_Ctrl1_FBSM1))
P42_Ctrl1_FBSM1 <- RunPCA(P42_Ctrl1_FBSM1, features = VariableFeatures(P42_Ctrl1_FBSM1))
ElbowPlot(P42_Ctrl1_FBSM1, ndims = 50)

P42_Ctrl1_FBSM1 <- FindNeighbors(P42_Ctrl1_FBSM1, reduction = "pca", dims = 1:20)
P42_Ctrl1_FBSM1 <- FindClusters(P42_Ctrl1_FBSM1, resolution = 0.8)
P42_Ctrl1_FBSM1 <- RunUMAP(P42_Ctrl1_FBSM1, reduction = "pca", dims = 1:20)

Idents(object = P42_Ctrl1_FBSM1) <- "seurat_clusters"
DimPlot(P42_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_Ctrl1_FBSM1) <- "seurat_clusters"
tiff(file = "P42_Ctrl1_FBSM1 UMAP dims20 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_Ctrl1_FBSM1) <- "Phase"
tiff(file = "P42_Ctrl1_FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "magenta2","orange"))
dev.off()

#P42_Ctrl1_FBSM1 celltype identification
DefaultAssay(P42_Ctrl1_FBSM1)<-"RNA"
tiff(file = "P42_Ctrl1_FBSM1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1_FBSM1, reduction = "umap", features = c(
  "Fbln1", "Dcn", "Vim", "Col1a1", 
  "Myh11", "Acta2", "Actg2", "Tagln"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = P42_Ctrl1_FBSM1) <- "seurat_clusters"
P42_Ctrl1_FBSM1 <- RenameIdents(object = P42_Ctrl1_FBSM1, 
                                '2' = "FB1", '1' = "FB2", '3' = "FB3", 
                                '4' = "FB4", '5' = "FB4",'0' = "FB5",
                                '6' = "SM1", '7' = "SM2", '8' = "SM2")
P42_Ctrl1_FBSM1[["FBSMCellTypes"]] <- Idents(object = P42_Ctrl1_FBSM1)

#CellTypes UMAP
Idents(object = P42_Ctrl1_FBSM1) <- "FBSMCellTypes"
tiff(file = "P42_Ctrl1_FBSM1 UMAP FBSMCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2",   "red",  "blue", "mediumorchid3",  "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

#Cell counts
Idents(object = P42_Ctrl1_FBSM1) <- "FBSMCellTypes"
P42_Ctrl1_FBSM1$EGFPExp.FBSMCellTypes <- paste(Idents(P42_Ctrl1_FBSM1), P42_Ctrl1_FBSM1$EGFPExp, sep = "_")
Idents(object = P42_Ctrl1_FBSM1) <- "EGFPExp.FBSMCellTypes"
table(Idents(P42_Ctrl1_FBSM1))

#Cell counts
Idents(object = P42_Ctrl1_FBSM1) <- "FBSMCellTypes"
P42_Ctrl1_FBSM1$Gli1Exp.FBSMCellTypes <- paste(Idents(P42_Ctrl1_FBSM1), P42_Ctrl1_FBSM1$Gli1Exp, sep = "_")
Idents(object = P42_Ctrl1_FBSM1) <- "Gli1Exp.FBSMCellTypes"
table(Idents(P42_Ctrl1_FBSM1))

#FeaturePlots
DefaultAssay(P42_Ctrl1_FBSM1)<-"RNA"
tiff(file = "P42_Ctrl1_FBSM1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1_FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_Ctrl1_FBSM1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1_FBSM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "P42_Ctrl1_FBSM1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(P42_Ctrl1_FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

