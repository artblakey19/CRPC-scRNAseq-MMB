#Gli1-lineage cells_E18.5_Ctrl

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
setwd("/Volumes/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/E18.5/Ctrl")

E18_Ctrl.unfiltered.data <- Read10X("/Volumes/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/CellRanger_Output/count_Project_COHP_36982_1_X3SC3/outs/filtered_feature_bc_matrix")
E18_Ctrl.unfiltered <- CreateSeuratObject(counts = E18_Ctrl.unfiltered.data,  min.cells = 3, min.features = 200, project = "E18_Ctrl")
E18_Ctrl.unfiltered <- NormalizeData(E18_Ctrl.unfiltered)

E18_Ctrl.unfiltered[["percent.mt"]] <- PercentageFeatureSet(E18_Ctrl.unfiltered, pattern = "^mt-")
tiff(file = "E18_Ctrl Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "E18_Ctrl Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_Ctrl.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E18_Ctrl Pre-filteration")
dev.off()
tiff(file = "E18_Ctrl Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_Ctrl.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E18_Ctrl Pre-filteration")
dev.off()

plot1 <- FeatureScatter(E18_Ctrl.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E18_Ctrl.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E18_Ctrl <- subset(E18_Ctrl.unfiltered, subset = nFeature_RNA > 900 & nFeature_RNA < 7400 & percent.mt < 10)
tiff(file = "E18_Ctrl Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(E18_Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "E18_Ctrl Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "skyblue", xlab = "nFeature_RNA", main = "E18_Ctrl Post-filteration")
dev.off()
tiff(file = "E18_Ctrl Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(E18_Ctrl@meta.data$percent.mt, breaks = 100, col = "skyblue", xlab = "percent.mt", main = "E18_Ctrl Post-filteration")
dev.off()

plot1 <- FeatureScatter(E18_Ctrl, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E18_Ctrl, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E18_Ctrl <- FindVariableFeatures(E18_Ctrl, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(E18_Ctrl)

#new filtering paramaters
#Clustering
DefaultAssay(E18_Ctrl) <- "RNA"
all.genes <- rownames(E18_Ctrl)
E18_Ctrl <- ScaleData(E18_Ctrl, features = all.genes)
E18_Ctrl <- RunPCA(E18_Ctrl, features = VariableFeatures(object = E18_Ctrl))
ElbowPlot(E18_Ctrl, ndims = 50)

E18_Ctrl <- FindNeighbors(E18_Ctrl, dims = 1:25)
E18_Ctrl <- FindClusters(E18_Ctrl, resolution = 0.5)
E18_Ctrl <- RunUMAP(E18_Ctrl, dims = 1:25)
DimPlot(E18_Ctrl, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_Ctrl) <- "stim"
tiff(file = "E18_Ctrl UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl, reduction = "umap", pt.size = 0.3, cols = c("grey")) 
dev.off()

#Cell cycle scoring
DefaultAssay(E18_Ctrl) <- "RNA"
all.genes <- rownames(E18_Ctrl)
E18_Ctrl <- ScaleData(E18_Ctrl, features = all.genes)
E18_Ctrl <- CellCycleScoring(E18_Ctrl, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_Ctrl) <- "Phase"
DimPlot(E18_Ctrl, reduction = "umap")
tiff(file = "E18_Ctrl Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell Cycle regression
E18_Ctrl1 <- E18_Ctrl
DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1 <- ScaleData(E18_Ctrl1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_Ctrl1))
E18_Ctrl1 <- RunPCA(E18_Ctrl1, features = VariableFeatures(E18_Ctrl1))
ElbowPlot(E18_Ctrl1, ndims = 50)

E18_Ctrl1 <- FindNeighbors(E18_Ctrl1, reduction = "pca", dims = 1:28)
E18_Ctrl1 <- FindClusters(E18_Ctrl1, resolution = 2)
E18_Ctrl1 <- RunUMAP(E18_Ctrl1, reduction = "pca", dims = 1:28)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_Ctrl1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1 UMAP dims28 res2.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_Ctrl1) <- "Phase"
tiff(file = "E18_Ctrl1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#UGS identification
DefaultAssay(E18_Ctrl1)<-"RNA"
tiff(file = "E18_Ctrl1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt14", 
                                                        "Upk1b", "Upk3b", "Pax8", 
                                                        "Vim", "Acta2", "Myh11", "Rgs5", "Plp1",
                                                        "Pecam1", "Tyrobp", "Rgs1", "Myog", "Syp"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#FeaturePlots
DefaultAssay(E18_Ctrl1)<-"RNA"
tiff(file = "E18_Ctrl1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
Idents(object = E18_Ctrl1) <- "seurat_clusters"
DefaultAssay(E18_Ctrl1) <- "RNA"
all.genes <- rownames(E18_Ctrl1)
E18_Ctrl1 <- ScaleData(E18_Ctrl1, features = all.genes)
E18_Ctrl1.seurat_cluster.markers <- FindAllMarkers(E18_Ctrl1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(E18_Ctrl1.seurat_cluster.markers, file = "E18_Ctrl1.seurat_cluster.markers.csv")

#Rename
Idents(object = E18_Ctrl1) <- "seurat_clusters"
E18_Ctrl1 <- RenameIdents(object = E18_Ctrl1, 
                          '31' = "UGE", '26' ="UGE",  '17' ="UGE",'19' ="UGE",
                          '22' = "Urothelium", '23' = "Urothelium",'29' = "Urothelium",
                          '8' = "WD",
                          '30' = "Mesenchyme",'15' = "Mesenchyme",'1' = "Mesenchyme",'9' = "Mesenchyme",
                          '4' = "Mesenchyme",'7' = "Mesenchyme",'6' = "Mesenchyme",'0' = "Mesenchyme",
                          '10' = "Mesenchyme",'13' = "Mesenchyme",'5' = "Mesenchyme",'3' = "Mesenchyme",
                          '5' = "Mesenchyme",'3' = "Mesenchyme",'2' = "Mesenchyme",'20' = "Mesenchyme",
                          '24' = "Mesenchyme",'11' = "Mesenchyme",'27' = "MyoBlast", 
                          '18' = "Pericyte",'25' = "Glia", '12' = "Glia", '16' = "VE",  
                          '14' = "NE", '21' = "NE",
                          '28'="Immune")
E18_Ctrl1[["CellTypes"]] <- Idents(object = E18_Ctrl1)

#CellTypes UMAP
Idents(object = E18_Ctrl1) <- "CellTypes"
tiff(file = "E18_Ctrl1 UMAP CellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red",  "green4", "bisque3", 
                                                                             "black"))
dev.off()

#EGFP expression
DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1EGFPPos <- subset(x=E18_Ctrl1,  subset = `EGFP` > 1.3)
E18_Ctrl1EGFPNeg <- subset(x=E18_Ctrl1,  subset = `EGFP` < 1.3)
Idents(object = E18_Ctrl1EGFPPos) <- "EGFPPos"
Idents(object = E18_Ctrl1EGFPNeg) <- "EGFPNeg"
E18_Ctrl1EGFPPos[["EGFPExp"]] <- Idents(object = E18_Ctrl1EGFPPos)
E18_Ctrl1EGFPNeg[["EGFPExp"]] <- Idents(object = E18_Ctrl1EGFPNeg)
E18_Ctrl1EGFP <- merge(x = E18_Ctrl1EGFPPos, y = E18_Ctrl1EGFPNeg)
Idents(object = E18_Ctrl1EGFP) <- "EGFPExp"
E18_Ctrl1$EGFPExp <- Idents(object = E18_Ctrl1EGFP)
Idents(object = E18_Ctrl1) <- "EGFPExp"
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Gli1 expression
DefaultAssay(E18_Ctrl1) <- "RNA"
E18_Ctrl1Gli1Pos <- subset(x=E18_Ctrl1,  subset = `Gli1` > 0)
E18_Ctrl1Gli1Neg <- subset(x=E18_Ctrl1,  subset = `Gli1` == 0)
Idents(object = E18_Ctrl1Gli1Pos) <- "Gli1Pos"
Idents(object = E18_Ctrl1Gli1Neg) <- "Gli1Neg"
E18_Ctrl1Gli1Pos[["Gli1Exp"]] <- Idents(object = E18_Ctrl1Gli1Pos)
E18_Ctrl1Gli1Neg[["Gli1Exp"]] <- Idents(object = E18_Ctrl1Gli1Neg)
E18_Ctrl1Gli1 <- merge(x = E18_Ctrl1Gli1Pos, y = E18_Ctrl1Gli1Neg)
Idents(object = E18_Ctrl1Gli1) <- "Gli1Exp"
E18_Ctrl1$Gli1Exp <- Idents(object = E18_Ctrl1Gli1)
Idents(object = E18_Ctrl1) <- "Gli1Exp"
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

#Cell counts
Idents(object = E18_Ctrl1) <- "CellTypes"
E18_Ctrl1$EGFPExp.CellTypes <- paste(Idents(E18_Ctrl1), E18_Ctrl1$EGFPExp, sep = "_")
Idents(object = E18_Ctrl1) <- "EGFPExp.CellTypes"
table(Idents(E18_Ctrl1))

#Cell counts
Idents(object = E18_Ctrl1) <- "CellTypes"
E18_Ctrl1$Gli1Exp.CellTypes <- paste(Idents(E18_Ctrl1), E18_Ctrl1$Gli1Exp, sep = "_")
Idents(object = E18_Ctrl1) <- "Gli1Exp.CellTypes"
table(Idents(E18_Ctrl1))

####Subclustering Stro####
Idents(object = E18_Ctrl1) <- "CellTypes"
E18_Ctrl1_Stro <- subset(E18_Ctrl1, idents = c("Mesenchyme", "MyoBlast", "Pericyte", "Glia", "VE", "Immune"))

DefaultAssay(E18_Ctrl1_Stro) <- "RNA"
#Run the standard workflow for visualization and clustering
E18_Ctrl1_Stro <- ScaleData(E18_Ctrl1_Stro, verbose = FALSE)
E18_Ctrl1_Stro <- RunPCA(E18_Ctrl1_Stro, npcs = 30, verbose = FALSE)
ElbowPlot(E18_Ctrl1_Stro, ndims = 50)
# UMAP and Clustering
E18_Ctrl1_Stro <- FindNeighbors(E18_Ctrl1_Stro, reduction = "pca", dims = 1:22)
E18_Ctrl1_Stro <- FindClusters(E18_Ctrl1_Stro, resolution = 0.8)
E18_Ctrl1_Stro <- RunUMAP(E18_Ctrl1_Stro, reduction = "pca", dims = 1:22)
DimPlot(E18_Ctrl1_Stro, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(E18_Ctrl1_Stro) <- "RNA"
all.genes <- rownames(E18_Ctrl1_Stro)
E18_Ctrl1_Stro <- ScaleData(E18_Ctrl1_Stro, features = all.genes)
E18_Ctrl1_Stro <- CellCycleScoring(E18_Ctrl1_Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_Ctrl1_Stro) <- "Phase"
DimPlot(E18_Ctrl1_Stro, reduction = "umap")
tiff(file = "E18_Ctrl1_Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_Stro, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen","orange", "magenta2"))
dev.off()

#Cell Cycle regression
E18_Ctrl1_Stro1 <- E18_Ctrl1_Stro
DefaultAssay(E18_Ctrl1_Stro1) <- "RNA"
E18_Ctrl1_Stro1 <- ScaleData(E18_Ctrl1_Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_Ctrl1_Stro1))
E18_Ctrl1_Stro1 <- RunPCA(E18_Ctrl1_Stro1, features = VariableFeatures(E18_Ctrl1_Stro1))
ElbowPlot(E18_Ctrl1_Stro1, ndims = 50)

E18_Ctrl1_Stro1 <- FindNeighbors(E18_Ctrl1_Stro1, reduction = "pca", dims = 1:25)
E18_Ctrl1_Stro1 <- FindClusters(E18_Ctrl1_Stro1, resolution = 1)
E18_Ctrl1_Stro1 <- RunUMAP(E18_Ctrl1_Stro1, reduction = "pca", dims = 1:25)
DimPlot(E18_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_Ctrl1_Stro1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1_Stro1 UMAP dims25 res1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_Ctrl1_Stro1) <- "Phase"
tiff(file = "E18_Ctrl1_Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#E18_Ctrl1_Stro1 celltype identification
DefaultAssay(E18_Ctrl1_Stro1)<-"RNA"
tiff(file = "E18_Ctrl1_Stro1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_Stro1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt14", 
                                                              "Upk1b", "Upk3b", "Pax8", 
                                                              "Vim", "Acta2", "Myh11", "Rgs5", "Plp1",
                                                              "Pecam1", "Tyrobp", "Crabp1", "Myog", "Syp"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = E18_Ctrl1_Stro1) <- "seurat_clusters"
E18_Ctrl1_Stro1 <- RenameIdents(object = E18_Ctrl1_Stro1, 
                                '16' = "FB", '6' = "FB",  '3' = "FB", '9' = "FB", '7' = "FB",
                                '5' = "FB", '0' = "FB",'2' = "FB",'18' = "FB",
                                '11' = "FB",'4' = "FB",'12' = "FB", '1'="WDM",
                                '17'="Myoblast",
                                '15' = "SM",'16' = "SM",'21' = "SM",'10' = "SM",
                                '14' = "Pericyte", '8' = "Glia",'19' = "Glia", '13' = "VE", 
                                '20'="Immune",
                                '22'="OS")
E18_Ctrl1_Stro1[["StroCellTypes"]] <- Idents(object = E18_Ctrl1_Stro1)

#CellTypes UMAP
Idents(object = E18_Ctrl1_Stro1) <- "StroCellTypes"
tiff(file = "E18_Ctrl1_Stro1 UMAP StroCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_Stro1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

####Subclustering FBSM####
Idents(object = E18_Ctrl1_Stro1) <- "StroCellTypes"
E18_Ctrl1_FBSM <- subset(E18_Ctrl1_Stro1, idents = c("FB", "SM"))

DefaultAssay(E18_Ctrl1_FBSM) <- "RNA"
#Run the standard workflow for visualization and clustering
E18_Ctrl1_FBSM <- ScaleData(E18_Ctrl1_FBSM, verbose = FALSE)
E18_Ctrl1_FBSM <- RunPCA(E18_Ctrl1_FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(E18_Ctrl1_FBSM, ndims = 50)
# UMAP and Clustering
E18_Ctrl1_FBSM <- FindNeighbors(E18_Ctrl1_FBSM, reduction = "pca", dims = 1:22)
E18_Ctrl1_FBSM <- FindClusters(E18_Ctrl1_FBSM, resolution = 0.5)
E18_Ctrl1_FBSM <- RunUMAP(E18_Ctrl1_FBSM, reduction = "pca", dims = 1:22)
DimPlot(E18_Ctrl1_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Cell cycle scoring
DefaultAssay(E18_Ctrl1_FBSM) <- "RNA"
all.genes <- rownames(E18_Ctrl1_FBSM)
E18_Ctrl1_FBSM <- ScaleData(E18_Ctrl1_FBSM, features = all.genes)
E18_Ctrl1_FBSM <- CellCycleScoring(E18_Ctrl1_FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = E18_Ctrl1_FBSM) <- "Phase"
DimPlot(E18_Ctrl1_FBSM, reduction = "umap")
tiff(file = "E18_Ctrl1_FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "orange","magenta2"))
dev.off()

#Cell Cycle regression
E18_Ctrl1_FBSM1 <- E18_Ctrl1_FBSM
DefaultAssay(E18_Ctrl1_FBSM1) <- "RNA"
E18_Ctrl1_FBSM1 <- ScaleData(E18_Ctrl1_FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(E18_Ctrl1_FBSM1))
E18_Ctrl1_FBSM1 <- RunPCA(E18_Ctrl1_FBSM1, features = VariableFeatures(E18_Ctrl1_FBSM1))
ElbowPlot(E18_Ctrl1_FBSM1, ndims = 50)

E18_Ctrl1_FBSM1 <- FindNeighbors(E18_Ctrl1_FBSM1, reduction = "pca", dims = 1:24)
E18_Ctrl1_FBSM1 <- FindClusters(E18_Ctrl1_FBSM1, resolution = 0.6)
E18_Ctrl1_FBSM1 <- RunUMAP(E18_Ctrl1_FBSM1, reduction = "pca", dims = 1:24)

Idents(object = E18_Ctrl1_FBSM1) <- "seurat_clusters"
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = E18_Ctrl1_FBSM1) <- "seurat_clusters"
tiff(file = "E18_Ctrl1_FBSM1 UMAP dims24 res0.6.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_Ctrl1_FBSM1) <- "Phase"
tiff(file = "E18_Ctrl1_FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen",  "orange","magenta2"))
dev.off()

#E18_Ctrl1_FBSM1 celltype identification
DefaultAssay(E18_Ctrl1_FBSM1)<-"RNA"
tiff(file = "E18_Ctrl1_FBSM1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c(
  "Fbln1", "Dcn", "Vim", "Col1a1", 
  "Myh11", "Acta2", "Actg2", "Tagln", "Thy1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = E18_Ctrl1_FBSM1) <- "seurat_clusters"
E18_Ctrl1_FBSM1 <- RenameIdents(object = E18_Ctrl1_FBSM1, 
                                '0' = "FB1", '3' = "FB2", '4' = "FB3", 
                                '2' = "FB4", '7' = "FB5",'1' = "FB6",
                                '8' = "FB7", '6' = "FB8",'10' = "FB9",
                                '9' = "SM1", '5' = "SM2" 
                                )
E18_Ctrl1_FBSM1[["FBSMCellTypes"]] <- Idents(object = E18_Ctrl1_FBSM1)

#CellTypes UMAP
Idents(object = E18_Ctrl1_FBSM1) <- "FBSMCellTypes"
tiff(file = "E18_Ctrl1_FBSM1 UMAP FBSMCellTypes.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl1_FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red", "turquoise3",
                                                                                   "slategray3", "green4", "bisque3", 
                                                                                   "black"))
dev.off()

#Cell counts
Idents(object = E18_Ctrl1_FBSM1) <- "FBSMCellTypes"
E18_Ctrl1_FBSM1$EGFPExp.FBSMCellTypes <- paste(Idents(E18_Ctrl1_FBSM1), E18_Ctrl1_FBSM1$EGFPExp, sep = "_")
Idents(object = E18_Ctrl1_FBSM1) <- "EGFPExp.FBSMCellTypes"
table(Idents(E18_Ctrl1_FBSM1))

#Cell counts
Idents(object = E18_Ctrl1_FBSM1) <- "FBSMCellTypes"
E18_Ctrl1_FBSM1$Gli1Exp.FBSMCellTypes <- paste(Idents(E18_Ctrl1_FBSM1), E18_Ctrl1_FBSM1$Gli1Exp, sep = "_")
Idents(object = E18_Ctrl1_FBSM1) <- "Gli1Exp.FBSMCellTypes"
table(Idents(E18_Ctrl1_FBSM1))

#FeaturePlots
DefaultAssay(E18_Ctrl1_FBSM1)<-"RNA"
tiff(file = "E18_Ctrl1_FBSM1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 mGFP UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "1.3", max.cutoff = "q90")
dev.off()
tiff(file = "E18_Ctrl1_FBSM1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(E18_Ctrl1_FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

