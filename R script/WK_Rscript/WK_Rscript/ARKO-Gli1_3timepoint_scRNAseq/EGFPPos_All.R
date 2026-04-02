#EGFP_All

#### Add necessary tools to library ####


devtools::install_github('satijalab/seurat-data')
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_version('speedglm', '0.3-4', repos = 'https://packagemanager.rstudio.com/cran/2023-03-31')
devtools::install_github('cole-trapnell-lab/monocle3')

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

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All")

#Individual UMAP
Idents(object = E18_Ctrl) <- "stim"
tiff(file = "E18_Ctrl UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_Ctrl, reduction = "umap", pt.size = 0.3, cols = c("salmon")) 
dev.off()
table(Idents(E18_Ctrl))

Idents(object = P11_Ctrl) <- "stim"
tiff(file = "P11_Ctrl UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P11_Ctrl, reduction = "umap", pt.size = 0.3, cols = c("green3")) 
dev.off()
table(Idents(P11_Ctrl))

Idents(object = P42_Ctrl) <- "stim"
tiff(file = "P42_Ctrl UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_Ctrl, reduction = "umap", pt.size = 0.3, cols = c("deeppink1")) 
dev.off()
table(Idents(P42_Ctrl))

Idents(object = E18_ARKO) <- "stim"
tiff(file = "E18_ARKO UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(E18_ARKO, reduction = "umap", pt.size = 0.3, cols = c("skyblue1")) 
dev.off()
table(Idents(E18_ARKO))

Idents(object = P11_ARKO) <- "stim"
tiff(file = "P11_ARKO UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P11_ARKO, reduction = "umap", pt.size = 0.3, cols = c("yellow2")) 
dev.off()
table(Idents(P11_ARKO))

Idents(object = P42_ARKO) <- "stim"
tiff(file = "P42_ARKO UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(P42_ARKO, reduction = "umap", pt.size = 0.3, cols = c("blue")) 
dev.off()
table(Idents(P42_ARKO))

#Stash old idents
Idents(object = E18_Ctrl) <- "seurat_clusters"
Idents(object = P11_Ctrl) <- "seurat_clusters"
Idents(object = P42_Ctrl) <- "seurat_clusters"
Idents(object = E18_ARKO) <- "seurat_clusters"
Idents(object = P11_ARKO) <- "seurat_clusters"
Idents(object = P42_ARKO) <- "seurat_clusters"
E18_Ctrl[["orig.clusters"]] <- Idents(object = E18_Ctrl)
P11_Ctrl[["orig.clusters"]] <- Idents(object = P11_Ctrl)
P42_Ctrl[["orig.clusters"]] <- Idents(object = P42_Ctrl)
E18_ARKO[["orig.clusters"]] <- Idents(object = E18_ARKO)
P11_ARKO[["orig.clusters"]] <- Idents(object = P11_ARKO)
P42_ARKO[["orig.clusters"]] <- Idents(object = P42_ARKO)

#Set Current idents
Idents(object = E18_Ctrl) <- "seurat_clusters"
Idents(object = P11_Ctrl) <- "seurat_clusters"
Idents(object = P42_Ctrl) <- "seurat_clusters"
Idents(object = E18_ARKO) <- "seurat_clusters"
Idents(object = P11_ARKO) <- "seurat_clusters"
Idents(object = P42_ARKO) <- "seurat_clusters"
E18_Ctrl$stim <- "E18_Ctrl"
P11_Ctrl$stim <- "P11_Ctrl"
P42_Ctrl$stim <- "P42_Ctrl"
E18_ARKO$stim <- "E18_ARKO"
P11_ARKO$stim <- "P11_ARKO"
P42_ARKO$stim <- "P42_ARKO"
features <- SelectIntegrationFeatures(object.list = list(E18_Ctrl, P11_Ctrl, P42_Ctrl, E18_ARKO, P11_ARKO, P42_ARKO))
All.combined.anchors <- FindIntegrationAnchors(object.list = list(E18_Ctrl, P11_Ctrl, P42_Ctrl, E18_ARKO, P11_ARKO,P42_ARKO), anchor.features = features, reduction = "rpca", k.anchor = 5)
All.combined <- IntegrateData(anchorset = All.combined.anchors)

DefaultAssay(All.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
All.combined <- ScaleData(All.combined, verbose = FALSE)
All.combined <- RunPCA(All.combined, npcs = 30, verbose = FALSE)
ElbowPlot(All.combined, ndims = 50)
# umap and Clustering
All.combined <- FindNeighbors(All.combined, reduction = "pca", dims = 1:19)
All.combined <- FindClusters(All.combined, resolution = 0.5)
All.combined <- RunUMAP(All.combined, reduction = "pca", dims = 1:19)
DimPlot(All.combined, reduction = "umap", pt.size = 0.3) 

Idents(object = All.combined) <- "stim"
All.combined <- RenameIdents(object = All.combined, 'E18_Ctrl' = "E18_Ctrl",'E18_ARKO' = "E18_ARKO", 
                             'P11_Ctrl' = "P11_Ctrl",'P11_ARKO' = "P11_ARKO",
                             'P42_Ctrl' = "P42_Ctrl",'P42_ARKO' = "P42_ARKO")
All.combined[["stim"]] <- Idents(object = All.combined)
DimPlot(All.combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c( "salmon", "skyblue1", "green3", "yellow2", "deeppink1", "blue")) 

tiff(file = "All.combined stim UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c( "salmon", "skyblue1", "green3", "yellow2", "deeppink1", "blue")) 
dev.off()

DefaultAssay(All.combined) <- "RNA"
VlnPlot(All.combined, features = "EGFP")

####Subclustering EGFPPos cells####
#EGFP expression
DefaultAssay(All.combined) <- "RNA"
All.combinedEGFPPos <- subset(x=All.combined,  subset = `EGFP` > 2)
All.combinedEGFPNeg <- subset(x=All.combined,  subset = `EGFP` < 2)
Idents(object = All.combinedEGFPPos) <- "EGFPPos"
Idents(object = All.combinedEGFPNeg) <- "EGFPNeg"
All.combinedEGFPPos[["EGFPExp"]] <- Idents(object = All.combinedEGFPPos)
All.combinedEGFPNeg[["EGFPExp"]] <- Idents(object = All.combinedEGFPNeg)
All.combinedEGFP <- merge(x = All.combinedEGFPPos, y = All.combinedEGFPNeg)
Idents(object = All.combinedEGFP) <- "EGFPExp"
All.combined$EGFPExp <- Idents(object = All.combinedEGFP)
Idents(object = All.combined) <- "EGFPExp"
DimPlot(All.combined, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

tiff(file = "All.combined EGFPExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))
dev.off()

save.image("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All/All.combined.RData")

#Subclustering EGFPPos cells
Idents(object = All.combined) <- "EGFPExp"
All.combined.EGFPPos <- subset(All.combined, idents = c("EGFPPos"))

DefaultAssay(All.combined.EGFPPos) <- "integrated"
#Run the standard workflow for visualization and clustering
All.combined.EGFPPos <- ScaleData(All.combined.EGFPPos, verbose = FALSE)
All.combined.EGFPPos <- RunPCA(All.combined.EGFPPos, npcs = 30, verbose = FALSE)
ElbowPlot(All.combined.EGFPPos, ndims = 50)
# UMAP and Clustering
All.combined.EGFPPos <- FindNeighbors(All.combined.EGFPPos, reduction = "pca", dims = 1:20)
All.combined.EGFPPos <- FindClusters(All.combined.EGFPPos, resolution = 0.8)
All.combined.EGFPPos <- RunUMAP(All.combined.EGFPPos, reduction = "pca", dims = 1:20)

Idents(object = All.combined.EGFPPos) <- "seurat_clusters"
DimPlot(All.combined.EGFPPos, reduction = "umap", pt.size = 0.3, label = TRUE) 
DimPlot(All.combined.EGFPPos, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim") 

#Cell Cycle scoring
DefaultAssay(All.combined.EGFPPos) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos)
All.combined.EGFPPos <- ScaleData(All.combined.EGFPPos, features = all.genes)
All.combined.EGFPPos <- CellCycleScoring(All.combined.EGFPPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = All.combined.EGFPPos) <- "Phase"
All.combined.EGFPPos <- RenameIdents(object = All.combined.EGFPPos, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
All.combined.EGFPPos[["Phase"]] <- Idents(object = All.combined.EGFPPos)
DimPlot(All.combined.EGFPPos, reduction = "umap")

tiff(file = "All.combined.EGFPPos Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos, reduction = "umap", pt.size = 0.5, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#Cell Cycle regression
All.combined.EGFPPos1 <- All.combined.EGFPPos
DefaultAssay(All.combined.EGFPPos1) <- "integrated"
All.combined.EGFPPos1 <- ScaleData(All.combined.EGFPPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(All.combined.EGFPPos1))
All.combined.EGFPPos1 <- RunPCA(All.combined.EGFPPos1, features = VariableFeatures(All.combined.EGFPPos1))
ElbowPlot(All.combined.EGFPPos1, ndims = 50)

All.combined.EGFPPos1 <- FindNeighbors(All.combined.EGFPPos1, reduction = "pca", dims = 1:32)
All.combined.EGFPPos1 <- FindClusters(All.combined.EGFPPos1, resolution = 2.5)
All.combined.EGFPPos1 <- RunUMAP(All.combined.EGFPPos1, reduction = "pca", dims = 1:34)
DimPlot(All.combined.EGFPPos1, reduction = "umap", pt.size = 0.3, label = TRUE)

DimPlot(All.combined.EGFPPos1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = All.combined.EGFPPos1) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1 UMAP dims34 res2.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = All.combined.EGFPPos1) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1 UMAP dims34 res2.5 stim.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
dev.off()

Idents(object = All.combined.EGFPPos1) <- "Phase"
tiff(file = "All.combined.EGFPPos1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1, reduction = "umap", pt.size = 0.3, cols = c("orange","lightseagreen", "magenta2"))
dev.off()

#Celltype identification
DefaultAssay(All.combined.EGFPPos1)<-"RNA"
tiff(file = "All.combined.EGFPPos1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Upk3b", "Syp", 
                                                                   "Vim",  "Fbln1", "Myh11", "Acta2", 
                                                                   "Pecam1", "Plp1", "Rgs5","Tyrobp", "Myod1","Crabp1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
Idents(object = All.combined.EGFPPos1) <- "seurat_clusters"
DefaultAssay(All.combined.EGFPPos1) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1)
All.combined.EGFPPos1 <- ScaleData(All.combined.EGFPPos1, features = all.genes)
All.combined.EGFPPos1.seurat_cluster.markers <- FindAllMarkers(All.combined.EGFPPos1, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(All.combined.EGFPPos1.seurat_cluster.markers, file = "All.combined.EGFPPos1.seurat_cluster.markers.csv")

#FeaturePlots
DefaultAssay(All.combined.EGFPPos1)<-"RNA"
tiff(file = "All.combined.EGFPPos1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1 Ar split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1 Gli1 split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = All.combined.EGFPPos1) <- "seurat_clusters"
All.combined.EGFPPos1 <- RenameIdents(object = All.combined.EGFPPos1, 
                          '8' = "FB", '16' ="FB", '5' ="FB",'13' ="FB",'0' ="FB",
                          '33' ="FB",'12' ="FB",'27' ="FB",'4' ="FB",'2' ="FB",'11' ="FB",
                          '31' ="FB",'7' ="FB",'3' ="FB",'18' ="FB",'23' ="FB",'9' ="FB",'15' ="FB",
                          '10' = "SM", '22' = "SM", '19' = "SM",'6' = "SM",
                          '14' = "MyoFB",'1' = "MyoFB",
                          '21' = "MyoE",'24' = "MyoE",'17' = "MyoE",'26' = "MyoE",
                          '36' = "Myoblast", 
                          '20' = "Pericyte",'34' = "Pericyte", '30' = "Glia",  '29' = "VE", '28' = "Bladder", '25'="SVM", 
                          '32'="Immune", '35' = "NE")
All.combined.EGFPPos1[["CellTypes"]] <- Idents(object = All.combined.EGFPPos1)

#CellTypes UMAP
Idents(object = All.combined.EGFPPos1) <- "CellTypes"
tiff(file = "All.combined.EGFPPos1 UMAP CellTypes.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "lightgrey", "olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red",  "green4", "bisque3", 
                                                                             "black"))
dev.off()
tiff(file = "All.combined.EGFPPos1 UMAP CellTypes split.tiff", width = 24, height =4, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "lightgrey","olivedrab2", "yellow2", "deeppink1",  "blue", "mediumorchid3",  "red",  "green4", "bisque3", 
                                                                                        "black"))
dev.off()

#Cell counts
Idents(object = All.combined.EGFPPos1) <- "CellTypes"
All.combined.EGFPPos1$stim.CellTypes <- paste(Idents(All.combined.EGFPPos1), All.combined.EGFPPos1$stim, sep = "_")
Idents(object = All.combined.EGFPPos1) <- "stim.CellTypes"
table(Idents(All.combined.EGFPPos1))

save.image("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All/All.combined.EGFPPos1.RData")

####Subclustering FBSM####
Idents(object = All.combined.EGFPPos1) <- "CellTypes"
All.combined.EGFPPos1.FBSM <- subset(All.combined.EGFPPos1, idents = c("FB", "SM", "MyoFB", "MyoE"))

DefaultAssay(All.combined.EGFPPos1.FBSM) <- "integrated"
#Run the standard workflow for visualization and clustering
All.combined.EGFPPos1.FBSM <- ScaleData(All.combined.EGFPPos1.FBSM, verbose = FALSE)
All.combined.EGFPPos1.FBSM <- RunPCA(All.combined.EGFPPos1.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(All.combined.EGFPPos1.FBSM, ndims = 50)
# UMAP and Clustering
All.combined.EGFPPos1.FBSM <- FindNeighbors(All.combined.EGFPPos1.FBSM, reduction = "pca", dims = 1:20)
All.combined.EGFPPos1.FBSM <- FindClusters(All.combined.EGFPPos1.FBSM, resolution = 0.8)
All.combined.EGFPPos1.FBSM <- RunUMAP(All.combined.EGFPPos1.FBSM, reduction = "pca", dims = 1:20)

Idents(object = All.combined.EGFPPos1.FBSM) <- "seurat_clusters"
DimPlot(All.combined.EGFPPos1.FBSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
DimPlot(All.combined.EGFPPos1.FBSM, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim") 

#Cell Cycle scoring
DefaultAssay(All.combined.EGFPPos1.FBSM) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM)
All.combined.EGFPPos1.FBSM <- ScaleData(All.combined.EGFPPos1.FBSM, features = all.genes)
All.combined.EGFPPos1.FBSM <- CellCycleScoring(All.combined.EGFPPos1.FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = All.combined.EGFPPos1.FBSM) <- "Phase"
All.combined.EGFPPos1.FBSM <- RenameIdents(object = All.combined.EGFPPos1.FBSM, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
All.combined.EGFPPos1.FBSM[["Phase"]] <- Idents(object = All.combined.EGFPPos1.FBSM)
DimPlot(All.combined.EGFPPos1.FBSM, reduction = "umap")

tiff(file = "All.combined.EGFPPos1.FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM, reduction = "umap", pt.size = 0.5, cols = c("orange", "lightseagreen", "magenta2"))
dev.off()

#Cell Cycle regression
All.combined.EGFPPos1.FBSM1 <- All.combined.EGFPPos1.FBSM
DefaultAssay(All.combined.EGFPPos1.FBSM1) <- "integrated"
All.combined.EGFPPos1.FBSM1 <- ScaleData(All.combined.EGFPPos1.FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(All.combined.EGFPPos1.FBSM1))
All.combined.EGFPPos1.FBSM1 <- RunPCA(All.combined.EGFPPos1.FBSM1, features = VariableFeatures(All.combined.EGFPPos1.FBSM1))
ElbowPlot(All.combined.EGFPPos1.FBSM1, ndims = 50)

All.combined.EGFPPos1.FBSM1 <- FindNeighbors(All.combined.EGFPPos1.FBSM1, reduction = "pca", dims = 1:25)
All.combined.EGFPPos1.FBSM1 <- FindClusters(All.combined.EGFPPos1.FBSM1, resolution = 2.5)
All.combined.EGFPPos1.FBSM1 <- RunUMAP(All.combined.EGFPPos1.FBSM1, reduction = "pca", dims = 1:25)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

#Optimize resolution
DefaultAssay(All.combined.EGFPPos1.FBSM1) <- "integrated"
All.combined.EGFPPos1.FBSM1 <- ScaleData(All.combined.EGFPPos1.FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(All.combined.EGFPPos1.FBSM1))
All.combined.EGFPPos1.FBSM1 <- RunPCA(All.combined.EGFPPos1.FBSM1, features = VariableFeatures(All.combined.EGFPPos1.FBSM1))
ElbowPlot(All.combined.EGFPPos1.FBSM1, ndims = 50)

All.combined.EGFPPos1.FBSM1 <- FindNeighbors(All.combined.EGFPPos1.FBSM1, reduction = "pca", dims = 1:25)
All.combined.EGFPPos1.FBSM1 <- FindClusters(All.combined.EGFPPos1.FBSM1, resolution = 0.8)
All.combined.EGFPPos1.FBSM1 <- RunUMAP(All.combined.EGFPPos1.FBSM1, reduction = "pca", dims = 1:25)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = All.combined.EGFPPos1.FBSM1) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM1 UMAP dims25 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
Idents(object = All.combined.EGFPPos1.FBSM1) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM1 split UMAP dims25 res0.8.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
dev.off()

Idents(object = All.combined.EGFPPos1.FBSM1) <- "Phase"
tiff(file = "All.combined.EGFPPos1.FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, cols = c("orange","lightseagreen", "magenta2"))
dev.off()

#Celltype identification
DefaultAssay(All.combined.EGFPPos1.FBSM1)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Upk3b", "Syp", 
                                                                   "Vim",  "Fbln1", "Myh11", "Acta2", 
                                                                   "Pecam1", "Plp1", "Rgs5","Tyrobp", "Myod1","Crabp1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#FeaturePlots
DefaultAssay(All.combined.EGFPPos1.FBSM1)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM1 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM1 Ar split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM1 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM1 Gli1 split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
DefaultAssay(All.combined.EGFPPos1.FBSM1)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM1 Cd34 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM1 Cd34 split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", split.by = "stim", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
DefaultAssay(All.combined.EGFPPos1.FBSM1)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM1 Alcam UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", features = c("Alcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM1 Alcam split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", split.by = "stim", features = c("Alcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = All.combined.EGFPPos1.FBSM1) <- "seurat_clusters"
All.combined.EGFPPos1.FBSM1 <- RenameIdents(object = All.combined.EGFPPos1.FBSM1, 
                                     '9' = "FB1", '7' ="FB2", '4' ="FB3",'1' ="FB4",'0' ="FB5",
                                     '11' ="FB6",'2' ="FB7",'3' ="FB8",'6' ="FB9",'18' ="FB10",'13' ="FB11",
                                     '14' = "MyoFB",  
                                     '10' = "SM1",'5' = "SM2",'8' = "SM3",'15' = "SM4",'17' = "SM5",
                                     '16' = "MyoE1",'19' = "MyoE2", '12' = "MyoE3")
All.combined.EGFPPos1.FBSM1[["FBSMCellTypes"]] <- Idents(object = All.combined.EGFPPos1.FBSM1)

#FBSMCellTypes UMAP
Idents(object = All.combined.EGFPPos1.FBSM1) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM1 UMAP FBSMCellTypes.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue1", "olivedrab2", "blue", "deeppink1",  "brown3", "darkorange", "yellow2", "red", "turquoise3",
                                                                                              "bisque3", "slategray3","palegreen2",   "palevioletred3", "khaki4", "goldenrod2", "royalblue1","mediumorchid3",
                                                                                               "green4", "black"))
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM1 UMAP FBSMCellTypes split.tiff", width = 15, height =5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM1, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("salmon", "skyblue1", "olivedrab2", "blue", "deeppink1",  "brown3", "darkorange", "yellow2", "red", "turquoise3",
                                                                                               "bisque3", "slategray3","palegreen2",   "palevioletred3", "khaki4", "goldenrod2", "royalblue1","mediumorchid3",
                                                                                               "green4", "black"))
dev.off()

#Cell counts
Idents(object = All.combined.EGFPPos1.FBSM1) <- "FBSMCellTypes"
All.combined.EGFPPos1.FBSM1$stim.FBSMCellTypes <- paste(Idents(All.combined.EGFPPos1.FBSM1), All.combined.EGFPPos1.FBSM1$stim, sep = "_")
Idents(object = All.combined.EGFPPos1.FBSM1) <- "stim.FBSMCellTypes"
table(Idents(All.combined.EGFPPos1.FBSM1))

#DEGs
Idents(object = All.combined.EGFPPos1.FBSM1) <- "FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM1) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM1)
All.combined.EGFPPos1.FBSM1 <- ScaleData(All.combined.EGFPPos1.FBSM1, features = all.genes)
All.combined.EGFPPos1.FBSM1.FBSMCellTypes.markers <- FindAllMarkers(All.combined.EGFPPos1.FBSM1, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.1)
write.csv(All.combined.EGFPPos1.FBSM1.FBSMCellTypes.markers, file = "All.combined.EGFPPos1.FBSM1.FBSMCellTypes.markers.csv")

#Dotplot
Idents(object = All.combined.EGFPPos1.FBSM1) <- "FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM1) <- "RNA"
tiff(file = "All.combined.EGFPPos1.FBSM1 FBSMCellTypes markers DotPlot.tiff", width =12 , height = 3.1, units = "in", compression = "lzw", res = 800)
DotPlot(All.combined.EGFPPos1.FBSM1, features = c("Dach2", "Nefm", "H19", "Fam19a2", "Sh2d5","Abcg2", "Nefl","Aard", "Fabp4",
                                              "Sox11", "Egr1", "Jun", "Fos", "Nr2f1", "Fosb", "Hspa1a", "Socs3",
                                              "Epha3", "Pdgfrb", "Aldh1a3", "Col26a1", "Rgcc", "Col13a1", "Stc1",
                                              "Inhba", "Scg2", "Wisp1", "Hmga2", "Tac1", "Fam150b", "Areg", "Sox4",
                                              "Aspn", "Ugdh", "Sult1e1", "Clec3b", "Lum", "Thbd", "Ogn", "Dpt", "Gpx3",
                                              "Ifit1", "C2", "Gm17416", "Wif1", "Rtp4", "Parp14", "Isg15", "Tor3a", "Pbsn",
                                              "Sbp", "Ly6c1", "Defb1", "Ecm1", "Entpd2", "Gm7714", "Rnase1", "Itgbl1", "Alpl", "Ly6a", "Casp4", 
                                              "Klf11", "Pbk", "Mki67", "Racgap1", "Hist1h2ap", "Shcbp1", "Incenp", "Spc24", "Prc1",
                                              "Lmod1", "Npnt", "Pgm5", "Tnnt1", "Fbxl22", "Plcb4", "Smoc2", "Nexn", "Cald1", "Hspb6", "Fermt2", "Hacd1", "Flnc", "Cdc42ep3",
                                              "Sh3bgr", "Tnnt2", "Csrp1", "Hspb7", "Tinagl1", "Rbfox3", "Flna", "Jph2", "Synm", "Synpo2", "Ctxn1",
                                              "Tpm2", "Pcp4", "Ckb", "Cav1", "H2afz", "Sfn", "Gm10709", "Tuba1c", "Ccdc107", "Fhl1",
                                              "Tacstd2", "Ly6d", "Fxyd3", "Cldn7", "Fgfbp1", "Kcnn4", "Plet1", "Cldn4", "Pkp3", 
                                              "Dsp", "Col17a1", "Cdcp1", "Fermt1", "Serpinb5", "Lamb3", "Irf6", "Ocln", "Lsr", "F11r",
                                              "Pnoc", "Hgf", "Fendrr", "Rab3b", "Ednrb", "Fam129a", "Mustn1", 
                                              "Tbx18", "Ndufa4l2", "Spon2", "Gata6", "Dkk3", "Tiparp", "Bnip3", "Col14a1", "Fam213a", "Htra1", "Thbs2"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

tiff(file = "All.combined.EGFPPos1.FBSM1 FBSMCellTypes markers DotPlot.tiff", width =15 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM1, features = c("Dach2", "Nefm", "Fam19a2", "Sh2d5","Aard",
                                                 "Egr1",  "Fos", "Fosb", "Hspa1a", "Socs3",
                                                 "Epha3", "Pdgfrb", "Aldh1a3", "Col26a1", "Rgcc", 
                                                 "Inhba", "Wisp1", "Fam150b", "Areg", "Sox4",
                                                 "Aspn", "Ugdh", "Sult1e1", "Lum", "Gpx3",
                                                 "Ifit1", "C2", "Wif1", "Parp14", "Tor3a",
                                                 "Sbp", "Ly6c1", "Entpd2","Itgbl1","Casp4", 
                                                 "Pbk", "Mki67", "Racgap1","Spc24", "Prc1",
                                                 "Pgm5", "Tnnt1", "Plcb4", "Smoc2",  "Fermt2",
                                                 "Sh3bgr", "Tnnt2", "Hspb7","Synm", "Synpo2", 
                                                 "Gadd45g", "Tuba1c", "Tuba1a", "Actb", "Cfl1",
                                                 "Tacstd2", "Ly6d", "Cldn7", "Fgfbp1", "Kcnn4",
                                                 "Dsp", "Col17a1", "Cdcp1", "Fermt1", "Serpinb5", 
                                                 "Hgf", "Fendrr", "Rab3b", "Ednrb", "Fam129a", 
                                                 "Tbx18", "Ndufa4l2", "Spon2", "Gata6",  "Fam213a"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

tiff(file = "All.combined.EGFPPos1.FBSM1 FBSM canonical markers DotPlot.tiff", width =10 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM1, features = c("Vim", "Fbln1", "Fbn1", "Apod", "Fbln2",
                                                 "Acta2", "Tagln", "Myh11", "Actg2", "Cnn1", 
                                                 "Cdh1", "Krt4", "Krt5", "Krt6a", "S100a14"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Option2
#Cell Cycle regression
All.combined.EGFPPos1.FBSM2 <- All.combined.EGFPPos1.FBSM
DefaultAssay(All.combined.EGFPPos1.FBSM2) <- "integrated"
All.combined.EGFPPos1.FBSM2 <- ScaleData(All.combined.EGFPPos1.FBSM2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(All.combined.EGFPPos1.FBSM2))
All.combined.EGFPPos1.FBSM2 <- RunPCA(All.combined.EGFPPos1.FBSM2, features = VariableFeatures(All.combined.EGFPPos1.FBSM2))
ElbowPlot(All.combined.EGFPPos1.FBSM2, ndims = 50)

All.combined.EGFPPos1.FBSM2 <- FindNeighbors(All.combined.EGFPPos1.FBSM2, reduction = "pca", dims = 1:25)
All.combined.EGFPPos1.FBSM2 <- FindClusters(All.combined.EGFPPos1.FBSM1, resolution = 2.5)
All.combined.EGFPPos1.FBSM2 <- RunUMAP(All.combined.EGFPPos1.FBSM2, reduction = "pca", dims = 1:25)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, label = TRUE)

#Optimize resolution
DefaultAssay(All.combined.EGFPPos1.FBSM2) <- "integrated"
All.combined.EGFPPos1.FBSM2 <- ScaleData(All.combined.EGFPPos1.FBSM2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(All.combined.EGFPPos1.FBSM2))
All.combined.EGFPPos1.FBSM2 <- RunPCA(All.combined.EGFPPos1.FBSM2, features = VariableFeatures(All.combined.EGFPPos1.FBSM2))
ElbowPlot(All.combined.EGFPPos1.FBSM2, ndims = 50)

tiff(file = "All.combined.EGFPPos1.FBSM2 elbowplot.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
ElbowPlot(All.combined.EGFPPos1.FBSM2, ndims = 50)
dev.off()

All.combined.EGFPPos1.FBSM2 <- FindNeighbors(All.combined.EGFPPos1.FBSM2, reduction = "pca", dims = 1:30)
All.combined.EGFPPos1.FBSM2 <- FindClusters(All.combined.EGFPPos1.FBSM2, resolution = 0.7)
All.combined.EGFPPos1.FBSM2 <- RunUMAP(All.combined.EGFPPos1.FBSM2, reduction = "pca", dims = 1:30)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = All.combined.EGFPPos1.FBSM2) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM2 UMAP dims30 res0.7.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
Idents(object = All.combined.EGFPPos1.FBSM2) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM2 split UMAP dims30 res0.7.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
dev.off()

Idents(object = All.combined.EGFPPos1.FBSM2) <- "Phase"
tiff(file = "All.combined.EGFPPos1.FBSM2 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, cols = c("orange","lightseagreen", "magenta2"))
dev.off()

#Celltype identification
#FeaturePlot
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM2 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", features = c("Epcam", "Cdh1", "Krt5", "Krt19", "Upk3b", "Syp", 
                                                                          "Vim",  "Fbln1", "Myh11", "Acta2", 
                                                                          "Pecam1", "Plp1", "Rgs5","Tyrobp", "Myod1","Crabp1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
#Dotplot
Idents(object = All.combined.EGFPPos1.FBSM2) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM2 FBSM canonical markers DotPlot.tiff", width =10 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM2, features = c("Vim", "Fbln1", "Fbn1", "Apod", "Fbln2",
                                                  "Acta2", "Tagln", "Myh11", "Actg2", "Cnn1", 
                                                  "Cdh1", "Krt4", "Krt5", "Krt6a", "S100a14"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#FeaturePlots
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM2 Ar UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Ar split UMAP.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Gli1 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Gli1 split UMAP.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM2 Cd34 UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Cd34 split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM2 Alcam UMAP.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", features = c("Alcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Alcam split UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Alcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = All.combined.EGFPPos1.FBSM2) <- "seurat_clusters"
All.combined.EGFPPos1.FBSM2 <- RenameIdents(object = All.combined.EGFPPos1.FBSM2, 
                                            '11' = "FB1", '9' ="FB2", '1' ="FB3",'0' ="FB4",'3' ="FB5",
                                            '2' ="FB6",'8' ="FB7",'10' ="FB8",'4' ="FB9",
                                            '14' = "MyoFB",  
                                            '5' = "SM1",'6' = "SM1",'7' = "SM2",'15' = "SM3",
                                            '13' = "MyoE1",'16' = "MyoE2", '12' = "MyoE3")
All.combined.EGFPPos1.FBSM2[["FBSMCellTypes"]] <- Idents(object = All.combined.EGFPPos1.FBSM2)

#FBSMCellTypes UMAP
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 UMAP FBSMCellTypes.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                               "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                               "black"))
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 UMAP FBSMCellTypes split.tiff", width = 15, height =5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                                    "black"))
dev.off()

#Cell counts
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
All.combined.EGFPPos1.FBSM2$stim.FBSMCellTypes <- paste(Idents(All.combined.EGFPPos1.FBSM2), All.combined.EGFPPos1.FBSM2$stim, sep = "_")
Idents(object = All.combined.EGFPPos1.FBSM2) <- "stim.FBSMCellTypes"
table(Idents(All.combined.EGFPPos1.FBSM2))

#DEGs
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM2) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM2)
All.combined.EGFPPos1.FBSM2 <- ScaleData(All.combined.EGFPPos1.FBSM2, features = all.genes)
All.combined.EGFPPos1.FBSM2.FBSMCellTypes.markers <- FindAllMarkers(All.combined.EGFPPos1.FBSM2, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.1)
write.csv(All.combined.EGFPPos1.FBSM2.FBSMCellTypes.markers, file = "All.combined.EGFPPos1.FBSM2.FBSMCellTypes.markers.csv")


#Dotplot
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM2) <- "RNA"
tiff(file = "All.combined.EGFPPos1.FBSM2 FBSMCellTypes markers DotPlot.tiff", width =12 , height = 3.1, units = "in", compression = "lzw", res = 800)
DotPlot(All.combined.EGFPPos1.FBSM2, features = c("Dach2", "Hmga2", "Nefm", "H19", "Cytip","Smox", "Tac1","Mest", "Nefl", "Dlk1",
                                                  "Scg2", "Wisp1", "Inhba", "Areg", "Penk", "Mmp3",
                                                  "Rpl35a", "Rps27a", "Rpl9", "Rps19", "Rpl13", "Rpl24", "Erdr1", "Pim1", "Marcksl1", "Igfbp3", "Malat1", "C77370", "Mad2l1",
                                                  "Dcn", "Aspn", "Ugdh", "Thbd", "Clec3b", "Gfpt2", "Ogn", "Lox", "Ifi203",
                                                  "Sox11", "Bmp2", "Nav3", "Plat", "Ier3", "Stc1", "Cited2", "Emb", "Fam46a",
                                                  "Fosb", "Fos", "Klf4", "Jun", "Tob1", "Socs3", "Hspa1a", "Nr4a1", "Gadd45b", "Mfap4", "Dkk2", "Sfrp2", "Col15a1",
                                                  "Six2", "Epha3", "Pdgfrb", "Rgcc", "Col13a1", "Aldh1a3", "Igfbp5", "Cthrc1",  "Nrk", "Igfbp2", "Tcf4", 
                                                  "Tbx18", "Ndufa4l2", "Rgs17", "Htra1", "Cyp1a1", "Ereg", "Spon2", "Col14a1", "Vcam1", "Itm2a",
                                                  "Sbp", "9530002B09Rik", "Defb50", "Sbpl", "Pbsn", "Isg15", "Ly6a", "Casp4", "Entpd2", 
                                                  "Wfdc15b", "Tacstd2", "Sfn", "Crb3", "Cldn7", "Cldn4", "Cldn3", "Pkp3", "Fxyd3", "Lcn2", "Plet1", "Sprr1a",
                                                  "Krt6a", "Krt16",  "Serpinb5", "Foxa1", "Spint2", "Ly6d",
                                                  "Dsp", "Irf6", "Kcnn4", "Col17a1", "F2rl1",  "Itga3", "Cdcp1", "Lsr", "Fermt1", "Ripk4", "Cbr2", "Bcam", "Eps8l2", "Lgals7", "F11r",
                                                  "Top2a", "Prc1", "Racgap1", "Pbk", "Spc24", "Tpx2", "Incenp", "Shcbp1", "Kif11", "Spc25",
                                                  "Tnnt2", "Sh3bgr", "Rbfox3", "Tpm2", "Fbxl22", "Lmod1", "Tinagl1",  "Ckb", "Slmap", "Pgm5", "Cdc42ep3", "Nexn",
                                                  "Fhl1", "Cav1", "Mgll", "Flna", "Npnt", "Ctxn1", "Tnnt1", "Dmpk", "Lpp", "Mustn1",
                                                  "Pcp4l1", "Inhbb", "Frem2", "Pcp4", "Wnt7b", "Dpp6", "Esr1", "Vps37b", "Entpd1",  "Tmem200a",
                                                  "Gm27187", "2310057J18Rik", "Esp8", "Mir143hg", "Dmd", "Aqp1", "Tuba4a", "Pgr"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

tiff(file = "All.combined.EGFPPos1.FBSM2 FBSMCellTypes markers DotPlot.tiff", width =15 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM2, features = c("Dach2", "Nefm", "H19",  "Tac1","Mest", 
                                                  "Scg2", "Wisp1", "Inhba", "Areg",  "Mmp3",
                                                  "Rpl35a", "Rpl9", "Rpl13", "Erdr1", "Ddx5", 
                                                 "Aspn",  "Thbd", "Clec3b",  "Ogn", "Lox",
                                                   "Bmp2", "Nav3", "Ier3", "Cited2", "Emb",
                                                  "Fosb", "Fos", "Klf4", "Tob1", "Hspa1a", 
                                                  "Epha3", "Pdgfrb", "Rgcc", "Col13a1", "Aldh1a3",
                                                 "Ndufa4l2", "Htra1", "Spon2", "Col14a1", "Vcam1",
                                                  "Sbp", "9530002B09Rik", "Defb50", "Sbpl", "Pbsn", 
                                                 "Top2a", "Prc1", "Racgap1", "Pbk", "Spc24", 
                                                 "Tnnt2", "Fhl1", "Cdc42ep3", "Fbxl22", "Hspb6", 
                                                 "Inhbb", "Frem2", "Dpp6", "Esr1", "Tmem200a",
                                                 "Gm27187", "2310057J18Rik", "Esp8", "Mir143hg", "Pgr",
                                                  "Wfdc15b", "Crb3", "Lcn2", "Cdh16", "Spink8",
                                                 "Sprr1a","Krt6a", "Krt16",  "Spint2", "Ly6d",
                                                  "Dsp", "Irf6", "Col17a1", "Pkp3",  "Itga3"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

tiff(file = "All.combined.EGFPPos1.FBSM2 FBSMCelltype canonical markers DotPlot.tiff", width =10 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM2, features = c("Vim", "Fbln1", "Fbn1", "Apod", "Fbln2",
                                                  "Acta2", "Tagln", "Myh11", "Actg2", "Cnn1", 
                                                  "Cdh1", "Krt4", "Krt5", "Krt6a", "S100a14"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#FeaturePlots
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM2 Stra6 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Stra6"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Rspo3 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Rspo3"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Vlnplots

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

DefaultAssay(All.combined.EGFPPos1.FBSM2) <- "RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Dach2 Vln.tiff", width = 10, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM2, features = "Dach2", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

####Finding marker for each cluster####
#Rename
Idents(object = All.combined.EGFPPos1.FBSM2) <- "stim"
DimPlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", pt.size = 0.3, label = TRUE)
All.combined.EGFPPos1.FBSM2 <- RenameIdents(object = All.combined.EGFPPos1.FBSM2, 
                                            'E18_Ctrl' = "Ctrl", 'P11_Ctrl' ="Ctrl", 'P42_Ctrl' ="Ctrl",
                                            'E18_ARKO' ="ARKO",'P11_ARKO' ="ARKO",'P42_ARKO' ="ARKO")
All.combined.EGFPPos1.FBSM2[["CtrlvARKO"]] <- Idents(object = All.combined.EGFPPos1.FBSM2)

#Cell counts
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
All.combined.EGFPPos1.FBSM2$CtrlvARKO.FBSMCellTypes <- paste(Idents(All.combined.EGFPPos1.FBSM2), All.combined.EGFPPos1.FBSM2$CtrlvARKO, sep = "_")
Idents(object = All.combined.EGFPPos1.FBSM2) <- "CtrlvARKO.FBSMCellTypes"
table(Idents(All.combined.EGFPPos1.FBSM2))

#DEGs
Idents(object = All.combined.EGFPPos1.FBSM2) <- "CtrlvARKO.FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM2) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM2)
All.combined.EGFPPos1.FBSM2 <- ScaleData(All.combined.EGFPPos1.FBSM2, features = all.genes)
All.combined.EGFPPos1.FBSM2.CtrlvARKO.FBSMCellTypes.markers <- FindAllMarkers(All.combined.EGFPPos1.FBSM2, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.1)
write.csv(All.combined.EGFPPos1.FBSM2.CtrlvARKO.FBSMCellTypes.markers, file = "All.combined.EGFPPos1.FBSM2.CtrlvARKO.FBSMCellTypes.markers.csv")

#All.combined.EGFPPos1.FBSM2.Ctrl
Idents(object = All.combined.EGFPPos1.FBSM2) <- "CtrlvARKO"
All.combined.EGFPPos1.FBSM2.Ctrl <- subset(All.combined.EGFPPos1.FBSM2, idents = c("Ctrl"))

#DEGs
Idents(object = All.combined.EGFPPos1.FBSM2.Ctrl) <- "CtrlvARKO.FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM2.Ctrl) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM2.Ctrl)
All.combined.EGFPPos1.FBSM2.Ctrl <- ScaleData(All.combined.EGFPPos1.FBSM2.Ctrl, features = all.genes)
All.combined.EGFPPos1.FBSM2.Ctrl.FBSMCellTypes.markers <- FindAllMarkers(All.combined.EGFPPos1.FBSM2.Ctrl, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.1)
write.csv(All.combined.EGFPPos1.FBSM2.Ctrl.FBSMCellTypes.markers, file = "All.combined.EGFPPos1.FBSM2.Ctrl.FBSMCellTypes.markers.csv")

#FeaturePlot for Telocytes
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Foxl1 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Foxf1 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Foxf1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Tcf4 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Tcf4"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Cd34 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Epha3 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Epha3"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Aldh1a3 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Aldh1a3"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Col13a1 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Col13a1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM2 Cthrc1 split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Cthrc1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FeaturePlot for FB1
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FeaturePlot for MyoE2
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FeaturePlot for FB9
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#All.combined.EGFPPos1.FBSM2.ARKO
Idents(object = All.combined.EGFPPos1.FBSM2) <- "CtrlvARKO"
All.combined.EGFPPos1.FBSM2.ARKO <- subset(All.combined.EGFPPos1.FBSM2, idents = c("ARKO"))

#DEGs
Idents(object = All.combined.EGFPPos1.FBSM2.ARKO) <- "CtrlvARKO.FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM2.ARKO) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM2.ARKO)
All.combined.EGFPPos1.FBSM2.ARKO <- ScaleData(All.combined.EGFPPos1.FBSM2.ARKO, features = all.genes)
All.combined.EGFPPos1.FBSM2.ARKO.FBSMCellTypes.markers <- FindAllMarkers(All.combined.EGFPPos1.FBSM2.ARKO, min.pct = 0.1, only.pos = TRUE, logfc.threshold = 0.1)
write.csv(AAll.combined.EGFPPos1.FBSM2.ARKO.FBSMCellTypes.markers, file = "All.combined.EGFPPos1.FBSM2.ARKO.FBSMCellTypes.markers.csv")

#FeaturePlot for SM2
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FeaturePlot for SM3
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FeaturePlot for FB3
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FeaturePlot for FB6
DefaultAssay(All.combined.EGFPPos1.FBSM2)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM2) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM2 Pdgfra split UMAP.tiff", width = 18, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM2, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

save.image("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All/All.combined.EGFPPos1.FBSM2.RData")

####All.combined.EGFPPos1.FBSM3####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All/All.combined.EGFPPos1.FBSM3")

#Cell Cycle regression
All.combined.EGFPPos1.FBSM3 <- All.combined.EGFPPos1.FBSM2
DefaultAssay(All.combined.EGFPPos1.FBSM3) <- "integrated"
All.combined.EGFPPos1.FBSM3 <- ScaleData(All.combined.EGFPPos1.FBSM3, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(All.combined.EGFPPos1.FBSM3))
All.combined.EGFPPos1.FBSM3 <- RunPCA(All.combined.EGFPPos1.FBSM3, features = VariableFeatures(All.combined.EGFPPos1.FBSM3))
ElbowPlot(All.combined.EGFPPos1.FBSM3, ndims = 50)

tiff(file = "All.combined.EGFPPos1.FBSM3 elbowplot.tiff", width = 4.5, height = 4, units = "in", compression = "lzw", res = 200)
ElbowPlot(All.combined.EGFPPos1.FBSM3, ndims = 50)
dev.off()

All.combined.EGFPPos1.FBSM3 <- FindNeighbors(All.combined.EGFPPos1.FBSM3, reduction = "pca", dims = 1:40)
All.combined.EGFPPos1.FBSM3 <- FindClusters(All.combined.EGFPPos1.FBSM3, resolution = 0.8)
All.combined.EGFPPos1.FBSM3 <- RunUMAP(All.combined.EGFPPos1.FBSM3, reduction = "pca", dims = 1:40)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = All.combined.EGFPPos1.FBSM3) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM3 UMAP dims40 res0.8.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
Idents(object = All.combined.EGFPPos1.FBSM3) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM3 split UMAP dims40 res0.8.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
dev.off()

Idents(object = All.combined.EGFPPos1.FBSM3) <- "Phase"
tiff(file = "All.combined.EGFPPos1.FBSM3 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, cols = c("orange","lightseagreen", "magenta2"))
dev.off()

#Celltype identification
#FeaturePlot
DefaultAssay(All.combined.EGFPPos1.FBSM3)<-"RNA"
tiff(file = "All.combined.EGFPPos1.FBSM3 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", features = c("Vim", "Fbln1", "Fbln2", "Fbn1", 
                                                                          "Pdgfra",  "Foxl1", "Tcf4", "Cd34", 
                                                                          "Acta2", "Tagln", "Myh11", "Actg2",
                                                                          "Cdh1", "Krt14", "Krt6a", "S100a14"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Dotplot
Idents(object = All.combined.EGFPPos1.FBSM3) <- "seurat_clusters"
tiff(file = "All.combined.EGFPPos1.FBSM3 FBSM canonical markers DotPlot.tiff", width =10 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM3, features = c("Vim", "Fbln1", "Fbn1", "Apod", "Fbln2",
                                                  "Acta2", "Tagln", "Myh11", "Actg2", "Cnn1", 
                                                  "Cdh1", "Krt4", "Krt5", "Krt6a", "S100a14"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Rename
Idents(object = All.combined.EGFPPos1.FBSM3) <- "seurat_clusters"
All.combined.EGFPPos1.FBSM3 <- RenameIdents(object = All.combined.EGFPPos1.FBSM3, 
                                            '5' = "FB1", '11' = "FB1",  '0' ="FB2",'19' ="FB2", '3' ="FB3",'1' ="FB4",
                                            '4' ="FB5",
                                            '2' ="FB6",'9' ="SEFB1",'18' ="SEFB2", '12' ="MM",
                                            '14' = "MyoFB",  
                                            '10' = "SM1",'8' = "SM1",'6' = "SM1",'7' = "SM2",'15' = "SM3",
                                            '16' = "MyoE1",'17' = "MyoE2", '13' = "MyoE2")
All.combined.EGFPPos1.FBSM3[["FBSMCellTypes"]] <- Idents(object = All.combined.EGFPPos1.FBSM3)

#FBSMCellTypes UMAP
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM3 UMAP FBSMCellTypes.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                               "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                               "black"))
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 UMAP FBSMCellTypes split.tiff", width = 15, height =5, units = "in", compression = "lzw", res = 200)
DimPlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                                    "black"))
dev.off()

#Cell counts
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
All.combined.EGFPPos1.FBSM3$stim.FBSMCellTypes <- paste(Idents(All.combined.EGFPPos1.FBSM3), All.combined.EGFPPos1.FBSM3$stim, sep = "_")
Idents(object = All.combined.EGFPPos1.FBSM3) <- "stim.FBSMCellTypes"
table(Idents(All.combined.EGFPPos1.FBSM3))

#Dotplot
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM3 FBSMCellTypes canonical markers DotPlot.tiff", width =10 , height = 5, units = "in", compression = "lzw", res = 200)
DotPlot(All.combined.EGFPPos1.FBSM3, features = c("Vim", "Fbln1", "Fbn1", "Apod", "Fbln2",
                                                  "Pdgfra", "Foxl1", "Tcf4", "Cd34",
                                                  "Acta2", "Tagln", "Myh11", "Actg2", "Cnn1", 
                                                   "Cdh1", "Krt4", "Krt5", "Krt6a", "S100a14"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#DEGs
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
DefaultAssay(All.combined.EGFPPos1.FBSM3) <- "RNA"
all.genes <- rownames(All.combined.EGFPPos1.FBSM3)
All.combined.EGFPPos1.FBSM3 <- ScaleData(All.combined.EGFPPos1.FBSM3, features = all.genes)
All.combined.EGFPPos1.FBSM3.FBSMCellTypes.markers <- FindAllMarkers(All.combined.EGFPPos1.FBSM3, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(All.combined.EGFPPos1.FBSM3.FBSMCellTypes.markers, file = "All.combined.EGFPPos1.FBSM3.FBSMCellTypes.markers.csv")

save.image("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All/All.combined.EGFPPos1.FBSM3/All.combined.EGFPPos1.FBSM3.Monocle3.RData")

#Telocyte1v2
#DEGs_hMETPosEpi_TriplevsDouble
DefaultAssay(All.combined.EGFPPos1.FBSM3) <- "RNA"
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
all.genes <- rownames(All.combined.EGFPPos1.FBSM3)
All.combined.EGFPPos1.FBSM3 <- ScaleData(All.combined.EGFPPos1.FBSM3, features = all.genes)
Telocyte1v2.0.Markers <- FindMarkers(All.combined.EGFPPos1.FBSM3, ident.1 = c("Telocyte1"), 
                                                      ident.2 = c("Telocyte2"), min.pct = 0, logfc.threshold = 0)
write.csv(Telocyte1v2.0.Markers, "Telocyte1v2.0.Markers.csv")

#p.adjust
DEG_Telocyte1v2 <- read.csv("Telocyte1v2.0.Markers.csv") 
DEG_Telocyte1v2_pvalue <- DEG_Telocyte1v2$p_val
DEG_Telocyte1v2_pvalue=as.numeric(DEG_Telocyte1v2_pvalue)
DEG_Telocyte1v2_BH = p.adjust(DEG_Telocyte1v2_pvalue, "BH")
write.csv(DEG_Telocyte1v2_BH, "DEG_Telocyte1v2_BH.csv")

#Telocyte1vothers
#DEGs_hMETPosEpi_TriplevsDouble
DefaultAssay(All.combined.EGFPPos1.FBSM3) <- "RNA"
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
all.genes <- rownames(All.combined.EGFPPos1.FBSM3)
All.combined.EGFPPos1.FBSM3 <- ScaleData(All.combined.EGFPPos1.FBSM3, features = all.genes)
Telocyte1vothers.0.Markers <- FindMarkers(All.combined.EGFPPos1.FBSM3, ident.1 = c("Telocyte1"), min.pct = 0, logfc.threshold = 0)
write.csv(Telocyte1vothers.0.Markers, "Telocyte1vothers.0.Markers.csv")

#p.adjust
DEG_Telocyte1vothers <- read.csv("Telocyte1vothers.0.Markers.csv") 
DEG_Telocyte1vothers_pvalue <- DEG_Telocyte1vothers$p_val
DEG_Telocyte1vothers_pvalue=as.numeric(DEG_Telocyte1vothers_pvalue)
DEG_Telocyte1vothers_BH = p.adjust(DEG_Telocyte1vothers_pvalue, "BH")
write.csv(DEG_Telocyte1vothers_BH, "DEG_Telocyte1vothers_BH.csv")

#FeaturePlot for Telocytes
DefaultAssay(All.combined.EGFPPos1.FBSM3)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM3 Ar split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Gli1 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Foxl1 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Pdgfra split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q99")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Cd34 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()

#FeaturePlot for Telocyte1
tiff(file = "All.combined.EGFPPos1.FBSM3 Aldh1a3 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Aldh1a3"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Col13a1 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Col13a1"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Epha3 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Epha3"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Cthrc1 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Cthrc1"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Tnc split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Tnc"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Igfbp5 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Igfbp5"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q99")
dev.off()

tiff(file = "All.combined.EGFPPos1.FBSM3 Igf1 split UMAP.tiff", width = 17, height = 3, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", split.by = "stim", features = c("Igf1"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()

#Volin plot
DefaultAssay(All.combined.EGFPPos1.FBSM3) <- "RNA"
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"

tiff(file = "All.combined.EGFPPos1.FBSM3 Pdgfra Vln.tiff", width = 12, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Pdgfra", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO") + NoLegend()
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Cd34 Vln.tiff", width = 12, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Cd34", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO") + NoLegend()
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Foxl1 Vln.tiff", width = 12, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Foxl1", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO") + NoLegend()
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Foxf1 Vln.tiff", width = 12, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Foxf1", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO") + NoLegend()
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Epha3 Vln.tiff", width = 12.5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Epha3", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO")
dev.off()

tiff(file = "All.combined.EGFPPos1.FBSM3 Igfbp3 Vln.tiff", width = 12, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Igfbp3", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO") + NoLegend()
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Igfbp5 Vln.tiff", width = 12, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(All.combined.EGFPPos1.FBSM3, features = "Igfbp5", pt.size = 0, split.plot = TRUE, split.by = "CtrlvARKO") + NoLegend()
dev.off()

#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named Symbol
genelist = DEG_Telocyte1vothers_preranked$symbol
ortho_DEGs_Telocyte1vothers_preranked<- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

DEGs_Telocyte1vothers_preranked <- inner_join(DEG_Telocyte1vothers_preranked,ortho_DEGs_Telocyte1vothers_preranked,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(DEGs_Telocyte1vothers_preranked,"DEGs_Telocyte1vothers_converted.txt")

####Monocle3_ALL####

BiocManager::install(version = "3.14")
BiocManager::install(c('BiocGenerics', 'DelayedArray', 'DelayedMatrixStats',
                       'limma', 'lme4', 'S4Vectors', 'SingleCellExperiment',
                       'SummarizedExperiment', 'batchelor', 'HDF5Array',
                       'terra', 'ggrastr'))
remotes::install_github('satijalab/seurat-wrappers')
devtools::install_github('cole-trapnell-lab/monocle3')

DefaultAssay(All.combined.EGFPPos1.FBSM3) <- "RNA"
cds <- as.cell_data_set(All.combined.EGFPPos1.FBSM3)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(All.combined.EGFPPos1.FBSM3)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "FBSMCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))

tiff(file = "All.combined.EGFPPos1.FBSM3 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)
plt

tiff(file = "All.combined.EGFPPos1.FBSM3 pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt
dev.off()

#plot genes in pseudotime
Solid_genes <- c("Bmpr1b")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$stim %in% c("E18_Ctrl")]

tiff(file = "All.combined.EGFPPos1.FBSM3 Bmpr1b in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 200)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="FBSMCellTypes"
) + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                  "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                  "black"))
dev.off()

#Cytotrace
EGFPPosCountsCytoTrace <- All.combined.EGFPPos1.FBSM3@assays$RNA@counts
write.table(EGFPPosCountsCytoTrace, file = "EGFPPosCountsCytoTrace_All.combined.EGFPPos1.FBSM3.txt")

phenotypedata = FetchData(All.combined.EGFPPos1.FBSM3, c("FBSMCellTypes"))
write.table(phenotypedata, file = "phenotypedata_All.combined.EGFPPos1.FBSM3.txt") #Change name "ATCT-1_1" to "ATCT.1_1"

#MSC markers
DefaultAssay(All.combined.EGFPPos1.FBSM3)<-"RNA"
Idents(object = All.combined.EGFPPos1.FBSM3) <- "FBSMCellTypes"
tiff(file = "All.combined.EGFPPos1.FBSM3 Thy1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", features = c("Thy1"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 CD34 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()
tiff(file = "All.combined.EGFPPos1.FBSM3 Alcam UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 200)
FeaturePlot(All.combined.EGFPPos1.FBSM3, reduction = "umap", features = c("Alcam"), cols = c("light grey", "red"), pt.size = 0.7, max.cutoff = "q90")
dev.off()

save.image("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARKOvCtrl_3timepoint_NEW/All/Monocle3/All.combined.EGFPPos1.FBSM2.RData")

####Monocle3_WT####
#without clustering
Idents(object = All.combined.EGFPPos1.FBSM3) <- "stim"
Ctrl.combined.EGFPPos1.FBSM3 <- subset(All.combined.EGFPPos1.FBSM3, idents = c("E18_Ctrl", "P11_Ctrl", "P42_Ctrl"))
DimPlot(Ctrl.combined.EGFPPos1.FBSM3, reduction = "umap", pt.size = 0.3, label = TRUE)
        
DefaultAssay(Ctrl.combined.EGFPPos1.FBSM3) <- "RNA"
cds <- as.cell_data_set(Ctrl.combined.EGFPPos1.FBSM3)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Ctrl.combined.EGFPPos1.FBSM3)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "FBSMCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))

tiff(file = "Ctrl.combined.EGFPPos1.FBSM3 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)
plt

tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt
dev.off()

#with clustering
#Cell Cycle regression
Ctrl.combined.EGFPPos1.FBSM4 <- Ctrl.combined.EGFPPos1.FBSM3
DefaultAssay(Ctrl.combined.EGFPPos1.FBSM4) <- "integrated"
Ctrl.combined.EGFPPos1.FBSM4 <- ScaleData(Ctrl.combined.EGFPPos1.FBSM4, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Ctrl.combined.EGFPPos1.FBSM4))
Ctrl.combined.EGFPPos1.FBSM4 <- RunPCA(Ctrl.combined.EGFPPos1.FBSM4, features = VariableFeatures(Ctrl.combined.EGFPPos1.FBSM4))
ElbowPlot(Ctrl.combined.EGFPPos1.FBSM4, ndims = 50)

tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 Elbow.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
ElbowPlot(Ctrl.combined.EGFPPos1.FBSM4, ndims = 50)
dev.off()

Ctrl.combined.EGFPPos1.FBSM4 <- FindNeighbors(Ctrl.combined.EGFPPos1.FBSM4, reduction = "pca", dims = 1:25)
Ctrl.combined.EGFPPos1.FBSM4 <- FindClusters(Ctrl.combined.EGFPPos1.FBSM4, resolution = 0.8)
Ctrl.combined.EGFPPos1.FBSM4 <- RunUMAP(Ctrl.combined.EGFPPos1.FBSM4, reduction = "pca", dims = 1:25)
Idents(object = Ctrl.combined.EGFPPos1.FBSM4) <- "FBSMCellTypes"
DimPlot(Ctrl.combined.EGFPPos1.FBSM4, reduction = "umap", pt.size = 0.3, label = TRUE)

#FBSMCellTypes UMAP
Idents(object = Ctrl.combined.EGFPPos1.FBSM4) <- "FBSMCellTypes"
tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 UMAP FBSMCellTypes.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
DimPlot(Ctrl.combined.EGFPPos1.FBSM4, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                               "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                               "black"))
dev.off()
tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 UMAP FBSMCellTypes split.tiff", width = 15, height =5, units = "in", compression = "lzw", res = 200)
DimPlot(Ctrl.combined.EGFPPos1.FBSM4, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                                    "black"))
dev.off()

##Convert Seurat to Monocle3 cell data set class
DefaultAssay(Ctrl.combined.EGFPPos1.FBSM4) <- "RNA"
cds <- as.cell_data_set(Ctrl.combined.EGFPPos1.FBSM4)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Ctrl.combined.EGFPPos1.FBSM4)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "FBSMCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))

tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#
plt <- plot_cells(cds,
                  color_cells_by = "stim",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 re-clustering stim trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)
plt

tiff(file = "Ctrl.combined.EGFPPos1.FBSM4 pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt
dev.off()

####Monocle3_ARKO####
#without clustering
Idents(object = All.combined.EGFPPos1.FBSM3) <- "stim"
ARKO.combined.EGFPPos1.FBSM3 <- subset(All.combined.EGFPPos1.FBSM3, idents = c("E18_ARKO", "P11_ARKO", "P42_ARKO"))

DefaultAssay(ARKO.combined.EGFPPos1.FBSM3) <- "RNA"
cds <- as.cell_data_set(ARKO.combined.EGFPPos1.FBSM3)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ARKO.combined.EGFPPos1.FBSM3)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "FBSMCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))

tiff(file = "ARKO.combined.EGFPPos1.FBSM3 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)
plt

tiff(file = "ARKO.combined.EGFPPos1.FBSM3 pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt
dev.off()

#with clustering
#Cell Cycle regression
ARKO.combined.EGFPPos1.FBSM4 <- ARKO.combined.EGFPPos1.FBSM3
DefaultAssay(ARKO.combined.EGFPPos1.FBSM4) <- "integrated"
ARKO.combined.EGFPPos1.FBSM4 <- ScaleData(ARKO.combined.EGFPPos1.FBSM4, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ARKO.combined.EGFPPos1.FBSM4))
ARKO.combined.EGFPPos1.FBSM4 <- RunPCA(ARKO.combined.EGFPPos1.FBSM4, features = VariableFeatures(ARKO.combined.EGFPPos1.FBSM4))
ElbowPlot(ARKO.combined.EGFPPos1.FBSM4, ndims = 50)

tiff(file = "ARKO.combined.EGFPPos1.FBSM4 Elbow.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
ElbowPlot(ARKO.combined.EGFPPos1.FBSM4, ndims = 50)
dev.off()

ARKO.combined.EGFPPos1.FBSM4 <- FindNeighbors(ARKO.combined.EGFPPos1.FBSM4, reduction = "pca", dims = 1:31)
ARKO.combined.EGFPPos1.FBSM4 <- FindClusters(ARKO.combined.EGFPPos1.FBSM4, resolution = 0.8)
ARKO.combined.EGFPPos1.FBSM4 <- RunUMAP(ARKO.combined.EGFPPos1.FBSM4, reduction = "pca", dims = 1:31)
Idents(object = ARKO.combined.EGFPPos1.FBSM4) <- "FBSMCellTypes"
DimPlot(ARKO.combined.EGFPPos1.FBSM4, reduction = "umap", pt.size = 0.3, label = TRUE)

#FBSMCellTypes UMAP
Idents(object = ARKO.combined.EGFPPos1.FBSM4) <- "FBSMCellTypes"
tiff(file = "ARKO.combined.EGFPPos1.FBSM4 UMAP FBSMCellTypes.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 200)
DimPlot(ARKO.combined.EGFPPos1.FBSM4, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                                "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                                "black"))
dev.off()
tiff(file = "ARKO.combined.EGFPPos1.FBSM4 UMAP FBSMCellTypes split.tiff", width = 15, height =5, units = "in", compression = "lzw", res = 200)
DimPlot(ARKO.combined.EGFPPos1.FBSM4, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                                                                                     "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                                                                                     "black"))
dev.off()

##Convert Seurat to Monocle3 cell data set class
DefaultAssay(ARKO.combined.EGFPPos1.FBSM4) <- "RNA"
cds <- as.cell_data_set(ARKO.combined.EGFPPos1.FBSM4)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(ARKO.combined.EGFPPos1.FBSM4)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "FBSMCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))

tiff(file = "ARKO.combined.EGFPPos1.FBSM4 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#
plt <- plot_cells(cds,
                  color_cells_by = "stim",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

tiff(file = "ARKO.combined.EGFPPos1.FBSM4 re-clustering stim trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt + scale_color_manual(values = c("salmon", "skyblue2", "olivedrab2", "brown4", "deeppink1",  "blue", "darkorange", "yellow2", "red", "turquoise1",
                                    "bisque2", "slategray3","mediumorchid3",  "khaki4", "goldenrod2", "green4", "royalblue2",
                                    "black"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1)
plt

tiff(file = "ARKO.combined.EGFPPos1.FBSM4 pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 200)
plt
dev.off()