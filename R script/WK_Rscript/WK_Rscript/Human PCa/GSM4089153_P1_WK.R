library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)

####GSM4089153####
setwd("//isi-dcnl/user_data/zjsun/group/Aly Buckley/human/GSM4089153 CRPC P3")

GSM4089153.data1 <- read.table(file = "GSM4089153_P3_gene_cell_exprs_table2.txt", sep = "", header = T, row.names=1, as.is=T)
GSM4089153_P3.1 <- CreateSeuratObject(counts = GSM4089153.data1, min.cells = 3, min.features = 500)
GSM4089153_P3.1 <- NormalizeData(GSM4089153_P3.1)
GSM4089153_P3.1 <- FindVariableFeatures(GSM4089153_P3.1, selection.method = "mean.var.plot", nfeatures = 5000)
#Run the standard workflow for visualization and clustering 
all.genes <- rownames(GSM4089153_P3.1)
GSM4089153_P3.1 <- ScaleData(GSM4089153_P3.1, features = all.genes)
GSM4089153_P3.1 <- RunPCA(GSM4089153_P3.1, npcs = 50, verbose = FALSE)
ElbowPlot(GSM4089153_P3.1, ndims = 30)
# UMAP and Clustering
GSM4089153_P3.1 <- FindNeighbors(GSM4089153_P3.1, reduction = "pca", dims = 1:20)
GSM4089153_P3.1 <- FindClusters(GSM4089153_P3.1, resolution = 0.5)
GSM4089153_P3.1 <- RunUMAP(GSM4089153_P3.1, reduction = "pca", dims = 1:20)
DimPlot(GSM4089153_P3.1, reduction = "umap", pt.size = 0.5, label = TRUE) 

Idents(object = GSM4089153_P3.1) <- "seurat_clusters"
tiff(file = "GSM4089153_P3.1 seurat_clusters UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089153_P3.1, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()

#DEGs
DefaultAssay(GSM4089153_P3.1)<-"RNA"
Idents(object = GSM4089153_P3.1) <- "seurat_clusters"
all.genes <- rownames(GSM4089153_P3.1)
GSM4089153_P3.1 <- ScaleData(GSM4089153_P3.1, features = all.genes)
GSM4089153_P3.1.seurat.markers <- FindAllMarkers(GSM4089153_P3.1, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(GSM4089153_P3.1.seurat.markers, file = "GSM4089153_P3.1.seurat.markers.csv")

tiff(file = "GSM4089153_P3.1 celltype marker expression plots.tiff", width = 16, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089153_P3.1, reduction = "umap", features = c("EPCAM", "VIM", "KRT5", "KRT14", "KRT8", "KRT19",
                                                              "FBLN1", "MYH11", "PECAM1",
                                                              "RGS1", "CCL5", "CHGA"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSM4089153_P3.1) <- "seurat_clusters"
DimPlot(GSM4089153_P3.1, reduction = "umap", pt.size = 0.3, label=TRUE)
GSM4089153_P3.1 <- RenameIdents(object = GSM4089153_P3.1, 
                              '0'="LE",'1'="BE",'2'="FB",'3'="Endo",'4'="Endo",
                              '5'="LE",'6'="BE",'7'="Leukocytes",'8'="LE",'9'="SM")  
GSM4089153_P3.1[["CellTypes"]] <- Idents(object = GSM4089153_P3.1)

#UMAP_CellTypes
Idents(object = GSM4089153_P3.1) <- "CellTypes"
tiff(file = "GSM4089153_P3.1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089153_P3.1, reduction = "umap", pt.size = 0.5, label = TRUE)
dev.off()
tiff(file = "GSM4089153_P3.1 MET WNT expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089153_P3.1, reduction = "umap", features = c("MET", "CTNNB1", "AR", "XPO1", "MYC", "TCF4", "EIF4A1", "RPL12", "RPS16"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
#Cell count
table(Idents(GSM4089153_P3.1))

GSM4089153_P3.1_epi <- subset(GSM4089153_P3.1, idents = c("BE","LE"))
Idents(object = GSM4089153_P3.1_epi) <- "seurat_clusters"

#Run the standard workflow for visualization and clustering
GSM4089153_P3.1_epi <- ScaleData(GSM4089153_P3.1_epi, verbose = FALSE)
GSM4089153_P3.1_epi <- RunPCA(GSM4089153_P3.1_epi, npcs = 50, verbose = FALSE)
ElbowPlot(GSM4089153_P3.1_epi, ndims = 50)

#Umap and Clustering
GSM4089153_P3.1_epi <- FindNeighbors(GSM4089153_P3.1_epi, reduction = "pca", dims = 1:19)
GSM4089153_P3.1_epi <- FindClusters(GSM4089153_P3.1_epi, resolution = 0.7)
GSM4089153_P3.1_epi <- RunUMAP(GSM4089153_P3.1_epi, reduction = "pca", dims = 1:19)
DimPlot(GSM4089153_P3.1_epi, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "GSM4089153_P3.1_epi EpiCellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089153_P3.1_epi, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = GSM4089153_P3.1_epi) <- "seurat_clusters"
FeaturePlot(GSM4089153_P3.1_epi, reduction = "umap", features = c("KRT8","KRT5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
GSM4089153_P3.1_epi <- RenameIdents(object = GSM4089153_P3.1_epi, 
                                '0'="BE1",'4'="BE2",'1'="LE1",'2'="LE2",'3'="LE3",
                                '5'="LE4",'6'="LE5")  
GSM4089153_P3.1_epi[["EpiCellTypes"]] <- Idents(object = GSM4089153_P3.1_epi)
Idents(object = GSM4089153_P3.1_epi) <- "EpiCellTypes"
DimPlot(GSM4089153_P3.1_epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#DEGs
DefaultAssay(GSM4089153_P3.1_epi)<-"RNA"
Idents(object = GSM4089153_P3.1_epi) <- "EpiCellTypes"
all.genes <- rownames(GSM4089153_P3.1_epi)
GSM4089153_P3.1_epi <- ScaleData(GSM4089153_P3.1_epi, features = all.genes)
GSM4089153_P3.1_epi.markers <- FindAllMarkers(GSM4089153_P3.1_epi, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(GSM4089153_P3.1_epi.markers, file = "GSM4089153_P3.1_epi.markers.csv")

tiff(file = "GSM4089153_P3.1_epi MET WNT expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089153_P3.1_epi, reduction = "umap", features = c("MET", "CTNNB1", "AR", "XPO1", "MYC", "TCF4", "EIF4A1", "RPL12", "RPS16"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####GSM4089151 P1####
setwd("//isi-dcnl/user_data/zjsun/group/Aly Buckley/human/GSM4089151 CRPC P1")

GSM4089151_P1.data <- read.table(file = "GSM4089151_P1_gene_cell_exprs_table2.txt", sep = "", header = T, row.names=1, as.is=T)
GSM4089151_P1 <- CreateSeuratObject(counts = GSM4089151_P1.data, min.cells = 3, min.features = 500)
GSM4089151_P1 <- NormalizeData(GSM4089151_P1)


tiff(file = "GSM4089151_P1 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "GSM4089151_P1 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(GSM4089151_P1@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "GSM4089151_P1 Pre-filteration")
dev.off()

GSM4089151_P1_1 <- subset(GSM4089151_P1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000)
GSM4089151_P1_1 <- FindVariableFeatures(GSM4089151_P1_1, selection.method = "vst", nfeatures = 5000)

all.genes <- rownames(GSM4089151_P1_1)
GSM4089151_P1_1 <- ScaleData(GSM4089151_P1_1, features = all.genes)
GSM4089151_P1_1 <- RunPCA(GSM4089151_P1_1, npcs = 50, verbose = FALSE)
ElbowPlot(GSM4089151_P1_1, ndims = 30)

# UMAP and Clustering
GSM4089151_P1_1 <- FindNeighbors(GSM4089151_P1_1, reduction = "pca", dims = 1:14)
GSM4089151_P1_1 <- FindClusters(GSM4089151_P1_1, resolution = 0.5)
GSM4089151_P1_1 <- RunUMAP(GSM4089151_P1_1, reduction = "pca", dims = 1:14)
DimPlot(GSM4089151_P1_1, reduction = "umap", pt.size = 0.5, label = TRUE) 

# Load Seurat and prepare cell cycle score data
library(Seurat)
data(cc.genes)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

DefaultAssay(GSM4089151_P1_1) <- "RNA"
all.genes <- rownames(GSM4089151_P1_1)
GSM4089151_P1_1 <- ScaleData(GSM4089151_P1_1, features = all.genes)
GSM4089151_P1_1 <- CellCycleScoring(GSM4089151_P1_1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = GSM4089151_P1_1) <- "Phase"
DimPlot(GSM4089151_P1_1, reduction = "umap")

tiff(file = "GSM4089151_P1_1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell Cycle regression
GSM4089151_P1_2 <- GSM4089151_P1_1
DefaultAssay(GSM4089151_P1_2) <- "RNA"
GSM4089151_P1_2 <- ScaleData(GSM4089151_P1_2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(GSM4089151_P1_2))
GSM4089151_P1_2 <- RunPCA(GSM4089151_P1_2, features = VariableFeatures(GSM4089151_P1_2))
ElbowPlot(GSM4089151_P1_2, ndims = 50)

GSM4089151_P1_2 <- FindNeighbors(GSM4089151_P1_2, reduction = "pca", dims = 1:15)
GSM4089151_P1_2 <- FindClusters(GSM4089151_P1_2, resolution = 0.5)
GSM4089151_P1_2 <- RunUMAP(GSM4089151_P1_2, reduction = "pca", dims = 1:15)
GSM4089151_P1_2 <- RunTSNE(GSM4089151_P1_2, reduction = "pca", dims = 1:15)
DimPlot(GSM4089151_P1_2, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = GSM4089151_P1_2) <- "Phase"
tiff(file = "GSM4089151_P1_2 after Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

tiff(file = "GSM4089151_P1_2 celltype marker expression plots.tiff", width = 16, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089151_P1_2, reduction = "umap", features = c("EPCAM", "VIM", "KRT5", "KRT14", "KRT8", "KRT19",
                                                                      "FBLN1", "MYH11", "PECAM1",
                                                                      "RGS1", "CCL5", "CHGA"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = GSM4089151_P1_2) <- "seurat_clusters"
DimPlot(GSM4089151_P1_2, reduction = "umap", pt.size = 0.3, label=TRUE)
GSM4089151_P1_2 <- RenameIdents(object = GSM4089151_P1_2, 
                              '2'="Epi",'4'="Epi",'7'="Epi",'1'="Epi",'0'="Epi",
                              '3'="Epi", '9'= "Epi", '5'="FB",'10'="VE",'8'="Immune",'11'="Immune",'12'="Immune",'6'="Immune")  
GSM4089151_P1_2[["CellTypes"]] <- Idents(object = GSM4089151_P1_2)

#UMAP_CellTypes
Idents(object = GSM4089151_P1_2) <- "CellTypes"
tiff(file = "GSM4089151_P1_2 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2, reduction = "umap", pt.size = 0.3)
dev.off()

####Subcluster Epi####
Idents(object = GSM4089151_P1_2) <- "CellTypes"
tiff(file = "GSM4089151_P1_2 Epi highlighted UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2, reduction = "umap", pt.size = 0.3, cols = c("salmon", "grey", "grey", "grey"))
dev.off()

GSM4089151_P1_2.epi <- subset(GSM4089151_P1_2, idents = c("Epi"))

#Run the standard workflow for visualization and clustering
GSM4089151_P1_2.epi <- ScaleData(GSM4089151_P1_2.epi, verbose = FALSE)
GSM4089151_P1_2.epi <- RunPCA(GSM4089151_P1_2.epi, npcs = 50, verbose = FALSE)
ElbowPlot(GSM4089151_P1_2.epi, ndims = 50)

#Umap and Clustering
GSM4089151_P1_2.epi <- FindNeighbors(GSM4089151_P1_2.epi, reduction = "pca", dims = 1:12)
GSM4089151_P1_2.epi <- FindClusters(GSM4089151_P1_2.epi, resolution = 0.4)
GSM4089151_P1_2.epi <- RunUMAP(GSM4089151_P1_2.epi, reduction = "pca", dims = 1:12)
DimPlot(GSM4089151_P1_2.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(GSM4089151_P1_2.epi) <- "RNA"
all.genes <- rownames(GSM4089151_P1_2.epi)
GSM4089151_P1_2.epi <- ScaleData(GSM4089151_P1_2.epi, features = all.genes)
GSM4089151_P1_2.epi <- CellCycleScoring(GSM4089151_P1_2.epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = GSM4089151_P1_2.epi) <- "Phase"
DimPlot(GSM4089151_P1_2.epi, reduction = "umap")

tiff(file = "GSM4089151_P1_2.epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2.epi, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell Cycle regression
GSM4089151_P1_2.epi2 <- GSM4089151_P1_2.epi
DefaultAssay(GSM4089151_P1_2.epi2) <- "RNA"
GSM4089151_P1_2.epi2 <- ScaleData(GSM4089151_P1_2.epi2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(GSM4089151_P1_2.epi2))
GSM4089151_P1_2.epi2 <- RunPCA(GSM4089151_P1_2.epi2, features = VariableFeatures(GSM4089151_P1_2.epi2))
ElbowPlot(GSM4089151_P1_2.epi2, ndims = 50)

GSM4089151_P1_2.epi2 <- FindNeighbors(GSM4089151_P1_2.epi2, reduction = "pca", dims = 1:16)
GSM4089151_P1_2.epi2 <- FindClusters(GSM4089151_P1_2.epi2, resolution = 0.4)
GSM4089151_P1_2.epi2 <- RunUMAP(GSM4089151_P1_2.epi2, reduction = "pca", dims = 1:16)
DimPlot(GSM4089151_P1_2.epi2, reduction = "umap", pt.size = 0.3, label = TRUE)
GSM4089151_P1_2.epi2 <- subset(GSM4089151_P1_2.epi2, idents = c("0", "1", "2", "3", "4", "5", "6", "7"))

Idents(object = GSM4089151_P1_2.epi2) <- "Phase"
tiff(file = "GSM4089151_P1_2.epi2 after Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2.epi2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Featureplots
tiff(file = "GSM4089151_P1_2.epi2 celltype marker expression plots.tiff", width = 16, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089151_P1_2.epi2, reduction = "umap", features = c("EPCAM", "VIM", "KRT5", "KRT14", "KRT8", "KRT19",
                                                              "FBLN1", "MYH11", "PECAM1",
                                                              "RGS1", "CCL5", "CHGA"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#EpiCellTypes UMAP
tiff(file = "GSM4089151_P1_2.epi2 seurat_clusters UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2.epi2, reduction = "umap", pt.size = 0.5, label = TRUE, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

#Rename
Idents(object = GSM4089151_P1_2.epi2) <- "seurat_clusters"
GSM4089151_P1_2.epi2 <- RenameIdents(object = GSM4089151_P1_2.epi2, 
                                        '5' = "LE1", '3' = "LE2", '6' = "LE3",
                                        '1' = "LE4", '0' = "LE5", '7' = "LE6", '4' = "LE7", 
                                        '2' = "LE8")  
GSM4089151_P1_2.epi2[["EpiCellTypes"]] <- Idents(object = GSM4089151_P1_2.epi2)
Idents(object = GSM4089151_P1_2.epi2) <- "EpiCellTypes"

#EpiCellTypes UMAP
tiff(file = "GSM4089151_P1_2.epi2 EpiCellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(GSM4089151_P1_2.epi2, reduction = "umap", pt.size = 0.5, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

#Featureplots for MET
tiff(file = "GSM4089151_P1_2.epi2 MET expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089151_P1_2.epi2, reduction = "umap", features = c("MET"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

#FeaturePlots for AR signature
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))

DefaultAssay(GSM4089151_P1_2.epi2) <- "RNA"
all.genes <- rownames(GSM4089151_P1_2.epi)
GSM4089151_P1_2.epi <- ScaleData(GSM4089151_P1_2.epi, features = all.genes)

GSM4089151_P1_2.epi2 <- AddModuleScore(object = GSM4089151_P1_2.epi2, features = ARlist, name = "AR_Downstream") 

tiff(file = "GSM4089151_P1_2.epi2 AR score plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089151_P1_2.epi2, reduction = "umap", features = "AR_Downstream1", cols = c("lightgrey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.3)
dev.off()

#FeaturePlots for NE signature
NElist <- list(c("SYP", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
GSM4089151_P1_2.epi2 <- AddModuleScore(object = GSM4089151_P1_2.epi2, features = NElist, name = "NE_module") 

DefaultAssay(GSM4089151_P1_2.epi2) <- "RNA"
tiff(file = "GSM4089151_P1_2.epi2 NE score plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(GSM4089151_P1_2.epi2, reduction = "umap", features = "NE_module1", cols = c("lightgrey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 0.3)
dev.off()

#AR signature vlnplot
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
GSM4089151_P1_2.epi2 <- AddModuleScore(object = GSM4089151_P1_2.epi2, features = ARlist, name = "ARscore")
Idents(object = GSM4089151_P1_2.epi2) <- "EpiCellTypes"
VlnPlot(GSM4089151_P1_2.epi2, features = "ARscore1", pt.size = 0.01) 

tiff(file = "GSM4089151_P1_2.epi2 AR signature Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "ARscore1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#NE signature vlnplot
NElist <- list(c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
GSM4089151_P1_2.epi2 <- AddModuleScore(object = GSM4089151_P1_2.epi2, features = NElist, name = "NEscore")
Idents(object = GSM4089151_P1_2.epi2) <- "EpiCellTypes"
VlnPlot(GSM4089151_P1_2.epi2, features = "NEscore1", pt.size = 0.01)

tiff(file = "GSM4089151_P1_2.epi2 NE signature Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "NEscore1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#DEGs EpiCellTypes
DefaultAssay(GSM4089151_P1_2.epi2) <- "RNA"
Idents(object = GSM4089151_P1_2.epi2) <- "EpiCellTypes"
all.genes <- rownames(GSM4089151_P1_2.epi2)
GSM4089151_P1_2.epi2 <- ScaleData(GSM4089151_P1_2.epi2, features = all.genes)
GSM4089151_P1_2.epi2.EpiCellTypes.markers <- FindAllMarkers(GSM4089151_P1_2.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(GSM4089151_P1_2.epi2.EpiCellTypes.markers, file = "GSM4089151_P1_2.epi2.EpiCellTypes.markers.csv")

#Violin for selected WNT downstream
DefaultAssay(GSM4089151_P1_2.epi2) <- "RNA"
Idents(object = GSM4089151_P1_2.epi2) <- "EpiCellTypes"
tiff(file = "GSM4089151_P1_2.epi2 AR Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "AR", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 MET Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "MET", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 CTNNB1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "CTNNB1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 CLDN1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "CLDN1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 DKK1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "DKK1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 PLAU Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "PLAU", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 SERPINB1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "SERPINB1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 MMP7 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "MMP7", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 CCND1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "CCND1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 MYC Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "MYC", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "GSM4089151_P1_2.epi2 XPO1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "XPO1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#StackedViolin for selected WNT downstream
Idents(object = GSM4089151_P1_2.epi2) <- "EpiCellTypes"

features <- c("AR", "MET", "CTNNB1", "CLDN1",
              "DKK1", "PLAU",  
              "MMP7", 'CCND1', "SERPINB1", "MYC",
              "XPO1")

b <- VlnPlot(GSM4089151_P1_2.epi2,  features, stack = TRUE, flip = TRUE) +
  theme(legend.position = "none")
tiff(file = "GSM4089151_P1_2.epi2 AR MET Wnt XPO1 StackedVln split.tiff", width = 6, height = 7, units = "in", compression = "lzw", res = 800)
plot_grid(b)
dev.off()

#Correlation_ARPC_final
Idents(object = GSM4089151_P1_2.epi2) <- "CRPC"
DimPlot(GSM4089151_P1_2.epi2, reduction = "umap", pt.size = 0.3)

GOI <- c("MET", "HGF", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4",
         "CTNNB1", "CLDN1", "DKK1", "PLAU", "SERPINB1", "MMP7", "CCND1", "XPO1")  
GOI <- c("MET", "KLK3", "KLK2", "PLPP1", "PDLIM5", "AR", 
         "DKK1", "PLAU", "CCND1", "MMP7", "CLDN1")  
GOI_index <- is.element(rownames(GSM4089151_P1_2.epi2),GOI)
Cell_index <- is.element(Idents(GSM4089151_P1_2.epi2), c('CRPC'))
expr_GOI <- GSM4089151_P1_2.epi2@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- GSM4089151_P1_2.epi2@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR WNT downstream spearman correlation CRPC.tiff", width = 7, height = 7, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()


