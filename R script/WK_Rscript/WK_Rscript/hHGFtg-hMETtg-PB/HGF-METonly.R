#hHGFtg-hMETtg

#Add necessary tools to library
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

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb")

####Loading data####
HGFMETunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/201120_TGen/RUN_mm10_hMETtg_hHGFtg/count_39553_mm10_hMETtg_hHGFtg/outs/filtered_feature_bc_matrix")
HGFMETunfiltered <- CreateSeuratObject(counts = HGFMETunfiltered.data,  min.cells = 3, min.features = 200, project = "HGF-MET")
HGFMETunfiltered <- NormalizeData(HGFMETunfiltered)

####Initial processing, Filtering and Clustering####
#HGFMETunfiltered
HGFMETunfiltered[["percent.mt"]] <- PercentageFeatureSet(HGFMETunfiltered, pattern = "^mt-")
tiff(file = "HGFMET Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGFMETunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "HGFMET Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMETunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "HGFMET Pre-filteration")
dev.off()
tiff(file = "HGFMET Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMETunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "HGFMET Pre-filteration")
dev.off()
#HGFMETfiltered
HGFMET <- subset(HGFMETunfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & percent.mt < 10)
tiff(file = "HGFMET Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGFMET, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "HGFMET Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMET@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "HGFMET Post-filteration")
dev.off()
tiff(file = "HGFMET Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGFMET@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "HGFMET Post-filteration")
dev.off()

#Genes and UMI counts per cell
mean(HGFMET$nCount_RNA)
mean(HGFMET$nFeature_RNA)

#Clustering
HGFMET <- FindVariableFeatures(HGFMET, selection.method = "vst", nfeatures = 5000)
tiff(file = "HGFMET Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(HGFMET)
dev.off()

HGFMET <- ScaleData(HGFMET, verbose = FALSE)
HGFMET <- RunPCA(HGFMET, npcs = 50, verbose = FALSE)
tiff(file = "HGFMET ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(HGFMET, ndims = 50)
dev.off()

HGFMET <- FindNeighbors(HGFMET, reduction = "pca", dims = 1:20)
HGFMET <- FindClusters(HGFMET, resolution = 0.5)
HGFMET <- RunTSNE(HGFMET, reduction = "pca", dims = 1:20)
HGFMET <- RunUMAP(HGFMET, reduction = "pca", dims = 1:20)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3) 

#Cell cycle assignment
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(HGFMET) <- "RNA"
all.genes <- rownames(HGFMET)
HGFMET <- ScaleData(HGFMET, features = all.genes)
HGFMET <- CellCycleScoring(HGFMET, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = HGFMET) <- "Phase"
DimPlot(HGFMET, reduction = "umap")
tiff(file = "HGFMET Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

Idents(object = HGFMET) <- "stim"
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3, cols = "grey")
tiff(file = "HGFMET grey UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.3, cols = "grey")
dev.off()

#
Idents(object = HGFMET) <- "stim"
DimPlot(HGFMET, reduction = "umap", pt.size = 0.5, cols = "grey")
tiff(file = "HGFMET grey 0.5 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMET, reduction = "umap", pt.size = 0.5, cols = "grey")
dev.off()

#Cell Cycle regression
HGF_MET1 <- HGFMET
DefaultAssay(HGF_MET1) <- "RNA"
HGF_MET1 <- ScaleData(HGF_MET1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(HGF_MET1))
HGF_MET1 <- RunPCA(HGF_MET1, features = VariableFeatures(HGF_MET1))
ElbowPlot(HGF_MET1, ndims = 50)

HGF_MET1 <- FindNeighbors(HGF_MET1, reduction = "pca", dims = 1:20)
HGF_MET1 <- FindClusters(HGF_MET1, resolution = 1.5)
HGF_MET1 <- RunUMAP(HGF_MET1, reduction = "pca", dims = 1:20)
HGF_MET1 <- RunTSNE(HGF_MET1, reduction = "pca", dims = 1:20)

DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = HGF_MET1) <- "Phase"
tiff(file = "HGF_MET1 Cell Cyle UMAP after Cell Cycle Regression.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#
Idents(object = HGF_MET1) <- "Phase"
tiff(file = "HGF_MET1 Cell Cyle 0.5 UMAP after Cell Cycle Regression.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.5, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#FeaturePlots for cell markers
DefaultAssay(HGF_MET1)<-"RNA"
tiff(file = "HGF_MET1 celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hHGFtg", "hMETtg", "Ar", "Pbsn",
                                                      "Krt5", "Trp63", "Krt8", "Cd24a", 
                                                      "Fbln1", "Myh11", "Plp1", "Pecam1",
                                                      "Rgs5", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell Type Identification
#Rename CellTypes
Idents(object = HGF_MET1) <- "seurat_clusters"
HGF_MET1 <- RenameIdents(object = HGF_MET1, 
                        '23'="BE", '24' = "BE", '0'= "BE", '4'="BE",
                         '21'="LE", '10'="LE",'5'="LE",'1'="LE",
                        '3'="LE",'2'="LE",'20'="LE", '7'="LE", '16'="LE", '6'="LE", '11'="LE", '13'="LE",'19'="SV",
                                        '15'="FB", '8'="FB",'17'="SM", '22'="Pericyte",
                        '25'="Glia",'14'="VE", '12'="Immune", '18'="Immune",
                                        '9'="Immune")  
HGF_MET1[["CellTypes"]] <- Idents(object = HGF_MET1)

tiff(file = "HGF_MET1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()

#
Idents(object = HGF_MET1) <- "CellTypes"
tiff(file = "HGF_MET1 CellTypes 0.5 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.5, cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()

#Cell counts
Idents(object = HGF_MET1) <- "CellTypes"
table(Idents(HGF_MET1))

#DEGs for Dotplot
#Degs celltype clusters
DefaultAssay(HGF_MET1) <- "RNA"
Idents(object = HGF_MET1) <- "CellTypes"
all.genes <- rownames(HGF_MET1)
HGF_MET1 <- ScaleData(HGF_MET1, features = all.genes)
HGF_MET1.allmarkers <- FindAllMarkers(HGF_MET1, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(HGF_MET1.allmarkers, file = "HGF_MET1.allmarkers.csv")

#Dotplots
Idents(object = HGF_MET1) <- "CellTypes"
tiff(file = "HGF_MET1 markers DotPlot.tiff", width =12 , height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(HGF_MET1, features = c("hMETtg", "hHGFtg",
  "Krt14", "Lgals7", "Krt5", "Aqp3", "Col17a1",
                                              "Gm5615", "Agr2", "5430419D17Rik", "Oit1", "Azgp1",
                                              "Defb42", "Elf5", "Hoxd8", "Svs4", "Pate4",
                                              "Apod", "Fbln1", "Dcn", "Crispld2", "Penk",
                                              "1500015O10Rik", "Dkk2", "Itgbl1", "Mfap5", "Actg2",
                                              "Rgs5", "Ndufa4l2", "Vtn", "Cox4i2", "Kcnj8",
                               "Plp1", "Fabp7", "Cdh19", "Kcna1", "Gfra3", 
                               "Flt1", "Plvap", "Aqp1", "Pecam1", "Cdh5", 
                               "Ccl5", "Rgs1", "Il1b", "Fcer1g", "Tyrobp"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Featureplots
DefaultAssay(HGF_MET1) <- "RNA"
tiff(file = "HGF_MET1 hHGFtg 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hHGFtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 hMETtg 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Pbsn 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Ar 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Epcam 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Vim 3.0 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()

#
DefaultAssay(HGF_MET1) <- "RNA"
tiff(file = "HGF_MET1 hHGFtg 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hHGFtg"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 hMETtg 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Pbsn 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Ar 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Epcam 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1 Vim 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()

####Re-clustering Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/Epi")

Idents(object = HGF_MET1) <- "CellTypes"
tiff(file = "HGF_MET1 Epi Highlight UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "light grey", "light grey", "light grey", "light grey"
                                                              , "light grey", "light grey", "light grey"))
dev.off()

Idents(object = HGF_MET1) <- "CellTypes"
HGF_MET1.epi <- subset(HGF_MET1, idents = c("BE","LE"))

#Run the standard workflow for visualization and clustering
Idents(object = HGF_MET1.epi) <- "seurat_clusters"
DefaultAssay(HGF_MET1.epi) <- "RNA"
HGF_MET1.epi <- ScaleData(HGF_MET1.epi, verbose = FALSE)
HGF_MET1.epi <- RunPCA(HGF_MET1.epi, npcs = 50, verbose = FALSE)
ElbowPlot(HGF_MET1.epi, ndims = 50)

#Umap and Clustering
HGF_MET1.epi <- FindNeighbors(HGF_MET1.epi, reduction = "pca", dims = 1:18)
HGF_MET1.epi <- FindClusters(HGF_MET1.epi, resolution = 0.5)
HGF_MET1.epi <- RunUMAP(HGF_MET1.epi, reduction = "pca", dims = 1:18)
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle assignment epi
DefaultAssay(HGF_MET1.epi) <- "RNA"
all.genes <- rownames(HGF_MET1.epi)
HGF_MET1.epi <- ScaleData(HGF_MET1.epi, features = all.genes)
HGF_MET1.epi <- CellCycleScoring(HGF_MET1.epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = HGF_MET1.epi) <- "Phase"
DimPlot(HGF_MET1.epi, reduction = "umap")
tiff(file = "HGF_MET1.epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

Idents(object = HGF_MET1.epi) <- "stim"
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, cols = "grey")
tiff(file = "HGF_MET1.epi grey UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.3, cols = "grey")
dev.off()

#
Idents(object = HGF_MET1.epi) <- "stim"
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.5, cols = "grey")
tiff(file = "HGF_MET1.epi grey 0.5 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi, reduction = "umap", pt.size = 0.5, cols = "grey")
dev.off()

#Cell Cycle Regression epi
HGF_MET1.epi1 <- HGF_MET1.epi
DefaultAssay(HGF_MET1.epi1) <- "RNA"
HGF_MET1.epi1 <- ScaleData(HGF_MET1.epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(HGF_MET1.epi1))
HGF_MET1.epi1 <- RunPCA(HGF_MET1.epi1, features = VariableFeatures(HGF_MET1.epi1))
ElbowPlot(HGF_MET1.epi1, ndims = 30)

HGF_MET1.epi1 <- FindNeighbors(HGF_MET1.epi1, reduction = "pca", dims = 1:15)
HGF_MET1.epi1 <- FindClusters(HGF_MET1.epi1, resolution = 0.3)
HGF_MET1.epi1 <- RunUMAP(HGF_MET1.epi1, reduction = "pca", dims = 1:15)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = HGF_MET1.epi1) <- "Phase"
tiff(file = "HGF_MET1.epi1 Cell Cyle UMAP after regression.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#
Idents(object = HGF_MET1.epi1) <- "Phase"
tiff(file = "HGF_MET1.epi1 0.5 Cell Cyle UMAP after regression.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.5, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Featureplots
DefaultAssay(HGF_MET1.epi1) <- "RNA"
tiff(file = "HGF_MET1.epi1 hHGFtg expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("hHGFtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 hMETtg expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Pbsn expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Ar expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Krt5 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Krt8 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()

#Featureplots for Wnt genes
tiff(file = "HGF_MET1.epi1 Ppp1r1b expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Ppp1r1b"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Ccnd1 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Ccnd1"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Tcf7l2 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Tcf7l2"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Mmp7 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Mmp7"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Sox9 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Sox9"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Plaur expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Plaur"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()
tiff(file = "HGF_MET1.epi1 Cd44 expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Cd44"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90")
dev.off()

#
DefaultAssay(HGF_MET1.epi1) <- "RNA"
tiff(file = "HGF_MET1.epi1 hHGFtg 3.0 expression plots size 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("hHGFtg"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 hMETtg 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Pbsn 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Ar 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Krt5 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()
tiff(file = "HGF_MET1.epi1 Krt8 3.0 expression plots 0.5.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HGF_MET1.epi1, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = 3.0)
dev.off()

#Rename
Idents(object = HGF_MET1.epi1) <- "seurat_clusters"
tiff(file = "HGF_MET1.epi1 UMAP seurat plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

HGF_MET1.epi1 <- RenameIdents(object = HGF_MET1.epi1, '0' = "BE", '5' = "LE1", '4' = "LE2", '3' = "LE3", '11' = "LE4", 
                                     '2' = "LE5", '10' = "LE5", '1' = "LE6", '8' = "LE7", '6' = "LE8", '7' = "UrLE", '9' = "OE")
HGF_MET1.epi1[["EpiCellTypes"]] <- Idents(object = HGF_MET1.epi1)

#Umap
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
tiff(file = "HGF_MET1.epi1 EpiCellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3, cols = c("red", "#FF9933", "#1D762E", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey"))
dev.off()

#
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
tiff(file = "HGF_MET1.epi1 EpiCellTypes 0.5 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.5, cols = c("red", "#FF9933", "#1D762E", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9","#FF5CFF", "grey"))
dev.off()

#Cell counts
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
table(Idents(HGF_MET1.epi1))

#DEGs_allclusters
DefaultAssay(HGF_MET1.epi1) <- "RNA"
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
HGF_MET1.epi1 <- ScaleData(HGF_MET1.epi1, features = rownames(HGF_MET1.epi1))
HGF_MET1.epi1.allMarkers <- FindAllMarkers(HGF_MET1.epi1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(HGF_MET1.epi1.allMarkers, "HGF_MET1.epi1.allMarkers.csv")

#Dotplots
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
tiff(file = "HGF_MET1.epi1 EpiCellType markers DotPlot.tiff", width =12 , height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(HGF_MET1.epi1, features = c("hMETtg", "hHGFtg", "Krt14", "Krt17", "Lgals7", "Krt5", "Col17a1", 
                                    "Gm42418", "Gm26917", "Lars2", "Atf3", "Slc7a5", 
                                    "Rps26", "Serf2", "Sec61g", "Cox8a", "Pfdn5",
                                    "Tac1", "Ank", "Kcnk3", "Crabp1", "Spink5",
                                    "Laptm5", "Tnfrsf9", "Ptprc", "Coro1a", "Rac2",
                                    "Msmb", "Chodl", "Cmbl", "Serpinb11", "Reg3g",
                                    "Areg", "Btc", "Ctse", "Podxl", "Nlrp10",
                                    "Cmpk2", "Apol9a", "Mx1", "Gbp2", "Igtp", 
                                    "Fgl1", "Chn2", "Glb1l3", "Gsdma", "Pgm2l1",
                                    "Gsdmc2", "Gsdmc3", "Naip5", "Ces1d", "Cxcl15", 
                                    "Igfbp6", "Col1a2", "Sparcl1", "Serping1", "Bgn"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

####hMETtg+ vs hMETtg-####

#Add hMETtg info
DefaultAssay(HGF_MET1.epi1) <- "RNA"
HGF_MET1.epi1.hMETtgPos <- subset(x=HGF_MET1.epi1, subset = hMETtg > 0)
HGF_MET1.epi1.hMETtgNeg <- subset(x=HGF_MET1.epi1, subset = hMETtg == 0)
Idents(object = HGF_MET1.epi1.hMETtgPos) <- "hMETtgPos"
Idents(object = HGF_MET1.epi1.hMETtgNeg) <- "hMETtgNeg"
HGF_MET1.epi1.hMETtgPos[["hMETtgExp"]] <- Idents(object = HGF_MET1.epi1.hMETtgPos)
HGF_MET1.epi1.hMETtgNeg[["hMETtgExp"]] <- Idents(object = HGF_MET1.epi1.hMETtgNeg)
HGF_MET1.epi1.hMETtg <- merge(x = HGF_MET1.epi1.hMETtgPos, y = HGF_MET1.epi1.hMETtgNeg)
Idents(object = HGF_MET1.epi1.hMETtg) <- "hMETtgExp"
HGF_MET1.epi1$hMETtgExp <- Idents(object = HGF_MET1.epi1.hMETtg)

#DEGs
Idents(object = HGF_MET1.epi1) <- "hMETtgExp"
tiff(file = "HGF_MET1.epi1 hMETtgExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1.epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
HGF_MET1.LE <- subset(HGF_MET1.epi1, idents = c("LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7", "LE8"))

DefaultAssay(HGF_MET1.LE) <- "RNA"
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
HGF_MET1.LE <- ScaleData(HGF_MET1.LE, features = rownames(HGF_MET1.LE))
HGF_MET1.LE.0.Markers <- FindMarkers(HGF_MET1.LE, ident.1 = "hMETtgPos", ident.2 = "hMETtgNeg", min.pct = 0, logfc.threshold = 0)
write.csv(HGF_MET1.LE.0.Markers, "HGF_MET1.LE.0.Markers.csv")

#p.adjust
DEG_hMETtgPosvhMETtgNeg <- read.csv("HGF_MET1.LE.0.Markers.csv") 
DEG_hMETtgPosvhMETtgNeg_pvalue <- DEG_hMETtgPosvhMETtgNeg$p_val
DEG_hMETtgPosvhMETtgNeg_pvalue=as.numeric(DEG_hMETtgPosvhMETtgNeg_pvalue)
DEG_hMETtgPosvhMETtgNeg_BH = p.adjust(DEG_hMETtgPosvhMETtgNeg_pvalue, "BH")
write.csv(DEG_hMETtgPosvhMETtgNeg_BH, "DEG_hMETtgPosvhMETtgNeg_BH-1.csv")

#Heatmap
#Heatmap
DefaultAssay(HGF_MET1.LE) <- "RNA"
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
HGF_MET1.LE <- ScaleData(HGF_MET1.LE, features = rownames(HGF_MET1.LE))
HGF_MET1.LE.all.markers <- FindAllMarkers(HGF_MET1.LE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HGF_MET1.LE.all.markers.Top50 <- HGF_MET1.LE.all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "HGF_MET1.LE Heatmap Top50 purple.tiff", width = 6, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(HGF_MET1.LE, features = c(HGF_MET1.LE.all.markers.Top50$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violin Plots
HGF_MET1.LE <- RenameIdents(object = HGF_MET1.LE, 'hMETtgNeg' = "hMETtgNeg", 'hMETtgPos' = "hMETtgPos")
HGF_MET1.LE[["hMETtgExp"]] <- Idents(object = HGF_MET1.LE)

Idents(object = HGF_MET1.LE) <- "hMETtgExp"
tiff(file = "HGF_MET1.LE hMETtg Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "hMETtg", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE hHGFtg Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "hHGFtg", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Plaur Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Plaur", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Sox9 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Sox9", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Mmp7 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Mmp7", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Cd44 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Cd44", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Tcf7l2 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.LE, features = "Tcf7l2", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()

#Boxplots
boxdata = FetchData(HGF_MET1.LE, c("hMETtgExp", "hMETtg", 'hHGFtg', "Mmp7", "Cd44", "Sox9", "Plaur", "Tcf7l2", "Ccnd1"))
tail(boxdata,9)

tiff(file = "HGF_MET1.LE hHGFtg Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=hHGFtg, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE hMETtg Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=hMETtg, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Mmp7 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=Mmp7, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Cd44 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=Cd44, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Sox9 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=Sox9, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Plaur Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=Plaur, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.LE Tcf7l2 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=hMETtgExp, y=Tcf7l2, fill = hMETtgExp)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()

#Gene-Gene Spearman correlation
Idents(object = HGF_MET1.LE) <- "hMETtgExp"
GOI <- c('hMETtg','hHGFtg', 'Plaur','Sox9', 'Mmp7', 'Cd44', 'Tcf7l2')  
GOI_index <- is.element(rownames(HGF_MET1.LE),GOI)
Cell_index <- is.element(Idents(HGF_MET1.LE), c('hMETtgPos','hMETtgNeg'))

expr_GOI <- HGF_MET1.LE@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- HGF_MET1.LE@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI))

tiff(file = "HGF_MET1.LE spearman correlation.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "HGF_MET1.LE spearman correlation NoLabel.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))
dev.off()

####hMETtg+ Epi only####
Idents(object = HGF_MET1.epi1) <- "hMETtgExp"
HGF_MET1.epi1.hMETtgPos <- subset(HGF_MET1.epi1, idents = c("hMETtgPos"))

####hMETtg+LE5-7 vs hMETtg-LE1-2####
#Dotplots
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
tiff(file = "HGF_MET1.epi1 EpiCellType Transgene & AR & Wnt DotPlot.tiff", width =7 , height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(HGF_MET1.epi1, features = c("hMETtg", "hHGFtg", "Ar", "Pbsn", "Krt5", "Krt14", "Krt8", "Krt18", "Mki67", "Pcna",
                                    "Plaur", "Cd44", "Sox9", "Tcf7l2", "Mmp7"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Cell counts
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes"
HGF_MET1.epi1$EpiCellTypes.hMETtgExp <- paste(Idents(HGF_MET1.epi1), HGF_MET1.epi1$hMETtgExp, sep = "_")
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes.hMETtgExp"
table(Idents(HGF_MET1.epi1))

#Subset
Idents(object = HGF_MET1.epi1) <- "EpiCellTypes.hMETtgExp"
HGF_MET1.PINvNormal <- subset(HGF_MET1.epi1, idents = c("LE1_hMETtgNeg", "LE2_hMETtgNeg", "LE5_hMETtgPos", "LE6_hMETtgPos", "LE7_hMETtgPos"))

#Rename
HGF_MET1.PINvNormal <- RenameIdents(object = HGF_MET1.PINvNormal, 'LE1_hMETtgNeg' = "Normal", 'LE2_hMETtgNeg' = "Normal", 
                              'LE5_hMETtgPos' = "PIN", 'LE6_hMETtgPos' = "PIN", 'LE7_hMETtgPos' = "PIN")
HGF_MET1.PINvNormal[["PINvNormal"]] <- Idents(object = HGF_MET1.PINvNormal)

#DEGs
DefaultAssay(HGF_MET1.PINvNormal) <- "RNA"
Idents(object = HGF_MET1.PINvNormal) <- "PINvNormal"
HGF_MET1.PINvNormal <- ScaleData(HGF_MET1.PINvNormal, features = rownames(HGF_MET1.PINvNormal))
HGF_MET1.PINvNormal.0.Markers <- FindMarkers(HGF_MET1.PINvNormal, ident.1 = "PIN", ident.2 = "Normal", min.pct = 0, logfc.threshold = 0)
write.csv(HGF_MET1.PINvNormal.0.Markers, "HGF_MET1.PINvNormal.0.Markers.csv")

#p.adjust
DEG_PINvNormal <- read.csv("HGF_MET1.PINvNormal.0.Markers.csv") 
DEG_PINvNormal_pvalue <- DEG_PINvNormal$p_val
DEG_PINvNormal_pvalue=as.numeric(DEG_PINvNormal_pvalue)
DEG_PINvNormal_BH = p.adjust(DEG_PINvNormal_pvalue, "BH")
write.csv(DEG_PINvNormal_BH, "DEG_PINvNormal_BH.csv")

#Heatmap
DefaultAssay(HGF_MET1.PINvNormal) <- "RNA"
Idents(object = HGF_MET1.PINvNormal) <- "PINvNormal"
HGF_MET1.PINvNormal <- ScaleData(HGF_MET1.PINvNormal, features = rownames(HGF_MET1.PINvNormal))
HGF_MET1.PINvNormal.all.markers <- FindAllMarkers(HGF_MET1.PINvNormal, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
HGF_MET1.PINvNormal.all.markers.Top50 <- HGF_MET1.PINvNormal.all.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "HGF_MET1.PINvNormal Heatmap Top50 purple.tiff", width = 6, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(HGF_MET1.PINvNormal, features = c(HGF_MET1.PINvNormal.all.markers.Top50$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violin Plots
Idents(object = HGF_MET1.PINvNormal) <- "PINvNormal"
tiff(file = "HGF_MET1.PINvNormal hMETtg Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "hMETtg", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal hHGFtg Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "hHGFtg", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Plaur Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "Plaur", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Sox9 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "Sox9", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Mmp7 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "Mmp7", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Cd44 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "Cd44", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Tcf7l2 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "Tcf7l2", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Ccnd1 Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET1.PINvNormal, features = "Ccnd1", pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()

#Boxplots
boxdata = FetchData(HGF_MET1.PINvNormal, c("PINvNormal", "hMETtg", 'hHGFtg', "Mmp7", "Cd44", "Sox9", "Plaur", "Tcf7l2", "Ccnd1"))
tail(boxdata,9)

tiff(file = "HGF_MET1.PINvNormal hHGFtg Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=hHGFtg, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal hMETtg Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=hMETtg, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Mmp7 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=Mmp7, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Cd44 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=Cd44, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Sox9 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=Sox9, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Plaur Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=Plaur, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Tcf7l2 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=Tcf7l2, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()
tiff(file = "HGF_MET1.PINvNormal Ccnd1 Boxplot with Dots and Median.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=PINvNormal, y=Ccnd1, fill = PINvNormal)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95) + scale_fill_manual(values=c("#3399FF", "#E06666"))
dev.off()

#Gene-Gene Spearman correlation
Idents(object = HGF_MET1.PINvNormal) <- "PINvNormal"
GOI <- c('hMETtg','hHGFtg', 'Plaur','Sox9', 'Mmp7', 'Cd44', 'Tcf7l2')  
GOI_index <- is.element(rownames(HGF_MET1.PINvNormal),GOI)
Cell_index <- is.element(Idents(HGF_MET1.PINvNormal), c('PIN','Normal'))

expr_GOI <- HGF_MET1.PINvNormal@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- HGF_MET1.PINvNormal@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI))

tiff(file = "HGF_MET1.PINvNormal spearman correlation.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)
dev.off()

tiff(file = "HGF_MET1.PINvNormal spearman correlation NoLabel.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"))
dev.off()
