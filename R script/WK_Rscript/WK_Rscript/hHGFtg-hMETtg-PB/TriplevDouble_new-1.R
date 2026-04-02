#hHGFtg-hMETtg & hHGFtg-hMETtg-Bcat

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

####Merging Dataset####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble")

#HGFMETBcat_Primary
Idents(object = PrimaryvWT.combined1) <- "stim"
HGFMETBcat_Primary <- subset(PrimaryvWT.combined1, idents = c("HGFMETBcat_Primary"))

HGFMETBcat_Primary <- FindVariableFeatures(HGFMETBcat_Primary, selection.method = "vst", nfeatures = 5000)
#Clustering
HGFMETBcat_Primary <- ScaleData(HGFMETBcat_Primary, verbose = FALSE)
HGFMETBcat_Primary <- RunPCA(HGFMETBcat_Primary, npcs = 50, verbose = FALSE)
ElbowPlot(HGFMETBcat_Primary, ndims = 50)

HGFMETBcat_Primary <- FindNeighbors(HGFMETBcat_Primary, reduction = "pca", dims = 1:20)
HGFMETBcat_Primary <- FindClusters(HGFMETBcat_Primary, resolution = 0.5)
HGFMETBcat_Primary <- RunTSNE(HGFMETBcat_Primary, reduction = "pca", dims = 1:20)
HGFMETBcat_Primary <- RunUMAP(HGFMETBcat_Primary, reduction = "pca", dims = 1:20)
tiff(file = "HGFMETBcat_Primary darkblue UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGFMETBcat_Primary, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("darkblue")) + NoLegend()
dev.off()

#HGFMET
Idents(object = HGF_MET1) <- "stim"
tiff(file = "HGF_MET1 gray UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET1, reduction = "umap", pt.size = 0.3, cols = c("gray")) + NoLegend()
dev.off()

#Stash old idents
HGFMETBcat_Primary[["orig.clusters"]] <- Idents(object = HGFMETBcat_Primary)
HGF_MET1[["orig.clusters"]] <- Idents(object = HGF_MET)

#Set Current idents
Idents(object = HGFMETBcat_Primary) <- "seurat_clusters"
Idents(object = HGF_MET1) <- "seurat_clusters"

HGFMETBcat_Primary$stim <- "Triple"
HGF_MET1$stim <- "Double"

TriplevDouble.anchors <- FindIntegrationAnchors(object.list = list(HGFMETBcat_Primary, HGF_MET1), dims = 1:20)
TriplevDouble.combined <- IntegrateData(anchorset = TriplevDouble.anchors, dims = 1:20)
DefaultAssay(TriplevDouble.combined) <- "integrated"

#clustering
TriplevDouble.combined <- ScaleData(TriplevDouble.combined, verbose = FALSE)
TriplevDouble.combined <- RunPCA(TriplevDouble.combined, npcs = 30, verbose = FALSE)
ElbowPlot(TriplevDouble.combined, ndims = 50)

#Umap and Clustering
TriplevDouble.combined <- FindNeighbors(TriplevDouble.combined, reduction = "pca", dims = 1:18)
TriplevDouble.combined <- FindClusters(TriplevDouble.combined, resolution = 0.5)
TriplevDouble.combined <- RunUMAP(TriplevDouble.combined, reduction = "pca", dims = 1:18)
DimPlot(TriplevDouble.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TriplevDouble.combined) <- "stim"
TriplevDouble.combined <- RenameIdents(object = TriplevDouble.combined, 'Double' = "Double", 'Triple' = "Triple")
TriplevDouble.combined[["stim"]] <- Idents(object = TriplevDouble.combined)

tiff(file = "TriplevDouble.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined, reduction = "umap", pt.size = 0.3, cols = c("grey", "darkblue"))
dev.off()

#Cell cycle assignment
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(TriplevDouble.combined) <- "RNA"
all.genes <- rownames(TriplevDouble.combined)
TriplevDouble.combined <- ScaleData(TriplevDouble.combined, features = all.genes)
TriplevDouble.combined <- CellCycleScoring(TriplevDouble.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TriplevDouble.combined) <- "Phase"
DimPlot(TriplevDouble.combined, reduction = "umap")
tiff(file = "TriplevDouble.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell Cycle regression
TriplevDouble.combined1 <- TriplevDouble.combined
DefaultAssay(TriplevDouble.combined1) <- "integrated"
TriplevDouble.combined1 <- ScaleData(TriplevDouble.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1))
TriplevDouble.combined1 <- RunPCA(TriplevDouble.combined1, features = VariableFeatures(TriplevDouble.combined1))
ElbowPlot(TriplevDouble.combined1, ndims = 50)

TriplevDouble.combined1 <- FindNeighbors(TriplevDouble.combined1, reduction = "pca", dims = 1:20)
TriplevDouble.combined1 <- FindClusters(TriplevDouble.combined1, resolution = 1.5)
TriplevDouble.combined1 <- RunUMAP(TriplevDouble.combined1, reduction = "pca", dims = 1:20)
TriplevDouble.combined1 <- RunTSNE(TriplevDouble.combined1, reduction = "pca", dims = 1:20)

Idents(object = TriplevDouble.combined1) <- "seurat_clusters"
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, split.by='stim',label = TRUE)

Idents(object = TriplevDouble.combined1) <- "Phase"
tiff(file = "TriplevDouble.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

Idents(object = TriplevDouble.combined1) <- "stim"
tiff(file = "TriplevDouble.combined1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("grey", "darkblue"))
dev.off()

#Cell count
Idents(object = TriplevDouble.combined1) <- "stim"
table(Idents(TriplevDouble.combined1))

#Cell type identification
DefaultAssay(TriplevDouble.combined1)<-"RNA"
tiff(file = "TriplevDouble.combined1 celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("hHGFtg", "hMETtg", "Ar", "Pbsn",
                                                                      "Krt5", "Krt14", "Krt8", "Cd24a", "Plp1",
                                                                      "Fbln1", "Myh11", "Pecam1",
                                                                      "Rgs5", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Featureplots
DefaultAssay(TriplevDouble.combined1)<-"RNA"
tiff(file = "TriplevDouble.combined1 hHGFtg split plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("hHGFtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 hMETtg split plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Ar split plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Pbsn split plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Krt5 split plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Krt8 split plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90", keep.scale = "all")
dev.off()

#Rename CellTypes
Idents(object = TriplevDouble.combined1) <- "seurat_clusters"
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, label=TRUE)
TriplevDouble.combined1 <- RenameIdents(object = TriplevDouble.combined1, 
                                        '29'="BE",'2'="BE",'4'="BE",'21'="BE", '25'="BE", 
                                        '8'="LE", '3'="LE", "12"="LE", '13'="LE",
                                        '20'="LE",'26'="LE",'24'="LE", '5'="LE",'7'="LE",
                                        '17'="LE",'6'="LE",'9'="LE",'10'="LE",'0'="LE",'1'="LE",'3'="LE",
                                        '15'="SV", '28'="SV",
                                        '16'="FB",'14'="FB",'23'="SM",'27'="Pericyte", '31'="Glia", '18'="VE",
                                        '11'="Immune",'22'="Immune", '19'="Immune", '30'="Immune")  
TriplevDouble.combined1[["CellTypes"]] <- Idents(object = TriplevDouble.combined1)

#UMAP
Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()
tiff(file = "TriplevDouble.combined1 CellTypes split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, split.by='stim', cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()

#DEGs
#Degs celltype clusters
DefaultAssay(TriplevDouble.combined1) <- "RNA"
Idents(object = TriplevDouble.combined1) <- "CellTypes"
all.genes <- rownames(TriplevDouble.combined1)
TriplevDouble.combined1 <- ScaleData(TriplevDouble.combined1, features = all.genes)
TriplevDouble.combined1.allmarkers <- FindAllMarkers(TriplevDouble.combined1, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(TriplevDouble.combined1.allmarkers, file = "TriplevDouble.combined1.allmarkers.csv")

#Dotplots
Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 markers DotPlot.tiff", width =12 , height = 5, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1, features = c("hMETtg", "hHGFtg",
                                              "Krt15", "Krt14", "Krt5", "Aqp3", "Col17a1", 
                                              "Prr9", "Slc12a2", "Krt8", "Krt19", "Arl14",
                                              "Svs2", "Svs5", "Svs4", "Pate4", "Wfdc15b",
                                              "Apod", "Penk", "Crispld2", "Sult1e1", "Inmt", 
                                              "Col3a1", "Col1a1", "Tagln", "Crlf1", "Ltbp2", 
                                              "Rgs5", "Ndufa4l2", "Cox4i2", "Abcc9", "Notch3",
                                              "Plp1", "Kcna1", "Cdh19", "S100b", "Abca8a", 
                                              "Aqp1", "Plvap", "Cdh5", "Pecam1", "Cd93", 
                                              "H2-Aa", "Rgs1", "Tyrobp", "Laptm5", "Bcl2a1b"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Cell Counts
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1$stim.CellTypes <- paste(Idents(TriplevDouble.combined1), TriplevDouble.combined1$stim, sep = "_")
Idents(object = TriplevDouble.combined1) <- "stim.CellTypes"
table(Idents(TriplevDouble.combined1))

####Subcluster epi DoublevTriple####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/Epi")

Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 Epi Highlight UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "light grey", "light grey", "light grey", "light grey"
                                                                             , "light grey", "light grey", "light grey", "light grey"))
dev.off()

Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1.epi <- subset(TriplevDouble.combined1, idents = c("BE","LE"))
Idents(object = TriplevDouble.combined1.epi) <- "seurat_clusters"

#Run the standard workflow for visualization and clustering
DefaultAssay(TriplevDouble.combined1.epi) <- "integrated"
TriplevDouble.combined1.epi <- ScaleData(TriplevDouble.combined1.epi, verbose = FALSE)
TriplevDouble.combined1.epi <- RunPCA(TriplevDouble.combined1.epi, npcs = 50, verbose = FALSE)
ElbowPlot(TriplevDouble.combined1.epi, ndims = 50)

#Umap and Clustering
TriplevDouble.combined1.epi <- FindNeighbors(TriplevDouble.combined1.epi, reduction = "pca", dims = 1:19)
TriplevDouble.combined1.epi <- FindClusters(TriplevDouble.combined1.epi, resolution = 0.7)
TriplevDouble.combined1.epi <- RunUMAP(TriplevDouble.combined1.epi, reduction = "pca", dims = 1:19)
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TriplevDouble.combined1.epi) <- "seurat_clusters"

#Cell cycle assignment epi
DefaultAssay(TriplevDouble.combined1.epi) <- "RNA"
all.genes <- rownames(TriplevDouble.combined1.epi)
TriplevDouble.combined1.epi <- ScaleData(TriplevDouble.combined1.epi, features = all.genes)
TriplevDouble.combined1.epi <- CellCycleScoring(TriplevDouble.combined1.epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TriplevDouble.combined1.epi) <- "Phase"
DimPlot(TriplevDouble.combined1.epi, reduction = "umap")
tiff(file = "TriplevDouble.combined1.epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell Cycle Regression epi
TriplevDouble.combined1.epi1 <- TriplevDouble.combined1.epi
DefaultAssay(TriplevDouble.combined1.epi1) <- "integrated"
TriplevDouble.combined1.epi1 <- ScaleData(TriplevDouble.combined1.epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi1))
TriplevDouble.combined1.epi1 <- RunPCA(TriplevDouble.combined1.epi1, features = VariableFeatures(TriplevDouble.combined1.epi1))
ElbowPlot(TriplevDouble.combined1.epi1, ndims = 50)

TriplevDouble.combined1.epi1 <- FindNeighbors(TriplevDouble.combined1.epi1, reduction = "pca", dims = 1:19)
TriplevDouble.combined1.epi1 <- FindClusters(TriplevDouble.combined1.epi1, resolution = 0.4)
TriplevDouble.combined1.epi1 <- RunUMAP(TriplevDouble.combined1.epi1, reduction = "pca", dims = 1:19)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TriplevDouble.combined1.epi1) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi1 after Cell Cyle Regression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Feature plots
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("hHGFtg","hMETtg","Ar","Krt5", "Krt8","Ppp1r1b", "Svs2", "Aqp3",
                                                                           "Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename Clusters
Idents(object = TriplevDouble.combined1.epi1) <- "seurat_clusters"
TriplevDouble.combined1.epi1 <- RenameIdents(object = TriplevDouble.combined1.epi1, 
                                             '6'="BE1",'5'="BE2",'8'="BE3", 
                                             '3'="LE1",'0'="LE2", '12'="LE2", '1'="LE3", '11' = "LE4",
                                             '2'="LE4",'7'="LE5",'4'="LE6",
                                             '9'="UrLE", 
                                             '10'="OE")  
TriplevDouble.combined1.epi1[["EpiCellTypes"]] <- Idents(object = TriplevDouble.combined1.epi1)

tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label.size = 6, cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "blueviolet", "blue", "bisque3", "steelblue1", "aquamarine3","deeppink2", "plum4",
                                                                                                  "darkorange1", "grey", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, split.by='stim', cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "blueviolet", "blue", "bisque3", "steelblue1", "aquamarine3","deeppink2", "plum4",
                                                                                                   "darkorange1", "grey", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes stim Label UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6, split.by='stim', cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "blueviolet", "blue", "bisque3", "steelblue1", "aquamarine3","deeppink2", "plum4",
                                                                                                                                 "darkorange1", "grey", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

#Featureplots
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi1 hHGFtg split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("hHGFtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 hMETtg split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Ar split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 2.5, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Tmprss2 split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Tmprss2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 2.5, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Pbsn split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Krt5 split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Krt8 split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90", keep.scale = "all")
dev.off()

#Cell Counts
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi1), TriplevDouble.combined1.epi1$stim, sep = "_")
Idents(object = TriplevDouble.combined1.epi1) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi1))

#Heatmap for Triple only
Idents(object = TriplevDouble.combined1.epi1) <- "stim"
Triple.combined1.epi <- subset(TriplevDouble.combined1.epi1, idents = c("Triple"))

DefaultAssay(Triple.combined1.epi) <- "RNA"
Idents(object = Triple.combined1.epi) <- "EpiCellTypes"
all.genes <- rownames(Triple.combined1.epi)
Triple.combined1.epi <- ScaleData(Triple.combined1.epi, features = all.genes)
Triple.combined1.epi.allmarkers <- FindAllMarkers(Triple.combined1.epi, min.pct = 0.1, logfc.threshold = 0.1, only.pos=TRUE)
write.csv(Triple.combined1.epi.allmarkers, file = "Triple.combined1.epi.allmarkers.csv")

Triple.combined1.epi.allmarkers.Top10 <- Triple.combined1.epi.allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
tiff(file = "Triple.combined1.epi.allmarkers.Top10 Heatmap.tiff", width =12, height = 12, units = "in", compression = "lzw", res = 200)
DoHeatmap(Triple.combined1.epi, features = c(Triple.combined1.epi.allmarkers.Top10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####Higher resolution for TriplevDouble Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/Epi/TriplevDouble.combined1.epi2")

TriplevDouble.combined1.epi2 <- TriplevDouble.combined1.epi1
DefaultAssay(TriplevDouble.combined1.epi2) <- "integrated"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2 <- RunPCA(TriplevDouble.combined1.epi2, features = VariableFeatures(TriplevDouble.combined1.epi2))
ElbowPlot(TriplevDouble.combined1.epi2, ndims = 50)

tiff(file = "TriplevDouble.combined1.epi2 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(TriplevDouble.combined1.epi2, ndims = 50)
dev.off()

#Res0.7
TriplevDouble.combined1.epi2 <- FindNeighbors(TriplevDouble.combined1.epi2, reduction = "pca", dims = 1:19)
TriplevDouble.combined1.epi2 <- FindClusters(TriplevDouble.combined1.epi2, resolution = 0.7)
TriplevDouble.combined1.epi2 <- RunUMAP(TriplevDouble.combined1.epi2, reduction = "pca", dims = 1:19)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

#seurat_clusters Dimplot
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
tiff(file = "TriplevDouble.combined1.epi2 Dim25 res0.7 UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6, split.by = "stim")
dev.off()

#Phase Dimplot
Idents(object = TriplevDouble.combined1.epi2) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi2 after Cell Cyle Regression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Featureplots
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi2 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt5", "Krt14", "Krt19", "Krt8", "Ppp1r1b", "Krt4", "Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2.allMarkers <- FindAllMarkers(TriplevDouble.combined1.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi2.allMarkers, "TriplevDouble.combined1.epi2.seurat.allMarkers.csv")

#Rename Clusters
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- RenameIdents(object = TriplevDouble.combined1.epi2, 
                                             '5'="BE1",'7'="BE2",'17'="BE3", '11' = "BE4", '13'="BE5", 
                                             '4'="LE1",'0'="LE2", '2'="LE3", '3' = "LE4", '10'="LE5",
                                             '9'="LE6",'12'="LE7",'6'="LE8", '18' = "LE8",
                                             '8'="LE9", '1'="LE10", '15'="LE11", '14'="UrLE", 
                                             '16'="OE")  
TriplevDouble.combined1.epi2[["EpiCellTypes"]] <- Idents(object = TriplevDouble.combined1.epi2)

tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label.size = 6, cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                  "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, split.by='stim', cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                   "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

#Cell Counts
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.epi2$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi2), TriplevDouble.combined1.epi2$stim, sep = "_")
Idents(object = TriplevDouble.combined1.epi2) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi2))

Idents(object = TriplevDouble.combined1.epi2) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi2 after Cell Cyle Regression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off() 

#Heatmap 
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.epi2)
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = all.genes)
TriplevDouble.combined1.epi2.allmarkers <- FindAllMarkers(TriplevDouble.combined1.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(TriplevDouble.combined1.epi2.allmarkers, file = "TriplevDouble.combined1.epi2.allmarkers.csv")

TriplevDouble.combined1.epi2.allmarkers.Top10 <- TriplevDouble.combined1.epi2.allmarkers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
tiff(file = "TriplevDouble.combined1.epi2.allmarkers.Top10 Heatmap.tiff", width =12, height = 12, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c(TriplevDouble.combined1.epi2.allmarkers.Top10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Dotplot
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes markers DotPlot.tiff", width =20 , height = 10, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1.epi2, features = c("hMETtg", "hHGFtg", 
                                                   "Apoe", "Dcn", "Pdpn", "Fhl1", "Axl", 
                                                   "Ets1", "Krt13","Inhba", "Nrg1", "Jag1", 
                                                   "Krt75", "Fmod", "Ncam1", "Smoc2", "S100a4", 
                                                   "Ivl", "Il33", "Gm42835", "Nccrp1", "Atp6v0d2", 
                                                   "Gpnmb", "Wnt10a", "Wif1", "Wnt6", "Notum",
                                                   "Defa22", "Defa5","Cryba4", "Lypd2", "Golga7b",
                                                   "Clu", "Ly6a", "Ly6c1","Ly6e", "Psca",
                                                   "Crip1", "Upk1b", "Mrps18a", "Areg", "Serpine1", 
                                                   "Coch", "Tgfb2", "Dkk2", "Ier3", "Zbtb20", 
                                                   "9530053A07Rik", "Tgm4", "C1rb", "C1s2", "Glb1l3",
                                                   "Msmb", "Defb50", "Mme", "Cldn10", "Apof",
                                                   "Rsad2", "Ifit1", "Zbp1", "Usp18", "Xaf1",
                                                   "Pigr", "Lrg1", "5330417C22Rik", "Gne", "Adh1",
                                                   "Sbp", "Spink1", "Sbpl","Kcnk3", "Gucy2g",
                                                   "Svs5", "Gm42418", "AY036118", "Lars2", "Svs2",
                                                   "Serping1", "Igfbp6", "Serpinf1", "Apod", "Sparcl1",
                                                   "Gsdmc3", "Scnn1a", "Krt4", "Ppp1r1b", "Cbr2",
                                                   "Srgn", "Rgs1", "H2-Eb1", "H2-Aa", "Cd74"),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Featureplots
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi2 hHGFtg split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("hHGFtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 hMETtg split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Ar split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 2.5, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Tmprss2 split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Tmprss2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 2.5, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Pbsn split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0, keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Krt5 split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Krt8 split expression plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90", keep.scale = "all")
dev.off()

##Cell type identification
#Heatmap for AR downstream
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi2 AR downstream Heatmap.tiff", width = 10, height = 20, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("hHGFtg", "hMETtg", "Ar", 
                                                     "Abcc4",	"Abhd2",	"Acsl3",	"Actn1",	"Adamts1",	"Adrm1",	"Akap12",	"Akt1",	"Aldh1a3",	"Ank",	"Appbp2",	"Arid5b",	"Azgp1",
                                                     "B2m",	"B4galt1",	"Bmpr1b",	"Camkk2",	"Ccnd1",	"Ccnd3",	"Cdc14b",	"Cdk6",	
                                                     "Cenpn",	"Dbi",	"Dhcr24",	"Dnajb9",	"Elk4",	"Ell2",	"Elovl5",	"Fads1",	"Fkbp5",
                                                     "Gnai3",	"Gpd1l",	"Gsr",	"H1f0",	"Herc3",	"Hmgcr",	"Hmgcs1",
                                                     "Homer2",	"Hpgd",	"Hsd17b14",	"Idi1",	"Inpp4b",	"Insig1",	"Iqgap2",	
                                                     "Itgav",	"Lifr",	"Lman1",	"Maf",	"Mak",	"Map7",	"Mertk",	"Myl12b",
                                                     "Ncoa4",	"Ndrg1",	"Ngly1",	"Nkx3-1",	"Pa2g4",	"Pdlim5",
                                                     "Pgm3",	"Pias1",	"Plpp1",	"Pmepa1",	"Ptk2b",	"Ptpn21",	
                                                     "Rab4a",	"Rps6ka3",	"Rrp12",	"Sat1",	"Scd1",	"Sec24d",	"Sgk1",
                                                     "Slc26a2",	"Slc38a2",	"Sms",	"Sord",	"Spcs3",	"Spdef",	"Srf",	"Srp19",	
                                                     "Steap4",	"Stk39",	"Tmem50a",	"Tmprss2",	"Tnfaip8",	"Tpd52",	"Tsc22d1",
                                                     "Uap1",	"Ube2i",	"Ube2j1",	"Vapa",	"Xrcc5",	"Xrcc6",	"Zbtb10",	"Zmiz1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 Wnt downstream Heatmap.tiff", width = 10, height = 15, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("hHGFtg", "hMETtg", "Ar", 
                                                     "Atoh1",	"Axin2",	"Birc5",	"Ccnd1",
                                                     "Cd44",	"Cdx4",	"Cldn1",	"Ctla4",	"Dkk1",	"Edn1",
                                                     "Efnb1",	"Egfr",	"En2",	"Fgf18",	"Fgf20",	"Fgf9",
                                                     "Fosl1",	"Fst",	"Fzd7",	"Gast",	"Gja1",	"Id2",	"Jag1",
                                                     "Jun",	"L1cam",	"Lef1",	"Lgr5",	"Met",	"Mmp7",	"Myc",
                                                     "Mycbp",	"Nrcam",	"Plaur",	"Ppard",	"Pttg1",	"Rhou",
                                                     "Sfrp2",	"Sox17",	"Stra6",	"Tcf20",	"Tcf4",	"Tcf7",
                                                     "Tcf7l2",	"Tiam1",	"Tnfrsf19",	"Twist1",	"Vcan",	"Vegfa")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 Biocarta MET downstream Heatmap.tiff", width = 10, height = 20, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("hHGFtg", "hMETtg", "Ar", 
                                                     "Arf6",	"Cbl",	"Col11a1",	"Col11a2",	"Col1a1",	"Col1a2",
                                                     "Col24a1",	"Col27a1",	"Col2a1",	"Col3a1",	"Col5a1",	"Col5a2",
                                                     "Col5a3",	"Crk",	"Crkl",	"Dock7",	"Eps15",	"Fn1",	"Gab1",
                                                     "Gga3",	"Grb2",	"Hgf",	"Hgfac",	"Hgs",	"Hpn",	
                                                     "Hras",	"Itga2",	"Itga3",	"Itgb1",	"Kras",
                                                     "Lama1",	"Lama2",	"Lama3",	"Lama4",	"Lama5",
                                                     "Lamb1",	"Lamb2",	"Lamb3",	"Lamc1",	"Lamc2",	
                                                     "Lamc3",	"Lrig1",	"Met",	"Muc20",	"Nras",	"Pik3ca",
                                                     "Pik3r1",	"Ptk2",	"Ptpn1",	"Ptpn11",	"Ptpn2",
                                                     "Ptprj",	"Rab4a",	"Rab4b",	"Rac1",
                                                     "Ranbp10",	"Ranbp9",	"Rap1a",	"Rap1b",	"Rapgef1",	"Rps27a",
                                                     "Sh3gl1",	"Sh3gl2",	"Sh3gl3",	"Sh3kbp1",	"Shc1",	"Sos1",
                                                     "Spint1",	"Spint2",	"Src",	"Stam",	"Stam2",	"Stat3",
                                                     "Tns3",	"Tns4",	"Uba52",	"Ubb",	"Ubc",	"Usp8")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.epi2)
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = all.genes)

tiff(file = "TriplevDouble.combined1.epi2 AR WNT MET downstream Heatmap.tiff", width = 10, height = 6, units = "in", compression = "lzw", res =200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("hMETtg", "Ar", 
                                                     "Ank", "Bmpr1b", "Fkbp5",
                                                     "H1f0", "Tmem50a", "Homer2",
                                                     "Hmgcr", "Hmgcs1", "Tmprss2",
                                                     "Axin2", "Lef1", "Tnfrsf19", "Fst", "Lgr5", "Tcf4",
                                                     "Grb2", "Kras", "Lamb2", "Pik3ca", "Rac1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 AR WNT MET downstream Heatmap-1.tiff", width = 10, height = 6, units = "in", compression = "lzw", res =200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("hMETtg", "Ar", 
                                                     "Ank", "Scd1", "Bmpr1b", "Fkbp5",
                                                     "Hpgd", "Homer2", "Aldh1a3",
                                                     "Hmgcr", "Hmgcs1", "Tmprss2", "Dhcr24",
                                                     "Axin2", "Lef1", "Tnfrsf19", "Stra6", "Birc5", 
                                                     "Fst", "Jun", "Mmp7", "Myc", "Tcf4",
                                                     "Grb2", "Kras", "Lamc1", "Pik3ca", "Rac1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Vlnplots for AR signature
ARlist <- list(c("Ank", "Scd1", "Bmpr1b", "Fkbp5", "Hpgd", "Homer2", "Aldh1a3", "Hmgcr", "Hmgcs1", "Tmprss2", "Dhcr24"))
TriplevDouble.combined1.epi2 <- AddModuleScore(object = TriplevDouble.combined1.epi2, features = ARlist, name = "AR_Downstream") 
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
VlnPlot(TriplevDouble.combined1.epi2, features = "AR_Downstream1", split.by = "stim", pt.size = 0.01)
tiff(file = "TriplevDouble.combined1.epi2 AR signature Vln.tiff", width = 8, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi2, features = "AR_Downstream1", split.by = "stim", pt.size = 0.01)
dev.off()

WNTlist <- list(c("Axin2", "Lef1", "Tnfrsf19", "Stra6", "Birc5", "Fst", "Jun", "Mmp7", "Myc", "Tcf4"))
TriplevDouble.combined1.epi2 <- AddModuleScore(object = TriplevDouble.combined1.epi2, features = ARlist, name = "AR_Downstream") 
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
VlnPlot(TriplevDouble.combined1.epi2, features = "AR_Downstream1", split.by = "stim", pt.size = 0.01)
tiff(file = "TriplevDouble.combined1.epi2 AR signature Vln.tiff", width = 8, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi2, features = "AR_Downstream1", split.by = "stim", pt.size = 0.01)
dev.off()

####epi3####
TriplevDouble.combined1.epi3 <- TriplevDouble.combined1.epi2
DefaultAssay(TriplevDouble.combined1.epi3) <- "integrated"
TriplevDouble.combined1.epi3 <- ScaleData(TriplevDouble.combined1.epi3, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi3))
TriplevDouble.combined1.epi3 <- RunPCA(TriplevDouble.combined1.epi3, features = VariableFeatures(TriplevDouble.combined1.epi3))
ElbowPlot(TriplevDouble.combined1.epi3, ndims = 50)
#Res0.7
TriplevDouble.combined1.epi3 <- FindNeighbors(TriplevDouble.combined1.epi3, reduction = "pca", dims = 1:26)
TriplevDouble.combined1.epi3 <- FindClusters(TriplevDouble.combined1.epi3, resolution = 0.7)
TriplevDouble.combined1.epi3 <- RunUMAP(TriplevDouble.combined1.epi3, reduction = "pca", dims = 1:26)
DimPlot(TriplevDouble.combined1.epi3, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

#seurat_clusters Dimplot
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
tiff(file = "TriplevDouble.combined1.epi2 Dim25 res0.7 UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6, split.by = "stim")
dev.off()

#Phase Dimplot
Idents(object = TriplevDouble.combined1.epi2) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi2 after Cell Cyle Regression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Featureplots
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi2 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt5", "Krt14", "Krt19", "Krt8", "Ppp1r1b", "Krt4", "Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2.allMarkers <- FindAllMarkers(TriplevDouble.combined1.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi2.allMarkers, "TriplevDouble.combined1.epi2.seurat.allMarkers.csv")


####Triple&Double LE####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/Epi/TriplevDouble.combined1.epi2/LE")

#Phase Dimplot
Idents(object = TriplevDouble.combined1.epi2) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi2 G2M highlight UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("red", "grey", "grey"))
dev.off()

Idents(object = TriplevDouble.combined1.epi2) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi2 cellcyleregression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey", "gold"))
dev.off()

Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.LE <- subset(TriplevDouble.combined1.epi2, idents = c("LE1", "LE2", "LE3", "LE4", "LE5",
                                                                              "LE6", "LE7", "LE8", "LE9", "LE10", "LE11"))

#Env
BiocManager::install("ComplexHeatmap")
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggplotify)

df <- read.csv("Heatmap_Triple_and_Double_LE_EpiCellTypes_1.csv", header = TRUE, sep = ",")
colnames(df)[1] <- "Clusters"
df <- as.data.frame(df)

#
df <- column_to_rownames(df, var = "Clusters")
df <- as.matrix(df)
d <- pheatmap::pheatmap(df, cluster_cols = F, cluster_rows = F)

##Heatmap1
RColorBrewer::brewer.pal(n=3,name = "PuBlye")
df <- column_to_rownames(df, var = "Clusters")
df <- as.matrix(df)

plt2 <- pheatmap(df, cluster_cols = F,cluster_rows = F,
                 color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                 cellwidth = 20,fontsize = 10,angle_col = 45)

tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes Heatmap G2M MET WNT AR Ribosome-3.tiff", width = 12, height = 8, units = "in", compression = "lzw", res = 800)
plt2
dev.off()

##Heatmap2
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
# Change to log10
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)

tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes Heatmap G2M MET WNT AR Ribosome-1.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)
dev.off()

#Heatmap for NE, AR, SQUAM markers in BELE
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.epi2)
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = all.genes)

tiff(file = "TriplevDouble.combined1.epi2 AR NE SQUAM Heatmap.tiff", width =12, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("Chga",	"Syp",	"Actl6b",	"Snap25",	"Insm1",	"Ascl1",	"Chrnb2",	"Srrm4",	"Celf3",
                                                   "Pcsk1",	"Sox2",	"Pou3f2",	"Lmo3",	"Nkx2-1",	"Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                   "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                   "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 AR SQUAM Heatmap.tiff", width =12, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                   "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                   "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap for NE, AR, SQUAM markers in BELE in triple only
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.epi2)
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = all.genes)

tiff(file = "TriplevDouble.combined1.epi2 AR NE SQUAM Heatmap.tiff", width =12, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("Chga",	"Syp",	"Actl6b",	"Snap25",	"Insm1",	"Ascl1",	"Chrnb2",	"Srrm4",	"Celf3",
                                                     "Pcsk1",	"Sox2",	"Pou3f2",	"Lmo3",	"Nkx2-1",	"Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                     "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                     "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 AR SQUAM Heatmap.tiff", width =12, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c("Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                     "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                     "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap for NE, AR, SQUAM markers in LE
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)

tiff(file = "TriplevDouble.combined1.LE AR NE SQUAM Heatmap.tiff", width =12, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.LE, features = c("Chga",	"Syp",	"Actl6b",	"Snap25",	"Insm1",	"Ascl1",	"Chrnb2",	"Srrm4",	"Celf3",
                                                   "Pcsk1",	"Sox2",	"Pou3f2",	"Lmo3",	"Nkx2-1",	"Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                   "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                   "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "TriplevDouble.combined1.LE AR SQUAM Heatmap.tiff", width =12, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.LE, features = c("Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                   "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                   "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap for NE, AR, SQUAM markers in triple only
Idents(object = TriplevDouble.combined1.LE) <- "stim"
Triple.combined1.LE <- subset(TriplevDouble.combined1.LE, idents = c("Triple"))

DefaultAssay(Triple.combined1.LE) <- "RNA"
Idents(object = Triple.combined1.LE) <- "EpiCellTypes"
all.genes <- rownames(Triple.combined1.LE)
Triple.combined1.LE <- ScaleData(Triple.combined1.LE, features = all.genes)

tiff(file = "Triple.combined1.LE AR NE SQUAM Heatmap.tiff", width =8, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(Triple.combined1.LE, features = c("Chga",	"Syp",	"Actl6b",	"Snap25",	"Insm1",	"Ascl1",	"Chrnb2",	"Srrm4",	"Celf3",
                                                   "Pcsk1",	"Sox2",	"Pou3f2",	"Lmo3",	"Nkx2-1",	"Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                   "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                   "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "Triple.combined1.LE AR SQUAM Heatmap.tiff", width =8, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(Triple.combined1.LE, features = c("Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                                   "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel",
                                                   "Muc4",	"Krt14",	"Anxa8",	"Sbsn"
)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#
epi2_df <- as.data.frame(t(TriplevDouble.combined1.epi2$RNA@scale.data))
epi2_df_selected <- cbind(epi2_df$Chga, epi2_df$Syp, epi2_df$Actl6b, epi2_df$Snap25, epi2_df$Insm1,
                epi2_df$Ascl1, epi2_df$Chrnb2, epi2_df$Srrm4, epi2_df$Celf3,
                epi2_df$Pcsk1, epi2_df$Sox2, epi2_df$Pou3f2, epi2_df$Lmo3, epi2_df$`Nkx2-1`,
                epi2_df$Ar, epi2_df$`Nkx3-1`, epi2_df$Klk1b24, epi2_df$Klk1b22, epi2_df$Chrna2,
                epi2_df$Slc45a3, epi2_df$Nap1l2, epi2_df$S100a14, epi2_df$Krt5, epi2_df$Krt6a,
                epi2_df$Krt6b, epi2_df$Fgfbp1, epi2_df$Dsg3, epi2_df$Scel, epi2_df$S100a7a, epi2_df$Muc4,
                epi2_df$Krt14, epi2_df$Anxa8, epi2_df$Sbsn)
epi2_df_selected <- as.data.frame(epi2_df_selected)
colnames(epi2_df_selected) <- c("Chga",	"Syp",	"Actl6b",	"Snap25",	"Insm1",	"Ascl1",	"Chrnb2",	"Srrm4",	"Celf3",
                                "Pcsk1",	"Sox2",	"Pou3f2",	"Lmo3",	"Nkx2-1",	"Ar",	"Nkx3-1",	"Klk1b24",	"Klk1b22",	"Chrna2",
                                "Slc45a3",	"Nap1l2",	"S100a14",	"Krt5",	"Krt6a",	"Krt6b",	"Fgfbp1",	"Dsg3",	"Scel", "S100a7a",
                                "Muc4",	"Krt14",	"Anxa8",	"Sbsn")
rownames(epi2_df_selected) <- row.names(epi2_df)
View(epi2_df_selected)

#Simplified Heatmap
df1 <- read.csv("df.csv", header = TRUE, sep = ",")
colnames(df1)[1] <- "Clusters"
df1 <- as.data.frame(df1)

##Heatmap1
df1 <- column_to_rownames(df1, var = "Clusters")
df1 <- as.matrix(df1)
d1 <- pheatmap::pheatmap(df1, cluster_cols = F, cluster_rows = F)

##Heatmap2
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
# Change to log10
pheatmap::pheatmap(log10(df1+1),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F) 

tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes Simplified Heatmap NE1 NE2 AR SQUAM-2.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
pheatmap::pheatmap(log2(df1+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)
dev.off()

#Simplified Heatmap-1
df2 <- read.csv("df_1.csv", header = TRUE, sep = ",")
colnames(df2)[1] <- "Clusters"
df2 <- as.data.frame(df2)

##Heatmap1
df2 <- column_to_rownames(df2, var = "Clusters")
df2 <- as.matrix(df2)
d2 <- pheatmap::pheatmap(df2, cluster_cols = F, cluster_rows = F)

##Heatmap2
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
# Change to log10
pheatmap::pheatmap(log10(df2+1),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = T, cluster_rows = F) 

tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes Simplified Heatmap NE1 NE2 AR SQUAM-2.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
pheatmap::pheatmap(log2(df2+1),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)
dev.off()



####Pseudotime analysis using hMETtg+ cells####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/Epi/TriplevDouble.combined1.epi2/hMETtgPos")

library(Seurat)
library(devtools)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 

#Add hMETtg info
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
TriplevDouble.combined1.epi2METPos <- subset(x=TriplevDouble.combined1.epi2,  subset = `hMETtg` > 0)
TriplevDouble.combined1.epi2METNeg <- subset(x=TriplevDouble.combined1.epi2,  subset = `hMETtg` == 0)
Idents(object = TriplevDouble.combined1.epi2METPos) <- "METPos"
Idents(object = TriplevDouble.combined1.epi2METNeg) <- "METNeg"
TriplevDouble.combined1.epi2METPos[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi2METPos)
TriplevDouble.combined1.epi2METNeg[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi2METNeg)
TriplevDouble.combined1.epi2MET <- merge(x = TriplevDouble.combined1.epi2METPos, y = TriplevDouble.combined1.epi2METNeg)
Idents(object = TriplevDouble.combined1.epi2MET) <- "METExp"
TriplevDouble.combined1.epi2$METExp <- Idents(object = TriplevDouble.combined1.epi2MET)
Idents(object = TriplevDouble.combined1.epi2) <- "METExp"
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label = TRUE)

TriplevDouble.combined1.epi2.METPos <- subset(TriplevDouble.combined1.epi2, idents = c("METPos"))
Idents(object = TriplevDouble.combined1.epi2.METPos) <- "EpiCellTypes"

TriplevDouble.combined1.epi2.METPos1 <- subset(TriplevDouble.combined1.epi2.METPos, idents = c("BE1", "BE2", "BE3", "BE4", "BE5",
                                                                                               "LE1", "LE2", "LE3", "LE4", "LE5",
                                                                                               "LE6", "LE7", "LE8", "LE9", "LE10",
                                                                                               "LE11"))

DefaultAssay(TriplevDouble.combined1.epi2.METPos1) <- "integrated"
TriplevDouble.combined1.epi2.METPos1 <- ScaleData(TriplevDouble.combined1.epi2.METPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi2.METPos1))
TriplevDouble.combined1.epi2.METPos1 <- RunPCA(TriplevDouble.combined1.epi2.METPos1, features = VariableFeatures(TriplevDouble.combined1.epi2.METPos1))
ElbowPlot(TriplevDouble.combined1.epi2.METPos1, ndims = 50)

#Res0.7
TriplevDouble.combined1.epi2.METPos1 <- FindNeighbors(TriplevDouble.combined1.epi2.METPos1, reduction = "pca", dims = 1:15)
TriplevDouble.combined1.epi2.METPos1 <- FindClusters(TriplevDouble.combined1.epi2.METPos1, resolution = 0.7)
TriplevDouble.combined1.epi2.METPos1 <- RunUMAP(TriplevDouble.combined1.epi2.METPos1, reduction = "pca", dims = 1:15)
DimPlot(TriplevDouble.combined1.epi2.METPos1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TriplevDouble.combined1.epi2.METPos1) <- "EpiCellTypes"
DimPlot(TriplevDouble.combined1.epi2.METPos1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim", 
        cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))

##Convert Seurat to Monocle3 cell data set class
cds <- as.cell_data_set(TriplevDouble.combined1.epi2.METPos1)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(TriplevDouble.combined1.epi2.METPos1[["RNA"]])

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 
cds <- learn_graph(cds, use_partition = FALSE)

#plot single cell trajectory
tiff(file = "TriplevDouble.combined1.epi2.METPos1 EpiCellTypes trajectory UMAP.tiff", width = 7, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = FALSE,
           label_cell_groups = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           cell_size = 0.5)+ scale_color_manual(values = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2.METPos1 BELE EpiCellTypes trajectory w Branch points UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = T, 
           label_leaves= T,
           label_roots = FALSE,
           group_label_size = 4,
           graph_label_size = 3)
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

tiff(file = "TriplevDouble.combined1.epi2.METPos1 pseudotime trajectory UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           label_roots = FALSE,
           cell_size = 0.5) 
dev.off()

tiff(file = "triplevdouble BELE pseudotime trajectory with branch points UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           graph_label_size=4)
dev.off()

###less resolution
TriplevDouble.combined1.epi2.METPos2 <- TriplevDouble.combined1.epi2.METPos1
DefaultAssay(TriplevDouble.combined1.epi2.METPos2) <- "integrated"
TriplevDouble.combined1.epi2.METPos2 <- ScaleData(TriplevDouble.combined1.epi2.METPos2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi2.METPos2))
TriplevDouble.combined1.epi2.METPos2 <- RunPCA(TriplevDouble.combined1.epi2.METPos2, features = VariableFeatures(TriplevDouble.combined1.epi2.METPos2))
ElbowPlot(TriplevDouble.combined1.epi2.METPos2, ndims = 50)

#Res0.5
TriplevDouble.combined1.epi2.METPos2 <- FindNeighbors(TriplevDouble.combined1.epi2.METPos2, reduction = "pca", dims = 1:23)
TriplevDouble.combined1.epi2.METPos2 <- FindClusters(TriplevDouble.combined1.epi2.METPos2, resolution = 0.5)
TriplevDouble.combined1.epi2.METPos2 <- RunUMAP(TriplevDouble.combined1.epi2.METPos2, reduction = "pca", dims = 1:23)
DimPlot(TriplevDouble.combined1.epi2.METPos2, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TriplevDouble.combined1.epi2.METPos2) <- "EpiCellTypes"
DimPlot(TriplevDouble.combined1.epi2.METPos2, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim", 
        cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))

##Convert Seurat to Monocle3 cell data set class
cds2 <- as.cell_data_set(TriplevDouble.combined1.epi2.METPos2)

#Construction single cell trajectories
cds2 <- cluster_cells(cds = cds2, reduction_method = "UMAP") 
cds2 <- learn_graph(cds2, use_partition = FALSE)

#plot single cell trajectory
tiff(file = "TriplevDouble.combined1.epi2.METPos1 EpiCellTypes trajectory UMAP-1.tiff", width = 7, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds2,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = FALSE,
           label_cell_groups = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           cell_size = 0.5)+ scale_color_manual(values = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2.METPos1 BELE EpiCellTypes trajectory w Branch points UMAP-1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds2,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = T, 
           label_leaves= T,
           label_roots = FALSE,
           group_label_size = 4,
           graph_label_size = 3)
dev.off()

#plot single cell trajectory by pseudotime
cds2 <- order_cells(cds2, reduction_method = "UMAP")

tiff(file = "TriplevDouble.combined1.epi2.METPos1 pseudotime trajectory UMAP-1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           label_roots = FALSE,
           cell_size = 0.5) 
dev.off()

tiff(file = "triplevdouble BELE pseudotime trajectory with branch points UMAP-1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds2,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           graph_label_size=4)
dev.off()


###Subset Triple
Idents(object = TriplevDouble.combined1.epi2.METPos1) <- "stim"
Triple.combined1.epi2.METPos <- subset(TriplevDouble.combined1.epi2.METPos1, idents = c("Triple"))

Triple.combined1.epi2.METPos1 <- Triple.combined1.epi2.METPos
DefaultAssay(Triple.combined1.epi2.METPos1) <- "integrated"
Triple.combined1.epi2.METPos1 <- ScaleData(Triple.combined1.epi2.METPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.epi2.METPos1))
Triple.combined1.epi2.METPos1 <- RunPCA(Triple.combined1.epi2.METPos1, features = VariableFeatures(Triple.combined1.epi2.METPos1))
ElbowPlot(Triple.combined1.epi2.METPos1, ndims = 50)

#Res0.7
Triple.combined1.epi2.METPos1 <- FindNeighbors(Triple.combined1.epi2.METPos1, reduction = "pca", dims = 1:16)
Triple.combined1.epi2.METPos1 <- FindClusters(Triple.combined1.epi2.METPos1, resolution = 0.7)
Triple.combined1.epi2.METPos1 <- RunUMAP(Triple.combined1.epi2.METPos1, reduction = "pca", dims = 1:16)
DimPlot(Triple.combined1.epi2.METPos1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = Triple.combined1.epi2.METPos1) <- "EpiCellTypes"
DimPlot(Triple.combined1.epi2.METPos1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))

##Convert Seurat to Monocle3 cell data set class
cds.Triple.combined1.epi2.METPos1 <- as.cell_data_set(Triple.combined1.epi2.METPos1)

#Construction single cell trajectories
cds.Triple.combined1.epi2.METPos1 <- cluster_cells(cds = cds.Triple.combined1.epi2.METPos1, reduction_method = "UMAP") 
cds.Triple.combined1.epi2.METPos1 <- learn_graph(cds.Triple.combined1.epi2.METPos1, use_partition = FALSE)

#plot single cell trajectory
tiff(file = "Triple.combined1.epi2.METPos1 EpiCellTypes trajectory UMAP.tiff", width = 7, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds = cds.Triple.combined1.epi2.METPos1,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = FALSE,
           label_cell_groups = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           cell_size = 0.5)+ scale_color_manual(values = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))
dev.off()

tiff(file = "Triple.combined1.epi2.METPos1 BELE EpiCellTypes trajectory w Branch points UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds = cds.Triple.combined1.epi2.METPos1,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = T, 
           label_leaves= T,
           label_roots = FALSE,
           group_label_size = 4,
           graph_label_size = 3)
dev.off()

#plot single cell trajectory by pseudotime
cds.Triple.combined1.epi2.METPos1 <- order_cells(cds.Triple.combined1.epi2.METPos1, reduction_method = "UMAP")

tiff(file = "Triple.combined1.epi2.METPos1 pseudotime trajectory UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds=cds.Triple.combined1.epi2.METPos1,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE, 
           label_roots = FALSE,
           cell_size = 0.5) 
dev.off()

tiff(file = "Triple.combined1.epi2.METPos1 pseudotime trajectory with branch points UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
plot_cells(cds=cds.Triple.combined1.epi2.METPos1,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=TRUE,
           label_branch_points=TRUE,
           label_roots = FALSE,
           graph_label_size=4)
dev.off()
#Subset Double
Idents(object = TriplevDouble.combined1.epi2.METPos1) <- "stim"
Double.combined1.epi2.METPos <- subset(TriplevDouble.combined1.epi2.METPos1, idents = c("Double"))

Double.combined1.epi2.METPos1 <- Double.combined1.epi2.METPos
DefaultAssay(Double.combined1.epi2.METPos1) <- "integrated"
Double.combined1.epi2.METPos1 <- ScaleData(Double.combined1.epi2.METPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Double.combined1.epi2.METPos1))
Double.combined1.epi2.METPos1 <- RunPCA(Double.combined1.epi2.METPos1, features = VariableFeatures(Double.combined1.epi2.METPos1))
ElbowPlot(Double.combined1.epi2.METPos1, ndims = 50)

#Res0.7
Double.combined1.epi2.METPos1 <- FindNeighbors(Double.combined1.epi2.METPos1, reduction = "pca", dims = 1:19)
Double.combined1.epi2.METPos1 <- FindClusters(Double.combined1.epi2.METPos1, resolution = 0.7)
Double.combined1.epi2.METPos1 <- RunUMAP(Double.combined1.epi2.METPos1, reduction = "pca", dims = 1:19)
DimPlot(Double.combined1.epi2.METPos1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = Double.combined1.epi2.METPos1) <- "EpiCellTypes"
DimPlot(Double.combined1.epi2.METPos1, reduction = "umap", pt.size = 0.3, label = TRUE, cols = c("chartreuse3", "chocolate4", "darkslategray3", "salmon", "brown3", "steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red"))    

####hMETtg+ Triple vs hMETtg+ Double####

#DEGs_hMETtg+ Triple vs hMETtg+ Double
Idents(object = TriplevDouble.combined1.epi2.METPos1) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.epi2.METPos1)
TriplevDouble.combined1.epi2.METPos1 <- ScaleData(TriplevDouble.combined1.epi2.METPos1, features = all.genes)
hMETPos_TriplevDouble.0.Markers <- FindMarkers(TriplevDouble.combined1.epi2.METPos1, ident.1 = "Triple", ident.2 = "Double", min.pct = 0, logfc.threshold = 0)
write.csv(hMETPos_TriplevDouble.0.Markers, "hMETPos_TriplevDouble.0.Markers.csv")

#p.adjust
DEG_hMETPos_TriplevDouble <- read.csv("hMETPos_TriplevDouble.0.Markers.csv") 
DEG_hMETPos_TriplevDouble_pvalue <- DEG_hMETPos_TriplevDouble$p_val
DEG_hMETPos_TriplevDouble_pvalue=as.numeric(DEG_hMETPos_TriplevDouble_pvalue)
DEG_hMETPos_TriplevDouble_BH = p.adjust(DEG_hMETPos_TriplevDouble_pvalue, "BH")
write.csv(DEG_hMETPos_TriplevDouble_BH, "DEG_hMETPos_TriplevDouble_BH.csv")



####LE1/LE4 vs others####

#DEGs_LE1vOthers
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)
LE1vsOthers.0.Markers <- FindMarkers(TriplevDouble.combined1.LE, ident.1 = "LE1", ident.2 = c("LE2", "LE3", "LE5", "LE6", "LE7", "LE8",
                                                                                              "LE9", "LE10", "LE11"), min.pct = 0, logfc.threshold = 0)
write.csv(LE1vsOthers.0.Markers, "LE1vsOthers.0.Markers.csv")

#p.adjust
DEG_LE1vsOthers <- read.csv("LE1vsOthers.0.Markers.csv") 
DEG_LE1vsOthers_pvalue <- DEG_LE1vsOthers$p_val
DEG_LE1vsOthers_pvalue=as.numeric(DEG_LE1vsOthers_pvalue)
DEG_LE1vsOthers_BH = p.adjust(DEG_LE1vsOthers_pvalue, "BH")
write.csv(DEG_LE1vsOthers_BH, "DEG_LE1vsOthers_BH.csv")

#DEGs_LE4vOthers
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
DimPlot(TriplevDouble.combined1.LE, reduction = "umap", pt.size = 0.3)
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)
LE4vsOthers.0.Markers <- FindMarkers(TriplevDouble.combined1.LE, ident.1 = "LE4", ident.2 = c("LE2", "LE3", "LE5", "LE6", "LE7", "LE8",
                                                                                              "LE9", "LE10", "LE11"), min.pct = 0, logfc.threshold = 0)
write.csv(LE4vsOthers.0.Markers, "LE4vsOthers.0.Markers.csv")

#p.adjust
DEG_LE4vsOthers <- read.csv("LE4vsOthers.0.Markers.csv") 
DEG_LE4vsOthers_pvalue <- DEG_LE4vsOthers$p_val
DEG_LE4vsOthers_pvalue=as.numeric(DEG_LE4vsOthers_pvalue)
DEG_LE4vsOthers_BH = p.adjust(DEG_LE4vsOthers_pvalue, "BH")
write.csv(DEG_LE4vsOthers_BH, "DEG_LE4vsOthers_BH.csv")

#Rename
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
TriplevDouble.combined1.LE <- RenameIdents(object = TriplevDouble.combined1.LE, 
                                             'LE1'="LE1",'LE4'="LE4",'LE2'="Others", 'LE3' = "Others", 'LE5'="Others", 
                                             'LE6'="Others",'LE7'="Others", 'LE8'="Others", 'LE9' = "Others", 'LE10'="Others",
                                             'LE11'="Others")  
TriplevDouble.combined1.LE[["LE14vOthers"]] <- Idents(object = TriplevDouble.combined1.LE)

#Vlnplots
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "LE14vOthers"
tiff(file = "TriplevDouble.combined1.LE Xpo1 Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Xpo1", pt.size = 0.01, cols = c("steelblue1", "blueviolet", "lightgrey")) + NoLegend()
dev.off()

#Cell Counts
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
TriplevDouble.combined1.LE$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.LE), TriplevDouble.combined1.LE$stim, sep = "_")
Idents(object = TriplevDouble.combined1.LE) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.LE))

#DEGs_LE4vOthers in Double
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "stim.EpiCellTypes"
DimPlot(TriplevDouble.combined1.LE, reduction = "umap", pt.size = 0.3)
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)
LE4_TriplevsOthers_Double.0.Markers <- FindMarkers(TriplevDouble.combined1.LE, ident.1 = "LE4_Triple", ident.2 = c("LE2_Double", "LE3_Double", "LE5_Double", 
                                                                                                     "LE6_Double", "LE7_Double", "LE8_Double",
                                                                                              "LE9_Double", "LE10_Double", "LE11_Double"), min.pct = 0, logfc.threshold = 0)
write.csv(LE4_TriplevsOthers_Double.0.Markers, "LE4_TriplevsOthers_Double.0.Markers.csv")

#p.adjust
DEG_LE4_TriplevsOthers_Double <- read.csv("LE4_TriplevsOthers_Double.0.Markers.csv") 
DEG_LE4_TriplevsOthers_Double_pvalue <- DEG_LE4_TriplevsOthers_Double$p_val
DEG_LE4_TriplevsOthers_Double_pvalue=as.numeric(DEG_LE4_TriplevsOthers_Double_pvalue)
DEG_LE4_TriplevsOthers_Double_BH = p.adjust(DEG_LE4_TriplevsOthers_Double_pvalue, "BH")
write.csv(DEG_LE4_TriplevsOthers_Double_BH, "DEG_LE4_TriplevsOthers_Double_BH.csv")

#DEGs_LE1vOthers in double
Idents(object = TriplevDouble.combined1.LE) <- "stim.EpiCellTypes"
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)
LE1_TriplevsOthers_Double.0.Markers <- FindMarkers(TriplevDouble.combined1.LE, ident.1 = "LE1_Triple", ident.2 = c("LE2_Double", "LE3_Double", "LE5_Double", 
                                                                                                                   "LE6_Double", "LE7_Double", "LE8_Double",
                                                                                                                   "LE9_Double", "LE10_Double", "LE11_Double"), min.pct = 0, logfc.threshold = 0)
write.csv(LE1_TriplevsOthers_Double.0.Markers, "LE1_TriplevsOthers_Double.0.Markers.csv")

#p.adjust
DEG_LE1_TriplevsOthers_Double <- read.csv("LE1_TriplevsOthers_Double.0.Markers.csv") 
DEG_LE1_TriplevsOthers_Double_pvalue <- DEG_LE1_TriplevsOthers_Double$p_val
DEG_LE1_TriplevsOthers_Double_pvalue=as.numeric(DEG_LE1_TriplevsOthers_Double_pvalue)
DEG_LE1_TriplevsOthers_Double_BH = p.adjust(DEG_LE1_TriplevsOthers_Double_pvalue, "BH")
write.csv(DEG_LE1_TriplevsOthers_Double_BH, "DEG_LE1_TriplevsOthers_Double_BH.csv")


####Other analysis####
##WNT
#Wnt downstream vlnplot
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
tiff(file = "TriplevDouble.combined1.LE Axin2 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Axin2", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Axin2 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Axin2", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.LE Tcf4 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Tcf4", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Tcf4 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Tcf4", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.LE Lef1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Lef1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Lef1 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Lef1", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.LE Stra6 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Stra6", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Stra6 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Stra6", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Ccnd1 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Ccnd1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Ccnd1 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Ccnd1", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

#StackedViolin for selected WNT downstream
features <- c("Axin2", "Tcf4", "Lef1", "Stra6")

b <- VlnPlot(TriplevDouble.combined1.LE, features, stack = TRUE, flip = TRUE) +
  theme(legend.position = "none") 
tiff(file = "TriplevDouble.combined1.LE Axin2, Tcf4, Lef1, Stra6 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
plot_grid(b)
dev.off()

features <- c("Axin2", "Tcf4", "Lef1", "Ccnd1", "Stra6")

b <- VlnPlot(TriplevDouble.combined1.LE, features, stack = TRUE, flip = TRUE) +
  theme(legend.position = "none") 
tiff(file = "TriplevDouble.combined1.LE Axin2, Tcf4, Lef1, Stra6 Ccnd1 Vln.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plot_grid(b)
dev.off()

#Stanford Wnt signature vlnplot
WNTlist <- list(c("Atoh1", "Axin2", "Birc5", "Ccnd1", "Cd44",
                  "Cdx4", "Cldn1", "Ctla4", "Dkk1", "Edn1", "Efnb1",
                  "Egfr", "En2", "Fgf18", "Fgf20", "Fgf9", "Fosl1",
                  "Fst", "Fzd7", "Gast", "Gja1", "Id2", "Jag1", "Jun",
                  "L1cam", "Lef1", "Lgr5", "Met", "Mmp7", "Myc", "Mycbp",
                  "Nrcam", "Plaur", "Ppard", "Pttg1", "Rhou", "Sfrp2",
                  "Sox17", "Stra6", "Tcf20", "Tcf4", "Tcf7", "Tcf7l2",
                  "Tiam1", "Tnfrsf19", "Twist1", "Vcan", "Vegfa"))
TriplevDouble.combined1.LE <- AddModuleScore(object = TriplevDouble.combined1.LE, features = WNTlist, name = "Stanford_WNT_Downstream") 
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
VlnPlot(TriplevDouble.combined1.LE, features = "Stanford_WNT_Downstream1", pt.size = 0.01)

tiff(file = "TriplevDouble.combined1.LE Standford Wnt signature Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Stanford_WNT_Downstream1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                             "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#GoBP Wnt signature vlnplot
GOBPWNTlist <- list(c("Lmbr1l",	 "Macf1",	 "Calcoco1",	 "Dlx3",	 "Tnks",	 "Cul3",	 "Ctr9",
                      "Lef1",	 "Ctnnd1",	 "Tcf7",	 "Tspan12",	 "Apcdd1",	 "Nkd1",	 "Wnt6",	 
                      "Znrf3",	 "Ccnd1",	 "Tnks2",	 "Wif1",	 "Lrrfip2",	 "Wnt4",	 "Mark1",
                      "Wls",	 "Rnf43",	 "Ryk",	 "Csnk1a1",	 "Wnt7b",	 "Csnk2a2",
                      "Wnt5a",	 "Csnk1d",	 "Axin2",	 "Csnk1e",	 "Dkk2",	 "Bcl9",	 "Bmp2",	 
                      "Apc",	 "Daam1",	 "Otulin",	 "Amotl2",	 "Hbp1",	 "Notum",	 "Tnik",	 "Cd44"))
TriplevDouble.combined1.LE <- AddModuleScore(object = TriplevDouble.combined1.LE, features = GOBPWNTlist, name = "GOBP_WNT_Downstream") 
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
VlnPlot(TriplevDouble.combined1.LE, features = "GOBP_WNT_Downstream1", pt.size = 0.01)

tiff(file = "TriplevDouble.combined1.LE GOBP_WNT_Downstream1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "GOBP_WNT_Downstream1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                    "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#KEGG signature vlnplot
KEGGWNTlist <- list(c("Rock2",	 "Ctnnd2",	 "Lef1",	 "Tcf7",	 "Cul1",	 "Prickle1",
                      "Apcdd1",	 "Nkd1",	 "Wnt6",	 "Mapk8",	 "Znrf3",	 "Ccnd1",
                      "Wif1",	 "Ep300",	 "Wnt4",	 "Rnf43",	 "Jun",	 "Crebbp",
                      "Smad3",	 "Ryk",	 "Csnk1a1",	 "Wnt7b",	 "Csnk2a2",	 "Wnt5a",
                      "Cacybp",	 "Axin2",	 "Csnk1e"	, "Dkk2",	 "Rbx1",	 "Vangl1",	 
                      "Apc",	 "Daam1",	 "Bambi",	 "Trp53",	 "Notum",	 "Lgr5",
                      "Ppard"))
TriplevDouble.combined1.LE <- AddModuleScore(object = TriplevDouble.combined1.LE, features = KEGGWNTlist, name = "KEGG_WNT_Downstream") 
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.LE KEGG_WNT_Downstream Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "KEGG_WNT_Downstream1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#StackedVolin for WNT downstream
features <- c("GOBP_WNT_Downstream1", "KEGG_WNT_Downstream1")

b <- VlnPlot(TriplevDouble.combined1.LE, features, stack = TRUE, flip = TRUE, pt.size = 0.01) +
  theme(legend.position = "none") 
tiff(file = "TriplevDouble.combined1.LE WNT signature Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_grid(b)
dev.off()

##RIOBOSOME
#Ribosome vlnplot
tiff(file = "TriplevDouble.combined1.LE Rpl12 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rpl12", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Rpl12 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rpl12", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.LE Rps16 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rps16", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Rps16 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rps16", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.LE Rpl3 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rpl3", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Rpl3 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rpl3", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.LE Rps5 Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rps5", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                 "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LE Rps5 split Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rps5", split.by = "stim", pt.size = 0.01) + NoLegend()
dev.off()

#StackedViolin for Ribosome
features <- c("Rpl12", "Rpl3", "Rps5", "Rps16", "Rps20")

b <- VlnPlot(TriplevDouble.combined1.LE, features, stack = TRUE, flip = TRUE, pt.size = 0.01) +
  theme(legend.position = "none") 
tiff(file = "TriplevDouble.combined1.LE Rpl12 Rpl3 Rps5 Rps16 Rps20 Vln.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plot_grid(b)
dev.off()

#Ribosome signature vlnplot
KEGG_Ribosomelist <- list(c("Rpl4", "Rpl5", "Rpl30", "Rpl3", "Rpl32", 
                            "Rpl34", "Rpl10a", "Rpl8", "Mrpl34", "Rpl9", 
                            "Rpl7", "Rps4x", "Rps15", "Rps14", "Mrpl3", "Rps17"
                            , "Rps16", "Rpl18a", "Rps19", "Rps18", "Rplp2", "Rpl37",
                            "Rps11", "Rps10", "Rps13", "Rps12", "Rps9", "Rps7",
                            "Rpl21", "Rps8", "Rps5", "Rpl23", "Rps6", "Rpl22", "Mrps18a"
                            ,  "Rpsa", "Mrps5", "Rpl24", "Rpl26", 
                            "Uba52", "Rpl28", "Rpl12", "Rpl11", "Rpl36a", "Mrpl14", "Rps15a"
                            , "Rpl14", "Rps3", "Rpl13", "Rps2", "Rpl15", "Rpl18", 
                            "Rps27a", "Rpl17", "Rpl19", "Rpl35a", "Rpl23a", "Rps26", 
                            "Rps3a1", "Rpl27a", "Rpl22l1", "Rps20", "Rsl24d1", "Rps24", "Rps23"))
TriplevDouble.combined1.LE <- AddModuleScore(object = TriplevDouble.combined1.LE, features = KEGG_Ribosomelist, ctrl = 100, name = "KEGG_Ribosomelist") 
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.LE KEGG_Ribosomelist Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "KEGG_Ribosomelist1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                                "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

GOBP_Ribosomelist <- list(c("Rpl4",	 "Rpl5",	 "Rpl3"	, "Rpl32",	 "Rpl34",	 "Rpl10a",
                            "Rpl8",	 "Mrpl34",	 "Rpl9",	 "Rpl7",	 "Eef1b2",
                            "Rbm3",	 "Rps15",	 "Rps14",	 "Mrpl3",	 "Rps17",	 "Rps16",
                            "Rpl18a",	 "Rps19",	 "Rps18",	 "Rars",	 "Rpl37",	 "Rps11",
                            "Rps13",	 "Rps12",	 "Eif5a",	 "Rps9",	 "Rps7",	 "Rpl21",	 
                            "Rps8",	 "Rps5",	 "Rpl23",	 "Rps6",	 "Rpl22",	 "Mrps18a",
                            "Rpl13a",	 "Rpsa",	 "Mrps5",	 "Mrpl43",	 "Eef1g",	 "Lars2",
                            "Eef1d",	 "Rpl24",	 "Rpl26",	 "Eif4e2",	 "Uba52",	 "Rpl28",
                            "Rpl12",	 "Rpl11",	 "Rpl36a",	 "Mrpl14",	 "Mrpl57",	 "Gspt1",
                            "Rps15a",	 "Rpl14",	 "Rps3",	 "Rpl13",	 "Rps2",	 "Rpl15",	 "Paip2",
                            "Rpl18",	 "Rps27a"	 ,"Rpl17"	, "Eif4b"	, "Rpl19"	, "Mcts1"	,
                            "Eif2b4",	 "Rpl35a"	, "Rpl23a",	 "Aimp1"	, "Rps26"	, "Eif3m",
                            "Eif6",	 "Eif3k",	 "Eif3l",	 "Rps3a1"	, "Rpl27a",	 "Eif3h",	 
                            "Rpl22l1",	 "Eif3e",	 "Rps20",	 "Eif3f",	 "Rsl24d1",	 "Rps24",	 
                            "Eif4g2",	 "Rps23"
))
TriplevDouble.combined1.LE <- AddModuleScore(object = TriplevDouble.combined1.LE, features = GOBP_Ribosomelist, ctrl = 100, name = "GOBP_Ribosomelist") 
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.LE GOBP_Ribosomelist Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "GOBP_Ribosomelist1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4",
                                                                                              "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

#StackedViolin Plot for Ribosome
features <- c("GOBP_Ribosomelist1", "KEGG_Ribosomelist1")

b <- VlnPlot(TriplevDouble.combined1.LE, features, stack = TRUE, flip = TRUE, pt.size = 0.01) +
  theme(legend.position = "none") 
tiff(file = "TriplevDouble.combined1.LE Ribosome signature Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_grid(b)
dev.off()

#G2M signature vlnplot
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
G2Mlist <- list(c("Sqle",	"Top1",	"Tubb4b",	"G3bp1",	"Gspt1",	"H2Afz",	"Hif1A",	"Nasp",	"Ncl",
                 "Nolc1",	"Sfpq",	"Snrpd1",	"Syncrip",	"Top1",	"Ube2S",	"Mki67",	"Cenpf",	"Cenpe",	"Cenpa"))
TriplevDouble.combined1.LE <- AddModuleScore(object = TriplevDouble.combined1.LE, features = G2Mlist, name = "G2Mscore")
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
VlnPlot(TriplevDouble.combined1.LE, features = "G2Mscore1", pt.size = 0.01)

tiff(file = "TriplevDouble.combined1.LE G2M signature Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "ARscore1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()





####Triple vs Double BE####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/Epi/BE")
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.BE <- subset(TriplevDouble.combined1.epi2, idents = c("BE1", "BE2", "BE3", "BE4", "BE5"))

#Heatmap 
DefaultAssay(TriplevDouble.combined1.BE) <- "RNA"
Idents(object = TriplevDouble.combined1.BE) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.BE)
TriplevDouble.combined1.BE <- ScaleData(TriplevDouble.combined1.BE, features = all.genes)
TriplevDouble.combined1.BE.allmarkers <- FindAllMarkers(TriplevDouble.combined1.BE, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(TriplevDouble.combined1.BE.allmarkers, file = "TriplevDouble.combined1.BE.allmarkers.csv")

TriplevDouble.combined1.BE.allmarkers.Top50 <- TriplevDouble.combined1.BE.allmarkers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TriplevDouble.combined1.BE.allmarkers.Top50 Heatmap.tiff", width =12, height = 12, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.BE, features = c(TriplevDouble.combined1.BE.allmarkers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
#DEGs_SolidvNormal
DefaultAssay(TriplevDouble.combined1.BE) <- "RNA"
Idents(object = TriplevDouble.combined1.BE) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.BE)
TriplevDouble.combined1.BE <- ScaleData(TriplevDouble.combined1.BE, features = all.genes)

BE_TriplevDouble.0.Markers <- FindMarkers(TriplevDouble.combined1.BE, ident.1 = "Triple", ident.2 = "Double", min.pct = 0, logfc.threshold = 0)
write.csv(BE_TriplevDouble.0.Markers, "BE_TriplevDouble.0.Markers.csv")

#p.adjust
DEG_BE_TriplevDouble <- read.csv("BE_TriplevDouble.0.Markers.csv") 
DEG_BE_TriplevDouble_pvalue <- DEG_BE_TriplevDouble$p_val
DEG_BE_TriplevDouble_pvalue=as.numeric(DEG_BE_TriplevDouble_pvalue)
DEG_BE_TriplevDouble_BH = p.adjust(DEG_BE_TriplevDouble_pvalue, "BH")
write.csv(DEG_BE_TriplevDouble_BH, "DEG_BE_TriplevDouble_BH.csv")

####Triple LE vs Double LE####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/Epi/LE")

Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.LE <- subset(TriplevDouble.combined1.epi2, idents = c("LE1", "LE2", "LE3", "LE4", "LE5",
                                                                              "LE6", "LE7", "LE8", "LE9", "LE10", "LE11"))

#Heatmap 
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)
TriplevDouble.combined1.LE.allmarkers <- FindAllMarkers(TriplevDouble.combined1.LE, min.pct = 0.25, logfc.threshold = 0.25, only.pos=TRUE)
write.csv(TriplevDouble.combined1.LE.allmarkers, file = "TriplevDouble.combined1.LE.allmarkers.csv")

TriplevDouble.combined1.LE.allmarkers.Top50 <- TriplevDouble.combined1.LE.allmarkers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TriplevDouble.combined1.LE.allmarkers.Top50 Heatmap.tiff", width =12, height = 12, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.LE, features = c(TriplevDouble.combined1.LE.allmarkers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
#DEGs_LE_TriplevDouble
DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.LE)
TriplevDouble.combined1.LE <- ScaleData(TriplevDouble.combined1.LE, features = all.genes)

LE_TriplevDouble.0.Markers <- FindMarkers(TriplevDouble.combined1.LE, ident.1 = "Triple", ident.2 = "Double", min.pct = 0, logfc.threshold = 0)
write.csv(LE_TriplevDouble.0.Markers, "LE_TriplevDouble.0.Markers.csv")

#p.adjust
DEG_LE_TriplevDouble <- read.csv("LE_TriplevDouble.0.Markers.csv") 
DEG_LE_TriplevDouble_pvalue <- DEG_LE_TriplevDouble$p_val
DEG_LE_TriplevDouble_pvalue=as.numeric(DEG_LE_TriplevDouble_pvalue)
DEG_LE_TriplevDouble_BH = p.adjust(DEG_LE_TriplevDouble_pvalue, "BH")
write.csv(DEG_LE_TriplevDouble_BH, "DEG_LE_TriplevDouble_BH.csv")



