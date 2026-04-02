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

TriplevDouble.combined1 <- FindNeighbors(TriplevDouble.combined1, reduction = "pca", dims = 1:26)
TriplevDouble.combined1 <- FindClusters(TriplevDouble.combined1, resolution = 1.5)
TriplevDouble.combined1 <- RunUMAP(TriplevDouble.combined1, reduction = "pca", dims = 1:26)

Idents(object = TriplevDouble.combined1) <- "seurat_clusters"
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, split.by='stim',label = TRUE)

Idents(object = TriplevDouble.combined1) <- "seurat_clusters"
tiff(file = "TriplevDouble.combined1 UMAP dims26.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
dev.off()

Idents(object = TriplevDouble.combined1) <- "Phase"
tiff(file = "TriplevDouble.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell type identification
DefaultAssay(TriplevDouble.combined1)<-"RNA"
tiff(file = "TriplevDouble.combined1 celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                      "Fbln1", "Myh11", "Pecam1",
                                                                      "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename CellTypes
Idents(object = TriplevDouble.combined1) <- "seurat_clusters"
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, label=TRUE)
TriplevDouble.combined1 <- RenameIdents(object = TriplevDouble.combined1, 
                                        '1'="BE",'3'="BE",'20'="BE",
                                        '6'="LE", '7'="LE", "0"="LE", '19'="LE",
                                        '2'="LE",'12'="LE",'5'="LE", '8'="LE",'11'="LE",
                                        '22'="LE",'21'="LE",'18'="LE",'13'="LE",'4'="LE",'25'="SV",
                                        '26'="SV", '16'="SV",
                                        '10'="FB",'14'="FB", '23'="SM",'27'="Pericyte", '17'="VE",
                                        '9'="Immune",'24'="Immune", '15'="Immune")  
TriplevDouble.combined1[["CellTypes"]] <- Idents(object = TriplevDouble.combined1)

#UMAP
Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 CellTypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()
tiff(file = "TriplevDouble.combined1 CellTypes split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1, reduction = "umap", pt.size = 0.3, split.by='stim', cols = c("chartreuse3", "salmon", "darkslategray3", "plum4", "brown3", "blueviolet", "blue", "bisque3", "steelblue1"))
dev.off()

#FeaturePlots
DefaultAssay(TriplevDouble.combined1) <- "RNA"
tiff(file = "TriplevDouble.combined1 hMETtg expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 2.5)
dev.off()
tiff(file = "TriplevDouble.combined1 hMETtg split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Ar split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Tmprss2 split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Tmprss2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Pbsn split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Krt5 split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Krt8 split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1 Vim split expression plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1, reduction = "umap", features = c("Vim"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()

#Cell counts
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1$stim.CellTypes <- paste(Idents(TriplevDouble.combined1), TriplevDouble.combined1$stim, sep = "_")
Idents(object = TriplevDouble.combined1) <- "stim.CellTypes"
table(Idents(TriplevDouble.combined1))

#DEGs
DefaultAssay(TriplevDouble.combined1) <- "RNA"
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1 <- ScaleData(TriplevDouble.combined1, features = rownames(TriplevDouble.combined1))
TriplevDouble.combined1.allMarkers <- FindAllMarkers(TriplevDouble.combined1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.allMarkers, "TriplevDouble.combined1.CellTypes.Markers.csv")

#Dotplot
Idents(object = TriplevDouble.combined1) <- "CellTypes"
tiff(file = "TriplevDouble.combined1 CellTypes markers DotPlot.tiff", width =12 , height = 3.1, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1, features = c("Krt15", "Krt14", "Krt5", "Aqp3", "Col17a1", 
                                           "Prr9", "Pigr","Slc12a2", "Arl14", "Bex1",
                                          "Svs4",  "Pate4", "A630095E13Rik", "D730048I06Rik", "Sprr2f", 
                                           "Apod", "Bgn", "Igfbp6", "Serping1", "Col1a2",
                                           "Ndufa4l2", "Adamts4", "Myh11", "Pdgfrb", "Des",
                                           "Plp1", "Kcna1", "Cdh19", "S100b", "Fxyd1",
                                           "Aqp1", "Plvap", "Cdh5", "Pecam1", "Cd93",
                                           "H2-Eb1", "C1qa", "Rgs1", "Tyrobp", "Fcer1g"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()


####Subcluster epi DoublevTriple####
Idents(object = TriplevDouble.combined1) <- "CellTypes"
TriplevDouble.combined1.epi <- subset(TriplevDouble.combined1, idents = c("BE","LE"))
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell counts
Idents(object = TriplevDouble.combined1.epi) <- "CellTypes"
TriplevDouble.combined1.epi$stim.CellTypes <- paste(Idents(TriplevDouble.combined1.epi), TriplevDouble.combined1.epi$stim, sep = "_")
Idents(object = TriplevDouble.combined1.epi) <- "stim.CellTypes"
table(Idents(TriplevDouble.combined1.epi))

Idents(object = TriplevDouble.combined1.epi) <- "seurat_clusters"

#Run the standard workflow for visualization and clustering
DefaultAssay(TriplevDouble.combined1.epi) <- "integrated"
TriplevDouble.combined1.epi <- ScaleData(TriplevDouble.combined1.epi, verbose = FALSE)
TriplevDouble.combined1.epi <- RunPCA(TriplevDouble.combined1.epi, npcs = 50, verbose = FALSE)
ElbowPlot(TriplevDouble.combined1.epi, ndims = 50)

#Umap and Clustering
TriplevDouble.combined1.epi <- FindNeighbors(TriplevDouble.combined1.epi, reduction = "pca", dims = 1:20)
TriplevDouble.combined1.epi <- FindClusters(TriplevDouble.combined1.epi, resolution = 0.7)
TriplevDouble.combined1.epi <- RunUMAP(TriplevDouble.combined1.epi, reduction = "pca", dims = 1:20)
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle assignment epi
DefaultAssay(TriplevDouble.combined1.epi) <- "RNA"
all.genes <- rownames(TriplevDouble.combined1.epi)
TriplevDouble.combined1.epi <- ScaleData(TriplevDouble.combined1.epi, features = all.genes)
TriplevDouble.combined1.epi <- CellCycleScoring(TriplevDouble.combined1.epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TriplevDouble.combined1.epi) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Cell Cycle Regression epi
TriplevDouble.combined1.epi1 <- TriplevDouble.combined1.epi
DefaultAssay(TriplevDouble.combined1.epi1) <- "integrated"
TriplevDouble.combined1.epi1 <- ScaleData(TriplevDouble.combined1.epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi1))
TriplevDouble.combined1.epi1 <- RunPCA(TriplevDouble.combined1.epi1, features = VariableFeatures(TriplevDouble.combined1.epi1))
ElbowPlot(TriplevDouble.combined1.epi1, ndims = 50)

TriplevDouble.combined1.epi1 <- FindNeighbors(TriplevDouble.combined1.epi1, reduction = "pca", dims = 1:23)
TriplevDouble.combined1.epi1 <- FindClusters(TriplevDouble.combined1.epi1, resolution = 0.8)
TriplevDouble.combined1.epi1 <- RunUMAP(TriplevDouble.combined1.epi1, reduction = "pca", dims = 1:23)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

tiff(file = "TriplevDouble.combined1.epi1 seurat dim23 UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
dev.off()

#DEGs
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1) <- "seurat_clusters"
TriplevDouble.combined1.epi1 <- ScaleData(TriplevDouble.combined1.epi1, features = rownames(TriplevDouble.combined1.epi1))
TriplevDouble.combined1.epi1.allMarkers <- FindAllMarkers(TriplevDouble.combined1.epi1, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi1.allMarkers, "TriplevDouble.combined1.epi1.seurat.allMarkers.csv")

Idents(object = TriplevDouble.combined1.epi1) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi1 after Cell Cyle Regression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Featureplots
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi1 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Krt5", "Krt14", "Krt19", "Krt8", "Ppp1r1b", "Vim", "hMETtg", "Trp63"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()


#Rename Clusters
Idents(object = TriplevDouble.combined1.epi1) <- "seurat_clusters"
TriplevDouble.combined1.epi1 <- RenameIdents(object = TriplevDouble.combined1.epi1, 
                                             '5'="BE1", '20'="BE1",'6'="BE2",'13'="BE3", '14' = "BE3", '18'="BE3", 
                                             '7'="LE1",'12'="LE1", '4'="LE2", '11' = "LE2", '3'="LE3",
                                             '0' = "LE4", '1'="LE4", '19'="LE4",
                                             '9'="LE5",
                                             '15'="LE6", '2' = "LE6",
                                             '10'="LE7", '8'="LE8", '16'="UrLE", 
                                             '17'="OE")  
TriplevDouble.combined1.epi1[["EpiCellTypes"]] <- Idents(object = TriplevDouble.combined1.epi1)

#Umap Epicelltype
tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes UMAP.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label.size = 6, 
        cols = c("skyblue1", "chartreuse3", "brown3", "deeppink1", "blue",  "red",  "khaki4","aquamarine3",
                 "bisque3", "salmon",  "plum4", "darkorange1", "grey"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes split UMAP.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label.size = 6, split.by = "stim",
        cols = c("skyblue1", "chartreuse3", "brown3", "deeppink1", "blue", "red", "khaki4", "aquamarine3", 
                 "bisque3", "salmon",  "plum4", "darkorange1", "grey"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes split label UMAP.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label.size = 6, split.by = "stim",
        cols = c("skyblue1", "chartreuse3", "brown3", "deeppink1", "blue", "red", "khaki4", "aquamarine3", 
                 "bisque3", "salmon",  "plum4", "darkorange1", "grey"), label = TRUE)
dev.off()

#Featureplot
DefaultAssay(TriplevDouble.combined1.epi1)<-"RNA"
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi1 hMETtg split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Tcf4 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Tcf4"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Axin2 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Axin2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Tmprss2 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Tmprss2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Nkx3-1 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Nkx3-1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", min.cutoff = "q5", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Pbsn split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Ar split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()

tiff(file = "TriplevDouble.combined1.epi1 Krt5 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Krt8 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", min.cutoff = "q5", keep.scale = "all")
dev.off()

DefaultAssay(TriplevDouble.combined1.epi1)<-"RNA"
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi1 hMETtg plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), pt.size = 0.5,  min.cutoff = "q5")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Tcf4 plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Tcf4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Axin2 plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Axin2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Tmprss2 plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Tmprss2"),cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Nkx3-1 plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Nkx3-1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", min.cutoff = "q5", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Pbsn plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()
tiff(file = "TriplevDouble.combined1.epi1 Ar plots.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()

#Cell counts
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi1), TriplevDouble.combined1.epi1$stim, sep = "_")
Idents(object = TriplevDouble.combined1.epi1) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi1))

Idents(object = TriplevDouble.combined1.epi1) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

Idents(object = TriplevDouble.combined1.epi1) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi1 Cell Cyle for count UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey", "gold"))
dev.off()

#Cell counts
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1$Phase.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi1), TriplevDouble.combined1.epi1$Phase, sep = "_")
Idents(object = TriplevDouble.combined1.epi1) <- "Phase.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi1))

#DEGs
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1 <- ScaleData(TriplevDouble.combined1.epi1, features = rownames(TriplevDouble.combined1.epi1))
TriplevDouble.combined1.epi1.allMarkers <- FindAllMarkers(TriplevDouble.combined1.epi1, min.pct = 0.1, logfc.threshold = 0.1, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi1.allMarkers, "TriplevDouble.combined1.epi1.EpiCellTypes.allMarkers.csv")

#Dotplot
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi1 EpiCellTypes markers DotPlot.tiff", width =16 , height = 4.3, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1.epi1, features = c("hMETtg", "hHGFtg",
                                                   "Apoe", "Ctsl", "Tpm1", "Ncl", "Fhl1", 
                                                   "Aqp3", "Col17a1", "Anxa8", "Lamb3", "Inhba", 
                                                   "Sncg", "Dapl1", "Ifi202b", "Gpnmb", "Fscn1",
                                                   "Lars2", "Gm42418", "Svs5", "AY036118", "Gm26917", 
                                                   "Defa22", "Defa21", "Defa5", "Cryba4", "Slc40a1",
                                                   "Coch", "Tgfb2", "Dkk2", "Ier3", "Apoc1",
                                                   "Clu", "Mmp7", "Areg", "Ly6a", "Basp1", 
                                                   "9530053A07Rik", "Tgm4", "Gm5615","Chn2","Spink8",
                                                   "Pigr", "Cldn3", "Lrg1", "Tnfrsf21", "Dcxr",
                                                   "Msmb", "Mme", "Apof", "Pcp4", "Reg3g", 
                                                   "Spink1", "Sbpl", "Timp4", "Crabp1", "Col6a3",
                                                   "Gsdmc2", "Gsdmc3", "Cxcl15", "Scnn1a", "Krt4",
                                                   "Serping1", "Igfbp6", "Fbln1", "Serpinf1", "Lum"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

##Heatmap
#scale.data
Epi_df <- as.data.frame(t(TriplevDouble.combined1.epi1$RNA@scale.data))
Epi_df_selected <- cbind( Epi_df$Cenpe, Epi_df$Cenpf, Epi_df$Top2a, Epi_df$Ube2s, Epi_df$Slc7a5,Epi_df$Pcna,Epi_df$Rangap1,
                          Epi_df$Tmprss2, Epi_df$Pbsn, Epi_df$`Nkx3-1`, Epi_df$Fkbp5, Epi_df$Azgp1,  
                          Epi_df$Rap1b, Epi_df$Grb2, Epi_df$Kras, Epi_df$Spint1, Epi_df$Uba52, Epi_df$Tns3, Epi_df$Bcar1, Epi_df$Hpn, 
                          Epi_df$Axin2, Epi_df$Cd44, Epi_df$Mmp7, Epi_df$Dkk2, Epi_df$Tcf4, Epi_df$Myc,
                          Epi_df$Rpl12, Epi_df$Rpl13a,Epi_df$Rpl3, Epi_df$Rps16, Epi_df$Rps20, Epi_df$Rps27a,Epi_df$Rps5)
Epi_df_selected <- as.data.frame(Epi_df_selected)
colnames(Epi_df_selected) <- c( "Cenpe",	"Cenpf",	"Top2a", "Ube2s", "Slc7a5","Pcna","Rangap1",
                                "Tmprss2", "Pbsn",	'Nkx3-1', "Fkbp5", "Azgp1",	
                                "Rap1b",	"Grb2",	"Kras",	"Spint1",	"Uba52",	"Tns3",	"Bcar1",	"Hpn",
                                "Axin2", "Cd44",	"Mmp7",	"Dkk2", "Tcf4", "Myc",
                                "Rpl12", "Rpl13a", "Rpl3", "Rps16", "Rps20", "Rps27a", "Rps5")
rownames(Epi_df_selected) <- row.names(Epi_df)
write.csv(Epi_df_selected, file = "Epi_df_selected_scaledata.csv")


#meta.data
write.csv(TriplevDouble.combined1.epi1@meta.data, file = "TriplevDouble.combined1.epi1_metadata.csv")

#Env
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggplotify)

df <- read.csv("Heatmap_Triple_and_Double_Epi_EpiCellTypes-1.csv", header = TRUE, sep = ",")
colnames(df)[1] <- "Clusters"
df <- as.data.frame(df)

#
df <- column_to_rownames(df, var = "Clusters")
df <- as.matrix(df)
d <- pheatmap::pheatmap(df, cluster_cols = F, cluster_rows = F)

##Heatmap2
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
# Change to log10
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)

tiff(file = "TriplevDouble.combined1.Epi EpiCellTypes Heatmap G2M AR MET WNT AR Ribosome.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)
dev.off()

####hMETtgPos vs hMETtgNeg####
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"

TriplevDouble.combined1.epi1METPos <- subset(x=TriplevDouble.combined1.epi1,  subset = `hMETtg` > 0)
TriplevDouble.combined1.epi1METNeg <- subset(x=TriplevDouble.combined1.epi1,  subset = `hMETtg` == 0)
Idents(object = TriplevDouble.combined1.epi1METPos) <- "METPos"
Idents(object = TriplevDouble.combined1.epi1METNeg) <- "METNeg"
TriplevDouble.combined1.epi1METPos[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi1METPos)
TriplevDouble.combined1.epi1METNeg[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi1METNeg)
TriplevDouble.combined1.epi1MET <- merge(x = TriplevDouble.combined1.epi1METPos, y = TriplevDouble.combined1.epi1METNeg)
Idents(object = TriplevDouble.combined1.epi1MET) <- "METExp"
TriplevDouble.combined1.epi1$METExp <- Idents(object = TriplevDouble.combined1.epi1MET)
Idents(object = TriplevDouble.combined1.epi1) <- "METExp"
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

tiff(file = "TriplevDouble.combined1.epi1 hMETPos highlighted stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"), split.by = "stim")
dev.off()

#DEGs_hMETtg+Triple vs hMETtg+Double
Idents(object = TriplevDouble.combined1.epi1METPos) <- "EpiCellTypes"
DimPlot(TriplevDouble.combined1.epi1METPos, reduction = "umap", pt.size = 0.3, label = TRUE)

TriplevDouble.combined1.BELEMETPos <- subset(TriplevDouble.combined1.epi1METPos, idents = c("BE1","BE2", "BE3",
                                                                                     "LE1", "LE2", "LE3", "LE4",
                                                                                     "LE5", "LE6", "LE7", "LE8"))
DimPlot(TriplevDouble.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TriplevDouble.combined1.BELEMETPos) <- "EpiCellTypes"
TriplevDouble.combined1.BELEMETPos$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.BELEMETPos), TriplevDouble.combined1.BELEMETPos$stim, sep = "_")
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.BELEMETPos))

all.genes <- rownames(TriplevDouble.combined1.BELEMETPos)
TriplevDouble.combined1.BELEMETPos <- ScaleData(TriplevDouble.combined1.BELEMETPos, features = all.genes)
hMETtgPos_LE_TriplevDouble.0.Markers <- FindMarkers(TriplevDouble.combined1.BELEMETPos, ident.1 = c("LE1_Triple", "LE2_Triple", "LE3_Triple",
                                                                                                    "LE4_Triple","LE5_Triple","LE6_Triple",
                                                                                                    "LE7_Triple","LE8_Triple"), 
                                                    ident.2 = c("LE1_Double", "LE2_Double","LE3_Double",
                                                                "LE4_Double","LE5_Double","LE6_Double",
                                                                "LE7_Double","LE8_Double"), min.pct = 0, logfc.threshold = 0)
write.csv(hMETtgPos_LE_TriplevDouble.0.Markers, "hMETtgPos_LE_TriplevDouble.0.Markers.csv")

#p.adjust
DEG_hMETtgPos_LE_TriplevDouble <- read.csv("hMETtgPos_LE_TriplevDouble.0.Markers.csv") 
DEG_hMETtgPos_LE_TriplevDouble_pvalue <- DEG_hMETtgPos_LE_TriplevDouble$p_val
DEG_hMETtgPos_LE_TriplevDouble_pvalue=as.numeric(DEG_hMETtgPos_LE_TriplevDouble_pvalue)
DEG_hMETtgPos_LE_TriplevDouble_BH = p.adjust(DEG_hMETtgPos_LE_TriplevDouble_pvalue, "BH")
write.csv(DEG_hMETtgPos_LE_TriplevDouble_BH, "DEG_hMETtgPos_LE_TriplevDouble_BH.csv")

#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named Symbol
genelist = DEG_hMETtgPos_LE_TriplevDouble_preranked$symbol
ortho_DEG_hMETtgPos_LE_TriplevDouble_preranked<- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

DEG_hMETtgPos_LE_TriplevDouble_preranked <- inner_join(DEG_hMETtgPos_LE_TriplevDouble_preranked,ortho_DEG_hMETtgPos_LE_TriplevDouble_preranked,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(DEG_hMETtgPos_LE_TriplevDouble_preranked,"converted_gsea_DEG_hMETtgPos_LE_TriplevDouble_preranked.txt")

#Subset hMETtgPos LE
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "EpiCellTypes"
TriplevDouble.combined1.LEMETPos <- subset(TriplevDouble.combined1.BELEMETPos, idents = c("LE1", "LE2", "LE3",
                                                                               "LE4", "LE5", "LE6", "LE7", "LE8"))

#Rename
Idents(object = TriplevDouble.combined1.LEMETPos) <- "stim"
TriplevDouble.combined1.LEMETPos <- RenameIdents(object = TriplevDouble.combined1.LEMETPos, 
                                               'Double'="Double",'Triple'="Triple")  
TriplevDouble.combined1.LEMETPos[["stim"]] <- Idents(object = TriplevDouble.combined1.LEMETPos)

#Vlnplots
DefaultAssay(TriplevDouble.combined1.LEMETPos) <- "RNA"
Idents(object = TriplevDouble.combined1.LEMETPos) <- "stim"
tiff(file = "TriplevDouble.combined1.LEMETPos Xpo1 Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LEMETPos, features = "Xpo1", pt.size = 0.01, cols = c("grey",  "red"), y.max = 1.2) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LEMETPos Ran Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LEMETPos, features = "Ran", pt.size = 0.01, cols = c("grey",  "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LEMETPos Rpl12 Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LEMETPos, features = "Rpl12", pt.size = 0.01, cols = c("grey",  "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LEMETPos Rps16 Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LEMETPos, features = "Rps16", pt.size = 0.01, cols = c("grey",  "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LEMETPos Myc Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LEMETPos, features = "Myc", pt.size = 0.01, cols = c("grey",  "red")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.LEMETPos Eif4e2 Vln.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LEMETPos, features = "Eif4e2", pt.size = 0.01, cols = c("grey",  "red")) + NoLegend()
dev.off()

####LE2 vs LE3 vs LE45678####
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1.LE <- subset(TriplevDouble.combined1.epi1, idents = c("LE1", "LE2", "LE3", "LE4","LE5", "LE6", "LE7", "LE8"))
DimPlot(TriplevDouble.combined1.epi1.LE, reduction = "umap", pt.size = 0.3, label = TRUE)

#Rename
Idents(object = TriplevDouble.combined1.epi1.LE) <- "EpiCellTypes"
TriplevDouble.combined1.epi1.LE <- RenameIdents(object = TriplevDouble.combined1.epi1.LE, 'LE2'="LE2",'LE3' = "LE3",
                                        'LE4'="LE145678", 'LE5'="LE145678", 'LE6' = "LE145678", 'LE7' = "LE145678", 'LE8' = "LE145678", 'LE1'="LE145678")  
TriplevDouble.combined1.epi1.LE[["OthersvLE23"]] <- Idents(object = TriplevDouble.combined1.epi1.LE)

#Vlnplots
DefaultAssay(TriplevDouble.combined1.epi1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1.LE) <- "OthersvLE23"

tiff(file = "TriplevDouble.combined1.epi1.LE Xpo1 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Xpo1", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey"), y.max = 1.5) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Ran Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Ran", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey"), y.max = 3.5) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Rpl12 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Rpl12", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Rps16 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Rps16", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Myc Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Myc", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey"), y.max = 3) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Eif4a1 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Eif4a1", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Eif4e2 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Eif4e2", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()

tiff(file = "TriplevDouble.combined1.epi1.LE Xpo1 Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Xpo1", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey"), y.max = 1.5) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Ran Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Ran", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey"), y.max = 3.5) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Rpl12 Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Rpl12", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Rps16 Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Rps16", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Myc Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Myc", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey"), y.max = 3) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Eif4a1 Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Eif4a1", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE Eif4e2 Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE, features = "Eif4e2", pt.size = 0, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
#DEGS
DefaultAssay(TriplevDouble.combined1.epi1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1.LE) <- "OthersvLE23"
all.genes <- rownames(TriplevDouble.combined1.epi1.LE)
TriplevDouble.combined1.epi1.LE <- ScaleData(TriplevDouble.combined1.epi1.LE, features = all.genes)
LE2v145678.0.Markers <- FindMarkers(TriplevDouble.combined1.epi1.LE, ident.1 = c("LE2"), 
                                                    ident.2 = c("LE145678"), min.pct = 0, logfc.threshold = 0)
write.csv(LE2v145678.0.Markers, "LE2v145678.0.Markers.csv")

#p.adjust
DEG_LE2v145678 <- read.csv("LE2v145678.0.Markers.csv") 
DEG_LE2v145678_pvalue <- DEG_LE2v45678$p_val
DEG_LE2v145678_pvalue=as.numeric(DEG_LE2v145678_pvalue)
DEG_LE2v145678_BH = p.adjust(DEG_LE2v145678_pvalue, "BH")
write.csv(DEG_LE2v145678_BH, "DEG_LE2v145678_BH.csv")

#DEGS
DefaultAssay(TriplevDouble.combined1.epi1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1.LE) <- "OthersvLE23"
all.genes <- rownames(TriplevDouble.combined1.epi1.LE)
TriplevDouble.combined1.epi1.LE <- ScaleData(TriplevDouble.combined1.epi1.LE, features = all.genes)
LE3v145678.0.Markers <- FindMarkers(TriplevDouble.combined1.epi1.LE, ident.1 = c("LE3"), 
                                   ident.2 = c("LE145678"), min.pct = 0, logfc.threshold = 0)
write.csv(LE3v145678.0.Markers, "LE3v145678.0.Markers.csv")

#p.adjust
DEG_LE3v145678 <- read.csv("LE3v145678.0.Markers.csv") 
DEG_LE3v145678_pvalue <- DEG_LE3v145678$p_val
DEG_LE3v145678_pvalue=as.numeric(DEG_LE3v145678_pvalue)
DEG_LE3v145678_BH = p.adjust(DEG_LE3v145678_pvalue, "BH")
write.csv(DEG_LE3v145678_BH, "DEG_LE3v145678_BH.csv")


####LE2 vs LE3 vs Others####
Idents(object = TriplevDouble.combined1.epi1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1.LE1 <- subset(TriplevDouble.combined1.epi1, idents = c("LE1", "LE2", "LE3", "LE5", "LE6", "LE7", "LE8"))
DimPlot(TriplevDouble.combined1.epi1.LE1, reduction = "umap", pt.size = 0.3, label = TRUE)

#Rename
Idents(object = TriplevDouble.combined1.epi1.LE1) <- "EpiCellTypes"
TriplevDouble.combined1.epi1.LE1 <- RenameIdents(object = TriplevDouble.combined1.epi1.LE1, 'LE2'="LE2",'LE3' = "LE3",
                                                'LE5' = "LE5678", 'LE6' = "LE5678", 'LE7' = "LE5678", 'LE8' = "LE5678", 'LE1'="LE5678")  
TriplevDouble.combined1.epi1.LE1[["OthersvLE123"]] <- Idents(object = TriplevDouble.combined1.epi1.LE1)

#Vlnplots
DefaultAssay(TriplevDouble.combined1.epi1.LE1) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1.LE1) <- "OthersvLE123"

tiff(file = "TriplevDouble.combined1.epi1.LE1 Xpo1 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Xpo1", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey"), y.max = 1.5) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE1 Ran Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Ran", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey"), y.max = 3.5) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE1 Rpl12 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Rpl12", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE1 Rps16 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Rps16", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE1 Myc Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Myc", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey"), y.max = 3) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE1 Eif4a1 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Eif4a1", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()
tiff(file = "TriplevDouble.combined1.epi1.LE1 Eif4e2 Vln-1.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.epi1.LE1, features = "Eif4e2", pt.size = 0.01, cols = c("deeppink1", "blue", "red", "grey")) + NoLegend()
dev.off()

#DEGS
DefaultAssay(TriplevDouble.combined1.epi1.LE1) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1.LE1) <- "OthersvLE123"
all.genes <- rownames(TriplevDouble.combined1.epi1.LE1)
TriplevDouble.combined1.epi1.LE1 <- ScaleData(TriplevDouble.combined1.epi1.LE1, features = all.genes)
LE2v5678.0.Markers <- FindMarkers(TriplevDouble.combined1.epi1.LE1, ident.1 = c("LE2"), 
                                  ident.2 = c("LE5678"), min.pct = 0, logfc.threshold = 0)
write.csv(LE2v5678.0.Markers, "LE2v5678.0.Markers.csv")

#p.adjust
DEG_LE2v5678 <- read.csv("LE2v5678.0.Markers.csv") 
DEG_LE2v5678_pvalue <- DEG_LE2v5678$p_val
DEG_LE2v5678_pvalue=as.numeric(DEG_LE2v5678_pvalue)
DEG_LE2v5678_BH = p.adjust(DEG_LE2v5678_pvalue, "BH")
write.csv(DEG_LE2v5678_BH, "DEG_LE2v5678_BH.csv")

#DEGS
DefaultAssay(TriplevDouble.combined1.epi1.LE1) <- "RNA"
Idents(object = TriplevDouble.combined1.epi1.LE1) <- "OthersvLE123"
all.genes <- rownames(TriplevDouble.combined1.epi1.LE1)
TriplevDouble.combined1.epi1.LE1 <- ScaleData(TriplevDouble.combined1.epi1.LE1, features = all.genes)
LE3v5678.0.Markers <- FindMarkers(TriplevDouble.combined1.epi1.LE1, ident.1 = c("LE3"), 
                                  ident.2 = c("LE5678"), min.pct = 0, logfc.threshold = 0)
write.csv(LE3v5678.0.Markers, "LE3v5678.0.Markers.csv")

#p.adjust
DEG_LE3v5678 <- read.csv("LE3v5678.0.Markers.csv") 
DEG_LE3v5678_pvalue <- DEG_LE3v5678$p_val
DEG_LE3v5678_pvalue=as.numeric(DEG_LE3v5678_pvalue)
DEG_LE3v5678_BH = p.adjust(DEG_LE3v5678_pvalue, "BH")
write.csv(DEG_LE3v5678_BH, "DEG_LE3v5678_BH.csv")

