#hHGFtg-hMETtg & hHGFtg-hMETtg-Bcat

#Add necessary tools to library
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(Seurat)
library(devtools)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/LabelledMulti")
list.files()
scRNA <- readRDS("scRNA.rds")

####Merging Dataset####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble")

#subset primary
Idents(object = HGFMETBcat_Multi1) <- "orig.ident"
DimPlot(HGFMETBcat_Multi1, reduction = "umap")
HGFMETBcat_Primary <- subset(HGFMETBcat_Multi1, idents = c("Primary"))
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


####Subcluster epi TriplevDouble####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/TriplevDouble.combined1.epi1")

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

#Cell Cycle Regression epi
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/TriplevDouble.combined1.epi1/TriplevDouble.combined1.epi2")

TriplevDouble.combined1.epi2 <- TriplevDouble.combined1.epi1
DefaultAssay(TriplevDouble.combined1.epi2) <- "integrated"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2 <- RunPCA(TriplevDouble.combined1.epi2, features = VariableFeatures(TriplevDouble.combined1.epi2))
ElbowPlot(TriplevDouble.combined1.epi2, ndims = 50)

TriplevDouble.combined1.epi2 <- FindNeighbors(TriplevDouble.combined1.epi2, reduction = "pca", dims = 1:22)
TriplevDouble.combined1.epi2 <- FindClusters(TriplevDouble.combined1.epi2, resolution = 2.5)
TriplevDouble.combined1.epi2 <- RunUMAP(TriplevDouble.combined1.epi2, reduction = "pca", dims = 1:22)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TriplevDouble.combined1.epi2) <- "Phase"
tiff(file = "TriplevDouble.combined1.epi2 after Cell Cyle Regression UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

#Featureplots
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
tiff(file = "TriplevDouble.combined1.epi2 celltype marker expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt5", "Krt14", "Krt19", "Krt8", "Ppp1r1b", "Vim", "hMETtg", "Trp63"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2.seurat.Markers <- FindAllMarkers(TriplevDouble.combined1.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi2.seurat.Markers, "TriplevDouble.combined1.epi2.seurat.Markers.csv")

#Heatmap
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
TriplevDouble.combined1.epi2.seuratTop30 <- TriplevDouble.combined1.epi2.seurat.Markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
tiff(file = "TriplevDouble.combined1.epi2 seurat Heatmap Top30.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c(TriplevDouble.combined1.epi2.seuratTop30$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Rename Clusters
Idents(object = TriplevDouble.combined1.epi2) <- "seurat_clusters"
TriplevDouble.combined1.epi2 <- RenameIdents(object = TriplevDouble.combined1.epi2, 
                                             '33'="BE1", 
                                             '6'="BE2",'16'="BE2",'29'="BE2",
                                             '15'="BE3", '5'="BE3", 
                                             '23' = "BE4", '30'="BE4", '28' = "BE4", '34'="BE4", 
                                             '35'="LE1",'21'="LE1", '12'="LE1", '13'="LE1",
                                             '25'="LE2", '10'="LE2",'17'="LE2",
                                             '18' = "LE3", '14' ="LE3",'32' = "LE3", '31' ="LE3",'4' ="LE3", '20'="LE3",
                                             '11' = "LE4", '1'="LE4",
                                               '0'="LE5",'9'="LE5",
                                             '2'="LE6",'36'="LE6",
                                              '3'="LE7", 
                                              '22'="LE8",'8'="LE8", '19'="LE8",
                                             '7'="LE9", '24'="LE9",
                                             '26'="UrLE", 
                                             '27'="OE")  
TriplevDouble.combined1.epi2[["EpiCellTypes"]] <- Idents(object = TriplevDouble.combined1.epi2)

#Assessing cluster separation
install.packages('ape')
library(ape)
TriplevDouble.combined1.epi3 <- TriplevDouble.combined1.epi2
DefaultAssay(TriplevDouble.combined1.epi3) <- "integrated"
Idents(object = TriplevDouble.combined1.epi3) <- "seurat_clusters"
TriplevDouble.combined1.epi3 <- BuildClusterTree(TriplevDouble.combined1.epi3,assay = "integrated")

BiocManager::install("ggtree")
library(ggtree)
myPhyTree <- Tool(object=TriplevDouble.combined1.epi3, slot = "BuildClusterTree")
ggtree(myPhyTree)+geom_tiplab()+theme_tree()+xlim(NA,400)

#Umap Epicelltype
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes UMAP.tiff", width = 5.5, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label.size = 6, 
        cols = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                  "bisque3", "slategray3", "mediumorchid3", 
                 "yellow2", "green4", "black"))
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes split UMAP.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, label.size = 6, split.by = "stim",
        cols = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                 "bisque3", "slategray3", "mediumorchid3", 
                 "yellow2", "green4", "black"))
dev.off()

#Featureplot
DefaultAssay(TriplevDouble.combined1.epi2)<-"RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi2 hMETtg split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Tcf4 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Tcf4"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Axin2 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Axin2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Tmprss2 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Tmprss2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q5", max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Homer2 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Homer2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3,max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Nkx3-1 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Nkx3-1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", min.cutoff = "q5", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Pbsn split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Ar split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()

tiff(file = "TriplevDouble.combined1.epi2 Krt5 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Krt8 split plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", min.cutoff = "q5", keep.scale = "all")
dev.off()
tiff(file = "TriplevDouble.combined1.epi2 Fkbp5 plots.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 800)
FeaturePlot(TriplevDouble.combined1.epi2, reduction = "umap", features = c("Fkbp5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95", keep.scale = "all")
dev.off()

#Cell counts
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.epi2$stim.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi2), TriplevDouble.combined1.epi2$stim, sep = "_")
Idents(object = TriplevDouble.combined1.epi2) <- "stim.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi2))

#Cell counts
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.epi2$Phase.EpiCellTypes <- paste(Idents(TriplevDouble.combined1.epi2), TriplevDouble.combined1.epi2$Phase, sep = "_")
Idents(object = TriplevDouble.combined1.epi2) <- "Phase.EpiCellTypes"
table(Idents(TriplevDouble.combined1.epi2))

#DEGs
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.epi2 <- ScaleData(TriplevDouble.combined1.epi2, features = rownames(TriplevDouble.combined1.epi2))
TriplevDouble.combined1.epi2.allMarkers <- FindAllMarkers(TriplevDouble.combined1.epi2, min.pct = 0.25, logfc.threshold = 0.25, only.pos = TRUE)
write.csv(TriplevDouble.combined1.epi2.allMarkers, "TriplevDouble.combined1.epi2.EpiCellTypes.allMarkers.csv")

#Heatmap
DefaultAssay(TriplevDouble.combined1.epi2) <- "RNA"
TriplevDouble.combined1.epi2Top20 <- TriplevDouble.combined1.epi2.allMarkers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
tiff(file = "TriplevDouble.combined1.epi2 Heatmap Top20.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TriplevDouble.combined1.epi2, features = c(TriplevDouble.combined1.epi2Top20$gene), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Dotplot
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
tiff(file = "TriplevDouble.combined1.epi2 EpiCellTypes markers DotPlot.tiff", width =17.5 , height = 4.8, units = "in", compression = "lzw", res = 800)
DotPlot(TriplevDouble.combined1.epi2, features = c("Tmem171", "Egfl6", "Ncam1", "Clca3a2", "Adm", 
                                                   "Krt15", "Palld", "Ctsl", "Tubb6", "Tpm1",
                                                   "Aqp3", "Lgals7", "Col17a1", "Lamb3", "Pvrl1",
                                                   "Sncg", "Ifi202b", "Gpnmb", "Dapl1", "Gpr87", 
                                                   "Lars2", "AY036118", "Gm42418", "Gm26917", "Hbb-bs",
                                                   "Defa21", "Gm15293", "Defa5", "Cryba4", "Ltf", 
                                                   "Rnf149", "Car2", "Lap3", "Plat", "Tm4sf1", 
                                                   "Coch", "Tgfb2", "Dkk2", "Zeb2", "Apoc1",
                                                   "Crip1", "Tspan8", "Ly6c1", "Btc", "B2m", 
                                                   "Tgm4", "9530053A07Rik", "Gm5615", "Man1a", "Spink8", 
                                                   "Msmb", "Mme", "Apof", "Pcp4", "Agtr1a",
                                                   "Pigr", "Tspan1", "Tnfrsf21", "Dcxr", "Cldn3",
                                                   "Spink1", "Sbpl", "Crabp1", "Col6a3", 'Gucy2g', 
                                                   "Gsdmc2", "Gsdmc3", "Barx2", "Cxcl15", "Krt4", 
                                                   "Serping1", "Igfbp6", "Fbln1", "Serpinf1", "Col1a2"
),cols = c("light grey", "red")) + RotatedAxis()
dev.off()

##Heatmap
#scale.data
Epi_df <- as.data.frame(t(TriplevDouble.combined1.epi2$RNA@scale.data)) 
Epi_df_selected <- cbind( Epi_df$Pcna, Epi_df$Cenpf, Epi_df$Cenpe, Epi_df$Top2a, Epi_df$Mki67,Epi_df$Ccnb2,Epi_df$Ccnb1,Epi_df$Pttg1, Epi_df$Tubb4b,Epi_df$Cdk4,
                         Epi_df$Pbsn, Epi_df$`Nkx3-1`, Epi_df$Fkbp5, Epi_df$Azgp1, Epi_df$Tmprss2, Epi_df$Homer2, 
                          Epi_df$Hpn, Epi_df$Uba52, Epi_df$Sh3kbp1, Epi_df$Stam, Epi_df$Nras, Epi_df$Mapk8, Epi_df$Sos1, Epi_df$Mtor, 
                          Epi_df$Axin2, Epi_df$Cd44, Epi_df$Dkk2, Epi_df$Tcf4, 
                          Epi_df$Rpl3, Epi_df$Rpl12,Epi_df$Rps5, Epi_df$Rps16)
Epi_df_selected <- as.data.frame(Epi_df_selected)
colnames(Epi_df_selected) <- c( "Pcna",	"Cenpf",	"Cenpe", "Top2a", "Mki67", "Ccnb2","Ccnb1", "Pttg1", "Tubb4b", "Cdk4",
                                "Pbsn",	'Nkx3-1', "Fkbp5", "Azgp1",	"Tmprss2", "Homer2",
                                "Hpn",	"Uba52",	"Sh3kbp1",	"Stam",	"Nras",	"Mapk8",	"Sos1",	"Mtor",
                                "Axin2", "Cd44",	"Dkk2", "Tcf4", 
                                 "Rpl3", "Rpl12", "Rps5", "Rps16")
rownames(Epi_df_selected) <- row.names(Epi_df)
write.csv(Epi_df_selected, file = "Epi_df_selected_scaledata.csv")

#meta.data
write.csv(TriplevDouble.combined1.epi2@meta.data, file = "TriplevDouble.combined1.epi2_metadata.csv")

#Env
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggplotify)

df <- read.csv("Heatmap_TriplevDouble.combined1.epi2.csv", header = TRUE, sep = ",")
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

##Heatmap-1
#scale.data
Epi_df <- as.data.frame(t(TriplevDouble.combined1.epi2$RNA@scale.data)) 
Epi_df_selected <- cbind( Epi_df$Pbsn, Epi_df$`Nkx3-1`, Epi_df$Fkbp5, Epi_df$Azgp1, 
                          Epi_df$Hpn, Epi_df$Mtor, Epi_df$Uba52, Epi_df$Mapk8, Epi_df$Spint1, 
                          Epi_df$Axin2, Epi_df$Cd44, Epi_df$Dkk2, Epi_df$Tcf4, Epi_df$Myc,
                          Epi_df$Rpl3, Epi_df$Rpl12,Epi_df$Rps5, Epi_df$Rps16,
                          Epi_df$Pcna, Epi_df$Cenpf, Epi_df$Pttg1, Epi_df$Cdk4 
                          )
Epi_df_selected <- as.data.frame(Epi_df_selected)
colnames(Epi_df_selected) <- c( "Pbsn",	'Nkx3-1', "Fkbp5", "Azgp1",
                                "Hpn", "Mtor", "Uba52", "Mapk8", "Spint1",
                                "Axin2", "Cd44",	"Dkk2", "Tcf4", "Myc",
                                "Rpl3", "Rpl12", "Rps5", "Rps16",
                                "Pcna",	"Cenpf", "Pttg1", "Cdk4")
rownames(Epi_df_selected) <- row.names(Epi_df)
write.csv(Epi_df_selected, file = "Epi_df_selected_scaledata-1.csv")

#meta.data
write.csv(TriplevDouble.combined1.epi2@meta.data, file = "TriplevDouble.combined1.epi2_metadata.csv")

#Env
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggplotify)

df <- read.csv("Heatmap_Triple_and_Double_Epi_EpiCellTypes-3.csv", header = TRUE, sep = ",")
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

####Triple Monocle3####
Idents(object = TriplevDouble.combined1.epi2) <- "stim"
Triple.combined1.epi2 <- subset(TriplevDouble.combined1.epi2, idents = c("Triple"))
Idents(object = Triple.combined1.epi2) <- "EpiCellTypes"
Triple.combined1.BELE <- subset(Triple.combined1.epi2, idents = c("BE1", "BE2", "BE3", "BE4",
                                                                  "LE1", "LE2", 'LE3', 'LE4',
                                                                  "LE5", "LE6", "LE7", "LE8", "LE9"))

##Triple.combined1.BELE re-clustering
DefaultAssay(Triple.combined1.BELE) <- "RNA"
Triple.combined1.BELE <- FindVariableFeatures(Triple.combined1.BELE, selection.method = "vst", nfeatures = 5000)
Triple.combined1.BELE <- ScaleData(Triple.combined1.BELE, verbose = FALSE)
Triple.combined1.BELE <- RunPCA(Triple.combined1.BELE, npcs = 50, verbose = FALSE)
ElbowPlot(Triple.combined1.BELE, ndims = 50)

Triple.combined1.BELE <- FindNeighbors(Triple.combined1.BELE, reduction = "pca", dims = 1:20)
Triple.combined1.BELE <- FindClusters(Triple.combined1.BELE, resolution = 0.5)
Triple.combined1.BELE <- RunUMAP(Triple.combined1.BELE, reduction = "pca", dims = 1:20)
DimPlot(Triple.combined1.BELE, reduction = "umap")

Idents(object = Triple.combined1.BELE) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELE, reduction = "umap", label = TRUE)

#Cell cycle scoring
DefaultAssay(Triple.combined1.BELE) <- "RNA"
all.genes <- rownames(Triple.combined1.BELE)
Triple.combined1.BELE <- ScaleData(Triple.combined1.BELE, features = all.genes)
Triple.combined1.BELE <- CellCycleScoring(Triple.combined1.BELE, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = Triple.combined1.BELE) <- "Phase"
DimPlot(Triple.combined1.BELE, reduction = "umap")

#Cell Cycle regression
Triple.combined1.BELE1 <- Triple.combined1.BELE
DefaultAssay(Triple.combined1.BELE1) <- "integrated"
Triple.combined1.BELE1 <- ScaleData(Triple.combined1.BELE1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.BELE1))
Triple.combined1.BELE1 <- RunPCA(Triple.combined1.BELE1, features = VariableFeatures(Triple.combined1.BELE1))
ElbowPlot(Triple.combined1.BELE1, ndims = 50)

Triple.combined1.BELE1 <- FindNeighbors(Triple.combined1.BELE1, reduction = "pca", dims = 1:12)
Triple.combined1.BELE1 <- FindClusters(Triple.combined1.BELE1, resolution = 0.5)
Triple.combined1.BELE1 <- RunUMAP(Triple.combined1.BELE1, reduction = "pca", dims = 1:12)

Idents(object = Triple.combined1.BELE1) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELE1, reduction = "umap", pt.size = 0.3, label = TRUE)

##Convert Seurat to Monocle3 cell data set class
DefaultAssay(Triple.combined1.BELE1) <- "RNA"
cds <- as.cell_data_set(Triple.combined1.BELE1)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Triple.combined1.BELE1)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.5,
                  alpha = 0.9) 

plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))

tiff(file = "Triple.combined1.BELE re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))
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
                  cell_size = 0.5,
                  alpha = 0.9)
plt

tiff(file = "Triple.combined1.BELE re-clustering pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

#plot genes in pseudotime
Solid_genes <- c("Eif4e2")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$stim %in% c("Triple")]

tiff(file = "Triple.combined1.BELE Eif4e2 in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes", min_expr = 0.5
                         ) + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                                                       "bisque3", "slategray3", "mediumorchid3", 
                                                                       "yellow2", "green4", "black"))
dev.off()

####hMETtgPos vs hMETtgNeg####
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
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

tiff(file = "TriplevDouble.combined1.epi2 hMETPos highlighted stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi2, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"), split.by = "stim")
dev.off()

#Subset hMETPos BELE
Idents(object = TriplevDouble.combined1.epi2) <- "METExp"
TriplevDouble.combined1.epi2METPos <- subset(TriplevDouble.combined1.epi2, idents = c("METPos"))

Idents(object = TriplevDouble.combined1.epi2METPos) <- "EpiCellTypes"
DimPlot(TriplevDouble.combined1.epi2METPos, reduction = "umap", pt.size = 0.3, label = TRUE)

TriplevDouble.combined1.BELEMETPos <- subset(TriplevDouble.combined1.epi2METPos, idents = c("BE1","BE2", "BE3", "BE4",
                                                                                     "LE1", "LE2", "LE3", "LE4",
                                                                                     "LE5", "LE6", "LE7", "LE8", "LE9"))
DimPlot(TriplevDouble.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
Triple.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELEMETPos, idents = c("Triple"))
Idents(object = Triple.combined1.BELEMETPos) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE)

#DEGs_hMETPosEpi_TriplevsDouble
DefaultAssay(TriplevDouble.combined1.BELEMETPos) <- "RNA"
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
all.genes <- rownames(TriplevDouble.combined1.BELEMETPos)
TriplevDouble.combined1.BELEMETPos <- ScaleData(TriplevDouble.combined1.BELEMETPos, features = all.genes)
hMETtgPos_BELE_TriplevDouble.0.Markers <- FindMarkers(TriplevDouble.combined1.BELEMETPos, ident.1 = c("Triple"), 
                                                    ident.2 = c("Double"), min.pct = 0, logfc.threshold = 0)
write.csv(hMETtgPos_BELE_TriplevDouble.0.Markers, "hMETtgPos_BELE_TriplevDouble.0.Markers.csv")

#p.adjust
DEG_hMETtgPos_BELE_TriplevDouble <- read.csv("hMETtgPos_BELE_TriplevDouble.0.Markers.csv") 
DEG_hMETtgPos_BELE_TriplevDouble_pvalue <- DEG_hMETtgPos_BELE_TriplevDouble$p_val
DEG_hMETtgPos_BELE_TriplevDouble_pvalue=as.numeric(DEG_hMETtgPos_BELE_TriplevDouble_pvalue)
DEG_hMETtgPos_BELE_TriplevDouble_BH = p.adjust(DEG_hMETtgPos_BELE_TriplevDouble_pvalue, "BH")
write.csv(DEG_hMETtgPos_BELE_TriplevDouble_BH, "DEG_hMETtgPos_BELE_TriplevDouble_BH.csv")

#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named Symbol
genelist = DEGs_hMETtgPos_BELE_preranked_2$symbol
ortho_DEGs_hMETtgPos_BELE_preranked_2<- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

DEGs_hMETtgPos_BELE_preranked_2 <- inner_join(DEGs_hMETtgPos_BELE_preranked_2,ortho_DEGs_hMETtgPos_BELE_preranked_2,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(DEGs_hMETtgPos_BELE_preranked_2,"DEGs_hMETtgPos_BELE_converted-2.txt")

#Rename
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
TriplevDouble.combined1.BELEMETPos <- RenameIdents(object = TriplevDouble.combined1.BELEMETPos, 
                                                 'Double'="Double",'Triple'="Triple")  
TriplevDouble.combined1.BELEMETPos[["stim"]] <- Idents(object = TriplevDouble.combined1.BELEMETPos)

#Vlnplots

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

DefaultAssay(TriplevDouble.combined1.BELEMETPos) <- "RNA"
Idents(object = TriplevDouble.combined1.BELEMETPos) <- "stim"
tiff(file = "TriplevDouble.combined1.BELEMETPos Xpo1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Xpo1", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Ran Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Ran", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Rpl12 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Rpl12", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Rps16 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Rps16", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Eif4e2 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Eif4e2", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Eif4a1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Eif4a1", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

tiff(file = "TriplevDouble.combined1.BELEMETPos Mapk8 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Mapk8", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1.5) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Hpn Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Hpn", pt.size = 0, cols = c("#3399FF", "#E06666"), y.max = 1.2) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Uba52 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Uba52", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Mtor Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Mtor", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1.2) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Spint1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Spint1", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()


tiff(file = "TriplevDouble.combined1.BELEMETPos Axin2 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Axin2", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Tcf4 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Tcf4", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Myc Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Myc", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Dkk2 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Dkk2", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Cd44 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Cd44", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 4) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

tiff(file = "TriplevDouble.combined1.BELEMETPos Pbsn Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Pbsn", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Fkbp5 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Fkbp5", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Nkx3-1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Nkx3-1", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Tmprss2 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Tmprss2", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Homer2 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Homer2", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Azgp1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Azgp1", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Slc45a3 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Slc45a3", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Plpp1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Plpp1", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

tiff(file = "TriplevDouble.combined1.BELEMETPos Pcna Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Pcna", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 1.5) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Cenpf Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Cenpf", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 0.5) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Pttg1 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Pttg1", pt.size = 0, cols = c("#3399FF",   "#E06666")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.BELEMETPos Cdk4 Vln.tiff", width = 3.7, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.BELEMETPos, features = "Cdk4", pt.size = 0, cols = c("#3399FF",   "#E06666"), y.max = 2.0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

####hMETtgPos Monocle3####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/TriplevDouble.combined1.epi1/BELEPos Monocle3")

Idents(object = TriplevDouble.combined1.BELE.METPos) <- "stim"
Double.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELE.METPos, idents = c("Double"))
Triple.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELE.METPos, idents = c("Triple"))

###Triple
##Triple.BELEMETPos re-clustering
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
Triple.combined1.BELEMETPos <- FindVariableFeatures(Triple.combined1.BELEMETPos, selection.method = "vst", nfeatures = 5000)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, verbose = FALSE)
Triple.combined1.BELEMETPos <- RunPCA(Triple.combined1.BELEMETPos, npcs = 50, verbose = FALSE)
ElbowPlot(Triple.combined1.BELEMETPos, ndims = 50)

Triple.combined1.BELEMETPos <- FindNeighbors(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:14)
Triple.combined1.BELEMETPos <- FindClusters(Triple.combined1.BELEMETPos, resolution = 0.5)
Triple.combined1.BELEMETPos <- RunUMAP(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:14)
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap")

Idents(object = Triple.combined1.BELEMETPos) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap", label = TRUE)

#Cell cycle scoring
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
all.genes <- rownames(Triple.combined1.BELEMETPos)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, features = all.genes)
Triple.combined1.BELEMETPos <- CellCycleScoring(Triple.combined1.BELEMETPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = Triple.combined1.BELEMETPos) <- "Phase"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap")

#Cell Cycle regression
Triple.combined1.BELEMETPos1 <- Triple.combined1.BELEMETPos
DefaultAssay(Triple.combined1.BELEMETPos1) <- "integrated"
Triple.combined1.BELEMETPos1 <- ScaleData(Triple.combined1.BELEMETPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.BELEMETPos1))
Triple.combined1.BELEMETPos1 <- RunPCA(Triple.combined1.BELEMETPos1, features = VariableFeatures(Triple.combined1.BELEMETPos1))
ElbowPlot(Triple.combined1.BELEMETPos1, ndims = 50)

Triple.combined1.BELEMETPos1 <- FindNeighbors(Triple.combined1.BELEMETPos1, reduction = "pca", dims = 1:16)
Triple.combined1.BELEMETPos1 <- FindClusters(Triple.combined1.BELEMETPos1, resolution = 0.5)
Triple.combined1.BELEMETPos1 <- RunUMAP(Triple.combined1.BELEMETPos1, reduction = "pca", dims = 1:16)

Idents(object = Triple.combined1.BELEMETPos1) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELEMETPos1, reduction = "umap", pt.size = 0.3, label = TRUE)

##Convert Seurat to Monocle3 cell data set class
DefaultAssay(Triple.combined1.BELEMETPos4) <- "RNA"
cds <- as.cell_data_set(Triple.combined1.BELEMETPos4)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Triple.combined1.BELEMETPos1)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))

tiff(file = "Triple.combined1.BELEMETPos4 re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))
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

tiff(file = "Triple.combined1.BELEMETPos1 re-clustering pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

#plot genes in pseudotime
Solid_genes <- c("Bmpr1b")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$stim %in% c("Triple")]

tiff(file = "Triple.combined1.BELEMETPos4 Bmpr1b in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes"
) + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                  "bisque3", "slategray3", "mediumorchid3", 
                                  "yellow2", "green4", "black"))
dev.off()

###Double
Idents(object = TriplevDouble.combined1.BELEMETPos1) <- "stim"
TriplevDouble.combined1.BELE.METPos.Double <- subset(TriplevDouble.combined1.BELEMETPos1, idents = c("Double"))

Idents(object = TriplevDouble.combined1.BELE.METPos.Double) <- "EpiCellTypes"
DimPlot(TriplevDouble.combined1.BELE.METPos.Double, reduction = "umap", pt.size = 0.3, label = TRUE)

##Convert Seurat to Monocle3 cell data set class
cds <- as.cell_data_set(TriplevDouble.combined1.BELE.METPos.Double)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(TriplevDouble.combined1.BELE.METPos.Double)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.6,
                  alpha = 1) 

plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))

tiff(file = "TriplevDouble.combined1.BELE.METPos.Double re-clustering trajectory UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                    "bisque3", "slategray3", "mediumorchid3", 
                                    "yellow2", "green4", "black"))
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

tiff(file = "TriplevDouble.combined1.BELE.METPos.Double re-clustering pseudotime UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

#plot genes in pseudotime
Solid_genes <- c("Dkk2")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$stim %in% c("Double")]

tiff(file = "TriplevDouble.combined1.BELE.METPos.Double Dkk2 in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes", min_expr = 0.1,
) + scale_color_manual(values = c("salmon", "skyblue1", "olivedrab2", "brown3", "deeppink1",  "blue", "darkorange", "red", "turquoise3",
                                  "bisque3", "slategray3", "mediumorchid3", 
                                  "yellow2", "green4", "black"))
dev.off()

####LE2 vs LE3 vs LE4 vs Others####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/TriplevDouble.combined1.epi1")

#LE4 gsea bubble plot
data <-  as.data.frame(GSEA_bubble_LE4)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE4v15678.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")
dev.off()

#LE4 gsea bubble plot-1
data <-  as.data.frame(GSEA_bubble_LE4_1)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`-Log10(FDR)`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE4v15678-1.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`-Log10(FDR)`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()

tiff(file = "GSEA bubble plot LE4v15678-2.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`NES`, size = `-Log10(FDR)`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()

#LE3 gsea bubble plot
data <-  as.data.frame(GSEA_bubble_LE3)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE3v15678.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")
dev.off()

#LE3 gsea bubble plot-1
data <-  as.data.frame(GSEA_bubble_LE3_1)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`-Log10(FDR)`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE3v15678-1.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`-Log10(FDR)`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()

tiff(file = "GSEA bubble plot LE3v15678-2.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`NES`, size = `-Log10(FDR)`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()

#LE2 gsea bubble plot
data <-  as.data.frame(GSEA_bubble_LE2)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE2v15678.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")
dev.off()

#LE2 gsea bubble plot-1
data <-  as.data.frame(GSEA_bubble_LE2_1)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`-Log10(FDR)`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE2v15678-1.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`-Log10(FDR)`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()

tiff(file = "GSEA bubble plot LE2v15678-2.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`NES`, size = `-Log10(FDR)`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()


#subset
Idents(object = TriplevDouble.combined1.epi2) <- "EpiCellTypes"
TriplevDouble.combined1.LE <- subset(TriplevDouble.combined1.epi2, idents = c("LE1","LE2", "LE3", "LE4",
                                                                               "LE5", "LE6", "LE7", "LE8", "LE9"))
DimPlot(TriplevDouble.combined1.LE, reduction = "umap", pt.size = 0.3, label = TRUE)

#Rename Clusters
Idents(object = TriplevDouble.combined1.LE) <- "EpiCellTypes"
TriplevDouble.combined1.LE <- RenameIdents(object = TriplevDouble.combined1.LE, 
                                             'LE1'="Others", 'LE5'="Others",'LE6'="Others",
                                             'LE7'="Others", 'LE8' = "Others", 'LE9' = "Others", 'LE2'="LE2", 
                                             'LE3'="LE3",'LE4'="LE4")  
TriplevDouble.combined1.LE[["OthersvLE234"]] <- Idents(object = TriplevDouble.combined1.LE)
DimPlot(TriplevDouble.combined1.LE, reduction = "umap", pt.size = 0.3, label = TRUE)

#Vlnplots

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

DefaultAssay(TriplevDouble.combined1.LE) <- "RNA"
Idents(object = TriplevDouble.combined1.LE) <- "OthersvLE234"
tiff(file = "TriplevDouble.combined1.LE Xpo1 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Xpo1", pt.size = 0, cols = c("grey",  "orange", "blue", "red"), y.max = 1) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Ran Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Ran", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Rpl12 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rpl12", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Rps16 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Rps16", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Myc Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Myc", pt.size = 0, cols = c("grey",  "orange", "blue", "red"), y.max = 3) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Eif4e2 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Eif4e2", pt.size = 0, cols = c("grey",  "orange", "blue", "red"), y.max = 2.5) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Eif4a1 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Eif4a1", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Axin2 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Axin2", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Tcf4 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Tcf4", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Pbsn Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Pbsn", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Fkbp5 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Fkbp5", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Tmprss2 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Tmprss2", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Homer2 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Homer2", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

tiff(file = "TriplevDouble.combined1.LE Hpn Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Hpn", pt.size = 0, cols = c("grey", "blue", "orange", "red"), y.max = 1.5) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Dkk2 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Dkk2", pt.size = 0, cols = c("grey", "blue", "orange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Pbsn Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Pbsn", pt.size = 0, cols = c("grey", "blue", "orange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "TriplevDouble.combined1.LE Fkbp5 Vln.tiff", width = 4, height = 3.5, units = "in", compression = "lzw", res = 800)
VlnPlot(TriplevDouble.combined1.LE, features = "Fkbp5", pt.size = 0, cols = c("grey", "blue", "orange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

####LE2 vs LE3 vs LE4 vs Others in Triple####
#subset
Idents(object = TriplevDouble.combined1.epi2) <- "stim.EpiCellTypes"
Triple.combined1.BELE <- subset(TriplevDouble.combined1.epi2, idents = c("LE1_Triple", "LE2_Triple",
                                                                         "LE3_Triple", "LE4_Triple",
                                                                         "LE5_Triple", "LE6_Triple",
                                                                         "LE7_Triple", "LE8_Triple", "LE9_Triple"))
DimPlot(Triple.combined1.BELE, reduction = "umap", pt.size = 0.3, label = TRUE)

#Rename Clusters
Idents(object = Triple.combined1.BELE) <- "EpiCellTypes"
Triple.combined1.BELE <- RenameIdents(object = Triple.combined1.BELE, 
                                           'LE1'="Others", 'LE5'="Others", 'LE6'="Others", 'LE7'="Others",'LE8'="Others",'LE9'="Others",
                                      'LE2'="LE2",'LE3'="LE3",
                                           'LE4'="LE4")  
Triple.combined1.BELE[["LE234vOthers"]] <- Idents(object = Triple.combined1.BELE)
DimPlot(Triple.combined1.BELE, reduction = "umap", pt.size = 0.3, label = TRUE)

#Vlnplots

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

DefaultAssay(Triple.combined1.BELE) <- "RNA"
Idents(object = Triple.combined1.BELE) <- "LE234vOthers"
tiff(file = "Triple.combined1.BELE Xpo1 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Xpo1", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red"), y.max = 1) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Ran Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Ran", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Rpl12 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Rpl12", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Rps16 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Rps16", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Myc Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Myc", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red"), y.max = 3) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Eif4e2 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Eif4e2", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red"), y.max = 2.5) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Eif4a1 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Eif4a1", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Axin2 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Axin2", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Tcf4 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Tcf4", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Pbsn Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Pbsn", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Fkbp5 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Fkbp5", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Tmprss2 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Tmprss2", pt.size = 0, cols = c("grey",  "orange", "blue", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Homer2 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Homer2", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Spint1 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Spint1", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Uba52 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Uba52", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Hpn Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Hpn", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()
tiff(file = "Triple.combined1.BELE Dkk2 Vln.tiff", width = 3.5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(Triple.combined1.BELE, features = "Dkk2", pt.size = 0, cols = c("grey",  "blue", "darkorange", "red")) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 23, colour = "green2", shape = 95)
dev.off()

#LE2vOthers
DefaultAssay(Triple.combined1.BELE) <- "RNA"
Idents(object = Triple.combined1.BELE) <- "LE234vOthers"
all.genes <- rownames(Triple.combined1.BELE)
Triple.combined1.BELE <- ScaleData(Triple.combined1.BELE, features = all.genes)
Triple_LE2vOthers.0.Markers <- FindMarkers(Triple.combined1.BELE, ident.1 = c("LE2"), 
                                                 ident.2 = c("Others"), min.pct = 0, logfc.threshold = 0)
write.csv(Triple_LE2vOthers.0.Markers, "Triple_LE2vOthers.0.Markers.csv")

#p.adjust
DEG_Triple_LE2vOthers <- read.csv("Triple_LE2vOthers.0.Markers.csv") 
DEG_Triple_LE2vOthers_pvalue <- DEG_Triple_LE2vOthers$p_val
DEG_Triple_LE2vOthers_pvalue=as.numeric(DEG_Triple_LE2vOthers_pvalue)
DEG_Triple_LE2vOthers_BH = p.adjust(DEG_Triple_LE2vOthers_pvalue, "BH")
write.csv(DEG_Triple_LE2vOthers_BH, "DEG_Triple_LE2vOthers_BH.csv")

#LE3vOthers
DefaultAssay(Triple.combined1.BELE) <- "RNA"
Idents(object = Triple.combined1.BELE) <- "LE234vOthers"
all.genes <- rownames(Triple.combined1.BELE)
Triple.combined1.BELE <- ScaleData(Triple.combined1.BELE, features = all.genes)
Triple_LE3vOthers.0.Markers <- FindMarkers(Triple.combined1.BELE, ident.1 = c("LE3"), 
                                           ident.2 = c("Others"), min.pct = 0, logfc.threshold = 0)
write.csv(Triple_LE3vOthers.0.Markers, "Triple_LE3vOthers.0.Markers.csv")

#p.adjust
DEG_Triple_LE3vOthers <- read.csv("Triple_LE3vOthers.0.Markers.csv") 
DEG_Triple_LE3vOthers_pvalue <- DEG_Triple_LE3vOthers$p_val
DEG_Triple_LE3vOthers_pvalue=as.numeric(DEG_Triple_LE3vOthers_pvalue)
DEG_Triple_LE3vOthers_BH = p.adjust(DEG_Triple_LE3vOthers_pvalue, "BH")
write.csv(DEG_Triple_LE3vOthers_BH, "DEG_Triple_LE3vOthers_BH.csv")

#LE4vOthers
DefaultAssay(Triple.combined1.BELE) <- "RNA"
Idents(object = Triple.combined1.BELE) <- "LE234vOthers"
all.genes <- rownames(Triple.combined1.BELE)
Triple.combined1.BELE <- ScaleData(Triple.combined1.BELE, features = all.genes)
Triple_LE4vOthers.0.Markers <- FindMarkers(Triple.combined1.BELE, ident.1 = c("LE4"), 
                                           ident.2 = c("Others"), min.pct = 0, logfc.threshold = 0)
write.csv(Triple_LE4vOthers.0.Markers, "Triple_LE4vOthers.0.Markers.csv")

#p.adjust
DEG_Triple_LE4vOthers <- read.csv("Triple_LE4vOthers.0.Markers.csv") 
DEG_Triple_LE4vOthers_pvalue <- DEG_Triple_LE4vOthers$p_val
DEG_Triple_LE4vOthers_pvalue=as.numeric(DEG_Triple_LE4vOthers_pvalue)
DEG_Triple_LE4vOthers_BH = p.adjust(DEG_Triple_LE4vOthers_pvalue, "BH")
write.csv(DEG_Triple_LE4vOthers_BH, "DEG_Triple_LE4vOthers_BH.csv")

#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named symbol
genelist = DEGs_Triple_LE2vOthers_preranked$symbol
ortho_DEGs_Triple_LE2vOthers_preranked <- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

DEGs_Triple_LE2vOthers_preranked <- inner_join(DEGs_Triple_LE2vOthers_preranked,ortho_DEGs_Triple_LE2vOthers_preranked,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(DEGs_Triple_LE2vOthers_preranked,"DEGs_Triple_LE2vOthers_converted.txt")

#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named symbol
genelist = DEGs_Triple_LE3vOthers_preranked$symbol
ortho_DEGs_Triple_LE3vOthers_preranked <- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

DEGs_Triple_LE3vOthers_preranked <- inner_join(DEGs_Triple_LE3vOthers_preranked,ortho_DEGs_Triple_LE3vOthers_preranked,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(DEGs_Triple_LE3vOthers_preranked,"DEGs_Triple_LE3vOthers_converted.txt")

#Mouse to human ortholog conversion
library(babelgene)
# Must have a column named symbol
genelist = DEGs_Triple_LE4vOthers_preranked$symbol
ortho_DEGs_Triple_LE4vOthers_preranked <- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

DEGs_Triple_LE4vOthers_preranked <- inner_join(DEGs_Triple_LE4vOthers_preranked,ortho_DEGs_Triple_LE4vOthers_preranked,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(DEGs_Triple_LE4vOthers_preranked,"DEGs_Triple_LE4vOthers_converted.txt")

#LE4 gsea bubble plot
data <-  as.data.frame(GSEA_bubble_plot_LE4vOthers)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE4vOthers.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")
dev.off()

#LE3 gsea bubble plot
data <-  as.data.frame(GSEA_bubble_plot_LE2vOthers)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot LE2vOthers.tiff", width = 2.4, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = `-Log10(FDR)`, y =`GSEA pathways`, 
                        color =`NES`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(-0.5,4.5)) +
  ggtitle("GSEA")
dev.off()

#LE4vOthers
DefaultAssay(Triple.combined1.BELE) <- "RNA"
Idents(object = Triple.combined1.BELE) <- "EpiCellTypes"
all.genes <- rownames(Triple.combined1.BELE)
Triple.combined1.BELE <- ScaleData(Triple.combined1.BELE, features = all.genes)
Triple_LE4vLE2.0.Markers <- FindMarkers(Triple.combined1.BELE, ident.1 = c("LE4"), 
                                           ident.2 = c("LE2"), min.pct = 0, logfc.threshold = 0)
write.csv(Triple_LE4vLE2.0.Markers, "Triple_LE4vLE2.0.Markers.csv")

#p.adjust
DEG_Triple_LE4vLE2 <- read.csv("Triple_LE4vLE2.0.Markers.csv") 
DEG_Triple_LE4vLE2_pvalue <- DEG_Triple_LE4vLE2$p_val
DEG_Triple_LE4vLE2_pvalue=as.numeric(DEG_Triple_LE4vLE2_pvalue)
DEG_Triple_LE4vLE2_BH = p.adjust(DEG_Triple_LE4vLE2_pvalue, "BH")
write.csv(DEG_Triple_LE4vLE2_BH, "DEG_Triple_LE4vLE2_BH.csv")



####RNAseq####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/RAJEEV/Heatmap")

####Plot####
#heatmap
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(pheatmap)
library(viridis)
library(RColorBrewer)
#This plot will use all DEGs, without filter
#create a dataframe of common genes
common_genes <- left_join(DEGs_IntactVsWT,DEGs_CasVsWT_significant,by = "gene", suffix = c("_IntactVsWT", "_CasVsWT")) %>%
  left_join(., DEGs_CasVsIntact, by='gene',suffix = c(".,", "_CasVsIntact"))
common_genes <- unique(common_genes)
common_genes <- common_genes %>% select(starts_with (c("gene","logFC")))
common_genes <- rename(common_genes, logFC_CasVsIntact = logFC) %>% column_to_rownames( var = "gene")
view(common_genes)
write.csv(common_genes, "common_genes_final.csv")
#create matrix file
cas_vs_intact <- unique(cas_vs_intact_degs)
cas_vs_intact1 <- cas_vs_intact %>% select(starts_with (c("gene","logFC","Group")))
cas_vs_intact <- cas_vs_intact %>% select(starts_with (c("gene","logFC")))
cas_vs_intact <-  column_to_rownames(cas_vs_intact, var = "gene")
view(cas_vs_intact)

mat <- as.matrix(cas_vs_intact) 
mat_breaks <- seq(min(mat), max(mat), length.out = 10)
phtmap <- pheatmap::pheatmap(cas_vs_intact, cluster_cols = F, cluster_rows = F, 
                             color = colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100),
                             legend = TRUE, annotation_names_row = TRUE, annotation_names_col = TRUE,
                            show_rownames = F, show_colnames = T,width = 2,height = 8)

tiff(file = "Triple_pheatmap.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
phtmap
dev.off()

rownames(mat) = cas_vs_intact$gene
colnames(mat) = c("Castrated","Intact")


htmap <- Heatmap(mat, cluster_columns = F, show_column_dend = F, show_row_dend = F,
                 show_column_names = T,show_row_names = F,
                 col = colorRamp2(c(-3,0,3), c("darkblue","white","darkred")),name = "Fold change",
                 width = unit (5,"cm"), height = unit(10,"cm"),
                 heatmap_legend_param = list(legend_direction = "horizontal",
                                             legend_positin = "bottom",legend_width = unit(2, "cm"))) 
htmap
tiff(file = "Triple_heatmap.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
htmap
dev.off()



mat <- as.matrix(common_genes_final) 
pheatmap::pheatmap(common_genes_final, cluster_cols = F, cluster_rows = T)
rownames(mat) = common_genes_final$gene
colnames(mat) = colnames(common_genes_final)
htmap <- Heatmap(mat, cluster_columns = FALSE, show_column_dend = F,show_column_names = T,show_row_names = F,
                 col = colorRamp2(c(-3,0,3), c("darkblue","white","darkred")),name = "Fold change",
                 width = unit (5,"cm"), height = unit(10,"cm"),
                 heatmap_legend_param = list(legend_direction = "horizontal",
                                             legend_positin = "bottom",legend_width = unit(2, "cm"))) 
htmap
tiff(file = "Triple_heatmap.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
htmap
dev.off()
# library
#venn diagram
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(eulerr)
#This plot will use significant DEG's. i,e filtered genes(p<0.05,logFC< -1 & logFC>1)
#Make the plot

x = list(Cas = Cas_vs_WT_3vs3_ranked$gene ,Pri = Primary_vs_WT_2vs3_ranked$gene,Difr = cas_vs_pri_3vs2_ranked$gene)

eulvenn <- plot(euler(x), shape = "ellipse",
                quantities = T,fintsize = 20,labels = F,imagetype="pdf",
                fill = "transparent",edges = c("red","green","purple"),
                height = 480 , width = 480 ,resolution = 300,compression = "lzw",
                lwd = 4,rotation = 10,cex = 0.5,cat.cex = 0.3,legend = T,
                cat.pos = c(-27, 27, 135),cat.dist = c(0.055, 0.055, 0.085))

tiff(file = "Triple_venn.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
eulvenn
dev.off()
#spiderplot
library(fmsb)
# To use the fmsb package, add 2 lines to the data frame: the max and min of each variable to show on the plot!
data <- rbind(rep(20,5) , rep(0,5) , data)

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart( data  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)

####Seurat object to 10X files####
BiocManager::install("DropletUtils")
library("DropletUtils")

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/LabelledMulti")
list.files()
data.seurat <- readRDS("scRNA.rds")

unique(data.seurat$orig.ident)

data.seurat.list  <- Seurat::SplitObject(data.seurat, split.by = "orig.ident")
sample.names <- unique(scRNA$orig.ident)
head(sample.names)

demultiplex_convert_to_10x <- function(obj, samples) {
  if(class(data.seurat) != "Seurat") {
    message("WARNING: this rds file does not contain a Seurat object! STOP RUNNING THIS SCRIPT")
    message("Check the data type by running:")
    message("class(data.seurat)")
    stop()
  }
  if(!dir.exists(file.path(getwd(), "demultiplexed"))) {
    dir.create(file.path(getwd(), "demultiplexed"))
  } else {
    print("WARNING! A demultiplexed directory already exists")
    return()
  }
  for (i in 1:length(samples)) {
    print(paste0("Converting sample ", samples[i]))
    obj.sub <- obj[[samples[i]]]
    DropletUtils::write10xCounts(path = paste0(getwd(),"/demultiplexed/",samples[i]), x = obj.sub[["RNA"]]@data, type = "sparse", version="3")
  }
}

demultiplex_convert_to_10x(obj = data.seurat.list, samples = sample.names)