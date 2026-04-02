#MycARKOvCtrl SCseq Work Flow

#### Add necessary tools to library ####

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Ctrl2vARKO1/New")

Idents(object = CtrlvARKO.combined) <- "stim"

tiff(file = "CtrlvARKO.combined stim UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 
dev.off()

Idents(object = CtrlvARKO.combined) <- "groups"

tiff(file = "CtrlvARKO.combined seurat UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3) 
dev.off()

tiff(file = "CtrlvARKO.combined seurat UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, split.by = "stim") 
dev.off()

Idents(object = CtrlvARKO.combined) <- "stim"
Ctrl.combined <- subset(CtrlvARKO.combined, idents = c("Ctrl"))
ARKO.combined <- subset(CtrlvARKO.combined, idents = c("ARKO"))

#Run the standard workflow for visualization and clustering
Ctrl <- ScaleData(Ctrl, verbose = FALSE)
Ctrl <- RunPCA(Ctrl, npcs = 30, verbose = FALSE)
ElbowPlot(Ctrl, ndims = 50)
# t-SNE and Clustering
Ctrl <- FindNeighbors(Ctrl, reduction = "pca", dims = 1:20)
Ctrl <- FindClusters(Ctrl, resolution = 0.5)
Ctrl <- RunTSNE(Ctrl, reduction = "pca", dims = 1:20)
Ctrl <- RunUMAP(Ctrl, reduction = "pca", dims = 1:20)
DimPlot(Ctrl, reduction = "umap", pt.size = 0.3) 

Idents(object = Ctrl) <- "stim"

tiff(file = "Ctrl UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Ctrl, reduction = "umap", pt.size = 0.3, cols = c("light grey")) 
dev.off()

#Run the standard workflow for visualization and clustering
MycARKO <- ScaleData(MycARKO, verbose = FALSE)
MycARKO <- RunPCA(MycARKO, npcs = 30, verbose = FALSE)
ElbowPlot(MycARKO, ndims = 50)
# t-SNE and Clustering
MycARKO <- FindNeighbors(MycARKO, reduction = "pca", dims = 1:20)
MycARKO <- FindClusters(MycARKO, resolution = 0.5)
MycARKO <- RunTSNE(MycARKO, reduction = "pca", dims = 1:20)
MycARKO <- RunUMAP(MycARKO, reduction = "pca", dims = 1:20)
DimPlot(MycARKO, reduction = "umap", pt.size = 0.3) 

Idents(object = MycARKO) <- "stim"

tiff(file = "MycARKO UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(MycARKO, reduction = "umap", pt.size = 0.3, cols = c("darkblue")) 
dev.off()

CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, '0' = "BE", '1' = "BE", '2' = "FB", '3' = "FB", '4' = "LE", '5' = "Leu", '6' = "LE", '7' = "Endo", '8' = "LE", '9' = "Leu", '10' = "LE", '11' = "Leu", '12' = "LE", '13' = "Peri", '14' ="Leu", '15' = "")


DefaultAssay(CtrlvARKO.combined) <- "RNA"
tiff(file = "combined Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Merging Datasets ####

#Stash old idents
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)

#Set Current idents
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = MycARKO) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
MycARKO$stim <- "ARKO"
CtrlvARKO.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, MycARKO), dims = 1:20)
CtrlvARKO.combined <- IntegrateData(anchorset = CtrlvARKO.anchors, dims = 1:20)
DefaultAssay(CtrlvARKO.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, verbose = FALSE)
CtrlvARKO.combined <- RunPCA(CtrlvARKO.combined, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined, ndims = 50)
# t-SNE and Clustering
CtrlvARKO.combined <- FindNeighbors(CtrlvARKO.combined, reduction = "pca", dims = 1:19)
CtrlvARKO.combined <- FindClusters(CtrlvARKO.combined, resolution = 0.5)
CtrlvARKO.combined <- RunTSNE(CtrlvARKO.combined, reduction = "pca", dims = 1:19)
CtrlvARKO.combined <- RunUMAP(CtrlvARKO.combined, reduction = "pca", dims = 1:19)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3) 
DimPlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", pt.size = 0.3) 

#FeaturePlot
DefaultAssay(CtrlvARKO.combined) <- "RNA"
tiff(file = "combined Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Epcam.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Vim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Myh11.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Fbln1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Reclustering Epithelial cells ####

Idents(object = CtrlvARKO.combined) <- "seurat_clusters"
CtrlvARKO.combined.Epi <- subset(CtrlvARKO.combined, idents = c("0", "3", "18", "6", "5", "7", "13"))
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"
DefaultAssay(CtrlvARKO.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, verbose = FALSE)
CtrlvARKO.combined.Epi <- RunPCA(CtrlvARKO.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi, ndims = 50)
#Umap and Clustering
CtrlvARKO.combined.Epi <- FindNeighbors(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:17)
CtrlvARKO.combined.Epi <- FindClusters(CtrlvARKO.combined.Epi, resolution = 0.5)
CtrlvARKO.combined.Epi <- RunTSNE(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:17)
CtrlvARKO.combined.Epi <- RunUMAP(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:17)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", label = TRUE)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", split.by = "stim")

#FeaturePlot
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"

tiff(file = "combined.Epi Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt14"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt18.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt18"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Myc.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("MYC-transgene"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "combined.Epi Myc-transgene split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("MYC-transgene"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Myc split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Myc"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Tcf7l2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Tcf7l2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Ccnd1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Ccnd1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Axin2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Axin2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "combined.Epi Fos split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Fos"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Igf1r split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Igf1r"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Jun split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Jun"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Junb split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Junb"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "combined.Epi Jund split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Jund"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Fosl2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Fosl2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Igf2r split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Igf2r"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined.Epi) <- "stim"
CtrlvARKO.combined.Epi$stim.seurat_clusters <- paste(Idents(CtrlvARKO.combined.Epi), CtrlvARKO.combined.Epi$seurat_clusters, sep = "_")
Idents(object = CtrlvARKO.combined.Epi) <- "stim.seurat_clusters"
table(Idents(CtrlvARKO.combined.Epi))

#DEGs
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, features = rownames(CtrlvARKO.combined.Epi))
combined.Epi.markers <- FindAllMarkers(CtrlvARKO.combined.Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.Epi.markers, "combined.Epi.Markers.csv")
Idents(object = CtrlvARKO.combined.Epi) <- "stim"
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, features = rownames(CtrlvARKO.combined.Epi))
CtrlvARKO.combined.Epi.markers <- FindAllMarkers(CtrlvARKO.combined.Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.Epi.markers, "CtrlvARKO.combined.Epi.Markers.csv")

#### Myc+ CtrlvARKO in Integrated ####
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"

#Add MYC info
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"

CtrlvARKO.combined.Epi2 <- subset(x=CtrlvARKO.combined.Epi, subset = `MYC-transgene` > 0)
DefaultAssay(CtrlvARKO.combined.Epi2) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi2 <- ScaleData(CtrlvARKO.combined.Epi2, verbose = FALSE)
CtrlvARKO.combined.Epi2 <- RunPCA(CtrlvARKO.combined.Epi2, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi2, ndims = 50)
# t-SNE and Clustering
CtrlvARKO.combined.Epi2 <- FindNeighbors(CtrlvARKO.combined.Epi2, reduction = "pca", dims = 1:19)
CtrlvARKO.combined.Epi2 <- FindClusters(CtrlvARKO.combined.Epi2, resolution = 0.5)
CtrlvARKO.combined.Epi2 <- RunTSNE(CtrlvARKO.combined.Epi2, reduction = "pca", dims = 1:19)
CtrlvARKO.combined.Epi2 <- RunUMAP(CtrlvARKO.combined.Epi2, reduction = "pca", dims = 1:19)
DimPlot(CtrlvARKO.combined.Epi2, reduction = "umap", pt.size = 0.3, label = TRUE) 
DimPlot(CtrlvARKO.combined.Epi2, reduction = "umap", split.by = "stim", pt.size = 0.3) 

#FeaturePlot
DefaultAssay(CtrlvARKO.combined.Epi2) <- "RNA"

tiff(file = "MycPos.Epi Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt14"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi Krt18.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt18"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi MYC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("MYC-transgene"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "MycPos.Epi Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi Pbsn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi Krt5 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MycPos.Epi Nkx3.1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c('Nkx3-1'), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined.Epi2) <- "stim"
CtrlvARKO.combined.Epi2$stim.seurat_clusters <- paste(Idents(CtrlvARKO.combined.Epi2), CtrlvARKO.combined.Epi2$seurat_clusters, sep = "_")
Idents(object = CtrlvARKO.combined.Epi2) <- "stim.seurat_clusters"
table(Idents(CtrlvARKO.combined.Epi2))

#DEGs
DefaultAssay(CtrlvARKO.combined.Epi2) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi2) <- "seurat_clusters"
CtrlvARKO.combined.Epi2 <- ScaleData(CtrlvARKO.combined.Epi2, features = rownames(CtrlvARKO.combined.Epi2))
combined.Epi2.markers <- FindAllMarkers(CtrlvARKO.combined.Epi2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.Epi2.markers, "combined.Epi2.Markers.csv")
Idents(object = CtrlvARKO.combined.Epi2) <- "stim"
CtrlvARKO.combined.Epi2 <- ScaleData(CtrlvARKO.combined.Epi2, features = rownames(CtrlvARKO.combined.Epi2))
CtrlvARKO.combined.Epi2.markers <- FindAllMarkers(CtrlvARKO.combined.Epi2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.Epi2.markers, "CtrlvARKO.combined.Epi2.Markers.csv")

#### Myc+ BE CtrlvARKO in Integrated ####
Idents(object = CtrlvARKO.combined.Epi2) <- "seurat_clusters"
CtrlvARKO.combined.BE2 <- subset(CtrlvARKO.combined.Epi2, idents = c("0", "3"))
DefaultAssay(CtrlvARKO.combined.BE2) <- "RNA"

Idents(object = CtrlvARKO.combined.BE2) <- "stim"
CtrlvARKO.combined.BE2 <- ScaleData(CtrlvARKO.combined.BE2, features = rownames(CtrlvARKO.combined.BE2))
CtrlvARKO.combined.BE2.markers <- FindAllMarkers(CtrlvARKO.combined.BE2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.BE2.markers, "CtrlvARKO.combined.BE2.Markers.csv")




#### Myc+ LE CtrlvARKO in Integrated ####
Idents(object = CtrlvARKO.combined.Epi2) <- "seurat_clusters"
CtrlvARKO.combined.Epi3 <- subset(CtrlvARKO.combined.Epi2, idents = c("1", "2", "4", "5", "6", "7", "8", "9", "10"))
DefaultAssay(CtrlvARKO.combined.Epi3) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi3 <- ScaleData(CtrlvARKO.combined.Epi3, verbose = FALSE)
CtrlvARKO.combined.Epi3 <- RunPCA(CtrlvARKO.combined.Epi3, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi3, ndims = 50)
# t-SNE and Clustering
CtrlvARKO.combined.Epi3 <- FindNeighbors(CtrlvARKO.combined.Epi3, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Epi3 <- FindClusters(CtrlvARKO.combined.Epi3, resolution = 0.5)
CtrlvARKO.combined.Epi3 <- RunTSNE(CtrlvARKO.combined.Epi3, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Epi3 <- RunUMAP(CtrlvARKO.combined.Epi3, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.Epi3, reduction = "umap", pt.size = 0.3, label = TRUE) 
DimPlot(CtrlvARKO.combined.Epi3, reduction = "umap", split.by = "stim", pt.size = 0.3) 

#Cell counts
Idents(object = CtrlvARKO.combined.Epi3) <- "stim"
CtrlvARKO.combined.Epi3$stim.seurat_clusters <- paste(Idents(CtrlvARKO.combined.Epi3), CtrlvARKO.combined.Epi3$seurat_clusters, sep = "_")
Idents(object = CtrlvARKO.combined.Epi3) <- "stim.seurat_clusters"
table(Idents(CtrlvARKO.combined.Epi3))

#DEGs
DefaultAssay(CtrlvARKO.combined.Epi3) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi3) <- "seurat_clusters"
CtrlvARKO.combined.Epi3 <- ScaleData(CtrlvARKO.combined.Epi3, features = rownames(CtrlvARKO.combined.Epi3))
combined.Epi3.markers <- FindAllMarkers(CtrlvARKO.combined.Epi3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.Epi3.markers, "combined.Epi3.Markers.csv")
Idents(object = CtrlvARKO.combined.Epi3) <- "stim"
CtrlvARKO.combined.Epi3 <- ScaleData(CtrlvARKO.combined.Epi3, features = rownames(CtrlvARKO.combined.Epi3))
CtrlvARKO.combined.Epi3.markers <- FindAllMarkers(CtrlvARKO.combined.Epi3, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.Epi3.markers, "CtrlvARKO.combined.Epi3.Markers.csv")

tiff(file = "MYCPosLE Cd44 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Cd44"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Ly6a split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Ly6a"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Plaur split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Plaur"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Aldh1a3 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Aldh1a3"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Ctnnb1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Ctnnb1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Psca split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Psca"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Fosl1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Fosl1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Ccnd1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Ccnd1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Tcf7l2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Tcf7l2"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Src split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Src"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Jak1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Jak1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Igf1r split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Igf1r"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPosLE Junb split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi3, reduction = "umap", features = c("Junb"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Pseudotime of Ctrl Epi ####
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"
CtrlvARKO.combined.Epi <- subset(CtrlvARKO.combined.Epi, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))


Idents(object = CtrlvARKO.combined.Epi) <- "stim"
CtrlEpi <- subset(CtrlvARKO.combined.Epi, idents = c("Ctrl"))
 
Idents(object = CtrlEpi) <- "seurat_clusters"
DefaultAssay(CtrlEpi) <- "RNA"
CtrlEpiPseudo <- as.CellDataSet(CtrlEpi)
CtrlEpiPseudo <- detectGenes(CtrlEpiPseudo, min_expr = 0.1)
print(head(fData(CtrlEpiPseudo)))

expressed_genes <- row.names(subset(fData(CtrlEpiPseudo),
                                    num_cells_expressed >= 10))

pData(CtrlEpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(CtrlEpiPseudo))
CtrlEpiPseudo <- CtrlEpiPseudo[,pData(CtrlEpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(CtrlEpiPseudo), geom =
        "density")

CtrlEpiPseudo <- estimateSizeFactors(CtrlEpiPseudo)
CtrlEpiPseudo <- estimateDispersions(CtrlEpiPseudo)

disp_table <- dispersionTable(CtrlEpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
CtrlEpiPseudo <- setOrderingFilter(CtrlEpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(CtrlEpiPseudo)

#EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(CtrlEpiPseudo, return_all = F) # norm_method='log'

CtrlEpiPseudo <- reduceDimension(CtrlEpiPseudo, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
CtrlEpiPseudo <- clusterCells(CtrlEpiPseudo, num_clusters = 2)

plot_cell_clusters(CtrlEpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(CtrlEpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
CtrlEpiPseudo <- setOrderingFilter(CtrlEpiPseudo, ordering_genes)
plot_ordering_genes(CtrlEpiPseudo)

CtrlEpiPseudo <- reduceDimension(CtrlEpiPseudo, max_components = 2,
                             method = 'DDRTree')

CtrlEpiPseudo <- orderCells(CtrlEpiPseudo)

GM_state <- function(CtrlEpiPseudo){
  if (length(unique(pData(CtrlEpiPseudo)$State)) > 1){
    T0_counts <- table(pData(CtrlEpiPseudo)$State, pData(CtrlEpiPseudo)$seurat_clusters)[,"2"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

CtrlEpiPseudo <- orderCells(CtrlEpiPseudo, root_state = GM_state(CtrlEpiPseudo))

#Visualization
plot_cell_trajectory(CtrlEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

#### Pseudotime of ARKO Epi ####
Idents(object = CtrlvARKO.combined.Epi) <- "stim"
ARKOEpi <- subset(CtrlvARKO.combined.Epi, idents = c("ARKO"))

Idents(object = ARKOEpi) <- "seurat_clusters"
DefaultAssay(ARKOEpi) <- "RNA"
ARKOEpiPseudo <- as.CellDataSet(ARKOEpi)
ARKOEpiPseudo <- detectGenes(ARKOEpiPseudo, min_expr = 0.1)
print(head(fData(ARKOEpiPseudo)))

expressed_genes <- row.names(subset(fData(ARKOEpiPseudo),
                                    num_cells_expressed >= 10))

pData(ARKOEpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARKOEpiPseudo))
ARKOEpiPseudo <- ARKOEpiPseudo[,pData(ARKOEpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARKOEpiPseudo), geom =
        "density")

ARKOEpiPseudo <- estimateSizeFactors(ARKOEpiPseudo)
ARKOEpiPseudo <- estimateDispersions(ARKOEpiPseudo)

disp_table <- dispersionTable(ARKOEpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARKOEpiPseudo <- setOrderingFilter(ARKOEpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(ARKOEpiPseudo)

#EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(ARKOEpiPseudo, return_all = F) # norm_method='log'

ARKOEpiPseudo <- reduceDimension(ARKOEpiPseudo, max_components = 2, num_dim = 20,
                                 reduction_method = 'tSNE', verbose = T)
ARKOEpiPseudo <- clusterCells(ARKOEpiPseudo, num_clusters = 2)

plot_cell_clusters(ARKOEpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(ARKOEpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARKOEpiPseudo <- setOrderingFilter(ARKOEpiPseudo, ordering_genes)
plot_ordering_genes(ARKOEpiPseudo)

ARKOEpiPseudo <- reduceDimension(ARKOEpiPseudo, max_components = 2,
                                 method = 'DDRTree')

ARKOEpiPseudo <- orderCells(ARKOEpiPseudo)

GM_state <- function(ARKOEpiPseudo){
  if (length(unique(pData(ARKOEpiPseudo)$State)) > 1){
    T0_counts <- table(pData(ARKOEpiPseudo)$State, pData(ARKOEpiPseudo)$seurat_clusters)[,"2"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARKOEpiPseudo <- orderCells(ARKOEpiPseudo, root_state = GM_state(ARKOEpiPseudo))

#Visualization
plot_cell_trajectory(ARKOEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

#### Pseudotime of Epi ####
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"
EpiPseudo <- as.CellDataSet(CtrlvARKO.combined.Epi)
EpiPseudo <- detectGenes(EpiPseudo, min_expr = 0.1)
print(head(fData(EpiPseudo)))

expressed_genes <- row.names(subset(fData(EpiPseudo),
                                    num_cells_expressed >= 10))

pData(EpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo))
EpiPseudo <- EpiPseudo[,pData(EpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo), geom =
        "density")

EpiPseudo <- estimateSizeFactors(EpiPseudo)
EpiPseudo <- estimateDispersions(EpiPseudo)

disp_table <- dispersionTable(EpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo <- setOrderingFilter(EpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo)

#EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo, return_all = F) # norm_method='log'

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2, num_dim = 20,
                                 reduction_method = 'tSNE', verbose = T)
EpiPseudo <- clusterCells(EpiPseudo, num_clusters = 2)

plot_cell_clusters(EpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(EpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
EpiPseudo <- setOrderingFilter(EpiPseudo, ordering_genes)
plot_ordering_genes(EpiPseudo)

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2,
                                 method = 'DDRTree')

EpiPseudo <- orderCells(EpiPseudo)

GM_state <- function(EpiPseudo){
  if (length(unique(pData(EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo)$State, pData(EpiPseudo)$seurat_clusters)[,"2"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo <- orderCells(EpiPseudo, root_state = GM_state(EpiPseudo))

#Visualization
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

tiff(file = "EpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "CtrlEpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "EpiPseudo Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo stim.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()

tiff(file = "EpiPseudo seurat.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "EpiPseudo seurat split tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()

tiff(file = "CtrlEpiPseudo seurat tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "CtrlEpiPseudo seurat split tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "ARKOEpiPseudo seurat.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "ARKOEpiPseudo seurat split tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()

tiff(file = "EpiPseudo CellMarkers purple.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("MYC-transgene", "Ar", "Pbsn", 'Nkx3-1', "Krt5", "Krt8", "Krt19", "Myc"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "CtrlEpiPseudo CellMarkerspurple.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, markers = c("MYC-transgene", "Ar", "Pbsn", 'Nkx3-1', "Krt5", "Krt8", "Krt19", "Myc"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOEpiPseudo CellMarkers purple.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, markers = c("MYC-transgene", "Ar", "Pbsn", 'Nkx3-1', "Krt5", "Krt8", "Krt19", "Myc"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()

tiff(file = "CtrlEpiPseudo MYC.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, markers = c("MYC-transgene"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOEpiPseudo MYC.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, markers = c("MYC-transgene"), cell_size = 1, use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo MYC.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("MYC-transgene"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()

tiff(file = "EpiPseudo Wntgenes red.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Myc", "Ccnd1", "Cd44", "Plaur", "Tcf712", "Ctnnb1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "CtrlEpiPseudo Wntgenes red.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, markers = c("Myc", "Ccnd1", "Cd44", "Plaur", "Tcf712", "Ctnnb1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOEpiPseudo Wntgenes red.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, markers = c("Myc", "Ccnd1", "Cd44", "Plaur", "Tcf712", "Ctnnb1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "EpiPseudo IGF signaling darkgreen.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Igf1r", "Src", "Jak1", "Fosl1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "darkgreen"))
dev.off()
tiff(file = "CtrlEpiPseudo IGF signaling darkgreen.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(CtrlEpiPseudo, markers = c("Igf1r", "Src", "Jak1", "Fosl1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "darkgreen"))
dev.off()
tiff(file = "ARKOEpiPseudo IGF signaling darkgreen.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOEpiPseudo, markers = c("Igf1r", "Src", "Jak1", "Fosl1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "darkgreen"))
dev.off()
