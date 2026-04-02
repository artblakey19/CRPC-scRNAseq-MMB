#WK HiMyc Ctrl2vARKO1 Workflow

#### Add necessary tools to library ####

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(patchwork)
library(reticulate)
library(monocle)
library(RColorBrewer)

library(GENIE3)
library(AUCell)
library(RcisTarget)
library(SingleCellSignalR)
Integrate

#### Initial Filtering and Clustering ####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final")

###HiMYC-Ctrl### 

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

Ctrlunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/2nd Run 6-4-2020/Sun Lab Analysis/Project_COHP_36595_1_XCTC1_count/outs/filtered_feature_bc_matrix")
Ctrlunfiltered <- CreateSeuratObject(counts = Ctrlunfiltered.data,  min.cells = 3, min.features = 500, project = "Ctrlunfiltered")
Ctrlunfiltered <- NormalizeData(Ctrlunfiltered)

#Initial processing & filtering
Ctrlunfiltered[["percent.mt"]] <- PercentageFeatureSet(Ctrlunfiltered, pattern = "^mt-")

tiff(file = "Ctrl Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Ctrlunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "Ctrl Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Ctrlunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Pre-filteration")
dev.off()
tiff(file = "Ctrl Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Ctrlunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Pre-filteration")
dev.off()

Ctrl <- subset(Ctrlunfiltered, subset = nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 15)
tiff(file = "Ctrl Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "Ctrl Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Post-filteration")
dev.off()
tiff(file = "Ctrl Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Ctrl@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Post-filteration")
dev.off()

Ctrl <- FindVariableFeatures(Ctrl, selection.method = "vst", nfeatures = 5000)
tiff(file = "Ctrl Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(Ctrl)
dev.off()

table(Idents(Ctrlunfiltered))
table(Idents(Ctrl))

#Run the standard workflow for visualization and clustering Ctrl2
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

###HiMYC-ARKO###

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

MycARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/190408_190425_combined/28345_count/outs/filtered_feature_bc_matrix")
MycARKOunfiltered <- CreateSeuratObject(counts = MycARKOunfiltered.data, project = "MycARKOSC", min.cells = 3, min.features = 200)

#Initial processing & filtering
MycARKOunfiltered[["percent.mt"]] <- PercentageFeatureSet(MycARKOunfiltered, pattern = "^mt-")

tiff(file = "MycARKO Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(MycARKOunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "MycARKO Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(MycARKOunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "MycARKO Pre-filteration")
dev.off()
tiff(file = "MycARKO Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(MycARKOunfiltered@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "MycARKO Pre-filteration")
dev.off()

MycARKO <- subset(MycARKOunfiltered, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mt < 15)
tiff(file = "MycARKO Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(MycARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "MycARKO Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(MycARKO@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "MycARKO Post-filteration")
dev.off()
tiff(file = "MycARKO Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(MycARKO@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "MycARKO Post-filteration")
dev.off()

MycARKO <- FindVariableFeatures(MycARKO, selection.method = "vst", nfeatures = 5000)
tiff(file = "MycARKO Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(MycARKO)
dev.off()

table(Idents(MycARKOunfiltered))
table(Idents(MycARKO))

#Run the standard workflow for visualization and clustering ARKO1
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

#Plot
tiff(file = "MycARKO UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(MycARKO, reduction = "umap", pt.size = 0.3, cols = c("darkblue")) 
dev.off()

#
mean(Ctrl$nCount_RNA)
mean(Ctrl$nFeature_RNA)

mean(MycARKO$nCount_RNA)
mean(MycARKO$nFeature_RNA)

#### Merging Datasets #### ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/CtrlvARKO.combined")

Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)

Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = MycARKO) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
MycARKO$stim <- "ARKO"
CtrlvARKO.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, MycARKO), dims = 1:20)
CtrlvARKO.combined <- IntegrateData(anchorset = CtrlvARKO.anchors, dims = 1:20)
DefaultAssay(CtrlvARKO.combined) <- "integrated"

#Run the standard workflow for visualization and clustering Ctrl2vARKO1.combined
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, verbose = FALSE)
CtrlvARKO.combined <- RunPCA(CtrlvARKO.combined, npcs = 30, verbose = FALSE)
CtrlvARKO.combined <- FindNeighbors(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
CtrlvARKO.combined <- FindClusters(CtrlvARKO.combined, resolution = 0.5)
CtrlvARKO.combined <- RunUMAP(CtrlvARKO.combined, reduction = "pca", dims = 1:20)
CtrlvARKO.combined <- RunTSNE(CtrlvARKO.combined, reduction = "pca", dims = 1:20)

Idents(object = CtrlvARKO.combined) <- "stim"

#Plot
tiff(file = "CtrlvARKO.combined stim UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 
dev.off()

#celltype identification
Idents(object = CtrlvARKO.combined) <- "seurat_clusters"

#DEGS
DefaultAssay(CtrlvARKO.combined) <- "RNA"
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, features = rownames(CtrlvARKO.combined))
CtrlvARKO.combined.all.Markers <- FindAllMarkers(CtrlvARKO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.all.Markers, "CtrlvARKO.combined.all.seurat.Markers.csv")

#Featureplots
DefaultAssay(CtrlvARKO.combined) <- "RNA"
tiff(file = "combined Epcam.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Vim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Fbln1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Acta2.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Acta2"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Cdh5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Rgs5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Plp1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Plp1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()


CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, '0' = "BE", '1' = "FB",
                                   '2' = "FB", '3' = "BE", '4' = "LE", '5' = "LE", 
                                   '6' = "FB", '7' = "Leu", '8' = "LE", '9' = "Endo",
                                   '10' = "Leu", '11' = "SM", '12' = "Leu", '13' = "Leu",
                                   '14' = "LE", '15' = "Peri", '16' = "FB", '17' = "Glia",
                                   '18' = "Leu")
CtrlvARKO.combined[["celltype"]] <- Idents(object = CtrlvARKO.combined)

Idents(object = CtrlvARKO.combined) <- "celltype"
tiff(file = "C2vA1 Celltype UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("Red", "Purple", "Blue", "Green", "Orange", "Brown", "Gold", "Black"))
dev.off()

tiff(file = "C2vA1 Celltype UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("Red", "Purple", "Blue", "Green", "Orange", "Brown", "Gold", "Black"))
dev.off()

Idents(object = CtrlvARKO.combined) <- "stim"
CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, 'Ctrl' = "Ctrl", 'ARKO' = "ARKO")
CtrlvARKO.combined[["stim"]] <- Idents(object = CtrlvARKO.combined)

#DEGS
DefaultAssay(CtrlvARKO.combined) <- "RNA"
Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, features = rownames(CtrlvARKO.combined))
CtrlvARKO.combined.all.CvA.Markers <- FindAllMarkers(CtrlvARKO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.all.Markers, "CtrlvARKO.combined.all.CvA.Markers_celltype.csv")


CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, 'BE' = "BE", 'LE' = "LE", 'FB' = "FB", 'SM' = "SM", 'Peri' = "Peri", 'Leu' = "Leu", 'Endo' = "Endo", 'Glia' = "Glia")
CtrlvARKO.combined[["celltype"]] <- Idents(object = CtrlvARKO.combined)

#Dotplot
Idents(object = CtrlvARKO.combined) <- "celltype"
tiff(file = "C2vA1 Celltype dotplot.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DotPlot(CtrlvARKO.combined, features = c("Ar", "Pbsn", "MYC-transgene", "EGFP", "Gli1", 
                                         "Krt5", "Krt14", "Krt15", "Krt17", "Col17a1", 
                                         "Tspan1", "Cldn3", "Agr2", "Bex1", "Tspan8",
                                         "Fbln1", "Apod", "Lum", "Igf1", "Col3a1",
                                         "Acta1", "Acta2", "Myh11", "Cnn1", "Actg2",
                                         "Rgs4", "Rgs5", "Gja4", "S1pr3", "Ndufa4l2",
                                         "Rgs1", "Ccl5", "Tyrobp", "Fcer1g", "Cd52",
                                         "Plvap", "Pecam1", "Aqp1", "Cd93", "Cdh5",
                                         "Kcna1", "Cd59a", "Mpz", "Plp1", "Fabp7"
), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Feature Plots Split and combined
DefaultAssay(CtrlvARKO.combined) <- "RNA"
tiff(file = "combined Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined MYCtg split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c('MYC-transgene'), cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q80")
dev.off()
tiff(file = "combined EGFP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Pbsn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Mki67 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Mki67"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Pcna split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Pcna"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Gli1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

tiff(file = "combined Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined MYCtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined EGFP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Pbsn.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Mki67.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Pcna.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Pcna"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Gli1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined$celltype.stim <- paste(Idents(CtrlvARKO.combined), CtrlvARKO.combined$stim, sep = "_")
Idents(object = CtrlvARKO.combined) <- "celltype.stim"
table(Idents(CtrlvARKO.combined))


#### Subclustering Epi ####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/CtrlvARKO.combined.Epi")

tiff(file = "Epi highlight UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("darkblue", "darkblue", "light grey", "light grey", "light grey", "light grey", "light grey", "light grey"))
dev.off()

Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined.Epi <- subset(CtrlvARKO.combined, idents = c("BE", "LE"))
DefaultAssay(CtrlvARKO.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, verbose = FALSE)
CtrlvARKO.combined.Epi <- RunPCA(CtrlvARKO.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi, ndims = 50)
#Umap and Clustering
CtrlvARKO.combined.Epi <- FindNeighbors(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKO.combined.Epi <- FindClusters(CtrlvARKO.combined.Epi, resolution = 0.5)
CtrlvARKO.combined.Epi <- RunTSNE(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKO.combined.Epi <- RunUMAP(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:16)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", label = TRUE)

#Cell cycle regression
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"
all.genes <- rownames(CtrlvARKO.combined.Epi)
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, features = all.genes)
CtrlvARKO.combined.Epi <- CellCycleScoring(CtrlvARKO.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = CtrlvARKO.combined.Epi) <- "Phase"
CtrlvARKO.combined.Epi <- RenameIdents(object = CtrlvARKO.combined.Epi, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
CtrlvARKO.combined.Epi[["Phase"]] <- Idents(object = CtrlvARKO.combined.Epi)

DimPlot(CtrlvARKO.combined.Epi, reduction = "umap")
tiff(file = "CtrlvARKO.combined.Epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#FFD966", "purple", "#1D762E"))
dev.off()

#Take Cell cycle out 
CtrlvARKO.combined.Epi1 <- CtrlvARKO.combined.Epi
DefaultAssay(CtrlvARKO.combined.Epi1) <- "integrated"
CtrlvARKO.combined.Epi1 <- ScaleData(CtrlvARKO.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(CtrlvARKO.combined.Epi1))
CtrlvARKO.combined.Epi1 <- RunPCA(CtrlvARKO.combined.Epi1, features = VariableFeatures(CtrlvARKO.combined.Epi1))
ElbowPlot(CtrlvARKO.combined.Epi1, ndims = 50)

CtrlvARKO.combined.Epi1 <- FindNeighbors(CtrlvARKO.combined.Epi1, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Epi1 <- FindClusters(CtrlvARKO.combined.Epi1, resolution = 0.4)
CtrlvARKO.combined.Epi1 <- RunUMAP(CtrlvARKO.combined.Epi1, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Epi1 <- RunTSNE(CtrlvARKO.combined.Epi1, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = CtrlvARKO.combined.Epi1) <- "stim"
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap")

CtrlvARKO.combined.Epi1 <- RenameIdents(object = CtrlvARKO.combined.Epi1, 'Ctrl' = "Ctrl", 'ARKO' = "ARKO")
CtrlvARKO.combined.Epi1[["stim"]] <- Idents(object = CtrlvARKO.combined.Epi1)
tiff(file = "CtrlvARKO.combined.Epi1 stim UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("grey", "darkblue"))
dev.off()

Idents(object = CtrlvARKO.combined.Epi1) <- "Phase"
CtrlvARKO.combined.Epi1 <- RenameIdents(object = CtrlvARKO.combined.Epi1, 'G1' = "G1", 'G2M' = "G2M", 'S' = "S")
CtrlvARKO.combined.Epi1[["Phase"]] <- Idents(object = CtrlvARKO.combined.Epi1)
tiff(file = "CtrlvARKO.combined.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#FFD966", "purple", "#1D762E"))
dev.off()

DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = 0, max.cutoff = "q90")

Idents(object = CtrlvARKO.combined.Epi1) <- "seurat_clusters"
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

#rename
Idents(object = CtrlvARKO.combined.Epi1) <- "seurat_clusters"
CtrlvARKO.combined.Epi1 <- RenameIdents(object = CtrlvARKO.combined.Epi1, '2' = "BE1", '0' = "BE2", '1' = "BE3", '4' = "LE1", '3' = "LE2", '5' = "LE3", '7' = "LE4", '8' = "LE5", '6' = "LE6", '9' = "OE")
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 1, label = TRUE)
CtrlvARKO.combined.Epi1[["Epicelltype"]] <- Idents(object = CtrlvARKO.combined.Epi1)

#Dimplot
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
tiff(file = "CtrlvARKO.combined.Epi1 Epicelltype split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Light grey")) 
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Epicelltype UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Light grey")) 
dev.off()

#DEGs_celltype
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
CtrlvARKO.combined.Epi1 <- ScaleData(CtrlvARKO.combined.Epi1, features = rownames(CtrlvARKO.combined.Epi1))
CtrlvARKO.combined.Epi1.celltype.marker <- FindAllMarkers(CtrlvARKO.combined.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.Epi1.celltype.marker, "CtrlvARKO.combined.Epi1.celltype.marker.csv")

#Dotplot
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
tiff(file = "CtrlvARKO.combined.Epi1 Celltype dotplot.tiff", width = 15, height = 6, units = "in", compression = "lzw", res = 800)
DotPlot(CtrlvARKO.combined.Epi1, features = c("Ar", "Pbsn", "MYC-transgene", 
                                              "Pim1", "Cyr61", "Ier3", "Tead1", "Slco2a1", 
                                             "Cp", "Lbp", "Col4a1", "Pltp", "Smoc2",
                                             "Calcb", "Anxa8", "Krt16", "Lypd3", "Adam8",
                                             "Pigr", "Cxcl17", "Pglyrp1", "Cldn10", "Spns2", 
                                             "Crabp1", "Asrgl1", "Defb29", "Gdpd1", "Sbpl",
                                             "Wfdc15b", "Serpinb11", "Fmo1",  "Fmo3", "Folh1",
                                             "Glb1l2", "Fuca1", "Mmp7", "Qsox1", "Tcaf2", 
                                             "Pnliprp1", "Gm37223", "Chn2", "Mt3", "Gsdma",
                                             "Gm26917", "Lars2", "mt-Co3", "AY036118", "Erdr1", 
                                             "Col6a1", "Lox", "Tcf21", "Col5a2", "Col1a1"
), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
CtrlvARKO.combined.Epi1$Epicelltype.stim <- paste(Idents(CtrlvARKO.combined.Epi1), CtrlvARKO.combined.Epi1$stim, sep = "_")
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim"
table(Idents(CtrlvARKO.combined.Epi1))

#Featureplot
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
tiff(file = "CtrlvARKO.combined.Epi1 Krt19 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Krt19"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 1, max.cutoff = 3.5)
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Krt14 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Krt14"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 3.5)
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 MycTg split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("MYC-transgene"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

tiff(file = "CtrlvARKO.combined.Epi1 Igf1r split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Igf1r"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Ctnnb1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Ctnnb1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Ccnd1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Ccnd1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Cd44 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", features = c("Cd44"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.7, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Add MYC info
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
CtrlvARKO.combined.Epi1MYCPos <- subset(x=CtrlvARKO.combined.Epi1,  subset = `MYC-transgene` > 0)
CtrlvARKO.combined.Epi1MYCNeg <- subset(x=CtrlvARKO.combined.Epi1,  subset = `MYC-transgene` == 0)
Idents(object = CtrlvARKO.combined.Epi1MYCPos) <- "MYCPos"
Idents(object = CtrlvARKO.combined.Epi1MYCNeg) <- "MYCNeg"
CtrlvARKO.combined.Epi1MYCPos[["MYCExp"]] <- Idents(object = CtrlvARKO.combined.Epi1MYCPos)
CtrlvARKO.combined.Epi1MYCNeg[["MYCExp"]] <- Idents(object = CtrlvARKO.combined.Epi1MYCNeg)
CtrlvARKO.combined.Epi1MYC <- merge(x = CtrlvARKO.combined.Epi1MYCPos, y = CtrlvARKO.combined.Epi1MYCNeg)
Idents(object = CtrlvARKO.combined.Epi1MYC) <- "MYCExp"
CtrlvARKO.combined.Epi1$MYCExp <- Idents(object = CtrlvARKO.combined.Epi1MYC)
Idents(object = CtrlvARKO.combined.Epi1) <- "MYCExp"
CtrlvARKO.combined.Epi1.MYC <- subset(CtrlvARKO.combined.Epi1, idents = c("MYCPos", "MYCNeg"))

#Cell counts
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim"
CtrlvARKO.combined.Epi1$Epicelltype.stim.MYCExp <- paste(Idents(CtrlvARKO.combined.Epi1), CtrlvARKO.combined.Epi1$MYCExp, sep = "_")
Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype.stim.MYCExp"
table(Idents(CtrlvARKO.combined.Epi1))

#Violin plots
Idents(object = CtrlvARKO.combined.Epi1) <- "celltype"
tiff(file = "CtrlvARKO.combined.Epi1 MYCtg split Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.Epi1, features = "MYC-transgene", group.by = "celltype", split.by = "stim",  pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()


#### Subclustering BE ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/CtrlvARKO.combined.Epi/BE")

Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
CtrlvARKO.combined.BE <- subset(CtrlvARKO.combined.Epi1, idents = c("BE1", "BE2", "BE3"))

Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype"
CtrlvARKO.combined.BE$Epicelltype.stim <- paste(Idents(CtrlvARKO.combined.BE), CtrlvARKO.combined.BE$stim, sep = "_")
Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype.stim"
CtrlvARKO.combined.BE$Epicelltype.stim.MYCExp <- paste(Idents(CtrlvARKO.combined.BE), CtrlvARKO.combined.BE$MYCExp, sep = "_")
Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype.stim.MYCExp"
table(Idents(CtrlvARKO.combined.BE))

#DEGs_celltype
DefaultAssay(CtrlvARKO.combined.BE) <- "RNA"
Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype.stim"
CtrlvARKO.combined.BE <- ScaleData(CtrlvARKO.combined.BE, features = rownames(CtrlvARKO.combined.BE))
CtrlvARKO.combined.BE.Epicelltype.stim.marker <- FindAllMarkers(CtrlvARKO.combined.BE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.BE.Epicelltype.stim.marker, "CtrlvARKO.combined.BE.Epicelltype.stim.marker.csv")

Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype.stim"
CtrlvARKO.combined.BE <- RenameIdents(object = CtrlvARKO.combined.BE, 'BE1_Ctrl' = "BE1_Ctrl", 'BE2_Ctrl' = "BE2_Ctrl", 'BE3_Ctrl' = "BE3_Ctrl", 
                                    'BE1_ARKO' = "BE1_ARKO", 'BE2_ARKO' = "BE2_ARKO", 'BE3_ARKO' = "BE3_ARKO")
CtrlvARKO.combined.BE[["Epicelltype.stim"]] <- Idents(object = CtrlvARKO.combined.BE)

CtrlvARKO.combined.BETop30 <- CtrlvARKO.combined.BE.Epicelltype.stim.marker %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
tiff(file = "CtrlvARKO.combined.BETop30 Heatmap purple.tiff", width =8, height = 12, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.BE, features = c(CtrlvARKO.combined.BETop30$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Dotplot
Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype.stim"
tiff(file = "CtrlvARKO.combined.BE Epicelltype.stim dotplot.tiff", width = 15, height = 6, units = "in", compression = "lzw", res = 800)
DotPlot(CtrlvARKO.combined.BE, features = c("Ar", "Pbsn", "MYC-transgene", "Krt5", "Krt14", "Krt7",
                                            "Igf1r", "Ctnnb1", "Ccnd1", "Cd44", "Tcf7l2",
                                              "Defb50", "Pnrc1", "Taf1d", "Wfdc12", "Rbm3", 
                                              "Lbp", "Dcn", "Apoe", "Col4a1", "Tagln",
                                              "Calcb", "Anxa1", "Anxa8", "Adam8", "Clu",
                                              "Sbp", "Ftl1", "Malat1", "Eif4a1", "Eif1", 
                                              "Hbegf", "Cyr61", "Fos", "Calm2", "Fosb",
                                              "Sfn", "Plaur", "Oaz1",  "Sub1", "Morf4l1"
), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#DEGs_celltype
DefaultAssay(CtrlvARKO.combined.BE) <- "RNA"
Idents(object = CtrlvARKO.combined.BE) <- "Epicelltype"
CtrlvARKO.combined.BE <- ScaleData(CtrlvARKO.combined.BE, features = rownames(CtrlvARKO.combined.BE))
CtrlvARKO.combined.BE.Epicelltype.marker <- FindAllMarkers(CtrlvARKO.combined.BE, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.BE.Epicelltype.marker, "CtrlvARKO.combined.BE.Epicelltype.marker.csv")


####subclustering MYC+ BE####
Idents(object = CtrlvARKO.combined.Epi) <- "MYCExp"
CtrlvARKO.combined.Epi$MYCExp.Epicelltype <- paste(Idents(CtrlvARKO.combined.Epi), CtrlvARKO.combined.Epi$Epicelltype, sep = "_")
Idents(object = CtrlvARKO.combined.Epi) <- "MYCExp.Epicelltype"
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap") 

CtrlvARKO.combined.MYCBE <- subset(CtrlvARKO.combined.Epi, idents = c("MYCPos_BE1", "MYCPos_BE2", "MYCPos_BE3"))

#DEGs
DefaultAssay(CtrlvARKO.combined.MYCBE) <- "RNA"
Idents(object = CtrlvARKO.combined.MYCBE) <- "stim"
CtrlvARKO.combined.MYCBE <- ScaleData(CtrlvARKO.combined.MYCBE, features = rownames(CtrlvARKO.combined.MYCBE))
CtrlvARKO.combined.MYCBE.marker <- FindAllMarkers(CtrlvARKO.combined.MYCBE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.MYCBE.marker, "CtrlvARKO.combined.MYCBE.marker.csv")

CtrlvARKO.combined.MYCBE.IPA.Markers <- FindMarkers(CtrlvARKO.combined.MYCBE, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.MYCBE.IPA.Markers, "CtrlvARKO.combined.MYCBE.IPA.Markers.csv")
CtrlvARKO.combined.MYCBE.GSEA.Markers <- FindMarkers(CtrlvARKO.combined.MYCBE, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(CtrlvARKO.combined.MYCBE.GSEA.Markers, "CtrlvARKO.combined.MYCBE.GSEA.Markers.csv")



#Heatmap
DefaultAssay(CtrlvARKO.combined.MYCBE) <- "RNA"
CtrlvARKO.combined.MYCBE.Top50 <- CtrlvARKO.combined.MYCBE.marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
tiff(file = "CtrlvARKO.combined.MYCBE Heatmap Top50 purple.tiff", width = 4, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.MYCBE, features = c(CtrlvARKO.combined.MYCBE.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violin plots
DefaultAssay(CtrlvARKO.combined.MYCBE) <- "RNA"
Idents(object = CtrlvARKO.combined.MYCBE) <- "stim"
tiff(file = "MYCBE Igf1r Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Igf1r", pt.size = 0)
dev.off()
tiff(file = "MYCBE Irs2 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Irs2", pt.size = 0)
dev.off()
tiff(file = "MYCBE Hras Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Hras", pt.size = 0)
dev.off()
tiff(file = "MYCBE Pik3ca Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Pik3ca", pt.size = 0)
dev.off()
tiff(file = "MYCBE Akt1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Akt1", pt.size = 0)
dev.off()
tiff(file = "MYCBE Grb2 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Grb2", pt.size = 0)
dev.off()

tiff(file = "MYCBE Ctnnb1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Ctnnb1", pt.size = 0)
dev.off()
tiff(file = "MYCBE Ccnd1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Ccnd1", pt.size = 0)
dev.off()
tiff(file = "MYCBE Cd44 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Cd44", pt.size = 0)
dev.off()
tiff(file = "MYCBE Tcf7l2 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Tcf7l2", pt.size = 0)
dev.off()
tiff(file = "MYCBE Mmp7 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.MYCBE, features = "Mmp7", pt.size = 0)
dev.off()

#Boxplot Generation
boxdata = FetchData(CtrlvARKO.combined.MYCBE, c("stim", "Igf1r", "Irs2", "Hras", "Pik3ca", "Akt1", "Grb2", "Ccnd1", "Myc", "Mmp7", "Cd44", "Tcf7l2", "Ctnnb1"))
tail(boxdata,5)

tiff(file = "CtrlvARKO.combined.MYCBE Igf1r Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Igf1r, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Hras Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Hras, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Irs2 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Irs2, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Pik3ca Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Pik3ca, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Akt1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Akt1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Grb2 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Grb2, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Ctnnb1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Ctnnb1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Myc Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Myc, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Cd44 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Cd44, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Mmp7 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Mmp7, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Tcf7l2 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Tcf7l2, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Ccnd1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Ccnd1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()

#Featureplot
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"

tiff(file = "combined Epi Ccnd1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Ccnd1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Ctnnb1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Ctnnb1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Cd44 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Cd44"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Tcf7l2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Tcf7l2"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Mmp7 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Mmp7"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()

tiff(file = "combined Epi Igf1r split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Igf1r"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Hras split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Hras"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Irs2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Irs2"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Pik3ca split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Pik3ca"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Akt1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Akt1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Grb2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Akt1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()

tiff(file = "combined Epi Tmem219 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Tmem219"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Lrp1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Lrp1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()
tiff(file = "combined Epi Ltbp1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Ltbp1"), split.by = "stim", cols = c("light grey", "purple"), pt.size = 0.3, min.cutoff = 0, max.cutoff = 2)
dev.off()

#Comparision with Adam's data

Idents(object = CtrlvARKO.combined.Epi) <- "celltype"
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap")

CtrlvARKO.combined.BE <- subset(CtrlvARKO.combined.Epi, idents = c("BE"))

#Boxplot Generation
DefaultAssay(CtrlvARKO.combined.BE) <- "RNA"

boxdata = FetchData(CtrlvARKO.combined.BE, c("stim", "Hoxb13", "Etv4", "Mmp14", "Bmp7", "Hes1", "Shh"))
tail(boxdata,5)

tiff(file = "CtrlvARKO.combined.BE Hoxb13 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Hoxb13, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()
tiff(file = "CtrlvARKO.combined.BE Etv4 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Etv4, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) 
dev.off()
tiff(file = "CtrlvARKO.combined.BE Mmp14 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Mmp14, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) 
dev.off()
tiff(file = "CtrlvARKO.combined.BE Bmp7 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Bmp7, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()
tiff(file = "CtrlvARKO.combined.BE Hes1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Hes1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()
tiff(file = "CtrlvARKO.combined.BE Shh Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Shh, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()

#Boxplot Generation
DefaultAssay(CtrlvARKO.combined.MYCBE) <- "RNA"

boxdata = FetchData(CtrlvARKO.combined.MYCBE, c("stim", "Hoxb13", "Etv4", "Mmp14", "Bmp7", "Hes1", "Shh"))
tail(boxdata,5)

tiff(file = "CtrlvARKO.combined.MYCBE Hoxb13 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Hoxb13, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Etv4 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Etv4, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) 
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Mmp14 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Mmp14, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) 
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Bmp7 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Bmp7, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Hes1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Hes1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()
tiff(file = "CtrlvARKO.combined.MYCBE Shh Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Shh, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25)
dev.off()

####subclustering MYC+ LE####
Idents(object = CtrlvARKO.combined.Epi) <- "MYCExp.Epicelltype"
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap") 

CtrlvARKO.combined.MYCLE <- subset(CtrlvARKO.combined.Epi, idents = c("MYCPos_LE1", "MYCPos_LE2", "MYCPos_LE3", "MYCPos_LE4", "MYCPos_LE5", "MYCPos_LE6"))

#DEGs
DefaultAssay(CtrlvARKO.combined.MYCLE) <- "RNA"
Idents(object = CtrlvARKO.combined.MYCLE) <- "stim"
CtrlvARKO.combined.MYCLE <- ScaleData(CtrlvARKO.combined.MYCLE, features = rownames(CtrlvARKO.combined.MYCLE))
CtrlvARKO.combined.MYCLE.marker <- FindAllMarkers(CtrlvARKO.combined.MYCLE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.MYCLE.marker, "CtrlvARKO.combined.MYCLE.marker.csv")

CtrlvARKO.combined.MYCLE.IPA.Markers <- FindMarkers(CtrlvARKO.combined.MYCLE, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.MYCLE.IPA.Markers, "CtrlvARKO.combined.MYCLE.IPA.Markers.csv")
CtrlvARKO.combined.MYCLE.GSEA.Markers <- FindMarkers(CtrlvARKO.combined.MYCLE, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(CtrlvARKO.combined.MYCLE.GSEA.Markers, "CtrlvARKO.combined.MYCLE.GSEA.Markers.csv")

#Heatmap
DefaultAssay(CtrlvARKO.combined.MYCLE) <- "RNA"
CtrlvARKO.combined.MYCLE.Top50 <- CtrlvARKO.combined.MYCLE.marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
tiff(file = "CtrlvARKO.combined.MYCLE Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.MYCLE, features = c(CtrlvARKO.combined.MYCLE.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####Pseudotime Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/CtrlvARKO.combined.Epi/Pseudotime/new")

Idents(object = CtrlvARKO.combined.Epi1) <- "Epicelltype"
combined.Epi.Pseudo <- subset(CtrlvARKO.combined.Epi1, idents = c("BE1", "BE2", "BE3", "LE1", "LE2", "LE3", "LE4", "LE5", "LE6"))
Idents(object = combined.Epi.Pseudo) <- "seurat_clusters"
DimPlot(combined.Epi.Pseudo, reduction = "umap")

DefaultAssay(combined.Epi.Pseudo) <- "RNA"
EpiPseudo <- as.CellDataSet(combined.Epi.Pseudo)
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
ordering_genes <- row.names (subset(diff_test_res, qval < 0.05))
EpiPseudo <- setOrderingFilter(EpiPseudo, ordering_genes)
plot_ordering_genes(EpiPseudo)

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2,
                                 method = 'DDRTree')

EpiPseudo <- orderCells(EpiPseudo)

GM_state <- function(EpiPseudo){
  if (length(unique(pData(EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo)$State, pData(EpiPseudo)$State)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo <- orderCells(EpiPseudo, root_state = GM_state(EpiPseudo))

#Visualization

tiff(file = "EpiPseudo Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "EpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

tiff(file = "EpiPseudo Epicelltype.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Epicelltype", show_branch_points = FALSE) + scale_color_manual(values = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue"))
dev.off()

tiff(file = "EpiPseudo Epicelltype stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Epicelltype", show_branch_points = FALSE) + scale_color_manual(values=c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue")) + facet_wrap(~stim, nrow = 1) 
dev.off()


tiff(file = "EpiPseudo Ar.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo MYCtg.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "MYC-transgene", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Krt5.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt5", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Krt5.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt5", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Krt19.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt19", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Krt8.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt8", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Pbsn.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Pbsn", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Igf1r.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Igf1r", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Ctnnb1.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ctnnb1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Ccnd1.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ccnd1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Cd44.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Cd44", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Krt14.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt14", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "EpiPseudo Pcna.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Pcna", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Mki67.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Mki67", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Trp63.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Trp63", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Ly6a.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ly6a", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "EpiPseudo Ar MYCtg Krt5 Krt19 Pbsn Pcna.tiff", width = 12, height = 8.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Ar", "MYC-transgene", "Krt5", "Krt19", "Pbsn", "Pcna"), use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "EpiPseudo Ar MYCtg Krt14 Igf1r Ctnnb1 Ccnd1.tiff", width = 12, height = 36, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Ar", "MYC-transgene", "Krt14", "Igf1r", "Ctnnb1", "Ccnd1"), use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 6) 
dev.off()

tiff(file = "EpiPseudo Igf1r.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Igf1r", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Hras.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Hras", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Irs2.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Irs2", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Akt1.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Akt1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "EpiPseudo Igf1r Hras Irs2 Akt1.tiff", width = 12, height = 8.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Igf1r", "Hras", "Irs2", "Akt1", "Ar", "MYC-transgene"), use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3"))
dev.off()

tiff(file = "EpiPseudo Ctnnb1.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ctnnb1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Ccnd1.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ccnd1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Tcf7l2.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Tcf7l2", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Mmp7.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Mmp7", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Cd44.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Cd44", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo Myc.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Myc", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3")) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "EpiPseudo Ctnnb1 Ccnd1 Tcf7l2 Mmp7 Cd44 Myc.tiff", width = 12, height = 8.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Ctnnb1", "Ccnd1", "Tcf7l2", "Mmp7", "Cd44", "Myc"), use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "purple3"))
dev.off()

#Expression level of Wnt downstreams in plot_genes_branched_pseudotime

EpiPseudo_genes <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Igf1r", "Hras")))
tiff(file = "PINvTumor plot_genes_branched_pseudotime Tcf4 Axin2.tiff", width = 4, height = 6, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(EpiPseudo[EpiPseudo_genes,],
                               branch_point = 3,
                               color_by = "Epicelltype",
                               ncol = 1)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#1D762E", "purple", "red", "#FF9933", "#E06666", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
dev.off()



#### Subclustering Stro ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/CtrlvARKO.combined.FBSM")

Idents(object = CtrlvARKO.combined) <- "celltype"
tiff(file = "Stro highlight UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "light grey", "darkblue", "darkblue", "darkblue", "darkblue", "darkblue", "darkblue"))
dev.off()

CtrlvARKO.combined.Stro <- subset(CtrlvARKO.combined, idents = c("FB", "SM", "Peri", "Leu", "Endo", "Glia"))
DefaultAssay(CtrlvARKO.combined.Stro) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Stro <- ScaleData(CtrlvARKO.combined.Stro, verbose = FALSE)
CtrlvARKO.combined.Stro <- RunPCA(CtrlvARKO.combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Stro, ndims = 50)
#Umap and Clustering
CtrlvARKO.combined.Stro <- FindNeighbors(CtrlvARKO.combined.Stro, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Stro <- FindClusters(CtrlvARKO.combined.Stro, resolution = 0.5)
CtrlvARKO.combined.Stro <- RunTSNE(CtrlvARKO.combined.Stro, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.Stro <- RunUMAP(CtrlvARKO.combined.Stro, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.Stro, reduction = "umap", label = TRUE)

#Featureplots
DefaultAssay(CtrlvARKO.combined.Stro) <- "RNA"
tiff(file = "combined Stro Epcam.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Vim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Fbln1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Acta2.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Acta2"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Cdh5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Rgs5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Stro Plp1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Stro, reduction = "umap", features = c("Plp1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Celltype identification
CtrlvARKO.combined.Stro <- RenameIdents(object = CtrlvARKO.combined.Stro, '1' = "FB1", '2' = "FB2",
                                   '0' = "FB3", '4' = "FB4", '3' = "FB5", '8' = "SM", 
                                   '12' = "Peri", '6' = "Endo", '13' = "Glia", '5' = "Leu",
                                   '11' = "Leu", '10' = "Leu", '7' = "Leu", '9' = "OS")
CtrlvARKO.combined.Stro[["Strocelltype"]] <- Idents(object = CtrlvARKO.combined.Stro)

Idents(object = CtrlvARKO.combined.Stro) <- "Strocelltype"
tiff(file = "FBSM highlight UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Stro, reduction = "umap", pt.size = 0.3, cols = c("red", "purple", "blue", "green", "orange", "brown", "light grey", "light grey", "light grey", "light grey", "light grey"))
dev.off()

####Subclustering FBSM####
Idents(object = CtrlvARKO.combined.Stro) <- "Strocelltype"

CtrlvARKO.combined.FBSM <- subset(CtrlvARKO.combined.Stro, idents = c("FB1", "FB2", "FB3", "FB4", "FB5", "SM"))
CtrlvARKO.combined.FBSM <- RunTSNE(CtrlvARKO.combined.FBSM, reduction = "pca", dims = 1:15)
CtrlvARKO.combined.FBSM <- RunUMAP(CtrlvARKO.combined.FBSM, reduction = "pca", dims = 1:15)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", label = TRUE)

#Rename
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", label = TRUE)
CtrlvARKO.combined.FBSM <- RenameIdents(object = CtrlvARKO.combined.FBSM, 'FB5' = "FB1", 'FB2' = "FB2",
                                        'FB3' = "FB3", 'FB1' = "FB4", 'FB4' = "FB5", 'SM' = "SM")
CtrlvARKO.combined.FBSM[["Strocelltype"]] <- Idents(object = CtrlvARKO.combined.FBSM)


tiff(file = "FBSM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", pt.size = 0.3, cols = c("orange", "purple", "blue", "red", "green", "brown", "light grey", "light grey", "light grey", "light grey", "light grey"))
dev.off()

Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
tiff(file = "FBSM split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("orange", "purple", "blue", "red", "green", "brown", "light grey", "light grey", "light grey", "light grey", "light grey"))
dev.off()

Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
table(Idents(CtrlvARKO.combined.FBSM))

#DEGs_celltype
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
CtrlvARKO.combined.FBSM.celltype.marker <- FindAllMarkers(CtrlvARKO.combined.FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.FBSM.celltype.marker, "CtrlvARKO.combined.FBSM.celltype.marker.csv")

#DEGs_FB5
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
CtrlvARKO.combined.FB5.marker <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "FB5", ident.2 = c("FB1", "FB2", "FB3", "FB4", "SM"), min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.FB5.marker, "CtrlvARKO.combined.FB5.marker.csv")

#DEGs_FB4vothers_Ctrl
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FBSM$Strocelltype.stim <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$stim, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype.stim"

DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype.stim"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
Ctrl.FB4vOthers.marker <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "FB4_Ctrl", ident.2 = c("FB1_Ctrl", "FB2_Ctrl", "FB3_Ctrl", "FB5_Ctrl"), min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Ctrl.FB4vOthers.marker, "Ctrl.FB4vOthers.marker.csv")

#DEGs_FB4_CtrlvARKO
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype.stim"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
FB4_CtrlvARKO.marker <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "FB4_Ctrl", ident.2 = "FB4_ARKO", min.pct = 0.25, logfc.threshold = 0.25)
write.csv(FB4_CtrlvARKO.marker, "FB4_CtrlvARKO.marker.csv")


Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype.stim"

table(Idents(CtrlvARKO.combined.FBSM))

#Dotplot
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
tiff(file = "CtrlvARKO.combined.FBSM Celltype dotplot.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
DotPlot(CtrlvARKO.combined.FBSM, features = c("Ar", "EGFP", "Gli1", 
                                              "Cxcl14", "F2r", "Itgbl1", "Hsd11b1", "Isg15",
                                              "Pcolce2", "Cd248", "Sfrp4", "Tmem100", "Mest",
                                              "Jun", "Klf2", "Hspa1a", "Nr1d1", "Tgfbi",
                                              "Ifi205", "Thbd", "Ccl11", "Ptx3", "Sult1e1",
                                              "Kdm6b", "Maff", "Pim1", "Egr1", "Nr4a1",
                                              "Acta2", "Tagln", "Myl9", "Actg2", "Myh11"
), cols = c("light grey", "red")) + RotatedAxis()
dev.off()


#Add Ar info
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
CtrlvARKO.combined.FBSM.ArPos <- subset(x=CtrlvARKO.combined.FBSM,  subset = `Ar` > 0)
CtrlvARKO.combined.FBSM.ArNeg <- subset(x=CtrlvARKO.combined.FBSM,  subset = `Ar` == 0)
Idents(object = CtrlvARKO.combined.FBSM.ArPos) <- "ArPos"
Idents(object = CtrlvARKO.combined.FBSM.ArNeg) <- "ArNeg"
CtrlvARKO.combined.FBSM.ArPos[["ArExp"]] <- Idents(object = CtrlvARKO.combined.FBSM.ArPos)
CtrlvARKO.combined.FBSM.ArNeg[["ArExp"]] <- Idents(object = CtrlvARKO.combined.FBSM.ArNeg)
CtrlvARKO.combined.FBSM.Ar <- merge(x = CtrlvARKO.combined.FBSM.ArPos, y = CtrlvARKO.combined.FBSM.ArNeg)
Idents(object = CtrlvARKO.combined.FBSM.Ar) <- "ArExp"
CtrlvARKO.combined.FBSM$ArExp <- Idents(object = CtrlvARKO.combined.FBSM.Ar)

Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp"
tiff(file = "CtrlvARKO.combined.FBSM ArExt split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("red", "light grey"))
dev.off()


#Add EGFP info
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
CtrlvARKO.combined.FBSM.EGFPPos <- subset(x=CtrlvARKO.combined.FBSM,  subset = `EGFP` > 0)
CtrlvARKO.combined.FBSM.EGFPNeg <- subset(x=CtrlvARKO.combined.FBSM,  subset = `EGFP` == 0)
Idents(object = CtrlvARKO.combined.FBSM.EGFPPos) <- "EGFPPos"
Idents(object = CtrlvARKO.combined.FBSM.EGFPNeg) <- "EGFPNeg"
CtrlvARKO.combined.FBSM.EGFPPos[["EGFPExp"]] <- Idents(object = CtrlvARKO.combined.FBSM.EGFPPos)
CtrlvARKO.combined.FBSM.EGFPNeg[["EGFPExp"]] <- Idents(object = CtrlvARKO.combined.FBSM.EGFPNeg)
CtrlvARKO.combined.FBSM.EGFP <- merge(x = CtrlvARKO.combined.FBSM.EGFPPos, y = CtrlvARKO.combined.FBSM.EGFPNeg)
Idents(object = CtrlvARKO.combined.FBSM.EGFP) <- "EGFPExp"
CtrlvARKO.combined.FBSM$EGFPExp <- Idents(object = CtrlvARKO.combined.FBSM.EGFP)


####Subclustering EGFP+ FBSM####
Idents(object = CtrlvARKO.combined.FBSM) <- "EGFPExp"

Idents(object = CtrlvARKO.combined.FBSM) <- "EGFPExp"
tiff(file = "CtrlvARKO.combined.FBSM EGFP split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("red", "light grey"))
dev.off()

combined.FBSM.EGFP <- subset(CtrlvARKO.combined.FBSM, idents = c("EGFPPos"))

#Heatmap_FBSM
DefaultAssay(combined.FBSM.EGFP) <- "RNA"
Idents(object = combined.FBSM.EGFP) <- "stim"
DimPlot(combined.FBSM.EGFP, reduction = "umap")
combined.FBSM.EGFP <- ScaleData(combined.FBSM.EGFP, features = rownames(combined.FBSM.EGFP))
combined.FBSM.EGFP.markers <- FindAllMarkers(combined.FBSM.EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.FBSM.EGFP.markers, "combined.FBSM.EGFP.markers.csv")
combined.FBSM.EGFPTop50 <- combined.FBSM.EGFP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "combined.FBSM.EGFP Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(combined.FBSM.EGFP, features = c(combined.FBSM.EGFPTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = combined.FBSM.EGFP) <- "celltype"
combined.FB.EGFP <- subset(combined.FBSM.EGFP, idents = c("FB"))
combined.SM.EGFP <- subset(combined.FBSM.EGFP, idents = c("SM"))

#Heatmap_FB
DefaultAssay(combined.FB.EGFP) <- "RNA"
Idents(object = combined.FB.EGFP) <- "stim"
combined.FB.EGFP <- ScaleData(combined.FB.EGFP, features = rownames(combined.FB.EGFP))
combined.FB.EGFP.markers <- FindAllMarkers(combined.FB.EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.FB.EGFP.markers, "combined.FB.EGFP.markers.csv")
combined.FB.EGFPTop50 <- combined.FB.EGFP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "combined.FB.EGFP Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(combined.FB.EGFP, features = c(combined.FB.EGFPTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_SM
DefaultAssay(combined.SM.EGFP) <- "RNA"
Idents(object = combined.SM.EGFP) <- "stim"
combined.SM.EGFP <- ScaleData(combined.SM.EGFP, features = rownames(combined.SM.EGFP))
combined.SM.EGFP.markers <- FindAllMarkers(combined.SM.EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.SM.EGFP.markers, "combined.SM.EGFP.markers.csv")
combined.SM.EGFPTop50 <- combined.SM.EGFP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "combined.SM.EGFP Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(combined.SM.EGFP, features = c(combined.SM.EGFPTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violinplot
DefaultAssay(combined.FB.EGFP) <- "RNA"
tiff(file = "combined.FB.EGFP Vln Igfbp3.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB.EGFP, features = "Igfbp3", group.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$EGFPExp, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp.Strocelltype <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$Strocelltype, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp.Strocelltype.stim <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$stim, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype.stim"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp.celltype <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$celltype, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.celltype"
CtrlvARKO.combined.FBSM$ArExp.EGFPExp.celltype.stim <- paste(Idents(CtrlvARKO.combined.FBSM), CtrlvARKO.combined.FBSM$stim, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.celltype.stim"

DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap")

table(Idents(CtrlvARKO.combined.FBSM))

#GFP+Ar+CtrlvGFP+Ar-ARKO
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.celltype.stim"
combined.FB.Ar.EGFP <- subset(CtrlvARKO.combined.FBSM, idents = c("ArPos_EGFPPos_FB_Ctrl", "ArNeg_EGFPPos_FB_ARKO"))
combined.SM.Ar.EGFP <- subset(CtrlvARKO.combined.FBSM, idents = c("ArPos_EGFPPos_SM_Ctrl", "ArNeg_EGFPPos_SM_ARKO"))
combined.Stro.Ar.EGFP <- subset(CtrlvARKO.combined.FBSM, idents = c("ArPos_EGFPPos_SM_Ctrl", "ArNeg_EGFPPos_SM_ARKO", "ArPos_EGFPPos_FB_Ctrl", "ArNeg_EGFPPos_FB_ARKO"))

#Heatmap_FB
DefaultAssay(combined.FB.Ar.EGFP) <- "RNA"
Idents(object = combined.FB.Ar.EGFP) <- "ArExp.EGFPExp.celltype.stim"
combined.FB.Ar.EGFP <- ScaleData(combined.FB.Ar.EGFP, features = rownames(combined.FB.Ar.EGFP))
combined.FB.Ar.EGFP.markers <- FindAllMarkers(combined.FB.Ar.EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.FB.Ar.EGFP.markers, "combined.FB.Ar.EGFP.markers.csv")
combined.FB.Ar.EGFPTop50 <- combined.FB.Ar.EGFP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "combined.FB.Ar.EGFP Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(combined.FB.Ar.EGFP, features = c(combined.FB.Ar.EGFPTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_SM
DefaultAssay(combined.SM.Ar.EGFP) <- "RNA"
Idents(object = combined.SM.Ar.EGFP) <- "ArExp.EGFPExp.celltype.stim"
combined.SM.Ar.EGFP <- ScaleData(combined.SM.Ar.EGFP, features = rownames(combined.SM.Ar.EGFP))
combined.SM.Ar.EGFP.markers <- FindAllMarkers(combined.SM.Ar.EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.SM.Ar.EGFP.markers, "combined.SM.Ar.EGFP.markers.csv")
combined.SM.Ar.EGFPTop50 <- combined.SM.Ar.EGFP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "combined.SM.Ar.EGFP Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(combined.SM.Ar.EGFP, features = c(combined.SM.Ar.EGFPTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_SM
DefaultAssay(combined.Stro.Ar.EGFP) <- "RNA"
Idents(object = combined.Stro.Ar.EGFP) <- "stim"
combined.Stro.Ar.EGFP <- ScaleData(combined.Stro.Ar.EGFP, features = rownames(combined.Stro.Ar.EGFP))
combined.Stro.Ar.EGFP.markers <- FindAllMarkers(combined.Stro.Ar.EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.Stro.Ar.EGFP.markers, "combined.Stro.Ar.EGFP.markers.csv")
combined.Stro.Ar.EGFPTop50 <- combined.Stro.Ar.EGFP.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "combined.Stro.Ar.EGFP Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(combined.Stro.Ar.EGFP, features = c(combined.Stro.Ar.EGFPTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violinplot
DefaultAssay(combined.FB.Ar.EGFP) <- "RNA"
tiff(file = "combined.FB.Ar.EGFP Vln Igfbp3.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB.Ar.EGFP, features = c("Igfbp3"), cols = c("#3399FF", "#E06666"), split.by = "ArExp.EGFPExp.celltype.stim", pt.size = 0)
dev.off()

#FeaturePlots
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
tiff(file = "combined FBSM Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined FBSM EGFP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined FBSM Gli1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined FBSM Igfbp3 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Igfbp3"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
tiff(file = "combined FBSM Ar split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM EGFP split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Gli1 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Igfbp3 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Igfbp3"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()

tiff(file = "combined FBSM Cxcl14 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Cxcl14"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Twist1 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Twist1"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Ly6e split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Ly6e"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Sgk1 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Sgk1"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Foxf1 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Foxf1"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Aldh1a1 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Aldh1a1"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Violinplot
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM, features = "Ar", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM, features = "EGFP", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Gli1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM, features = "Gli1", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Igfbp3 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM, features = "Igfbp3", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

tiff(file = "CtrlvARKO.combined.FB Twist1 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Twist1", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Foxf1 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Foxf1", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Il11 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Il11", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Cxcl10 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Cxcl10", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Sox9 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Sox9", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Pdgfrb Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Pdgfrb", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Pdgfra Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Pdgfra", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Fap Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Fap", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB S100a4 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "S100a4", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Acta2 Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Acta2", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FB Vim Vln.tiff", width = 6, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB, features = "Vim", group.by = "Strocelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

#DEGs_celltype
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
CtrlvARKO.combined.FBSM.celltype.marker <- FindAllMarkers(CtrlvARKO.combined.FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.FBSM.celltype.marker, "CtrlvARKO.combined.FBSM.celltype.marker.csv")

#DEGs_ctrlvARKO
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))
CtrlvARKO.combined.FBSM.marker <- FindAllMarkers(CtrlvARKO.combined.FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.FBSM.marker, "CtrlvARKO.combined.FBSM.marker.csv")

#Heatmap
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
CtrlvARKO.combined.FBSM.Top50 <- CtrlvARKO.combined.FBSM.marker %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
tiff(file = "CtrlvARKO.combined.FBSM Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.FBSM, features = c(CtrlvARKO.combined.FBSM.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs_FB1-5
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FB <- subset(CtrlvARKO.combined.FBSM, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
DefaultAssay(CtrlvARKO.combined.FB) <- "RNA"
Idents(object = CtrlvARKO.combined.FB) <- "Strocelltype"
CtrlvARKO.combined.FB <- ScaleData(CtrlvARKO.combined.FB, features = rownames(CtrlvARKO.combined.FB))
CtrlvARKO.combined.FB.marker <- FindAllMarkers(CtrlvARKO.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(CtrlvARKO.combined.FB, "CtrlvARKO.combined.FB.marker.csv")

#Heatmap_FB1-5
DefaultAssay(CtrlvARKO.combined.FB) <- "RNA"
tiff(file = "CtrlvARKO.combined.FB Heatmap CAF markers-1.tiff", width = 10, height = 2, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_FB1-5_Ctrl
Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
CtrlvARKO.combined.FB <- subset(CtrlvARKO.combined.FBSM, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
Idents(object = CtrlvARKO.combined.FB) <- "stim"
Ctrl.FB <- subset(CtrlvARKO.combined.FB, idents = c("Ctrl"))
DefaultAssay(Ctrl.FB) <- "RNA"
Idents(object = Ctrl.FB) <- "Strocelltype"
Ctrl.FB <- ScaleData(Ctrl.FB, features = rownames(Ctrl.FB))
Ctrl.FB.marker <- FindAllMarkers(Ctrl.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

tiff(file = "CtrlvARKO.combined.FB stim Heatmap CAF markers.tiff", width = 10, height = 2, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "Ctrl.FB Heatmap CAF markers.tiff", width = 10, height = 2, units = "in", compression = "lzw", res = 200)
DoHeatmap(Ctrl.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = CtrlvARKO.combined.FB) <- "stim"
ARKO.FB <- subset(CtrlvARKO.combined.FB, idents = c("ARKO"))
DefaultAssay(ARKO.FB) <- "RNA"
Idents(object = ARKO.FB) <- "Strocelltype"
ARKO.FB <- ScaleData(ARKO.FB, features = rownames(ARKO.FB))
ARKO.FB.marker <- FindAllMarkers(ARKO.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

tiff(file = "ARKO.FB Heatmap CAF markers.tiff", width = 10, height = 2, units = "in", compression = "lzw", res = 200)
DoHeatmap(ARKO.FB, features = c("Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs_ArPosEGFPPos_CtrlFB5vArNegEGFPPos_ARKOFB5
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype.stim"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))

ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.IPA.Markers <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "ArPos_EGFPPos_FB5_Ctrl", ident.2 = "ArNeg_EGFPPos_FB5_ARKO", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.IPA.Markers, "ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.IPA.Markers.csv")
ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.GSEA.Markers <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "ArPos_EGFPPos_FB5_Ctrl", ident.2 = "ArNeg_EGFPPos_FB5_ARKO", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.GSEA.Markers, "ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.GSEA.Markers.csv")


#FB5 markers
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
tiff(file = "combined FBSM Nrg1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Nrg1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Wnt2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Wnt2"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Sp1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Sp1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

tiff(file = "combined FBSM Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM EGFP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Igfbp3 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Igfbp3"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Wnt2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Wnt2"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Sp1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Sp1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Fos split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Fos"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM Cdkn1a split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Cdkn1a"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM P53 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Trp53"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Coexpression plot
DefaultAssay(CtrlvARKO.combined.FBSM.Ctrl) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM Ctrl EGFP Igfbp3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.Ctrl, reduction = "umap", features = c("EGFP", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Ctrl Ar Igfbp3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.Ctrl, reduction = "umap", features = c("Ar", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
DefaultAssay(CtrlvARKO.combined.FBSM.ARKO) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM ARKO EGFP Igfbp3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.ARKO, reduction = "umap", features = c("EGFP", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM ARKO Ar Igfbp3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.ARKO, reduction = "umap", features = c("Ar", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

#FB5 CAF markers
tiff(file = "combined FBSM Ogn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Ogn"), cols = c("light grey", "blue"), pt.size = 1, max.cutoff = "q90")
dev.off()

#FB Heatmap
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype.stim"
CtrlvARKO.combined.FB1 <- subset(CtrlvARKO.combined.FBSM, idents = c("ArPos_EGFPPos_FB1_Ctrl", "ArNeg_EGFPPos_FB1_ARKO"))
DefaultAssay(CtrlvARKO.combined.FB1) <- "RNA"
CtrlvARKO.combined.FB1 <- ScaleData(CtrlvARKO.combined.FB1, features = rownames(CtrlvARKO.combined.FB1))

CtrlvARKO.combined.FB1.Markers <- FindMarkers(CtrlvARKO.combined.FB1, ident.1 = "ArPos_EGFPPos_FB1_Ctrl", ident.2 = "ArNeg_EGFPPos_FB1_ARKO", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.FB1.Markers, "ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB1.IPA.Markers.csv")


tiff(file = "CtrlvARKO.combined.FB1 Heatmap-1.tiff", width = 4, height = 9, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.FB1, features = c("Ar", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "CtrlvARKO.combined.FB1 Heatmap-2.tiff", width = 4, height = 9, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.FB1, features = c("Ar", "Fkbp5", "Aldh1a1", "Ly6a", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

tiff(file = "CtrlvARKO.combined.FB1 Heatmap-3.tiff", width = 4, height = 9, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.FB1, features = c("Ar", "Fkbp5", "Aldh1a1", "Ly6a", "Fos", "Fosb", "Jun", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violinplot
DefaultAssay(CtrlvARKO.combined.FB5) <- "RNA"
tiff(file = "CtrlvARKO.combined.FB5 Vln all.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FB5, features = c("Ar", "Fkbp5", "Sgk1", "Aldh1a1", "Ly6a",  "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Gadd45a", "Fas", "Bcl2", "Igfbp3", "Igfbp7"), split.by = "ArExp.EGFPExp.Strocelltype.stim", pt.size = 0)
dev.off()


####comparison with Adam's data####
#Violinplot
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
DimPlot(CtrlvARKO.combined.FBSM, reduction = "umap")

Idents(object = CtrlvARKO.combined.FBSM) <- "celltype"
CtrlvARKO.combined.FBSM1 <- subset(CtrlvARKO.combined.FBSM, idents = c("FB", "SM"))


DefaultAssay(CtrlvARKO.combined.FBSM1) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM1 Igfbp3 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Igfbp3", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Dcn celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Dcn", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Mknk2 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Mknk2", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Cdkn1a celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Cdkn1a", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Lum celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Lum", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Apod celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Apod", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Filip1l celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Filip1l", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 Cystm1 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Cystm1", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

tiff(file = "CtrlvARKO.combined.FBSM1 Gli1 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "Gli1", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM1 mGFP celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM1, features = "EGFP", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

DefaultAssay(combined.FBSM.EGFP) <- "RNA"
Idents(object = combined.FBSM.EGFP) <- "celltype"
combined.FBSM.EGFP1 <- subset(combined.FBSM.EGFP, idents = c("FB", "SM"))

DefaultAssay(combined.FBSM.EGFP1) <- "RNA"
tiff(file = "combined.FBSM.EGFP1 Igfbp3 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Igfbp3", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Dcn celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Dcn", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Mknk2 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Mknk2", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Cdkn1a celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Cdkn1a", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Lum celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Lum", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Apod celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Apod", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Filip1l celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Filip1l", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 Cystm1 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Cystm1", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

tiff(file = "combined.FBSM.EGFP1 Gli1 celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "Gli1", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "combined.FBSM.EGFP1 mGFP celltype Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.EGFP1, features = "EGFP", group.by = "celltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()


####Pseudotime Stro####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/CtrlvARKO.combined.FBSM/Pseudotime")

Idents(object = CtrlvARKO.combined.FBSM) <- "Strocelltype"
combined.FB <- subset(CtrlvARKO.combined.FBSM, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
Idents(object = combined.FB) <- "seurat_clusters"


DefaultAssay(combined.FB) <- "RNA"
tiff(file = "combined.FB Thy1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB, reduction = "umap", features = c("Thy1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = combined.FB) <- "Strocelltype"

combined.FB <- RenameIdents(object = combined.FB, 'FB5' = "FB1", 'FB2' = "FB2",
                                        'FB3' = "FB3", 'FB1' = "FB4", 'FB4' = "FB5")
combined.FB[["Strocelltype"]] <- Idents(object = combined.FB)

Idents(object = combined.FB) <- "Strocelltype"
DimPlot(combined.FB, reduction = "umap")


#Combined
DefaultAssay(combined.FB) <- "RNA"
FBPseudo <- as.CellDataSet(combined.FB)
FBPseudo <- detectGenes(FBPseudo, min_expr = 0.1)
print(head(fData(FBPseudo)))

expressed_genes <- row.names(subset(fData(FBPseudo),
                                    num_cells_expressed >= 10))

pData(FBPseudo)$Total_mRNAs <- Matrix::colSums(exprs(FBPseudo))
FBPseudo <- FBPseudo[,pData(FBPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(FBPseudo), geom =
        "density")

FBPseudo <- estimateSizeFactors(FBPseudo)
FBPseudo <- estimateDispersions(FBPseudo)

disp_table <- dispersionTable(FBPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
FBPseudo <- setOrderingFilter(FBPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(FBPseudo)

#FBPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(FBPseudo, return_all = F) # norm_method='log'

FBPseudo <- reduceDimension(FBPseudo, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
FBPseudo <- clusterCells(FBPseudo, num_clusters = 2)

plot_cell_clusters(FBPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(FBPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
FBPseudo <- setOrderingFilter(FBPseudo, ordering_genes)
plot_ordering_genes(FBPseudo)

FBPseudo <- reduceDimension(FBPseudo, max_components = 2,
                             method = 'DDRTree')

FBPseudo <- orderCells(FBPseudo)

GM_state <- function(FBPseudo){
  if (length(unique(pData(FBPseudo)$State)) > 1){
    T0_counts <- table(pData(FBPseudo)$State, pData(FBPseudo)$seurat_clusters)[,"3"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

FBPseudo <- orderCells(FBPseudo, root_state = GM_state(FBPseudo))

#Visualization

tiff(file = "FBPseudo Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "FBPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

tiff(file = "FBPseudo Strocelltype.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("red", "purple", "blue", "green", "orange"))
dev.off()

tiff(file = "FBPseudo Strocelltype stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("red", "purple", "blue", "green", "orange")) + facet_wrap(~stim, nrow = 1) 
dev.off()

tiff(file = "FBPseudo Strocelltype split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("red", "purple", "blue", "green", "orange")) + facet_wrap(~Strocelltype, nrow = 2) 
dev.off()

tiff(file = "FBPseudo Ar.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "FBPseudo EGFP.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "FBPseudo Igfbp3.tiff", width = 8, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo, markers = "Igfbp3", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) + facet_wrap(~stim, nrow = 1) 
dev.off()


tiff(file = "EpiPseudo Ar MYCtg Krt5 Krt19 Pbsn Pcna.tiff", width = 12, height = 8.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Ar", "MYC-transgene", "Krt5", "Krt19", "Pbsn", "Pcna"), use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#Ctrl
Idents(object = combined.FB) <- "stim"
combined.FB.Ctrl <- subset(combined.FB, idents = c("Ctrl"))

DefaultAssay(combined.FB.Ctrl) <- "RNA"
FBPseudo.Ctrl <- as.CellDataSet(combined.FB.Ctrl)
FBPseudo.Ctrl <- detectGenes(FBPseudo.Ctrl, min_expr = 0.1)
print(head(fData(FBPseudo.Ctrl)))

expressed_genes <- row.names(subset(fData(FBPseudo.Ctrl),
                                    num_cells_expressed >= 10))

pData(FBPseudo.Ctrl)$Total_mRNAs <- Matrix::colSums(exprs(FBPseudo.Ctrl))
FBPseudo.Ctrl <- FBPseudo.Ctrl[,pData(FBPseudo.Ctrl)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(FBPseudo.Ctrl), geom =
        "density")

FBPseudo.Ctrl <- estimateSizeFactors(FBPseudo.Ctrl)
FBPseudo.Ctrl <- estimateDispersions(FBPseudo.Ctrl)

disp_table <- dispersionTable(FBPseudo.Ctrl)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
FBPseudo.Ctrl <- setOrderingFilter(FBPseudo.Ctrl, unsup_clustering_genes$gene_id)
plot_ordering_genes(FBPseudo.Ctrl)

#FBPseudo.Ctrl@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(FBPseudo.Ctrl, return_all = F) # norm_method='log'

FBPseudo.Ctrl <- reduceDimension(FBPseudo.Ctrl, max_components = 2, num_dim = 20,
                            reduction_method = 'tSNE', verbose = T)
FBPseudo.Ctrl <- clusterCells(FBPseudo.Ctrl, num_clusters = 2)

plot_cell_clusters(FBPseudo.Ctrl, color_by = "Strocelltype")

diff_test_res <- differentialGeneTest(FBPseudo.Ctrl[expressed_genes,],
                                      fullModelFormulaStr = "~Strocelltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
FBPseudo.Ctrl <- setOrderingFilter(FBPseudo.Ctrl, ordering_genes)
plot_ordering_genes(FBPseudo.Ctrl)

FBPseudo.Ctrl <- reduceDimension(FBPseudo.Ctrl, max_components = 2,
                            method = 'DDRTree')

FBPseudo.Ctrl <- orderCells(FBPseudo.Ctrl)

GM_state <- function(FBPseudo.Ctrl){
  if (length(unique(pData(FBPseudo.Ctrl)$State)) > 1){
    T0_counts <- table(pData(FBPseudo.Ctrl)$State, pData(FBPseudo.Ctrl)$State)[,"9"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

FBPseudo.Ctrl <- orderCells(FBPseudo.Ctrl, root_state = GM_state(FBPseudo.Ctrl))

plot_cell_trajectory(FBPseudo.Ctrl, color_by = "State", show_branch_points = FALSE)

#Visualization

tiff(file = "FBPseudo.Ctrl Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

tiff(file = "FBPseudo.Ctrl Strocelltype.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("orange", "purple", "blue", "red", "green"))
dev.off()

tiff(file = "FBPseudo.Ctrl Strocelltype split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("orange", "purple", "blue", "red", "green")) + facet_wrap(~Strocelltype, nrow = 2) 
dev.off()

tiff(file = "FBPseudo.Ctrl Ar.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl EGFP.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl Igfbp3.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Igfbp3", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl Thy1.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Thy1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#
ordering_genes <- row.names (subset(diff_test_res, qval < 0.05))
FBPseudo.Ctrl <- setOrderingFilter(FBPseudo.Ctrl, ordering_genes)
plot_ordering_genes(FBPseudo.Ctrl)

FBPseudo.Ctrl <- reduceDimension(FBPseudo.Ctrl, max_components = 2,
                                 method = 'DDRTree')

FBPseudo.Ctrl <- orderCells(FBPseudo.Ctrl)

GM_state <- function(FBPseudo.Ctrl){
  if (length(unique(pData(FBPseudo.Ctrl)$State)) > 1){
    T0_counts <- table(pData(FBPseudo.Ctrl)$State, pData(FBPseudo.Ctrl)$State)[,"5"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

FBPseudo.Ctrl <- orderCells(FBPseudo.Ctrl, root_state = GM_state(FBPseudo.Ctrl))

plot_cell_trajectory(FBPseudo.Ctrl, color_by = "State", show_branch_points = FALSE)

#Visualization

tiff(file = "FBPseudo.Ctrl Pseudotime 0.05.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

tiff(file = "FBPseudo.Ctrl Strocelltype 0.05.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("orange", "purple", "blue", "red", "green"))
dev.off()

tiff(file = "FBPseudo.Ctrl Strocelltype split 0.05.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("orange", "purple", "blue", "red", "green")) + facet_wrap(~Strocelltype, nrow = 2) 
dev.off()

tiff(file = "FBPseudo.Ctrl Ar 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl EGFP 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl Igfbp3 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Igfbp3", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl Thy1 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Thy1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "FBPseudo.Ctrl Pdgfrb 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Pdgfrb", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl Il11 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Il11", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.Ctrl Twist1 0.05.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.Ctrl, markers = "Twist1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#ARKO
Idents(object = combined.FB) <- "stim"
combined.FB.ARKO <- subset(combined.FB, idents = c("ARKO"))

DefaultAssay(combined.FB.ARKO) <- "RNA"
FBPseudo.ARKO <- as.CellDataSet(combined.FB.ARKO)
FBPseudo.ARKO <- detectGenes(FBPseudo.ARKO, min_expr = 0.1)
print(head(fData(FBPseudo.ARKO)))

expressed_genes <- row.names(subset(fData(FBPseudo.ARKO),
                                    num_cells_expressed >= 10))

pData(FBPseudo.ARKO)$Total_mRNAs <- Matrix::colSums(exprs(FBPseudo.ARKO))
FBPseudo.ARKO <- FBPseudo.ARKO[,pData(FBPseudo.ARKO)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(FBPseudo.ARKO), geom =
        "density")

FBPseudo.ARKO <- estimateSizeFactors(FBPseudo.ARKO)
FBPseudo.ARKO <- estimateDispersions(FBPseudo.ARKO)

disp_table <- dispersionTable(FBPseudo.ARKO)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
FBPseudo.ARKO <- setOrderingFilter(FBPseudo.ARKO, unsup_clustering_genes$gene_id)
plot_ordering_genes(FBPseudo.ARKO)

#FBPseudo.ARKO@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(FBPseudo.ARKO, return_all = F) # norm_method='log'

FBPseudo.ARKO <- reduceDimension(FBPseudo.ARKO, max_components = 2, num_dim = 20,
                                 reduction_method = 'tSNE', verbose = T)
FBPseudo.ARKO <- clusterCells(FBPseudo.ARKO, num_clusters = 2)

plot_cell_clusters(FBPseudo.ARKO, color_by = "Strocelltype")

diff_test_res <- differentialGeneTest(FBPseudo.ARKO[expressed_genes,],
                                      fullModelFormulaStr = "~Strocelltype")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
FBPseudo.ARKO <- setOrderingFilter(FBPseudo.ARKO, ordering_genes)
plot_ordering_genes(FBPseudo.ARKO)

FBPseudo.ARKO <- reduceDimension(FBPseudo.ARKO, max_components = 2,
                                 method = 'DDRTree')

FBPseudo.ARKO <- orderCells(FBPseudo.ARKO)

GM_state <- function(FBPseudo.ARKO){
  if (length(unique(pData(FBPseudo.ARKO)$State)) > 1){
    T0_counts <- table(pData(FBPseudo.ARKO)$State, pData(FBPseudo.ARKO)$Strocelltype)[,"FB1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

FBPseudo.ARKO <- orderCells(FBPseudo.ARKO, root_state = GM_state(FBPseudo.ARKO))

#Visualization

tiff(file = "FBPseudo.ARKO Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

tiff(file = "FBPseudo.ARKO Strocelltype.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("orange", "purple", "blue", "red", "green"))
dev.off()

tiff(file = "FBPseudo.ARKO Strocelltype split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, color_by = "Strocelltype", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("orange", "purple", "blue", "red", "green")) + facet_wrap(~Strocelltype, nrow = 2) 
dev.off()

tiff(file = "FBPseudo.ARKO Ar.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3")) 
dev.off()
tiff(file = "FBPseudo.ARKO EGFP.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.ARKO Igfbp3.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, markers = "Igfbp3", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "FBPseudo.ARKO Pdgfrb.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, markers = "Pdgfrb", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.ARKO Il11.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, markers = "Il11", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "FBPseudo.ARKO Twist1.tiff", width = 4, height = 4.5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(FBPseudo.ARKO, markers = "Twist1", use_color_gradient = TRUE, show_branch_points = FALSE, cell_size = 1) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

####Subclustering FBSM-1####
Idents(object = CtrlvARKO.combined.Stro) <- "Strocelltype"

CtrlvARKO.combined.FBSM2 <- subset(CtrlvARKO.combined.Stro, idents = c("FB1", "FB2", "FB3", "FB4", "FB5", "SM"))

DefaultAssay(CtrlvARKO.combined.FBSM2) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.FBSM2 <- ScaleData(CtrlvARKO.combined.FBSM2, verbose = FALSE)
CtrlvARKO.combined.FBSM2 <- RunPCA(CtrlvARKO.combined.FBSM2, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.FBSM2, ndims = 50)
#Umap and Clustering
CtrlvARKO.combined.FBSM2 <- FindNeighbors(CtrlvARKO.combined.FBSM2, reduction = "pca", dims = 1:16)
CtrlvARKO.combined.FBSM2 <- FindClusters(CtrlvARKO.combined.FBSM2, resolution = 0.5)
CtrlvARKO.combined.FBSM2 <- RunTSNE(CtrlvARKO.combined.FBSM2, reduction = "pca", dims = 1:16)
CtrlvARKO.combined.FBSM2 <- RunUMAP(CtrlvARKO.combined.FBSM2, reduction = "pca", dims = 1:16)
DimPlot(CtrlvARKO.combined.FBSM2, reduction = "umap", label = TRUE)
DimPlot(CtrlvARKO.combined.FBSM2, split.by = "stim", reduction = "umap", label = TRUE)

#Rename
Idents(object = CtrlvARKO.combined.FBSM2) <- "seurat_clusters"
CtrlvARKO.combined.FBSM2 <- RenameIdents(object = CtrlvARKO.combined.FBSM2, '3' = "FB1", '0' = "FB2",
                                        '2' = "FB3", '8' = "FB3", '4' = "FB4", '1' = "FB5", '5' = "FB6", 
                                        '7' = "FB7", '6' = "SM")
CtrlvARKO.combined.FBSM2[["FBSMcelltype"]] <- Idents(object = CtrlvARKO.combined.FBSM2)


Idents(object = CtrlvARKO.combined.FBSM2) <- "FBSMcelltype"
tiff(file = "FBSM2 UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.FBSM2, reduction = "umap", pt.size = 0.5)
dev.off()
tiff(file = "FBSM2 split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", pt.size = 0.5)
dev.off()


#Add Ar info
DefaultAssay(CtrlvARKO.combined.FBSM2) <- "RNA"
CtrlvARKO.combined.FBSM2.ArPos <- subset(x=CtrlvARKO.combined.FBSM2,  subset = `Ar` > 0)
CtrlvARKO.combined.FBSM2.ArNeg <- subset(x=CtrlvARKO.combined.FBSM2,  subset = `Ar` == 0)
Idents(object = CtrlvARKO.combined.FBSM2.ArPos) <- "ArPos"
Idents(object = CtrlvARKO.combined.FBSM2.ArNeg) <- "ArNeg"
CtrlvARKO.combined.FBSM2.ArPos[["ArExp"]] <- Idents(object = CtrlvARKO.combined.FBSM2.ArPos)
CtrlvARKO.combined.FBSM2.ArNeg[["ArExp"]] <- Idents(object = CtrlvARKO.combined.FBSM2.ArNeg)
CtrlvARKO.combined.FBSM2.Ar <- merge(x = CtrlvARKO.combined.FBSM2.ArPos, y = CtrlvARKO.combined.FBSM2.ArNeg)
Idents(object = CtrlvARKO.combined.FBSM2.Ar) <- "ArExp"
CtrlvARKO.combined.FBSM2$ArExp <- Idents(object = CtrlvARKO.combined.FBSM2.Ar)

#Add EGFP info
DefaultAssay(CtrlvARKO.combined.FBSM2) <- "RNA"
CtrlvARKO.combined.FBSM2.EGFPPos <- subset(x=CtrlvARKO.combined.FBSM2,  subset = `EGFP` > 0)
CtrlvARKO.combined.FBSM2.EGFPNeg <- subset(x=CtrlvARKO.combined.FBSM2,  subset = `EGFP` == 0)
Idents(object = CtrlvARKO.combined.FBSM2.EGFPPos) <- "EGFPPos"
Idents(object = CtrlvARKO.combined.FBSM2.EGFPNeg) <- "EGFPNeg"
CtrlvARKO.combined.FBSM2.EGFPPos[["EGFPExp"]] <- Idents(object = CtrlvARKO.combined.FBSM2.EGFPPos)
CtrlvARKO.combined.FBSM2.EGFPNeg[["EGFPExp"]] <- Idents(object = CtrlvARKO.combined.FBSM2.EGFPNeg)
CtrlvARKO.combined.FBSM2.EGFP <- merge(x = CtrlvARKO.combined.FBSM2.EGFPPos, y = CtrlvARKO.combined.FBSM2.EGFPNeg)
Idents(object = CtrlvARKO.combined.FBSM2.EGFP) <- "EGFPExp"
CtrlvARKO.combined.FBSM2$EGFPExp <- Idents(object = CtrlvARKO.combined.FBSM2.EGFP)

#Cell counts
Idents(object = CtrlvARKO.combined.FBSM2) <- "ArExp"
CtrlvARKO.combined.FBSM2$ArExp.EGFPExp <- paste(Idents(CtrlvARKO.combined.FBSM2), CtrlvARKO.combined.FBSM2$EGFPExp, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM2) <- "ArExp.EGFPExp"
CtrlvARKO.combined.FBSM2$ArExp.EGFPExp.FBSMcelltype <- paste(Idents(CtrlvARKO.combined.FBSM2), CtrlvARKO.combined.FBSM2$FBSMcelltype, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM2) <- "ArExp.EGFPExp.FBSMcelltype"
CtrlvARKO.combined.FBSM2$ArExp.EGFPExp.FBSMcelltype.stim <- paste(Idents(CtrlvARKO.combined.FBSM2), CtrlvARKO.combined.FBSM2$stim, sep = "_")
Idents(object = CtrlvARKO.combined.FBSM2) <- "ArExp.EGFPExp.FBSMcelltype.stim"

DimPlot(CtrlvARKO.combined.FBSM2, reduction = "umap")

table(Idents(CtrlvARKO.combined.FBSM2))

#FeaturePlots
DefaultAssay(CtrlvARKO.combined.FBSM2) <- "RNA"
tiff(file = "combined FBSM2 Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined FBSM2 EGFP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined FBSM2 Gli1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(CtrlvARKO.combined.FBSM2) <- "RNA"
tiff(file = "combined FBSM2 Ar split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM2 EGFP split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "combined FBSM2 Gli1 split-1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM2, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 1, max.cutoff = "q90")
dev.off()


#Violinplot
DefaultAssay(CtrlvARKO.combined.FBSM2) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM2 Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM2, features = "Ar", group.by = "FBSMcelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM2 EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM2, features = "EGFP", group.by = "FBSMcelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()


#DEGs_ArPosEGFPPos_CtrlFB5vArNegEGFPPos_ARKOFB5
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
Idents(object = CtrlvARKO.combined.FBSM) <- "ArExp.EGFPExp.Strocelltype.stim"
CtrlvARKO.combined.FBSM <- ScaleData(CtrlvARKO.combined.FBSM, features = rownames(CtrlvARKO.combined.FBSM))

ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.IPA.Markers <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "ArPos_EGFPPos_FB5_Ctrl", ident.2 = "ArNeg_EGFPPos_FB5_ARKO", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.IPA.Markers, "ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.IPA.Markers.csv")
ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.GSEA.Markers <- FindMarkers(CtrlvARKO.combined.FBSM, ident.1 = "ArPos_EGFPPos_FB5_Ctrl", ident.2 = "ArNeg_EGFPPos_FB5_ARKO", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.GSEA.Markers, "ArPosEGFPPosCtrlvArNegEGFPPosARKO.FB5.GSEA.Markers.csv")

####Subclustering Ctrl & ARKO Stro####
Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
CtrlvARKO.combined.FBSM.Ctrl <- subset(CtrlvARKO.combined.FBSM, idents = c("Ctrl"))

DefaultAssay(CtrlvARKO.combined.FBSM.Ctrl) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM.Ctrl EGFP Ar Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.Ctrl, reduction = "umap", features = c("EGFP", "Ar"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM.Ctrl EGFP Igfbp3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.Ctrl, reduction = "umap", features = c("EGFP", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
CtrlvARKO.combined.FBSM.ARKO <- subset(CtrlvARKO.combined.FBSM, idents = c("ARKO"))

DefaultAssay(CtrlvARKO.combined.FBSM.ARKO) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM.ARKO EGFP Ar Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.ARKO, reduction = "umap", features = c("EGFP", "Ar"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM.ARKO EGFP Igfbp3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM.ARKO, reduction = "umap", features = c("EGFP", "Igfbp3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

#### Running SingleCellSignalR  ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/Ligand_Receptor_Interaction")

Idents(object = CtrlvARKO.combined.MYCBE) <- "stim"
Ctrlbasal <- subset(CtrlvARKO.combined.MYCBE, idents = c("Ctrl"))
Ctrlbasal <- RenameIdents(object = Ctrlbasal, 'Ctrl' = "Ctrlbasal")
ARKObasal <- subset(CtrlvARKO.combined.MYCBE, idents = c("ARKO"))
ARKObasal <- RenameIdents(object = ARKObasal, 'ARKO' = "ARKObasal")

Idents(object = CtrlvARKO.combined.FB5) <- "stim"
Ctrlfb <- subset(CtrlvARKO.combined.FB5, idents = c("Ctrl"))
Ctrlfb <- RenameIdents(object = Ctrlfb, 'Ctrl' = "Ctrlfb")
ARKOfb <- subset(CtrlvARKO.combined.FB5, idents = c("ARKO"))
ARKOfb <- RenameIdents(object = ARKOfb, 'ARKO' = "ARKOfb")

FBandBasalCtrl <- merge(x = Ctrlfb, y = Ctrlbasal)
FBandBasalARKO <- merge(x = ARKOfb, y = ARKObasal)

InteractomeObject <- merge(FBandBasalCtrl, FBandBasalARKO)
InteractomeObject[["CtrlvARKO"]] <- Idents(object = InteractomeObject)
Idents(object = InteractomeObject) <- "CtrlvARKO"
table(Idents(InteractomeObject))

#
data = data.frame(InteractomeObject[["RNA"]]@data)
write.csv(data, file = "InteractomeObjectGenes.csv")

cluster = as.numeric(Idents(InteractomeObject))
genes <- rownames(data)
cluster_analysis(data = data, genes = genes, cluster = cluster)
signal = cell_signaling(data = data, genes = genes, cluster = cluster, s.score = 0.1, species = "mus musculus", tol = 1)

# to select specifc genes
#import excel file with LRint column as numeric
Ctrl_FBtoBE <- as.data.frame(read.csv("Ctrl_FBtoBE.csv"))
signal$int13 <- Ctrl_FBtoBE
visualize_interactions(signal, show.in = 13)

tiff(file = "FBtoBE_Ctrl_LigandR_Selected-1.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
visualize_interactions(signal, show.in = 13)
dev.off()

tiff(file = "FBtoBE_Ctrl_LigandR_Selected-2.tiff", width = 5.5, height = 4.5, units = "in", compression = "lzw", res = 800)
visualize_interactions(signal, show.in = 13)
dev.off()

ARKO_FBtoBE <- as.data.frame(read.csv("ARKO_FBtoBE-1.csv"))
signal$int14 <- ARKO_FBtoBE
visualize_interactions(signal, show.in = 14)

tiff(file = "FBtoBE_ARKO_LigandR_Selected-1.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
visualize_interactions(signal, show.in = 14, limit = 25)
dev.off()

tiff(file = "FBtoBE_ARKO_LigandR_Selected-2.tiff", width = 5.5, height = 4.5, units = "in", compression = "lzw", res = 800)
visualize_interactions(signal, show.in = 14, limit = 25)
dev.off()

####GEO####
CtrlvARKO.combinedExpdata = data.frame(CtrlvARKO.combined[["RNA"]]@data)
write.csv(CtrlvARKO.combinedExpdata , file = "CtrlvARKO.combinedExpdata.csv")


####Comparison with WT HiMYC HiMYC-ARKO####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/NCB revision")

WTunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191214_32848_33580_test/counts/33580_WT/outs/filtered_feature_bc_matrix")
WTunfiltered <- CreateSeuratObject(counts = WTunfiltered.data,  min.cells = 3, min.features = 500, project = "WTunfiltered")
WTunfiltered <- NormalizeData(WTunfiltered)

#Initial processing & filtering
WTunfiltered[["percent.mt"]] <- PercentageFeatureSet(WTunfiltered, pattern = "^mt-")

tiff(file = "WT Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WTunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "darkgreen", xlab = "nFeature_RNA", main = "WT Pre-filteration")
dev.off()
tiff(file = "WT Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$percent.mt, breaks = 100, col = "darkgreen", xlab = "percent.mt", main = "WT Pre-filteration")
dev.off()

WT <- subset(WTunfiltered, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 15)
tiff(file = "WT Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$nFeature_RNA, breaks = 100, col = "darkgreen", xlab = "nFeature_RNA", main = "WT Post-filteration")
dev.off()
tiff(file = "WT Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$percent.mt, breaks = 100, col = "darkgreen", xlab = "percent.mt", main = "WT Post-filteration")
dev.off()

WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 5000)
tiff(file = "WT Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT)
dev.off()

table(Idents(WTunfiltered))
table(Idents(WT))

#Run the standard workflow for visualization and clustering Ctrl2
WT <- ScaleData(WT, verbose = FALSE)
WT <- RunPCA(WT, npcs = 30, verbose = FALSE)
ElbowPlot(WT, ndims = 50)
# t-SNE and Clustering
WT <- FindNeighbors(WT, reduction = "pca", dims = 1:20)
WT <- FindClusters(WT, resolution = 0.5)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:20)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:20)
DimPlot(WT, reduction = "umap", pt.size = 0.3) 

Idents(object = WT) <- "stim"

tiff(file = "WT UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, cols = c("darkred")) 
dev.off()

#### Merging Datasets #### ####
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)
WT[["orig.clusters"]] <- Idents(object = WT)

Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = MycARKO) <- "seurat_clusters"
Idents(object = WT) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
MycARKO$stim <- "ARKO"
WT$stim <- "WT"

CtrlvARKOvWT.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, ARKO, WT), dims = 1:20)
CtrlvARKOvWT.combined <- IntegrateData(anchorset = CtrlvARKOvWT.anchors, dims = 1:20)
DefaultAssay(CtrlvARKOvWT.combined) <- "integrated"

#Run the standard workflow for visualization and clustering Ctrl2vARKO1.combined
CtrlvARKOvWT.combined <- ScaleData(CtrlvARKOvWT.combined, verbose = FALSE)
CtrlvARKOvWT.combined <- RunPCA(CtrlvARKOvWT.combined, npcs = 30, verbose = FALSE)
CtrlvARKOvWT.combined <- FindNeighbors(CtrlvARKOvWT.combined, reduction = "pca", dims = 1:20)
CtrlvARKOvWT.combined <- FindClusters(CtrlvARKOvWT.combined, resolution = 0.5)
CtrlvARKOvWT.combined <- RunUMAP(CtrlvARKOvWT.combined, reduction = "pca", dims = 1:20)
CtrlvARKOvWT.combined <- RunTSNE(CtrlvARKOvWT.combined, reduction = "pca", dims = 1:20)

Idents(object = CtrlvARKOvWT.combined) <- "stim"
DimPlot(CtrlvARKOvWT.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 

#celltype identification
Idents(object = CtrlvARKOvWT.combined) <- "seurat_clusters"
DimPlot(CtrlvARKOvWT.combined, reduction = "umap", pt.size = 0.3, label = TRUE) 

#Featureplots
DefaultAssay(CtrlvARKOvWT.combined) <- "RNA"
tiff(file = "CtrlvARKOvWT.combined Epcam.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Vim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Fbln1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Acta2.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Acta2"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Cdh5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Rgs5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined Plp1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined, reduction = "umap", features = c("Plp1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Epi reclustering
Idents(object = CtrlvARKOvWT.combined) <- "seurat_clusters"
CtrlvARKOvWT.combined.Epi <- subset(CtrlvARKOvWT.combined, idents = c("1", "2", "20", "3", "17", "7", "13", "5", "9"))
DefaultAssay(CtrlvARKOvWT.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKOvWT.combined.Epi <- ScaleData(CtrlvARKOvWT.combined.Epi, verbose = FALSE)
CtrlvARKOvWT.combined.Epi <- RunPCA(CtrlvARKOvWT.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKOvWT.combined.Epi, ndims = 50)
#Umap and Clustering
CtrlvARKOvWT.combined.Epi <- FindNeighbors(CtrlvARKOvWT.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKOvWT.combined.Epi <- FindClusters(CtrlvARKOvWT.combined.Epi, resolution = 0.5)
CtrlvARKOvWT.combined.Epi <- RunTSNE(CtrlvARKOvWT.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKOvWT.combined.Epi <- RunUMAP(CtrlvARKOvWT.combined.Epi, reduction = "pca", dims = 1:16)
DimPlot(CtrlvARKOvWT.combined.Epi, reduction = "umap", label = TRUE)

DefaultAssay(CtrlvARKOvWT.combined.Epi) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvARKOvWT.combined.Epi <- ScaleData(CtrlvARKOvWT.combined.Epi, verbose = FALSE)
CtrlvARKOvWT.combined.Epi <- RunPCA(CtrlvARKOvWT.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKOvWT.combined.Epi, ndims = 50)
#Umap and Clustering
CtrlvARKOvWT.combined.Epi <- FindNeighbors(CtrlvARKOvWT.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKOvWT.combined.Epi <- FindClusters(CtrlvARKOvWT.combined.Epi, resolution = 0.5)
CtrlvARKOvWT.combined.Epi <- RunTSNE(CtrlvARKOvWT.combined.Epi, reduction = "pca", dims = 1:16)
CtrlvARKOvWT.combined.Epi <- RunUMAP(CtrlvARKOvWT.combined.Epi, reduction = "pca", dims = 1:16)
DimPlot(CtrlvARKOvWT.combined.Epi, reduction = "umap", label = TRUE)

CtrlvARKOvWT.combined.Epi1 <- subset(CtrlvARKOvWT.combined.Epi, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13"))

CtrlvARKOvWT.combined.Epi1 <- RenameIdents(object = CtrlvARKOvWT.combined.Epi1, '2' = "BE1", '0' = "BE2",
                                   '1' = "BE3", '10' = "BE3", '6' = "LE1", '7' = "LE2", 
                                   '3' = "LE3", '13' = "LE3", '5' = "LE4", '8' = "LE4",
                                   '9' = "LE5", '4' = "LE6", '12' = "ProE", '11' = "OE")
CtrlvARKOvWT.combined.Epi1[["EpiCellType"]] <- Idents(object = CtrlvARKOvWT.combined.Epi1)
DimPlot(CtrlvARKOvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Black", "Light grey")) 

Idents(object = CtrlvARKOvWT.combined.Epi1) <- "stim"
CtrlvWT.combined.Epi1 <- subset(CtrlvARKOvWT.combined.Epi1, idents = c("Ctrl", "WT"))

Idents(object = CtrlvWT.combined.Epi1) <- "stim"
tiff(file = "CtrlvWT.combined.Epi1 UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown")) 
dev.off()

Idents(object = CtrlvWT.combined.Epi1) <- "EpiCellType"

tiff(file = "CtrlvWT.combined.Epi1 EpiCellType UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Black", "Light grey")) 
dev.off()

tiff(file = "CtrlvWT.combined.Epi1 split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvWT.combined.Epi1, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Black", "Light grey")) 
dev.off()

Idents(object = CtrlvWT.combined.Epi1) <- "stim"
CtrlvWT.combined.Epi1_Ctrl <- subset(CtrlvWT.combined.Epi1, idents = c("Ctrl"))
CtrlvWT.combined.Epi1_WT <- subset(CtrlvWT.combined.Epi1, idents = c("WT"))

DefaultAssay(CtrlvWT.combined.Epi1_Ctrl) <- "RNA"
tiff(file = "CtrlvWT.combined.Epi1_Ctrl hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_Ctrl, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvWT.combined.Epi1_Ctrl Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_Ctrl, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()
tiff(file = "CtrlvWT.combined.Epi1_Ctrl Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_Ctrl, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvWT.combined.Epi1_Ctrl Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_Ctrl, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(CtrlvWT.combined.Epi1_WT) <- "RNA"
tiff(file = "CtrlvWT.combined.Epi1_WT hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_WT, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "light grey"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvWT.combined.Epi1_WT Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_WT, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()
tiff(file = "CtrlvWT.combined.Epi1_WT Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_WT, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvWT.combined.Epi1_WT Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined.Epi1_WT, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#BE reclustering
Idents(object = CtrlvARKOvWT.combined.Epi) <- "seurat_clusters"
CtrlvARKOvWT.combined.BE <- subset(CtrlvARKOvWT.combined.Epi, idents = c("0", "1", "2", "11"))
DefaultAssay(CtrlvARKOvWT.combined.BE) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKOvWT.combined.BE <- ScaleData(CtrlvARKOvWT.combined.BE, verbose = FALSE)
CtrlvARKOvWT.combined.BE <- RunPCA(CtrlvARKOvWT.combined.BE, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKOvWT.combined.BE, ndims = 50)
#Umap and Clustering
CtrlvARKOvWT.combined.BE <- FindNeighbors(CtrlvARKOvWT.combined.BE, reduction = "pca", dims = 1:16)
CtrlvARKOvWT.combined.BE <- FindClusters(CtrlvARKOvWT.combined.BE, resolution = 0.3)
CtrlvARKOvWT.combined.BE <- RunTSNE(CtrlvARKOvWT.combined.BE, reduction = "pca", dims = 1:16)
CtrlvARKOvWT.combined.BE <- RunUMAP(CtrlvARKOvWT.combined.BE, reduction = "pca", dims = 1:16)
DimPlot(CtrlvARKOvWT.combined.BE, reduction = "umap", label = TRUE, split.by = "stim")
#Rename
CtrlvARKOvWT.combined.BE <- RenameIdents(object = CtrlvARKOvWT.combined.BE, '3' = "BE1", '2' = "BE2",
                                           '0' = "BE3", '1' = "BE4", '4' = "BE5")
CtrlvARKOvWT.combined.BE[["BECellType"]] <- Idents(object = CtrlvARKOvWT.combined.BE)

tiff(file = "CtrlvARKOvWT.combined.BE UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKOvWT.combined.BE, reduction = "umap", pt.size = 0.3) 
dev.off()

Idents(object = CtrlvARKOvWT.combined.BE) <- "stim"
CtrlvARKOvWT.combined.BE_Ctrl <- subset(CtrlvARKOvWT.combined.BE, idents = c("Ctrl"))
CtrlvARKOvWT.combined.BE_ARKO <- subset(CtrlvARKOvWT.combined.BE, idents = c("ARKO"))
CtrlvARKOvWT.combined.BE_WT <- subset(CtrlvARKOvWT.combined.BE, idents = c("WT"))

DefaultAssay(CtrlvARKOvWT.combined.BE_Ctrl) <- "RNA"
tiff(file = "CtrlvARKOvWT.combined.BE_Ctrl hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_Ctrl, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_Ctrl Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_Ctrl, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_Ctrl Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_Ctrl, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_Ctrl Ccnd1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_Ctrl, reduction = "umap", features = c("Ccnd1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(CtrlvARKOvWT.combined.BE_ARKO) <- "RNA"
tiff(file = "CtrlvARKOvWT.combined.BE_ARKO hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_ARKO, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_ARKO Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_ARKO, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_ARKO Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_ARKO, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_ARKO Ccnd1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_ARKO, reduction = "umap", features = c("Ccnd1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(CtrlvARKOvWT.combined.BE_WT) <- "RNA"
tiff(file = "CtrlvARKOvWT.combined.BE_WT hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_WT, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "light grey"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_WT Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_WT, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_WT Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_WT, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKOvWT.combined.BE_WT Ccnd1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKOvWT.combined.BE_WT, reduction = "umap", features = c("Ccnd1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

####Comparison with WTEpi HiMYCEpi ARKOEpi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/NCB revision")

#### Merging Datasets Epi#### ####
Idents(object = WT) <- "seurat_clusters"
tiff(file = "WT UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()

DefaultAssay(WT) <- "RNA"
tiff(file = "WT Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "WT Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

WT.Epi <- subset(WT, idents = c("0", "4", "2", "7", "11", "6", "3", "9"))

Idents(object = CtrlvARKO.combined.Epi) <- "stim"
Ctrl.combined.Epi <- subset(CtrlvARKO.combined.Epi, idents = c("Ctrl"))
ARKO.combined.Epi <- subset(CtrlvARKO.combined.Epi, idents = c("ARKO"))

Ctrl.combined.Epi[["orig.clusters"]] <- Idents(object = Ctrl.combined.Epi)
ARKO.combined.Epi[["orig.clusters"]] <- Idents(object = ARKO.combined.Epi)
WT.Epi[["orig.clusters"]] <- Idents(object = WT.Epi)

Idents(object = Ctrl.combined.Epi) <- "Epicelltype"
Idents(object = ARKO.combined.Epi) <- "Epicelltype"
Idents(object = WT.Epi) <- "seurat_clusters"
Ctrl.combined.Epi$stim <- "HiMyc"
ARKO.combined.Epi$stim <- "ARKO"
WT.Epi$stim <- "WT"

HiMYCvWT.Epi.anchors <- FindIntegrationAnchors(object.list = list(Ctrl.combined.Epi, ARKO.combined.Epi, WT.Epi), dims = 1:20)
HiMYCvWT.Epi.combined <- IntegrateData(anchorset = HiMYCvWT.Epi.anchors, dims = 1:20)
DefaultAssay(HiMYCvWT.Epi.combined) <- "integrated"

#Run the standard workflow for visualization and clustering Ctrl2vARKO1.combined
HiMYCvWT.Epi.combined <- ScaleData(HiMYCvWT.Epi.combined, verbose = FALSE)
HiMYCvWT.Epi.combined <- RunPCA(HiMYCvWT.Epi.combined, npcs = 30, verbose = FALSE)
HiMYCvWT.Epi.combined <- RunUMAP(HiMYCvWT.Epi.combined, reduction = "pca", dims = 1:16)
DimPlot(HiMYCvWT.Epi.combined, reduction = "umap", pt.size = 0.3, label = TRUE) 

tiff(file = "HiMYCvWT.Epi.combined UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HiMYCvWT.Epi.combined, reduction = "umap", pt.size = 0.3) 
dev.off()
tiff(file = "HiMYCvWT.Epi.combined split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HiMYCvWT.Epi.combined, reduction = "umap", pt.size = 0.3, split.by= "stim", label = TRUE) 
dev.off()

Idents(object = HiMYCvWT.Epi.combined) <- "stim"
HiMYCvWT.Epi.combined_HiMyc <- subset(HiMYCvWT.Epi.combined, idents = c("HiMyc"))
HiMYCvWT.Epi.combined_ARKO <- subset(HiMYCvWT.Epi.combined, idents = c("ARKO"))
HiMYCvWT.Epi.combined_WT <- subset(HiMYCvWT.Epi.combined, idents = c("WT"))

DefaultAssay(HiMYCvWT.Epi.combined_HiMyc) <- "RNA"
tiff(file = "HiMYCvWT.Epi.combined_HiMyc hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HiMYCvWT.Epi.combined_HiMyc, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "HiMYCvWT.Epi.combined_HiMyc Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HiMYCvWT.Epi.combined_HiMyc, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(HiMYCvWT.Epi.combined_ARKO) <- "RNA"
tiff(file = "HiMYCvWT.Epi.combined_ARKO hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HiMYCvWT.Epi.combined_ARKO, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "HiMYCvWT.Epi.combined_ARKO Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HiMYCvWT.Epi.combined_ARKO, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()


DefaultAssay(HiMYCvWT.Epi.combined_WT) <- "RNA"
tiff(file = "HiMYCvWT.Epi.combined_WT hMycTg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HiMYCvWT.Epi.combined_WT, reduction = "umap", features = c('MYC-transgene'), cols = c("light grey", "light grey"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "HiMYCvWT.Epi.combined_WT Ar.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(HiMYCvWT.Epi.combined_WT, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Rename
Idents(object = HiMYCvWT.Epi.combined) <- "Epicelltype"

HiMYCvWT.Epi.combined <- RenameIdents(object = HiMYCvWT.Epi.combined, 'BE1' = "BE1", 'BE2' = "BE2",
                                           'BE3' = "BE3", 'BE4' = "BE4", 'LE1' = "LE1", 'LE2' = "LE2", 
                                           'LE3' = "LE3", 'LE4' = "LE4", 'LE5' = "LE5", 'LE6' = "LE6",
                                           'ProE' = "ProE", 'OE' = "OE")
HiMYCvWT.Epi.combined[["stimcelltype"]] <- Idents(object = HiMYCvWT.Epi.combined)
DimPlot(HiMYCvWT.Epi.combined, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Black", "Light grey")) 

tiff(file = "HiMYCvWT.Epi.combined WT grey UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HiMYCvWT.Epi.combined, reduction = "umap", pt.size = 0.3, cols = c("darkcyan", "saddlebrown", "red1", "Purple2", "Gold", "Green3", "hotpink", "#8494FF", "Blue", "Black", "Light grey")) 
dev.off()

###Merge WT, HiMYC, HiMYC-ARKO final version####
####CAF FB1####

combined.FB1 <- subset(combined.FB, idents = c("FB1"))
DimPlot(combined.FB1, reduction = "umap", split.by = "stim", pt.size = 0.3) 

DefaultAssay(combined.FB1) <- "RNA"
tiff(file = "combined.FB1 Ar.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Igf1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Igf1"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Igfbp3.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Igfbp3"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()

tiff(file = "combined.FB1 Pdgfrb.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Pdgfrb"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Twist1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Twist1"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Foxl1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Il11.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Il11"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Cxcl10.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Cxcl10"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB1 Sox9.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB1, reduction = "umap", features = c("Sox9"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()

combined.FB11 <- combined.FB1
DefaultAssay(combined.FB11) <- "integrated"

#Run the standard workflow for visualization and clustering
combined.FB11 <- ScaleData(combined.FB11, verbose = FALSE)
combined.FB11 <- RunPCA(combined.FB11, npcs = 30, verbose = FALSE)
ElbowPlot(combined.FB11, ndims = 50)
#Umap and Clustering
combined.FB11 <- FindNeighbors(combined.FB11, reduction = "pca", dims = 1:20)
combined.FB11 <- FindClusters(combined.FB11, resolution = 0.5)
combined.FB11 <- RunTSNE(combined.FB11, reduction = "pca", dims = 1:20)
combined.FB11 <- RunUMAP(combined.FB11, reduction = "pca", dims = 1:20)
DimPlot(combined.FB11, reduction = "umap", label = TRUE, split.by = "stim")

tiff(file = "combined.FB11 stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.FB11, reduction = "umap", split.by = "stim")
dev.off()

DefaultAssay(combined.FB11) <- "RNA"
tiff(file = "combined.FB11 Ar.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Igf1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Igf1"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Igfbp3.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Igfbp3"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()

tiff(file = "combined.FB11 Pdgfrb.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Pdgfrb"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Twist1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Twist1"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Foxl1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Il11.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Il11"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Cxcl10.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Cxcl10"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FB11 Sox9.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FB11, reduction = "umap", features = c("Sox9"), cols = c("light grey", "red"), split.by = "stim", pt.size = 1, min.cutoff = 0, max.cutoff = "q90")
dev.off()

####NC revision####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/NC revision")

#IGF1 expression
DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM IGF1 Strocelltype Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM, features = "Igf1", group.by = "Strocelltype", split.by = "stim", pt.size = 0.3, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM IGF1 stim Vln.tiff", width = 4, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined.FBSM, features = "Igf1", group.by = "stim",  pt.size = 0.3, cols = c("#3399FF", "#E06666"))
dev.off()

tiff(file = "CtrlvARKO.combined.FBSM Igf1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, reduction = "umap", split.by = "stim", features = c("Igf1"), cols = c("light grey", "red"), pt.size = 0.3, keep.scale = "all")
dev.off()

#AR downstream gene in BE/LE from HiMYC/HiMYC-ARKO
#Cell counts
Idents(object = CtrlvARKO.combined.Epi1) <- "celltype"
CtrlvARKO.combined.Epi1$celltype.stim <- paste(Idents(CtrlvARKO.combined.Epi1), CtrlvARKO.combined.Epi1$stim, sep = "_")
Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
table(Idents(CtrlvARKO.combined.Epi1))

Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
CtrlvARKO.combined.Epi1 <- RenameIdents(object = CtrlvARKO.combined.Epi1, 'BE_Ctrl' = "BE_Ctrl", 'LE_Ctrl' = "LE_Ctrl",
                                        'BE_ARKO' = "BE_ARKO", 'LE_ARKO' = "LE_ARKO")
DimPlot(CtrlvARKO.combined.Epi1, reduction = "umap", pt.size = 1, label = TRUE)
CtrlvARKO.combined.Epi1[["celltype.stim"]] <- Idents(object = CtrlvARKO.combined.Epi1)

DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
CtrlvARKO.combined.Epi1 <- ScaleData(CtrlvARKO.combined.Epi1, features = rownames(CtrlvARKO.combined.Epi1))
tiff(file = "CtrlvARKO.combined.Epi1 AR downstream Heatmap.tiff", width = 5, height = 15, units = "in", compression = "lzw", res = 200)
DoHeatmap(CtrlvARKO.combined.Epi1, features = c("Ar", "Myc",
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

#Violinplot for Gli1
Idents(object = CtrlvARKO.combined) <- "celltype"
CtrlvARKO.combined <- RenameIdents(object = CtrlvARKO.combined, 'BE' = "Epi", 'LE' = "Epi", 'FB' = "Stro", 'SM' = "Stro", 'Peri' = "Stro", 'Leu' = "Stro", 'Endo' = "Stro", 'Glia' = "Stro")
CtrlvARKO.combined[["EpiStro"]] <- Idents(object = CtrlvARKO.combined)

DefaultAssay(CtrlvARKO.combined) <- "RNA"
Idents(object = CtrlvARKO.combined) <- "EpiStro"
tiff(file = "CtrlvARKO.combined Gli1 EpiStro Vln.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined, features = "Gli1", group.by = "EpiStro", split.by = "stim", pt.size = 0.3, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined mGFP EpiStro Vln.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined, features = "EGFP", group.by = "EpiStro", split.by = "stim", pt.size = 0.3, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

DefaultAssay(CtrlvARKO.combined) <- "RNA"
Idents(object = CtrlvARKO.combined) <- "celltype"
tiff(file = "CtrlvARKO.combined Gli1 celltype Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined, features = "Gli1", group.by = "celltype", split.by = "stim", pt.size = 0.3, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()
tiff(file = "CtrlvARKO.combined mGFP celltype Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(CtrlvARKO.combined, features = "EGFP", group.by = "celltype", split.by = "stim", pt.size = 0.3, cols = c("#3399FF", "#E06666"), split.plot = TRUE)
dev.off()

#Expressionplot for Gli1
tiff(file = "CtrlvARKO.combined Gli1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q85", keep.scale = "all")
dev.off()

#UMI
Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
table(CtrlvARKO.combined.Epi1$nFeature_RNA)
write.csv(CtrlvARKO.combined.Epi1@assays[["RNA"]]@counts, file = "CtrlvARKO.combined.Epi1@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
BE_Ctrl <- subset(CtrlvARKO.combined.Epi1, idents = c("BE_Ctrl"))
write.csv(BE_Ctrl@assays[["RNA"]]@counts, file = "BE_Ctrl@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
LE_Ctrl <- subset(CtrlvARKO.combined.Epi1, idents = c("LE_Ctrl"))
write.csv(LE_Ctrl@assays[["RNA"]]@counts, file = "LE_Ctrl@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
BE_ARKO <- subset(CtrlvARKO.combined.Epi1, idents = c("BE_ARKO"))
write.csv(BE_ARKO@assays[["RNA"]]@counts, file = "BE_ARKO@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.Epi1) <- "celltype.stim"
LE_ARKO <- subset(CtrlvARKO.combined.Epi1, idents = c("LE_ARKO"))
write.csv(LE_ARKO@assays[["RNA"]]@counts, file = "LE_ARKO@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
FBSM_Ctrl <- subset(CtrlvARKO.combined.FBSM, idents = c("Ctrl"))
write.csv(FBSM_Ctrl@assays[["RNA"]]@counts, file = "FBSM_Ctrl@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.FBSM) <- "stim"
FBSM_ARKO <- subset(CtrlvARKO.combined.FBSM, idents = c("ARKO"))
write.csv(FBSM_ARKO@assays[["RNA"]]@counts, file = "FBSM_ARKO@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.Stro) <- "stim"
Stro_Ctrl <- subset(CtrlvARKO.combined.Stro, idents = c("Ctrl"))
write.csv(Stro_Ctrl@assays[["RNA"]]@counts, file = "Stro_Ctrl@assays@RNAcounts.csv")

Idents(object = CtrlvARKO.combined.Stro) <- "stim"
Stro_ARKO <- subset(CtrlvARKO.combined.Stro, idents = c("ARKO"))
write.csv(Stro_ARKO@assays[["RNA"]]@counts, file = "Stro_ARKO@assays@RNAcounts.csv")


#Expressionplot for Wnt downstream
DefaultAssay(CtrlvARKO.combined.Epi1) <- "RNA"
tiff(file = "CtrlvARKO.combined.Epi1 Ccnd1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", split.by = "stim", features = c("Ccnd1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q85", keep.scale = "all")
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Cd44 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", split.by = "stim", features = c("Cd44"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q85", keep.scale = "all")
dev.off()
tiff(file = "CtrlvARKO.combined.Epi1 Tcf7l2 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi1, reduction = "umap", split.by = "stim", features = c("Tcf7l2"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q85", keep.scale = "all")
dev.off()

#
sce.sling2 <- slingshot(HGFMETBcat_Multi1.epi1.1, clusterLabels = HGFMETBcat_Multi1.epi1$EpiCellTypes, reducedDim='UMAP', start.clus = "BE2", end.clus = "LE10")
pseudo.paths <- slingPseudotime(sce.sling2)
head(pseudo.paths)

reducedDim(sce.sling2, "UMAP") <- reducedDim(HGFMETBcat_Multi1.epi1.1, "UMAP")
shared.pseudo <- rowMeans(pseudo.paths, na.rm=TRUE)

gg <- plotUMAP(sce.sling2, colour_by=I(shared.pseudo))

embedded <- embedCurves(sce.sling2, "UMAP")
embedded <- slingCurves(embedded)
for (path in embedded) {
  embedded <- data.frame(path$s[path$ord,])
  gg <- gg + geom_path(data=embedded, aes(x=UMAP_1, y=UMAP_2), size=1.2)
}

tiff(file = "HGFMETBcat_Multi1.epi1 slingshot pseudotime final-1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
gg
dev.off()

#Slingshot
CtrlvARKO.combined.Epi1_sds <- slingshot(Embeddings(CtrlvARKO.combined.Epi1, "umap"), clusterLabels = CtrlvARKO.combined.Epi1$Epicelltype, 
                                        stretch = 2, start.clus = "BE2", end.clus = "LE5")

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors_clust <- cell_pal(CtrlvARKO.combined.Epi1$Epicelltype, hue_pal()) 

plot(reducedDim(CtrlvARKO.combined.Epi1_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(CtrlvARKO.combined.Epi1_sds, lwd = 2, col = 'black')

tiff(file = "CtrlvARKO.combined.Epi1 slingshot Pseudotime-1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot(reducedDim(CtrlvARKO.combined.Epi1_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(CtrlvARKO.combined.Epi1_sds, lwd = 2, col = 'black')
dev.off()


