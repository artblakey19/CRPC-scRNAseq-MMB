#ARQ9-mTmG-Gli1 at 3 timepoints workflow

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
library(ggpubr)
library(GGally)
library(scDblFinder)

####2M########

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M")

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

TwoM.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/v3_count__ARQ9_length_3191/v3_count_len3191_Project_COHP_46245_1_X3SC4/outs/filtered_feature_bc_matrix")
TwoM.unfiltered <- CreateSeuratObject(counts = TwoM.data,  min.cells = 3, min.features = 500, project = "ARQ9-Gli1_2M")
TwoM.unfiltered <- NormalizeData(TwoM.unfiltered)

#Initial processing & filtering

TwoM.unfiltered[["percent.mt"]] <- PercentageFeatureSet(TwoM.unfiltered, pattern = "^mt-")

tiff(file = "TwoM.unfiltered Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoM.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "TwoM.unfiltered Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(TwoM.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "TwoM.unfiltered Pre-filteration")
dev.off()
tiff(file = "TwoM.unfiltered Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(TwoM.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "TwoM.unfiltered Pre-filteration")
dev.off()

plot1 <- FeatureScatter(TwoM.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TwoM.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

TwoM <- subset(TwoM.unfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)

tiff(file = "TwoM Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "TwoM Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(TwoM@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "TwoM Post-filteration")
dev.off()
tiff(file = "TwoM Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(TwoM@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "TwoM Post-filteration")
dev.off()

plot1 <- FeatureScatter(TwoM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TwoM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

TwoM <- FindVariableFeatures(TwoM, selection.method = "vst", nfeatures = 5000)
tiff(file = "TwoM Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(TwoM)
dev.off()

#Clustering
all.genes <- rownames(TwoM)
TwoM <- ScaleData(TwoM, features = all.genes)
TwoM <- RunPCA(TwoM, features = VariableFeatures(object = TwoM))
print(TwoM[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(TwoM, dims = 1:2, reduction = "pca")
DimPlot(TwoM, reduction = "pca")

tiff(file = "TwoM ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(TwoM, ndims = 50)
dev.off()

TwoM <- FindNeighbors(TwoM, dims = 1:25)
TwoM <- FindClusters(TwoM, resolution = 0.5)
head(Idents(TwoM), 5)
TwoM <- RunTSNE(TwoM, reduction = "pca", dims = 1:25)
TwoM <- RunUMAP(TwoM, reduction = "pca", dims = 1:25)

Idents(object = TwoM) <- "seurat_clusters"
tiff(file = "TwoM seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3)
dev.off()

DefaultAssay(TwoM) <- "RNA"
FeaturePlot(TwoM, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(TwoM, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(TwoM, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(TwoM, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)

tiff(file = "TwoM hARtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM mGFP expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Ar expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Gli1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell Cycle Regression
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")

s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(TwoM) <- "RNA"
all.genes <- rownames(TwoM)
TwoM <- ScaleData(TwoM, features = all.genes)
TwoM <- CellCycleScoring(TwoM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoM) <- "Phase"
tiff(file = "TwoM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoM1 <- TwoM
DefaultAssay(TwoM1) <- "RNA"
TwoM1 <- ScaleData(TwoM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoM1))
TwoM1 <- RunPCA(TwoM1, features = VariableFeatures(TwoM1))
print(TwoM1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(TwoM1, dims = 1:2, reduction = "pca")
DimPlot(TwoM1, reduction = "pca")
ElbowPlot(TwoM1, ndims = 30)

TwoM1 <- FindNeighbors(TwoM1, reduction = "pca", dims = 1:20)
TwoM1 <- FindClusters(TwoM1, resolution = 0.5)
TwoM1 <- RunUMAP(TwoM1, reduction = "pca", dims = 1:20)
TwoM1 <- RunTSNE(TwoM1, reduction = "pca", dims = 1:20)
DimPlot(TwoM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoM1) <- "Phase"
tiff(file = "TwoM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoM1) <- "seurat_clusters"
tiff(file = "TwoM1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM1, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell type identification
DefaultAssay(TwoM1) <- "RNA"
Idents(object = TwoM1) <- "seurat_clusters"
TwoM1 <- ScaleData(TwoM1, features = rownames(TwoM1))
TwoM1.markers <- FindAllMarkers(TwoM1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoM1.markers, "TwoM1.seurat.markers.csv")

DefaultAssay(TwoM1) <- "RNA"
tiff(file = "TwoM1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", 
                                                    "Fbln1", "Myh11",  "Pecam1", 
                                                    "Mpz", "Rgs5", "Spink8", "Upk3a", "Htr3a", 
                                                    "Gcg", "Ccl5", "C1qa", "Myog"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = TwoM1) <- "seurat_clusters"
TwoM1 <- RenameIdents(object = TwoM1, '4' = "BE1", '6' = "BE2", '0' = "LE1", '2' = "LE2", '1' = "LE3", '9' = "LE4", '3' = "LE5", '5' = "LE6", '8' = "LE7", '7' = "UrLE", '10' = "OE")
TwoM1[["celltype"]] <- Idents(object = TwoM1)

#Umap
Idents(object = TwoM1) <- "celltype"
tiff(file = "TwoM1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM1, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell counts
Idents(object = TwoM1) <- "celltype"
table(Idents(TwoM1))

####P35 WT########

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35")

load("//isi-dcnl/user_data/zjsun/group/Adam Olson/AO Integrated Merging/Old files/AO ARKO-Gli1 Integrated Merge Environment.R.RData")

WTunfiltered[["percent.mt"]] <- PercentageFeatureSet(WTunfiltered, pattern = "^mt-")

tiff(file = "WTunfiltered Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WTunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WTunfiltered Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "WTunfiltered Pre-filteration")
dev.off()
tiff(file = "WTunfiltered Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "WTunfiltered Pre-filteration")
dev.off()

plot1 <- FeatureScatter(WTunfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WTunfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

WT <- subset(WTunfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)

tiff(file = "WT Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "WT Post-filteration")
dev.off()
tiff(file = "WT Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "WT Post-filteration")
dev.off()

plot1 <- FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2500)
tiff(file = "WT Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT)
dev.off()

#Clustering
all.genes <- rownames(WT)
WT <- ScaleData(WT, features = all.genes)
WT <- RunPCA(WT, features = VariableFeatures(object = WT))
print(WT[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WT, dims = 1:2, reduction = "pca")
DimPlot(WT, reduction = "pca")

tiff(file = "WT ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(WT, ndims = 50)
dev.off()

WT <- FindNeighbors(WT, dims = 1:25)
WT <- FindClusters(WT, resolution = 0.5)
head(Idents(WT), 5)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:25)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:25)

Idents(object = WT) <- "seurat_clusters"
tiff(file = "WT UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(WT) <- "RNA"
all.genes <- rownames(WT)
WT <- ScaleData(WT, features = all.genes)
WT <- CellCycleScoring(WT, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = WT) <- "Phase"
tiff(file = "WT Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
WT1 <- WT
DefaultAssay(WT1) <- "RNA"
WT1 <- ScaleData(WT1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(WT1))
WT1 <- RunPCA(WT1, features = VariableFeatures(WT1))
ElbowPlot(WT1, ndims = 30)

WT1 <- FindNeighbors(WT1, reduction = "pca", dims = 1:20)
WT1 <- FindClusters(WT1, resolution = 0.5)
WT1 <- RunUMAP(WT1, reduction = "pca", dims = 1:20)
WT1 <- RunTSNE(WT1, reduction = "pca", dims = 1:20)
DimPlot(WT1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = WT1) <- "Phase"
tiff(file = "WT1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = WT1) <- "seurat_clusters"
tiff(file = "WT1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT1, reduction = "umap", pt.size = 0.3)
dev.off()

#### Merging Datasets TwoMvWT ####

#Stash old idents
TwoM1[["orig.clusters"]] <- Idents(object = TwoM1)
WT1[["orig.clusters"]] <- Idents(object = WT1)

#Set Current idents
Idents(object = TwoM1) <- "seurat_clusters"
Idents(object = WT1) <- "seurat_clusters"
TwoM1$stim <- "ARQ9_2M"
WT1$stim <- "WT_P35"
TwoMvWT.anchors <- FindIntegrationAnchors(object.list = list(TwoM1, WT1), dims = 1:20)
TwoMvWT.combined <- IntegrateData(anchorset = TwoMvWT.anchors, dims = 1:20)

DefaultAssay(TwoMvWT.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined <- ScaleData(TwoMvWT.combined, verbose = FALSE)
TwoMvWT.combined <- RunPCA(TwoMvWT.combined, npcs = 30, verbose = FALSE)

#Umap and Clustering
TwoMvWT.combined <- FindNeighbors(TwoMvWT.combined, reduction = "pca", dims = 1:20)
TwoMvWT.combined <- FindClusters(TwoMvWT.combined, resolution = 0.5)
TwoMvWT.combined <- RunTSNE(TwoMvWT.combined, reduction = "pca", dims = 1:20)
TwoMvWT.combined <- RunUMAP(TwoMvWT.combined, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvWT.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined) <- "RNA"
all.genes <- rownames(TwoMvWT.combined)
TwoMvWT.combined <- ScaleData(TwoMvWT.combined, features = all.genes)
TwoMvWT.combined <- CellCycleScoring(TwoMvWT.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined) <- "Phase"
DimPlot(TwoMvWT.combined, reduction = "umap")
tiff(file = "TwoMvWT.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined1 <- TwoMvWT.combined
DefaultAssay(TwoMvWT.combined1) <- "integrated"
TwoMvWT.combined1 <- ScaleData(TwoMvWT.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined1))
TwoMvWT.combined1 <- RunPCA(TwoMvWT.combined1, features = VariableFeatures(TwoMvWT.combined1))
ElbowPlot(TwoMvWT.combined1, ndims = 50)

TwoMvWT.combined1 <- FindNeighbors(TwoMvWT.combined1, reduction = "pca", dims = 1:25)
TwoMvWT.combined1 <- FindClusters(TwoMvWT.combined1, resolution = 0.5)
TwoMvWT.combined1 <- RunUMAP(TwoMvWT.combined1, reduction = "pca", dims = 1:25)
TwoMvWT.combined1 <- RunTSNE(TwoMvWT.combined1, reduction = "pca", dims = 1:25)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWT.combined1) <- "Phase"
tiff(file = "TwoMvWT.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvWT.combined1) <- "RNA"
Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
TwoMvWT.combined1 <- ScaleData(TwoMvWT.combined1, features = rownames(TwoMvWT.combined1))
TwoMvWT.combined1.markers <- FindAllMarkers(TwoMvWT.combined1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoMvWT.combined1.markers, "TwoMvWT.combined1.seurat.markers.csv")

DefaultAssay(TwoMvWT.combined1) <- "RNA"
tiff(file = "TwoMvWT.combined1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2",
                                                    "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                    "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
TwoMvWT.combined1 <- RenameIdents(object = TwoMvWT.combined1, '0' = "BE", '10' = "BE", '12' = "BE", '1' = "LE", '3' = "LE", '6' = "LE", '8' = "LE", '9' = "LE", '4' = "SV", '2' = "FB", '7' = "FB",
                                  '5' = "SM", '13' = "VE", '15' = "Pericyte", '17' = "Glia", '11' = "Leu", '14' = "Lym", '16' = "Lym")
TwoMvWT.combined1[["CellType"]] <- Idents(object = TwoMvWT.combined1)

#Umap
Idents(object = TwoMvWT.combined1) <- "CellType"
tiff(file = "TwoMvWT.combined1 CellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined1 CellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined1) <- "stim"
TwoMvWT.combined1$stim.CellType <- paste(Idents(TwoMvWT.combined1), TwoMvWT.combined1$CellType, sep = "_")
Idents(object = TwoMvWT.combined1) <- "stim.CellType"
table(Idents(TwoMvWT.combined1))

####Sub-clustering FBSM####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35/FBSM")

Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.Stro <- subset(TwoMvWT.combined1, idents = c("FB", "SM", "VE", "Pericyte", "Glia", "Leu", "Lym"))
Idents(object = TwoMvWT.combined.Stro) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.Stro) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.Stro <- ScaleData(TwoMvWT.combined.Stro, verbose = FALSE)
TwoMvWT.combined.Stro <- RunPCA(TwoMvWT.combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.Stro, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.Stro <- FindNeighbors(TwoMvWT.combined.Stro, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Stro <- FindClusters(TwoMvWT.combined.Stro, resolution = 0.5)
TwoMvWT.combined.Stro <- RunTSNE(TwoMvWT.combined.Stro, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Stro <- RunUMAP(TwoMvWT.combined.Stro, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.Stro, reduction = "umap", pt.size = 0.3)

tiff(file = "TwoMvWT.combined.Stro UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Stro, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.Stro) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.Stro)
TwoMvWT.combined.Stro <- ScaleData(TwoMvWT.combined.Stro, features = all.genes)
TwoMvWT.combined.Stro <- CellCycleScoring(TwoMvWT.combined.Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.Stro) <- "Phase"
DimPlot(TwoMvWT.combined.Stro, reduction = "umap")
tiff(file = "TwoMvWT.combined.Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Stro, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.Stro1 <- TwoMvWT.combined.Stro
DefaultAssay(TwoMvWT.combined.Stro1) <- "integrated"
TwoMvWT.combined.Stro1 <- ScaleData(TwoMvWT.combined.Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.Stro1))
TwoMvWT.combined.Stro1 <- RunPCA(TwoMvWT.combined.Stro1, features = VariableFeatures(TwoMvWT.combined.Stro1))
ElbowPlot(TwoMvWT.combined.Stro1, ndims = 50)

TwoMvWT.combined.Stro1 <- FindNeighbors(TwoMvWT.combined.Stro1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Stro1 <- FindClusters(TwoMvWT.combined.Stro1, resolution = 0.5)
TwoMvWT.combined.Stro1 <- RunUMAP(TwoMvWT.combined.Stro1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Stro1 <- RunTSNE(TwoMvWT.combined.Stro1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.Stro1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvWT.combined.Stro1) <- "Phase"
tiff(file = "TwoMvWT.combined.Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Stro1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.Stro1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Stro1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

#Rename
DefaultAssay(TwoMvWT.combined.Stro1) <- "RNA"
tiff(file = "TwoMvWT.combined.Stro1 expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Stro1, reduction = "umap", features = c("Fbln1", "Myh11", "Pecam1", "Plp1", "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 1)
dev.off()

Idents(object = TwoMvWT.combined.Stro1) <- "seurat_clusters"
TwoMvWT.combined.Stro1 <- RenameIdents(object = TwoMvWT.combined.Stro1, '1' = "FB1", '3' = "FB2", '4' = "FB3", '0' = "FB4", '9' = "FB5", 
                                       '2' = "SM1", '5' = "SM2", '6' = "VE", '15' = "Pericyte", '12' = "Glia", '16' = "Glia", '11' = "Leu",
                                       '10' = "Leu", '14' = "Leu", '13' = "Lym", '7' = "Lym")
TwoMvWT.combined.Stro1[["StroCellType"]] <- Idents(object = TwoMvWT.combined.Stro1)

#Umap
Idents(object = TwoMvWT.combined.Stro1) <- "StroCellType"
tiff(file = "TwoMvWT.combined.Stro1 StroCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

#subset FBSM
Idents(object = TwoMvWT.combined.Stro1) <- "StroCellType"
TwoMvWT.combined.FBSM <- subset(TwoMvWT.combined.Stro1, idents = c("FB1", "FB2", "FB3", "FB4", "FB5", "SM1", "SM2"))
TwoMvWT.combined.FBSM <- RunTSNE(TwoMvWT.combined.FBSM, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FBSM <- RunUMAP(TwoMvWT.combined.FBSM, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.FBSM, reduction = "umap", pt.size = 0.3, split.by = "stim")

tiff(file = "TwoMvWT.combined.FBSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM, reduction = "umap", pt.size = 1)
dev.off()
tiff(file = "TwoMvWT.combined.FBSM split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.FBSM) <- "stim"
TwoMvWT.combined.FBSM$stim.StroCellType <- paste(Idents(TwoMvWT.combined.FBSM), TwoMvWT.combined.FBSM$StroCellType, sep = "_")
Idents(object = TwoMvWT.combined.FBSM) <- "stim.StroCellType"
table(Idents(TwoMvWT.combined.FBSM))

FeaturePlot(TwoMvWT.combined.FBSM, reduction = "umap", features = c("Fbln1", "Myh11"), cols = c("light grey", "red"), pt.size = 1)

#Violinplot
DefaultAssay(TwoMvWT.combined.FBSM) <- "RNA"
tiff(file = "TwoMvWT.combined.FBSM StroClusters EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FBSM, group.by = "StroCellType", split.by = "stim", features = c("EGFP"), pt.size = 0, split.plot = TRUE)
dev.off()
tiff(file = "TwoMvWT.combined.FBSM StroClusters Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FBSM, group.by = "StroCellType", split.by = "stim", features = c("Ar"), pt.size = 0, split.plot = TRUE)
dev.off()
tiff(file = "TwoMvWT.combined.FBSM StroClusters ARQ Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FBSM, group.by = "StroCellType", split.by = "stim", features = c("ARQ"), pt.size = 0, split.plot = TRUE)
dev.off()
tiff(file = "TwoMvWT.combined.FBSM StroClusters Gli1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FBSM, group.by = "StroCellType", split.by = "stim", features = c("Gli1"), pt.size = 0, split.plot = TRUE)
dev.off()

#EGFP expression
Idents(object = TwoMvWT.combined.FBSM) <- "stim"
TwoM.combined.FBSM <- subset(TwoMvWT.combined.FBSM, idents = c("ARQ9_2M"))
WT.combined.FBSM <- subset(TwoMvWT.combined.FBSM, idents = c("WT_P35"))

DefaultAssay(TwoM.combined.FBSM) <- "RNA"
tiff(file = "TwoM.combined.FBSM EGFP expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM Ar expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM ARQ expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

DefaultAssay(WT.combined.FBSM) <- "RNA"
tiff(file = "WT.combined.FBSM EGFP expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FBSM Ar expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FBSM ARQ expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM, reduction = "umap", features = c("ARQ"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()
####Sub-clustering FB-1####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35/FB-1")

Idents(object = TwoMvWT.combined.FBSM) <- "StroCellType"
TwoMvWT.combined.FB <- subset(TwoMvWT.combined.FBSM, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
TwoMvWT.combined.FB <- RunTSNE(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB <- RunUMAP(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3)

tiff(file = "TwoMvWT.combined.FB UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 1)
dev.off()

tiff(file = "TwoMvWT.combined.FB split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoMvWT.combined.FB$stim.StroCellType <- paste(Idents(TwoMvWT.combined.FB), TwoMvWT.combined.FB$StroCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
table(Idents(TwoMvWT.combined.FB))

#Violin Plots
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
tiff(file = "TwoMvWT.combined.FB StroClusters EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FB, group.by = "StroCellType", split.by = "stim", features = c("EGFP"), pt.size = 0, split.plot = TRUE)
dev.off()
tiff(file = "TwoMvWT.combined.FB StroClusters Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FB, group.by = "StroCellType", split.by = "stim", features = c("Ar"), pt.size = 0, split.plot = TRUE)
dev.off()
tiff(file = "TwoMvWT.combined.FB StroClusters ARQ Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TwoMvWT.combined.FB, group.by = "StroCellType", split.by = "stim", features = c("ARQ"), pt.size = 0, split.plot = TRUE)
dev.off()

#Feature plots
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoM.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("ARQ9_2M"))
WT.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("WT_P35"))

DefaultAssay(TwoM.combined.FB) <- "RNA"
tiff(file = "TwoM.combined.FB EGFP expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FB Ar expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FB ARQ expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

DefaultAssay(WT.combined.FB) <- "RNA"
tiff(file = "WT.combined.FB EGFP expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FB Ar expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FB ARQ expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB, reduction = "umap", features = c("ARQ"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Add EGFP info
DefaultAssay(TwoM.combined.FB) <- "RNA"
TwoM.combined.FBEGFPNeg <- subset(x=TwoM.combined.FB, subset = EGFP == 0)
TwoM.combined.FBEGFPPos <- subset(x=TwoM.combined.FB, subset = EGFP > 0)
Idents(object = TwoM.combined.FBEGFPNeg) <- "EGFPNeg"
Idents(object = TwoM.combined.FBEGFPPos) <- "EGFPPos"
TwoM.combined.FBEGFPNeg[["EGFPExp"]] <- Idents(object = TwoM.combined.FBEGFPNeg)
TwoM.combined.FBEGFPPos[["EGFPExp"]] <- Idents(object = TwoM.combined.FBEGFPPos)
TwoM.combined.FBEGFP <- merge(x = TwoM.combined.FBEGFPNeg, y = TwoM.combined.FBEGFPPos)
Idents(object = TwoM.combined.FBEGFP) <- "EGFPExp"
TwoM.combined.FB$EGFPExp <- Idents(object = TwoM.combined.FBEGFP)

#Add ARQ info
DefaultAssay(TwoM.combined.FB) <- "RNA"
TwoM.combined.FBARQNeg <- subset(x=TwoM.combined.FB, subset = ARQ == 0)
TwoM.combined.FBARQPos <- subset(x=TwoM.combined.FB, subset = ARQ > 0)
Idents(object = TwoM.combined.FBARQNeg) <- "ARQNeg"
Idents(object = TwoM.combined.FBARQPos) <- "ARQPos"
TwoM.combined.FBARQNeg[["ARQExp"]] <- Idents(object = TwoM.combined.FBARQNeg)
TwoM.combined.FBARQPos[["ARQExp"]] <- Idents(object = TwoM.combined.FBARQPos)
TwoM.combined.FBARQ <- merge(x = TwoM.combined.FBARQNeg, y = TwoM.combined.FBARQPos)
Idents(object = TwoM.combined.FBARQ) <- "ARQExp"
TwoM.combined.FB$ARQExp <- Idents(object = TwoM.combined.FBARQ)

Idents(object = TwoM.combined.FB) <- "EGFPExp"
TwoM.combined.FB$EGFPExp.ARQExp <- paste(Idents(TwoM.combined.FB), TwoM.combined.FB$ARQExp, sep = "_")
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp"
table(Idents(TwoM.combined.FB))

TwoM.combined.FB$EGFPExp.ARQExp.FBCellType <- paste(Idents(TwoM.combined.FB), TwoM.combined.FB$FBCellType, sep = "_")
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp.FBCellType"
table(Idents(TwoM.combined.FB))

#DEGs_Combined.FB_TwoMvWT
Idents(object = TwoMvWT.combined.FB) <- "stim"
DimPlot(TwoMvWT.combined.FB, reduction = "umap")
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB.0.Markers <- FindMarkers(TwoMvWT.combined.FB, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.FB.0.Markers, "TwoMvWT.combined.FB.0.Markers.csv")

#p.adjust
Combined.FB.TwoMvWT <- read.csv("TwoMvWT.combined.FB.0.Markers.csv") 
Combined.FB.TwoMvWT_pvalue <- Combined.FB.TwoMvWT$p_val
Combined.FB.TwoMvWT_pvalue=as.numeric(Combined.FB.TwoMvWT_pvalue)
Combined.FB.TwoMvWT_BH = p.adjust(Combined.FB.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB.TwoMvWT_BH, "Combined.FB.TwoMvWT_BH.csv")

#DEGs_TwoMvWT.combined.FB1
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB1 <- subset(TwoMvWT.combined.FB, idents = c("FB1"))

Idents(object = TwoMvWT.combined.FB1) <- "stim"
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
TwoMvWT.combined.FB1.0.Markers <- FindMarkers(TwoMvWT.combined.FB1, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.FB1.0.Markers, "TwoMvWT.combined.FB1.0.Markers.csv")

#p.adjust
Combined.FB1.TwoMvWT <- read.csv("TwoMvWT.combined.FB1.0.Markers.csv") 
Combined.FB1.TwoMvWT_pvalue <- Combined.FB1.TwoMvWT$p_val
Combined.FB1.TwoMvWT_pvalue=as.numeric(Combined.FB1.TwoMvWT_pvalue)
Combined.FB1.TwoMvWT_BH = p.adjust(Combined.FB1.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB1.TwoMvWT_BH, "Combined.FB1.TwoMvWT_BH.csv")

#DEGs_TwoMvWT.combined.FB2
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB2 <- subset(TwoMvWT.combined.FB, idents = c("FB2"))

Idents(object = TwoMvWT.combined.FB2) <- "stim"
DefaultAssay(TwoMvWT.combined.FB2) <- "RNA"
TwoMvWT.combined.FB2.0.Markers <- FindMarkers(TwoMvWT.combined.FB2, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.FB2.0.Markers, "TwoMvWT.combined.FB2.0.Markers.csv")

#p.adjust
Combined.FB2.TwoMvWT <- read.csv("TwoMvWT.combined.FB2.0.Markers.csv") 
Combined.FB2.TwoMvWT_pvalue <- Combined.FB2.TwoMvWT$p_val
Combined.FB2.TwoMvWT_pvalue=as.numeric(Combined.FB2.TwoMvWT_pvalue)
Combined.FB2.TwoMvWT_BH = p.adjust(Combined.FB2.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB2.TwoMvWT_BH, "Combined.FB2.TwoMvWT_BH.csv")

#DEGs_TwoMvWT.combined.FB3
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB3 <- subset(TwoMvWT.combined.FB, idents = c("FB3"))

Idents(object = TwoMvWT.combined.FB3) <- "stim"
DefaultAssay(TwoMvWT.combined.FB3) <- "RNA"
TwoMvWT.combined.FB3.0.Markers <- FindMarkers(TwoMvWT.combined.FB3, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.FB3.0.Markers, "TwoMvWT.combined.FB3.0.Markers.csv")

#p.adjust
Combined.FB3.TwoMvWT <- read.csv("TwoMvWT.combined.FB3.0.Markers.csv") 
Combined.FB3.TwoMvWT_pvalue <- Combined.FB3.TwoMvWT$p_val
Combined.FB3.TwoMvWT_pvalue=as.numeric(Combined.FB3.TwoMvWT_pvalue)
Combined.FB3.TwoMvWT_BH = p.adjust(Combined.FB3.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB3.TwoMvWT_BH, "Combined.FB3.TwoMvWT_BH.csv")

#DEGs_TwoMvWT.combined.FB4
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB4 <- subset(TwoMvWT.combined.FB, idents = c("FB4"))

Idents(object = TwoMvWT.combined.FB4) <- "stim"
DefaultAssay(TwoMvWT.combined.FB4) <- "RNA"
TwoMvWT.combined.FB4.0.Markers <- FindMarkers(TwoMvWT.combined.FB4, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.FB4.0.Markers, "TwoMvWT.combined.FB4.0.Markers.csv")

#p.adjust
Combined.FB4.TwoMvWT <- read.csv("TwoMvWT.combined.FB4.0.Markers.csv") 
Combined.FB4.TwoMvWT_pvalue <- Combined.FB4.TwoMvWT$p_val
Combined.FB4.TwoMvWT_pvalue=as.numeric(Combined.FB4.TwoMvWT_pvalue)
Combined.FB4.TwoMvWT_BH = p.adjust(Combined.FB4.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB4.TwoMvWT_BH, "Combined.FB4.TwoMvWT_BH.csv")

#DEGs_TwoMvWT.combined.FB5
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB5 <- subset(TwoMvWT.combined.FB, idents = c("FB5"))

Idents(object = TwoMvWT.combined.FB5) <- "stim"
DefaultAssay(TwoMvWT.combined.FB5) <- "RNA"
TwoMvWT.combined.FB5.0.Markers <- FindMarkers(TwoMvWT.combined.FB5, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.FB5.0.Markers, "TwoMvWT.combined.FB5.0.Markers.csv")

#p.adjust
Combined.FB5.TwoMvWT <- read.csv("TwoMvWT.combined.FB5.0.Markers.csv") 
Combined.FB5.TwoMvWT_pvalue <- Combined.FB5.TwoMvWT$p_val
Combined.FB5.TwoMvWT_pvalue=as.numeric(Combined.FB5.TwoMvWT_pvalue)
Combined.FB5.TwoMvWT_BH = p.adjust(Combined.FB5.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB5.TwoMvWT_BH, "Combined.FB5.TwoMvWT_BH.csv")

#Heatmap
Idents(object = TwoMvWT.combined.FB4) <- "stim"
tiff(file = "TwoMvWT.combined.FB4 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB4, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB4 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB4, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt11", "Wnt5a", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Igfbp3", "Igfbp4", "Igfbp5", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
TwoMvWT.combined.FB.markers <- FindAllMarkers(TwoMvWT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FBTop50 <- TwoMvWT.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoMvWT.combined.FBTop50, "TwoMvWT.combined.FBTop50.csv")

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3)
TwoMvWT.combined.FB <- RenameIdents(object = TwoMvWT.combined.FB, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5",
                                     'WT_P35_FB1' = "WT_P35_FB1", 'WT_P35_FB2' = "WT_P35_FB2", 'WT_P35_FB3' = "WT_P35_FB3", 'WT_P35_FB4' = "WT_P35_FB4", 'WT_P35_FB5' = "WT_P35_FB5")
TwoMvWT.combined.FB[["stim.StroCellType"]] <- Idents(object = TwoMvWT.combined.FB)
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
tiff(file = "TwoMvWT.combined.FB stim Heatmap Top50 purple.tiff", width = 16, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
TwoMvWT.combined.FB.markers.1 <- FindAllMarkers(TwoMvWT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FBTop50.1 <- TwoMvWT.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB stim individual Heatmap Top50 purple.tiff", width = 16, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoMvWT.combined.FBTop50.1, "TwoMvWT.combined.FBTop50.1.csv")

#Heatmap_AllFB_TwoM
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoM.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("ARQ9_2M"))
WT.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("WT_P35"))

DefaultAssay(TwoM.combined.FB) <- "RNA"
Idents(object = TwoM.combined.FB) <- "StroCellType"
TwoM.combined.FB <- ScaleData(TwoM.combined.FB, features = rownames(TwoM.combined.FB))
TwoM.combined.FB.markers <- FindAllMarkers(TwoM.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoM.combined.FBTop50 <- TwoM.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoM.combined.FB Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.combined.FB, features = c("ARQ", "EGFP", TwoM.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoM.combined.FBTop50, "TwoM.combined.FBTop50.csv")

#Heatmap_AllFB_WT
DefaultAssay(WT.combined.FB) <- "RNA"
Idents(object = WT.combined.FB) <- "StroCellType"
WT.combined.FB <- ScaleData(WT.combined.FB, features = rownames(WT.combined.FB))
WT.combined.FB.markers <- FindAllMarkers(WT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
WT.combined.FBTop50 <- WT.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "WT.combined.FB Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(WT.combined.FB, features = c("ARQ", "EGFP", WT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(WT.combined.FBTop50, "WT.combined.FBTop50.csv")

####Sub-clustering Epi-1####

Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.Epi <- subset(TwoMvWT.combined1, idents = c("BE", "LE", "SV"))
Idents(object = TwoMvWT.combined.Epi) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.Epi <- ScaleData(TwoMvWT.combined.Epi, verbose = FALSE)
TwoMvWT.combined.Epi <- RunPCA(TwoMvWT.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.Epi, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.Epi <- FindNeighbors(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:16)
TwoMvWT.combined.Epi <- FindClusters(TwoMvWT.combined.Epi, resolution = 0.5)
TwoMvWT.combined.Epi <- RunTSNE(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:16)
TwoMvWT.combined.Epi <- RunUMAP(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:16)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, split.by = "stim")

tiff(file = "TwoMvWT.combined.Epi UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.Epi) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.Epi)
TwoMvWT.combined.Epi <- ScaleData(TwoMvWT.combined.Epi, features = all.genes)
TwoMvWT.combined.Epi <- CellCycleScoring(TwoMvWT.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.Epi) <- "Phase"
DimPlot(TwoMvWT.combined.Epi, reduction = "umap")
tiff(file = "TwoMvWT.combined.Epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.Epi1 <- TwoMvWT.combined.Epi
DefaultAssay(TwoMvWT.combined.Epi1) <- "integrated"
TwoMvWT.combined.Epi1 <- ScaleData(TwoMvWT.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.Epi1))
TwoMvWT.combined.Epi1 <- RunPCA(TwoMvWT.combined.Epi1, features = VariableFeatures(TwoMvWT.combined.Epi1))
ElbowPlot(TwoMvWT.combined.Epi1, ndims = 50)

TwoMvWT.combined.Epi1 <- FindNeighbors(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Epi1 <- FindClusters(TwoMvWT.combined.Epi1, resolution = 0.5)
TwoMvWT.combined.Epi1 <- RunUMAP(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Epi1 <- RunTSNE(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvWT.combined.Epi1) <- "Phase"
tiff(file = "TwoMvWT.combined.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
tiff(file = "TwoMvWT.combined.Epi1 Epicelltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Epi1, reduction = "umap", features = c("Krt5", "Krt19", "Fbln1", "Tagln",
                                                                    "Gsdma", "Msmb", "Lrrc26", "Krt4",
                                                                    "Clu", "Svs5", "Ly6a", "Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
TwoMvWT.combined.Epi1 <- RenameIdents(object = TwoMvWT.combined.Epi1, '0' = "BE1", '6' = "BE2", '12' = "BE3", '1' = "LE1", '4' = "LE2", '10' = "LE3", '9' = "LE4", '7' = "LE5", '2' = "LE6", 
                                      '14' = "LE7",
                                      '5' = "LumP", '11' = "UrLE", '3' = "SV", '13' = "SV", '8' = "OE", '15' = "OE")
TwoMvWT.combined.Epi1[["EpiCellType"]] <- Idents(object = TwoMvWT.combined.Epi1)

Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
tiff(file = "TwoMvWT.combined.Epi1 EpiCellType  UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
tiff(file = "TwoMvWT.combined.Epi1 EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.Epi1) <- "stim"
TwoMvWT.combined.Epi1$stim.EpiCellType <- paste(Idents(TwoMvWT.combined.Epi1), TwoMvWT.combined.Epi1$EpiCellType, sep = "_")
Idents(object = TwoMvWT.combined.Epi1) <- "stim.EpiCellType"
table(Idents(TwoMvWT.combined.Epi1))

#DEGs_Combined.BE_TwoMvWT
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
TwoMvWT.combined.BE <- subset(TwoMvWT.combined.Epi1, idents = c("BE1", "BE2", "BE3"))

Idents(object = TwoMvWT.combined.BE) <- "stim"

DefaultAssay(TwoMvWT.combined.BE) <- "RNA"
TwoMvWT.combined.BE.0.Markers <- FindMarkers(TwoMvWT.combined.BE, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.BE.0.Markers, "TwoMvWT.combined.BE.0.Markers.csv")

#p.adjust
Combined.BE.TwoMvWT <- read.csv("TwoMvWT.combined.BE.0.Markers.csv") 
Combined.BE.TwoMvWT_pvalue <- Combined.BE.TwoMvWT$p_val
Combined.BE.TwoMvWT_pvalue=as.numeric(Combined.BE.TwoMvWT_pvalue)
Combined.BE.TwoMvWT_BH = p.adjust(Combined.BE.TwoMvWT_pvalue, "BH")
write.csv(Combined.BE.TwoMvWT_BH, "Combined.BE.TwoMvWT_BH.csv")

#DEGs_Combined.BE1_TwoMvWT
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
TwoMvWT.combined.BE1 <- subset(TwoMvWT.combined.Epi1, idents = c("BE1"))

Idents(object = TwoMvWT.combined.BE1) <- "stim"
DefaultAssay(TwoMvWT.combined.BE1) <- "RNA"
TwoMvWT.combined.BE1.0.Markers <- FindMarkers(TwoMvWT.combined.BE1, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.BE1.0.Markers, "TwoMvWT.combined.BE1.0.Markers.csv")

#p.adjust
Combined.BE1.TwoMvWT <- read.csv("TwoMvWT.combined.BE1.0.Markers.csv") 
Combined.BE1.TwoMvWT_pvalue <- Combined.BE1.TwoMvWT$p_val
Combined.BE1.TwoMvWT_pvalue=as.numeric(Combined.BE1.TwoMvWT_pvalue)
Combined.BE1.TwoMvWT_BH = p.adjust(Combined.BE1.TwoMvWT_pvalue, "BH")
write.csv(Combined.BE1.TwoMvWT_BH, "Combined.BE1.TwoMvWT_BH.csv")

#DEGs_Combined.BE2_TwoMvWT
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
TwoMvWT.combined.BE2 <- subset(TwoMvWT.combined.Epi1, idents = c("BE2"))

Idents(object = TwoMvWT.combined.BE2) <- "stim"
DefaultAssay(TwoMvWT.combined.BE2) <- "RNA"
TwoMvWT.combined.BE2.0.Markers <- FindMarkers(TwoMvWT.combined.BE2, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.BE2.0.Markers, "TwoMvWT.combined.BE2.0.Markers.csv")

#p.adjust
Combined.BE2.TwoMvWT <- read.csv("TwoMvWT.combined.BE2.0.Markers.csv") 
Combined.BE2.TwoMvWT_pvalue <- Combined.BE2.TwoMvWT$p_val
Combined.BE2.TwoMvWT_pvalue=as.numeric(Combined.BE2.TwoMvWT_pvalue)
Combined.BE2.TwoMvWT_BH = p.adjust(Combined.BE2.TwoMvWT_pvalue, "BH")
write.csv(Combined.BE2.TwoMvWT_BH, "Combined.BE2.TwoMvWT_BH.csv")

#DEGs_Combined.BE3_TwoMvWT
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
TwoMvWT.combined.BE3 <- subset(TwoMvWT.combined.Epi1, idents = c("BE3"))

Idents(object = TwoMvWT.combined.BE3) <- "stim"
DefaultAssay(TwoMvWT.combined.BE3) <- "RNA"
TwoMvWT.combined.BE3.0.Markers <- FindMarkers(TwoMvWT.combined.BE3, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(TwoMvWT.combined.BE3.0.Markers, "TwoMvWT.combined.BE3.0.Markers.csv")

#p.adjust
Combined.BE3.TwoMvWT <- read.csv("TwoMvWT.combined.BE3.0.Markers.csv") 
Combined.BE3.TwoMvWT_pvalue <- Combined.BE3.TwoMvWT$p_val
Combined.BE3.TwoMvWT_pvalue=as.numeric(Combined.BE3.TwoMvWT_pvalue)
Combined.BE3.TwoMvWT_BH = p.adjust(Combined.BE3.TwoMvWT_pvalue, "BH")
write.csv(Combined.BE3.TwoMvWT_BH, "Combined.BE3.TwoMvWT_BH.csv")

#Heatmap_BE1
DefaultAssay(TwoMvWT.combined.BE1) <- "RNA"
Idents(object = TwoMvWT.combined.BE1) <- "stim"
TwoMvWT.combined.BE1 <- ScaleData(TwoMvWT.combined.BE1, features = rownames(TwoMvWT.combined.BE1))
TwoMvWT.combined.BE1.markers <- FindAllMarkers(TwoMvWT.combined.BE1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.BE1Top100 <- TwoMvWT.combined.BE1.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.BE1 Heatmap Top100 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.BE1, features = c(TwoMvWT.combined.BE1Top100$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
#Heatmap_BE2
DefaultAssay(TwoMvWT.combined.BE2) <- "RNA"
Idents(object = TwoMvWT.combined.BE2) <- "stim"
TwoMvWT.combined.BE2 <- ScaleData(TwoMvWT.combined.BE2, features = rownames(TwoMvWT.combined.BE2))
TwoMvWT.combined.BE2.markers <- FindAllMarkers(TwoMvWT.combined.BE2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.BE2Top100 <- TwoMvWT.combined.BE2.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.BE2 Heatmap Top100 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.BE2, features = c(TwoMvWT.combined.BE2Top100$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
#Heatmap_BE3
DefaultAssay(TwoMvWT.combined.BE3) <- "RNA"
Idents(object = TwoMvWT.combined.BE3) <- "stim"
TwoMvWT.combined.BE3 <- ScaleData(TwoMvWT.combined.BE3, features = rownames(TwoMvWT.combined.BE3))
TwoMvWT.combined.BE3.markers <- FindAllMarkers(TwoMvWT.combined.BE3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.BE3Top100 <- TwoMvWT.combined.BE3.markers %>% group_by(cluster) %>% top_n(n = 100, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.BE3 Heatmap Top100 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.BE3, features = c(TwoMvWT.combined.BE3Top100$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoMvWT.combined.BE1Top50, "TwoMvWT.combined.BE1Top50.csv")
write.csv(TwoMvWT.combined.BE2Top50, "TwoMvWT.combined.BE2Top50.csv")
write.csv(TwoMvWT.combined.BE3Top50, "TwoMvWT.combined.BE3Top50.csv")
write.csv(TwoMvWT.combined.BE1Top100, "TwoMvWT.combined.BE1Top100.csv")
write.csv(TwoMvWT.combined.BE2Top100, "TwoMvWT.combined.BE2Top100.csv")
write.csv(TwoMvWT.combined.BE3Top100, "TwoMvWT.combined.BE3Top100.csv")

####Sub-clustering FB-2####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35/FB")

Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.FB <- subset(TwoMvWT.combined1, idents = c("FB"))
Idents(object = TwoMvWT.combined.FB) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.FB) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, verbose = FALSE)
TwoMvWT.combined.FB <- RunPCA(TwoMvWT.combined.FB, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.FB, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.FB <- FindNeighbors(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB <- FindClusters(TwoMvWT.combined.FB, resolution = 0.5)
TwoMvWT.combined.FB <- RunTSNE(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB <- RunUMAP(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvWT.combined.FB UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.FB)
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = all.genes)
TwoMvWT.combined.FB <- CellCycleScoring(TwoMvWT.combined.FB, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.FB) <- "Phase"
DimPlot(TwoMvWT.combined.FB, reduction = "umap")
tiff(file = "TwoMvWT.combined.FB Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.FB1 <- TwoMvWT.combined.FB
DefaultAssay(TwoMvWT.combined.FB1) <- "integrated"
TwoMvWT.combined.FB1 <- ScaleData(TwoMvWT.combined.FB1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.FB1))
TwoMvWT.combined.FB1 <- RunPCA(TwoMvWT.combined.FB1, features = VariableFeatures(TwoMvWT.combined.FB1))
ElbowPlot(TwoMvWT.combined.FB1, ndims = 50)

TwoMvWT.combined.FB1 <- FindNeighbors(TwoMvWT.combined.FB1, reduction = "pca", dims = 1:22)
TwoMvWT.combined.FB1 <- FindClusters(TwoMvWT.combined.FB1, resolution = 0.4)
TwoMvWT.combined.FB1 <- RunUMAP(TwoMvWT.combined.FB1, reduction = "pca", dims = 1:22)
TwoMvWT.combined.FB1 <- RunTSNE(TwoMvWT.combined.FB1, reduction = "pca", dims = 1:22)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvWT.combined.FB1) <- "Phase"
tiff(file = "TwoMvWT.combined.FB1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FB1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FB1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Rename
Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
TwoMvWT.combined.FB1 <- RenameIdents(object = TwoMvWT.combined.FB1, '0' = "FB1", '3' = "FB2", '1' = "FB3", '2' = "FB4", '5' = "FB5", '4' = "FB6", '6' = "FB7")
TwoMvWT.combined.FB1[["FBCellType"]] <- Idents(object = TwoMvWT.combined.FB1)

#Umap
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
tiff(file = "TwoMvWT.combined.FB1 FBCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined.FB1 FBCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 1, split.by = "stim", label = TRUE)
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.FB1) <- "stim"
TwoMvWT.combined.FB1$stim.FBCellType <- paste(Idents(TwoMvWT.combined.FB1), TwoMvWT.combined.FB1$FBCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB1) <- "stim.FBCellType"
table(Idents(TwoMvWT.combined.FB1))

#EGFP expression
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
tiff(file = "TwoMvWT.combined.FB1 EGFP expression stim plots.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, split.by = "stim")
dev.off()

#Add EGFP info
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
TwoMvWT.combined.FB1EGFPNeg <- subset(x=TwoMvWT.combined.FB1, subset = EGFP < 1)
TwoMvWT.combined.FB1EGFPPos <- subset(x=TwoMvWT.combined.FB1, subset = EGFP >= 1)
Idents(object = TwoMvWT.combined.FB1EGFPNeg) <- "EGFPNeg"
Idents(object = TwoMvWT.combined.FB1EGFPPos) <- "EGFPPos"
TwoMvWT.combined.FB1EGFPNeg[["EGFPExp"]] <- Idents(object = TwoMvWT.combined.FB1EGFPNeg)
TwoMvWT.combined.FB1EGFPPos[["EGFPExp"]] <- Idents(object = TwoMvWT.combined.FB1EGFPPos)
TwoMvWT.combined.FB1EGFP <- merge(x = TwoMvWT.combined.FB1EGFPNeg, y = TwoMvWT.combined.FB1EGFPPos)
Idents(object = TwoMvWT.combined.FB1EGFP) <- "EGFPExp"
TwoMvWT.combined.FB1$EGFPExp <- Idents(object = TwoMvWT.combined.FB1EGFP)

Idents(object = TwoMvWT.combined.FB1) <- "EGFPExp"

tiff(file = "TwoMvWT.combined.FB1 EGFPExp split UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 1)
dev.off()

TwoMvWT.combined.EGFP.FB <- subset(TwoMvWT.combined.FB1, idents = c("EGFPPos"))
Idents(object = TwoMvWT.combined.EGFP.FB) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.EGFP.FB) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.EGFP.FB <- ScaleData(TwoMvWT.combined.EGFP.FB, verbose = FALSE)
TwoMvWT.combined.EGFP.FB <- RunPCA(TwoMvWT.combined.EGFP.FB, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.EGFP.FB, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.EGFP.FB <- FindNeighbors(TwoMvWT.combined.EGFP.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.EGFP.FB <- FindClusters(TwoMvWT.combined.EGFP.FB, resolution = 0.5)
TwoMvWT.combined.EGFP.FB <- RunTSNE(TwoMvWT.combined.EGFP.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.EGFP.FB <- RunUMAP(TwoMvWT.combined.EGFP.FB, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvWT.combined.EGFP.FB UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.EGFP.FB) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.EGFP.FB)
TwoMvWT.combined.EGFP.FB <- ScaleData(TwoMvWT.combined.EGFP.FB, features = all.genes)
TwoMvWT.combined.EGFP.FB <- CellCycleScoring(TwoMvWT.combined.EGFP.FB, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.EGFP.FB) <- "Phase"
DimPlot(TwoMvWT.combined.EGFP.FB, reduction = "umap")
tiff(file = "TwoMvWT.combined.EGFP.FB Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.EGFP.FB1 <- TwoMvWT.combined.EGFP.FB
DefaultAssay(TwoMvWT.combined.EGFP.FB1) <- "integrated"
TwoMvWT.combined.EGFP.FB1 <- ScaleData(TwoMvWT.combined.EGFP.FB1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.EGFP.FB1))
TwoMvWT.combined.EGFP.FB1 <- RunPCA(TwoMvWT.combined.EGFP.FB1, features = VariableFeatures(TwoMvWT.combined.EGFP.FB1))
ElbowPlot(TwoMvWT.combined.EGFP.FB1, ndims = 50)

TwoMvWT.combined.EGFP.FB1 <- FindNeighbors(TwoMvWT.combined.EGFP.FB1, reduction = "pca", dims = 1:22)
TwoMvWT.combined.EGFP.FB1 <- FindClusters(TwoMvWT.combined.EGFP.FB1, resolution = 0.4)
TwoMvWT.combined.EGFP.FB1 <- RunUMAP(TwoMvWT.combined.EGFP.FB1, reduction = "pca", dims = 1:22)
TwoMvWT.combined.EGFP.FB1 <- RunTSNE(TwoMvWT.combined.EGFP.FB1, reduction = "pca", dims = 1:22)
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvWT.combined.EGFP.FB1) <- "Phase"
tiff(file = "TwoMvWT.combined.EGFP.FB1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.EGFP.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.EGFP.FB1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.EGFP.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.EGFP.FB1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Rename
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "seurat_clusters"
TwoMvWT.combined.EGFP.FB1 <- RenameIdents(object = TwoMvWT.combined.EGFP.FB1, '0' = "FB1", '2' = "FB2", '3' = "FB3", '1' = "FB4", '4' = "FB5", '5' = "FB6", '6' = "FB7")
TwoMvWT.combined.EGFP.FB1[["EGFPFBCellType"]] <- Idents(object = TwoMvWT.combined.EGFP.FB1)

#Umap
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
tiff(file = "TwoMvWT.combined.EGFP.FB1 EGFPFBCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined.EGFP.FB1 EGFPFBCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap", pt.size = 1, split.by = "stim", label = TRUE)
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "stim"
TwoMvWT.combined.EGFP.FB1$stim.EGFPFBCellType <- paste(Idents(TwoMvWT.combined.EGFP.FB1), TwoMvWT.combined.EGFP.FB1$EGFPFBCellType, sep = "_")
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "stim.EGFPFBCellType"
table(Idents(TwoMvWT.combined.EGFP.FB1))

#DEGs_Combined.EGFPFB_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "stim"
DimPlot(TwoMvWT.combined.EGFP.FB1, reduction = "umap")
DefaultAssay(TwoMvWT.combined.EGFP.FB1) <- "RNA"
Combined.EGFPFB.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB1, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB.TwoMvWT.0.Markers, "Combined.EGFPFB.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB.TwoMvWT <- read.csv("Combined.EGFPFB.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB.TwoMvWT_pvalue <- Combined.EGFPFB.TwoMvWT$p_val
Combined.EGFPFB.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB.TwoMvWT_pvalue)
Combined.EGFPFB.TwoMvWT_BH = p.adjust(Combined.EGFPFB.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB.TwoMvWT_BH, "Combined.EGFPFB.TwoMvWT_BH.csv")


#DEGs_Combined.EGFPFB1_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB11 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB1"))

Idents(object = TwoMvWT.combined.EGFP.FB11) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB11) <- "RNA"
Combined.EGFPFB1.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB11, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB1.TwoMvWT.0.Markers, "Combined.EGFPFB1.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB1.TwoMvWT <- read.csv("Combined.EGFPFB1.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB1.TwoMvWT_pvalue <- Combined.EGFPFB1.TwoMvWT$p_val
Combined.EGFPFB1.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB1.TwoMvWT_pvalue)
Combined.EGFPFB1.TwoMvWT_BH = p.adjust(Combined.EGFPFB1.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB1.TwoMvWT_BH, "Combined.EGFPFB1.TwoMvWT_BH.csv")

#DEGs_Combined.EGFPFB2_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB2 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB2"))

Idents(object = TwoMvWT.combined.EGFP.FB2) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB2) <- "RNA"
Combined.EGFPFB2.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB2, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB2.TwoMvWT.0.Markers, "Combined.EGFPFB2.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB2.TwoMvWT <- read.csv("Combined.EGFPFB2.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB2.TwoMvWT_pvalue <- Combined.EGFPFB2.TwoMvWT$p_val
Combined.EGFPFB2.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB2.TwoMvWT_pvalue)
Combined.EGFPFB2.TwoMvWT_BH = p.adjust(Combined.EGFPFB2.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB2.TwoMvWT_BH, "Combined.EGFPFB2.TwoMvWT_BH.csv")

#DEGs_Combined.EGFPFB3_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB3 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB3"))

Idents(object = TwoMvWT.combined.EGFP.FB3) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB3) <- "RNA"
Combined.EGFPFB3.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB3, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB3.TwoMvWT.0.Markers, "Combined.EGFPFB3.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB3.TwoMvWT <- read.csv("Combined.EGFPFB3.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB3.TwoMvWT_pvalue <- Combined.EGFPFB3.TwoMvWT$p_val
Combined.EGFPFB3.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB3.TwoMvWT_pvalue)
Combined.EGFPFB3.TwoMvWT_BH = p.adjust(Combined.EGFPFB3.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB3.TwoMvWT_BH, "Combined.EGFPFB3.TwoMvWT_BH.csv")

#DEGs_Combined.EGFPFB4_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB4 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB4"))

Idents(object = TwoMvWT.combined.EGFP.FB4) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB4) <- "RNA"
Combined.EGFPFB4.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB4, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB4.TwoMvWT.0.Markers, "Combined.EGFPFB4.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB4.TwoMvWT <- read.csv("Combined.EGFPFB4.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB4.TwoMvWT_pvalue <- Combined.EGFPFB4.TwoMvWT$p_val
Combined.EGFPFB4.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB4.TwoMvWT_pvalue)
Combined.EGFPFB4.TwoMvWT_BH = p.adjust(Combined.EGFPFB4.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB4.TwoMvWT_BH, "Combined.EGFPFB4.TwoMvWT_BH.csv")

#DEGs_Combined.EGFPFB5_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB5 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB5"))

Idents(object = TwoMvWT.combined.EGFP.FB5) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB5) <- "RNA"
Combined.EGFPFB5.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB5, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB5.TwoMvWT.0.Markers, "Combined.EGFPFB5.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB5.TwoMvWT <- read.csv("Combined.EGFPFB5.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB5.TwoMvWT_pvalue <- Combined.EGFPFB5.TwoMvWT$p_val
Combined.EGFPFB5.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB5.TwoMvWT_pvalue)
Combined.EGFPFB5.TwoMvWT_BH = p.adjust(Combined.EGFPFB5.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB5.TwoMvWT_BH, "Combined.EGFPFB5.TwoMvWT_BH.csv")

#DEGs_Combined.EGFPFB6_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB6 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB6"))

Idents(object = TwoMvWT.combined.EGFP.FB6) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB6) <- "RNA"
Combined.EGFPFB6.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB6, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB6.TwoMvWT.0.Markers, "Combined.EGFPFB6.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB6.TwoMvWT <- read.csv("Combined.EGFPFB6.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB6.TwoMvWT_pvalue <- Combined.EGFPFB6.TwoMvWT$p_val
Combined.EGFPFB6.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB6.TwoMvWT_pvalue)
Combined.EGFPFB6.TwoMvWT_BH = p.adjust(Combined.EGFPFB6.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB6.TwoMvWT_BH, "Combined.EGFPFB6.TwoMvWT_BH.csv")

#DEGs_Combined.EGFPFB7_TwoMvWT
Idents(object = TwoMvWT.combined.EGFP.FB1) <- "EGFPFBCellType"
TwoMvWT.combined.EGFP.FB7 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("FB7"))

Idents(object = TwoMvWT.combined.EGFP.FB7) <- "stim"
DefaultAssay(TwoMvWT.combined.EGFP.FB7) <- "RNA"
Combined.EGFPFB7.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.EGFP.FB7, ident.1 = "ARQ9_2M", ident.2 = "WT_P35", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.EGFPFB7.TwoMvWT.0.Markers, "Combined.EGFPFB7.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.EGFPFB7.TwoMvWT <- read.csv("Combined.EGFPFB7.TwoMvWT.0.Markers.csv") 
Combined.EGFPFB7.TwoMvWT_pvalue <- Combined.EGFPFB7.TwoMvWT$p_val
Combined.EGFPFB7.TwoMvWT_pvalue=as.numeric(Combined.EGFPFB7.TwoMvWT_pvalue)
Combined.EGFPFB7.TwoMvWT_BH = p.adjust(Combined.EGFPFB7.TwoMvWT_pvalue, "BH")
write.csv(Combined.EGFPFB7.TwoMvWT_BH, "Combined.EGFPFB7.TwoMvWT_BH.csv")

#ARQ expression

Idents(object = TwoMvWT.combined.EGFP.FB1) <- "stim"
TwoM.combined.EGFP.FB1 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("ARQ9_2M"))
WT.combined.EGFP.FB1 <- subset(TwoMvWT.combined.EGFP.FB1, idents = c("WT_P35"))

DefaultAssay(TwoM.combined.EGFP.FB1) <- "RNA"
tiff(file = "TwoM.combined.EGFP.FB1 ARQ expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.EGFP.FB1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
DefaultAssay(WT.combined.EGFP.FB1) <- "RNA"
tiff(file = "WT.combined.EGFP.FB1 ARQ expression stim plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.EGFP.FB1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = TwoMvWT.combined.EGFP.FB11) <- "stim"
tiff(file = "TwoMvWT.combined.EGFP.FB11 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFP.FB11, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.EGFP.FB11 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFP.FB11, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = TwoMvWT.combined.EGFP.FB3) <- "stim"
tiff(file = "TwoMvWT.combined.EGFP.FB3 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFP.FB3, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.EGFP.FB3 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFP.FB3, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = TwoMvWT.combined.EGFP.FB5) <- "stim"
tiff(file = "TwoMvWT.combined.EGFP.FB5 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFP.FB5, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.EGFP.FB5 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 6, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFP.FB5, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####Sub-clustering Epi-2####

Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.Epi <- subset(TwoMvWT.combined1, idents = c("BE", "LE"))
Idents(object = TwoMvWT.combined.Epi) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.Epi <- ScaleData(TwoMvWT.combined.Epi, verbose = FALSE)
TwoMvWT.combined.Epi <- RunPCA(TwoMvWT.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.Epi, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.Epi <- FindNeighbors(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Epi <- FindClusters(TwoMvWT.combined.Epi, resolution = 0.5)
TwoMvWT.combined.Epi <- RunTSNE(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Epi <- RunUMAP(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvWT.combined.Epi UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.Epi) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.Epi)
TwoMvWT.combined.Epi <- ScaleData(TwoMvWT.combined.Epi, features = all.genes)
TwoMvWT.combined.Epi <- CellCycleScoring(TwoMvWT.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.Epi) <- "Phase"
DimPlot(TwoMvWT.combined.Epi, reduction = "umap")
tiff(file = "TwoMvWT.combined.Epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.Epi1 <- TwoMvWT.combined.Epi
DefaultAssay(TwoMvWT.combined.Epi1) <- "integrated"
TwoMvWT.combined.Epi1 <- ScaleData(TwoMvWT.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.Epi1))
TwoMvWT.combined.Epi1 <- RunPCA(TwoMvWT.combined.Epi1, features = VariableFeatures(TwoMvWT.combined.Epi1))
ElbowPlot(TwoMvWT.combined.Epi1, ndims = 50)

TwoMvWT.combined.Epi1 <- FindNeighbors(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:25)
TwoMvWT.combined.Epi1 <- FindClusters(TwoMvWT.combined.Epi1, resolution = 0.5)
TwoMvWT.combined.Epi1 <- RunUMAP(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:25)
TwoMvWT.combined.Epi1 <- RunTSNE(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:25)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWT.combined.Epi1) <- "Phase"
tiff(file = "TwoMvWT.combined.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
TwoMvWT.combined.Epi1 <- ScaleData(TwoMvWT.combined.Epi1, features = rownames(TwoMvWT.combined.Epi1))
TwoMvWT.combined.Epi1.markers <- FindAllMarkers(TwoMvWT.combined.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoMvWT.combined.Epi1.markers, "TwoMvWT.combined.Epi1.seurat.markers.csv")

DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
tiff(file = "TwoMvWT.combined.Epi1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Epi1, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()


####P60 WT####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P60")

WT.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191214_32848_33580_test/counts/33580_WT/outs/filtered_feature_bc_matrix")
WTunfiltered <- CreateSeuratObject(counts = WT.data,  min.cells = 3, min.features = 200, project = "WT_P60")
WTunfiltered <- NormalizeData(WTunfiltered)

WTunfiltered[["percent.mt"]] <- PercentageFeatureSet(WTunfiltered, pattern = "^mt-")

tiff(file = "WTunfiltered Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WTunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WTunfiltered Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "WTunfiltered Pre-filteration")
dev.off()
tiff(file = "WTunfiltered Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "WTunfiltered Pre-filteration")
dev.off()

plot1 <- FeatureScatter(WTunfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WTunfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

WT <- subset(WTunfiltered, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 15)

tiff(file = "WT Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "WT Post-filteration")
dev.off()
tiff(file = "WT Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "WT Post-filteration")
dev.off()

plot1 <- FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 2500)
tiff(file = "WT Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT)
dev.off()

#Clustering
all.genes <- rownames(WT)
WT <- ScaleData(WT, features = all.genes)
WT <- RunPCA(WT, features = VariableFeatures(object = WT))
print(WT[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WT, dims = 1:2, reduction = "pca")
DimPlot(WT, reduction = "pca")

tiff(file = "WT ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(WT, ndims = 50)
dev.off()

WT <- FindNeighbors(WT, dims = 1:24)
WT <- FindClusters(WT, resolution = 0.5)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:24)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:24)

Idents(object = WT) <- "seurat_clusters"
tiff(file = "WT UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(WT) <- "RNA"
all.genes <- rownames(WT)
WT <- ScaleData(WT, features = all.genes)
WT <- CellCycleScoring(WT, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = WT) <- "Phase"
tiff(file = "WT Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
WT1 <- WT
DefaultAssay(WT1) <- "RNA"
WT1 <- ScaleData(WT1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(WT1))
WT1 <- RunPCA(WT1, features = VariableFeatures(WT1))
ElbowPlot(WT1, ndims = 30)

WT1 <- FindNeighbors(WT1, reduction = "pca", dims = 1:20)
WT1 <- FindClusters(WT1, resolution = 0.5)
WT1 <- RunUMAP(WT1, reduction = "pca", dims = 1:20)
WT1 <- RunTSNE(WT1, reduction = "pca", dims = 1:20)
DimPlot(WT1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = WT1) <- "Phase"
tiff(file = "WT1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = WT1) <- "seurat_clusters"
tiff(file = "WT1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT1, reduction = "umap", pt.size = 0.3)
dev.off()

#### Merging Datasets TwoMvWT ####

#Stash old idents
TwoM1[["orig.clusters"]] <- Idents(object = TwoM1)
WT1[["orig.clusters"]] <- Idents(object = WT1)

#Set Current idents
Idents(object = TwoM1) <- "seurat_clusters"
Idents(object = WT1) <- "seurat_clusters"
TwoM1$stim <- "ARQ9_2M"
WT1$stim <- "WT_P60"
TwoMvWT.anchors <- FindIntegrationAnchors(object.list = list(TwoM1, WT1), dims = 1:20)
TwoMvWT.combined <- IntegrateData(anchorset = TwoMvWT.anchors, dims = 1:20)

DefaultAssay(TwoMvWT.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined <- ScaleData(TwoMvWT.combined, verbose = FALSE)
TwoMvWT.combined <- RunPCA(TwoMvWT.combined, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined, ndims = 50)
#Umap and Clustering
TwoMvWT.combined <- FindNeighbors(TwoMvWT.combined, reduction = "pca", dims = 1:20)
TwoMvWT.combined <- FindClusters(TwoMvWT.combined, resolution = 0.5)
TwoMvWT.combined <- RunTSNE(TwoMvWT.combined, reduction = "pca", dims = 1:20)
TwoMvWT.combined <- RunUMAP(TwoMvWT.combined, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvWT.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(TwoMvWT.combined) <- "RNA"
all.genes <- rownames(TwoMvWT.combined)
TwoMvWT.combined <- ScaleData(TwoMvWT.combined, features = all.genes)
TwoMvWT.combined <- CellCycleScoring(TwoMvWT.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined) <- "Phase"
DimPlot(TwoMvWT.combined, reduction = "umap")
tiff(file = "TwoMvWT.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined1 <- TwoMvWT.combined
DefaultAssay(TwoMvWT.combined1) <- "integrated"
TwoMvWT.combined1 <- ScaleData(TwoMvWT.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined1))
TwoMvWT.combined1 <- RunPCA(TwoMvWT.combined1, features = VariableFeatures(TwoMvWT.combined1))
ElbowPlot(TwoMvWT.combined1, ndims = 50)

TwoMvWT.combined1 <- FindNeighbors(TwoMvWT.combined1, reduction = "pca", dims = 1:25)
TwoMvWT.combined1 <- FindClusters(TwoMvWT.combined1, resolution = 0.5)
TwoMvWT.combined1 <- RunUMAP(TwoMvWT.combined1, reduction = "pca", dims = 1:25)
TwoMvWT.combined1 <- RunTSNE(TwoMvWT.combined1, reduction = "pca", dims = 1:25)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWT.combined1) <- "Phase"
tiff(file = "TwoMvWT.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvWT.combined1) <- "RNA"
Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
TwoMvWT.combined1 <- ScaleData(TwoMvWT.combined1, features = rownames(TwoMvWT.combined1))
TwoMvWT.combined1.markers <- FindAllMarkers(TwoMvWT.combined1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoMvWT.combined1.markers, "TwoMvWT.combined1.seurat.markers.csv")

DefaultAssay(TwoMvWT.combined1) <- "RNA"
tiff(file = "TwoMvWT.combined1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2",
                                                                "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                                "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
TwoMvWT.combined1 <- RenameIdents(object = TwoMvWT.combined1, '0' = "BE", '22' = "BE", '8' = "BE", '2' = "LE", '7' = "LE", '11' = "LE", '5' = "LE", '16' = "LE", '3' = "LE",
                                  '9' = "LE", '10' = "LE", '6' = "SV", '19' = "SV", '18' = "FB", '1' = "FB", '4' = "FB",
                                  '14' = "SM", '13' = "VE", '17' = "Pericyte", '20' = "Glia", '12' = "Leu", '21' = "Leu", '15' = "Lym")
TwoMvWT.combined1[["CellType"]] <- Idents(object = TwoMvWT.combined1)

#Umap
Idents(object = TwoMvWT.combined1) <- "CellType"
tiff(file = "TwoMvWT.combined1 CellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined1 CellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined1) <- "stim"
TwoMvWT.combined1$stim.CellType <- paste(Idents(TwoMvWT.combined1), TwoMvWT.combined1$CellType, sep = "_")
Idents(object = TwoMvWT.combined1) <- "stim.CellType"
table(Idents(TwoMvWT.combined1))

#PIN markers for SV
Idents(object = TwoMvWT.combined1) <- "stim"
TwoM.combined1 <- subset(TwoMvWT.combined1, idents = c("ARQ9_2M"))

DefaultAssay(TwoM.combined1) <- "RNA"
tiff(file = "TwoM.combined1 PIN marker Svs5 expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined1, reduction = "umap", features = c("Svs5"), cols = c("light grey", "red"), pt.size = 0.3)
dev.off()
tiff(file = "TwoM.combined1 PIN marker Slc31a2 expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined1, reduction = "umap", features = c("Slc31a2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q95")
dev.off()

####Re-clustering Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P60/Epi")


Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.Epi <- subset(TwoMvWT.combined1, idents = c("BE", "LE", "SV"))
Idents(object = TwoMvWT.combined.Epi) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.Epi <- ScaleData(TwoMvWT.combined.Epi, verbose = FALSE)
TwoMvWT.combined.Epi <- RunPCA(TwoMvWT.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.Epi, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.Epi <- FindNeighbors(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:21)
TwoMvWT.combined.Epi <- FindClusters(TwoMvWT.combined.Epi, resolution = 0.5)
TwoMvWT.combined.Epi <- RunTSNE(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:21)
TwoMvWT.combined.Epi <- RunUMAP(TwoMvWT.combined.Epi, reduction = "pca", dims = 1:21)

tiff(file = "TwoMvWT.combined.Epi UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.Epi) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.Epi)
TwoMvWT.combined.Epi <- ScaleData(TwoMvWT.combined.Epi, features = all.genes)
TwoMvWT.combined.Epi <- CellCycleScoring(TwoMvWT.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.Epi) <- "Phase"
DimPlot(TwoMvWT.combined.Epi, reduction = "umap")
tiff(file = "TwoMvWT.combined.Epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.Epi1 <- TwoMvWT.combined.Epi
DefaultAssay(TwoMvWT.combined.Epi1) <- "integrated"
TwoMvWT.combined.Epi1 <- ScaleData(TwoMvWT.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.Epi1))
TwoMvWT.combined.Epi1 <- RunPCA(TwoMvWT.combined.Epi1, features = VariableFeatures(TwoMvWT.combined.Epi1))
ElbowPlot(TwoMvWT.combined.Epi1, ndims = 50)

TwoMvWT.combined.Epi1 <- FindNeighbors(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Epi1 <- FindClusters(TwoMvWT.combined.Epi1, resolution = 0.5)
TwoMvWT.combined.Epi1 <- RunUMAP(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.Epi1 <- RunTSNE(TwoMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWT.combined.Epi1) <- "Phase"
tiff(file = "TwoMvWT.combined.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
TwoMvWT.combined.Epi1 <- ScaleData(TwoMvWT.combined.Epi1, features = rownames(TwoMvWT.combined.Epi1))
TwoMvWT.combined.Epi1.markers <- FindAllMarkers(TwoMvWT.combined.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoMvWT.combined.Epi1.markers, "TwoMvWT.combined.Epi1.seurat.markers.csv")

DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
tiff(file = "TwoMvWT.combined.Epi1 Epicelltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Epi1, reduction = "umap", features = c("Krt5", "Krt19", "Fbln1", "Tagln",
                                                                    "Gsdma", "Msmb", "Lrrc26", "Krt4",
                                                                    "Clu", "Svs5", "Ly6a", "Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
tiff(file = "TwoMvWT.combined.Epi1 Epicelltype marker expression plots-1.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Epi1, reduction = "umap", features = c("Krt5", "Krt14", "Trp63", "Krt4", "Ppp1r1b", "Krt7", "Krt19", "Tgm4", "Msmb", "Svs2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
TwoMvWT.combined.Epi1 <- RenameIdents(object = TwoMvWT.combined.Epi1, '0' = "BE1", '8' = "BE2", '11' = "BE3", '17' = "BE4", '14' = "BE5", '12' = "UrLE", '13' = "Prog1", '16' = "Prog2", 
                                      '7' = "VP1",
                                  '1' = "VP2", '15' = "VP3", '6' = "LP", '3' = "AP1", '4' = "AP2", '18' = "AP3", '9' = "AP4",
                                  '2' = "AP5", '5' = "SV", '10' = "OE")
TwoMvWT.combined.Epi1[["EpiCellType"]] <- Idents(object = TwoMvWT.combined.Epi1)

#Umap
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
tiff(file = "TwoMvWT.combined.Epi1 EpiCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined.Epi1 EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.Epi1) <- "stim"
TwoMvWT.combined.Epi1$stim.EpiCellType <- paste(Idents(TwoMvWT.combined.Epi1), TwoMvWT.combined.Epi1$EpiCellType, sep = "_")
Idents(object = TwoMvWT.combined.Epi1) <- "stim.EpiCellType"
table(Idents(TwoMvWT.combined.Epi1))


####Re-clustering FBSM####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P60/FBSM")

Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.FBSM <- subset(TwoMvWT.combined1, idents = c("FB", "SM"))
Idents(object = TwoMvWT.combined.FBSM) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.FBSM) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.FBSM <- ScaleData(TwoMvWT.combined.FBSM, verbose = FALSE)
TwoMvWT.combined.FBSM <- RunPCA(TwoMvWT.combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.FBSM, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.FBSM <- FindNeighbors(TwoMvWT.combined.FBSM, reduction = "pca", dims = 1:21)
TwoMvWT.combined.FBSM <- FindClusters(TwoMvWT.combined.FBSM, resolution = 0.5)
TwoMvWT.combined.FBSM <- RunTSNE(TwoMvWT.combined.FBSM, reduction = "pca", dims = 1:21)
TwoMvWT.combined.FBSM <- RunUMAP(TwoMvWT.combined.FBSM, reduction = "pca", dims = 1:21)

tiff(file = "TwoMvWT.combined.FBSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.FBSM) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.FBSM)
TwoMvWT.combined.FBSM <- ScaleData(TwoMvWT.combined.FBSM, features = all.genes)
TwoMvWT.combined.FBSM <- CellCycleScoring(TwoMvWT.combined.FBSM, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.FBSM) <- "Phase"
DimPlot(TwoMvWT.combined.FBSM, reduction = "umap")
tiff(file = "TwoMvWT.combined.FBSM Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.FBSM1 <- TwoMvWT.combined.FBSM
DefaultAssay(TwoMvWT.combined.FBSM1) <- "integrated"
TwoMvWT.combined.FBSM1 <- ScaleData(TwoMvWT.combined.FBSM1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.FBSM1))
TwoMvWT.combined.FBSM1 <- RunPCA(TwoMvWT.combined.FBSM1, features = VariableFeatures(TwoMvWT.combined.FBSM1))
ElbowPlot(TwoMvWT.combined.FBSM1, ndims = 50)

TwoMvWT.combined.FBSM1 <- FindNeighbors(TwoMvWT.combined.FBSM1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FBSM1 <- FindClusters(TwoMvWT.combined.FBSM1, resolution = 0.5)
TwoMvWT.combined.FBSM1 <- RunUMAP(TwoMvWT.combined.FBSM1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FBSM1 <- RunTSNE(TwoMvWT.combined.FBSM1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWT.combined.FBSM1) <- "Phase"
tiff(file = "TwoMvWT.combined.FBSM1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.FBSM1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FBSM1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.FBSM1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FBSM1 seurat UMAP LABEL.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = TwoMvWT.combined.FBSM1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FBSM1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvWT.combined.FBSM1) <- "RNA"
Idents(object = TwoMvWT.combined.FBSM1) <- "seurat_clusters"
TwoMvWT.combined.FBSM1 <- ScaleData(TwoMvWT.combined.FBSM1, features = rownames(TwoMvWT.combined.FBSM1))
TwoMvWT.combined.FBSM1.markers <- FindAllMarkers(TwoMvWT.combined.FBSM1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoMvWT.combined.FBSM1.markers, "TwoMvWT.combined.FBSM1.seurat.markers.csv")

#Dotplot
Idents(object = TwoMvWT.combined.FBSM1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FBSM1 seurat DotPlot.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(TwoMvWT.combined.FBSM1, features = c("Myh11", "Acta2", "Bmp5", "Pitx1", "Pitx2", "Prdm6", "Itga4", "Angptl1", "Jph2", "Cyp26a1", "Foxp2", "Slc36a2", "Crabp1", "Robo2", "Clstn2", "Mkx", "Efemp1", "Tfap2a", "Crlf1", "Bmp7", "Foxf1", "Foxf2", "Hand2", "Hoxd13", "Fgfr2", "Wif1", "Cxcl14", "Mfap5", "Col5a3", "Anpep", "Clec3b", "Col14a1", "Emcn", "Mef2c", "Fxyd7", "Fam150b", "Inhba", "Thy1", "Tac2", "Dach2", "Gsta3", "Ramp1", "Ar", "EGFP", "ARQ"), cols = c("white", "red")) + RotatedAxis()
dev.off()

#FeaturePlot
Idents(object = TwoMvWT.combined.FBSM1) <- "stim"
TwoM.combined.FBSM1 <- subset(TwoMvWT.combined.FBSM1, idents = c("ARQ9_2M"))
WT.combined.FBSM1 <- subset(TwoMvWT.combined.FBSM1, idents = c("WT_P60"))

DefaultAssay(TwoM.combined.FBSM1) <- "RNA"
tiff(file = "TwoM.combined.FBSM1 ARQ marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM1 EGFP marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM1 Ar marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FBSM1 Gli1 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

DefaultAssay(WT.combined.FBSM1) <- "RNA"
tiff(file = "WT.combined.FBSM1 ARQ marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FBSM1 EGFP marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FBSM1 Ar marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FBSM1 Gli1 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FBSM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.FBSM1) <- "stim"
TwoMvWT.combined.FBSM1$stim.seurat_clusters <- paste(Idents(TwoMvWT.combined.FBSM1), TwoMvWT.combined.FBSM1$seurat_clusters, sep = "_")
Idents(object = TwoMvWT.combined.FBSM1) <- "stim.seurat_clusters"
table(Idents(TwoMvWT.combined.FBSM1))

#

#Rename
Idents(object = TwoMvWT.combined.FBSM1) <- "seurat_clusters"
TwoMvWT.combined.FBSM1 <- RenameIdents(object = TwoMvWT.combined.FBSM1, '0' = "BE1", '8' = "BE2", '11' = "BE3", '17' = "BE4", '14' = "BE5", '12' = "UrLE", '13' = "Prog1", '16' = "Prog2", 
                                      '7' = "VP1",
                                      '1' = "VP2", '15' = "VP3", '6' = "LP", '3' = "AP1", '4' = "AP2", '18' = "AP3", '9' = "AP4",
                                      '2' = "AP5", '5' = "SV", '10' = "OE")
TwoMvWT.combined.FBSM1[["FBSMCellType"]] <- Idents(object = TwoMvWT.combined.FBSM1)

#Umap
Idents(object = TwoMvWT.combined.FBSM1) <- "EpiCellType"
tiff(file = "TwoMvWT.combined.FBSM1 EpiCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined.FBSM1 EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FBSM1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.FBSM1) <- "stim"
TwoMvWT.combined.FBSM1$stim.FBSMCellType <- paste(Idents(TwoMvWT.combined.FBSM1), TwoMvWT.combined.FBSM1$FBSMCellType, sep = "_")
Idents(object = TwoMvWT.combined.FBSM1) <- "stim.EpiCellType"
table(Idents(TwoMvWT.combined.FBSM1))

#Heatmap
DefaultAssay(onlyBE2) <- "RNA"
onlyBE2 <- ScaleData(onlyBE2, features = rownames(onlyBE2))
onlyBE2.markers <- FindAllMarkers(onlyBE2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
onlyBE2Top50 <- onlyBE2.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "onlyBE2 Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(onlyBE2, features = c(onlyBE2Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()



####Re-clustering FB####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P60/FB")

Idents(object = TwoMvWT.combined1) <- "CellType"
TwoMvWT.combined.FB <- subset(TwoMvWT.combined1, idents = c("FB"))
Idents(object = TwoMvWT.combined.FB) <- "seurat_clusters"
DefaultAssay(TwoMvWT.combined.FB) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, verbose = FALSE)
TwoMvWT.combined.FB <- RunPCA(TwoMvWT.combined.FB, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined.FB, ndims = 50)

#Umap and Clustering
TwoMvWT.combined.FB <- FindNeighbors(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB <- FindClusters(TwoMvWT.combined.FB, resolution = 0.5)
TwoMvWT.combined.FB <- RunTSNE(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB <- RunUMAP(TwoMvWT.combined.FB, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvWT.combined.FB UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
all.genes <- rownames(TwoMvWT.combined.FB)
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = all.genes)
TwoMvWT.combined.FB <- CellCycleScoring(TwoMvWT.combined.FB, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvWT.combined.FB) <- "Phase"
DimPlot(TwoMvWT.combined.FB, reduction = "umap")
tiff(file = "TwoMvWT.combined.FB Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvWT.combined.FB1 <- TwoMvWT.combined.FB
DefaultAssay(TwoMvWT.combined.FB1) <- "integrated"
TwoMvWT.combined.FB1 <- ScaleData(TwoMvWT.combined.FB1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvWT.combined.FB1))
TwoMvWT.combined.FB1 <- RunPCA(TwoMvWT.combined.FB1, features = VariableFeatures(TwoMvWT.combined.FB1))
ElbowPlot(TwoMvWT.combined.FB1, ndims = 50)

TwoMvWT.combined.FB1 <- FindNeighbors(TwoMvWT.combined.FB1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB1 <- FindClusters(TwoMvWT.combined.FB1, resolution = 0.4)
TwoMvWT.combined.FB1 <- RunUMAP(TwoMvWT.combined.FB1, reduction = "pca", dims = 1:20)
TwoMvWT.combined.FB1 <- RunTSNE(TwoMvWT.combined.FB1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoMvWT.combined.FB1) <- "Phase"
tiff(file = "TwoMvWT.combined.FB1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FB1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FB1 seurat UMAP LABEL.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.FB1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Rename
Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
TwoMvWT.combined.FB1 <- RenameIdents(object = TwoMvWT.combined.FB1, '1' = "FB1", '2' = "FB2", '0' = "FB3", '3' = "FB4", '5' = "FB5", '6' = "FB6", '7' = "FB7", '4' = "FB8")
TwoMvWT.combined.FB1[["FBCellType"]] <- Idents(object = TwoMvWT.combined.FB1)

#Umap
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
tiff(file = "TwoMvWT.combined.FB1 FBCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined.FB1 EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.FB1) <- "stim"
TwoMvWT.combined.FB1$stim.FBCellType <- paste(Idents(TwoMvWT.combined.FB1), TwoMvWT.combined.FB1$FBCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB1) <- "stim.FBCellType"
table(Idents(TwoMvWT.combined.FB1))

#Add ARQ info
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
TwoMvWT.combined.FB1ARQNeg <- subset(x=TwoMvWT.combined.FB1, subset = ARQ == 0)
TwoMvWT.combined.FB1ARQPos <- subset(x=TwoMvWT.combined.FB1, subset = ARQ > 0)
Idents(object = TwoMvWT.combined.FB1ARQNeg) <- "ARQNeg"
Idents(object = TwoMvWT.combined.FB1ARQPos) <- "ARQPos"
TwoMvWT.combined.FB1ARQNeg[["ARQExp"]] <- Idents(object = TwoMvWT.combined.FB1ARQNeg)
TwoMvWT.combined.FB1ARQPos[["ARQExp"]] <- Idents(object = TwoMvWT.combined.FB1ARQPos)
TwoMvWT.combined.FB1ARQ <- merge(x = TwoMvWT.combined.FB1ARQNeg, y = TwoMvWT.combined.FB1ARQPos)
Idents(object = TwoMvWT.combined.FB1ARQ) <- "ARQExp"
TwoMvWT.combined.FB1$ARQExp <- Idents(object = TwoMvWT.combined.FB1ARQ)

Idents(object = TwoMvWT.combined.FB1) <- "ARQExp"
TwoMvWT.combined.FB1$ARQExp.FBCellType <- paste(Idents(TwoMvWT.combined.FB1), TwoMvWT.combined.FB1$FBCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB1) <- "ARQExp.FBCellType"
TwoMvWT.combined.FB1$ARQExp.FBCellType.stim <- paste(Idents(TwoMvWT.combined.FB1), TwoMvWT.combined.FB1$stim, sep = "_")
Idents(object = TwoMvWT.combined.FB1) <- "ARQExp.FBCellType.stim"
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)

#Heatmap_All
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB1 <- ScaleData(TwoMvWT.combined.FB1, features = rownames(TwoMvWT.combined.FB1))
TwoMvWT.combined.FB1.markers <- FindAllMarkers(TwoMvWT.combined.FB1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FB1Top50 <- TwoMvWT.combined.FB1.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB1 Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB1, features = c(TwoMvWT.combined.FB1Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoMvWT.combined.FB1Top50, "TwoMvWT.combined.FB.Top50.csv")
write.csv(TwoMvWT.combined.FB1.markers, "TwoMvWT.combined.FB1.FBCellType.marker.csv")

#Heatmap_FB1-8_combined
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB1 <- ScaleData(TwoMvWT.combined.FB1, features = rownames(TwoMvWT.combined.FB1))

tiff(file = "TwoMvWT.combined.FB1 Heatmap CAF markers.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB1, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB1 Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB1, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_FB1-8_combined_stim
Idents(object = TwoMvWT.combined.FB1) <- "stim.FBCellType"
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)
TwoMvWT.combined.FB1 <- RenameIdents(object = TwoMvWT.combined.FB1, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5", 'ARQ9_2M_FB6' = "ARQ9_2M_FB6", 'ARQ9_2M_FB7' = "ARQ9_2M_FB7", 'ARQ9_2M_FB8' = "ARQ9_2M_FB8",
                                     'WT_P60_FB1' = "WT_P60_FB1", 'WT_P60_FB2' = "WT_P60_FB2", 'WT_P60_FB3' = "WT_P60_FB3", 'WT_P60_FB4' = "WT_P60_FB4", 'WT_P60_FB5' = "WT_P60_FB5", 'WT_P60_FB6' = "WT_P60_FB6", 'WT_P60_FB7' = "WT_P60_FB7", 'WT_P60_FB8' = "WT_P60_FB8")
TwoMvWT.combined.FB1[["stim.FBCellType"]] <- Idents(object = TwoMvWT.combined.FB1)
DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
TwoMvWT.combined.FB1 <- ScaleData(TwoMvWT.combined.FB1, features = rownames(TwoMvWT.combined.FB1))

tiff(file = "TwoMvWT.combined.FB1 stim Heatmap CAF markers.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB1, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB1 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 10, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB1, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_FB1-8_TwoM
Idents(object = TwoMvWT.combined.FB1) <- "stim"
TwoM.FB <- subset(TwoMvWT.combined.FB1, idents = c("ARQ9_2M"))
DefaultAssay(TwoM.FB) <- "RNA"
Idents(object = TwoM.FB) <- "FBCellType"
TwoM.FB <- ScaleData(TwoM.FB, features = rownames(TwoM.FB))
TwoM.FB.marker <- FindAllMarkers(TwoM.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(TwoM.FB.marker, "TwoM.FB.FBCellType.marker.csv")
tiff(file = "TwoM.FB Heatmap CAF markers.tiff", width = 10, height = 2, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.FB, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_FB5/6/2_combined_stim
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
DimPlot(TwoMvWT.combined.FB1, reduction = "umap", pt.size = 0.3)
TwoMvWT.combined.FB5 <- subset(TwoMvWT.combined.FB1, idents = c("FB5"))
TwoMvWT.combined.FB6 <- subset(TwoMvWT.combined.FB1, idents = c("FB6"))
TwoMvWT.combined.FB2 <- subset(TwoMvWT.combined.FB1, idents = c("FB2"))
TwoMvWT.combined.FB256 <- subset(TwoMvWT.combined.FB1, idents = c("FB2", "FB5", "FB6"))
TwoMvWT.combined.FB56 <- subset(TwoMvWT.combined.FB1, idents = c("FB5", "FB6"))

Idents(object = TwoMvWT.combined.FB5) <- "stim.FBCellType"
DefaultAssay(TwoMvWT.combined.FB5) <- "RNA"
TwoMvWT.combined.FB5 <- ScaleData(TwoMvWT.combined.FB5, features = rownames(TwoMvWT.combined.FB5))

tiff(file = "TwoMvWT.combined.FB5 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB5, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB5 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 7, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB5, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt4", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = TwoMvWT.combined.FB6) <- "stim.FBCellType"
DefaultAssay(TwoMvWT.combined.FB6) <- "RNA"
TwoMvWT.combined.FB6 <- ScaleData(TwoMvWT.combined.FB6, features = rownames(TwoMvWT.combined.FB6))

tiff(file = "TwoMvWT.combined.FB6 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB6, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB6 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 7, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB6, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt4", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = TwoMvWT.combined.FB2) <- "stim.FBCellType"
DefaultAssay(TwoMvWT.combined.FB2) <- "RNA"
TwoMvWT.combined.FB2 <- ScaleData(TwoMvWT.combined.FB2, features = rownames(TwoMvWT.combined.FB2))

tiff(file = "TwoMvWT.combined.FB2 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB2, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB2 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 7, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB2, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt4", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = TwoMvWT.combined.FB256) <- "stim.FBCellType"
DefaultAssay(TwoMvWT.combined.FB256) <- "RNA"
TwoMvWT.combined.FB256 <- ScaleData(TwoMvWT.combined.FB256, features = rownames(TwoMvWT.combined.FB256))

tiff(file = "TwoMvWT.combined.FB256 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB256, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB256 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 7, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB256, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt4", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

Idents(object = TwoMvWT.combined.FB56) <- "stim.FBCellType"
DefaultAssay(TwoMvWT.combined.FB56) <- "RNA"
TwoMvWT.combined.FB56 <- ScaleData(TwoMvWT.combined.FB56, features = rownames(TwoMvWT.combined.FB56))

tiff(file = "TwoMvWT.combined.FB56 stim Heatmap CAF markers.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB56, features = c("ARQ", "EGFP", "Ar", "Vim", "Pdgfrb", "Fap", "S100a4", "Twist1", "Foxf1", "Il11",  "Sox9", "Cxcl10")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
tiff(file = "TwoMvWT.combined.FB56 stim Heatmap AR WNT SP1 P53 IGF markers.tiff", width = 5, height = 7, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB56, features = c("ARQ", "EGFP", "Fkbp5", "Aldh1a1", "Ly6a", "Pbsn", "Wnt2", "Wnt4", "Wnt6", "Fos", "Fosb", "Jun", "Junb", "Egr1", "Serpine1", "Cdkn1a", "Dnajb1", "Bcl2", "Igfbp3", "Igfbp7", "Igf1")) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()
        
#FeaturePlot
Idents(object = TwoMvWT.combined.FB1) <- "stim"
TwoM.combined.FB1 <- subset(TwoMvWT.combined.FB1, idents = c("ARQ9_2M"))
WT.combined.FB1 <- subset(TwoMvWT.combined.FB1, idents = c("WT_P60"))

DefaultAssay(TwoM.combined.FB1) <- "RNA"
tiff(file = "TwoM.combined.FB1 ARQ marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FB1 EGFP marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FB1 Ar marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM.combined.FB1 Gli1 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

DefaultAssay(WT.combined.FB1) <- "RNA"
tiff(file = "WT.combined.FB1 ARQ marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FB1 EGFP marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "light grey"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FB1 Ar marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "WT.combined.FB1 Gli1 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(WT.combined.FB1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
tiff(file = "TwoMvWT.combined.FB1 Apod marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("Apod"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoMvWT.combined.FB1 C7 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("C7"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

tiff(file = "TwoMvWT.combined.FB1 Thy1 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("Thy1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoMvWT.combined.FB1 Prom1 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("Prom1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoMvWT.combined.FB1 Ly6a marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("Ly6a"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TwoMvWT.combined.FB1 Cd44 marker expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.FB1, reduction = "umap", features = c("Cd44"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

####Pseudotime using Slingslot####

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

DefaultAssay(TwoMvWT.combined.FB1) <- "RNA"
Idents(object = TwoMvWT.combined.FB1) <- "seurat_clusters"
DimPlot(TwoMvWT.combined.FB1, reduction = "umap")

FB_sds <- slingshot(Embeddings(TwoMvWT.combined.FB1, "umap"), clusterLabels = TwoMvWT.combined.FB1$seurat_clusters, 
                     start.clus = 1, stretch = 2)

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

cell_colors_clust <- cell_pal(TwoMvWT.combined.FB1$seurat_clusters, hue_pal())

tiff(file = "TwoMvWT.combined.FB1 slingshot Pseudotime.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot(reducedDim(FB_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(FB_sds, lwd = 2, type = 'lineages', col = 'black')
dev.off()

####DEGs in FB####

#DEGs_Combined.FB_TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "stim"
DimPlot(TwoMvWT.combined.FB1, reduction = "umap")
Combined.FB.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB1, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB.TwoMvWT.0.Markers, "Combined.FB.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB.TwoMvWT <- read.csv("Combined.FB.TwoMvWT.0.Markers.csv") 
Combined.FB.TwoMvWT_pvalue <- Combined.FB.TwoMvWT$p_val
Combined.FB.TwoMvWT_pvalue=as.numeric(Combined.FB.TwoMvWT_pvalue)
Combined.FB.TwoMvWT_BH = p.adjust(Combined.FB.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB.TwoMvWT_BH, "Combined.FB.TwoMvWT_BH.csv")

#DEGs_TwoM.FB.ARQPosvARQNeg
Idents(object = TwoMvWT.combined.FB1) <- "stim"
TwoM.combined.FB <- subset(TwoMvWT.combined.FB1, idents = c("ARQ9_2M"))

Idents(object = TwoM.combined.FB) <- "ARQExp"
TwoM.FB.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB.ARQPosvARQNeg.0.Markers, "TwoM.FB.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB.ARQPosvARQNeg <- read.csv("TwoM.FB.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB.ARQPosvARQNeg_pvalue <- TwoM.FB.ARQPosvARQNeg$p_val
TwoM.FB.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB.ARQPosvARQNeg_pvalue)
TwoM.FB.ARQPosvARQNeg_BH = p.adjust(TwoM.FB.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB.ARQPosvARQNeg_BH, "TwoM.FB.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB1.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB1 <- subset(TwoM.combined.FB, idents = c("FB1"))

Idents(object = TwoM.combined.FB1) <- "ARQExp"
TwoM.FB1.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB1, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB1.ARQPosvARQNeg.0.Markers, "TwoM.FB1.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB1.ARQPosvARQNeg <- read.csv("TwoM.FB1.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB1.ARQPosvARQNeg_pvalue <- TwoM.FB1.ARQPosvARQNeg$p_val
TwoM.FB1.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB1.ARQPosvARQNeg_pvalue)
TwoM.FB1.ARQPosvARQNeg_BH = p.adjust(TwoM.FB1.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB1.ARQPosvARQNeg_BH, "TwoM.FB1.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB2.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB2 <- subset(TwoM.combined.FB, idents = c("FB2"))

Idents(object = TwoM.combined.FB2) <- "ARQExp"
TwoM.FB2.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB2, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB2.ARQPosvARQNeg.0.Markers, "TwoM.FB2.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB2.ARQPosvARQNeg <- read.csv("TwoM.FB2.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB2.ARQPosvARQNeg_pvalue <- TwoM.FB2.ARQPosvARQNeg$p_val
TwoM.FB2.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB2.ARQPosvARQNeg_pvalue)
TwoM.FB2.ARQPosvARQNeg_BH = p.adjust(TwoM.FB2.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB2.ARQPosvARQNeg_BH, "TwoM.FB2.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB3.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB3 <- subset(TwoM.combined.FB, idents = c("FB3"))

Idents(object = TwoM.combined.FB3) <- "ARQExp"
TwoM.FB3.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB3, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB3.ARQPosvARQNeg.0.Markers, "TwoM.FB3.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB3.ARQPosvARQNeg <- read.csv("TwoM.FB3.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB3.ARQPosvARQNeg_pvalue <- TwoM.FB3.ARQPosvARQNeg$p_val
TwoM.FB3.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB3.ARQPosvARQNeg_pvalue)
TwoM.FB3.ARQPosvARQNeg_BH = p.adjust(TwoM.FB3.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB3.ARQPosvARQNeg_BH, "TwoM.FB3.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB4.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB4 <- subset(TwoM.combined.FB, idents = c("FB4"))

Idents(object = TwoM.combined.FB4) <- "ARQExp"
TwoM.FB4.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB4, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB4.ARQPosvARQNeg.0.Markers, "TwoM.FB4.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB4.ARQPosvARQNeg <- read.csv("TwoM.FB4.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB4.ARQPosvARQNeg_pvalue <- TwoM.FB4.ARQPosvARQNeg$p_val
TwoM.FB4.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB4.ARQPosvARQNeg_pvalue)
TwoM.FB4.ARQPosvARQNeg_BH = p.adjust(TwoM.FB4.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB4.ARQPosvARQNeg_BH, "TwoM.FB4.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB5.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB5 <- subset(TwoM.combined.FB, idents = c("FB5"))

Idents(object = TwoM.combined.FB5) <- "ARQExp"
TwoM.FB5.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB5, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB5.ARQPosvARQNeg.0.Markers, "TwoM.FB5.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB5.ARQPosvARQNeg <- read.csv("TwoM.FB5.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB5.ARQPosvARQNeg_pvalue <- TwoM.FB5.ARQPosvARQNeg$p_val
TwoM.FB5.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB5.ARQPosvARQNeg_pvalue)
TwoM.FB5.ARQPosvARQNeg_BH = p.adjust(TwoM.FB5.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB5.ARQPosvARQNeg_BH, "TwoM.FB5.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB6.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB6 <- subset(TwoM.combined.FB, idents = c("FB6"))

Idents(object = TwoM.combined.FB6) <- "ARQExp"
TwoM.FB6.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB6, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB6.ARQPosvARQNeg.0.Markers, "TwoM.FB6.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB6.ARQPosvARQNeg <- read.csv("TwoM.FB6.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB6.ARQPosvARQNeg_pvalue <- TwoM.FB6.ARQPosvARQNeg$p_val
TwoM.FB6.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB6.ARQPosvARQNeg_pvalue)
TwoM.FB6.ARQPosvARQNeg_BH = p.adjust(TwoM.FB6.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB6.ARQPosvARQNeg_BH, "TwoM.FB6.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB7.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB7 <- subset(TwoM.combined.FB, idents = c("FB7"))

Idents(object = TwoM.combined.FB7) <- "ARQExp"
TwoM.FB7.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB7, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB7.ARQPosvARQNeg.0.Markers, "TwoM.FB7.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB7.ARQPosvARQNeg <- read.csv("TwoM.FB7.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB7.ARQPosvARQNeg_pvalue <- TwoM.FB7.ARQPosvARQNeg$p_val
TwoM.FB7.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB7.ARQPosvARQNeg_pvalue)
TwoM.FB7.ARQPosvARQNeg_BH = p.adjust(TwoM.FB7.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB7.ARQPosvARQNeg_BH, "TwoM.FB7.ARQPosvARQNeg_BH.csv")

#DEGs_TwoM.FB8.ARQPosvARQNeg
Idents(object = TwoM.combined.FB) <- "FBCellType"
TwoM.combined.FB8 <- subset(TwoM.combined.FB, idents = c("FB8"))

Idents(object = TwoM.combined.FB8) <- "ARQExp"
TwoM.FB8.ARQPosvARQNeg.0.Markers <- FindMarkers(TwoM.combined.FB8, ident.1 = "ARQPos", ident.2 = "ARQNeg", min.pct = 0, logfc.threshold = 0)
write.csv(TwoM.FB8.ARQPosvARQNeg.0.Markers, "TwoM.FB8.ARQPosvARQNeg.0.Markers.csv")

#p.adjust
TwoM.FB8.ARQPosvARQNeg <- read.csv("TwoM.FB8.ARQPosvARQNeg.0.Markers.csv") 
TwoM.FB8.ARQPosvARQNeg_pvalue <- TwoM.FB8.ARQPosvARQNeg$p_val
TwoM.FB8.ARQPosvARQNeg_pvalue=as.numeric(TwoM.FB8.ARQPosvARQNeg_pvalue)
TwoM.FB8.ARQPosvARQNeg_BH = p.adjust(TwoM.FB8.ARQPosvARQNeg_pvalue, "BH")
write.csv(TwoM.FB8.ARQPosvARQNeg_BH, "TwoM.FB8.ARQPosvARQNeg_BH.csv")



#DEGs_combined.FB1.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB11 <- subset(TwoMvWT.combined.FB1, idents = c("FB1"))

Idents(object = TwoMvWT.combined.FB11) <- "stim"
Combined.FB1.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB11, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB1.TwoMvWT.0.Markers, "Combined.FB1.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB1.TwoMvWT <- read.csv("Combined.FB1.TwoMvWT.0.Markers.csv") 
Combined.FB1.TwoMvWT_pvalue <- Combined.FB1.TwoMvWT$p_val
Combined.FB1.TwoMvWT_pvalue=as.numeric(Combined.FB1.TwoMvWT_pvalue)
Combined.FB1.TwoMvWT_BH = p.adjust(Combined.FB1.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB1.TwoMvWT_BH, "Combined.FB1.TwoMvWT_BH.csv")

#DEGs_combined.FB2.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB2 <- subset(TwoMvWT.combined.FB1, idents = c("FB2"))

Idents(object = TwoMvWT.combined.FB2) <- "stim"
Combined.FB2.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB2, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB2.TwoMvWT.0.Markers, "Combined.FB2.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB2.TwoMvWT <- read.csv("Combined.FB2.TwoMvWT.0.Markers.csv") 
Combined.FB2.TwoMvWT_pvalue <- Combined.FB2.TwoMvWT$p_val
Combined.FB2.TwoMvWT_pvalue=as.numeric(Combined.FB2.TwoMvWT_pvalue)
Combined.FB2.TwoMvWT_BH = p.adjust(Combined.FB2.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB2.TwoMvWT_BH, "Combined.FB2.TwoMvWT_BH.csv")

#DEGs_combined.FB3.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB3 <- subset(TwoMvWT.combined.FB1, idents = c("FB3"))

Idents(object = TwoMvWT.combined.FB3) <- "stim"
Combined.FB3.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB3, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB3.TwoMvWT.0.Markers, "Combined.FB3.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB3.TwoMvWT <- read.csv("Combined.FB3.TwoMvWT.0.Markers.csv") 
Combined.FB3.TwoMvWT_pvalue <- Combined.FB3.TwoMvWT$p_val
Combined.FB3.TwoMvWT_pvalue=as.numeric(Combined.FB3.TwoMvWT_pvalue)
Combined.FB3.TwoMvWT_BH = p.adjust(Combined.FB3.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB3.TwoMvWT_BH, "Combined.FB3.TwoMvWT_BH.csv")

#DEGs_combined.FB4.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB4 <- subset(TwoMvWT.combined.FB1, idents = c("FB4"))

Idents(object = TwoMvWT.combined.FB4) <- "stim"
Combined.FB4.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB4, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB4.TwoMvWT.0.Markers, "Combined.FB4.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB4.TwoMvWT <- read.csv("Combined.FB4.TwoMvWT.0.Markers.csv") 
Combined.FB4.TwoMvWT_pvalue <- Combined.FB4.TwoMvWT$p_val
Combined.FB4.TwoMvWT_pvalue=as.numeric(Combined.FB4.TwoMvWT_pvalue)
Combined.FB4.TwoMvWT_BH = p.adjust(Combined.FB4.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB4.TwoMvWT_BH, "Combined.FB4.TwoMvWT_BH.csv")

#DEGs_combined.FB5.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB5 <- subset(TwoMvWT.combined.FB1, idents = c("FB5"))

Idents(object = TwoMvWT.combined.FB5) <- "stim"
TwoMvWT.combined.FB5 <- ScaleData(TwoMvWT.combined.FB5, features = rownames(TwoMvWT.combined.FB5))
Combined.FB5.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB5, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB5.TwoMvWT.0.Markers, "Combined.FB5.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB5.TwoMvWT <- read.csv("Combined.FB5.TwoMvWT.0.Markers.csv") 
Combined.FB5.TwoMvWT_pvalue <- Combined.FB5.TwoMvWT$p_val
Combined.FB5.TwoMvWT_pvalue=as.numeric(Combined.FB5.TwoMvWT_pvalue)
Combined.FB5.TwoMvWT_BH = p.adjust(Combined.FB5.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB5.TwoMvWT_BH, "Combined.FB5.TwoMvWT_BH.csv")

#DEGs_combined.FB6.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB6 <- subset(TwoMvWT.combined.FB1, idents = c("FB6"))

Idents(object = TwoMvWT.combined.FB6) <- "stim"
Combined.FB6.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB6, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB6.TwoMvWT.0.Markers, "Combined.FB6.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB6.TwoMvWT <- read.csv("Combined.FB6.TwoMvWT.0.Markers.csv") 
Combined.FB6.TwoMvWT_pvalue <- Combined.FB6.TwoMvWT$p_val
Combined.FB6.TwoMvWT_pvalue=as.numeric(Combined.FB6.TwoMvWT_pvalue)
Combined.FB6.TwoMvWT_BH = p.adjust(Combined.FB6.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB6.TwoMvWT_BH, "Combined.FB6.TwoMvWT_BH.csv")

#DEGs_combined.FB7.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB7 <- subset(TwoMvWT.combined.FB1, idents = c("FB7"))

Idents(object = TwoMvWT.combined.FB7) <- "stim"
Combined.FB7.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB7, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB7.TwoMvWT.0.Markers, "Combined.FB7.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB7.TwoMvWT <- read.csv("Combined.FB7.TwoMvWT.0.Markers.csv") 
Combined.FB7.TwoMvWT_pvalue <- Combined.FB7.TwoMvWT$p_val
Combined.FB7.TwoMvWT_pvalue=as.numeric(Combined.FB7.TwoMvWT_pvalue)
Combined.FB7.TwoMvWT_BH = p.adjust(Combined.FB7.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB7.TwoMvWT_BH, "Combined.FB7.TwoMvWT_BH.csv")

#DEGs_combined.FB8.TwoMvWT
Idents(object = TwoMvWT.combined.FB1) <- "FBCellType"
TwoMvWT.combined.FB8 <- subset(TwoMvWT.combined.FB1, idents = c("FB8"))

Idents(object = TwoMvWT.combined.FB8) <- "stim"
Combined.FB8.TwoMvWT.0.Markers <- FindMarkers(TwoMvWT.combined.FB8, ident.1 = "ARQ9_2M", ident.2 = "WT_P60", min.pct = 0, logfc.threshold = 0)
write.csv(Combined.FB8.TwoMvWT.0.Markers, "Combined.FB8.TwoMvWT.0.Markers.csv")

#p.adjust
Combined.FB8.TwoMvWT <- read.csv("Combined.FB8.TwoMvWT.0.Markers.csv") 
Combined.FB8.TwoMvWT_pvalue <- Combined.FB8.TwoMvWT$p_val
Combined.FB8.TwoMvWT_pvalue=as.numeric(Combined.FB8.TwoMvWT_pvalue)
Combined.FB8.TwoMvWT_BH = p.adjust(Combined.FB8.TwoMvWT_pvalue, "BH")
write.csv(Combined.FB8.TwoMvWT_BH, "Combined.FB8.TwoMvWT_BH.csv")

####Cell Count####
Idents(object = TwoM.combined.FB) <- "FBCellType"

tiff(file = "TwoM.combined.FB FBCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM.combined.FB, reduction = "umap", pt.size = 1, label = TRUE,cols = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF"))
dev.off()

Idents(object = TwoM.combined.FB) <- "ARQExp"

tiff(file = "TwoM.combined.FB ARQExp split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM.combined.FB, reduction = "umap", pt.size = 1, split.by = "ARQExp", cols = c("skyblue", "darkorange"))
dev.off()

TwoM.combined.FB$ARQExp.FBCellType <- paste(Idents(TwoM.combined.FB), TwoM.combined.FB$FBCellType, sep = "_")
Idents(object = TwoM.combined.FB) <- "ARQExp.FBCellType"
table(Idents(TwoM.combined.FB))


tiff(file = "TwoM.combined.FB EGFP expression plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1)
dev.off()
tiff(file = "TwoM.combined.FB EGFP expression over0.5 plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, min.cutoff = 0.5)
dev.off()
tiff(file = "TwoM.combined.FB EGFP expression over1 plots.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM.combined.FB, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, min.cutoff = 1.0)
dev.off()

#Add EGFP info
DefaultAssay(TwoM.combined.FB) <- "RNA"
TwoM.combined.FBEGFPNeg <- subset(x=TwoM.combined.FB, subset = EGFP < 1)
TwoM.combined.FBEGFPPos <- subset(x=TwoM.combined.FB, subset = EGFP >= 1)
Idents(object = TwoM.combined.FBEGFPNeg) <- "EGFPNeg"
Idents(object = TwoM.combined.FBEGFPPos) <- "EGFPPos"
TwoM.combined.FBEGFPNeg[["EGFPExp"]] <- Idents(object = TwoM.combined.FBEGFPNeg)
TwoM.combined.FBEGFPPos[["EGFPExp"]] <- Idents(object = TwoM.combined.FBEGFPPos)
TwoM.combined.FBEGFP <- merge(x = TwoM.combined.FBEGFPNeg, y = TwoM.combined.FBEGFPPos)
Idents(object = TwoM.combined.FBEGFP) <- "EGFPExp"
TwoM.combined.FB$EGFPExp <- Idents(object = TwoM.combined.FBEGFP)

Idents(object = TwoM.combined.FB) <- "EGFPExp"
TwoM.combined.FB$EGFPExp.ARQExp <- paste(Idents(TwoM.combined.FB), TwoM.combined.FB$ARQExp, sep = "_")
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp"
table(Idents(TwoM.combined.FB))

TwoM.combined.FB$EGFPExp.ARQExp.FBCellType <- paste(Idents(TwoM.combined.FB), TwoM.combined.FB$FBCellType, sep = "_")
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp.FBCellType"
table(Idents(TwoM.combined.FB))

####6M####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/6M")

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

SixM.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/v3_count__ARQ9_length_3191/v3_count_len3191_Project_COHP_46247_2_X3SC4/outs/filtered_feature_bc_matrix")
SixM.unfiltered <- CreateSeuratObject(counts = SixM.data,  min.cells = 3, min.features = 500, project = "ARQ9-Gli1_6M")
SixM.unfiltered <- NormalizeData(SixM.unfiltered)

#Initial processing & filtering

SixM.unfiltered[["percent.mt"]] <- PercentageFeatureSet(SixM.unfiltered, pattern = "^mt-")

tiff(file = "SixM.unfiltered Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(SixM.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "SixM.unfiltered Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(SixM.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "SixM.unfiltered Pre-filteration")
dev.off()
tiff(file = "SixM.unfiltered Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(SixM.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "SixM.unfiltered Pre-filteration")
dev.off()

plot1 <- FeatureScatter(SixM.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SixM.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

SixM <- subset(SixM.unfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 & percent.mt < 15)

tiff(file = "SixM Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(SixM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "SixM Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(SixM@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "SixM Post-filteration")
dev.off()
tiff(file = "SixM Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(SixM@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "SixM Post-filteration")
dev.off()

plot1 <- FeatureScatter(SixM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(SixM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

SixM <- FindVariableFeatures(SixM, selection.method = "vst", nfeatures = 2500)
tiff(file = "SixM Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(SixM)
dev.off()

#Clustering
all.genes <- rownames(SixM)
SixM <- ScaleData(SixM, features = all.genes)
SixM <- RunPCA(SixM, features = VariableFeatures(object = SixM))
print(SixM[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(SixM, dims = 1:2, reduction = "pca")
DimPlot(SixM, reduction = "pca")

tiff(file = "SixM ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(SixM, ndims = 50)
dev.off()

SixM <- FindNeighbors(SixM, dims = 1:25)
SixM <- FindClusters(SixM, resolution = 0.5)
head(Idents(SixM), 5)
SixM <- RunTSNE(SixM, reduction = "pca", dims = 1:25)
SixM <- RunUMAP(SixM, reduction = "pca", dims = 1:25)
DimPlot(SixM, reduction = "umap", pt.size = 1, label = TRUE)

DefaultAssay(SixM) <- "RNA"
FeaturePlot(SixM, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(SixM, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(SixM, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(SixM, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)


tiff(file = "SixM hARtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM mGFP expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Ar expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Gli1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####Merging Datasets TwoMvSixM####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/6M/TwoMvSixM")

#Stash old idents
TwoM[["orig.clusters"]] <- Idents(object = TwoM)
SixM[["orig.clusters"]] <- Idents(object = SixM)

#Set Current idents
Idents(object = TwoM) <- "seurat_clusters"
Idents(object = SixM) <- "seurat_clusters"
tiff(file = "TwoM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "SixM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixM, reduction = "umap", pt.size = 0.3)
dev.off()

TwoM$stim <- "ARQ9_2M"
SixM$stim <- "ARQ9_6M"
TwoMvSixM.anchors <- FindIntegrationAnchors(object.list = list(TwoM, SixM), dims = 1:20)
TwoMvSixM.combined <- IntegrateData(anchorset = TwoMvSixM.anchors, dims = 1:20)

DefaultAssay(TwoMvSixM.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvSixM.combined <- ScaleData(TwoMvSixM.combined, verbose = FALSE)
TwoMvSixM.combined <- RunPCA(TwoMvSixM.combined, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvSixM.combined, ndims = 50)

#Umap and Clustering
TwoMvSixM.combined <- FindNeighbors(TwoMvSixM.combined, reduction = "pca", dims = 1:20)
TwoMvSixM.combined <- FindClusters(TwoMvSixM.combined, resolution = 0.5)
TwoMvSixM.combined <- RunTSNE(TwoMvSixM.combined, reduction = "pca", dims = 1:20)
TwoMvSixM.combined <- RunUMAP(TwoMvSixM.combined, reduction = "pca", dims = 1:20)

tiff(file = "TwoMvSixM.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvSixM.combined) <- "RNA"
all.genes <- rownames(TwoMvSixM.combined)
TwoMvSixM.combined <- ScaleData(TwoMvSixM.combined, features = all.genes)
TwoMvSixM.combined <- CellCycleScoring(TwoMvSixM.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvSixM.combined) <- "Phase"
DimPlot(TwoMvSixM.combined, reduction = "umap")
tiff(file = "TwoMvSixM.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvSixM.combined1 <- TwoMvSixM.combined
DefaultAssay(TwoMvSixM.combined1) <- "integrated"
TwoMvSixM.combined1 <- ScaleData(TwoMvSixM.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvSixM.combined1))
TwoMvSixM.combined1 <- RunPCA(TwoMvSixM.combined1, features = VariableFeatures(TwoMvSixM.combined1))
ElbowPlot(TwoMvSixM.combined1, ndims = 50)

TwoMvSixM.combined1 <- FindNeighbors(TwoMvSixM.combined1, reduction = "pca", dims = 1:20)
TwoMvSixM.combined1 <- FindClusters(TwoMvSixM.combined1, resolution = 0.5)
TwoMvSixM.combined1 <- RunUMAP(TwoMvSixM.combined1, reduction = "pca", dims = 1:20)
TwoMvSixM.combined1 <- RunTSNE(TwoMvSixM.combined1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvSixM.combined1) <- "Phase"
tiff(file = "TwoMvSixM.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = TwoMvSixM.combined1) <- "seurat_clusters"
tiff(file = "TwoMvSixM.combined1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = TwoMvSixM.combined1) <- "seurat_clusters"
tiff(file = "TwoMvSixM.combined1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(TwoMvSixM.combined1) <- "RNA"
tiff(file = "TwoMvSixM.combined1 celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvSixM.combined1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2",
                                                                "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                                "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

FeaturePlot(TwoMvSixM.combined1, reduction = "umap", features = c("Chga"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")


#Rename
Idents(object = TwoMvSixM.combined1) <- "seurat_clusters"
TwoMvSixM.combined1 <- RenameIdents(object = TwoMvSixM.combined1, '0' = "BE", '3' = "BE", '15' = "BE", '8' = "LE", '2' = "LE", '13' = "LE", '7' = "LE", '11' = "LE",
                                  '5' = "LE", '19' = "LE",  '6' = "SV", '23' = "SV", '1' = "FB", '18' = "FB", '16' = "FB", '4' = "FB",
                                  '12' = "SM", '9' = "VE", '22' = "VE", '14' = "Pericyte", '20' = "Glia", '21' = "Leu", '17' = "Leu", '10' = "Leu")
TwoMvSixM.combined1[["CellType"]] <- Idents(object = TwoMvSixM.combined1)

#Umap
Idents(object = TwoMvSixM.combined1) <- "CellType"
tiff(file = "TwoMvSixM.combined1 CellType LABEL UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()

tiff(file = "TwoMvSixM.combined1 CellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvSixM.combined1) <- "stim"
TwoMvSixM.combined1$stim.CellType <- paste(Idents(TwoMvSixM.combined1), TwoMvSixM.combined1$CellType, sep = "_")
Idents(object = TwoMvSixM.combined1) <- "stim.CellType"
table(Idents(TwoMvSixM.combined1))

####Reclustering Stro####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/6M/TwoMvSixM/FBSM")

Idents(object = TwoMvSixM.combined1) <- "CellType"
TwoMvSixM.combined.Stro <- subset(TwoMvSixM.combined1, idents = c("FB", "SM", "VE", "Pericyte", "Glia", "Leu"))
Idents(object = TwoMvSixM.combined.Stro) <- "seurat_clusters"
DefaultAssay(TwoMvSixM.combined.Stro) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvSixM.combined.Stro <- ScaleData(TwoMvSixM.combined.Stro, verbose = FALSE)
TwoMvSixM.combined.Stro <- RunPCA(TwoMvSixM.combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvSixM.combined.Stro, ndims = 50)

#Umap and Clustering
TwoMvSixM.combined.Stro <- FindNeighbors(TwoMvSixM.combined.Stro, reduction = "pca", dims = 1:16)
TwoMvSixM.combined.Stro <- FindClusters(TwoMvSixM.combined.Stro, resolution = 0.5)
TwoMvSixM.combined.Stro <- RunTSNE(TwoMvSixM.combined.Stro, reduction = "pca", dims = 1:16)
TwoMvSixM.combined.Stro <- RunUMAP(TwoMvSixM.combined.Stro, reduction = "pca", dims = 1:16)
DimPlot(TwoMvSixM.combined.Stro, reduction = "umap", pt.size = 0.3)

tiff(file = "TwoMvSixM.combined.Stro UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.Stro, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(TwoMvSixM.combined.Stro) <- "RNA"
all.genes <- rownames(TwoMvSixM.combined.Stro)
TwoMvSixM.combined.Stro <- ScaleData(TwoMvSixM.combined.Stro, features = all.genes)
TwoMvSixM.combined.Stro <- CellCycleScoring(TwoMvSixM.combined.Stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = TwoMvSixM.combined.Stro) <- "Phase"
DimPlot(TwoMvSixM.combined.Stro, reduction = "umap")
tiff(file = "TwoMvSixM.combined.Stro Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.Stro, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
TwoMvSixM.combined.Stro1 <- TwoMvSixM.combined.Stro
DefaultAssay(TwoMvSixM.combined.Stro1) <- "integrated"
TwoMvSixM.combined.Stro1 <- ScaleData(TwoMvSixM.combined.Stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(TwoMvSixM.combined.Stro1))
TwoMvSixM.combined.Stro1 <- RunPCA(TwoMvSixM.combined.Stro1, features = VariableFeatures(TwoMvSixM.combined.Stro1))
ElbowPlot(TwoMvSixM.combined.Stro1, ndims = 50)

TwoMvSixM.combined.Stro1 <- FindNeighbors(TwoMvSixM.combined.Stro1, reduction = "pca", dims = 1:20)
TwoMvSixM.combined.Stro1 <- FindClusters(TwoMvSixM.combined.Stro1, resolution = 0.5)
TwoMvSixM.combined.Stro1 <- RunUMAP(TwoMvSixM.combined.Stro1, reduction = "pca", dims = 1:20)
TwoMvSixM.combined.Stro1 <- RunTSNE(TwoMvSixM.combined.Stro1, reduction = "pca", dims = 1:20)
DimPlot(TwoMvSixM.combined.Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

DimPlot(TwoMvSixM.combined.Stro1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvSixM.combined.Stro1) <- "Phase"
tiff(file = "TwoMvSixM.combined.Stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.Stro1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
Idents(object = TwoMvSixM.combined.Stro1) <- "seurat_clusters"
tiff(file = "TwoMvSixM.combined.Stro1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

#Rename
DefaultAssay(TwoMvSixM.combined.Stro1) <- "RNA"
tiff(file = "TwoMvSixM.combined.Stro1 expression plots.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvSixM.combined.Stro1, reduction = "umap", features = c("Fbln1", "Myh11", "Pecam1", "Plp1", "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 1)
dev.off()

Idents(object = TwoMvSixM.combined.Stro1) <- "seurat_clusters"
TwoMvSixM.combined.Stro1 <- RenameIdents(object = TwoMvSixM.combined.Stro1, '1' = "FB1", '0' = "FB2", '5' = "FB3", '2' = "FB4", '8' = "FB5", 
                                       '7' = "SM", '3' = "VE1", '14' = "VE2", '10' = "Glia", '11' = "Pericyte1", '13' = "Pericyte2",
                                       '6' = "Leu1", '12' = "Leu2", '9' = "Leu3", '4' = "OF")
TwoMvSixM.combined.Stro1[["StroCellType"]] <- Idents(object = TwoMvSixM.combined.Stro1)

#Umap
Idents(object = TwoMvSixM.combined.Stro1) <- "StroCellType"
tiff(file = "TwoMvSixM.combined.Stro1 StroCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvSixM.combined.Stro1 StroCellType stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.Stro1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvSixM.combined.Stro1) <- "stim"
TwoMvSixM.combined.Stro1$stim.StroCellType <- paste(Idents(TwoMvSixM.combined.Stro1), TwoMvSixM.combined.Stro1$StroCellType, sep = "_")
Idents(object = TwoMvSixM.combined.Stro1) <- "stim.StroCellType"
table(Idents(TwoMvSixM.combined.Stro1))

####Subset FB####
Idents(object = TwoMvSixM.combined.Stro1) <- "StroCellType"

TwoMvSixM.combined.FBSM <- subset(TwoMvSixM.combined.Stro1, idents = c("FB1", "FB2", "FB3", "FB4", "FB5", "SM"))
TwoMvSixM.combined.FBSM <- RunTSNE(TwoMvSixM.combined.FBSM, reduction = "pca", dims = 1:20)
TwoMvSixM.combined.FBSM <- RunUMAP(TwoMvSixM.combined.FBSM, reduction = "pca", dims = 1:20)
DimPlot(TwoMvSixM.combined.FBSM, reduction = "umap", pt.size = 0.3)

tiff(file = "TwoMvSixM.combined.FBSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.FBSM, reduction = "umap", pt.size = 1)
dev.off()

tiff(file = "TwoMvSixM.combined.FBSM split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvSixM.combined.FBSM, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvSixM.combined.FBSM) <- "stim"
TwoMvSixM.combined.FBSM$stim.StroCellType <- paste(Idents(TwoMvSixM.combined.FBSM), TwoMvSixM.combined.FBSM$StroCellType, sep = "_")
Idents(object = TwoMvSixM.combined.FBSM) <- "stim.StroCellType"
table(Idents(TwoMvSixM.combined.FBSM))

#Add EGFP info
DefaultAssay(TwoMvSixM.combined.FBSM) <- "RNA"
TwoMvSixM.combined.FBSMEGFPNeg <- subset(x=TwoMvSixM.combined.FBSM, subset = EGFP == 0)
TwoMvSixM.combined.FBSMEGFPPos <- subset(x=TwoMvSixM.combined.FBSM, subset = EGFP > 0)
Idents(object = TwoMvSixM.combined.FBSMEGFPNeg) <- "EGFPNeg"
Idents(object = TwoMvSixM.combined.FBSMEGFPPos) <- "EGFPPos"
TwoMvSixM.combined.FBSMEGFPNeg[["EGFPExp"]] <- Idents(object = TwoMvSixM.combined.FBSMEGFPNeg)
TwoMvSixM.combined.FBSMEGFPPos[["EGFPExp"]] <- Idents(object = TwoMvSixM.combined.FBSMEGFPPos)
TwoMvSixM.combined.FBSMEGFP <- merge(x = TwoMvSixM.combined.FBSMEGFPNeg, y = TwoMvSixM.combined.FBSMEGFPPos)
Idents(object = TwoMvSixM.combined.FBSMEGFP) <- "EGFPExp"
TwoMvSixM.combined.FBSM$EGFPExp <- Idents(object = TwoMvSixM.combined.FBSMEGFP)

#Add ARQ info
DefaultAssay(TwoMvSixM.combined.FBSM) <- "RNA"
TwoMvSixM.combined.FBSMARQNeg <- subset(x=TwoMvSixM.combined.FBSM, subset = ARQ == 0)
TwoMvSixM.combined.FBSMARQPos <- subset(x=TwoMvSixM.combined.FBSM, subset = ARQ > 0)
Idents(object = TwoMvSixM.combined.FBSMARQNeg) <- "ARQNeg"
Idents(object = TwoMvSixM.combined.FBSMARQPos) <- "ARQPos"
TwoMvSixM.combined.FBSMARQNeg[["ARQExp"]] <- Idents(object = TwoMvSixM.combined.FBSMARQNeg)
TwoMvSixM.combined.FBSMARQPos[["ARQExp"]] <- Idents(object = TwoMvSixM.combined.FBSMARQPos)
TwoMvSixM.combined.FBSMARQ <- merge(x = TwoMvSixM.combined.FBSMARQNeg, y = TwoMvSixM.combined.FBSMARQPos)
Idents(object = TwoMvSixM.combined.FBSMARQ) <- "ARQExp"
TwoMvSixM.combined.FBSM$ARQExp <- Idents(object = TwoMvSixM.combined.FBSMARQ)

#Cellcounts
Idents(object = TwoMvSixM.combined.FBSM) <- "EGFPExp"
TwoMvSixM.combined.FBSM$EGFPExp.ARQExp <- paste(Idents(TwoMvSixM.combined.FBSM), TwoMvSixM.combined.FBSM$ARQExp, sep = "_")
Idents(object = TwoMvSixM.combined.FBSM) <- "EGFPExp.ARQExp"
TwoMvSixM.combined.FBSM$EGFPExp.ARQExp.StroCellType <- paste(Idents(TwoMvSixM.combined.FBSM), TwoMvSixM.combined.FBSM$StroCellType, sep = "_")
Idents(object = TwoMvSixM.combined.FBSM) <- "EGFPExp.ARQExp.StroCellType"
TwoMvSixM.combined.FBSM$EGFPExp.ARQExp.StroCellType.stim <- paste(Idents(TwoMvSixM.combined.FBSM), TwoMvSixM.combined.FBSM$stim, sep = "_")
Idents(object = TwoMvSixM.combined.FBSM) <- "EGFPExp.ARQExp.StroCellType.stim"
table(Idents(TwoMvSixM.combined.FBSM))


#Heatmap_AllFB
Idents(object = TwoMvSixM.combined.FBSM) <- "EGFPExp.ARQExp"
TwoMvSixM.combined.FBSM.EGFPARQ <- subset(x=TwoMvSixM.combined.FBSM, idents = c("EGFPPos_ARQPos"))


DefaultAssay(TwoMvSixM.combined.FBSM.EGFPARQ) <- "RNA"
Idents(object = TwoMvSixM.combined.FBSM.EGFPARQ) <- "StroCellType"
TwoMvSixM.combined.FBSM.EGFPARQ <- ScaleData(TwoMvSixM.combined.FBSM.EGFPARQ, features = rownames(TwoMvSixM.combined.FBSM.EGFPARQ))
TwoMvSixM.combined.FBSM.EGFPARQ.markers <- FindAllMarkers(TwoMvSixM.combined.FBSM.EGFPARQ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvSixM.combined.FBSM.EGFPARQTop50 <- TwoMvSixM.combined.FBSM.EGFPARQ.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

Idents(object = TwoMvSixM.combined.FBSM.EGFPARQ) <- "stim.StroCellType"
DimPlot(TwoMvSixM.combined.FBSM.EGFPARQ, reduction = "umap", pt.size = 0.3)

TwoMvSixM.combined.FBSM.EGFPARQ <- RenameIdents(object = TwoMvSixM.combined.FBSM.EGFPARQ, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5", 'ARQ9_2M_SM' = "ARQ9_2M_SM",
                                                'ARQ9_6M_FB1' = "ARQ9_6M_FB1", 'ARQ9_6M_FB2' = "ARQ9_6M_FB2", 'ARQ9_6M_FB3' = "ARQ9_6M_FB3", 'ARQ9_6M_FB4' = "ARQ9_6M_FB4", 'ARQ9_6M_FB5' = "ARQ9_6M_FB5", 'ARQ9_6M_SM' = "ARQ9_6M_SM")
TwoMvSixM.combined.FBSM.EGFPARQ[["stim.StroCellType"]] <- Idents(object = TwoMvSixM.combined.FBSM.EGFPARQ)
DimPlot(TwoMvSixM.combined.FBSM.EGFPARQ, reduction = "umap", pt.size = 0.3)

DefaultAssay(TwoMvSixM.combined.FBSM.EGFPARQ) <- "RNA"
TwoMvSixM.combined.FBSM.EGFPARQ <- ScaleData(TwoMvSixM.combined.FBSM.EGFPARQ, features = rownames(TwoMvSixM.combined.FBSM.EGFPARQ))
tiff(file = "TwoMvSixM.combined.FBSM.EGFPARQ stim Heatmap Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvSixM.combined.FBSM.EGFPARQ, features = c("ARQ", "EGFP", TwoMvSixM.combined.FBSM.EGFPARQTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))+ theme(axis.text.y = element_text(size = 5))
dev.off()

#Heatmap_AllFBSM
Idents(object = TwoMvSixM.combined.FBSM.EGFPARQ) <- "stim.StroCellType"

DefaultAssay(TwoMvSixM.combined.FBSM.EGFPARQ) <- "RNA"
TwoMvSixM.combined.FBSM.EGFPARQ <- ScaleData(TwoMvSixM.combined.FBSM.EGFPARQ, features = rownames(TwoMvSixM.combined.FBSM.EGFPARQ))
TwoMvSixM.combined.FBSM.EGFPARQ.stim.markers <- FindAllMarkers(TwoMvSixM.combined.FBSM.EGFPARQ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvSixM.combined.FBSM.EGFPARQ.stim.Top50 <- TwoMvSixM.combined.FBSM.EGFPARQ.stim.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

tiff(file = "TwoMvSixM.combined.FBSM.EGFPARQ.stim.Top50 heatmap purple.tiff", width = 12, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvSixM.combined.FBSM.EGFPARQ, features = c("ARQ", "EGFP", TwoMvSixM.combined.FBSM.EGFPARQ.stim.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) 
dev.off()

write.csv(TwoMvSixM.combined.FBSM.EGFPARQ.stim.Top50, "TwoMvSixM.combined.FBSM.EGFPARQ.stim.Top50.csv")

#Heatmap_AllFB
Idents(object = TwoMvSixM.combined.FBSM.EGFPARQ) <- "StroCellType"
TwoMvSixM.combined.FB.EGFPARQ <- subset(x=TwoMvSixM.combined.FBSM.EGFPARQ, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))

DefaultAssay(TwoMvSixM.combined.FB.EGFPARQ) <- "RNA"
Idents(object = TwoMvSixM.combined.FB.EGFPARQ) <- "stim"
TwoMvSixM.combined.FB.EGFPARQ <- ScaleData(TwoMvSixM.combined.FB.EGFPARQ, features = rownames(TwoMvSixM.combined.FB.EGFPARQ))
TwoMvSixM.combined.FB.EGFPARQ.markers <- FindAllMarkers(TwoMvSixM.combined.FB.EGFPARQ, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvSixM.combined.FB.EGFPARQ.Top50 <- TwoMvSixM.combined.FB.EGFPARQ.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)

tiff(file = "TwoMvSixM.combined.FB.EGFPARQ heatmap purple.tiff", width = 10, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvSixM.combined.FB.EGFPARQ, features = c("ARQ", "EGFP", TwoMvSixM.combined.FB.EGFPARQ.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####11M-1####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/11M/A0843")

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

ElevenM1.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/v3_count__ARQ9_length_3191/v3_count_len3191_Project_COHP_44263_1_X3SC4/outs/filtered_feature_bc_matrix")
ElevenM1.unfiltered <- CreateSeuratObject(counts = ElevenM1.data,  min.cells = 3, min.features = 500, project = "ARQ9-Gli1_11M")
ElevenM1.unfiltered <- NormalizeData(ElevenM1.unfiltered)

#Initial processing & filtering

ElevenM1.unfiltered[["percent.mt"]] <- PercentageFeatureSet(ElevenM1.unfiltered, pattern = "^mt-")

tiff(file = "ElevenM1.unfiltered Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(ElevenM1.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "ElevenM1.unfiltered Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM1.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ElevenM1.unfiltered Pre-filteration")
dev.off()
tiff(file = "ElevenM1.unfiltered Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM1.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ElevenM1.unfiltered Pre-filteration")
dev.off()

plot1 <- FeatureScatter(ElevenM1.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ElevenM1.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ElevenM1 <- subset(ElevenM1.unfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 6000 & percent.mt < 10)

tiff(file = "ElevenM1 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(ElevenM1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "ElevenM1 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM1@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "ElevenM1 Post-filteration")
dev.off()
tiff(file = "ElevenM1 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM1@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "ElevenM1 Post-filteration")
dev.off()

plot1 <- FeatureScatter(ElevenM1, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ElevenM1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ElevenM1 <- FindVariableFeatures(ElevenM1, selection.method = "vst", nfeatures = 2500)
tiff(file = "ElevenM1 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(ElevenM1)
dev.off()

#Clustering
all.genes <- rownames(ElevenM1)
ElevenM1 <- ScaleData(ElevenM1, features = all.genes)
ElevenM1 <- RunPCA(ElevenM1, features = VariableFeatures(object = ElevenM1))
print(ElevenM1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ElevenM1, dims = 1:2, reduction = "pca")
DimPlot(ElevenM1, reduction = "pca")

tiff(file = "ElevenM1 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(ElevenM1, ndims = 50)
dev.off()

ElevenM1 <- FindNeighbors(ElevenM1, dims = 1:24)
ElevenM1 <- FindClusters(ElevenM1, resolution = 0.5)
head(Idents(ElevenM1), 5)
ElevenM1 <- RunTSNE(ElevenM1, reduction = "pca", dims = 1:24)
ElevenM1 <- RunUMAP(ElevenM1, reduction = "pca", dims = 1:24)
DimPlot(ElevenM1, reduction = "umap", pt.size = 1, label = TRUE)

DefaultAssay(ElevenM1) <- "RNA"
FeaturePlot(ElevenM1, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(ElevenM1, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(ElevenM1, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(ElevenM1, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)

tiff(file = "ElevenM1 hARtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM1, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ElevenM1 mGFP expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ElevenM1 Ar expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ElevenM1 Gli1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off

####11M-2####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/11M/A1350")

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

ElevenM2.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/v3_count__ARQ9_length_3191/v3_count_len3191_Project_COHP_46246_2_X3SC4/outs/filtered_feature_bc_matrix")
ElevenM2.unfiltered <- CreateSeuratObject(counts = ElevenM2.data,  min.cells = 3, min.features = 500, project = "ARQ9-Gli1_11M")
ElevenM2.unfiltered <- NormalizeData(ElevenM2.unfiltered)

#Initial processing & filtering

ElevenM2.unfiltered[["percent.mt"]] <- PercentageFeatureSet(ElevenM2.unfiltered, pattern = "^mt-")

tiff(file = "ElevenM2.unfiltered Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(ElevenM2.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "ElevenM2.unfiltered Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM2.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ElevenM2.unfiltered Pre-filteration")
dev.off()
tiff(file = "ElevenM2.unfiltered Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM2.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ElevenM2.unfiltered Pre-filteration")
dev.off()

plot1 <- FeatureScatter(ElevenM2.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ElevenM2.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ElevenM2 <- subset(ElevenM2.unfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 15)

tiff(file = "ElevenM2 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(ElevenM2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "ElevenM2 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "ElevenM2 Post-filteration")
dev.off()
tiff(file = "ElevenM2 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(ElevenM2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "ElevenM2 Post-filteration")
dev.off()

plot1 <- FeatureScatter(ElevenM2, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(ElevenM2, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

ElevenM2 <- FindVariableFeatures(ElevenM2, selection.method = "vst", nfeatures = 2500)
tiff(file = "ElevenM2 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(ElevenM2)
dev.off()

#Clustering
all.genes <- rownames(ElevenM2)
ElevenM2 <- ScaleData(ElevenM2, features = all.genes)
ElevenM2 <- RunPCA(ElevenM2, features = VariableFeatures(object = ElevenM2))
print(ElevenM2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(ElevenM2, dims = 1:2, reduction = "pca")
DimPlot(ElevenM2, reduction = "pca")

tiff(file = "ElevenM2 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(ElevenM2, ndims = 50)
dev.off()

ElevenM2 <- FindNeighbors(ElevenM2, dims = 1:20)
ElevenM2 <- FindClusters(ElevenM2, resolution = 0.5)
head(Idents(ElevenM2), 5)
ElevenM2 <- RunTSNE(ElevenM2, reduction = "pca", dims = 1:20)
ElevenM2 <- RunUMAP(ElevenM2, reduction = "pca", dims = 1:20)
DimPlot(ElevenM2, reduction = "umap", pt.size = 1, label = TRUE)

DefaultAssay(ElevenM2) <- "RNA"
FeaturePlot(ElevenM2, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(ElevenM2, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(ElevenM2, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(ElevenM2, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)

tiff(file = "ElevenM2 hARtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM2, reduction = "umap", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ElevenM2 mGFP expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM2, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ElevenM2 Ar expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM2, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ElevenM2 Gli1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM2, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()