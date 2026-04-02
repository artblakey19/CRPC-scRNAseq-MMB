#WK HiMyc WTvARKOvMYCvMYCARKO Workflow

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
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)
library(SingleCellSignalR)

#### Initial Filtering and Clustering ####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO")

###HiMYC-Ctrl### 

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

Mycunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/2nd Run 6-4-2020/Sun Lab Analysis/Project_COHP_36595_1_XCTC1_count/outs/filtered_feature_bc_matrix")
Mycunfiltered <- CreateSeuratObject(counts = Mycunfiltered.data,  min.cells = 3, min.features = 500, project = "Mycunfiltered")
Mycunfiltered<- NormalizeData(Mycunfiltered)

#Initial processing & filtering
Mycunfiltered[["percent.mt"]] <- PercentageFeatureSet(Mycunfiltered, pattern = "^mt-")

tiff(file = "MYC Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Mycunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "MYC Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Mycunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "MYC Pre-filteration")
dev.off()
tiff(file = "MYC Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Mycunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "MYC Pre-filteration")
dev.off()

Myc <- subset(Mycunfiltered, subset = nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 15)
tiff(file = "MYC Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Myc, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "MYC Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Myc@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "MYC Post-filteration")
dev.off()
tiff(file = "MYC Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Myc@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "MYC Post-filteration")
dev.off()

Myc <- FindVariableFeatures(Myc, selection.method = "vst", nfeatures = 5000)
tiff(file = "Myc Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(Myc)
dev.off()

table(Idents(Mycunfiltered))
table(Idents(Myc))

#Run the standard workflow for visualization and clustering Myc
Myc <- ScaleData(Myc, verbose = FALSE)
Myc <- RunPCA(Myc, npcs = 30, verbose = FALSE)
ElbowPlot(Myc, ndims = 50)
# t-SNE and Clustering
Myc <- FindNeighbors(Myc, reduction = "pca", dims = 1:20)
Myc <- FindClusters(Myc, resolution = 0.5)
Myc <- RunTSNE(Myc, reduction = "pca", dims = 1:20)
Myc <- RunUMAP(Myc, reduction = "pca", dims = 1:20)
DimPlot(Myc, reduction = "umap", pt.size = 0.3) 

Idents(object = Myc) <- "stim"

tiff(file = "MYC UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Myc, reduction = "umap", pt.size = 0.3, cols = c("dark blue")) 
dev.off()

DefaultAssay(Myc) <- "RNA"
tiff(file = "Myc MYCtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Myc, reduction = "umap", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

###HiMYC-ARKO###

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

MycARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/190408_190425_combined/28345_count/outs/filtered_feature_bc_matrix")
MycARKOunfiltered <- CreateSeuratObject(counts = MycARKOunfiltered.data, project = "MycARKO", min.cells = 3, min.features = 200)
MycARKOunfiltered <- NormalizeData(MycARKOunfiltered)

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
DimPlot(MycARKO, reduction = "umap", pt.size = 0.3, cols = c("darkred")) 
dev.off()

DefaultAssay(MycARKO) <- "RNA"
tiff(file = "MycARKO MYCtg.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(MycARKO, reduction = "umap", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

###WT### 

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

WTunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/ARKO-Gli1/Single_Cell_Seq/190219_scRNA/27347_WT_count_EGFPmm10/outs/filtered_feature_bc_matrix")
WTunfiltered <- CreateSeuratObject(counts = WTunfiltered.data,  min.cells = 3, min.features = 500, project = "WTunfiltered")
WTunfiltered <- NormalizeData(WTunfiltered)

#Initial processing & filtering
WTunfiltered[["percent.mt"]] <- PercentageFeatureSet(WTunfiltered, pattern = "^mt-")

tiff(file = "WT Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WTunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "grey50", xlab = "nFeature_RNA", main = "WT Pre-filteration")
dev.off()
tiff(file = "WT Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WTunfiltered@meta.data$percent.mt, breaks = 100, col = "grey50", xlab = "percent.mt", main = "WT Pre-filteration")
dev.off()

WT <- subset(WTunfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)
tiff(file = "WT Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "WT Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$nFeature_RNA, breaks = 100, col = "grey50", xlab = "nFeature_RNA", main = "WT Post-filteration")
dev.off()
tiff(file = "WT Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT@meta.data$percent.mt, breaks = 100, col = "grey50", xlab = "percent.mt", main = "WT Post-filteration")
dev.off()

WT <- FindVariableFeatures(WT, selection.method = "vst", nfeatures = 5000)
tiff(file = "WT Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT)
dev.off()

table(Idents(WTunfiltered))
table(Idents(WT))

#Run the standard workflow for visualization and clustering WT
WT <- ScaleData(WT, verbose = FALSE)
WT <- RunPCA(WT, npcs = 30, verbose = FALSE)
ElbowPlot(WT, ndims = 50)
# t-SNE and Clustering
WT <- FindNeighbors(WT, reduction = "pca", dims = 1:25)
WT <- FindClusters(WT, resolution = 0.5)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:25)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:25)
DimPlot(WT, reduction = "umap", pt.size = 0.3) 

Idents(object = WT) <- "stim"

tiff(file = "WT UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, cols = c("light grey")) 
dev.off()


#### Merging Datasets #### ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO/combined")

Myc[["orig.clusters"]] <- Idents(object = Myc)
MycARKO[["orig.clusters"]] <- Idents(object = MycARKO)
WT[["orig.clusters"]] <- Idents(object = WT)

Idents(object = Myc) <- "seurat_clusters"
Idents(object = MycARKO) <- "seurat_clusters"
Idents(object = WT) <- "seurat_clusters"
Myc$stim <- "Myc"
MycARKO$stim <- "MycARKO"
WT$stim <- "WT"

combined.anchors <- FindIntegrationAnchors(object.list = list(MycARKO, Myc, WT), dims = 1:20)
combined <- IntegrateData(anchorset = combined.anchors, dims = 1:20)
DefaultAssay(combined) <- "integrated"

#Run the standard workflow for visualization and clustering Ctrl2vARKO1.combined
combined <- ScaleData(combined, verbose = FALSE)
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)
ElbowPlot(combined, ndims = 50)

combined <- FindNeighbors(combined, reduction = "pca", dims = 1:20)
combined <- FindClusters(combined, resolution = 0.5)
combined <- RunUMAP(combined, reduction = "pca", dims = 1:20)
combined <- RunTSNE(combined, reduction = "pca", dims = 1:20)

Idents(object = combined) <- "stim"
combined <- RenameIdents(object = combined, 'WT' = "WT", 'Myc' = "Myc", 'MycARKO' = "MycARKO")
combined[["stim"]] <- Idents(object = combined)

#Plot
tiff(file = "combined stim UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue", "darkred")) 
dev.off()

#Featureplots
DefaultAssay(combined) <- "RNA"
tiff(file = "combined Epcam.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Vim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Fbln1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Myh11.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Myh11"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Cdh5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Rgs5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Plp1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Plp1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Mki67.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Ccl5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Ccl5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Tyrobp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Tyrobp"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Cdo1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", features = c("Cdo1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()


#Celltype identification
Idents(object = combined) <- "seurat_clusters"
combined <- RenameIdents(object = combined, '0' = "BE", '2' = "BE",
                                   '6' = "LE", '7' = "LE", '8' = "LE", '10' = "LE", 
                                   '12' = "LE", '17' = "ProE", '19' = "OE", '15' = "SV",
                                   '1' = "FB", '3' = "FB", '4' = "FB", '5' = "SM",
                                   '21' = "ProS", '18' = "Pericyte", '20' = "Glia", '11' = "Endo",
                                   '9' = "Leu", '13' = "Leu", '14' = "Leu", '16' = "Leu")
DimPlot(combined, reduction = "umap", pt.size = 1)
combined[["celltype"]] <- Idents(object = combined)

Idents(object = combined) <- "celltype"
combined <- RenameIdents(object = combined, 'BE' = "Epithelia", 'LE' = "Epithelia", 'ProE' = "Epithelia", 'OE' = "Epithelia", 'SV' = "Epithelia", 'FB' = "Stroma", 'SM' = "Stroma", 'ProS' = "Stroma", 'Pericyte' = "Stroma", 'Glia' = "Stroma", 'Endo' = "Stroma", 'Leu' = "Stroma")
DimPlot(combined, reduction = "umap", pt.size = 1)
combined[["EpiStro"]] <- Idents(object = combined)


#Dimplot
Idents(object = combined) <- "celltype"
tiff(file = "combined celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = combined) <- "stim"
combined.WT <- subset(combined, idents = c("WT"))
Idents(object = combined.WT) <- "celltype"
tiff(file = "combined.WT celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.WT, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = combined) <- "stim"
combined.Myc <- subset(combined, idents = c("Myc"))
Idents(object = combined.Myc) <- "celltype"
tiff(file = "combined.Myc celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.Myc, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = combined) <- "stim"
combined.MycARKO <- subset(combined, idents = c("MycARKO"))
Idents(object = combined.MycARKO) <- "celltype"
tiff(file = "combined.MycARKO celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.MycARKO, reduction = "umap", pt.size = 0.3)
dev.off()

#DEGS
DefaultAssay(combined) <- "RNA"
Idents(object = combined) <- "celltype"
combined <- ScaleData(combined, features = rownames(CtrlvARKO.combined))
combined.celltype.markers <- FindAllMarkers(CtrlvARKO.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.celltype.markers, "combined.celltype.markers.csv")

#Dotplot

#Feature Plots Split
DefaultAssay(combined) <- "RNA"
tiff(file = "combined Ar split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Pcna split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Pcna"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined MYCtg split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Pbsn split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined EGFP split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined Gli1 split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Cell counts
Idents(object = combined) <- "stim"
combined$stim.celltype <- paste(Idents(combined), combined$celltype, sep = "_")
Idents(object = combined) <- "stim.celltype"
table(Idents(combined))

####Subclustering FB/SM####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO/combined.FBSM")

#Dimplot
Idents(object = combined) <- "EpiStro"
tiff(file = "combined Stro highlight UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue"))
dev.off()

#subclustering Stroma
Idents(object = combined) <- "EpiStro"
combined.Stro <- subset(combined, idents = c("Stroma"))
DefaultAssay(combined.Stro) <- "integrated"

#Run the standard workflow for visualization and clustering
combined.Stro <- ScaleData(combined.Stro, verbose = FALSE)
combined.Stro <- RunPCA(combined.Stro, npcs = 30, verbose = FALSE)
ElbowPlot(combined.Stro, ndims = 50)
#Umap and Clustering
combined.Stro <- FindNeighbors(combined.Stro, reduction = "pca", dims = 1:16)
combined.Stro <- FindClusters(combined.Stro, resolution = 0.5)
combined.Stro <- RunTSNE(combined.Stro, reduction = "pca", dims = 1:16)
combined.Stro <- RunUMAP(combined.Stro, reduction = "pca", dims = 1:16)
DimPlot(combined.Stro, reduction = "umap", label = TRUE)

#Featureplots
DefaultAssay(combined.Stro) <- "RNA"

tiff(file = "combined.Stro Fbln1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Myh11.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Myh11"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Cdh5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Rgs5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Plp1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Plp1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Mki67.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Ccl5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Ccl5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Stro Tyrobp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Stro, reduction = "umap", features = c("Tyrobp"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

#Celltype identification
Idents(object = combined.Stro) <- "seurat_clusters"
combined.Stro <- RenameIdents(object = combined.Stro, '0' = "FB", '1' = "FB",
                         '2' = "FB", '4' = "FB", '5' = "FB", '3' = "SM", 
                         '12' = "ProS", '14' = "Glia", '11' = "Pericyte", '7' = "Endo",
                         '6' = "Leu", '8' = "Leu", '9' = "Leu", '10' = "Leu",
                         '15' = "Leu", '13' = "OS")
DimPlot(combined.Stro, reduction = "umap", pt.size = 1)
combined.Stro[["Strocelltype"]] <- Idents(object = combined.Stro)

Idents(object = combined.Stro) <- "Strocelltype"
tiff(file = "combined.Stro celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.Stro, reduction = "umap", pt.size = 0.3)
dev.off()

####Subclustering FBSM####
Idents(object = combined.Stro) <- "Strocelltype"

combined.FBSM <- subset(combined.Stro, idents = c("FB", "SM"))
combined.FBSM <- RunTSNE(combined.FBSM, reduction = "pca", dims = 1:15)
combined.FBSM <- RunUMAP(combined.FBSM, reduction = "pca", dims = 1:15)

Idents(object = combined.FBSM) <- "seurat_clusters"
DimPlot(combined.FBSM, reduction = "umap", label = TRUE)

#Reclustering
DefaultAssay(combined.FBSM) <- "integrated"
combined.FBSM <- ScaleData(combined.FBSM, verbose = FALSE)
combined.FBSM <- RunPCA(combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(combined.FBSM, ndims = 50)
#Umap and Clustering
combined.FBSM <- FindNeighbors(combined.FBSM, reduction = "pca", dims = 1:13)
combined.FBSM <- FindClusters(combined.FBSM, resolution = 0.5)
combined.FBSM <- RunTSNE(combined.FBSM, reduction = "pca", dims = 1:13)
combined.FBSM <- RunUMAP(combined.FBSM, reduction = "pca", dims = 1:13)
DimPlot(combined.FBSM, reduction = "umap", label = TRUE)
DimPlot(combined.FBSM, split.by = "stim", reduction = "umap", label = TRUE)

#Rename
Idents(object = combined.FBSM) <- "seurat_clusters"
combined.FBSM <- RenameIdents(object = combined.FBSM, '0' = "FB1", '1' = "FB2",
                              '8' = "FB3", '2' = "FB4", '6' = "FB5", '4' = "FB6", 
                              '3' = "FB7", '5' = "SM1", '7' = "SM2")
DimPlot(combined.FBSM, reduction = "umap", pt.size = 1)
combined.FBSM[["FBSMcelltype"]] <- Idents(object = combined.FBSM)

#Dimplot
tiff(file = "combined.FBSM FBSMcelltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.FBSM, reduction = "umap", pt.size = 0.7)
dev.off()
tiff(file = "combined.FBSM FBSMcelltype UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.FBSM, reduction = "umap", split.by = "stim", pt.size = 0.7)
dev.off()

#Violinplot
Idents(object = combined.FBSM) <- "stim"
combined.FBSM.WT <- subset(combined.FBSM, idents = c("WT"))

DefaultAssay(combined.FBSM.WT) <- "RNA"
Idents(object = combined.FBSM.WT) <- "FBSMcelltype"

tiff(file = "combined.FBSM.WT Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.WT, features = "Ar", pt.size = 0)
dev.off()
tiff(file = "combined.FBSM.WT EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.WT, features = "EGFP", pt.size = 0)
dev.off()

Idents(object = combined.FBSM) <- "stim"
combined.FBSM.Myc <- subset(combined.FBSM, idents = c("Myc"))
DefaultAssay(combined.FBSM.Myc) <- "RNA"
Idents(object = combined.FBSM.Myc) <- "FBSMcelltype"
tiff(file = "combined.FBSM.Myc Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.Myc, features = "Ar", pt.size = 0)
dev.off()
tiff(file = "combined.FBSM.Myc EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.Myc, features = "EGFP", pt.size = 0)
dev.off()

Idents(object = combined.FBSM) <- "stim"
combined.FBSM.MycARKO <- subset(combined.FBSM, idents = c("MycARKO"))
DefaultAssay(combined.FBSM.MycARKO) <- "RNA"
Idents(object = combined.FBSM.MycARKO) <- "FBSMcelltype"
tiff(file = "combined.FBSM.MycARKO Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.MycARKO, features = "Ar", pt.size = 0)
dev.off()
tiff(file = "combined.FBSM.MycARKO EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM.MycARKO, features = "EGFP", pt.size = 0)
dev.off()

#Cell counts
Idents(object = combined.FBSM) <- "stim"
combined.FBSM$stim.FBSMcelltype <- paste(Idents(combined.FBSM), combined.FBSM$FBSMcelltype, sep = "_")
Idents(object = combined.FBSM) <- "stim.FBSMcelltype"
table(Idents(combined.FBSM))

#FeaturePlot
DefaultAssay(combined.FBSM) <- "RNA"
tiff(file = "combined.FBSM Ar Split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FBSM EGFP Split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FBSM Igfbp3 Split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Igfbp3"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.FBSM Cxcl14 Split.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Cxcl14"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

#DEGs
DefaultAssay(combined.FBSM) <- "RNA"
Idents(object = combined.FBSM) <- "FBSMcelltype"
combined.FBSM <- ScaleData(combined.FBSM, features = rownames(combined.FBSM))
combined.FBSM.FBSMcelltype.marker <- FindAllMarkers(combined.FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.FBSM.FBSMcelltype.marker, "combined.FBSM.FBSMcelltype.marker.csv")

