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

#### Initial Filtering and Clustering ####

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO_P42")

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

###WT### 

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

WTunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/CellRanger_Output/count_Project_COHP_36980_1_X3SC3/outs/filtered_feature_bc_matrix")
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

WT <- subset(WTunfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & percent.mt < 15)
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

tiff(file = "WT ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(WT, ndims = 50)
dev.off()

# t-SNE and Clustering
WT <- FindNeighbors(WT, reduction = "pca", dims = 1:20)
WT <- FindClusters(WT, resolution = 0.5)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:20)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:20)
DimPlot(WT, reduction = "umap", pt.size = 0.3) 

Idents(object = WT) <- "stim"

tiff(file = "WT UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, cols = c("light grey")) 
dev.off()


#### Merging Datasets #### ####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO_P42/combined")

Myc[["orig.clusters"]] <- Idents(object = Myc)
WT[["orig.clusters"]] <- Idents(object = WT)

Idents(object = Myc) <- "seurat_clusters"
Idents(object = WT) <- "seurat_clusters"
Myc$stim <- "Myc"
WT$stim <- "WT"

combined.anchors <- FindIntegrationAnchors(object.list = list(Myc, WT), dims = 1:20)
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
combined <- RenameIdents(object = combined, 'WT' = "WT", 'Myc' = "Myc")
combined[["stim"]] <- Idents(object = combined)

#Plot
tiff(file = "combined stim UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 
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

tiff(file = "combined Myctg.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()


#Celltype identification
Idents(object = combined) <- "seurat_clusters"
DimPlot(combined, reduction = "umap", pt.size = 0.3, label = TRUE) 

combined <- RenameIdents(object = combined, '0' = "BE", '1' = "BE", '19' = "BE",
                                   '3' = "LE", '10' = "LE", '10' = "LE", 
                                   '5' = "LE", '6' = "LE", '12' = "LE", '20' = "OE", '15' = "OE", 
                                   '17' = "ProE", '11' = "SV", '2' = "FB", '4' = "FB",
                                   '16' = "SM", '13' = "Pericyte", '14' = "Glia", '9' = "Endo",
                                   '7' = "Leu", '8' = "Leu", '18' = "Leu")
DimPlot(combined, reduction = "umap", pt.size = 1)
combined[["celltype"]] <- Idents(object = combined)

Idents(object = combined) <- "celltype"
combined <- RenameIdents(object = combined, 'BE' = "Epithelia", 'LE' = "Epithelia", 'ProE' = "Epithelia", 'OE' = "Epithelia", 'SV' = "Epithelia", 'FB' = "Stroma", 'SM' = "Stroma", 'Pericyte' = "Stroma", 'Glia' = "Stroma", 'Endo' = "Stroma", 'Leu' = "Stroma")
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

#DEGS
DefaultAssay(combined) <- "RNA"
Idents(object = combined) <- "celltype"
combined <- ScaleData(combined, features = rownames(combined))
combined.celltype.markers <- FindAllMarkers(combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.celltype.markers, "combined.celltype.markers.csv")

#Dotplot

#Feature Plots Split
DefaultAssay(combined) <- "RNA"
tiff(file = "combined Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined Pcna split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Pcna"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined MYCtg split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined Pbsn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined EGFP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined Gli1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined, reduction = "umap", split.by = "stim", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

#Cell counts
Idents(object = combined) <- "stim"
combined$stim.celltype <- paste(Idents(combined), combined$celltype, sep = "_")
Idents(object = combined) <- "stim.celltype"
table(Idents(combined))

####Subclustering FB/SM####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO_P42/combined.FBSM")

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
combined.Stro <- FindNeighbors(combined.Stro, reduction = "pca", dims = 1:20)
combined.Stro <- FindClusters(combined.Stro, resolution = 0.5)
combined.Stro <- RunTSNE(combined.Stro, reduction = "pca", dims = 1:20)
combined.Stro <- RunUMAP(combined.Stro, reduction = "pca", dims = 1:20)
DimPlot(combined.Stro, reduction = "umap", label = TRUE)
DimPlot(combined.Stro, reduction = "umap", split.by = "stim", label = TRUE)

#Featureplots
DefaultAssay(combined.Stro) <- "RNA"
FeaturePlot(combined.Stro, reduction = "umap", split.by = "stim", features = c("Cxcl14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)



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
combined.Stro <- RenameIdents(object = combined.Stro, '0' = "FB", '1' = "FB", '2' = "FB", '6' = "FB", '8' = "FB",
                         '9' = "SM", 
                         '11' = "Pericyte", '14' = "Pericyte", '13' = "Glia", '4' = "Endo", '17' = "Endo", 
                         '3' = "Leu", '15' = "Leu", '12' = "Leu", '16' = "Leu",
                         '10' = "Leu", '5' = "Leu", '7' = "OS")
DimPlot(combined.Stro, reduction = "umap", pt.size = 1)
combined.Stro[["Strocelltype"]] <- Idents(object = combined.Stro)

Idents(object = combined.Stro) <- "Strocelltype"
tiff(file = "combined.Stro celltype UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.Stro, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "combined.Stro FBSM Highlight UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.Stro, reduction = "umap", pt.size = 0.3, cols = c("darkred", "darkred", "light grey", "light grey", "light grey", "light grey", "light grey"))
dev.off()

#Subclustering FBSM
Idents(object = combined.Stro) <- "Strocelltype"
combined.FBSM <- subset(combined.Stro, idents = c("FB", "SM"))

#Reclustering
DefaultAssay(combined.FBSM) <- "integrated"
combined.FBSM <- ScaleData(combined.FBSM, verbose = FALSE)
combined.FBSM <- RunPCA(combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(combined.FBSM, ndims = 50)
#Umap and Clustering
combined.FBSM <- FindNeighbors(combined.FBSM, reduction = "pca", dims = 1:15)
combined.FBSM <- FindClusters(combined.FBSM, resolution = 0.5)
combined.FBSM <- RunTSNE(combined.FBSM, reduction = "pca", dims = 1:15)
combined.FBSM <- RunUMAP(combined.FBSM, reduction = "pca", dims = 1:15)
DimPlot(combined.FBSM, reduction = "umap", label = TRUE)
DimPlot(combined.FBSM, split.by = "stim", reduction = "umap", label = TRUE)

#Rename
Idents(object = combined.FBSM) <- "seurat_clusters"
combined.FBSM <- RenameIdents(object = combined.FBSM, '2' = "FB1", '0' = "FB2",
                              '3' = "FB3", '6' = "FB4", '8' = "FB5", '5' = "FB6", 
                              '1' = "FB7", '4' = "FB8", '9' = "FB9", '7' = "SM")
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
DefaultAssay(combined.FBSM) <- "RNA"
tiff(file = "combined.FBSM Ar Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM, features = "Ar", group.by = "FBSMcelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "combined.FBSM EGFP Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FBSM, features = "EGFP", group.by = "FBSMcelltype", split.by = "stim", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()

#Cell counts
Idents(object = combined.FBSM) <- "stim"
combined.FBSM$stim.FBSMcelltype <- paste(Idents(combined.FBSM), combined.FBSM$FBSMcelltype, sep = "_")
Idents(object = combined.FBSM) <- "stim.FBSMcelltype"
table(Idents(combined.FBSM))

#DEGs
DefaultAssay(combined.FBSM) <- "RNA"
Idents(object = combined.FBSM) <- "FBSMcelltype"
combined.FBSM <- ScaleData(combined.FBSM, features = rownames(combined.FBSM))
combined.FBSM.FBSMcelltype.marker <- FindAllMarkers(combined.FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(combined.FBSM.FBSMcelltype.marker, "combined.FBSM.FBSMcelltype.marker.csv")

#DEGs
Idents(object = combined.FBSM) <- "stim"
combined.FBSM.Myc <- subset(combined.FBSM, idents = c("Myc"))
Idents(object = combined.FBSM.Myc) <- "FBSMcelltype"
DefaultAssay(combined.FBSM.Myc) <- "RNA"
combined.FBSM.Myc <- ScaleData(combined.FBSM.Myc, features = rownames(combined.FBSM.Myc))

combined.FBSM.Myc.FB6.IPA.Markers <- FindMarkers(combined.FBSM.Myc, ident.1 = "FB6", ident.2 = c("FB1", "FB2", "FB3", "FB4", "FB5", "FB7", "FB8", "FB9"), min.pct = 0.1, logfc.threshold = 0.1)
write.csv(combined.FBSM.Myc.FB6.IPA.Markers, "combined.FBSM.Myc.FB6.IPA.Markers.csv")
combined.FBSM.Myc.FB6.GSEA.Markers <- FindMarkers(combined.FBSM.Myc, ident.1 = "FB6", ident.2 = c("FB1", "FB2", "FB3", "FB4", "FB5", "FB7", "FB8", "FB9"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(combined.FBSM.Myc.FB6.GSEA.Markers, "combined.FBSM.Myc.FB6.GSEA.Markers.csv")

#DEGs
Idents(object = combined.FBSM) <- "stim.FBSMcelltype"
DefaultAssay(combined.FBSM) <- "RNA"
combined.FBSM <- ScaleData(combined.FBSM, features = rownames(combined.FBSM.Myc))

HiMYCvWT.FBSM.FB6.IPA.Markers <- FindMarkers(combined.FBSM, ident.1 = "Myc_FB6", ident.2 = "WT_FB6", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(HiMYCvWT.FBSM.FB6.IPA.Markers, "HiMYCvWT.FBSM.FB6.IPA.Markers.csv")
HiMYCvWT.FBSM.FB6.GSEA.Markers <- FindMarkers(combined.FBSM, ident.1 = "Myc_FB6", ident.2 = "WT_FB6", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(HiMYCvWT.FBSM.FB6.GSEA.Markers, "HiMYCvWT.FBSM.FB6.GSEA.Markers.csv")

#Violin Plot
Idents(object = combined.FBSM) <- "FBSMcelltype"
combined.FB6 <- subset(combined.FBSM, idents = c("FB6"))

DefaultAssay(combined.FB6) <- "RNA"
Idents(object = combined.FB6) <- "stim"
tiff(file = "combined.FB6 Wnt2 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Wnt2", group.by = "stim", pt.size = 0)
dev.off()
tiff(file = "combined.FB6 Igfbp3 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Igfbp3", group.by = "stim", pt.size = 0)
dev.off()
tiff(file = "combined.FB6 Fn1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Fn1", group.by = "stim", pt.size = 0)
dev.off()
tiff(file = "combined.FB6 Il11 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Il11", group.by = "stim", pt.size = 0)
dev.off()

tiff(file = "combined.FB6 Fn1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Fn1", group.by = "stim", pt.size = 0)
dev.off()
tiff(file = "combined.FB6 Foxf1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Foxf1", group.by = "stim", pt.size = 0)
dev.off()
tiff(file = "combined.FB6 Twist1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Twist1", group.by = "stim", pt.size = 0)
dev.off()
tiff(file = "combined.FB6 Il11 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.FB6, features = "Il11", group.by = "stim", pt.size = 0)
dev.off()


#Boxplot Generation
boxdata = FetchData(combined.FB6, c("stim", "Il11", "Wnt2", "Igfbp3", "Fn1"))
tail(boxdata,5)

tiff(file = "combined.FB6 Il11 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Il11, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "combined.FB6 Wnt2 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Il11, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "combined.FB6 Igfbp3 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Igfbp3, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()
tiff(file = "combined.FB6 Fn1 Boxplot with Dots and Median.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
ggplot(boxdata, aes(x=stim, y=Fn1, fill = stim)) + geom_boxplot(outlier.color = NA)  + geom_jitter(size = 0.25) + stat_summary(fun.y = median, geom='point', size = 15, colour = "red", shape = 95)
dev.off()

#FeaturePlot
DefaultAssay(combined.FBSM) <- "RNA"

tiff(file = "combined.FBSM Fn1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Fn1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Il11 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Il11"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Igfbp3 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Igfbp3"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Wnt2 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Wnt2"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

####CAF markers####
DefaultAssay(combined.FBSM) <- "RNA"
tiff(file = "combined.FBSM Il11 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Il11"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Cxcl10 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Cxcl10"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Sox9 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Sox9"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Twist1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Twist1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Foxf1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Foxf1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Aldh1a1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Aldh1a1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.FBSM Rbp4 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.FBSM, split.by = "stim", reduction = "umap", features = c("Rbp4"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()

DefaultAssay(CtrlvARKO.combined.FBSM) <- "RNA"
tiff(file = "CtrlvARKO.combined.FBSM Il11 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Il11"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Cxcl10 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Cxcl10"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Sox9 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Sox9"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Twist1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Twist1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Foxf1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Foxf1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Aldh1a1 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Aldh1a1"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "CtrlvARKO.combined.FBSM Rbp4 Split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.FBSM, split.by = "stim", reduction = "umap", features = c("Rbp4"), cols = c("light grey", "purple"), pt.size = 1, min.cutoff = 0, max.cutoff = 2.5)
dev.off()


####Subclustering Epi####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/WTvMYCvMYCARKO_P42/combined.Epi")

#Dimplot
Idents(object = combined) <- "EpiStro"
tiff(file = "combined Epi highlight UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined, reduction = "umap", pt.size = 0.3, cols = c("darkblue", "light grey"))
dev.off()

#subclustering Stroma
Idents(object = combined) <- "EpiStro"
combined.Epi <- subset(combined, idents = c("Epithelia"))
DefaultAssay(combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
combined.Epi <- ScaleData(combined.Epi, verbose = FALSE)
combined.Epi <- RunPCA(combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(combined.Epi, ndims = 50)
#Umap and Clustering
combined.Epi <- FindNeighbors(combined.Epi, reduction = "pca", dims = 1:15)
combined.Epi <- FindClusters(combined.Epi, resolution = 0.5)
combined.Epi <- RunTSNE(combined.Epi, reduction = "pca", dims = 1:15)
combined.Epi <- RunUMAP(combined.Epi, reduction = "pca", dims = 1:15)
DimPlot(combined.Epi, reduction = "umap", label = TRUE)

#Featureplots
DefaultAssay(combined.Epi) <- "RNA"

tiff(file = "combined.Epi Krt5.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Epi Krt19.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt8.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt14.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Krt14"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Epi Vim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Epi Mki67.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()
tiff(file = "combined.Epi Cdo1.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.Epi, reduction = "umap", features = c("Cdo1"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = 0, max.cutoff = 2.5)
dev.off()


#Celltype identification
Idents(object = combined.Epi) <- "seurat_clusters"
combined.Epi <- RenameIdents(object = combined.Epi, '0' = "BE", '2' = "BE",
                              '3' = "BE", '9' = "BE", '15' = "BE", '11' = "LE", '5' = "LE", 
                              '1' = "LE", '14' = "LE", '6' = "LE", '7' = "LE",
                              '15' = "LE", '13' = "LE", '4' = "LE", '12' = "LE",
                              '10' = "OE", '8' = "SV")
DimPlot(combined.Epi, reduction = "umap", pt.size = 1)
combined.Epi[["Epicelltype"]] <- Idents(object = combined.Epi)

Idents(object = combined.Epi) <- "Epicelltype"
tiff(file = "combined.Epi celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = combined.Epi) <- "Epicelltype"
tiff(file = "combined.Epi BELE highlight UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("dark red", "dark red", "light grey", "light grey", "light grey"))
dev.off()

#Add MYC info
DefaultAssay(combined.Epi) <- "RNA"
combined.Epi.MYCPos <- subset(x=combined.Epi,  subset = `MYC-transgene` > 0)
combined.Epi.MYCNeg <- subset(x=combined.Epi,  subset = `MYC-transgene` == 0)
Idents(object = combined.Epi.MYCPos) <- "MYCPos"
Idents(object = combined.Epi.MYCNeg) <- "MYCNeg"
combined.Epi.MYCPos[["MYCExp"]] <- Idents(object = combined.Epi.MYCPos)
combined.Epi.MYCNeg[["MYCExp"]] <- Idents(object = combined.Epi.MYCNeg)
combined.EpiMYC <- merge(x = combined.Epi.MYCPos, y = combined.Epi.MYCNeg)
Idents(object = combined.EpiMYC) <- "MYCExp"
combined.Epi$MYCExp <- Idents(object = combined.EpiMYC)
Idents(object = combined.Epi) <- "MYCExp"

Idents(object = combined.Epi) <- "stim"
combined.Epi$stim.MYCExp <- paste(Idents(combined.Epi), combined.Epi$MYCExp, sep = "_")
Idents(object = combined.Epi) <- "stim.MYCExp"
DimPlot(combined.Epi, reduction = "umap")

Idents(object = combined.Epi) <- "Epicelltype"
combined.Epi.BE <- subset(combined.Epi, idents = c("BE"))
combined.Epi.LE <- subset(combined.Epi, idents = c("LE"))
Idents(object = combined.Epi.BE) <- "stim.MYCExp"
DimPlot(combined.Epi.BE, reduction = "umap")

#DEGs
DefaultAssay(combined.Epi.BE) <- "RNA"
Idents(object = combined.Epi.BE) <- "stim.MYCExp"
combined.Epi.BE <- ScaleData(combined.Epi.BE, features = rownames(combined.Epi.BE))
combined.Epi.BE.GSEA.Markers <- FindMarkers(combined.Epi.BE, ident.1 = "Myc_MYCPos", ident.2 = "WT_MYCNeg", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(combined.Epi.BE.GSEA.Markers, " combined.Epi.BE.GSEA.Markers.csv")

#ViolinPlots
DefaultAssay(combined.Epi.BE) <- "RNA"
Idents(object = combined.Epi.BE) <- "stim.MYCExp"
combined.BE.MYC <- subset(combined.Epi.BE, idents = c("Myc_MYCPos", "WT_MYCNeg"))
DefaultAssay(combined.BE.MYC) <- "RNA"
tiff(file = "combined.BE.MYC Igf1r Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.BE.MYC, features = "Igf1r", pt.size = 0)
dev.off()

#DEGs
DefaultAssay(combined.Epi.LE) <- "RNA"
Idents(object = combined.Epi.LE) <- "stim.MYCExp"
combined.Epi.LE <- ScaleData(combined.Epi.LE, features = rownames(combined.Epi.LE))
combined.Epi.LE.GSEA.Markers <- FindMarkers(combined.Epi.LE, ident.1 = "Myc_MYCPos", ident.2 = "WT_MYCNeg", min.pct = 0.01, logfc.threshold = 0.01)
write.csv(combined.Epi.LE.GSEA.Markers, "combined.Epi.LE.GSEA.Markers.csv")

###Subclustering BELE###
Idents(object = combined.Epi) <- "Epicelltype"
combined.BELE <- subset(combined.Epi, idents = c("BE", "LE"))

#Reclustering
DefaultAssay(combined.BELE) <- "integrated"
combined.BELE <- ScaleData(combined.BELE, verbose = FALSE)
combined.BELE <- RunPCA(combined.BELE, npcs = 30, verbose = FALSE)
ElbowPlot(combined.BELE, ndims = 50)
#Umap and Clustering
combined.BELE <- FindNeighbors(combined.BELE, reduction = "pca", dims = 1:15)
combined.BELE <- FindClusters(combined.BELE, resolution = 0.5)
combined.BELE <- RunTSNE(combined.BELE, reduction = "pca", dims = 1:15)
combined.BELE <- RunUMAP(combined.BELE, reduction = "pca", dims = 1:15)
DimPlot(combined.BELE, reduction = "umap", label = TRUE)
DimPlot(combined.BELE, split.by = "stim", reduction = "umap", label = TRUE)

#Featureplots
DefaultAssay(combined.BELE) <- "RNA"
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.7, max.cutoff = "q80")

#Rename
Idents(object = combined.BELE) <- "seurat_clusters"
combined.BELE <- RenameIdents(object = combined.BELE, '10' = "BE1", '0' = "BE2",
                             '2' = "BE3", '13' = "BE3", '6' = "BE4", '3' = "BE5", '15' = "ProBE",  
                             '9' = "LE1", '8' = "LE2", '1' = "LE3", '7' = "LE4", '4' = "LE5",
                             '12' = "LE6", '5' = "LE7", '11' = "ProLE", '14' = "OE")
DimPlot(combined.BELE, reduction = "umap", pt.size = 1)
combined.BELE[["BELEcelltype"]] <- Idents(object = combined.BELE)

#Dimplot
Idents(object = combined.BELE) <- "BELEcelltype"
tiff(file = "combined.BELE celltype UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.BELE, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "combined.BELE celltype split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(combined.BELE, reduction = "umap", split.by = "stim", pt.size = 0.3)
dev.off()

#Cell counts
Idents(object = combined.BELE) <- "stim"
combined.BELE$stim.BELEcelltype <- paste(Idents(combined.BELE), combined.BELE$BELEcelltype, sep = "_")
Idents(object = combined.BELE) <- "stim.BELEcelltype"
table(Idents(combined.BELE))

#Featureplots
DefaultAssay(combined.BELE) <- "RNA"
tiff(file = "combined.BELE Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.BELE Pbsn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.BELE Myctg split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("MYC-transgene"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.BELE Krt5 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.5, max.cutoff = "q90")
dev.off()
tiff(file = "combined.BELE Krt19 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("Krt19"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

tiff(file = "combined.BELE Krt19 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(combined.BELE, reduction = "umap", split.by = "stim", features = c("Igf1r"), cols = c("light grey", "purple"), pt.size = 0.5, min.cutoff = "q10", max.cutoff = "q90")
dev.off()

DefaultAssay(combined.BELE) <- "RNA"
tiff(file = "combined.BELE MYCtg Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(combined.BELE, features = "MYC-transgene", group.by = "BELEcelltype", split.by = "stim", pt.size = 0)
dev.off()

