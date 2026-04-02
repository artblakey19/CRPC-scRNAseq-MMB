#HGF-MET SCseq Work Flow

#### Download librarys ####

install.packages('Seurat')
install.packages('devtools')
install.packages('dplyr')
install.packages('Matrix')
install.packages('cowplot')
install.packages('ggplot2')
install.packages('reticulate')
install.packages("BiocManager")
BiocManager::install("monocle")
install.packages('RColorBrewer')

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

#### Initial Filtering and Clustering ####

setwd("//isi-dcnl/user_data/zjsun/group/Christian Nenninger/HGF")

#HGF_MET

HGF_METunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/201120_TGen/new_counts/count_for_163bp_hMETtg2_mm10/count_39553_HGF_hMET_Pb_8M/outs/filtered_feature_bc_matrix")
HGF_METunfiltered <- CreateSeuratObject(counts = HGF_METunfiltered.data,  min.cells = 3, min.features = 500, project = "HGF_METunfiltered")
HGF_METunfiltered <- NormalizeData(HGF_METunfiltered)

HGF_METunfiltered[["percent.mt"]] <- PercentageFeatureSet(HGF_METunfiltered, pattern = "^mt-")

#new filtering paramaters

table(Idents(HGF_METunfiltered))

tiff(file = "HGF_MET Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_METunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()

tiff(file = "HGF_MET Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGF_METunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "HGF_MET pre filteration")
dev.off()

tiff(file = "HGF_MET Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGF_METunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "HGF_MET pre filteration")
dev.off()

HGF_MET <- subset(HGF_METunfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 15)

table(Idents(HGF_MET))
tiff(file = "HGF_MET Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(HGF_MET, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()
tiff(file = "HGF_MET Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGF_MET@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "HGF_MET post filteration")
dev.off()
tiff(file = "HGF_MET Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(HGF_MET@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "HGF_MET post filteration")
dev.off()

HGF_MET <- FindVariableFeatures(HGF_MET, selection.method = "vst", nfeatures = 2500)
tiff(file = "HGF_MET Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(HGF_MET)
dev.off()


all.genes <- rownames(HGF_MET)
HGF_MET <- ScaleData(HGF_MET, features = all.genes)
HGF_MET <- RunPCA(HGF_MET, features = VariableFeatures(object = HGF_MET))
print(HGF_MET[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(HGF_MET, dims = 1:2, reduction = "pca")
DimPlot(HGF_MET, reduction = "pca")
ElbowPlot(HGF_MET, ndims = 50)

tiff(file = "HGF_MET ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(HGF_MET, ndims = 50)
dev.off()

HGF_MET <- FindNeighbors(HGF_MET, dims = 1:25)
HGF_MET <- FindClusters(HGF_MET, resolution = 0.5)
head(Idents(HGF_MET), 5)
HGF_MET <- RunTSNE(HGF_MET, dims = 1:25)
HGF_MET <- RunUMAP(HGF_MET, dims = 1:25)
DimPlot(HGF_MET, reduction = "umap", pt.size = 0.5, label = TRUE)

Idents(object = HGF_MET) <- "stim"
tiff(file = "HGF_MET UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(HGF_MET, reduction = "umap", pt.size = 0.3, cols = c("darkblue")) 
dev.off()

#WT1

WT1unfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/190219_scRNA/fastq_path/Zijie_Sun/27347_count_ref5/outs/filtered_feature_bc_matrix")
WT1unfiltered <- CreateSeuratObject(counts = WT1unfiltered.data,  min.cells = 3, min.features = 500, project = "WT1unfiltered")
WT1unfiltered <- NormalizeData(WT1unfiltered)

WT1unfiltered[["percent.mt"]] <- PercentageFeatureSet(WT1unfiltered, pattern = "^mt-")

#new filtering paramaters

table(Idents(WT1unfiltered))
tiff(file = "WT1 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT1unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()
tiff(file = "WT1 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT1unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "skyblue", xlab = "nFeature_RNA", main = "WT pre filteration")
dev.off()
tiff(file = "WT1 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT1unfiltered@meta.data$percent.mt, breaks = 100, col = "skyblue", xlab = "percent.mt", main = "WT pre filteration")
dev.off()

WT1 <- subset(WT1unfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 8000 & percent.mt < 15)

table(Idents(WT1))
tiff(file = "WT1 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()
tiff(file = "WT1 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT1@meta.data$nFeature_RNA, breaks = 100, col = "skyblue", xlab = "nFeature_RNA", main = "WT post filteration")
dev.off()
tiff(file = "WT1 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT1@meta.data$percent.mt, breaks = 100, col = "skyblue", xlab = "percent.mt", main = "WT post filteration")
dev.off()

WT1 <- FindVariableFeatures(WT1, selection.method = "vst", nfeatures = 2500)
tiff(file = "WT1 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT1)
dev.off()

all.genes <- rownames(WT1)
WT1 <- ScaleData(WT1, features = all.genes)
WT1 <- RunPCA(WT1, features = VariableFeatures(object = WT1))
print(WT1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WT1, dims = 1:2, reduction = "pca")
DimPlot(WT1, reduction = "pca")
ElbowPlot(WT1, ndims = 50)
tiff(file = "WT1 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(WT1, ndims = 50)
dev.off()

WT1 <- FindNeighbors(WT1, dims = 1:26)
WT1 <- FindClusters(WT1, resolution = 0.5)
head(Idents(WT1), 5)
WT1 <- RunTSNE(WT1, dims = 1:26)
WT1 <- RunUMAP(WT1, dims = 1:26)
DimPlot(WT1, reduction = "umap", pt.size = 1, label = TRUE)

#WT2

WT2unfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/36980_4-P42-WT-Gli1/36980_count_ref5/outs/filtered_feature_bc_matrix")
WT2unfiltered <- CreateSeuratObject(counts = WT2unfiltered.data,  min.cells = 3, min.features = 500, project = "WT2unfiltered")
WT2unfiltered <- NormalizeData(WT2unfiltered)

WT2unfiltered[["percent.mt"]] <- PercentageFeatureSet(WT2unfiltered, pattern = "^mt-")

#new filtering paramaters

table(Idents(WTunfiltered))
tiff(file = "WT2 Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT2unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()
tiff(file = "WT2 Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT2unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "purple", xlab = "nFeature_RNA", main = "WT2 pre filteration")
dev.off()
tiff(file = "WT2 Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT2unfiltered@meta.data$percent.mt, breaks = 100, col = "purple", xlab = "percent.mt", main = "WT2 pre filteration")
dev.off()

WT2 <- subset(WT2unfiltered, subset = nFeature_RNA > 200 & nFeature_RNA < 8000 & percent.mt < 15)

table(Idents(WT2))
tiff(file = "WT2 Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(WT2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
dev.off()
tiff(file = "WT2 Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT2@meta.data$nFeature_RNA, breaks = 100, col = "purple", xlab = "nFeature_RNA", main = "WT2 post filteration")
dev.off()
tiff(file = "WT2 Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(WT2@meta.data$percent.mt, breaks = 100, col = "purple", xlab = "percent.mt", main = "WT2 post filteration")
dev.off()

WT2 <- FindVariableFeatures(WT2, selection.method = "vst", nfeatures = 2500)
tiff(file = "WT2 Variable Genes.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VariableFeaturePlot(WT2)
dev.off()

all.genes <- rownames(WT2)
WT2 <- ScaleData(WT2, features = all.genes)
WT2 <- RunPCA(WT1, features = VariableFeatures(object = WT2))
print(WT2[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WT2, dims = 1:2, reduction = "pca")
DimPlot(WT2, reduction = "pca")
ElbowPlot(WT2, ndims = 50)
tiff(file = "WT2 ElbowPlot.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
ElbowPlot(WT2, ndims = 50)
dev.off()

WT2 <- FindNeighbors(WT2, dims = 1:25)
WT2 <- FindClusters(WT2, resolution = 0.5)
head(Idents(WT2), 5)
WT2 <- RunTSNE(WT2, dims = 1:25)
WT2 <- RunUMAP(WT2, dims = 1:25)
DimPlot(WT2, reduction = "umap", pt.size = 1, label = TRUE)

#### WT Integration ####

Idents(object = WT1) <- "orig.ident"
Idents(object = WT2) <- "orig.ident"

WT1$stim <- "WT"
WT2$stim <- "WT"

WT.anchors <- FindIntegrationAnchors(object.list = list(WT1, WT2), dims = 1:20)

WT <- IntegrateData(anchorset = WT.anchors, dims = 1:20)

DefaultAssay(WT) <- "integrated"

# Run the standard workflow for visualization and clustering
WT <- ScaleData(WT, verbose = FALSE)
WT <- RunPCA(WT, npcs = 30, verbose = FALSE)
ElbowPlot(WT, ndims = 30)
WT <- FindNeighbors(WT, reduction = "pca", dims = 1:26)
WT <- FindClusters(WT, resolution = 0.5)
WT <- RunUMAP(WT, reduction = "pca", dims = 1:26)
WT <- RunTSNE(WT, reduction = "pca", dims = 1:26)

DimPlot(WT, reduction = "umap", pt.size = .5, label = TRUE)

tiff(file = "WT UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("grey")) 
dev.off()

#### HGF_MET WT Integration ####
Idents(object = HGF_MET) <- "orig.ident"
Idents(object = WT) <- "orig.ident"

HGF_MET$stim <- "HGF_MET"
WT$stim <- "WT"

Combined.anchors <- FindIntegrationAnchors(object.list = list(HGF_MET, WT), dims = 1:20)

Combined <- IntegrateData(anchorset = Combined.anchors, dims = 1:20)

DefaultAssay(Combined) <- "integrated"

# Run the standard workflow for visualization and clustering
Combined  <- ScaleData(Combined , verbose = FALSE)
Combined  <- RunPCA(Combined , npcs = 30, verbose = FALSE)
ElbowPlot(Combined , ndims = 30)

Combined  <- FindNeighbors(Combined , reduction = "pca", dims = 1:25)
Combined  <- FindClusters(Combined , resolution = 0.5)
Combined  <- RunUMAP(Combined , reduction = "pca", dims = 1:25)
Combined  <- RunTSNE(Combined , reduction = "pca", dims = 1:25)

DimPlot(Combined, reduction = "umap", pt.size = .5, label = TRUE)

tiff(file = "Combined UMAP.tiff", width = 12, height = 12, units = "in", compression = "lzw", res = 800)
DimPlot(Combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("darkblue", "grey")) 
dev.off()

tiff(file = "Combined UMAP label.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined, reduction = "umap", pt.size = .5, label = TRUE)
dev.off()

DefaultAssay(Combined) <- "RNA"
FeaturePlot(Combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Krt19"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Rgs1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("Pecam1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined, reduction = "umap", features = c("hMETtg"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

tiff(file = "Combined hMETtg Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Hgf Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Hgf"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Ar Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Pbsn Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Epcam Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Epcam"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Vim Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Vim"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Myh11 Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Myh11"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Fbln1 Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Fbln1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Krt8.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined Krt19.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined, reduction = "umap", features = c("Krt19"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, min.cutoff = "q20", max.cutoff = "q90")
dev.off()

Idents(object = Combined) <- "seurat_clusters"
Combined <- RenameIdents(object = Combined,  '0' = "BE", '17' = "BE", '1' = "LE", '4' = "LE", '5' = "LE", '6' = "LE", '8' = "LE", '9' = "LE", '15' = "LE", '3' = "FB", '7' = "FB", '16' = "FB", '19' = "FB", '2' = "SM", '11' = "Immune", '12' = "Immune", '14' = "Immune", '20' = "Immune", '10' = "SV", '13' = "Endothelial", '18' = "Pericyte")
DimPlot(Combined, reduction = "umap", pt.size = 0.3, label = TRUE)
Combined[["CellTypes"]] <- Idents(object = Combined)

tiff(file = "Combined Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = Combined) <- "stim"

Combined <- RenameIdents(object = Combined, 'WT' = "WT", 'HGF_MET' = "HGF_MET")
Combined[["stim"]] <- Idents(object = Combined)

Idents(object = Combined) <- "CellTypes"

tiff(file = "Combined Celltypes UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

Idents(object = Combined) <- "CellTypes"
Combined$stim.CellTypes <- paste(Idents(Combined), Combined$stim, sep = "_")
Idents(object = Combined) <- "stim.CellTypes"
table(Idents(Combined))

Idents(object = Combined) <- "seurat_clusters"
Combined$stim.seurat_clusters <- paste(Idents(Combined), Combined$stim, sep = "_")
Idents(object = Combined) <- "stim.seurat_clusters"
table(Idents(Combined))

Idents(object = Combined) <- "CellTypes"
DefaultAssay(Combined) <- "RNA"
all.genes <- rownames(Combined)
Combined <- ScaleData(Combined, features = all.genes)
Combined.celltype.markers <- FindAllMarkers(Combined, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(Combined.celltype.markers, file = "Combined.celltype.markers.csv")

#Dotplot
Idents(object = Combined) <- "CellTypes"
DotPlot(Combined, features = c("Vtn", "Mustn1", "Cox4i2", "Ndufa4l2", "Rgs5", "Cd93", "Pecam1", "Ctla2a", "Plvap", "Flt1", "Sprr2f", "Cd52", "Wfdc15b", "Svs6", "Pate4", "Tyrobp", "Fcer1g", "Srgn", "Rgs1", "Ccl5", "Myh11", "Cnn1", "Actg2", "Acta1", "Tagln", "Fbln1", "Col1a2", "Apod", "Col1a1", "Col3a1", "Agr2", "Ceacam1", "Tgm4", "Oit1", "Msmb", "Col17a1", "Lgals7", "Krt5", "Krt17", "Krt14", "Pbsn", "hMETtg", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

tiff(file = "Combined CellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(Combined, features = c("Vtn", "Mustn1", "Cox4i2", "Ndufa4l2", "Rgs5", "Cd93", "Pecam1", "Ctla2a", "Plvap", "Flt1", "Sprr2f", "Cd52", "Wfdc15b", "Svs6", "Pate4", "Tyrobp", "Fcer1g", "Srgn", "Rgs1", "Ccl5", "Myh11", "Cnn1", "Actg2", "Acta1", "Tagln", "Fbln1", "Col1a2", "Apod", "Col1a1", "Col3a1", "Agr2", "Ceacam1", "Tgm4", "Oit1", "Msmb", "Col17a1", "Lgals7", "Krt5", "Krt17", "Krt14", "Pbsn", "hMETtg", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()


####Subclustering All_Pro_Epi####

Idents(object = Combined) <- "CellTypes"
Combined_Epi <- subset(Combined, idents = c("BE", "LE"))

DefaultAssay(Combined_Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
Combined_Epi <- ScaleData(Combined_Epi, verbose = FALSE)
Combined_Epi <- RunPCA(Combined_Epi, npcs = 30, verbose = FALSE)
ElbowPlot(Combined_Epi, ndims = 50)
# t-SNE and Clustering
Combined_Epi <- FindNeighbors(Combined_Epi, reduction = "pca", dims = 1:15)
Combined_Epi <- FindClusters(Combined_Epi, resolution = 0.5)
Combined_Epi <- RunTSNE(Combined_Epi, reduction = "pca", dims = 1:15)
Combined_Epi <- RunUMAP(Combined_Epi, reduction = "pca", dims = 1:15)
DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = Combined_Epi) <- "stim"

Combined_Epi <- RenameIdents(object = Combined_Epi, 'WT' = "WT", 'HGF_MET' = "HGF_MET")
Combined_Epi[["stim"]] <- Idents(object = Combined_Epi)

DefaultAssay(Combined_Epi) <- "RNA"
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Krt19"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

tiff(file = "Combined_Epi Ar Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_Epi Krt5 Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_Epi Krt19 Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Krt19"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_Epi Pbsn Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_Epi hMETtg Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("hMETtg"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_Epi Met Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Met"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Combined_Epi Mki67 Exp.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Combined_Epi, reduction = "umap", features = c("Mki67"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()



Idents(object = Combined_Epi) <- "seurat_clusters"
Combined_Epi <- RenameIdents(object = Combined_Epi,  '0' = "BE1", '9' = "BE2", '13' = "BE3", '5' = "LE1", '6' = "LE2", '2' = "LE3", '1' = "LE4", '3' = "LE5", '8' = "LE6", '14' = "LE6", '4' = "LE7", '10' = "UrLE", '7' = "OE1", '11' = "OE2", '12' = "OE3")
DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)
Combined_Epi[["CellTypes"]] <- Idents(object = Combined_Epi)

tiff(file = "Combined_Epi Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "Combined_Epi Celltypes UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

Idents(object = Combined_Epi) <- "CellTypes"
Combined_Epi$stim.CellTypes <- paste(Idents(Combined_Epi), Combined_Epi$stim, sep = "_")
Idents(object = Combined_Epi) <- "stim.CellTypes"
table(Idents(Combined_Epi))

Idents(object = Combined_Epi) <- "CellTypes"
DefaultAssay(Combined_Epi) <- "RNA"
all.genes <- rownames(Combined_Epi)
Combined_Epi <- ScaleData(Combined_Epi, features = all.genes)
Combined_Epi.celltype.markers <- FindAllMarkers(Combined_Epi, min.pct = 0.25, only.pos = TRUE, logfc.threshold = 0.25)
write.csv(Combined_Epi.celltype.markers, file = "Combined_Epi.celltype.markers.csv")

#Dotplot
Idents(object = Combined_Epi) <- "CellTypes"
DotPlot(Combined_Epi, features = c("Fscn1", "Dll1", "Cxcl12", "Egfl6", "Clca3a2", "Lgals1", "Vim", "Col1a2", "Col1a1", "Col3a1", "Cd52", "Wfdc15b", "Acta2", "Hba-a1", "Hbb-bs", "Ppp1r1b", "Barx2", "Gsdmc3", "Clca1", "Gsdmc2", "Slc26a4", "Serpinb11", "Gjb2", "Slc5a8", "Gm7714", "Fam25c", "Gucy2g", "Kcnk3", "Crabp1", "Col6a3", "Pcp4", "Ccdc80", "Chodl", "Derl3", "Bmpr1b", "Ngf", "Cyp2b10", "Fxyd4", "Aldh1a3", "Ly6c1", "Gsdma", "Spink8", "Tcaf2", "A630095E13Rik", "Gm5615", "Fam71a", "Gdf15", "Gm26802", 'Nkx3-1', "Habp2", "Mki67", "Birc5", "Top2a", "Ube2c", "Cdk1", "Gins2", "Mcm3", "Chaf1b", "Hells", "Mcm5", "Krt13", "Ltbp1", "Lama3", "Anxa8", "Aqp3", "Fmo6", "Col16a1", "C2cd4b", "Fmo2", "Smoc2", "Ar" ,"Pbsn", "hMETtg"), cols = c("light grey", "red")) + RotatedAxis()

tiff(file = "Combined_Epi CellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(Combined_Epi, features = c("Fscn1", "Dll1", "Cxcl12", "Egfl6", "Clca3a2", "Lgals1", "Vim", "Col1a2", "Col1a1", "Col3a1", "Cd52", "Wfdc15b", "Acta2", "Hba-a1", "Hbb-bs", "Ppp1r1b", "Barx2", "Gsdmc3", "Clca1", "Gsdmc2", "Slc26a4", "Serpinb11", "Gjb2", "Slc5a8", "Gm7714", "Fam25c", "Gucy2g", "Kcnk3", "Crabp1", "Col6a3", "Pcp4", "Ccdc80", "Chodl", "Derl3", "Bmpr1b", "Ngf", "Cyp2b10", "Fxyd4", "Aldh1a3", "Ly6c1", "Gsdma", "Spink8", "Tcaf2", "A630095E13Rik", "Gm5615", "Fam71a", "Gdf15", "Gm26802", 'Nkx3-1', "Habp2", "Mki67", "Birc5", "Top2a", "Ube2c", "Cdk1", "Gins2", "Mcm3", "Chaf1b", "Hells", "Mcm5", "Krt13", "Ltbp1", "Lama3", "Anxa8", "Aqp3", "Fmo6", "Col16a1", "C2cd4b", "Fmo2", "Smoc2", "Ar" ,"Pbsn", "hMETtg"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#### hMETtg+vs- LE in Combined_Epi ####
Idents(object = Combined_Epi) <- "CellTypes"

Combined_Epi <- RenameIdents(object = Combined_Epi,  'BE1' = "BE", 'BE2' = "BE", 'BE3' = "BE", 'LE1' = "LE", 'LE2' = "LE", 'LE3' = "LE", 'LE4' = "LE", 'LE5' = "LE", 'LE6' = "LE", 'LE7' = "LE", 'LE8' = "LE", 'OE1' = "OE", 'OE2' = "OE", 'UrLE' = "UrLE")
DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)
Combined_Epi[["celltype"]] <- Idents(object = Combined_Epi)


#Add HMETtg info
DefaultAssay(Combined_Epi) <- "RNA"
Combined_EpihMETtgPos <- subset(x=Combined_Epi, subset = hMETtg > 0)
Combined_EpihMETtgNeg <- subset(x=Combined_Epi, subset = hMETtg == 0)
Idents(object = Combined_EpihMETtgPos) <- "hMETtgPos"
Idents(object = Combined_EpihMETtgNeg) <- "hMETtgNeg"
Combined_EpihMETtgPos[["hMETtgExp"]] <- Idents(object = Combined_EpihMETtgPos)
Combined_EpihMETtgNeg[["hMETtgExp"]] <- Idents(object = Combined_EpihMETtgNeg)
Combined_EpihMETtg <- merge(x = Combined_EpihMETtgPos, y = Combined_EpihMETtgNeg)
Idents(object = Combined_EpihMETtg) <- "hMETtgExp"
Combined_Epi$hMETtgExp <- Idents(object = Combined_EpihMETtg)
Idents(object = Combined_Epi) <- "hMETtgExp"

#
Idents(object = Combined_Epi) <- "hMETtgExp"
Combined_Epi$hMETtgExp.celltype <- paste(Idents(Combined_Epi), Combined_Epi$celltype, sep = "_")
Idents(object = Combined_Epi) <- "hMETtgExp.celltype"

DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "#E06666", "#3399FF", "grey", "grey", "grey", "grey"))

tiff(file = "Combined_Epi hMETtg+v- LE UMAP.tiff", width = 8, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Combined_Epi, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "#E06666", "#3399FF", "grey", "grey", "grey", "grey"))
dev.off()

#
Combined_LE <- subset(Combined_Epi, idents = c("hMETtgPos_LE", "hMETtgNeg_LE"))
DefaultAssay(Combined_LE) <- "RNA"
all.genes <- rownames(Combined_LE)
Combined_LE <- ScaleData(Combined_LE, features = all.genes)
Combined_LE.markers <- FindAllMarkers(Combined_LE, min.pct = 0.01, only.pos = TRUE, logfc.threshold = 0.01)
write.csv(Combined_LE.markers, file = "Combined_LE.markers.csv")

#Heatmap
Combined_LE.markersTop50 <- Combined_LE.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
tiff(file = "Combined_LE.markersTop50 Heatmap.tiff", width = 6, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(Combined_LE, features = c(Combined_LE.markersTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Violin
Combined_LE <- RenameIdents(object = Combined_LE, 'hMETtgNeg_LE' = "hMETtgNeg_LE", 'hMETtgPos_LE' = "hMETtgPos_LE")



tiff(file = "Combined_LE hMETtg Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "hMETtg", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Hgf Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Hgf", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Egf Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Egf", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Egfr Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Egfr", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Src Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Src", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Bcar1 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Bcar1", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Igf1r Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Igf1r", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Sox9 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Sox9", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()
tiff(file = "Combined_LE Stat3 Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Combined_LE, features = "Stat3", pt.size = 0, cols = c("#3399FF", "#E06666"))
dev.off()