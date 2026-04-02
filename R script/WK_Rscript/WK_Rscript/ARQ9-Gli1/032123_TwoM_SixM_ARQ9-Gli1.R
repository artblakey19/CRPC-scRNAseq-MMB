#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/ARQ9-Gli1/ARQ9-Gli1_new_NEW")

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

WT_P60.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200108_Tgen_scRNA/Project_COHP_33580_1_X3SC3/outs/filtered_feature_bc_matrix")
WT_P60.unfiltered <- CreateSeuratObject(counts = WT_P60.data,  min.cells = 3, min.features = 500, project = "WT_P60")
WT_P60.unfiltered <- NormalizeData(WT_P60.unfiltered)

####TwoM-ARQ9####
TwoM.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/count_new_3ref/count_Project_COHP_46245_1_X3SC4/outs/filtered_feature_bc_matrix")
TwoM.unfiltered <- CreateSeuratObject(counts = TwoM.data,  min.cells = 3, min.features = 500, project = "ARQ9-Gli1-2M")
TwoM.unfiltered <- NormalizeData(TwoM.unfiltered)

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

TwoM <- subset(TwoM.unfiltered, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 15)

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

TwoM <- FindVariableFeatures(TwoM, selection.method = "vst", nfeatures = 2500)
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

TwoM <- FindNeighbors(TwoM, dims = 1:20)
TwoM <- FindClusters(TwoM, resolution = 0.5)
head(Idents(TwoM), 5)
TwoM <- RunTSNE(TwoM, reduction = "pca", dims = 1:20)
TwoM <- RunUMAP(TwoM, reduction = "pca", dims = 1:20)

Idents(object = TwoM) <- "seurat_clusters"
tiff(file = "TwoM seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3)
dev.off()

DefaultAssay(TwoM) <- "RNA"
tiff(file = "TwoM hARtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("new-ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM mGFP expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM hMETtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("MET"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Myh11 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Fbln1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Cd34 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Cnn1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Cnn1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Des expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Des"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Vim expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Acta2 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "TwoM Lama2 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoM, reduction = "umap", features = c("Lama2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####SixM-ARQ9####
SixM.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/count_new_3ref/count_Project_COHP_46247_2_X3SC4/outs/filtered_feature_bc_matrix")
SixM.unfiltered <- CreateSeuratObject(counts = SixM.data,  min.cells = 3, min.features = 500, project = "ARQ9-Gli1_6M")
SixM.unfiltered <- NormalizeData(SixM.unfiltered)

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

SixM <- subset(SixM.unfiltered, subset = nFeature_RNA > 500 & nFeature_RNA < 7000 & percent.mt < 15)

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
tiff(file = "SixM hARtg expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("new-ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
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
tiff(file = "SixM Myh11 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Fbln1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Cd34 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Cnn1 expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Cnn1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Des expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Des"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "SixM Vim expression plot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixM, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

save.image("mG-Gli1_P35.RData")

write.csv(SixM@assays$RNA@counts, "SixM.counts.csv")