####2M####

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

####P35 WT####

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

####Sub-clustering Epi####

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
tiff(file = "TwoMvWT.combined.Epi1 Epicelltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Epi1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Fbln1", "Myh11", "Ppp1r1b"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()