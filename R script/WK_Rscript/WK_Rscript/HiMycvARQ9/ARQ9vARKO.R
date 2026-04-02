#ARQ9vARKO

setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/HiMycvARQ9")

#Setup the seurat object

TwoM.unfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/count_new_3ref/count_Project_COHP_46245_1_X3SC4/outs/filtered_feature_bc_matrix")
TwoM.unfiltered <- CreateSeuratObject(counts = TwoM.unfiltered.data,  min.cells = 3, min.features = 200, project = "ARQ9-2M")

#Initial processing & filtering

TwoM.unfiltered[["percent.mt"]] <- PercentageFeatureSet(TwoM.unfiltered, pattern = "^mt-")
tiff(file = "TwoM Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(TwoM.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "TwoM Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "TwoM Pre-filteration")
dev.off()
tiff(file = "TwoM Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "TwoM Pre-filteration")
dev.off()

plot1 <- FeatureScatter(TwoM.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TwoM.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

TwoM <- subset(TwoM.unfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & percent.mt < 15)
tiff(file = "TwoM Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 200)
VlnPlot(TwoM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
dev.off()
tiff(file = "TwoM Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "TwoM Post-filteration")
dev.off()
tiff(file = "TwoM Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 200)
hist(TwoM@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "TwoM Post-filteration")
dev.off()

plot1 <- FeatureScatter(TwoM, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(TwoM, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

#Normalizing the data
TwoM <- NormalizeData(TwoM)

#Identification of highly variable features
TwoM <- FindVariableFeatures(TwoM, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(TwoM)

#Scaling the data
all.genes <- rownames(TwoM)
TwoM <- ScaleData(TwoM, features = all.genes)

#Perform linear dimensional reduction
TwoM <- RunPCA(TwoM, features = VariableFeatures(object = TwoM))

#Determine the dimensionality of the dataset
ElbowPlot(TwoM, ndims = 50)

#Cluster the cells
TwoM <- FindNeighbors(TwoM, dims = 1:30)
TwoM <- FindClusters(TwoM, resolution = 0.8)

#Run non-linear dimensional reduction
TwoM <- RunUMAP(TwoM, dims = 1:30)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TwoM) <- "stim"
tiff(file = "TwoM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3, cols = c("darkblue")) 
dev.off()

table(Idents(TwoM))

#Load ARKO
ARKOunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/ARKO-Gli1/Single_Cell_Seq/190219_scRNA/27348_ARKO1_count_EGFPmm10/outs/filtered_feature_bc_matrix/")
ARKOunfiltered <- CreateSeuratObject(counts = ARKOunfiltered.data,  min.cells = 3, min.features = 500, project = "ARKOunfiltered")
ARKOunfiltered <- NormalizeData(ARKOunfiltered)
ARKOunfiltered[["percent.mt"]] <- PercentageFeatureSet(ARKOunfiltered, pattern = "^mt-")
ARKO <- subset(ARKOunfiltered, subset = nFeature_RNA > 750 & nFeature_RNA < 7500 & percent.mt < 15)
VlnPlot(ARKO, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)
hist(ARKO@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "ARKO post filteration")
hist(ARKO@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "ARKO post filteration")
ARKO <- FindVariableFeatures(ARKO, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(ARKO)
ARKO <- ScaleData(ARKO, features = all.genes)
ARKO <- RunPCA(ARKO, features = VariableFeatures(object = ARKO))
ElbowPlot(ARKO, ndims = 50)
ARKO <- FindNeighbors(ARKO, reduction = "pca", dims = 1:25)
ARKO <- FindClusters(ARKO, resolution = 0.5)
ARKO <- RunUMAP(ARKO, reduction = "pca", dims = 1:25)
DimPlot(ARKO, reduction = "umap", pt.size = 1, label = TRUE)

####Integration####
#Set Current idents
TwoM[["orig.clusters"]] <- Idents(object = TwoM)
ARKO[["orig.clusters"]] <- Idents(object = ARKO)
Idents(object = TwoM) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"
TwoM$stim <- "ARQ9"
ARKO$stim <- "ARKO"
ARQ9vARKO.anchors <- FindIntegrationAnchors(object.list = list(TwoM, ARKO), dims = 1:20)
ARQ9vARKO.combined <- IntegrateData(anchorset = ARQ9vARKO.anchors, dims = 1:20)
DefaultAssay(ARQ9vARKO.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
ARQ9vARKO.combined <- ScaleData(ARQ9vARKO.combined, verbose = FALSE)
ARQ9vARKO.combined <- RunPCA(ARQ9vARKO.combined, npcs = 30, verbose = FALSE)
ARQ9vARKO.combined <- FindNeighbors(ARQ9vARKO.combined, reduction = "pca", dims = 1:20)
ARQ9vARKO.combined <- FindClusters(ARQ9vARKO.combined, resolution = 0.5)
ARQ9vARKO.combined <- RunUMAP(ARQ9vARKO.combined, reduction = "pca", dims = 1:20)
Idents(object = ARQ9vARKO.combined) <- "stim"
DimPlot(ARQ9vARKO.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 

DefaultAssay(ARQ9vARKO.combined)<-"RNA"
FeaturePlot(ARQ9vARKO.combined, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Cell type identification
DefaultAssay(ARQ9vARKO.combined)<-"RNA"
tiff(file = "ARQ9vARKO.combined celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQ9vARKO.combined, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                 "Fbln1", "Myh11", "Pecam1",
                                                                 "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

Idents(object = ARQ9vARKO.combined) <- "seurat_clusters"
tiff(file = "ARQ9vARKO.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQ9vARKO.combined, reduction = "umap", pt.size = 0.3, label = T) 
dev.off()

####Subclustering FBSM####
Idents(object = ARQ9vARKO.combined) <- "seurat_clusters"
ARQ9vARKO.combined.FBSM <- subset(ARQ9vARKO.combined, idents = c("1", "2"))
ARQ9vARKO.combined.FBSM <- RunUMAP(ARQ9vARKO.combined.FBSM, reduction = "pca", dims = 1:15)
DimPlot(ARQ9vARKO.combined.FBSM, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(ARQ9vARKO.combined.FBSM) <- "integrated"
#Run the standard workflow for visualization and clustering
ARQ9vARKO.combined.FBSM <- ScaleData(ARQ9vARKO.combined.FBSM, verbose = FALSE)
ARQ9vARKO.combined.FBSM <- RunPCA(ARQ9vARKO.combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(ARQ9vARKO.combined.FBSM, ndims = 50)
#Umap and Clustering
DefaultAssay(ARQ9vARKO.combined.FBSM) <- "integrated"
ARQ9vARKO.combined.FBSM <- FindNeighbors(ARQ9vARKO.combined.FBSM, reduction = "pca", dims = 1:10)
ARQ9vARKO.combined.FBSM <- FindClusters(ARQ9vARKO.combined.FBSM, resolution = 2)
ARQ9vARKO.combined.FBSM <- RunUMAP(ARQ9vARKO.combined.FBSM, reduction = "pca", dims = 1:10)

Idents(object = ARQ9vARKO.combined.FBSM) <- "seurat_clusters"
DimPlot(ARQ9vARKO.combined.FBSM, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(ARQ9vARKO.combined.FBSM)<-"RNA"
FeaturePlot(ARQ9vARKO.combined.FBSM, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#CtrlvARKO.combined.FBSM1
Idents(object = ARQ9vARKO.combined.FBSM) <- "seurat_clusters"
ARQ9vARKO.combined.FBSM1 <- subset(ARQ9vARKO.combined.FBSM, idents = c("0", "1", "2", "4", "5", "6", "7", 
                                                                       "8", "9", "10", "12", "13", "14",
                                                                       "15", "16", "18", "19", "20",
                                                                       "21", "22"))
ARQ9vARKO.combined.FBSM1 <- RunUMAP(ARQ9vARKO.combined.FBSM1, reduction = "pca", dims = 1:9)
DimPlot(ARQ9vARKO.combined.FBSM1, reduction = "umap", label = TRUE, split.by = "stim")

DefaultAssay(ARQ9vARKO.combined.FBSM1)<-"RNA"
FeaturePlot(ARQ9vARKO.combined.FBSM1, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM1, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#higher resolution
ARQ9vARKO.combined.FBSM2 <- ARQ9vARKO.combined.FBSM1
DefaultAssay(ARQ9vARKO.combined.FBSM2) <- "integrated"
#Run the standard workflow for visualization and clustering
ARQ9vARKO.combined.FBSM2 <- ScaleData(ARQ9vARKO.combined.FBSM2, verbose = FALSE)
ARQ9vARKO.combined.FBSM2 <- RunPCA(ARQ9vARKO.combined.FBSM2, npcs = 30, verbose = FALSE)
ElbowPlot(ARQ9vARKO.combined.FBSM2, ndims = 50)
#Umap and Clustering
DefaultAssay(ARQ9vARKO.combined.FBSM2) <- "integrated"
ARQ9vARKO.combined.FBSM2 <- FindNeighbors(ARQ9vARKO.combined.FBSM2, reduction = "pca", dims = 1:5)
ARQ9vARKO.combined.FBSM2 <- FindClusters(ARQ9vARKO.combined.FBSM2, resolution = 3)
ARQ9vARKO.combined.FBSM2 <- RunUMAP(ARQ9vARKO.combined.FBSM2, reduction = "pca", dims = 1:5)

DimPlot(ARQ9vARKO.combined.FBSM2, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(ARQ9vARKO.combined.FBSM2)<-"RNA"
FeaturePlot(ARQ9vARKO.combined.FBSM2, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q99")
FeaturePlot(ARQ9vARKO.combined.FBSM2, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM2, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM2, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(ARQ9vARKO.combined.FBSM2, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")


####SCT####
#Load ##ARQ9
TwoM.unfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/220106_WonKyungKim/count_new_3ref/count_Project_COHP_46245_1_X3SC4/outs/filtered_feature_bc_matrix")
TwoM.unfiltered <- CreateSeuratObject(counts = TwoM.unfiltered.data,  min.cells = 3, min.features = 200, project = "ARQ9-2M")
TwoM.unfiltered[["percent.mt"]] <- PercentageFeatureSet(TwoM.unfiltered, pattern = "^mt-")
VlnPlot(TwoM.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(TwoM.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "TwoM Pre-filteration")
hist(TwoM.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "TwoM Pre-filteration")
TwoM <- subset(TwoM.unfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & percent.mt < 15)
VlnPlot(TwoM, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(TwoM@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "TwoM Post-filteration")
hist(TwoM@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "TwoM Post-filteration")
#Apply sctransform normalization
TwoM <- SCTransform(TwoM, vars.to.regress = "percent.mt", verbose = FALSE)

#Perform dimensionality reduction by PCA and UMAP embedding
TwoM <- RunPCA(TwoM, verbose = FALSE)
ElbowPlot(TwoM, ndims = 50)
TwoM <- FindNeighbors(TwoM, dims = 1:25, verbose = FALSE)
TwoM <- FindClusters(TwoM, resolution = 0.5, verbose = FALSE)
TwoM <- RunUMAP(TwoM, dims = 1:25, verbose = FALSE)
DimPlot(TwoM, reduction = "umap", pt.size = 0.3, label = TRUE)

#Stash old idents
TwoM[["orig.clusters"]] <- Idents(object = TwoM)
ARKO[["orig.clusters"]] <- Idents(object = ARKO)

#Set Current idents
Idents(object = TwoM) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"

TwoM$stim <- "ARQ9"
ARKO$stim <- "ARKO"

#Perform integration using pearson residuals
ifnb.list <- list(TwoM = TwoM, ARKO = ARKO)
features <- SelectIntegrationFeatures(object.list = ifnb.list, nfeatures = 3000)
ifnb.list <- PrepSCTIntegration(object.list = ifnb.list, anchor.features = features)
TwoMvARKO.anchors <- FindIntegrationAnchors(object.list = ifnb.list, normalization.method = "SCT",
                                            anchor.features = features)
TWoMvARKO.combined.sct <- IntegrateData(anchorset = TwoMvARKO.anchors, normalization.method = "SCT")

TWoMvARKO.combined.sct <- RunPCA(TWoMvARKO.combined.sct, verbose = FALSE)
ElbowPlot(TWoMvARKO.combined.sct, ndims = 50)
TWoMvARKO.combined.sct <- RunUMAP(TWoMvARKO.combined.sct, reduction = "pca", dims = 1:30, verbose = FALSE)
TWoMvARKO.combined.sct <- FindNeighbors(TWoMvARKO.combined.sct, reduction = "pca", dims = 1:30)
TWoMvARKO.combined.sct <- FindClusters(TWoMvARKO.combined.sct, resolution = 0.5)
DimPlot(TWoMvARKO.combined.sct, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TWoMvARKO.combined.sct) <- "stim"
TWoMvARKO.combined.sct <- RenameIdents(object = TWoMvARKO.combined.sct, 'ARQ9' = "ARQ9", 'ARKO' = "ARKO")
TWoMvARKO.combined.sct[["stim"]] <- Idents(object = TWoMvARKO.combined.sct)
DimPlot(TWoMvARKO.combined.sct, reduction = "umap", pt.size = 0.3, cols = c("darkblue", "grey"))

#Cell type identification
DefaultAssay(TWoMvARKO.combined.sct)<-"SCT"
tiff(file = "TWoMvARKO.combined.sct celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoMvARKO.combined.sct, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                     "Fbln1", "Myh11", "Pecam1",
                                                                     "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

Idents(object = TWoMvARKO.combined.sct) <- "seurat_clusters"
DimPlot(TWoMvARKO.combined.sct, reduction = "umap", pt.size = 0.3, label = T) 


####Subclustering FBSM####
Idents(object = TWoMvARKO.combined.sct) <- "seurat_clusters"
TWoMvARKO.combined.sct.FB <- subset(TWoMvARKO.combined.sct, idents = c("1", "3"))
TWoMvARKO.combined.sct.FB <- RunUMAP(TWoMvARKO.combined.sct.FB, reduction = "pca", dims = 1:10)
DimPlot(TWoMvARKO.combined.sct.FB, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(TWoMvARKO.combined.sct.FB) <- "integrated"
#Run the standard workflow for visualization and clustering
TWoMvARKO.combined.sct.FB <- ScaleData(TWoMvARKO.combined.sct.FB, verbose = FALSE)
TWoMvARKO.combined.sct.FB <- RunPCA(TWoMvARKO.combined.sct.FB, npcs = 30, verbose = FALSE)
ElbowPlot(TWoMvARKO.combined.sct.FB, ndims = 30)
#Umap and Clustering
DefaultAssay(TWoMvARKO.combined.sct.FB) <- "integrated"
TWoMvARKO.combined.sct.FB <- FindNeighbors(TWoMvARKO.combined.sct.FB, reduction = "pca", dims = 1:10)
TWoMvARKO.combined.sct.FB <- FindClusters(TWoMvARKO.combined.sct.FB, resolution = 2)
TWoMvARKO.combined.sct.FB <- RunUMAP(TWoMvARKO.combined.sct.FB, reduction = "pca", dims = 1:10)

Idents(object = TWoMvARKO.combined.sct.FB) <- "seurat_clusters"
DimPlot(TWoMvARKO.combined.sct.FB, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(TWoMvARKO.combined.sct.FB)<-"SCT"
FeaturePlot(TWoMvARKO.combined.sct.FB, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#CtrlvARKO.combined.sct.FB1
Idents(object = TWoMvARKO.combined.sct.FB) <- "seurat_clusters"
TWoMvARKO.combined.sct.FB1 <- subset(TWoMvARKO.combined.sct.FB, idents = c("0", "1", "2", "3", "4", "5", "6", "7", 
                                                                           "8", "9"))
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(TWoMvARKO.combined.sct.FB1) <- "integrated"
#Run the standard workflow for visualization and clustering
TWoMvARKO.combined.sct.FB1 <- ScaleData(TWoMvARKO.combined.sct.FB1, verbose = FALSE)
TWoMvARKO.combined.sct.FB1 <- RunPCA(TWoMvARKO.combined.sct.FB1, npcs = 30, verbose = FALSE)
ElbowPlot(TWoMvARKO.combined.sct.FB1, ndims = 30)
#Umap and Clustering
DefaultAssay(TWoMvARKO.combined.sct.FB1) <- "integrated"
TWoMvARKO.combined.sct.FB1 <- FindNeighbors(TWoMvARKO.combined.sct.FB1, reduction = "pca", dims = 1:10)
TWoMvARKO.combined.sct.FB1 <- FindClusters(TWoMvARKO.combined.sct.FB1, resolution = 2)
TWoMvARKO.combined.sct.FB1 <- RunUMAP(TWoMvARKO.combined.sct.FB1, reduction = "pca", dims = 1:10)
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", label = TRUE, split.by = "stim")

DefaultAssay(TWoMvARKO.combined.sct.FB1)<-"SCT"
FeaturePlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Rename
Idents(object = TWoMvARKO.combined.sct.FB1) <- "seurat_clusters"
TWoMvARKO.combined.sct.FB1 <- RenameIdents(object = TWoMvARKO.combined.sct.FB1, 
                                           '4' = "Telocyte", 
                                           '26' = "FB1",'6' = "FB1",'13' = "FB1",'3' = "FB1",'5' = "FB1",
                                           '24' = "FB2",'20' = "FB2",'23' = "FB2",'8' = "FB2",
                                           '0' = "FB3",'2' = "FB3",'7' = "FB3",'18' = "FB3",
                                           '10' = "FB4",'25' = "FB4",
                                           '21' = "FB4",'1' = "FB4",'11' = "FB4",
                                           '14' = "FB5",'19' = "FB5",'12' = "FB5",'9' = "FB5"
                                           ,'16' = "FB6",'15' = "FB6",'22' = "FB6")
TWoMvARKO.combined.sct.FB1[["FBSMCellTypes1"]] <- Idents(object = TWoMvARKO.combined.sct.FB1)
DimPlot(TWoMvARKO.combined.sct.FB1, split.by = "stim", reduction = "umap", cols = c("red","skyblue1","darkorange","green4", "bisque3","palevioletred3",
                                                                                    "steelblue3","grey"))

Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBSMCellTypes1"
tiff(file = "TWoMvARKO.combined.sct.FB1 FBSMCellTypes1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", cols = c("red","skyblue1","darkorange","green4", "bisque3","palevioletred3",
                                                                 "steelblue3","grey"))
dev.off()
tiff(file = "TWoMvARKO.combined.sct.FB1 FBSMCellTypes1 split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap", split.by="stim", cols = c("red","skyblue1","darkorange","green4", "bisque3","palevioletred3",
                                                                                  "steelblue3","grey"))
dev.off()

Idents(object = TWoMvARKO.combined.sct.FB1) <- "stim"
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap")
TWoM.combined.sct.FB1 <- subset(TWoMvARKO.combined.sct.FB1, idents = c("Ctrl"))
ARKO.combined.sct.FB1 <- subset(TWoMvARKO.combined.sct.FB1, idents = c("ARKO"))


DefaultAssay(TWoM.combined.sct.FB1) <- "SCT"
tiff(file = "TWoM.combined.sct.FB1 PDGFRa.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoM.combined.sct.FB1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TWoM.combined.sct.FB1 Cd34.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoM.combined.sct.FB1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TWoM.combined.sct.FB1 Foxl1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoM.combined.sct.FB1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q50")
dev.off()
tiff(file = "TWoM.combined.sct.FB1 Ar.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoM.combined.sct.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "TWoM.combined.sct.FB1 Gli1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoM.combined.sct.FB1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q70")
dev.off()
tiff(file = "TWoM.combined.sct.FB1 EGFP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TWoM.combined.sct.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()

DefaultAssay(ARKO.combined.sct.FB1) <- "SCT"
tiff(file = "ARKO.combined.sct.FB1 PDGFRa.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Cd34.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Foxl1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Ar.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 Gli1.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "ARKO.combined.sct.FB1 EGFP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKO.combined.sct.FB1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q99", min.cutoff = 2)
dev.off()

#split
Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBSMCellTypes1"
TWoMvARKO.combined.sct.FB1 <- RenameIdents(object = TWoMvARKO.combined.sct.FB1, 
                                           'Telocyte' = "Telocyte", 'FB1' = "FB",'FB2' = "FB",'FB3' = "FB",
                                           'FB4' = "FB",'FB5' = "FB", "FB6" = "FB")
TWoMvARKO.combined.sct.FB1[["FBvTC"]] <- Idents(object = TWoMvARKO.combined.sct.FB1)
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap")

Idents(object = TWoMvARKO.combined.sct.FB1) <- "stim"
TWoMvARKO.combined.sct.FB1 <- RenameIdents(object = TWoMvARKO.combined.sct.FB1, 
                                           'Ctrl' = "Ctrl", 'ARKO' = "ARKO")
TWoMvARKO.combined.sct.FB1[["stim"]] <- Idents(object = TWoMvARKO.combined.sct.FB1)
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap")

DefaultAssay(TWoMvARKO.combined.sct.FB1) <- "SCT"
Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBvTC"
tiff(file = "TWoMvARKO.combined.sct.FB1 vlnplot EGFP split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TWoMvARKO.combined.sct.FB1, features = "EGFP", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "TWoMvARKO.combined.sct.FB1 vlnplot Foxl1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TWoMvARKO.combined.sct.FB1, features = "Foxl1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "TWoMvARKO.combined.sct.FB1 vlnplot Pdgfra split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TWoMvARKO.combined.sct.FB1, features = "Pdgfra", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "TWoMvARKO.combined.sct.FB1 vlnplot Gli1 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TWoMvARKO.combined.sct.FB1, features = "Gli1", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "TWoMvARKO.combined.sct.FB1 vlnplot CD34 split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TWoMvARKO.combined.sct.FB1, features = "Cd34", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()
tiff(file = "TWoMvARKO.combined.sct.FB1 vlnplot Ar split.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(TWoMvARKO.combined.sct.FB1, features = "Ar", pt.size = 0, split.by = "stim", split.plot = TRUE)
dev.off()

#Cell counts
Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBvTC"
TWoMvARKO.combined.sct.FB1$FBvTC.stim <- paste(Idents(TWoMvARKO.combined.sct.FB1), TWoMvARKO.combined.sct.FB1$stim, sep = "_")
Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBvTC.stim"
table(Idents(TWoMvARKO.combined.sct.FB1))


Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBvTC.stim"
TWoMvARKO.combined.sct.FB1 <- RenameIdents(object = TWoMvARKO.combined.sct.FB1, 
                                           'Telocyte_Ctrl' = "Telocyte_Ctrl", 'FB_Ctrl' = "FB_Ctrl",
                                           'Telocyte_ARKO' = "Telocyte_ARKO",'FB_ARKO' = "FB_ARKO")
TWoMvARKO.combined.sct.FB1[["FBvTC.stim"]] <- Idents(object = TWoMvARKO.combined.sct.FB1)
DimPlot(TWoMvARKO.combined.sct.FB1, reduction = "umap")


#Dotplot
Idents(object = TWoMvARKO.combined.sct.FB1) <- "FBvTC.stim"
DefaultAssay(TWoMvARKO.combined.sct.FB1) <- "SCT"