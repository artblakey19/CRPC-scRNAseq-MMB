#Setup the seurat object

WT_P60.unfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200108_Tgen_scRNA/33580_count_EGFPmm10_v3/outs/filtered_feature_bc_matrix")
WT_P60.unfiltered <- CreateSeuratObject(counts = WT_P60.unfiltered.data,  min.cells = 3, min.features = 200, project = "WT-P60")
WT_P60.unfiltered[["percent.mt"]] <- PercentageFeatureSet(WT_P60.unfiltered, pattern = "^mt-")
VlnPlot(WT_P60.unfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(WT_P60.unfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "WT_P60 Pre-filteration")
hist(WT_P60.unfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "WT_P60 Pre-filteration")

plot1 <- FeatureScatter(WT_P60.unfiltered, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(WT_P60.unfiltered, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

WT_P60 <- subset(WT_P60.unfiltered, subset = nFeature_RNA > 700 & nFeature_RNA < 7000 & percent.mt < 15)

#Normalizing the data
WT_P60 <- NormalizeData(WT_P60)

#Identification of highly variable features
WT_P60 <- FindVariableFeatures(WT_P60, selection.method = "vst", nfeatures = 2000)
VariableFeaturePlot(WT_P60)

#Scaling the data
all.genes <- rownames(WT_P60)
WT_P60 <- ScaleData(WT_P60, features = all.genes)

#Perform linear dimensional reduction
WT_P60 <- RunPCA(WT_P60, features = VariableFeatures(object = WT_P60))

#Determine the dimensionality of the dataset
ElbowPlot(WT_P60, ndims = 50)

#Cluster the cells
WT_P60 <- FindNeighbors(WT_P60, dims = 1:27)
WT_P60 <- FindClusters(WT_P60, resolution = 0.8)

#Run non-linear dimensional reduction
WT_P60 <- RunUMAP(WT_P60, dims = 1:27)
DimPlot(WT_P60, reduction = "umap", pt.size = 0.3, label = TRUE)

#Load HiMyc
##HiMYC
Ctrlunfiltered.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Hi-myc-ARKO/2nd Run 6-4-2020/Sun Lab Analysis/Project_COHP_36595_1_XCTC1_count/outs/filtered_feature_bc_matrix")
Ctrlunfiltered <- CreateSeuratObject(counts = Ctrlunfiltered.data,  min.cells = 3, min.features = 500, project = "Ctrlunfiltered")
Ctrlunfiltered <- NormalizeData(Ctrlunfiltered)
Ctrlunfiltered[["percent.mt"]] <- PercentageFeatureSet(Ctrlunfiltered, pattern = "^mt-")
VlnPlot(Ctrlunfiltered, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrlunfiltered@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Pre-filteration")
hist(Ctrlunfiltered@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Pre-filteration")
Ctrl <- subset(Ctrlunfiltered, subset = nFeature_RNA > 600 & nFeature_RNA < 6500 & percent.mt < 15)
VlnPlot(Ctrl, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
hist(Ctrl@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "Ctrl Post-filteration")
hist(Ctrl@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "Ctrl Post-filteration")
Ctrl <- FindVariableFeatures(Ctrl, selection.method = "vst", nfeatures = 5000)
VariableFeaturePlot(Ctrl)
Ctrl <- ScaleData(Ctrl, verbose = FALSE)
Ctrl <- RunPCA(Ctrl, npcs = 30, verbose = FALSE)
ElbowPlot(Ctrl, ndims = 50)
Ctrl <- FindNeighbors(Ctrl, reduction = "pca", dims = 1:20)
Ctrl <- FindClusters(Ctrl, resolution = 0.5)
Ctrl <- RunUMAP(Ctrl, reduction = "pca", dims = 1:20)
DimPlot(Ctrl, reduction = "umap", pt.size = 0.3)

####Integration####
#Set Current idents
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
WT_P60[["orig.clusters"]] <- Idents(object = WT_P60)
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = WT_P60) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
WT_P60$stim <- "WT_P60"
CtrlvWT.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, WT_P60), dims = 1:20)
CtrlvWT.combined <- IntegrateData(anchorset = CtrlvWT.anchors, dims = 1:20)
DefaultAssay(CtrlvWT.combined) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvWT.combined <- ScaleData(CtrlvWT.combined, verbose = FALSE)
CtrlvWT.combined <- RunPCA(CtrlvWT.combined, npcs = 30, verbose = FALSE)
CtrlvWT.combined <- FindNeighbors(CtrlvWT.combined, reduction = "pca", dims = 1:20)
CtrlvWT.combined <- FindClusters(CtrlvWT.combined, resolution = 0.5)
CtrlvWT.combined <- RunUMAP(CtrlvWT.combined, reduction = "pca", dims = 1:20)
Idents(object = CtrlvWT.combined) <- "stim"
DimPlot(CtrlvWT.combined, reduction = "umap", pt.size = 0.3, cols = c("light grey", "darkblue")) 

#Cell type identification
DefaultAssay(CtrlvWT.combined)<-"RNA"
tiff(file = "CtrlvWT.combined celltype marker expression plots.tiff", width = 20, height = 20, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvWT.combined, reduction = "umap", features = c("Krt5", "Krt4", "Krt8", "Cd24a", "Plp1",
                                                                 "Fbln1", "Myh11", "Pecam1",
                                                                 "Vim", "Tyrobp", "Ccl5", "Pate4"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

Idents(object = CtrlvWT.combined) <- "seurat_clusters"
tiff(file = "CtrlvWT.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvWT.combined, reduction = "umap", pt.size = 0.3, label = T) 
dev.off()

####Subclustering FBSM####
Idents(object = CtrlvWT.combined) <- "seurat_clusters"
CtrlvWT.combined.FBSM <- subset(CtrlvWT.combined, idents = c("1", "3"))
CtrlvWT.combined.FBSM <- RunUMAP(CtrlvWT.combined.FBSM, reduction = "pca", dims = 1:10)
DimPlot(CtrlvWT.combined.FBSM, reduction = "umap", label = TRUE, split.by = "stim")

#higher resolution
DefaultAssay(CtrlvWT.combined.FBSM) <- "integrated"
#Run the standard workflow for visualization and clustering
CtrlvWT.combined.FBSM <- ScaleData(CtrlvWT.combined.FBSM, verbose = FALSE)
CtrlvWT.combined.FBSM <- RunPCA(CtrlvWT.combined.FBSM, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvWT.combined.FBSM, ndims = 50)
#Umap and Clustering
DefaultAssay(CtrlvWT.combined.FBSM) <- "integrated"
CtrlvWT.combined.FBSM <- FindNeighbors(CtrlvWT.combined.FBSM, reduction = "pca", dims = 1:12)
CtrlvWT.combined.FBSM <- FindClusters(CtrlvWT.combined.FBSM, resolution = 2)
CtrlvWT.combined.FBSM <- RunUMAP(CtrlvWT.combined.FBSM, reduction = "pca", dims = 1:12)

Idents(object = CtrlvWT.combined.FBSM) <- "seurat_clusters"
DimPlot(CtrlvWT.combined.FBSM, reduction = "umap", split.by="stim", label = TRUE)

DefaultAssay(CtrlvWT.combined.FBSM)<-"RNA"
FeaturePlot(CtrlvWT.combined.FBSM, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvWT.combined.FBSM, reduction = "umap", features = c("Cd34"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvWT.combined.FBSM, reduction = "umap", features = c("Foxl1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvWT.combined.FBSM, reduction = "umap", features = c("Gli1"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(CtrlvWT.combined.FBSM, reduction = "umap", features = c("EGFP"),split.by = "stim",  cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

Idents(object = CtrlvWT.combined.FBSM) <- "seurat_clusters"
CtrlvWT.combined.FBSM1 <- subset(CtrlvWT.combined.FBSM, idents = c("1", "2", "3", "4", "5", "6", "7", 
                                                                        "8", "9", "10", "11", "12", "13", "14",
                                                                        "15", "16", "17", "18", "19", 
                                                                        "21", "22", "24", "25", "26", "27", "28", "29"))
DimPlot(CtrlvWT.combined.FBSM1, reduction = "umap", label = TRUE, split.by = "stim")




#Rename
Idents(object = CtrlvWT.combined.FBSM) <- "seurat_clusters"
CtrlvWT.combined.FBSM <- RenameIdents(object = CtrlvWT.combined.FBSM, 
                                         '7' = "Telocyte", '12' = "Telocyte",
                                         '1' = "FB",'2' = "FB",'3' = "FB",'4' = "FB",
                                         '5' = "FB",'6' = "FB",'7' = "FB",'8' = "FB",'9' = "FB",
                                         '10' = "FB",'11' = "FB",'12' = "FB",'13' = "FB",'14' = "FB",'15' = "FB",
                                         '16' = "FB",'17' = "FB",'18' = "FB",'19' = "FB",'21' = "FB",'22' = "FB",
                                         '24' = "FB", '25' = "FB",'26' = "FB",'27' = "FB",'28' = "FB",
                                         '29' = "FB")
CtrlvWT.combined.FBSM[["FBSMCellTypes1"]] <- Idents(object = CtrlvWT.combined.FBSM)
DimPlot(CtrlvWT.combined.FBSM, split.by = "stim", reduction = "umap", cols = c("red","blue","darkorange","skyblue1", "bisque3","palevioletred3",
                                                                                  "green4"))
#Dotplot
Idents(object = CtrlvWT.combined.FBSM) <- "celltype"
tiff(file = "CtrlvWT.combined.FBSM TC markers dotplot.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DotPlot(CtrlvWT.combined.FBSM, features = c("Pdgfra", "Foxl1", "Gli1", "Cd34"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()