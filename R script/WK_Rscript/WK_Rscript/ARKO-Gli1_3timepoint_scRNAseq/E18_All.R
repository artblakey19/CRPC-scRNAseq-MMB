setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E17.5")

#Set Current idents
Idents(object = WTE17) <- "seurat_clusters"
Idents(object = ARKOE17) <- "seurat_clusters"

WTE17$stim <- "E17_Ctrl"
ARKOE17$stim <- "E17_ARKO"

E17.Combined.anchors <- FindIntegrationAnchors(object.list = list(WTE17, ARKOE17), dims = 1:20)
E17.Combined <- IntegrateData(anchorset = E17.Combined.anchors, dims = 1:20)
DefaultAssay(E17.Combined) <- "integrated"

#Run the standard workflow for visualization and clustering
E17.Combined <- ScaleData(E17.Combined, verbose = FALSE)
E17.Combined <- RunPCA(E17.Combined, npcs = 50, verbose = FALSE)
ElbowPlot(E17.Combined)

# umap and Clustering
E17.Combined <- FindNeighbors(E17.Combined, reduction = "pca", dims = 1:15)
E17.Combined <- FindClusters(E17.Combined, resolution = 0.5)
E17.Combined <- RunUMAP(E17.Combined, reduction = "pca", dims = 1:15)
E17.Combined <- RunTSNE(E17.Combined, reduction = "pca", dims = 1:15)

DimPlot(E17.Combined, reduction = "umap", pt.size = 0.3, label = TRUE) 


DefaultAssay(E17.Combined) <- "RNA"
FeaturePlot(E17.Combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E17.Combined, reduction = "umap", features = c("Krt15"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E17.Combined, reduction = "umap", features = c("Upk3a"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E17.Combined, reduction = "umap", features = c("Upk3b"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E17.Combined, reduction = "umap", features = c("Krt19"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

E17.combined_UGE <- subset(E17.Combined, idents = c("3", "6"))
Idents(object = E17.combined_UGE) <- "seurat_clusters"

DefaultAssay(E17.combined_UGE) <- "integrated"
E17.combined_UGE <- ScaleData(E17.combined_UGE, verbose = FALSE)
E17.combined_UGE <- RunPCA(E17.combined_UGE, npcs = 50, verbose = FALSE)
ElbowPlot(E17.combined_UGE, ndims = 50)

E17.combined_UGE <- FindNeighbors(E17.combined_UGE, reduction = "pca", dims = 1:16)
E17.combined_UGE <- FindClusters(E17.combined_UGE, resolution = 0.5)
E17.combined_UGE <- RunUMAP(E17.combined_UGE, reduction = "pca", dims = 1:16)
E17.combined_UGE <- RunTSNE(E17.combined_UGE, reduction = "pca", dims = 1:16)
DimPlot(E17.combined_UGE, reduction = "umap", pt.size = 0.3, label = TRUE)

DimPlot(E17.combined_UGE, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)