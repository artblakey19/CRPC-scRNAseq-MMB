setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E17.5")

#clustering
Idents(object = E17.combined) <- "seurat_clusters"
DefaultAssay(E17.combined) <- "integrated"
E17.combined <- ScaleData(E17.combined, verbose = FALSE)
E17.combined <- RunPCA(E17.combined, npcs = 50, verbose = FALSE)
ElbowPlot(E17.combined, ndims = 50)
E17.combined <- FindNeighbors(E17.combined, reduction = "pca", dims = 1:20)
E17.combined <- FindClusters(E17.combined, resolution = 0.5)
E17.combined <- RunUMAP(E17.combined, reduction = "pca", dims = 1:20)
E17.combined <- RunTSNE(E17.combined, reduction = "pca", dims = 1:20)
DimPlot(E17.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E17.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E17.combined, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

DefaultAssay(E17.combined) <- "RNA"
tiff(file = "E17.combined Cdh1 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(E17.combined, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### subclustering Epi ####

E17.combined_Epi <- subset(E17.combined, idents = c("2", "9", "12", "14", "17"))

Idents(object = E17.combined_Epi) <- "seurat_clusters"
DefaultAssay(E17.combined_Epi) <- "integrated"
E17.combined_Epi <- ScaleData(E17.combined_Epi, verbose = FALSE)
E17.combined_Epi <- RunPCA(E17.combined_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(E17.combined_Epi, ndims = 50)
E17.combined_Epi <- FindNeighbors(E17.combined_Epi, reduction = "pca", dims = 1:16)
E17.combined_Epi <- FindClusters(E17.combined_Epi, resolution = 0.5)
E17.combined_Epi <- RunUMAP(E17.combined_Epi, reduction = "pca", dims = 1:16)
E17.combined_Epi <- RunTSNE(E17.combined_Epi, reduction = "pca", dims = 1:16)
DimPlot(E17.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E17.combined_Epi UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E17.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

#DEGs
Idents(object = E17.combined_Epi) <- "seurat_clusters"
DefaultAssay(E17.combined_Epi) <- "RNA"
E17.combined_Epi <- ScaleData(E17.combined_Epi, features = rownames(E17.combined_Epi))
E17.combined_Epi.markers <- FindAllMarkers(E17.combined_Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(E17.combined_Epi.markers, "E17.combined_Epi.markers.csv")


