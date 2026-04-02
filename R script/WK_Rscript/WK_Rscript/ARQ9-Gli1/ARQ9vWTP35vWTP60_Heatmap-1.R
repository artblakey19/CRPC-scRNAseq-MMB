####WT_P60####

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

#Setup workspace to make file calling & saving easy
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35_P60")

WTP60.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191214_32848_33580_test/counts/33580_WT/outs/filtered_feature_bc_matrix")
WTP60unfiltered <- CreateSeuratObject(counts = WTP60.data,  min.cells = 3, min.features = 200, project = "WT_P60")
WTP60unfiltered <- NormalizeData(WTP60unfiltered)

WTP60unfiltered[["percent.mt"]] <- PercentageFeatureSet(WTP60unfiltered, pattern = "^mt-")

WT_P60 <- subset(WTP60unfiltered, subset = nFeature_RNA > 250 & nFeature_RNA < 4000 & percent.mt < 15)

WT_P60 <- FindVariableFeatures(WT_P60, selection.method = "vst", nfeatures = 2500)
VariableFeaturePlot(WT_P60)

#Clustering
all.genes <- rownames(WT_P60)
WT_P60 <- ScaleData(WT_P60, features = all.genes)
WT_P60 <- RunPCA(WT_P60, features = VariableFeatures(object = WT_P60))
print(WT_P60[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(WT_P60, dims = 1:2, reduction = "pca")
DimPlot(WT_P60, reduction = "pca")
ElbowPlot(WT_P60, ndims = 50)


WT_P60 <- FindNeighbors(WT_P60, dims = 1:24)
WT_P60 <- FindClusters(WT_P60, resolution = 0.5)
WT_P60 <- RunTSNE(WT_P60, reduction = "pca", dims = 1:24)
WT_P60 <- RunUMAP(WT_P60, reduction = "pca", dims = 1:24)

Idents(object = WT_P60) <- "seurat_clusters"
DimPlot(WT_P60, reduction = "umap", pt.size = 0.3)

#### Merging Datasets TwoMvWT-P35vWT-P60 ####

#Stash old idents
TwoM[["orig.clusters"]] <- Idents(object = TwoM)
WT[["orig.clusters"]] <- Idents(object = WT)
WT_P60[["orig.clusters"]] <- Idents(object = WT_P60)

#Set Current idents
Idents(object = TwoM) <- "seurat_clusters"
Idents(object = WT) <- "seurat_clusters"
Idents(object = WT_P60) <- "seurat_clusters"

TwoM$stim <- "ARQ9_2M"
WT$stim <- "WT_P35"
WT_P60$stim <- "WT_P60"

TwoMvWT.anchors <- FindIntegrationAnchors(object.list = list(TwoM, WT, WT_P60), dims = 1:20)
TwoMvWT.combined <- IntegrateData(anchorset = TwoMvWT.anchors, dims = 1:20)

DefaultAssay(TwoMvWT.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
TwoMvWT.combined <- ScaleData(TwoMvWT.combined, verbose = FALSE)
TwoMvWT.combined <- RunPCA(TwoMvWT.combined, npcs = 30, verbose = FALSE)
ElbowPlot(TwoMvWT.combined, ndims = 50)

#Umap and Clustering
TwoMvWT.combined <- FindNeighbors(TwoMvWT.combined, reduction = "pca", dims = 1:16)
TwoMvWT.combined <- FindClusters(TwoMvWT.combined, resolution = 0.5)
TwoMvWT.combined <- RunTSNE(TwoMvWT.combined, reduction = "pca", dims = 1:16)
TwoMvWT.combined <- RunUMAP(TwoMvWT.combined, reduction = "pca", dims = 1:16)
DimPlot(TwoMvWT.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

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

TwoMvWT.combined1 <- FindNeighbors(TwoMvWT.combined1, reduction = "pca", dims = 1:20)
TwoMvWT.combined1 <- FindClusters(TwoMvWT.combined1, resolution = 0.5)
TwoMvWT.combined1 <- RunUMAP(TwoMvWT.combined1, reduction = "pca", dims = 1:20)
TwoMvWT.combined1 <- RunTSNE(TwoMvWT.combined1, reduction = "pca", dims = 1:20)
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

DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell type identification
DefaultAssay(TwoMvWT.combined1) <- "RNA"
tiff(file = "TwoMvWT.combined1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2",
                                                                "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                                "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

####Re-clustering Stro####

Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
TwoMvWT.combined.Stro <- subset(TwoMvWT.combined1, idents = c("1", "5", "6", "14", "16", "15", "11", "13"))
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
DimPlot(TwoMvWT.combined.Stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

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
TwoMvWT.combined.Stro1 <- RenameIdents(object = TwoMvWT.combined.Stro1, '1' = "FB1", '4' = "FB2", '0' = "FB3", '2' = "FB4", '15' = "FB5", 
                                       '3' = "SM1", '6' = "SM2", '7' = "VE", '13' = "Pericyte", '10' = "Glia", '14' = "Glia", '5' = "Leu",
                                       '12' = "Leu", '8' = "Lym", '11' = "Lym", '9' = "OF")
TwoMvWT.combined.Stro1[["StroCellType"]] <- Idents(object = TwoMvWT.combined.Stro1)

#Umap
Idents(object = TwoMvWT.combined.Stro1) <- "StroCellType"
tiff(file = "TwoMvWT.combined.Stro1 StroCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Stro1, reduction = "umap", pt.size = 0.3)
dev.off()

####Subset FB####
Idents(object = TwoMvWT.combined.Stro1) <- "StroCellType"

Idents(object = TwoMvWT.combined.Stro1) <- "StroCellType"
TwoMvWT.combined.Stro1 <- RenameIdents(object = TwoMvWT.combined.Stro1, 'FB1' = "FB1", 'FB4' = "FB2", 'FB2' = "FB3", 'FB3' = "FB4", 'FB5' = "FB5", 
                                       'SM1' = "SM1", 'SM2' = "SM2", 'VE' = "VE", 'Pericyte' = "Pericyte", 'Glia' = "Glia", 'Leu' = "Leu",
                                       'Lym' = "Lym",'OF' = "OF")
TwoMvWT.combined.Stro1[["StroCellType"]] <- Idents(object = TwoMvWT.combined.Stro1)

TwoMvWT.combined.FB <- subset(TwoMvWT.combined.Stro1, idents = c("FB1", "FB2", "FB3", "FB4", "FB5"))
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
TwoMvWT.combined.FB <- RenameIdents(object = TwoMvWT.combined.FB, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5", 'WT_P35_FB1' = "WT_P35_FB1", 'WT_P35_FB2' = "WT_P35_FB2", 'WT_P35_FB3' = "WT_P35_FB3", 'WT_P35_FB4' = "WT_P35_FB4", 'WT_P35_FB5' = "WT_P35_FB5",
                                    'WT_P60_FB1' = "WT_P60_FB1", 'WT_P60_FB2' = "WT_P60_FB2", 'WT_P60_FB3' = "WT_P60_FB3", 'WT_P60_FB4' = "WT_P60_FB4", 'WT_P60_FB5' = "WT_P60_FB5")
TwoMvWT.combined.FB[["stim.StroCellType"]] <- Idents(object = TwoMvWT.combined.FB)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3)

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

#Heatmap_AllFB
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
TwoMvWT.combined.FB.markers <- FindAllMarkers(TwoMvWT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FBTop50 <- TwoMvWT.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB Heatmap Top50 purple-1.tiff", width = 10, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
tiff(file = "TwoMvWT.combined.FB stim Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
TwoMvWT.combined.FB.markers.1 <- FindAllMarkers(TwoMvWT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FBTop50.1 <- TwoMvWT.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB stim individual Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()


####TwoMvWT_P60####
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoMvWTP60.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("ARQ9_2M", "WT_P60"))


#Heatmap_AllFB
DefaultAssay(TwoMvWTP60.combined.FB) <- "RNA"
Idents(object = TwoMvWTP60.combined.FB) <- "StroCellType"
TwoMvWTP60.combined.FB <- ScaleData(TwoMvWTP60.combined.FB, features = rownames(TwoMvWTP60.combined.FB))
TwoMvWTP60.combined.FB.markers <- FindAllMarkers(TwoMvWTP60.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP60.combined.FBTop50 <- TwoMvWTP60.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP60.combined.FB Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP60.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP60.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) 
dev.off()

write.csv(TwoMvWTP60.combined.FBTop50, "TwoMvWTP60.combined.FBTop50.csv")

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWTP60.combined.FB) <- "stim.StroCellType"
DimPlot(TwoMvWTP60.combined.FB, reduction = "umap", pt.size = 0.3)
TwoMvWTP60.combined.FB <- RenameIdents(object = TwoMvWTP60.combined.FB, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5", 'WT_P35_FB1' = "WT_P35_FB1", 'WT_P35_FB2' = "WT_P35_FB2", 'WT_P35_FB3' = "WT_P35_FB3", 'WT_P35_FB4' = "WT_P35_FB4", 'WT_P35_FB5' = "WT_P35_FB5",
                                    'WT_P60_FB1' = "WT_P60_FB1", 'WT_P60_FB2' = "WT_P60_FB2", 'WT_P60_FB3' = "WT_P60_FB3", 'WT_P60_FB4' = "WT_P60_FB4", 'WT_P60_FB5' = "WT_P60_FB5")
TwoMvWTP60.combined.FB[["stim.StroCellType"]] <- Idents(object = TwoMvWTP60.combined.FB)
DimPlot(TwoMvWTP60.combined.FB, reduction = "umap", pt.size = 0.3)

DefaultAssay(TwoMvWTP60.combined.FB) <- "RNA"
TwoMvWTP60.combined.FB <- ScaleData(TwoMvWTP60.combined.FB, features = rownames(TwoMvWTP60.combined.FB))
tiff(file = "TwoMvWTP60.combined.FB stim Heatmap Top50 purple.tiff", width = 16, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP60.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP60.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWTP60.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWTP60.combined.FB) <- "RNA"
TwoMvWTP60.combined.FB <- ScaleData(TwoMvWTP60.combined.FB, features = rownames(TwoMvWTP60.combined.FB))
TwoMvWTP60.combined.FB.markers.1 <- FindAllMarkers(TwoMvWTP60.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP60.combined.FBTop50.1 <- TwoMvWTP60.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP60.combined.FB stim individual Heatmap Top50 purple.tiff", width = 16, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP60.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP60.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) 
dev.off()

write.csv(TwoMvWTP60.combined.FBTop50.1, "TwoMvWTP60.combined.FBTop50.1.csv")

#Heatmap_AllFB
DefaultAssay(TwoMvWTP60.combined.FB) <- "RNA"
Idents(object = TwoMvWTP60.combined.FB) <- "StroCellType"
TwoMvWTP60.combined.FB <- ScaleData(TwoMvWTP60.combined.FB, features = rownames(TwoMvWTP60.combined.FB))
TwoMvWTP60.combined.FB.markers <- FindAllMarkers(TwoMvWTP60.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP60.combined.FBTop50 <- TwoMvWTP60.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP60.combined.FB Heatmap Top50 purple-1.tiff", width = 10, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP60.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP60.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 20))
dev.off()

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWTP60.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWTP60.combined.FB) <- "RNA"
TwoMvWTP60.combined.FB <- ScaleData(TwoMvWTP60.combined.FB, features = rownames(TwoMvWTP60.combined.FB))
tiff(file = "TwoMvWTP60.combined.FB stim Heatmap Top50 purple-1.tiff", width = 20, height = 70, units = "in", compression = "lzw", res = 300)
DoHeatmap(TwoMvWTP60.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP60.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 18))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWTP60.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWTP60.combined.FB) <- "RNA"
TwoMvWTP60.combined.FB <- ScaleData(TwoMvWTP60.combined.FB, features = rownames(TwoMvWTP60.combined.FB))
TwoMvWTP60.combined.FB.markers.1 <- FindAllMarkers(TwoMvWTP60.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP60.combined.FBTop50.1 <- TwoMvWTP60.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP60.combined.FB stim individual Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 500)
DoHeatmap(TwoMvWTP60.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP60.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

####TwoMvWT_P35####
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoMvWTP35.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("ARQ9_2M", "WT_P35"))

#Heatmap_AllFB
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "StroCellType"
TwoMvWTP35.combined.FB <- ScaleData(TwoMvWTP35.combined.FB, features = rownames(TwoMvWTP35.combined.FB))
TwoMvWTP35.combined.FB.markers <- FindAllMarkers(TwoMvWTP35.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP35.combined.FBTop50 <- TwoMvWTP35.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP35.combined.FB Heatmap Top50 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP35.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP35.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoMvWTP35.combined.FBTop50, "TwoMvWTP35.combined.FBTop50.csv")

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWTP35.combined.FB) <- "stim.StroCellType"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 0.3)
TwoMvWTP35.combined.FB <- RenameIdents(object = TwoMvWTP35.combined.FB, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5", 'WT_P35_FB1' = "WT_P35_FB1", 'WT_P35_FB2' = "WT_P35_FB2", 'WT_P35_FB3' = "WT_P35_FB3", 'WT_P35_FB4' = "WT_P35_FB4", 'WT_P35_FB5' = "WT_P35_FB5")
TwoMvWTP35.combined.FB[["stim.StroCellType"]] <- Idents(object = TwoMvWTP35.combined.FB)
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 0.3)

DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
TwoMvWTP35.combined.FB <- ScaleData(TwoMvWTP35.combined.FB, features = rownames(TwoMvWTP35.combined.FB))
tiff(file = "TwoMvWTP35.combined.FB stim Heatmap Top50 purple.tiff", width = 16, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP35.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP35.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWTP35.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
TwoMvWTP35.combined.FB <- ScaleData(TwoMvWTP35.combined.FB, features = rownames(TwoMvWTP35.combined.FB))
TwoMvWTP35.combined.FB.markers.1 <- FindAllMarkers(TwoMvWTP35.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP35.combined.FBTop50.1 <- TwoMvWTP35.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP35.combined.FB stim individual Heatmap Top50 purple.tiff", width = 16, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP35.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP35.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

write.csv(TwoMvWTP35.combined.FBTop50.1, "TwoMvWTP35.combined.FBTop50.1.csv")

#Heatmap_AllFB
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "StroCellType"
TwoMvWTP35.combined.FB <- ScaleData(TwoMvWTP35.combined.FB, features = rownames(TwoMvWTP35.combined.FB))
TwoMvWTP35.combined.FB.markers <- FindAllMarkers(TwoMvWTP35.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP35.combined.FBTop50 <- TwoMvWTP35.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP35.combined.FB Heatmap Top50 purple-1.tiff", width = 10, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP35.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP35.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWTP35.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
TwoMvWTP35.combined.FB <- ScaleData(TwoMvWTP35.combined.FB, features = rownames(TwoMvWTP35.combined.FB))
tiff(file = "TwoMvWTP35.combined.FB stim Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP35.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP35.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWTP35.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
TwoMvWTP35.combined.FB <- ScaleData(TwoMvWTP35.combined.FB, features = rownames(TwoMvWTP35.combined.FB))
TwoMvWTP35.combined.FB.markers.1 <- FindAllMarkers(TwoMvWTP35.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWTP35.combined.FBTop50.1 <- TwoMvWTP35.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWTP35.combined.FB stim individual Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWTP35.combined.FB, features = c("ARQ", "EGFP", TwoMvWTP35.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

####mGFP+ARQ+vsothers####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35_P60/mGFP+")

#Add EGFP info
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FBEGFPNeg <- subset(x=TwoMvWT.combined.FB, subset = EGFP < 1)
TwoMvWT.combined.FBEGFPPos <- subset(x=TwoMvWT.combined.FB, subset = EGFP >= 1)
Idents(object = TwoMvWT.combined.FBEGFPNeg) <- "EGFPNeg"
Idents(object = TwoMvWT.combined.FBEGFPPos) <- "EGFPPos"
TwoMvWT.combined.FBEGFPNeg[["EGFPExp"]] <- Idents(object = TwoMvWT.combined.FBEGFPNeg)
TwoMvWT.combined.FBEGFPPos[["EGFPExp"]] <- Idents(object = TwoMvWT.combined.FBEGFPPos)
TwoMvWT.combined.FBEGFP <- merge(x = TwoMvWT.combined.FBEGFPNeg, y = TwoMvWT.combined.FBEGFPPos)
Idents(object = TwoMvWT.combined.FBEGFP) <- "EGFPExp"
TwoMvWT.combined.FB$EGFPExp <- Idents(object = TwoMvWT.combined.FBEGFP)

#Add ARQ info
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FBARQNeg <- subset(x=TwoMvWT.combined.FB, subset = ARQ == 0)
TwoMvWT.combined.FBARQPos <- subset(x=TwoMvWT.combined.FB, subset = ARQ > 0)
Idents(object = TwoMvWT.combined.FBARQNeg) <- "ARQNeg"
Idents(object = TwoMvWT.combined.FBARQPos) <- "ARQPos"
TwoMvWT.combined.FBARQNeg[["ARQExp"]] <- Idents(object = TwoMvWT.combined.FBARQNeg)
TwoMvWT.combined.FBARQPos[["ARQExp"]] <- Idents(object = TwoMvWT.combined.FBARQPos)
TwoMvWT.combined.FBARQ <- merge(x = TwoMvWT.combined.FBARQNeg, y = TwoMvWT.combined.FBARQPos)
Idents(object = TwoMvWT.combined.FBARQ) <- "ARQExp"
TwoMvWT.combined.FB$ARQExp <- Idents(object = TwoMvWT.combined.FBARQ)

#Cellcounts
Idents(object = TwoMvWT.combined.FB) <- "EGFPExp."
TwoMvWT.combined.FB$EGFPExp.StroCellType <- paste(Idents(TwoMvWT.combined.FB), TwoMvWT.combined.FB$StroCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB) <- "EGFPExp"
TwoMvWT.combined.FB$EGFPExp.ARQExp <- paste(Idents(TwoMvWT.combined.FB), TwoMvWT.combined.FB$ARQExp, sep = "_")
Idents(object = TwoMvWT.combined.FB) <- "EGFPExp.ARQExp"
TwoMvWT.combined.FB$EGFPExp.ARQExp.StroCellType <- paste(Idents(TwoMvWT.combined.FB), TwoMvWT.combined.FB$StroCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB) <- "EGFPExp.ARQExp.StroCellType"
TwoMvWT.combined.FB$EGFPExp.ARQExp.StroCellType.stim <- paste(Idents(TwoMvWT.combined.FB), TwoMvWT.combined.FB$stim, sep = "_")
Idents(object = TwoMvWT.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"

#DEGs_TwoM_mGFP+ARQ+FB1vsmGFP+ARQ+FB2-5
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoM.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("ARQ9_2M"))

DefaultAssay(TwoM.combined.FB) <- "RNA"
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp.StroCellType"
DimPlot(TwoM.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB1vFB2_5.0.Markers <- FindMarkers(TwoM.combined.FB, ident.1 = "EGFPPos_ARQPos_FB1", ident.2 = c("EGFPPos_ARQPos_FB2", "EGFPPos_ARQPos_FB3",
                                                                                                            "EGFPPos_ARQPos_FB4", "EGFPPos_ARQPos_FB5"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB1vFB2_5.0.Markers, "TwoM_EGFPPosARQPos_FB1vFB2_5.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB1vFB2_5 <- read.csv("TwoM_EGFPPosARQPos_FB1vFB2_5.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB1vFB2_5_pvalue <- TwoM_EGFPPosARQPos_FB1vFB2_5$p_val
TwoM_EGFPPosARQPos_FB1vFB2_5_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB1vFB2_5_pvalue)
TwoM_EGFPPosARQPos_FB1vFB2_5_BH = p.adjust(TwoM_EGFPPosARQPos_FB1vFB2_5_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB1vFB2_5_BH, "TwoM_EGFPPosARQPos_FB1vFB2_5_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB1vsmGFP+ARQ+FB2-5
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp.StroCellType"
DimPlot(TwoM.combined.FB, reduction = "umap", pt.size = 1)

TwoM.combined.EGFPARQFB <- subset(TwoM.combined.FB, idents = c("EGFPPos_ARQPos_FB1", "EGFPPos_ARQPos_FB2", "EGFPPos_ARQPos_FB3", "EGFPPos_ARQPos_FB4", "EGFPPos_ARQPos_FB5"))
TwoM.combined.EGFPARQFB <- RenameIdents(object = TwoM.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB1' = "EGFPPos_ARQPos_FB1", 'EGFPPos_ARQPos_FB2' = "EGFPPos_ARQPos_FB2_3_4_5", 'EGFPPos_ARQPos_FB3' = "EGFPPos_ARQPos_FB2_3_4_5", 'EGFPPos_ARQPos_FB4' = "EGFPPos_ARQPos_FB2_3_4_5", 'EGFPPos_ARQPos_FB5' = "EGFPPos_ARQPos_FB2_3_4_5")
TwoM.combined.EGFPARQFB[["FB1vFB2345"]] <- Idents(object = TwoM.combined.EGFPARQFB)

DefaultAssay(TwoM.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoM.combined.EGFPARQFB) <- "FB1vFB2345"
TwoM.combined.EGFPARQFB <- ScaleData(TwoM.combined.EGFPARQFB, features = rownames(TwoM.combined.EGFPARQFB))
TwoM.combined.EGFPARQFB.markers <- FindAllMarkers(TwoM.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoM.combined.EGFPARQFB.markers.Top50 <- TwoM.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB1vsmGFP+ARQ+FB2-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.combined.EGFPARQFB, features = c(TwoM.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#DEGs_TwoM_mGFP+ARQ+FB2vsmGFP+ARQ+FB1-5
TwoM_EGFPPosARQPos_FB2vFB1_5.0.Markers <- FindMarkers(TwoM.combined.FB, ident.1 = "EGFPPos_ARQPos_FB2", ident.2 = c("EGFPPos_ARQPos_FB1", "EGFPPos_ARQPos_FB3",
                                                                                                                    "EGFPPos_ARQPos_FB4", "EGFPPos_ARQPos_FB5"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB2vFB1_5.0.Markers, "TwoM_EGFPPosARQPos_FB2vFB1_5.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB2vFB1_5 <- read.csv("TwoM_EGFPPosARQPos_FB2vFB1_5.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB2vFB1_5_pvalue <- TwoM_EGFPPosARQPos_FB2vFB1_5$p_val
TwoM_EGFPPosARQPos_FB2vFB1_5_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB2vFB1_5_pvalue)
TwoM_EGFPPosARQPos_FB2vFB1_5_BH = p.adjust(TwoM_EGFPPosARQPos_FB2vFB1_5_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB2vFB1_5_BH, "TwoM_EGFPPosARQPos_FB2vFB1_5_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB2vsmGFP+ARQ+FB1-5
Idents(object = TwoM.combined.EGFPARQFB) <- "EGFPExp.ARQExp.StroCellType"
TwoM.combined.EGFPARQFB <- RenameIdents(object = TwoM.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB2' = "EGFPPos_ARQPos_FB2", 'EGFPPos_ARQPos_FB1' = "EGFPPos_ARQPos_FB1_3_4_5", 'EGFPPos_ARQPos_FB3' = "EGFPPos_ARQPos_FB1_3_4_5", 'EGFPPos_ARQPos_FB4' = "EGFPPos_ARQPos_FB1_3_4_5", 'EGFPPos_ARQPos_FB5' = "EGFPPos_ARQPos_FB1_3_4_5")
TwoM.combined.EGFPARQFB[["FB2vFB1345"]] <- Idents(object = TwoM.combined.EGFPARQFB)

DefaultAssay(TwoM.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoM.combined.EGFPARQFB) <- "FB2vFB1345"
TwoM.combined.EGFPARQFB <- ScaleData(TwoM.combined.EGFPARQFB, features = rownames(TwoM.combined.EGFPARQFB))
TwoM.combined.EGFPARQFB.markers <- FindAllMarkers(TwoM.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoM.combined.EGFPARQFB.markers.Top50 <- TwoM.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB2vsmGFP+ARQ+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.combined.EGFPARQFB, features = c(TwoM.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#DEGs_TwoM_mGFP+ARQ+FB3vsmGFP+ARQ+FB1-5
TwoM_EGFPPosARQPos_FB3vFB1_5.0.Markers <- FindMarkers(TwoM.combined.FB, ident.1 = "EGFPPos_ARQPos_FB3", ident.2 = c("EGFPPos_ARQPos_FB1", "EGFPPos_ARQPos_FB2",
                                                                                                                    "EGFPPos_ARQPos_FB4", "EGFPPos_ARQPos_FB5"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB3vFB1_5.0.Markers, "TwoM_EGFPPosARQPos_FB3vFB1_5.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB3vFB1_5 <- read.csv("TwoM_EGFPPosARQPos_FB3vFB1_5.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB3vFB1_5_pvalue <- TwoM_EGFPPosARQPos_FB3vFB1_5$p_val
TwoM_EGFPPosARQPos_FB3vFB1_5_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB3vFB1_5_pvalue)
TwoM_EGFPPosARQPos_FB3vFB1_5_BH = p.adjust(TwoM_EGFPPosARQPos_FB3vFB1_5_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB3vFB1_5_BH, "TwoM_EGFPPosARQPos_FB3vFB1_5_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB3vsmGFP+ARQ+FB1-5
Idents(object = TwoM.combined.EGFPARQFB) <- "EGFPExp.ARQExp.StroCellType"
TwoM.combined.EGFPARQFB <- RenameIdents(object = TwoM.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB3' = "EGFPPos_ARQPos_FB3", 'EGFPPos_ARQPos_FB1' = "EGFPPos_ARQPos_FB1_2_4_5", 'EGFPPos_ARQPos_FB2' = "EGFPPos_ARQPos_FB1_2_4_5", 'EGFPPos_ARQPos_FB4' = "EGFPPos_ARQPos_FB1_2_4_5", 'EGFPPos_ARQPos_FB5' = "EGFPPos_ARQPos_FB1_2_4_5")
TwoM.combined.EGFPARQFB[["FB3vFB1245"]] <- Idents(object = TwoM.combined.EGFPARQFB)

DefaultAssay(TwoM.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoM.combined.EGFPARQFB) <- "FB3vFB1245"
TwoM.combined.EGFPARQFB <- ScaleData(TwoM.combined.EGFPARQFB, features = rownames(TwoM.combined.EGFPARQFB))
TwoM.combined.EGFPARQFB.markers <- FindAllMarkers(TwoM.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoM.combined.EGFPARQFB.markers.Top50 <- TwoM.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB3vsmGFP+ARQ+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.combined.EGFPARQFB, features = c(TwoM.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#DEGs_TwoM_mGFP+ARQ+FB4vsmGFP+ARQ+FB1-5
TwoM_EGFPPosARQPos_FB4vFB1_5.0.Markers <- FindMarkers(TwoM.combined.FB, ident.1 = "EGFPPos_ARQPos_FB4", ident.2 = c("EGFPPos_ARQPos_FB1", "EGFPPos_ARQPos_FB3",
                                                                                                                    "EGFPPos_ARQPos_FB2", "EGFPPos_ARQPos_FB5"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB4vFB1_5.0.Markers, "TwoM_EGFPPosARQPos_FB4vFB1_5.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB4vFB1_5 <- read.csv("TwoM_EGFPPosARQPos_FB4vFB1_5.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB4vFB1_5_pvalue <- TwoM_EGFPPosARQPos_FB4vFB1_5$p_val
TwoM_EGFPPosARQPos_FB4vFB1_5_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB4vFB1_5_pvalue)
TwoM_EGFPPosARQPos_FB4vFB1_5_BH = p.adjust(TwoM_EGFPPosARQPos_FB4vFB1_5_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB4vFB1_5_BH, "TwoM_EGFPPosARQPos_FB4vFB1_5_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB4vsmGFP+ARQ+FB1-5
Idents(object = TwoM.combined.EGFPARQFB) <- "EGFPExp.ARQExp.StroCellType"
TwoM.combined.EGFPARQFB <- RenameIdents(object = TwoM.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB4' = "EGFPPos_ARQPos_FB4", 'EGFPPos_ARQPos_FB1' = "EGFPPos_ARQPos_FB1_2_3_5", 'EGFPPos_ARQPos_FB2' = "EGFPPos_ARQPos_FB1_2_3_5", 'EGFPPos_ARQPos_FB3' = "EGFPPos_ARQPos_FB1_2_3_5", 'EGFPPos_ARQPos_FB5' = "EGFPPos_ARQPos_FB1_2_3_5")
TwoM.combined.EGFPARQFB[["FB4vFB1235"]] <- Idents(object = TwoM.combined.EGFPARQFB)

DefaultAssay(TwoM.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoM.combined.EGFPARQFB) <- "FB4vFB1235"
TwoM.combined.EGFPARQFB <- ScaleData(TwoM.combined.EGFPARQFB, features = rownames(TwoM.combined.EGFPARQFB))
TwoM.combined.EGFPARQFB.markers <- FindAllMarkers(TwoM.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoM.combined.EGFPARQFB.markers.Top50 <- TwoM.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB4vsmGFP+ARQ+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.combined.EGFPARQFB, features = c(TwoM.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#DEGs_TwoM_mGFP+ARQ+FB5vsmGFP+ARQ+FB1-4
TwoM_EGFPPosARQPos_FB5vFB1_4.0.Markers <- FindMarkers(TwoM.combined.FB, ident.1 = "EGFPPos_ARQPos_FB5", ident.2 = c("EGFPPos_ARQPos_FB1", "EGFPPos_ARQPos_FB3",
                                                                                                                    "EGFPPos_ARQPos_FB2", "EGFPPos_ARQPos_FB4"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB5vFB1_4.0.Markers, "TwoM_EGFPPosARQPos_FB5vFB1_4.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB5vFB1_4 <- read.csv("TwoM_EGFPPosARQPos_FB5vFB1_4.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB5vFB1_4_pvalue <- TwoM_EGFPPosARQPos_FB5vFB1_4$p_val
TwoM_EGFPPosARQPos_FB5vFB1_4_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB5vFB1_4_pvalue)
TwoM_EGFPPosARQPos_FB5vFB1_4_BH = p.adjust(TwoM_EGFPPosARQPos_FB5vFB1_4_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB5vFB1_4_BH, "TwoM_EGFPPosARQPos_FB5vFB1_4_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB5vsmGFP+ARQ+FB1-4
Idents(object = TwoM.combined.EGFPARQFB) <- "EGFPExp.ARQExp.StroCellType"
TwoM.combined.EGFPARQFB <- RenameIdents(object = TwoM.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB5' = "EGFPPos_ARQPos_FB5", 'EGFPPos_ARQPos_FB1' = "EGFPPos_ARQPos_FB1_2_3_4", 'EGFPPos_ARQPos_FB2' = "EGFPPos_ARQPos_FB1_2_3_4", 'EGFPPos_ARQPos_FB3' = "EGFPPos_ARQPos_FB1_2_3_4", 'EGFPPos_ARQPos_FB4' = "EGFPPos_ARQPos_FB1_2_3_4")
TwoM.combined.EGFPARQFB[["FB5vFB1234"]] <- Idents(object = TwoM.combined.EGFPARQFB)

DefaultAssay(TwoM.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoM.combined.EGFPARQFB) <- "FB5vFB1234"
TwoM.combined.EGFPARQFB <- ScaleData(TwoM.combined.EGFPARQFB, features = rownames(TwoM.combined.EGFPARQFB))
TwoM.combined.EGFPARQFB.markers <- FindAllMarkers(TwoM.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoM.combined.EGFPARQFB.markers.Top50 <- TwoM.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB5vsmGFP+ARQ+FB1-4 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoM.combined.EGFPARQFB, features = c(TwoM.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#Cell count
Idents(object = TwoM.combined.FB) <- "EGFPExp.ARQExp.StroCellType"
table(Idents(TwoM.combined.FB))


#DEGs_TwoM_mGFP+ARQ+FB1vsWTP35_mGFP+FB1-5
Idents(object = TwoMvWT.combined.FB) <- "stim"
TwoMvWTP35.combined.FB <- subset(TwoMvWT.combined.FB, idents = c("ARQ9_2M", "WT_P35"))

DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB.0.Markers <- FindMarkers(TwoMvWTP35.combined.FB, ident.1 = "EGFPPos_ARQPos_FB1_ARQ9_2M", ident.2 = c("EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35",
                                                                                                                    "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB.0.Markers, "TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB <- read.csv("TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_pvalue <- TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB$p_val
TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_pvalue)
TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_BH = p.adjust(TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_BH, "TwoM_EGFPPosARQPos_FB1vWT_EGFPPos_FB_BH.csv")

#DEGs_TwoM_mGFP+ARQ+FB1vsWTP35_mGFP+FB1-5
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB.0.Markers <- FindMarkers(TwoMvWTP35.combined.FB, ident.1 = "EGFPPos_ARQPos_FB2_ARQ9_2M", ident.2 = c("EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35",
                                                                                                                                          "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB.0.Markers, "TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB <- read.csv("TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_pvalue <- TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB$p_val
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_pvalue)
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_BH = p.adjust(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_BH, "TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB_BH.csv")

#DEGs_TwoM_mGFP+ARQ+FB3vsWTP35_mGFP+FB1-5
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB.0.Markers <- FindMarkers(TwoMvWTP35.combined.FB, ident.1 = "EGFPPos_ARQPos_FB3_ARQ9_2M", ident.2 = c("EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35",
                                                                                                                                          "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB.0.Markers, "TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB <- read.csv("TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_pvalue <- TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB$p_val
TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_pvalue)
TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_BH = p.adjust(TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_BH, "TwoM_EGFPPosARQPos_FB3vWT_EGFPPos_FB_BH.csv")

#DEGs_TwoM_mGFP+ARQ+FB4vsWTP35_mGFP+FB1-5
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB.0.Markers <- FindMarkers(TwoMvWTP35.combined.FB, ident.1 = "EGFPPos_ARQPos_FB4_ARQ9_2M", ident.2 = c("EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35",
                                                                                                                                          "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB.0.Markers, "TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB <- read.csv("TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_pvalue <- TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB$p_val
TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_pvalue)
TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_BH = p.adjust(TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_BH, "TwoM_EGFPPosARQPos_FB4vWT_EGFPPos_FB_BH.csv")

#DEGs_TwoM_mGFP+ARQ+FB5vsWTP35_mGFP+FB1-5
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB.0.Markers <- FindMarkers(TwoMvWTP35.combined.FB, ident.1 = "EGFPPos_ARQPos_FB5_ARQ9_2M", ident.2 = c("EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35",
                                                                                                                                          "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB.0.Markers, "TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB <- read.csv("TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_pvalue <- TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB$p_val
TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_pvalue)
TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_BH = p.adjust(TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_BH, "TwoM_EGFPPosARQPos_FB5vWT_EGFPPos_FB_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB1vsWTP35_mGFP+FB1-5


Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)
TwoMvWT.combined.EGFPARQFB <- subset(TwoMvWTP35.combined.FB, idents = c("EGFPPos_ARQPos_FB1_ARQ9_2M", "EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35", "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"))
TwoMvWT.combined.EGFPARQFB <- RenameIdents(object = TwoMvWT.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB1_ARQ9_2M' = "EGFPPos_ARQPos_FB1_ARQ9_2M", 'EGFPPos_ARQNeg_FB1_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB2_WT_P35' = "EGFPPos_ARQNeg_FB_WT"
                                           , 'EGFPPos_ARQNeg_FB3_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB4_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB5_WT_P35' = "EGFPPos_ARQNeg_FB_WT")
TwoMvWT.combined.EGFPARQFB[["TwoMvWT"]] <- Idents(object = TwoMvWT.combined.EGFPARQFB)

DefaultAssay(TwoMvWT.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoMvWT.combined.EGFPARQFB) <- "TwoMvWT"
TwoMvWT.combined.EGFPARQFB <- ScaleData(TwoMvWT.combined.EGFPARQFB, features = rownames(TwoMvWT.combined.EGFPARQFB))
TwoMvWT.combined.EGFPARQFB.markers <- FindAllMarkers(TwoMvWT.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.EGFPARQFB.markers.Top50 <- TwoMvWT.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB1vsWTP35_mGFP+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFPARQFB, features = c(TwoMvWT.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#DEGs_TwoM_mGFP+ARQ+FB2vsWTP35_mGFP+FB2
DefaultAssay(TwoMvWTP35.combined.FB) <- "RNA"
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)

TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2.0.Markers <- FindMarkers(TwoMvWTP35.combined.FB, ident.1 = "EGFPPos_ARQPos_FB2_ARQ9_2M", ident.2 = c("EGFPPos_ARQNeg_FB2_WT_P35"
                                                                                                                                          ), min.pct = 0, logfc.threshold = 0)
write.csv(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2.0.Markers, "TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2.0.Markers.csv")

#p.adjust
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2 <- read.csv("TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2.0.Markers.csv") 
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_pvalue <- TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2$p_val
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_pvalue=as.numeric(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_pvalue)
TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_BH = p.adjust(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_pvalue, "BH")
write.csv(TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_BH, "TwoM_EGFPPosARQPos_FB2vWT_EGFPPos_FB2_BH.csv")

#Heatmap_TwoM_mGFP+ARQ+FB2vsWTP35_mGFP+FB1-5
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)
TwoMvWT.combined.EGFPARQFB <- subset(TwoMvWTP35.combined.FB, idents = c("EGFPPos_ARQPos_FB2_ARQ9_2M", "EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35", "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"))
TwoMvWT.combined.EGFPARQFB <- RenameIdents(object = TwoMvWT.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB2_ARQ9_2M' = "EGFPPos_ARQPos_FB2_ARQ9_2M", 'EGFPPos_ARQNeg_FB1_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB2_WT_P35' = "EGFPPos_ARQNeg_FB_WT"
                                           , 'EGFPPos_ARQNeg_FB3_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB4_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB5_WT_P35' = "EGFPPos_ARQNeg_FB_WT")
TwoMvWT.combined.EGFPARQFB[["TwoMvWT"]] <- Idents(object = TwoMvWT.combined.EGFPARQFB)

DefaultAssay(TwoMvWT.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoMvWT.combined.EGFPARQFB) <- "TwoMvWT"
TwoMvWT.combined.EGFPARQFB <- ScaleData(TwoMvWT.combined.EGFPARQFB, features = rownames(TwoMvWT.combined.EGFPARQFB))
TwoMvWT.combined.EGFPARQFB.markers <- FindAllMarkers(TwoMvWT.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.EGFPARQFB.markers.Top50 <- TwoMvWT.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB2vsWTP35_mGFP+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFPARQFB, features = c(TwoMvWT.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#Heatmap_TwoM_mGFP+ARQ+FB3vsWTP35_mGFP+FB1-5
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)
TwoMvWT.combined.EGFPARQFB <- subset(TwoMvWTP35.combined.FB, idents = c("EGFPPos_ARQPos_FB3_ARQ9_2M", "EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35", "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"))
TwoMvWT.combined.EGFPARQFB <- RenameIdents(object = TwoMvWT.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB3_ARQ9_2M' = "EGFPPos_ARQPos_FB3_ARQ9_2M", 'EGFPPos_ARQNeg_FB1_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB2_WT_P35' = "EGFPPos_ARQNeg_FB_WT"
                                           , 'EGFPPos_ARQNeg_FB3_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB4_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB5_WT_P35' = "EGFPPos_ARQNeg_FB_WT")
TwoMvWT.combined.EGFPARQFB[["TwoMvWT"]] <- Idents(object = TwoMvWT.combined.EGFPARQFB)

DefaultAssay(TwoMvWT.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoMvWT.combined.EGFPARQFB) <- "TwoMvWT"
TwoMvWT.combined.EGFPARQFB <- ScaleData(TwoMvWT.combined.EGFPARQFB, features = rownames(TwoMvWT.combined.EGFPARQFB))
TwoMvWT.combined.EGFPARQFB.markers <- FindAllMarkers(TwoMvWT.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.EGFPARQFB.markers.Top50 <- TwoMvWT.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB3vsWTP35_mGFP+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFPARQFB, features = c(TwoMvWT.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#Heatmap_TwoM_mGFP+ARQ+FB4vsWTP35_mGFP+FB1-5
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)
TwoMvWT.combined.EGFPARQFB <- subset(TwoMvWTP35.combined.FB, idents = c("EGFPPos_ARQPos_FB4_ARQ9_2M", "EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35", "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"))
TwoMvWT.combined.EGFPARQFB <- RenameIdents(object = TwoMvWT.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB4_ARQ9_2M' = "EGFPPos_ARQPos_FB3_ARQ9_2M", 'EGFPPos_ARQNeg_FB1_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB2_WT_P35' = "EGFPPos_ARQNeg_FB_WT"
                                           , 'EGFPPos_ARQNeg_FB3_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB4_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB5_WT_P35' = "EGFPPos_ARQNeg_FB_WT")
TwoMvWT.combined.EGFPARQFB[["TwoMvWT"]] <- Idents(object = TwoMvWT.combined.EGFPARQFB)

DefaultAssay(TwoMvWT.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoMvWT.combined.EGFPARQFB) <- "TwoMvWT"
TwoMvWT.combined.EGFPARQFB <- ScaleData(TwoMvWT.combined.EGFPARQFB, features = rownames(TwoMvWT.combined.EGFPARQFB))
TwoMvWT.combined.EGFPARQFB.markers <- FindAllMarkers(TwoMvWT.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.EGFPARQFB.markers.Top50 <- TwoMvWT.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB4vsWTP35_mGFP+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFPARQFB, features = c(TwoMvWT.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

#Heatmap_TwoM_mGFP+ARQ+FB5vsWTP35_mGFP+FB1-5
Idents(object = TwoMvWTP35.combined.FB) <- "EGFPExp.ARQExp.StroCellType.stim"
DimPlot(TwoMvWTP35.combined.FB, reduction = "umap", pt.size = 1)
TwoMvWT.combined.EGFPARQFB <- subset(TwoMvWTP35.combined.FB, idents = c("EGFPPos_ARQPos_FB5_ARQ9_2M", "EGFPPos_ARQNeg_FB1_WT_P35", "EGFPPos_ARQNeg_FB2_WT_P35", "EGFPPos_ARQNeg_FB3_WT_P35", "EGFPPos_ARQNeg_FB4_WT_P35", "EGFPPos_ARQNeg_FB5_WT_P35"))
TwoMvWT.combined.EGFPARQFB <- RenameIdents(object = TwoMvWT.combined.EGFPARQFB, 'EGFPPos_ARQPos_FB5_ARQ9_2M' = "EGFPPos_ARQPos_FB5_ARQ9_2M", 'EGFPPos_ARQNeg_FB1_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB2_WT_P35' = "EGFPPos_ARQNeg_FB_WT"
                                           , 'EGFPPos_ARQNeg_FB3_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB4_WT_P35' = "EGFPPos_ARQNeg_FB_WT", 'EGFPPos_ARQNeg_FB5_WT_P35' = "EGFPPos_ARQNeg_FB_WT")
TwoMvWT.combined.EGFPARQFB[["TwoMvWT"]] <- Idents(object = TwoMvWT.combined.EGFPARQFB)

DefaultAssay(TwoMvWT.combined.EGFPARQFB) <- "RNA"
Idents(object = TwoMvWT.combined.EGFPARQFB) <- "TwoMvWT"
TwoMvWT.combined.EGFPARQFB <- ScaleData(TwoMvWT.combined.EGFPARQFB, features = rownames(TwoMvWT.combined.EGFPARQFB))
TwoMvWT.combined.EGFPARQFB.markers <- FindAllMarkers(TwoMvWT.combined.EGFPARQFB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.EGFPARQFB.markers.Top50 <- TwoMvWT.combined.EGFPARQFB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "Heatmap_TwoM_mGFP+ARQ+FB5vsWTP35_mGFP+FB1-5 Top50 purple.tiff", width = 10, height = 25, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.EGFPARQFB, features = c(TwoMvWT.combined.EGFPARQFB.markers.Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 12))
dev.off()

####Re-clustering Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/2M/WT_P35_P60/Epi")

Idents(object = TwoMvWT.combined1) <- "seurat_clusters"
DimPlot(TwoMvWT.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
DefaultAssay(TwoMvWT.combined1) <- "RNA"
FeaturePlot(TwoMvWT.combined1, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

TwoMvWT.combined.Epi <- subset(TwoMvWT.combined1, idents = c("0", "10", "12", "7", "9", "2", "4", "17", "8", "3"))
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
DimPlot(TwoMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)

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
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = TwoMvWT.combined.Epi1) <- "Phase"
tiff(file = "TwoMvWT.combined.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "TwoMvWT.combined.Epi1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "TwoMvWT.combined.Epi1 seurat LABEL UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()

#Rename
DefaultAssay(TwoMvWT.combined.Epi1) <- "RNA"
tiff(file = "TwoMvWT.combined.Epi1 expression plots.tiff", width = 18, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(TwoMvWT.combined.Epi1, reduction = "umap", features = c("Krt5", "Krt14", "Trp63", "Krt4", "Ppp1r1b", "Krt7", "Krt19", "Tgm4", "Msmb", "Svs2"), cols = c("light grey", "red"), pt.size = 0.3)
dev.off()

Idents(object = TwoMvWT.combined.Epi1) <- "seurat_clusters"
TwoMvWT.combined.Epi1 <- RenameIdents(object = TwoMvWT.combined.Epi1, '0' = "BE1", '10' = "BE2", '6' = "BE3", '12' = "BE4", '13' = "UrLE", 
                                      '7' = "LumP", '2' = "LE1", '3' = "LE2", '8' = "LE3", '15' = "LE3", '5' = "LE4", '1' = "LE5", '11' = "LE6",
                                      '4' = "SV", '9' = "OE1", '14' = "OE2")
TwoMvWT.combined.Epi1[["EpiCellType"]] <- Idents(object = TwoMvWT.combined.Epi1)

#Umap
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
tiff(file = "TwoMvWT.combined.Epi1 EpiCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()

tiff(file = "TwoMvWT.combined.Epi1 EpiCellType stim UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.Epi1) <- "stim"
TwoMvWT.combined.Epi1$stim.EpiCellType <- paste(Idents(TwoMvWT.combined.Epi1), TwoMvWT.combined.Epi1$EpiCellType, sep = "_")
Idents(object = TwoMvWT.combined.Epi1) <- "stim.EpiCellType"
table(Idents(TwoMvWT.combined.Epi1))


####Subset BE####
Idents(object = TwoMvWT.combined.Epi1) <- "EpiCellType"
TwoMvWT.combined.BE <- subset(TwoMvWT.combined.Epi1, idents = c("BE1", "BE2", "BE3"))
TwoMvWT.combined.BE <- RunTSNE(TwoMvWT.combined.BE, reduction = "pca", dims = 1:20)
TwoMvWT.combined.BE <- RunUMAP(TwoMvWT.combined.BE, reduction = "pca", dims = 1:20)
DimPlot(TwoMvWT.combined.BE, reduction = "umap", pt.size = 0.3)

tiff(file = "TwoMvWT.combined.BE UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.BE, reduction = "umap", pt.size = 1)
dev.off()

tiff(file = "TwoMvWT.combined.BE split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TwoMvWT.combined.BE, reduction = "umap", pt.size = 1, split.by = "stim")
dev.off()

#Cell counts
Idents(object = TwoMvWT.combined.BE) <- "stim"
TwoMvWT.combined.BE$stim.EpiCellType <- paste(Idents(TwoMvWT.combined.BE), TwoMvWT.combined.BE$EpiCellType, sep = "_")
Idents(object = TwoMvWT.combined.FB) <- "stim.EpiCellType"
table(Idents(TwoMvWT.combined.FB))

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
TwoMvWT.combined.FB <- RenameIdents(object = TwoMvWT.combined.FB, 'ARQ9_2M_FB1' = "ARQ9_2M_FB1", 'ARQ9_2M_FB2' = "ARQ9_2M_FB2", 'ARQ9_2M_FB3' = "ARQ9_2M_FB3", 'ARQ9_2M_FB4' = "ARQ9_2M_FB4", 'ARQ9_2M_FB5' = "ARQ9_2M_FB5", 'WT_P35_FB1' = "WT_P35_FB1", 'WT_P35_FB2' = "WT_P35_FB2", 'WT_P35_FB3' = "WT_P35_FB3", 'WT_P35_FB4' = "WT_P35_FB4", 'WT_P35_FB5' = "WT_P35_FB5",
                                    'WT_P60_FB1' = "WT_P60_FB1", 'WT_P60_FB2' = "WT_P60_FB2", 'WT_P60_FB3' = "WT_P60_FB3", 'WT_P60_FB4' = "WT_P60_FB4", 'WT_P60_FB5' = "WT_P60_FB5")
TwoMvWT.combined.FB[["stim.StroCellType"]] <- Idents(object = TwoMvWT.combined.FB)
DimPlot(TwoMvWT.combined.FB, reduction = "umap", pt.size = 0.3)

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

#Heatmap_AllFB
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
Idents(object = TwoMvWT.combined.FB) <- "StroCellType"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
TwoMvWT.combined.FB.markers <- FindAllMarkers(TwoMvWT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FBTop50 <- TwoMvWT.combined.FB.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB Heatmap Top50 purple-1.tiff", width = 10, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
tiff(file = "TwoMvWT.combined.FB stim Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Heatmap_AllFB_combined_stim_individual
Idents(object = TwoMvWT.combined.FB) <- "stim.StroCellType"
DefaultAssay(TwoMvWT.combined.FB) <- "RNA"
TwoMvWT.combined.FB <- ScaleData(TwoMvWT.combined.FB, features = rownames(TwoMvWT.combined.FB))
TwoMvWT.combined.FB.markers.1 <- FindAllMarkers(TwoMvWT.combined.FB, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
TwoMvWT.combined.FBTop50.1 <- TwoMvWT.combined.FB.markers.1 %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "TwoMvWT.combined.FB stim individual Heatmap Top50 purple-1.tiff", width = 20, height = 50, units = "in", compression = "lzw", res = 200)
DoHeatmap(TwoMvWT.combined.FB, features = c("ARQ", "EGFP", TwoMvWT.combined.FBTop50.1$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()


