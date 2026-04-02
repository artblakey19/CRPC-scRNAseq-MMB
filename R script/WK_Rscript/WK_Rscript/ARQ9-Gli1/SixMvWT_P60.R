####WT_P60####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/6M/SixMvWT_P60")

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
SixM[["orig.clusters"]] <- Idents(object = SixM)
WT_P60[["orig.clusters"]] <- Idents(object = WT_P60)

#Set Current idents
Idents(object = SixM) <- "seurat_clusters"
Idents(object = WT_P60) <- "seurat_clusters"

tiff(file = "SixM seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixM, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "WT_P60 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(WT_P60, reduction = "umap", pt.size = 0.3)
dev.off()

SixM$stim <- "ARQ9_6M"
WT_P60$stim <- "WT"

SixMvWT.anchors <- FindIntegrationAnchors(object.list = list(SixM, WT_P60), dims = 1:20)
SixMvWT.combined <- IntegrateData(anchorset = SixMvWT.anchors, dims = 1:20)

DefaultAssay(SixMvWT.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
SixMvWT.combined <- ScaleData(SixMvWT.combined, verbose = FALSE)
SixMvWT.combined <- RunPCA(SixMvWT.combined, npcs = 30, verbose = FALSE)
ElbowPlot(SixMvWT.combined, ndims = 50)

#Umap and Clustering
SixMvWT.combined <- FindNeighbors(SixMvWT.combined, reduction = "pca", dims = 1:20)
SixMvWT.combined <- FindClusters(SixMvWT.combined, resolution = 0.5)
SixMvWT.combined <- RunTSNE(SixMvWT.combined, reduction = "pca", dims = 1:20)
SixMvWT.combined <- RunUMAP(SixMvWT.combined, reduction = "pca", dims = 1:20)
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "SixMvWT.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
DefaultAssay(SixMvWT.combined) <- "RNA"
all.genes <- rownames(SixMvWT.combined)
SixMvWT.combined <- ScaleData(SixMvWT.combined, features = all.genes)
SixMvWT.combined <- CellCycleScoring(SixMvWT.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = SixMvWT.combined) <- "Phase"
DimPlot(SixMvWT.combined, reduction = "umap")
tiff(file = "SixMvWT.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
SixMvWT.combined1 <- SixMvWT.combined
DefaultAssay(SixMvWT.combined1) <- "integrated"
SixMvWT.combined1 <- ScaleData(SixMvWT.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SixMvWT.combined1))
SixMvWT.combined1 <- RunPCA(SixMvWT.combined1, features = VariableFeatures(SixMvWT.combined1))
ElbowPlot(SixMvWT.combined1, ndims = 50)

SixMvWT.combined1 <- FindNeighbors(SixMvWT.combined1, reduction = "pca", dims = 1:20)
SixMvWT.combined1 <- FindClusters(SixMvWT.combined1, resolution = 0.5)
SixMvWT.combined1 <- RunUMAP(SixMvWT.combined1, reduction = "pca", dims = 1:20)
SixMvWT.combined1 <- RunTSNE(SixMvWT.combined1, reduction = "pca", dims = 1:20)
DimPlot(SixMvWT.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(SixMvWT.combined1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")
#No need to do cell cycle regression...

Idents(object = SixMvWT.combined1) <- "Phase"
tiff(file = "SixMvWT.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
Idents(object = SixMvWT.combined1) <- "seurat_clusters"
tiff(file = "SixMvWT.combined1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = SixMvWT.combined) <- "seurat_clusters"
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)

#Cell type identification
DefaultAssay(SixMvWT.combined) <- "RNA"
tiff(file = "SixMvWT.combined celltype marker expression plots.tiff", width = 20, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixMvWT.combined, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2",
                                                                "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                                "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

FeaturePlot(SixMvWT.combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(SixMvWT.combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Rename
Idents(object = SixMvWT.combined) <- "seurat_clusters"
SixMvWT.combined <- RenameIdents(object = SixMvWT.combined, '1' = "BE", '2' = "BE", '19' = "BE", '3' = "LE", '7' = "LE", 
                                      '6' = "LE", '10' = "LE", '13' = "LE", '5' = "LE", '22' = "LE", '11' = "SV", '0' = "FB", '4' = "FB",
                                      '18' = "FB", '15' = "FB", '21' = "SM", '9' = "VE", '20' = "VE", '16' = "Glia", '12' = "Pericyte", '8' = "Leu", '17' = "Leu", '14' = "Leu")
SixMvWT.combined[["CellType"]] <- Idents(object = SixMvWT.combined)

Idents(object = SixMvWT.combined) <- "CellType"
tiff(file = "SixMvWT.combined CellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3)
dev.off()
tiff(file = "SixMvWT.combined CellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()
tiff(file = "SixMvWT.combined CellType LABEL UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()

#Cell counts
Idents(object = SixMvWT.combined) <- "stim"
SixMvWT.combined$stim.CellType <- paste(Idents(SixMvWT.combined), SixMvWT.combined$CellType, sep = "_")
Idents(object = SixMvWT.combined) <- "stim.CellType"
table(Idents(SixMvWT.combined))


####Re-clustering Epi####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Gli1/6M/SixMvWT_P60/Epi")

Idents(object = SixMvWT.combined) <- "CellType"
SixMvWT.combined.Epi <- subset(SixMvWT.combined, idents = c("BE", "LE", "SV"))
Idents(object = SixMvWT.combined.Epi) <- "seurat_clusters"
DefaultAssay(SixMvWT.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
SixMvWT.combined.Epi <- ScaleData(SixMvWT.combined.Epi, verbose = FALSE)
SixMvWT.combined.Epi <- RunPCA(SixMvWT.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(SixMvWT.combined.Epi, ndims = 50)

#Umap and Clustering
SixMvWT.combined.Epi <- FindNeighbors(SixMvWT.combined.Epi, reduction = "pca", dims = 1:20)
SixMvWT.combined.Epi <- FindClusters(SixMvWT.combined.Epi, resolution = 0.5)
SixMvWT.combined.Epi <- RunTSNE(SixMvWT.combined.Epi, reduction = "pca", dims = 1:20)
SixMvWT.combined.Epi <- RunUMAP(SixMvWT.combined.Epi, reduction = "pca", dims = 1:20)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, split.by = "stim")

tiff(file = "SixMvWT.combined.Epi UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "SixMvWT.combined.Epi Label UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()

#Cell cycle regression
DefaultAssay(SixMvWT.combined.Epi) <- "RNA"
all.genes <- rownames(SixMvWT.combined.Epi)
SixMvWT.combined.Epi <- ScaleData(SixMvWT.combined.Epi, features = all.genes)
SixMvWT.combined.Epi <- CellCycleScoring(SixMvWT.combined.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = SixMvWT.combined.Epi) <- "Phase"
DimPlot(SixMvWT.combined.Epi, reduction = "umap")
tiff(file = "SixMvWT.combined.Epi Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
SixMvWT.combined.Epi1 <- SixMvWT.combined.Epi
DefaultAssay(SixMvWT.combined.Epi1) <- "integrated"
SixMvWT.combined.Epi1 <- ScaleData(SixMvWT.combined.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(SixMvWT.combined.Epi1))
SixMvWT.combined.Epi1 <- RunPCA(SixMvWT.combined.Epi1, features = VariableFeatures(SixMvWT.combined.Epi1))
ElbowPlot(SixMvWT.combined.Epi1, ndims = 50)

SixMvWT.combined.Epi1 <- FindNeighbors(SixMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
SixMvWT.combined.Epi1 <- FindClusters(SixMvWT.combined.Epi1, resolution = 0.5)
SixMvWT.combined.Epi1 <- RunUMAP(SixMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
SixMvWT.combined.Epi1 <- RunTSNE(SixMvWT.combined.Epi1, reduction = "pca", dims = 1:20)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE, split.by = "stim")

Idents(object = SixMvWT.combined.Epi1) <- "Phase"
tiff(file = "SixMvWT.combined.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = SixMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "SixMvWT.combined.Epi1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = SixMvWT.combined.Epi1) <- "seurat_clusters"
tiff(file = "SixMvWT.combined.Epi1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(SixMvWT.combined.Epi) <- "RNA"
Idents(object = SixMvWT.combined.Epi) <- "seurat_clusters"
SixMvWT.combined.Epi <- ScaleData(SixMvWT.combined.Epi, features = rownames(SixMvWT.combined.Epi1))
SixMvWT.combined.Epi.markers <- FindAllMarkers(SixMvWT.combined.Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(SixMvWT.combined.Epi.markers, "SixMvWT.combined.Epi.seurat.markers.csv")

DefaultAssay(SixMvWT.combined.Epi) <- "RNA"
tiff(file = "SixMvWT.combined.Epi Epicelltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixMvWT.combined.Epi, reduction = "umap", features = c("Krt5", "Krt19", "Fbln1", "Tagln",
                                                                    "Gsdma", "Msmb", "Lrrc26", "Krt4",
                                                                    "Clu", "Svs5", "Ly6a", "Pbsn"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

DefaultAssay(SixMvWT.combined.Epi) <- "RNA"
tiff(file = "SixMvWT.combined.Epi Epicelltype marker expression plots-1.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(SixMvWT.combined.Epi, reduction = "umap", features = c("Krt5", "Krt14", "Trp63", "Krt4", "Ppp1r1b", "Krt7", "Krt19", "Tgm4", "Msmb", "Svs2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Rename
Idents(object = SixMvWT.combined.Epi) <- "seurat_clusters"
SixMvWT.combined.Epi <- RenameIdents(object = SixMvWT.combined.Epi, '0' = "BE1", '1' = "BE2", '10' = "BE3", '8' = "LE1", '2' = "LE2", '4' = "LE3", '9' = "LE4", '3' = "LE5", 
                                      '5' = "LE6", '7' = "SV", '6' = "OE", '11' = "OE", '13' = "OE", '12' = "OE", '14' = "OE")
SixMvWT.combined.Epi[["EpiCellType"]] <- Idents(object = SixMvWT.combined.Epi)

#Umap
Idents(object = SixMvWT.combined.Epi) <- "EpiCellType"
tiff(file = "SixMvWT.combined.Epi EpiCellType Label UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()
tiff(file = "SixMvWT.combined.Epi EpiCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "SixMvWT.combined.Epi EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = SixMvWT.combined.Epi) <- "stim"
SixMvWT.combined.Epi$stim.EpiCellType <- paste(Idents(SixMvWT.combined.Epi), SixMvWT.combined.Epi$EpiCellType, sep = "_")
Idents(object = SixMvWT.combined.Epi) <- "stim.EpiCellType"
table(Idents(SixMvWT.combined.Epi))

#W/O LP & VP
Idents(object = SixMvWT.combined.Epi) <- "EpiCellType"
SixMvWT.combined.Epi1 <- subset(SixMvWT.combined.Epi, idents = c("BE1", "BE2", "BE3", "LE1", "LE2", "LE3", "LE6", "SV", "OE"))

Idents(object = SixMvWT.combined.Epi1) <- "EpiCellType"
SixMvWT.combined.Epi1 <- RenameIdents(object = SixMvWT.combined.Epi1, 'BE1' = "BE1", 'BE2' = "BE2", 'BE3' = "BE3", 'LE1' = "LE1", 'LE2' = "LE2", 'LE3' = "LE3", 
                                     'LE6' = "LE4", 'SV' = "SV", 'OE' = "OE")
SixMvWT.combined.Epi1[["EpiCellType1"]] <- Idents(object = SixMvWT.combined.Epi1)

#Umap
Idents(object = SixMvWT.combined.Epi1) <- "EpiCellType1"
tiff(file = "SixMvWT.combined.Epi1 EpiCellType Label UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE, label.size = 6)
dev.off()
tiff(file = "SixMvWT.combined.Epi1 EpiCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

tiff(file = "SixMvWT.combined.Epi1 EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(SixMvWT.combined.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = SixMvWT.combined.Epi1) <- "stim"
SixMvWT.combined.Epi1$stim.EpiCellType1 <- paste(Idents(SixMvWT.combined.Epi1), SixMvWT.combined.Epi1$EpiCellType1, sep = "_")
Idents(object = SixMvWT.combined.Epi1) <- "stim.EpiCellType1"
table(Idents(SixMvWT.combined.Epi1))

#Heatmap_AllBE_SixMvWT
Idents(object = SixMvWT.combined.Epi) <- "EpiCellType"
SixMvWT.combined.BE <- subset(SixMvWT.combined.Epi, idents = c("BE1", "BE2", "BE3"))

DefaultAssay(SixMvWT.combined.BE) <- "RNA"
Idents(object = SixMvWT.combined.BE) <- "stim"
SixMvWT.combined.BE <- ScaleData(SixMvWT.combined.BE, features = rownames(SixMvWT.combined.BE))
SixMvWT.combined.BE.markers <- FindAllMarkers(SixMvWT.combined.BE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SixMvWT.combined.BETop50 <- SixMvWT.combined.BE.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "SixMvWT.combined.BE Heatmap Top50 purple.tiff", width = 20, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(SixMvWT.combined.BE, features = c(SixMvWT.combined.BETop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 15))
dev.off()

#DEGs
DefaultAssay(SixMvWT.combined.BE) <- "RNA"
Idents(object = SixMvWT.combined.BE) <- "stim"
SixMvWT.combined.BE <- ScaleData(SixMvWT.combined.BE, features = rownames(SixMvWT.combined.BE))
SixMvWT.combined.BE.0.Markers <- FindMarkers(SixMvWT.combined.BE, ident.1 = "ARQ9_6M", ident.2 = "WT", min.pct = 0, logfc.threshold = 0)
write.csv(SixMvWT.combined.BE.0.Markers, "SixMvWT.combined.BE.0.Markers.csv")

#p.adjust
SixMvWT.combined.BE <- read.csv("SixMvWT.combined.BE.0.Markers.csv") 
SixMvWT.combined.BE_pvalue <- SixMvWT.combined.BE$p_val
SixMvWT.combined.BE_pvalue=as.numeric(SixMvWT.combined.BE_pvalue)
SixMvWT.combined.BE_BH = p.adjust(SixMvWT.combined.BE_pvalue, "BH")
write.csv(SixMvWT.combined.BE_BH, "SixMvWT.combined.BE_BH.csv")

#Heatmap_BE1_SixMvWT
Idents(object = SixMvWT.combined.Epi) <- "EpiCellType"
SixMvWT.combined.BE1 <- subset(SixMvWT.combined.Epi, idents = c("BE1"))

DefaultAssay(SixMvWT.combined.BE1) <- "RNA"
Idents(object = SixMvWT.combined.BE1) <- "stim"
SixMvWT.combined.BE1 <- ScaleData(SixMvWT.combined.BE1, features = rownames(SixMvWT.combined.BE1))
SixMvWT.combined.BE1.markers <- FindAllMarkers(SixMvWT.combined.BE1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SixMvWT.combined.BE1Top50 <- SixMvWT.combined.BE1.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "SixMvWT.combined.BE1 Heatmap Top50 purple.tiff", width = 20, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(SixMvWT.combined.BE1, features = c(SixMvWT.combined.BE1Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 15))
dev.off()

#DEGs
DefaultAssay(SixMvWT.combined.BE1) <- "RNA"
Idents(object = SixMvWT.combined.BE1) <- "stim"
SixMvWT.combined.BE1 <- ScaleData(SixMvWT.combined.BE1, features = rownames(SixMvWT.combined.BE1))
SixMvWT.combined.BE1.0.Markers <- FindMarkers(SixMvWT.combined.BE1, ident.1 = "ARQ9_6M", ident.2 = "WT", min.pct = 0, logfc.threshold = 0)
write.csv(SixMvWT.combined.BE1.0.Markers, "SixMvWT.combined.BE1.0.Markers.csv")

#p.adjust
SixMvWT.combined.BE1 <- read.csv("SixMvWT.combined.BE1.0.Markers.csv") 
SixMvWT.combined.BE1_pvalue <- SixMvWT.combined.BE1$p_val
SixMvWT.combined.BE1_pvalue=as.numeric(SixMvWT.combined.BE1_pvalue)
SixMvWT.combined.BE1_BH = p.adjust(SixMvWT.combined.BE1_pvalue, "BH")
write.csv(SixMvWT.combined.BE1_BH, "SixMvWT.combined.BE1_BH.csv")

#Heatmap_BE2_SixMvWT
Idents(object = SixMvWT.combined.Epi) <- "EpiCellType"
SixMvWT.combined.BE2 <- subset(SixMvWT.combined.Epi, idents = c("BE2"))

DefaultAssay(SixMvWT.combined.BE2) <- "RNA"
Idents(object = SixMvWT.combined.BE2) <- "stim"
SixMvWT.combined.BE2 <- ScaleData(SixMvWT.combined.BE2, features = rownames(SixMvWT.combined.BE2))
SixMvWT.combined.BE2.markers <- FindAllMarkers(SixMvWT.combined.BE2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SixMvWT.combined.BE2Top50 <- SixMvWT.combined.BE2.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "SixMvWT.combined.BE2 Heatmap Top50 purple.tiff", width = 20, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(SixMvWT.combined.BE2, features = c(SixMvWT.combined.BE2Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 15))
dev.off()

#DEGs
DefaultAssay(SixMvWT.combined.BE2) <- "RNA"
Idents(object = SixMvWT.combined.BE2) <- "stim"
SixMvWT.combined.BE2 <- ScaleData(SixMvWT.combined.BE2, features = rownames(SixMvWT.combined.BE2))
SixMvWT.combined.BE2.0.Markers <- FindMarkers(SixMvWT.combined.BE2, ident.1 = "ARQ9_6M", ident.2 = "WT", min.pct = 0, logfc.threshold = 0)
write.csv(SixMvWT.combined.BE2.0.Markers, "SixMvWT.combined.BE2.0.Markers.csv")

#p.adjust
SixMvWT.combined.BE2 <- read.csv("SixMvWT.combined.BE2.0.Markers.csv") 
SixMvWT.combined.BE2_pvalue <- SixMvWT.combined.BE2$p_val
SixMvWT.combined.BE2_pvalue=as.numeric(SixMvWT.combined.BE2_pvalue)
SixMvWT.combined.BE2_BH = p.adjust(SixMvWT.combined.BE2_pvalue, "BH")
write.csv(SixMvWT.combined.BE2_BH, "SixMvWT.combined.BE2_BH.csv")

#Heatmap_BE3_SixMvWT
Idents(object = SixMvWT.combined.Epi) <- "EpiCellType"
SixMvWT.combined.BE3 <- subset(SixMvWT.combined.Epi, idents = c("BE3"))

DefaultAssay(SixMvWT.combined.BE3) <- "RNA"
Idents(object = SixMvWT.combined.BE3) <- "stim"
SixMvWT.combined.BE3 <- ScaleData(SixMvWT.combined.BE3, features = rownames(SixMvWT.combined.BE3))
SixMvWT.combined.BE3.markers <- FindAllMarkers(SixMvWT.combined.BE3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SixMvWT.combined.BE3Top50 <- SixMvWT.combined.BE3.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "SixMvWT.combined.BE3 Heatmap Top50 purple.tiff", width = 20, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(SixMvWT.combined.BE3, features = c(SixMvWT.combined.BE3Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 15))
dev.off()

#DEGs
DefaultAssay(SixMvWT.combined.BE3) <- "RNA"
Idents(object = SixMvWT.combined.BE3) <- "stim"
SixMvWT.combined.BE3 <- ScaleData(SixMvWT.combined.BE3, features = rownames(SixMvWT.combined.BE3))
SixMvWT.combined.BE3.0.Markers <- FindMarkers(SixMvWT.combined.BE3, ident.1 = "ARQ9_6M", ident.2 = "WT", min.pct = 0, logfc.threshold = 0)
write.csv(SixMvWT.combined.BE3.0.Markers, "SixMvWT.combined.BE3.0.Markers.csv")

#p.adjust
SixMvWT.combined.BE3 <- read.csv("SixMvWT.combined.BE3.0.Markers.csv") 
SixMvWT.combined.BE3_pvalue <- SixMvWT.combined.BE3$p_val
SixMvWT.combined.BE3_pvalue=as.numeric(SixMvWT.combined.BE3_pvalue)
SixMvWT.combined.BE3_BH = p.adjust(SixMvWT.combined.BE3_pvalue, "BH")
write.csv(SixMvWT.combined.BE3_BH, "SixMvWT.combined.BE3_BH.csv")

#Heatmap_SixM_BE1vBE2VBE3
Idents(object = SixMvWT.combined.Epi) <- "stim.EpiCellType"
SixM.BE <- subset(SixMvWT.combined.Epi, idents = c("ARQ9_6M_BE1", "ARQ9_6M_BE2", "ARQ9_6M_BE3"))

DefaultAssay(SixM.BE) <- "RNA"
Idents(object = SixM.BE) <- "EpiCellType"
SixM.BE <- ScaleData(SixM.BE, features = rownames(SixM.BE))
SixM.BE.markers <- FindAllMarkers(SixM.BE, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
SixM.BETop50 <- SixM.BE.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_log2FC)
tiff(file = "SixM.BE Heatmap Top50 purple.tiff", width = 20, height = 30, units = "in", compression = "lzw", res = 200)
DoHeatmap(SixM.BE, features = c(SixM.BETop50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow")) + theme(axis.text.y = element_text(size = 15))
dev.off()

#DEGs
DefaultAssay(SixM.BE) <- "RNA"
Idents(object = SixM.BE) <- "EpiCellType"
SixM.BE <- ScaleData(SixM.BE, features = rownames(SixM.BE))
SixM.BE3vBE1_2.0.Markers <- FindMarkers(SixM.BE, ident.1 = "BE3", ident.2 = c("BE1", "BE2"), min.pct = 0, logfc.threshold = 0)
write.csv(SixM.BE3vBE1_2.0.Markers, "SixM.BE3vBE1_2.0.Markers.csv")

#p.adjust
SixM.BE3vBE1_2 <- read.csv("SixM.BE3vBE1_2.0.Markers.csv") 
SixM.BE3vBE1_2_pvalue <- SixM.BE3vBE1_2$p_val
SixM.BE3vBE1_2_pvalue=as.numeric(SixM.BE3vBE1_2_pvalue)
SixM.BE3vBE1_2_BH = p.adjust(SixM.BE3vBE1_2_pvalue, "BH")
write.csv(SixM.BE3vBE1_2_BH, "SixM.BE3vBE1_2_BH.csv")