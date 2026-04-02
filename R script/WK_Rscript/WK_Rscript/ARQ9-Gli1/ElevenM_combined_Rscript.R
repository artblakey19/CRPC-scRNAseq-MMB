#### Merging Datasets TwoMvWT ####

#Stash old idents
ElevenM1[["orig.clusters"]] <- Idents(object = ElevenM1)
ElevenM2[["orig.clusters"]] <- Idents(object = ElevenM2)

#Set Current idents
Idents(object = ElevenM1) <- "seurat_clusters"
Idents(object = ElevenM2) <- "seurat_clusters"
ElevenM1$stim <- "ARQ9_11M_A0843"
ElevenM2$stim <- "ARQ9_11M_A1350"
ElevenM.anchors <- FindIntegrationAnchors(object.list = list(ElevenM1, ElevenM2), dims = 1:20)
ElevenM.combined <- IntegrateData(anchorset = ElevenM.anchors, dims = 1:20)

DefaultAssay(ElevenM.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
ElevenM.combined <- ScaleData(ElevenM.combined, verbose = FALSE)
ElevenM.combined <- RunPCA(ElevenM.combined, npcs = 30, verbose = FALSE)

#Umap and Clustering
ElevenM.combined <- FindNeighbors(ElevenM.combined, reduction = "pca", dims = 1:20)
ElevenM.combined <- FindClusters(ElevenM.combined, resolution = 0.5)
ElevenM.combined <- RunTSNE(ElevenM.combined, reduction = "pca", dims = 1:20)
ElevenM.combined <- RunUMAP(ElevenM.combined, reduction = "pca", dims = 1:20)

tiff(file = "ElevenM.combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ElevenM.combined, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell cycle regression
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(ElevenM.combined) <- "RNA"
all.genes <- rownames(ElevenM.combined)
ElevenM.combined <- ScaleData(ElevenM.combined, features = all.genes)
ElevenM.combined <- CellCycleScoring(ElevenM.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = ElevenM.combined) <- "Phase"
DimPlot(ElevenM.combined, reduction = "umap")
tiff(file = "ElevenM.combined Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ElevenM.combined, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
ElevenM.combined1 <- ElevenM.combined
DefaultAssay(ElevenM.combined1) <- "integrated"
ElevenM.combined1 <- ScaleData(ElevenM.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(ElevenM.combined1))
ElevenM.combined1 <- RunPCA(ElevenM.combined1, features = VariableFeatures(ElevenM.combined1))
ElbowPlot(ElevenM.combined1, ndims = 50)

ElevenM.combined1 <- FindNeighbors(ElevenM.combined1, reduction = "pca", dims = 1:25)
ElevenM.combined1 <- FindClusters(ElevenM.combined1, resolution = 0.5)
ElevenM.combined1 <- RunUMAP(ElevenM.combined1, reduction = "pca", dims = 1:25)
ElevenM.combined1 <- RunTSNE(ElevenM.combined1, reduction = "pca", dims = 1:25)
DimPlot(ElevenM.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = ElevenM.combined1) <- "Phase"
tiff(file = "ElevenM.combined1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ElevenM.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = ElevenM.combined1) <- "seurat_clusters"
tiff(file = "ElevenM.combined1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ElevenM.combined1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = ElevenM.combined1) <- "seurat_clusters"
tiff(file = "ElevenM.combined1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ElevenM.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell type identification
DefaultAssay(ElevenM.combined1) <- "RNA"
Idents(object = ElevenM.combined1) <- "seurat_clusters"
ElevenM.combined1 <- ScaleData(ElevenM.combined1, features = rownames(ElevenM.combined1))
ElevenM.combined1.markers <- FindAllMarkers(ElevenM.combined1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(ElevenM.combined1.markers, "ElevenM.combined1.seurat.markers.csv")

DefaultAssay(ElevenM.combined1) <- "RNA"
tiff(file = "ElevenM.combined1 celltype marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(ElevenM.combined1, reduction = "umap", features = c("Krt5", "Trp63", "Krt19", "Pbsn", "Svs2",
                                                                "Fbln1", "Myh11",  "Pecam1", "Plp1",
                                                                "Rgs5", "Tyrobp", "Ccl5"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()






#Rename
Idents(object = ElevenM.combined1) <- "seurat_clusters"
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
