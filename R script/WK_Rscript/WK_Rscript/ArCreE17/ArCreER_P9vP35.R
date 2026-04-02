###ARCreER P9 P35 Workflow###

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
library(ggpubr)
library(GGally)

setwd("//isi-dcnl/user_data/zjsun/group/Aly Buckley/ARCreER/ARCreER 20220120")

#### Merging Datasets ####

#Stash old idents
P9[["orig.clusters"]] <- Idents(object = P9)
P35[["orig.clusters"]] <- Idents(object = P35)

#Set Current idents
Idents(object = P9) <- "seurat_clusters"
Idents(object = P35) <- "seurat_clusters"
P9$stim <- "P9"
P35$stim <- "P35"
P9vP35.anchors <- FindIntegrationAnchors(object.list = list(P9, P35), dims = 1:20)
P9vP35.combined <- IntegrateData(anchorset = P9vP35.anchors, dims = 1:20)
DefaultAssay(P9vP35.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
P9vP35.combined <- ScaleData(P9vP35.combined, verbose = FALSE)
P9vP35.combined <- RunPCA(P9vP35.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
P9vP35.combined <- FindNeighbors(P9vP35.combined, reduction = "pca", dims = 1:20)
P9vP35.combined <- FindClusters(P9vP35.combined, resolution = 0.5)
P9vP35.combined <- RunTSNE(P9vP35.combined, reduction = "pca", dims = 1:20)
P9vP35.combined <- RunUMAP(P9vP35.combined, reduction = "pca", dims = 1:20)

####Cell Cycle Regression####
mouse_cell_cycle_genes <- readRDS("mouse_cell_cycle_genes.Rds")
s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(P9vP35.combined) <- "RNA"
all.genes <- rownames(P9vP35.combined)
P9vP35.combined <- ScaleData(P9vP35.combined, features = all.genes)
P9vP35.combined <- CellCycleScoring(P9vP35.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P9vP35.combined) <- "Phase"
DimPlot(P9vP35.combined, reduction = "umap")
tiff(file = "P9vP35.combined Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()
DimPlot(P9vP35.combined1)
#Take Cell cycle out 
P9vP35.combined1 <- P9vP35.combined
DefaultAssay(P9vP35.combined1) <- "integrated"
P9vP35.combined1 <- ScaleData(P9vP35.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P9vP35.combined1))
P9vP35.combined1 <- RunPCA(P9vP35.combined1, features = VariableFeatures(P9vP35.combined1))
ElbowPlot(P9vP35.combined1, ndims = 30)
P9vP35.combined1 <- FindNeighbors(P9vP35.combined1, reduction = "pca", dims = 1:18)
P9vP35.combined1 <- FindClusters(P9vP35.combined1, resolution = 0.5)
P9vP35.combined1 <- RunUMAP(P9vP35.combined1, reduction = "pca", dims = 1:18)
P9vP35.combined1 <- RunTSNE(P9vP35.combined1, reduction = "pca", dims = 1:18)
DimPlot(P9vP35.combined1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P9vP35.combined1) <- "Phase"
tiff(file = "P9vP35.combined1 Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = P9vP35.combined1) <- "stim"
P9vP35.combined1 <- RenameIdents(object = P9vP35.combined1, 'P9' = "P9", 'P35' = "P35")

tiff(file = "P9vP35.combined1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1, reduction = "umap", pt.size = 0.3, cols = c("gray", "#B71800")) 
dev.off()

P9.combined1 <- subset(P9vP35.combined1, idents = c("P9"))
Idents(object = P9.combined1) <- "stim"
tiff(file = "P9.combined1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9.combined1, reduction = "umap", pt.size = 0.3, cols = c("gray"))
dev.off()

P35.combined1 <- subset(P9vP35.combined1, idents = c("P35"))
Idents(object = P35.combined1) <- "stim"
tiff(file = "P35.combined1 UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P35.combined1, reduction = "umap", pt.size = 0.3, cols = c("#B71800"))
dev.off()

table(Idents(P9vP35.combined1))

#Featureplots
DefaultAssay(P9vP35.combined1) <- "RNA"
tiff(file = "P9vP35.combined1 expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1, reduction = "umap", features = c("Fbln1", "Myh11", "Krt5", "Krt19", "Pecam1", 
                                                             "Mpz", "Cd52", "Hoxb9", "Rgs5", "C1qa",
                                                             "Gcg", "Chgb", "Upk3a", "Ccl5", "Mki67"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "P9vP35.combined1 Ar expression plots split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", split.by = "stim")
dev.off()
tiff(file = "P9vP35.combined1 EGFP expression plots split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90", split.by = "stim")
dev.off()

#New labeling
Idents(object = P9vP35.combined1) <- "seurat_clusters"
P9vP35.combined1 <- RenameIdents(object = P9vP35.combined1, '2' = "BE", '8' = "BE", '11' = "Ductus Deferen", '7' = "LE", '5' = "LE", '12' = "LE", '1' = "FB", '3' = "FB", '4' = "FB", '0' = "SM", '13' = "SM", '16' = "SM", '6' = "Glia", '9' = "Endo", '10' = "Peri", '14' = "Leu", '17' = "Leu", '18' = "Leu", '15' = "Bladder")
P9vP35.combined1[["CellType"]] <- Idents(object = P9vP35.combined1)

tiff(file = "P9vP35.combined1 CellType UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()

tiff(file = "P9vP35.combined1 CellType UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

#Cell counts
Idents(object = P9vP35.combined1) <- "stim"
P9vP35.combined1$stim.CellType <- paste(Idents(P9vP35.combined1), P9vP35.combined1$CellType, sep = "_")
Idents(object = P9vP35.combined1) <- "stim.CellType"
table(Idents(P9vP35.combined1))

#Add GFP info
DefaultAssay(P9vP35.combined1) <- "RNA"
P9vP35.combined1EGFPPos <- subset(x=P9vP35.combined1, subset = EGFP > 0)
P9vP35.combined1EGFPNeg <- subset(x=P9vP35.combined1, subset = EGFP == 0)
Idents(object = P9vP35.combined1EGFPPos) <- "EGFPPos"
Idents(object = P9vP35.combined1EGFPNeg) <- "EGFPNeg"
P9vP35.combined1EGFPPos[["EGFPExp"]] <- Idents(object = P9vP35.combined1EGFPPos)
P9vP35.combined1EGFPNeg[["EGFPExp"]] <- Idents(object = P9vP35.combined1EGFPNeg)
P9vP35.combined1EGFP <- merge(x = P9vP35.combined1EGFPPos, y = P9vP35.combined1EGFPNeg)
Idents(object = P9vP35.combined1EGFP) <- "EGFPExp"
P9vP35.combined1$EGFPExp <- Idents(object = P9vP35.combined1EGFP)
Idents(object = P9vP35.combined1) <- "EGFPExp"
P9vP35.combined1$EGFPExp.CellType <- paste(Idents(P9vP35.combined1), P9vP35.combined1$CellType, sep = "_")

#Cell counts
Idents(object = P9vP35.combined1) <- "EGFPExp.CellType"
P9vP35.combined1$stim.EGFPExp.CellType <- paste(Idents(P9vP35.combined1), P9vP35.combined1$stim, sep = "_")
Idents(object = P9vP35.combined1) <- "stim.EGFPExp.CellType"
table(Idents(P9vP35.combined1))

#Stem cell marker
DefaultAssay(P9vP35.combined1) <- "RNA"
tiff(file = "P9vP35.combined1 stem cell marker expression plots.tiff", width = 15, height = 15, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1, reduction = "umap", features = c("EGFP", "Ar", "Zeb1", "Prom1", "Itga6", "Ly6a", "Tacstd2", "Trp63", "Psca"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Reclustering Epithelial cells ####
Idents(object = P9vP35.combined1) <- "CellType"
P9vP35.combined1.Epi <- subset(P9vP35.combined1, idents = c("BE", "LE"))
Idents(object = P9vP35.combined1.Epi) <- "seurat_clusters"
DefaultAssay(P9vP35.combined1.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
P9vP35.combined1.Epi <- ScaleData(P9vP35.combined1.Epi, verbose = FALSE)
P9vP35.combined1.Epi <- RunPCA(P9vP35.combined1.Epi, npcs = 30, verbose = FALSE)

#Umap and Clustering
ElbowPlot(P9vP35.combined1.Epi, ndims = 30)
P9vP35.combined1.Epi <- FindNeighbors(P9vP35.combined1.Epi, reduction = "pca", dims = 1:18)
P9vP35.combined1.Epi <- FindClusters(P9vP35.combined1.Epi, resolution = 0.5)
P9vP35.combined1.Epi <- RunTSNE(P9vP35.combined1.Epi, reduction = "pca", dims = 1:18)
P9vP35.combined1.Epi <- RunUMAP(P9vP35.combined1.Epi, reduction = "pca", dims = 1:18)

tiff(file = "P9vP35.combined1.Epi UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell Cycle
DefaultAssay(P9vP35.combined1.Epi) <- "RNA"
all.genes <- rownames(P9vP35.combined1.Epi)
P9vP35.combined1.Epi <- ScaleData(P9vP35.combined1.Epi, features = all.genes)
P9vP35.combined1.Epi <- CellCycleScoring(P9vP35.combined1.Epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P9vP35.combined1.Epi) <- "Phase"
DimPlot(P9vP35.combined1.Epi, reduction = "umap")
tiff(file = "P9vP35.combined1.Epi Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
P9vP35.combined1.Epi1 <- P9vP35.combined1.Epi
DefaultAssay(P9vP35.combined1.Epi1) <- "integrated"
P9vP35.combined1.Epi1 <- ScaleData(P9vP35.combined1.Epi1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P9vP35.combined1.Epi1))
P9vP35.combined1.Epi1 <- RunPCA(P9vP35.combined1.Epi1, features = VariableFeatures(P9vP35.combined1.Epi1))
ElbowPlot(P9vP35.combined1.Epi1, ndims = 30)
P9vP35.combined1.Epi1 <- FindNeighbors(P9vP35.combined1.Epi1, reduction = "pca", dims = 1:17)
P9vP35.combined1.Epi1 <- FindClusters(P9vP35.combined1.Epi1, resolution = 0.5)
P9vP35.combined1.Epi1 <- RunUMAP(P9vP35.combined1.Epi1, reduction = "pca", dims = 1:17)
P9vP35.combined1.Epi1 <- RunTSNE(P9vP35.combined1.Epi1, reduction = "pca", dims = 1:17)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P9vP35.combined1.Epi1) <- "Phase"
tiff(file = "P9vP35.combined1.Epi1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = P9vP35.combined1.Epi1) <- "stim"
P9vP35.combined1.Epi1 <- RenameIdents(object = P9vP35.combined1.Epi1, 'P9' = "P9", 'P35' = "P35")
P9vP35.combined1.Epi1[["stim"]] <- Idents(object = P9vP35.combined1.Epi1)

Idents(object = P9vP35.combined1.Epi1) <- "seurat_clusters"
tiff(file = "P9vP35.combined1.Epi1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", pt.size = 0.3, label= TRUE)
dev.off()

Idents(object = P9vP35.combined1.Epi1) <- "seurat_clusters"
tiff(file = "P9vP35.combined1.Epi1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

Idents(object = P9vP35.combined1.Epi1) <- "seurat_clusters"
P9vP35.combined1.Epi1 <- subset(P9vP35.combined1.Epi1, idents = c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))


#Cell counts
Idents(object = P9vP35.combined1.Epi1) <- "stim"
table(Idents(P9vP35.combined1.Epi1))

#FeaturePlots
DefaultAssay(P9vP35.combined1.Epi1) <- "RNA"
tiff(file = "P9vP35.combined1.Epi1 expression plots.tiff", width = 16, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1.Epi1, reduction = "umap", features = c("Ar", "EGFP", "Epcam", "Krt5", "Trp63", "Krt8", "Krt19", "Pbsn", "Vim", "Fbln1", "Acta2", "Cd34"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()
FeaturePlot(P9vP35.combined1.Epi1, reduction = "umap", features = c("Tgm4","Mmp7", "Cldn10","Msmb", "Trpv6", "Ppp1r1b"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
FeaturePlot(P9vP35.combined1.Epi1, reduction = "umap", features = c("Ube2c","Top2a","Birc5","Mki67","Cdk1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
tiff(file = "P9vP35.combined1.Epi1 expression plots2.tiff", width = 16, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1.Epi1, reduction = "umap", features = c("Krt5", "Krt4", "Aqp3", "Ppp1r1b","Tgm4", "Msmb", "Trpv6", "Trp63"))
dev.off()

FeaturePlot(P9vP35.combined1.Epi1, reduction = "umap", features = c("Wfdc2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P9vP35.combined1.Epi1, reduction = "umap", features = c("Atp10b"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Rename
Idents(object = P9vP35.combined1.Epi1) <- "seurat_clusters"
P9vP35.combined1.Epi1 <- RenameIdents(object = P9vP35.combined1.Epi1, '0' = "BE1", '4' = "BE2", '1' = "BE3", '3' = "Intermediate", '2' = "LE1", '8' = "LE2", '7' = "LE3", '9' = "UrLE", '6' = "FB-BE", '5' = "SMA-BE")
P9vP35.combined1.Epi1[["EpiCellType"]] <- Idents(object = P9vP35.combined1.Epi1)

#Umap
Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
tiff(file = "P9vP35.combined1.Epi1 EpiCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
tiff(file = "P9vP35.combined1.Epi1 EpiCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
P9vP35.combined1.Epi1$stim.EpiCellType <- paste(Idents(P9vP35.combined1.Epi1), P9vP35.combined1.Epi1$stim, sep = "_")
Idents(object = P9vP35.combined1.Epi1) <- "stim.EpiCellType"
table(Idents(P9vP35.combined1.Epi1))

Idents(object = P9vP35.combined1.Epi1) <- "EGFPExp"
P9vP35.combined1.Epi1$EGFPExp.stim.EpiCellType <- paste(Idents(P9vP35.combined1.Epi1), P9vP35.combined1.Epi1$stim.EpiCellType, sep = "_")
Idents(object = P9vP35.combined1.Epi1) <- "EGFPExp.stim.EpiCellType"
table(Idents(P9vP35.combined1.Epi1))

##ALY

#### Reclustering Stromal cells ####
Idents(object = P9vP35.combined1) <- "CellType"
P9vP35.combined1.stro <- subset(P9vP35.combined1, idents = c("FB", "SM"))
Idents(object = P9vP35.combined1.stro) <- "seurat_clusters"
DefaultAssay(P9vP35.combined1.stro) <- "integrated"

#Run the standard workflow for visualization and clustering
P9vP35.combined1.stro <- ScaleData(P9vP35.combined1.stro, verbose = FALSE)
P9vP35.combined1.stro <- RunPCA(P9vP35.combined1.stro, npcs = 30, verbose = FALSE)

#Umap and Clustering
ElbowPlot(P9vP35.combined1.stro, ndims = 30)
P9vP35.combined1.stro <- FindNeighbors(P9vP35.combined1.stro, reduction = "pca", dims = 1:18)
P9vP35.combined1.stro <- FindClusters(P9vP35.combined1.stro, resolution = 0.5)
P9vP35.combined1.stro <- RunTSNE(P9vP35.combined1.stro, reduction = "pca", dims = 1:18)
P9vP35.combined1.stro <- RunUMAP(P9vP35.combined1.stro, reduction = "pca", dims = 1:18)

tiff(file = "P9vP35.combined1.stro UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro, reduction = "umap", pt.size = 0.3)
dev.off()

#Cell Cycle
DefaultAssay(P9vP35.combined1.stro) <- "RNA"
all.genes <- rownames(P9vP35.combined1.stro)
P9vP35.combined1.stro <- ScaleData(P9vP35.combined1.stro, features = all.genes)
P9vP35.combined1.stro <- CellCycleScoring(P9vP35.combined1.stro, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = P9vP35.combined1.stro) <- "Phase"
DimPlot(P9vP35.combined1.stro, reduction = "umap")
tiff(file = "P9vP35.combined1.stro Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

#Take Cell cycle out 
P9vP35.combined1.stro1 <- P9vP35.combined1.stro
DefaultAssay(P9vP35.combined1.stro1) <- "integrated"
P9vP35.combined1.stro1 <- ScaleData(P9vP35.combined1.stro1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(P9vP35.combined1.stro1))
P9vP35.combined1.stro1 <- RunPCA(P9vP35.combined1.stro1, features = VariableFeatures(P9vP35.combined1.stro1))
ElbowPlot(P9vP35.combined1.stro1, ndims = 30)
P9vP35.combined1.stro1 <- FindNeighbors(P9vP35.combined1.stro1, reduction = "pca", dims = 1:17)
P9vP35.combined1.stro1 <- FindClusters(P9vP35.combined1.stro1, resolution = 0.5)
P9vP35.combined1.stro1 <- RunUMAP(P9vP35.combined1.stro1, reduction = "pca", dims = 1:17)
P9vP35.combined1.stro1 <- RunTSNE(P9vP35.combined1.stro1, reduction = "pca", dims = 1:17)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P9vP35.combined1.stro1) <- "seurat_clusters"
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, label = TRUE)
P9vP35.combined1.stro1 <- subset(P9vP35.combined1.stro1, idents = c("0", "1", "2","3","4","5","6","7","8","9"))
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, label = TRUE)


Idents(object = P9vP35.combined1.stro1) <- "Phase"
tiff(file = "P9vP35.combined1.stro1 Cell Cyle UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, cols = c("#1D762E", "#FFD966", "purple"))
dev.off()

Idents(object = P9vP35.combined1.stro1) <- "stim"
P9vP35.combined1.stro1 <- RenameIdents(object = P9vP35.combined1.stro1, 'P9' = "P9", 'P35' = "P35")
P9vP35.combined1.stro1[["stim"]] <- Idents(object = P9vP35.combined1.stro1)

Idents(object = P9vP35.combined1.stro1) <- "seurat_clusters"
tiff(file = "P9vP35.combined1.stro1 seurat UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, label= TRUE)
dev.off()

Idents(object = P9vP35.combined1.stro1) <- "seurat_clusters"
tiff(file = "P9vP35.combined1.stro1 seurat split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = P9vP35.combined1.stro1) <- "stim"
table(Idents(P9vP35.combined1.stro1))

#FeaturePlots
DefaultAssay(P9vP35.combined1.stro1) <- "RNA"
tiff(file = "P9vP35.combined1.stro1 expression plots.tiff", width = 16, height = 12, units = "in", compression = "lzw", res = 800)
FeaturePlot(P9vP35.combined1.stro1, reduction = "umap", features = c("Ar", "EGFP", "Lum","Apod", "Igfbp3", "Fbln1","Col15a1","Actg2", "Myh11", "Tagln", "Myl9", "Acta2"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = 3.0)
dev.off()


DotPlot(P9vP35.combined1.stro1, features = c("Lum", "Apod", "Igfbp3", "Fbln1", "Col15a1", "Actg2", "Myh11", "Acta2", "Myl9", "Tagln", "Ar", "EGFP"), cols = c("light grey", "red")) + RotatedAxis()

#Rename
Idents(object = P9vP35.combined1.stro1) <- "seurat_clusters"
P9vP35.combined1.stro1 <- RenameIdents(object = P9vP35.combined1.stro1, '0' = "SM1", '2'="SM2", '7'="SM3", '8'= "SM4", '1'= "FB1", '3'="FB2", '4'="FB3", '5'="FB4", '6'="FB5", '9'="FB6")
P9vP35.combined1.stro1[["StroCellType"]] <- Idents(object = P9vP35.combined1.stro1)

#Umap
Idents(object = P9vP35.combined1.stro1) <- "StroCellType"
tiff(file = "P9vP35.combined1.Stro1 StroCellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3)
dev.off()

Idents(object = P9vP35.combined1.stro1) <- "StroCellType"
tiff(file = "P9vP35.combined1.Stro1 STroCellType split UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell counts
Idents(object = P9vP35.combined1.stro1) <- "StroCellType"
P9vP35.combined1.stro1$stim.stroCellType <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$stim, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "stim.stroCellType"
table(Idents(P9vP35.combined1.stro1))

Idents(object = P9vP35.combined1.stro1) <- "EGFPExp"
P9vP35.combined1.stro1$EGFPExp.stim.StroCellType <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$stim.StroCellType, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "EGFPExp.stim.StroCellType"
table(Idents(P9vP35.combined1.stro1))

Idents(object = P9vP35.combined1.stro1) <- "ArExp"
P9vP35.combined1.stro1$ArExp.stim.StroCellType <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$stim.StroCellType, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "ArExp.stim.StroCellType"
table(Idents(P9vP35.combined1.stro1))


#### AR EGFP Exp Analysis #### P9vP35

DoublePos <- subset(x=P9vP35.combined1, subset = Ar > 0.01 & EGFP > 1)
Aronly <- subset(x=P9vP35.combined1, subset = Ar > 0.01 & EGFP < 1)
EGFPonly <- subset(x=P9vP35.combined1, subset = Ar < 0.01 & EGFP > 1)
DoubleNeg <- subset(x=P9vP35.combined1, subset = Ar < 0.01 & EGFP < 1)
Idents(object = DoublePos) <- "DoublePos"
Idents(object = Aronly) <- "Aronly"
Idents(object = EGFPonly) <- "EGFPonly"
Idents(object = DoubleNeg) <- "DoubleNeg"
DoublePos[["ArEGFPExp"]] <- Idents(object = DoublePos)
Aronly[["ArEGFPExp"]] <- Idents(object = Aronly)
EGFPonly[["ArEGFPExp"]] <- Idents(object = EGFPonly)
DoubleNeg[["ArEGFPExp"]] <- Idents(object = DoubleNeg)
Temp4 <- list(DoublePos, DoubleNeg, Aronly, EGFPonly)
ArEGFP1 <- merge(x = DoublePos, y = Aronly)
ArEGFP2 <- merge(x = EGFPonly, y = DoubleNeg)
ArEGFP <- merge(x = ArEGFP1, y = ArEGFP2)
P9vP35.combined1[["ArEGFPExp"]] <- Idents(object = ArEGFP)
tiff(file = "P9vP35.combined1 ArEGFPExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1, reduction = "umap", group.by = "ArEGFPExp", split.by='stim')
dev.off()

Idents(object = P9vP35.combined1) <- "CellClass"
P9vP35.combined1$CellClass.ArEGFPExp <- paste(Idents(P9vP35.combined1), P9vP35.combined1$ArEGFPExp, sep = "_")
Idents(object = P9vP35.combined1) <- "CellClass.ArEGFPExp"
table(Idents(P9vP35.combined1))

Idents(object = P9vP35.combined1) <- "OverAllCellType"
P9vP35.combined1$OverAllCellType.ArEGFPExp <- paste(Idents(P9vP35.combined1), P9vP35.combined1$ArEGFPExp, sep = "_")
Idents(object = P9vP35.combined1) <- "OverAllCellType.ArEGFPExp"
table(Idents(P9vP35.combined1))

#### AR EGFP Exp Analysis #### P9only

DoublePos <- subset(x=P9.combined1, subset = Ar > 0.01 & EGFP > 1)
Aronly <- subset(x=P9.combined1, subset = Ar > 0.01 & EGFP < 1)
EGFPonly <- subset(x=P9.combined1, subset = Ar < 0.01 & EGFP > 1)
DoubleNeg <- subset(x=P9.combined1, subset = Ar < 0.01 & EGFP < 1)
Idents(object = DoublePos) <- "DoublePos"
Idents(object = Aronly) <- "Aronly"
Idents(object = EGFPonly) <- "EGFPonly"
Idents(object = DoubleNeg) <- "DoubleNeg"
DoublePos[["ArEGFPExp"]] <- Idents(object = DoublePos)
Aronly[["ArEGFPExp"]] <- Idents(object = Aronly)
EGFPonly[["ArEGFPExp"]] <- Idents(object = EGFPonly)
DoubleNeg[["ArEGFPExp"]] <- Idents(object = DoubleNeg)
Temp4 <- list(DoublePos, DoubleNeg, Aronly, EGFPonly)
ArEGFP1 <- merge(x = DoublePos, y = Aronly)
ArEGFP2 <- merge(x = EGFPonly, y = DoubleNeg)
ArEGFP <- merge(x = ArEGFP1, y = ArEGFP2)
P9.combined1[["ArEGFPExp"]] <- Idents(object = ArEGFP)
DimPlot(P9.combined1, reduction = "umap", group.by = "ArEGFPExp")

Idents(object = P9.combined1) <- "CellClass"
P9.combined1$CellClass.ArEGFPExp <- paste(Idents(P9.combined1), P9.combined1$ArEGFPExp, sep = "_")
Idents(object = P9.combined1) <- "CellClass.ArEGFPExp"
table(Idents(P9.combined1))

Idents(object = PP9.combined1) <- "OverAllCellType"
P9.combined1$OverAllCellType.ArEGFPExp <- paste(Idents(P9.combined1), P9.combined1$ArEGFPExp, sep = "_")
Idents(object = P9.combined1) <- "OverAllCellType.ArEGFPExp"
table(Idents(P9.combined1))

#### AR EGFP Exp Analysis #### P35only

DoublePos <- subset(x=P35.combined1, subset = Ar > 0.01 & EGFP > 1)
Aronly <- subset(x=P35.combined1, subset = Ar > 0.01 & EGFP < 1)
EGFPonly <- subset(x=P35.combined1, subset = Ar < 0.01 & EGFP > 1)
DoubleNeg <- subset(x=P35.combined1, subset = Ar < 0.01 & EGFP < 1)
Idents(object = DoublePos) <- "DoublePos"
Idents(object = Aronly) <- "Aronly"
Idents(object = EGFPonly) <- "EGFPonly"
Idents(object = DoubleNeg) <- "DoubleNeg"
DoublePos[["ArEGFPExp"]] <- Idents(object = DoublePos)
Aronly[["ArEGFPExp"]] <- Idents(object = Aronly)
EGFPonly[["ArEGFPExp"]] <- Idents(object = EGFPonly)
DoubleNeg[["ArEGFPExp"]] <- Idents(object = DoubleNeg)
Temp4 <- list(DoublePos, DoubleNeg, Aronly, EGFPonly)
ArEGFP1 <- merge(x = DoublePos, y = Aronly)
ArEGFP2 <- merge(x = EGFPonly, y = DoubleNeg)
ArEGFP <- merge(x = ArEGFP1, y = ArEGFP2)
P35.combined1[["ArEGFPExp"]] <- Idents(object = ArEGFP)
DimPlot(P35.combined1, reduction = "umap", group.by = "ArEGFPExp")

Idents(object = P35.combined1) <- "CellClass"
P35.combined1$CellClass.ArEGFPExp <- paste(Idents(P35.combined1), P35.combined1$ArEGFPExp, sep = "_")
Idents(object = P35.combined1) <- "CellClass.ArEGFPExp"
table(Idents(P35.combined1))

Idents(object = PP9.combined1) <- "OverAllCellType"
P35.combined1$OverAllCellType.ArEGFPExp <- paste(Idents(P35.combined1),P35.combined1$ArEGFPExp, sep = "_")
Idents(object = P35.combined1) <- "OverAllCellType.ArEGFPExp"
table(Idents(P35.combined1))

#### AR EGFP Exp Analysis #### P9vP35 Epithelium

DoublePos <- subset(x=P9vP35.combined1.Epi1, subset = Ar > 0.01 & EGFP > 1)
Aronly <- subset(x=P9vP35.combined1.Epi1, subset = Ar > 0.01 & EGFP < 1)
EGFPonly <- subset(x=P9vP35.combined1.Epi1, subset = Ar < 0.01 & EGFP > 1)
DoubleNeg <- subset(x=P9vP35.combined1.Epi1, subset = Ar < 0.01 & EGFP < 1)
Idents(object = DoublePos) <- "DoublePos"
Idents(object = Aronly) <- "Aronly"
Idents(object = EGFPonly) <- "EGFPonly"
Idents(object = DoubleNeg) <- "DoubleNeg"
DoublePos[["ArEGFPExp"]] <- Idents(object = DoublePos)
Aronly[["ArEGFPExp"]] <- Idents(object = Aronly)
EGFPonly[["ArEGFPExp"]] <- Idents(object = EGFPonly)
DoubleNeg[["ArEGFPExp"]] <- Idents(object = DoubleNeg)
Temp4 <- list(DoublePos, DoubleNeg, Aronly, EGFPonly)
ArEGFP1 <- merge(x = DoublePos, y = Aronly)
ArEGFP2 <- merge(x = EGFPonly, y = DoubleNeg)
ArEGFP <- merge(x = ArEGFP1, y = ArEGFP2)
P9vP35.combined1.Epi1[["ArEGFPExp"]] <- Idents(object = ArEGFP)
tiff(file = "P9vP35.combined1.Epi1 ArEGFPExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", group.by = "ArEGFPExp", split.by = 'stim')
dev.off()

Idents(object = P9vP35.combined1.Epi1) <- "ArEGFPExp"
P9vP35.combined1.Epi1$ArEGFPExp.stim.EpiCellType <- paste(Idents(P9vP35.combined1.Epi1), P9vP35.combined1.Epi1$stim.EpiCellType, sep = "_")
Idents(object = P9vP35.combined1.Epi1) <- "ArEGFPExp.stim.EpiCellType"
table(Idents(P9vP35.combined1.Epi1))

Idents(object = P9vP35.combined1.Epi1) <- "OverAllCellType"
P9vP35.combined1.Epi1$OverAllCellType.stim.ArEGFPExp <- paste(Idents(P9vP35.combined1.Epi1), P9vP35.combined1.Epi1$stim.ArEGFPExp, sep = "_")
Idents(object = P9vP35.combined1.Epi1) <- "OverAllCellType.stim.ArEGFPExp"
table(Idents(P9vP35.combined1.Epi1))

#### AR EGFP Exp Analysis #### P9vP35 Stromal

DoublePos <- subset(x=P9vP35.combined1.stro1, subset = Ar > 0.01 & EGFP > 1)
Aronly <- subset(x=P9vP35.combined1.stro1, subset = Ar > 0.01 & EGFP < 1)
EGFPonly <- subset(x=P9vP35.combined1.stro1, subset = Ar < 0.01 & EGFP > 1)
DoubleNeg <- subset(x=P9vP35.combined1.stro1, subset = Ar < 0.01 & EGFP < 1)
Idents(object = DoublePos) <- "DoublePos"
Idents(object = Aronly) <- "Aronly"
Idents(object = EGFPonly) <- "EGFPonly"
Idents(object = DoubleNeg) <- "DoubleNeg"
DoublePos[["ArEGFPExp"]] <- Idents(object = DoublePos)
Aronly[["ArEGFPExp"]] <- Idents(object = Aronly)
EGFPonly[["ArEGFPExp"]] <- Idents(object = EGFPonly)
DoubleNeg[["ArEGFPExp"]] <- Idents(object = DoubleNeg)
Temp4 <- list(DoublePos, DoubleNeg, Aronly, EGFPonly)
ArEGFP1 <- merge(x = DoublePos, y = Aronly)
ArEGFP2 <- merge(x = EGFPonly, y = DoubleNeg)
ArEGFP <- merge(x = ArEGFP1, y = ArEGFP2)
P9vP35.combined1.stro1[["ArEGFPExp"]] <- Idents(object = ArEGFP)
tiff(file = "P9vP35.combined1.stro1 ArEGFPExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", group.by = "ArEGFPExp")
dev.off()

Idents(object = P9vP35.combined1.stro1) <- "stim"
P9vP35.combined1.stro1$stim.ArEGFPExp <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$ArEGFPExp, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "stim.ArEGFPExp"
table(Idents(P9vP35.combined1.stro1))

Idents(object = P9vP35.combined1.stro1) <- "OverAllCellType"
P9vP35.combined1.stro1$OverAllCellType.ArEGFPExp <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$ArEGFPExp, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "OverAllCellType.ArEGFPExp"
table(Idents(P9vP35.combined1.stro1))

#DEGs Epithelial Cells 

Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
DimPlot(P9vP35.combined1.Epi1, reduction = "umap")
P9vP35.combined1 <- subset(P9vP35.combined1.Epi1)
DefaultAssay(P9vP35.combined1) <- "RNA"

P9vP35.combined1.Epi1degs.markers <- FindAllMarkers(P9vP35.combined1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P9vP35.combined1.Epi1degs.markers, file = "P9vP35combined1Epi1markers.csv")

#DEGs Stromal Cells
Idents(object = P9vP35.combined1.stro1) <- "StroCellType"
DimPlot(P9vP35.combined1.stro1, reduction = "umap")
P9vP35.combined1 <- subset(P9vP35.combined1.stro1)
DefaultAssay(P9vP35.combined1) <- "RNA"

P9vP35.combined1.stro1degs.markers <- FindAllMarkers(P9vP35.combined1, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P9vP35.combined1.stro1degs.markers, file = "P9vP35.combined1.stro1markers.csv")

#Set up GFP+ & GFP- idents in Epithelial Cells + DEGs

Idents(object = P9vP35.combined1.Epi1) <- "EGFPExp"
EGFPPosP9vP35.combined1.Epi1 <- subset(x=P9vP35.combined1.Epi1, subset = EGFP > 1)
EGFPNegP9vP35.combined1.Epi1 <- subset(x=P9vP35.combined1.Epi1, subset = EGFP < 1)
Idents(object = EGFPPosP9vP35.combined1.Epi1) <- "EGFPpos"
Idents(object = EGFPNegP9vP35.combined1.Epi1) <- "EGFPneg"
EGFPPosP9vP35.combined1.Epi1[["EGFPExp"]] <- Idents(object = EGFPPosP9vP35.combined1.Epi1)
EGFPNegP9vP35.combined1.Epi1[["EGFPExp"]] <- Idents(object = EGFPNegP9vP35.combined1.Epi1)
P9vP35.combined1.Epi1EGFP <- merge(x = EGFPPosP9vP35.combined1.Epi1, y = EGFPNegP9vP35.combined1.Epi1)
P9vP35.combined1.Epi1[["EGFPExp"]] <- Idents(object = P9vP35.combined1.Epi1EGFP)
tiff(file = "P9vP35.combined1.Epi1 EGFPExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", group.by = "EGFPExp")
dev.off()

P9vP35.combined1.Epi1EGFP.markers <- FindAllMarkers(P9vP35.combined1.Epi1EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(P9vP35.combined1.Epi1EGFP.markers, file = "P9vP35combined1Epi1EGFPmarkers.csv")

P9vP35.combined1.Epi1EGFPdegs.markers <- FindMarkers(P9vP35.combined1.Epi1EGFP, ident.1 = "EGFPpos", ident.2 = "EGFPneg", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(P9vP35.combined1.Epi1EGFPdegs.markers, file = "P9vP35.combined1.Epi1EGFPdegsmarkers.csv")

Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
DimPlot(P9vP35.combined1.Epi1, reduction = "umap")
P9vP35.combined1.BE3 <- subset(P9vP35.combined1.Epi1, idents = c("BE3"))
DefaultAssay(P9vP35.combined1.BE3) <- "RNA"
Idents(object = P9vP35.combined1.BE3) <- "EGFPExp"
P9vP35.combined1.BE3EGFPdegs.markers <- FindMarkers(P9vP35.combined1.BE3, ident.1 = "EGFPpos", ident.2 = "EGFPneg", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(P9vP35.combined1.BE3EGFPdegs.markers, file = "P9vP35.combined1.BE3EGFPdegs.markers.csv")

Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
DimPlot(P9vP35.combined1.Epi1, reduction = "umap")
P9vP35.combined1.LE2 <- subset(P9vP35.combined1.Epi1, ident  = c("LE2"))
DefaultAssay(P9vP35.combined1.LE2) <- "RNA"
Idents(object = P9vP35.combined1.LE2) <- "EGFPExp"
P9vP35.combined1.LE2EGFPdegs.markers <- FindMarkers(P9vP35.combined1.LE2, ident.1 = "EGFPpos", ident.2 = "EGFPneg", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(P9vP35.combined1.LE2EGFPdegs.markers, file = "P9vP35.combined1.LE2EGFPdegs.markers.csv")

Idents(object = P9vP35.combined1.Epi1) <- "EpiCellType"
DimPlot(P9vP35.combined1.Epi1, reduction = "umap")
P9vP35.combined1.LE3 <- subset(P9vP35.combined1.Epi1, idents = c("LE3"))
DefaultAssay(P9vP35.combined1.LE3) <- "RNA"
Idents(object = P9vP35.combined1.LE3) <- "EGFPExp"
P9vP35.combined1.LE3EGFPdegs.markers <- FindMarkers(P9vP35.combined1.LE3, ident.1 = "EGFPpos", ident.2 = "EGFPneg", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(P9vP35.combined1.LE3EGFPdegs.markers, file = "P9vP35.combined1.LE3EGFPdegs.markers.csv")

#Set up GFP+ & GFP- idents in Stromal Cells + DEGs

Idents(object = P9vP35.combined1.stro1) <- "EGFPExp"
EGFPPosP9vP35.combined1.stro1 <- subset(x=P9vP35.combined1.stro1, subset = EGFP > 1)
EGFPNegP9vP35.combined1.stro1 <- subset(x=P9vP35.combined1.stro1, subset = EGFP < 1)
Idents(object = EGFPPosP9vP35.combined1.stro1) <- "EGFPpos"
Idents(object = EGFPNegP9vP35.combined1.stro1) <- "EGFPneg"
EGFPPosP9vP35.combined1.stro1[["EGFPExp"]] <- Idents(object = EGFPPosP9vP35.combined1.stro1)
EGFPNegP9vP35.combined1.stro1[["EGFPExp"]] <- Idents(object = EGFPNegP9vP35.combined1.stro1)
P9vP35.combined1.stro1EGFP <- merge(x = EGFPPosP9vP35.combined1.stro1, y = EGFPNegP9vP35.combined1.stro1)
P9vP35.combined1.stro1[["EGFPExp"]] <- Idents(object = P9vP35.combined1.stro1EGFP)
tiff(file = "P9vP35.combined1.stro1 EGFPExp UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", group.by = "EGFPExp")
dev.off()
P9vP35.combined1.stro1EGFP.markers <- FindAllMarkers(P9vP35.combined1.stro1EGFP, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(P9vP35.combined1.stro1EGFP.markers, file = "P9vP35combined1stro1EGFPmarkers.csv")

P9vP35.combined1.stro1EGFPdegs.markers <- FindMarkers(P9vP35.combined1.stro1EGFP, ident.1 = "EGFPpos", ident.2 = "EGFPneg", logfc.threshold = 0.1, min.pct = 0.1)
write.csv(P9vP35.combined1.stro1EGFPdegs.markers, file = "P9vP35.combined1.stro1EGFPdegsmarkers.csv")

#Set up AR+ & AR- idents in Epithelial Cells
DefaultAssay(P9vP35.combined1.Epi1) <- "RNA"
Idents(object = P9vP35.combined1.Epi1) <- "CellType"
P9vP35.combined1.Epi1$CellType.ArExp <- paste(Idents(P9vP35.combined1.Epi1), P9vP35.combined1.Epi1$ArExp, sep = "_")
Idents(object = P9vP35.combined1.Epi1) <- "CellType.ArExp"

ArPosP9vP35.combined1.Epi1 <- subset(x=P9vP35.combined1.Epi1, subset = Ar > 0.01)
ArNegP9vP35.combined1.Epi1 <- subset(x=P9vP35.combined1.Epi1, subset = Ar < 0.01)
Idents(object = ArPosP9vP35.combined1.Epi1) <- "Arpos"
Idents(object = ArNegP9vP35.combined1.Epi1) <- "Arneg"
ArPosP9vP35.combined1.Epi1[["ArExp"]] <- Idents(object = ArPosP9vP35.combined1.Epi1)
ArNegP9vP35.combined1.Epi1[["ArExp"]] <- Idents(object = ArNegP9vP35.combined1.Epi1)
P9vP35.combined1.Epi1AR <- merge(x = ArPosP9vP35.combined1.Epi1, y = ArNegP9vP35.combined1.Epi1)
P9vP35.combined1.Epi1[["ArExp"]] <- Idents(object = P9vP35.combined1.Epi1AR)
tiff(file = "P9vP35.combined1.Epi1 ArExpEpistim UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.Epi1, reduction = "umap", group.by = "ArExp", split.by = 'stim')
dev.off()

Idents(object = P9vP35.combined1.Epi1) <- "ArExp"
P9vP35.combined1.Epi1AR.markers <- FindAllMarkers(P9vP35.combined1.Epi1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P9vP35.combined1.Epi1AR.markers %>% group_by(ArExp) %>% top_n(n = 2, wt = avg_logFC)
P9vP35.combined1.Epi1ARtop10 <- P9vP35.combined1.Epi1AR.markers %>% group_by(ArExp) %>% top_n(n = 10, wt = avg_logFC)

Idents(object = P9vP35.combined1.Epi1) <- "ArExp"
P9vP35.combined1.Epi1$ArExp.stim.EpiCellType <- paste(Idents(P9vP35.combined1.Epi1), P9vP35.combined1.Epi1$stim.EpiCellType, sep = "_")
Idents(object = P9vP35.combined1.Epi1) <- "ArExp.stim.EpiCellType"
table(Idents(P9vP35.combined1.Epi1))

#Set up AR+ & AR- idents in stromal cells
DefaultAssay(P9vP35.combined1.stro1) <- "RNA"
Idents(object = P9vP35.combined1.stro1) <- "CellType"
P9vP35.combined1.stro1$CellType.ArExp <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$ArExp, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "CellType.ArExp"

ArPosP9vP35.combined1.stro1 <- subset(x=P9vP35.combined1.stro1, subset = Ar > 0.01)
ArNegP9vP35.combined1.stro1 <- subset(x=P9vP35.combined1.stro1, subset = Ar < 0.01)
Idents(object = ArPosP9vP35.combined1.stro1) <- "Arpos"
Idents(object = ArNegP9vP35.combined1.stro1) <- "Arneg"
ArPosP9vP35.combined1.stro1[["ArExp"]] <- Idents(object = ArPosP9vP35.combined1.stro1)
ArNegP9vP35.combined1.stro1[["ArExp"]] <- Idents(object = ArNegP9vP35.combined1.stro1)
P9vP35.combined1.stro1AR <- merge(x = ArPosP9vP35.combined1.stro1, y = ArNegP9vP35.combined1.stro1)
P9vP35.combined1.stro1[["ArExp"]] <- Idents(object = P9vP35.combined1.stro1AR)
tiff(file = "P9vP35.combined1.stro1 ArExpStro UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P9vP35.combined1.stro1, reduction = "umap", group.by = "ArExp")
dev.off()
P9vP35.combined1.stro1AR.markers <- FindAllMarkers(P9vP35.combined1.stro1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P9vP35.combined1.stro1AR.markers %>% group_by(ArExp) %>% top_n(n = 2, wt = avg_logFC)
P9vP35.combined1.stro1ARtop10 <- P9vP35.combined1.stro1AR.markers %>% group_by(ArExp) %>% top_n(n = 10, wt = avg_logFC)

Idents(object = P9vP35.combined1.stro1) <- "ArExp"
P9vP35.combined1.stro1$ArExp.stim.StroCellType <- paste(Idents(P9vP35.combined1.stro1), P9vP35.combined1.stro1$stim.StroCellType, sep = "_")
Idents(object = P9vP35.combined1.stro1) <- "ArExp.stim.StroCellType"
table(Idents(P9vP35.combined1.stro1))

