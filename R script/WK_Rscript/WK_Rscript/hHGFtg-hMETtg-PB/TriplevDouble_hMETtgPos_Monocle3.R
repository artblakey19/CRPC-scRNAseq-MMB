####Pseudotime analysis using hMETtg+ cells####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble")

library(Seurat)
library(devtools)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/TriplevDouble.combined1.epi1")

#Rename Clusters
Idents(object = TriplevDouble.combined1.epi1) <- "seurat_clusters"
TriplevDouble.combined1.epi1 <- RenameIdents(object = TriplevDouble.combined1.epi1, 
                                             '18'="BE1", '5'="BE2",'20'="BE2",
                                             '6'="BE3", '13' = "BE4", '14'="BE4", 
                                             '7'="LE1",'12'="LE1", '0'="LE2", '1' = "LE2", '19' ="LE2",
                                             '4' = "LE3", '11'="LE3",
                                             '3' = "LE4", '9'="LE5", '2'="LE6",'15'="LE6",
                                             '10'="LE7", '8'="LE8", '16'="UrLE", 
                                             '17'="OE")  
TriplevDouble.combined1.epi1[["EpiCellTypes"]] <- Idents(object = TriplevDouble.combined1.epi1)

####Add METExp####
DefaultAssay(TriplevDouble.combined1.epi1) <- "RNA"

TriplevDouble.combined1.epi1METPos <- subset(x=TriplevDouble.combined1.epi1,  subset = `hMETtg` > 0)
TriplevDouble.combined1.epi1METNeg <- subset(x=TriplevDouble.combined1.epi1,  subset = `hMETtg` == 0)
Idents(object = TriplevDouble.combined1.epi1METPos) <- "METPos"
Idents(object = TriplevDouble.combined1.epi1METNeg) <- "METNeg"
TriplevDouble.combined1.epi1METPos[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi1METPos)
TriplevDouble.combined1.epi1METNeg[["METExp"]] <- Idents(object = TriplevDouble.combined1.epi1METNeg)
TriplevDouble.combined1.epi1MET <- merge(x = TriplevDouble.combined1.epi1METPos, y = TriplevDouble.combined1.epi1METNeg)
Idents(object = TriplevDouble.combined1.epi1MET) <- "METExp"
TriplevDouble.combined1.epi1$METExp <- Idents(object = TriplevDouble.combined1.epi1MET)
Idents(object = TriplevDouble.combined1.epi1) <- "METExp"
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"))

tiff(file = "TriplevDouble.combined1.epi1 hMETPos highlighted stim UMAP.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, cols = c("red", "lightgrey"), split.by = "stim")
dev.off()

#Subset hMETPos BELE
Idents(object = TriplevDouble.combined1.epi1) <- "METExp"
DimPlot(TriplevDouble.combined1.epi1, reduction = "umap", pt.size = 0.3, label = TRUE)

TriplevDouble.combined1.epi1.METPos <- subset(TriplevDouble.combined1.epi1, idents = c("METPos"))
DimPlot(TriplevDouble.combined1.epi1.METPos, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = TriplevDouble.combined1.epi1.METPos) <- "EpiCellTypes"
TriplevDouble.combined1.BELE.METPos <- subset(TriplevDouble.combined1.epi1.METPos, idents = c("BE1", "BE2", "BE3", "BE4",
                                                                                              "LE1", "LE2", "LE3", "LE4",
                                                                                              "LE5", "LE6", "LE7", "LE8"))
DimPlot(TriplevDouble.combined1.BELE.METPos, reduction = "umap", pt.size = 0.3, label = TRUE)

####hMETtgPos Double####
Idents(object = TriplevDouble.combined1.BELE.METPos) <- "stim"
Double.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELE.METPos, idents = c("Double"))
Triple.combined1.BELEMETPos <- subset(TriplevDouble.combined1.BELE.METPos, idents = c("Triple"))

DefaultAssay(Double.combined1.BELEMETPos) <- "integrated"

#Run the standard workflow for visualization and clustering
Double.combined1.BELEMETPos <- ScaleData(Double.combined1.BELEMETPos, verbose = FALSE)
Double.combined1.BELEMETPos <- RunPCA(Double.combined1.BELEMETPos, npcs = 30, verbose = FALSE)
ElbowPlot(Double.combined1.BELEMETPos, ndims = 50)

#Umap and Clustering
Double.combined1.BELEMETPos <- FindNeighbors(Double.combined1.BELEMETPos, reduction = "pca", dims = 1:20)
Double.combined1.BELEMETPos <- FindClusters(Double.combined1.BELEMETPos, resolution = 0.5)
Double.combined1.BELEMETPos <- RunUMAP(Double.combined1.BELEMETPos, reduction = "pca", dims = 1:20)
DimPlot(Double.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell cycle regression
DefaultAssay(Double.combined1.BELEMETPos) <- "RNA"
all.genes <- rownames(Double.combined1.BELEMETPos)
Double.combined1.BELEMETPos <- ScaleData(Double.combined1.BELEMETPos, features = all.genes)
Double.combined1.BELEMETPos <- CellCycleScoring(Double.combined1.BELEMETPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = Double.combined1.BELEMETPos) <- "Phase"
DimPlot(Double.combined1.BELEMETPos, reduction = "umap")

#Take Cell cycle out 
DefaultAssay(Double.combined1.BELEMETPos) <- "RNA"
Double.combined1.BELEMETPos <- ScaleData(Double.combined1.BELEMETPos, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Double.combined1.BELEMETPos))
Double.combined1.BELEMETPos <- RunPCA(Double.combined1.BELEMETPos, features = VariableFeatures(Double.combined1.BELEMETPos))
ElbowPlot(Double.combined1.BELEMETPos, ndims = 50)

Double.combined1.BELEMETPos <- FindNeighbors(Double.combined1.BELEMETPos, reduction = "pca", dims = 1:20)
Double.combined1.BELEMETPos <- FindClusters(Double.combined1.BELEMETPos, resolution = 0.5)
Double.combined1.BELEMETPos <- RunUMAP(Double.combined1.BELEMETPos, reduction = "pca", dims = 1:20)
DimPlot(Double.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Double.combined1.BELEMETPos) <- "EpiCellTypes"
DimPlot(Double.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE,  
        cols = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                 "aquamarine3", "bisque3", "plum4", 
                 "salmon"))

##Convert Seurat to Monocle3 cell data set class
cds <- as.cell_data_set(Double.combined1.BELEMETPos)

## Calculate size factors using built-in function in monocle3
cds <- estimate_size_factors(cds)

## Add gene names into CDS
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Double.combined1.BELEMETPos)

#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") 

#plot single cell trajectory
cds <- learn_graph(cds, use_partition = FALSE)

#
plt <- plot_cells(cds,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.7,
                  alpha = 1) 

plt + scale_color_manual(values = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                                    "aquamarine3", "bisque3", "plum4", 
                                    "salmon"))

tiff(file = "Double.combined1.BELEMETPos1 EpiCellTypes trajectory black UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                                    "aquamarine3", "bisque3", "plum4", 
                                    "salmon"))
dev.off()

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP")

plt <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.7,
                  alpha = 1)

tiff(file = "Double.combined1.BELEMETPos1 pseudotime trajectory black UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

####plot genes in pseudotime####
Solid_genes <- c("Eif4a1")
Solid_lineage_cds <- cds[rowData(cds)$gene_short_name %in% Solid_genes,
                         colData(cds)$EpiCellTypes %in% c("BE1", "BE2", "BE3", "BE4",
                                                          "LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7", "LE8")]

tiff(file = "Double.combined1.BELEPos1 Eif4a1 in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes",
                         min_expr=10) + scale_color_manual(values = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                                                                      "aquamarine3", "bisque3", "plum4", 
                                                                      "salmon"))
dev.off()

####hMETtgPos Triple####
#Cell cycle regression
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
all.genes <- rownames(Triple.combined1.BELEMETPos)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, features = all.genes)
Triple.combined1.BELEMETPos <- CellCycleScoring(Triple.combined1.BELEMETPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = Triple.combined1.BELEMETPos) <- "Phase"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap")

#Take Cell cycle out 
Triple.combined1.BELEMETPos1 <- Triple.combined1.BELEMETPos
DefaultAssay(Triple.combined1.BELEMETPos1) <- "RNA"
Triple.combined1.BELEMETPos1 <- ScaleData(Triple.combined1.BELEMETPos1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.BELEMETPos1))
Triple.combined1.BELEMETPos1 <- RunPCA(Triple.combined1.BELEMETPos1, features = VariableFeatures(Triple.combined1.BELEMETPos1))
ElbowPlot(Triple.combined1.BELEMETPos1, ndims = 50)

#Res0.5
Triple.combined1.BELEMETPos1 <- FindNeighbors(Triple.combined1.BELEMETPos1, reduction = "pca", dims = 1:15)
Triple.combined1.BELEMETPos1 <- FindClusters(Triple.combined1.BELEMETPos1, resolution = 0.5)
Triple.combined1.BELEMETPos1 <- RunUMAP(Triple.combined1.BELEMETPos1, reduction = "pca", dims = 1:15)
DimPlot(Triple.combined1.BELEMETPos1, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Triple.combined1.BELEMETPos1) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELEMETPos1, reduction = "umap", pt.size = 0.3, label = TRUE,  
        cols = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                 "aquamarine3", "bisque3", "plum4", 
                 "salmon"))

##Convert Seurat to Monocle3 cell data set class
cds1 <- as.cell_data_set(Triple.combined1.BELEMETPos1)

## Calculate size factors using built-in function in monocle3
cds1 <- estimate_size_factors(cds1)

## Add gene names into CDS
cds1@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(Triple.combined1.BELEMETPos1)

#Construction single cell trajectories
cds1 <- cluster_cells(cds1 = cds1, reduction_method = "UMAP") 

#plot single cell trajectory
cds1 <- learn_graph(cds1, use_partition = FALSE)

#
plt <- plot_cells(cds1,
                  color_cells_by = "EpiCellTypes",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.7,
                  alpha = 1) 

plt + scale_color_manual(values = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                                    "aquamarine3", "bisque3", "plum4", 
                                    "salmon"))

tiff(file = "Triple.combined1.BELEMETPos1 EpiCellTypes trajectory black UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt + scale_color_manual(values = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                                    "aquamarine3", "bisque3", "plum4", 
                                    "salmon"))
dev.off()

#plot single cell trajectory by pseudotime
cds1 <- order_cells(cds1, reduction_method = "UMAP")

plt <- plot_cells(cds1,
                  color_cells_by = "pseudotime",
                  show_trajectory_graph = TRUE,
                  label_branch_points = FALSE,
                  trajectory_graph_color = "black",
                  trajectory_graph_segment_size = 1.8,
                  label_cell_groups = FALSE,
                  label_leaves=FALSE,
                  label_roots = FALSE,
                  cell_size = 0.7,
                  alpha = 1)

tiff(file = "Triple.combined1.BELEMETPos1 pseudotime trajectory black UMAP.tiff", width = 6, height = 5, units = "in", compression = "lzw", res = 800)
plt
dev.off()

####plot genes in pseudotime####
Solid_genes <- c("Eif4a1")
Solid_lineage_cds <- cds1[rowData(cds1)$gene_short_name %in% Solid_genes,
                         colData(cds1)$EpiCellTypes %in% c("BE1", "BE2", "BE3", "BE4",
                                                          "LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7", "LE8")]

tiff(file = "Triple.combined1.BELEMETPos1 Eif4a1 in pseudotime.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 1000)
plot_genes_in_pseudotime(Solid_lineage_cds, cell_size = 1.5,
                         color_cells_by="EpiCellTypes",
                         min_expr=10) + scale_color_manual(values = c("orchid", "skyblue1", "chartreuse3", "brown3", "deeppink1", "khaki4",  "blue",  "red", 
                                                                      "aquamarine3", "bisque3", "plum4", 
                                                                      "salmon"))
dev.off()

