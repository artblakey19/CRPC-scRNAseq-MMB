#MycCtrlvMycARKO SCseq Work Flow

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

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Ctrl2vARKO2")

#### Merging Datasets ####

#Stash old idents
Ctrl[["orig.clusters"]] <- Idents(object = Ctrl)
ARKO[["orig.clusters"]] <- Idents(object = ARKO)

#Set Current idents
Idents(object = Ctrl) <- "seurat_clusters"
Idents(object = ARKO) <- "seurat_clusters"
Ctrl$stim <- "Ctrl"
ARKO$stim <- "ARKO"
CtrlvARKO.anchors <- FindIntegrationAnchors(object.list = list(Ctrl, ARKO), dims = 1:20)
CtrlvARKO.combined <- IntegrateData(anchorset = CtrlvARKO.anchors, dims = 1:20)
DefaultAssay(CtrlvARKO.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined <- ScaleData(CtrlvARKO.combined, verbose = FALSE)
CtrlvARKO.combined <- RunPCA(CtrlvARKO.combined, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined, ndims = 50)

# t-SNE and Clustering
CtrlvARKO.combined <- FindNeighbors(CtrlvARKO.combined, reduction = "pca", dims = 1:25)
CtrlvARKO.combined <- FindClusters(CtrlvARKO.combined, resolution = 0.5)
CtrlvARKO.combined <- RunTSNE(CtrlvARKO.combined, reduction = "pca", dims = 1:25)
CtrlvARKO.combined <- RunUMAP(CtrlvARKO.combined, reduction = "pca", dims = 1:25)

DimPlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", pt.size = 0.3) 
tiff(file = "Combined UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)  
dev.off()
tiff(file = "Combined UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)  
dev.off()

#FeaturePlot
DefaultAssay(CtrlvARKO.combined) <- "RNA"
tiff(file = "combined Krt5.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Krt8.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Epcam.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Epcam"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Vim.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Vim"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Myh11.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Myh11"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined Fbln1.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined, reduction = "umap", features = c("Fbln1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Reclustering Epithelial cells ####

Idents(object = CtrlvARKO.combined) <- "seurat_clusters"
CtrlvARKO.combined.Epi <- subset(CtrlvARKO.combined, idents = c("0", "1", "4", "8", "11", "15", "17", "6", "10", "19"))
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"
DefaultAssay(CtrlvARKO.combined.Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi <- ScaleData(CtrlvARKO.combined.Epi, verbose = FALSE)
CtrlvARKO.combined.Epi <- RunPCA(CtrlvARKO.combined.Epi, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi, ndims = 50)

#Umap and Clustering
CtrlvARKO.combined.Epi <- FindNeighbors(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:17)
CtrlvARKO.combined.Epi <- FindClusters(CtrlvARKO.combined.Epi, resolution = 0.5)
CtrlvARKO.combined.Epi <- RunTSNE(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:17)
CtrlvARKO.combined.Epi <- RunUMAP(CtrlvARKO.combined.Epi, reduction = "pca", dims = 1:17)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", label = TRUE)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", split.by = "stim")
tiff(file = "Combined.Epi UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", pt.size = 0.3, label = TRUE)  
dev.off()
tiff(file = "Combined.Epi UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)  
dev.off()

#FeaturePlot
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"
tiff(file = "combined.Epi MYC-transgene split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c('MYC-transgene'), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt5 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt8 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Pbsn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Nkx3.1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c('Nkx3-1'), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "combined.Epi Krt6a split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi, reduction = "umap", features = c("Krt6a"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell counts
Idents(object = CtrlvARKO.combined.Epi) <- "stim"
CtrlvARKO.combined.Epi$stim.seurat_clusters <- paste(Idents(CtrlvARKO.combined.Epi), CtrlvARKO.combined.Epi$seurat_clusters, sep = "_")
Idents(object = CtrlvARKO.combined.Epi) <- "stim.seurat_clusters"
table(Idents(CtrlvARKO.combined.Epi))

#### Myc+ CtrlvARKO in Integrated ####
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"

#Add MYC info
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"

CtrlvARKO.combined.Epi2 <- subset(x=CtrlvARKO.combined.Epi, subset = `MYC-transgene` > 0)
DefaultAssay(CtrlvARKO.combined.Epi2) <- "integrated"

#Run the standard workflow for visualization and clustering
CtrlvARKO.combined.Epi2 <- ScaleData(CtrlvARKO.combined.Epi2, verbose = FALSE)
CtrlvARKO.combined.Epi2 <- RunPCA(CtrlvARKO.combined.Epi2, npcs = 30, verbose = FALSE)
ElbowPlot(CtrlvARKO.combined.Epi2, ndims = 50)
# t-SNE and Clustering
CtrlvARKO.combined.Epi2 <- FindNeighbors(CtrlvARKO.combined.Epi2, reduction = "pca", dims = 1:22)
CtrlvARKO.combined.Epi2 <- FindClusters(CtrlvARKO.combined.Epi2, resolution = 0.5)
CtrlvARKO.combined.Epi2 <- RunTSNE(CtrlvARKO.combined.Epi2, reduction = "pca", dims = 1:22)
CtrlvARKO.combined.Epi2 <- RunUMAP(CtrlvARKO.combined.Epi2, reduction = "pca", dims = 1:22)

DimPlot(CtrlvARKO.combined.Epi2, reduction = "umap", split.by = "stim", pt.size = 0.3) 

tiff(file = "MYCPos.Combined.Epi UMAP split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(CtrlvARKO.combined.Epi2, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)  
dev.off()

#FeaturePlot
DefaultAssay(CtrlvARKO.combined.Epi2) <- "RNA"
tiff(file = "MYCPos.combined.Epi MYC-transgene split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c('MYC-transgene'), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPos.combined.Epi Ar split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Ar"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPos.combined.Epi Krt5 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt5"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPos.combined.Epi Krt8 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Krt8"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPos.combined.Epi Pbsn split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c("Pbsn"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()
tiff(file = "MYCPos.combined.Epi Nkx3.1 split.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(CtrlvARKO.combined.Epi2, reduction = "umap", features = c('Nkx3-1'), split.by = "stim", cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

#DEGs
DefaultAssay(CtrlvARKO.combined.Epi2) <- "RNA"
Idents(object = CtrlvARKO.combined.Epi2) <- "stim"
CtrlvARKO.combined.Epi2 <- ScaleData(CtrlvARKO.combined.Epi2, features = rownames(CtrlvARKO.combined.Epi2))
CtrlvARKO.combined.Epi2.markers <- FindAllMarkers(CtrlvARKO.combined.Epi2, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(CtrlvARKO.combined.Epi2.markers, "CtrlvARKO.combined.Epi2.Markers.csv")

#### Pseudotime of Epi ####
Idents(object = CtrlvARKO.combined.Epi) <- "seurat_clusters"
DefaultAssay(CtrlvARKO.combined.Epi) <- "RNA"
EpiPseudo <- as.CellDataSet(CtrlvARKO.combined.Epi)
EpiPseudo <- detectGenes(EpiPseudo, min_expr = 0.1)
print(head(fData(EpiPseudo)))

expressed_genes <- row.names(subset(fData(EpiPseudo),
                                    num_cells_expressed >= 10))

pData(EpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo))
EpiPseudo <- EpiPseudo[,pData(EpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo), geom =
        "density")

EpiPseudo <- estimateSizeFactors(EpiPseudo)
EpiPseudo <- estimateDispersions(EpiPseudo)

disp_table <- dispersionTable(EpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo <- setOrderingFilter(EpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo)

#EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo, return_all = F) # norm_method='log'

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
EpiPseudo <- clusterCells(EpiPseudo, num_clusters = 2)

plot_cell_clusters(EpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(EpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
EpiPseudo <- setOrderingFilter(EpiPseudo, ordering_genes)
plot_ordering_genes(EpiPseudo)

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2,
                             method = 'DDRTree')

EpiPseudo <- orderCells(EpiPseudo)

GM_state <- function(EpiPseudo){
  if (length(unique(pData(EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo)$State, pData(EpiPseudo)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo <- orderCells(EpiPseudo, root_state = GM_state(EpiPseudo))

#Visualization
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

tiff(file = "EpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "EpiPseudo Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1) 
dev.off()
tiff(file = "EpiPseudo stim.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()

tiff(file = "EpiPseudo seurat.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "EpiPseudo seurat split. tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()

tiff(file = "EpiPseudo MYC.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c('MYC-transgene'), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Ar.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Ar"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Krt5.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Krt5"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Krt8.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Krt8"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Pbsn.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Pbsn"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Nkx3.1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c('Nkx3-1'), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()

tiff(file = "EpiPseudo Ccnd1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Ccnd1"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Cd44.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Cd44"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red")) + facet_wrap(~stim, nrow = 1)
dev.off()

tiff(file = "EpiPseudo Src.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Src"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "darkgreen")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo Fosl1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Fosl1"), cell_size = 1,  use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "darkgreen")) + facet_wrap(~stim, nrow = 1)
dev.off()


tiff(file = "EpiPseudo Wntgenes red.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Myc", "Ccnd1", "Cd44", "Plaur", "Tcf712", "Ctnnb1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "EpiPseudo IGF signaling darkgreen.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Igf1r", "Src", "Jak1", "Fosl1"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "darkgreen"))
dev.off()
