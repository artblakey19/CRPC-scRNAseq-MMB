#### load each samples ####
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/E18.5_only_Analysis/E18only_Rworkspace.RData")
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P11_only_Analysis/P11only_Rworkspace.RData")
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P42_only_Analysis/P42only_Rworkspace.RData")

#### Subclustering E18.5 UGE ####

Idents(object = E18_CtrlvARKO.combined) <- "CellTypes"
Combined_E18_UGE <- subset(E18_CtrlvARKO.combined, idents = c("UGE"))

Idents(object = Combined_E18_UGE) <- "stim"
DimPlot(Combined_E18_UGE, reduction = "umap", pt.size = 0.3, label = TRUE) 
E18_Ctrl_UGE <- subset(Combined_E18_UGE, idents = c("E18_Ctrl"))
E18_ARKO_UGE <- subset(Combined_E18_UGE, idents = c("E18_ARKO"))

#### Subclustering P11 Epi ####

Idents(object = P11_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P11_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)
Combined_P11_Epi <- subset(P11_CtrlvARKO.combined, idents = c("Basal", "Intermediate", "Luminal"))

Idents(object = Combined_P11_Epi) <- "stim"
DimPlot(Combined_P11_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 
P11_Ctrl_Epi <- subset(Combined_P11_Epi, idents = c("P11_Ctrl"))
P11_ARKO_Epi <- subset(Combined_P11_Epi, idents = c("P11_ARKO"))

#### Subclustering P42 Epi ####

Idents(object = P42_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P42_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)
Combined_P42_Epi <- subset(P42_CtrlvARKO.combined, idents = c("Basal", "Prolif Epi", "Luminal"))

Idents(object = Combined_P42_Epi) <- "stim"
DimPlot(Combined_P42_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 
P42_Ctrl_Epi <- subset(Combined_P42_Epi, idents = c("P42_Ctrl"))
P42_ARKO_Epi <- subset(Combined_P42_Epi, idents = c("P42_ARKO"))

####Integrate EGFP+ FB/SM of E18.5 P11 P42####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 Epi")

#Set Current idents
Idents(object = E18_Ctrl_UGE) <- "seurat_clusters"
Idents(object = E18_ARKO_UGE) <- "seurat_clusters"
Idents(object = P11_Ctrl_Epi) <- "seurat_clusters"
Idents(object = P11_ARKO_Epi) <- "seurat_clusters"
Idents(object = P42_Ctrl_Epi) <- "seurat_clusters"
Idents(object = P42_ARKO_Epi) <- "seurat_clusters"
E18_Ctrl_UGE$stim <- "E18_Ctrl"
E18_ARKO_UGE$stim <- "E18_ARKO"
P11_Ctrl_Epi$stim <- "P11_Ctrl"
P11_ARKO_Epi$stim <- "P11_ARKO"
P42_Ctrl_Epi$stim <- "P42_Ctrl"
P42_ARKO_Epi$stim <- "P42_ARKO"
ARKOvCtrl_Epi.anchors <- FindIntegrationAnchors(object.list = list(E18_Ctrl_UGE, E18_ARKO_UGE, P11_Ctrl_Epi, P11_ARKO_Epi, P42_Ctrl_Epi, P42_ARKO_Epi), dims = 1:20)
ARKOvCtrl_Epi <- IntegrateData(anchorset = ARKOvCtrl_Epi.anchors, dims = 1:20)
DefaultAssay(ARKOvCtrl_Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
ARKOvCtrl_Epi <- ScaleData(ARKOvCtrl_Epi, verbose = FALSE)
ARKOvCtrl_Epi <- RunPCA(ARKOvCtrl_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_Epi)

# umap and Clustering
ARKOvCtrl_Epi <- FindNeighbors(ARKOvCtrl_Epi, reduction = "pca", dims = 1:15)
ARKOvCtrl_Epi <- FindClusters(ARKOvCtrl_Epi, resolution = 0.5)
ARKOvCtrl_Epi <- RunUMAP(ARKOvCtrl_Epi, reduction = "pca", dims = 1:15)
ARKOvCtrl_Epi <- RunTSNE(ARKOvCtrl_Epi, reduction = "pca", dims = 1:15)

tiff(file = "ARKOvCtrl_Epi UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_Epi stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

Idents(object = ARKOvCtrl_Epi) <- "CellTypes"
tiff(file = "ARKOvCtrl_Epi Celltype UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_Epi celltype split UMAP.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#Cell count (split)
Idents(object = ARKOvCtrl_Epi) <- "stim"
ARKOvCtrl_Epi$stim.seurat_clusters <- paste(Idents(ARKOvCtrl_Epi), ARKOvCtrl_Epi$seurat_clusters, sep = "_")
Idents(object = ARKOvCtrl_Epi) <- "stim.seurat_clusters"
table(Idents(ARKOvCtrl_Epi))

#Featureplot
DefaultAssay(ARKOvCtrl_Epi) <- "RNA"
tiff(file = "ARKOvCtrl_Epi split Krt5.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi split Krt8.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi split Pbsn.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi split Mki67.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "ARKOvCtrl_Epi  Krt5.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Krt5"), cols = c("light grey", "red"),  min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi Krt8.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Krt8"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi  Pbsn.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Pbsn"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi Mki67.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()




#Heatmap
Idents(object = ARKOvCtrl_Epi) <- "seurat_clusters"
DefaultAssay(ARKOvCtrl_Epi) <- "RNA"
ARKOvCtrl_Epi <- ScaleData(ARKOvCtrl_Epi, features = rownames(ARKOvCtrl_Epi))
ARKOvCtrl_Epi.markers <- FindAllMarkers(ARKOvCtrl_Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ARKOvCtrl_EpiTop10 <- ARKOvCtrl_Epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
tiff(file = "ARKOvCtrl_Epi Heatmap Top10.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(ARKOvCtrl_Epi, features = c(ARKOvCtrl_EpiTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
ARKOvCtrl_Epi.0.1markers <- FindAllMarkers(ARKOvCtrl_Epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(ARKOvCtrl_Epi.0.1markers, "ARKOvCtrl_Epi.0.1markers.csv")

tiff(file = "ARKOvCtrl_Epi split Vim.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi split Myh11.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi Vim.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi Myh11.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#cluster 1 markers
tiff(file = "ARKOvCtrl_Epi Shh.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Shh"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi Pdgfa.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Pdgfa"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi Fgfr2.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Fgfr2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "ARKOvCtrl_Epi stim Shh.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Shh"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi stim Pdgfa.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Pdgfa"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi stim Fgfr2.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Fgfr2"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#cluster 3 markers
tiff(file = "ARKOvCtrl_Epi Sdc3.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Sdc3"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi  Itgav.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Itgav"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "ARKOvCtrl_Epi stim Sdc3.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Sdc3"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi stim Itgav.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Itgav"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()


#squamous markers
tiff(file = "ARKOvCtrl_Epi Krt10.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi, reduction = "umap", features = c("Krt10"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Pseudotime of 6 EGFP_Epi ####
Idents(object = ARKOvCtrl_Epi) <- "stim"
DefaultAssay(ARKOvCtrl_Epi) <- "RNA"
ARKOvCtrl_EpiPseudo <- as.CellDataSet(ARKOvCtrl_Epi)
ARKOvCtrl_EpiPseudo  <- detectGenes(ARKOvCtrl_EpiPseudo, min_expr = 0.1)
print(head(fData(ARKOvCtrl_EpiPseudo)))

expressed_genes <- row.names(subset(fData(ARKOvCtrl_EpiPseudo),
                                    num_cells_expressed >= 10))

pData(ARKOvCtrl_EpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARKOvCtrl_EpiPseudo))
ARKOvCtrl_EpiPseudo <- ARKOvCtrl_EpiPseudo[,pData(ARKOvCtrl_EpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARKOvCtrl_EpiPseudo), geom =
        "density")

ARKOvCtrl_EpiPseudo <- estimateSizeFactors(ARKOvCtrl_EpiPseudo)
ARKOvCtrl_EpiPseudo <- estimateDispersions(ARKOvCtrl_EpiPseudo)

disp_table <- dispersionTable(ARKOvCtrl_EpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARKOvCtrl_EpiPseudo <- setOrderingFilter(ARKOvCtrl_EpiPseudo, unsup_clustering_genes$gene_id)
#plot_ordering_genes(ARKOvCtrl_EpiPseudo)

#ARKOvCtrl_EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(ARKOvCtrl_EpiPseudo, return_all = F) # norm_method='log'

ARKOvCtrl_EpiPseudo <- reduceDimension(ARKOvCtrl_EpiPseudo, max_components = 2, num_dim = 15,
                                         reduction_method = 'tSNE', verbose = T)
ARKOvCtrl_EpiPseudo <- clusterCells(ARKOvCtrl_EpiPseudo, num_clusters = 2)

plot_cell_clusters(ARKOvCtrl_EpiPseudo, color_by = "stim")

diff_test_res <- differentialGeneTest(ARKOvCtrl_EpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~stim")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARKOvCtrl_EpiPseudo <- setOrderingFilter(ARKOvCtrl_EpiPseudo, ordering_genes)
plot_ordering_genes(ARKOvCtrl_EpiPseudo)

ARKOvCtrl_EpiPseudo <- reduceDimension(ARKOvCtrl_EpiPseudo, max_components = 2,
                                         method = 'DDRTree')
ARKOvCtrl_EpiPseudo <- orderCells(ARKOvCtrl_EpiPseudo)

GM_state <- function(ARKOvCtrl_EpiPseudo){
  if (length(unique(pData(ARKOvCtrl_EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(ARKOvCtrl_EpiPseudo)$State, pData(ARKOvCtrl_EpiPseudo)$stim)[,"E18_ARKO"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARKOvCtrl_EpiPseudo <- orderCells(ARKOvCtrl_EpiPseudo, root_state = GM_state(ARKOvCtrl_EpiPseudo))

plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)


tiff(file = "ARKOvCtrl_EpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_EpiPseudo seurat_clusters Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "ARKOvCtrl_EpiPseudo stim Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_EpiPseudo CellTypes Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()


tiff(file = "ARKOvCtrl_EpiPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "ARKOvCtrl_EpiPseudo stim split Pseudotime.tiff", width = 6, height = 9, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 3)
dev.off()

plot_cell_trajectory(ARKOvCtrl_EpiPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))


####Eliminate Cdh1-Vim+ cells####

Idents(object = ARKOvCtrl_Epi) <- "seurat_clusters"
DefaultAssay(ARKOvCtrl_Epi) <- "RNA"

ARKOvCtrl_Epi <- subset(x=ARKOvCtrl_Epi, subset = Cdh1 > 0)

DefaultAssay(ARKOvCtrl_Epi) <- "integrated"

#Run the standard workflow for visualization and clustering
ARKOvCtrl_Epi <- ScaleData(ARKOvCtrl_Epi, verbose = FALSE)
ARKOvCtrl_Epi <- RunPCA(ARKOvCtrl_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_Epi)

# umap and Clustering
ARKOvCtrl_Epi <- FindNeighbors(ARKOvCtrl_Epi, reduction = "pca", dims = 1:15)
ARKOvCtrl_Epi <- FindClusters(ARKOvCtrl_Epi, resolution = 0.5)
ARKOvCtrl_Epi <- RunUMAP(ARKOvCtrl_Epi, reduction = "pca", dims = 1:15)
ARKOvCtrl_Epi <- RunTSNE(ARKOvCtrl_Epi, reduction = "pca", dims = 1:15)

tiff(file = "ARKOvCtrl_Epi UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_Epi stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

Idents(object = ARKOvCtrl_Epi) <- "CellTypes"
tiff(file = "ARKOvCtrl_Epi Celltype UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_Epi celltype split UMAP.tiff", width = 24, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Epi, reduction = "umap", pt.size = 0.3, split.by = "stim")
dev.off()

#### Add State info to Seurat object ####
Combined_ARKO_Epi$PseudoState <- ARKO_Epi_Pseudo$State
DimPlot(Combined_ARKO_Epi, reduction = "umap", group.by = "PseudoState")
Idents(object = Combined_ARKO_Epi) <- "PseudoState"

#DEGs
setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/Combined_ARKO")
Combined_ARKO_Epi_State.0.1markers <- FindAllMarkers(Combined_ARKO_Epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(Combined_ARKO_Epi_State.0.1markers, "Combined_ARKO_Epi_State.0.1markers.csv")
