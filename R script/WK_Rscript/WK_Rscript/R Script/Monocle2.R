library(monocle)
library(RColorBrewer)

Idents(object = PINvTumor.combined.Epi2) <- "ARQExp"
ARQPosEpi <- subset(PINvTumor.combined.Epi2, idents = c("ARQPos"))
ARQNegEpi <- subset(PINvTumor.combined.Epi2, idents = c("ARQNeg"))

#ARQPosEpi
Idents(object = ARQPosEpi) <- "seurat_clusters"
DimPlot(ARQPosEpi, reduction = "tsne")


DefaultAssay(ARQPosEpi) <- "RNA"
EpiPseudo <- as.CellDataSet(ARQPosEpi)
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
ordering_genes <- row.names (subset(diff_test_res, qval < 0.1))
EpiPseudo <- setOrderingFilter(EpiPseudo, ordering_genes)
plot_ordering_genes(EpiPseudo)

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2,
                             method = 'DDRTree')

EpiPseudo <- orderCells(EpiPseudo)

GM_state <- function(EpiPseudo){
  if (length(unique(pData(EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo)$State, pData(EpiPseudo)$seurat_clusters)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo <- orderCells(EpiPseudo, root_state = GM_state(EpiPseudo))

#Visualization
tiff(file = "EpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "EpiPseudo Pseudotime split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo EpiCellType.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values = c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) 
dev.off()
tiff(file = "EpiPseudo stim EpiCellType stim split.tiff", width = 16, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "EpiCellType", show_branch_points = FALSE) + scale_color_manual(values=c("#1D762E", "purple", "red", "#FF9933", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "EpiPseudo CellMarkers Wnt targets purple.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt5", "Krt14", "Krt8", "Krt18", "Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo CellMarkers Wnt targets.tiff", width = 16, height = 16, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Pbsn", "Krt5", "Krt14", "Krt8", "Krt18", "Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "EpiPseudo Wnt targets.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "EpiPseudo Wnt targets redcolor.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Tcf4", "Ccnd1", "Axin2", "Lgr5"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "EpiPseudo IGF1R Wnt stem cell targets.tiff", width = 16, height = 21, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("Igf1r", "Mapk13", "Jak2", "Tcf4", "Ccnd1", "Axin2", "Lgr5", "Trp63", "Ly6a", "Psca", "Tacstd2", "Itga6", "Osr1", "Pbsn"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "EpiPseudo celltype IGF1R Wnt stem cell targets.tiff", width = 16, height = 21, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = c("ARQ", "EGFP", "Ar", "Krt5", "Krt8", "Igf1r", "Mapk13", "Jak2", "Tcf4", "Ccnd1", "Axin2", "Lgr5", "Trp63", "Ly6a", "Psca", "Tacstd2", "Itga6", "Osr1", "Pbsn"), use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
