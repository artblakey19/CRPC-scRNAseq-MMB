library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

####E18.5 EGFP+ FB/SM w/o ProS > Separate Ctrl & ARKO ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/GFP+/without ProS")

#subclustering EGFP+
Idents(object = E18_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined) <- "RNA"

E18_CtrlvARKO.combined_EGFP <- subset(x=E18_CtrlvARKO.combined, subset = EGFP > 0)

#Add Ar info
DefaultAssay(E18_CtrlvARKO.combined_EGFP) <- "RNA"
E18_CtrlvARKO.combined_EGFP_ArPos <- subset(x=E18_CtrlvARKO.combined_EGFP, subset = Ar > 0)
E18_CtrlvARKO.combined_EGFP_ArNeg <- subset(x=E18_CtrlvARKO.combined_EGFP, subset = Ar == 0)
Idents(object = E18_CtrlvARKO.combined_EGFP_ArPos) <- "ArPos"
Idents(object = E18_CtrlvARKO.combined_EGFP_ArNeg) <- "ArNeg"
E18_CtrlvARKO.combined_EGFP_ArPos[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_EGFP_ArPos)
E18_CtrlvARKO.combined_EGFP_ArNeg[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_EGFP_ArNeg)
E18_CtrlvARKO.combined_EGFP_Ar <- merge(x = E18_CtrlvARKO.combined_EGFP_ArPos, y = E18_CtrlvARKO.combined_EGFP_ArNeg)
Idents(object = E18_CtrlvARKO.combined_EGFP_Ar) <- "ArExp"
E18_CtrlvARKO.combined_EGFP$ArExp <- Idents(object = E18_CtrlvARKO.combined_EGFP_Ar)

Idents(object = E18_CtrlvARKO.combined_EGFP) <- "ArExp"
E18_CtrlvARKO.combined_EGFP <- subset(E18_CtrlvARKO.combined_EGFP, idents = c("ArPos", "ArNeg"))
E18_CtrlvARKO.combined_EGFP <- RenameIdents(object = E18_CtrlvARKO.combined_EGFP, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(E18_CtrlvARKO.combined_EGFP, reduction = "umap", pt.size = 0.3)

Idents(object = E18_CtrlvARKO.combined_EGFP) <- "stim"
E18_CtrlvARKO.combined_EGFP$stim.ArExp <- paste(Idents(E18_CtrlvARKO.combined_EGFP), E18_CtrlvARKO.combined_EGFP$ArExp, sep = "_")
Idents(object = E18_CtrlvARKO.combined_EGFP) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = E18_CtrlvARKO.combined_EGFP) <- "CellTypes"
E18_CtrlvARKO.combined_EGFP_FibSM <- subset(E18_CtrlvARKO.combined_EGFP, idents = c("Fibroblast", "SM"))
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"
E18_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(E18_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
E18_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(E18_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(E18_CtrlvARKO.combined_EGFP_FibSM, ndims = 50)

E18_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
E18_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(E18_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
E18_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
E18_CtrlvARKO.combined_EGFP_FibSM <- RunTSNE(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Arstim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
E18_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(E18_CtrlvARKO.combined_EGFP_FibSM), E18_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(E18_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(E18_CtrlvARKO.combined_EGFP_FibSM, features = rownames(E18_CtrlvARKO.combined_EGFP_FibSM))
E18_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(E18_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E18_CtrlvARKO.combined_EGFP_FibSMTop10 <- E18_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "E18_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(E18_CtrlvARKO.combined_EGFP_FibSM, features = c(E18_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
E18_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(E18_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E18_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "E18_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")

#Separating Ctrl & ARKO
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "stim"
E18_Ctrl_FBSM <- subset(E18_CtrlvARKO.combined_EGFP_FibSM, idents = c("E18_Ctrl"))
E18_ARKO_FBSM <- subset(E18_CtrlvARKO.combined_EGFP_FibSM, idents = c("E18_ARKO"))

#### Cytotrace ####
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18_CtrlvARKO.combined_EGFP_FibSMExpdata = data.frame(E18_CtrlvARKO.combined_EGFP_FibSM[["RNA"]]@data)
E18_CtrlvARKO.combined_EGFP_FibSMphenotypedata = FetchData(E18_CtrlvARKO.combined_EGFP_FibSM, c("seurat_clusters"))
write.table(E18_CtrlvARKO.combined_EGFP_FibSMExpdata, file = "FibSMExpdata.txt")
write.table(E18_CtrlvARKO.combined_EGFP_FibSMphenotypedata, file = "FibSMphenotypedata.txt")

#### Pseudotime of combined_E18_EGFP_FibSM w/o ProS ####

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18EGFPFibSMPseudo <- as.CellDataSet(E18_CtrlvARKO.combined_EGFP_FibSM)
E18EGFPFibSMPseudo  <- detectGenes(E18EGFPFibSMPseudo, min_expr = 0.1)
print(head(fData(E18EGFPFibSMPseudo)))

expressed_genes <- row.names(subset(fData(E18EGFPFibSMPseudo),
                                    num_cells_expressed >= 10))

pData(E18EGFPFibSMPseudo)$Total_mRNAs <- Matrix::colSums(exprs(E18EGFPFibSMPseudo))
pData(E18EGFPFibSMPseudo)$Pseudotime <- max(pData(E18EGFPFibSMPseudo)$Pseudotime) - pData(E18EGFPFibSMPseudo)$Pseudotime
E18EGFPFibSMPseudo <- E18EGFPFibSMPseudo[,pData(E18EGFPFibSMPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(E18EGFPFibSMPseudo), geom =
        "density")

E18EGFPFibSMPseudo <- estimateSizeFactors(E18EGFPFibSMPseudo)
E18EGFPFibSMPseudo <- estimateDispersions(E18EGFPFibSMPseudo)

disp_table <- dispersionTable(E18EGFPFibSMPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, unsup_clustering_genes$gene_id)
#plot_ordering_genes(E18EGFPFibSMPseudo)

#E18EGFPFibSMPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(E18EGFPFibSMPseudo, return_all = F) # norm_method='log'

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2, num_dim = 20,
                                      reduction_method = 'tSNE', verbose = T)
E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, num_clusters = 2)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "seurat_clusters")
plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "stim")
plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "CellTypes")

diff_test_res <- differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, ordering_genes)
plot_ordering_genes(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2,
                                      method = 'DDRTree')

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo)

GM_state <- function(E18EGFPFibSMPseudo){
  if (length(unique(pData(E18EGFPFibSMPseudo)$State)) > 1){
    T0_counts <- table(pData(E18EGFPFibSMPseudo)$State, pData(E18EGFPFibSMPseudo)$seurat_clusters)[,"7"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo, root_state = GM_state(E18EGFPFibSMPseudo), reverse = TRUE)

plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)


tiff(file = "E18EGFPFibSMPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo seurat_clusters Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim.ArExp Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo CellTypes Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()


tiff(file = "E18EGFPFibSMPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim.ArExp split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)
dev.off()

####P11 EGFP+ FB/SM w/o ProS > Separate Ctrl & ARKO ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/P11/GFP+/without ProS")

#subclustering EGFP+
Idents(object = P11_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined) <- "RNA"

P11_CtrlvARKO.combined_EGFP <- subset(x=P11_CtrlvARKO.combined, subset = EGFP > 0)

#Add Ar info
DefaultAssay(P11_CtrlvARKO.combined_EGFP) <- "RNA"
P11_CtrlvARKO.combined_EGFP_ArPos <- subset(x=P11_CtrlvARKO.combined_EGFP, subset = Ar > 0)
P11_CtrlvARKO.combined_EGFP_ArNeg <- subset(x=P11_CtrlvARKO.combined_EGFP, subset = Ar == 0)
Idents(object = P11_CtrlvARKO.combined_EGFP_ArPos) <- "ArPos"
Idents(object = P11_CtrlvARKO.combined_EGFP_ArNeg) <- "ArNeg"
P11_CtrlvARKO.combined_EGFP_ArPos[["ArExp"]] <- Idents(object = P11_CtrlvARKO.combined_EGFP_ArPos)
P11_CtrlvARKO.combined_EGFP_ArNeg[["ArExp"]] <- Idents(object = P11_CtrlvARKO.combined_EGFP_ArNeg)
P11_CtrlvARKO.combined_EGFP_Ar <- merge(x = P11_CtrlvARKO.combined_EGFP_ArPos, y = P11_CtrlvARKO.combined_EGFP_ArNeg)
Idents(object = P11_CtrlvARKO.combined_EGFP_Ar) <- "ArExp"
P11_CtrlvARKO.combined_EGFP$ArExp <- Idents(object = P11_CtrlvARKO.combined_EGFP_Ar)

Idents(object = P11_CtrlvARKO.combined_EGFP) <- "ArExp"
P11_CtrlvARKO.combined_EGFP <- subset(P11_CtrlvARKO.combined_EGFP, idents = c("ArPos", "ArNeg"))
P11_CtrlvARKO.combined_EGFP <- RenameIdents(object = P11_CtrlvARKO.combined_EGFP, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(P11_CtrlvARKO.combined_EGFP, reduction = "umap", pt.size = 0.3)

Idents(object = P11_CtrlvARKO.combined_EGFP) <- "stim"
P11_CtrlvARKO.combined_EGFP$stim.ArExp <- paste(Idents(P11_CtrlvARKO.combined_EGFP), P11_CtrlvARKO.combined_EGFP$ArExp, sep = "_")
Idents(object = P11_CtrlvARKO.combined_EGFP) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = P11_CtrlvARKO.combined_EGFP) <- "CellTypes"
P11_CtrlvARKO.combined_EGFP_FibSM <- subset(P11_CtrlvARKO.combined_EGFP, idents = c("Fibroblast", "SM"))
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"
P11_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P11_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(P11_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(P11_CtrlvARKO.combined_EGFP_FibSM, ndims = 50)

P11_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
P11_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(P11_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunTSNE(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Arstim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
P11_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(P11_CtrlvARKO.combined_EGFP_FibSM), P11_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(P11_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P11_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P11_CtrlvARKO.combined_EGFP_FibSM, features = rownames(P11_CtrlvARKO.combined_EGFP_FibSM))
P11_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(P11_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P11_CtrlvARKO.combined_EGFP_FibSMTop10 <- P11_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P11_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P11_CtrlvARKO.combined_EGFP_FibSM, features = c(P11_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(P11_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")

#Separating Ctrl & ARKO
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim"
P11_Ctrl_FBSM <- subset(P11_CtrlvARKO.combined_EGFP_FibSM, idents = c("P11_Ctrl"))
P11_ARKO_FBSM <- subset(P11_CtrlvARKO.combined_EGFP_FibSM, idents = c("P11_ARKO"))

#### eliminate Epi cells ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/P11/GFP+/without ProS")

DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Add Cdh1 info
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P11_CtrlvARKO.combined_EGFP_FibSMCdh1Pos <- subset(x=P11_CtrlvARKO.combined_EGFP_FibSM, subset = Cdh1 > 0)
P11_CtrlvARKO.combined_EGFP_FibSMCdh1Neg <- subset(x=P11_CtrlvARKO.combined_EGFP_FibSM, subset = Cdh1 == 0)
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSMCdh1Pos) <- "Cdh1Pos"
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSMCdh1Neg) <- "Cdh1Neg"
P11_CtrlvARKO.combined_EGFP_FibSMCdh1Pos[["Cdh1Exp"]] <- Idents(object = P11_CtrlvARKO.combined_EGFP_FibSMCdh1Pos)
P11_CtrlvARKO.combined_EGFP_FibSMCdh1Neg[["Cdh1Exp"]] <- Idents(object = P11_CtrlvARKO.combined_EGFP_FibSMCdh1Neg)
P11_CtrlvARKO.combined_EGFP_FibSMCdh1  <- merge(x = P11_CtrlvARKO.combined_EGFP_FibSMCdh1Pos, y = P11_CtrlvARKO.combined_EGFP_FibSMCdh1Neg)
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSMCdh1) <- "Cdh1Exp"
P11_CtrlvARKO.combined_EGFP_FibSM$Cdh1Exp <- Idents(object = P11_CtrlvARKO.combined_EGFP_FibSMCdh1)

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "Cdh1Exp"
P11_CtrlvARKO.combined_EGFP_FibSM <- subset(P11_CtrlvARKO.combined_EGFP_FibSM, idents = c("Cdh1Neg"))
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Run the standard workflow for visualization and clustering
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"

P11_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P11_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(P11_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(P11_CtrlvARKO.combined_EGFP_FibSM)

# umap and Clustering
P11_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:15)
P11_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(P11_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:15)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM_umap.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM_split_umap.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE) 
dev.off()

#Cell count (split)
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
P11_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(P11_CtrlvARKO.combined_EGFP_FibSM), P11_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(P11_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P11_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P11_CtrlvARKO.combined_EGFP_FibSM, features = rownames(P11_CtrlvARKO.combined_EGFP_FibSM))
P11_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(P11_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P11_CtrlvARKO.combined_EGFP_FibSMTop10 <- P11_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P11_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P11_CtrlvARKO.combined_EGFP_FibSM, features = c(P11_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(P11_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")

#Separating Ctrl & ARKO
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim"
P11_Ctrl_FBSM <- subset(P11_CtrlvARKO.combined_EGFP_FibSM, idents = c("P11_Ctrl"))
P11_ARKO_FBSM <- subset(P11_CtrlvARKO.combined_EGFP_FibSM, idents = c("P11_ARKO"))




####P42 EGFP+ FB/SM w/o ProS > Separate Ctrl & ARKO ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/P42/GFP+/without ProS")

#subclustering EGFP+
Idents(object = P42_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined) <- "RNA"

P42_CtrlvARKO.combined_EGFP <- subset(x=P42_CtrlvARKO.combined, subset = EGFP > 0)

#Add Ar info
DefaultAssay(P42_CtrlvARKO.combined_EGFP) <- "RNA"
P42_CtrlvARKO.combined_EGFP_ArPos <- subset(x=P42_CtrlvARKO.combined_EGFP, subset = Ar > 0)
P42_CtrlvARKO.combined_EGFP_ArNeg <- subset(x=P42_CtrlvARKO.combined_EGFP, subset = Ar == 0)
Idents(object = P42_CtrlvARKO.combined_EGFP_ArPos) <- "ArPos"
Idents(object = P42_CtrlvARKO.combined_EGFP_ArNeg) <- "ArNeg"
P42_CtrlvARKO.combined_EGFP_ArPos[["ArExp"]] <- Idents(object = P42_CtrlvARKO.combined_EGFP_ArPos)
P42_CtrlvARKO.combined_EGFP_ArNeg[["ArExp"]] <- Idents(object = P42_CtrlvARKO.combined_EGFP_ArNeg)
P42_CtrlvARKO.combined_EGFP_Ar <- merge(x = P42_CtrlvARKO.combined_EGFP_ArPos, y = P42_CtrlvARKO.combined_EGFP_ArNeg)
Idents(object = P42_CtrlvARKO.combined_EGFP_Ar) <- "ArExp"
P42_CtrlvARKO.combined_EGFP$ArExp <- Idents(object = P42_CtrlvARKO.combined_EGFP_Ar)

Idents(object = P42_CtrlvARKO.combined_EGFP) <- "ArExp"
P42_CtrlvARKO.combined_EGFP <- subset(P42_CtrlvARKO.combined_EGFP, idents = c("ArPos", "ArNeg"))
P42_CtrlvARKO.combined_EGFP <- RenameIdents(object = P42_CtrlvARKO.combined_EGFP, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(P42_CtrlvARKO.combined_EGFP, reduction = "umap", pt.size = 0.3)

Idents(object = P42_CtrlvARKO.combined_EGFP) <- "stim"
P42_CtrlvARKO.combined_EGFP$stim.ArExp <- paste(Idents(P42_CtrlvARKO.combined_EGFP), P42_CtrlvARKO.combined_EGFP$ArExp, sep = "_")
Idents(object = P42_CtrlvARKO.combined_EGFP) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = P42_CtrlvARKO.combined_EGFP) <- "CellTypes"
P42_CtrlvARKO.combined_EGFP_FibSM <- subset(P42_CtrlvARKO.combined_EGFP, idents = c("Fibroblast", "SM"))
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"
P42_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P42_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(P42_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(P42_CtrlvARKO.combined_EGFP_FibSM, ndims = 50)

P42_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
P42_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(P42_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunTSNE(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Arstim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
P42_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(P11_CtrlvARKO.combined_EGFP_FibSM), P11_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(P42_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P42_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P42_CtrlvARKO.combined_EGFP_FibSM, features = rownames(P42_CtrlvARKO.combined_EGFP_FibSM))
P42_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(P42_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P42_CtrlvARKO.combined_EGFP_FibSMTop10 <- P42_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P42_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P42_CtrlvARKO.combined_EGFP_FibSM, features = c(P42_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(P42_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")

#Separating Ctrl & ARKO
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim"
P42_Ctrl_FBSM <- subset(P42_CtrlvARKO.combined_EGFP_FibSM, idents = c("P42_Ctrl"))
P42_ARKO_FBSM <- subset(P42_CtrlvARKO.combined_EGFP_FibSM, idents = c("P42_ARKO"))

#### eliminate Epi cells ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/P42/GFP+/without ProS")

DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Add Cdh1 info
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P42_CtrlvARKO.combined_EGFP_FibSMCdh1Pos <- subset(x=P42_CtrlvARKO.combined_EGFP_FibSM, subset = Cdh1 > 0)
P42_CtrlvARKO.combined_EGFP_FibSMCdh1Neg <- subset(x=P42_CtrlvARKO.combined_EGFP_FibSM, subset = Cdh1 == 0)
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSMCdh1Pos) <- "Cdh1Pos"
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSMCdh1Neg) <- "Cdh1Neg"
P42_CtrlvARKO.combined_EGFP_FibSMCdh1Pos[["Cdh1Exp"]] <- Idents(object = P42_CtrlvARKO.combined_EGFP_FibSMCdh1Pos)
P42_CtrlvARKO.combined_EGFP_FibSMCdh1Neg[["Cdh1Exp"]] <- Idents(object = P42_CtrlvARKO.combined_EGFP_FibSMCdh1Neg)
P42_CtrlvARKO.combined_EGFP_FibSMCdh1  <- merge(x = P42_CtrlvARKO.combined_EGFP_FibSMCdh1Pos, y = P42_CtrlvARKO.combined_EGFP_FibSMCdh1Neg)
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSMCdh1) <- "Cdh1Exp"
P42_CtrlvARKO.combined_EGFP_FibSM$Cdh1Exp <- Idents(object = P42_CtrlvARKO.combined_EGFP_FibSMCdh1)

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "Cdh1Exp"
P42_CtrlvARKO.combined_EGFP_FibSM <- subset(P42_CtrlvARKO.combined_EGFP_FibSM, idents = c("Cdh1Neg"))
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Run the standard workflow for visualization and clustering
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"

P42_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P42_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(P42_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(P42_CtrlvARKO.combined_EGFP_FibSM)

# umap and Clustering
P42_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
P42_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(P42_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM_umap.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM_split_umap.tiff", width = 20, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE) 
dev.off()

#Cell count (split)
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
P42_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(P11_CtrlvARKO.combined_EGFP_FibSM), P11_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(P42_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P42_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P42_CtrlvARKO.combined_EGFP_FibSM, features = rownames(P42_CtrlvARKO.combined_EGFP_FibSM))
P42_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(P42_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P42_CtrlvARKO.combined_EGFP_FibSMTop10 <- P42_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P42_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P42_CtrlvARKO.combined_EGFP_FibSM, features = c(P42_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(P42_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")

#Separating Ctrl & ARKO
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim"
P42_Ctrl_FBSM <- subset(P42_CtrlvARKO.combined_EGFP_FibSM, idents = c("P42_Ctrl"))
P42_ARKO_FBSM <- subset(P42_CtrlvARKO.combined_EGFP_FibSM, idents = c("P42_ARKO"))





####Integrate EGFP+ FB/SM of E18.5 P11 P42####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 FibSM/EGFP+")

#Set Current idents
Idents(object = E18_Ctrl_FBSM) <- "seurat_clusters"
Idents(object = E18_ARKO_FBSM) <- "seurat_clusters"
Idents(object = P11_Ctrl_FBSM) <- "seurat_clusters"
Idents(object = P11_ARKO_FBSM) <- "seurat_clusters"
Idents(object = P42_Ctrl_FBSM) <- "seurat_clusters"
Idents(object = P42_ARKO_FBSM) <- "seurat_clusters"
E18_Ctrl_FBSM$stim <- "E18_Ctrl"
E18_ARKO_FBSM$stim <- "E18_ARKO"
P11_Ctrl_FBSM$stim <- "P11_Ctrl"
P11_ARKO_FBSM$stim <- "P11_ARKO"
P42_Ctrl_FBSM$stim <- "P42_Ctrl"
P42_ARKO_FBSM$stim <- "P42_ARKO"
ARKOvCtrl_FibSM.anchors <- FindIntegrationAnchors(object.list = list(E18_Ctrl_FBSM, E18_ARKO_FBSM, P11_Ctrl_FBSM, P11_ARKO_FBSM, P42_Ctrl_FBSM, P42_ARKO_FBSM), dims = 1:20)
ARKOvCtrl_FibSM<- IntegrateData(anchorset = ARKOvCtrl_FibSM.anchors, dims = 1:20)
DefaultAssay(ARKOvCtrl_FibSM) <- "integrated"

#Run the standard workflow for visualization and clustering
ARKOvCtrl_FibSM <- ScaleData(ARKOvCtrl_FibSM, verbose = FALSE)
ARKOvCtrl_FibSM <- RunPCA(ARKOvCtrl_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_FibSM)

# umap and Clustering
ARKOvCtrl_FibSM <- FindNeighbors(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:20)
ARKOvCtrl_FibSM <- FindClusters(ARKOvCtrl_FibSM, resolution = 0.5)
ARKOvCtrl_FibSM <- RunUMAP(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:20)
ARKOvCtrl_FibSM <- RunTSNE(ARKOvCtrl_FibSM, reduction = "pca", dims = 1:20)

tiff(file = "ARKOvCtrl_FibSM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_FibSM stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

Idents(object = ARKOvCtrl_FibSM) <- "CellTypes"
tiff(file = "ARKOvCtrl_FibSM Celltype UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_FibSM celltype split UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim")dev.off()

#Cell count (split)
Idents(object = ARKOvCtrl_FibSM) <- "stim"
ARKOvCtrl_FibSM$stim.seurat_clusters <- paste(Idents(ARKOvCtrl_FibSM), ARKOvCtrl_FibSM$seurat_clusters, sep = "_")
Idents(object = ARKOvCtrl_FibSM) <- "stim.seurat_clusters"
table(Idents(ARKOvCtrl_FibSM))

#Featureplot
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
tiff(file = "ARKOvCtrl_FibSM split EGFP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM split Ar.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM split Gli1.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "ARKOvCtrl_FibSM EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM Ar.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_FibSM Gli1.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = ARKOvCtrl_FibSM) <- "seurat_clusters"
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
ARKOvCtrl_FibSM <- ScaleData(ARKOvCtrl_FibSM, features = rownames(ARKOvCtrl_FibSM))
ARKOvCtrl_FibSM.markers <- FindAllMarkers(ARKOvCtrl_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ARKOvCtrl_FibSMTop10 <- ARKOvCtrl_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
tiff(file = "ARKOvCtrl_FibSM Heatmap Top10.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(ARKOvCtrl_FibSM, features = c(ARKOvCtrl_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
ARKOvCtrl_FibSM.0.1markers <- FindAllMarkers(ARKOvCtrl_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(ARKOvCtrl_FibSM.0.1markers, "ARKOvCtrl_FibSM.0.1markers.csv")



#### Pseudotime of 6 EGFP_FBSM ####
Idents(object = ARKOvCtrl_FibSM) <- "stim"
DefaultAssay(ARKOvCtrl_FibSM) <- "RNA"
ARKOvCtrl_FibSMPseudo <- as.CellDataSet(ARKOvCtrl_FibSM)
ARKOvCtrl_FibSMPseudo  <- detectGenes(ARKOvCtrl_FibSMPseudo, min_expr = 0.1)
print(head(fData(ARKOvCtrl_FibSMPseudo)))

expressed_genes <- row.names(subset(fData(ARKOvCtrl_FibSMPseudo),
                                    num_cells_expressed >= 10))

pData(ARKOvCtrl_FibSMPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARKOvCtrl_FibSMPseudo))
ARKOvCtrl_FibSMPseudo <- ARKOvCtrl_FibSMPseudo[,pData(ARKOvCtrl_FibSMPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARKOvCtrl_FibSMPseudo), geom =
        "density")

ARKOvCtrl_FibSMPseudo <- estimateSizeFactors(ARKOvCtrl_FibSMPseudo)
ARKOvCtrl_FibSMPseudo <- estimateDispersions(ARKOvCtrl_FibSMPseudo)

disp_table <- dispersionTable(ARKOvCtrl_FibSMPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARKOvCtrl_FibSMPseudo <- setOrderingFilter(ARKOvCtrl_FibSMPseudo, unsup_clustering_genes$gene_id)
#plot_ordering_genes(ARKOvCtrl_FibSMPseudo)

#ARKOvCtrl_FibSMPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(ARKOvCtrl_FibSMPseudo, return_all = F) # norm_method='log'

ARKOvCtrl_FibSMPseudo <- reduceDimension(ARKOvCtrl_FibSMPseudo, max_components = 2, num_dim = 20,
                                      reduction_method = 'tSNE', verbose = T)
ARKOvCtrl_FibSMPseudo <- clusterCells(ARKOvCtrl_FibSMPseudo, num_clusters = 2)

plot_cell_clusters(ARKOvCtrl_FibSMPseudo, color_by = "stim")

diff_test_res <- differentialGeneTest(ARKOvCtrl_FibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~stim")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.005))
ARKOvCtrl_FibSMPseudo <- setOrderingFilter(ARKOvCtrl_FibSMPseudo, ordering_genes)
plot_ordering_genes(ARKOvCtrl_FibSMPseudo)

ARKOvCtrl_FibSMPseudo <- reduceDimension(ARKOvCtrl_FibSMPseudo, max_components = 2,
                                      method = 'DDRTree')

ARKOvCtrl_FibSMPseudo <- orderCells(ARKOvCtrl_FibSMPseudo)

GM_state <- function(ARKOvCtrl_FibSMPseudo){
  if (length(unique(pData(ARKOvCtrl_FibSMPseudo)$State)) > 1){
    T0_counts <- table(pData(ARKOvCtrl_FibSMPseudo)$State, pData(ARKOvCtrl_FibSMPseudo)$stim)[,"E18_ARKO"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARKOvCtrl_FibSMPseudo <- orderCells(ARKOvCtrl_FibSMPseudo, root_state = GM_state(ARKOvCtrl_FibSMPseudo))


tiff(file = "ARKOvCtrl_FibSMPseudo Pseudotime.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo seurat_clusters.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo stim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Celltypes.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()

tiff(file = "ARKOvCtrl_FibSMPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo stim split Pseudotime.tiff", width = 6, height = 9, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 3)
dev.off()

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 FibSM/EGFP+/6 FBSM")

tiff(file = "ARKOvCtrl_FibSMPseudo Ar.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Gli1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Gli1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Rspo3.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Rspo3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Mki67.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Mki67", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Sfrp1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Sfrp1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Ptn.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Ptn", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Igf2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Igf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Ngf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Ngf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Csf1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Csf1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Fgf2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Fgf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Fgf18.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Fgf18", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Ctgf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Ctgf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Bmp2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Bmp2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Bmp4.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Bmp4", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Bmp7.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Bmp7", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Dkk2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Dkk2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Nkd1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Nkd1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()





####Subclustering EGFP+ FB & SM separately of E18.5 P11 P42####

Idents(object = ARKOvCtrl_FibSM) <- "CellTypes"
ARKOvCtrl_Fib <- subset(ARKOvCtrl_FibSM, idents = c("Fibroblast"))
ARKOvCtrl_SM <- subset(ARKOvCtrl_FibSM, idents = c("SM"))

DefaultAssay(ARKOvCtrl_Fib) <- "integrated"

#Run the standard workflow for visualization and clustering
ARKOvCtrl_Fib <- ScaleData(ARKOvCtrl_Fib, verbose = FALSE)
ARKOvCtrl_Fib <- RunPCA(ARKOvCtrl_Fib, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_Fib)

# umap and Clustering
ARKOvCtrl_Fib <- FindNeighbors(ARKOvCtrl_Fib, reduction = "pca", dims = 1:20)
ARKOvCtrl_Fib <- FindClusters(ARKOvCtrl_Fib, resolution = 0.5)
ARKOvCtrl_Fib <- RunUMAP(ARKOvCtrl_Fib, reduction = "pca", dims = 1:20)
ARKOvCtrl_Fib <- RunTSNE(ARKOvCtrl_Fib, reduction = "pca", dims = 1:20)

tiff(file = "ARKOvCtrl_Fib UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Fib, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_Fib stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_Fib, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

DefaultAssay(ARKOvCtrl_SM) <- "integrated"

#Run the standard workflow for visualization and clustering
ARKOvCtrl_SM <- ScaleData(ARKOvCtrl_SM, verbose = FALSE)
ARKOvCtrl_SM <- RunPCA(ARKOvCtrl_SM, npcs = 50, verbose = FALSE)
ElbowPlot(ARKOvCtrl_SM)

# umap and Clustering
ARKOvCtrl_SM <- FindNeighbors(ARKOvCtrl_SM, reduction = "pca", dims = 1:20)
ARKOvCtrl_SM <- FindClusters(ARKOvCtrl_SM, resolution = 0.5)
ARKOvCtrl_SM <- RunUMAP(ARKOvCtrl_SM, reduction = "pca", dims = 1:20)
ARKOvCtrl_SM <- RunTSNE(ARKOvCtrl_SM, reduction = "pca", dims = 1:20)

tiff(file = "ARKOvCtrl_SM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_SM, reduction = "umap", pt.size = 0.3, label = TRUE) 
dev.off()
tiff(file = "ARKOvCtrl_SM stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(ARKOvCtrl_SM, reduction = "umap", pt.size = 0.3, split.by = "stim", label = TRUE) 
dev.off()

#### Pseudotime of 6 EGFP_FB ####
setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 FibSM/EGFP+")

Idents(object = ARKOvCtrl_Fib) <- "stim"
DefaultAssay(ARKOvCtrl_Fib) <- "RNA"
ARKOvCtrl_FibPseudo <- as.CellDataSet(ARKOvCtrl_Fib)
ARKOvCtrl_FibPseudo  <- detectGenes(ARKOvCtrl_FibPseudo, min_expr = 0.1)
print(head(fData(ARKOvCtrl_FibPseudo)))

expressed_genes <- row.names(subset(fData(ARKOvCtrl_FibPseudo),
                                    num_cells_expressed >= 10))

pData(ARKOvCtrl_FibPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARKOvCtrl_FibPseudo))
ARKOvCtrl_FibPseudo <- ARKOvCtrl_FibPseudo[,pData(ARKOvCtrl_FibPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARKOvCtrl_FibPseudo), geom =
        "density")

ARKOvCtrl_FibPseudo <- estimateSizeFactors(ARKOvCtrl_FibPseudo)
ARKOvCtrl_FibPseudo <- estimateDispersions(ARKOvCtrl_FibPseudo)

disp_table <- dispersionTable(ARKOvCtrl_FibPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARKOvCtrl_FibPseudo <- setOrderingFilter(ARKOvCtrl_FibPseudo, unsup_clustering_genes$gene_id)
#plot_ordering_genes(ARKOvCtrl_FibSMPseudo)

#ARKOvCtrl_FibSMPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(ARKOvCtrl_FibSMPseudo, return_all = F) # norm_method='log'

ARKOvCtrl_FibPseudo <- reduceDimension(ARKOvCtrl_FibPseudo, max_components = 2, num_dim = 20,
                                         reduction_method = 'tSNE', verbose = T)
ARKOvCtrl_FibPseudo <- clusterCells(ARKOvCtrl_FibPseudo, num_clusters = 2)

plot_cell_clusters(ARKOvCtrl_FibPseudo, color_by = "stim")

diff_test_res <- differentialGeneTest(ARKOvCtrl_FibPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~stim")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARKOvCtrl_FibPseudo <- setOrderingFilter(ARKOvCtrl_FibPseudo, ordering_genes)
plot_ordering_genes(ARKOvCtrl_FibPseudo)

ARKOvCtrl_FibPseudo <- reduceDimension(ARKOvCtrl_FibPseudo, max_components = 2,
                                         method = 'DDRTree')

ARKOvCtrl_FibPseudo <- orderCells(ARKOvCtrl_FibPseudo)

GM_state <- function(ARKOvCtrl_FibPseudo){
  if (length(unique(pData(ARKOvCtrl_FibPseudo)$State)) > 1){
    T0_counts <- table(pData(ARKOvCtrl_FibPseudo)$State, pData(ARKOvCtrl_FibPseudo)$stim)[,"E18_ARKO"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARKOvCtrl_FibPseudo <- orderCells(ARKOvCtrl_FibPseudo, root_state = GM_state(ARKOvCtrl_FibPseudo))

tiff(file = "ARKOvCtrl_FibPseudo Pseudotime.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo seurat_clusters.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, color_by = "seurat_clusters", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo stim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Celltypes.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()

tiff(file = "ARKOvCtrl_FibPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo stim split Pseudotime.tiff", width = 6, height = 9, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 3)
dev.off()


setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 FibSM/EGFP+/6 FB")

tiff(file = "ARKOvCtrl_FibPseudo Ar.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Gli1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Gli1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Ptn.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ptn", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Fgf2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Fgf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Igf2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Igf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Ctgf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ctgf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Rspo3.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Rspo3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Fgf10.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Fgf10", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibSMPseudo Csf1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibSMPseudo, markers = "Csf1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Ngf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ngf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Fst.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Fst", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Aldh1a3.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Aldh1a3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Ly6e.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ly6e", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Pdgfra.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Pdgfra", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Aldh1a1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Aldh1a1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Pdgfrb.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Pdgfrb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Egfr.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Egfr", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()

tiff(file = "ARKOvCtrl_FibPseudo Bmp2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Bmp2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Bmp4.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Bmp4", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Bmp7.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Bmp7", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Wnt5a.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Wnt5a", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Dkk2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Dkk2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Nkd1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Nkd1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Lum.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Lum", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Ptch1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ptch1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Fzd2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Fzd2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Epha4.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Epha4", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_FibPseudo Fgfr2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Fgfr2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

tiff(file = "ARKOvCtrl_FibPseudo Ly6e.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_FibPseudo, markers = "Ly6e", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()

#Cell count (split)
Idents(object = ARKOvCtrl_Fib) <- "stim"
ARKOvCtrl_Fib$stim.seurat_clusters <- paste(Idents(ARKOvCtrl_Fib), ARKOvCtrl_Fib$seurat_clusters, sep = "_")
Idents(object = ARKOvCtrl_Fib) <- "stim.seurat_clusters"
table(Idents(ARKOvCtrl_Fib))



#### Pseudotime of 6 EGFP_SM ####
Idents(object = ARKOvCtrl_SM) <- "stim"
DefaultAssay(ARKOvCtrl_SM) <- "RNA"
ARKOvCtrl_SMPseudo <- as.CellDataSet(ARKOvCtrl_SM)
ARKOvCtrl_SMPseudo  <- detectGenes(ARKOvCtrl_SMPseudo, min_expr = 0.1)
print(head(fData(ARKOvCtrl_SMPseudo)))

expressed_genes <- row.names(subset(fData(ARKOvCtrl_SMPseudo),
                                    num_cells_expressed >= 10))

pData(ARKOvCtrl_SMPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARKOvCtrl_SMPseudo))
ARKOvCtrl_SMPseudo <- ARKOvCtrl_SMPseudo[,pData(ARKOvCtrl_SMPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARKOvCtrl_SMPseudo), geom =
        "density")

ARKOvCtrl_SMPseudo <- estimateSizeFactors(ARKOvCtrl_SMPseudo)
ARKOvCtrl_SMPseudo <- estimateDispersions(ARKOvCtrl_SMPseudo)

disp_table <- dispersionTable(ARKOvCtrl_SMPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARKOvCtrl_SMPseudo <- setOrderingFilter(ARKOvCtrl_SMPseudo, unsup_clustering_genes$gene_id)
#plot_ordering_genes(ARKOvCtrl_FibSMPseudo)

#ARKOvCtrl_FibSMPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(ARKOvCtrl_FibSMPseudo, return_all = F) # norm_method='log'

ARKOvCtrl_SMPseudo <- reduceDimension(ARKOvCtrl_SMPseudo, max_components = 2, num_dim = 20,
                                         reduction_method = 'tSNE', verbose = T)
ARKOvCtrl_SMPseudo <- clusterCells(ARKOvCtrl_SMPseudo, num_clusters = 2)

plot_cell_clusters(ARKOvCtrl_SMPseudo, color_by = "stim")

diff_test_res <- differentialGeneTest(ARKOvCtrl_SMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~stim")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARKOvCtrl_SMPseudo <- setOrderingFilter(ARKOvCtrl_SMPseudo, ordering_genes)
plot_ordering_genes(ARKOvCtrl_SMPseudo)

ARKOvCtrl_SMPseudo <- reduceDimension(ARKOvCtrl_SMPseudo, max_components = 2,
                                         method = 'DDRTree', ncenter = 50)

ARKOvCtrl_SMPseudo <- orderCells(ARKOvCtrl_SMPseudo)

GM_state <- function(ARKOvCtrl_SMPseudo){
  if (length(unique(pData(ARKOvCtrl_SMPseudo)$State)) > 1){
    T0_counts <- table(pData(ARKOvCtrl_SMPseudo)$State, pData(ARKOvCtrl_SMPseudo)$stim)[,"E18_ARKO"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARKOvCtrl_SMPseudo <- orderCells(ARKOvCtrl_SMPseudo, root_state = GM_state(ARKOvCtrl_SMPseudo))
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

tiff(file = "ARKOvCtrl_SMPseudo Pseudotime.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo seurat_clusters.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo stim.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Celltypes.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()

tiff(file = "ARKOvCtrl_SMPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo stim split Pseudotime.tiff", width = 6, height = 9, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 3)
dev.off()

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/6 FibSM/EGFP+/6 SM")

tiff(file = "ARKOvCtrl_SMPseudo Ar.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Gli1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Gli1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Wnt5a.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Wnt5a", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Wnt7a.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Wnt7a", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off
tiff(file = "ARKOvCtrl_SMPseudo Inhbb.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Inhbb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Tgfb3.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Tgfb3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Tgfb2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Tgfb2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Esr1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Esr1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Pgr.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Pgr", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Pcp4.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Pcp4", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Fgfr2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Fgfr2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off

tiff(file = "ARKOvCtrl_SMPseudo Fst.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Fst", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Ngf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Ngf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Ctgf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Ctgf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()
tiff(file = "ARKOvCtrl_SMPseudo Csf1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARKOvCtrl_SMPseudo, markers = "Csf1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "blue"))
dev.off()

#Cell count (split)
Idents(object = ARKOvCtrl_SM) <- "stim"
ARKOvCtrl_SM$stim.seurat_clusters <- paste(Idents(ARKOvCtrl_SM), ARKOvCtrl_SM$seurat_clusters, sep = "_")
Idents(object = ARKOvCtrl_SM) <- "stim.seurat_clusters"
table(Idents(ARKOvCtrl_SM))






####P42 EGFP+ FB/SM w/o ProS ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/P42/GFP+/without ProS")

#subclustering EGFP+
Idents(object = P42_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined) <- "RNA"

P42_CtrlvARKO.combined_EGFP <- subset(x=P42_CtrlvARKO.combined, subset = EGFP > 0)
P42_CtrlvARKO.combined_EGFP <- subset(x=P42_CtrlvARKO.combined_EGFP, subset = Cdh1 == 0)

#Add Ar info
DefaultAssay(P42_CtrlvARKO.combined_EGFP) <- "RNA"
P42_CtrlvARKO.combined_EGFP_ArPos <- subset(x=P42_CtrlvARKO.combined_EGFP, subset = Ar > 0)
P42_CtrlvARKO.combined_EGFP_ArNeg <- subset(x=P42_CtrlvARKO.combined_EGFP, subset = Ar == 0)
Idents(object = P42_CtrlvARKO.combined_EGFP_ArPos) <- "ArPos"
Idents(object = P42_CtrlvARKO.combined_EGFP_ArNeg) <- "ArNeg"
P42_CtrlvARKO.combined_EGFP_ArPos[["ArExp"]] <- Idents(object = P42_CtrlvARKO.combined_EGFP_ArPos)
P42_CtrlvARKO.combined_EGFP_ArNeg[["ArExp"]] <- Idents(object = P42_CtrlvARKO.combined_EGFP_ArNeg)
P42_CtrlvARKO.combined_EGFP_Ar <- merge(x = P42_CtrlvARKO.combined_EGFP_ArPos, y = P42_CtrlvARKO.combined_EGFP_ArNeg)
Idents(object = P42_CtrlvARKO.combined_EGFP_Ar) <- "ArExp"
P42_CtrlvARKO.combined_EGFP$ArExp <- Idents(object = P42_CtrlvARKO.combined_EGFP_Ar)

Idents(object = P42_CtrlvARKO.combined_EGFP) <- "ArExp"
P42_CtrlvARKO.combined_EGFP <- subset(P42_CtrlvARKO.combined_EGFP, idents = c("ArPos", "ArNeg"))
P42_CtrlvARKO.combined_EGFP <- RenameIdents(object = P42_CtrlvARKO.combined_EGFP, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(P42_CtrlvARKO.combined_EGFP, reduction = "umap", pt.size = 0.3)

Idents(object = P42_CtrlvARKO.combined_EGFP) <- "stim"
P42_CtrlvARKO.combined_EGFP$stim.ArExp <- paste(Idents(P42_CtrlvARKO.combined_EGFP), P42_CtrlvARKO.combined_EGFP$ArExp, sep = "_")
Idents(object = P42_CtrlvARKO.combined_EGFP) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = P42_CtrlvARKO.combined_EGFP) <- "CellTypes"
P42_CtrlvARKO.combined_EGFP_FibSM <- subset(P42_CtrlvARKO.combined_EGFP, idents = c("Fibroblast", "SM"))
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"
P42_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P42_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(P42_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(P42_CtrlvARKO.combined_EGFP_FibSM, ndims = 50)

P42_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:25)
P42_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(P42_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.7)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:25)
P42_CtrlvARKO.combined_EGFP_FibSM <- RunTSNE(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:25)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Arstim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
P42_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(P42_CtrlvARKO.combined_EGFP_FibSM), P42_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(P42_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P42_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = P42_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P42_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P42_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P42_CtrlvARKO.combined_EGFP_FibSM, features = rownames(P42_CtrlvARKO.combined_EGFP_FibSM))
P42_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(P42_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P42_CtrlvARKO.combined_EGFP_FibSMTop10 <- P42_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P42_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P42_CtrlvARKO.combined_EGFP_FibSM, features = c(P42_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(P42_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "P42_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")
