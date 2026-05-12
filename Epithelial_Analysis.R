library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(clustree)

# Multi-core processing
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scRNA_utils.R")

dir.create("Results/Epithelial", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load: epithelial subset ----
# ============================================================
combined_CRPC <- readRDS("Results/Integrated/combined_CRPC.rds")

epi <- subset(combined_CRPC, subset = celltype == "Epithelial")

rm(combined_CRPC)
gc()

# Reset assay state (drop inherited SCT, re-split RNA by patient for Harmony)
DefaultAssay(epi) <- "RNA"
if ("SCT" %in% names(epi@assays)) epi[["SCT"]] <- NULL
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)

# ============================================================
# SCTransform + PCA + Clustering ----
# ============================================================
epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)

# Harmony integration on epithelial subset — corrects patient-level batch
# effects on the PCA embedding using the split RNA layer structure.
epi <- IntegrateLayers(
    object = epi,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    normalization.method = "SCT",
    verbose = FALSE
)
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])

ggsave("Results/Epithelial/ElbowPlot.png",
    plot = ElbowPlot(epi, ndims = 50), width = 12, height = 8
)

epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:50)

# Multi-resolution clustering for clustree stability assessment.
# All resolutions stored as SCT_snn_res.X columns; final ident set to 0.6.
res_vec <- c(0.1, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2)
epi <- FindClusters(epi, resolution = res_vec)
Idents(epi) <- "SCT_snn_res.0.6"
epi$seurat_clusters <- epi$SCT_snn_res.0.6

epi <- RunUMAP(epi,
    reduction = "harmony", dims = 1:50,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

# clustree: cluster stability tree across resolutions
# - node size  = cell count in that cluster
# - edge width = # cells flowing between clusters across resolutions
# - in_prop    = fraction of cells entering a cluster from the dominant parent
#                (low in_prop edges = unstable cluster pulling from multiple parents)
p_tree <- clustree(epi, prefix = "SCT_snn_res.")
ggsave("Results/Epithelial/clustree.png",
    plot = p_tree, width = 12, height = 14, dpi = 200
)

# UMAPs at each resolution for side-by-side comparison
umap_list <- lapply(res_vec, function(r) {
    DimPlot(epi, group.by = paste0("SCT_snn_res.", r), label = TRUE, pt.size = 0.2) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
p_umap_grid <- wrap_plots(umap_list, ncol = 3)
ggsave("Results/Epithelial/UMAP_multires.png",
    plot = p_umap_grid, width = 18, height = 14, dpi = 200
)

p1 <- DimPlot(epi, label = TRUE, pt.size = 0.3) + ggtitle("Epithelial clusters")
p2 <- DimPlot(epi, group.by = "orig.ident", pt.size = 0.3) + ggtitle("By patient")
ggsave("Results/Epithelial/UMAP_cluster.png", plot = p1, width = 10, height = 8)
ggsave("Results/Epithelial/UMAP_patient.png", plot = p2, width = 10, height = 8)

# Patient-split UMAP — patient-specific clusters are tumor candidates
p <- DimPlot(epi, split.by = "orig.ident", label = TRUE, pt.size = 0.3)
ggsave("Results/Epithelial/UMAP_split_by_patient.png",
    plot = p, width = 24, height = 8
)

# ============================================================
# copyKAT malignant overlay ----
# ============================================================
# Overlay copyKAT predictions (aneuploid/diploid/not.defined) inherited
# from the parent integrated object onto the reclustered epithelial UMAP.
# Clusters with high aneuploid fraction are malignant subclusters.
copykat_cols <- c(
    "aneuploid" = "#D7261E", "diploid" = "#1F77B4", "not.defined" = "grey80"
)

p_pred <- DimPlot(epi,
    group.by = "copykat_prediction", pt.size = 0.3, cols = copykat_cols
) + ggtitle("copyKAT prediction")
ggsave("Results/Epithelial/UMAP_copykat_prediction.png",
    plot = p_pred, width = 10, height = 8
)

p_pred_split <- DimPlot(epi,
    group.by = "copykat_prediction", split.by = "orig.ident",
    pt.size = 0.3, cols = copykat_cols
)
ggsave("Results/Epithelial/UMAP_copykat_prediction_by_patient.png",
    plot = p_pred_split, width = 24, height = 8
)

# Per-cluster composition table — sort by aneuploid fraction to surface
# malignant-enriched clusters.
copykat_cluster_table <- epi@meta.data %>%
    dplyr::count(seurat_clusters, copykat_prediction, name = "n") %>%
    tidyr::pivot_wider(
        names_from = copykat_prediction, values_from = n, values_fill = 0
    )
num_cols <- setdiff(names(copykat_cluster_table), "seurat_clusters")
copykat_cluster_table$total <-
    rowSums(copykat_cluster_table[, num_cols, drop = FALSE])
if ("aneuploid" %in% num_cols) {
    copykat_cluster_table$aneuploid_frac <-
        copykat_cluster_table$aneuploid / copykat_cluster_table$total
    copykat_cluster_table <- copykat_cluster_table[
        order(-copykat_cluster_table$aneuploid_frac),
    ]
}
write.csv(copykat_cluster_table,
    "Results/Epithelial/copykat_cluster_composition.csv", row.names = FALSE
)

# Stacked barplot of aneuploid/diploid fractions per cluster
p_bar <- ggplot(epi@meta.data,
    aes(x = seurat_clusters, fill = copykat_prediction)
) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = copykat_cols) +
    labs(x = "Cluster", y = "Proportion", fill = "copyKAT") +
    theme_classic()
ggsave("Results/Epithelial/copykat_cluster_barplot.png",
    plot = p_bar, width = 10, height = 6
)

# ============================================================
# Cell cycle scoring ----
# ============================================================
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, nbin = 12)

p <- DimPlot(epi, group.by = "Phase", pt.size = 0.3)
ggsave("Results/Epithelial/CellCycle_UMAP.tiff",
    plot = p, width = 10, height = 7, dpi = 300
)

# ============================================================
# Marker visualization ----
# ============================================================
luminal_basal_club <- c(
    "FOXA1", "HOXB13", "NKX3-1", "KLK3", "TOP2A", "MKI67", # Luminal
    "KRT5", "KRT14", "TP63", # Basal
    "MMP7", "WFDC2" # Club
)
dnpc_markers <- c("CHD7", "MYC", "KMT2C", "KRT7", "SOX2", "SYP", "AR")

p <- FeaturePlot(epi, features = luminal_basal_club, ncol = 4, pt.size = 0.1)
ggsave("Results/Epithelial/Luminal_Basal_Club_FeaturePlot.png",
    plot = p, width = 20, height = 16
)

p <- DotPlot(epi, features = luminal_basal_club) + RotatedAxis()
ggsave("Results/Epithelial/Luminal_Basal_Club_DotPlot.png",
    plot = p, width = 14, height = 8
)

p <- FeaturePlot(epi, features = dnpc_markers, ncol = 4, pt.size = 0.1)
ggsave("Results/Epithelial/DNPC_FeaturePlot.png",
    plot = p, width = 20, height = 16
)

p <- DotPlot(epi, features = dnpc_markers) + RotatedAxis()
ggsave("Results/Epithelial/DNPC_DotPlot.png",
    plot = p, width = 12, height = 8
)

# ============================================================
# Find all markers ----
# ============================================================
utils_save_all_markers(epi, "Results/Epithelial/all_markers.csv")

saveRDS(epi, "Results/Epithelial/epi_clustered.rds")