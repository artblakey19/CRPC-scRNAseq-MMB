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

source("scripts/00_utils/scRNA_utils.R")

dir.create("Results/02_Epithelial_Initial", showWarnings = FALSE, recursive = TRUE)

# Load-or-compute: when the clustered object already exists, skip the heavy
# SCT/Harmony/clustering pipeline and just reload it to regenerate figures.
EPI_RDS  <- "Results/02_Epithelial_Initial/epi_clustered.rds"
computed <- !file.exists(EPI_RDS)

# Resolutions swept for clustree (also used by the figure section below).
res_vec <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)

if (computed) {

# ============================================================
# Load: epithelial subset ----
# ============================================================
combined_CRPC <- readRDS("Results/01_Integrated/combined_CRPC.rds")

epi <- subset(combined_CRPC, subset = celltype == "Epithelial")

rm(combined_CRPC)
gc()

# Reset assay state (drop inherited SCT, re-split RNA by patient for Harmony)
DefaultAssay(epi) <- "RNA"
if ("SCT" %in% names(epi@assays)) epi[["SCT"]] <- NULL
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)

# Drop inherited *_snn_res.* / seurat_clusters columns from the parent object
# so clustree() doesn't pick up stale resolutions (e.g. 0.5 from Integrated).
epi@meta.data <- dplyr::select(epi@meta.data,
    -matches("_snn_res\\."), -any_of("seurat_clusters"))

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

ggsave("Results/02_Epithelial_Initial/ElbowPlot.png",
    plot = ElbowPlot(epi, ndims = 50), width = 12, height = 8, bg = "white"
)

epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)

# Multi-resolution clustering for clustree stability assessment.
# All resolutions stored as SCT_snn_res.X columns; final ident set to 0.6.
epi <- FindClusters(epi, resolution = res_vec)

# Curate small clusters at res 0.4 based on clustree (0.2 parent inheritance):
#   12 + 13 -> 12  (both are fragments of 0.2 cluster 10; ~80 cells total,
#                   mostly diploid/not.defined, treated as low-quality residue)
#   14 -> 10       (14 is small split-off of 0.2 cluster 8; 10 is its main body)
#   15 ->  9       (15 is small split-off of 0.2 cluster 7;  9 is its main body)
#   11 kept       (stable across 0.2 and 0.4, ~206 cells, ~62% aneuploid mixed)
merge_map <- c("13" = "12", "14" = "10", "15" = "9")
curated <- dplyr::recode(as.character(epi$SCT_snn_res.0.4), !!!merge_map)
epi$seurat_clusters <- factor(curated,
    levels = as.character(sort(unique(as.integer(curated)))))
Idents(epi) <- "seurat_clusters"

epi <- RunUMAP(epi,
    reduction = "harmony", dims = 1:30,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

} else {
    message("Loading cached ", EPI_RDS,
            " — skipping compute, regenerating figures")
    epi <- readRDS(EPI_RDS)
    Idents(epi) <- "seurat_clusters"
}

# clustree: cluster stability tree across resolutions
# - node size  = cell count in that cluster
# - edge width = # cells flowing between clusters across resolutions
# - in_prop    = fraction of cells entering a cluster from the dominant parent
#                (low in_prop edges = unstable cluster pulling from multiple parents)
p_tree <- clustree(epi, prefix = "SCT_snn_res.")
ggsave("Results/02_Epithelial_Initial/clustree.png",
    plot = p_tree, width = 12, height = 14, dpi = 200, bg = "white"
)

# UMAPs at each resolution for side-by-side comparison
umap_list <- lapply(res_vec, function(r) {
    rcol <- paste0("SCT_snn_res.", r)
    DimPlot(epi, group.by = rcol, label = TRUE, pt.size = 0.2,
            cols = utils_cb_palette(dplyr::n_distinct(epi@meta.data[[rcol]]))) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
p_umap_grid <- wrap_plots(umap_list, ncol = 3)
ggsave("Results/02_Epithelial_Initial/UMAP_multires.png",
    plot = p_umap_grid, width = 18, height = 14, dpi = 200, bg = "white"
)

p1 <- DimPlot(epi, label = TRUE, pt.size = 0.3,
              cols = utils_cb_palette(nlevels(factor(epi$seurat_clusters)))) +
    ggtitle("Epithelial clusters")
p2 <- DimPlot(epi, group.by = "orig.ident", pt.size = 0.3,
              cols = utils_cb_palette(dplyr::n_distinct(epi$orig.ident))) +
    ggtitle("By patient")
ggsave("Results/02_Epithelial_Initial/UMAP_cluster.png", plot = p1, width = 10, height = 8, bg = "white")
ggsave("Results/02_Epithelial_Initial/UMAP_patient.png", plot = p2, width = 10, height = 8, bg = "white")

# Patient-split UMAP — patient-specific clusters are tumor candidates
p <- DimPlot(epi, split.by = "orig.ident", label = TRUE, pt.size = 0.3,
             cols = utils_cb_palette(nlevels(factor(epi$seurat_clusters))))
ggsave("Results/02_Epithelial_Initial/UMAP_split_by_patient.png",
    plot = p, width = 24, height = 8, bg = "white"
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
ggsave("Results/02_Epithelial_Initial/UMAP_copykat_prediction.png",
    plot = p_pred, width = 10, height = 8, bg = "white"
)

p_pred_split <- DimPlot(epi,
    group.by = "copykat_prediction", split.by = "orig.ident",
    pt.size = 0.3, cols = copykat_cols
)
ggsave("Results/02_Epithelial_Initial/UMAP_copykat_prediction_by_patient.png",
    plot = p_pred_split, width = 24, height = 8, bg = "white"
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
    "Results/02_Epithelial_Initial/copykat_cluster_composition.csv", row.names = FALSE
)

# Stacked barplot of aneuploid/diploid fractions per cluster
p_bar <- ggplot(epi@meta.data,
    aes(x = seurat_clusters, fill = copykat_prediction)
) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = copykat_cols) +
    labs(x = "Cluster", y = "Proportion", fill = "copyKAT") +
    theme_classic()
ggsave("Results/02_Epithelial_Initial/copykat_cluster_barplot.png",
    plot = p_bar, width = 10, height = 6, bg = "white"
)

# ============================================================
# Cell cycle scoring ----
# ============================================================
if (computed) {
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, nbin = 12)
}

p <- DimPlot(epi, group.by = "Phase", pt.size = 0.3,
             cols = utils_cb_palette(dplyr::n_distinct(epi$Phase)))
ggsave("Results/02_Epithelial_Initial/CellCycle_UMAP.png", plot = p, width = 10, height = 7, bg = "white")

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
ggsave("Results/02_Epithelial_Initial/Luminal_Basal_Club_FeaturePlot.png",
    plot = p, width = 20, height = 16, bg = "white"
)

p <- DotPlot(epi, features = luminal_basal_club) + RotatedAxis()
ggsave("Results/02_Epithelial_Initial/Luminal_Basal_Club_DotPlot.png",
    plot = p, width = 14, height = 8, bg = "white"
)

p <- FeaturePlot(epi, features = dnpc_markers, ncol = 4, pt.size = 0.1)
ggsave("Results/02_Epithelial_Initial/DNPC_FeaturePlot.png",
    plot = p, width = 20, height = 16, bg = "white"
)

p <- DotPlot(epi, features = dnpc_markers) + RotatedAxis()
ggsave("Results/02_Epithelial_Initial/DNPC_DotPlot.png",
    plot = p, width = 12, height = 8, bg = "white"
)

# CRPC subtype panel: ARPC / DNPC-FGF / DNPC-KRT7 / NEPC
# Genes grouped by subtype axis (column order = subtype order in the table).
subtype_markers <- list(
    ARPC       = c("AR", "KLK3", "FOLH1"),
    `DNPC-FGF` = c("FGF8", "FGFR1", "CCL2", "DKK1"),
    `DNPC-KRT7` = c("KRT7", "SOX2", "FOXA2"),
    NEPC       = c("CHGA", "SYP", "PROX1")
)

p <- DotPlot(epi, features = subtype_markers, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "CRPC subtype markers (ARPC / DNPC-FGF / DNPC-KRT7 / NEPC)")
ggsave("Results/02_Epithelial_Initial/Subtype_DotPlot.png",
    plot = p, width = 14, height = 8, bg = "white"
)

# Hillock / club-hillock verification panel for cluster 5
# True hillock = KRT13+ KRT4+ KRT5+ TP63+(weak) AQP3+ LY6D+
# Club-hillock intermediate = KRT13+ AQP3+ LY6D+ PSCA+ LYPD3+ but KRT5/TP63 low
hillock_markers <- list(
    Hillock_core    = c("KRT13", "KRT4", "AQP3", "LY6D", "PSCA", "LYPD3"),
    Basal_axis      = c("KRT5", "KRT14", "TP63"),
    Squamous_diff   = c("SERPINB5", "SERPINB1", "ZNF750", "S100A2"),
    Club            = c("MMP7", "WFDC2", "PIGR", "LCN2"),
    Luminal_ref     = c("KRT8", "KRT18", "AR", "KLK3")
)

p <- DotPlot(epi, features = hillock_markers, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "Hillock / club-hillock verification (cluster 5)")
ggsave("Results/02_Epithelial_Initial/Hillock_DotPlot.png",
    plot = p, width = 16, height = 8, bg = "white"
)

# ============================================================
# Find all markers + save (compute mode only) ----
# ============================================================
if (computed) {
    utils_save_all_markers(epi, "Results/02_Epithelial_Initial/all_markers.csv")
    saveRDS(epi, EPI_RDS)
}