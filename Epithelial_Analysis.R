library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
# Multi-core processing
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scRNA_utils.R")

# Load data & Epithelial subset ----
combined_CRPC <- readRDS("Results/Integrated/combined_CRPC.rds")
epi <- subset(combined_CRPC, subset = celltype == "Epithelial")
rm(combined_CRPC)

# Remove inherited SCT assay from parent integration to avoid conflicts
DefaultAssay(epi) <- "RNA"
epi[["SCT"]] <- NULL

# Split by patient for batch correction
epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)

# Shared SCTransform + PCA ----
epi_sct <- SCTransform(epi, verbose = FALSE)
epi_sct <- RunPCA(epi_sct, verbose = FALSE)
ggsave("Results/Epithelial/SCT_ElbowPlot.png",
    plot = ElbowPlot(epi_sct, ndims = 50), width = 12, height = 8
)

# ============================================================
# Method 1: SCTransform + CCA ----
# ============================================================
epi_sct_cca <- IntegrateLayers(
    object = epi_sct,
    method = CCAIntegration,
    normalization.method = "SCT",
    verbose = FALSE
)
epi_sct_cca[["RNA"]] <- JoinLayers(epi_sct_cca[["RNA"]])

epi_sct_cca <- FindNeighbors(epi_sct_cca, reduction = "integrated.dr", dims = 1:30)
epi_sct_cca <- FindClusters(epi_sct_cca, resolution = 0.6)
epi_sct_cca <- RunUMAP(epi_sct_cca,
    reduction = "integrated.dr", dims = 1:30,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

p1 <- DimPlot(epi_sct_cca, label = TRUE, pt.size = 0.3) +
    ggtitle("SCT + CCA")
p2 <- DimPlot(epi_sct_cca, group.by = "orig.ident", pt.size = 0.3) +
    ggtitle("SCT + CCA (by patient)")
ggsave("Results/Epithelial/SCT_CCA_UMAP_cluster.png", plot = p1, width = 10, height = 8)
ggsave("Results/Epithelial/SCT_CCA_UMAP_patient.png", plot = p2, width = 10, height = 8)

# ============================================================
# Method 2: LogNormalize + Harmony ----
# ============================================================
epi_log_harmony <- NormalizeData(epi, verbose = FALSE)
epi_log_harmony <- FindVariableFeatures(epi_log_harmony, verbose = FALSE)
epi_log_harmony <- ScaleData(epi_log_harmony, verbose = FALSE)
epi_log_harmony <- RunPCA(epi_log_harmony, verbose = FALSE)

epi_log_harmony <- IntegrateLayers(
    object = epi_log_harmony,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    verbose = FALSE
)
epi_log_harmony[["RNA"]] <- JoinLayers(epi_log_harmony[["RNA"]])

epi_log_harmony <- FindNeighbors(epi_log_harmony, reduction = "harmony", dims = 1:30)
epi_log_harmony <- FindClusters(epi_log_harmony, resolution = 0.6)
epi_log_harmony <- RunUMAP(epi_log_harmony,
    reduction = "harmony", dims = 1:30,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

p1 <- DimPlot(epi_log_harmony, label = TRUE, pt.size = 0.3) +
    ggtitle("LogNorm + Harmony")
p2 <- DimPlot(epi_log_harmony, group.by = "orig.ident", pt.size = 0.3) +
    ggtitle("LogNorm + Harmony (by patient)")
ggsave("Results/Epithelial/LogNorm_Harmony_UMAP_cluster.png", plot = p1, width = 10, height = 8)
ggsave("Results/Epithelial/LogNorm_Harmony_UMAP_patient.png", plot = p2, width = 10, height = 8)

# ============================================================
# Method 3: LogNormalize + CCA ----
# ============================================================
epi_log_cca <- NormalizeData(epi, verbose = FALSE)
epi_log_cca <- FindVariableFeatures(epi_log_cca, verbose = FALSE)
epi_log_cca <- ScaleData(epi_log_cca, verbose = FALSE)
epi_log_cca <- RunPCA(epi_log_cca, verbose = FALSE)
ggsave("Results/Epithelial/LogNorm_CCA_ElbowPlot.png",
    plot = ElbowPlot(epi_log_cca, ndims = 50), width = 12, height = 8
)

epi_log_cca <- IntegrateLayers(
    object = epi_log_cca,
    method = CCAIntegration,
    verbose = FALSE
)
epi_log_cca[["RNA"]] <- JoinLayers(epi_log_cca[["RNA"]])

epi_log_cca <- FindNeighbors(epi_log_cca, reduction = "integrated.dr", dims = 1:30)
epi_log_cca <- FindClusters(epi_log_cca, resolution = 0.6)
epi_log_cca <- RunUMAP(epi_log_cca,
    reduction = "integrated.dr", dims = 1:30,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

p1 <- DimPlot(epi_log_cca, label = TRUE, pt.size = 0.3) +
    ggtitle("LogNorm + CCA")
p2 <- DimPlot(epi_log_cca, group.by = "orig.ident", pt.size = 0.3) +
    ggtitle("LogNorm + CCA (by patient)")
ggsave("Results/Epithelial/LogNorm_CCA_UMAP_cluster.png", plot = p1, width = 10, height = 8)
ggsave("Results/Epithelial/LogNorm_CCA_UMAP_patient.png", plot = p2, width = 10, height = 8)

# ============================================================
# Comparison Plot ----
# ============================================================
p_compare <- (
    DimPlot(epi_sct_cca, label = TRUE, pt.size = 0.1) +
        ggtitle("SCT + CCA") + NoLegend() |
        DimPlot(epi_log_harmony, label = TRUE, pt.size = 0.1) +
            ggtitle("LogNorm + Harmony") + NoLegend() |
        DimPlot(epi_log_cca, label = TRUE, pt.size = 0.1) +
            ggtitle("LogNorm + CCA") + NoLegend()
)
ggsave("Results/Epithelial/Comparison_UMAP.png",
    plot = p_compare, width = 24, height = 8
)

p_compare_patient <- (
    DimPlot(epi_sct_cca, group.by = "orig.ident", pt.size = 0.1) +
        ggtitle("SCT + CCA") |
        DimPlot(epi_log_harmony, group.by = "orig.ident", pt.size = 0.1) +
            ggtitle("LogNorm + Harmony") |
        DimPlot(epi_log_cca, group.by = "orig.ident", pt.size = 0.1) +
            ggtitle("LogNorm + CCA")
)
ggsave("Results/Epithelial/Comparison_UMAP_by_patient.png",
    plot = p_compare_patient, width = 24, height = 8
)

# ============================================================
# Downstream Analysis (all methods) ----
# ============================================================
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes

methods <- list(
    list(obj = epi_sct_cca, tag = "SCT_CCA"),
    list(obj = epi_log_harmony, tag = "LogNorm_Harmony"),
    list(obj = epi_log_cca, tag = "LogNorm_CCA")
)

for (m in methods) {
    obj <- m$obj
    tag <- m$tag

    # Cell Cycle Identification
    obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, nbin = 12)
    p <- DimPlot(obj, group.by = "Phase", pt.size = 0.3)
    ggsave(paste0("Results/Epithelial/", tag, "_CellCycle_UMAP.tiff"),
        plot = p, width = 10, height = 7, dpi = 300
    )

    # UMAP split by patient
    p <- DimPlot(obj, split.by = "orig.ident", label = TRUE, pt.size = 0.3)
    ggsave(paste0("Results/Epithelial/", tag, "_UMAP_by_patient.png"),
        plot = p, width = 24, height = 8
    )

    # Marker gene sets
    luminal_basal_club <- c(
        "FOXA1", "HOXB13", "NKX3-1", "KLK3", "TOP2A", "MKI67", # Luminal
        "KRT5", "KRT14", "TP63", # Basal
        "MMP7", "WFDC2" # Club
    )
    dnpc_markers <- c("CHD7", "MYC", "KMT2C", "KRT7", "SOX2", "SYP", "AR")

    # Luminal/Basal, Club — FeaturePlot
    p <- FeaturePlot(obj, features = luminal_basal_club, ncol = 4, pt.size = 0.1)
    ggsave(paste0("Results/Epithelial/", tag, "_Luminal_Basal_Club_FeaturePlot.png"),
        plot = p, width = 20, height = 16
    )

    # Luminal/Basal, Club — DotPlot
    p <- DotPlot(obj, features = luminal_basal_club) + RotatedAxis()
    ggsave(paste0("Results/Epithelial/", tag, "_Luminal_Basal_Club_DotPlot.png"),
        plot = p, width = 14, height = 8
    )

    # DNPC — FeaturePlot
    p <- FeaturePlot(obj, features = dnpc_markers, ncol = 4, pt.size = 0.1)
    ggsave(paste0("Results/Epithelial/", tag, "_DNPC_FeaturePlot.png"),
        plot = p, width = 20, height = 16
    )

    # DNPC — DotPlot
    p <- DotPlot(obj, features = dnpc_markers) + RotatedAxis()
    ggsave(paste0("Results/Epithelial/", tag, "_DNPC_DotPlot.png"),
        plot = p, width = 12, height = 8
    )

    # Find all markers
    utils_save_all_markers(obj, paste0("Results/Epithelial/", tag, "_all_markers.csv"))
}

# Save RDS ----
saveRDS(epi_sct_cca, "Results/Epithelial/epi_SCT_CCA.rds")
saveRDS(epi_log_harmony, "Results/Epithelial/epi_LogNorm_Harmony.rds")
saveRDS(epi_log_cca, "Results/Epithelial/epi_LogNorm_CCA.rds")

# Read RDS
# epi_sct_cca <- readRDS("Results/Epithelial/epi_SCT_CCA.rds")
# epi_log_harmony <- readRDS("Results/Epithelial/epi_LogNorm_Harmony.rds")
# epi_log_cca <- readRDS("Results/Epithelial/epi_LogNorm_CCA.rds")
