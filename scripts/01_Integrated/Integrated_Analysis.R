library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(SingleR)
library(celldex)
# Multi-core processing
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)
# Fixed worker count + RNGseed makes scDblFinder reproducible across machines
# (BiocParallel allocates one deterministic L'Ecuyer RNG stream per worker), so
# the doublet calls — and the integrated clustering — are stable as long as
# `workers` is not changed. Also used by SingleR below.
bpp <- MulticoreParam(workers = 30, RNGseed = 42)

source("scripts/00_utils/scRNA_utils.R")

# Load-or-compute: when the integrated object already exists, skip the heavy
# QC / scDblFinder / SCT / Harmony / SingleR pipeline and just reload
# it to regenerate figures. Each compute block below is guarded by `computed`;
# the figure/ggsave blocks run in both modes.
INT_RDS  <- "Results/01_Integrated/combined_CRPC.rds"
computed <- !file.exists(INT_RDS)
if (!computed) {
    message("Loading cached ", INT_RDS,
            " — skipping compute, regenerating figures")
    combined_CRPC <- readRDS(INT_RDS)
}

if (computed) {

# Load 3 Samples ----
p1 <- Read10X(data.dir = "Raw_data/CRPC1/filtered_feature_bc_matrix/")
p2 <- Read10X(data.dir = "Raw_data/CRPC2/filtered_feature_bc_matrix/")
p3 <- Read10X(data.dir = "Raw_data/CRPC3/filtered_feature_bc_matrix/")
# Create Seurat Objects
s1 <- CreateSeuratObject(counts = p1, project = "CRPC1", min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = p2, project = "CRPC2", min.cells = 3, min.features = 200)
s3 <- CreateSeuratObject(counts = p3, project = "CRPC3", min.cells = 3, min.features = 200)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")
s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")
s3[["percent.mt"]] <- PercentageFeatureSet(s3, pattern = "^MT-")

# Doublet Removal (per-sample, before merge) ----
run_scDblFinder <- function(seu) {
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce, BPPARAM = bpp)
    seu$scDblFinder.class <- sce$scDblFinder.class
    seu$scDblFinder.score <- sce$scDblFinder.score
    seu <- subset(seu, subset = scDblFinder.class == "singlet")
    return(seu)
}
set.seed(42)
s1 <- run_scDblFinder(s1)
s2 <- run_scDblFinder(s2)
s3 <- run_scDblFinder(s3)

# Integration
combined_CRPC_raw <- merge(s1,
    y = list(s2, s3),
    add.cell.ids = c("P1", "P2", "P3"), # To avoid barcode collision
    project = "CRPC"
)

# QC ----
# Check MT gene ratio
combined_CRPC_raw[["percent.mt"]] <- PercentageFeatureSet(combined_CRPC_raw, pattern = "^MT-")
patient_cols <- utils_cb_palette(dplyr::n_distinct(combined_CRPC_raw$orig.ident))
p <- VlnPlot(combined_CRPC_raw,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "orig.ident", ncol = 3, pt.size = 0, cols = patient_cols
)
ggsave("Results/01_Integrated/QC/mt_by_patient_VlnPlot.png", plot = p, width = 12, height = 8, bg = "white")

# VlnPlot with threshold lines to verify QC cutoffs
qc_thresholds <- list(
    nFeature_RNA = c(200, 9000),
    nCount_RNA   = NULL,
    percent.mt   = c(NA, 20)
)
qc_vln <- lapply(names(qc_thresholds), function(feat) {
    vp <- VlnPlot(combined_CRPC_raw,
        features = feat,
        group.by = "orig.ident", pt.size = 0, cols = patient_cols
    ) + NoLegend()
    th <- qc_thresholds[[feat]]
    if (!is.null(th)) {
        for (y in th) if (!is.na(y)) vp <- vp + geom_hline(yintercept = y, linetype = "dashed", color = "red")
    }
    vp
})
qc_vln_combined <- wrap_plots(qc_vln, ncol = 3)
ggsave("Results/01_Integrated/QC/QC_VlnPlot_with_thresholds.png", plot = qc_vln_combined, width = 14, height = 6, bg = "white")

# FeatureScatter plot - separate panels for each patient (P1, P2, P3)
split_obj <- SplitObject(combined_CRPC_raw, split.by = "orig.ident")
scatter_plots <- lapply(names(split_obj), function(sample) {
    p1 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(sample)
    p2 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(sample)
    p1 + p2
})
combined_scatter <- wrap_plots(scatter_plots, ncol = 1)
ggsave("Results/01_Integrated/QC/FeatureScatter_by_patient.png", plot = combined_scatter, width = 12, height = 15, bg = "white")

# Filter
combined_CRPC <- subset(combined_CRPC_raw,
    subset = nFeature_RNA > 200 &
        nFeature_RNA < 9000 &
        percent.mt < 20
)

# Normalization(SCTransform) & PCA ----
combined_CRPC <- SCTransform(combined_CRPC, verbose = FALSE)
combined_CRPC <- RunPCA(combined_CRPC, verbose = FALSE)
# Draw elbow plot
ggsave("Results/01_Integrated/QC/ElbowPlot.png", plot = ElbowPlot(combined_CRPC, ndims = 50), width = 12, height = 15, bg = "white")

# Integration ----
# Harmony Integration — preserves biological heterogeneity better than CCA
# for tumor data; iterative correction on PCA embeddings.
combined_CRPC <- IntegrateLayers(
    object = combined_CRPC,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    normalization.method = "SCT",
    verbose = FALSE
)
# Join layers
combined_CRPC[["RNA"]] <- JoinLayers(combined_CRPC[["RNA"]])

# Clustering(Harmony) ----
combined_CRPC <- FindNeighbors(combined_CRPC, reduction = "harmony", dims = 1:30)
combined_CRPC <- FindClusters(combined_CRPC, resolution = 0.5)
combined_CRPC <- RunUMAP(combined_CRPC, reduction = "harmony", dims = 1:30)

# Cell Cycle Identification ----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
# Cell cycle scoring
combined_CRPC <- CellCycleScoring(combined_CRPC, s.features = s.genes, g2m.features = g2m.genes)

}  # end if (computed) — compute pipeline

# Draw on UMAP & Save
p <- DimPlot(combined_CRPC, group.by = "Phase", pt.size = 0.3,
    cols = utils_cb_palette(dplyr::n_distinct(combined_CRPC$Phase)))
ggsave("Results/01_Integrated/CellCycle_Phase_UMAP.png", plot = p, width = 10, height = 7, bg = "white")

# Annotation ----
# Check marker genes — 각 패널 제목에 cell-type 마커 정보 함께 표기
marker_genes <- c(
    EPCAM  = "Epithelial",
    KRT8   = "Epithelial (luminal)",
    AR     = "Luminal / AR pathway",
    PTPRC  = "Immune",
    VIM    = "Mesenchymal",
    PECAM1 = "Endothelial",
    ACTA2  = "Smooth muscle",
    SYP    = "Neuroendocrine",
    COL1A1 = "Fibroblast",
    DCN    = "Fibroblast / stromal",
    TPSAB1 = "Mast cell"
)
fp <- FeaturePlot(combined_CRPC,
    features = names(marker_genes),
    ncol = 4, pt.size = 0.1, combine = FALSE
)
# combine = FALSE → features 순서대로 패널 리스트 반환; 제목을 "타입 — 유전자"로 덮어씀
fp <- lapply(seq_along(fp), function(i) {
    fp[[i]] + ggtitle(bquote(.(marker_genes[[i]]) ~ "—" ~ italic(.(names(marker_genes)[i]))))
})
p <- wrap_plots(fp, ncol = 4)
ggsave("Results/01_Integrated/FeaturePlot_markers.png", plot = p, width = 20, height = 19, bg = "white")

# Find all markers
if (computed) {
    utils_save_all_markers(combined_CRPC, "Results/01_Integrated/all_markers.csv")
}

# Save cluster UMAP (annotation 전 단계 — cluster 번호로 확인)
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "seurat_clusters",
    label = TRUE,
    cols = utils_cb_palette(nlevels(factor(combined_CRPC$seurat_clusters))))
ggsave("Results/01_Integrated/Cluster_UMAP_integrated.png", plot = p, width = 15, height = 15, bg = "white")

# Save cluster UMAP by patient
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "seurat_clusters",
    split.by = "orig.ident", label = TRUE,
    cols = utils_cb_palette(nlevels(factor(combined_CRPC$seurat_clusters))))
ggsave("Results/01_Integrated/Cluster_UMAP_by_patient.png", plot = p, width = 24, height = 15, bg = "white")

if (computed) {

# SingleR Annotation ----
# Reference: Human Primary Cell Atlas (bulk RNA-seq, broad cell types)
ref_hpca <- celldex::HumanPrimaryCellAtlasData()

# Seurat -> SCE (use SCT-normalized data for SingleR input)
sce_CRPC <- as.SingleCellExperiment(combined_CRPC, assay = "SCT")

# Per-cell annotation (SingleR vignette default — safer, preserves ambiguity)
singler_pred <- SingleR(
    test = sce_CRPC,
    ref = ref_hpca,
    labels = ref_hpca$label.main,
    BPPARAM = bpp
)

# Assign labels to Seurat metadata (row order matches cell order)
combined_CRPC$SingleR <- singler_pred$labels
# pruned.labels = labels with low-confidence calls set to NA
combined_CRPC$SingleR_pruned <- singler_pred$pruned.labels

# Save annotation table
write.csv(
    data.frame(
        cell = rownames(singler_pred),
        label = singler_pred$labels,
        pruned_label = singler_pred$pruned.labels
    ),
    "Results/01_Integrated/SingleR_cell_labels.csv",
    row.names = FALSE
)

# Diagnostic plots — score heatmap & delta distribution
png("Results/01_Integrated/SingleR_score_heatmap.png", width = 1200, height = 900, bg = "white")
plotScoreHeatmap(singler_pred)
dev.off()

png("Results/01_Integrated/SingleR_delta_distribution.png",
    width = 1200, height = 900, bg = "white"
)
plotDeltaDistribution(singler_pred, ncol = 4)
dev.off()

}  # end if (computed) — SingleR

# UMAP with SingleR labels (pruned — uncertain cells appear as NA/grey)
p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "SingleR_pruned",
    label = TRUE, repel = TRUE,
    cols = utils_cb_palette(length(unique(na.omit(combined_CRPC$SingleR_pruned))))
)
ggsave("Results/01_Integrated/SingleR_UMAP.png",
    plot = p, width = 15, height = 15, bg = "white"
)

# UMAP with SingleR labels by patient
p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "SingleR",
    split.by = "orig.ident", label = TRUE, repel = TRUE,
    cols = utils_cb_palette(length(unique(na.omit(combined_CRPC$SingleR))))
)
ggsave("Results/01_Integrated/SingleR_UMAP_by_patient.png",
    plot = p, width = 24, height = 15, bg = "white"
)

if (computed) {

# Label Annotation
combined_CRPC <- RenameIdents(
    object = combined_CRPC,
    `0` = "Epithelial",
    `1` = "Fibroblast",
    `2` = "Epithelial",
    `3` = "Epithelial",
    `4` = "Epithelial",
    `5` = "Epithelial",
    `6` = "Epithelial",
    `7` = "Endothelial",
    `8` = "T/NK cells",
    `9` = "Epithelial",
    `10` = "Smooth muscle cells",
    `11` = "Epithelial",
    `12` = "Epithelial",
    `13` = "Epithelial",
    `14` = "Stromal",
    `15` = "Phagocytes",
    `16` = "Immune cells",
    `17` = "Mast cells"
)
# Save annotation to metadata
combined_CRPC$celltype <- Idents(combined_CRPC)

}  # end if (computed) — celltype annotation

# Save labelled UMAP
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", label = TRUE,
    cols = utils_cb_palette(nlevels(factor(combined_CRPC$celltype))))
ggsave("Results/01_Integrated/Labelled_UMAP_integrated.png", plot = p, width = 15, height = 15, bg = "white")

# Save labelled UMAP by patient
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", split.by = "orig.ident", label = TRUE,
    cols = utils_cb_palette(nlevels(factor(combined_CRPC$celltype))))
ggsave("Results/01_Integrated/Labelled_UMAP_by_patient.png", plot = p, width = 24, height = 15, bg = "white")

# Save RDS ----
if (computed) {
    saveRDS(combined_CRPC, INT_RDS)
}
# Read RDS
# combined_CRPC <- readRDS("Results/01_Integrated/combined_CRPC.rds")
