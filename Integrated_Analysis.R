library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
# Multi-core processing
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scRNA_utils.R")

# Load 3 Samples ----
p1 <- Read10X(data.dir = "~/Human CRPC scRNAseq/Raw data/CRPC1/filtered_feature_bc_matrix/")
p2 <- Read10X(data.dir = "~/Human CRPC scRNAseq/Raw data/CRPC2/filtered_feature_bc_matrix/")
p3 <- Read10X(data.dir = "~/Human CRPC scRNAseq/Raw data/CRPC3/filtered_feature_bc_matrix/")
# Create Seurat Objects
s1 <- CreateSeuratObject(counts = p1, project = "CRPC1", min.cells = 3, min.features = 200)
s2 <- CreateSeuratObject(counts = p2, project = "CRPC2", min.cells = 3, min.features = 200)
s3 <- CreateSeuratObject(counts = p3, project = "CRPC3", min.cells = 3, min.features = 200)
s1[["percent.mt"]] <- PercentageFeatureSet(s1, pattern = "^MT-")
s2[["percent.mt"]] <- PercentageFeatureSet(s2, pattern = "^MT-")
s3[["percent.mt"]] <- PercentageFeatureSet(s3, pattern = "^MT-")
# Integration
combined_CRPC_raw <- merge(s1,
    y = list(s2, s3),
    add.cell.ids = c("P1", "P2", "P3"), # To avoid barcode collision
    project = "CRPC"
)

# QC ----
# Check MT gene ratio
combined_CRPC_raw[["percent.mt"]] <- PercentageFeatureSet(combined_CRPC_raw, pattern = "^MT-")
p <- VlnPlot(combined_CRPC_raw,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "orig.ident", ncol = 3, pt.size = 0
)
ggsave("Results/mt_by_patient_VlnPlot.png", plot = p, width = 12, height = 8)

# FeatureScatter plot - separate panels for each patient (P1, P2, P3)
split_obj <- SplitObject(combined_CRPC_raw, split.by = "orig.ident")
scatter_plots <- lapply(names(split_obj), function(sample) {
    p1 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(sample)
    p2 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(sample)
    p1 + p2
})
combined_scatter <- wrap_plots(scatter_plots, ncol = 1)
ggsave("Results/FeatureScatter_by_patient.png", plot = combined_scatter, width = 12, height = 15)

# Filter
combined_CRPC <- subset(combined_CRPC_raw,
    subset = nFeature_RNA > 800 &
        nFeature_RNA < 9000 &
        percent.mt < 20
)

# Normalization(SCTransform) & PCA ----
combined_CRPC <- SCTransform(combined_CRPC, verbose = FALSE)
combined_CRPC <- RunPCA(combined_CRPC, verbose = FALSE)
# Draw elbow plot
ggsave("Results/ElbowPlot.png", plot = ElbowPlot(combined_CRPC, ndims = 50), width = 12, height = 15)

# Integration ----
# CCA Integration
combined_CRPC <- IntegrateLayers(
    object = combined_CRPC,
    method = CCAIntegration,
    normalization.method = "SCT",
    verbose = FALSE
)
# Join layers
combined_CRPC[["RNA"]] <- JoinLayers(combined_CRPC[["RNA"]])

# Clustering(CCA) ----
combined_CRPC <- FindNeighbors(combined_CRPC, reduction = "integrated.dr", dims = 1:30)
combined_CRPC <- FindClusters(combined_CRPC, resolution = 0.5)
combined_CRPC <- RunUMAP(combined_CRPC, reduction = "integrated.dr", dims = 1:30)

# Cell Cycle Identification ----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
# Cell cycle scoring
combined_CRPC <- CellCycleScoring(combined_CRPC, s.features = s.genes, g2m.features = g2m.genes)
# Draw on UMAP & Save
p <- DimPlot(combined_CRPC, group.by = "Phase", pt.size = 0.3)
ggsave(
    "Results/CellCycle_Phase_UMAP.tiff",
    plot = p,
    width = 10, height = 7, dpi = 300
)

# Annotation ----
# Check marker genes
p <- FeaturePlot(combined_CRPC,
    features = c(
        "EPCAM", "KRT8", "AR", "PTPRC",
        "VIM", "PECAM1", "ACTA2", "SYP"
    ),
    ncol = 4, pt.size = 0.1
)
ggsave("Results/FeaturePlot_markers.png", plot = p, width = 20, height = 15)

# Find all markers
utils_save_all_markers(combined_CRPC, "Results/all_markers.csv")

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
    `7` = "CD8⁺ T",
    `8` = "Endothelial",
    `9` = "Epithelial",
    `10` = "Epithelial",
    `11` = "Smooth muscle",
    `12` = "Epithelial",
    `13` = "Epithelial",
    `14` = "Macrophage",
    `15` = "Neuronal",
    `16` = "Mast"
)
# Save annotation to metadata
combined_CRPC$celltype <- Idents(combined_CRPC)

# Save labelled UMAP
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave("Results/Labelled_UMAP_integrated.png", plot = p, width = 15, height = 15)

# Save labelled UMAP by patient
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", split.by = "orig.ident", label = TRUE)
ggsave("Results/Labelled_UMAP_by_patient.png", plot = p, width = 24, height = 15)

# Epithelial Reclustering ----
epi <- subset(combined_CRPC, subset = celltype == "Epithelial")

# Split by patient for re-batch correction
epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)

# Re-run SCTransform on the subset
epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)

# Elbow plot for epithelial subset
ggsave("Results/Epithelial/Epithelial_ElbowPlot.png", plot = ElbowPlot(epi, ndims = 50), width = 12, height = 8)

# Re-integrate across patients
epi <- IntegrateLayers(
    object = epi,
    method = CCAIntegration,
    normalization.method = "SCT",
    verbose = FALSE
)
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])

# Clustering
epi <- FindNeighbors(epi, reduction = "integrated.dr", dims = 1:30)
epi <- FindClusters(epi, resolution = 0.6)
epi <- RunUMAP(epi,
    reduction = "integrated.dr", dims = 1:30,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

# Cell Cycle Identification(Epi) ----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
# Cell cycle scoring
epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes)
# Draw on UMAP & Save
p <- DimPlot(epi, group.by = "Phase", pt.size = 0.3)
ggsave(
    "Results/Epithelial/Epithelial_CellCycle_UMAP.tiff",
    plot = p,
    width = 10, height = 7, dpi = 300
)

# UMAP visualization
p <- DimPlot(epi, reduction = "umap", label = TRUE, pt.size = 0.3)
ggsave("Results/Epithelial/Epithelial_UMAP.png", plot = p, width = 10, height = 8)

# UMAP split by patient
p <- DimPlot(epi, reduction = "umap", split.by = "orig.ident", label = TRUE, pt.size = 0.3)
ggsave("Results/Epithelial/Epithelial_UMAP_by_patient.png", plot = p, width = 24, height = 8)

# Epithelial marker FeaturePlot

# Luminal/Basal, Club markers
p <- FeaturePlot(epi,
    features = c(
        "FOXA1", "HOXB13", "NKX3-1", "KLK3", "TOP2A", "MKI67", # Luminal
        "KRT5", "KRT14", "TP63", # Basal
        "MMP7", "WFDC2" # Club
    ),
    ncol = 4, pt.size = 0.1
)
ggsave("Results/Epithelial/Luminal_Basal_Club_FeaturePlot_markers.png", plot = p, width = 20, height = 16)

# DNPC-related markers
p <- FeaturePlot(epi,
    features = c("CHD7", "MYC", "KMT2C", "KRT7", "SOX2", "SYP", "AR"),
    ncol = 4, pt.size = 0.1
)
ggsave("Results/Epithelial/DNPC_FeaturePlot_markers.png", plot = p, width = 20, height = 16)

# Find all markers for epithelial subclusters
utils_save_all_markers(epi, "Results/Epithelial/Epithelial_all_markers.csv")

# Save RDS ----
saveRDS(combined_CRPC, "Results/combined_CRPC.rds")
saveRDS(epi, "Results/epi.rds")
# Read RDS
combined_CRPC <- readRDS("Results/combined_CRPC.rds")
epi <- readRDS("Results/epi.rds")
