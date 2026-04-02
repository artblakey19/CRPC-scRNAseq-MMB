library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
# Multi-core processing
library(future)
plan("multisession", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scRNA_utils.R")

# Load & QC ----
# Load 3 Samples
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

# Check MT gene ratio & visualization
combined_CRPC_raw[["percent.mt"]] <- PercentageFeatureSet(combined_CRPC_raw, pattern = "^MT-")
VlnPlot(combined_CRPC_raw,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "orig.ident", ncol = 3, pt.size = 0
)

# FeatureScatter plot - separate panels for each patient (P1, P2, P3)
split_obj <- SplitObject(combined_CRPC_raw, split.by = "orig.ident")
scatter_plots <- lapply(names(split_obj), function(sample) {
    p1 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(sample)
    p2 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(sample)
    p1 + p2
})
combined_scatter <- wrap_plots(scatter_plots, ncol = 1)
ggsave("Results/FeatureScatter_by_patient.png", plot = combined_scatter, width = 12, height = 15)

# QC
combined_CRPC <- subset(combined_CRPC_raw,
    subset = nFeature_RNA > 800 &
        nFeature_RNA < 9000 &
        percent.mt < 20
)

# Normalization(SCTransform) & PCA ----
combined_CRPC <- SCTransform(combined_CRPC, verbose = FALSE)
combined_CRPC <- RunPCA(combined_CRPC, verbose = FALSE)

# Integration ----
# CCA Integration
combined_CRPC <- IntegrateLayers(
    object = combined_CRPC,
    method = CCAIntegration,
    normalization.method = "SCT",
    verbose = FALSE
)

# Harmony Integration
combined_CRPC <- IntegrateLayers(
    object = combined_CRPC,
    method = HarmonyIntegration,
    normalization.method = "SCT",
    verbose = FALSE
)
