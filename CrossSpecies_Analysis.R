# Cross-Species Harmony Integration
# Mouse DNPC Model (DoubleTg + TripleTg) vs Human CRPC Biopsy (CRPC1/2/3)

# Hypothesis:
#   CRPC1 (mostly normal/benign) → clusters with DoubleTg (PIN)
#   CRPC2, CRPC3 (cancer)        → clusters with TripleTg (DNPC-like)

library(Seurat)
library(harmony)
library(dplyr)
library(ggplot2)
library(patchwork)
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scRNA_utils.R")

dir.create("Results/Cross_Species", recursive = TRUE, showWarnings = FALSE)

# Load Data ----
# Human epithelial (already reclustered)
human_epi <- readRDS("Results/epi.rds")

# Mouse integrated epithelial (DoubleTg + TripleTg, cell cycle regressed)
mouse_epi <- readRDS("TG Mouse/TriplevDouble.combined1.epi.rds")

# Ortholog Gene Name Conversion (Mouse → Human) ----
DefaultAssay(human_epi) <- "RNA"
DefaultAssay(mouse_epi) <- "RNA"

# Load ortholog table from Ensembl BioMart (release 115, GRCm39 / GRCh38)
orth_raw <- read.delim("TG Mouse/mart_release115_GRCm39_GRCh38.txt")
colnames(orth_raw) <- c("mouse_gene", "human_gene", "homology_type", "confidence")

# Remove empty values and low confidence orthologs
orth_1to1 <- orth_raw %>%
    filter(
        mouse_gene != "",
        human_gene != "",
        homology_type == "ortholog_one2one",
        confidence == 1
    ) %>%
    distinct(mouse_gene, human_gene)

# Create mouse gene counts matrix for manipulation
mouse_counts <- LayerData(mouse_epi, assay = "RNA", layer = "counts")

# Keep only genes with orthologs
mouse_counts <- mouse_counts[rownames(mouse_counts) %in% orth_1to1$mouse_gene, ]

# Match human gene symbols to mouse gene order
gene_map <- orth_1to1[match(rownames(mouse_counts), orth_1to1$mouse_gene), ]
rownames(mouse_counts) <- gene_map$human_gene

# Create new mouse Seurat object with human gene symbols
mouse_epi_human <- CreateSeuratObject(
    counts = mouse_counts,
    meta.data = mouse_epi@meta.data[colnames(mouse_counts), ],
    project = "Mouse_as_Human"
)

# species label
mouse_epi_human$species <- "Mouse"
human_epi$species <- "Human"

# Keep only common genes
common_genes <- intersect(rownames(human_epi), rownames(mouse_epi_human))

human_sub <- subset(human_epi, features = common_genes)
mouse_sub <- subset(mouse_epi_human, features = common_genes)

# Per-species normalization and HVG selection for shared features
obj_list <- list(
    human = NormalizeData(human_sub),
    mouse = NormalizeData(mouse_sub)
)
obj_list <- lapply(obj_list, FindVariableFeatures)
shared_hvg <- SelectIntegrationFeatures(obj_list)

# Merge human and mouse subsets
CombinedEpi <- merge(human_sub, mouse_sub, project = "CombinedEpi")

# Normalize and use shared HVG for PCA
CombinedEpi <- NormalizeData(CombinedEpi)
VariableFeatures(CombinedEpi) <- shared_hvg
CombinedEpi <- ScaleData(CombinedEpi, features = rownames(CombinedEpi))
CombinedEpi <- RunPCA(CombinedEpi, features = shared_hvg)

# Harmony Integration
CombinedEpi <- RunHarmony(CombinedEpi,
    group.by.vars = "species",
    theta = 4, max.iter.harmony = 30
)

# UMAP
CombinedEpi <- RunUMAP(CombinedEpi, reduction = "harmony", dims = 1:30)

# Clustering
CombinedEpi <- FindNeighbors(CombinedEpi, reduction = "harmony", dims = 1:30)
CombinedEpi <- FindClusters(CombinedEpi, resolution = 0.5)

# Save
# saveRDS(CombinedEpi, "Results/Cross_Species/CombinedEpi.rds")

# UMAP Visualization
# Combined UMAP colored by species
p1 <- DimPlot(CombinedEpi, group.by = "species", cols = c("Human" = "#E64B35", "Mouse" = "#4DBBD5")) +
    ggtitle("Cross-Species Integration (Harmony)") +
    theme(plot.title = element_text(hjust = 0.5))

# Combined UMAP colored by cluster
p2 <- DimPlot(CombinedEpi, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle("Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

# Split by species
p3 <- DimPlot(CombinedEpi,
    split.by = "species", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, ncol = 2
) +
    ggtitle("Clusters Split by Species") +
    theme(plot.title = element_text(hjust = 0.5))

# Save plots
ggsave("Results/Cross_Species/UMAP_species_harmony.png",
    p1 + p2,
    width = 16, height = 7, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_split_by_species_harmony.png",
    p3,
    width = 16, height = 7, dpi = 300
)

# CCA Integration ----
anchors <- FindIntegrationAnchors(obj_list,
    anchor.features = shared_hvg,
    reduction = "cca", dims = 1:30
)
CombinedEpi_CCA <- IntegrateData(anchors, dims = 1:30)

DefaultAssay(CombinedEpi_CCA) <- "integrated"
CombinedEpi_CCA <- ScaleData(CombinedEpi_CCA)
CombinedEpi_CCA <- RunPCA(CombinedEpi_CCA)
CombinedEpi_CCA <- RunUMAP(CombinedEpi_CCA, reduction = "pca", dims = 1:30)
CombinedEpi_CCA <- FindNeighbors(CombinedEpi_CCA, reduction = "pca", dims = 1:30)
CombinedEpi_CCA <- FindClusters(CombinedEpi_CCA, resolution = 0.5)

# CCA UMAP Visualization
p4 <- DimPlot(CombinedEpi_CCA, group.by = "species", cols = c("Human" = "#E64B35", "Mouse" = "#4DBBD5")) +
    ggtitle("Cross-Species Integration (CCA)") +
    theme(plot.title = element_text(hjust = 0.5))

p5 <- DimPlot(CombinedEpi_CCA, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle("Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

p6 <- DimPlot(CombinedEpi_CCA,
    split.by = "species", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, ncol = 2
) +
    ggtitle("Clusters Split by Species") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("Results/Cross_Species/UMAP_species_CCA.png",
    p4 + p5,
    width = 16, height = 7, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_split_by_species_CCA.png",
    p6,
    width = 16, height = 7, dpi = 300
)
