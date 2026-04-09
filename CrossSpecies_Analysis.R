# Cross-Species Integration (Harmony / LogNorm-CCA / SCT-CCA)
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
# Human integrated epithelial (Three patients: CRPC1/2/3)
human_epi <- readRDS("Results/Epithelial/epi_SCT_CCA.rds")

# Mouse integrated epithelial (DoubleTg + TripleTg)
# Rebuild epi2 with EpiCellTypes if not already saved
if (!file.exists("TG Mouse/TriplevDouble.combined1.epi2.rds")) {
    source("TG Mouse/rebuild_epi2.R")
}
mouse_epi <- readRDS("TG Mouse/TriplevDouble.combined1.epi2.rds")

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
CombinedEpi <- ScaleData(CombinedEpi,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(CombinedEpi)
)
CombinedEpi <- RunPCA(CombinedEpi, features = shared_hvg)

# Evaluate PCA dimensions
p_elbow <- ElbowPlot(CombinedEpi, ndims = 50) +
    ggtitle("LogNorm+Harmony Elbow Plot")
ggsave("Results/Cross_Species/ElbowPlot_LogNorm_Harmony.png",
    p_elbow,
    width = 7, height = 5, dpi = 300
)

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

# Create unified sample label: CRPC1/2/3 from orig.ident, DoubleTg/TripleTg from stim
CombinedEpi$sample <- ifelse(CombinedEpi$species == "Human",
    CombinedEpi$orig.ident,
    ifelse(CombinedEpi$stim == "Double", "DoubleTg", "TripleTg")
)
CombinedEpi$sample <- factor(CombinedEpi$sample,
    levels = c("CRPC1", "CRPC2", "CRPC3", "DoubleTg", "TripleTg")
)

sample_cols <- c(
    "CRPC1" = "#E64B35", "CRPC2" = "#F39B7F", "CRPC3" = "#B09C85",
    "DoubleTg" = "#4DBBD5", "TripleTg" = "#3C5488"
)

# UMAP Visualization ----
# Combined UMAP colored by species
p1 <- DimPlot(CombinedEpi, group.by = "species", cols = c("Human" = "#E64B35", "Mouse" = "#4DBBD5")) +
    ggtitle("LogNorm+Harmony: Species") +
    theme(plot.title = element_text(hjust = 0.5))

# Combined UMAP colored by cluster
p2 <- DimPlot(CombinedEpi, group.by = "seurat_clusters", label = TRUE, repel = TRUE) +
    ggtitle("LogNorm+Harmony: Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

# Combined UMAP colored by sample
p3 <- DimPlot(CombinedEpi, group.by = "sample", cols = sample_cols) +
    ggtitle("LogNorm+Harmony: Samples") +
    theme(plot.title = element_text(hjust = 0.5))

# Split by sample
p4 <- DimPlot(CombinedEpi,
    split.by = "sample", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, ncol = 5
) +
    ggtitle("LogNorm+Harmony: Clusters by Sample") +
    theme(plot.title = element_text(hjust = 0.5))

# Save plots
ggsave("Results/Cross_Species/UMAP_LogNorm_Harmony_overview.png",
    p1 + p2 + p3 + plot_layout(ncol = 3),
    width = 24, height = 7, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_LogNorm_Harmony_split_by_sample.png",
    p4,
    width = 28, height = 7, dpi = 300
)

# LogNorm + CCA Integration ----
anchors <- FindIntegrationAnchors(obj_list,
    anchor.features = shared_hvg,
    reduction = "cca", dims = 1:30
)
CombinedEpi_CCA <- IntegrateData(anchors, dims = 1:30)

DefaultAssay(CombinedEpi_CCA) <- "integrated"
CombinedEpi_CCA <- ScaleData(CombinedEpi_CCA)
CombinedEpi_CCA <- RunPCA(CombinedEpi_CCA)

# Evaluate LogNorm+CCA PCA dimensions
p_elbow_cca <- ElbowPlot(CombinedEpi_CCA, ndims = 50) + ggtitle("LogNorm+CCA Elbow Plot")
ggsave("Results/Cross_Species/ElbowPlot_LogNorm_CCA.png", p_elbow_cca, width = 7, height = 5, dpi = 300)

CombinedEpi_CCA <- RunUMAP(CombinedEpi_CCA, reduction = "pca", dims = 1:30)
CombinedEpi_CCA <- FindNeighbors(CombinedEpi_CCA, reduction = "pca", dims = 1:30)
CombinedEpi_CCA <- FindClusters(CombinedEpi_CCA, resolution = 0.5)

# Create unified sample label for CCA object
CombinedEpi_CCA$sample <- ifelse(CombinedEpi_CCA$species == "Human",
    CombinedEpi_CCA$orig.ident,
    ifelse(CombinedEpi_CCA$stim == "Double", "DoubleTg", "TripleTg")
)
CombinedEpi_CCA$sample <- factor(CombinedEpi_CCA$sample,
    levels = c("CRPC1", "CRPC2", "CRPC3", "DoubleTg", "TripleTg")
)

# LogNorm+CCA UMAP Visualization
species_cols <- c("Human" = "#E64B35", "Mouse" = "#4DBBD5")

p5 <- DimPlot(CombinedEpi_CCA,
    group.by = "species", cols = species_cols
) +
    ggtitle("LogNorm+CCA: Species") +
    theme(plot.title = element_text(hjust = 0.5))

p6 <- DimPlot(CombinedEpi_CCA,
    group.by = "seurat_clusters",
    label = TRUE, repel = TRUE
) +
    ggtitle("LogNorm+CCA: Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

p7 <- DimPlot(CombinedEpi_CCA,
    group.by = "sample", cols = sample_cols
) +
    ggtitle("LogNorm+CCA: Samples") +
    theme(plot.title = element_text(hjust = 0.5))

p8 <- DimPlot(CombinedEpi_CCA,
    split.by = "sample", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, ncol = 5
) +
    ggtitle("LogNorm+CCA: Clusters by Sample") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("Results/Cross_Species/UMAP_LogNorm_CCA_overview.png",
    p5 + p6 + p7 + plot_layout(ncol = 3),
    width = 24, height = 7, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_LogNorm_CCA_split_by_sample.png",
    p8,
    width = 28, height = 7, dpi = 300
)

# SCTransform + Harmony Integration ----
sct_harmony_list <- list(
    human = SCTransform(human_sub,
        vars.to.regress = c("S.Score", "G2M.Score"),
        verbose = FALSE
    ),
    mouse = SCTransform(mouse_sub,
        vars.to.regress = c("S.Score", "G2M.Score"),
        verbose = FALSE
    )
)
sct_harmony_hvg <- SelectIntegrationFeatures(sct_harmony_list)

CombinedEpi_SCT_Harmony <- merge(sct_harmony_list$human, sct_harmony_list$mouse,
    project = "CombinedEpi_SCT_Harmony"
)
VariableFeatures(CombinedEpi_SCT_Harmony) <- sct_harmony_hvg
CombinedEpi_SCT_Harmony <- RunPCA(CombinedEpi_SCT_Harmony, features = sct_harmony_hvg)

p_elbow_sct_harmony <- ElbowPlot(CombinedEpi_SCT_Harmony, ndims = 50) +
    ggtitle("SCT+Harmony Elbow Plot")
ggsave("Results/Cross_Species/ElbowPlot_SCT_Harmony.png",
    p_elbow_sct_harmony,
    width = 7, height = 5, dpi = 300
)

CombinedEpi_SCT_Harmony <- RunHarmony(CombinedEpi_SCT_Harmony,
    group.by.vars = "species",
    theta = 4, max.iter.harmony = 30
)

CombinedEpi_SCT_Harmony <- RunUMAP(CombinedEpi_SCT_Harmony,
    reduction = "harmony", dims = 1:30
)
CombinedEpi_SCT_Harmony <- FindNeighbors(CombinedEpi_SCT_Harmony,
    reduction = "harmony", dims = 1:30
)
CombinedEpi_SCT_Harmony <- FindClusters(CombinedEpi_SCT_Harmony, resolution = 0.5)

# Create unified sample label
CombinedEpi_SCT_Harmony$sample <- ifelse(
    CombinedEpi_SCT_Harmony$species == "Human",
    CombinedEpi_SCT_Harmony$orig.ident,
    ifelse(CombinedEpi_SCT_Harmony$stim == "Double", "DoubleTg", "TripleTg")
)
CombinedEpi_SCT_Harmony$sample <- factor(
    CombinedEpi_SCT_Harmony$sample,
    levels = c("CRPC1", "CRPC2", "CRPC3", "DoubleTg", "TripleTg")
)

# SCT+Harmony UMAP Visualization
p_sh1 <- DimPlot(CombinedEpi_SCT_Harmony,
    group.by = "species", cols = c("Human" = "#E64B35", "Mouse" = "#4DBBD5")
) +
    ggtitle("SCT+Harmony: Species") +
    theme(plot.title = element_text(hjust = 0.5))

p_sh2 <- DimPlot(CombinedEpi_SCT_Harmony,
    group.by = "seurat_clusters",
    label = TRUE, repel = TRUE
) +
    ggtitle("SCT+Harmony: Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

p_sh3 <- DimPlot(CombinedEpi_SCT_Harmony,
    group.by = "sample", cols = sample_cols
) +
    ggtitle("SCT+Harmony: Samples") +
    theme(plot.title = element_text(hjust = 0.5))

p_sh4 <- DimPlot(CombinedEpi_SCT_Harmony,
    split.by = "sample", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, ncol = 5
) +
    ggtitle("SCT+Harmony: Clusters by Sample") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("Results/Cross_Species/UMAP_SCT_Harmony_overview.png",
    p_sh1 + p_sh2 + p_sh3 + plot_layout(ncol = 3),
    width = 24, height = 7, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_SCT_Harmony_split_by_sample.png",
    p_sh4,
    width = 28, height = 7, dpi = 300
)

# SCTransform + CCA Integration ----
sct_list <- list(
    human = SCTransform(human_sub,
        vars.to.regress = c("S.Score", "G2M.Score"),
        verbose = FALSE
    ),
    mouse = SCTransform(mouse_sub,
        vars.to.regress = c("S.Score", "G2M.Score"),
        verbose = FALSE
    )
)
sct_hvg <- SelectIntegrationFeatures(sct_list)
sct_list <- PrepSCTIntegration(sct_list, anchor.features = sct_hvg)

sct_anchors <- FindIntegrationAnchors(sct_list,
    normalization.method = "SCT",
    anchor.features = sct_hvg,
    reduction = "cca", dims = 1:30
)
CombinedEpi_SCT <- IntegrateData(sct_anchors,
    normalization.method = "SCT", dims = 1:30
)

DefaultAssay(CombinedEpi_SCT) <- "integrated"
CombinedEpi_SCT <- ScaleData(CombinedEpi_SCT)
CombinedEpi_SCT <- RunPCA(CombinedEpi_SCT)

p_elbow_sct <- ElbowPlot(CombinedEpi_SCT, ndims = 50) +
    ggtitle("SCT+CCA Elbow Plot")
ggsave("Results/Cross_Species/ElbowPlot_SCT_CCA.png",
    p_elbow_sct,
    width = 7, height = 5, dpi = 300
)

CombinedEpi_SCT <- RunUMAP(CombinedEpi_SCT,
    reduction = "pca", dims = 1:30
)
CombinedEpi_SCT <- FindNeighbors(CombinedEpi_SCT,
    reduction = "pca", dims = 1:30
)
CombinedEpi_SCT <- FindClusters(CombinedEpi_SCT, resolution = 0.3)

# Create unified sample label for SCT object
CombinedEpi_SCT$sample <- ifelse(
    CombinedEpi_SCT$species == "Human",
    CombinedEpi_SCT$orig.ident,
    ifelse(CombinedEpi_SCT$stim == "Double",
        "DoubleTg", "TripleTg"
    )
)
CombinedEpi_SCT$sample <- factor(
    CombinedEpi_SCT$sample,
    levels = c(
        "CRPC1", "CRPC2", "CRPC3",
        "DoubleTg", "TripleTg"
    )
)

# SCT+CCA UMAP Visualization
p9 <- DimPlot(CombinedEpi_SCT,
    group.by = "species", cols = species_cols
) +
    ggtitle("SCT+CCA: Species") +
    theme(plot.title = element_text(hjust = 0.5))

p10 <- DimPlot(CombinedEpi_SCT,
    group.by = "seurat_clusters",
    label = TRUE, repel = TRUE
) +
    ggtitle("SCT+CCA: Clusters") +
    theme(plot.title = element_text(hjust = 0.5))

p11 <- DimPlot(CombinedEpi_SCT,
    group.by = "sample", cols = sample_cols
) +
    ggtitle("SCT+CCA: Samples") +
    theme(plot.title = element_text(hjust = 0.5))

p12 <- DimPlot(CombinedEpi_SCT,
    split.by = "sample", group.by = "seurat_clusters",
    label = TRUE, repel = TRUE, ncol = 5
) +
    ggtitle("SCT+CCA: Clusters by Sample") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("Results/Cross_Species/UMAP_SCT_CCA_overview.png",
    p9 + p10 + p11 + plot_layout(ncol = 3),
    width = 24, height = 7, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_SCT_CCA_split_by_sample.png",
    p12,
    width = 28, height = 7, dpi = 300
)

# Mouse EpiCellTypes + Human sample overlay on SCT+CCA UMAP ----
# Create label: Mouse → EpiCellTypes (BE1-4, LE1-9, UrLE, OE), Human → CRPC1/2/3
CombinedEpi_SCT$celltype_label <- ifelse(
    CombinedEpi_SCT$species == "Mouse",
    as.character(CombinedEpi_SCT$EpiCellTypes),
    as.character(CombinedEpi_SCT$orig.ident)
)

epi_levels <- c(
    "CRPC1", "CRPC2", "CRPC3",
    "BE1", "BE2", "BE3", "BE4",
    "LE1", "LE2", "LE3", "LE4", "LE5", "LE6", "LE7", "LE8", "LE9",
    "UrLE", "OE"
)
epi_levels <- epi_levels[epi_levels %in% CombinedEpi_SCT$celltype_label]
CombinedEpi_SCT$celltype_label <- factor(
    CombinedEpi_SCT$celltype_label,
    levels = epi_levels
)

celltype_cols <- c(
    "CRPC1" = "#E64B35", "CRPC2" = "#F39B7F", "CRPC3" = "#B09C85",
    "BE1" = "salmon", "BE2" = "olivedrab2", "BE3" = "#98FB98", "BE4" = "brown3",
    "LE1" = "deeppink1", "LE2" = "blue", "LE3" = "darkorange",
    "LE4" = "red", "LE5" = "turquoise3", "LE6" = "bisque3",
    "LE7" = "slategray3", "LE8" = "mediumorchid3", "LE9" = "skyblue1",
    "UrLE" = "yellow2", "OE" = "green4"
)
celltype_cols <- celltype_cols[names(celltype_cols) %in% epi_levels]

# Panel 1: Mouse cells only, colored by EpiCellTypes (reference map)
mouse_cells <- Cells(CombinedEpi_SCT)[CombinedEpi_SCT$species == "Mouse"]
p_mouse_ref <- DimPlot(CombinedEpi_SCT,
    cells = mouse_cells,
    group.by = "celltype_label", cols = celltype_cols,
    label = TRUE, repel = TRUE, label.size = 3.5
) +
    ggtitle("Mouse EpiCellTypes (Reference)") +
    theme(plot.title = element_text(hjust = 0.5))

# Panel 2-4: Each CRPC patient highlighted on grey background
plot_highlight <- function(obj, cells_highlight, title, col) {
    DimPlot(obj,
        cells.highlight = cells_highlight,
        cols.highlight = col, cols = "grey85",
        sizes.highlight = 0.3, pt.size = 0.1
    ) +
        ggtitle(title) +
        theme(plot.title = element_text(hjust = 0.5)) +
        NoLegend()
}

p_crpc1 <- plot_highlight(
    CombinedEpi_SCT,
    Cells(CombinedEpi_SCT)[CombinedEpi_SCT$orig.ident == "CRPC1"],
    "CRPC1", "#E64B35"
)
p_crpc2 <- plot_highlight(
    CombinedEpi_SCT,
    Cells(CombinedEpi_SCT)[CombinedEpi_SCT$orig.ident == "CRPC2"],
    "CRPC2", "#F39B7F"
)
p_crpc3 <- plot_highlight(
    CombinedEpi_SCT,
    Cells(CombinedEpi_SCT)[CombinedEpi_SCT$orig.ident == "CRPC3"],
    "CRPC3", "#B09C85"
)

# Panel 5-6: DoubleTg / TripleTg highlighted separately
p_double <- plot_highlight(
    CombinedEpi_SCT,
    Cells(CombinedEpi_SCT)[CombinedEpi_SCT$species == "Mouse" &
        CombinedEpi_SCT$stim == "Double"],
    "DoubleTg", "#4DBBD5"
)
p_triple <- plot_highlight(
    CombinedEpi_SCT,
    Cells(CombinedEpi_SCT)[CombinedEpi_SCT$species == "Mouse" &
        CombinedEpi_SCT$stim == "Triple"],
    "TripleTg", "#3C5488"
)

# Save: reference + all 5 samples highlighted
ggsave("Results/Cross_Species/UMAP_SCT_CCA_EpiCellTypes_reference.png",
    p_mouse_ref,
    width = 10, height = 8, dpi = 300
)
ggsave("Results/Cross_Species/UMAP_SCT_CCA_samples_on_reference.png",
    p_mouse_ref + p_crpc1 + p_crpc2 + p_crpc3 + p_double + p_triple +
        plot_layout(ncol = 3),
    width = 24, height = 14, dpi = 300
)

# Export for scANVI ----
# Merge raw counts with metadata for Python scANVI integration
combined_for_scanvi <- merge(human_sub, mouse_sub, project = "CrossSpecies")
combined_for_scanvi$sample <- ifelse(
    combined_for_scanvi$species == "Human",
    combined_for_scanvi$orig.ident,
    ifelse(combined_for_scanvi$stim == "Double", "DoubleTg", "TripleTg")
)
# scANVI label: Mouse → EpiCellTypes (BE1-4, LE1-9, UrLE, OE), Human → "Unknown"
combined_for_scanvi$scanvi_label <- ifelse(
    combined_for_scanvi$species == "Mouse",
    as.character(combined_for_scanvi$EpiCellTypes),
    "Unknown"
)

# Export counts + metadata for Python
combined_for_scanvi[["RNA"]] <- JoinLayers(combined_for_scanvi[["RNA"]])
counts_mat <- LayerData(combined_for_scanvi, assay = "RNA", layer = "counts")
Matrix::writeMM(counts_mat, "Results/Cross_Species/scanvi_counts.mtx")
writeLines(rownames(counts_mat), "Results/Cross_Species/scanvi_genes.txt")
writeLines(colnames(counts_mat), "Results/Cross_Species/scanvi_barcodes.txt")
write.csv(combined_for_scanvi@meta.data, "Results/Cross_Species/scanvi_metadata.csv")
