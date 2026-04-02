library(dplyr)
library(Seurat)
library(patchwork)
library(SingleR)
library(celldex)
library(SingleCellExperiment)
library(ggplot2)

# Load the PBMC dataset
CRPC1.raw.data <- Read10X(data.dir = "~/Human CRPC scRNAseq/Raw data/CRPC1/filtered_feature_bc_matrix")
# Initialize the Seurat object with the raw (non-normalized data).
CRPC1.raw <- CreateSeuratObject(counts = CRPC1.raw.data, project = "CRPC1", min.cells = 3, min.features = 200)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
CRPC1.raw[["percent.mt"]] <- PercentageFeatureSet(CRPC1.raw, pattern = "^MT-")
# Visualize QC metrics as a violin plot
VlnPlot(CRPC1.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

# # FeatureScatter is typically used to visualize feature-feature relationships, but can be used
# # for anything calculated by the object, i.e. columns in object metadata, PC scores etc.
# plot1 <- FeatureScatter(CRPC1.raw, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(CRPC1.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1+plot2

# QC ----
CRPC1 <- subset(CRPC1.raw, subset = nFeature_RNA > 800 & nFeature_RNA < 9000 & percent.mt < 15)

# plot1 <- FeatureScatter(CRPC1, feature1 = "nCount_RNA", feature2 = "percent.mt")
# plot2 <- FeatureScatter(CRPC1, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
# plot1+plot2

# Normalization ----
CRPC1 <- NormalizeData(CRPC1)
# VST
CRPC1 <- FindVariableFeatures(CRPC1, selection.method = "vst", nfeatures = 2000)

# # Identify the 10 most highly variable genes
# top10 <- head(VariableFeatures(CRPC1), 10)
# # plot variable features with and without labels
# plot1 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)

# Scaling ----
CRPC1 <- ScaleData(CRPC1, features = rownames(CRPC1))
# PCA
CRPC1 <- RunPCA(CRPC1, features = VariableFeatures(object = CRPC1))
ElbowPlot(CRPC1, ndims = 50)
# Clustering ----
CRPC1 <- FindNeighbors(CRPC1, dims = 1:30)
CRPC1 <- FindClusters(CRPC1, resolution = 0.2)
# Dimensional reduction by UMAP
CRPC1 <- RunUMAP(CRPC1, dims = 1:30)
DimPlot(CRPC1, reduction = "umap", label = TRUE)
# GSEA each PCs to check biological feasibility
# source("PCA_GSEA.R")
# gsea<-run_pc_gsea(CRPC1, num_pcs = 50)

# Cell Cycle Identification ----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
# Cell cycle scoring
CRPC1 <- CellCycleScoring(CRPC1, s.features = s.genes, g2m.features = g2m.genes)
# Draw on UMAP & Save
p <- DimPlot(CRPC1, group.by = "Phase", pt.size = 0.3)
ggsave("CRPC1/CellCycle_Phase_UMAP.tiff", plot = p, width = 10, height = 7, dpi = 300)

# Cell Annotation ----
# Use SingleR to basic annotation
# Load reference Data
ref <- HumanPrimaryCellAtlasData()
# SingleR uses Bioconductor(SingleCellExperiment, SCE) objects, not Seurat objects
# Convert Seurat -> SCE by as.SingleCellExperiment() function(from Seurat lib)
sce <- as.SingleCellExperiment(CRPC1)
predictions <- SingleR(test = sce, ref = ref, labels = ref$label.main)
# Add annotations to Seurat object
CRPC1$celltype <- predictions$labels
# Draw UMAP & Save
p <- DimPlot(CRPC1, group.by = "celltype", label = TRUE, repel = TRUE)
ggsave("CRPC1/SingleR_Annotations.tiff", plot = p, width = 15, height = 7, dpi = 300)

# Check Markers & Save
p <- FeaturePlot(CRPC1, features = c("EPCAM", "KRT8", "AR", "PTPRC", 
                                "VIM", "PECAM1", "ACTA2", "SYP"), 
            ncol = 4, pt.size = 0.1)
ggsave("CRPC1/Annotation_Markers.tiff", plot = p, width = 10, height = 7, dpi = 300)
# Find Markers & Save
all_markers <- FindAllMarkers(CRPC1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5)
write.csv(all_markers, file = "CRPC1/CRPC1_all_markers.csv", row.names = FALSE)

# Label Annotation
CRPC1 <- RenameIdents(object = CRPC1, 
                      `0` = "Epithelial",
                      `1` = "Epithelial",
                      `2` = "Epithelial",
                      `3` = "Epithelial",
                      `4` = "B/T",
                      `5` = "Fibroblast", 
                      `6` = "Endothelial",
                      `7` = "Smooth Muscle", 
                      `8` = "Leukocyte", 
                      `9` = "Epithelial", 
                      `10` = "Immune", 
                      `11` = "Epithelial", 
                      `12` = "Mast")
p <- DimPlot(CRPC1, reduction = "umap", label = TRUE)
ggsave("CRPC1/SingleR_Annotations.tiff", plot = p, width = 10, height = 7, dpi = 300)



