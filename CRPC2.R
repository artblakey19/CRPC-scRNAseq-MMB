library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(ggplot2)

source("scRNA_utils.R")

# Load & QC ----
CRPC2 <- utils_load_and_qc(
    data_dir =
        "~/Human CRPC scRNAseq/Raw data/CRPC2/filtered_feature_bc_matrix",
    project_name = "CRPC2",
    min_features = 800, max_features = 9000, max_mt = 15
)

# Normalization, VST, Scaling, PCA ----
CRPC2 <- utils_run_standard_preprocessing(CRPC2, nfeatures = 2000)
ElbowPlot(CRPC2, ndims = 50)
# Clustering ----
CRPC2 <- FindNeighbors(CRPC2, dims = 1:30)
CRPC2 <- FindClusters(CRPC2, resolution = 0.2)
# Dimensional reduction by UMAP
CRPC2 <- RunUMAP(CRPC2, dims = 1:30)
DimPlot(CRPC2, reduction = "umap", label = TRUE)
# GSEA each PCs to check biological feasibility
# source("PCA_GSEA.R")
# gsea<-run_pc_gsea(CRPC2, num_pcs = 50)

# Cell Cycle Identification ----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
# Cell cycle scoring
CRPC2 <- CellCycleScoring(CRPC2, s.features = s.genes, g2m.features = g2m.genes)
# Draw on UMAP & Save
p <- DimPlot(CRPC2, group.by = "Phase", pt.size = 0.3)
ggsave(
    "CRPC2/CellCycle_Phase_UMAP.tiff",
    plot = p,
    width = 10, height = 7, dpi = 300
)

# Cell Annotation ----
ref <- celldex::HumanPrimaryCellAtlasData()
CRPC2 <- utils_run_singleR_annotation(CRPC2, ref_data = ref)
# Draw UMAP & Save
p <- DimPlot(CRPC2, group.by = "celltype", label = TRUE, repel = TRUE)
ggsave(
    "CRPC2/SingleR_Annotations.tiff",
    plot = p,
    width = 15, height = 7, dpi = 300
)

# Check Markers & Save
p <- FeaturePlot(CRPC2,
    features = c(
        "EPCAM", "KRT8", "AR", "PTPRC",
        "VIM", "PECAM1", "ACTA2", "SYP"
    ),
    ncol = 4, pt.size = 0.1
)
ggsave(
    "CRPC2/Annotation_Markers.tiff",
    plot = p,
    width = 10, height = 7, dpi = 300
)
# Find Markers & Save
all_markers <- utils_save_all_markers(
    CRPC2,
    output_csv = "CRPC2/CRPC2_all_markers.csv",
    logfc_threshold = 0.5
)

# Label Annotation
CRPC2 <- RenameIdents(
    object = CRPC2,
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
    `12` = "Mast"
)
p <- DimPlot(CRPC2, reduction = "umap", label = TRUE)
ggsave(
    "CRPC2/SingleR_Annotations.tiff",
    plot = p,
    width = 10, height = 7, dpi = 300
)
