library(dplyr)
library(Seurat)
library(patchwork)
library(celldex)
library(ggplot2)
library(presto)

source("scRNA_utils.R")

# Load & QC ----
CRPC1 <- utils_load_and_qc(
  data_dir =
    "~/Human CRPC scRNAseq/Raw data/CRPC1/filtered_feature_bc_matrix",
  project_name = "CRPC1",
  min_features = 800, max_features = 9000, max_mt = 15
)

# Normalization, VST, Scaling, PCA ----
CRPC1 <- utils_run_standard_preprocessing(CRPC1, nfeatures = 2000)
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
ggsave(
  "Results/CRPC1/CellCycle_Phase_UMAP.tiff",
  plot = p,
  width = 10, height = 7, dpi = 300
)

# Cell Annotation ----
ref <- celldex::HumanPrimaryCellAtlasData()
CRPC1 <- utils_run_singleR_annotation(CRPC1, ref_data = ref)
# Draw UMAP & Save
p <- DimPlot(CRPC1, group.by = "celltype", label = TRUE, repel = TRUE)
ggsave(
  "Results/CRPC1/SingleR_Annotations.tiff",
  plot = p,
  width = 15, height = 7, dpi = 300
)

# Check Markers & Save
p <- FeaturePlot(CRPC1,
  features = c(
    "EPCAM", "KRT8", "AR", "PTPRC",
    "VIM", "PECAM1", "ACTA2", "SYP"
  ),
  ncol = 4, pt.size = 0.1
)
ggsave(
  "Results/CRPC1/Annotation_Markers.tiff",
  plot = p,
  width = 10, height = 7, dpi = 300
)
# Find Markers & Save
all_markers <- utils_save_all_markers(
  CRPC1,
  output_csv = "Results/CRPC1/CRPC1_all_markers.csv",
  logfc_threshold = 0.5
)

# Label Annotation
CRPC1 <- RenameIdents(
  object = CRPC1,
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
p <- DimPlot(CRPC1, reduction = "umap", label = TRUE)
ggsave(
  "Results/CRPC1/SingleR_Annotations.tiff",
  plot = p,
  width = 10, height = 7, dpi = 300
)
