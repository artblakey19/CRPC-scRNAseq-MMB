# scRNA_utils.R
# Useful functions for scRNAseq analysis
# Loading data & QC: utils_load_and_qc()
# Normalization, VST, Scale, PCA: utils_run_standard_preprocessing()
# SingleR Annotation: utils_run_singleR_annotation()
# Find and Save All Markers: utils_save_all_markers()

library(Seurat)
library(SingleR)
library(SingleCellExperiment)

# Load & QC ----
#' 1. Load & QC
#' @param data_dir Path to the 10X Genomics dataset (filtered_feature_bc_matrix)
#' @param project_name Name of the project
#' @param min_features Minimum number of features(genes) to retain a cell
#' @param max_features Maximum number of features(genes) to filter out doublets
#' @param max_mt Maximum allowed mitochondrial percentage
#' @return A QC-filtered Seurat object
utils_load_and_qc <- function(
  data_dir, project_name, min_features = 800, max_features = 9000,
  max_mt = 15
) {
  raw_data <- Read10X(data.dir = data_dir)
  seurat_obj <- CreateSeuratObject(
    counts = raw_data, project = project_name,
    min.cells = 3, min.features = 200
  )
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^MT-")

  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > min_features &
      nFeature_RNA < max_features &
      percent.mt < max_mt
  )
  return(seurat_obj)
}

# Normalize, VST, Scale, PCA ----
#' 2. Normalize, VST, Scale, PCA
#' @param seurat_obj A Seurat object
#' @param nfeatures Number of highly variable features to select (default: 2000)
#' @return A Seurat object with PCA dimensionality reduction
utils_run_standard_preprocessing <- function(seurat_obj, nfeatures = 2000) {
  seurat_obj <- NormalizeData(seurat_obj)
  seurat_obj <- FindVariableFeatures(
    seurat_obj,
    selection.method = "vst", nfeatures = nfeatures
  )
  seurat_obj <- ScaleData(seurat_obj, features = rownames(seurat_obj))
  seurat_obj <- RunPCA(
    seurat_obj,
    features = VariableFeatures(object = seurat_obj)
  )
  return(seurat_obj)
}

# SingleR Annotation ----
#' 3. SingleR Annotation
#' @param seurat_obj A Seurat object
#' @param ref_data Reference dataset via celldex (e.g., HumanPrimaryCellAtlasData())
#' @return A Seurat object with basic celltype predictions in the metadata
utils_run_singleR_annotation <- function(seurat_obj, ref_data) {
  sce <- as.SingleCellExperiment(seurat_obj)
  predictions <- SingleR(
    test = sce, ref = ref_data, labels = ref_data$label.main
  )

  seurat_obj$celltype <- predictions$labels
  return(seurat_obj)
}

# Find and Save All Markers ----
#' 4. Find and Save All Markers
#' @param seurat_obj A Seurat object
#' @param output_csv Path to save the markers CSV
#' @param logfc_threshold Log2 fold change threshold (default: 0.5)
#' @return A dataframe of all markers
utils_save_all_markers <- function(
  seurat_obj, output_csv, logfc_threshold = 0.5
) {
  all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = logfc_threshold
  )
  write.csv(all_markers, file = output_csv, row.names = FALSE)
  return(all_markers)
}
