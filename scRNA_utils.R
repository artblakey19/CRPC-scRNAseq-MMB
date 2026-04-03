# scRNA_utils.R
# Useful functions for scRNAseq analysis
# Loading data & QC: utils_load_and_qc()
# Normalization, VST, Scale, PCA: utils_run_standard_preprocessing()
# SingleR Annotation: utils_run_singleR_annotation()
# Find and Save All Markers: utils_save_all_markers()
# GSEA by cluster: utils_run_cluster_gsea()

library(Seurat)
library(SingleR)
library(SingleCellExperiment)
library(clusterProfiler)
library(msigdbr)
library(openxlsx)
library(future.apply)

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
  # PrepSCTFindMarkers is required for SCTransform normalized data
  if (DefaultAssay(seurat_obj) == "SCT") {
    seurat_obj <- PrepSCTFindMarkers(seurat_obj)
  }

  all_markers <- FindAllMarkers(
    seurat_obj,
    only.pos = TRUE,
    min.pct = 0.25, logfc.threshold = logfc_threshold
  )
  write.csv(all_markers, file = output_csv, row.names = FALSE)
  return(all_markers)
}

# GSEA by cluster ----
#' 5. GSEA by cluster (Hallmark, KEGG_MEDICUS, GO:BP)
#' @param seurat_obj A Seurat object with clustering results
#' @param output_xlsx Path to save the GSEA results Excel file
#' @param ident Column name for cluster identity (default: active ident)
#' @param species "Homo sapiens" or "Mus musculus"
#' @param pval_cutoff P-value cutoff for GSEA results (default: 0.05)
#' @return Invisibly returns a list of GSEA result data frames per collection
utils_run_cluster_gsea <- function(
  seurat_obj,
  output_xlsx,
  ident = NULL,
  species = "Homo sapiens",
  pval_cutoff = 0.05
) {
  # Set identity if specified
  if (!is.null(ident)) {
    Idents(seurat_obj) <- ident
  }

  # PrepSCTFindMarkers for SCT assay
  if (DefaultAssay(seurat_obj) == "SCT") {
    seurat_obj <- PrepSCTFindMarkers(seurat_obj)
  }

  # Retrieve gene set collections from msigdbr
  hallmark_df <- msigdbr(species = species, category = "H")
  kegg_med_df <- msigdbr(species = species, category = "C2", subcategory = "CP:KEGG_MEDICUS")
  gobp_df <- msigdbr(species = species, category = "C5", subcategory = "GO:BP")

  collections <- list(
    Hallmark     = hallmark_df[, c("gs_name", "gene_symbol")],
    KEGG_MEDICUS = kegg_med_df[, c("gs_name", "gene_symbol")],
    GO_BP        = gobp_df[, c("gs_name", "gene_symbol")]
  )

  clusters <- sort(unique(Idents(seurat_obj)))
  message(
    "Running GSEA for ", length(clusters), " clusters x ",
    length(collections), " collections..."
  )

  # Parallel GSEA per cluster via future.apply
  # (inherits plan() set by the user, e.g. plan("multicore", workers = 30))
  cluster_results <- future_lapply(clusters, function(cl) {
    message("  Cluster: ", cl)

    # DEG: cluster vs all others → ranked gene list by avg_log2FC
    markers <- FindMarkers(
      seurat_obj,
      ident.1 = cl,
      min.pct = 0.1,
      logfc.threshold = 0,
      only.pos = FALSE
    )
    gene_rank <- setNames(markers$avg_log2FC, rownames(markers))
    gene_rank <- sort(gene_rank, decreasing = TRUE)

    res_per_coll <- list()
    for (coll_name in names(collections)) {
      gs <- collections[[coll_name]]

      gsea_res <- tryCatch(
        {
          res <- GSEA(
            geneList = gene_rank,
            TERM2GENE = gs,
            pvalueCutoff = pval_cutoff,
            pAdjustMethod = "BH",
            minGSSize = 10,
            maxGSSize = 500,
            verbose = FALSE
          )
          as.data.frame(res)
        },
        error = function(e) {
          message(
            "    [", coll_name, "] No significant results or error: ",
            conditionMessage(e)
          )
          data.frame()
        }
      )

      if (nrow(gsea_res) > 0) {
        gsea_res$cluster <- cl
        gsea_res <- gsea_res[, c("cluster", setdiff(names(gsea_res), "cluster"))]
      }
      res_per_coll[[coll_name]] <- gsea_res
    }
    return(res_per_coll)
  }, future.seed = TRUE)

  # Merge parallel results
  results_list <- list(
    Hallmark     = data.frame(),
    KEGG_MEDICUS = data.frame(),
    GO_BP        = data.frame()
  )
  for (cr in cluster_results) {
    for (coll_name in names(results_list)) {
      if (nrow(cr[[coll_name]]) > 0) {
        results_list[[coll_name]] <- rbind(results_list[[coll_name]], cr[[coll_name]])
      }
    }
  }

  # Write to Excel: one sheet per collection
  wb <- createWorkbook()
  for (sheet_name in names(results_list)) {
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, results_list[[sheet_name]])
  }
  saveWorkbook(wb, output_xlsx, overwrite = TRUE)
  message("GSEA results saved to: ", output_xlsx)

  invisible(results_list)
}
