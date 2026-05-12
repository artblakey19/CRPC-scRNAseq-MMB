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
bpp <- MulticoreParam(workers = 30, RNGseed = 42)

source("scRNA_utils.R")

run_copykat <- TRUE
copykat_cores <- 30
copykat_reference_celltypes <- c(
    "T/NK cells", "Phagocytes", "Immune cells"
)

if (run_copykat) {
    if (!requireNamespace("copykat", quietly = TRUE)) {
        stop(
            "copykat package is required. Install it with: remotes::install_github('navinlabcode/copykat')",
            call. = FALSE
        )
    }
    library(copykat)
}

# Load 3 Samples ----
p1 <- Read10X(data.dir = "Raw data/CRPC1/filtered_feature_bc_matrix/")
p2 <- Read10X(data.dir = "Raw data/CRPC2/filtered_feature_bc_matrix/")
p3 <- Read10X(data.dir = "Raw data/CRPC3/filtered_feature_bc_matrix/")
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
p <- VlnPlot(combined_CRPC_raw,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "orig.ident", ncol = 3, pt.size = 0
)
ggsave("Results/QC/mt_by_patient_VlnPlot.png", plot = p, width = 12, height = 8)

# VlnPlot with threshold lines to verify QC cutoffs
qc_thresholds <- list(
    nFeature_RNA = c(200, 9000),
    nCount_RNA   = NULL,
    percent.mt   = c(NA, 20)
)
qc_vln <- lapply(names(qc_thresholds), function(feat) {
    vp <- VlnPlot(combined_CRPC_raw,
        features = feat,
        group.by = "orig.ident", pt.size = 0
    ) + NoLegend()
    th <- qc_thresholds[[feat]]
    if (!is.null(th)) {
        for (y in th) if (!is.na(y)) vp <- vp + geom_hline(yintercept = y, linetype = "dashed", color = "red")
    }
    vp
})
qc_vln_combined <- wrap_plots(qc_vln, ncol = 3)
ggsave("Results/QC/QC_VlnPlot_with_thresholds.png", plot = qc_vln_combined, width = 14, height = 6)

# FeatureScatter plot - separate panels for each patient (P1, P2, P3)
split_obj <- SplitObject(combined_CRPC_raw, split.by = "orig.ident")
scatter_plots <- lapply(names(split_obj), function(sample) {
    p1 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(sample)
    p2 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(sample)
    p1 + p2
})
combined_scatter <- wrap_plots(scatter_plots, ncol = 1)
ggsave("Results/QC/FeatureScatter_by_patient.png", plot = combined_scatter, width = 12, height = 15)

# Filter
combined_CRPC <- subset(combined_CRPC_raw,
    subset = nFeature_RNA > 200 &
        nFeature_RNA < 9000 &
        percent.mt < 20
)

# copyKAT helper functions ----
# copyKAT is executed after cell type annotation so each patient can use
# same-sample immune cells as an explicit diploid reference.
dir.create("Results/copyKAT", showWarnings = FALSE, recursive = TRUE)

extract_copykat_prediction <- function(copykat_res, sample_id) {
    if (is.null(copykat_res$prediction)) {
        stop("copyKAT returned no prediction table for ", sample_id, call. = FALSE)
    }

    pred <- as.data.frame(copykat_res$prediction)
    cell_col <- intersect(c("cell.names", "cell.name", "cell", "barcode"), colnames(pred))
    pred_col <- intersect(c("copykat.pred", "copykat_prediction", "prediction"), colnames(pred))

    if (length(cell_col) == 0 || length(pred_col) == 0) {
        stop("copyKAT prediction columns were not recognized for ", sample_id, call. = FALSE)
    }

    pred <- pred[, c(cell_col[1], pred_col[1]), drop = FALSE]
    colnames(pred) <- c("cell", "copykat_prediction")
    pred$copykat_sample <- sample_id
    pred
}

select_copykat_reference_cells <- function(
    seu,
    sample_id,
    reference_celltypes
) {
    sample_meta <- seu@meta.data[seu$orig.ident == sample_id, , drop = FALSE]
    ref_cells <- rownames(sample_meta)[sample_meta$celltype %in% reference_celltypes]
    reference_source <- "immune"

    if (length(ref_cells) == 0) {
        stop("No immune reference cells found for ", sample_id, call. = FALSE)
    }

    list(cells = ref_cells, source = reference_source)
}

run_copykat_sample <- function(
    seu,
    sample_id,
    output_dir,
    normal_cells = character(),
    reference_source = "immune",
    n_cores = 30
) {
    message("Running copyKAT for ", sample_id)

    sample_name <- gsub("[^A-Za-z0-9_.-]+", "_", sample_id)
    sample_dir <- file.path(output_dir, sample_name)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)

    sample_cells <- colnames(seu)[seu$orig.ident == sample_id]
    normal_cells <- intersect(normal_cells, sample_cells)
    norm_cell_arg <- if (length(normal_cells) > 0) normal_cells else ""
    message("  normal reference mode: ", reference_source, " (n = ", length(normal_cells), ")")

    sample_obj <- subset(seu, cells = sample_cells)
    DefaultAssay(sample_obj) <- "RNA"
    sample_obj[["RNA"]] <- JoinLayers(sample_obj[["RNA"]])

    raw_counts <- GetAssayData(sample_obj, assay = "RNA", layer = "counts")
    raw_counts <- as.matrix(raw_counts)
    storage.mode(raw_counts) <- "numeric"

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(sample_dir)
    copykat_res <- copykat::copykat(
        rawmat = raw_counts,
        id.type = "S",
        cell.line = "no",
        ngene.chr = 5,
        win.size = 25,
        norm.cell.names = norm_cell_arg,
        KS.cut = 0.1,
        sam.name = sample_name,
        distance = "euclidean",
        genome = "hg20",
        n.cores = n_cores
    )
    setwd(old_wd)

    saveRDS(copykat_res, file.path(sample_dir, paste0(sample_name, "_copykat.rds")))
    write.csv(
        data.frame(cell = normal_cells, copykat_reference_source = reference_source),
        file.path(sample_dir, paste0(sample_name, "_copykat_reference_cells.csv")),
        row.names = FALSE
    )
    pred <- extract_copykat_prediction(copykat_res, sample_id)
    pred$copykat_reference_source <- reference_source
    pred$copykat_reference_n <- length(normal_cells)
    write.csv(pred, file.path(sample_dir, paste0(sample_name, "_copykat_predictions.csv")), row.names = FALSE)

    rm(sample_obj, raw_counts, copykat_res)
    gc()
    pred
}

# Normalization(SCTransform) & PCA ----
combined_CRPC <- SCTransform(combined_CRPC, verbose = FALSE)
combined_CRPC <- RunPCA(combined_CRPC, verbose = FALSE)
# Draw elbow plot
ggsave("Results/QC/ElbowPlot.png", plot = ElbowPlot(combined_CRPC, ndims = 50), width = 12, height = 15)

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
# Draw on UMAP & Save
p <- DimPlot(combined_CRPC, group.by = "Phase", pt.size = 0.3)
ggsave(
    "Results/Integrated/CellCycle_Phase_UMAP.tiff",
    plot = p,
    width = 10, height = 7, dpi = 300
)

# Annotation ----
# Check marker genes
p <- FeaturePlot(combined_CRPC,
    features = c(
        "EPCAM", "KRT8", "AR", "PTPRC",
        "VIM", "PECAM1", "ACTA2", "SYP",
        "COL1A1", "DCN", "TPSAB1"
    ),
    ncol = 4, pt.size = 0.1
)
ggsave("Results/Integrated/FeaturePlot_markers.png", plot = p, width = 20, height = 19)

# Find all markers
utils_save_all_markers(combined_CRPC, "Results/Integrated/all_markers.csv")

# Save cluster UMAP (annotation 전 단계 — cluster 번호로 확인)
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "seurat_clusters", label = TRUE)
ggsave("Results/Integrated/Cluster_UMAP_integrated.png", plot = p, width = 15, height = 15)

# Save cluster UMAP by patient
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "seurat_clusters", split.by = "orig.ident", label = TRUE)
ggsave("Results/Integrated/Cluster_UMAP_by_patient.png", plot = p, width = 24, height = 15)

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
    "Results/Integrated/SingleR_cell_labels.csv",
    row.names = FALSE
)

# Diagnostic plots — score heatmap & delta distribution
png("Results/Integrated/SingleR_score_heatmap.png", width = 1200, height = 900)
plotScoreHeatmap(singler_pred)
dev.off()

png("Results/Integrated/SingleR_delta_distribution.png",
    width = 1200, height = 900
)
plotDeltaDistribution(singler_pred, ncol = 4)
dev.off()

# UMAP with SingleR labels (pruned — uncertain cells appear as NA/grey)
p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "SingleR_pruned",
    label = TRUE, repel = TRUE
)
ggsave("Results/Integrated/SingleR_UMAP.png",
    plot = p, width = 15, height = 15
)

# UMAP with SingleR labels by patient
p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "SingleR",
    split.by = "orig.ident", label = TRUE, repel = TRUE
)
ggsave("Results/Integrated/SingleR_UMAP_by_patient.png",
    plot = p, width = 24, height = 15
)

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

# copyKAT CNV inference ----
# Run per patient on QC-filtered raw RNA counts. Immune cells are used as the
# same-sample diploid reference.
if (run_copykat) {
    copykat_samples <- sort(unique(combined_CRPC$orig.ident))
    copykat_reference_info <- setNames(
        lapply(copykat_samples, function(sample_id) {
            select_copykat_reference_cells(
                seu = combined_CRPC,
                sample_id = sample_id,
                reference_celltypes = copykat_reference_celltypes
            )
        }),
        copykat_samples
    )

    copykat_reference_summary <- bind_rows(lapply(copykat_samples, function(sample_id) {
        ref_info <- copykat_reference_info[[sample_id]]
        data.frame(
            orig.ident = sample_id,
            copykat_reference_source = ref_info$source,
            copykat_reference_n = length(ref_info$cells)
        )
    }))
    write.csv(
        copykat_reference_summary,
        "Results/copyKAT/copykat_reference_summary.csv",
        row.names = FALSE
    )

    copykat_pred <- bind_rows(lapply(copykat_samples, function(sample_id) {
        ref_info <- copykat_reference_info[[sample_id]]
        run_copykat_sample(
            seu = combined_CRPC,
            sample_id = sample_id,
            output_dir = "Results/copyKAT",
            normal_cells = ref_info$cells,
            reference_source = ref_info$source,
            n_cores = copykat_cores
        )
    }))

    write.csv(copykat_pred, "Results/copyKAT/copykat_predictions_all_samples.csv", row.names = FALSE)

    copykat_idx <- match(colnames(combined_CRPC), copykat_pred$cell)
    combined_CRPC$copykat_prediction <- copykat_pred$copykat_prediction[copykat_idx]
    combined_CRPC$copykat_sample <- copykat_pred$copykat_sample[copykat_idx]
    combined_CRPC$copykat_reference_source <- copykat_pred$copykat_reference_source[copykat_idx]
    combined_CRPC$copykat_reference_n <- copykat_pred$copykat_reference_n[copykat_idx]
} else {
    combined_CRPC$copykat_prediction <- NA_character_
    combined_CRPC$copykat_sample <- NA_character_
    combined_CRPC$copykat_reference_source <- NA_character_
    combined_CRPC$copykat_reference_n <- NA_integer_
}

copykat_pred_lower <- tolower(combined_CRPC$copykat_prediction)
combined_CRPC$copykat_malignant <- case_when(
    copykat_pred_lower == "aneuploid" ~ "malignant",
    copykat_pred_lower == "diploid" ~ "non-malignant",
    TRUE ~ "not.defined"
)

copykat_patient_summary <- combined_CRPC@meta.data %>%
    count(orig.ident, copykat_prediction, copykat_malignant, name = "n")
write.csv(
    copykat_patient_summary,
    "Results/copyKAT/copykat_patient_summary.csv",
    row.names = FALSE
)

copykat_cluster_summary <- combined_CRPC@meta.data %>%
    count(orig.ident, seurat_clusters, copykat_prediction, copykat_malignant, name = "n")
write.csv(
    copykat_cluster_summary,
    "Results/copyKAT/copykat_cluster_summary.csv",
    row.names = FALSE
)

p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "copykat_malignant",
    pt.size = 0.3
)
ggsave("Results/copyKAT/copykat_malignant_UMAP.png",
    plot = p, width = 10, height = 8
)

p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "copykat_malignant",
    split.by = "orig.ident", pt.size = 0.3
)
ggsave("Results/copyKAT/copykat_malignant_UMAP_by_patient.png",
    plot = p, width = 24, height = 8
)

copykat_pred_lower <- tolower(combined_CRPC$copykat_prediction)
combined_CRPC$malignant_status <- case_when(
    combined_CRPC$celltype == "Epithelial" & copykat_pred_lower == "aneuploid" ~ "Malignant epithelial",
    combined_CRPC$celltype == "Epithelial" & copykat_pred_lower == "diploid" ~ "Diploid epithelial",
    combined_CRPC$celltype != "Epithelial" & copykat_pred_lower == "aneuploid" ~ "Aneuploid non-epithelial",
    combined_CRPC$celltype != "Epithelial" & copykat_pred_lower == "diploid" ~ "Diploid non-epithelial",
    TRUE ~ "Not defined"
)

copykat_celltype_summary <- combined_CRPC@meta.data %>%
    count(orig.ident, celltype, malignant_status, name = "n")
write.csv(
    copykat_celltype_summary,
    "Results/copyKAT/copykat_celltype_summary.csv",
    row.names = FALSE
)

copykat_meta <- combined_CRPC@meta.data
copykat_meta$cell <- rownames(copykat_meta)
copykat_meta <- copykat_meta[, c("cell", setdiff(colnames(copykat_meta), "cell"))]
write.csv(
    copykat_meta[
        !is.na(copykat_meta$malignant_status) &
            copykat_meta$malignant_status == "Malignant epithelial",
    ],
    "Results/copyKAT/malignant_epithelial_cells.csv",
    row.names = FALSE
)

p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "malignant_status",
    pt.size = 0.3
)
ggsave("Results/copyKAT/malignant_status_UMAP.png",
    plot = p, width = 12, height = 9
)

p <- DimPlot(combined_CRPC,
    reduction = "umap", group.by = "malignant_status",
    split.by = "orig.ident", pt.size = 0.3
)
ggsave("Results/copyKAT/malignant_status_UMAP_by_patient.png",
    plot = p, width = 24, height = 8
)

# Save labelled UMAP
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", label = TRUE)
ggsave("Results/Integrated/Labelled_UMAP_integrated.png", plot = p, width = 15, height = 15)

# Save labelled UMAP by patient
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", split.by = "orig.ident", label = TRUE)
ggsave("Results/Integrated/Labelled_UMAP_by_patient.png", plot = p, width = 24, height = 15)

# Save RDS ----
saveRDS(combined_CRPC, "Results/Integrated/combined_CRPC.rds")
# Read RDS
# combined_CRPC <- readRDS("Results/Integrated/combined_CRPC.rds")
