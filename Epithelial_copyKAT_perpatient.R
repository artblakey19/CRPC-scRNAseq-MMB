# ============================================================
# Epithelial_copyKAT_perpatient.R
#
# Per-patient copyKAT CNV inference on the epithelial subset, followed by
# subclone calling and subclone-tree visualisation.
#
# Reference (diploid anchor): cluster 8 — a benign-appearing epithelial
# population transcriptionally close to the DNPC malignant states. copyKAT is
# run once PER PATIENT so each tumour's acquired CNVs are called against that
# patient's own cluster-8 baseline rather than a cross-patient pool (pooling
# folds patient-specific germline/baseline differences into the CNV signal).
#
# Pipeline (per patient):
#   copyKAT  ->  cut hclust tree into subclones  ->  subclone trees + CNA heatmap
# ============================================================

library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)
library(data.table)
library(ComplexHeatmap)
library(circlize)
library(patchwork)
library(cluster)

library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)
Sys.setenv(OMP_NUM_THREADS = "30", OPENBLAS_NUM_THREADS = "30", MKL_NUM_THREADS = "30")

# ============================================================
# Config ----
# ============================================================
REF_CLUSTER <- "8"          # diploid anchor cluster
K_RANGE     <- 2:5          # candidate subclone counts scanned per patient
# The primary `copykat_subclone` label is no longer a fixed K. Per patient we
# pick the K in K_RANGE that maximises the average silhouette width of the
# CNA-profile clustering (Euclidean distance + Ward.D2 linkage), following
# Wang et al., Brief. Bioinform. 26(2):bbaf076 (2025). bbaf076 also merges
# similar subclones post hoc; we omit that step — it reintroduces an arbitrary
# similarity threshold, and silhouette-optimal K already penalises over-splitting.
N_CORES     <- 30
ROOT_OUT    <- "Results/Epithelial/copyKAT_perpatient"
dir.create(ROOT_OUT, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load epithelial subset ----
# ============================================================
epi <- readRDS("Results/Epithelial/epi_clustered.rds")
DefaultAssay(epi) <- "RNA"
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
epi$seurat_clusters <- as.character(epi$seurat_clusters)

patients <- sort(unique(epi$orig.ident))
message("Patients: ", paste(patients, collapse = ", "))

ck_cols  <- c(aneuploid = "#D7261E", diploid = "#1F77B4", not.defined = "grey80")
cna_cols <- colorRamp2(c(-0.10, 0, 0.10), c("#3B4CC0", "white", "#B40426"))

# ============================================================
# Helper: run copyKAT for one patient ----
# ============================================================
# copyKAT writes intermediate files to the working directory, so we chdir into
# a per-patient subdir and restore the wd on exit.
run_copykat_patient <- function(sample_id) {
    message("\n==== copyKAT: ", sample_id, " ====")
    sub <- epi[, epi$orig.ident == sample_id]
    ref_cells <- colnames(sub)[sub$seurat_clusters == REF_CLUSTER]
    message("  cells = ", ncol(sub),
            " | cluster-", REF_CLUSTER, " reference = ", length(ref_cells))
    if (length(ref_cells) < 50)
        warning("  reference < 50 cells for ", sample_id,
                " — CNV calls may be unstable", call. = FALSE)

    raw <- as.matrix(GetAssayData(sub, assay = "RNA", layer = "counts"))
    storage.mode(raw) <- "numeric"

    sample_dir <- file.path(ROOT_OUT, sample_id)
    dir.create(sample_dir, showWarnings = FALSE, recursive = TRUE)
    write.csv(data.frame(cell = ref_cells),
              file.path(sample_dir, paste0(sample_id, "_reference_cells.csv")),
              row.names = FALSE)

    old_wd <- getwd()
    on.exit(setwd(old_wd), add = TRUE)
    setwd(sample_dir)
    res <- copykat::copykat(
        rawmat          = raw,
        id.type         = "S",
        cell.line       = "no",
        ngene.chr       = 5,
        win.size        = 25,
        norm.cell.names = ref_cells,
        KS.cut          = 0.1,
        sam.name        = sample_id,
        distance        = "euclidean",
        genome          = "hg20",
        n.cores         = N_CORES
    )
    setwd(old_wd)

    saveRDS(res, file.path(sample_dir, paste0(sample_id, "_copykat.rds")))
    rm(raw, sub); gc()
    res
}

# ============================================================
# Helper: cell dendrogram with subclone / cluster / copyKAT ribbons ----
# ============================================================
plot_cell_dendrogram <- function(sample_id, hc, ann, sample_dir) {
    ann <- ann[hc$order, ]
    clu_levels <- sort(unique(as.character(ann$cluster)))
    clu_cols   <- setNames(scales::hue_pal()(length(clu_levels)), clu_levels)
    sub_levels <- sort(unique(as.character(ann$subclone)))
    sub_cols   <- setNames(
        RColorBrewer::brewer.pal(max(3, length(sub_levels)), "Set1")[seq_along(sub_levels)],
        sub_levels)

    png(file.path(sample_dir, paste0(sample_id, "_subclone_dendrogram.png")),
        width = 14, height = 7, units = "in", res = 200, bg = "white")
    layout(matrix(1:4, nrow = 4), heights = c(5, 0.6, 0.6, 0.6))
    par(mar = c(0, 6, 3, 1))
    plot(as.dendrogram(hc), leaflab = "none",
         main = paste0(sample_id, " — copyKAT hclust (cells), ribbons = subclone / cluster / ploidy"))
    par(mar = c(0, 6, 0, 1))
    image(seq_len(nrow(ann)), 1,
          matrix(match(as.character(ann$subclone), sub_levels), ncol = 1),
          col = sub_cols, axes = FALSE, xlab = "", ylab = "subclone")
    par(mar = c(0, 6, 0, 1))
    image(seq_len(nrow(ann)), 1,
          matrix(match(as.character(ann$cluster), clu_levels), ncol = 1),
          col = clu_cols, axes = FALSE, xlab = "", ylab = "cluster")
    par(mar = c(2, 6, 0, 1))
    image(seq_len(nrow(ann)), 1,
          matrix(match(as.character(ann$copykat), names(ck_cols)), ncol = 1),
          col = ck_cols, axes = FALSE, xlab = "", ylab = "ploidy")
    dev.off()

    png(file.path(sample_dir, paste0(sample_id, "_subclone_dendrogram_legend.png")),
        width = 7, height = 4, units = "in", res = 200, bg = "white")
    par(mar = c(0, 0, 0, 0)); plot.new()
    legend("topleft",   legend = sub_levels, fill = sub_cols,
           title = "subclone", bty = "n")
    legend("top",       legend = clu_levels, fill = clu_cols,
           title = "seurat_cluster", ncol = 2, bty = "n")
    legend("topright",  legend = names(ck_cols), fill = ck_cols,
           title = "copyKAT", bty = "n")
    dev.off()
}

# ============================================================
# Per-patient run ----
# ============================================================
pred_list  <- list()
sub_list   <- list()
bestk_list <- list()

for (sid in patients) {
    res <- run_copykat_patient(sid)
    sample_dir <- file.path(ROOT_OUT, sid)

    # ---- prediction table ----
    pred <- as.data.frame(res$prediction)
    cell_col <- intersect(c("cell.names", "cell_names", "cells"), names(pred))[1]
    pred_col <- intersect(c("copykat.pred", "copykat_pred", "prediction"), names(pred))[1]
    pred_df <- data.frame(
        cell            = sub("\\.1$", "-1", pred[[cell_col]]),
        copykat_clu8ref = pred[[pred_col]],
        stringsAsFactors = FALSE)

    # ---- CNA matrix written by copyKAT ----
    cna <- data.table::fread(
        file.path(sample_dir, paste0(sid, "_copykat_CNA_results.txt")),
        showProgress = FALSE)
    bin_chr <- cna$chrom
    mat <- as.matrix(cna[, -(1:3)])
    colnames(mat) <- sub("\\.1$", "-1", colnames(mat))
    rm(cna); gc()

    # ---- hierarchical clustering of CNA profiles (Euclidean + Ward.D2, bbaf076) ----
    cna_dist <- dist(t(mat))
    hc <- hclust(cna_dist, method = "ward.D2")
    stopifnot(all(hc$labels %in% colnames(epi)))

    # ---- pick K by maximum average silhouette width over K_RANGE ----
    sil_avg <- vapply(K_RANGE, function(k)
        mean(cluster::silhouette(cutree(hc, k = k), cna_dist)[, "sil_width"]),
        numeric(1))
    names(sil_avg) <- K_RANGE
    best_k <- as.integer(names(sil_avg)[which.max(sil_avg)])
    message("  avg silhouette: ",
            paste(sprintf("K%d=%.3f", K_RANGE, sil_avg), collapse = "  "),
            "  -> K=", best_k)

    png(file.path(sample_dir, paste0(sid, "_silhouette_K.png")),
        width = 5, height = 4, units = "in", res = 200, bg = "white")
    par(mar = c(4, 4, 3, 1))
    plot(K_RANGE, sil_avg, type = "b", pch = 19,
         xlab = "K (subclones)", ylab = "mean silhouette width",
         main = paste0(sid, " — silhouette-optimal K = ", best_k))
    abline(v = best_k, col = "#D7261E", lty = 2)
    dev.off()

    # ---- cut hclust tree into subclones (patient-prefixed, not cross-comparable) ----
    sub_df <- data.frame(cell = hc$labels, stringsAsFactors = FALSE)
    for (k in K_RANGE) {
        v <- cutree(hc, k = k)
        sub_df[[paste0("subclone_k", k)]] <- paste0(sid, "_S", v[hc$labels])
    }
    sub_df$subclone <- sub_df[[paste0("subclone_k", best_k)]]
    primary_col <- paste0("subclone_k", best_k)

    # ---- subclone x seurat_cluster crosstab ----
    ann <- data.frame(
        cell     = hc$labels,
        cluster  = epi@meta.data[hc$labels, "seurat_clusters"],
        subclone = sub_df[[primary_col]],
        copykat  = pred_df$copykat_clu8ref[match(hc$labels, pred_df$cell)],
        stringsAsFactors = FALSE)
    write.csv(
        as.data.frame.matrix(table(subclone = ann$subclone, cluster = ann$cluster)),
        file.path(sample_dir, paste0(sid, "_subclone_x_cluster_k", best_k, ".csv")))
    write.csv(
        as.data.frame.matrix(table(subclone = ann$subclone, ploidy = ann$copykat)),
        file.path(sample_dir, paste0(sid, "_subclone_x_ploidy_k", best_k, ".csv")))

    # ---- tree view 1: cell dendrogram with ribbons ----
    plot_cell_dendrogram(sid, hc, ann, sample_dir)

    # ---- aggregate mean CNA per subclone (silhouette-optimal K) ----
    prim <- cutree(hc, k = best_k); names(prim) <- hc$labels
    sub_levels <- sort(unique(prim))
    agg <- vapply(sub_levels, function(s)
        rowMeans(mat[, names(prim)[prim == s], drop = FALSE]),
        numeric(nrow(mat)))
    colnames(agg) <- paste0(sid, "_S", sub_levels)

    # ---- tree view 2: subclone relationship tree (clonal tree) ----
    sub_hc <- hclust(dist(t(agg)), method = "ward.D2")
    png(file.path(sample_dir, paste0(sid, "_subclone_tree_k", best_k, ".png")),
        width = 6, height = 5, units = "in", res = 200, bg = "white")
    par(mar = c(2, 4, 3, 1))
    plot(sub_hc, main = paste0(sid, " — subclone relationship tree (K=", best_k, ")"),
         xlab = "", sub = "", ylab = "CNA distance")
    dev.off()

    # ---- per-subclone mean CNA heatmap ----
    chr_factor <- droplevels(factor(paste0("chr", bin_chr),
                                    levels = paste0("chr", c(1:22, "X", "Y"))))
    n_per_sub <- as.integer(table(prim)[as.character(sub_levels)])
    ht <- Heatmap(
        t(agg), name = "mean CNA", col = cna_cols,
        cluster_rows = TRUE, cluster_columns = FALSE,
        show_column_names = FALSE, row_names_side = "left",
        column_split = chr_factor, column_gap = unit(0.5, "mm"),
        column_title_gp = gpar(fontsize = 8), border = TRUE,
        use_raster = TRUE, raster_quality = 4,
        left_annotation = rowAnnotation(
            n_cells = anno_barplot(n_per_sub, gp = gpar(fill = "grey40"),
                                   width = unit(2, "cm"))),
        row_title = "copyKAT subclone",
        heatmap_legend_param = list(direction = "horizontal",
                                    title_position = "lefttop"))
    png(file.path(sample_dir, paste0(sid, "_subclone_CNA_heatmap_k", best_k, ".png")),
        width = 18, height = 4, units = "in", res = 300, bg = "white")
    draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE,
         column_title = paste0(sid, " — mean copyKAT CNA per subclone (K=", best_k, ")"))
    dev.off()

    pred_list[[sid]]  <- pred_df
    sub_list[[sid]]   <- sub_df
    bestk_list[[sid]] <- best_k
    rm(mat, hc, agg, cna_dist); gc()
    message("  done: ", sid)
}

# ============================================================
# Merge results back into the epithelial object ----
# ============================================================
pred_all <- bind_rows(pred_list)
sub_all  <- bind_rows(sub_list)

epi$copykat_clu8ref <- pred_all$copykat_clu8ref[match(colnames(epi), pred_all$cell)]
for (k in K_RANGE) {
    col <- paste0("subclone_k", k)
    epi[[paste0("copykat_", col)]] <- sub_all[[col]][match(colnames(epi), sub_all$cell)]
}
epi$copykat_subclone <- sub_all$subclone[match(colnames(epi), sub_all$cell)]

saveRDS(epi, file.path(ROOT_OUT, "epi_copykat_perpatient.rds"))

# ============================================================
# Combined summary plots ----
# ============================================================
p_pred <- DimPlot(epi, group.by = "copykat_clu8ref", reduction = "umap",
                  cols = ck_cols, pt.size = 0.3, raster = FALSE) +
    ggtitle("Per-patient copyKAT (cluster 8 reference)")
ggsave(file.path(ROOT_OUT, "UMAP_copykat_clu8ref.png"),
       plot = p_pred, width = 9, height = 7, bg = "white")

p_pred_pat <- DimPlot(epi, group.by = "copykat_clu8ref", split.by = "orig.ident",
                      reduction = "umap", cols = ck_cols, pt.size = 0.3,
                      raster = FALSE)
ggsave(file.path(ROOT_OUT, "UMAP_copykat_clu8ref_by_patient.png"),
       plot = p_pred_pat, width = 20, height = 7, bg = "white")

# subclone labels are patient-specific and NOT cross-comparable, so a single
# 12-colour legend is misleading and the per-panel hues end up poorly separated.
# Strip the patient prefix to the bare S1..S5 id and reuse one fixed palette
# per patient → each facet shows only its own 4 well-separated colours.
epi$subclone_id <- factor(sub("^.*_(S[0-9]+)$", "\\1", epi$copykat_subclone),
                          levels = paste0("S", sort(unique(
                              as.integer(sub("^.*_S", "",
                                  na.omit(epi$copykat_subclone)))))))
sub_id_cols <- setNames(
    RColorBrewer::brewer.pal(max(3, nlevels(epi$subclone_id)), "Set1")[
        seq_len(nlevels(epi$subclone_id))],
    levels(epi$subclone_id))

p_sub <- DimPlot(epi, group.by = "subclone_id", split.by = "orig.ident",
                 reduction = "umap", cols = sub_id_cols, pt.size = 0.3,
                 ncol = length(patients), raster = FALSE) &
    ggtitle(paste0("copyKAT subclone (silhouette-optimal K per patient)",
                   " — S-ids are patient-specific, not cross-comparable")) &
    theme(legend.position = "right")
ggsave(file.path(ROOT_OUT, "UMAP_subclone_by_patient.png"),
       plot = p_sub, width = 7 * length(patients), height = 6, bg = "white")

# subclone composition per seurat cluster, per patient
comp_df <- epi@meta.data |>
    dplyr::filter(!is.na(copykat_subclone)) |>
    dplyr::count(orig.ident, seurat_clusters, copykat_subclone)
p_comp <- ggplot(comp_df, aes(seurat_clusters, n, fill = copykat_subclone)) +
    geom_col(position = "fill") +
    facet_wrap(~orig.ident, ncol = 1) +
    labs(x = "Seurat cluster", y = "Subclone proportion",
         title = "Subclone composition per cluster (silhouette-optimal K)") +
    theme_classic()
ggsave(file.path(ROOT_OUT, "Subclone_composition_by_cluster.png"),
       plot = p_comp, width = 10, height = 9, bg = "white")

# ============================================================
# Report ----
# ============================================================
cat("\n===== copyKAT prediction by patient =====\n")
print(as.data.frame.matrix(
    table(patient = epi$orig.ident, epi$copykat_clu8ref, useNA = "ifany")))
cat("\n===== silhouette-optimal K per patient =====\n")
print(data.frame(patient = names(bestk_list),
                 best_k  = unlist(bestk_list), row.names = NULL))
cat("\n===== subclones per patient (silhouette-optimal K) =====\n")
print(epi@meta.data |>
      dplyr::filter(!is.na(copykat_subclone)) |>
      dplyr::count(orig.ident, copykat_subclone) |>
      as.data.frame())
cat("\nSaved outputs to: ", ROOT_OUT, "\n")
cat("Merged object:    ", file.path(ROOT_OUT, "epi_copykat_perpatient.rds"), "\n")
