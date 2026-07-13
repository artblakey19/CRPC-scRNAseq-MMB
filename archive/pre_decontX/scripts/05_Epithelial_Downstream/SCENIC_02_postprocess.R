#!/usr/bin/env Rscript
# Read the AUCell loom from pySCENIC back into R, attach to epi_annotated,
# compute per-cluster (and per-patient) regulon activity, and plot top regulons.

suppressPackageStartupMessages({
  library(Seurat)
  library(SCopeLoomR)
  library(SCENIC)        # for AUCell helpers (optional)
  library(tidyverse)
  library(pheatmap)
  library(AUCell)
})

PROJ <- "/home/MMB/projects/Human CRPC scRNAseq"
OUT_DIR <- file.path(PROJ, "Results/05_Epithelial_Downstream/SCENIC")
AUC_LOOM <- file.path(OUT_DIR, "epi_pyscenic_aucell.loom")
ANNOT_LEVELS <- file.path(PROJ, "Results/05_Epithelial_Downstream/Annotation/label_levels.txt")
stopifnot(file.exists(AUC_LOOM))

message("[1/4] Reading AUC loom")
loom <- open_loom(AUC_LOOM)
regulons_incidMat <- get_regulons(loom, column.attr.name = "Regulons")
regulons          <- regulonsToGeneLists(regulons_incidMat)
regulon_auc       <- get_regulons_AUC(loom, column.attr.name = "RegulonsAUC")
close_loom(loom)

auc_mat <- AUCell::getAUC(regulon_auc)   # regulon x cell
message("regulons: ", nrow(auc_mat), "  cells: ", ncol(auc_mat))

message("[2/4] Attaching to epi_annotated")
epi <- readRDS(file.path(PROJ, "Results/05_Epithelial_Downstream/epi_annotated.rds"))
stopifnot("annotation column missing" = "annotation" %in% colnames(epi@meta.data))
common_cells <- intersect(colnames(epi), colnames(auc_mat))
message("  shared cells: ", length(common_cells))
auc_mat <- auc_mat[, common_cells]
epi <- epi[, common_cells]
annotation <- as.character(epi$annotation)
if (file.exists(ANNOT_LEVELS)) {
  cluster_order <- readLines(ANNOT_LEVELS, warn = FALSE)
} else {
  cluster_order <- levels(epi$annotation)
}
if (is.null(cluster_order)) cluster_order <- unique(annotation)
cluster_order <- cluster_order[cluster_order %in% unique(annotation)]
clusters <- factor(annotation, levels = cluster_order)
message("  annotation levels: ", paste(levels(clusters), collapse = " / "))

epi[["AUC"]] <- CreateAssayObject(counts = auc_mat)
DefaultAssay(epi) <- "AUC"

message("[3/4] Per-cluster mean regulon activity")
mean_by_cluster <- sapply(split(seq_along(clusters), clusters), function(idx) {
  rowMeans(auc_mat[, idx, drop = FALSE])
})
mean_by_cluster <- mean_by_cluster[, levels(clusters), drop = FALSE]
write.csv(mean_by_cluster, file.path(OUT_DIR, "regulon_mean_by_cluster.csv"))

# Regulon Specificity Score: how cluster-specific each regulon is
rss <- function(auc_mat, clusters) {
  cluster_levels <- levels(droplevels(clusters))
  freq <- table(clusters)/length(clusters)
  out <- matrix(NA, nrow = nrow(auc_mat), ncol = length(cluster_levels),
                dimnames = list(rownames(auc_mat), cluster_levels))
  # Jensen-Shannon-based specificity (PCA-MarkerEnrichr / Suo 2018)
  for (r in rownames(auc_mat)) {
    p <- auc_mat[r, ] / sum(auc_mat[r, ])
    if (any(is.na(p)) || sum(p) == 0) next
    for (cl in cluster_levels) {
      q <- as.numeric(clusters == cl); q <- q / sum(q)
      m <- 0.5 * (p + q)
      jsd <- 0.5 * sum(ifelse(p>0, p*log2(p/m), 0)) +
             0.5 * sum(ifelse(q>0, q*log2(q/m), 0))
      out[r, cl] <- 1 - sqrt(jsd)
    }
  }
  out
}
rss_mat <- rss(auc_mat, clusters)
rss_mat <- rss_mat[, levels(clusters), drop = FALSE]
write.csv(rss_mat, file.path(OUT_DIR, "regulon_specificity_score.csv"))

message("[4/4] Plots")
# Top-10 regulons by RSS per cluster
top_per_cluster <- apply(rss_mat, 2, function(x) names(sort(x, decreasing = TRUE))[1:10])
top_regs <- unique(as.vector(top_per_cluster))

pheatmap::pheatmap(
  mean_by_cluster[top_regs, ],
  scale       = "row",
  color       = colorRampPalette(c("steelblue","white","firebrick"))(50),
  main        = "Top RSS regulons × cluster (z-scored mean AUC)",
  filename    = file.path(OUT_DIR, "top_regulons_by_cluster_heatmap.pdf"),
  width       = 8, height = max(8, length(top_regs)*0.18),
  fontsize_row = 6
)

write.csv(data.frame(top_per_cluster),
          file.path(OUT_DIR, "top10_regulons_per_cluster.csv"))

saveRDS(epi, file.path(OUT_DIR, "epi_with_scenic_auc.rds"))
message("Done.  outputs: ", OUT_DIR)
