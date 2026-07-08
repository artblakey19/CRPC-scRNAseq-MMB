# 02b_RegressOut_Denisenko17_BaselineControl.R
# Corrected regress-out test: use the STORED Stage-04 baseline (the embedding /
# clustering that the annotation actually rests on) as the control, instead of a
# fresh JoinLayers+RunHarmony "control" that does not reproduce the baseline.
#
# Why: script 02's control arm used JoinLayers + RunHarmony(group.by.vars), but the
# baseline (scripts/04_Epithelial_Filtered/Epithelial_FilteredAnalysis.R) used
# split-by-sample SCT + IntegrateLayers(HarmonyIntegration). Different pipeline ->
# the fresh control merged OE 1/BE 1 (not the regression's doing), confounding the
# read-out. Here every recompute replicates the 04 pipeline exactly.
#
# Three arms (annotation overlaid as labels on all):
#   (1) BASELINE  : stored seurat_clusters (== annotation 1:1). The real control.
#   (2) RECOMPUTE : 04 pipeline, NO regression. Validates that the pipeline
#                   reproduces the baseline (does the recompute alone keep OE1/BE1
#                   and OE2/OE3 separate?).
#   (3) REGRESS   : 04 pipeline + SCTransform(vars.to.regress = stress_score).
# Clean isolation of the regression effect = arm (2) vs arm (3) (identical pipeline,
# differ only by vars.to.regress). Faithfulness of recompute = baseline vs (2).
#
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds   (read-only)
# Writes: Results/06_QC_Artifact_Diagnostics/02b_RegressOut_Denisenko17_BaselineControl/
#           arm_summary.csv             ARI(baseline,recompute), ARI(recompute,regress), R2
#           pair_merge_OE2OE3_OE1BE1.csv  co-clustering of the two pairs across 3 arms
#           original_annotation_fate.csv  where each annotation lands (recompute & regress)
#           crosstab_{recompute,regress}_vs_annotation.csv
#           UMAP_baseline_by_annotation.png / UMAP_regress_by_annotation.png
# Existing objects are NOT modified.

suppressMessages({
    library(Seurat); library(dplyr); library(ggplot2)
    library(harmony); library(future)
})
plan("sequential")
options(future.globals.maxSize = 128 * 1024^3)
source("scripts/00_utils/scRNA_utils.R")

set.seed(42)
IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/06_QC_Artifact_Diagnostics/02b_RegressOut_Denisenko17_BaselineControl"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

denisenko17 <- c("FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1",
                 "HSPA1A","HSPA1B","HSP90AB1","HSPA8","HSPB1",
                 "IER3","IER2","BTG1","BTG2","DUSP1")

stopifnot(file.exists(IN_RDS))
epi <- readRDS(IN_RDS)
message("Loaded ", IN_RDS, " (", ncol(epi), " cells)")
orig   <- as.character(epi$annotation)
base_cl <- as.character(epi$seurat_clusters)   # stored Stage-04 baseline (== annotation)

# Stress score on joined, normalized RNA (persists as a per-cell meta column).
DefaultAssay(epi) <- "RNA"
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
epi <- NormalizeData(epi, verbose = FALSE)
den_used <- intersect(denisenko17, rownames(epi[["RNA"]]))
set.seed(42)
epi <- AddModuleScore(epi, features = list(den_used), name = "den17_sc", assay = "RNA")
epi$stress_score <- epi$den17_sc1; epi$den17_sc1 <- NULL
message(sprintf("Denisenko17 stress_score from %d genes", length(den_used)))

# ------------------------------------------------------------------
# Replicate the Stage-04 pipeline exactly; only vars.to.regress differs.
#   split(RNA, orig.ident) -> SCTransform -> RunPCA -> IntegrateLayers(Harmony,SCT)
#   -> JoinLayers -> FindNeighbors(harmony,1:30) -> FindClusters(0.4)
# ------------------------------------------------------------------
run04 <- function(obj, regress, want_umap = FALSE) {
    DefaultAssay(obj) <- "RNA"
    for (r in intersect(c("harmony","pca","umap"), Reductions(obj))) obj[[r]] <- NULL
    if ("SCT" %in% Assays(obj)) obj[["SCT"]] <- NULL
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    obj[["RNA"]] <- split(obj[["RNA"]], f = obj$orig.ident)
    obj <- if (regress) SCTransform(obj, vars.to.regress = "stress_score", verbose = FALSE)
           else         SCTransform(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    set.seed(42)
    obj <- IntegrateLayers(object = obj, method = HarmonyIntegration,
                           orig.reduction = "pca", new.reduction = "harmony",
                           normalization.method = "SCT", verbose = FALSE)
    obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
    obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
    set.seed(42)
    obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
    H  <- Embeddings(obj, "harmony")[, 1:30]
    r2 <- summary(lm(obj$stress_score ~ H))$r.squared
    um <- NULL
    if (want_umap) {
        set.seed(42)
        obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, n.neighbors = 40,
                       min.dist = 0.3, spread = 1.5, seed.use = 42, verbose = FALSE)
        um <- Embeddings(obj, "umap")
    }
    list(cluster = as.character(obj$seurat_clusters), r2 = r2, umap = um)
}

message("ARM 2 — RECOMPUTE (04 pipeline, no regression) ...")
rec <- run04(epi, regress = FALSE)
message("ARM 3 — REGRESS (04 pipeline + stress regressed) ...")
reg <- run04(epi, regress = TRUE, want_umap = TRUE)

# ------------------------------------------------------------------
# ARI + manipulation check
# ------------------------------------------------------------------
ari <- function(a, b) {
    tab <- table(a, b); n <- sum(tab); s_ij <- sum(choose(tab, 2))
    a_i <- sum(choose(rowSums(tab), 2)); b_j <- sum(choose(colSums(tab), 2))
    e <- a_i*b_j/choose(n,2); (s_ij - e)/((a_i+b_j)/2 - e)
}
arm_summary <- data.frame(
    metric = c("ARI(baseline, recompute)  [pipeline faithfulness]",
               "ARI(recompute, regress)   [regression effect]",
               "ARI(baseline, regress)",
               "stress_R2 recompute", "stress_R2 regress",
               "n_clusters baseline", "n_clusters recompute", "n_clusters regress"),
    value = c(round(ari(base_cl, rec$cluster), 3),
              round(ari(rec$cluster, reg$cluster), 3),
              round(ari(base_cl, reg$cluster), 3),
              round(rec$r2, 3), round(reg$r2, 3),
              length(unique(base_cl)), length(unique(rec$cluster)),
              length(unique(reg$cluster))))
write.csv(arm_summary, file.path(OUT_DIR, "arm_summary.csv"), row.names = FALSE)
cat("\n=== Arm summary ===\n"); print(arm_summary, row.names = FALSE)

# ------------------------------------------------------------------
# Pairwise co-clustering: do OE2/OE3 and OE1/BE1 share a cluster, per arm?
# ------------------------------------------------------------------
pair_row <- function(clusters, arm, A, B) {
    isA <- orig == A; isB <- orig == B
    dA <- names(sort(table(clusters[isA]), TRUE))[1]
    dB <- names(sort(table(clusters[isB]), TRUE))[1]
    data.frame(arm = arm, pair = paste(A, "/", B),
               dom_A = dA, dom_B = dB, same_cluster = (dA == dB),
               A_in_domA = round(mean(clusters[isA] == dA), 3),
               B_in_domA = round(mean(clusters[isB] == dA), 3),
               B_in_domB = round(mean(clusters[isB] == dB), 3),
               stringsAsFactors = FALSE)
}
pair_tab <- rbind(
    pair_row(base_cl,    "1_baseline",  "OE 2", "OE 3"),
    pair_row(rec$cluster,"2_recompute", "OE 2", "OE 3"),
    pair_row(reg$cluster,"3_regress",   "OE 2", "OE 3"),
    pair_row(base_cl,    "1_baseline",  "OE 1", "BE 1"),
    pair_row(rec$cluster,"2_recompute", "OE 1", "BE 1"),
    pair_row(reg$cluster,"3_regress",   "OE 1", "BE 1"))
write.csv(pair_tab, file.path(OUT_DIR, "pair_merge_OE2OE3_OE1BE1.csv"), row.names = FALSE)
cat("\n=== Pair co-clustering (same_cluster = pair merged into one cluster) ===\n")
print(pair_tab, row.names = FALSE)

# ------------------------------------------------------------------
# Full fate table (recompute & regress) + crosstabs
# ------------------------------------------------------------------
ann_levels <- levels(factor(orig))
fate_one <- function(clusters, a) {
    isa <- orig == a; tab <- table(clusters[isa]); dom <- names(tab)[which.max(tab)]
    in_dom <- clusters == dom; other <- orig[in_dom & !isa]
    c(dom = dom, cohesion = round(max(tab)/sum(isa), 3),
      purity = round(sum(isa & in_dom)/sum(in_dom), 3),
      top_other = if (length(other)) names(sort(table(other), TRUE))[1] else NA,
      top_other_frac = if (length(other)) round(max(table(other))/sum(in_dom), 3) else 0)
}
fate <- do.call(rbind, lapply(ann_levels, function(a) {
    r <- fate_one(rec$cluster, a); g <- fate_one(reg$cluster, a)
    data.frame(annotation = a, n = sum(orig == a),
               rec_dom = r["dom"], rec_cohesion = r["cohesion"], rec_purity = r["purity"],
               rec_top_other = r["top_other"], rec_top_other_frac = r["top_other_frac"],
               reg_dom = g["dom"], reg_cohesion = g["cohesion"], reg_purity = g["purity"],
               reg_top_merge = g["top_other"], reg_top_merge_frac = g["top_other_frac"],
               stringsAsFactors = FALSE)
}))
write.csv(fate, file.path(OUT_DIR, "original_annotation_fate.csv"), row.names = FALSE)
cat("\n=== Annotation fate: recompute (no reg) vs regress ===\n"); print(fate, row.names = FALSE)

write.csv(as.data.frame.matrix(table(annotation = orig, recompute = rec$cluster)),
          file.path(OUT_DIR, "crosstab_recompute_vs_annotation.csv"))
write.csv(as.data.frame.matrix(table(annotation = orig, regress = reg$cluster)),
          file.path(OUT_DIR, "crosstab_regress_vs_annotation.csv"))

# ------------------------------------------------------------------
# Figures: baseline (stored umap) vs regress (new umap), by annotation
# ------------------------------------------------------------------
pal <- utils_cb_palette(length(ann_levels)); names(pal) <- ann_levels
mk <- function(emb, title) {
    df <- data.frame(UMAP1 = emb[,1], UMAP2 = emb[,2],
                     annotation = factor(orig, levels = ann_levels))
    ggplot(df, aes(UMAP1, UMAP2, colour = annotation)) + geom_point(size = .3) +
        scale_colour_manual(values = pal) + ggtitle(title) + theme_bw() +
        guides(colour = guide_legend(override.aes = list(size = 3)))
}
ggsave(file.path(OUT_DIR, "UMAP_baseline_by_annotation.png"),
       mk(Embeddings(epi, "umap"), "BASELINE (stored Stage-04) — by annotation"),
       width = 10, height = 8, dpi = 200, bg = "white")
if (!is.null(reg$umap))
    ggsave(file.path(OUT_DIR, "UMAP_regress_by_annotation.png"),
           mk(reg$umap, "REGRESS (04 pipeline + stress regressed) — by annotation"),
           width = 10, height = 8, dpi = 200, bg = "white")

cat("\nDONE — outputs in ", OUT_DIR, "\n")
