# 04_OE2_Reclustering_RegressOut.R
# QC-artifact diagnostic #4 — STRICTER confirmatory test for OE 2.
#
# Script 03 only dropped 32 stress genes from the PCA feature list; correlated
# genes could still carry the stress signal. Here we regress the dissociation
# stress SCORE out of every gene (Seurat SCTransform(vars.to.regress=...)), which
# removes the stress-correlated component genome-wide, then re-cluster. If OE 2 is
# merely a stress overlay, regressing the program out should dissolve it; if OE 2
# survives this too, its identity is genuinely stress-independent.
#
# Controlled design (only vars.to.regress differs):
#   (A) CONTROL : SCTransform(epi)                       -> reproduce OE 2
#   (B) TEST    : SCTransform(epi, vars.to.regress=stress) -> does OE 2 dissolve?
# Manipulation check: R^2 of (stress_score ~ harmony 1:30) — high in control,
# should drop in test if the regression actually removed the stress axis.
#
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds   (read-only)
# Writes: Results/06_QC_Artifact_Diagnostics/04_OE2_Reclustering_RegressOut/
#           cell_clusterings.csv / OE2_fate_summary.csv
#           crosstab_{control,test}_vs_annotation.csv
#           regression_effectiveness.csv
#           UMAP_{control,test}_by_annotation.png / UMAP_test_OE2_highlight.png
# Existing Results/04 and Results/05 objects are NOT modified.

suppressMessages({
    library(Seurat); library(dplyr); library(ggplot2)
    library(patchwork); library(harmony); library(future)
})
plan("sequential")
options(future.globals.maxSize = 128 * 1024^3)
source("scripts/00_utils/scRNA_utils.R")

set.seed(42)
IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/06_QC_Artifact_Diagnostics/04_OE2_Reclustering_RegressOut"
FOCUS   <- "OE 2"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Broad dissociation/stress panel -> single module score that gets regressed out.
stress_genes <- unique(c(
    "FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1","HSPA1A","HSPA1B","HSP90AB1",
    "HSPA8","HSPB1","IER3","IER2","BTG1","BTG2","DUSP1","FOSL1","FOSL2","EGR2",
    "EGR3","IER5","IER5L","NR4A1","NR4A2","NR4A3","DUSP2","DUSP5","ZFP36",
    "GADD45B","GADD45G","SOCS3","CCN1","KLF2","KLF4","KLF6","MCL1","RHOB","MAFB",
    "HSPA6","HSP90AA1","DNAJB1","DNAJA1","BAG3","HSPH1"))

stopifnot(file.exists(IN_RDS))
epi <- readRDS(IN_RDS)
message("Loaded ", IN_RDS, " (", ncol(epi), " cells)")
orig_annotation <- epi$annotation

# ------------------------------------------------------------------
# Stress module score on RNA (the regression covariate)
# ------------------------------------------------------------------
DefaultAssay(epi) <- "RNA"
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
epi <- NormalizeData(epi, verbose = FALSE)
set.seed(42)
sg <- intersect(stress_genes, rownames(epi[["RNA"]]))
epi <- AddModuleScore(epi, features = list(sg), name = "stress_sc", assay = "RNA")
epi$stress_score <- epi$stress_sc1; epi$stress_sc1 <- NULL
message(sprintf("stress_score from %d genes; range %.2f..%.2f",
                length(sg), min(epi$stress_score), max(epi$stress_score)))

# ------------------------------------------------------------------
# Pipeline helper: PCA -> Harmony -> clusters(0.4) -> UMAP
# ------------------------------------------------------------------
run_harmony <- function(obj, dims = 1:30) {
    set.seed(42)
    tryCatch(
        RunHarmony(obj, group.by.vars = "orig.ident", reduction.use = "pca",
                   dims.use = dims, reduction.save = "harmony_tmp", verbose = FALSE),
        error = function(e)
            RunHarmony(obj, group.by.vars = "orig.ident", reduction = "pca",
                       dims = dims, reduction.save = "harmony_tmp", verbose = FALSE))
}
run_pipeline <- function(obj, tag, dims = 1:30) {
    DefaultAssay(obj) <- "SCT"
    set.seed(42)
    obj <- RunPCA(obj, features = VariableFeatures(obj, assay = "SCT"),
                  npcs = 50, verbose = FALSE)
    obj <- run_harmony(obj, dims = dims)
    obj <- FindNeighbors(obj, reduction = "harmony_tmp", dims = dims, verbose = FALSE)
    set.seed(42)
    obj <- FindClusters(obj, resolution = 0.4, verbose = FALSE)
    obj <- RunUMAP(obj, reduction = "harmony_tmp", dims = dims, n.neighbors = 40,
                   min.dist = 0.3, spread = 1.5,
                   reduction.name = paste0("umap_", tag), verbose = FALSE)
    H <- Embeddings(obj, "harmony_tmp")[, dims]
    r2 <- summary(lm(obj$stress_score ~ H))$r.squared     # manipulation check
    list(cluster = as.character(obj$seurat_clusters),
         umap = Embeddings(obj, paste0("umap_", tag)), stress_r2 = r2)
}

# CONTROL: plain SCT
message("CONTROL: SCTransform (no regression) ...")
epi <- SCTransform(epi, verbose = FALSE)
ctrl <- run_pipeline(epi, "control")

# TEST: SCT with stress score regressed out genome-wide
message("TEST: SCTransform(vars.to.regress = stress_score) ...")
epi <- SCTransform(epi, vars.to.regress = "stress_score", verbose = FALSE)
test <- run_pipeline(epi, "test")

clust_df <- data.frame(
    cell = colnames(epi), annotation = as.character(orig_annotation),
    control_cluster = ctrl$cluster, test_cluster = test$cluster,
    stringsAsFactors = FALSE)
write.csv(clust_df, file.path(OUT_DIR, "cell_clusterings.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Manipulation check + cross-tabs
# ------------------------------------------------------------------
reg_eff <- data.frame(run = c("control","test"),
                      stress_var_explained_by_harmony_R2 =
                          round(c(ctrl$stress_r2, test$stress_r2), 3))
write.csv(reg_eff, file.path(OUT_DIR, "regression_effectiveness.csv"), row.names = FALSE)
cat("\n=== Manipulation check: R^2 of stress_score ~ harmony(1:30) ===\n")
print(reg_eff, row.names = FALSE)

ct_ctrl <- table(annotation = clust_df$annotation, control = clust_df$control_cluster)
ct_test <- table(annotation = clust_df$annotation, test = clust_df$test_cluster)
write.csv(as.data.frame.matrix(ct_ctrl), file.path(OUT_DIR, "crosstab_control_vs_annotation.csv"))
write.csv(as.data.frame.matrix(ct_test), file.path(OUT_DIR, "crosstab_test_vs_annotation.csv"))

# ------------------------------------------------------------------
# OE 2 fate (same metrics as script 03)
# ------------------------------------------------------------------
fate <- function(new_cluster) {
    is_f <- clust_df$annotation == FOCUS
    tab  <- table(new_cluster[is_f]); dom <- names(tab)[which.max(tab)]
    in_dom <- new_cluster == dom
    other  <- clust_df$annotation[in_dom & !is_f]
    data.frame(dom_cluster = dom, n_focus_cells = sum(is_f),
               cohesion = round(max(tab)/sum(is_f), 3),
               dom_cluster_purity_focus = round(sum(is_f & in_dom)/sum(in_dom), 3),
               dom_cluster_size = sum(in_dom),
               absorbed_into = if (length(other)) names(sort(table(other), TRUE))[1] else NA,
               absorbed_frac = if (length(other)) round(max(table(other))/sum(in_dom), 3) else 0)
}
fate_tab <- rbind(cbind(run = "control", fate(clust_df$control_cluster)),
                  cbind(run = "test",    fate(clust_df$test_cluster)))
write.csv(fate_tab, file.path(OUT_DIR, "OE2_fate_summary.csv"), row.names = FALSE)
cat("\n=== OE 2 fate (control = no regression; test = stress regressed out) ===\n")
print(fate_tab, row.names = FALSE)

ari <- function(a, b) {
    tab <- table(a, b); n <- sum(tab); s_ij <- sum(choose(tab, 2))
    a_i <- sum(choose(rowSums(tab), 2)); b_j <- sum(choose(colSums(tab), 2))
    exp <- a_i*b_j/choose(n,2); (s_ij - exp)/((a_i+b_j)/2 - exp)
}
cat(sprintf("\nARI(control, test) overall = %.3f\n", ari(clust_df$control_cluster, clust_df$test_cluster)))

# ------------------------------------------------------------------
# Figures
# ------------------------------------------------------------------
ann_levels <- levels(orig_annotation)
pal <- utils_cb_palette(length(ann_levels)); names(pal) <- ann_levels
mk_umap <- function(emb, title) {
    df <- data.frame(UMAP1 = emb[,1], UMAP2 = emb[,2],
                     annotation = factor(clust_df$annotation, levels = ann_levels))
    ggplot(df, aes(UMAP1, UMAP2, colour = annotation)) +
        geom_point(size = .3) + scale_colour_manual(values = pal) +
        ggtitle(title) + theme_bw() +
        guides(colour = guide_legend(override.aes = list(size = 3)))
}
ggsave(file.path(OUT_DIR, "UMAP_control_by_annotation.png"),
       mk_umap(ctrl$umap, "CONTROL (no regression) — by original annotation"),
       width = 10, height = 8, dpi = 200, bg = "white")
ggsave(file.path(OUT_DIR, "UMAP_test_by_annotation.png"),
       mk_umap(test$umap, "TEST (stress score regressed out) — by original annotation"),
       width = 10, height = 8, dpi = 200, bg = "white")
df_h <- data.frame(UMAP1 = test$umap[,1], UMAP2 = test$umap[,2],
                   focus = clust_df$annotation == FOCUS)
ggsave(file.path(OUT_DIR, "UMAP_test_OE2_highlight.png"),
       ggplot() +
           geom_point(data = subset(df_h, !focus), aes(UMAP1, UMAP2), colour = "grey80", size = .3) +
           geom_point(data = subset(df_h, focus), aes(UMAP1, UMAP2), colour = "#D55E00", size = .5) +
           ggtitle(paste0("TEST (stress regressed) — original ", FOCUS, " cells (orange)")) +
           theme_bw(),
       width = 10, height = 8, dpi = 200, bg = "white")

cat("\nDONE — outputs in ", OUT_DIR, "\n")
