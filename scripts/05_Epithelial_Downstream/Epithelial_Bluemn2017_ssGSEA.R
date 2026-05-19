# Epithelial_Bluemn2017_ssGSEA.R
# Bluemn et al. 2017 Cancer Cell 재현 — AR-null neuroendocrine-low CRPC의 FGF/MAPK
# driver phenotype을 epithelial 단일세포 수준에서 ssGSEA로 검증.
# Bluemn 논문의 핵심 finding: AR signaling 소실 + FGFR1/MAPK/ERK/MEK 활성화 +
# EMT 동반 = "double-negative" PCa. 본 cohort는 DNPC late-stage이므로 직접 비교 대상.
#
# Per-cell ssGSEA scoring for 6 Bluemn-signature gene sets.
# Gene sets sourced from MSigDB v6.0 local GMTs (Resources/msigdb_v6.0_GMTs/).
# v6.0 is the last release that retained ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_UP
# with mapped HGNC symbols; later releases dropped it (UniGene deprecation).
#
# Signatures (6, Bluemn 2017 panel):
#   HALLMARK_ANDROGEN_RESPONSE                         -> 예상 DOWN in DNPC clusters
#   HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION         -> UP — plasticity
#   GO_POSITIVE_REGULATION_OF_MAPK_CASCADE             -> UP — RTK downstream
#   GO_REGULATION_OF_ERK1_AND_ERK2_CASCADE             -> UP — ERK MAPK axis
#   ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_UP  -> UP — FGFR1-driven PCa program
#   MEK_UP.V1_UP                                       -> UP — oncogenic MEK activation (C6)
#
# Method: GSVA 2.x ssgseaParam() + gsva(). Per-cell scores written to meta.data.
# One-way ANOVA + TukeyHSD across seurat_clusters per signature.

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(GSVA)
    library(BiocParallel)
})
suppressMessages(source("scripts/00_utils/scRNA_utils.R"))

IN_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Bluemn2017_ssGSEA"
OUT_RDS <- "Results/05_Epithelial_Downstream/epi_bluemn2017_ssgsea.rds"
GMT_DIR <- "Resources/msigdb_v6.0_GMTs"

TARGETS <- c(
    "HALLMARK_ANDROGEN_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "GO_POSITIVE_REGULATION_OF_MAPK_CASCADE",
    "GO_REGULATION_OF_ERK1_AND_ERK2_CASCADE",
    "ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_UP",
    "MEK_UP.V1_UP"
)
GMT_FILES <- c(
    "h.all.v6.0.symbols.gmt",
    "c2.cgp.v6.0.symbols.gmt",
    "c5.bp.v6.0.symbols.gmt",
    "c6.all.v6.0.symbols.gmt"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

# ============================================================
# Load-or-compute: skip heavy ssGSEA when OUT_RDS already exists ----
# ============================================================
if (file.exists(OUT_RDS)) {
    message("Loading cached ", OUT_RDS, " — skipping ssGSEA, regenerating figures")
    epi <- readRDS(OUT_RDS)
    score_cols <- grep("_ssGSEA$", colnames(epi@meta.data), value = TRUE)
} else {

# ============================================================
# Load gene sets from GMT files ----
# ============================================================
read_gmt <- function(path) {
    lines <- readLines(path, warn = FALSE)
    sets <- lapply(lines, function(l) {
        f <- strsplit(l, "\t", fixed = TRUE)[[1]]
        list(name = f[1], genes = f[-(1:2)])
    })
    nms <- vapply(sets, `[[`, character(1), "name")
    out <- lapply(sets, `[[`, "genes")
    names(out) <- nms
    out
}

all_sets <- list()
for (gmt in GMT_FILES) {
    p <- file.path(GMT_DIR, gmt)
    if (!file.exists(p)) stop("Missing GMT: ", p)
    all_sets <- c(all_sets, read_gmt(p))
}

missing <- setdiff(TARGETS, names(all_sets))
if (length(missing) > 0) {
    stop("Sets not found in v6.0 GMTs: ", paste(missing, collapse = ", "))
}
gene_sets <- all_sets[TARGETS]
for (nm in names(gene_sets)) {
    message(sprintf("  %s: %d genes (v6.0)", nm, length(gene_sets[[nm]])))
}

# ============================================================
# Load Seurat object, normalize ----
# ============================================================
epi <- readRDS(IN_RDS)
message("Loaded: ", IN_RDS, " | cells=", ncol(epi), " genes=", nrow(epi))

DefaultAssay(epi) <- "RNA"
if (!"data" %in% Layers(epi[["RNA"]])) {
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)
}
expr <- GetAssayData(epi, assay = "RNA", layer = "data")

# Restrict gene sets to genes present in the expression matrix
present <- rownames(expr)
for (nm in names(gene_sets)) {
    n0 <- length(gene_sets[[nm]])
    gene_sets[[nm]] <- intersect(gene_sets[[nm]], present)
    message(sprintf("  %s: %d / %d genes present", nm, length(gene_sets[[nm]]), n0))
}
keep <- vapply(gene_sets, length, integer(1)) >= 5
if (any(!keep)) {
    warning("Dropping sets with <5 genes: ", paste(names(gene_sets)[!keep], collapse = ", "))
    gene_sets <- gene_sets[keep]
}

# ============================================================
# ssGSEA via GSVA 2.x ----
# ============================================================
# Sparse dgCMatrix -> dense matrix for GSVA. Cells × 6 sets is the output.
expr_dense <- as.matrix(expr)
param <- ssgseaParam(
    exprData   = expr_dense,
    geneSets   = gene_sets,
    minSize    = 5,
    maxSize    = Inf,
    alpha      = 0.25,
    normalize  = TRUE
)

bpp <- tryCatch(MulticoreParam(workers = max(1, parallel::detectCores() - 1)),
                error = function(e) SerialParam())
message("Running ssGSEA on ", ncol(expr_dense), " cells × ",
        length(gene_sets), " sets ...")
t0 <- Sys.time()
ssgsea_mat <- gsva(param, BPPARAM = bpp, verbose = TRUE)
message("ssGSEA done in ", round(as.numeric(Sys.time() - t0, units = "mins"), 2), " min")

# Persist raw matrix (sets × cells)
saveRDS(ssgsea_mat, file.path(OUT_DIR, "ssgsea_matrix.rds"))
write.csv(t(ssgsea_mat),
    file.path(OUT_DIR, "ssgsea_scores_per_cell.csv"),
    row.names = TRUE
)

# Inject scores into meta.data: column name = <SET>_ssGSEA
score_cols <- character(0)
for (nm in rownames(ssgsea_mat)) {
    col <- paste0(nm, "_ssGSEA")
    epi@meta.data[[col]] <- as.numeric(ssgsea_mat[nm, colnames(epi)])
    score_cols <- c(score_cols, col)
}

# Persist annotated object (compute branch only)
saveRDS(epi, OUT_RDS)

}  # end load-or-compute

# Group cells by curated annotation (Epithelial_Annotation.R output) rather
# than raw seurat_clusters, so violin/boxplot/ANOVA panels use the cell-type
# label the user assigned. Falls back to seurat_clusters only if the
# annotation CSV is missing.
ANNOT_CSV    <- "Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv"
ANNOT_LEVELS <- "Results/05_Epithelial_Downstream/Annotation/label_levels.txt"
if (file.exists(ANNOT_CSV) && file.exists(ANNOT_LEVELS)) {
    .a <- read.csv(ANNOT_CSV, stringsAsFactors = FALSE)
    epi$annotation <- .a$annotation[match(colnames(epi), .a$cell)]
    .lvls <- readLines(ANNOT_LEVELS)
    clusters <- factor(epi$annotation, levels = .lvls)
} else {
    warning("annotation files not found — falling back to seurat_clusters")
    clusters <- factor(epi$seurat_clusters)
}

# ============================================================
# UMAP feature plots (viridis, colorblind-safe continuous scale) ----
# ============================================================
paper_grad <- function(score_vec) {
    rng <- range(score_vec, na.rm = TRUE)
    scale_color_viridis_c(option = "viridis", limits = rng)
}
for (sc in score_cols) {
    scores <- epi@meta.data[[sc]]
    p <- FeaturePlot(epi, features = sc, pt.size = 0.3, order = TRUE) +
        paper_grad(scores) +
        ggtitle(sub("_ssGSEA$", "", sc))
    ggsave(file.path(OUT_DIR, paste0("UMAP_", sub("_ssGSEA$", "", sc), ".png")),
        plot = p, width = 20, height = 20, units = "cm", dpi = 200, bg = "white"
    )
}

# ============================================================
# Per-cluster violin + boxplot ----
# ============================================================
for (sc in score_cols) {
    df <- data.frame(cluster = clusters, score = epi@meta.data[[sc]])
    ymax <- max(df$score, na.rm = TRUE)
    title <- sub("_ssGSEA$", "", sc)

    p_v <- ggplot(df, aes(cluster, score, fill = cluster)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = utils_cb_palette(nlevels(clusters))) +
        stat_summary(fun = mean, geom = "point",
                     size = 6, colour = "blue", shape = 95) +
        stat_compare_means(method = "anova",
                           label.x = 3, label.y = ymax + 0.02) +
        ggtitle(title) + NoLegend() +
        theme(text = element_text(size = 11),
              axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, paste0(title, "_violin.png")),
        plot = p_v, width = 20, height = 15, units = "cm", dpi = 200, bg = "white"
    )

    p_b <- ggplot(df, aes(cluster, score, fill = cluster)) +
        geom_boxplot(notch = FALSE, outlier.size = 0.3) +
        scale_fill_manual(values = utils_cb_palette(nlevels(clusters))) +
        stat_summary(fun = mean, geom = "point",
                     size = 6, colour = "blue", shape = 95) +
        ggtitle(title) + NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, paste0(title, "_boxplot.png")),
        plot = p_b, width = 20, height = 15, units = "cm", dpi = 200, bg = "white"
    )
}

# ============================================================
# One-way ANOVA + TukeyHSD per signature ----
# ============================================================
anova_rows <- list()
tukey_rows <- list()
for (sc in score_cols) {
    df <- data.frame(cluster = clusters, score = epi@meta.data[[sc]])
    fit <- aov(score ~ cluster, data = df)
    s <- summary(fit)[[1]]
    anova_rows[[length(anova_rows) + 1]] <- data.frame(
        signature  = sub("_ssGSEA$", "", sc),
        df_between = s["cluster", "Df"],
        df_within  = s["Residuals", "Df"],
        F_stat     = s["cluster", "F value"],
        p_value    = s["cluster", "Pr(>F)"]
    )
    tuk <- TukeyHSD(fit)$cluster
    tuk_df <- as.data.frame(tuk)
    tuk_df$signature  <- sub("_ssGSEA$", "", sc)
    tuk_df$comparison <- rownames(tuk_df)
    rownames(tuk_df) <- NULL
    tukey_rows[[length(tukey_rows) + 1]] <- tuk_df
}
anova_df <- do.call(rbind, anova_rows)
tukey_df <- do.call(rbind, tukey_rows)
tukey_df <- tukey_df[, c("signature", "comparison", "diff", "lwr", "upr", "p adj")]
names(tukey_df)[names(tukey_df) == "p adj"] <- "p_adj"

write.csv(anova_df,
    file.path(OUT_DIR, "OneWay_ANOVA_per_signature.csv"), row.names = FALSE)
write.csv(tukey_df,
    file.path(OUT_DIR, "TukeyHSD_signature_cluster_pairs.csv"), row.names = FALSE)

# Per-cluster mean score summary (compact "what's lit up where")
mean_by_clu <- sapply(score_cols, function(sc)
    tapply(epi@meta.data[[sc]], clusters, mean, na.rm = TRUE))
colnames(mean_by_clu) <- sub("_ssGSEA$", "", colnames(mean_by_clu))
write.csv(round(mean_by_clu, 4),
    file.path(OUT_DIR, "mean_ssGSEA_per_cluster.csv"))

# ============================================================
# Save ----
# ============================================================
# Note: saveRDS(epi, OUT_RDS) happens inside the compute branch above.
message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
