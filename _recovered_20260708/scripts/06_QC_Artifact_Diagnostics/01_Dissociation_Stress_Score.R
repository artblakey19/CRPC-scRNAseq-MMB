# 01_Dissociation_Stress_Score.R
# QC-artifact diagnostic #1 — is the OE 2 cluster a tissue-dissociation / AP-1·IEG
# stress technical state rather than a biological cell type?
#
# Rationale (literature, deep-research verified 2026-06-29):
#   - van den Brink 2017 (Nat Methods 14:935) PROVED warm (37C) enzymatic
#     dissociation induces an artifactual immediate-early / AP-1 stress program
#     (FOS/JUN/... + heat-shock); FOS 0/80 in intact tissue vs 36% in dissociated.
#     Recommended handling: SCORE cells, then flag/filter or regress — do NOT
#     annotate the cluster as a cell type.
#   - Denisenko 2020 (Genome Biol) per-cell "stress score" = 17 genes via Seurat
#     AddModuleScore. <-- reproduced exactly below.
# OE 2 top markers (FOS/JUN/EGR1/IER2/IER5L/DDIT3/GADD45G) match this program.
#
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds   (annotation + QC, untouched)
# Writes: Results/06_QC_Artifact_Diagnostics/01_Dissociation_Stress_Score/
#           per_cluster_stress_scores.csv      mean/median module score per cluster
#           ANOVA_summary.csv                   F / eta^2 per signature
#           ANOVA_TukeyHSD_OE2_vs_rest.csv      stats: is OE 2 the significant outlier?
#           VlnPlot_Denisenko17_byCluster.png
#           UMAP_Denisenko17_stress.png         (viridis continuous)
#
# Source RDS left untouched. Diagnostic only — applies no filtering.

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
})
source("scripts/00_utils/scRNA_utils.R")

set.seed(42)  # AddModuleScore samples control genes -> seed for reproducibility
              # (see project note: reproducibility_harmony_seed)

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/06_QC_Artifact_Diagnostics/01_Dissociation_Stress_Score"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ------------------------------------------------------------------
# Denisenko 2020 EXACT 17-gene dissociation stress score (their Methods).
# https://doi.org/10.1186/s13059-020-02048-6  — scored via Seurat AddModuleScore.
# ------------------------------------------------------------------
denisenko17 <- c("FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1",
                 "HSPA1A","HSPA1B","HSP90AB1","HSPA8","HSPB1",
                 "IER3","IER2","BTG1","BTG2","DUSP1")

signatures <- list(Denisenko17 = denisenko17)

# ------------------------------------------------------------------
# Load annotated epithelial object
# ------------------------------------------------------------------
stopifnot("Input rds not found" = file.exists(IN_RDS))
epi <- readRDS(IN_RDS)
message("Loaded: ", IN_RDS, "  (", ncol(epi), " cells)")

DefaultAssay(epi) <- "RNA"
if (!"data" %in% Layers(epi[["RNA"]])) {
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)
}
Idents(epi) <- "annotation"
ann_levels <- levels(epi$annotation)

# Restrict each signature to genes actually present; report coverage.
present <- rownames(epi[["RNA"]])
for (s in names(signatures)) {
    g  <- intersect(signatures[[s]], present)
    miss <- setdiff(signatures[[s]], present)
    message(sprintf("  %-12s %d/%d genes present%s", s, length(g),
                    length(signatures[[s]]),
                    if (length(miss)) paste0(" (missing: ", paste(miss, collapse=","), ")") else ""))
    signatures[[s]] <- g
}

# ------------------------------------------------------------------
# AddModuleScore per signature
# ------------------------------------------------------------------
set.seed(42)
for (s in names(signatures)) {
    epi <- tryCatch(
        AddModuleScore(epi, features = list(signatures[[s]]),
                       name = paste0(s, "_sc"), assay = "RNA"),
        error = function(e)
            AddModuleScore(epi, features = list(signatures[[s]]),
                           name = paste0(s, "_sc"), nbin = 12, assay = "RNA"))
    # AddModuleScore appends a trailing "1"; normalise the column name.
    epi@meta.data[[paste0(s, "_sc")]] <- epi@meta.data[[paste0(s, "_sc1")]]
    epi@meta.data[[paste0(s, "_sc1")]] <- NULL
}
score_cols <- paste0(names(signatures), "_sc")

# ------------------------------------------------------------------
# (A) Per-cluster score summary
# ------------------------------------------------------------------
md <- epi@meta.data
summ <- md %>%
    group_by(annotation) %>%
    summarise(n = n(),
              across(all_of(score_cols),
                     list(mean = ~round(mean(.x), 4), median = ~round(median(.x), 4)),
                     .names = "{.col}_{.fn}"),
              .groups = "drop") %>%
    arrange(desc(.data[[paste0(score_cols[1], "_mean")]]))
write.csv(summ, file.path(OUT_DIR, "per_cluster_stress_scores.csv"), row.names = FALSE)
cat("\n=== Per-cluster stress scores (sorted by Denisenko17 mean) ===\n")
print(as.data.frame(summ), row.names = FALSE)

# ------------------------------------------------------------------
# (B) ANOVA + TukeyHSD: is OE 2 the significant outlier?
# ------------------------------------------------------------------
aov_rows <- list(); tuk_rows <- list()
for (sc in score_cols) {
    df  <- data.frame(score = md[[sc]], grp = md$annotation)
    fit <- aov(score ~ grp, data = df)
    fs  <- summary(fit)[[1]]
    ss  <- fs[["Sum Sq"]]; eta2 <- ss[1] / sum(ss)         # effect size
    aov_rows[[sc]] <- data.frame(signature = sc,
        F = round(fs[["F value"]][1], 1),
        p = signif(fs[["Pr(>F)"]][1], 3),
        eta2 = round(eta2, 3),
        grp_mean_rank1 = names(sort(tapply(df$score, df$grp, mean), TRUE))[1],
        grp_mean_rank2 = names(sort(tapply(df$score, df$grp, mean), TRUE))[2])
    # OE 2 vs every other cluster
    tu <- as.data.frame(TukeyHSD(fit)$grp)
    tu$pair <- rownames(tu)
    oe2 <- tu[grepl("(^|-)OE 2($|-)", tu$pair), ]
    oe2$signature <- sc
    tuk_rows[[sc]] <- oe2[, c("signature","pair","diff","lwr","upr","p adj")]
}
aov_tab <- do.call(rbind, aov_rows)
tuk_tab <- do.call(rbind, tuk_rows)
write.csv(rbind(
    cbind(section = "ANOVA", aov_tab[, 1:4],
          detail = paste0("top1=", aov_tab$grp_mean_rank1, "; top2=", aov_tab$grp_mean_rank2),
          stringsAsFactors = FALSE) |> setNames(c("section","signature","stat_F","stat_p","stat_eta2","detail")),
    NULL), file.path(OUT_DIR, "ANOVA_summary.csv"), row.names = FALSE)
write.csv(tuk_tab, file.path(OUT_DIR, "ANOVA_TukeyHSD_OE2_vs_rest.csv"), row.names = FALSE)
cat("\n=== ANOVA (score ~ cluster) ===\n"); print(aov_tab, row.names = FALSE)
# How many of the 11 OE2-vs-other contrasts are positive (OE2 higher) & significant?
oe2_den <- tuk_tab[tuk_tab$signature == "Denisenko17_sc", ]
# diff sign depends on pair ordering "A-B"; normalise so positive = OE2 higher
oe2_den$oe2_higher <- ifelse(grepl("^OE 2-", oe2_den$pair), oe2_den$diff > 0, oe2_den$diff < 0)
cat(sprintf("\nDenisenko17: OE 2 higher than the other cluster in %d/%d contrasts; significant (padj<0.05) in %d/%d.\n",
            sum(oe2_den$oe2_higher), nrow(oe2_den),
            sum(oe2_den$oe2_higher & oe2_den$`p adj` < 0.05), nrow(oe2_den)))

# ------------------------------------------------------------------
# (C) Plots
# ------------------------------------------------------------------
pal <- utils_cb_palette(length(ann_levels))

p_den <- VlnPlot(epi, features = "Denisenko17_sc", group.by = "annotation",
                 cols = pal, pt.size = 0) +
    ggtitle("Denisenko 2020 dissociation stress score (17 genes)") +
    theme(legend.position = "none", axis.title.x = element_blank())
ggsave(file.path(OUT_DIR, "VlnPlot_Denisenko17_byCluster.png"),
       p_den, width = 11, height = 6, dpi = 200, bg = "white")

p_umap <- FeaturePlot(epi, features = "Denisenko17_sc", pt.size = 0.3, order = TRUE) +
    scale_color_viridis_c(option = "viridis") +
    ggtitle("Denisenko17 dissociation stress score")
ggsave(file.path(OUT_DIR, "UMAP_Denisenko17_stress.png"),
       p_umap, width = 10, height = 8, dpi = 200, bg = "white")

cat("\nDONE — outputs in ", OUT_DIR, "\n")
