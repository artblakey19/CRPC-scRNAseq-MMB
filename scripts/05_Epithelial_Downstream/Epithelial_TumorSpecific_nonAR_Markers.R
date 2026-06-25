# Epithelial_TumorSpecific_nonAR_Markers.R
# Question: are the "ARPC" cells (cluster 9) actually malignant, or just
#   AR-pathway-active luminal cells (which benign prostate luminal also are)?
#
# AR-response / PCa-AR signatures CANNOT separate benign vs malignant luminal —
# both are AR+. This script contrasts:
#   (A) canonical AR-luminal markers (magnitude-driven; benign AND tumor share)
#   (B) TUMOR-SPECIFIC, non-androgen markers that are aberrant in PCa but ~absent
#       in benign luminal (AMACR, ERG/ETS fusion readout, SPINK1, EZH2, MYC,
#       GOLM1, lncRNAs PCA3/SCHLAP1)
#   (C) down-in-cancer contrast markers (MSMB/PSP94, GSTP1) — benign-enriched
#   (D) Song et al. 2022 PCa-DEG tumor signatures (ERGpos_Tumor / ERGneg_Tumor)
#       vs Normal LE, scored with AddModuleScore.
#
# Logic: if ARPC is distinctly HIGH on (B)/(D) and LOW on (C) relative to the
# other luminal-ish clusters, that SUPPORTS a malignant call beyond mere AR
# activity. If ARPC is indistinguishable on (B)/(C)/(D), the malignant call
# rests only on AR activity + clinical context (not provable from transcriptome).
#
# CNV-free by design: copyKAT/inferCNV deliberately abandoned in this project;
# this is an orthogonal, expression-only line of evidence.
#
# Caveats (apply when interpreting):
#   - 10x 3' chemistry under-captures lncRNAs (PCA3, SCHLAP1): a zero is NOT
#     evidence of absence.
#   - TMPRSS2-ERG fusion is a patient-level event -> ERG is reported per patient
#     within ARPC, not just pooled.
#   - AMACR is the single most robust non-AR PCa marker here; benign luminal is
#     AMACR-low but not strictly zero.

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(openxlsx)
})

IN_RDS       <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
ANNOT_CSV    <- "Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv"
ANNOT_LEVELS <- "Results/05_Epithelial_Downstream/Annotation/label_levels.txt"
SIG_XLSX     <- "Resources/Song2022_SuppData2.xlsx"
OUT_DIR      <- "Results/05_Epithelial_Downstream/TumorSpecific_nonAR_Markers"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

stopifnot("Input rds not found" = file.exists(IN_RDS))

# ------------------------------------------------------------------
# Marker panels --------------------------------------------------
# ------------------------------------------------------------------
panels <- list(
    "A_AR_luminal"      = c("KLK3", "KLK2", "KLK4", "NKX3-1", "AR",
                            "FOLH1", "TMPRSS2", "ACPP", "STEAP2"),
    "B_tumor_nonAR_up"  = c("AMACR", "ERG", "ETV1", "ETV4", "SPINK1",
                            "EZH2", "MYC", "GOLM1", "SCHLAP1", "PCA3"),
    "C_benign_down"     = c("MSMB", "GSTP1"),
    "D_proliferation"   = c("MKI67", "TOP2A")
)
all_genes <- unique(unlist(panels))

# ------------------------------------------------------------------
# Load object + annotation --------------------------------------
# ------------------------------------------------------------------
epi <- readRDS(IN_RDS)
message("Loaded: ", IN_RDS, " | cells=", ncol(epi))

DefaultAssay(epi) <- "RNA"
if (!"data" %in% Layers(epi[["RNA"]])) {
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)
}

stopifnot("annotation CSV missing"    = file.exists(ANNOT_CSV),
          "annotation levels missing" = file.exists(ANNOT_LEVELS))
.a <- read.csv(ANNOT_CSV, stringsAsFactors = FALSE)
epi$annotation <- .a$annotation[match(colnames(epi), .a$cell)]
lvls <- readLines(ANNOT_LEVELS)
epi$annotation <- factor(epi$annotation, levels = lvls)
clu <- epi$annotation

# ------------------------------------------------------------------
# Per-cluster mean (normalized) + pct expressing ----------------
# ------------------------------------------------------------------
present <- intersect(all_genes, rownames(epi))
missing <- setdiff(all_genes, present)
if (length(missing) > 0) message("Genes NOT in object (cannot assess): ",
                                 paste(missing, collapse = ", "))

mat <- as.matrix(GetAssayData(epi, assay = "RNA", layer = "data")[present, , drop = FALSE])

mean_by <- t(apply(mat, 1, function(x) tapply(x, clu, mean)))
pct_by  <- t(apply(mat, 1, function(x) tapply(x > 0, clu, function(z) 100 * mean(z))))
mean_by <- mean_by[, lvls, drop = FALSE]
pct_by  <- pct_by[,  lvls, drop = FALSE]

panel_of <- setNames(rep(names(panels), lengths(panels)), unlist(panels))
mean_df <- data.frame(panel = panel_of[rownames(mean_by)], gene = rownames(mean_by),
                      round(mean_by, 4), check.names = FALSE)
pct_df  <- data.frame(panel = panel_of[rownames(pct_by)],  gene = rownames(pct_by),
                      round(pct_by, 1),  check.names = FALSE)
write.csv(mean_df, file.path(OUT_DIR, "mean_expr_per_cluster.csv"), row.names = FALSE)
write.csv(pct_df,  file.path(OUT_DIR, "pct_expressing_per_cluster.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# Dotplot: gene x cluster, size = pct, color = mean expr ---------
# ------------------------------------------------------------------
long <- merge(
    as.data.frame(as.table(mean_by)) |> setNames(c("gene", "cluster", "mean")),
    as.data.frame(as.table(pct_by))  |> setNames(c("gene", "cluster", "pct")),
    by = c("gene", "cluster")
)
long$panel   <- panel_of[as.character(long$gene)]
long$gene    <- factor(long$gene, levels = rev(unlist(panels)[unlist(panels) %in% present]))
long$cluster <- factor(long$cluster, levels = lvls)

p <- ggplot(long, aes(cluster, gene)) +
    geom_point(aes(size = pct, color = mean)) +
    scale_color_viridis_c(option = "viridis", name = "mean\nexpr") +
    scale_size_continuous(range = c(0, 7), name = "% expr") +
    facet_grid(panel ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Tumor-specific (non-AR) vs AR-luminal markers by cluster",
         x = NULL, y = NULL) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          strip.text.y = element_text(angle = 0))
ggsave(file.path(OUT_DIR, "DotPlot_TumorSpecific_nonAR.png"),
       plot = p, width = 9, height = 9, dpi = 200, bg = "white")

# ------------------------------------------------------------------
# ERG by patient within each cluster (fusion is patient-level) ---
# ------------------------------------------------------------------
if ("ERG" %in% present) {
    erg <- mat["ERG", ]
    erg_df <- data.frame(cluster = clu, patient = epi$orig.ident, erg = erg) |>
        group_by(cluster, patient) |>
        summarise(n = n(), mean_ERG = round(mean(erg), 4),
                  pct_ERG = round(100 * mean(erg > 0), 1), .groups = "drop")
    write.csv(erg_df, file.path(OUT_DIR, "ERG_by_cluster_patient.csv"), row.names = FALSE)
}

# ------------------------------------------------------------------
# (D) Song 2022 PCa-DEG tumor signatures via AddModuleScore ------
# ------------------------------------------------------------------
if (file.exists(SIG_XLSX)) {
    pca_sig <- tryCatch(read.xlsx(SIG_XLSX, sheet = "PCa signature"),
                        error = function(e) NULL)
    if (!is.null(pca_sig)) {
        wanted <- grep("ERGpos|ERGneg|^LE$|luminal", colnames(pca_sig),
                       ignore.case = TRUE, value = TRUE)
        sig_list <- list()
        for (cn in wanted) {
            g <- unique(na.omit(as.character(pca_sig[[cn]])))
            g <- intersect(g[nzchar(g)], rownames(epi))
            if (length(g) >= 5) sig_list[[cn]] <- g
        }
        message("Song PCa-signature sets scored: ",
                paste(names(sig_list), collapse = ", "))
        score_cols <- character(0)
        for (cn in names(sig_list)) {
            nm <- paste0("Song_", gsub("[^A-Za-z0-9]", "_", cn))
            epi <- AddModuleScore(epi, features = list(sig_list[[cn]]),
                                  name = nm, assay = "RNA")
            epi@meta.data[[nm]] <- epi@meta.data[[paste0(nm, "1")]]
            epi@meta.data[[paste0(nm, "1")]] <- NULL
            score_cols <- c(score_cols, nm)
        }
        if (length(score_cols) > 0) {
            song_mean <- sapply(score_cols, function(sc)
                tapply(epi@meta.data[[sc]], clu, mean))[lvls, , drop = FALSE]
            write.csv(round(song_mean, 4),
                      file.path(OUT_DIR, "Song_PCa_tumor_signature_mean_per_cluster.csv"))
        }
    }
}

message("Done. Outputs in: ", OUT_DIR)
