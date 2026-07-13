# ============================================================
# Epithelial_Signature_Scoring_Normal_vs_PCa.R
# Normal Signature vs PCa signature를 구분해서 AddModuleScore + plot
# ============================================================

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(openxlsx)
library(ggpubr)
library(msigdbr)
library(patchwork)
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
SIG_XLSX <- "Resources/Song2022_SuppData2.xlsx"

OUT_DIR <- "Results/05_Epithelial_Downstream/Song2022_signature_scores"
OUT_RDS <- "Results/05_Epithelial_Downstream/Song2022_signature_scores.rds"

FORCE_RECOMPUTE <- FALSE
ADD_AR_SIGNATURE <- TRUE

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

SIGNATURE_SHEETS <- c(
    Normal = "Normal Signature",
    PCa    = "PCa signature"
)

# ============================================================
# Signature parser
# ============================================================

canonical_match <- function(nm) {
    nm_l <- tolower(nm)

    if (grepl("basal|^be$", nm_l)) return("BE")
    if (grepl("luminal|^le$", nm_l)) return("LE")
    if (grepl("hillock", nm_l)) return("Hillock")
    if (grepl("club", nm_l)) return("Club")
    if (grepl("neuroendocrine|^ne$", nm_l)) return("NE")

    # PCa signature sheet
    if (grepl("ergpos|erg_pos|erg\\+|ergpositive", nm_l)) return("ERGpos_Tumor")
    if (grepl("ergneg|erg_neg|erg\\-|ergnegative", nm_l)) return("ERGneg_Tumor")
    if (grepl("^tumou?r$", nm_l)) return("Tumor")

    NA_character_
}

read_sigs_from_sheet <- function(path, sheet,
                                 wanted = c(
                                     "BE", "LE", "Hillock", "Club", "NE",
                                     "ERGpos_Tumor", "ERGneg_Tumor", "Tumor"
                                 )) {
    df <- read.xlsx(path, sheet = sheet)

    out <- list()

    for (col_nm in colnames(df)) {
        canon <- canonical_match(col_nm)

        if (is.na(canon) || !canon %in% wanted) next

        genes <- unique(na.omit(as.character(df[[col_nm]])))
        genes <- genes[nzchar(genes)]

        if (length(genes) > 0) {
            out[[canon]] <- genes
        }
    }

    out
}

read_all_sigs <- function(path, sheets) {
    all_sigs <- list()

    for (src in names(sheets)) {
        sheet <- sheets[[src]]
        sigs <- read_sigs_from_sheet(path, sheet)

        if (length(sigs) == 0) {
            warning("No signatures parsed from sheet: ", sheet)
            next
        }

        message("Loaded ", src, " signatures from sheet: ", sheet)
        message("  ", paste(names(sigs), collapse = ", "))

        for (sig_nm in names(sigs)) {
            # column base name: Normal_BE, PCa_Club, PCa_ERGpos_Tumor ...
            all_sigs[[paste(src, sig_nm, sep = "_")]] <- sigs[[sig_nm]]
        }
    }

    all_sigs
}

# ============================================================
# Load or compute
# ============================================================

if (file.exists(OUT_RDS) && !FORCE_RECOMPUTE) {

    message("Loading cached object: ", OUT_RDS)
    epi <- readRDS(OUT_RDS)

} else {

    stopifnot(
        "epi_filtered_clustered.rds not found — run Epithelial_Analysis.R first" =
            file.exists(EPI_RDS),
        "Song2022_SuppData2.xlsx not found" =
            file.exists(SIG_XLSX)
    )

    epi <- readRDS(EPI_RDS)

    DefaultAssay(epi) <- "RNA"
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)

    sigs <- read_all_sigs(SIG_XLSX, SIGNATURE_SHEETS)

    # Optional: Hallmark AR score. Normal/PCa 소속이 아니라 공통 reference score로 둠.
    if (ADD_AR_SIGNATURE) {
        hallmark_df <- tryCatch(
            msigdbr(species = "Homo sapiens", collection = "H"),
            error = function(e) msigdbr(species = "Homo sapiens", category = "H")
        )

        hallmark_ar <- unique(
            hallmark_df$gene_symbol[
                hallmark_df$gs_name == "HALLMARK_ANDROGEN_RESPONSE"
            ]
        )

        if (length(hallmark_ar) > 0) {
            sigs[["Hallmark_AR"]] <- hallmark_ar
            message("Loaded Hallmark_AR: ", length(hallmark_ar), " genes")
        } else {
            warning("HALLMARK_ANDROGEN_RESPONSE not found in msigdbr")
        }
    }

    # Keep genes present in object
    present <- rownames(epi)

    for (s in names(sigs)) {
        n0 <- length(sigs[[s]])
        sigs[[s]] <- intersect(sigs[[s]], present)
        message(sprintf("  %s: %d / %d genes present", s, length(sigs[[s]]), n0))
    }

    # ========================================================
    # AddModuleScore
    # ========================================================

    for (s in names(sigs)) {

        if (length(sigs[[s]]) < 5) {
            message("Skip ", s, " — fewer than 5 valid genes")
            next
        }

        score_base <- paste0(s, "_Score")

        epi <- tryCatch(
            AddModuleScore(
                object = epi,
                features = list(sigs[[s]]),
                name = score_base,
                assay = "RNA"
            ),
            error = function(e) {
                message("nbin=24 failed for ", s, " — retry with nbin=12: ", e$message)

                AddModuleScore(
                    object = epi,
                    features = list(sigs[[s]]),
                    name = score_base,
                    nbin = 12,
                    assay = "RNA"
                )
            }
        )

        src_col <- paste0(score_base, "1")
        dst_col <- score_base

        epi@meta.data[[dst_col]] <- epi@meta.data[[src_col]]
        epi@meta.data[[src_col]] <- NULL
    }

    saveRDS(epi, OUT_RDS)
}

# ============================================================
# Grouping: curated annotation first, otherwise seurat_clusters
# ============================================================

ANNOT_CSV    <- "Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv"
ANNOT_LEVELS <- "Results/05_Epithelial_Downstream/Annotation/label_levels.txt"

if (file.exists(ANNOT_CSV) && file.exists(ANNOT_LEVELS)) {

    ann <- read.csv(ANNOT_CSV, stringsAsFactors = FALSE)

    epi$annotation <- ann$annotation[match(colnames(epi), ann$cell)]

    if (anyNA(epi$annotation)) {
        warning(sum(is.na(epi$annotation)), " cells have no annotation")
    }

    lvls <- readLines(ANNOT_LEVELS)
    epi$plot_group <- droplevels(factor(epi$annotation, levels = lvls))

} else {

    warning("annotation files not found — using seurat_clusters")
    epi$plot_group <- factor(epi$seurat_clusters)
}

# ============================================================
# Score columns
# ============================================================

normal_score_cols <- grep("^Normal_.*_Score$", colnames(epi@meta.data), value = TRUE)
pca_score_cols    <- grep("^PCa_.*_Score$", colnames(epi@meta.data), value = TRUE)
ar_score_cols     <- grep("^Hallmark_AR_Score$", colnames(epi@meta.data), value = TRUE)

score_cols <- c(normal_score_cols, pca_score_cols, ar_score_cols)

message("Normal score columns: ", paste(normal_score_cols, collapse = ", "))
message("PCa score columns: ", paste(pca_score_cols, collapse = ", "))
message("Common score columns: ", paste(ar_score_cols, collapse = ", "))

# ============================================================
# Helpers
# ============================================================

score_meta <- function(sc) {
    source <- sub("^([^_]+)_.*_Score$", "\\1", sc)
    signature <- sub("^[^_]+_(.*)_Score$", "\\1", sc)

    if (grepl("^Hallmark_", sc)) {
        source <- "Hallmark"
        signature <- sub("^Hallmark_(.*)_Score$", "\\1", sc)
    }

    list(source = source, signature = signature)
}

score_range <- function(x) {
    c(min(x, na.rm = TRUE), max(x, na.rm = TRUE))
}

make_umap <- function(sc, out_dir) {
    meta <- score_meta(sc)
    x <- epi@meta.data[[sc]]

    p <- FeaturePlot(
        epi,
        features = sc,
        pt.size = 0.3,
        order = TRUE
    ) +
        scale_color_viridis_c(
            option = "viridis",
            limits = score_range(x)
        ) +
        ggtitle(paste(meta$source, meta$signature, sep = " - "))

    ggsave(
        file.path(out_dir, paste0("UMAP_", sc, ".png")),
        plot = p,
        width = 20,
        height = 20,
        units = "cm",
        dpi = 200,
        bg = "white"
    )
}

make_split_umap <- function(sc, out_dir) {
    meta <- score_meta(sc)
    x <- epi@meta.data[[sc]]
    n_samples <- length(unique(epi$orig.ident))

    p <- FeaturePlot(
        epi,
        features = sc,
        split.by = "orig.ident",
        pt.size = 0.3,
        order = TRUE
    ) &
        scale_color_viridis_c(
            option = "viridis",
            limits = score_range(x)
        )

    p <- p + plot_annotation(
        title = paste(meta$source, meta$signature, "split by sample", sep = " - ")
    )

    ggsave(
        file.path(out_dir, paste0("UMAP_", sc, "_splitBySample.png")),
        plot = p,
        width = 8 * n_samples,
        height = 8,
        dpi = 200,
        bg = "white"
    )
}

make_boxplot <- function(sc, out_dir) {
    meta <- score_meta(sc)

    df <- data.frame(
        group = epi$plot_group,
        score = epi@meta.data[[sc]]
    )

    p <- ggplot(df, aes(group, score, fill = group)) +
        geom_boxplot(notch = FALSE, outlier.size = 0.1) +
        scale_fill_manual(values = utils_cb_palette(nlevels(df$group))) +
        stat_summary(
            fun = mean,
            geom = "point",
            size = 12,
            colour = "blue",
            shape = 95
        ) +
        ggtitle(paste(meta$source, meta$signature, sep = " - ")) +
        NoLegend() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    ggsave(
        file.path(out_dir, paste0("Boxplot_", sc, ".png")),
        plot = p,
        width = 20,
        height = 20,
        units = "cm",
        dpi = 200,
        bg = "white"
    )
}

make_violin <- function(sc, out_dir) {
    meta <- score_meta(sc)

    df <- data.frame(
        group = epi$plot_group,
        score = epi@meta.data[[sc]]
    )

    ymax <- max(df$score, na.rm = TRUE)

    p <- ggplot(df, aes(group, score, fill = group)) +
        geom_violin(scale = "width") +
        scale_fill_manual(values = utils_cb_palette(nlevels(df$group))) +
        stat_compare_means(
            method = "anova",
            label.x = 1,
            label.y = ymax + 0.05
        ) +
        ggtitle(paste(meta$source, meta$signature, sep = " - ")) +
        NoLegend() +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1)
        )

    ggsave(
        file.path(out_dir, paste0("Violin_", sc, ".png")),
        plot = p,
        width = 20,
        height = 15,
        units = "cm",
        dpi = 200,
        bg = "white"
    )
}

# ============================================================
# Separate plots: Normal / PCa / Hallmark
# ============================================================

for (src in c("Normal", "PCa", "Hallmark")) {
    dir.create(file.path(OUT_DIR, src), showWarnings = FALSE, recursive = TRUE)
}

for (sc in normal_score_cols) {
    out <- file.path(OUT_DIR, "Normal")
    make_umap(sc, out)
    make_split_umap(sc, out)
    make_boxplot(sc, out)
    make_violin(sc, out)
}

for (sc in pca_score_cols) {
    out <- file.path(OUT_DIR, "PCa")
    make_umap(sc, out)
    make_split_umap(sc, out)
    make_boxplot(sc, out)
    make_violin(sc, out)
}

for (sc in ar_score_cols) {
    out <- file.path(OUT_DIR, "Hallmark")
    make_umap(sc, out)
    make_split_umap(sc, out)
    make_boxplot(sc, out)
    make_violin(sc, out)
}

# ============================================================
# Combined comparison plot: Normal vs PCa signatures
# ============================================================

long_df <- epi@meta.data %>%
    mutate(
        cell = colnames(epi),
        group = epi$plot_group
    ) %>%
    select(cell, group, all_of(c(normal_score_cols, pca_score_cols))) %>%
    pivot_longer(
        cols = all_of(c(normal_score_cols, pca_score_cols)),
        names_to = "score_name",
        values_to = "score"
    ) %>%
    mutate(
        source = sub("^([^_]+)_.*_Score$", "\\1", score_name),
        signature = sub("^[^_]+_(.*)_Score$", "\\1", score_name)
    )

p_all <- ggplot(long_df, aes(group, score, fill = group)) +
    geom_violin(scale = "width") +
    facet_grid(source ~ signature, scales = "free_y") +
    scale_fill_manual(values = utils_cb_palette(nlevels(epi$plot_group))) +
    NoLegend() +
    theme(
        axis.text.x = element_text(angle = 45, hjust = 1),
        strip.text = element_text(size = 10)
    ) +
    ggtitle("Normal vs PCa signature scores")

ggsave(
    file.path(OUT_DIR, "Violin_Normal_vs_PCa_all_signatures.png"),
    plot = p_all,
    width = 34,
    height = 18,
    units = "cm",
    dpi = 200,
    bg = "white"
)

# ============================================================
# Delta plots for shared signatures
# PCa_x - Normal_x
# ============================================================

normal_sigs <- sub("^Normal_(.*)_Score$", "\\1", normal_score_cols)
pca_sigs    <- sub("^PCa_(.*)_Score$", "\\1", pca_score_cols)

shared_sigs <- intersect(normal_sigs, pca_sigs)

dir.create(file.path(OUT_DIR, "Delta_PCa_minus_Normal"),
           showWarnings = FALSE, recursive = TRUE)

for (sig in shared_sigs) {
    normal_col <- paste0("Normal_", sig, "_Score")
    pca_col    <- paste0("PCa_", sig, "_Score")
    delta_col  <- paste0("Delta_PCa_minus_Normal_", sig, "_Score")

    epi@meta.data[[delta_col]] <- epi@meta.data[[pca_col]] - epi@meta.data[[normal_col]]

    x <- epi@meta.data[[delta_col]]
    lim <- max(abs(x), na.rm = TRUE)

    p <- FeaturePlot(
        epi,
        features = delta_col,
        pt.size = 0.3,
        order = TRUE
    ) +
        scale_color_gradient2(
            low = "blue",
            mid = "white",
            high = "red",
            midpoint = 0,
            limits = c(-lim, lim)
        ) +
        ggtitle(paste0("PCa - Normal: ", sig))

    ggsave(
        file.path(
            OUT_DIR,
            "Delta_PCa_minus_Normal",
            paste0("UMAP_Delta_PCa_minus_Normal_", sig, ".png")
        ),
        plot = p,
        width = 20,
        height = 20,
        units = "cm",
        dpi = 200,
        bg = "white"
    )

    df <- data.frame(
        group = epi$plot_group,
        delta = epi@meta.data[[delta_col]]
    )

    p2 <- ggplot(df, aes(group, delta, fill = group)) +
        geom_violin(scale = "width") +
        geom_hline(yintercept = 0, linetype = "dashed") +
        scale_fill_manual(values = utils_cb_palette(nlevels(df$group))) +
        NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
        ggtitle(paste0("PCa - Normal delta: ", sig))

    ggsave(
        file.path(
            OUT_DIR,
            "Delta_PCa_minus_Normal",
            paste0("Violin_Delta_PCa_minus_Normal_", sig, ".png")
        ),
        plot = p2,
        width = 20,
        height = 15,
        units = "cm",
        dpi = 200,
        bg = "white"
    )
}

# Save object with delta scores included
saveRDS(epi, OUT_RDS)

message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
