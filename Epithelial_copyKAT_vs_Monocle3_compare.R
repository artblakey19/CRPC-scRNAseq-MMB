# ============================================================
# Epithelial_copyKAT_vs_Monocle3_compare.R
#
# Compare the malignisation process inferred two ways on the SAME epithelial
# cells, per patient:
#   - monocle3 : transcriptional pseudotime (graph distance from the cluster-8
#                root) — Epithelial_Monocle3_persample.R
#   - copyKAT  : genomic CNA progression (CNA-space distance from the cluster-8
#                diploid centroid) + CNV-defined malignant "forms"
#
# Both methods share the SAME root — cluster 8, the diploid anchor / copyKAT
# norm.cell.names population — so each yields a "distance from the benign root".
# That shared origin is what makes the two progressions directly comparable.
#
# Hypothesis under test (DNPC late-stage cohort):
#   H1  aneuploid cells split into 2-3 "forms" by CNV pattern
#   H2  those CNV forms map onto distinct transcriptional terminal states
#   H3  CNA accumulation is monotonic along transcriptional pseudotime
#
# Per-patient outputs:
#   {sid}_silhouette_aneuploid_K.png  - avg silhouette vs K for form calling (H1)
#   {sid}_DualProgression_scatter.png - pseudotime vs CNA distance, loess + rho (H3)
#   {sid}_CNV_pseudotime_heatmap.png  - per-cell CNA, rows ordered by pseudotime,
#                                       split by CNV form    <- main plot (H1+H3)
#   {sid}_Form_vs_cluster_alluvial.png- CNV form <-> transcriptional cluster (H2)
# Combined:
#   DualProgression_scatter_allpatients.png
#   copykat_monocle3_percell.csv
#   Form_spearman_summary.csv
# ============================================================

library(Seurat)
library(data.table)
library(cluster)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggalluvial)
library(dplyr)

# ============================================================
# Config ----
# ============================================================
REF_CLUSTER <- "8"          # diploid anchor cluster — shared root of both methods
K_RANGE     <- 2:5          # candidate CNV-form counts scanned per patient
N_CORES     <- 30           # threads for parallelDist
CK_ROOT     <- "Results/Epithelial/copyKAT_perpatient"
M3_CSV      <- "Results/Epithelial/Monocle3_persample/pseudotime_persample.csv"
OUT_DIR     <- "Results/Epithelial/copyKAT_vs_Monocle3"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

cna_cols <- colorRamp2(c(-0.10, 0, 0.10), c("#3B4CC0", "white", "#B40426"))
ck_cols  <- c(aneuploid = "#D7261E", diploid = "#1F77B4", not.defined = "grey80")

# ============================================================
# Load shared inputs ----
# ============================================================
epi <- readRDS(file.path(CK_ROOT, "epi_copykat_perpatient.rds"))
epi$seurat_clusters <- as.character(epi$seurat_clusters)

meta <- data.frame(
    cell            = colnames(epi),
    patient         = as.character(epi$orig.ident),
    seurat_cluster  = epi$seurat_clusters,
    copykat_clu8ref = epi$copykat_clu8ref,
    copykat_subclone = epi$copykat_subclone,
    stringsAsFactors = FALSE)

# monocle3 per-cell pseudotime; pst_norm is per-patient min-max rescaled.
m3 <- data.table::fread(M3_CSV)[, .(cell, pseudotime, pst_norm)]
meta <- dplyr::left_join(meta, as.data.frame(m3), by = "cell")

patients <- sort(unique(meta$patient))
message("Patients: ", paste(patients, collapse = ", "))

# ============================================================
# Helper: per-patient comparison ----
# ============================================================
run_patient <- function(sid) {
    message("\n==== compare: ", sid, " ====")
    sample_dir <- file.path(CK_ROOT, sid)
    out_dir    <- file.path(OUT_DIR, sid)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

    # ---- copyKAT CNA matrix (bins x cells) ----
    cna <- data.table::fread(
        file.path(sample_dir, paste0(sid, "_copykat_CNA_results.txt")),
        showProgress = FALSE)
    bin_chr <- cna$chrom
    mat <- as.matrix(cna[, -(1:3)])
    colnames(mat) <- sub("\\.1$", "-1", colnames(mat))
    rm(cna); gc()

    pm <- meta[meta$patient == sid, ]
    pm <- pm[pm$cell %in% colnames(mat), ]
    mat <- mat[, pm$cell]

    # ---- genomic "pseudotime": CNA distance from the cluster-8 diploid centroid
    ref_cells <- pm$cell[pm$seurat_cluster == REF_CLUSTER]
    if (length(ref_cells) < 20)
        warning("  cluster-", REF_CLUSTER, " reference < 20 cells for ", sid,
                call. = FALSE)
    ref_centroid <- rowMeans(mat[, ref_cells, drop = FALSE])
    pm$cna_dist   <- sqrt(colSums((mat - ref_centroid)^2))
    pm$cna_burden <- colMeans(abs(mat))   # secondary, scale-free CNA summary

    # ---- define CNV "forms" on ANEUPLOID cells only (the malignant forms) ----
    # bbaf076 method: Euclidean + Ward.D2 on CNA profiles, K by max avg silhouette.
    aneu <- pm$cell[pm$copykat_clu8ref == "aneuploid"]
    pm$cnv_form <- NA_character_
    best_k <- NA_integer_
    if (length(aneu) >= 2 * max(K_RANGE)) {
        d_aneu  <- parallelDist::parDist(t(mat[, aneu]), method = "euclidean",
                                         threads = N_CORES)
        hc_aneu <- hclust(d_aneu, method = "ward.D2")
        sil_avg <- vapply(K_RANGE, function(k)
            mean(cluster::silhouette(cutree(hc_aneu, k = k), d_aneu)[, "sil_width"]),
            numeric(1))
        names(sil_avg) <- K_RANGE
        best_k <- as.integer(names(sil_avg)[which.max(sil_avg)])
        message("  aneuploid cells = ", length(aneu),
                " | avg silhouette: ",
                paste(sprintf("K%d=%.3f", K_RANGE, sil_avg), collapse = "  "),
                " -> K=", best_k)

        png(file.path(out_dir, paste0(sid, "_silhouette_aneuploid_K.png")),
            width = 5, height = 4, units = "in", res = 200, bg = "white")
        par(mar = c(4, 4, 3, 1))
        plot(K_RANGE, sil_avg, type = "b", pch = 19,
             xlab = "K (CNV forms)", ylab = "mean silhouette width",
             main = paste0(sid, " — aneuploid CNV forms, optimal K = ", best_k))
        abline(v = best_k, col = "#D7261E", lty = 2)
        dev.off()

        form_id <- cutree(hc_aneu, k = best_k)
        pm$cnv_form[match(aneu, pm$cell)] <- paste0(sid, "_F", form_id)
    } else {
        warning("  too few aneuploid cells for form calling in ", sid,
                call. = FALSE)
    }

    # ---- row grouping for the heatmap: forms + non-aneuploid baseline ----
    pm$group <- ifelse(is.na(pm$cnv_form), pm$copykat_clu8ref, pm$cnv_form)
    grp_levels <- c(intersect(c("diploid", "not.defined"), pm$group),
                    sort(unique(na.omit(pm$cnv_form))))
    pm$group <- factor(pm$group, levels = grp_levels)

    # ---- Plot: pseudotime-ordered per-cell CNV heatmap (main, H1 + H3) ----
    hm <- pm[is.finite(pm$pseudotime), ]
    hm <- hm[order(hm$group, hm$pseudotime), ]
    M  <- t(mat[, hm$cell])
    chr_factor <- droplevels(factor(paste0("chr", bin_chr),
                                    levels = paste0("chr", c(1:22, "X", "Y"))))
    pst_col <- colorRamp2(range(hm$pseudotime), c("grey92", "#08306B"))
    ht <- Heatmap(
        M, name = "CNA", col = cna_cols,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_split = hm$group, row_gap = unit(1, "mm"),
        show_row_names = FALSE, show_column_names = FALSE,
        column_split = chr_factor, column_gap = unit(0.5, "mm"),
        column_title_gp = gpar(fontsize = 8),
        row_title_gp = gpar(fontsize = 9), border = TRUE,
        use_raster = TRUE, raster_quality = 4,
        left_annotation = rowAnnotation(
            pseudotime = hm$pseudotime,
            col = list(pseudotime = pst_col),
            annotation_name_gp = gpar(fontsize = 8)),
        heatmap_legend_param = list(direction = "horizontal",
                                    title_position = "lefttop"))
    png(file.path(out_dir, paste0(sid, "_CNV_pseudotime_heatmap.png")),
        width = 18, height = 9, units = "in", res = 300, bg = "white")
    draw(ht, heatmap_legend_side = "bottom", merge_legend = TRUE,
         column_title = paste0(
             sid, " — per-cell copyKAT CNA, rows split by CNV form, ",
             "ordered by monocle3 pseudotime within each form"))
    dev.off()

    # ---- Plot: dual-progression scatter (H3) ----
    sc <- pm[is.finite(pm$pseudotime), ]
    sc_aneu <- sc[!is.na(sc$cnv_form), ]
    rho_all <- if (nrow(sc_aneu) > 5)
        suppressWarnings(cor(sc_aneu$pseudotime, sc_aneu$cna_dist,
                             method = "spearman")) else NA_real_
    p_sc <- ggplot() +
        geom_point(data = sc[is.na(sc$cnv_form), ],
                   aes(pseudotime, cna_dist), colour = "grey80", size = 0.4) +
        geom_point(data = sc_aneu,
                   aes(pseudotime, cna_dist, colour = cnv_form), size = 0.5) +
        geom_smooth(data = sc_aneu,
                    aes(pseudotime, cna_dist, colour = cnv_form),
                    method = "loess", se = FALSE, linewidth = 0.9) +
        labs(x = "monocle3 pseudotime (transcriptional progression)",
             y = "CNA distance from cluster-8 centroid (genomic progression)",
             colour = "CNV form",
             title = paste0(sid, " — transcriptional vs genomic malignisation"),
             subtitle = sprintf("Spearman rho (aneuploid) = %.2f", rho_all)) +
        theme_classic()
    ggsave(file.path(out_dir, paste0(sid, "_DualProgression_scatter.png")),
           plot = p_sc, width = 8, height = 6, bg = "white")

    # ---- Plot: CNV form <-> transcriptional cluster alluvial (H2) ----
    if (any(!is.na(pm$cnv_form))) {
        al <- pm[!is.na(pm$cnv_form), ] |>
            dplyr::count(cnv_form, seurat_cluster)
        p_al <- ggplot(al, aes(axis1 = cnv_form, axis2 = seurat_cluster, y = n)) +
            geom_alluvium(aes(fill = cnv_form), alpha = 0.7) +
            geom_stratum(width = 0.3) +
            geom_text(stat = "stratum",
                      aes(label = after_stat(stratum)), size = 3) +
            scale_x_discrete(limits = c("CNV form", "Transcriptional cluster"),
                             expand = c(0.1, 0.1)) +
            labs(y = "cells", fill = "CNV form",
                 title = paste0(sid, " — CNV form vs transcriptional cluster")) +
            theme_classic()
        ggsave(file.path(out_dir, paste0(sid, "_Form_vs_cluster_alluvial.png")),
               plot = p_al, width = 7, height = 6, bg = "white")
    }

    rm(mat, M); gc()
    pm$best_k <- best_k
    pm
}

# ============================================================
# Run all patients ----
# ============================================================
percell <- dplyr::bind_rows(lapply(patients, run_patient))
write.csv(percell[, c("cell", "patient", "seurat_cluster", "copykat_clu8ref",
                      "copykat_subclone", "cnv_form", "cna_dist", "cna_burden",
                      "pseudotime", "pst_norm")],
          file.path(OUT_DIR, "copykat_monocle3_percell.csv"), row.names = FALSE)

# ---- combined dual-progression scatter, faceted by patient ----
# Strip the patient prefix so F1/F2/.. reuse one palette per facet (forms are
# patient-private and not cross-comparable).
sc_all <- percell[is.finite(percell$pseudotime), ]
sc_all$form_id <- sub("^.*_(F[0-9]+)$", "\\1", sc_all$cnv_form)
p_all <- ggplot() +
    geom_point(data = sc_all[is.na(sc_all$cnv_form), ],
               aes(pst_norm, cna_dist), colour = "grey80", size = 0.3) +
    geom_point(data = sc_all[!is.na(sc_all$cnv_form), ],
               aes(pst_norm, cna_dist, colour = form_id), size = 0.4) +
    geom_smooth(data = sc_all[!is.na(sc_all$cnv_form), ],
                aes(pst_norm, cna_dist, colour = form_id),
                method = "loess", se = FALSE, linewidth = 0.8) +
    facet_wrap(~patient, scales = "free", ncol = length(patients)) +
    labs(x = "monocle3 pseudotime (per-patient min-max normalised)",
         y = "CNA distance from cluster-8 centroid",
         colour = "CNV form",
         title = "Transcriptional vs genomic malignisation — per patient") +
    theme_classic()
ggsave(file.path(OUT_DIR, "DualProgression_scatter_allpatients.png"),
       plot = p_all, width = 7 * length(patients), height = 6, bg = "white")

# ---- per-form monotonicity summary (H3) ----
summary_df <- percell |>
    dplyr::filter(!is.na(cnv_form), is.finite(pseudotime)) |>
    dplyr::group_by(patient, cnv_form) |>
    dplyr::summarise(
        n               = dplyr::n(),
        spearman_rho    = suppressWarnings(
            cor(pseudotime, cna_dist, method = "spearman")),
        median_pst      = median(pseudotime),
        median_cna_dist = median(cna_dist),
        .groups = "drop")
write.csv(summary_df, file.path(OUT_DIR, "Form_spearman_summary.csv"),
          row.names = FALSE)

# ============================================================
# Report ----
# ============================================================
cat("\n===== CNV forms per patient (silhouette-optimal K on aneuploid cells) =====\n")
print(percell |>
      dplyr::filter(!is.na(cnv_form)) |>
      dplyr::count(patient, cnv_form) |>
      as.data.frame())
cat("\n===== per-form monotonicity (Spearman rho: pseudotime vs CNA distance) =====\n")
print(as.data.frame(summary_df))
cat("\nSaved outputs to:", OUT_DIR, "\n")
