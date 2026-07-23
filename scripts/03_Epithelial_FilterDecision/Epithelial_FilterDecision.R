# Stage 2.5 (decontX variant) — Filter Decision for Epithelial Reclustering
#
# 원본 Epithelial_FilterDecision.R 와 동일한 진단이되 입력/출력만 decontX 경로로.
#   Input : Results/02_Epithelial_Initial/epi_clustered.rds
#   Output: Results/03_Epithelial_FilterDecision/
# filter_decision.csv 는 모든 cluster keep=TRUE 템플릿으로 나오며, figure/summary
# 검토 후 keep/qc_flag/lineage_call/rationale 를 채운 뒤 stage 04(decontX) 실행.

library(Seurat)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

source("scripts/00_utils/scRNA_utils.R")

IN_RDS  <- "Results/02_Epithelial_Initial/epi_clustered.rds"
OUT_DIR <- "Results/03_Epithelial_FilterDecision"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

epi <- readRDS(IN_RDS)
Idents(epi) <- "seurat_clusters"

cluster_levels <- levels(factor(epi$seurat_clusters))
n_clusters     <- length(cluster_levels)
cluster_cols   <- utils_cb_palette(n_clusters)

# 1. Cluster-level QC metrics ----
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
qc_vln_list <- lapply(qc_features, function(feat) {
    VlnPlot(epi, features = feat, pt.size = 0, cols = cluster_cols) +
        NoLegend() + labs(x = "Cluster", title = feat)
})
ggsave(file.path(OUT_DIR, "cluster_QC_violin.png"),
    plot = wrap_plots(qc_vln_list, ncol = 1),
    width = 12, height = 12, bg = "white"
)

# 2. Lineage purity — epi vs stromal/immune/endothelial ----
lineage_markers <- list(
    Epithelial    = c("EPCAM", "KRT8", "KRT18", "KRT19"),
    Basal_Hillock = c("KRT5", "KRT14", "TP63", "KRT13", "KRT4"),
    Fibroblast    = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
    Smooth_muscle = c("ACTA2", "MYH11", "RGS5", "TAGLN"),
    Endothelial   = c("PECAM1", "VWF", "CDH5"),
    Immune        = c("PTPRC", "CD3D", "CD68", "TPSAB1")
)
p_dot <- DotPlot(epi, features = lineage_markers, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "Lineage purity — cluster별 epi/stromal/immune/endo marker (decontX)")
ggsave(file.path(OUT_DIR, "cluster_lineage_DotPlot.png"),
    plot = p_dot, width = 18, height = 8, bg = "white"
)

key_markers <- c(
    "EPCAM", "KRT8", "KRT5", "TP63",
    "COL1A1", "DCN", "ACTA2", "RGS5",
    "PECAM1", "PTPRC"
)
p_feat <- FeaturePlot(epi, features = key_markers, ncol = 4, pt.size = 0.1)
ggsave(file.path(OUT_DIR, "cluster_lineage_FeaturePlot.png"),
    plot = p_feat, width = 20, height = 14, bg = "white"
)

# 3. SingleR (carry-over) label distribution per cluster ----
if ("SingleR_pruned" %in% colnames(epi@meta.data)) {
    sr_col <- "SingleR_pruned"
} else if ("SingleR" %in% colnames(epi@meta.data)) {
    sr_col <- "SingleR"
} else {
    sr_col <- NULL
}

dominant_singleR <- NULL
if (!is.null(sr_col)) {
    sr_df <- epi@meta.data %>%
        dplyr::mutate(.lab = ifelse(is.na(.data[[sr_col]]),
                                    "NA_pruned", .data[[sr_col]])) %>%
        dplyr::count(seurat_clusters, .lab, name = "n")

    n_labs <- dplyr::n_distinct(sr_df$.lab)
    p_sr <- ggplot(sr_df,
        aes(x = seurat_clusters, y = n, fill = .lab)
    ) +
        geom_col(position = "fill") +
        scale_fill_manual(values = utils_cb_palette(n_labs)) +
        labs(x = "Cluster", y = "Proportion", fill = sr_col,
             title = paste0(sr_col, " distribution per cluster")) +
        theme_classic() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, "cluster_SingleR_distribution.png"),
        plot = p_sr, width = 12, height = 6, bg = "white"
    )

    sr_table <- sr_df %>%
        tidyr::pivot_wider(names_from = .lab, values_from = n, values_fill = 0)
    write.csv(sr_table,
        file.path(OUT_DIR, "cluster_SingleR_table.csv"), row.names = FALSE
    )

    dominant_singleR <- sr_df %>%
        dplyr::group_by(seurat_clusters) %>%
        dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
        dplyr::ungroup() %>%
        dplyr::transmute(seurat_clusters,
                         dominant_SingleR   = .lab,
                         dominant_SingleR_n = n)
}

# 4. Cluster × patient cross-table ----
patient_table <- epi@meta.data %>%
    dplyr::count(seurat_clusters, orig.ident, name = "n") %>%
    tidyr::pivot_wider(names_from = orig.ident, values_from = n, values_fill = 0)
patient_cols <- setdiff(names(patient_table), "seurat_clusters")
patient_table$total <- rowSums(
    patient_table[, patient_cols, drop = FALSE]
)
write.csv(patient_table,
    file.path(OUT_DIR, "cluster_patient_table.csv"), row.names = FALSE
)

dominant_patient <- epi@meta.data %>%
    dplyr::count(seurat_clusters, orig.ident, name = "n") %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::slice_max(n, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(seurat_clusters,
                     dominant_patient   = orig.ident,
                     dominant_patient_n = n)

# 5. Per-cluster summary + decontX contamination ----
qc_summary <- epi@meta.data %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
        n_cells         = dplyr::n(),
        median_nFeature = stats::median(nFeature_RNA),
        median_nCount   = stats::median(nCount_RNA),
        median_pct_mt   = round(stats::median(percent.mt), 2),
        mean_contam     = round(mean(decontX_contamination), 3),
        .groups = "drop"
    )

cluster_summary <- qc_summary %>%
    dplyr::left_join(dominant_patient, by = "seurat_clusters")

if (!is.null(dominant_singleR)) {
    cluster_summary <- cluster_summary %>%
        dplyr::left_join(dominant_singleR, by = "seurat_clusters")
} else {
    cluster_summary$dominant_SingleR   <- NA_character_
    cluster_summary$dominant_SingleR_n <- NA_integer_
}

cluster_summary <- cluster_summary %>%
    dplyr::arrange(suppressWarnings(as.integer(as.character(seurat_clusters))))
write.csv(cluster_summary,
    file.path(OUT_DIR, "cluster_summary.csv"), row.names = FALSE
)

# 6. filter_decision.csv (수동 편집 템플릿) ----
filter_decision <- cluster_summary %>%
    dplyr::transmute(
        cluster          = as.character(seurat_clusters),
        n_cells,
        median_nFeature,
        median_pct_mt,
        mean_contam,
        dominant_patient,
        dominant_SingleR,
        qc_flag      = NA_character_,
        lineage_call = NA_character_,
        keep         = TRUE,
        rationale    = NA_character_
    )
write.csv(filter_decision,
    file.path(OUT_DIR, "filter_decision.csv"), row.names = FALSE
)

message("Stage 2.5 (decontX) figures + tables → ", OUT_DIR)
