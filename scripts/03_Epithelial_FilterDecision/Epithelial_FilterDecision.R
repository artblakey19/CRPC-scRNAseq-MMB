# Stage 2.5 — Filter Decision for Epithelial Reclustering
#
# 1차 epithelial reclustering (Epithelial_Analysis.R) 결과를 진단하여, stromal
# contamination · low-quality cluster · single-patient noise cluster 중 어느
# 것을 배제할지 판단하기 위한 figures/tables 를 생성한다. 실제 cell 배제와
# reclustering은 다음 단계 (scripts/04_Epithelial_Filtered/) 에서 수행한다.
#
# Input : Results/02_Epithelial_Initial/epi_clustered.rds
# Output: Results/03_Epithelial_FilterDecision/
#   - cluster_QC_violin.png            cluster별 nFeature/nCount/percent.mt
#   - cluster_lineage_DotPlot.png      epi / stromal / immune / endo marker
#   - cluster_lineage_FeaturePlot.png  핵심 lineage marker FeaturePlot
#   - cluster_copykat_barplot.png      cluster별 aneuploid/diploid 비율
#   - cluster_SingleR_distribution.png HPCA SingleR_pruned 분포 (carry-over)
#   - cluster_SingleR_table.csv
#   - cluster_patient_table.csv        cluster × patient 교차표
#   - cluster_summary.csv              종합 진단 테이블
#   - filter_decision.csv              keep/drop 결정용 템플릿 (수동 편집)

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

# =============================================================================
# 1. Cluster-level QC metrics ----
# =============================================================================
# Low-quality cluster 식별 — 다른 cluster 대비 nFeature/nCount 가 현저히 낮거나
# percent.mt 가 높은 cluster는 배제 후보.
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
qc_vln_list <- lapply(qc_features, function(feat) {
    VlnPlot(epi, features = feat, pt.size = 0, cols = cluster_cols) +
        NoLegend() + labs(x = "Cluster", title = feat)
})
ggsave(file.path(OUT_DIR, "cluster_QC_violin.png"),
    plot = wrap_plots(qc_vln_list, ncol = 1),
    width = 12, height = 12, bg = "white"
)

# =============================================================================
# 2. Lineage purity — epi vs stromal/immune/endothelial ----
# =============================================================================
# Cluster 에서 epithelial marker (EPCAM/KRT) 가 약하고 stromal/immune/endo
# marker 가 강하면 contamination 후보. Basal/Hillock 축 (KRT5/14/TP63/KRT13/KRT4)
# 을 같이 그려 squamous-like vs myoepi/basal stroma 경계 cluster 도 본다.
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
    labs(title = "Lineage purity — cluster별 epi/stromal/immune/endo marker")
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

# =============================================================================
# 3. copyKAT cluster composition ----
# =============================================================================
# Stromal contamination cluster 는 diploid 우세인 경향. Epi 진성 cluster 는
# ARPC vs club/hillock 에 따라 aneuploid 비율이 다양 — diploid 우세 자체가
# 자동 배제 사유는 아님. 다른 진단과 교차 확인 용도.
copykat_cols <- c(
    "aneuploid" = "#D7261E", "diploid" = "#1F77B4", "not.defined" = "grey80"
)
p_ck <- ggplot(epi@meta.data,
    aes(x = seurat_clusters, fill = copykat_prediction)
) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = copykat_cols, na.value = "grey90") +
    labs(x = "Cluster", y = "Proportion", fill = "copyKAT",
         title = "copyKAT prediction per cluster") +
    theme_classic()
ggsave(file.path(OUT_DIR, "cluster_copykat_barplot.png"),
    plot = p_ck, width = 10, height = 6, bg = "white"
)

# =============================================================================
# 4. SingleR (carry-over) label distribution per cluster ----
# =============================================================================
# Integrated 단계에서 부여한 HPCA-based SingleR_pruned label 분포. Stromal/immune
# label 이 우세한 cluster 는 1차 celltype 정의 단계에서 잘못 epi 로 넘어왔을 가능성.
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

# =============================================================================
# 5. Cluster × patient cross-table ----
# =============================================================================
# Single-patient cluster (한 환자에 거의 몰림) 는 tumor-specific subclone 일 수도
# 있고 noise/doublet 일 수도. cell count 와 함께 봐서 판단.
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

# =============================================================================
# 6. Per-cluster summary (배제 결정 종합 테이블) ----
# =============================================================================
qc_summary <- epi@meta.data %>%
    dplyr::group_by(seurat_clusters) %>%
    dplyr::summarise(
        n_cells         = dplyr::n(),
        median_nFeature = stats::median(nFeature_RNA),
        median_nCount   = stats::median(nCount_RNA),
        median_pct_mt   = round(stats::median(percent.mt), 2),
        .groups = "drop"
    )

copykat_long <- epi@meta.data %>%
    dplyr::count(seurat_clusters, copykat_prediction, name = "n") %>%
    tidyr::pivot_wider(names_from = copykat_prediction, values_from = n,
                       values_fill = 0)
ck_cols <- setdiff(names(copykat_long), "seurat_clusters")
copykat_long$total_ck <- rowSums(copykat_long[, ck_cols, drop = FALSE])
copykat_long$aneuploid_frac <- if ("aneuploid" %in% ck_cols) {
    round(copykat_long$aneuploid / copykat_long$total_ck, 3)
} else {
    NA_real_
}

cluster_summary <- qc_summary %>%
    dplyr::left_join(
        copykat_long %>% dplyr::select(seurat_clusters, aneuploid_frac),
        by = "seurat_clusters"
    ) %>%
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

# =============================================================================
# 7. filter_decision.csv (수동 편집 템플릿) ----
# =============================================================================
# Stage 3 (Epithelial_FilteredAnalysis.R) 가 이 파일을 읽어 keep == TRUE 인
# cluster 만 보존. 우선 모든 cluster 를 keep = TRUE 로 출력 — figure/summary
# 검토 후 수동으로 FALSE 와 rationale 을 채워넣은 뒤 Stage 3 를 실행.
filter_decision <- cluster_summary %>%
    dplyr::transmute(
        cluster          = as.character(seurat_clusters),
        n_cells,
        median_nFeature,
        median_pct_mt,
        aneuploid_frac,
        dominant_patient,
        dominant_SingleR,
        qc_flag      = NA_character_,
        lineage_call = NA_character_,
        copykat_call = NA_character_,
        keep         = TRUE,
        rationale    = NA_character_
    )
write.csv(filter_decision,
    file.path(OUT_DIR, "filter_decision.csv"), row.names = FALSE
)

message("Stage 2.5 figures + tables → ", OUT_DIR)
message(
    "다음 단계: ", file.path(OUT_DIR, "filter_decision.csv"),
    " 를 검토하고 keep/qc_flag/lineage_call/rationale 채운 뒤 ",
    "scripts/04_Epithelial_Filtered/Epithelial_FilteredAnalysis.R 실행"
)
