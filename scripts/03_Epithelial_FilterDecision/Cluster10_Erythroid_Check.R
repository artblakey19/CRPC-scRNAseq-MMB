# Cluster 10 erythroid marker check — drop rationale
#
# FilterDecision 단계에서 cluster 10을 "low_quality / ambient" 로 배제했는데
# (filter_decision.csv L12), HBB/HBA1/HBA2 ambient hemoglobin contamination
# 가설이 실제로 맞는지(=진성 erythroblast/RBC cluster 가 아닌지)를
# 적혈구 lineage 마커 패널로 재확인한다.
#
#   Erythroid lineage panel
#     - HBB / HBA1 / HBA2 : hemoglobin chains (ambient에서도 흔히 검출)
#     - ALAS2             : heme synthesis (erythroid-specific)
#     - GYPA              : RBC membrane glycophorin A
#     - KLF1              : erythroid TF
#
#   진성 erythroid cluster 라면 → 위 마커 모두에서 % expressing 이 매우 높고
#   세포 대부분이 동시에 양성이어야 함. Ambient 오염이면 → HBB/HBA1/HBA2 의
#   raw count 만 일부 droplet에서 폭증하고 (% expressing 은 낮음), ALAS2/GYPA/KLF1
#   은 거의 부재.
#
# Input : Results/02_Epithelial_Initial/epi_clustered.rds
# Output: Results/03_Epithelial_FilterDecision/
#   - cluster10_erythroid_FeaturePlot.png
#   - cluster10_erythroid_VlnPlot.png
#   - cluster10_erythroid_DotPlot.png
#   - cluster10_erythroid_summary.csv  (per-cluster % expressing & mean counts)

suppressPackageStartupMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
})

source("scripts/00_utils/scRNA_utils.R")

IN_RDS  <- "Results/02_Epithelial_Initial/epi_clustered.rds"
OUT_DIR <- "Results/03_Epithelial_FilterDecision"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

epi <- readRDS(IN_RDS)
Idents(epi) <- "seurat_clusters"

ery_genes <- c("HBB", "HBA1", "HBA2", "ALAS2", "GYPA", "KLF1")

# =============================================================================
# 1. Per-cluster summary table — % expressing + mean raw counts ----
# =============================================================================
# RNA assay 의 raw counts 기준. Ambient 오염은 mean count 가 폭증해도
# % expressing 은 낮은 특징적인 패턴을 보임.
DefaultAssay(epi) <- "RNA"
present <- intersect(ery_genes, rownames(epi[["RNA"]]))
missing <- setdiff(ery_genes, rownames(epi[["RNA"]]))

cnts <- LayerData(epi, assay = "RNA", layer = "counts")[present, , drop = FALSE]
raw  <- as.data.frame(t(as.matrix(cnts)))
raw$cluster <- as.character(epi$seurat_clusters)

pct_expr <- raw %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
        n = dplyr::n(),
        dplyr::across(dplyr::all_of(present),
                      ~ round(mean(.x > 0) * 100, 2),
                      .names = "pct_{.col}"),
        .groups = "drop"
    )
mean_cnt <- raw %>%
    dplyr::group_by(cluster) %>%
    dplyr::summarise(
        dplyr::across(dplyr::all_of(present),
                      ~ round(mean(.x), 2),
                      .names = "mean_{.col}"),
        .groups = "drop"
    )

summary_tbl <- dplyr::left_join(pct_expr, mean_cnt, by = "cluster") %>%
    dplyr::arrange(suppressWarnings(as.integer(cluster)))

# Mark genes absent from feature matrix so the CSV records that explicitly.
for (g in missing) {
    summary_tbl[[paste0("pct_", g)]]  <- NA_real_
    summary_tbl[[paste0("mean_", g)]] <- NA_real_
}
summary_tbl <- summary_tbl[, c(
    "cluster", "n",
    paste0("pct_",  ery_genes),
    paste0("mean_", ery_genes)
)]
write.csv(summary_tbl,
    file.path(OUT_DIR, "cluster10_erythroid_summary.csv"), row.names = FALSE
)

# =============================================================================
# 2. FeaturePlot — UMAP에서 hemoglobin이 cluster 10에 국한되는지 ----
# =============================================================================
# SCT data layer 의 normalized expression 으로 시각화. ALAS2/GYPA/KLF1 은
# 거의 0 일 것이라 plot이 비어도 그것 자체가 ambient 가설 증거.
DefaultAssay(epi) <- "SCT"
sct_present <- intersect(ery_genes, rownames(epi[["SCT"]]))

p_feat <- FeaturePlot(
    epi, features = sct_present, ncol = 3, pt.size = 0.1,
    order = TRUE
) & scale_color_viridis_c(option = "magma", direction = -1)
ggsave(file.path(OUT_DIR, "cluster10_erythroid_FeaturePlot.png"),
    plot = p_feat, width = 18, height = 12, bg = "white"
)

# =============================================================================
# 3. VlnPlot — cluster 별 분포 (cluster 10 만 튀는지) ----
# =============================================================================
n_clusters   <- nlevels(factor(epi$seurat_clusters))
cluster_cols <- utils_cb_palette(n_clusters)

p_vln <- VlnPlot(
    epi, features = sct_present, pt.size = 0,
    cols = cluster_cols, ncol = 3
) & NoLegend() & labs(x = "Cluster")
ggsave(file.path(OUT_DIR, "cluster10_erythroid_VlnPlot.png"),
    plot = p_vln, width = 18, height = 12, bg = "white"
)

# =============================================================================
# 4. DotPlot — % expressing × avg expression 동시 비교 ----
# =============================================================================
p_dot <- DotPlot(epi, features = sct_present, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "Erythroid lineage markers per cluster",
         subtitle = "Drop rationale for cluster 10 (ambient hemoglobin vs true erythroid)")
ggsave(file.path(OUT_DIR, "cluster10_erythroid_DotPlot.png"),
    plot = p_dot, width = 12, height = 7, bg = "white"
)

# =============================================================================
# 5. Console summary — cluster 10 vs the rest ----
# =============================================================================
c10 <- summary_tbl %>% dplyr::filter(cluster == "10")
others <- summary_tbl %>% dplyr::filter(cluster != "10")

message("\n--- Cluster 10 (n = ", c10$n, ") ---")
print(as.data.frame(c10), row.names = FALSE)

message("\n--- Other clusters: max % expressing per gene (for reference) ---")
ref <- sapply(present, function(g) max(others[[paste0("pct_", g)]], na.rm = TRUE))
print(round(ref, 2))

message("\nFigures + summary CSV → ", OUT_DIR)
if (length(missing) > 0) {
    message("NOTE: feature matrix 에 없는 gene (전 세포 0) → ",
            paste(missing, collapse = ", "))
}
