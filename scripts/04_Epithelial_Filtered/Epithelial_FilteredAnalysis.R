# Stage 3 — Filtered Reclustering
#
# Stage 2.5 에서 stromal contamination / low-quality 로 판정된 cluster 의 cell 을
# 제외하고 epithelial subset 을 다시 SCT + Harmony + cluster + UMAP 한다.
# 본격 annotation (Hillock, ARPC/DNPC subtype, module score 등) 은
# scripts/05_Epithelial_Downstream/ 에서 이 단계의 출력 RDS 를 입력으로 받아 진행.
#
# Input:
#   - Results/02_Epithelial_Initial/epi_clustered.rds        (Stage 2 결과)
#   - Results/03_Epithelial_FilterDecision/filter_decision.csv   (Stage 2.5 의 keep/drop)
# Output:
#   - Results/04_Epithelial_Filtered/
#       epi_filtered_clustered.rds          최종 2차 reclustering 객체
#       ElbowPlot.png / clustree.png / UMAP_multires.png
#       UMAP_cluster.png / UMAP_patient.png / UMAP_split_by_patient.png
#       CellCycle_UMAP.png
#       cluster_QC_violin.png / cluster_lineage_DotPlot.png   ← post-filter 검증
#       all_markers.csv

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(clustree)

library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scripts/00_utils/scRNA_utils.R")

IN_RDS       <- "Results/02_Epithelial_Initial/epi_clustered.rds"
DECISION_CSV <- "Results/03_Epithelial_FilterDecision/filter_decision.csv"
OUT_DIR      <- "Results/04_Epithelial_Filtered"
OUT_RDS      <- file.path(OUT_DIR, "epi_filtered_clustered.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Resolutions swept for clustree (Stage 2 와 동일)
res_vec <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)

# Post-clustering cluster merge map. Initially empty — clustree.png 검토 후,
# 작은 fragment cluster 를 부모로 합쳐야 하면 여기 채울 것.
# 예: c("12" = "10", "13" = "9")
merge_map <- character()

computed <- !file.exists(OUT_RDS)

if (computed) {

    # =========================================================================
    # 1. Load + filter ----
    # =========================================================================
    epi <- readRDS(IN_RDS)

    if (!file.exists(DECISION_CSV)) {
        stop(
            "filter_decision.csv 가 없음: ", DECISION_CSV,
            ". 먼저 scripts/03_Epithelial_FilterDecision/Epithelial_FilterDecision.R ",
            "실행 후 keep 컬럼을 채워주세요.",
            call. = FALSE
        )
    }
    decision <- read.csv(DECISION_CSV, stringsAsFactors = FALSE)
    if (!all(c("cluster", "keep") %in% names(decision))) {
        stop("filter_decision.csv 에 cluster, keep 컬럼 필요", call. = FALSE)
    }

    keep_clusters <- as.character(
        decision$cluster[decision$keep %in% c(TRUE, "TRUE", "true")]
    )
    drop_clusters <- setdiff(as.character(decision$cluster), keep_clusters)
    message("Keep ", length(keep_clusters), " clusters; drop ",
            length(drop_clusters))
    if (length(drop_clusters) > 0) {
        message("  dropped: ", paste(drop_clusters, collapse = ", "))
    }

    epi$keep_after_qc <- as.character(epi$seurat_clusters) %in% keep_clusters
    n_before <- ncol(epi)
    epi <- subset(epi, subset = keep_after_qc)
    message("Cells: ", n_before, " → ", ncol(epi),
            " (", n_before - ncol(epi), " dropped)")

    # =========================================================================
    # 2. Reset assay/reduction state ----
    # =========================================================================
    DefaultAssay(epi) <- "RNA"
    if ("SCT"     %in% names(epi@assays))     epi[["SCT"]]     <- NULL
    if ("harmony" %in% names(epi@reductions)) epi[["harmony"]] <- NULL
    if ("pca"     %in% names(epi@reductions)) epi[["pca"]]     <- NULL
    if ("umap"    %in% names(epi@reductions)) epi[["umap"]]    <- NULL

    epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
    epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)

    # Stage 2 의 stale clustering metadata 제거 — clustree 가 옛 SCT_snn_res.* 를
    # 같이 잡아 그려지지 않도록.
    epi@meta.data <- dplyr::select(epi@meta.data,
        -matches("_snn_res\\."),
        -any_of(c("seurat_clusters", "keep_after_qc"))
    )

    # =========================================================================
    # 3. SCT + Harmony + Cluster + UMAP (Stage 2 와 동일 파라미터) ----
    # =========================================================================
    epi <- SCTransform(epi, verbose = FALSE)
    epi <- RunPCA(epi, verbose = FALSE)

    epi <- IntegrateLayers(
        object = epi,
        method = HarmonyIntegration,
        orig.reduction = "pca",
        new.reduction = "harmony",
        normalization.method = "SCT",
        verbose = FALSE
    )
    epi[["RNA"]] <- JoinLayers(epi[["RNA"]])

    ggsave(file.path(OUT_DIR, "ElbowPlot.png"),
        plot = ElbowPlot(epi, ndims = 50),
        width = 12, height = 8, bg = "white"
    )

    epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)
    epi <- FindClusters(epi, resolution = res_vec)

    # Default ident = res 0.4 + 사용자 merge_map (위에서 채울 것)
    base_ident <- as.character(epi$SCT_snn_res.0.4)
    curated <- if (length(merge_map) > 0) {
        dplyr::recode(base_ident, !!!merge_map)
    } else {
        base_ident
    }
    epi$seurat_clusters <- factor(curated,
        levels = as.character(sort(unique(
            suppressWarnings(as.integer(curated))
        )))
    )
    Idents(epi) <- "seurat_clusters"

    epi <- RunUMAP(epi,
        reduction = "harmony", dims = 1:30,
        n.neighbors = 20, min.dist = 0.1, spread = 4.0
    )

    # =========================================================================
    # 4. Cell cycle re-scoring (SCT 모델이 바뀌었으므로 재계산) ----
    # =========================================================================
    s.genes   <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    epi <- CellCycleScoring(epi,
        s.features = s.genes, g2m.features = g2m.genes, nbin = 12)

} else {
    message("Loading cached ", OUT_RDS, " — regenerating figures only")
    epi <- readRDS(OUT_RDS)
    Idents(epi) <- "seurat_clusters"
}

# =============================================================================
# 5. Diagnostic figures ----
# =============================================================================
p_tree <- clustree(epi, prefix = "SCT_snn_res.")
ggsave(file.path(OUT_DIR, "clustree.png"),
    plot = p_tree, width = 12, height = 14, dpi = 200, bg = "white"
)

umap_list <- lapply(res_vec, function(r) {
    rcol <- paste0("SCT_snn_res.", r)
    DimPlot(epi, group.by = rcol, label = TRUE, pt.size = 0.2,
            cols = utils_cb_palette(dplyr::n_distinct(epi@meta.data[[rcol]]))) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
ggsave(file.path(OUT_DIR, "UMAP_multires.png"),
    plot = wrap_plots(umap_list, ncol = 3),
    width = 18, height = 14, dpi = 200, bg = "white"
)

n_final <- nlevels(factor(epi$seurat_clusters))
p1 <- DimPlot(epi, label = TRUE, pt.size = 0.3,
              cols = utils_cb_palette(n_final)) +
    ggtitle("Filtered epithelial clusters")
p2 <- DimPlot(epi, group.by = "orig.ident", pt.size = 0.3,
              cols = utils_cb_palette(dplyr::n_distinct(epi$orig.ident))) +
    ggtitle("By patient")
ggsave(file.path(OUT_DIR, "UMAP_cluster.png"),
    plot = p1, width = 10, height = 8, bg = "white")
ggsave(file.path(OUT_DIR, "UMAP_patient.png"),
    plot = p2, width = 10, height = 8, bg = "white")

p_split <- DimPlot(epi, split.by = "orig.ident", label = TRUE, pt.size = 0.3,
              cols = utils_cb_palette(n_final))
ggsave(file.path(OUT_DIR, "UMAP_split_by_patient.png"),
    plot = p_split, width = 24, height = 8, bg = "white")

p_phase <- DimPlot(epi, group.by = "Phase", pt.size = 0.3,
              cols = utils_cb_palette(dplyr::n_distinct(epi$Phase)))
ggsave(file.path(OUT_DIR, "CellCycle_UMAP.png"),
    plot = p_phase, width = 10, height = 7, bg = "white")

# =============================================================================
# 6. Post-filter QC verification ----
# =============================================================================
# Stage 2.5 와 동일 형식의 QC + lineage purity figure 를 filter 후 객체에
# 다시 그려, stromal residue 가 의도대로 제거되었는지 시각 확인.
qc_features <- c("nFeature_RNA", "nCount_RNA", "percent.mt")
qc_vln_list <- lapply(qc_features, function(feat) {
    VlnPlot(epi, features = feat, pt.size = 0,
            cols = utils_cb_palette(n_final)) +
        NoLegend() + labs(x = "Cluster", title = feat)
})
ggsave(file.path(OUT_DIR, "cluster_QC_violin.png"),
    plot = wrap_plots(qc_vln_list, ncol = 1),
    width = 12, height = 12, bg = "white"
)

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
    labs(title = "Post-filter lineage purity check")
ggsave(file.path(OUT_DIR, "cluster_lineage_DotPlot.png"),
    plot = p_dot, width = 18, height = 8, bg = "white"
)

# =============================================================================
# 7. Find markers + save ----
# =============================================================================
if (computed) {
    utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(epi, OUT_RDS)
}

message("Stage 3 결과 → ", OUT_DIR)
message("다음: scripts/05_Epithelial_Downstream/*.R 의 입력 RDS 경로를 ",
        OUT_RDS, " 로 교체 후 재실행")
