# Stage 1 (decontX variant) — Integrated analysis with ambient-RNA decontamination
#
# 기존 scripts/01_Integrated/Integrated_Analysis.R 와 동일한 파이프라인이되,
# per-sample 로 decontX (ambient RNA 제거) 를 scDblFinder 앞에 삽입한 버전.
# 원본 결과 객체를 보존하기 위해 출력은 별도 경로(Results/01_Integrated/)로
# 나가고, 원본 스크립트는 그대로 둔다.
#
# 파이프라인 순서:
#   Read10X(filtered) → CreateSeuratObject(min.cells=3, min.features=200)
#     → decontX (sample별, seed 고정) → round(decontXcounts) 로 RNA counts 교체
#     → nCount/nFeature/percent.mt 재계산
#     → scDblFinder (decontaminated counts) → merge → QC filter
#     → SCT → Harmony(seed) → cluster res 0.5 → CellCycle → SingleR
#     → [checkpoint 저장]  ── celltype RenameIdents 는 새 cluster 검토 후 아래 map 으로 부여
#     → stromal refinement: Fibroblast+SMC 만 재클러스터링(SCT→Harmony→res 0.3) 후
#       subcluster 맵으로 celltype 갱신 (Pericyte / Schwann cell 추가) → 최종 저장
#
# 배경 매트릭스(raw/empty droplet)는 이 데이터에 없으므로 decontX default 모드
# (filtered matrix + 내부 cluster 로 ambient 추정)로 실행 — decontX 는 SoupX 와 달리
# empty-droplet 프로파일이 필수가 아니다.

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(scDblFinder)
library(SingleCellExperiment)
library(BiocParallel)
library(SingleR)
library(celldex)
library(decontX)
# Multi-core processing
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)
# Fixed worker count + RNGseed makes scDblFinder reproducible across machines
bpp <- MulticoreParam(workers = 30, RNGseed = 42)

source("scripts/00_utils/scRNA_utils.R")

OUT_DIR <- "Results/01_Integrated"
QC_DIR  <- file.path(OUT_DIR, "QC")
DX_DIR  <- file.path(OUT_DIR, "DecontX")
dir.create(QC_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(DX_DIR, showWarnings = FALSE, recursive = TRUE)

INT_RDS  <- file.path(OUT_DIR, "combined_CRPC.rds")
computed <- !file.exists(INT_RDS)
if (!computed) {
    message("Loading cached ", INT_RDS,
            " — skipping compute, regenerating figures / applying celltype map")
    combined_CRPC <- readRDS(INT_RDS)
}

# ---------------------------------------------------------------------------
# Per-sample decontX: load → filter → decontaminate → replace counts
# ---------------------------------------------------------------------------
# decontX 는 내부적으로 quick clustering + UMAP 에 randomness 를 쓰므로 seed 고정
# (Harmony seed 이슈와 동일 맥락). 반환 SCE 의 decontXcounts / decontX_contamination
# 를 사용한다.
run_decontX_sample <- function(data_dir, project) {
    raw <- Read10X(data.dir = data_dir)
    seu <- CreateSeuratObject(counts = raw, project = project,
                              min.cells = 3, min.features = 200)

    counts <- GetAssayData(seu, assay = "RNA", layer = "counts")
    sce <- SingleCellExperiment(assays = list(counts = counts))
    set.seed(42)
    sce <- decontX(sce)

    # --- per-sample diagnostics ---
    # contamination on decontX's own UMAP
    p_cont <- tryCatch(
        plotDecontXContamination(sce) +
            ggtitle(paste0(project, " — decontX contamination")),
        error = function(e) NULL
    )
    if (!is.null(p_cont)) {
        ggsave(file.path(DX_DIR, paste0("contamination_UMAP_", project, ".png")),
               plot = p_cont, width = 8, height = 6, bg = "white")
    }
    # marker leakage before/after decontamination (lineage markers that should
    # not co-express across clusters — ambient RNA is the usual cause)
    marker_leak <- list(
        Epithelial  = c("EPCAM", "KRT8", "KRT18"),
        Tcell       = c("CD3D", "CD3E", "TRAC"),
        Myeloid     = c("LYZ", "CD68", "C1QA"),
        Endothelial = c("PECAM1", "VWF"),
        Fibroblast  = c("DCN", "LUM", "COL1A1")
    )
    marker_leak <- lapply(marker_leak, function(g) intersect(g, rownames(sce)))
    marker_leak <- marker_leak[lengths(marker_leak) > 0]
    p_mk <- tryCatch({
        sce_ln <- scater::logNormCounts(sce, assay.type = "counts")
        plotDecontXMarkerPercentage(sce_ln, markers = marker_leak,
                                    groupClusters = NULL) +
            ggtitle(paste0(project, " — marker %: raw vs decontaminated"))
    }, error = function(e) NULL)
    if (!is.null(p_mk)) {
        ggsave(file.path(DX_DIR, paste0("marker_percentage_", project, ".png")),
               plot = p_mk, width = 10, height = 6, bg = "white")
    }

    # --- replace counts with rounded decontaminated matrix ---
    decont <- round(decontXcounts(sce))
    decont <- as(decont, "CsparseMatrix")
    seu[["RNA"]] <- CreateAssay5Object(counts = decont)
    # decontamination only removes counts → recompute QC metrics on new counts
    seu$decontX_contamination <- sce$decontX_contamination
    seu$nCount_RNA   <- colSums(decont)
    seu$nFeature_RNA <- colSums(decont > 0)
    seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

    message(sprintf(
        "  %s: %d cells, contamination mean=%.3f median=%.3f",
        project, ncol(seu),
        mean(seu$decontX_contamination), median(seu$decontX_contamination)))
    seu
}

if (computed) {

# Load + decontaminate 3 samples ----
s1 <- run_decontX_sample("Raw_data/CRPC1/filtered_feature_bc_matrix/", "CRPC1")
s2 <- run_decontX_sample("Raw_data/CRPC2/filtered_feature_bc_matrix/", "CRPC2")
s3 <- run_decontX_sample("Raw_data/CRPC3/filtered_feature_bc_matrix/", "CRPC3")

# Per-sample contamination summary (pre-QC) ----
contam_df <- dplyr::bind_rows(lapply(list(s1, s2, s3), function(x) {
    data.frame(orig.ident = x$orig.ident,
               decontX_contamination = x$decontX_contamination)
}))
contam_summary <- contam_df |>
    dplyr::group_by(orig.ident) |>
    dplyr::summarise(n = dplyr::n(),
                     mean = mean(decontX_contamination),
                     median = median(decontX_contamination),
                     q90 = quantile(decontX_contamination, 0.90),
                     frac_gt_0.2 = mean(decontX_contamination > 0.2),
                     .groups = "drop")
write.csv(contam_summary,
          file.path(DX_DIR, "contamination_summary_by_sample.csv"),
          row.names = FALSE)

patient_cols_dx <- utils_cb_palette(dplyr::n_distinct(contam_df$orig.ident))
p_dx_vln <- ggplot(contam_df, aes(orig.ident, decontX_contamination,
                                  fill = orig.ident)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.12, outlier.size = 0.2, fill = "white") +
    scale_fill_manual(values = patient_cols_dx) +
    labs(x = NULL, y = "decontX contamination fraction",
         title = "Estimated ambient RNA contamination per sample") +
    theme_bw() + theme(legend.position = "none")
ggsave(file.path(DX_DIR, "contamination_violin_by_sample.png"),
       plot = p_dx_vln, width = 8, height = 6, bg = "white")

# Doublet Removal (per-sample, before merge, on decontaminated counts) ----
run_scDblFinder <- function(seu) {
    sce <- as.SingleCellExperiment(seu)
    sce <- scDblFinder(sce, BPPARAM = bpp)
    seu$scDblFinder.class <- sce$scDblFinder.class
    seu$scDblFinder.score <- sce$scDblFinder.score
    seu <- subset(seu, subset = scDblFinder.class == "singlet")
    return(seu)
}
set.seed(42)
s1 <- run_scDblFinder(s1)
s2 <- run_scDblFinder(s2)
s3 <- run_scDblFinder(s3)

# Integration
combined_CRPC_raw <- merge(s1,
    y = list(s2, s3),
    add.cell.ids = c("P1", "P2", "P3"), # To avoid barcode collision
    project = "CRPC"
)

# QC ----
combined_CRPC_raw[["percent.mt"]] <- PercentageFeatureSet(combined_CRPC_raw, pattern = "^MT-")
patient_cols <- utils_cb_palette(dplyr::n_distinct(combined_CRPC_raw$orig.ident))
p <- VlnPlot(combined_CRPC_raw,
    features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
    group.by = "orig.ident", ncol = 3, pt.size = 0, cols = patient_cols
)
ggsave(file.path(QC_DIR, "mt_by_patient_VlnPlot.png"), plot = p, width = 12, height = 8, bg = "white")

# decontX contamination distribution after doublet removal, by patient
p <- VlnPlot(combined_CRPC_raw, features = "decontX_contamination",
    group.by = "orig.ident", pt.size = 0, cols = patient_cols) + NoLegend()
ggsave(file.path(QC_DIR, "decontX_contamination_VlnPlot.png"), plot = p, width = 8, height = 6, bg = "white")

# VlnPlot with threshold lines to verify QC cutoffs (same thresholds as original)
qc_thresholds <- list(
    nFeature_RNA = c(200, 9000),
    nCount_RNA   = NULL,
    percent.mt   = c(NA, 20)
)
qc_vln <- lapply(names(qc_thresholds), function(feat) {
    vp <- VlnPlot(combined_CRPC_raw,
        features = feat,
        group.by = "orig.ident", pt.size = 0, cols = patient_cols
    ) + NoLegend()
    th <- qc_thresholds[[feat]]
    if (!is.null(th)) {
        for (y in th) if (!is.na(y)) vp <- vp + geom_hline(yintercept = y, linetype = "dashed", color = "red")
    }
    vp
})
qc_vln_combined <- wrap_plots(qc_vln, ncol = 3)
ggsave(file.path(QC_DIR, "QC_VlnPlot_with_thresholds.png"), plot = qc_vln_combined, width = 14, height = 6, bg = "white")

# FeatureScatter plot - separate panels for each patient (P1, P2, P3)
split_obj <- SplitObject(combined_CRPC_raw, split.by = "orig.ident")
scatter_plots <- lapply(names(split_obj), function(sample) {
    p1 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "percent.mt") + ggtitle(sample)
    p2 <- FeatureScatter(split_obj[[sample]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + ggtitle(sample)
    p1 + p2
})
combined_scatter <- wrap_plots(scatter_plots, ncol = 1)
ggsave(file.path(QC_DIR, "FeatureScatter_by_patient.png"), plot = combined_scatter, width = 12, height = 15, bg = "white")

# Filter (same thresholds as original stage 01)
combined_CRPC <- subset(combined_CRPC_raw,
    subset = nFeature_RNA > 200 &
        nFeature_RNA < 9000 &
        percent.mt < 20
)

# Normalization(SCTransform) & PCA ----
combined_CRPC <- SCTransform(combined_CRPC, verbose = FALSE)
combined_CRPC <- RunPCA(combined_CRPC, verbose = FALSE)
ggsave(file.path(QC_DIR, "ElbowPlot.png"), plot = ElbowPlot(combined_CRPC, ndims = 50), width = 12, height = 15, bg = "white")

# Integration (Harmony) ----
# Anchor the RNG before Harmony — its k-means init consumes the ambient RNG
# stream. Makes the embedding, and therefore the broad clusters, stable.
set.seed(42)
combined_CRPC <- IntegrateLayers(
    object = combined_CRPC,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony",
    normalization.method = "SCT",
    verbose = FALSE
)
combined_CRPC[["RNA"]] <- JoinLayers(combined_CRPC[["RNA"]])

# Clustering(Harmony) ----
combined_CRPC <- FindNeighbors(combined_CRPC, reduction = "harmony", dims = 1:30)
combined_CRPC <- FindClusters(combined_CRPC, resolution = 0.5)
combined_CRPC <- RunUMAP(combined_CRPC, reduction = "harmony", dims = 1:30)

# Cell Cycle Identification ----
s.genes <- cc.genes.updated.2019$s.genes
g2m.genes <- cc.genes.updated.2019$g2m.genes
combined_CRPC <- CellCycleScoring(combined_CRPC, s.features = s.genes, g2m.features = g2m.genes)

# SingleR Annotation ----
ref_hpca <- celldex::HumanPrimaryCellAtlasData()
sce_CRPC <- as.SingleCellExperiment(combined_CRPC, assay = "SCT")
singler_pred <- SingleR(
    test = sce_CRPC,
    ref = ref_hpca,
    labels = ref_hpca$label.main,
    BPPARAM = bpp
)
combined_CRPC$SingleR <- singler_pred$labels
combined_CRPC$SingleR_pruned <- singler_pred$pruned.labels
write.csv(
    data.frame(
        cell = rownames(singler_pred),
        label = singler_pred$labels,
        pruned_label = singler_pred$pruned.labels
    ),
    file.path(OUT_DIR, "SingleR_cell_labels.csv"),
    row.names = FALSE
)
png(file.path(OUT_DIR, "SingleR_score_heatmap.png"), width = 1200, height = 900, bg = "white")
plotScoreHeatmap(singler_pred)
dev.off()
png(file.path(OUT_DIR, "SingleR_delta_distribution.png"), width = 1200, height = 900, bg = "white")
plotDeltaDistribution(singler_pred, ncol = 4)
dev.off()

# Checkpoint save (seurat_clusters + SingleR + contamination; celltype 미부여) ----
utils_save_all_markers(combined_CRPC, file.path(OUT_DIR, "all_markers.csv"))
saveRDS(combined_CRPC, INT_RDS)

}  # end if (computed)

# ===========================================================================
# Figures (both modes) ----
# ===========================================================================
p <- DimPlot(combined_CRPC, group.by = "Phase", pt.size = 0.3,
    cols = utils_cb_palette(dplyr::n_distinct(combined_CRPC$Phase)))
ggsave(file.path(OUT_DIR, "CellCycle_Phase_UMAP.png"), plot = p, width = 10, height = 7, bg = "white")

# Marker genes for lineage assignment (same panel as original)
marker_genes <- c(
    EPCAM  = "Epithelial",
    KRT8   = "Epithelial (luminal)",
    AR     = "Luminal / AR pathway",
    PTPRC  = "Immune",
    VIM    = "Mesenchymal",
    PECAM1 = "Endothelial",
    ACTA2  = "Smooth muscle (ACTA2)",
    MYH11  = "Smooth muscle (MYH11)",
    SYP    = "Neuroendocrine",
    COL1A1 = "Fibroblast",
    DCN    = "Fibroblast / stromal",
    TPSAB1 = "Mast cell"
)
fp <- FeaturePlot(combined_CRPC, features = names(marker_genes),
    ncol = 4, pt.size = 0.1, combine = FALSE)
fp <- lapply(seq_along(fp), function(i) {
    fp[[i]] + ggtitle(bquote(.(marker_genes[[i]]) ~ "—" ~ italic(.(names(marker_genes)[i]))))
})
p <- wrap_plots(fp, ncol = 4)
ggsave(file.path(OUT_DIR, "FeaturePlot_markers.png"), plot = p, width = 20, height = 19, bg = "white")

# Cluster UMAP (번호 확인용 — 새 clustering 이므로 원본과 cluster ID 다름)
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "seurat_clusters", label = TRUE,
    cols = utils_cb_palette(nlevels(factor(combined_CRPC$seurat_clusters))))
ggsave(file.path(OUT_DIR, "Cluster_UMAP_integrated.png"), plot = p, width = 15, height = 15, bg = "white")

p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "seurat_clusters",
    split.by = "orig.ident", label = TRUE,
    cols = utils_cb_palette(nlevels(factor(combined_CRPC$seurat_clusters))))
ggsave(file.path(OUT_DIR, "Cluster_UMAP_by_patient.png"), plot = p, width = 24, height = 15, bg = "white")

# decontX contamination on the integrated UMAP (continuous → viridis)
p <- FeaturePlot(combined_CRPC, features = "decontX_contamination", pt.size = 0.2) +
    viridis::scale_color_viridis(option = "viridis") +
    ggtitle("decontX contamination (integrated UMAP)")
ggsave(file.path(OUT_DIR, "decontX_contamination_UMAP.png"), plot = p, width = 10, height = 8, bg = "white")

# SingleR UMAPs
p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "SingleR_pruned",
    label = TRUE, repel = TRUE,
    cols = utils_cb_palette(length(unique(na.omit(combined_CRPC$SingleR_pruned)))))
ggsave(file.path(OUT_DIR, "SingleR_UMAP.png"), plot = p, width = 15, height = 15, bg = "white")

p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "SingleR",
    split.by = "orig.ident", label = TRUE, repel = TRUE,
    cols = utils_cb_palette(length(unique(na.omit(combined_CRPC$SingleR)))))
ggsave(file.path(OUT_DIR, "SingleR_UMAP_by_patient.png"), plot = p, width = 24, height = 15, bg = "white")

# ===========================================================================
# Celltype assignment ----
# ===========================================================================
# ⚠️ decontX 재클러스터링으로 cluster 정수 ID 가 원본과 다르다. 아래 맵은
# all_markers.csv + FeaturePlot_markers.png + SingleR (majority label) 검토로
# 도출한 res-0.5 (19 clusters) lineage 배정.
#   Epithelial : 0,2,3,4,5,7,9,10,15,16,18  (EPCAM/KRT8/18/19 + SingleR Epithelial)
#     - cl0=basal(KRT14/COL17A1), cl5=분비성(CFTR/SLC4A4), cl10=TP63,
#       cl15=luminal/ARPC(AR/KLK3 최고), cl16=club(SCGB2A1/SPDEF)
#   Fibroblast : 1   (DCN/LUM/COL1A1)
#   Endothelial: 8   (PECAM1/VWF)
#   Smooth muscle cells: 11(CARMN/DMD/KCNMA1), 12(ACTA2/MYH11)
#     - cl11 은 SingleR 이 Epi/SM 로 갈렸으나 top 마커가 CARMN/DMD → 평활근(Epi=ambient)
#   T/NK cells : 6(CD3D/NKG7), 14(THEMIS/BCL11B/SKAP1 → T 계열)
#   Phagocytes : 13  (LYZ/CD68/C1QA)
#   Mast cells : 17  (TPSAB1/CPA3)
# 비어 있으면 celltype 미부여 + checkpoint 만 저장.
celltype_map <- c(
    `0`  = "Epithelial",          `1`  = "Fibroblast",
    `2`  = "Epithelial",          `3`  = "Epithelial",
    `4`  = "Epithelial",          `5`  = "Epithelial",
    `6`  = "T/NK cells",          `7`  = "Epithelial",
    `8`  = "Endothelial",         `9`  = "Epithelial",
    `10` = "Epithelial",          `11` = "Smooth muscle cells",
    `12` = "Smooth muscle cells", `13` = "Phagocytes",
    `14` = "T/NK cells",          `15` = "Epithelial",
    `16` = "Epithelial",          `17` = "Mast cells",
    `18` = "Epithelial"
)

if (length(celltype_map) > 0) {
    stopifnot(all(levels(factor(combined_CRPC$seurat_clusters)) %in% names(celltype_map)))
    # named-vector lookup (RenameIdents !!! splice 대신 견고한 방식)
    combined_CRPC$celltype <- factor(
        unname(celltype_map[as.character(combined_CRPC$seurat_clusters)])
    )
    Idents(combined_CRPC) <- "celltype"

    # =======================================================================
    # Stromal refinement ----
    #   위 coarse map 의 'Smooth muscle cells'(cl 11+12) 는 pericyte / vascular SMC /
    #   fibroblast / Schwann / 저품질 mesenchyme 의 lump 이고, 'Fibroblast'(cl 1) 에도
    #   SMC 가 섞여 있다 (2026-08 외부 CRPC1 재주석 비교 + 기질 subclustering).
    #   기질은 전체의 6% 라 mural 축이 공용 PCA 에 실려 있지 않아 전체 해상도를 올려도
    #   갈리지 않는다 → 상피(stage 02)와 같은 방식으로 기질만 떼어 재클러스터링하고
    #   subcluster → 라벨 맵으로 celltype 을 갱신한다. 진단 그림은 OUT_DIR/Stromal/.
    # =======================================================================
    STR_DIR <- file.path(OUT_DIR, "Stromal")
    dir.create(STR_DIR, showWarnings = FALSE, recursive = TRUE)
    str_res_vec <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)

    str_obj <- subset(combined_CRPC,
        subset = celltype %in% c("Fibroblast", "Smooth muscle cells"))
    message("Stromal subset: ", ncol(str_obj), " cells")
    str_obj$parent_celltype <- as.character(str_obj$celltype)
    str_obj$parent_cluster  <- as.character(str_obj$seurat_clusters)

    DefaultAssay(str_obj) <- "RNA"
    if ("SCT" %in% names(str_obj@assays)) str_obj[["SCT"]] <- NULL
    str_obj[["RNA"]] <- JoinLayers(str_obj[["RNA"]])
    str_obj[["RNA"]] <- split(str_obj[["RNA"]], f = str_obj$orig.ident)
    str_obj@meta.data <- dplyr::select(str_obj@meta.data,
        -matches("_snn_res\\."), -any_of("seurat_clusters"))

    str_obj <- SCTransform(str_obj, verbose = FALSE)
    str_obj <- RunPCA(str_obj, verbose = FALSE)
    set.seed(42)
    str_obj <- IntegrateLayers(
        object = str_obj, method = HarmonyIntegration,
        orig.reduction = "pca", new.reduction = "harmony",
        normalization.method = "SCT", verbose = FALSE
    )
    str_obj[["RNA"]] <- JoinLayers(str_obj[["RNA"]])
    str_obj <- FindNeighbors(str_obj, reduction = "harmony", dims = 1:30)
    str_obj <- FindClusters(str_obj, resolution = str_res_vec)
    str_obj$stromal_subcluster <- factor(as.character(str_obj$SCT_snn_res.0.3),
        levels = as.character(sort(unique(as.integer(as.character(str_obj$SCT_snn_res.0.3))))))
    str_obj <- RunUMAP(str_obj, reduction = "harmony", dims = 1:30,
        n.neighbors = 20, min.dist = 0.1, spread = 4.0)

    # res-0.3 subcluster → 라벨 (마커 + 외부 CRPC1 재주석 교차검증으로 확정)
    #   c6  pericyte (HIGD1B/NDUFA4L2/GJA4/ADGRF5)
    #   c4  vascular SMC (PLN/RERGL/CASQ2, MYH11+ 99.8%, NOTCH3+ 65%) → Pericyte 로 통합
    #   c9  분화 SM (DES/ACTG2/MYLK/CNN1),  c8  MYH11+ 저복잡도 SM
    #   c10 Schwann (PLP1/S100B/SOX10/CDH19)
    #   c0 APOD+ / c1 ASPN·OGN·C7 / c3 PENK·PTGDS / c7 CA3·APOE fibroblast
    #   c2·c5·c8 은 부모 cl 11 유래 저복잡도(nCount 중앙값 1.2k~2.2k) — 계통 힌트대로 배정
    # ⚠️ cluster 정수 ID 는 rerun 간 불안정할 수 있다. ID 가 맵에 없으면 아래 stopifnot 이
    #    멈추므로 Stromal/ 진단 그림을 보고 맵을 다시 채운다 (조용한 오배정 방지).
    stromal_map <- c(
        `0` = "Fibroblast", `1` = "Fibroblast", `2` = "Fibroblast",
        `3` = "Fibroblast", `7` = "Fibroblast",
        `4` = "Pericyte",   `5` = "Pericyte",   `6` = "Pericyte",
        `8` = "Smooth muscle cells", `9` = "Smooth muscle cells",
        `10` = "Schwann cell"
    )
    str_sub <- as.character(str_obj$stromal_subcluster)
    stopifnot(all(str_sub %in% names(stromal_map)))
    str_obj$celltype <- factor(unname(stromal_map[str_sub]),
        levels = c("Fibroblast", "Pericyte", "Smooth muscle cells", "Schwann cell"))
    Idents(str_obj) <- "celltype"

    # canonical celltype 갱신 (기질 세포만 바뀜)
    ct <- as.character(combined_CRPC$celltype)
    ct[match(colnames(str_obj), colnames(combined_CRPC))] <- as.character(str_obj$celltype)
    celltype_levels <- c("Epithelial", "Fibroblast", "Pericyte", "Smooth muscle cells",
                         "Schwann cell", "Endothelial", "T/NK cells", "Phagocytes", "Mast cells")
    stopifnot(all(ct %in% celltype_levels))
    combined_CRPC$celltype <- factor(ct, levels = celltype_levels)
    Idents(combined_CRPC) <- "celltype"
    message("Stromal refinement 적용:")
    print(table(str_obj$parent_celltype, str_obj$celltype))

    # --- Stromal 진단 출력 --------------------------------------------------
    n_sub <- nlevels(str_obj$stromal_subcluster)
    p <- DimPlot(str_obj, group.by = "stromal_subcluster", label = TRUE, pt.size = 0.4,
                 cols = utils_cb_palette(n_sub)) + ggtitle("Stromal subclusters (res 0.3)")
    ggsave(file.path(STR_DIR, "UMAP_subcluster.png"), plot = p, width = 10, height = 8, bg = "white")
    p <- DimPlot(str_obj, group.by = "celltype", label = TRUE, pt.size = 0.4,
                 cols = utils_cb_palette(4)) + ggtitle("Stromal cell types")
    ggsave(file.path(STR_DIR, "UMAP_celltype.png"), plot = p, width = 10, height = 8, bg = "white")
    p <- DimPlot(str_obj, group.by = "celltype", split.by = "orig.ident", label = TRUE,
                 pt.size = 0.4, cols = utils_cb_palette(4))
    ggsave(file.path(STR_DIR, "UMAP_celltype_by_patient.png"), plot = p, width = 24, height = 8, bg = "white")
    p <- DimPlot(str_obj, group.by = "parent_cluster", label = TRUE, pt.size = 0.4,
                 cols = utils_cb_palette(dplyr::n_distinct(str_obj$parent_cluster))) +
        ggtitle("Parent res-0.5 cluster (1 = Fibroblast, 11/12 = coarse 'SMC')")
    ggsave(file.path(STR_DIR, "UMAP_parent_cluster.png"), plot = p, width = 10, height = 8, bg = "white")
    umap_list <- lapply(str_res_vec, function(r) {
        rcol <- paste0("SCT_snn_res.", r)
        DimPlot(str_obj, group.by = rcol, label = TRUE, pt.size = 0.3,
                cols = utils_cb_palette(dplyr::n_distinct(str_obj@meta.data[[rcol]]))) +
            ggtitle(paste0("res = ", r)) + NoLegend()
    })
    ggsave(file.path(STR_DIR, "UMAP_multires.png"), plot = wrap_plots(umap_list, ncol = 3),
           width = 18, height = 14, dpi = 200, bg = "white")

    p <- VlnPlot(str_obj, group.by = "stromal_subcluster", pt.size = 0, ncol = 2,
                 features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination"),
                 cols = utils_cb_palette(n_sub))
    ggsave(file.path(STR_DIR, "QC_violin_by_subcluster.png"), plot = p, width = 16, height = 10, bg = "white")

    stromal_markers <- list(
        Fibroblast    = c("DCN", "LUM", "PDGFRA", "COL1A1", "FBLN1", "APOD"),
        Pericyte      = c("RGS5", "NOTCH3", "HIGD1B", "NDUFA4L2", "KCNJ8", "GJA4"),
        `Smooth muscle cells` = c("MYH11", "ACTG2", "DES", "MYLK", "CNN1", "ACTA2"),
        Schwann       = c("PLP1", "S100B", "MPZ", "SOX10", "CDH19"),
        Contamination = c("PECAM1", "PTPRC", "EPCAM", "KRT8")
    )
    p <- DotPlot(str_obj, features = stromal_markers, group.by = "stromal_subcluster",
                 cluster.idents = FALSE) + RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        labs(title = "Stromal markers by subcluster")
    ggsave(file.path(STR_DIR, "DotPlot_by_subcluster.png"), plot = p, width = 20, height = 8, bg = "white")
    p <- DotPlot(str_obj, features = stromal_markers, group.by = "celltype",
                 cluster.idents = FALSE) + RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        labs(title = "Stromal markers by cell type")
    ggsave(file.path(STR_DIR, "DotPlot_by_celltype.png"), plot = p, width = 18, height = 6, bg = "white")

    write.csv(as.data.frame.matrix(table(parent = str_obj$parent_cluster,
                                         subcluster = str_obj$stromal_subcluster)),
              file.path(STR_DIR, "parent_vs_subcluster.csv"))
    write.csv(as.data.frame.matrix(table(str_obj$stromal_subcluster, str_obj$celltype)),
              file.path(STR_DIR, "subcluster_to_celltype.csv"))
    write.csv(as.data.frame.matrix(table(str_obj$celltype, str_obj$orig.ident)),
              file.path(STR_DIR, "celltype_by_patient.csv"))
    qc_tab <- str_obj@meta.data %>%
        group_by(subcluster = stromal_subcluster) %>%
        summarise(n = n(), median_nCount = median(nCount_RNA),
                  median_nFeature = median(nFeature_RNA),
                  median_mt = round(median(percent.mt), 2),
                  median_decontX = round(median(decontX_contamination), 3), .groups = "drop")
    write.csv(qc_tab, file.path(STR_DIR, "subcluster_QC_summary.csv"), row.names = FALSE)

    Idents(str_obj) <- "stromal_subcluster"
    utils_save_all_markers(str_obj, file.path(STR_DIR, "all_markers_by_subcluster.csv"))
    Idents(str_obj) <- "celltype"
    saveRDS(str_obj, file.path(STR_DIR, "stromal_annotated.rds"))
    rm(str_obj); gc()

    # --- 최종 celltype 그림 ---------------------------------------------------
    n_ct <- nlevels(combined_CRPC$celltype)
    p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", label = TRUE,
        cols = utils_cb_palette(n_ct))
    ggsave(file.path(OUT_DIR, "Labelled_UMAP_integrated.png"), plot = p, width = 15, height = 15, bg = "white")

    p <- DimPlot(combined_CRPC, reduction = "umap", group.by = "celltype", split.by = "orig.ident", label = TRUE,
        cols = utils_cb_palette(n_ct))
    ggsave(file.path(OUT_DIR, "Labelled_UMAP_by_patient.png"), plot = p, width = 24, height = 15, bg = "white")

    comp <- as.data.frame(table(celltype = combined_CRPC$celltype, patient = combined_CRPC$orig.ident))
    p <- ggplot(comp, aes(x = patient, y = Freq, fill = celltype)) +
        geom_col(position = "fill") +
        scale_fill_manual(values = utils_cb_palette(n_ct)) +
        labs(y = "fraction of cells", x = NULL, title = "Cell type composition by patient") +
        theme_classic(base_size = 14)
    ggsave(file.path(OUT_DIR, "Celltype_composition_by_patient.png"), plot = p, width = 8, height = 6, bg = "white")
    write.csv(as.data.frame.matrix(table(combined_CRPC$celltype, combined_CRPC$orig.ident)),
              file.path(OUT_DIR, "celltype_by_patient.csv"))

    saveRDS(combined_CRPC, INT_RDS)   # persist celltype for stage 02 subset
    message("celltype 부여 완료 (stromal refinement 포함) → ", INT_RDS)
    print(table(combined_CRPC$celltype))
} else {
    message("celltype_map 비어 있음 — cluster/marker 검토 후 채우고 재실행하면 ",
            "reload 모드로 celltype 부여 + 저장")
}
