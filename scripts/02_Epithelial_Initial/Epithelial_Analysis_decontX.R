# Stage 2 (decontX variant) — Epithelial initial reclustering
#
# 원본 scripts/02_Epithelial_Initial/Epithelial_Analysis.R 와 동일한 방법론이되,
# 입력을 decontX 파이프라인의 stage-01 객체로 바꾸고 출력을 별도 경로로 낸다.
#   입력: Results/01_Integrated_decontX/combined_CRPC_decontX.rds  (celltype 부여됨)
#   출력: Results/02_Epithelial_Initial_decontX/
# 파라미터(res_vec, Harmony seed=42, UMAP 설정, 마커 패널)는 원본과 동일하게 유지.

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

IN_RDS  <- "Results/01_Integrated_decontX/combined_CRPC_decontX.rds"
OUT_DIR <- "Results/02_Epithelial_Initial_decontX"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

EPI_RDS  <- file.path(OUT_DIR, "epi_clustered.rds")
computed <- !file.exists(EPI_RDS)

# Resolutions swept for clustree (원본과 동일)
res_vec <- c(0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2)

if (computed) {

# ============================================================
# Load: epithelial subset ----
# ============================================================
combined_CRPC <- readRDS(IN_RDS)

epi <- subset(combined_CRPC, subset = celltype == "Epithelial")
message("Epithelial subset: ", ncol(epi), " cells")

rm(combined_CRPC)
gc()

# Reset assay state (drop inherited SCT, re-split RNA by patient for Harmony)
DefaultAssay(epi) <- "RNA"
if ("SCT" %in% names(epi@assays)) epi[["SCT"]] <- NULL
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)

# Drop inherited *_snn_res.* / seurat_clusters columns from the parent object
epi@meta.data <- dplyr::select(epi@meta.data,
    -matches("_snn_res\\."), -any_of("seurat_clusters"))

# ============================================================
# SCTransform + PCA + Clustering ----
# ============================================================
epi <- SCTransform(epi, verbose = FALSE)
epi <- RunPCA(epi, verbose = FALSE)

# Anchor the RNG before Harmony (원본과 동일 — k-means init 이 ambient RNG 사용)
set.seed(42)
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
    plot = ElbowPlot(epi, ndims = 50), width = 12, height = 8, bg = "white"
)

epi <- FindNeighbors(epi, reduction = "harmony", dims = 1:30)
epi <- FindClusters(epi, resolution = res_vec)

# Curate small fragment clusters at res 0.4 into their parent. clustree.png 검토 후
# 채운다 — cluster 정수 ID 는 rerun 간 불안정하므로 하드코딩 map 은 조용히 오병합.
# 비어 있으면 raw res-0.4 클러스터 유지.
merge_map <- character()
curated <- if (length(merge_map) > 0) {
    dplyr::recode(as.character(epi$SCT_snn_res.0.4), !!!merge_map)
} else {
    as.character(epi$SCT_snn_res.0.4)
}
epi$seurat_clusters <- factor(curated,
    levels = as.character(sort(unique(as.integer(curated)))))
Idents(epi) <- "seurat_clusters"

epi <- RunUMAP(epi,
    reduction = "harmony", dims = 1:30,
    n.neighbors = 20, min.dist = 0.1, spread = 4.0
)

} else {
    message("Loading cached ", EPI_RDS,
            " — skipping compute, regenerating figures")
    epi <- readRDS(EPI_RDS)
    Idents(epi) <- "seurat_clusters"
}

# clustree: cluster stability tree across resolutions
p_tree <- clustree(epi, prefix = "SCT_snn_res.")
ggsave(file.path(OUT_DIR, "clustree.png"),
    plot = p_tree, width = 12, height = 14, dpi = 200, bg = "white"
)

# UMAPs at each resolution for side-by-side comparison
umap_list <- lapply(res_vec, function(r) {
    rcol <- paste0("SCT_snn_res.", r)
    DimPlot(epi, group.by = rcol, label = TRUE, pt.size = 0.2,
            cols = utils_cb_palette(dplyr::n_distinct(epi@meta.data[[rcol]]))) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
p_umap_grid <- wrap_plots(umap_list, ncol = 3)
ggsave(file.path(OUT_DIR, "UMAP_multires.png"),
    plot = p_umap_grid, width = 18, height = 14, dpi = 200, bg = "white"
)

p1 <- DimPlot(epi, label = TRUE, pt.size = 0.3,
              cols = utils_cb_palette(nlevels(factor(epi$seurat_clusters)))) +
    ggtitle("Epithelial clusters (decontX)")
p2 <- DimPlot(epi, group.by = "orig.ident", pt.size = 0.3,
              cols = utils_cb_palette(dplyr::n_distinct(epi$orig.ident))) +
    ggtitle("By patient")
ggsave(file.path(OUT_DIR, "UMAP_cluster.png"), plot = p1, width = 10, height = 8, bg = "white")
ggsave(file.path(OUT_DIR, "UMAP_patient.png"), plot = p2, width = 10, height = 8, bg = "white")

# Patient-split UMAP — patient-specific clusters are tumor candidates
p <- DimPlot(epi, split.by = "orig.ident", label = TRUE, pt.size = 0.3,
             cols = utils_cb_palette(nlevels(factor(epi$seurat_clusters))))
ggsave(file.path(OUT_DIR, "UMAP_split_by_patient.png"),
    plot = p, width = 24, height = 8, bg = "white"
)

# ============================================================
# Cell cycle scoring ----
# ============================================================
if (computed) {
    s.genes <- cc.genes.updated.2019$s.genes
    g2m.genes <- cc.genes.updated.2019$g2m.genes
    epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, nbin = 12)
}

p <- DimPlot(epi, group.by = "Phase", pt.size = 0.3,
             cols = utils_cb_palette(dplyr::n_distinct(epi$Phase)))
ggsave(file.path(OUT_DIR, "CellCycle_UMAP.png"), plot = p, width = 10, height = 7, bg = "white")

# ============================================================
# Marker visualization (원본과 동일 패널) ----
# ============================================================
luminal_basal_club <- c(
    "FOXA1", "HOXB13", "NKX3-1", "KLK3", "TOP2A", "MKI67", # Luminal
    "KRT5", "KRT14", "TP63", # Basal
    "MMP7", "WFDC2" # Club
)
dnpc_markers <- c("CHD7", "MYC", "KMT2C", "KRT7", "SOX2", "SYP", "AR")

p <- FeaturePlot(epi, features = luminal_basal_club, ncol = 4, pt.size = 0.1)
ggsave(file.path(OUT_DIR, "Luminal_Basal_Club_FeaturePlot.png"),
    plot = p, width = 20, height = 16, bg = "white"
)
p <- DotPlot(epi, features = luminal_basal_club) + RotatedAxis()
ggsave(file.path(OUT_DIR, "Luminal_Basal_Club_DotPlot.png"),
    plot = p, width = 14, height = 8, bg = "white"
)
p <- FeaturePlot(epi, features = dnpc_markers, ncol = 4, pt.size = 0.1)
ggsave(file.path(OUT_DIR, "DNPC_FeaturePlot.png"),
    plot = p, width = 20, height = 16, bg = "white"
)
p <- DotPlot(epi, features = dnpc_markers) + RotatedAxis()
ggsave(file.path(OUT_DIR, "DNPC_DotPlot.png"),
    plot = p, width = 12, height = 8, bg = "white"
)

# CRPC subtype panel: ARPC / DNPC-FGF / DNPC-KRT7 / NEPC
subtype_markers <- list(
    ARPC       = c("AR", "KLK3", "FOLH1"),
    `DNPC-FGF` = c("FGF8", "FGFR1", "CCL2", "DKK1"),
    `DNPC-KRT7` = c("KRT7", "SOX2", "FOXA2"),
    NEPC       = c("CHGA", "SYP", "PROX1")
)
p <- DotPlot(epi, features = subtype_markers, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "CRPC subtype markers (ARPC / DNPC-FGF / DNPC-KRT7 / NEPC)")
ggsave(file.path(OUT_DIR, "Subtype_DotPlot.png"),
    plot = p, width = 14, height = 8, bg = "white"
)

# Hillock / club-hillock verification panel
hillock_markers <- list(
    Hillock_core    = c("KRT13", "KRT4", "AQP3", "LY6D", "PSCA", "LYPD3"),
    Basal_axis      = c("KRT5", "KRT14", "TP63"),
    Squamous_diff   = c("SERPINB5", "SERPINB1", "ZNF750", "S100A2"),
    Club            = c("MMP7", "WFDC2", "PIGR", "LCN2"),
    Luminal_ref     = c("KRT8", "KRT18", "AR", "KLK3")
)
p <- DotPlot(epi, features = hillock_markers, cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    labs(title = "Hillock / club-hillock verification")
ggsave(file.path(OUT_DIR, "Hillock_DotPlot.png"),
    plot = p, width = 16, height = 8, bg = "white"
)

# ============================================================
# Find all markers + save (compute mode only) ----
# ============================================================
if (computed) {
    utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(epi, EPI_RDS)
}

message("Stage 2 (decontX) 결과 → ", OUT_DIR)
