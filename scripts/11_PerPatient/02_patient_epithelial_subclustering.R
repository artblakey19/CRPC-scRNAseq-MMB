#!/usr/bin/env Rscript
# Stage 11-02 — Per-patient epithelial subclustering + annotation (no integration)
#
# 입력: Results/11_PerPatient/<PATIENT>/01_Major/<PATIENT>_major.rds (celltype 부여됨)
#   → celltype == "Epithelial" (+ "Neuroendocrine" 이 있으면 포함) 만 subset
# 파이프라인: SCTransform → PCA → FindNeighbors(dims 1:30) → FindClusters(다중 res)
#   → UMAP(stage-02 canonical 파라미터 nn20/md0.1/spread4) → CellCycle → markers
#   → 상피 subtype 마커 패널 (luminal/basal/club/hillock/OE/ionocyte/NE/ARPC/DNPC)
#   → 통합 상피 annotation (Results/05_Epithelial_Downstream/epi_annotated.rds) 과 교차표
#   → annotation_map (환자별, 검토 후 채움) → 최종 저장
#
# 실행 (프로젝트 루트): Rscript scripts/11_PerPatient/02_patient_epithelial_subclustering.R CRPC2
# 출력: Results/11_PerPatient/<PATIENT>/02_Epithelial/

args <- commandArgs(trailingOnly = TRUE)
PATIENT <- if (length(args) >= 1) args[1] else "CRPC2"
stopifnot(PATIENT %in% c("CRPC1", "CRPC2", "CRPC3"))

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
    library(glmGamPoi)
    library(clustree)
    library(future)
})
N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)

source("scripts/00_utils/scRNA_utils.R")

IN_RDS   <- file.path("Results/11_PerPatient", PATIENT, "01_Major", paste0(PATIENT, "_major.rds"))
EPI_INT  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"   # 통합 상피 annotation (비교용)
OUT_DIR  <- file.path("Results/11_PerPatient", PATIENT, "02_Epithelial")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OBJ_RDS  <- file.path(OUT_DIR, paste0(PATIENT, "_epi.rds"))

RES_VEC <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
RES_EPI <- list(CRPC1 = "0.4", CRPC2 = "0.4", CRPC3 = "0.4")[[PATIENT]]

computed <- !file.exists(OBJ_RDS)

if (computed) {
    major <- readRDS(IN_RDS)
    stopifnot("celltype" %in% colnames(major@meta.data))
    epi <- subset(major, subset = celltype %in% c("Epithelial", "Neuroendocrine"))
    message(PATIENT, " epithelial: ", ncol(epi), " cells")
    rm(major); gc()

    epi$major_celltype <- as.character(epi$celltype)
    epi$major_cluster  <- as.character(epi$seurat_clusters)
    DefaultAssay(epi) <- "RNA"
    epi[["SCT"]] <- NULL
    for (rd in names(epi@reductions)) epi[[rd]] <- NULL
    epi@meta.data <- dplyr::select(epi@meta.data, -matches("_snn_res\\."), -any_of(c("seurat_clusters", "celltype")))

    # 통합 상피 annotation 붙이기 (같은 barcode; 통합 상피에 없던 셀은 NA)
    epi_int <- readRDS(EPI_INT)
    int_ann <- setNames(as.character(epi_int$annotation), colnames(epi_int))
    epi$integrated_epi_annotation <- unname(int_ann[colnames(epi)])
    rm(epi_int); gc()

    epi <- SCTransform(epi, verbose = FALSE)
    epi <- RunPCA(epi, verbose = FALSE)
    ggsave(file.path(OUT_DIR, "ElbowPlot.png"), plot = ElbowPlot(epi, ndims = 50),
           width = 10, height = 6, bg = "white")
    epi <- FindNeighbors(epi, reduction = "pca", dims = 1:30)
    epi <- FindClusters(epi, resolution = RES_VEC, random.seed = 42)
    cl <- as.character(epi@meta.data[[paste0("SCT_snn_res.", RES_EPI)]])
    epi$seurat_clusters <- factor(cl, levels = as.character(sort(unique(as.integer(cl)))))
    Idents(epi) <- "seurat_clusters"
    set.seed(42)
    epi <- RunUMAP(epi, reduction = "pca", dims = 1:30, seed.use = 42,
                   n.neighbors = 20, min.dist = 0.1, spread = 4.0)
    epi <- CellCycleScoring(epi, s.features = cc.genes.updated.2019$s.genes,
                            g2m.features = cc.genes.updated.2019$g2m.genes, nbin = 12)
    utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(epi, OBJ_RDS)
} else {
    message("Loading cached ", OBJ_RDS, " — figures / annotation map only")
    epi <- readRDS(OBJ_RDS)
    Idents(epi) <- "seurat_clusters"
}

# ======================================================================
# Diagnostics
# ======================================================================
n_cl <- nlevels(epi$seurat_clusters)
ggsave(file.path(OUT_DIR, "clustree.png"), plot = clustree(epi, prefix = "SCT_snn_res."),
       width = 12, height = 14, dpi = 200, bg = "white")
umap_list <- lapply(RES_VEC, function(r) {
    rcol <- paste0("SCT_snn_res.", r)
    DimPlot(epi, group.by = rcol, label = TRUE, pt.size = 0.2,
            cols = utils_cb_palette(dplyr::n_distinct(epi@meta.data[[rcol]]))) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
ggsave(file.path(OUT_DIR, "UMAP_multires.png"), plot = wrap_plots(umap_list, ncol = 4),
       width = 24, height = 12, dpi = 200, bg = "white")

p <- DimPlot(epi, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3, cols = utils_cb_palette(n_cl)) +
    ggtitle(paste0(PATIENT, " alone — epithelial clusters (res ", RES_EPI, ")"))
ggsave(file.path(OUT_DIR, "UMAP_cluster.png"), plot = p, width = 10, height = 8, bg = "white")
p <- DimPlot(epi, group.by = "integrated_epi_annotation", label = TRUE, repel = TRUE, pt.size = 0.3,
             cols = utils_cb_palette(length(unique(na.omit(epi$integrated_epi_annotation))))) +
    ggtitle("Integrated epithelial annotation on patient-alone UMAP")
ggsave(file.path(OUT_DIR, "UMAP_integrated_annotation.png"), plot = p, width = 11, height = 8, bg = "white")
p <- DimPlot(epi, group.by = "Phase", pt.size = 0.3, cols = utils_cb_palette(dplyr::n_distinct(epi$Phase)))
ggsave(file.path(OUT_DIR, "UMAP_CellCycle.png"), plot = p, width = 10, height = 7, bg = "white")
p <- VlnPlot(epi, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination"),
             group.by = "seurat_clusters", pt.size = 0, ncol = 2, cols = utils_cb_palette(n_cl))
ggsave(file.path(OUT_DIR, "QC_violin_by_cluster.png"), plot = p, width = 16, height = 10, bg = "white")

epi_markers <- list(
    Luminal_AR   = c("AR", "KLK3", "KLK2", "NKX3-1", "FOLH1", "ACP3", "TMPRSS2", "HOXB13", "FOXA1"),
    Luminal_pan  = c("KRT8", "KRT18", "CD24"),
    Basal        = c("KRT5", "KRT14", "KRT15", "TP63", "COL17A1"),
    Club         = c("MMP7", "WFDC2", "PIGR", "SCGB3A1", "SCGB1A1", "LCN2", "LTF"),
    Hillock      = c("KRT13", "KRT4", "AQP3", "LY6D", "PSCA", "S100A2", "KRT6A"),
    OE_CFTR      = c("CFTR", "SLC4A4", "MUC5B", "OLFM4"),
    Ionocyte     = c("FOXI1", "ATP6V1B1", "ASCL3"),
    NE           = c("CHGA", "CHGB", "SYP", "ENO2", "ASCL1", "INSM1", "NCAM1"),
    DNPC         = c("KRT7", "SOX2", "FOXA2", "FGFR1", "FGF8", "DKK1", "CHD7", "MYC"),
    Tumor_nonAR  = c("PCA3", "SCHLAP1", "GOLM1", "AMACR", "ERG"),
    Stress_AP1   = c("FOS", "JUN", "JUNB", "ATF3", "HSPA1A", "DDIT3"),
    Cycling      = c("MKI67", "TOP2A"),
    Squamous     = c("KRT1", "KRT10", "SPRR1B", "ZNF750"),
    Urothelial   = c("UPK1B", "UPK3A"),
    Contam       = c("PTPRC", "DCN", "PECAM1")
)
epi_markers <- lapply(epi_markers, function(g) intersect(g, rownames(epi)))
p <- DotPlot(epi, features = epi_markers, group.by = "seurat_clusters", cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    theme(axis.text.x = element_text(size = 7)) +
    labs(title = paste0(PATIENT, " alone — epithelial subtype markers by cluster"))
ggsave(file.path(OUT_DIR, "DotPlot_epi_markers_by_cluster.png"), plot = p, width = 30, height = 8, bg = "white")

fp_genes <- c("AR", "KLK3", "NKX3-1", "KRT8", "KRT5", "KRT14", "TP63", "MMP7", "PIGR",
              "KRT13", "KRT4", "CFTR", "FOXI1", "SYP", "CHGA", "KRT7", "SOX2", "PCA3",
              "AMACR", "FOS", "MKI67", "TOP2A", "S100A2", "UPK1B")
p <- FeaturePlot(epi, features = fp_genes, ncol = 6, pt.size = 0.1)
ggsave(file.path(OUT_DIR, "FeaturePlot_epi_markers.png"), plot = p, width = 30, height = 20, bg = "white")

write.csv(as.data.frame.matrix(table(cluster = epi$seurat_clusters,
                                     integrated = addNA(factor(epi$integrated_epi_annotation)))),
          file.path(OUT_DIR, "cluster_vs_integrated_annotation.csv"))
write.csv(as.data.frame.matrix(table(cluster = epi$seurat_clusters, major_cluster = epi$major_cluster)),
          file.path(OUT_DIR, "cluster_vs_major_cluster.csv"))
avg <- AverageExpression(epi, features = unlist(epi_markers), group.by = "seurat_clusters",
                         assay = "SCT", layer = "data")$SCT
write.csv(round(as.matrix(avg), 3), file.path(OUT_DIR, "epi_marker_avg_by_cluster.csv"))
pct <- sapply(levels(epi$seurat_clusters), function(k) {
    m <- GetAssayData(epi, assay = "SCT", layer = "data")[unlist(epi_markers), epi$seurat_clusters == k]
    Matrix::rowMeans(m > 0)
})
write.csv(round(pct, 3), file.path(OUT_DIR, "epi_marker_pct_by_cluster.csv"))
qc_tab <- epi@meta.data %>% group_by(cluster = seurat_clusters) %>%
    summarise(n = n(), median_nCount = median(nCount_RNA), median_nFeature = median(nFeature_RNA),
              median_mt = round(median(percent.mt), 2),
              median_decontX = round(median(decontX_contamination), 3),
              frac_cycling = round(mean(Phase != "G1"), 3), .groups = "drop")
write.csv(qc_tab, file.path(OUT_DIR, "cluster_QC_summary.csv"), row.names = FALSE)

# ======================================================================
# Annotation — 환자별 맵 (마커 DotPlot / all_markers / 교차표 검토 후 채움)
# ======================================================================
# 라벨은 마커 기반 서술형 (환자 단독 구조를 그대로 드러내기 위해 통합 BE1–6 번호 체계를 강제하지 않음).
# CRPC2 (res 0.4, 11 clusters; 2026-09-07 검토):
#   0 KRT5/KRT14/KRT15 basal + DKK1/CCL2/SFRP1 (통합 BE1/2/3)
#   1 저count(nFeature 1.4k, mt 10.6%) 핵 lncRNA (TALAM1/MAGI2/ZBTB20) — 통합 BE4/BE6 = 같은 low-count 상태
#   2 FOS/JUN/ATF3/DDIT3 AP-1 stress basal (통합 BE5)
#   3 KRT13/OLFM4/GABRP/FOXC1/AREG — basal→hillock 이행 (통합 BE1/BE2/Hillock1/2 혼재)
#   4 KRT13 76/KRT6A/SPRR1B/S100A8/KRT1 squamous hillock (통합 Hillock2)
#   5 저count(nFeature 1.1k) PPARG/BMPR1B/TMPRSS2 luminal-like (통합 OE 242/LE 74)
#   6 SCGB1A1/SCGB3A1/TFF1/PIGR/MMP7 club (+ACP3/KLK3 약함) (통합 Club/LE)
#   7 FOXI1/ATP6V0A4/ATP6V0D2/TMEM213 ionocyte 151셀 (통합 Ionocyte 148)
#   8 DSG3/DSG1/TMPRSS11D/SERPINB13/ESR1, TP63+ squamous (mt 15.6%)
#   9 AR 6.4/TMPRSS2 12/FOLH1/ACP3/KLK2 + PCA3(38%)/SCHLAP1/GOLM1/AMACR, KRT8 낮음 — AR+ 종양성 luminal (통합 LE(ARPC))
#   10 KIT/TFCP2L1/ATP6V0D2 + LAMA2/COL4A2, FOXI1 17% — ionocyte-like 36셀 (major cl18)
# CRPC3 (res 0.4, 8 clusters):
#   0 KRT5/KRT15 basal + DKK1 20/TNC/CCL2/S100A2 (통합 BE1/BE2/BE4)
#   1 KRT13 29/KRT6A/LY6D/OLFM4/SERPINB3 hillock (통합 Hillock1/2 + BE1)
#   2 FOS/JUN/ATF3/HSPA1A AP-1 stress basal (통합 BE5 857/900)
#   3 KRT14 23/POSTN/VCAN/DKK2 basal (통합 BE3)
#   4 MMP7/WFDC2/PIGR/SCGB3A1/LTF/OLFM4, KRT5- club (통합 Club 532/561)
#   5 KRT13 35/ZNF750/PADI3/SPRR1B/ESR1/OVOL1 squamous-differentiated hillock (통합 Hillock1)
#   6 TOP2A(26%)/MKI67 cycling basal, decontX 0.26 (통합 BE1 286)
#   7 KLK3/KLK2/FOLH1/ACP3/MSMB AR+ luminal 61셀 (통합 LE(ARPC) 61/61)
#   ※ CRPC3 단독에는 OE(CFTR+)·Ionocyte·NE 클러스터 없음 (통합 Ionocyte 36셀은 산재)
# CRPC1 (res 0.4, 12 clusters; 2026-09-07 검토; 우리 QC 셀 7,852):
#   0 MMP7/WFDC2/SCGB3A1/OLFM4/FOLR1/HLA-DR club (통합 Club 1271 + LE 184)
#   1 CFTR/SLC4A4/PDE8B, nFeature 972 저count (통합 OE 758 + LE 286) — 통합 OE 의 실체
#   2 KRT13 29/LY6D 15/UPK1B/PSCA/S100P/H19 hillock-urothelial (통합 Hillock1 963)
#   3 TP63+ 저count(nFeature 790) basal (통합 BE6 804) — 통합 BE6 = 이 상태
#   4 KRT5/KRT15/COL17A1/CRYAB/SFRP1 basal, FOS/JUN 최고, decontX 0.35 (통합 BE1 504 + BE5 158)
#   5 DSG3/DSC3/KRT13 17/OLFM4/SNAI2 squamous-basal (통합 BE1 289 + Hillock2 349)
#   6 저count(1.2k, mt 9.6%) 핵 lncRNA, TP63 1.0 (통합 BE4 526) — CRPC2 'Basal low-count' 와 같은 상태
#   7 ALDH3A1/SHH/S100A2/DKK1/AKR1C3/MMP2/TGFBI basal (통합 BE2 498)
#   8 저count(1.0k) SHROOM3/PDE4C/TMPRSS2 luminal-like (통합 OE 275 + BE4 124) — CRPC2 cl5 와 같은 축
#   9 mt 12%, PPARG/SCHLAP1/lncRNA 저품질 171셀 (통합 BE4)
#   10 GNG4/CALML3/GABRP/CLDN8/OLFM4 11/KRT13 club-hillock 중간, decontX 0.41 (통합 BE5 120 + Club 24)
#   11 KLK3/KLK2/ACP3/MSMB/FOLH1 AR+ luminal 95셀 (통합 LE(ARPC) 95/95)
#   ※ CRPC1 단독에 Ionocyte(통합 43셀 산재)·NE 클러스터 없음
# 라벨 규칙 (2026-09-07 개정): "계통 + 잘 알려진 마커 1–2개 (+ 상태 수식어)". lncRNA·비유명 유전자는 라벨에 넣지 않음
# (근거 유전자 전체는 위 주석 / all_markers.csv / DotPlot 참조).
annotation_maps <- list(
    CRPC1 = c(
        `0`  = "Club MMP7/SCGB3A1",
        `1`  = "OE CFTR (low-count)",
        `2`  = "Hillock KRT13/PSCA",
        `3`  = "Basal TP63 (low-count)",
        `4`  = "Basal KRT5/KRT15 (AP-1 high)",
        `5`  = "Basal-squamous DSG3/KRT13",
        `6`  = "Basal (low-count)",
        `7`  = "Basal DKK1/SHH",
        `8`  = "Luminal TMPRSS2 (low-count)",
        `9`  = "Low-quality PPARG (mt high)",
        `10` = "Club-Hillock OLFM4/MMP7",
        `11` = "Luminal AR+ (KLK3)"
    ),
    CRPC2 = c(
        `0`  = "Basal DKK1/CCL2",
        `1`  = "Basal (low-count)",
        `2`  = "Basal stress (FOS/JUN)",
        `3`  = "Hillock KRT13/OLFM4",
        `4`  = "Hillock squamous (KRT1/SPRR1B)",
        `5`  = "Luminal TMPRSS2/PPARG (low-count)",
        `6`  = "Club SCGB1A1/TFF1",
        `7`  = "Ionocyte FOXI1",
        `8`  = "Squamous DSG3/KRT6A",
        `9`  = "Luminal AR+ (KLK3)",
        `10` = "Ionocyte-like KIT+"
    ),
    CRPC3 = c(
        `0` = "Basal DKK1/TNC",
        `1` = "Hillock KRT13/KRT6A",
        `2` = "Basal stress (FOS/JUN)",
        `3` = "Basal KRT14/POSTN",
        `4` = "Club MMP7/PIGR",
        `5` = "Hillock squamous (KRT13/ESR1)",
        `6` = "Basal cycling (TOP2A)",
        `7` = "Luminal AR+ (KLK3)"
    )
)
annotation_map <- annotation_maps[[PATIENT]]

if (length(annotation_map) > 0) {
    stopifnot(all(levels(epi$seurat_clusters) %in% names(annotation_map)))
    ann <- unname(annotation_map[as.character(epi$seurat_clusters)])
    epi$annotation <- factor(ann, levels = unique(unname(annotation_map)))
    Idents(epi) <- "annotation"
    n_ann <- nlevels(epi$annotation)

    p <- DimPlot(epi, group.by = "annotation", label = TRUE, repel = TRUE, pt.size = 0.3,
                 cols = utils_cb_palette(n_ann)) +
        ggtitle(paste0(PATIENT, " alone — epithelial annotation"))
    ggsave(file.path(OUT_DIR, "UMAP_annotation.png"), plot = p, width = 11, height = 8, bg = "white")
    p1 <- DimPlot(epi, group.by = "annotation", label = TRUE, repel = TRUE, pt.size = 0.2,
                  cols = utils_cb_palette(n_ann)) + ggtitle("Patient-alone annotation")
    p2 <- DimPlot(epi, group.by = "integrated_epi_annotation", label = TRUE, repel = TRUE, pt.size = 0.2,
                  cols = utils_cb_palette(length(unique(na.omit(epi$integrated_epi_annotation))))) +
        ggtitle("Integrated annotation")
    ggsave(file.path(OUT_DIR, "UMAP_annotation_vs_integrated.png"), plot = p1 + p2,
           width = 22, height = 8, bg = "white")
    p <- DotPlot(epi, features = epi_markers, group.by = "annotation", cluster.idents = FALSE) +
        RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        theme(axis.text.x = element_text(size = 7))
    ggsave(file.path(OUT_DIR, "DotPlot_epi_markers_by_annotation.png"), plot = p, width = 30, height = 7, bg = "white")

    tab <- as.data.frame.matrix(table(patient_alone = epi$annotation,
                                      integrated = addNA(factor(epi$integrated_epi_annotation))))
    write.csv(tab, file.path(OUT_DIR, "annotation_vs_integrated_annotation.csv"))
    write.csv(as.data.frame(table(annotation = epi$annotation)),
              file.path(OUT_DIR, "annotation_counts.csv"), row.names = FALSE)
    write.csv(data.frame(cell = colnames(epi), cluster = as.character(epi$seurat_clusters),
                         annotation = as.character(epi$annotation),
                         integrated_epi_annotation = epi$integrated_epi_annotation),
              file.path(OUT_DIR, "annotation_per_cell.csv"), row.names = FALSE)
    Idents(epi) <- "annotation"
    utils_save_all_markers(epi, file.path(OUT_DIR, "all_markers_by_annotation.csv"))
    saveRDS(epi, OBJ_RDS)
    message("annotation 부여 완료 → ", OBJ_RDS)
    print(table(epi$annotation))
    print(tab)
} else {
    message("annotation_map[", PATIENT, "] 비어 있음 — 진단 출력 검토 후 채우고 재실행 (reload 모드)")
}
