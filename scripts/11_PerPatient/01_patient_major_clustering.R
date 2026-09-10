#!/usr/bin/env Rscript
# Stage 11-01 — Per-patient (no-integration) major-lineage clustering + annotation
#
# 배경: 기존 파이프라인은 세 환자를 Harmony 로 통합한 뒤 클러스터링/annotation 했다.
# 세 원발암이 서로 크게 달라 통합 과정에서 환자 고유 구조가 묻혔을 가능성이 있어
# 환자별로 **단독** 클러스터링·annotation 을 수행한다 (CRPC1 은 외부 재주석
# reannotated_CRPC1.rds 가 이미 있으므로 CRPC2 / CRPC3 대상).
#
# 입력 셀 집합: Results/01_Integrated/combined_CRPC.rds 에서 orig.ident == PATIENT 만
#   subset. stage-01 의 decontX·scDblFinder·QC(nFeature 200–9000, mt<20) 는 전부
#   per-sample 로 수행됐으므로 이 subset 은 단독 재처리와 동일한 셀·count 집합이다.
#   통합 객체에서 가져오는 것은 raw(decontX) count 와 QC 메타데이터뿐이며, SCT /
#   PCA / Harmony / cluster / UMAP 은 모두 버리고 환자 단독으로 새로 계산한다.
#   통합 라벨(celltype)은 `integrated_celltype` 로만 보관해 비교용으로 쓴다.
#
# 파이프라인: SCTransform → PCA → FindNeighbors(dims 1:30) → FindClusters(다중 res)
#   → UMAP → CellCycle → SingleR(HPCA) → markers → [checkpoint]
#   → celltype_map (환자별, 마커/SingleR/통합라벨 교차표 검토 후 채움) → 최종 저장
#
# 실행 (프로젝트 루트): Rscript scripts/11_PerPatient/01_patient_major_clustering.R CRPC2
# 출력: Results/11_PerPatient/<PATIENT>/01_Major/

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
    library(SingleCellExperiment)
    library(BiocParallel)
    library(SingleR)
    library(celldex)
    library(future)
})
N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)
bpp <- MulticoreParam(workers = N_CORES, RNGseed = 42)

source("scripts/00_utils/scRNA_utils.R")

IN_RDS  <- "Results/01_Integrated/combined_CRPC.rds"
OUT_DIR <- file.path("Results/11_PerPatient", PATIENT, "01_Major")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
OBJ_RDS <- file.path(OUT_DIR, paste0(PATIENT, "_major.rds"))

RES_VEC   <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0)
RES_MAJOR <- "0.5"   # canonical 통합 stage-01 과 동일

computed <- !file.exists(OBJ_RDS)

if (computed) {
    # ------------------------------------------------------------------
    # Subset patient from canonical object → fresh Seurat object
    # ------------------------------------------------------------------
    combined <- readRDS(IN_RDS)
    keep <- colnames(combined)[combined$orig.ident == PATIENT]
    counts <- GetAssayData(combined[["RNA"]], layer = "counts")[, keep]
    counts <- counts[Matrix::rowSums(counts) > 0, ]   # 환자에서 0 인 유전자 제거
    md <- combined@meta.data[keep, ]
    md <- md[, intersect(c("orig.ident", "nCount_RNA", "nFeature_RNA", "percent.mt",
                           "decontX_contamination", "scDblFinder.score",
                           "S.Score", "G2M.Score", "Phase",
                           "SingleR", "SingleR_pruned",
                           "seurat_clusters", "celltype"), colnames(md))]
    md <- dplyr::rename(md,
        integrated_cluster  = seurat_clusters,
        integrated_celltype = celltype,
        integrated_SingleR  = SingleR,
        integrated_SingleR_pruned = SingleR_pruned,
        integrated_S.Score = S.Score, integrated_G2M.Score = G2M.Score,
        integrated_Phase = Phase)
    rm(combined); gc()

    obj <- CreateSeuratObject(counts = counts, meta.data = md, project = PATIENT)
    message(PATIENT, ": ", ncol(obj), " cells x ", nrow(obj), " genes")

    # ------------------------------------------------------------------
    # SCT → PCA → clustering → UMAP (환자 단독, integration 없음)
    # ------------------------------------------------------------------
    obj <- SCTransform(obj, verbose = FALSE)
    obj <- RunPCA(obj, verbose = FALSE)
    ggsave(file.path(OUT_DIR, "ElbowPlot.png"),
           plot = ElbowPlot(obj, ndims = 50), width = 10, height = 6, bg = "white")

    obj <- FindNeighbors(obj, reduction = "pca", dims = 1:30)
    obj <- FindClusters(obj, resolution = RES_VEC, random.seed = 42)
    obj$seurat_clusters <- obj@meta.data[[paste0("SCT_snn_res.", RES_MAJOR)]]
    obj$seurat_clusters <- factor(as.character(obj$seurat_clusters),
        levels = as.character(sort(unique(as.integer(as.character(obj$seurat_clusters))))))
    Idents(obj) <- "seurat_clusters"
    set.seed(42)
    obj <- RunUMAP(obj, reduction = "pca", dims = 1:30, seed.use = 42)

    obj <- CellCycleScoring(obj,
        s.features = cc.genes.updated.2019$s.genes,
        g2m.features = cc.genes.updated.2019$g2m.genes, nbin = 12)

    # ------------------------------------------------------------------
    # SingleR (HPCA main) — 환자 단독 SCT 데이터 기준으로 재실행
    # ------------------------------------------------------------------
    ref_hpca <- celldex::HumanPrimaryCellAtlasData()
    sce <- as.SingleCellExperiment(obj, assay = "SCT")
    pred <- SingleR(test = sce, ref = ref_hpca, labels = ref_hpca$label.main, BPPARAM = bpp)
    obj$SingleR <- pred$labels
    obj$SingleR_pruned <- pred$pruned.labels
    rm(sce, pred); gc()

    # ------------------------------------------------------------------
    # Markers at RES_MAJOR + checkpoint
    # ------------------------------------------------------------------
    utils_save_all_markers(obj, file.path(OUT_DIR, "all_markers.csv"))
    saveRDS(obj, OBJ_RDS)
} else {
    message("Loading cached ", OBJ_RDS, " — figures / celltype map only")
    obj <- readRDS(OBJ_RDS)
    Idents(obj) <- "seurat_clusters"
}

# ======================================================================
# Diagnostics (both modes)
# ======================================================================
n_cl <- nlevels(obj$seurat_clusters)

p_tree <- clustree(obj, prefix = "SCT_snn_res.")
ggsave(file.path(OUT_DIR, "clustree.png"), plot = p_tree, width = 12, height = 14, dpi = 200, bg = "white")

umap_list <- lapply(RES_VEC, function(r) {
    rcol <- paste0("SCT_snn_res.", r)
    DimPlot(obj, group.by = rcol, label = TRUE, pt.size = 0.2,
            cols = utils_cb_palette(dplyr::n_distinct(obj@meta.data[[rcol]]))) +
        ggtitle(paste0("res = ", r)) + NoLegend()
})
ggsave(file.path(OUT_DIR, "UMAP_multires.png"), plot = wrap_plots(umap_list, ncol = 4),
       width = 24, height = 12, dpi = 200, bg = "white")

p <- DimPlot(obj, group.by = "seurat_clusters", label = TRUE, pt.size = 0.3,
             cols = utils_cb_palette(n_cl)) +
    ggtitle(paste0(PATIENT, " alone — clusters (res ", RES_MAJOR, ")"))
ggsave(file.path(OUT_DIR, "UMAP_cluster.png"), plot = p, width = 10, height = 8, bg = "white")

p <- DimPlot(obj, group.by = "integrated_celltype", label = TRUE, repel = TRUE, pt.size = 0.3,
             cols = utils_cb_palette(dplyr::n_distinct(obj$integrated_celltype))) +
    ggtitle("Integrated (Harmony) celltype projected on patient-alone UMAP")
ggsave(file.path(OUT_DIR, "UMAP_integrated_celltype.png"), plot = p, width = 11, height = 8, bg = "white")

p <- DimPlot(obj, group.by = "SingleR_pruned", label = TRUE, repel = TRUE, pt.size = 0.3,
             cols = utils_cb_palette(length(unique(na.omit(obj$SingleR_pruned))))) +
    ggtitle("SingleR (HPCA main, patient alone)")
ggsave(file.path(OUT_DIR, "UMAP_SingleR.png"), plot = p, width = 12, height = 8, bg = "white")

p <- DimPlot(obj, group.by = "Phase", pt.size = 0.3,
             cols = utils_cb_palette(dplyr::n_distinct(obj$Phase)))
ggsave(file.path(OUT_DIR, "UMAP_CellCycle.png"), plot = p, width = 10, height = 7, bg = "white")

p <- FeaturePlot(obj, features = "decontX_contamination", pt.size = 0.2) +
    viridis::scale_color_viridis()
ggsave(file.path(OUT_DIR, "UMAP_decontX.png"), plot = p, width = 10, height = 8, bg = "white")

p <- VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt", "decontX_contamination"),
             group.by = "seurat_clusters", pt.size = 0, ncol = 2, cols = utils_cb_palette(n_cl))
ggsave(file.path(OUT_DIR, "QC_violin_by_cluster.png"), plot = p, width = 16, height = 10, bg = "white")

lineage_markers <- list(
    Epithelial  = c("EPCAM", "KRT8", "KRT18", "KRT5", "KRT15", "AR", "KLK3"),
    NE          = c("CHGA", "SYP", "ENO2"),
    Fibroblast  = c("DCN", "LUM", "COL1A1", "PDGFRA"),
    Pericyte    = c("RGS5", "NOTCH3", "HIGD1B"),
    SMC         = c("MYH11", "ACTG2", "DES", "CNN1"),
    Schwann     = c("PLP1", "S100B", "SOX10"),
    Endothelial = c("PECAM1", "VWF", "CLDN5"),
    Lymphatic   = c("PROX1", "LYVE1"),
    T           = c("CD3D", "CD3E", "CD8A", "CD4", "IL7R"),
    NK          = c("NKG7", "GNLY", "KLRD1"),
    B           = c("MS4A1", "CD79A"),
    Plasma      = c("JCHAIN", "MZB1", "IGKC"),
    Myeloid     = c("LYZ", "CD68", "C1QA", "CD14", "FCGR3A"),
    DC          = c("CLEC9A", "CD1C", "LAMP3"),
    Mast        = c("TPSAB1", "CPA3", "KIT"),
    Cycling     = c("MKI67", "TOP2A")
)
lineage_markers <- lapply(lineage_markers, function(g) intersect(g, rownames(obj)))
p <- DotPlot(obj, features = lineage_markers, group.by = "seurat_clusters", cluster.idents = FALSE) +
    RotatedAxis() +
    scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
    theme(axis.text.x = element_text(size = 8)) +
    labs(title = paste0(PATIENT, " alone — lineage markers by cluster"))
ggsave(file.path(OUT_DIR, "DotPlot_lineage_by_cluster.png"), plot = p, width = 26, height = 8, bg = "white")

fp_genes <- c("EPCAM", "KRT8", "KRT5", "AR", "KLK3", "PTPRC", "CD3E", "LYZ",
              "PECAM1", "DCN", "RGS5", "MYH11", "TPSAB1", "JCHAIN", "SYP", "MKI67")
p <- FeaturePlot(obj, features = fp_genes, ncol = 4, pt.size = 0.1)
ggsave(file.path(OUT_DIR, "FeaturePlot_lineage.png"), plot = p, width = 20, height = 20, bg = "white")

# 교차표: 환자 단독 cluster × 통합 celltype / SingleR
write.csv(as.data.frame.matrix(table(cluster = obj$seurat_clusters, integrated = obj$integrated_celltype)),
          file.path(OUT_DIR, "cluster_vs_integrated_celltype.csv"))
write.csv(as.data.frame.matrix(table(cluster = obj$seurat_clusters, SingleR = obj$SingleR)),
          file.path(OUT_DIR, "cluster_vs_SingleR.csv"))
# 마커 패널 cluster 평균 (annotation 근거용)
avg <- AverageExpression(obj, features = unlist(lineage_markers), group.by = "seurat_clusters",
                         assay = "SCT", layer = "data")$SCT
write.csv(round(as.matrix(avg), 3), file.path(OUT_DIR, "lineage_marker_avg_by_cluster.csv"))
qc_tab <- obj@meta.data %>% group_by(cluster = seurat_clusters) %>%
    summarise(n = n(), median_nCount = median(nCount_RNA), median_nFeature = median(nFeature_RNA),
              median_mt = round(median(percent.mt), 2),
              median_decontX = round(median(decontX_contamination), 3),
              frac_cycling = round(mean(Phase != "G1"), 3), .groups = "drop")
write.csv(qc_tab, file.path(OUT_DIR, "cluster_QC_summary.csv"), row.names = FALSE)

# ======================================================================
# Celltype assignment — 환자별 맵 (all_markers / DotPlot / 교차표 검토 후 채움)
# ======================================================================
# ⚠️ cluster 정수 ID 는 rerun 간 불안정할 수 있음. ID 가 맵에 없으면 stopifnot 으로 멈춤.
# CRPC3 (res 0.5, 16 clusters; 2026-09-07 검토):
#   Epithelial: 0,1,3,5,9 (KRT5/KRT15 basal-like; cl1=KRT6A/KRT13 hillock-like, cl3=DKK1/CCL2/TNC),
#               2 (KRT15 최고·lncRNA), 6 (KRT8/18+ KRT5- luminal; MMP7/PIGR/LTF club-like),
#               15 (KLK3/KLK2/FOLH1/MSMB = AR+ luminal, 62셀)
#   Fibroblast: 4 (ADH1B/CFD/WNT2), 7 (SFRP4/ASPN/THBS2)
#   Pericyte: 13 (RGS5/HIGD1B/NDUFA4L2 + MYH11 → mural; 통합 Pericyte 121/SMC 9)
#   Endothelial: 10 (VWF/SELE/ACKR1)
#   T/NK cells: 8 (CD3D/CD2/NKG7; JCHAIN+ plasma·TPSAB1+ mast 소수 포함, res1.0 에서도 분리 안 됨)
#   Phagocytes: 12 (IL1B/MMP9/EREG 염증성 단핵구), 14 (F13A1/FOLR2/MRC1/LYVE1 조직 대식세포)
#   Low-quality: 11 (nFeature 중앙값 1322, MALAT1/NEAT1 핵 RNA, 통합라벨 Epi179/Fibro30/Peri8/SMC6 혼재)
# CRPC2 (res 0.5, 19 clusters; 2026-09-07 검토):
#   Epithelial: 0 (KRT14/COL17A1/DKK1/TNC basal), 2 (저count lncRNA, SingleR Epi 99%), 3 (SPRR1B/KRT1/KRT13 squamous),
#               4 (SOX2/KRT15/HOXB13/FGFR2), 8 (저count PPARG/TMPRSS2/TMC5 luminal-like), 9 (TFF1/SCGB3A1/PIGR club),
#               13 (FOXI1/ATP6V ionocyte 153셀), 15 (DSG3/SERPINB13 squamous, mt 16%), 17 (AR 최고·XACT/POTE, 47셀),
#               18 (TFCP2L1/ATP6V0D2/KIT ionocyte-like 40셀)
#   Fibroblast: 1 (LUM/DCN/SFRP2), 7 (저count GLI2/GLIS1/EYA1; 통합 Fibro 391/SMC 36)
#   Pericyte: 10 (RGS5/PLN/RERGL/NDUFA4L2 + MYH11 = 혈관 mural)
#   Smooth muscle cells: 11 (MYH11/CARMN/MYOCD, RGS5 낮음; 통합은 Pericyte 로 묶었던 셀)
#   Schwann cell: 16 (PLP1/S100B/CDH19),  Endothelial: 5,  T/NK cells: 6,  Mast cells: 12,  Phagocytes: 14
# CRPC1 (res 0.5, 21 clusters; 2026-09-07 검토):
#   Epithelial: 0 (SCGB2A1/MMP7/HOXB13 club-luminal), 1 (CFTR/SLC4A4, nFeature 976 저count = 통합 OE 축),
#               2 (UPK1B/PSCA/LY6D/S100P), 5 (TP63+ 저count 769), 7 (저count 1.3k, mt 10%), 8 (KRT5/KRT15/COL17A1 basal),
#               9 (DSG3/DSC3/SERPINB5 squamous-basal), 10 (ALDH3A1/SHH/DKK1/S100A2 basal), 11 (저count 1.1k, SHROOM3)
#   Fibroblast: 3 (SFRP2/DPT/APOD/PDGFRA)
#   Pericyte: 13 (RGS5/HGF/PDE3A), 14 (RGS5/MYH11/PLN/RERGL 혈관 mural)
#   Smooth muscle cells: 17 (CARMN/ACTG2/MYH11, 저count)
#   Endothelial: 6,  T/NK cells: 4 (CD8 T), 12 (저count THEMIS/SKAP1 T), 15 (NK: GNLY/KLRF1/FGFBP2)
#   Phagocytes: 16 (MMP9/IL1B/CD163), 19 (APOC1/C1Q/APOE 저count),  Mast cells: 18,  Erythrocytes: 20 (HBB/HBA, 25셀)
celltype_maps <- list(
    CRPC1 = c(
        `0` = "Epithelial", `1` = "Epithelial", `2` = "Epithelial", `3` = "Fibroblast",
        `4` = "T/NK cells", `5` = "Epithelial", `6` = "Endothelial", `7` = "Epithelial",
        `8` = "Epithelial", `9` = "Epithelial", `10` = "Epithelial", `11` = "Epithelial",
        `12` = "T/NK cells", `13` = "Pericyte", `14` = "Pericyte", `15` = "T/NK cells",
        `16` = "Phagocytes", `17` = "Smooth muscle cells", `18` = "Mast cells",
        `19` = "Phagocytes", `20` = "Erythrocytes"
    ),
    CRPC2 = c(
        `0` = "Epithelial", `1` = "Fibroblast", `2` = "Epithelial", `3` = "Epithelial",
        `4` = "Epithelial", `5` = "Endothelial", `6` = "T/NK cells", `7` = "Fibroblast",
        `8` = "Epithelial", `9` = "Epithelial", `10` = "Pericyte", `11` = "Smooth muscle cells",
        `12` = "Mast cells", `13` = "Epithelial", `14` = "Phagocytes", `15` = "Epithelial",
        `16` = "Schwann cell", `17` = "Epithelial", `18` = "Epithelial"
    ),
    CRPC3 = c(
        `0` = "Epithelial", `1` = "Epithelial", `2` = "Epithelial", `3` = "Epithelial",
        `4` = "Fibroblast", `5` = "Epithelial", `6` = "Epithelial", `7` = "Fibroblast",
        `8` = "T/NK cells", `9` = "Epithelial", `10` = "Endothelial", `11` = "Low-quality",
        `12` = "Phagocytes", `13` = "Pericyte", `14` = "Phagocytes", `15` = "Epithelial"
    )
)
celltype_map <- celltype_maps[[PATIENT]]
celltype_levels <- c("Epithelial", "Neuroendocrine", "Fibroblast", "Pericyte", "Smooth muscle cells",
                     "Schwann cell", "Endothelial", "Lymphatic endothelial",
                     "T/NK cells", "B cells", "Plasma cells",
                     "Phagocytes", "Dendritic cells", "Mast cells", "Erythrocytes", "Cycling", "Low-quality")

if (length(celltype_map) > 0) {
    stopifnot(all(levels(obj$seurat_clusters) %in% names(celltype_map)))
    ct <- unname(celltype_map[as.character(obj$seurat_clusters)])
    stopifnot(all(ct %in% celltype_levels))
    obj$celltype <- factor(ct, levels = intersect(celltype_levels, unique(ct)))
    Idents(obj) <- "celltype"
    n_ct <- nlevels(obj$celltype)

    p <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.3,
                 cols = utils_cb_palette(n_ct)) +
        ggtitle(paste0(PATIENT, " alone — cell types"))
    ggsave(file.path(OUT_DIR, "UMAP_celltype.png"), plot = p, width = 11, height = 8, bg = "white")

    p1 <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.2,
                  cols = utils_cb_palette(n_ct)) + ggtitle("Patient-alone annotation")
    p2 <- DimPlot(obj, group.by = "integrated_celltype", label = TRUE, repel = TRUE, pt.size = 0.2,
                  cols = utils_cb_palette(dplyr::n_distinct(obj$integrated_celltype))) +
        ggtitle("Integrated annotation")
    ggsave(file.path(OUT_DIR, "UMAP_celltype_vs_integrated.png"), plot = p1 + p2,
           width = 22, height = 8, bg = "white")

    p <- DotPlot(obj, features = lineage_markers, group.by = "celltype", cluster.idents = FALSE) +
        RotatedAxis() +
        scale_color_gradient2(low = "#1F77B4", mid = "grey90", high = "#D7261E") +
        theme(axis.text.x = element_text(size = 8))
    ggsave(file.path(OUT_DIR, "DotPlot_lineage_by_celltype.png"), plot = p, width = 26, height = 7, bg = "white")

    tab <- as.data.frame.matrix(table(patient_alone = obj$celltype, integrated = obj$integrated_celltype))
    write.csv(tab, file.path(OUT_DIR, "celltype_vs_integrated_celltype.csv"))
    write.csv(as.data.frame(table(celltype = obj$celltype)),
              file.path(OUT_DIR, "celltype_counts.csv"), row.names = FALSE)
    write.csv(data.frame(cell = colnames(obj), cluster = as.character(obj$seurat_clusters),
                         celltype = as.character(obj$celltype),
                         integrated_celltype = as.character(obj$integrated_celltype)),
              file.path(OUT_DIR, "celltype_per_cell.csv"), row.names = FALSE)

    saveRDS(obj, OBJ_RDS)
    message("celltype 부여 완료 → ", OBJ_RDS)
    print(table(obj$celltype))
    print(tab)
} else {
    message("celltype_map[", PATIENT, "] 비어 있음 — 진단 출력 검토 후 채우고 재실행 (reload 모드)")
}
