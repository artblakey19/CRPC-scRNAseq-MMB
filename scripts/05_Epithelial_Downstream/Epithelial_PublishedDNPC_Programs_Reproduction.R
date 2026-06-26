# Epithelial_PublishedDNPC_Programs_Reproduction.R
# ------------------------------------------------------------------------------
# 목적: "여러 DNPC 연구의 transcriptional program이 우리 DNPC late-stage cohort의
#        epithelial 단일세포에서 재현되는가?"를 한 번의 객체 로드로 통합 검정.
#
# 이미 별도 스크립트로 검정한 연구(중복 방지):
#   - Bluemn 2017 (FGF/MAPK)          -> Epithelial_Bluemn2017_ssGSEA.R / _pseudobulk_GSEA.R
#   - Kim 2024  (AR/HGF/Wnt, XPO1)    -> Epithelial_Kim2024_AR_HGF_Wnt_*.R
#   - Song 2022 (BE/LE/Hillock/Club)  -> Epithelial_Song2022_*.R
#
# 새로 검정하는 미검정 DNPC program (10-study 문헌 sweep + 합성 카탈로그 기반):
#
#   [A] He et al. 2024, npj Prec Oncol (HuPSA): KRT7 subtype / SOX2·FOXA2 progenitor
#   [B] Tang et al. 2022, Science: CRPC-AR / -NE / -WNT / -SCL(stem-cell-like)
#       * CRPC-SCL은 AP-1(FOS/JUN) + YAP/TAZ-target + stem/basal로 교정.
#         (TWIST2/SOX4/STAT1/2는 Tang이 아니라 Beltran 2024 PNAS임 — 초기 버전 오귀속 수정)
#   [C] Han et al. 2025, Cancer Cell: KMT2C->ΔNp63(basal)->SREBP1c/FASN lipogenesis
#   [D] Labrecque et al. 2019, JCI: 5-phenotype mCRPC 분류의 DNPC 정의.
#       * DNPC-특이 POSITIVE 마커(CXCL8/CXCR1/RUNX2/TGFB) — AR-/NE- "gate"가 아닌
#         양성 readout. CXCL8 축은 cohort의 BE2 CXCL1/CCL2 finding과 연결.
#       * SQUAM panel (KRT5/6A/6B/14/15/17/FGFBP1/DSG3/IVL/SBSN/S100A7).
#       * NEURO master-TF panel = NE-null 음성대조.
#   [E] POU2F3 tuft-cell chemosensory program (Oncogene 2023; SCN 4-TF framework):
#       AR-/classical-NE- lineage. cohort의 FOXI1+ Ionocyte(cl10)의 sibling.
#       기존 어떤 시그니처로도 안 잡힌 최대 gap.
#   [F] PNAS 2025 / CCS 2026 review: KLF5/RUNX1-coordinated mixed basal/club/hillock 정체성.
#   [G] Bluemn FGFR/MAPK output (gene-level): DUSP6/ETV4/SPRY2/4/ID1-4 — MSigDB proxy보다
#       sharp한 cross-check.
#   [H] Thienger 2024/26 (Cancer Research): SMARCA2/4-degrader-vulnerable TF panel +
#       CRC_WNT target-output score (WNT-component 축과 분리 — [[project_wnt_signature_divergence]]).
#   [I] Quigley 2018 double-negative MYC-target program: HALLMARK_MYC_TARGETS_V1 (msigdbr).
#
#   [QC] AP-1/IEG(FOS/JUN/EGR1/DUSP1/NR4A...)는 droplet scRNA의 대표적 dissociation
#        artifact. 이를 검증하기 위해 heat-shock dissociation 대조 program을 함께 스코어링하고
#        AP-1 ↔ heat-shock ↔ QC(percent.mt/nCount/nFeature)의 상관 + per-sample 분포를 산출.
#        AP-1이 heat-shock/QC와 강하게 동행하면 biology가 아닌 prep artifact로 해석해야 함.
#
# 방법(기존 module-score 스크립트 idiom 준수):
#   - 입력: Results/04_Epithelial_Filtered/epi_filtered_clustered.rds (+ annotation CSV 재부착)
#   - Seurat::AddModuleScore (assay=RNA, default ctrl=100), 증분 캐시(미계산 program만 계산)
#   - 연속형 score UMAP은 viridis ([[feedback_featureplot_color]])
#   - discrete fill은 색약 친화 utils_cb_palette ([[feedback_colorblind_palettes]])
#   - per-cluster mean / one-way ANOVA / marker DotPlot / 통합 verdict heatmap
#
# Caveat (memory 준수):
#   - benign vs malignant 미확정([[project_epi_benign_malignant_unknown]]) — "어느 cluster가
#     program을 발현하는가"만 기술. malignancy/방향성 trajectory 단정 금지.
#   - CNV/genomic-driver argument 재도입 금지([[feedback_no_cnv_arguments]]).
# ------------------------------------------------------------------------------

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(ggplot2)
    library(ggpubr)
    library(patchwork)
    library(msigdbr)
})
suppressMessages(source("scripts/00_utils/scRNA_utils.R"))

IN_RDS       <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR      <- "Results/05_Epithelial_Downstream/PublishedDNPC_Programs"
OUT_RDS      <- "Results/05_Epithelial_Downstream/epi_publishedDNPC_programs.rds"
ANNOT_CSV    <- "Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv"
ANNOT_LEVELS <- "Results/05_Epithelial_Downstream/Annotation/label_levels.txt"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

# ==============================================================================
# Multi-gene PROGRAMS (AddModuleScore). gene 별칭은 둘 다 넣고 intersect로 거름.
# ==============================================================================
PROGRAMS <- list(
    # --- [A] He et al. 2024 (HuPSA) DNPC subtypes ---
    He2024_KRT7_DNPC        = c("KRT7", "S100A9", "S100A8", "S100A4", "S100A6",
                                "S100A10", "S100A11", "S100A14", "S100A16", "S100P"),
    He2024_ProgenitorLike   = c("SOX2", "FOXA2", "SRRM4", "ONECUT2", "POU3F2", "PEG10"),

    # --- [B] Tang et al. 2022 (Science) CRPC subtypes ---
    Tang2022_CRPC_AR        = c("AR", "KLK3", "KLK2", "NKX3-1", "FOLH1", "TMPRSS2", "STEAP2"),
    Tang2022_CRPC_NE        = c("CHGA", "CHGB", "SYP", "NCAM1", "INSM1", "ASCL1", "SOX2"),
    Tang2022_CRPC_WNT       = c("TCF7L2", "LEF1", "AXIN2", "ASCL2", "NKD1", "RNF43",
                                "LGR5", "SOX9"),
    # CRPC-SCL 교정판: AP-1 + YAP/TAZ-target + stem/basal (CCN1=CYR61, CCN2=CTGF 별칭 포함)
    Tang2022_CRPC_SCL       = c("FOSL1", "FOS", "FOSB", "JUNB", "JDP2", "BATF", "MAFF",
                                "ATF3", "JUN", "JUND", "FOSL2", "CD44", "TACSTD2",
                                "CCN1", "CYR61", "CCN2", "CTGF", "ANKRD1", "AJUBA",
                                "AXL", "CCND1", "KRT5", "KRT14", "KRT17", "TP63"),

    # --- [C] Han et al. 2025 (Cancer Cell) KMT2C/ΔNp63/lipogenesis ---
    Han2025_dNp63_basal     = c("TP63", "KRT5", "KRT14", "KRT6A", "KRT15"),
    Han2025_SREBP_lipogen   = c("SREBF1", "FASN", "SCD", "ACACA", "ACLY", "ELOVL6",
                                "HMGCR", "INSIG1", "LDLR", "MVD", "SQLE"),

    # --- [D] Labrecque et al. 2019 (JCI) DNPC definition ---
    Labrecque2019_DNPC_pos  = c("CXCL8", "IL8", "CXCR1", "RUNX2", "TGFB1", "TGFB2"),
    Labrecque2019_SQUAM     = c("KRT5", "KRT6A", "KRT6B", "FGFBP1", "DSG3", "KRT14",
                                "KRT15", "KRT17", "TP63", "IVL", "SBSN", "S100A7"),
    Labrecque2019_NEURO_ctl = c("ASCL1", "NEUROD1", "POU3F2", "FOXA2", "PDX1", "INSM1",
                                "SOX2", "CHGA", "SYP", "NCAM1", "ENO2", "MYCN", "DLX5",
                                "SRRM4", "SNAP25", "VGF", "SCG3", "CHGB", "CHRNB2",
                                "NKX2-1", "LMO3"),

    # --- [E] POU2F3 tuft-cell program (최대 gap; Ionocyte sibling) ---
    POU2F3_tuft             = c("POU2F3", "GFI1B", "ASCL2", "TRPM5", "AVIL", "PLCB2",
                                "GNG13", "LRMP", "BMX", "HCK", "PLCG2", "RGS13",
                                "SOX9", "PTGS1"),

    # --- [F] PNAS 2025 KLF5/RUNX1 mixed basal/club/hillock identity ---
    PNAS2025_KLF5_mixed     = c("KLF5", "RUNX1", "RARG", "LTF", "MMP7", "PIGR",
                                "SCGB1A1", "SCGB3A1", "KRT13", "KRT4", "KRT5",
                                "KRT14", "KRT15", "CP", "OLFM4"),

    # --- [G] Bluemn FGFR/MAPK output, gene-level cross-check ---
    Bluemn2017_FGFR_output  = c("DUSP6", "ETV4", "EGR1", "SPRY2", "SPRY4", "ETV5",
                                "PHLDA1", "SPRED1", "DUSP4", "ID1", "ID2", "ID3", "ID4"),

    # --- [H] Thienger WNT-target output + SWI/SNF-degrader TF panel ---
    Thienger_WNT_target     = c("LGR5", "AXIN2", "ASCL2", "OLFM4", "SLC12A2", "NKD1",
                                "WIF1", "RNF43", "TCF7L2", "LEF1"),
    Thienger_SWISNF_TF      = c("TCF7L2", "TCF7L1", "LEF1", "KLF2", "KLF5", "SOX4",
                                "SOX13", "RUNX3", "SIX2", "TCF7"),

    # --- [QC] AP-1/IEG stress-persister + heat-shock dissociation 대조 ---
    AP1_IEG_stress          = c("FOS", "FOSB", "FOSL1", "FOSL2", "JUN", "JUNB", "JUND",
                                "ATF3", "EGR1", "DUSP1", "NR4A1", "NR4A2", "NR4A3",
                                "ZFP36", "IER2", "KLF6", "MAFF"),
    HeatShock_dissoc_ctl    = c("HSPA1A", "HSPA1B", "HSPB1", "DNAJB1", "HSPH1",
                                "HSP90AA1", "BAG3", "HSPA6", "DNAJA1")
)

# --- [I] Quigley double-negative MYC program = HALLMARK_MYC_TARGETS_V1 (msigdbr) ---
myc <- tryCatch({
    h <- msigdbr(species = "Homo sapiens", collection = "H")
    unique(h$gene_symbol[h$gs_name == "HALLMARK_MYC_TARGETS_V1"])
}, error = function(e) { warning("msigdbr MYC fetch failed: ", conditionMessage(e)); character(0) })
if (length(myc) >= 10) PROGRAMS[["Quigley2018_MYC_targets"]] <- myc

# ==============================================================================
# MARKER PANELS (DotPlot) — subtype/lineage 판별 1차 근거
# ==============================================================================
MARKER_PANELS <- list(
    `AR/Luminal`     = c("AR", "KLK3", "NKX3-1", "HOXB13", "FOXA1"),
    Basal            = c("KRT5", "KRT14", "TP63"),
    Squamous         = c("KRT6A", "DSG3", "FGFBP1", "IVL", "S100A7"),
    Club             = c("MMP7", "WFDC2", "SCGB1A1", "LTF", "PIGR"),
    NE               = c("CHGA", "SYP", "NCAM1", "INSM1", "ASCL1", "NEUROD1"),
    `He:KRT7-DNPC`   = c("KRT7", "S100A9", "S100A8", "S100A4"),
    `He:Progenitor`  = c("SOX2", "FOXA2", "ONECUT2", "POU3F2", "PEG10"),
    `Lab:DNPC+`      = c("CXCL8", "CXCR1", "RUNX2", "TGFB1"),
    Tuft             = c("POU2F3", "GFI1B", "TRPM5", "AVIL", "LRMP"),
    Ionocyte         = c("FOXI1", "ASCL3", "CFTR", "ATP6V0B"),
    `Tang:SCL`       = c("FOSL1", "JUNB", "CD44", "TACSTD2", "AXL", "CCND1"),
    `WNT-target`     = c("TCF7L2", "LGR5", "AXIN2", "OLFM4", "SOX9"),
    `KLF5-mixed`     = c("KLF5", "RUNX1", "KRT13", "KRT4"),
    `FGFR-output`    = c("DUSP6", "ETV4", "ID1", "ID3"),
    `AP1/IEG`        = c("FOS", "JUN", "EGR1", "DUSP1"),
    `HeatShock(QC)`  = c("HSPA1A", "HSPA1B", "DNAJB1")
)

# ==============================================================================
# Load-or-compute (증분): OUT_RDS 있으면 로드 후 미계산 program만 추가 스코어링
# ==============================================================================
if (file.exists(OUT_RDS)) {
    message("Loading cached ", OUT_RDS)
    epi <- readRDS(OUT_RDS)
    DefaultAssay(epi) <- "RNA"
    if (!"data" %in% Layers(epi[["RNA"]]))
        epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)
} else {
    epi <- readRDS(IN_RDS)
    message("Loaded: ", IN_RDS, " | cells=", ncol(epi), " genes=", nrow(epi))
    DefaultAssay(epi) <- "RNA"
    if (!"data" %in% Layers(epi[["RNA"]]))
        epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)
}
present <- rownames(epi)

# 현재 PROGRAMS에 없는 stale score 컬럼 제거 (예: 초기 버그 Tang2022_CRPC_SCL_AP1)
current_cols <- paste0(names(PROGRAMS), "_score")
stale <- setdiff(grep("_score$", colnames(epi@meta.data), value = TRUE), current_cols)
for (s in stale) { message("Dropping stale score column: ", s); epi@meta.data[[s]] <- NULL }

# program별 gene presence 보고(전체) + 미계산분만 AddModuleScore
presence_rows <- list()
n_new <- 0
for (nm in names(PROGRAMS)) {
    g0 <- PROGRAMS[[nm]]
    g1 <- intersect(g0, present)
    presence_rows[[nm]] <- data.frame(program = nm, n_total = length(g0),
        n_present = length(g1), absent = paste(setdiff(g0, present), collapse = ", "))
    col <- paste0(nm, "_score")
    if (col %in% colnames(epi@meta.data)) next            # 이미 스코어링됨 -> skip
    if (length(g1) < 2) { warning("Program dropped (<2 genes): ", nm); next }
    message(sprintf("  scoring %-26s %d/%d genes", nm, length(g1), length(g0)))
    epi <- tryCatch(
        AddModuleScore(epi, features = list(g1), name = paste0(nm, "_score"), assay = "RNA"),
        error = function(e)
            AddModuleScore(epi, features = list(g1), name = paste0(nm, "_score"),
                           nbin = 12, assay = "RNA"))
    src <- paste0(nm, "_score1")
    epi@meta.data[[col]] <- epi@meta.data[[src]]; epi@meta.data[[src]] <- NULL
    n_new <- n_new + 1
}
write.csv(do.call(rbind, presence_rows),
          file.path(OUT_DIR, "program_gene_presence.csv"), row.names = FALSE)
if (n_new > 0) { saveRDS(epi, OUT_RDS); message("Saved (", n_new, " new programs): ", OUT_RDS) }

score_cols <- grep("_score$", colnames(epi@meta.data), value = TRUE)
prog_label <- function(col) sub("_score$", "", col)

# ==============================================================================
# annotation 재부착 + 정렬된 factor (Epithelial_Bluemn2017_ssGSEA.R 패턴)
# ==============================================================================
if (file.exists(ANNOT_CSV) && file.exists(ANNOT_LEVELS)) {
    .a <- read.csv(ANNOT_CSV, stringsAsFactors = FALSE)
    epi$annotation <- factor(.a$annotation[match(colnames(epi), .a$cell)],
                             levels = readLines(ANNOT_LEVELS))
    clusters <- epi$annotation
} else {
    warning("annotation files not found — falling back to seurat_clusters")
    clusters <- factor(epi$seurat_clusters)
}
Idents(epi) <- clusters
pal <- utils_cb_palette(nlevels(clusters))

# ==============================================================================
# (1) MASTER DotPlot — subtype/lineage 판별 1차 근거
# ==============================================================================
mpp <- lapply(MARKER_PANELS, function(g) intersect(g, rownames(epi)))
.seen <- character(0)                          # panel 간 중복 gene 제거(중복 factor level 방지)
for (nm in names(mpp)) { g <- setdiff(mpp[[nm]], .seen); mpp[[nm]] <- g; .seen <- c(.seen, g) }
mpp <- mpp[vapply(mpp, length, 1L) > 0]
dp <- DotPlot(epi, features = mpp, cluster.idents = FALSE) +
    scale_color_viridis_c(option = "viridis") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          strip.text.x = element_text(angle = 90, size = 7)) +
    ggtitle("Published DNPC markers across epithelial clusters")
ggsave(file.path(OUT_DIR, "DotPlot_DNPC_markers_by_cluster.png"),
       plot = dp, width = 50, height = 16, units = "cm", dpi = 220, bg = "white")

# ==============================================================================
# (2) Per-program score UMAP (viridis) + per-cluster violin (cb palette)
# ==============================================================================
for (sc in score_cols) {
    title <- prog_label(sc); rng <- range(epi@meta.data[[sc]], na.rm = TRUE)
    p_umap <- FeaturePlot(epi, features = sc, pt.size = 0.3, order = TRUE) +
        scale_color_viridis_c(option = "viridis", limits = rng) + ggtitle(title)
    ggsave(file.path(OUT_DIR, paste0("UMAP_", title, ".png")),
           plot = p_umap, width = 20, height = 18, units = "cm", dpi = 200, bg = "white")

    df <- data.frame(cluster = clusters, score = epi@meta.data[[sc]])
    ymax <- max(df$score, na.rm = TRUE)
    p_v <- ggplot(df, aes(cluster, score, fill = cluster)) +
        geom_violin(scale = "width") + scale_fill_manual(values = pal) +
        stat_summary(fun = mean, geom = "point", size = 6, colour = "blue", shape = 95) +
        stat_compare_means(method = "anova", label.x = 3, label.y = ymax + 0.02) +
        ggtitle(title) + NoLegend() +
        theme(text = element_text(size = 11), axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, paste0(title, "_violin.png")),
           plot = p_v, width = 20, height = 15, units = "cm", dpi = 200, bg = "white")
}

# 핵심 DNPC program: split-by-sample UMAP (CRPC1/2/3 분포; [[project_pi_hypothesis_be12_is_dnpc]])
n_samples <- length(unique(epi$orig.ident))
split_targets <- intersect(c("He2024_KRT7_DNPC_score", "Tang2022_CRPC_SCL_score",
    "Labrecque2019_DNPC_pos_score", "POU2F3_tuft_score", "AP1_IEG_stress_score"), score_cols)
for (sc in split_targets) {
    title <- prog_label(sc); rng <- range(epi@meta.data[[sc]], na.rm = TRUE)
    p <- FeaturePlot(epi, features = sc, split.by = "orig.ident", pt.size = 0.3, order = TRUE) &
        scale_color_viridis_c(option = "viridis", limits = rng)
    ggsave(file.path(OUT_DIR, paste0("UMAP_", title, "_splitBySample.png")),
           plot = p + plot_annotation(title = title),
           width = 8 * n_samples, height = 8, dpi = 200, bg = "white")
}

# ==============================================================================
# (3) One-way ANOVA + per-cluster mean
# ==============================================================================
anova_rows <- lapply(score_cols, function(sc) {
    fit <- aov(epi@meta.data[[sc]] ~ clusters); s <- summary(fit)[[1]]
    data.frame(program = prog_label(sc), F_stat = s[1, "F value"], p_value = s[1, "Pr(>F)"])
})
write.csv(do.call(rbind, anova_rows),
          file.path(OUT_DIR, "OneWay_ANOVA_per_program.csv"), row.names = FALSE)

mean_by_clu <- sapply(score_cols, function(sc) tapply(epi@meta.data[[sc]], clusters, mean, na.rm = TRUE))
colnames(mean_by_clu) <- prog_label(colnames(mean_by_clu))
write.csv(round(mean_by_clu, 4), file.path(OUT_DIR, "mean_program_score_per_cluster.csv"))

# ==============================================================================
# (4) 통합 verdict heatmap — program × cluster, cluster 간 z-score
# ==============================================================================
z <- t(scale(mean_by_clu))
hm_df <- expand.grid(program = rownames(z), cluster = colnames(z),
                     KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
hm_df$z <- as.vector(z)
hm_df$cluster <- factor(hm_df$cluster, levels = colnames(z))
hm_df$program <- factor(hm_df$program, levels = rev(rownames(z)))
p_hm <- ggplot(hm_df, aes(cluster, program, fill = z)) +
    geom_tile(color = "grey85") + geom_text(aes(label = sprintf("%.1f", z)), size = 2.4) +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", midpoint = 0,
                         name = "z (across\nclusters)") +
    labs(title = "Published DNPC programs — relative reproduction per cluster",
         subtitle = "AddModuleScore mean per cluster, z-scored across clusters", x = NULL, y = NULL) +
    theme_minimal(base_size = 10) + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "Verdict_program_x_cluster_heatmap.png"),
       plot = p_hm, width = 26, height = 22, units = "cm", dpi = 220, bg = "white")

top_cluster <- data.frame(
    program     = colnames(mean_by_clu),
    top_cluster = rownames(mean_by_clu)[apply(mean_by_clu, 2, which.max)],
    top_mean    = round(apply(mean_by_clu, 2, max), 4),
    runnerup    = apply(mean_by_clu, 2, function(v) rownames(mean_by_clu)[order(v, decreasing = TRUE)[2]]))
write.csv(top_cluster, file.path(OUT_DIR, "Verdict_top_cluster_per_program.csv"), row.names = FALSE)

# ==============================================================================
# (5) DISSOCIATION QC 교차검증 — AP-1/IEG가 prep artifact인가?
#   AP-1 ↔ heat-shock ↔ QC(percent.mt/nCount/nFeature) Spearman 상관 + per-sample 분포.
# ==============================================================================
qc_candidates <- c("nCount_RNA", "nFeature_RNA",
                   grep("^percent", colnames(epi@meta.data), value = TRUE))
focus_scores <- intersect(c("AP1_IEG_stress_score", "HeatShock_dissoc_ctl_score",
    "Tang2022_CRPC_SCL_score", "He2024_KRT7_DNPC_score", "Labrecque2019_DNPC_pos_score",
    "POU2F3_tuft_score"), score_cols)
cor_vars <- c(focus_scores, intersect(qc_candidates, colnames(epi@meta.data)))
M <- as.data.frame(epi@meta.data[, cor_vars, drop = FALSE])
cmat <- suppressWarnings(cor(M, use = "pairwise.complete.obs", method = "spearman"))
write.csv(round(cmat, 3), file.path(OUT_DIR, "QC_dissociation_correlation.csv"))

if ("AP1_IEG_stress_score" %in% score_cols) {
    ap1_ps <- tapply(epi$AP1_IEG_stress_score, epi$orig.ident, mean, na.rm = TRUE)
    write.csv(data.frame(sample = names(ap1_ps), AP1_IEG_mean = round(as.numeric(ap1_ps), 4)),
              file.path(OUT_DIR, "AP1_per_sample_mean.csv"), row.names = FALSE)
}

message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
cat("\n=== top cluster per program ===\n"); print(top_cluster)
cat("\n=== AP-1/IEG dissociation cross-check (Spearman) ===\n")
print(round(cmat[intersect(c("AP1_IEG_stress_score","Tang2022_CRPC_SCL_score","He2024_KRT7_DNPC_score"),
            rownames(cmat)),
          intersect(c("HeatShock_dissoc_ctl_score", qc_candidates), colnames(cmat)), drop = FALSE], 3))
