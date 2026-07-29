# Epithelial_LuminalBasal_LiteratureMarkers_FeaturePlot.R
#
# Song 2022 이외의 문헌에서 "명시적으로 표로 제시된" prostate luminal / basal
# (+ club / hillock) 마커를 우리 상피 데이터에 FeaturePlot 으로 얹어 확인한다.
#
# 마커 출처 (모두 human prostate scRNA-seq, 본문/그림에 유전자명이 직접 표기됨)
#   [Henry2018]  Henry GH et al. Cell Rep 2018;25:3530-3542.e5  (PMID 30566875)
#                normal human prostate + prostatic urethra atlas.
#                "luminal KLK3+, basal KRT14+, OE1(club) SCGB1A1+, OE2(hillock) KRT13+"
#                BE: KRT5/KRT14/TP63 + NOTCH4/LTBP2/DKK1
#                LE: KLK3/ACPP(=ACP3)/MSMB + GP2/NEFH/NPY
#                NE: CHGA/CHGB (+ LY6D/SCGB3A1/PSCA 가 'other epithelia' 에 동반)
#                → 이 프로젝트의 BE / LE / Club / Hillock 라벨 체계의 원출처.
#   [Joseph2020] Joseph DB et al. Prostate 2020;80:872-884       (PMID 32497356)
#                urethral luminal(=human club/hillock 등가) 세포. progenitor 마커
#                TACSTD2(TROP2)/PSCA/KRT4/LY6D/KRT19 를 dot plot 으로 제시.
#   [Hirz2023]   Hirz T et al. Nat Commun 2023;14:663            (PMID 36750562)
#                종양 조직에서도 basal/luminal/club/hillock 4개 상피 아형이
#                동일한 key-marker 로 재현됨을 보인 human PCa TME 아틀라스.
#   [Pitzen2025] Pitzen SP et al. PNAS 2025;122:e2415308122      (PMID 39913208)
#                Henry2018 마커를 서로 겹치지 않게 정리한 BPECT gene set 으로
#                CRPC 에서 basal/club/hillock 정체성이 '동시에' 올라감을 보고.
#                readout 유전자: basal KRT14, club PI3, hillock KRT13.
#                DNPC 코호트인 우리 데이터와 가장 직접적으로 맞물리는 논문.
#   [Canonical]  교과서적 basal/luminal 마커(위 논문들이 공통 인용):
#                basal KRT5/KRT14/TP63/KRT15/ITGA6/NGFR(CD271)/PDPN,
#                luminal KRT8/KRT18/AR/NKX3-1/KLK2/KLK3/TMPRSS2/FOLH1/DPP4(CD26).
#
# 입력 : Results/05_Epithelial_Downstream/epi_annotated.rds
# 출력 : Results/05_Epithelial_Downstream/LuminalBasal_LiteratureMarkers/
#
# assay : 파이프라인 관례대로 SCT(정규화 data) 사용.
# 주의  : 이 그림은 lineage identity(luminal/basal 계열) 확인용이며,
#         benign vs malignant 판정 근거로 쓰지 않는다.

suppressMessages({
    library(Seurat); library(ggplot2); library(patchwork); library(dplyr)
})
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/LuminalBasal_LiteratureMarkers"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# 마커 테이블 (gene, cell_type, source)
# ACPP 는 현재 HGNC 심볼이 ACP3 이라 ACP3 로 적는다.
# ============================================================

marker_tbl <- rbind(
    # ---- Basal epithelia ----
    data.frame(cell_type = "Basal",   source = "Henry2018",
               gene = c("KRT5", "KRT14", "TP63", "NOTCH4", "LTBP2", "DKK1")),
    data.frame(cell_type = "Basal",   source = "Canonical",
               gene = c("KRT15", "DST", "COL17A1", "ITGA6", "ITGB4",
                        "NGFR", "PDPN", "BNC1", "MIR205HG", "CDH3")),

    # ---- Luminal epithelia ----
    data.frame(cell_type = "Luminal", source = "Henry2018",
               gene = c("KLK3", "ACP3", "MSMB", "GP2", "NEFH", "NPY", "DHRS7")),
    data.frame(cell_type = "Luminal", source = "Canonical",
               gene = c("KRT8", "KRT18", "AR", "NKX3-1", "KLK2", "TMPRSS2",
                        "FOLH1", "TGM4", "DPP4", "CD38", "PLA2G2A", "RDH11")),

    # ---- Club (Henry OE1) ----
    data.frame(cell_type = "Club",    source = "Henry2018",
               gene = c("SCGB1A1", "SCGB3A1", "LY6D", "PSCA")),
    data.frame(cell_type = "Club",    source = "Joseph2020",
               gene = c("TACSTD2", "KRT4", "PIGR", "MMP7", "LTF", "CP")),
    data.frame(cell_type = "Club",    source = "Hirz2023",
               gene = c("WFDC2", "LCN2", "OLFM4", "CRABP2")),
    data.frame(cell_type = "Club",    source = "Pitzen2025",
               gene = "PI3"),

    # ---- Hillock (Henry OE2) ----
    data.frame(cell_type = "Hillock", source = "Henry2018",
               gene = c("KRT13", "AKR1C1", "AKR1C2", "KRT19")),
    data.frame(cell_type = "Hillock", source = "Joseph2020",
               gene = c("KRT7", "RARRES1", "LYPD3", "AQP3")),
    data.frame(cell_type = "Hillock", source = "Hirz2023",
               gene = c("S100A2", "SERPINB3", "SERPINB4", "APOBEC3A", "CLDN4"))
)

CELL_TYPES <- c("Basal", "Luminal", "Club", "Hillock")
marker_tbl$cell_type <- factor(marker_tbl$cell_type, levels = CELL_TYPES)

# 클러스터 판별에 가장 자주 쓰이는 core 마커 (환자별 split 용)
CORE_GENES <- c("KRT5", "TP63", "KRT14",          # basal
                "KLK3", "NKX3-1", "AR",           # luminal
                "SCGB1A1", "PIGR",                # club
                "KRT13", "KRT19")                 # hillock

# ============================================================
# Load
# ============================================================

stopifnot("epi_annotated.rds not found" = file.exists(EPI_RDS))
epi <- readRDS(EPI_RDS)
epi$orig.ident <- factor(epi$orig.ident, levels = SAMPLES)
DefaultAssay(epi) <- "SCT"

lvls <- readLines("Results/05_Epithelial_Downstream/Annotation/label_levels.txt")
epi$plot_group <- factor(epi$annotation, levels = lvls)

# SCT 에 없는 유전자는 표에서 제외하고 상태를 기록
marker_tbl$in_SCT <- marker_tbl$gene %in% rownames(epi[["SCT"]])
marker_tbl$in_RNA <- marker_tbl$gene %in% rownames(epi[["RNA"]])

write.table(marker_tbl, file.path(OUT_DIR, "marker_gene_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

dropped <- marker_tbl[!marker_tbl$in_SCT, ]
if (nrow(dropped)) {
    message("SCT 에 없어 플롯에서 제외: ",
            paste(sprintf("%s(%s)", dropped$gene, dropped$cell_type), collapse = ", "))
}
mk <- marker_tbl[marker_tbl$in_SCT, ]

# ============================================================
# 공통 UMAP 프레임 / 테마
# ============================================================

emb <- Embeddings(epi, "umap")
xr  <- range(emb[, 1]) + 0.03 * c(-1, 1) * diff(range(emb[, 1]))
yr  <- range(emb[, 2]) + 0.03 * c(-1, 1) * diff(range(emb[, 2]))

base_theme <- theme(plot.title  = element_text(size = 11, face = "bold"),
                    axis.title  = element_text(size = 7),
                    axis.text   = element_text(size = 6),
                    legend.text = element_text(size = 6),
                    legend.key.width = unit(0.25, "cm"))

p_anno <- DimPlot(epi, group.by = "plot_group", label = TRUE, repel = TRUE,
                  pt.size = 0.25, label.size = 3,
                  cols = utils_cb_palette(nlevels(epi$plot_group))) +
    ggtitle("epithelial annotation") +
    coord_equal(xlim = xr, ylim = yr) + NoLegend() + base_theme

ggsave(file.path(OUT_DIR, "UMAP_annotation_reference.png"), p_anno,
       width = 16, height = 16, units = "cm", dpi = 200, bg = "white")

# 유전자 1개 FeaturePlot (viridis, 연속형 스케일)
gene_umap <- function(g, src) {
    FeaturePlot(epi, features = g, order = TRUE, pt.size = 0.25) +
        scale_color_viridis_c(option = "viridis", name = "expr") +
        ggtitle(sprintf("%s  [%s]", g, src)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}

# ============================================================
# 1) cell type 별 FeaturePlot 그리드
# ============================================================

for (ct in CELL_TYPES) {
    sub <- mk[mk$cell_type == ct, ]
    if (!nrow(sub)) next

    panels <- c(list(p_anno), Map(gene_umap, sub$gene, sub$source))
    ncol   <- 4
    nrow   <- ceiling(length(panels) / ncol)

    comb <- wrap_plots(panels, ncol = ncol) +
        plot_annotation(
            title    = sprintf("%s markers — literature panel (Song 2022 제외)", ct),
            subtitle = paste0("좌상단=annotation 참조 · 각 패널 [출처] 표기 · ",
                              "출처: Henry 2018 Cell Rep / Joseph 2020 Prostate / ",
                              "Hirz 2023 Nat Commun / Pitzen 2025 PNAS"),
            caption  = "assay=SCT(data), 색=viridis. lineage identity 확인용 — malignancy 판정 근거 아님.",
            theme    = theme(plot.title    = element_text(size = 16, face = "bold"),
                             plot.subtitle = element_text(size = 9),
                             plot.caption  = element_text(size = 8, colour = "grey35"))
        )

    out <- file.path(OUT_DIR, sprintf("FeaturePlot_%s.png", ct))
    ggsave(out, comb, width = 4.2 * ncol, height = 4.4 * nrow,
           dpi = 200, bg = "white", limitsize = FALSE)
    message("Wrote: ", out)
}

# ============================================================
# 2) core 마커 환자별 split FeaturePlot
#    (유전자별로 3환자 공통 색스케일 → 환자 간 비교 가능)
# ============================================================

core <- CORE_GENES[CORE_GENES %in% rownames(epi[["SCT"]])]
expr_core <- FetchData(epi, vars = c(core, "orig.ident"), layer = "data")

gene_panel_patient <- function(g, s, lim) {
    cells <- colnames(epi)[epi$orig.ident == s]
    FeaturePlot(epi, features = g, cells = cells, order = TRUE, pt.size = 0.3) +
        scale_color_viridis_c(option = "viridis", limits = lim, name = "expr") +
        ggtitle(sprintf("%s — %s", g, s)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}

panels <- list(p_anno, patchwork::plot_spacer(), patchwork::plot_spacer())
for (g in core) {
    mx  <- suppressWarnings(max(expr_core[[g]], na.rm = TRUE))
    lim <- if (!is.finite(mx) || mx <= 0) c(0, 1) else c(0, mx)
    panels <- c(panels, lapply(SAMPLES, function(s) gene_panel_patient(g, s, lim)))
}

comb <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
        title    = "Core luminal / basal / club / hillock markers — patient-split",
        subtitle = "행=마커(유전자별 3환자 공통 색스케일), 열=CRPC1/2/3 · 좌상단=annotation 참조",
        caption  = "assay=SCT(data). 출처: Henry 2018 / Joseph 2020 / Hirz 2023 / Pitzen 2025.",
        theme    = theme(plot.title    = element_text(size = 17, face = "bold"),
                         plot.subtitle = element_text(size = 10),
                         plot.caption  = element_text(size = 8, colour = "grey35"))
    )

out <- file.path(OUT_DIR, "FeaturePlot_CoreMarkers_byPatient.png")
ggsave(out, comb, width = 14, height = 4.4 * (length(core) + 1),
       dpi = 200, bg = "white", limitsize = FALSE)
message("Wrote: ", out)

# ============================================================
# 3) DotPlot — 클러스터별 요약 (어느 클러스터가 luminal/basal 인지 한 장에)
# ============================================================

dot_genes <- mk$gene[!duplicated(mk$gene)]
dot_ct    <- mk$cell_type[!duplicated(mk$gene)]

p_dot <- DotPlot(epi, features = dot_genes, group.by = "plot_group",
                 assay = "SCT", dot.scale = 5) +
    scale_color_viridis_c(option = "viridis", name = "avg expr\n(scaled)") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
          axis.text.y = element_text(size = 9)) +
    labs(x = NULL, y = NULL,
         title = "Literature luminal/basal/club/hillock markers by epithelial cluster",
         subtitle = paste("유전자 순서 =", paste(CELL_TYPES, collapse = " → "),
                          "· 출처: Henry 2018 / Joseph 2020 / Hirz 2023 / Pitzen 2025"))

ggsave(file.path(OUT_DIR, "DotPlot_AllMarkers_byCluster.png"), p_dot,
       width = 0.28 * length(dot_genes) + 6, height = 10,
       dpi = 200, bg = "white", limitsize = FALSE)

# ============================================================
# 4) cell type 별 module score (AddModuleScore) + UMAP + violin
#    개별 유전자 dropout 에 덜 흔들리는 요약 지표
# ============================================================

score_cols <- character(0)

for (ct in CELL_TYPES) {
    genes <- unique(mk$gene[mk$cell_type == ct])
    if (length(genes) < 5) { message("Skip module score: ", ct); next }

    nm <- paste0("Lit_", ct, "_Score")
    epi <- AddModuleScore(epi, features = list(genes), name = nm,
                          assay = "SCT", ctrl = 50)
    epi@meta.data[[nm]] <- epi@meta.data[[paste0(nm, "1")]]
    epi@meta.data[[paste0(nm, "1")]] <- NULL
    score_cols <- c(score_cols, nm)
}

score_umaps <- lapply(score_cols, function(sc) {
    FeaturePlot(epi, features = sc, order = TRUE, pt.size = 0.25) +
        scale_color_viridis_c(option = "viridis", name = "score") +
        ggtitle(sub("^Lit_(.*)_Score$", "\\1 signature score", sc)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
})

comb <- wrap_plots(c(list(p_anno), score_umaps), ncol = 3) +
    plot_annotation(
        title    = "Literature cell-type signature scores (AddModuleScore)",
        subtitle = "Henry 2018 / Joseph 2020 / Hirz 2023 / Pitzen 2025 마커 union, cell type 별",
        caption  = "assay=SCT(data), ctrl=50. 개별 유전자 dropout 에 덜 민감한 요약 지표.",
        theme    = theme(plot.title    = element_text(size = 16, face = "bold"),
                         plot.subtitle = element_text(size = 9),
                         plot.caption  = element_text(size = 8, colour = "grey35"))
    )

ggsave(file.path(OUT_DIR, "FeaturePlot_SignatureScores.png"), comb,
       width = 13, height = 4.4 * ceiling((length(score_cols) + 1) / 3),
       dpi = 200, bg = "white", limitsize = FALSE)

vln <- lapply(score_cols, function(sc) {
    df <- data.frame(group = epi$plot_group, score = epi@meta.data[[sc]])
    ggplot(df, aes(group, score, fill = group)) +
        geom_violin(scale = "width") +
        geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white", alpha = 0.6) +
        scale_fill_manual(values = utils_cb_palette(nlevels(df$group))) +
        ggtitle(sub("^Lit_(.*)_Score$", "\\1 signature score", sc)) +
        labs(x = NULL, y = "module score") + NoLegend() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              plot.title  = element_text(size = 12, face = "bold"))
})

ggsave(file.path(OUT_DIR, "Violin_SignatureScores_byCluster.png"),
       wrap_plots(vln, ncol = 2), width = 26, height = 9 * ceiling(length(vln) / 2),
       units = "cm", dpi = 200, bg = "white", limitsize = FALSE)

# 클러스터별 평균 점수 표 (숫자 근거)
score_summary <- epi@meta.data %>%
    group_by(plot_group) %>%
    summarise(n = n(), across(all_of(score_cols), ~ round(mean(.x), 4)), .groups = "drop")

write.csv(score_summary, file.path(OUT_DIR, "signature_score_mean_byCluster.csv"),
          row.names = FALSE)
print(as.data.frame(score_summary))

message("Done. Outputs in: ", OUT_DIR)
