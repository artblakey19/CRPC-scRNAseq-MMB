# Epithelial_NEmarkers_FeaturePlot.R
# 신경내분비(NE) 마커를 환자별 split FeaturePlot 으로 시각화 → DNPC 의
# "double-negative" 중 NE-negativity 를 확증하는 용도.
#
# 입력 : Results/05_Epithelial_Downstream/epi_annotated.rds
# 출력 : Results/05_Epithelial_Downstream/Annotation/FeaturePlot_NEmarkers_byPatient.png
#
# assay 주의 : 파이프라인 관례대로 SCT(정규화 data)를 그대로 사용. INSM1 은 SCT model
#   에서 제외돼(발현 극소, RNA counts>0 = 3세포) SCT 에 없으므로 이 그림에서 제외한다.

suppressMessages({
    library(Seurat); library(ggplot2); library(patchwork); library(viridis)
})
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# NE 마커 (표준 neuroendocrine/NEPC 패널) — INSM1 은 SCT 에 없어 제외
NE_GENES <- c("CHGA", "CHGB", "SYP", "NCAM1", "ASCL1")

epi <- readRDS(EPI_RDS)
epi$orig.ident <- factor(epi$orig.ident, levels = SAMPLES)
DefaultAssay(epi) <- "SCT"

missing <- setdiff(NE_GENES, rownames(epi[["SCT"]]))
if (length(missing)) message("SCT 에 없는 유전자(건너뜀): ", paste(missing, collapse = ", "))
NE_GENES <- NE_GENES[!NE_GENES %in% missing]

# 정량 요약(콘솔) — assay 무관 근거
expr_all <- FetchData(epi, vars = c(NE_GENES, "orig.ident"), layer = "data")
message("=== NE 마커 발현 세포수(>0) / ", ncol(epi), " ===")
for (g in NE_GENES) {
    tb <- tapply(expr_all[[g]] > 0, expr_all$orig.ident, sum)
    message(sprintf("  %-6s total=%4d (%.2f%%)  CRPC1=%d CRPC2=%d CRPC3=%d",
                    g, sum(expr_all[[g]] > 0), 100 * mean(expr_all[[g]] > 0),
                    tb[["CRPC1"]], tb[["CRPC2"]], tb[["CRPC3"]]))
}

# 모든 패널 공통 UMAP 프레임
emb <- Embeddings(epi, "umap")
xr  <- range(emb[, 1]) + 0.03 * c(-1, 1) * diff(range(emb[, 1]))
yr  <- range(emb[, 2]) + 0.03 * c(-1, 1) * diff(range(emb[, 2]))

base_theme <- theme(plot.title = element_text(size = 12, face = "bold"),
                    axis.title = element_text(size = 8), axis.text = element_text(size = 7))

# 클러스터 annotation 참조 패널(색약 친화 팔레트)
p_anno <- DimPlot(epi, group.by = "annotation", label = TRUE, repel = TRUE,
                  pt.size = 0.25, label.size = 3,
                  cols = utils_cb_palette(nlevels(factor(epi$annotation)))) +
    ggtitle("epithelial annotation") +
    coord_equal(xlim = xr, ylim = yr) + NoLegend() + base_theme

# 유전자별 3환자 패널 (유전자마다 3환자 공통 색스케일)
gene_panel <- function(g, s, lim) {
    cells <- colnames(epi)[epi$orig.ident == s]
    FeaturePlot(epi, features = g, cells = cells, order = TRUE, pt.size = 0.3) +
        scale_color_viridis_c(option = "viridis", limits = lim, name = "expr") +
        ggtitle(sprintf("%s — %s", g, s)) +
        coord_equal(xlim = xr, ylim = yr) + base_theme
}

# 행 배치: [annotation | spacer | spacer] 다음 유전자마다 [CRPC1 | CRPC2 | CRPC3]
panels <- list(p_anno, patchwork::plot_spacer(), patchwork::plot_spacer())
for (g in NE_GENES) {
    ev  <- expr_all[[g]]
    mx  <- suppressWarnings(max(ev, na.rm = TRUE))
    lim <- if (!is.finite(mx) || mx <= 0) c(0, 1) else c(0, mx)
    panels <- c(panels, lapply(SAMPLES, function(s) gene_panel(g, s, lim)))
}

combined <- wrap_plots(panels, ncol = 3) +
    plot_annotation(
        title    = "Neuroendocrine markers — patient-split (DNPC NE-negativity check)",
        subtitle = "행=NE 마커(유전자별 3환자 공통 색스케일), 열=CRPC1/2/3 · 좌상단=클러스터 참조",
        caption  = "assay=SCT (data); INSM1 은 SCT 제외 유전자(발현 극소)라 미표시. malignancy argument 근거로 사용 금지.",
        theme    = theme(plot.title    = element_text(size = 18, face = "bold"),
                         plot.subtitle = element_text(size = 11),
                         plot.caption  = element_text(size = 9, colour = "grey35"))
    )

out <- file.path(OUT_DIR, "FeaturePlot_NEmarkers_byPatient.png")
ggsave(out, combined, width = 15, height = 4.3 * (length(NE_GENES) + 1),
       dpi = 200, bg = "white", limitsize = FALSE)
message("Wrote: ", out)
