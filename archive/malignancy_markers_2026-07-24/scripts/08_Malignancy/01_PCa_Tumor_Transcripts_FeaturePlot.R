# 08_Malignancy/01_PCa_Tumor_Transcripts_FeaturePlot.R
# 전립선암 특이 전사체를 환자별 split FeaturePlot 으로 시각화.
# 유전자 1개당 PNG 1장 = [클러스터 annotation | CRPC1 | CRPC2 | CRPC3] 4패널.
#
# 입력 : Results/05_Epithelial_Downstream/epi_annotated.rds
# 출력 : Results/08_Malignancy/Figures/Fig_<GENE>_by_patient.png
#
# 주의 : CNV 결과와 마찬가지로 이 플롯도 malignancy 를 단정하는 argument 근거로
#        쓰지 않는다. 클러스터별 수렴 증거(환자 공유 패턴·정상 아틀라스 대응물)와
#        함께 읽을 것.

suppressMessages({
    library(Seurat); library(ggplot2); library(patchwork); library(viridis)
})
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
FIG_DIR <- "Results/08_Malignancy/Figures"
ASSAY   <- "SCT"        # RNA 의 data 레이어가 비어 있어(counts 전용) SCT 사용
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

# 전립선암 특이/연관 전사체 ----
# (종양 ↑) 암에서 발현 증가  /  (종양 ↓) 종양억제자, 암에서 소실
GENES <- c(
    PCA3     = "PCa 특이 lncRNA — 임상 소변 바이오마커 (종양 ↑)",
    AMACR    = "P504S — 병리 IHC 표준 PCa 마커 (종양 ↑)",
    SCHLAP1  = "공격성 PCa lncRNA (종양 ↑)",
    ERG      = "TMPRSS2-ERG fusion 산물 (fusion+ 종양 ↑)",
    TMEFF2   = "전립선 제한 발현, PCa 연관 (종양 ↑)",
    GOLM1    = "PCa 연관 당단백 (종양 ↑)",
    MYC      = "8q24 증폭 흔함 (종양 ↑)",
    `NKX3-1` = "전립선 종양억제자 — PCa 에서 소실 (종양 ↓)",
    PTEN     = "종양억제자 — PCa 에서 소실 (종양 ↓)"
)

o <- readRDS(EPI_RDS)
DefaultAssay(o) <- ASSAY
o$orig.ident <- factor(o$orig.ident, levels = SAMPLES)

missing <- setdiff(names(GENES), rownames(o[[ASSAY]]))
if (length(missing)) message("assay 에 없는 유전자(건너뜀): ", paste(missing, collapse = ", "))
GENES <- GENES[!names(GENES) %in% missing]

# 모든 패널이 동일한 UMAP 프레임을 쓰도록 임베딩 범위 고정
emb <- Embeddings(o, "umap")
xr  <- range(emb[, 1]); yr <- range(emb[, 2])
pad <- 0.03 * c(-diff(xr), diff(xr))
xr  <- xr + pad; yr <- yr + 0.03 * c(-diff(yr), diff(yr))

base_theme <- theme(
    plot.title   = element_text(size = 13, face = "bold"),
    axis.title   = element_text(size = 9),
    axis.text    = element_text(size = 8)
)

# 좌측 패널: 클러스터 annotation (색약 친화 팔레트)
p_anno <- DimPlot(o, group.by = "annotation", label = TRUE, repel = TRUE,
                  pt.size = 0.3, label.size = 3.4,
                  cols = utils_cb_palette(nlevels(factor(o$annotation)))) +
    ggtitle("epithelial annotation") +
    coord_equal(xlim = xr, ylim = yr) + NoLegend() + base_theme

for (g in names(GENES)) {
    expr <- FetchData(o, vars = g, layer = "data")[, 1]
    lim  <- range(expr, na.rm = TRUE)
    if (!is.finite(diff(lim)) || diff(lim) <= 0) {
        message("[", g, "] 발현 변이 없음 — 건너뜀"); next
    }

    # 환자별 패널: 전체 대비 동일 색스케일(lim) + 동일 프레임(xr/yr)
    panels <- lapply(SAMPLES, function(s) {
        cells <- colnames(o)[o$orig.ident == s]
        FeaturePlot(o, features = g, cells = cells, order = TRUE, pt.size = 0.35) +
            scale_color_viridis_c(option = "viridis", limits = lim, name = "expr") +
            ggtitle(sprintf("%s  (n=%s)", s, format(length(cells), big.mark = ","))) +
            coord_equal(xlim = xr, ylim = yr) + base_theme
    })

    combined <- wrap_plots(c(list(p_anno), panels), nrow = 1,
                           widths = c(1.25, 1, 1, 1)) +
        plot_layout(guides = "collect") +
        plot_annotation(
            title    = g,
            subtitle = unname(GENES[[g]]),
            caption  = sprintf("assay=%s (data), 색스케일 3환자 공통 [%.2f, %.2f]",
                               ASSAY, lim[1], lim[2]),
            theme    = theme(plot.title    = element_text(size = 20, face = "bold"),
                             plot.subtitle = element_text(size = 12),
                             plot.caption  = element_text(size = 9, colour = "grey35"))
        )

    out <- file.path(FIG_DIR, sprintf("Fig_%s_by_patient.png", g))
    ggsave(out, combined, width = 22, height = 6.2, dpi = 200, bg = "white")
    message("saved → ", out)
}

message("완료 → ", FIG_DIR)
