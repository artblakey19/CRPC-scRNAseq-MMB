# 08_CopyKAT/02_CopyKAT_UMAP_Epithelial.R
# copyKAT 예측을 상피 UMAP(stage 05 annotation)에 투사.
# 라벨 annotation(좌) + 환자별 예측(우) 조합 플롯 포함.
#
# 입력 : Results/05_Epithelial_Downstream/epi_annotated.rds
#        Results/08_CopyKAT/copykat_prediction_all_samples.csv
# 출력 : Results/08_CopyKAT/Figures/
#          UMAP_epi_copykat_pred_by_patient.png
#          UMAP_epi_annotation_and_pred_by_patient.png   (라벨 좌 + 환자별 우)

suppressMessages({
    library(Seurat); library(dplyr); library(ggplot2); library(patchwork)
})
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
PRED_CSV <- "Results/08_CopyKAT/copykat_prediction_all_samples.csv"
FIG_DIR  <- "Results/08_CopyKAT/Figures"
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)

PRED_LVL <- c("aneuploid", "diploid", "not.defined")
PRED_COL <- c(aneuploid = "#D55E00", diploid = "#0072B2", not.defined = "#999999")
ord <- c("not.defined", "diploid", "aneuploid")

o <- readRDS(EPI_RDS)
pred <- read.csv(PRED_CSV, stringsAsFactors = FALSE)
pv <- setNames(pred$copykat.pred, pred$cell.names)[colnames(o)]
pv[is.na(pv)] <- "not.defined"
o$copykat.pred <- factor(pv, levels = PRED_LVL)

# 환자별 예측
p2 <- DimPlot(o, group.by = "copykat.pred", split.by = "orig.ident",
              cols = PRED_COL, pt.size = 0.3, order = ord) +
    ggtitle("copyKAT prediction on epithelial UMAP, by patient")
ggsave(file.path(FIG_DIR, "UMAP_epi_copykat_pred_by_patient.png"), p2, width = 18, height = 6.5, bg = "white")

# 라벨 annotation(좌) + 환자별 예측(우)
pa <- DimPlot(o, group.by = "annotation", pt.size = 0.3, label = TRUE, repel = TRUE,
              cols = utils_cb_palette(nlevels(factor(o$annotation)))) +
    ggtitle("epithelial annotation") + coord_equal() + NoLegend()
combined <- pa + p2 + patchwork::plot_layout(widths = c(1.3, 3))
ggsave(file.path(FIG_DIR, "UMAP_epi_annotation_and_pred_by_patient.png"),
       combined, width = 24, height = 6.8, bg = "white")

message("완료 → ", FIG_DIR)
