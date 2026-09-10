#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 06 : inferCNV per-cell 지표를 상피 UMAP 에 투영
# -----------------------------------------------------------------------------
# 05 객체(epi_annotated.rds)의 통합 UMAP 에 inferCNV per-cell CNV 지표를 얹는다.
#   - inferCNV_PGA   : 세포 소속 subcluster 의 PGA (05, HMM Bayesian 기반)
#   - inferCNV_score : mean((denoise expr - 1)^2) (04)
# 연속형 → viridis. 환자 특이적 CNV 이므로 split.by=orig.ident 도 생성.
# ⚠ PGA/score 는 reference 오염으로 노이즈 큼(README_CAVEAT) → 해석 주의.
# =============================================================================
suppressPackageStartupMessages({
    library(Seurat); library(infercnv); library(ggplot2); library(patchwork)
})

OUT <- "Results/09_inferCNV"
VIS <- file.path(OUT, "Visualization")
dir.create(VIS, showWarnings = FALSE, recursive = TRUE)
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")

# Okabe-Ito 계열 색약 친화 팔레트 (discrete, 13 clusters)
cb <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
        "#999999","#000000","#332288","#117733","#88CCEE","#AA4499")

epi <- readRDS("Results/05_Epithelial_Downstream/epi_annotated.rds")
red <- if ("umap" %in% Reductions(epi)) "umap" else Reductions(epi)[1]
epi_key <- paste(as.character(epi$orig.ident),
                 sub("^P[0-9]+_", "", colnames(epi)), sep = "|")

# --- per-cell PGA : subcluster PGA(csv) × 세포 인덱스(object) -----------------
pga_sc <- read.csv(file.path(OUT, "pga_by_subcluster.csv"), stringsAsFactors = FALSE)
cellpga <- setNames(numeric(0), character(0))
for (S in SAMPLES) {
    o  <- readRDS(file.path(OUT, S, "infercnv", "run.final.infercnv_obj"))
    cn <- colnames(o@expr.data)
    sc <- o@tumor_subclusters$subclusters
    ps <- pga_sc[pga_sc$sample == S, ]
    pmap <- setNames(ps$pga, paste(ps$group, ps$subcl, sep = "\r"))
    for (g in names(sc)) for (s in names(sc[[g]])) {
        v  <- pmap[[paste(g, s, sep = "\r")]]
        bc <- cn[sc[[g]][[s]]]
        cellpga[paste(S, bc, sep = "|")] <- v
    }
    rm(o); gc()
}

# --- per-cell mean-sq score (04) --------------------------------------------
cs  <- read.csv(file.path(OUT, "cnv_score_percell.csv"), stringsAsFactors = FALSE)
csv <- setNames(cs$cnv_score, paste(cs$sample, cs$cell, sep = "|"))

epi$inferCNV_PGA   <- unname(cellpga[epi_key])
epi$inferCNV_score <- unname(csv[epi_key])
cat(sprintf("mapped PGA: %d/%d cells (NA %d) | score: %d/%d\n",
            sum(!is.na(epi$inferCNV_PGA)), ncol(epi), sum(is.na(epi$inferCNV_PGA)),
            sum(!is.na(epi$inferCNV_score)), ncol(epi)))

# 대비를 위해 상위 1% 는 clip (노이즈 극단값이 색 스케일 잠식 방지)
cap <- as.numeric(quantile(epi$inferCNV_PGA, 0.99, na.rm = TRUE))
epi$inferCNV_PGA_capped <- pmin(epi$inferCNV_PGA, cap)

vir <- function(name) scale_color_viridis_c(option = "viridis", na.value = "grey88", name = name)

# --- 1) overview: annotation + PGA + score ----------------------------------
p_ann <- DimPlot(epi, reduction = red, group.by = "annotation", cols = cb,
                 label = TRUE, repel = TRUE, label.size = 3) +
         ggtitle("Epithelial clusters") + theme(legend.position = "none")
p_pga <- FeaturePlot(epi, "inferCNV_PGA_capped", reduction = red, order = TRUE, pt.size = 0.2) +
         vir("PGA") + ggtitle("inferCNV PGA (per cell)")
p_sco <- FeaturePlot(epi, "inferCNV_score", reduction = red, order = TRUE, pt.size = 0.2) +
         vir("score") + ggtitle("inferCNV mean-sq score")
ggsave(file.path(VIS, "umap_infercnv_overview.png"),
       p_ann + p_pga + p_sco + plot_layout(nrow = 1), width = 18, height = 5.5, dpi = 150)

# --- 2) PGA, 환자별 split ----------------------------------------------------
# ⚠ FeaturePlot(split.by=) 는 legend 를 지우고 feature 명을 우측 y축 라벨로만 남긴다.
#   patchwork guides="collect" + legend.position 으로 공통 colorbar 를 되살린다.
p_split <- FeaturePlot(epi, "inferCNV_PGA_capped", reduction = red,
                       split.by = "orig.ident", order = TRUE, pt.size = 0.25) &
           vir("PGA")
p_split <- p_split + plot_layout(guides = "collect") &
           theme(legend.position = "right")
ggsave(file.path(VIS, "umap_pga_by_patient.png"), p_split,
       width = 17, height = 5.2, dpi = 150)

# --- 3) score, 환자별 split --------------------------------------------------
p_split2 <- FeaturePlot(epi, "inferCNV_score", reduction = red,
                        split.by = "orig.ident", order = TRUE, pt.size = 0.25) &
            vir("score")
p_split2 <- p_split2 + plot_layout(guides = "collect") &
            theme(legend.position = "right")
ggsave(file.path(VIS, "umap_score_by_patient.png"), p_split2,
       width = 17, height = 5.2, dpi = 150)

cat("[06] done ->", VIS, "\n")
