#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 10 : Henry-ref inferCNV 결과를 상피 UMAP 에 투영
# -----------------------------------------------------------------------------
# 06(stroma-ref)과 동일 방식으로, 08 Henry-ref 결과의 per-cell 지표를 상피 UMAP 에 얹는다.
#   - inferCNV_PGA_henry   : 세포 소속 subcluster 의 PGA (09: pga_by_subcluster_henryref)
#   - inferCNV_score_henry : mean((denoise expr - 1)^2)  (henryref final obj)
# observations(우리 상피)만 매핑; Henry reference 세포(HENRY_* )는 UMAP 에 없음.
# ⚠ Henry-ref 결과는 chemistry batch 로 confound(≈100% aneuploid) → 배치 아티팩트 확인용.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(infercnv); library(ggplot2); library(patchwork) })

BASE <- "Results/09_inferCNV"
OUT  <- file.path(BASE, "henryref")
VIS  <- file.path(OUT, "Visualization")
dir.create(VIS, showWarnings = FALSE, recursive = TRUE)
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
cb <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
        "#999999","#000000","#332288","#117733","#88CCEE","#AA4499")

epi <- readRDS("Results/05_Epithelial_Downstream/epi_annotated.rds")
red <- if ("umap" %in% Reductions(epi)) "umap" else Reductions(epi)[1]
epi_key <- paste(as.character(epi$orig.ident), sub("^P[0-9]+_", "", colnames(epi)), sep = "|")

pga_sc <- read.csv(file.path(OUT, "pga_by_subcluster_henryref.csv"), stringsAsFactors = FALSE)
cellpga <- setNames(numeric(0), character(0)); cellscore <- setNames(numeric(0), character(0))
for (S in SAMPLES) {
    o  <- readRDS(file.path(OUT, S, "infercnv", "run.final.infercnv_obj"))
    cn <- colnames(o@expr.data)
    # per-cell PGA : subcluster → cells
    sc <- o@tumor_subclusters$subclusters
    ps <- pga_sc[pga_sc$sample == S, ]
    pmap <- setNames(ps$pga, paste(ps$group, ps$subcl, sep = "\r"))
    for (g in names(sc)) for (s in names(sc[[g]])) {
        v <- pmap[[paste(g, s, sep = "\r")]]; bc <- cn[sc[[g]][[s]]]
        if (!is.null(v)) cellpga[paste(S, bc, sep = "|")] <- v
    }
    # per-cell mean-sq score (관측+ref 모두 계산; ref 는 epi 에 없어 자동 제외)
    sco <- colMeans((o@expr.data - 1)^2)
    cellscore[paste(S, cn, sep = "|")] <- sco
    rm(o); gc()
}
epi$inferCNV_PGA_henry   <- unname(cellpga[epi_key])
epi$inferCNV_score_henry <- unname(cellscore[epi_key])
cat(sprintf("mapped PGA %d/%d | score %d/%d cells\n",
            sum(!is.na(epi$inferCNV_PGA_henry)), ncol(epi),
            sum(!is.na(epi$inferCNV_score_henry)), ncol(epi)))

cap <- as.numeric(quantile(epi$inferCNV_PGA_henry, 0.99, na.rm = TRUE))
epi$PGA_henry_capped <- pmin(epi$inferCNV_PGA_henry, cap)
vir <- function(nm) scale_color_viridis_c(option = "viridis", na.value = "grey88", name = nm)

# 1) overview: annotation + PGA + score
p_ann <- DimPlot(epi, reduction = red, group.by = "annotation", cols = cb,
                 label = TRUE, repel = TRUE, label.size = 3) +
         ggtitle("Epithelial clusters") + theme(legend.position = "none")
p_pga <- FeaturePlot(epi, "PGA_henry_capped", reduction = red, order = TRUE, pt.size = 0.2) +
         vir("PGA") + ggtitle("Henry-ref inferCNV PGA")
p_sco <- FeaturePlot(epi, "inferCNV_score_henry", reduction = red, order = TRUE, pt.size = 0.2) +
         vir("score") + ggtitle("Henry-ref mean-sq score")
ggsave(file.path(VIS, "umap_henryref_overview.png"),
       p_ann + p_pga + p_sco + plot_layout(nrow = 1), width = 18, height = 5.5, dpi = 150)

# 2) PGA, 환자별 split
p_split <- FeaturePlot(epi, "PGA_henry_capped", reduction = red,
                       split.by = "orig.ident", order = TRUE, pt.size = 0.25) & vir("PGA")
ggsave(file.path(VIS, "umap_henryref_pga_by_patient.png"), p_split, width = 16, height = 5.2, dpi = 150)
cat("[10] done ->", VIS, "\n")
