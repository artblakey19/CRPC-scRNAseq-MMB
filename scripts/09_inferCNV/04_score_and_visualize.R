#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 04 : per-cluster CNV burden 정량 + 시각화
# -----------------------------------------------------------------------------
# run.final.infercnv_obj 의 denoise 잔차발현(≈1 중심)에서
#   세포별 CNV score = mean((expr - 1)^2)  (유전체 전반 이탈량).
# 상피 fine cluster / 환자별 요약표 + boxplot. reference 그룹은 baseline 으로 표시.
# =============================================================================
suppressPackageStartupMessages({ library(infercnv); library(ggplot2) })

OUT <- "Results/09_inferCNV"
VIS <- file.path(OUT, "Visualization")
dir.create(VIS, recursive = TRUE, showWarnings = FALSE)
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")

# Okabe-Ito 색약 친화 팔레트 (discrete)
cb <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
        "#999999", "#000000", "#332288", "#117733", "#88CCEE", "#DDCC77", "#AA4499")

rows <- list()
for (S in SAMPLES) {
    f <- file.path(OUT, S, "infercnv", "run.final.infercnv_obj")
    if (!file.exists(f)) { cat("skip (no obj):", S, "\n"); next }
    o    <- readRDS(f)
    expr <- o@expr.data                          # gene x cell, denoise 잔차
    score <- colMeans((expr - 1)^2)              # 세포별 CNV score

    grp <- setNames(rep(NA_character_, ncol(expr)), colnames(expr))
    for (g in names(o@observation_grouped_cell_indices))
        grp[o@observation_grouped_cell_indices[[g]]] <- g
    for (g in names(o@reference_grouped_cell_indices))
        grp[o@reference_grouped_cell_indices[[g]]] <- paste0("[ref] ", g)

    rows[[S]] <- data.frame(sample = S, cell = colnames(expr),
                            cluster = unname(grp), cnv_score = score, row.names = NULL)
}
df <- do.call(rbind, rows)
write.csv(df, file.path(OUT, "cnv_score_percell.csv"), row.names = FALSE)

# per-cluster 요약 (mean / median / n)
agg  <- aggregate(cnv_score ~ sample + cluster, df,
                  function(x) c(mean = mean(x), median = median(x), n = length(x)))
summ <- do.call(data.frame, agg)
summ <- summ[order(-summ$cnv_score.mean), ]
write.csv(summ, file.path(OUT, "cnv_score_by_cluster.csv"), row.names = FALSE)
cat("\n=== per-cluster CNV score (top) ===\n"); print(head(summ, 25))

# boxplot: reference baseline 은 아래쪽으로 모아서 표시
is_ref <- grepl("^\\[ref\\]", df$cluster)
lev <- c(sort(unique(df$cluster[!is_ref])), sort(unique(df$cluster[is_ref])))
df$cluster <- factor(df$cluster, levels = rev(lev))
p <- ggplot(df, aes(cluster, cnv_score, fill = sample)) +
    geom_boxplot(outlier.size = 0.2, lwd = 0.3, position = position_dodge(preserve = "single")) +
    scale_fill_manual(values = cb) +
    coord_flip() +
    labs(title = "inferCNV per-cell CNV score by epithelial cluster",
         subtitle = "score = mean((denoised expr - 1)^2);  [ref] = normal reference baseline",
         y = "CNV score", x = NULL) +
    theme_bw(base_size = 9)
ggsave(file.path(VIS, "cnv_score_by_cluster.pdf"), p, width = 8, height = 7)
cat("\n[04] done ->", VIS, "\n")
