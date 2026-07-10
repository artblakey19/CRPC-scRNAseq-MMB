# Stage 8b — Is BE2 (cl7) the same entity as BE1 (cl0)? (rigorous test)
#
# seed 강건성 분석에서 BE1/BE2/BE3 가 취약 트리오로 나왔고, one-vs-rest 마커상
# BE2 는 약함(FC≤1.41, DEG 최소, ribosome-biogenesis 다수)이라 "BE1 의 over-split
# 조각"으로 의심됨(memory: filtered_epi_annotation). 여기서는 one-vs-rest 가 아니라
# BE1 vs BE2 를 직접 5 각도로 검정한다. 증거가 반대(=별개 실체)면 그대로 보고.
#
#   1. 직접 pairwise DE: BE1 vs BE2 (+ 대조 BE1-vs-BE3/BE4/BE6, BE2-vs-BE3)
#   2. BE1-vs-BE2 DEG 성격: cell-cycle / ribosome-biogenesis(state) vs lineage
#   3. Centroid 유사도: BE2 의 최근접 이웃이 BE1 인가 (발현 상관 + harmony 거리)
#   4. 해상도 계층: coarser res 에서 BE2 가 BE1 로 먼저 합쳐지는가
#   5. Seed co-assignment: 21 seed 에서 BE2 세포가 BE1 과 같은 클러스터로 가는 빈도
#
# Input : Results/04_Epithelial_Filtered/epi_filtered_clustered.rds  (canonical)
#         Results/08_Cluster_Stability/seed_labels.csv               (21 seed 라벨)
# Output: Results/08_Cluster_Stability/BE1_BE2_equivalence/

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(pheatmap)
library(viridisLite)

library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scripts/00_utils/scRNA_utils.R")

CANON_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
LABELS_CSV <- "Results/08_Cluster_Stability/seed_labels.csv"
OUT_DIR    <- "Results/08_Cluster_Stability/BE1_BE2_equivalence"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
sink_con <- file(file.path(OUT_DIR, "REPORT.txt"), open = "wt")
say <- function(...) { msg <- sprintf(...); message(msg); writeLines(msg, sink_con) }

CLUSTER_LABELS <- c(
  "0" = "BE1", "7" = "BE2", "10" = "BE3", "1" = "BE4", "2" = "BE5", "8" = "BE6",
  "6" = "OE", "9" = "LE(ARPC)", "3" = "Club", "5" = "Hillock1", "4" = "Hillock2",
  "11" = "Ionocyte"
)
lab_of <- function(cl) ifelse(as.character(cl) %in% names(CLUSTER_LABELS),
                              CLUSTER_LABELS[as.character(cl)], paste0("cl", cl))

epi <- readRDS(CANON_RDS)
epi$be <- factor(lab_of(epi$SCT_snn_res.0.4),
                 levels = c("BE1","BE2","BE3","BE4","BE5","BE6","OE",
                            "LE(ARPC)","Club","Hillock1","Hillock2","Ionocyte"))
Idents(epi) <- "be"
DefaultAssay(epi) <- "SCT"
epi <- PrepSCTFindMarkers(epi)

say("################ BE2 == BE1 ? rigorous evidence ################")
say("Cells: BE1=%d  BE2=%d  BE3=%d  BE4=%d  BE6=%d",
    sum(epi$be=="BE1"), sum(epi$be=="BE2"), sum(epi$be=="BE3"),
    sum(epi$be=="BE4"), sum(epi$be=="BE6"))

# =============================================================================
# 1. Direct pairwise DE (BE1 vs BE2) + controls ----
# =============================================================================
say("\n===== [1] Direct pairwise DE  (|log2FC|>0.5 & padj<0.05) =====")
pair_de <- function(a, b) {
    d <- FindMarkers(epi, ident.1 = a, ident.2 = b,
                     logfc.threshold = 0.1, min.pct = 0.1, verbose = FALSE)
    d$gene <- rownames(d)
    sig <- subset(d, p_val_adj < 0.05 & abs(avg_log2FC) > 0.5)
    list(all = d, sig = sig,
         n = nrow(sig),
         up_a = sum(sig$avg_log2FC > 0), up_b = sum(sig$avg_log2FC < 0))
}
pairs <- list(
    "BE1_vs_BE2" = c("BE1","BE2"),
    "BE1_vs_BE3" = c("BE1","BE3"),
    "BE1_vs_BE4" = c("BE1","BE4"),
    "BE1_vs_BE6" = c("BE1","BE6"),
    "BE2_vs_BE3" = c("BE2","BE3")
)
de_res <- lapply(pairs, function(p) pair_de(p[1], p[2]))
say("  %-12s  %6s  %8s  %8s", "pair", "nDEG", "up_left", "up_right")
for (nm in names(de_res)) {
    r <- de_res[[nm]]
    say("  %-12s  %6d  %8d  %8d", nm, r$n, r$up_a, r$up_b)
}
de_bar <- data.frame(pair = names(de_res),
                     nDEG = sapply(de_res, `[[`, "n"))
ggsave(file.path(OUT_DIR, "01_pairwise_nDEG.png"),
    ggplot(de_bar, aes(reorder(pair, nDEG), nDEG, fill = pair == "BE1_vs_BE2")) +
        geom_col() + coord_flip() +
        geom_text(aes(label = nDEG), hjust = -0.2, size = 4) +
        scale_fill_manual(values = c("grey60", "#D55E00"), guide = "none") +
        labs(x = NULL, y = "# DEG (|log2FC|>0.5, padj<0.05)",
             title = "직접 pairwise DEG 수",
             subtitle = "BE1-vs-BE2(주황)가 대조 쌍보다 훨씬 적으면 = 같은 실체") +
        expand_limits(y = max(de_bar$nDEG) * 1.1) + theme_bw(),
    width = 8, height = 5, bg = "white")

# =============================================================================
# 2. Nature of BE1-vs-BE2 DEGs: state (cell-cycle / ribosome-biogenesis) vs lineage
# =============================================================================
say("\n===== [2] BE1-vs-BE2 DEG 성격 분류 =====")
cc <- c(cc.genes.updated.2019$s.genes, cc.genes.updated.2019$g2m.genes)
ribo_pat <- "^(RPL|RPS|MRPL|MRPS|RRP|RRS|NOP|POLR1|UTP|BYSL|DDX2[0-9]|GNL[0-9]|WDR12|BOP1|PES1|FBL|NCL|NPM1|NAT10|LYAR|MYBBP1A|GAR1|NHP2|NOL[0-9])"
d12 <- de_res[["BE1_vs_BE2"]]$sig
be2_up <- subset(d12, avg_log2FC < 0)   # BE2 에서 상향
be1_up <- subset(d12, avg_log2FC > 0)   # BE1 에서 상향
classify <- function(g) {
    ifelse(g %in% cc, "cell_cycle",
    ifelse(grepl(ribo_pat, g), "ribosome_biogenesis", "other"))
}
if (nrow(be2_up) > 0) {
    cl2 <- table(classify(be2_up$gene))
    say("  BE2-up DEG (n=%d) 구성: %s", nrow(be2_up),
        paste(sprintf("%s=%d(%.0f%%)", names(cl2), cl2, 100*cl2/sum(cl2)), collapse=", "))
    say("  BE2-up top12: %s", paste(head(be2_up$gene[order(be2_up$avg_log2FC)], 12), collapse=", "))
} else say("  BE2-up DEG: 0")
if (nrow(be1_up) > 0) {
    say("  BE1-up top12: %s", paste(head(be1_up$gene[order(-be1_up$avg_log2FC)], 12), collapse=", "))
}
say("  최대 |log2FC| (BE1 vs BE2) = %.2f  vs  BE1-vs-BE3 = %.2f",
    if (nrow(d12)) max(abs(d12$avg_log2FC)) else 0,
    max(abs(de_res[["BE1_vs_BE3"]]$sig$avg_log2FC)))

# BE1 대표 마커를 BE2 도 공유하는가 (DotPlot across BE1-6)
be1_markers <- c("CCL2","TNC","CCN2","EDN1","DKK1","LAMC2","S100A2","PHLDA1","GLIPR1","NRG1")
be1_markers <- intersect(be1_markers, rownames(epi))
p_dot <- DotPlot(epi, features = be1_markers,
                 idents = c("BE1","BE2","BE3","BE4","BE5","BE6")) +
    RotatedAxis() +
    scale_color_viridis_c() +
    labs(title = "BE1 one-vs-rest 대표마커의 BE1–6 발현",
         subtitle = "BE2 가 BE1 프로그램을 공유(높음)하면 = 같은 계열")
ggsave(file.path(OUT_DIR, "02_BE1_markers_across_BE.png"),
       p_dot, width = 10, height = 5, bg = "white")

# =============================================================================
# 3. Centroid similarity (expression corr + harmony distance) ----
# =============================================================================
say("\n===== [3] Centroid 유사도: BE2 의 최근접 이웃 =====")
hvg <- head(VariableFeatures(epi), 2000)
hvg <- intersect(hvg, rownames(epi))
av <- AverageExpression(epi, assays = "SCT", features = hvg,
                        group.by = "be", verbose = FALSE)$SCT
av <- log1p(as.matrix(av))
cormat <- cor(av)
png(file.path(OUT_DIR, "03_centroid_expr_corr.png"), width = 1400, height = 1200, res = 150)
pheatmap(cormat, color = viridis(100), display_numbers = TRUE,
         number_format = "%.2f", main = "Cluster centroid 발현 상관 (top2000 HVG)")
dev.off()
nn_be2 <- sort(cormat["BE2", setdiff(colnames(cormat), "BE2")], decreasing = TRUE)
say("  발현상관 기준 BE2 최근접 3: %s",
    paste(sprintf("%s(%.3f)", names(nn_be2)[1:3], nn_be2[1:3]), collapse=", "))

# harmony 공간 centroid 거리
he <- Embeddings(epi, "harmony")[, 1:30]
cent <- t(sapply(levels(epi$be), function(g) colMeans(he[epi$be == g, , drop = FALSE])))
dmat <- as.matrix(dist(cent))
nn_be2_h <- sort(dmat["BE2", setdiff(rownames(dmat), "BE2")])
say("  harmony거리 기준 BE2 최근접 3: %s",
    paste(sprintf("%s(%.2f)", names(nn_be2_h)[1:3], nn_be2_h[1:3]), collapse=", "))

# =============================================================================
# 4. Resolution hierarchy: BE2 가 BE1 로 먼저 합쳐지는가 ----
# =============================================================================
say("\n===== [4] 해상도 계층 (coarser res 에서 cl0=BE1, cl7=BE2 의 소속) =====")
for (r in c("SCT_snn_res.0.1","SCT_snn_res.0.2","SCT_snn_res.0.3")) {
    if (!r %in% colnames(epi@meta.data)) next
    be1_par <- names(which.max(table(epi@meta.data[[r]][epi$SCT_snn_res.0.4 == "0"])))
    be2_par <- names(which.max(table(epi@meta.data[[r]][epi$SCT_snn_res.0.4 == "7"])))
    frac_be2 <- mean(epi@meta.data[[r]][epi$SCT_snn_res.0.4 == "7"] == be1_par)
    say("  %s: BE1→cluster %s, BE2→cluster %s  | BE2 세포의 %.0f%% 가 BE1 부모에 속함 %s",
        sub("SCT_snn_res.","res ",r), be1_par, be2_par, 100*frac_be2,
        if (be1_par == be2_par) "[MERGED]" else "[split]")
}

# =============================================================================
# 5. Seed co-assignment: 21 seed 에서 BE2→BE1 병합률 ----
# =============================================================================
say("\n===== [5] Seed co-assignment (21 random+ref seed) =====")
if (file.exists(LABELS_CSV)) {
    lab <- read.csv(LABELS_CSV, stringsAsFactors = FALSE, check.names = FALSE)
    seed_cols <- grep("^seed_", names(lab), value = TRUE)
    r04 <- setNames(as.character(epi$SCT_snn_res.0.4), colnames(epi))
    idx <- match(lab$cell, names(r04))
    cell_lab <- r04[idx]                      # canonical cluster id per row of lab
    # 각 seed 에서: query 클러스터 세포가 target 의 plurality seed-cluster 로 가는 비율
    merge_rate <- function(query_id, target_id) {
        mean(sapply(seed_cols, function(sc) {
            tg <- lab[[sc]][cell_lab == target_id]
            qy <- lab[[sc]][cell_lab == query_id]
            if (!length(tg) || !length(qy)) return(NA)
            cmaj <- names(which.max(table(tg)))
            mean(qy == cmaj)
        }), na.rm = TRUE)
    }
    r_be2 <- merge_rate("7","0")   # BE2 → BE1
    r_be3 <- merge_rate("10","0")  # BE3 → BE1 (대조)
    r_be4 <- merge_rate("1","0")   # BE4 → BE1 (대조)
    r_self<- merge_rate("0","0")   # BE1 자기응집 (상한)
    say("  BE1 자기응집(상한)          : %.2f", r_self)
    say("  BE2 세포가 BE1 majority-cluster로 가는 비율 : %.2f", r_be2)
    say("  대조 BE3→BE1               : %.2f", r_be3)
    say("  대조 BE4→BE1               : %.2f", r_be4)
} else say("  seed_labels.csv 없음 → skip")

# =============================================================================
# 6. Proliferation state check (BE2 ribosome-biogenesis 근거) ----
# =============================================================================
say("\n===== [6] Proliferation/state (BE1 vs BE2) =====")
ph <- prop.table(table(epi$be, epi$Phase), 1)
for (g in c("BE1","BE2")) {
    say("  %s Phase: %s | S=%.2f G2M=%.2f",
        g, paste(sprintf("%s %.0f%%", colnames(ph), 100*ph[g,]), collapse=" "),
        mean(epi$S.Score[epi$be==g]), mean(epi$G2M.Score[epi$be==g]))
}
prolif <- intersect(c("MKI67","TOP2A","UBE2C","PCNA"), rownames(epi))
p_vln <- VlnPlot(epi, features = prolif, idents = c("BE1","BE2","BE3","BE4","BE5","BE6"),
                 pt.size = 0, cols = utils_cb_palette(6), stack = FALSE, ncol = 2)
ggsave(file.path(OUT_DIR, "06_proliferation_vln.png"),
       p_vln, width = 11, height = 8, bg = "white")

say("\n################ END ################")
close(sink_con)
message("Report → ", file.path(OUT_DIR, "REPORT.txt"))
