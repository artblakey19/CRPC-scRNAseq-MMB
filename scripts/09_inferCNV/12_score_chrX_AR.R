#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 12 : chrX/AR 신호 판독 (11 결과)
# -----------------------------------------------------------------------------
# 11(chrX 포함) 결과에서 AR amplification 흔적을 본다.
#   - chrX heatmap 렌더 (관측 vs reference; AR = Xq12 ~67.5Mb)
#   - per-cell denoise 잔차발현 평균: 전 chrX / AR 영역(Xq12, chrX 60–75Mb) / AR 단독
#   - 상피 cluster × 환자, 그리고 Numbat tumor/normal 별로 요약 → 환자 간 교차비교
#     (lineage AR 발현은 환자 간 비슷 → 암많은 CRPC2/CRPC3 tumor 초과분 = amp 시사)
# =============================================================================
suppressPackageStartupMessages({ library(infercnv) })

BASE <- "Results/09_inferCNV"
OUT  <- file.path(BASE, "withX")
VIS  <- file.path(OUT, "Visualization")
NB   <- "Results/08_Numbat"
dir.create(VIS, recursive = TRUE, showWarnings = FALSE)
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
REF_GRP <- c("Endothelial","Fibroblast","Ionocyte","Mast cells",
             "Phagocytes","Smooth muscle cells","T/NK cells")

pick_final <- function(S) {
    fs <- list.files(file.path(NB, S, "numbat"), pattern = "^clone_post_[0-9]+\\.tsv$", full.names = TRUE)
    fs[which.max(as.integer(sub(".*clone_post_([0-9]+)\\.tsv$", "\\1", fs)))]
}

rows <- list()
for (S in SAMPLES) {
    f <- file.path(OUT, S, "infercnv", "run.final.infercnv_obj")
    if (!file.exists(f)) { cat("skip (no obj):", S, "\n"); next }
    o  <- readRDS(f)

    # chrX heatmap 렌더
    infercnv::plot_cnv(o, out_dir = VIS, output_filename = paste0("heatmap_", S, "_withX"),
                       output_format = "png", png_res = 140, cluster_by_groups = TRUE,
                       x.center = 1, title = paste(S, "(chrX 포함)"),
                       color_safe_pal = FALSE, write_expr_matrix = FALSE)

    gu <- o@gene_order
    E  <- o@expr.data
    chrX  <- rownames(gu)[gu$chr == "chrX"]
    arreg <- rownames(gu)[gu$chr == "chrX" & gu$start >= 60e6 & gu$stop <= 75e6]  # Xq12 근방
    chrX  <- intersect(chrX, rownames(E)); arreg <- intersect(arreg, rownames(E))

    d <- data.frame(cell = colnames(E),
                    chrX_sig = colMeans(E[chrX, , drop = FALSE]),
                    ARreg_sig = colMeans(E[arreg, , drop = FALSE]),
                    AR_sig = if ("AR" %in% rownames(E)) E["AR", ] else NA_real_,
                    row.names = NULL, stringsAsFactors = FALSE)
    # annotation (obs) + reference 라벨
    ann <- read.delim(file.path(BASE, S, "annotations.txt"), header = FALSE, stringsAsFactors = FALSE)
    d$group <- ann$V2[match(d$cell, ann$V1)]
    # Numbat compartment (tumor/normal)
    cp <- read.delim(pick_final(S), stringsAsFactors = FALSE)
    d$numbat <- cp$compartment_opt[match(d$cell, cp$cell)]
    d$sample <- S
    d$is_ref <- d$group %in% REF_GRP
    rows[[S]] <- d
    cat(S, sprintf(": chrX genes %d | AR-region genes %d | rendered\n", length(chrX), length(arreg)))
}
df <- do.call(rbind, rows)
write.csv(df, file.path(OUT, "chrX_AR_percell.csv"), row.names = FALSE)

cat("\n=== [환자별] reference(정상) 대비 상피 tumor/normal 의 AR 영역(Xq12) 신호 ===\n")
cat("   (denoise 잔차; 1=중립, >1=gain. reference=비상피 정상, 환자내 baseline)\n")
agg <- function(sub) round(c(chrX = mean(sub$chrX_sig), ARregion = mean(sub$ARreg_sig),
                             AR = mean(sub$AR_sig, na.rm = TRUE)), 4)
for (S in SAMPLES) {
    d <- df[df$sample == S, ]
    cat(sprintf("\n%s:\n", S))
    cat("  reference(비상피)      "); print(agg(d[d$is_ref, ]))
    cat("  epi Numbat-normal      "); print(agg(d[!d$is_ref & d$numbat == "normal", ]))
    cat("  epi Numbat-tumor       "); print(agg(d[!d$is_ref & d$numbat == "tumor", ]))
}

cat("\n=== [상피 cluster × 환자] AR 영역(Xq12) 신호 ===\n")
ep <- df[!df$is_ref, ]
tab <- aggregate(ARreg_sig ~ sample + group, ep, mean)
tab$ARreg_sig <- round(tab$ARreg_sig, 4)
tab <- tab[order(tab$sample, -tab$ARreg_sig), ]
print(tab, row.names = FALSE)
write.csv(tab, file.path(OUT, "chrX_AR_by_cluster.csv"), row.names = FALSE)
cat("\n[12] done ->", OUT, "\n")
