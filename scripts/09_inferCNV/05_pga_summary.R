#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 05 : subcluster별 유전체 변이 비율(PGA) + aneuploid fraction
# -----------------------------------------------------------------------------
# mean-squared score(04)는 전유전체 평균이라 국소 arm-level CNV 를 희석한다.
# 제대로 된 malignancy 지표 = HMM(i6, Bayesian Pnorm_0.5)이 콜한 CNV region 기반
# Proportion of Genome Altered(PGA).
#   - state 3 = 중립(Bayesian 필터로 정상 복귀), 1·2 = loss, 4·5·6 = gain
#   - 분모 D = 평가된 상염색체(chr1–22) 길이(gene_order 기반). chrX/Y 는 run 기본 제외.
#   - subcluster별 PGA 를 세포 수로 가중 → 상피 cluster/환자별 요약.
#   - reference 그룹 subcluster PGA = 노이즈 바닥 → aneuploid threshold 로 사용.
# =============================================================================
suppressPackageStartupMessages({ library(infercnv) })

OUT <- "Results/09_inferCNV"
GO  <- file.path(OUT, "gene_order_GRCh38.txt")
SAMPLES  <- c("CRPC1", "CRPC2", "CRPC3")
REF_GRP  <- c("Endothelial","Fibroblast","Ionocyte","Mast cells",
              "Phagocytes","Smooth muscle cells","T/NK cells")
LOSS <- c(1,2); GAIN <- c(4,5,6)                 # i6 HMM states (3 = neutral)

# 분모 D : chr1–22 평가 길이(Mb) ---------------------------------------------
go <- read.table(GO, sep = "\t", header = FALSE, stringsAsFactors = FALSE,
                 col.names = c("gene","chr","start","stop"))
auto <- paste0("chr", 1:22)
span <- sapply(auto, function(c) { g <- go[go$chr == c, ]; max(g$stop) - min(g$start) })
D_mb <- sum(span) / 1e6
cat(sprintf("[05] assessed autosomal span D = %.0f Mb (chr1-22)\n", D_mb))

per_cell <- list(); per_grp <- list()
for (S in SAMPLES) {
    o  <- readRDS(file.path(OUT, S, "infercnv", "run.final.infercnv_obj"))
    sc <- o@tumor_subclusters$subclusters       # group -> (subcluster -> cell idx)

    # subcluster -> n_cells, group
    scdf <- do.call(rbind, lapply(names(sc), function(g)
        data.frame(group = g, subcl = names(sc[[g]]),
                   n = sapply(sc[[g]], length), stringsAsFactors = FALSE)))

    # regions (Bayesian Pnorm_0.5) -> subcluster별 loss/gain Mb
    rf <- file.path(OUT, S, "infercnv",
                    "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
    r  <- read.delim(rf, stringsAsFactors = FALSE)          # cell_group_name,cnv_name,state,chr,start,end
    r  <- r[r$chr %in% auto, ]
    r$subcl <- sub("^[^.]*\\.", "", r$cell_group_name)      # 첫 "." 이후 = subcluster 키
    r$len_mb <- (r$end - r$start) / 1e6
    agg_mb <- function(states) {
        x <- tapply(r$len_mb[r$state %in% states], r$subcl[r$state %in% states], sum)
        setNames(as.numeric(x), names(x))
    }
    loss_mb <- agg_mb(LOSS); gain_mb <- agg_mb(GAIN)
    scdf$loss_mb <- ifelse(is.na(loss_mb[scdf$subcl]), 0, loss_mb[scdf$subcl])
    scdf$gain_mb <- ifelse(is.na(gain_mb[scdf$subcl]), 0, gain_mb[scdf$subcl])
    scdf$pga <- (scdf$loss_mb + scdf$gain_mb) / D_mb
    scdf$sample <- S
    per_grp[[S]] <- scdf

    # 세포 단위로 확장(각 세포 = 소속 subcluster 의 PGA)
    pc <- scdf[rep(seq_len(nrow(scdf)), scdf$n),
               c("sample","group","subcl","pga","loss_mb","gain_mb")]
    per_cell[[S]] <- pc
}
grp <- do.call(rbind, per_grp)
cell <- do.call(rbind, per_cell)
write.csv(grp,  file.path(OUT, "pga_by_subcluster.csv"), row.names = FALSE)

# reference 노이즈 바닥 → aneuploid threshold ---------------------------------
ref_pga <- cell$pga[cell$group %in% REF_GRP]
thr <- as.numeric(quantile(ref_pga, 0.99))       # 정상세포 99분위
cat(sprintf("[05] reference PGA: median=%.3f  99pct(threshold)=%.3f  max=%.3f\n",
            median(ref_pga), thr, max(ref_pga)))

# 상피 cluster × 환자 요약 ----------------------------------------------------
cell$is_epi <- !cell$group %in% REF_GRP
summ <- aggregate(cbind(pga, loss_mb, gain_mb, aneuploid = as.integer(pga > thr)) ~ sample + group,
                  cell, mean)
summ$n_cells <- aggregate(pga ~ sample + group, cell, length)$pga
summ$type <- ifelse(summ$group %in% REF_GRP, "reference", "epithelial")
summ <- summ[order(summ$type, summ$sample, -summ$pga), ]
names(summ)[names(summ)=="aneuploid"] <- "aneuploid_frac"
write.csv(summ, file.path(OUT, "pga_by_cluster.csv"), row.names = FALSE)

cat("\n=== 상피 cluster PGA / aneuploid fraction (환자별) ===\n")
ep <- summ[summ$type == "epithelial", c("sample","group","n_cells","pga","loss_mb","gain_mb","aneuploid_frac")]
ep[,c("pga","aneuploid_frac")] <- round(ep[,c("pga","aneuploid_frac")], 3)
ep[,c("loss_mb","gain_mb")] <- round(ep[,c("loss_mb","gain_mb")], 0)
print(ep, row.names = FALSE)

cat("\n=== 환자별 상피 전체 aneuploid fraction (Numbat 비교용) ===\n")
for (S in SAMPLES) {
    e <- cell[cell$sample == S & cell$is_epi, ]
    cat(sprintf("  %s: aneuploid %.1f%% (%d/%d epi cells > thr %.3f) | mean PGA %.3f\n",
                S, 100*mean(e$pga > thr), sum(e$pga > thr), nrow(e), thr, mean(e$pga)))
}
cat("\n[05] done -> pga_by_cluster.csv / pga_by_subcluster.csv\n")
