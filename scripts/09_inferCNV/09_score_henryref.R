#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 09 : Henry-ref 결과 채점 + 시각화 + 이전 결과 비교
# -----------------------------------------------------------------------------
# 08(Henry 외부 정상 상피 reference) 결과에 대해:
#   (a) 환자별 heatmap PNG 렌더 (관측 vs Henry reference)
#   (b) subcluster PGA(HMM Bayesian Pnorm0.5) → 상피 cluster/환자별 요약 (05 와 동일식)
#   (c) Henry reference 그룹 PGA = 노이즈 바닥 → threshold
#   (d) 이전 stroma/immune-ref PGA(05: pga_by_cluster.csv) 와 나란히 비교
# =============================================================================
suppressPackageStartupMessages({ library(infercnv) })

BASE <- "Results/09_inferCNV"
OUT  <- file.path(BASE, "henryref")
VIS  <- file.path(OUT, "Visualization")
dir.create(VIS, recursive = TRUE, showWarnings = FALSE)
GO   <- file.path(BASE, "gene_order_GRCh38.txt")
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
REF_GRP <- c("Henry_basal", "Henry_luminal", "Henry_secretory")
LOSS <- c(1,2); GAIN <- c(4,5,6)

# 분모 D = chr1-22 평가 span (05 와 동일, 비교 위해 고정)
go <- read.table(GO, sep="\t", header=FALSE, col.names=c("gene","chr","start","stop"))
auto <- paste0("chr", 1:22)
D_mb <- sum(sapply(auto, function(c){g<-go[go$chr==c,]; max(g$stop)-min(g$start)}))/1e6

per_cell <- list(); grp_all <- list()
for (S in SAMPLES) {
    f <- file.path(OUT, S, "infercnv", "run.final.infercnv_obj")
    if (!file.exists(f)) { cat("skip (no obj):", S, "\n"); next }
    o  <- readRDS(f)

    # (a) heatmap 렌더
    infercnv::plot_cnv(o, out_dir = VIS, output_filename = paste0("heatmap_", S, "_henryref"),
                       output_format = "png", png_res = 140, cluster_by_groups = TRUE,
                       x.center = 1, title = paste(S, "(Henry ref)"),
                       color_safe_pal = FALSE, write_expr_matrix = FALSE)

    # (b) subcluster PGA
    sc <- o@tumor_subclusters$subclusters
    scdf <- do.call(rbind, lapply(names(sc), function(g)
        data.frame(group=g, subcl=names(sc[[g]]), n=sapply(sc[[g]], length), stringsAsFactors=FALSE)))
    rf <- file.path(OUT, S, "infercnv",
                    "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
    r <- read.delim(rf, stringsAsFactors=FALSE); r <- r[r$chr %in% auto, ]
    r$subcl <- sub("^[^.]*\\.", "", r$cell_group_name); r$len <- (r$end - r$start)/1e6
    mb <- function(st){x<-tapply(r$len[r$state%in%st], r$subcl[r$state%in%st], sum); setNames(as.numeric(x),names(x))}
    lo<-mb(LOSS); ga<-mb(GAIN)
    scdf$loss_mb <- ifelse(is.na(lo[scdf$subcl]),0,lo[scdf$subcl])
    scdf$gain_mb <- ifelse(is.na(ga[scdf$subcl]),0,ga[scdf$subcl])
    scdf$pga <- (scdf$loss_mb+scdf$gain_mb)/D_mb; scdf$sample <- S
    grp_all[[S]] <- scdf
    per_cell[[S]] <- scdf[rep(seq_len(nrow(scdf)), scdf$n), c("sample","group","pga")]
    cat(S, "rendered + scored\n")
}
grp  <- do.call(rbind, grp_all)
cell <- do.call(rbind, per_cell)
write.csv(grp, file.path(OUT, "pga_by_subcluster_henryref.csv"), row.names=FALSE)

# (c) Henry reference threshold
ref_pga <- cell$pga[cell$group %in% REF_GRP]
thr <- as.numeric(quantile(ref_pga, 0.99))
cat(sprintf("\n[Henry ref] PGA median=%.3f 99pct(thr)=%.3f max=%.3f\n",
            median(ref_pga), thr, max(ref_pga)))

cell$is_epi <- !cell$group %in% REF_GRP
summ <- aggregate(cbind(pga, aneuploid=as.integer(pga>thr)) ~ sample+group, cell, mean)
summ$n <- aggregate(pga~sample+group, cell, length)$pga
summ$type <- ifelse(summ$group %in% REF_GRP, "reference", "epithelial")
summ <- summ[order(summ$type, summ$sample, -summ$pga), ]
write.csv(summ, file.path(OUT, "pga_by_cluster_henryref.csv"), row.names=FALSE)
cat("\n=== 상피 cluster PGA / aneuploid frac (Henry ref) ===\n")
ep <- summ[summ$type=="epithelial",]; ep$pga<-round(ep$pga,3); ep$aneuploid<-round(ep$aneuploid,3)
print(ep[,c("sample","group","n","pga","aneuploid")], row.names=FALSE)
cat("\n=== 환자별 상피 aneuploid fraction ===\n")
for (S in SAMPLES){ e<-cell[cell$sample==S & cell$is_epi,]
    cat(sprintf("  %s: %.1f%% aneuploid | mean PGA %.3f\n", S, 100*mean(e$pga>thr), mean(e$pga))) }

# (d) 이전 stroma/immune-ref 비교
old <- file.path(BASE, "pga_by_cluster.csv")
if (file.exists(old)) {
    o5 <- read.csv(old, stringsAsFactors=FALSE)
    o5 <- o5[o5$type=="epithelial", c("sample","group","pga")]; names(o5)[3]<-"pga_stroma_ref"
    nw <- ep[,c("sample","group","pga")]; names(nw)[3]<-"pga_henry_ref"
    cmp <- merge(nw, o5, by=c("sample","group"))
    cmp$delta <- round(cmp$pga_henry_ref - cmp$pga_stroma_ref, 3)
    cmp <- cmp[order(cmp$sample, -cmp$pga_henry_ref),]
    write.csv(cmp, file.path(OUT, "pga_compare_stroma_vs_henry.csv"), row.names=FALSE)
    cat("\n=== stroma-ref vs Henry-ref PGA 비교 ===\n"); print(cmp, row.names=FALSE)
}
cat("\n[09] done ->", OUT, "\n")
