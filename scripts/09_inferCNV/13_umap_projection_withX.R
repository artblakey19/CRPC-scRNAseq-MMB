#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 13 : withX(chrX 포함) 결과를 상피 UMAP 에 투영
# -----------------------------------------------------------------------------
# 06(stroma-ref)/10(henryref)과 동일 형식 + withX 전용 AR 패널.
#   - inferCNV_PGA_withX   : subcluster PGA (autosome chr1-22; 이전과 비교 위해 고정)
#   - inferCNV_score_withX : mean((denoise expr - 1)^2)  (chrX 포함)
#   - inferCNV_AR_withX    : AR 영역(Xq12) denoise 잔차 신호 (12 chrX_AR_percell.csv)
# observations(우리 상피)만 매핑. 연속형 → viridis.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(infercnv); library(ggplot2); library(patchwork) })

BASE <- "Results/09_inferCNV"
OUT  <- file.path(BASE, "withX")
VIS  <- file.path(OUT, "Visualization")
GO   <- file.path(BASE, "gene_order_GRCh38.txt")
dir.create(VIS, showWarnings = FALSE, recursive = TRUE)
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
LOSS <- c(1,2); GAIN <- c(4,5,6)
cb <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
        "#999999","#000000","#332288","#117733","#88CCEE","#AA4499")

# 분모 D: chr1-22 span (05 와 동일)
go <- read.table(GO, sep="\t", header=FALSE, col.names=c("gene","chr","start","stop"))
auto <- paste0("chr", 1:22)
D_mb <- sum(sapply(auto, function(c){g<-go[go$chr==c,]; max(g$stop)-min(g$start)}))/1e6

epi <- readRDS("Results/05_Epithelial_Downstream/epi_annotated.rds")
red <- if ("umap" %in% Reductions(epi)) "umap" else Reductions(epi)[1]
epi_key <- paste(as.character(epi$orig.ident), sub("^P[0-9]+_", "", colnames(epi)), sep = "|")

cellpga <- setNames(numeric(0), character(0)); cellscore <- setNames(numeric(0), character(0))
for (S in SAMPLES) {
    o  <- readRDS(file.path(OUT, S, "infercnv", "run.final.infercnv_obj"))
    cn <- colnames(o@expr.data)
    cellscore[paste(S, cn, sep="|")] <- colMeans((o@expr.data - 1)^2)   # mean-sq (chrX 포함)

    # subcluster PGA (autosome only)
    sc <- o@tumor_subclusters$subclusters
    rf <- file.path(OUT, S, "infercnv",
                    "HMM_CNV_predictions.HMMi6.leiden.hmm_mode-subclusters.Pnorm_0.5.pred_cnv_regions.dat")
    r <- read.delim(rf, stringsAsFactors=FALSE); r <- r[r$chr %in% auto, ]
    r$subcl <- sub("^[^.]*\\.", "", r$cell_group_name); r$len <- (r$end - r$start)/1e6
    nn <- tapply(r$len[r$state %in% c(LOSS,GAIN)], r$subcl[r$state %in% c(LOSS,GAIN)], sum)
    pga_sc <- setNames(as.numeric(nn)/D_mb, names(nn))
    for (g in names(sc)) for (s in names(sc[[g]])) {
        v <- pga_sc[[s]]; bc <- cn[sc[[g]][[s]]]
        cellpga[paste(S, bc, sep="|")] <- if (is.null(v)) 0 else v
    }
    rm(o); gc()
}

# AR 영역 (12 결과)
ar <- read.csv(file.path(OUT, "chrX_AR_percell.csv"), stringsAsFactors=FALSE)
arv <- setNames(ar$ARreg_sig, paste(ar$sample, ar$cell, sep="|"))

epi$inferCNV_PGA_withX   <- unname(cellpga[epi_key])
epi$inferCNV_score_withX <- unname(cellscore[epi_key])
epi$inferCNV_AR_withX    <- unname(arv[epi_key])
cat(sprintf("mapped PGA %d | score %d | AR %d / %d cells\n",
            sum(!is.na(epi$inferCNV_PGA_withX)), sum(!is.na(epi$inferCNV_score_withX)),
            sum(!is.na(epi$inferCNV_AR_withX)), ncol(epi)))

cap <- as.numeric(quantile(epi$inferCNV_PGA_withX, 0.99, na.rm=TRUE))
epi$PGA_withX_capped <- pmin(epi$inferCNV_PGA_withX, cap)
vir <- function(nm) scale_color_viridis_c(option="viridis", na.value="grey88", name=nm)

# 1) overview
p_ann <- DimPlot(epi, reduction=red, group.by="annotation", cols=cb, label=TRUE, repel=TRUE, label.size=3) +
         ggtitle("Epithelial clusters") + theme(legend.position="none")
p_pga <- FeaturePlot(epi, "PGA_withX_capped", reduction=red, order=TRUE, pt.size=0.2) + vir("PGA") + ggtitle("withX inferCNV PGA (autosome)")
p_sco <- FeaturePlot(epi, "inferCNV_score_withX", reduction=red, order=TRUE, pt.size=0.2) + vir("score") + ggtitle("withX mean-sq score")
ggsave(file.path(VIS, "umap_withX_overview.png"), p_ann + p_pga + p_sco + plot_layout(nrow=1), width=18, height=5.5, dpi=150)

# 2) PGA by patient
ggsave(file.path(VIS, "umap_withX_pga_by_patient.png"),
       FeaturePlot(epi, "PGA_withX_capped", reduction=red, split.by="orig.ident", order=TRUE, pt.size=0.25) & vir("PGA"),
       width=16, height=5.2, dpi=150)

# 3) AR 영역(Xq12) by patient — withX 전용. AR amp 있으면 tumor 영역이 밝아야.
arcap <- quantile(epi$inferCNV_AR_withX, c(0.01,0.99), na.rm=TRUE)
epi$AR_withX_clip <- pmin(pmax(epi$inferCNV_AR_withX, arcap[1]), arcap[2])
ggsave(file.path(VIS, "umap_withX_AR_by_patient.png"),
       FeaturePlot(epi, "AR_withX_clip", reduction=red, split.by="orig.ident", order=TRUE, pt.size=0.25) &
         scale_color_viridis_c(option="viridis", na.value="grey88", name="AR region\n(Xq12)"),
       width=16, height=5.2, dpi=150)
cat("[13] done ->", VIS, "\n")
