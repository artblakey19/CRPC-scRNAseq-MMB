#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 03 : 환자별 inferCNV 실행 (HMM i6, subclusters)
# -----------------------------------------------------------------------------
# 입력: 08 raw count_mat.rds + 02 annotations.txt + 01 gene_order.
#   reference = 비-상피 7종 + Ionocyte anchor / observations = 상피 fine cluster.
#   cutoff=0.1 (10x droplet), denoise, HMM i6, analysis_mode=subclusters,
#   BayesMaxPNormal=0.5 (JAGS 사후필터로 정상세포 오탐 억제).
# 실행: Rscript scripts/09_inferCNV/03_run_infercnv.R   (프로젝트 루트에서)
# =============================================================================
suppressPackageStartupMessages({ library(infercnv) })
options(scipen = 100)                            # analysis_mode="subclusters" hclust 에러 예방(inferCNV 권고)

N_CORES <- 30                                   # 하드코딩 (프로젝트 규약; detectCores 금지)
NB  <- "Results/08_Numbat"
OUT <- "Results/09_inferCNV"
GO  <- file.path(OUT, "gene_order_GRCh38.txt")
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")

for (S in SAMPLES) {
    cat(sprintf("\n==================== %s ====================\n", S))
    cm  <- readRDS(file.path(NB, S, "count_mat.rds"))
    ref <- read.csv(file.path(NB, S, "cell_annot_ref.csv"), stringsAsFactors = FALSE)
    ref_groups <- sort(unique(ref$group))
    cat("reference groups:", paste(ref_groups, collapse = ", "), "\n")

    obj <- CreateInfercnvObject(
        raw_counts_matrix = cm,
        annotations_file  = file.path(OUT, S, "annotations.txt"),
        gene_order_file   = GO,
        ref_group_names   = ref_groups)

    infercnv::run(
        obj,
        cutoff            = 0.1,                 # 10x droplet 권장값
        out_dir           = file.path(OUT, S, "infercnv"),
        cluster_by_groups = TRUE,                # 상피 cluster 별 관측 그룹화
        denoise           = TRUE,
        HMM               = TRUE,
        HMM_type          = "i6",
        analysis_mode     = "subclusters",       # cluster 내 subclone CNV state 콜링
        BayesMaxPNormal   = 0.5,                 # JAGS 사후필터
        num_threads       = N_CORES,
        resume_mode       = TRUE,                # 재실행 시 완료단계 skip
        no_prelim_plot    = TRUE,
        output_format     = "pdf")

    cat(sprintf("[03] %s done -> %s\n", S, file.path(OUT, S, "infercnv")))
}
cat("\n[03] all samples done.\n")
