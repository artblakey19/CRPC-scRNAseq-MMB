#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 11 : chrX 포함 재실행 (AR amp 확인용)
# -----------------------------------------------------------------------------
# 03/08 은 inferCNV 기본 chr_exclude=chrX,chrY,chrM 라 AR(Xq12)이 안 보였다.
# CRPC 대표 CNV 인 AR amplification 을 놓쳤는지 확인하려고 chrX 를 포함해 재실행.
#   reference = stroma/immune + Ionocyte (02 annotations, 환자 내 → batch 없음)
#   gene_order(01)는 이미 chrX 1208개 포함(AR chrX:67.5Mb).
#   chr_exclude = c("chrY","chrM")  → chrX 만 추가로 포함.
# ⚠ AR 은 luminal lineage 유전자라 stroma-ref 대비 상피 전체가 높게 보임 → AR amp 는
#   "환자 간 tumor 세포 비교"로 해석(lineage 발현은 환자 간 비슷 → 암많은 환자의 초과분=amp).
# ⚠ 남성 chrX 는 hemizygous → allele(Numbat) 불가, expression(inferCNV)만 가능.
# 실행: Rscript scripts/09_inferCNV/11_run_infercnv_withX.R [SAMPLE]
# =============================================================================
suppressPackageStartupMessages({ library(infercnv) })
options(scipen = 100)
N_CORES <- 30

NB   <- "Results/08_Numbat"
BASE <- "Results/09_inferCNV"
OUT  <- file.path(BASE, "withX")
GO   <- file.path(BASE, "gene_order_GRCh38.txt")
.args <- commandArgs(trailingOnly = TRUE)
SAMPLES <- if (length(.args) > 0) .args else c("CRPC1", "CRPC2", "CRPC3")

for (S in SAMPLES) {
    cat(sprintf("\n==================== %s (chrX 포함) ====================\n", S))
    dir.create(file.path(OUT, S, "infercnv"), recursive = TRUE, showWarnings = FALSE)
    cm  <- readRDS(file.path(NB, S, "count_mat.rds"))
    ref <- read.csv(file.path(NB, S, "cell_annot_ref.csv"), stringsAsFactors = FALSE)
    ref_groups <- sort(unique(ref$group))

    obj <- CreateInfercnvObject(
        raw_counts_matrix = cm,
        annotations_file  = file.path(BASE, S, "annotations.txt"),
        gene_order_file   = GO,
        ref_group_names   = ref_groups,
        chr_exclude       = c("chrY", "chrM"))    # chrX 포함(기본은 chrX,chrY,chrM 제외)

    infercnv::run(obj,
        cutoff = 0.1, out_dir = file.path(OUT, S, "infercnv"),
        cluster_by_groups = TRUE, denoise = TRUE,
        HMM = TRUE, HMM_type = "i6", analysis_mode = "subclusters",
        BayesMaxPNormal = 0.5, num_threads = N_CORES,
        resume_mode = TRUE, no_prelim_plot = TRUE, output_format = "pdf")
    cat(sprintf("[11] %s done -> %s\n", S, file.path(OUT, S, "infercnv")))
}
cat("\n[11] done.\n")
