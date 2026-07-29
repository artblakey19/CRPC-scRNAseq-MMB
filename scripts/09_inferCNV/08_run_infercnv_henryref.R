#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 08 : Henry 외부 정상 상피 reference 로 inferCNV 재실행
# -----------------------------------------------------------------------------
# 09(01~05) 의 confound = reference 가 비-상피(stroma/immune) → 상피발현 CNV leak.
# 여기서는 필드 표준(Kfoury 2021·Dong 2020)대로 **외부 정상 전립선 상피(Henry 2018,
#   07 추출: basal/luminal/secretory)**를 reference 로 써서 상피-vs-상피 비교로 만든다.
#   observations = 우리 상피 fine cluster(11개; Ionocyte 제외, 09 와 동일).
#   reference     = Henry_basal / Henry_luminal / Henry_secretory (3 groups).
# ⚠ 외부 데이터라 study/chemistry batch 는 남는 위험(결과는 비판적으로 평가).
# 실행: Rscript scripts/09_inferCNV/08_run_infercnv_henryref.R
# =============================================================================
suppressPackageStartupMessages({ library(infercnv); library(Matrix) })
options(scipen = 100)
N_CORES <- 30

NB   <- "Results/08_Numbat"
BASE <- "Results/09_inferCNV"
OUT  <- file.path(BASE, "henryref")
GO   <- file.path(BASE, "gene_order_GRCh38.txt")
REFD <- "Raw_data/Henry_GSE120716/ref_matrix"
# 환자 인자 지정 시 해당 환자만(동시 병렬 실행용); 없으면 3환자 순차
.args <- commandArgs(trailingOnly = TRUE)
SAMPLES <- if (length(.args) > 0) .args else c("CRPC1", "CRPC2", "CRPC3")
EPI <- c("LE(ARPC)","Club","Hillock 1","Hillock 2",
         "BE 1","BE 2","BE 3","BE 4","BE 5","BE 6","OE")     # observations (09 와 동일)

# --- Henry reference matrix (한 번만 로드) ----------------------------------
hf   <- read.delim(file.path(REFD, "features.tsv.gz"), header = FALSE)   # ensembl, symbol
hbc  <- readLines(gzfile(file.path(REFD, "barcodes.tsv.gz")))
hmat <- as(Matrix::readMM(gzfile(file.path(REFD, "matrix.mtx.gz"))), "CsparseMatrix")
rownames(hmat) <- make.unique(hf$V2)          # 심볼(우리 count_mat 명명과 동일 규약)
colnames(hmat) <- hbc
hgrp <- read.csv(file.path(REFD, "cell_groups.csv"), stringsAsFactors = FALSE)  # cell, group
ref_groups <- sort(unique(hgrp$group))
cat(sprintf("Henry reference: %d genes x %d cells | groups: %s\n",
            nrow(hmat), ncol(hmat), paste(ref_groups, collapse = ", ")))

for (S in SAMPLES) {
    cat(sprintf("\n==================== %s ====================\n", S))
    cm  <- readRDS(file.path(NB, S, "count_mat.rds"))
    ann <- read.delim(file.path(BASE, S, "annotations.txt"), header = FALSE,
                      stringsAsFactors = FALSE)                 # V1=cell, V2=group
    obs <- ann[ann$V2 %in% EPI, ]                              # 우리 상피 관측세포
    our <- cm[, obs$V1, drop = FALSE]

    common <- intersect(rownames(our), rownames(hmat))
    cat(sprintf("obs %d cells | common genes (ours∩Henry) = %d\n", ncol(our), length(common)))

    comb <- cbind(our[common, , drop = FALSE], hmat[common, , drop = FALSE])
    annot <- rbind(data.frame(cell = obs$V1,  group = obs$V2),
                   data.frame(cell = hgrp$cell, group = hgrp$group))
    stopifnot(identical(annot$cell, colnames(comb)))

    dir.create(file.path(OUT, S), recursive = TRUE, showWarnings = FALSE)
    af <- file.path(OUT, S, "annotations_henryref.txt")
    write.table(annot, af, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

    obj <- CreateInfercnvObject(raw_counts_matrix = comb,
                                annotations_file = af,
                                gene_order_file  = GO,
                                ref_group_names  = ref_groups)
    infercnv::run(obj,
        cutoff = 0.1, out_dir = file.path(OUT, S, "infercnv"),
        cluster_by_groups = TRUE, denoise = TRUE,
        HMM = TRUE, HMM_type = "i6", analysis_mode = "subclusters",
        BayesMaxPNormal = 0.5, num_threads = N_CORES,
        resume_mode = TRUE, no_prelim_plot = TRUE, output_format = "pdf")
    cat(sprintf("[08] %s done -> %s\n", S, file.path(OUT, S, "infercnv")))
}
cat("\n[08] all samples done.\n")
