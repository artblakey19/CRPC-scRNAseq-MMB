#!/usr/bin/env Rscript
# =============================================================================
# Numbat step 2b  —  CONTAINER / pkharchenkolab/numbat-rbase  (numbat side)
# -----------------------------------------------------------------------------
# 02a 가 내보낸 환자별 count_mat.rds + cell_annot_ref.csv + step1 allele tsv 로
# normal 레퍼런스를 만들고 run_numbat 을 환자별로 실행한다.
# 컨테이너에서 cwd = /mnt/project 로 실행되므로 상대경로 사용.
# 실행은 02_run_numbat.sh 가 docker run 으로 진행.
# =============================================================================
suppressPackageStartupMessages({ library(numbat); library(Matrix); library(data.table) })

N_CORES <- 30
INIT_K  <- 12         # 기본 3. epi-only query 를 이 수로 초기 hclust 분할 → subclonal 집단을
                      #   독립 pseudobulk 로 분리, 희석으로 묻힌 CNV 검출 (numbat 권장 레버).
                      #   ⚠️ 반드시 run_numbat(init_k=INIT_K) 로 전달할 것 (과거 호출부 누락 → default 3).
NB_DIR  <- "Results/08_Numbat"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")

for (S in SAMPLES) {
    message("\n==================  ", S, "  ==================")
    s_dir   <- file.path(NB_DIR, S)
    out_dir <- file.path(s_dir, "numbat")
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    count_mat <- readRDS(file.path(s_dir, "count_mat.rds"))                 # gene x cell 정수 (전 celltype)
    ref_annot <- read.csv(file.path(s_dir, "cell_annot_ref.csv"),          # non-epi + Ionocyte
                          stringsAsFactors = FALSE)

    # normal 레퍼런스 (gene x cell type, normalized) — non-epi + Ionocyte(정상 상피)로 build
    ref_internal <- aggregate_counts(count_mat[, ref_annot$cell, drop = FALSE], ref_annot)

    # query = Epithelial 세포 (단 Ionocyte 제외). count_mat 은 전 세포, ref_annot 은
    #   non-epi + Ionocyte 이므로 여집합 = Epithelial \ Ionocyte. 레퍼런스 세포는 query 에서
    #   제외하고 normal anchor 로만 사용 (epi-only CNV inference).
    epi_cells <- setdiff(colnames(count_mat), ref_annot$cell)
    stopifnot(length(epi_cells) > 0)
    query_mat <- count_mat[, epi_cells, drop = FALSE]
    message(S, ": query ", length(epi_cells), " epithelial cells | ref ",
            length(unique(ref_annot$group)), " reference types (non-epi + Ionocyte)")

    # step1 산출물(allele counts)
    df_allele <- fread(file.path(s_dir, paste0(S, "_allele_counts.tsv.gz")))

    run_numbat(
        query_mat,          # query: Epithelial \ Ionocyte (benign epi = query 내부 diploid anchor)
        ref_internal,       # non-Epithelial + Ionocyte 로 만든 matched-normal 레퍼런스
        df_allele,          # step1 allele counts
        genome  = "hg38",
        t       = 1e-6,
        init_k  = INIT_K,   # epi 를 잘게 분할 기본값 3
        ncores  = N_CORES,
        plot    = TRUE,
        out_dir = out_dir,
        call_clonal_loh = TRUE
    )
    message(S, ": run_numbat done -> ", out_dir)
}
message("\n[02b] all samples done.")
