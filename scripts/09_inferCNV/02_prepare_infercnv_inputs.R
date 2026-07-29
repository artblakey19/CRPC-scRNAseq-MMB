#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 02 : 환자별 annotations 파일 생성
# -----------------------------------------------------------------------------
# 08 의 count_mat.rds(원 CellRanger raw UMI) 모든 세포에 group 라벨 부여:
#   - reference 세포(비-상피 7종 + Ionocyte anchor) = cell_annot_ref.csv 의 group
#   - query 상피 세포                               = epi_annotated.rds fine annotation
#                                                     (LE(ARPC)/Club/Hillock/BE1-6/OE)
# → Numbat 과 완전히 동일한 query/reference 구성.
# annotations.txt 포맷: <cell>\t<group> (헤더 없음, tab 구분 → 라벨 내 공백 허용)
# =============================================================================
suppressPackageStartupMessages({ library(Seurat) })
NB  <- "Results/08_Numbat"
OUT <- "Results/09_inferCNV"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")

epi <- readRDS("Results/05_Epithelial_Downstream/epi_annotated.rds")
epi_key   <- paste(as.character(epi$orig.ident),
                   sub("^P[0-9]+_", "", colnames(epi)), sep = "|")  # sample|barcode
epi_annot <- setNames(as.character(epi$annotation), epi_key)
rm(epi); gc()

for (S in SAMPLES) {
    cm    <- readRDS(file.path(NB, S, "count_mat.rds"))
    cells <- colnames(cm)
    ref   <- read.csv(file.path(NB, S, "cell_annot_ref.csv"), stringsAsFactors = FALSE)  # cell, group

    grp <- setNames(rep(NA_character_, length(cells)), cells)
    grp[ref$cell] <- ref$group                       # reference (Ionocyte 포함)
    q <- cells[is.na(grp)]                            # 나머지 = query 상피
    grp[q] <- epi_annot[paste(S, q, sep = "|")]
    stopifnot(!anyNA(grp))                            # 모든 세포에 라벨

    dir.create(file.path(OUT, S), recursive = TRUE, showWarnings = FALSE)
    ann <- data.frame(cell = cells, group = unname(grp))
    write.table(ann, file.path(OUT, S, "annotations.txt"),
                sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
    cat(sprintf("\n%s: %d cells (query %d / ref %d)\n",
                S, length(cells), length(q), nrow(ref)))
    print(table(grp))
}
cat("\n[02] done.\n")
