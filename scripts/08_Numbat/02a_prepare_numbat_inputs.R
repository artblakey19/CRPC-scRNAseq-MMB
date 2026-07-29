#!/usr/bin/env Rscript
# =============================================================================
# Numbat step 2a  —  HOST / renv  (Seurat side)
# -----------------------------------------------------------------------------
# 컨테이너에 Seurat이 없어 host에서 numbat 입력을 준비하는 스크립트
# 이 스크립트는 Seurat 으로 combined_CRPC.rds를 읽어, 환자별로 numbat 입력을 파일로 내보낸다.
#   - count_mat.rds      : gene x cell 정수 UMI (CellRanger raw count), QC 통과 세포
#   - cell_annot_ref.csv : 레퍼런스로 쓸 세포 (cell, group=celltype)
#
# 레퍼런스 = 01_Integrated 주석에서 Epithelial 을 제외한 모든 세포
#   Fibroblast / Endothelial / Smooth muscle cells / T/NK cells / Phagocytes / Mast cells
#   + Ionocyte (상피 세부주석의 정상 상피 서브타입; epi_annotated.rds 에서 가져옴)
#     → 상피 계통의 diploid normal anchor 를 레퍼런스에 추가 (stromal/immune 만으로는
#       상피 특이 발현이 레퍼런스에 결핍 → epithelial baseline 보강). query 에서는 제외.
# barcode: combined_CRPC.rds는 P1_/P2_/P3_ 접두사 → allele df 의 cell(AAAC..-1)과 맞추려 접두사 제거
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

INT_RDS <- "Results/01_Integrated/combined_CRPC.rds"
EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"   # 상피 세부주석(Ionocyte 라벨)
NB_DIR  <- "Results/08_Numbat"
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")

obj <- readRDS(INT_RDS)

# combined_CRPC.rds의 주석(barcode + celltype)만 annot에 저장
# count값은 Read10X()로 raw count를 가져옴.
annot <- data.frame(
    cell     = sub("^P[0-9]+_", "", colnames(obj)),
    sample   = as.character(obj$orig.ident),
    celltype = as.character(obj$celltype),
    stringsAsFactors = FALSE
)

# --- Ionocyte 를 레퍼런스로 승격 --------------------------------------------
# coarse celltype(combined_CRPC)은 Ionocyte 를 'Epithelial' 로 뭉뚱그린다. 상피 세부주석
# (epi_annotated.rds)에서 Ionocyte 바코드를 가져와 annot$celltype 을 'Ionocyte' 로 재지정.
# 그러면 아래 non-Epithelial 레퍼런스 선택에 자연히 포함되고, query(=Epithelial)에서는 빠진다.
# 매칭 키 = sample|stripped-barcode (환자 간 10x 바코드 충돌 방지).
epi <- readRDS(EPI_RDS)
iono_key <- paste(as.character(epi$orig.ident),
                  sub("^P[0-9]+_", "", colnames(epi)), sep = "|")[epi$annotation == "Ionocyte"]
rm(epi); gc()

annot_key <- paste(annot$sample, annot$cell, sep = "|")
is_iono   <- annot_key %in% iono_key
stopifnot(sum(is_iono) == length(iono_key))          # 모든 Ionocyte 가 annot 에 매칭되어야 함
annot$celltype[is_iono] <- "Ionocyte"
cat(sprintf("Ionocyte promoted to reference: %d cells\n", sum(is_iono)))

for (S in SAMPLES) {
    dir.create(file.path(NB_DIR, S), recursive = TRUE, showWarnings = FALSE)
    a <- annot[annot$sample == S, ]

    # count_mat = 환자 전체 세포(raw UMI). 단 02b 에서 query 는 Epithelial 만 subset,
    #   non-epi 는 레퍼런스 build 용으로만 사용 (epi-only CNV inference).
    # count 는 원 CellRanger raw UMI (decontX 아님). 객체 세포 ⊆ raw filtered 세포 (검증됨).
    raw <- Read10X(file.path("Raw_data", S, "filtered_feature_bc_matrix"))
    stopifnot(all(a$cell %in% colnames(raw)))
    cm <- raw[, a$cell, drop = FALSE]
    saveRDS(cm, file.path(NB_DIR, S, "count_mat.rds"))

    # reference = 비-Epithelial 세포 (matched normal) + Ionocyte(정상 상피 anchor), celltype 별 group
    ref <- a[a$celltype != "Epithelial", c("cell", "celltype")]
    colnames(ref) <- c("cell", "group")
    write.csv(ref, file.path(NB_DIR, S, "cell_annot_ref.csv"), row.names = FALSE)

    cat(sprintf("%s: query %d cells (raw counts) | reference %d cells across %d types (incl. Ionocyte)\n",
                S, ncol(cm), nrow(ref), length(unique(ref$group))))
    print(table(ref$group))
}
cat("\n[02a] done. Next: 02b (container) via 02_run_numbat.sh\n")
