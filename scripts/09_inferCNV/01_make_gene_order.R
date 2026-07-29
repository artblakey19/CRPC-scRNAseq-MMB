#!/usr/bin/env Rscript
# =============================================================================
# 09 inferCNV — step 01 : GRCh38 gene ordering file 생성 (환자 공통)
# -----------------------------------------------------------------------------
# inferCNV gene_order_file 포맷: <gene>\t<chr>\t<start>\t<stop> (헤더 없음)
#   gene = count_mat rownames 와 동일해야 함 (=Read10X make.unique 심볼).
# 전략: count_mat 행순서 = features.tsv.gz 행순서 (Read10X 가 보존; 검증됨).
#   각 행의 ENSEMBL id 를 features.tsv.gz 에서 위치기반으로 복원하고, 좌표는
#   Ensembl GTF(GRCh38, AnnotationHub)에서 ENSG 기준으로 가져온다.
#   ENSG 는 버전 간 안정 → GTF 릴리스와 무관하게 좌표 정확(심볼 drift 회피).
# =============================================================================
suppressPackageStartupMessages({ library(AnnotationHub); library(GenomicRanges) })

FEAT <- "Raw_data/CRPC1/filtered_feature_bc_matrix/features.tsv.gz"  # 세 환자 동일 레퍼런스
CM   <- "Results/08_Numbat/CRPC1/count_mat.rds"
OUT  <- "Results/09_inferCNV/gene_order_GRCh38.txt"
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)

# 1) features (ENSG, symbol) — matrix 행순서 그대로 --------------------------
feat  <- read.delim(FEAT, header = FALSE, stringsAsFactors = FALSE)
cm_rn <- rownames(readRDS(CM))
stopifnot(nrow(feat) == length(cm_rn),
          identical(make.unique(feat$V2), cm_rn))   # 위치기반 ENSG 복원 유효성
map <- data.frame(row_name = cm_rn, ensg = feat$V1, stringsAsFactors = FALSE)

# 2) Ensembl GTF (GRCh38) 좌표 — AnnotationHub, release 110 우선 -------------
#    (CellRanger 2024-A 레퍼런스가 Ensembl 110/GENCODE v44 기반)
ah  <- AnnotationHub()
q   <- query(ah, c("Homo sapiens", "Ensembl", "GRCh38", "gtf"))
rel <- suppressWarnings(as.integer(sub(".*GRCh38\\.(\\d+)\\.gtf.*", "\\1", q$title)))
sel <- if (any(rel == 110, na.rm = TRUE)) which(rel == 110)[1] else which.max(rel)
cat(sprintf("[01] AnnotationHub %s : %s\n", names(q)[sel], q$title[sel]))
gtf <- ah[[ names(q)[sel] ]]
g   <- gtf[gtf$type == "gene"]
gd  <- data.frame(ensg  = g$gene_id,
                  chr   = as.character(seqnames(g)),
                  start = start(g), stop = end(g), stringsAsFactors = FALSE)

# 3) join → 표준 염색체만 → 유전체 순서 정렬 ---------------------------------
std <- c(1:22, "X", "Y")
m <- merge(map, gd, by = "ensg")
m <- m[m$chr %in% std, ]
m$chr <- factor(m$chr, levels = std)
m <- m[order(m$chr, m$start), ]
go <- data.frame(chr = paste0("chr", as.character(m$chr)), start = m$start, stop = m$stop)
rownames(go) <- m$row_name
write.table(go, OUT, sep = "\t", quote = FALSE, col.names = FALSE)

cat(sprintf("[01] gene_order rows: %d  (dropped %d unmapped/non-std)\n",
            nrow(go), length(cm_rn) - nrow(go)))
print(table(factor(go$chr, levels = paste0("chr", std))))
cat("[01] done ->", OUT, "\n")
