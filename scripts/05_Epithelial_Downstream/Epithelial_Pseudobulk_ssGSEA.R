# Epithelial_Pseudobulk_ssGSEA.R
# Filtered epithelial을 두 가지 방식으로 pseudobulk 한 뒤
# 네 collection (Hallmark / WikiPathways / C6 oncogenic / C5 GO:BP) ssGSEA.
#
#   Mode A) Sample × Cluster  (patient-별 효과 보존)
#   Mode B) Cluster only      (3 patient 합산, 10 cluster)
#
# Note: Hillock-like 1+2는 Bluemn pseudobulk 스크립트 관례를 따라 통합 → 10 그룹.
# 입력은 RNA raw counts → sum aggregate → log2(CPM+1) → ssGSEA(GSVA, alpha=0.25).

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(GSVA)
    library(msigdbr)
    library(dplyr)
    library(ggplot2)
    library(pheatmap)
})
suppressMessages(source("scripts/00_utils/scRNA_utils.R"))

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Pseudobulk_ssGSEA"

MIN_CELLS_PER_PB <- 20    # Sample×Cluster에서 셀 수 적은 조합 제거
MIN_GENE_TOTAL   <- 10    # 저발현 유전자 제거
MIN_SET_SIZE     <- 10
MAX_SET_SIZE     <- 500
SSGSEA_ALPHA     <- 0.25

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

# =============================================================
# Load + group construction
# =============================================================
epi <- readRDS(IN_RDS)
DefaultAssay(epi) <- "RNA"
message("Loaded: ", IN_RDS, " | cells=", ncol(epi), " genes=", nrow(epi))

md <- epi@meta.data
stopifnot("annotation column missing"  = "annotation"  %in% colnames(md))
stopifnot("orig.ident column missing"  = "orig.ident"  %in% colnames(md))

grp <- as.character(md$annotation)
grp[grp %in% c("Hillock-like 1", "Hillock-like 2")] <- "Hillock-like"
md$group   <- grp
md$patient <- as.character(md$orig.ident)

cluster_order <- c(
    "ARPC", "Club-like", "Hillock-like",
    "BE 1", "BE 2",
    "OE 1", "OE 2", "OE 3", "OE 4",
    "Ionocyte-like"
)
patient_order <- sort(unique(md$patient))   # CRPC1/2/3

# =============================================================
# Pseudobulk aggregation (sum of RNA counts)
# =============================================================
counts <- GetAssayData(epi, assay = "RNA", layer = "counts")

aggregate_pb <- function(counts_mat, key_vec) {
    f <- factor(key_vec)
    ind <- Matrix::sparse.model.matrix(~ 0 + f)
    colnames(ind) <- levels(f)
    pb <- as.matrix(counts_mat %*% ind)
    storage.mode(pb) <- "integer"
    pb
}

# Mode A: patient__group
pb_key_A <- paste(md$patient, md$group, sep = "__")
pb_A <- aggregate_pb(counts, pb_key_A)
nc_A <- as.integer(table(pb_key_A)[colnames(pb_A)])
keep_A <- nc_A >= MIN_CELLS_PER_PB
pb_A   <- pb_A[, keep_A, drop = FALSE]
coldata_A <- data.frame(
    pb_key  = colnames(pb_A),
    patient = sub("__.*$", "", colnames(pb_A)),
    group   = sub("^.*__", "", colnames(pb_A)),
    n_cells = nc_A[keep_A]
)
write.csv(coldata_A, file.path(OUT_DIR, "pb_SampleXCluster_design.csv"), row.names = FALSE)

# Mode B: group only
pb_key_B <- md$group
pb_B <- aggregate_pb(counts, pb_key_B)
nc_B <- as.integer(table(pb_key_B)[colnames(pb_B)])
keep_B <- nc_B >= MIN_CELLS_PER_PB
pb_B   <- pb_B[, keep_B, drop = FALSE]
coldata_B <- data.frame(
    pb_key  = colnames(pb_B),
    group   = colnames(pb_B),
    n_cells = nc_B[keep_B]
)
write.csv(coldata_B, file.path(OUT_DIR, "pb_Cluster_design.csv"), row.names = FALSE)

message("Mode A profiles (kept ≥", MIN_CELLS_PER_PB, " cells): ", ncol(pb_A))
print(coldata_A)
message("Mode B profiles: ", ncol(pb_B))
print(coldata_B)

saveRDS(list(pb_A = pb_A, pb_B = pb_B,
             coldata_A = coldata_A, coldata_B = coldata_B),
        file.path(OUT_DIR, "pseudobulk_counts.rds"))

# =============================================================
# Normalize: log2(CPM+1)  (ssGSEA 입력)
# =============================================================
normalize_pb <- function(pb) {
    pb <- pb[rowSums(pb) >= MIN_GENE_TOTAL, , drop = FALSE]
    libsize <- colSums(pb)
    cpm <- t(t(pb) / libsize) * 1e6
    log2(cpm + 1)
}
expr_A <- normalize_pb(pb_A)
expr_B <- normalize_pb(pb_B)
message("Expr matrix A: ", nrow(expr_A), " × ", ncol(expr_A))
message("Expr matrix B: ", nrow(expr_B), " × ", ncol(expr_B))

# =============================================================
# Gene sets (msigdbr v2026.1.Hs)
# =============================================================
hallmark <- msigdbr(species = "Homo sapiens", collection = "H")
wikipath <- msigdbr(species = "Homo sapiens", collection = "C2",
                    subcollection = "CP:WIKIPATHWAYS")
c6_onc   <- msigdbr(species = "Homo sapiens", collection = "C6")
gobp     <- msigdbr(species = "Homo sapiens", collection = "C5",
                    subcollection = "GO:BP")

to_list <- function(df) split(df$gene_symbol, df$gs_name)
COLLECTIONS <- list(
    Hallmark     = to_list(hallmark),
    WikiPathways = to_list(wikipath),
    C6_Oncogenic = to_list(c6_onc),
    C5_GOBP      = to_list(gobp)
)
for (nm in names(COLLECTIONS)) {
    message(sprintf("Collection %s: %d sets loaded", nm, length(COLLECTIONS[[nm]])))
}

# =============================================================
# ssGSEA driver
# =============================================================
run_ssgsea <- function(expr_mat, sets) {
    sets_f <- lapply(sets, function(g) intersect(g, rownames(expr_mat)))
    sizes  <- vapply(sets_f, length, integer(1))
    sets_f <- sets_f[sizes >= MIN_SET_SIZE & sizes <= MAX_SET_SIZE]
    param <- ssgseaParam(
        exprData  = expr_mat,
        geneSets  = sets_f,
        minSize   = MIN_SET_SIZE,
        maxSize   = MAX_SET_SIZE,
        alpha     = SSGSEA_ALPHA,
        normalize = TRUE
    )
    suppressMessages(gsva(param, verbose = FALSE))
}

modes <- list(SampleXCluster = expr_A, Cluster = expr_B)
score_matrices <- list()
for (mode_name in names(modes)) {
    for (coll_name in names(COLLECTIONS)) {
        message(sprintf("  ssGSEA: %-14s × %-12s ...", mode_name, coll_name))
        t0 <- Sys.time()
        sc <- run_ssgsea(modes[[mode_name]], COLLECTIONS[[coll_name]])
        message(sprintf("    done. %d sets × %d profiles  [%.1fs]",
                        nrow(sc), ncol(sc),
                        as.numeric(difftime(Sys.time(), t0, units = "secs"))))
        tag <- paste(mode_name, coll_name, sep = "_")
        score_matrices[[tag]] <- sc
        write.csv(sc, file.path(OUT_DIR, sprintf("ssGSEA_%s.csv", tag)))
    }
}
saveRDS(score_matrices, file.path(OUT_DIR, "ssGSEA_scores.rds"))

# =============================================================
# Heatmaps (per mode × collection)
# =============================================================
# 컬럼 정렬 헬퍼
order_cols_A <- function(cols) {
    # cluster를 outer로 두어 같은 annotation의 patient들이 인접하게 묶이도록.
    want <- as.vector(outer(patient_order, cluster_order,
                            function(p, g) paste(p, g, sep = "__")))
    intersect(want, cols)
}
order_cols_B <- function(cols) intersect(cluster_order, cols)

# annotation bar (Mode A 용)
ann_col_A <- data.frame(
    patient = coldata_A$patient,
    cluster = coldata_A$group,
    row.names = coldata_A$pb_key
)
ann_colors_A <- list(
    patient = setNames(utils_cb_palette(length(patient_order)), patient_order),
    cluster = setNames(utils_cb_palette(length(cluster_order)), cluster_order)
)

plot_heatmap <- function(sc, fname, title, col_order,
                         ann_col = NA, ann_colors = NA, top_n = NULL) {
    if (!is.null(top_n) && nrow(sc) > top_n) {
        v <- apply(sc, 1, var)
        sc <- sc[order(-v)[seq_len(top_n)], , drop = FALSE]
    }
    sc <- sc[, col_order, drop = FALSE]
    # 행별 z-score (diverging 표시 위해 row scale)
    rng <- max(abs(sc), na.rm = TRUE)
    n_row <- nrow(sc)
    fs_row <- if (n_row > 200) 4 else if (n_row > 100) 5 else 7
    pheatmap(
        sc,
        color  = colorRampPalette(c("steelblue", "white", "firebrick"))(50),
        cluster_cols  = FALSE,
        cluster_rows  = TRUE,
        scale         = "row",
        show_rownames = n_row <= 400,
        show_colnames = TRUE,
        fontsize_row  = fs_row,
        fontsize_col  = 8,
        annotation_col   = ann_col,
        annotation_colors = ann_colors,
        main     = title,
        filename = fname,
        width    = max(8, ncol(sc) * 0.45 + 4),
        height   = max(6, n_row * 0.14 + 2)
    )
}

# Hallmark/C6는 전체, WikiPathways/GO:BP는 top variable로 압축(여전히 충분히 많이).
top_n_per_collection <- list(
    Hallmark     = NULL,
    WikiPathways = 200,
    C6_Oncogenic = NULL,
    C5_GOBP      = 300
)

for (mode_name in names(modes)) {
    is_A <- mode_name == "SampleXCluster"
    for (coll_name in names(COLLECTIONS)) {
        tag <- paste(mode_name, coll_name, sep = "_")
        sc  <- score_matrices[[tag]]
        col_ord <- if (is_A) order_cols_A(colnames(sc)) else order_cols_B(colnames(sc))
        plot_heatmap(
            sc,
            fname  = file.path(OUT_DIR, sprintf("heatmap_%s.png", tag)),
            title  = sprintf("ssGSEA — %s (%s, row z-score)", coll_name, mode_name),
            col_order = col_ord,
            ann_col   = if (is_A) ann_col_A     else NA,
            ann_colors = if (is_A) ann_colors_A else NA,
            top_n     = top_n_per_collection[[coll_name]]
        )
    }
}

message("\nDone. Outputs in: ", OUT_DIR)
