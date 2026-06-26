# Epithelial_PROGENy_pseudobulk_activity.R
# PROGENy activity scoring (contrast 없이 전체 cluster) — saezlab pw_sc.html
# vignette의 pseudobulk variant.
#   patient × cluster (Mode A)  와  cluster only (Mode B) 두 가지로 pseudobulk →
#   log2(CPM+1) → decoupleR::run_mlm() against PROGENy human network →
#   14 pathways × N profiles activity matrix → heatmap.
#
# Contrast 스크립트(Epithelial_PROGENy_pseudobulk_contrasts.R)와 같은 net을 사용해
# "ARPC 기준 contrast"와 "전체 cluster 활성 지도"를 한 짝으로 본다.

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(decoupleR)
    library(dplyr)
    library(tidyr)
    library(pheatmap)
})
suppressMessages(source("scripts/00_utils/scRNA_utils.R"))

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/PROGENy_pseudobulk_activity"

MIN_CELLS_PER_PB <- 20
MIN_GENE_TOTAL   <- 10
PROGENY_TOP      <- 500

CLUSTER_ORDER <- c(
    "ARPC", "Club-like", "Hillock-like",
    "BE 1", "BE 2",
    "OE 1", "OE 2", "OE 3", "OE 4", "OE 5",
    "Ionocyte"
)

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

# ============================================================
# Load + group ----
# ============================================================
epi <- readRDS(IN_RDS)
DefaultAssay(epi) <- "RNA"
message("Loaded: ", IN_RDS, " | cells=", ncol(epi))

md <- epi@meta.data
grp <- as.character(md$annotation)
grp[grp %in% c("Hillock-like 1", "Hillock-like 2")] <- "Hillock-like"
md$group   <- grp
md$patient <- as.character(md$orig.ident)

counts <- GetAssayData(epi, assay = "RNA", layer = "counts")
patient_order <- sort(unique(md$patient))

# ============================================================
# Pseudobulk aggregation + log2(CPM+1) ----
# ============================================================
aggregate_pb <- function(counts_mat, key_vec) {
    f <- factor(key_vec)
    ind <- Matrix::sparse.model.matrix(~ 0 + f)
    colnames(ind) <- levels(f)
    pb <- as.matrix(counts_mat %*% ind)
    storage.mode(pb) <- "integer"
    pb
}
normalize_pb <- function(pb) {
    pb <- pb[rowSums(pb) >= MIN_GENE_TOTAL, , drop = FALSE]
    libsize <- colSums(pb)
    cpm <- t(t(pb) / libsize) * 1e6
    log2(cpm + 1)
}

# Mode A: patient × cluster
pb_key_A <- paste(md$patient, md$group, sep = "__")
pb_A <- aggregate_pb(counts, pb_key_A)
nc_A <- as.integer(table(pb_key_A)[colnames(pb_A)])
pb_A   <- pb_A[, nc_A >= MIN_CELLS_PER_PB, drop = FALSE]
expr_A <- normalize_pb(pb_A)
message("Mode A (patient × cluster): ", ncol(expr_A), " profiles, ",
        nrow(expr_A), " genes")

# Mode B: cluster only
pb_key_B <- md$group
pb_B <- aggregate_pb(counts, pb_key_B)
nc_B <- as.integer(table(pb_key_B)[colnames(pb_B)])
pb_B   <- pb_B[, nc_B >= MIN_CELLS_PER_PB, drop = FALSE]
expr_B <- normalize_pb(pb_B)
message("Mode B (cluster only): ", ncol(expr_B), " profiles, ",
        nrow(expr_B), " genes")

# ============================================================
# PROGENy network + decoupleR run_mlm ----
# ============================================================
net <- get_progeny(organism = "human", top = PROGENY_TOP)
message(sprintf("PROGENy network: %d pathways, %d interactions",
                length(unique(net$source)), nrow(net)))

run_progeny <- function(expr_mat, tag) {
    set.seed(1)
    res <- run_mlm(mat = expr_mat, net = net,
                   .source = "source", .target = "target", .mor = "weight",
                   minsize = 5)
    res <- res %>% filter(statistic == "mlm")
    to_mat <- function(value_col) {
        wide <- res %>%
            select(source, condition, all_of(value_col)) %>%
            pivot_wider(names_from = condition, values_from = all_of(value_col)) %>%
            as.data.frame()
        rownames(wide) <- wide$source
        wide$source <- NULL
        as.matrix(wide)
    }
    score_mat <- to_mat("score")
    pval_mat  <- to_mat("p_value")
    write.csv(score_mat, file.path(OUT_DIR, sprintf("PROGENy_%s_score.csv",  tag)))
    write.csv(pval_mat,  file.path(OUT_DIR, sprintf("PROGENy_%s_pvalue.csv", tag)))
    list(score = score_mat, pvalue = pval_mat)
}

modeA_out <- run_progeny(expr_A, "SampleXCluster")
modeB_out <- run_progeny(expr_B, "Cluster")

# ============================================================
# Heatmaps ----
# ============================================================
order_cols_A <- function(cols) {
    want <- as.vector(outer(CLUSTER_ORDER, patient_order,
                            function(g, p) paste(p, g, sep = "__")))
    intersect(want, cols)
}
order_cols_B <- function(cols) intersect(CLUSTER_ORDER, cols)

plot_heatmap <- function(score_mat, fname, title, col_order,
                         ann_col = NA, ann_colors = NA) {
    score_mat <- score_mat[, col_order, drop = FALSE]
    rng <- max(abs(score_mat), na.rm = TRUE)
    pheatmap(
        score_mat,
        color  = colorRampPalette(c("steelblue", "white", "firebrick"))(50),
        breaks = seq(-rng, rng, length.out = 51),  # 0 = white (diverging)
        cluster_cols = FALSE,
        cluster_rows = TRUE,
        scale  = "none",                           # MLM score는 z-score 안 함
        show_rownames = TRUE, show_colnames = TRUE,
        fontsize_row  = 9, fontsize_col = 9,
        annotation_col   = ann_col,
        annotation_colors = ann_colors,
        main     = title,
        filename = fname,
        width    = max(8, ncol(score_mat) * 0.5 + 4),
        height   = max(6, nrow(score_mat) * 0.32 + 2)
    )
}

coldata_A <- data.frame(
    patient = sub("__.*$", "", colnames(modeA_out$score)),
    cluster = sub("^.*__", "", colnames(modeA_out$score)),
    row.names = colnames(modeA_out$score)
)
ann_colors_A <- list(
    patient = setNames(utils_cb_palette(length(patient_order)), patient_order),
    cluster = setNames(utils_cb_palette(length(CLUSTER_ORDER)), CLUSTER_ORDER)
)

plot_heatmap(modeA_out$score,
             file.path(OUT_DIR, "PROGENy_SampleXCluster_heatmap.png"),
             "PROGENy activity (MLM) — patient × cluster",
             order_cols_A(colnames(modeA_out$score)),
             ann_col = coldata_A, ann_colors = ann_colors_A)

plot_heatmap(modeB_out$score,
             file.path(OUT_DIR, "PROGENy_Cluster_heatmap.png"),
             "PROGENy activity (MLM) — cluster (patients pooled)",
             order_cols_B(colnames(modeB_out$score)))

message("\nDone. Outputs in: ", OUT_DIR)
