suppressMessages({
    library(Seurat)
    library(copykat)
})
source("scripts/00_utils/scRNA_utils.R")
set.seed(42)

OBJ_RDS  <- "Results/01_Integrated/combined_CRPC.rds"
OUT_ROOT <- "Results/08_CopyKAT"
N_CORES  <- 30
dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

obj <- readRDS(OBJ_RDS)
RAW_DIRS <- c(
    CRPC1 = "Raw_data/CRPC1/filtered_feature_bc_matrix",
    CRPC2 = "Raw_data/CRPC2/filtered_feature_bc_matrix",
    CRPC3 = "Raw_data/CRPC3/filtered_feature_bc_matrix"
)
samples <- names(RAW_DIRS)

run_copykat <- function(sample_id) {
    obj_cells <- colnames(obj)[obj$orig.ident == sample_id]   # QC-통과 세포(접두어 P#_)
    bc_raw    <- sub("^[^_]*_", "", obj_cells)                # 접두어 제거 → CellRanger 바코드

    raw  <- Read10X(RAW_DIRS[[sample_id]], gene.column = 1)   # Ensembl ID rowname
    stopifnot(all(bc_raw %in% colnames(raw)))
    mat  <- as.matrix(raw[, bc_raw])
    colnames(mat) <- obj_cells                               # 바코드로 복원
    
ANCHOR_TYPES <- c("T/NK cells", "Phagocytes", "Mast cells")
norm_cells <- obj_cells[obj@meta.data[obj_cells, "celltype"] %in% ANCHOR_TYPES]
    message(sprintf("[%s] cells=%d  anchor(non-epi)=%d  genes=%d",
                    sample_id, ncol(mat), length(norm_cells), nrow(mat)))

    out_dir <- file.path(OUT_ROOT, sample_id)
    dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
    old_wd <- getwd(); on.exit(setwd(old_wd), add = TRUE); setwd(out_dir)

    res <- copykat(rawmat = mat, id.type = "E", ngene.chr = 5, win.size = 25, KS.cut = 0.1,
                   norm.cell.names = norm_cells, sam.name = sample_id, distance = "euclidean",
                   genome = "hg20", n.cores = N_CORES)
    setwd(old_wd)

    saveRDS(res, file.path(out_dir, paste0(sample_id, "_copykat.rds")))
    pred <- res$prediction
    if (!is.null(pred)) {
        pred$sample <- sample_id
        write.csv(pred, file.path(out_dir, paste0(sample_id, "_copykat_prediction.csv")), row.names = FALSE)
    }
    pred
}

pred_list <- lapply(samples, function(s)
    tryCatch(run_copykat(s), error = function(e) {
        message(sprintf("[%s] failed: %s", s, conditionMessage(e))); NULL }))
pred_all <- do.call(rbind, Filter(Negate(is.null), pred_list))
if (!is.null(pred_all))
    write.csv(pred_all, file.path(OUT_ROOT, "copykat_prediction_all_samples.csv"), row.names = FALSE)
message("copyKAT 완료 → ", OUT_ROOT)
