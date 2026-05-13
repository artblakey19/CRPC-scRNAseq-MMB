library(Seurat)
library(copykat)
library(dplyr)
# library(copykat) loads plyr, which masks dplyr::count() — qualify with
# dplyr:: when needed (e.g. dplyr::count).

# Multi-core processing (matches other scripts in this project)
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)
Sys.setenv(OMP_NUM_THREADS = "30", OPENBLAS_NUM_THREADS = "30", MKL_NUM_THREADS = "30")

out_dir <- file.path(getwd(), "Results/Epithelial/copyKAT_clu8ref")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load: epithelial subset ----
# ============================================================
epi <- readRDS("Results/Epithelial/epi_clustered.rds")
DefaultAssay(epi) <- "RNA"
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])

# Cluster 8 cells as the diploid reference (predominantly diploid per the
# integrated copyKAT run — see Results/Epithelial/copykat_cluster_composition.csv).
ref_cells <- colnames(epi)[as.character(epi$seurat_clusters) == "8"]
message("Reference cells (cluster 8): n = ", length(ref_cells))

raw_counts <- GetAssayData(epi, assay = "RNA", layer = "counts")
raw_counts <- as.matrix(raw_counts)
storage.mode(raw_counts) <- "numeric"

# ============================================================
# Run copyKAT ----
# ============================================================
# copyKAT writes intermediate files to the working directory, so chdir in.
old_wd <- getwd()
on.exit(setwd(old_wd), add = TRUE)
setwd(out_dir)

copykat_res <- copykat::copykat(
    rawmat          = raw_counts,
    id.type         = "S",
    cell.line       = "no",
    ngene.chr       = 5,
    win.size        = 25,
    norm.cell.names = ref_cells,
    KS.cut          = 0.1,
    sam.name        = "epi_clu8ref",
    distance        = "euclidean",
    genome          = "hg20",
    n.cores         = 30
)

setwd(old_wd)

saveRDS(copykat_res, file.path(out_dir, "epi_clu8ref_copykat.rds"))

# ============================================================
# Overlay new predictions on the epithelial UMAP ----
# ============================================================
pred <- copykat_res$prediction
# copyKAT prediction column naming varies; pick whichever is present.
cell_col <- intersect(c("cell.names", "cell_names", "cells"), names(pred))[1]
pred_col <- intersect(c("copykat.pred", "copykat_pred", "prediction"), names(pred))[1]
pred_df <- data.frame(
    cell = pred[[cell_col]],
    copykat_clu8ref = pred[[pred_col]],
    stringsAsFactors = FALSE
)
write.csv(pred_df, file.path(out_dir, "epi_clu8ref_predictions.csv"), row.names = FALSE)

epi$copykat_clu8ref <- pred_df$copykat_clu8ref[match(colnames(epi), pred_df$cell)]

copykat_cols <- c(
    "aneuploid" = "#D7261E", "diploid" = "#1F77B4", "not.defined" = "grey80"
)

library(ggplot2)
p_pred <- DimPlot(epi,
    group.by = "copykat_clu8ref", pt.size = 0.3, cols = copykat_cols
) + ggtitle("copyKAT (cluster 8 reference)")
ggsave(file.path(out_dir, "UMAP_copykat_clu8ref.png"),
    plot = p_pred, width = 10, height = 8, bg = "white"
)

p_pred_split <- DimPlot(epi,
    group.by = "copykat_clu8ref", split.by = "orig.ident",
    pt.size = 0.3, cols = copykat_cols
)
ggsave(file.path(out_dir, "UMAP_copykat_clu8ref_by_patient.png"),
    plot = p_pred_split, width = 24, height = 8, bg = "white"
)

# Per-cluster composition under the new reference choice.
clu8_table <- epi@meta.data %>%
    dplyr::count(seurat_clusters, copykat_clu8ref, name = "n") %>%
    tidyr::pivot_wider(
        names_from = copykat_clu8ref, values_from = n, values_fill = 0
    )
num_cols <- setdiff(names(clu8_table), "seurat_clusters")
clu8_table$total <- rowSums(clu8_table[, num_cols, drop = FALSE])
if ("aneuploid" %in% num_cols) {
    clu8_table$aneuploid_frac <- clu8_table$aneuploid / clu8_table$total
    clu8_table <- clu8_table[order(-clu8_table$aneuploid_frac), ]
}
write.csv(clu8_table,
    file.path(out_dir, "epi_clu8ref_cluster_composition.csv"), row.names = FALSE
)

p_bar <- ggplot(epi@meta.data,
    aes(x = seurat_clusters, fill = copykat_clu8ref)
) +
    geom_bar(position = "fill") +
    scale_fill_manual(values = copykat_cols) +
    labs(x = "Cluster", y = "Proportion", fill = "copyKAT (clu8 ref)") +
    theme_classic()
ggsave(file.path(out_dir, "epi_clu8ref_cluster_barplot.png"),
    plot = p_bar, width = 10, height = 6, bg = "white"
)
