library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(dplyr)
library(ggplot2)
library(patchwork)

# Multi-core processing
library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

dir.create("Results/Epithelial/Monocle3", showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load: clustered epithelial subset ----
# ============================================================
epi <- readRDS("Results/Epithelial/epi_clustered.rds")

# Curated res 0.4 clusters are stored in `seurat_clusters`.
Idents(epi) <- "seurat_clusters"
DefaultAssay(epi) <- "SCT"

# ============================================================
# Seurat → cell_data_set ----
# ============================================================
# as.cell_data_set carries over counts + UMAP. We then overwrite the PCA slot
# with the Harmony embedding so any downstream step that consumes "PCA"
# operates on the patient-batch-corrected basis used in Seurat's FindNeighbors.
cds <- as.cell_data_set(epi)
reducedDims(cds)$PCA <- Embeddings(epi, "harmony")[, 1:30]
rowData(cds)$gene_short_name <- rownames(cds)

# ============================================================
# Cluster cells + learn principal graph ----
# ============================================================
# monocle3 needs its own clusters/partitions before learn_graph. We cluster on
# the inherited Seurat UMAP so the partitioning aligns with the visualization.
cds <- cluster_cells(cds, reduction_method = "UMAP")

# use_partition = FALSE → single connected graph spanning all partitions, so a
# trajectory can run from the cluster 8 (diploid anchor) into malignant DNPC
# clusters that monocle3 might otherwise split into a separate partition.
cds <- learn_graph(cds, use_partition = FALSE, close_loop = FALSE)

# ============================================================
# Root: cluster 8 ----
# ============================================================
# Cluster 8 is the diploid/normal anchor used as norm.cell.names in the
# epithelial-only copyKAT rerun, so it is the natural origin for pseudotime.
root_cells <- colnames(epi)[epi$seurat_clusters == "8"]
cds <- order_cells(cds, root_cells = root_cells)

epi$monocle3_pseudotime <- pseudotime(cds)[colnames(epi)]

# ============================================================
# Trajectory plots ----
# ============================================================
p_pseudo <- plot_cells(cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
    cell_size = 0.4
) + ggtitle("Monocle3 pseudotime (root = cluster 8)")
ggsave("Results/Epithelial/Monocle3/UMAP_pseudotime.png",
    plot = p_pseudo, width = 10, height = 8, bg = "white"
)

p_clu <- plot_cells(cds,
    color_cells_by = "seurat_clusters",
    label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = TRUE,
    cell_size = 0.4, group_label_size = 5
) + ggtitle("Monocle3 graph by cluster")
ggsave("Results/Epithelial/Monocle3/UMAP_graph_clusters.png",
    plot = p_clu, width = 10, height = 8, bg = "white"
)

p_part <- plot_cells(cds,
    color_cells_by = "partition",
    label_cell_groups = TRUE, label_leaves = FALSE, label_branch_points = FALSE,
    cell_size = 0.4
) + ggtitle("Monocle3 partitions")
ggsave("Results/Epithelial/Monocle3/UMAP_partitions.png",
    plot = p_part, width = 10, height = 8, bg = "white"
)

p_pat <- plot_cells(cds,
    color_cells_by = "orig.ident",
    label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
    cell_size = 0.4
) + ggtitle("Monocle3 graph by patient")
ggsave("Results/Epithelial/Monocle3/UMAP_graph_patient.png",
    plot = p_pat, width = 10, height = 8, bg = "white"
)

# copyKAT overlay — diploid cluster 8 should sit at low pseudotime,
# aneuploid mass should occupy higher pseudotime.
copykat_cols <- c(
    "aneuploid" = "#D7261E", "diploid" = "#1F77B4", "not.defined" = "grey80"
)
p_ck <- plot_cells(cds,
    color_cells_by = "copykat_prediction",
    label_cell_groups = FALSE, label_leaves = FALSE, label_branch_points = FALSE,
    cell_size = 0.4
) + scale_color_manual(values = copykat_cols) +
    ggtitle("Monocle3 graph by copyKAT prediction")
ggsave("Results/Epithelial/Monocle3/UMAP_graph_copykat.png",
    plot = p_ck, width = 10, height = 8, bg = "white"
)

# ============================================================
# Pseudotime distribution by cluster / patient ----
# ============================================================
p_box_clu <- ggplot(epi@meta.data,
    aes(x = reorder(seurat_clusters, monocle3_pseudotime, FUN = median),
        y = monocle3_pseudotime, fill = seurat_clusters)
) +
    geom_boxplot(outlier.size = 0.3) +
    labs(x = "Cluster (sorted by median pseudotime)",
         y = "Pseudotime (root = cluster 8)") +
    theme_classic() + NoLegend()
ggsave("Results/Epithelial/Monocle3/Pseudotime_by_cluster.png",
    plot = p_box_clu, width = 10, height = 6, bg = "white"
)

p_box_pat <- ggplot(epi@meta.data,
    aes(x = orig.ident, y = monocle3_pseudotime, fill = orig.ident)
) +
    geom_boxplot(outlier.size = 0.3) +
    labs(x = "Patient", y = "Pseudotime (root = cluster 8)") +
    theme_classic() + NoLegend()
ggsave("Results/Epithelial/Monocle3/Pseudotime_by_patient.png",
    plot = p_box_pat, width = 8, height = 6, bg = "white"
)

# ============================================================
# Marker dynamics along pseudotime ----
# ============================================================
# Track ARPC → DNPC axis as a function of pseudotime.
traj_markers <- c(
    "AR", "KLK3", "FOLH1", "NKX3-1",        # ARPC
    "FGF8", "FGFR1", "CCL2", "DKK1",        # DNPC-FGF
    "KRT7", "SOX2", "FOXA2",                # DNPC-KRT7
    "CHGA", "SYP", "PROX1"                  # NEPC
)
traj_markers <- intersect(traj_markers, rownames(cds))

p_genes <- plot_genes_in_pseudotime(
    cds[traj_markers, !is.infinite(pseudotime(cds))],
    color_cells_by = "seurat_clusters",
    min_expr = 0.1, ncol = 3
)
ggsave("Results/Epithelial/Monocle3/Genes_in_pseudotime.png",
    plot = p_genes, width = 14, height = 18, bg = "white"
)

# ============================================================
# Save ----
# ============================================================
saveRDS(cds, "Results/Epithelial/Monocle3/cds_monocle3.rds")
saveRDS(epi, "Results/Epithelial/Monocle3/epi_with_pseudotime.rds")
