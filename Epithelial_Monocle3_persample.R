# ============================================================
# Epithelial_Monocle3_persample.R
#
# Per-patient monocle3 trajectory / pseudotime on the epithelial subset.
#
# Rationale: the cohort is 3 independent late-stage DNPC patients. Learning a
# single principal graph on the Harmony-integrated embedding implicitly assumes
# all 3 tumours malignised along one shared path. Instead we run monocle3
# INDEPENDENTLY per patient on UNINTEGRATED data (fresh per-patient
# SCT/PCA/UMAP), so each trajectory reflects only that patient's own structure.
#
# The Harmony-integrated clustering is used ONLY to define the root: cluster 8
# (diploid anchor, the copyKAT norm.cell.names population). It is carried in as
# a cell-level label (`integrated_cluster`); order_cells() takes root cells as
# explicit barcodes, so the per-patient monocle3 clustering can number its own
# clusters freely without conflict.
#
# Input: epi_copykat_perpatient.rds — the epi_clustered object with the
# per-patient cluster-8-referenced copyKAT columns merged in
# (Epithelial_copyKAT_perpatient.R). The copyKAT overlay uses `copykat_clu8ref`
# from that rerun, NOT the stale immune-referenced `copykat_prediction` that
# was inherited from the integrated object.
# ============================================================

library(Seurat)
library(SeuratWrappers)
library(monocle3)
library(glmGamPoi)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

# ============================================================
# Config ----
# ============================================================
ROOT_CLUSTER <- "8"        # integrated diploid anchor cluster used as root
UMAP_DIMS    <- 1:30
OUT_DIR      <- "Results/Epithelial/Monocle3_persample"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

copykat_cols <- c(
    "aneuploid" = "#D7261E", "diploid" = "#1F77B4", "not.defined" = "grey80"
)

# Track ARPC -> DNPC -> NEPC axis as a function of pseudotime.
traj_markers <- c(
    "AR", "KLK3", "FOLH1", "NKX3-1", # ARPC
    "FGF8", "FGFR1", "CCL2", "DKK1", # DNPC-FGF
    "KRT7", "SOX2", "FOXA2",         # DNPC-KRT7
    "CHGA", "SYP", "PROX1"           # NEPC
)

# ============================================================
# Load epithelial subset ----
# ============================================================
# epi_clustered + per-patient cluster-8-referenced copyKAT columns.
epi <- readRDS("Results/Epithelial/copyKAT_perpatient/epi_copykat_perpatient.rds")

# Preserve the integrated res-0.4 clustering as a carried label. It is used for
# root definition and as a colour overlay only — NOT for the per-patient graph.
epi$integrated_cluster <- as.character(epi$seurat_clusters)

DefaultAssay(epi) <- "RNA"
epi[["RNA"]] <- JoinLayers(epi[["RNA"]])

patients <- sort(unique(epi$orig.ident))
message("Patients: ", paste(patients, collapse = ", "))

# ============================================================
# Per-patient trajectory ----
# ============================================================
run_patient <- function(sample_id) {
    message("\n==== monocle3: ", sample_id, " ====")
    out <- file.path(OUT_DIR, sample_id)
    dir.create(out, showWarnings = FALSE, recursive = TRUE)

    # --- subset + fresh unintegrated SCT/PCA/UMAP ---
    sub <- subset(epi, subset = orig.ident == sample_id)
    sub[["SCT"]] <- NULL
    DefaultAssay(sub) <- "RNA"
    sub <- SCTransform(sub, assay = "RNA", verbose = FALSE)
    sub <- RunPCA(sub, verbose = FALSE)
    sub <- RunUMAP(sub,
        reduction = "pca", dims = UMAP_DIMS,
        n.neighbors = 20, min.dist = 0.1, spread = 4.0
    )

    root_cells <- colnames(sub)[sub$integrated_cluster == ROOT_CLUSTER]
    message("  cells = ", ncol(sub),
            " | cluster-", ROOT_CLUSTER, " root = ", length(root_cells))
    if (length(root_cells) < 50)
        warning("  root cluster < 50 cells for ", sample_id,
                " — pseudotime origin may be unstable", call. = FALSE)

    # --- Seurat -> cell_data_set ---
    cds <- as.cell_data_set(sub)
    rowData(cds)$gene_short_name <- rownames(cds)

    # monocle3 needs its own clusters/partitions; cluster on the per-patient
    # UMAP. use_partition = TRUE -> learn a separate graph per partition.
    # Forcing FALSE (one connected graph) would structurally assume EVERY
    # epithelial subtype originates from the cluster-8 diploid anchor — an
    # unjustified claim for a cohort of patient-private aneuploid subclones
    # (Hillock/Club are separate subclones). With TRUE, partitions not
    # connected to the cluster-8 root correctly get Inf pseudotime rather than
    # a forced ordering; downstream code already filters on is.finite().
    cds <- cluster_cells(cds, reduction_method = "UMAP")
    cds <- learn_graph(cds, use_partition = TRUE, close_loop = FALSE)
    cds <- order_cells(cds, root_cells = root_cells)

    sub$pseudotime <- pseudotime(cds)[colnames(sub)]

    # --- trajectory plots ---
    p_pseudo <- plot_cells(cds,
        color_cells_by = "pseudotime",
        label_cell_groups = FALSE, label_leaves = FALSE,
        label_branch_points = FALSE, cell_size = 0.5
    ) + ggtitle(paste0(sample_id, " — monocle3 pseudotime (root = cluster ",
                       ROOT_CLUSTER, ")"))
    ggsave(file.path(out, "UMAP_pseudotime.png"),
        plot = p_pseudo, width = 10, height = 8, bg = "white")

    # label_cell_groups = FALSE: integrated clusters are often spatially split
    # on the per-patient UMAP, so median-position labels land in empty space or
    # on top of other clusters. A discrete colour legend reads cleaner.
    p_clu <- plot_cells(cds,
        color_cells_by = "integrated_cluster",
        label_cell_groups = FALSE, label_leaves = FALSE,
        label_branch_points = TRUE, cell_size = 0.5
    ) + ggtitle(paste0(sample_id, " — graph by integrated cluster"))
    ggsave(file.path(out, "UMAP_graph_integrated_clusters.png"),
        plot = p_clu, width = 10, height = 8, bg = "white")

    p_part <- plot_cells(cds,
        color_cells_by = "partition",
        label_cell_groups = TRUE, label_leaves = FALSE,
        label_branch_points = FALSE, cell_size = 0.5
    ) + ggtitle(paste0(sample_id, " — monocle3 partitions"))
    ggsave(file.path(out, "UMAP_partitions.png"),
        plot = p_part, width = 10, height = 8, bg = "white")

    p_ck <- plot_cells(cds,
        color_cells_by = "copykat_clu8ref",
        label_cell_groups = FALSE, label_leaves = FALSE,
        label_branch_points = FALSE, cell_size = 0.5
    ) + scale_color_manual(values = copykat_cols) +
        ggtitle(paste0(sample_id, " — graph by copyKAT (cluster-8 ref)"))
    ggsave(file.path(out, "UMAP_graph_copykat.png"),
        plot = p_ck, width = 10, height = 8, bg = "white")

    # --- pseudotime distribution by integrated cluster ---
    p_box <- ggplot(sub@meta.data,
        aes(x = reorder(integrated_cluster, pseudotime, FUN = median),
            y = pseudotime, fill = integrated_cluster)
    ) +
        geom_boxplot(outlier.size = 0.3) +
        labs(x = "Integrated cluster (sorted by median pseudotime)",
             y = paste0("Pseudotime (root = cluster ", ROOT_CLUSTER, ")"),
             title = sample_id) +
        theme_classic() + NoLegend()
    ggsave(file.path(out, "Pseudotime_by_cluster.png"),
        plot = p_box, width = 10, height = 6, bg = "white")

    # --- marker dynamics along pseudotime ---
    markers <- intersect(traj_markers, rownames(cds))
    finite  <- !is.infinite(pseudotime(cds))
    if (length(markers) > 0 && sum(finite) > 0) {
        p_genes <- plot_genes_in_pseudotime(
            cds[markers, finite],
            color_cells_by = "integrated_cluster",
            min_expr = 0.1, ncol = 3
        ) + ggtitle(sample_id)
        ggsave(file.path(out, "Genes_in_pseudotime.png"),
            plot = p_genes, width = 14, height = 18, bg = "white")
    }

    # --- report ---
    tab <- sub@meta.data |>
        dplyr::group_by(integrated_cluster) |>
        dplyr::summarise(
            n          = dplyr::n(),
            median_pst = median(pseudotime[is.finite(pseudotime)], na.rm = TRUE),
            n_inf      = sum(!is.finite(pseudotime)),
            .groups    = "drop"
        ) |>
        dplyr::arrange(median_pst) |>
        as.data.frame()
    cat("\n----- ", sample_id,
        ": median pseudotime by integrated cluster -----\n", sep = "")
    print(tab)

    saveRDS(cds, file.path(out, paste0("cds_monocle3_", sample_id, ".rds")))
    saveRDS(sub, file.path(out, paste0("epi_", sample_id, "_pseudotime.rds")))

    # Per-cell pseudotime + trajectory-marker expression, for cross-patient
    # comparison in pseudotime space (see below). SCT `data` layer = log-norm.
    meta <- sub@meta.data |>
        dplyr::transmute(
            cell = colnames(sub),
            patient = sample_id,
            integrated_cluster,
            copykat_clu8ref,
            pseudotime
        )
    expr <- FetchData(sub, vars = markers, layer = "data")
    cbind(meta, expr[meta$cell, , drop = FALSE])
}

# bind_rows (not rbind): patients may retain slightly different marker sets
# after per-patient SCT gene filtering; missing columns are filled with NA.
pst_all <- dplyr::bind_rows(lapply(patients, run_patient))

# Per-patient min-max normalised pseudotime — each patient is rooted in its own
# cluster-8 cells on an independent scale, so [0,1] rescaling is required before
# any cross-patient overlay.
pst_all <- pst_all |>
    dplyr::group_by(patient) |>
    dplyr::mutate(pst_norm = {
        v <- ifelse(is.finite(pseudotime), pseudotime, NA_real_)
        (v - min(v, na.rm = TRUE)) /
            (max(v, na.rm = TRUE) - min(v, na.rm = TRUE))
    }) |>
    dplyr::ungroup()

# ============================================================
# Cross-patient comparison ----
# ============================================================
# Pseudotime is on an independent scale per patient (each rooted in its own
# cluster-8 cells), so this faceted view compares the *ordering* of clusters
# across patients, not absolute values.
write.csv(pst_all, file.path(OUT_DIR, "pseudotime_persample.csv"),
    row.names = FALSE)

p_cmp <- ggplot(
    dplyr::filter(pst_all, is.finite(pseudotime)),
    aes(x = integrated_cluster, y = pseudotime, fill = integrated_cluster)
) +
    geom_boxplot(outlier.size = 0.2) +
    facet_wrap(~patient, ncol = 1, scales = "free_y") +
    labs(x = "Integrated cluster",
         y = paste0("Pseudotime (root = cluster ", ROOT_CLUSTER, ", per-patient)")) +
    theme_classic() + NoLegend()
ggsave(file.path(OUT_DIR, "Pseudotime_by_cluster_allpatients.png"),
    plot = p_cmp, width = 10, height = 10, bg = "white")

# --- marker dynamics vs pseudotime, patients overlaid ---
# The core "do the 3 patients follow a similar path?" view: each trajectory
# marker as a function of normalised pseudotime, one loess curve per patient on
# a shared panel. Similar paths -> similar curve shapes (e.g. AR down, KRT7 up).
marker_long <- pst_all |>
    dplyr::filter(is.finite(pst_norm)) |>
    tidyr::pivot_longer(
        cols = dplyr::any_of(traj_markers),
        names_to = "gene", values_to = "expr"
    ) |>
    dplyr::filter(!is.na(expr)) |>
    dplyr::mutate(gene = factor(gene, levels = traj_markers))

p_marker <- ggplot(marker_long, aes(pst_norm, expr, colour = patient)) +
    geom_smooth(method = "loess", se = FALSE, span = 0.75, linewidth = 0.9) +
    facet_wrap(~gene, scales = "free_y", ncol = 3) +
    labs(x = "Pseudotime (per-patient min-max normalised)",
         y = "SCT log-normalised expression",
         colour = "Patient",
         title = "Trajectory-marker dynamics across patients") +
    theme_classic()
ggsave(file.path(OUT_DIR, "Markers_vs_pseudotime_allpatients.png"),
    plot = p_marker, width = 14, height = 18, bg = "white")

# --- integrated-cluster pseudotime ordering across patients (slope chart) ---
# Rank each integrated cluster by its median pseudotime within each patient.
# Flat lines across patients -> consistent cluster ordering (similar path);
# crossing lines -> the trajectory diverges between patients.
rank_df <- pst_all |>
    dplyr::filter(is.finite(pseudotime)) |>
    dplyr::group_by(patient, integrated_cluster) |>
    dplyr::summarise(median_pst = median(pseudotime), n = dplyr::n(),
                     .groups = "drop") |>
    dplyr::group_by(patient) |>
    dplyr::mutate(pst_rank = rank(median_pst)) |>
    dplyr::ungroup()

p_rank <- ggplot(rank_df, aes(patient, pst_rank,
                              group = integrated_cluster,
                              colour = integrated_cluster)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    geom_text(
        data = dplyr::filter(rank_df, patient == max(patient)),
        aes(label = integrated_cluster), hjust = -0.4, size = 3,
        show.legend = FALSE
    ) +
    scale_y_reverse(breaks = seq_len(max(rank_df$pst_rank))) +
    labs(x = "Patient", y = "Pseudotime rank (1 = earliest)",
         colour = "Integrated\ncluster",
         title = "Integrated-cluster pseudotime ordering across patients") +
    theme_classic()
ggsave(file.path(OUT_DIR, "Cluster_pseudotime_rank_allpatients.png"),
    plot = p_rank, width = 8, height = 7, bg = "white")

cat("\nSaved outputs to:", OUT_DIR, "\n")
