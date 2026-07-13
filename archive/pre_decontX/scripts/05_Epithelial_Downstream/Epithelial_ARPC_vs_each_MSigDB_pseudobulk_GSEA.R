# Epithelial_ARPC_vs_each_MSigDB_pseudobulk_GSEA.R
# Discovery-style pseudobulk GSEA of LE(ARPC) vs every other epithelial cluster,
# run against whole MSigDB collections (not a curated target list):
#   - Hallmark (H)
#   - WikiPathways (C2 CP:WIKIPATHWAYS)
#   - C6 oncogenic signatures (C6)
#   - C5 GO Biological Process (C5 GO:BP)
#
# Behavior:
#   - reference arm: the annotation label matching LE(ARPC), ARPC, or uniquely
#     containing ARPC.
#   - test arms: every other annotation cluster, independently vs ARPC.
#   - contrast-specific pseudobulk: RNA raw counts summed by patient x arm.
#   - DESeq2 paired design: ~ patient + group; ranking metric = Wald stat.
#   - fgsea run separately per collection so FDR is collection-specific
#     (GO:BP's thousands of sets do not dilute Hallmark's FDR).
#   - gene sets pulled live from the installed msigdbr (MSigDB 2026.1.Hs) so
#     WikiPathways is available (absent from the local v6.0 GMTs).

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(DESeq2)
    library(fgsea)
    library(msigdbr)
    library(dplyr)
    library(ggplot2)
})
suppressMessages(source("scripts/00_utils/scRNA_utils.R"))

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/ARPC_vs_each_MSigDB_pseudobulk_GSEA"

REF_GROUP <- "ARPC"   # internal token for the reference arm (DESeq2-safe name)

# MSigDB collections to sweep (msigdbr new API: collection / subcollection).
COLLECTIONS <- list(
    Hallmark     = list(collection = "H",  subcollection = NULL,              label = "Hallmark"),
    WikiPathways = list(collection = "C2", subcollection = "CP:WIKIPATHWAYS", label = "WikiPathways"),
    C6_Oncogenic = list(collection = "C6", subcollection = NULL,              label = "C6 oncogenic"),
    GOBP         = list(collection = "C5", subcollection = "GO:BP",           label = "C5 GO:BP")
)

# fgsea gene-set size window (standard for whole-collection discovery GSEA).
MIN_SIZE <- 15
MAX_SIZE <- 500

# Dotplot: how many top pathways (by |NES| among significant) per cluster to
# take into the cross-cluster summary for each collection.
TOP_PER_CLUSTER <- 6
PADJ_CUTOFF     <- 0.05

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

cleanup_patterns <- c(
    "DESeq2_*.csv",
    "fgsea_*.csv",
    "NES_dotplot_*.png",
    "contrast_definitions.csv",
    "pseudobulk_counts*.rds",
    "pseudobulk_design*.csv"
)
old_outputs <- unlist(lapply(cleanup_patterns, function(x) {
    Sys.glob(file.path(OUT_DIR, x))
}), use.names = FALSE)
if (length(old_outputs) > 0) unlink(old_outputs, force = TRUE)

# ============================================================
# Helpers ----
# ============================================================
sanitize_key <- function(label) {
    key <- gsub("[^A-Za-z0-9]+", "", label)
    if (!grepl("^[A-Za-z]", key)) key <- paste0("Cluster", key)
    paste0(key, "_vs_ARPC")
}

find_ref_annotation <- function(annotation_levels) {
    exact <- intersect(c("LE(ARPC)", "ARPC"), annotation_levels)
    if (length(exact) > 0) return(exact[[1]])
    hits <- annotation_levels[grepl("ARPC", annotation_levels, ignore.case = TRUE)]
    if (length(hits) != 1) {
        stop("Could not uniquely identify ARPC annotation. Candidates: ",
             paste(hits, collapse = ", "))
    }
    hits[[1]]
}

make_contrasts <- function(annotation_levels, ref_annotation, ref_display) {
    test_annotations <- setdiff(annotation_levels, ref_annotation)
    keys <- make.unique(vapply(test_annotations, sanitize_key, character(1)), sep = "_")
    out <- lapply(test_annotations, function(ann) {
        list(label = paste(ann, "vs", ref_display), annotations = ann, short = ann)
    })
    names(out) <- keys
    out
}

stars_of <- function(p) {
    ifelse(is.na(p), "ns",
           ifelse(p < 1e-3, "***",
                  ifelse(p < 1e-2, "**",
                         ifelse(p < 0.05, "*", "ns"))))
}

# Human-readable pathway label: drop collection prefix, underscores -> spaces.
pretty_pathway <- function(x, width = 55) {
    x2 <- sub("^(HALLMARK|GOBP|GOCC|GOMF|WP|REACTOME|KEGG|BIOCARTA|PID)_", "", x)
    x2 <- gsub("_", " ", x2)
    ifelse(nchar(x2) > width, paste0(substr(x2, 1, width - 1), "…"), x2)
}

# ============================================================
# Gene sets from msigdbr ----
# ============================================================
message("Pulling MSigDB collections via msigdbr (", as.character(packageVersion("msigdbr")), ")")
collection_sets <- lapply(COLLECTIONS, function(cfg) {
    df <- msigdbr(species = "human",
                  collection = cfg$collection,
                  subcollection = cfg$subcollection)
    sets <- split(df$gene_symbol, df$gs_name)
    sets <- lapply(sets, function(g) unique(g))
    message(sprintf("  %-13s %s: %d sets", cfg$label, cfg$collection,
                    length(sets)))
    sets
})
names(collection_sets) <- names(COLLECTIONS)

# ============================================================
# Load annotated object and define current contrasts ----
# ============================================================
epi <- readRDS(IN_RDS)
DefaultAssay(epi) <- "RNA"
message("Loaded: ", IN_RDS, " | cells=", ncol(epi), " genes=", nrow(epi))

md <- epi@meta.data
md$annotation <- as.character(md$annotation)
md$patient <- as.character(md$orig.ident)

annotation_levels <- levels(epi$annotation)
if (is.null(annotation_levels)) annotation_levels <- unique(md$annotation)
annotation_levels <- annotation_levels[annotation_levels %in% unique(md$annotation)]

REF_ANNOTATION <- find_ref_annotation(annotation_levels)
REF_DISPLAY    <- REF_ANNOTATION
CONTRASTS      <- make_contrasts(annotation_levels, REF_ANNOTATION, REF_DISPLAY)
CONTRAST_KEYS  <- names(CONTRASTS)
DISPLAY        <- vapply(CONTRASTS, `[[`, character(1), "label")
SHORT          <- vapply(CONTRASTS, `[[`, character(1), "short")

message("Reference annotation: ", REF_ANNOTATION)
message("Test clusters (", length(CONTRAST_KEYS), "): ",
        paste(SHORT, collapse = ", "))
message("Annotation counts:")
print(sort(table(md$annotation), decreasing = TRUE))

write.csv(
    data.frame(
        contrast_key = CONTRAST_KEYS,
        contrast_label = unname(DISPLAY),
        reference_annotation = REF_ANNOTATION,
        test_annotations = unname(SHORT),
        stringsAsFactors = FALSE
    ),
    file.path(OUT_DIR, "contrast_definitions.csv"),
    row.names = FALSE
)

counts <- tryCatch(
    GetAssayData(epi, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(epi, assay = "RNA", slot = "counts")
)

# ============================================================
# Per-contrast pseudobulk + DESeq2 (ranking shared across collections) ----
# ============================================================
make_pseudobulk <- function(contrast_key, spec) {
    test_annotations <- spec$annotations
    keep_cell <- md$annotation %in% c(REF_ANNOTATION, test_annotations)
    md_c <- md[keep_cell, , drop = FALSE]
    md_c$group <- ifelse(md_c$annotation == REF_ANNOTATION, REF_GROUP, contrast_key)

    patient_group_tab <- table(md_c$patient, md_c$group)
    paired_patients <- rownames(patient_group_tab)[
        patient_group_tab[, REF_GROUP] > 0 & patient_group_tab[, contrast_key] > 0
    ]
    if (length(paired_patients) < 2) {
        stop("Contrast ", contrast_key,
             " has fewer than 2 paired patients with both ARPC and target cells.")
    }
    dropped <- setdiff(unique(md_c$patient), paired_patients)
    if (length(dropped) > 0) {
        message("  ", contrast_key, ": dropping unpaired patient(s): ",
                paste(dropped, collapse = ", "))
    }
    md_c <- md_c[md_c$patient %in% paired_patients, , drop = FALSE]

    pb_key <- paste(md_c$patient, md_c$group, sep = "__")
    ind <- Matrix::sparse.model.matrix(~ 0 + factor(pb_key))
    colnames(ind) <- levels(factor(pb_key))

    pb <- as.matrix(counts[, rownames(md_c), drop = FALSE] %*% ind)
    storage.mode(pb) <- "integer"

    coldata <- data.frame(
        contrast = contrast_key,
        contrast_label = spec$label,
        reference_annotation = REF_ANNOTATION,
        test_annotations = paste(test_annotations, collapse = ";"),
        pb_key = colnames(pb),
        patient = sub("__.*$", "", colnames(pb)),
        group = sub("^.*__", "", colnames(pb)),
        stringsAsFactors = FALSE
    )
    coldata$n_cells <- as.integer(table(pb_key)[coldata$pb_key])
    list(pb = pb, coldata = coldata)
}

run_deseq <- function(contrast_key, spec) {
    message("\nDESeq2 contrast: ", spec$label)
    pb_obj <- make_pseudobulk(contrast_key, spec)
    pb <- pb_obj$pb
    cd <- pb_obj$coldata
    message("Pseudobulk profiles:")
    print(cd[, c("patient", "group", "n_cells")])

    cd$group <- factor(cd$group, levels = c(REF_GROUP, contrast_key))
    cd$patient <- droplevels(factor(cd$patient))

    m <- pb[rowSums(pb) >= 10, , drop = FALSE]
    dds <- DESeqDataSetFromMatrix(m, colData = cd, design = ~ patient + group)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("group", contrast_key, REF_GROUP))
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$contrast <- contrast_key
    res_df$contrast_label <- spec$label
    res_df$reference_annotation <- REF_ANNOTATION
    res_df$test_annotations <- paste(spec$annotations, collapse = ";")

    ranks <- res_df$stat
    names(ranks) <- res_df$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

    list(
        res = res_df[order(res_df$stat, decreasing = TRUE), ],
        ranks = ranks,
        pb = pb,
        coldata = cd
    )
}

deseq_list <- lapply(CONTRAST_KEYS, function(cn) run_deseq(cn, CONTRASTS[[cn]]))
names(deseq_list) <- CONTRAST_KEYS

# Save pseudobulk + DESeq2 tables (DESeq2 result shared by all collections).
pb_by_contrast <- lapply(deseq_list, `[[`, "pb")
coldata_by_contrast <- lapply(deseq_list, `[[`, "coldata")
saveRDS(
    list(pb_by_contrast = pb_by_contrast, coldata_by_contrast = coldata_by_contrast),
    file.path(OUT_DIR, "pseudobulk_counts_by_contrast.rds")
)
write.csv(bind_rows(coldata_by_contrast),
          file.path(OUT_DIR, "pseudobulk_design_by_contrast.csv"),
          row.names = FALSE)
for (cn in names(deseq_list)) {
    write.csv(deseq_list[[cn]]$res,
              file.path(OUT_DIR, paste0("DESeq2_", cn, ".csv")),
              row.names = FALSE)
}

# ============================================================
# fgsea per (contrast x collection) ----
# ============================================================
run_fgsea <- function(ranks, gene_sets) {
    set.seed(1)
    fgsea(
        pathways = gene_sets,
        stats = ranks,
        eps = 0,
        minSize = MIN_SIZE,
        maxSize = MAX_SIZE,
        BPPARAM = BiocParallel::SerialParam()
    )
}

fg_records <- list()
for (coll_name in names(collection_sets)) {
    gene_sets <- collection_sets[[coll_name]]
    coll_label <- COLLECTIONS[[coll_name]]$label
    message("\n=== Collection: ", coll_label, " (", length(gene_sets), " sets) ===")
    coll_rows <- list()
    for (cn in CONTRAST_KEYS) {
        message("  fgsea ", coll_label, " | ", DISPLAY[[cn]])
        fg <- run_fgsea(deseq_list[[cn]]$ranks, gene_sets)
        fg <- fg[order(fg$padj, -abs(fg$NES)), ]
        rec <- data.frame(
            collection = coll_name,
            collection_label = coll_label,
            contrast_key = cn,
            contrast_label = DISPLAY[[cn]],
            cluster = SHORT[[cn]],
            reference_annotation = REF_ANNOTATION,
            pathway = fg$pathway,
            pathway_label = pretty_pathway(fg$pathway),
            NES = fg$NES,
            ES = fg$ES,
            pval = fg$pval,
            padj = fg$padj,
            size = fg$size,
            leadingEdge = vapply(fg$leadingEdge, paste, character(1), collapse = ";"),
            stringsAsFactors = FALSE
        )
        coll_rows[[cn]] <- rec
    }
    coll_df <- bind_rows(coll_rows)
    write.csv(coll_df,
              file.path(OUT_DIR, paste0("fgsea_", coll_name, ".csv")),
              row.names = FALSE)
    fg_records[[coll_name]] <- coll_df
}
fg_all <- bind_rows(fg_records)
write.csv(fg_all, file.path(OUT_DIR, "fgsea_results_all.csv"), row.names = FALSE)

# ============================================================
# Cross-cluster NES dotplot per collection ----
# ============================================================
cluster_order <- SHORT  # annotation-level order (ref excluded)

make_dotplot <- function(coll_name) {
    coll_label <- COLLECTIONS[[coll_name]]$label
    df <- fg_records[[coll_name]]

    sig <- df[!is.na(df$padj) & df$padj < PADJ_CUTOFF, ]
    if (nrow(sig) == 0) {
        message("  [", coll_label, "] no pathway with padj < ", PADJ_CUTOFF,
                " in any contrast; skipping dotplot.")
        return(invisible(NULL))
    }
    # union of top-N (by |NES|) significant pathways per cluster
    top_paths <- sig %>%
        group_by(cluster) %>%
        arrange(desc(abs(NES)), .by_group = TRUE) %>%
        slice_head(n = TOP_PER_CLUSTER) %>%
        ungroup() %>%
        pull(pathway) %>%
        unique()

    plot_df <- df[df$pathway %in% top_paths, , drop = FALSE]
    plot_df$neglog10padj <- -log10(pmax(plot_df$padj, .Machine$double.xmin))
    plot_df$is_sig <- !is.na(plot_df$padj) & plot_df$padj < PADJ_CUTOFF
    plot_df$cluster <- factor(plot_df$cluster,
                              levels = cluster_order[cluster_order %in% plot_df$cluster])

    # order pathways by mean NES across clusters (top = most up vs ARPC)
    path_ord <- plot_df %>%
        group_by(pathway, pathway_label) %>%
        summarise(mean_nes = mean(NES, na.rm = TRUE), .groups = "drop") %>%
        arrange(mean_nes)
    plot_df$pathway_label <- factor(plot_df$pathway_label,
                                    levels = path_ord$pathway_label)

    nes_lim <- max(abs(plot_df$NES), na.rm = TRUE)
    if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1

    p <- ggplot(plot_df, aes(x = cluster, y = pathway_label)) +
        geom_point(aes(size = neglog10padj, fill = NES, colour = is_sig),
                   shape = 21, stroke = 0.6) +
        scale_fill_gradient2(
            low = "#2166AC", mid = "grey92", high = "#B2182B",
            midpoint = 0, limits = c(-nes_lim, nes_lim),
            name = "NES\n(+ up vs ARPC)"
        ) +
        scale_colour_manual(
            values = c("TRUE" = "black", "FALSE" = "grey75"),
            labels = c("TRUE" = "FDR < 0.05", "FALSE" = "n.s."),
            name = NULL
        ) +
        scale_size_continuous(range = c(1, 5), name = "-log10(FDR)") +
        # pad the axes so edge dots are not clipped and rows/cols get breathing room
        scale_x_discrete(expand = expansion(add = 0.7)) +
        scale_y_discrete(expand = expansion(add = 0.7)) +
        labs(
            x = paste0("Epithelial cluster vs ", REF_DISPLAY),
            y = NULL,
            title = paste0(coll_label, ": each cluster vs ", REF_DISPLAY),
            subtitle = paste0(
                "Pseudobulk DESeq2 (paired ~patient+group), fgsea ranked by Wald stat | ",
                "union of top ", TOP_PER_CLUSTER, " FDR<0.05 pathways / cluster"
            )
        ) +
        theme_bw(base_size = 11) +
        theme(
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
            axis.text.y = element_text(size = 8),
            panel.grid.minor = element_blank(),
            plot.subtitle = element_text(size = 8),
            legend.key.size = unit(0.4, "cm")
        )

    n_paths <- nlevels(plot_df$pathway_label)
    n_clust <- length(levels(plot_df$cluster))
    # ~0.75 cm per pathway row and ~1.6 cm per cluster column keeps the largest
    # dot (size range max 5) clear of its neighbours; generous base offsets leave
    # room for the long y-axis pathway labels and the right-hand legend.
    h <- max(12, 3 + n_paths * 0.75)
    w <- max(22, 14 + n_clust * 1.6)
    ggsave(
        file.path(OUT_DIR, paste0("NES_dotplot_", coll_name, ".png")),
        plot = p, width = w, height = h, units = "cm",
        dpi = 220, bg = "white", limitsize = FALSE
    )
    message("  [", coll_label, "] dotplot: ", n_paths, " pathways x ",
            length(levels(plot_df$cluster)), " clusters")
}

for (coll_name in names(collection_sets)) make_dotplot(coll_name)

# ============================================================
# Console summary: significant-pathway counts per contrast x collection ----
# ============================================================
message("\n=== # pathways FDR < ", PADJ_CUTOFF, " (up / down vs ARPC) ===")
summ <- fg_all %>%
    filter(!is.na(padj) & padj < PADJ_CUTOFF) %>%
    mutate(dir = ifelse(NES >= 0, "up", "down")) %>%
    count(collection_label, cluster, dir) %>%
    tidyr::pivot_wider(names_from = dir, values_from = n, values_fill = 0)
if (!"up" %in% names(summ))   summ$up <- 0L
if (!"down" %in% names(summ)) summ$down <- 0L
summ <- summ %>%
    mutate(cell = paste0(up, " / ", down)) %>%
    select(collection_label, cluster, cell) %>%
    tidyr::pivot_wider(names_from = collection_label, values_from = cell,
                       values_fill = "0 / 0")
print(as.data.frame(summ), row.names = FALSE)

message("\nDone. Outputs in: ", OUT_DIR)
