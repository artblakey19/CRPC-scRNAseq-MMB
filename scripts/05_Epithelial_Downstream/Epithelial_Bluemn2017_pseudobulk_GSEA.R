# Epithelial_Bluemn2017_pseudobulk_GSEA.R
# Bluemn et al. 2017 Cancer Cell Fig 4A/4B style panel, adapted to the current
# epithelial annotation.
#
# Current behavior:
#   - reference arm: the annotation label matching LE(ARPC), ARPC, or uniquely
#     containing ARPC.
#   - test arms: every other annotation cluster, independently vs ARPC.
#   - contrast-specific pseudobulk: RNA raw counts summed by patient x arm.
#   - DESeq2 paired design: ~ patient + group.
#   - fgsea ranking metric: DESeq2 Wald stat.

suppressMessages({
    library(Seurat)
    library(Matrix)
    library(DESeq2)
    library(fgsea)
    library(dplyr)
    library(ggplot2)
    library(patchwork)
})
suppressMessages(source("scripts/00_utils/scRNA_utils.R"))

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Bluemn2017_pseudobulk_GSEA"
GMT_DIR <- "Resources/msigdb_v6.0_GMTs"

TARGETS <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
    "MEK_UP.V1_UP",
    "ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_UP",
    "GO_REGULATION_OF_ERK1_AND_ERK2_CASCADE",
    "GO_POSITIVE_REGULATION_OF_MAPK_CASCADE",
    "HALLMARK_ANDROGEN_RESPONSE"
)
GMT_FILES <- c(
    "h.all.v6.0.symbols.gmt",
    "c2.cgp.v6.0.symbols.gmt",
    "c5.bp.v6.0.symbols.gmt",
    "c6.all.v6.0.symbols.gmt"
)
LABELS <- c(
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"        = "EMT (Hallmark)",
    "MEK_UP.V1_UP"                                      = "MEK activation (C6)",
    "ACEVEDO_FGFR1_TARGETS_IN_PROSTATE_CANCER_MODEL_UP" = "FGFR1 targets (Acevedo)",
    "GO_REGULATION_OF_ERK1_AND_ERK2_CASCADE"            = "ERK1/2 cascade (GO)",
    "GO_POSITIVE_REGULATION_OF_MAPK_CASCADE"            = "MAPK cascade+ (GO)",
    "HALLMARK_ANDROGEN_RESPONSE"                        = "Androgen response (Hallmark)"
)

REF_GROUP <- "ARPC"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

cleanup_patterns <- c(
    "DESeq2_*.csv",
    "NES_barplot_*.png",
    "NES_vertical_barplot_*.png",
    "NES_dotplot_*.png",
    "NES_heatmap_*.png",
    "RunningEnrichment_*.png",
    "contrast_definitions.csv",
    "fgsea_results.csv",
    "pseudobulk_counts*.rds",
    "pseudobulk_design*.csv"
)
old_outputs <- unlist(lapply(cleanup_patterns, function(x) {
    Sys.glob(file.path(OUT_DIR, x))
}), use.names = FALSE)
if (length(old_outputs) > 0) {
    unlink(old_outputs, force = TRUE)
}

# ============================================================
# Helpers ----
# ============================================================
read_gmt <- function(path) {
    lines <- readLines(path, warn = FALSE)
    sets <- lapply(lines, function(l) {
        f <- strsplit(l, "\t", fixed = TRUE)[[1]]
        list(name = f[1], genes = f[-(1:2)])
    })
    out <- lapply(sets, `[[`, "genes")
    names(out) <- vapply(sets, `[[`, character(1), "name")
    out
}

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

make_contrasts <- function(annotation_levels, ref_annotation) {
    test_annotations <- setdiff(annotation_levels, ref_annotation)
    keys <- vapply(test_annotations, sanitize_key, character(1))
    keys <- make.unique(keys, sep = "_")
    out <- lapply(test_annotations, function(ann) {
        list(label = paste(ann, "vs ARPC"), annotations = ann)
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

# ============================================================
# Gene sets from local v6.0 GMTs ----
# ============================================================
all_sets <- list()
for (gmt in GMT_FILES) {
    p <- file.path(GMT_DIR, gmt)
    if (!file.exists(p)) stop("Missing GMT: ", p)
    all_sets <- c(all_sets, read_gmt(p))
}
missing_sets <- setdiff(TARGETS, names(all_sets))
if (length(missing_sets) > 0) {
    stop("Sets not found in v6.0 GMTs: ", paste(missing_sets, collapse = ", "))
}
gene_sets <- all_sets[TARGETS]
for (nm in names(gene_sets)) {
    message(sprintf("  %s: %d genes (v6.0)", nm, length(gene_sets[[nm]])))
}

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
if (is.null(annotation_levels)) {
    annotation_levels <- unique(md$annotation)
}
annotation_levels <- annotation_levels[annotation_levels %in% unique(md$annotation)]
REF_ANNOTATION <- find_ref_annotation(annotation_levels)
CONTRASTS <- make_contrasts(annotation_levels, REF_ANNOTATION)
CONTRAST_KEYS <- names(CONTRASTS)
DISPLAY <- vapply(CONTRASTS, `[[`, character(1), "label")

message("Reference annotation: ", REF_ANNOTATION)
message("Annotation counts:")
print(sort(table(md$annotation), decreasing = TRUE))

write.csv(
    data.frame(
        contrast_key = CONTRAST_KEYS,
        contrast_label = unname(DISPLAY),
        reference_annotation = REF_ANNOTATION,
        test_annotations = vapply(
            CONTRASTS,
            function(x) paste(x$annotations, collapse = ";"),
            character(1)
        ),
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
# Per-contrast pseudobulk + DESeq2 + fgsea preranked ----
# ============================================================
make_pseudobulk <- function(contrast_key, spec) {
    test_annotations <- spec$annotations
    keep_cell <- md$annotation %in% c(REF_ANNOTATION, test_annotations)
    md_c <- md[keep_cell, , drop = FALSE]
    md_c$group <- ifelse(md_c$annotation == REF_ANNOTATION,
                         REF_GROUP, contrast_key)

    patient_group_tab <- table(md_c$patient, md_c$group)
    paired_patients <- rownames(patient_group_tab)[
        patient_group_tab[, REF_GROUP] > 0 &
            patient_group_tab[, contrast_key] > 0
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

run_contrast <- function(contrast_key, spec) {
    message("\nRunning contrast: ", spec$label)
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

    set.seed(1)
    fg <- fgsea(
        pathways = gene_sets,
        stats = ranks,
        eps = 0,
        minSize = 5,
        maxSize = Inf,
        BPPARAM = BiocParallel::SerialParam()
    )

    list(
        res = res_df[order(res_df$stat, decreasing = TRUE), ],
        fgsea = fg,
        ranks = ranks,
        pb = pb,
        coldata = cd
    )
}

results_list <- lapply(CONTRAST_KEYS, function(cn) run_contrast(cn, CONTRASTS[[cn]]))
names(results_list) <- CONTRAST_KEYS

pb_by_contrast <- lapply(results_list, `[[`, "pb")
coldata_by_contrast <- lapply(results_list, `[[`, "coldata")
saveRDS(
    list(pb_by_contrast = pb_by_contrast,
         coldata_by_contrast = coldata_by_contrast),
    file.path(OUT_DIR, "pseudobulk_counts_by_contrast.rds")
)
write.csv(
    bind_rows(coldata_by_contrast),
    file.path(OUT_DIR, "pseudobulk_design_by_contrast.csv"),
    row.names = FALSE
)

fg_all <- bind_rows(lapply(names(results_list), function(cn) {
    fg <- results_list[[cn]]$fgsea
    data.frame(
        contrast_key = cn,
        contrast_label = DISPLAY[[cn]],
        reference_annotation = REF_ANNOTATION,
        test_annotations = paste(CONTRASTS[[cn]]$annotations, collapse = ";"),
        pathway = fg$pathway,
        pathway_label = unname(LABELS[fg$pathway]),
        NES = fg$NES,
        ES = fg$ES,
        pval = fg$pval,
        padj = fg$padj,
        size = fg$size,
        leadingEdge = vapply(fg$leadingEdge, paste, character(1), collapse = ";")
    )
}))
write.csv(fg_all, file.path(OUT_DIR, "fgsea_results.csv"), row.names = FALSE)

for (cn in names(results_list)) {
    write.csv(
        results_list[[cn]]$res,
        file.path(OUT_DIR, paste0("DESeq2_", cn, ".csv")),
        row.names = FALSE
    )
}

# ============================================================
# Figure 4B style: NES bar plot ----
# ============================================================
fg_all$pathway_label <- factor(
    fg_all$pathway_label,
    levels = rev(unname(LABELS[TARGETS]))
)
fg_all$contrast_label <- factor(fg_all$contrast_label, levels = unname(DISPLAY))
fg_all$stars <- stars_of(fg_all$padj)
fg_all$hjust <- ifelse(fg_all$NES >= 0, -0.2, 1.2)
fg_all$star_label <- ifelse(fg_all$stars == "ns", "", fg_all$stars)
fg_all$contrast_short <- factor(
    sub(" vs ARPC$", "", as.character(fg_all$contrast_label)),
    levels = sub(" vs ARPC$", "", unname(DISPLAY))
)
fg_all$pathway_facet <- factor(
    as.character(fg_all$pathway_label),
    levels = unname(LABELS[TARGETS])
)
fg_all$NES_direction <- factor(
    ifelse(fg_all$NES >= 0, "Up vs ARPC", "Down vs ARPC"),
    levels = c("Up vs ARPC", "Down vs ARPC")
)

pal <- utils_cb_palette(length(CONTRAST_KEYS))
xpad <- max(abs(fg_all$NES), na.rm = TRUE) * 0.18
if (!is.finite(xpad) || xpad == 0) xpad <- 0.1
nes_lim <- max(abs(fg_all$NES), na.rm = TRUE)
if (!is.finite(nes_lim) || nes_lim == 0) nes_lim <- 1
nes_lim <- ceiling(nes_lim * 10) / 10

# Vertical bar plot layout: keep bars tall, not flat.
# Width scales with bars-per-facet (clusters); height scales with the number of
# facet rows so the NES axis stays tall relative to bar width even when there
# are many clusters. The old formula grew width ~linearly with cluster count but
# fixed the height small, so each facet was wide and short and bars looked squat.
VBAR_NCOL      <- 2
VBAR_PER_BAR   <- 0.95  # cm of facet width allotted per cluster bar
VBAR_FACET_H   <- 6.0   # cm of plotting height per facet row
vbar_nrow      <- ceiling(length(TARGETS) / VBAR_NCOL)
vbar_width_cm  <- VBAR_NCOL * (length(CONTRAST_KEYS) * VBAR_PER_BAR) + 4.5
vbar_height_cm <- vbar_nrow * VBAR_FACET_H + 6

p_bar <- ggplot(fg_all, aes(x = NES, y = pathway_label, fill = contrast_label)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_text(
        aes(label = stars, hjust = hjust),
        position = position_dodge(width = 0.85),
        size = 2.6
    ) +
    scale_fill_manual(values = pal, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    coord_cartesian(xlim = c(min(fg_all$NES) - xpad, max(fg_all$NES) + xpad)) +
    labs(
        x = "Normalized Enrichment Score (NES)\n(+ = up vs ARPC)",
        y = NULL,
        title = "Bluemn 2017 panel - each epithelial cluster vs ARPC",
        subtitle = paste0(
            "Reference annotation: ", REF_ANNOTATION,
            " | DESeq2 paired, ranked by Wald stat"
        )
    ) +
    theme_bw(base_size = 12) +
    theme(
        legend.position = "top",
        legend.text = element_text(size = 8),
        plot.subtitle = element_text(size = 8.5),
        panel.grid.major.y = element_blank()
    ) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))

ggsave(
    file.path(OUT_DIR, "NES_barplot_Fig4B_style_each_cluster.png"),
    plot = p_bar,
    width = 30,
    height = 20,
    units = "cm",
    dpi = 220,
    bg = "white",
    limitsize = FALSE
)

# Readable vertical summary: pathway facets ----
fg_all$star_y <- ifelse(
    fg_all$NES >= 0,
    fg_all$NES + nes_lim * 0.08,
    fg_all$NES - nes_lim * 0.08
)
fg_all$star_vjust <- ifelse(fg_all$NES >= 0, 0, 1)

p_vertical <- ggplot(fg_all, aes(x = contrast_short, y = NES, fill = NES_direction)) +
    geom_col(width = 0.72, colour = "grey35", linewidth = 0.15) +
    geom_hline(yintercept = 0, linewidth = 0.35, colour = "grey30") +
    geom_text(
        aes(y = star_y, label = star_label, vjust = star_vjust),
        size = 2.7,
        fontface = "bold"
    ) +
    facet_wrap(~ pathway_facet, ncol = VBAR_NCOL) +
    scale_fill_manual(
        values = c("Up vs ARPC" = "#B2182B", "Down vs ARPC" = "#2166AC"),
        name = NULL
    ) +
    scale_y_continuous(
        limits = c(-nes_lim * 1.18, nes_lim * 1.18),
        expand = expansion(mult = c(0.02, 0.02))
    ) +
    labs(
        x = "Epithelial cluster vs ARPC",
        y = "Normalized Enrichment Score (NES)",
        title = "Bluemn 2017 panel - vertical NES summary",
        subtitle = paste0(
            "Bars: NES (+ = up vs ARPC) | stars: FDR < 0.05 | Reference: ",
            REF_ANNOTATION
        )
    ) +
    theme_bw(base_size = 11) +
    theme(
        legend.position = "top",
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5),
        panel.grid.major.x = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill = "grey92", colour = "grey70"),
        strip.text = element_text(face = "bold", size = 9),
        plot.subtitle = element_text(size = 8.5)
    )

ggsave(
    file.path(OUT_DIR, "NES_vertical_barplot_Fig4B_style_each_cluster.png"),
    plot = p_vertical,
    width = vbar_width_cm,
    height = vbar_height_cm,
    units = "cm",
    dpi = 220,
    bg = "white",
    limitsize = FALSE
)

# ============================================================
# Figure 4A style: running enrichment plot ----
# ============================================================
HIGHLIGHT <- c(
    "HALLMARK_ANDROGEN_RESPONSE",
    "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)

enr_panels <- list()
for (cn in names(results_list)) {
    ranks <- results_list[[cn]]$ranks
    fg <- results_list[[cn]]$fgsea
    for (pw in HIGHLIGHT) {
        st <- fg[fg$pathway == pw, ]
        ttl <- sprintf(
            "%s\n%s | NES=%.2f, FDR=%.1e",
            DISPLAY[[cn]], LABELS[[pw]], st$NES, st$padj
        )
        enr_panels[[paste(cn, pw, sep = "__")]] <-
            plotEnrichment(gene_sets[[pw]], ranks) +
            labs(title = ttl) +
            theme(plot.title = element_text(size = 9))
    }
}

panel_keys <- as.vector(t(outer(CONTRAST_KEYS, HIGHLIGHT, paste, sep = "__")))
p_enr <- wrap_plots(enr_panels[panel_keys], ncol = length(HIGHLIGHT), byrow = TRUE) +
    plot_annotation(
        title = "Running enrichment - each epithelial cluster vs ARPC",
        subtitle = "Columns: Androgen response / EMT | ranked by DESeq2 Wald stat",
        theme = theme(plot.subtitle = element_text(size = 10))
    )
ggsave(
    file.path(OUT_DIR, "RunningEnrichment_Fig4A_style_each_cluster.png"),
    plot = p_enr,
    width = 26,
    height = max(24, 5.5 * length(CONTRAST_KEYS)),
    units = "cm",
    dpi = 200,
    bg = "white",
    limitsize = FALSE
)

# ============================================================
# Compact console summary ----
# ============================================================
message("\n=== fgsea NES (padj) - each cluster vs ARPC ===")
summ <- fg_all %>%
    mutate(cell = sprintf("%+.2f (%.1e)%s", NES, padj,
                          ifelse(stars == "ns", "", stars))) %>%
    select(pathway_label, contrast_label, cell) %>%
    tidyr::pivot_wider(names_from = contrast_label, values_from = cell)
print(as.data.frame(summ), row.names = FALSE)

message("\nDone. Outputs in: ", OUT_DIR)
