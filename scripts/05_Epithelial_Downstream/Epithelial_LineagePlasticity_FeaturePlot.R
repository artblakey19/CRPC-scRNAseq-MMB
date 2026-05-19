# Epithelial_LineagePlasticity_FeaturePlot.R
# FeaturePlot of lineage-plasticity / NEPC / tumor-suppressor genes on the
# filtered epithelial UMAP.
#
# Caveats (do not over-interpret):
#   - 10x 3' chemistry: SNV/indel calls are unreliable -> "expressed" does NOT
#     rule out LOF in PTEN/TP53/RB1/KMT2C/CHD7 etc. Cross-check with copyKAT
#     CNV (chr10 for PTEN, chr17p for TP53, chr13q for RB1).
#   - dNp63: isoform-specific; 3' tag reads cannot separate TAp63 vs dNp63.
#     Plotted as TP63 only.
#   - "FOX2" in the request is interpreted as FOXA2 (NEPC pioneer TF). FOXA1
#     also plotted for ARPC contrast.
#
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds
# Writes: Results/05_Epithelial_Downstream/Annotation/
#           FeaturePlot_LineagePlasticity.png
#           FeaturePlot_LineagePlasticity_byPatient.png
#           VlnPlot_LineagePlasticity.png
#           VlnPlot_LineagePlasticity_byPatient.png
#           LineagePlasticity_gene_status.tsv

suppressMessages({
    library(Seurat)
    library(ggplot2)
    library(patchwork)
})

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Requested symbols -> HGNC-resolved symbols used for plotting
gene_map <- c(
    "CHD7"   = "CHD7",
    "MYC"    = "MYC",
    "KMT2C"  = "KMT2C",
    "KRT7"   = "KRT7",
    "SOX2"   = "SOX2",
    "FOXA1"  = "FOXA1",   # ARPC pioneer TF (added for contrast)
    "SYP"    = "SYP",
    "AR"     = "AR",
    "TP53"   = "TP53",
    "RB1"    = "RB1",
    "PTEN"   = "PTEN",
    "REST"   = "REST",
    "YAP1"   = "YAP1",
    "TP63"   = "TP63"     # dNp63 isoform not resolvable; gene-level only
)

epi <- readRDS(IN_RDS)
DefaultAssay(epi) <- "SCT"

present <- gene_map %in% rownames(epi)
status  <- data.frame(
    requested = names(gene_map),
    resolved  = unname(gene_map),
    in_object = present,
    stringsAsFactors = FALSE
)
write.table(status,
            file = file.path(OUT_DIR, "LineagePlasticity_gene_status.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
if (any(!present)) {
    message("Not in object: ",
            paste(gene_map[!present], collapse = ", "))
}
features <- gene_map[present]

fp_list <- lapply(features, function(g) {
    FeaturePlot(epi, features = g, order = TRUE, pt.size = 0.25) +
        ggtitle(g) +
        theme(plot.title = element_text(face = "bold"))
})
fp <- wrap_plots(fp_list, ncol = 4)
ggsave(file.path(OUT_DIR, "FeaturePlot_LineagePlasticity.png"),
       plot = fp, width = 20, height = 18, bg = "white", dpi = 200)

group_var <- if ("annotation" %in% colnames(epi@meta.data)) "annotation" else "seurat_clusters"
vp <- VlnPlot(epi, features = features, group.by = group_var,
              pt.size = 0, stack = TRUE, flip = TRUE) +
    NoLegend() +
    ggtitle(paste0("Lineage-plasticity / NEPC / TSG markers by ", group_var))
ggsave(file.path(OUT_DIR, "VlnPlot_LineagePlasticity.png"),
       plot = vp, width = 10, height = 11, bg = "white", dpi = 200)

message("Wrote: ", file.path(OUT_DIR, "FeaturePlot_LineagePlasticity.png"))
message("Wrote: ", file.path(OUT_DIR, "VlnPlot_LineagePlasticity.png"))
message("Wrote: ", file.path(OUT_DIR, "LineagePlasticity_gene_status.tsv"))

# ---- Per-patient split ----------------------------------------------------
# Layout: rows = genes, cols = patients. To make patients visually comparable
# we fix the colorbar ceiling per gene to the GLOBAL q99 expression across
# all epithelial cells (not per-panel). Caveat: patient == batch in this
# cohort, so absolute levels are not quantitatively comparable; use these
# plots for pattern/where-is-it-expressed comparison, not magnitude.
patients <- sort(unique(as.character(epi$orig.ident)))

expr_mat <- GetAssayData(epi, assay = "SCT", layer = "data")[features, , drop = FALSE]
gene_q99 <- apply(expr_mat, 1, function(x) {
    v <- as.numeric(x)
    q <- stats::quantile(v[v > 0], 0.99, na.rm = TRUE)
    if (!is.finite(q) || q <= 0) max(v, na.rm = TRUE) else q
})

fp_split_rows <- lapply(features, function(g) {
    cap <- gene_q99[[g]]
    row_plots <- lapply(patients, function(p) {
        cells_p <- colnames(epi)[epi$orig.ident == p]
        FeaturePlot(epi, features = g, cells = cells_p,
                    order = TRUE, pt.size = 0.25,
                    min.cutoff = 0, max.cutoff = cap) +
            ggtitle(paste0(g, " | ", p)) +
            theme(plot.title = element_text(face = "bold", size = 10),
                  legend.position = "right")
    })
    wrap_plots(row_plots, nrow = 1)
})
fp_split <- wrap_plots(fp_split_rows, ncol = 1)
ggsave(file.path(OUT_DIR, "FeaturePlot_LineagePlasticity_byPatient.png"),
       plot = fp_split,
       width = 5 * length(patients),
       height = 4 * length(features),
       bg = "white", dpi = 180, limitsize = FALSE)

vp_split <- VlnPlot(epi, features = features, group.by = group_var,
                    split.by = "orig.ident",
                    pt.size = 0, stack = TRUE, flip = TRUE) +
    ggtitle(paste0("Markers by ", group_var, ", split by patient"))
ggsave(file.path(OUT_DIR, "VlnPlot_LineagePlasticity_byPatient.png"),
       plot = vp_split, width = 12, height = 12, bg = "white", dpi = 200)

message("Wrote: ", file.path(OUT_DIR, "FeaturePlot_LineagePlasticity_byPatient.png"))
message("Wrote: ", file.path(OUT_DIR, "VlnPlot_LineagePlasticity_byPatient.png"))
