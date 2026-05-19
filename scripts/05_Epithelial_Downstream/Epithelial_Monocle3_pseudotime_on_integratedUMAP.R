# ============================================================
# Epithelial_Monocle3_pseudotime_on_integratedUMAP.R
#
# Overlay the PER-PATIENT monocle3 pseudotime onto the SHARED
# Harmony-integrated epithelial UMAP, faceted by sample.
#
# Rationale: Epithelial_Monocle3_persample.R learns each patient's trajectory
# independently on that patient's own unintegrated UMAP — so its native
# pseudotime plots live in 3 separate per-patient embeddings and can't be
# eyeballed side-by-side. Here we keep the pseudotime *values* from that run
# but paint them back onto the integrated UMAP coordinates (epi[["umap"]],
# built on Harmony), so all 3 patients are viewed in one shared embedding.
#
# Because each patient's pseudotime is rooted in its own cluster-8 cells on an
# independent scale, the cross-sample faceted overlay uses `pst_norm` (the
# per-patient min-max [0,1] rescaling from the persample script). Raw
# `pseudotime` is also plotted, but only as a per-facet free-scale view.
#
# Inputs:
#   Results/Epithelial/epi_clustered.rds
#       — integrated epithelial object; provides the Harmony-integrated `umap`
#         reduction + orig.ident. Only the embedding, barcodes and sample label
#         are used here — no copyKAT / pseudotime columns needed off this object.
#   Results/Epithelial/Monocle3_persample/pseudotime_persample.csv
#       — per-cell `pseudotime` / `pst_norm` from Epithelial_Monocle3_persample.R
# ============================================================

library(Seurat)
library(data.table)
library(dplyr)
library(ggplot2)

# ============================================================
# Config ----
# ============================================================
EPI_RDS  <- "Results/Epithelial/epi_clustered.rds"
M3_CSV   <- "Results/Epithelial/Monocle3_persample/pseudotime_persample.csv"
OUT_DIR  <- "Results/Epithelial/Monocle3_persample"   # alongside the persample run
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load + join ----
# ============================================================
epi <- readRDS(EPI_RDS)
stopifnot("umap" %in% Reductions(epi))

m3 <- data.table::fread(M3_CSV)[, .(cell, pseudotime, pst_norm)]

umap_xy <- as.data.frame(Embeddings(epi, "umap"))
colnames(umap_xy) <- c("UMAP_1", "UMAP_2")
umap_xy$cell    <- rownames(umap_xy)
umap_xy$patient <- as.character(epi$orig.ident)[match(umap_xy$cell, colnames(epi))]

df <- umap_xy |>
    dplyr::left_join(as.data.frame(m3), by = "cell")

# monocle3 leaves cells the principal graph can't reach at Inf pseudotime;
# treat those (and any cell missing from the CSV) as unassigned -> grey.
df$pseudotime[!is.finite(df$pseudotime)] <- NA_real_
df$pst_norm[!is.finite(df$pst_norm)]     <- NA_real_

n_assigned <- sum(!is.na(df$pst_norm))
message(sprintf("Cells: %d total | %d with finite pseudotime | %d unassigned",
                nrow(df), n_assigned, nrow(df) - n_assigned))

# ============================================================
# Plot helpers ----
# ============================================================
# Faint full-cohort backdrop reused in every facet, so each sample's cells are
# seen against the whole integrated epithelial UMAP, not floating in space.
backdrop <- df[, c("UMAP_1", "UMAP_2")]

base_umap <- function(point_size = 0.45) {
    list(
        geom_point(data = backdrop, aes(UMAP_1, UMAP_2),
                   colour = "grey88", size = point_size, inherit.aes = FALSE),
        theme_classic(),
        theme(strip.background = element_rect(fill = "grey95", colour = NA),
              strip.text = element_text(face = "bold")),
        coord_fixed()
    )
}

# ============================================================
# 1. Faceted overlay — per-patient min-max normalised pseudotime ----
#    Shared [0,1] colour scale: comparable across the 3 samples.
# ============================================================
p_facet_norm <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    base_umap() +
    geom_point(aes(colour = pst_norm), size = 0.45) +
    scale_colour_viridis_c(option = "plasma", na.value = "grey88",
                           name = "Pseudotime\n(per-patient\nmin-max norm.)") +
    facet_wrap(~patient) +
    labs(title = "Per-patient monocle3 pseudotime on the integrated epithelial UMAP",
         subtitle = "Normalised [0,1] within each patient — comparable across samples",
         x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "IntegratedUMAP_pseudotime_norm_facet.png"),
       plot = p_facet_norm, width = 16, height = 6, bg = "white")

# ============================================================
# 2. Faceted overlay — raw pseudotime, per-facet free colour scale ----
#    Raw values live on independent per-patient scales, so this is a
#    within-sample view only (one ggplot per patient, then patchwork).
# ============================================================
patients <- sort(unique(df$patient))
raw_panels <- lapply(patients, function(p) {
    d <- df[df$patient == p, ]
    ggplot(d, aes(UMAP_1, UMAP_2)) +
        base_umap() +
        geom_point(aes(colour = pseudotime), size = 0.45) +
        scale_colour_viridis_c(option = "viridis", na.value = "grey88",
                               name = "Pseudotime") +
        labs(title = p, x = "UMAP 1", y = "UMAP 2") +
        theme_classic() +
        theme(plot.title = element_text(face = "bold"))
})
p_facet_raw <- patchwork::wrap_plots(raw_panels, nrow = 1)
ggsave(file.path(OUT_DIR, "IntegratedUMAP_pseudotime_raw_facet.png"),
       plot = p_facet_raw, width = 16, height = 6, bg = "white")

# ============================================================
# 3. Single combined UMAP — all cells, normalised pseudotime ----
# ============================================================
p_combined <- ggplot(df, aes(UMAP_1, UMAP_2)) +
    geom_point(aes(colour = pst_norm), size = 0.4) +
    scale_colour_viridis_c(option = "plasma", na.value = "grey88",
                           name = "Pseudotime\n(per-patient\nmin-max norm.)") +
    theme_classic() + coord_fixed() +
    labs(title = "Per-patient monocle3 pseudotime on the integrated epithelial UMAP",
         subtitle = "All samples overlaid",
         x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(OUT_DIR, "IntegratedUMAP_pseudotime_norm_combined.png"),
       plot = p_combined, width = 9, height = 8, bg = "white")

# ============================================================
# Persist the joined table ----
# ============================================================
write.csv(df, file.path(OUT_DIR, "pseudotime_on_integratedUMAP.csv"),
          row.names = FALSE)

cat("\nSaved outputs to:", OUT_DIR, "\n")
cat("  IntegratedUMAP_pseudotime_norm_facet.png    — faceted, shared [0,1] scale\n")
cat("  IntegratedUMAP_pseudotime_raw_facet.png     — faceted, per-facet raw scale\n")
cat("  IntegratedUMAP_pseudotime_norm_combined.png — single overlaid UMAP\n")
cat("  pseudotime_on_integratedUMAP.csv            — joined coords + pseudotime\n")
