# Epithelial_Song2022_PaperPCa_ModuleScore.R
# Song et al. 2022 Nat Comm 재현 — franklinhuanglab repo Seqwell_Combined_E.R 패턴.
# Song 2022가 tumor cluster annotation용으로 사용한 paper-named PCa signatures
# (Liu/Tomlins/Wallace/Setlur 4종, C2CGP) — Seurat AddModuleScore로 스코어링.
#
# PCa-known gene set들의 epithelial UMAP 상 activity 표시 (참고용).
#
# DNPC cohort caveat:
#   - AR-driven LIU_PROSTATE_CANCER_UP, SETLUR_TMPRSS2_ERG_FUSION_UP은 muted 예상.
#   - WALLACE_PROSTATE_CANCER_UP / TOMLINS_PROSTATE_CANCER_UP이 ARPC
#     (cluster 9)에서 가장 강하게 잡힐 것.

library(Seurat)
library(dplyr)
library(msigdbr)
library(ggplot2)
library(patchwork)   # plot_annotation / `&` for split.by FeaturePlot panels
source("scripts/00_utils/scRNA_utils.R")

IN_RDS  <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Song2022_PaperPCa_ModuleScore"
OUT_RDS <- "Results/05_Epithelial_Downstream/epi_song2022_paperpca_modulescore.rds"

# Paper named PCa gene sets (Seqwell_Combined_E.R line 896-906에서 사용한 4종)
PAPER_PCA_SETS <- c(
    "LIU_PROSTATE_CANCER_UP",
    "TOMLINS_PROSTATE_CANCER_UP",
    "WALLACE_PROSTATE_CANCER_UP",
    "SETLUR_PROSTATE_CANCER_TMPRSS2_ERG_FUSION_UP"
)

stopifnot("Input rds not found" = file.exists(IN_RDS))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Load-or-compute: skip heavy module scoring when OUT_RDS exists ----
# ============================================================
if (file.exists(OUT_RDS)) {
    message("Loading cached ", OUT_RDS, " — skipping module scoring, regenerating figures")
    epi <- readRDS(OUT_RDS)
} else {

epi <- readRDS(IN_RDS)
message("Loaded: ", IN_RDS)

DefaultAssay(epi) <- "RNA"
if (!"data" %in% Layers(epi[["RNA"]])) {
    epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)
}

# ============================================================
# PCa gene sets from C2CGP ----
# ============================================================
c2cgp_df <- msigdbr(
    species = "Homo sapiens",
    collection = "C2", subcollection = "CGP"
)
c2cgp <- split(c2cgp_df$gene_symbol, c2cgp_df$gs_name)

missing_named <- setdiff(PAPER_PCA_SETS, names(c2cgp))
if (length(missing_named) > 0) {
    warning("Paper PCa sets missing in msigdbr: ", paste(missing_named, collapse = ", "))
}
pca_sets <- c2cgp[intersect(PAPER_PCA_SETS, names(c2cgp))]
message("Paper PCa sets matched: ", length(pca_sets), " / ", length(PAPER_PCA_SETS))

# ============================================================
# AddModuleScore on PCa gene sets ----
# ============================================================
# 저자 repo Seqwell_Combined_E.R:600 패턴 — assay="RNA" + ctrl=default(100).
for (set_name in names(pca_sets)) {
    epi <- tryCatch(
        AddModuleScore(epi,
            features = list(pca_sets[[set_name]]),
            name = paste0(set_name, "_sig"),
            assay = "RNA"
        ),
        error = function(e) {
            AddModuleScore(epi,
                features = list(pca_sets[[set_name]]),
                name = paste0(set_name, "_sig"),
                nbin = 12,
                assay = "RNA"
            )
        }
    )
    src <- paste0(set_name, "_sig1")
    dst <- paste0(set_name, "_sig")
    epi@meta.data[[dst]] <- epi@meta.data[[src]]
    epi@meta.data[[src]] <- NULL
}

# Persist annotated object (compute branch only)
saveRDS(epi, OUT_RDS)

}  # end load-or-compute

# Signature set names usable by the figure loop in both modes.
# In load mode pca_sets does not exist — derive from the *_sig meta.data columns.
set_names <- if (exists("pca_sets")) names(pca_sets) else
    sub("_sig$", "", grep("_sig$", colnames(epi@meta.data), value = TRUE))

# ============================================================
# UMAP feature plots (viridis continuous scale) ----
# ============================================================
paper_grad <- function() scale_color_viridis_c(option = "viridis")
for (set_name in set_names) {
    sc_col <- paste0(set_name, "_sig")
    p <- FeaturePlot(epi, features = sc_col, pt.size = 0.3, order = TRUE) +
        paper_grad() + ggtitle(set_name)
    ggsave(file.path(OUT_DIR, paste0("UMAP_", set_name, ".png")),
        plot = p, width = 10, height = 8, dpi = 200, bg = "white"
    )
}

# ============================================================
# Per-sample split UMAP feature plots (viridis continuous scale) ----
# ============================================================
# split.by = orig.ident -> one panel per sample (CRPC1/2/3). Shared viridis
# limits (global score range) applied to every panel via `&` so samples are
# directly comparable.
n_samples <- length(unique(epi$orig.ident))
for (set_name in set_names) {
    sc_col <- paste0(set_name, "_sig")
    rng <- range(epi@meta.data[[sc_col]], na.rm = TRUE)
    p <- FeaturePlot(epi, features = sc_col, split.by = "orig.ident",
                     pt.size = 0.3, order = TRUE) &
        scale_color_viridis_c(option = "viridis", limits = rng)
    p <- p + plot_annotation(title = set_name)
    ggsave(file.path(OUT_DIR, paste0("UMAP_", set_name, "_splitBySample.png")),
        plot = p, width = 8 * n_samples, height = 8, dpi = 200, bg = "white"
    )
}

# ============================================================
# Save ----
# ============================================================
# Note: saveRDS(epi, OUT_RDS) happens inside the compute branch above.
message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
