# Epithelial_Song2022_ssGSEA.R
# Song et al. 2022 Nat Comm 재현 — franklinhuanglab repo Seqwell_Combined_E.R 패턴.
# Liu/Tomlins/Wallace/Setlur 4종 paper PCa gene sets에 AddModuleScore 적용.
#
# 자매 스크립트: Epithelial_Bluemn2017_ssGSEA.R (Bluemn 2017 AR-null/FGF-MAPK
# phenotype 재현, GSVA 기반 진짜 ssGSEA). 둘은 출처 논문/목적 다름; 합치지 말 것.
#
# PCa-known gene set들의 epithelial UMAP 상 activity 표시 (참고용).
#
# Note:
#   원래 Song et al. 2022 page 4 two-stage workflow (Stage 1 AddModuleScore +
#   Stage 2 ssGSEA top 1% C2CGP)으로 tumor cluster를 자동 호출했으나, top 1%
#   cutoff와 PCa set 정의가 명확하지 않아 classification은 제거. 단순히 PCa
#   signature가 어느 cluster에서 활성화되는지 시각화만 수행.
#
# DNPC cohort caveat:
#   - AR-driven LIU_PROSTATE_CANCER_UP, SETLUR_TMPRSS2_ERG_FUSION_UP은 muted 예상.
#   - WALLACE_PROSTATE_CANCER_UP / TOMLINS_PROSTATE_CANCER_UP이 cluster 7
#     (ARPC remnant)에서 가장 강하게 잡힐 것.

library(Seurat)
library(dplyr)
library(msigdbr)
library(ggplot2)

IN_RDS  <- "Results/Epithelial/epi_clustered.rds"
OUT_DIR <- "Results/Epithelial/Song2022_PaperPCa"
OUT_RDS <- "Results/Epithelial/epi_song2022_paperpca.rds"

# Paper named PCa gene sets (Seqwell_Combined_E.R line 896-906에서 사용한 4종)
PAPER_PCA_SETS <- c(
    "LIU_PROSTATE_CANCER_UP",
    "TOMLINS_PROSTATE_CANCER_UP",
    "WALLACE_PROSTATE_CANCER_UP",
    "SETLUR_PROSTATE_CANCER_TMPRSS2_ERG_FUSION_UP"
)

stopifnot("Input rds not found" = file.exists(IN_RDS))
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

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

# ============================================================
# UMAP feature plots ----
# ============================================================
paper_grad <- function() scale_color_gradientn(
    colours = c("blue", "green", "yellow", "red")
)
for (set_name in names(pca_sets)) {
    sc_col <- paste0(set_name, "_sig")
    p <- FeaturePlot(epi, features = sc_col, pt.size = 0.3, order = TRUE) +
        paper_grad() + ggtitle(set_name)
    ggsave(file.path(OUT_DIR, paste0("UMAP_", set_name, ".png")),
        plot = p, width = 10, height = 8, dpi = 200, bg = "white"
    )
}

# ============================================================
# Save ----
# ============================================================
saveRDS(epi, OUT_RDS)
message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
