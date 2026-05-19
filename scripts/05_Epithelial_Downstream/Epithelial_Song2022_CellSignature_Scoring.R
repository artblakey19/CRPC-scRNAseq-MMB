# Epithelial_Signature_Scoring.R
# Song et al. 2022 (Nat Comm) BE/LE/Hillock/Club signature scoring 재현:
# 저자 GitHub repo (franklinhuanglab/scRNA-seq-Analysis-of-Prostate-Cancer-Samples)
# Seqwell_Combined_E.R 의 두 핵심 패턴을 그대로 따름:
#   - Line 600 패턴: AddModuleScore on dge_E (전체 epithelial), assay=RNA, ctrl=default(100)
#   - Line 603-614: boxplot + stat_summary(mean, shape=95, blue bar)
#   - Line 921:     scale_color_gradientn(c("blue","green","yellow","red"), limits=c(0, max))
#   - Line 927-934: violin + ggpubr::stat_compare_means(method="anova")
# Paper Methods (page 17): "A One-way ANOVA test was then conducted to determine
# if the signature score of each cluster was significantly different from the rest."
#
# 입력:
#   Results/Epithelial/epi_clustered.rds   (Epithelial_Analysis.R 출력)
#   Resources/Song2022_SuppData2.xlsx      (https://www.nature.com/articles/s41467-021-27322-4)
#
# DNPC cohort caveat:
#   - LE score는 AR-loss로 muted; ARPC remnant cluster에서만 강하게 발현 예상.
#   - Hillock signature는 primary PCa에는 없음(Song et al.). DNPC squamous/hillock-like
#     plasticity state 검출 가능성 때문에 계산은 하되 해석 주의.
#   - Club signature가 intermediate / progenitor-like state에 가장 정보량 큼.

library(Seurat)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(ggpubr)   # stat_compare_means — 저자 repo Seqwell_Combined_E.R:929 사용
library(msigdbr)  # HALLMARK_ANDROGEN_RESPONSE fetch (Fig 4h / Fig 5c)

EPI_RDS   <- "Results/Epithelial/epi_clustered.rds"
SIG_XLSX  <- "Resources/Song2022_SuppData2.xlsx"
SIG_SHEET <- "Normal Signature"   # Henry et al. 2018 derived (논문 기본).
                                  # 대안: "PCa signature" — Song et al.의 PCa DEG-derived
                                  # (BE/Club/LE + ERGpos_Tumor/ERGneg_Tumor).
OUT_DIR  <- "Results/Epithelial/Signatures"
OUT_RDS  <- "Results/Epithelial/epi_signature_scored.rds"

stopifnot(
    "epi_clustered.rds not found — run Epithelial_Analysis.R first" = file.exists(EPI_RDS),
    "Song2022_SuppData2.xlsx not found — download Supplementary Data 2 from Nat Comm article" =
        file.exists(SIG_XLSX)
)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

epi <- readRDS(EPI_RDS)

# 저자 repo Seqwell_Combined_E.R:600 — assay = "RNA" 사용.
# Epithelial_Analysis.R에서 RNA layer는 JoinLayers된 raw counts 상태이므로
# NormalizeData 적용 필요.
DefaultAssay(epi) <- "RNA"
epi <- NormalizeData(epi, assay = "RNA", verbose = FALSE)

# ============================================================
# Parse Supplementary Data 2 ----
# ============================================================
# 실제 구조 (확인됨):
#   Sheet "Normal Signature":   BE | Club | Hillock | LE | NE | Endothelial.cells | ...
#   Sheet "PCa signature":      BE | Club | LE | ERGpos_Tumor | ERGneg_Tumor | ...
#   Sheet "Organoid signature": BE | Club | Hillock | hMSC | MKI67+.Cells | Tumor
#   Sheet "PCa Epithelial cell DEG": long-format DEG table (사용하지 않음)
# 각 sheet 내 컬럼 = cell type, 컬럼 안의 값 = gene symbol (variable length, NA padded).

canonical_match <- function(nm) {
    nm_l <- tolower(nm)
    if (grepl("basal|^be$", nm_l)) return("BE")
    if (grepl("luminal|^le$", nm_l)) return("LE")
    if (grepl("hillock", nm_l)) return("Hillock")
    if (grepl("club", nm_l)) return("Club")
    if (grepl("neuroendocrine|^ne$", nm_l)) return("NE")
    NA_character_
}

read_sigs_from_sheet <- function(path, sheet,
                                 wanted = c("BE", "LE", "Hillock", "Club", "NE")) {
    df <- read.xlsx(path, sheet = sheet)
    out <- list()
    for (col_nm in colnames(df)) {
        canon <- canonical_match(col_nm)
        if (is.na(canon) || !canon %in% wanted) next
        genes <- unique(na.omit(as.character(df[[col_nm]])))
        genes <- genes[nzchar(genes)]
        if (length(genes) > 0) out[[canon]] <- genes
    }
    out
}

sigs <- read_sigs_from_sheet(SIG_XLSX, SIG_SHEET)
message("Loaded sheet: ", SIG_SHEET)
message("Signatures parsed: ", paste(names(sigs), collapse = ", "))

sig_order <- intersect(c("BE", "LE", "Hillock", "Club", "NE"), names(sigs))
sigs <- sigs[sig_order]

if (length(sig_order) < 4) {
    warning(
        "Expected BE/LE/Hillock/Club. Got: ",
        paste(sig_order, collapse = ", "),
        ". Check SIG_SHEET configuration."
    )
}

# ============================================================
# Hallmark Androgen Response signature ----
# ============================================================
# Paper 인용 (Methods + Fig 5c bottom, Fig 4h):
#   "we tested for androgen responsiveness among the epithelial cell populations
#    and identified LE cells and tumor cells as the most androgen-responsive as
#    they scored significantly higher than other epithelial cell types in AR
#    signature scores"
#   "Comparison of Hallmark AR pathway signature scores of each epithelial cell type"
#
# 저자 repo는 GenePattern ssGSEA로 외부 실행했지만, per-cell violin/boxplot용
# Fig 4h / Fig 5c bottom은 AddModuleScore + Wilcoxon/ANOVA로 재현 가능.
# (paper Fig 5c bottom: ***P<0.001, Wilcoxon rank-sum test — pairwise Tukey HSD
#  와 동등한 신호 검출.)
#
# DNPC cohort caveat:
#   - AR signaling 전반적으로 muted 예상. 그러나 cluster 7 (ARPC remnant)에서
#     상대적으로 강하게 검출되어야 cluster annotation 검증에 유용.
hallmark_df <- msigdbr(species = "Homo sapiens", collection = "H")
hallmark_ar <- unique(
    hallmark_df$gene_symbol[hallmark_df$gs_name == "HALLMARK_ANDROGEN_RESPONSE"]
)
if (length(hallmark_ar) == 0) {
    warning("HALLMARK_ANDROGEN_RESPONSE not found in msigdbr — skipping AR scoring.")
} else {
    sigs[["AR"]] <- hallmark_ar
    message(sprintf("  AR (HALLMARK_ANDROGEN_RESPONSE): %d genes fetched", length(hallmark_ar)))
}

# Object에 실제로 존재하는 gene만 유지
present <- rownames(epi)
for (s in names(sigs)) {
    n0 <- length(sigs[[s]])
    sigs[[s]] <- intersect(sigs[[s]], present)
    message(sprintf("  %s: %d / %d genes present in object", s, length(sigs[[s]]), n0))
}

# ============================================================
# AddModuleScore ----
# ============================================================
# 저자 repo의 "score all E population" 패턴 (Seqwell_Combined_E.R:600) 정확히 재현:
#   dge <- AddModuleScore(object = dge, features = list(Features),
#                         name = name_feature, assay = "RNA")
# → ctrl 명시 안 함 = Seurat default ctrl=100.
# (저자가 ctrl=5를 쓴 곳은 Fig 3 club integration / Fig 8 organoid subset 분석.
#  Fig 2e/f의 전체-E scoring은 default ctrl 사용.)
score_cols <- character(0)
for (s in names(sigs)) {
    if (length(sigs[[s]]) < 5) {
        message("Skip ", s, " — fewer than 5 valid genes")
        next
    }
    epi <- tryCatch(
        AddModuleScore(epi,
            features = list(sigs[[s]]),
            name = paste0(s, "_Score"),
            assay = "RNA"
        ),
        error = function(e) {
            message("nbin=24 failed for ", s, " — retry with nbin=12: ", e$message)
            AddModuleScore(epi,
                features = list(sigs[[s]]),
                name = paste0(s, "_Score"),
                nbin = 12,
                assay = "RNA"
            )
        }
    )
    # AddModuleScore appends "1" — rename for clarity
    src <- paste0(s, "_Score1")
    dst <- paste0(s, "_Score")
    epi@meta.data[[dst]] <- epi@meta.data[[src]]
    epi@meta.data[[src]] <- NULL
    score_cols <- c(score_cols, dst)
}

# ============================================================
# Boxplot per cluster (저자 repo Seqwell_Combined_E.R:603-614 정확 재현) ----
# ============================================================
# 저자 코드 그대로:
#   V1<-as.vector(FetchData(object = dge,vars = paste0(name_feature,"1")))
#   BOX_df<-NULL
#   BOX_df$id<-dge@active.ident
#   BOX_df$value<-as.vector(V1[,1])
#   BOX_df<-data.frame(BOX_df)
#   ggplot(BOX_df, aes(id,value,fill=id)) + geom_boxplot(notch = F) +
#     stat_summary(fun=mean,geom="point",size=20,colour="blue",shape=95)+
#     ggtitle(name_feature)+NoLegend()
#   ggsave(file=paste0(name_feature,"_PCA_boxplot.eps"),width = 20,height = 20,units = "cm")
make_boxplot <- function(sc) {
    BOX_df <- data.frame(
        id = factor(epi$seurat_clusters),
        value = epi@meta.data[[sc]]
    )
    ggplot(BOX_df, aes(id, value, fill = id)) +
        geom_boxplot(notch = FALSE) +
        stat_summary(
            fun = mean, geom = "point",
            size = 20, colour = "blue", shape = 95
        ) +
        ggtitle(sub("_Score$", "", sc)) +
        NoLegend()
}

for (sc in score_cols) {
    p <- make_boxplot(sc)
    # 저자 spec: width=20cm, height=20cm. ggsave에서 cm 단위 그대로.
    ggsave(file.path(OUT_DIR, paste0(sub("_Score$", "", sc), "_PCA_boxplot.png")),
        plot = p, width = 20, height = 20, units = "cm", dpi = 200, bg = "white"
    )
}

# ============================================================
# Violin per cluster (저자 repo Seqwell_Combined_E.R:927-934 정확 재현) ----
# ============================================================
# Fig 3 club integration 스타일. ggpubr::stat_compare_means(method="anova").
# Paper Fig 2e/f 본문에 "signature score violin plots across all clusters" 명시.
make_violin <- function(sc) {
    BOX_df <- data.frame(
        id = factor(epi$seurat_clusters),
        value = epi@meta.data[[sc]]
    )
    ymax <- max(BOX_df$value, na.rm = TRUE)
    ggplot(BOX_df, aes(id, value, fill = id)) +
        geom_violin() +
        stat_compare_means(method = "anova",
                           label.x = 3, label.y = ymax + 0.05) +
        ggtitle(sub("_Score$", "", sc)) +
        NoLegend() +
        theme(text = element_text(size = 11))
}

for (sc in score_cols) {
    p <- make_violin(sc)
    ggsave(file.path(OUT_DIR, paste0(sub("_Score$", "", sc), "_violin.png")),
        plot = p, width = 20, height = 15, units = "cm", dpi = 200, bg = "white"
    )
}

# ============================================================
# UMAP feature plots (저자 repo line 921 color gradient) ----
# ============================================================
# scale_color_gradientn(colours = c("blue","green","yellow","red"),
#                        limits = c(0, max(V1)))
paper_grad <- function(score_vec) {
    scale_color_gradientn(
        colours = c("blue", "green", "yellow", "red"),
        limits = c(0, max(score_vec, na.rm = TRUE))
    )
}

for (sc in score_cols) {
    scores <- epi@meta.data[[sc]]
    p <- FeaturePlot(epi, features = sc, pt.size = 0.3, order = TRUE) +
        paper_grad(scores) +
        ggtitle(sub("_Score$", "", sc))
    ggsave(file.path(OUT_DIR, paste0("UMAP_", sub("_Score$", "", sc), ".png")),
        plot = p, width = 20, height = 20, units = "cm", dpi = 200, bg = "white"
    )
}

# ============================================================
# One-way ANOVA per signature (paper Methods 문구 정확 재현) ----
# ============================================================
# Paper Methods: "Mean basal, luminal, hillock, and club signature scores were
# calculated for each cluster ... A One-way ANOVA test was then conducted to
# determine if the signature score of each cluster was significantly different
# from the rest."
# → cluster를 factor로 한 one-way ANOVA + Tukey HSD post-hoc.
anova_rows <- list()
tukey_rows <- list()
for (sc in score_cols) {
    df <- data.frame(
        cluster = factor(epi$seurat_clusters),
        score = epi@meta.data[[sc]]
    )
    fit <- aov(score ~ cluster, data = df)
    aov_summary <- summary(fit)[[1]]
    anova_rows[[length(anova_rows) + 1]] <- data.frame(
        signature = sub("_Score$", "", sc),
        df_between = aov_summary["cluster", "Df"],
        df_within = aov_summary["Residuals", "Df"],
        F_stat = aov_summary["cluster", "F value"],
        p_value = aov_summary["cluster", "Pr(>F)"]
    )
    # Tukey HSD — cluster vs cluster pairwise
    tuk <- TukeyHSD(fit)$cluster
    tuk_df <- as.data.frame(tuk)
    tuk_df$signature <- sub("_Score$", "", sc)
    tuk_df$comparison <- rownames(tuk_df)
    rownames(tuk_df) <- NULL
    tukey_rows[[length(tukey_rows) + 1]] <- tuk_df
}
anova_df <- do.call(rbind, anova_rows)
tukey_df <- do.call(rbind, tukey_rows)
tukey_df <- tukey_df[, c("signature", "comparison", "diff", "lwr", "upr", "p adj")]
names(tukey_df)[names(tukey_df) == "p adj"] <- "p_adj"

write.csv(anova_df,
    file.path(OUT_DIR, "OneWay_ANOVA_per_signature.csv"),
    row.names = FALSE
)
write.csv(tukey_df,
    file.path(OUT_DIR, "TukeyHSD_signature_cluster_pairs.csv"),
    row.names = FALSE
)

# ============================================================
# Save ----
# ============================================================
saveRDS(epi, OUT_RDS)
message("Done. Outputs in: ", OUT_DIR)
message("Annotated object: ", OUT_RDS)
