# Epithelial_Kim2024_AR_HGF_Wnt_pseudobulk_GSEA.R
# Kim et al. 2024 Nat Commun (s41467-024-45489-4) Fig 1b 재현 — 단, 원 논문이 SU2C/PCF
# bulk mCRPC에서 "DNPC tumor vs ARPC tumor" preranked GSEA를 그렸다면, 여기서는 본
# cohort의 filtered epithelial cluster 간 contrast로 변환:
#       Club-like      vs ARPC
#       Hillock-like 1 vs ARPC
#       Hillock-like 2 vs ARPC
#       BE 1           vs ARPC
#       BE 2           vs ARPC
#   3종 Kim2024 패널(Androgen response down, HGF up, β-catenin targets up)에 대해
#   pseudobulk(DESeq2) → fgsea preranked → Fig 1b 스타일 running enrichment plot 생성.
#   Hillock-like / BE 모두 sub-cluster를 합치지 않고 각각 contrast.
#
# Cohort caveat (memory):
#   - 본 cohort는 DNPC late-stage이고 ARPC는 transcriptional AR-positive "remnant".
#   - 이 그림은 cross-sectional cluster contrast일 뿐, ARPC→DNPC 방향성 trajectory를
#     주장하지 않는다. ARPC를 reference로 둔 것은 부호 해석 편의일 뿐.
#   - BE 1/2 benign vs malignant 미확정 (project_epi_benign_malignant_unknown).
#
# Design 선택:
#   - pseudobulk: RNA raw counts를 patient × group으로 합산 (3 patient = 3 replicate).
#   - DESeq2 paired design ~ patient + group.
#   - ranking metric: sign(log2FC) × -log10(pvalue)  ← 논문 Methods 표기 그대로
#       "DEGs ... were pre-ranked based on the log2 fold change and adjusted P-value"
#       (Kim 2024 Methods, Human samples and data analysis / RNA-seq and data analysis)
#       padj 대신 raw pvalue 사용: padj 0 → Inf 회피 + 순위 안정. 부호는 log2FC,
#       크기는 -log10(p). 잔여 Inf는 finite max × 1.1로 캡.
#   - gene sets: Epithelial_Kim2024_AR_HGF_Wnt_ssGSEA.R와 동일한 MSigDB v6.0 3종 패널.
#
# 기존 Epithelial_Kim2024_AR_HGF_Wnt_ssGSEA.R(per-cell ssGSEA)는 건드리지 않는 별도 스크립트.

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
OUT_DIR <- "Results/05_Epithelial_Downstream/Kim2024_AR_HGF_Wnt_pseudobulk_GSEA"
GMT_DIR <- "Resources/msigdb_v6.0_GMTs"

# 3종 Kim2024 패널 (Fig 1b) — ssGSEA 스크립트와 동일.
TARGETS <- c(
    "HALLMARK_ANDROGEN_RESPONSE",
    "RUTELLA_RESPONSE_TO_HGF_UP",
    "FEVR_CTNNB1_TARGETS_UP"
)
GMT_FILES <- c(
    "h.all.v6.0.symbols.gmt",
    "c2.cgp.v6.0.symbols.gmt"
)
# 그림용 짧은 라벨 + 표시 순서 (Fig 1b: Androgen response 좌, HGF 중앙, CTNNB1 우).
LABELS <- c(
    "HALLMARK_ANDROGEN_RESPONSE" = "HALLMARK Androgen Response",
    "RUTELLA_RESPONSE_TO_HGF_UP" = "RUTELLA Response to HGF UP",
    "FEVR_CTNNB1_TARGETS_UP"     = "FEVR CTNNB1 Targets UP"
)
# 논문 Fig 1b 부호 기대값 — 색 띠/검토용
EXPECTED_SIGN <- c(
    "HALLMARK_ANDROGEN_RESPONSE" = "down vs ARPC",
    "RUTELLA_RESPONSE_TO_HGF_UP" = "up vs ARPC",
    "FEVR_CTNNB1_TARGETS_UP"     = "up vs ARPC"
)

# Core 5종 contrast (기존 plot 유지용) — 모두 ARPC reference.
# group 값은 공백/구두점 없이 통일 (DESeq2 안전).
CONTRASTS <- list(
    "Club_vs_ARPC"     = "Club",
    "Hillock1_vs_ARPC" = "Hillock1",
    "Hillock2_vs_ARPC" = "Hillock2",
    "BE1_vs_ARPC"      = "BE1",
    "BE2_vs_ARPC"      = "BE2"
)
# Extra 5종 contrast — 전체 epithelial cluster 포함 plot용 (OE 1-4 + Ionocyte-like)
EXTRA_CONTRASTS <- list(
    "OE1_vs_ARPC"          = "OE1",
    "OE2_vs_ARPC"          = "OE2",
    "OE3_vs_ARPC"          = "OE3",
    "OE4_vs_ARPC"          = "OE4",
    "IonocyteLike_vs_ARPC" = "IonocyteLike"
)
ALL_CONTRASTS <- c(CONTRASTS, EXTRA_CONTRASTS)
# 그림 범례용 표시명 (공백 살림, 표시 순서)
DISPLAY <- c(
    "Club_vs_ARPC"         = "Club-like vs ARPC",
    "Hillock1_vs_ARPC"     = "Hillock-like 1 vs ARPC",
    "Hillock2_vs_ARPC"     = "Hillock-like 2 vs ARPC",
    "BE1_vs_ARPC"          = "BE 1 vs ARPC",
    "BE2_vs_ARPC"          = "BE 2 vs ARPC",
    "OE1_vs_ARPC"          = "OE 1 vs ARPC",
    "OE2_vs_ARPC"          = "OE 2 vs ARPC",
    "OE3_vs_ARPC"          = "OE 3 vs ARPC",
    "OE4_vs_ARPC"          = "OE 4 vs ARPC",
    "IonocyteLike_vs_ARPC" = "Ionocyte-like vs ARPC"
)
REF_GROUP <- "ARPC"

dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
stopifnot("Input rds not found" = file.exists(IN_RDS))

# ============================================================
# Gene sets from local v6.0 GMTs ----
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
all_sets <- list()
for (gmt in GMT_FILES) {
    p <- file.path(GMT_DIR, gmt)
    if (!file.exists(p)) stop("Missing GMT: ", p)
    all_sets <- c(all_sets, read_gmt(p))
}
missing <- setdiff(TARGETS, names(all_sets))
if (length(missing) > 0) stop("Sets not found in v6.0 GMTs: ", paste(missing, collapse = ", "))
gene_sets <- all_sets[TARGETS]
for (nm in names(gene_sets)) message(sprintf("  %s: %d genes (v6.0)", nm, length(gene_sets[[nm]])))

# ============================================================
# Pseudobulk: sum RNA counts by patient × group ----
# ============================================================
epi <- readRDS(IN_RDS)
DefaultAssay(epi) <- "RNA"
message("Loaded: ", IN_RDS, " | cells=", ncol(epi), " genes=", nrow(epi))

md <- epi@meta.data
ann <- as.character(md$annotation)
# annotation → CONTRASTS의 공백 없는 group 값으로 매핑.
grp <- dplyr::case_when(
    ann == "ARPC"            ~ "ARPC",
    ann == "Club-like"       ~ "Club",
    ann == "Hillock-like 1"  ~ "Hillock1",
    ann == "Hillock-like 2"  ~ "Hillock2",
    ann == "BE 1"            ~ "BE1",
    ann == "BE 2"            ~ "BE2",
    ann == "OE 1"            ~ "OE1",
    ann == "OE 2"            ~ "OE2",
    ann == "OE 3"            ~ "OE3",
    ann == "OE 4"            ~ "OE4",
    ann == "Ionocyte-like"   ~ "IonocyteLike",
    TRUE                     ~ NA_character_
)
md$group   <- grp
md$patient <- as.character(md$orig.ident)

needed <- c(REF_GROUP, unlist(ALL_CONTRASTS, use.names = FALSE))
keep_cell <- !is.na(md$group) & md$group %in% needed
message("Cells used (10 contrasts + ARPC): ", sum(keep_cell), " / ", nrow(md))

counts <- GetAssayData(epi, assay = "RNA", layer = "counts")[, keep_cell, drop = FALSE]
md_k   <- md[keep_cell, ]
pb_key <- paste(md_k$patient, md_k$group, sep = "__")

# cells × pseudobulk indicator → genes × pseudobulk count matrix
ind <- Matrix::sparse.model.matrix(~ 0 + factor(pb_key))
colnames(ind) <- levels(factor(pb_key))
pb <- as.matrix(counts %*% ind)
storage.mode(pb) <- "integer"

coldata <- data.frame(
    pb_key  = colnames(pb),
    patient = sub("__.*$", "", colnames(pb)),
    group   = sub("^.*__", "", colnames(pb)),
    stringsAsFactors = FALSE
)
ncell <- as.integer(table(pb_key)[coldata$pb_key])
coldata$n_cells <- ncell
message("Pseudobulk profiles (patient × group):")
print(coldata[, c("patient", "group", "n_cells")])

saveRDS(list(pb = pb, coldata = coldata), file.path(OUT_DIR, "pseudobulk_counts.rds"))
write.csv(coldata, file.path(OUT_DIR, "pseudobulk_design.csv"), row.names = FALSE)

# ============================================================
# Per-contrast DESeq2 (paired) + fgsea preranked ----
# ============================================================
run_contrast <- function(test_group) {
    sel <- coldata$group %in% c(REF_GROUP, test_group)
    cd  <- coldata[sel, ]
    cd$group   <- factor(cd$group, levels = c(REF_GROUP, test_group))  # ARPC = reference
    cd$patient <- factor(cd$patient)
    m <- pb[, sel, drop = FALSE]
    m <- m[rowSums(m) >= 10, ]  # 저발현 유전자 제거 → ranking 안정화

    dds <- DESeqDataSetFromMatrix(m, colData = cd, design = ~ patient + group)
    dds <- DESeq(dds, quiet = TRUE)
    res <- results(dds, contrast = c("group", test_group, REF_GROUP))  # +LFC = up vs ARPC
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)

    # ranking metric: sign(log2FC) × -log10(pvalue) — Kim 2024 Methods 표기 그대로.
    # NA pvalue/log2FC 제거 후, Inf(p≈0)는 유한값 max × 1.1로 캡.
    rk <- sign(res_df$log2FoldChange) * (-log10(res_df$pvalue))
    keep <- is.finite(res_df$log2FoldChange) & !is.na(res_df$pvalue) & res_df$pvalue > 0
    rk[!keep] <- NA_real_
    finite_max <- max(abs(rk[is.finite(rk)]), na.rm = TRUE)
    rk[is.infinite(rk) & rk > 0] <-  finite_max * 1.1
    rk[is.infinite(rk) & rk < 0] <- -finite_max * 1.1
    res_df$rank_metric <- rk

    ranks <- rk
    names(ranks) <- res_df$gene
    ranks <- sort(ranks[is.finite(ranks)], decreasing = TRUE)

    set.seed(1)
    fg <- fgsea(pathways = gene_sets, stats = ranks, eps = 0,
                minSize = 5, maxSize = Inf)
    list(res = res_df[order(res_df$rank_metric, decreasing = TRUE), ],
         fgsea = fg, ranks = ranks)
}

results_list <- lapply(ALL_CONTRASTS, function(g) run_contrast(g))
names(results_list) <- names(ALL_CONTRASTS)

# fgsea 결과 취합
fg_all <- bind_rows(lapply(names(results_list), function(cn) {
    fg <- results_list[[cn]]$fgsea
    data.frame(
        contrast = cn,
        pathway  = fg$pathway,
        label    = unname(LABELS[fg$pathway]),
        NES      = fg$NES,
        ES       = fg$ES,
        pval     = fg$pval,
        padj     = fg$padj,
        size     = fg$size,
        leadingEdge = vapply(fg$leadingEdge, paste, character(1), collapse = ";")
    )
}))
write.csv(fg_all, file.path(OUT_DIR, "fgsea_results.csv"), row.names = FALSE)

# DESeq2 결과 저장
for (cn in names(results_list)) {
    write.csv(results_list[[cn]]$res,
              file.path(OUT_DIR, paste0("DESeq2_", cn, ".csv")), row.names = FALSE)
}

# ============================================================
# Figure 1b 스타일: running enrichment plot (3 pathway × contrast 격자) ----
# ============================================================
# 행 = contrast, 열 = pathway (Androgen / HGF / CTNNB1) — Fig 1b 패널 순서.
enr_panels <- list()
for (cn in names(results_list)) {
    ranks <- results_list[[cn]]$ranks
    fg    <- results_list[[cn]]$fgsea
    for (pw in TARGETS) {
        st <- fg[fg$pathway == pw, ]
        ttl <- sprintf("%s\n%s\nNES=%.2f, FDR=%.2g",
                       DISPLAY[[cn]], LABELS[[pw]], st$NES, st$padj)
        enr_panels[[paste(cn, pw, sep = "__")]] <-
            plotEnrichment(gene_sets[[pw]], ranks) +
            labs(title = ttl,
                 x = "Rank in ordered DEG list — sign(log2FC) × -log10(p) (+ = up vs ARPC)",
                 y = "Enrichment score") +
            theme(plot.title = element_text(size = 9))
    }
}

# --- Core 5 contrast 격자 (기존 file 유지: Club/Hillock1/Hillock2/BE1/BE2 vs ARPC) ---
core_keys <- as.vector(t(outer(names(CONTRASTS), TARGETS, paste, sep = "__")))
enr_core  <- enr_panels[core_keys]
p_enr <- wrap_plots(enr_core, ncol = length(TARGETS), byrow = TRUE) +
    plot_annotation(
        title = "Kim 2024 Fig 1b panel — pseudobulk preranked GSEA",
        subtitle = "Rows: Club / Hillock1 / Hillock2 / BE1 / BE2 vs ARPC (DESeq2 paired, ranked by sign(log2FC) × -log10(p)) | Cols: Androgen / HGF / CTNNB1"
    )
ggsave(file.path(OUT_DIR, "RunningEnrichment_Fig1b_style.png"),
       plot = p_enr,
       width  = 10 * length(TARGETS),
       height = 7  * length(CONTRASTS),
       units  = "cm", dpi = 200, bg = "white")

# --- All 10 contrast 격자 (신규: + OE 1-4, Ionocyte-like) ---
all_keys <- as.vector(t(outer(names(ALL_CONTRASTS), TARGETS, paste, sep = "__")))
enr_all  <- enr_panels[all_keys]
p_enr_all <- wrap_plots(enr_all, ncol = length(TARGETS), byrow = TRUE) +
    plot_annotation(
        title = "Kim 2024 Fig 1b panel — pseudobulk preranked GSEA (all epithelial clusters)",
        subtitle = "Rows: Club / Hillock 1/2 / BE 1/2 / OE 1-4 / Ionocyte-like vs ARPC | Cols: Androgen / HGF / CTNNB1"
    )
ggsave(file.path(OUT_DIR, "RunningEnrichment_Fig1b_style_all_clusters.png"),
       plot = p_enr_all,
       width  = 10 * length(TARGETS),
       height = 7  * length(ALL_CONTRASTS),
       units  = "cm", dpi = 200, bg = "white", limitsize = FALSE)

# 단일 contrast별 3패널 (Fig 1b를 그대로 한 줄로) — pathway별 분리 보고용
for (cn in names(results_list)) {
    panels_cn <- enr_panels[grep(paste0("^", cn, "__"), names(enr_panels))]
    pcn <- wrap_plots(panels_cn, nrow = 1) +
        plot_annotation(title = paste0("Kim 2024 Fig 1b — ", DISPLAY[[cn]]))
    ggsave(file.path(OUT_DIR, paste0("RunningEnrichment_Fig1b_", cn, ".png")),
           plot = pcn,
           width = 10 * length(TARGETS), height = 8,
           units = "cm", dpi = 200, bg = "white")
}

# ============================================================
# Supplementary: NES 막대그래프 (3 pathway × 5 contrast 한 장 요약) ----
# ============================================================
stars_of <- function(p) ifelse(is.na(p), "ns",
    ifelse(p < 5e-4, "***", ifelse(p < 5e-3, "**", ifelse(p < 0.05, "*", "ns"))))

# fg_all에 표시용 컬럼 (raw contrast key는 contrast_key로 보존)
fg_all$contrast_key <- fg_all$contrast
fg_all$label        <- factor(fg_all$label, levels = rev(unname(LABELS[TARGETS])))
fg_all$contrast     <- factor(unname(DISPLAY[fg_all$contrast_key]),
                              levels = unname(DISPLAY[names(ALL_CONTRASTS)]))
fg_all$stars        <- stars_of(fg_all$padj)
fg_all$hjust        <- ifelse(fg_all$NES >= 0, -0.2, 1.2)

# --- Core (Hillock/Club/BE) 5개 NES summary (기존 file 유지) ---
fg_core <- fg_all[fg_all$contrast_key %in% names(CONTRASTS), ]
fg_core$contrast <- droplevels(fg_core$contrast)

pal_n <- utils_cb_palette(length(CONTRASTS))  # Okabe-Ito 5색
xpad <- max(abs(fg_core$NES), na.rm = TRUE) * 0.18

p_bar <- ggplot(fg_core, aes(x = NES, y = label, fill = contrast)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_text(aes(label = stars, hjust = hjust),
              position = position_dodge(width = 0.8), size = 3.2) +
    scale_fill_manual(values = pal_n, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    coord_cartesian(xlim = c(min(fg_core$NES) - xpad, max(fg_core$NES) + xpad)) +
    labs(x = "Normalized Enrichment Score (NES)\n(+ = up vs ARPC)", y = NULL,
         title = "Kim 2024 Fig 1b panel — pseudobulk preranked GSEA NES summary",
         subtitle = "Club-like / Hillock-like 1/2 / BE 1/2 vs ARPC (DESeq2 paired, ranked by sign(log2FC) × -log10(p))") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          legend.text = element_text(size = 9),
          panel.grid.major.y = element_blank()) +
    guides(fill = guide_legend(nrow = 2, byrow = TRUE))

ggsave(file.path(OUT_DIR, "NES_barplot_summary.png"),
       plot = p_bar, width = 24, height = 12, units = "cm", dpi = 220, bg = "white")

# --- All clusters (OE 1-4 + Ionocyte-like 추가) NES summary ---
pal_all <- utils_cb_palette(length(ALL_CONTRASTS))  # 10색 (Okabe-Ito interpolated)
xpad_all <- max(abs(fg_all$NES), na.rm = TRUE) * 0.18

p_bar_all <- ggplot(fg_all, aes(x = NES, y = label, fill = contrast)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.8) +
    geom_vline(xintercept = 0, linewidth = 0.4, colour = "grey40") +
    geom_text(aes(label = stars, hjust = hjust),
              position = position_dodge(width = 0.85), size = 2.6) +
    scale_fill_manual(values = pal_all, name = NULL) +
    scale_x_continuous(expand = expansion(mult = c(0.12, 0.12))) +
    coord_cartesian(xlim = c(min(fg_all$NES) - xpad_all, max(fg_all$NES) + xpad_all)) +
    labs(x = "Normalized Enrichment Score (NES)\n(+ = up vs ARPC)", y = NULL,
         title = "Kim 2024 Fig 1b panel — pseudobulk preranked GSEA (all epithelial clusters)",
         subtitle = "Club-like / Hillock-like 1/2 / BE 1/2 / OE 1-4 / Ionocyte-like vs ARPC (DESeq2 paired, ranked by sign(log2FC) × -log10(p))") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          legend.text = element_text(size = 8),
          panel.grid.major.y = element_blank()) +
    guides(fill = guide_legend(nrow = 3, byrow = TRUE))

ggsave(file.path(OUT_DIR, "NES_barplot_summary_all_clusters.png"),
       plot = p_bar_all, width = 30, height = 16, units = "cm", dpi = 220, bg = "white")

# ============================================================
# Compact 콘솔 요약 ----
# ============================================================
message("\n=== fgsea NES (padj) — Kim 2024 Fig 1b panel ===")
summ <- fg_all %>%
    mutate(cell = sprintf("%+.2f (%.1e)%s", NES, padj,
                          ifelse(stars == "ns", "", stars))) %>%
    select(label, contrast, cell) %>%
    tidyr::pivot_wider(names_from = contrast, values_from = cell)
print(as.data.frame(summ), row.names = FALSE)

message("\nExpected sign vs ARPC (논문 Fig 1b):")
for (pw in TARGETS) message(sprintf("  %s: %s", LABELS[[pw]], EXPECTED_SIGN[[pw]]))

message("\nDone. Outputs in: ", OUT_DIR)
