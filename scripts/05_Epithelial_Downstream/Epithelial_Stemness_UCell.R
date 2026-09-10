#!/usr/bin/env Rscript
# =============================================================================
# Epithelial_Stemness_UCell.R
# -----------------------------------------------------------------------------
# 질문: 이 데이터에서 stemness 가 강한 세포는 무엇인가? (암생물학적 stemness)
#
# 접근: 문헌 기반 stemness gene set 패널을 UCell(rank 기반, 셀당 독립·depth 비의존)
#   로 스코어링하고, Seurat AddModuleScore 로 교차검증한다. 여기에 Malta 2018 Cell
#   mRNAsi(one-class LR 가중치 → 셀별 Spearman) 를 별도 축으로 계산한다.
#
# Gene set 패널 (msigdbr, live MSigDB — 프로젝트 컨벤션):
#   [A] ESC/pluripotency-like "oncogenic dedifferentiation" stemness
#       BENPORATH_ES_1, BENPORATH_ES_2        Ben-Porath 2008 Nat Genet
#       WONG_EMBRYONIC_STEM_CELL_CORE         Wong 2008 Cell Stem Cell
#       MUELLER_PLURINET                      Müller 2008 Nature
#       BHATTACHARYA_EMBRYONIC_STEM_CELL      Bhattacharya 2004 Blood
#       MALTA_CURATED_STEMNESS_MARKERS        Malta 2018 Cell (21 curated CSC markers)
#       WP_EMBRYONIC_STEM_CELL_PLURIPOTENCY_PATHWAYS
#       + mRNAsi (Malta 2018; Resources/Malta2018_mRNAsi/SC-pcbc-stemsig.tsv)
#   [B] Somatic / tissue stem-progenitor
#       GOBP_SOMATIC_STEM_CELL_POPULATION_MAINTENANCE
#       LIM_MAMMARY_STEM_CELL_UP              Lim 2009 Nat Med (basal/MaSC; 전립선 basal 유사)
#       RAMALHO_STEMNESS_UP                   Ramalho-Santos 2002 Science ("stemness" 공통 유전자)
#   [C] Prostate CSC 마커 패널 (큐레이션): CD44/PROM1/ITGA6/ALDH1A1·A3/KLF4/SOX2/
#       POU5F1/NANOG/MYC/BMI1/EZH2/LGR5/ABCG2/NES 등 — 전립선 CSC 문헌(Collins 2005,
#       Patrawala 2006, Li 2010, Qin 2012 등) 에서 반복 보고되는 마커
#   [D] 교란 통제: BENPORATH_PROLIFERATION + 기존 S.Score/G2M.Score, nFeature_RNA
#       (ESC-like signature 는 증식과 강하게 공변 → 증식 보정 잔차 점수도 산출)
#
# 해석 caveat (메모리):
#   - benign/malignant 미확정 → stemness 높음 = "암 줄기세포" 단정 금지.
#   - Numbat aneuploid(=malignant 확정) 세포 vs 나머지 비교로 보조 근거만 제시.
#   - 연속형 = viridis, discrete = Okabe-Ito(utils_cb_palette).
#
# 입력: Results/05_Epithelial_Downstream/epi_annotated.rds (상피, annotation)
#       Results/01_Integrated/combined_CRPC.rds              (전체 세포, celltype)
#       Results/08_Numbat/<S>/numbat/clone_post_*.tsv        (aneuploid 라벨, 선택)
# 출력: Results/05_Epithelial_Downstream/Stemness_UCell/
# =============================================================================
suppressMessages({
    library(Seurat); library(UCell); library(msigdbr)
    library(dplyr); library(tidyr); library(ggplot2); library(patchwork)
    library(viridis)
})
source("scripts/00_utils/scRNA_utils.R")

N_CORES <- 30
set.seed(42)

IN_EPI   <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
IN_ALL   <- "Results/01_Integrated/combined_CRPC.rds"
NB_DIR   <- "Results/08_Numbat"
MRNASI_W <- "Resources/Malta2018_mRNAsi/SC-pcbc-stemsig.tsv"
OUT_DIR  <- "Results/05_Epithelial_Downstream/Stemness_UCell"
OUT_EPI  <- file.path(OUT_DIR, "epi_stemness_scores.rds")
OUT_ALL  <- file.path(OUT_DIR, "allcells_stemness_scores.rds")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(file.path(OUT_DIR, "UMAP"), showWarnings = FALSE)
dir.create(file.path(OUT_DIR, "AllCells"), showWarnings = FALSE)

SAMPLES      <- c("CRPC1", "CRPC2", "CRPC3")
LABEL_LEVELS <- c("LE(ARPC)", "Club", "Hillock 1", "Hillock 2",
                  "BE 1", "BE 2", "BE 3", "BE 4", "BE 5", "BE 6", "OE", "Ionocyte")
ANNO_COLS <- setNames(utils_cb_palette(length(LABEL_LEVELS)), LABEL_LEVELS)

# ============================================================
# 1. Gene sets ----
# ============================================================
ESC_SETS <- c("BENPORATH_ES_1", "BENPORATH_ES_2", "WONG_EMBRYONIC_STEM_CELL_CORE",
              "MUELLER_PLURINET", "BHATTACHARYA_EMBRYONIC_STEM_CELL",
              "MALTA_CURATED_STEMNESS_MARKERS",
              "WP_EMBRYONIC_STEM_CELL_PLURIPOTENCY_PATHWAYS")
SOMATIC_SETS <- c("GOBP_SOMATIC_STEM_CELL_POPULATION_MAINTENANCE",
                  "LIM_MAMMARY_STEM_CELL_UP", "RAMALHO_STEMNESS_UP")
CTRL_SETS <- c("BENPORATH_PROLIFERATION")
MSIG_SETS <- c(ESC_SETS, SOMATIC_SETS, CTRL_SETS)

msig <- msigdbr(species = "Homo sapiens")
msig <- msig[msig$gs_name %in% MSIG_SETS, c("gs_name", "gene_symbol")]
gene_sets <- split(msig$gene_symbol, msig$gs_name)
gene_sets <- lapply(gene_sets, unique)
missing <- setdiff(MSIG_SETS, names(gene_sets))
if (length(missing)) stop("Missing in msigdbr: ", paste(missing, collapse = ", "))

# [C] Prostate CSC marker panel (curated)
gene_sets[["PCa_CSC_markers_curated"]] <- c(
    "CD44", "PROM1", "ITGA6", "ITGB1", "ALDH1A1", "ALDH1A3", "ALDH7A1",
    "KLF4", "SOX2", "POU5F1", "NANOG", "MYC", "BMI1", "EZH2", "LGR5",
    "ABCG2", "NES", "TROP2" , "PSCA", "KIT", "LIN28B", "NOTCH1"
)
gene_sets[["PCa_CSC_markers_curated"]] <-
    sub("^TROP2$", "TACSTD2", gene_sets[["PCa_CSC_markers_curated"]])

SET_GROUP <- c(setNames(rep("ESC-like", length(ESC_SETS)), ESC_SETS),
               setNames(rep("Somatic stem", length(SOMATIC_SETS)), SOMATIC_SETS),
               PCa_CSC_markers_curated = "Prostate CSC",
               setNames(rep("Control", length(CTRL_SETS)), CTRL_SETS))
ALL_SETS <- names(SET_GROUP)
gene_sets <- gene_sets[ALL_SETS]

# gene set 정의 저장 (재현성)
gs_df <- bind_rows(lapply(names(gene_sets), function(n)
    data.frame(gene_set = n, group = SET_GROUP[[n]], gene = gene_sets[[n]])))
write.csv(gs_df, file.path(OUT_DIR, "gene_sets_used.csv"), row.names = FALSE)

# mRNAsi weights
stemsig <- read.delim(MRNASI_W, header = FALSE, col.names = c("gene", "w"))
stemsig <- stemsig[!duplicated(stemsig$gene), ]
w_vec   <- setNames(stemsig$w, stemsig$gene)

# ============================================================
# 2. Scoring helpers ----
# ============================================================
short_name <- function(x) {
    x <- sub("^BENPORATH_", "BP_", x)
    x <- sub("^WONG_EMBRYONIC_STEM_CELL_CORE$", "WONG_ESC_CORE", x)
    x <- sub("^BHATTACHARYA_EMBRYONIC_STEM_CELL$", "BHATTACHARYA_ESC", x)
    x <- sub("^MALTA_CURATED_STEMNESS_MARKERS$", "MALTA_CSC_MARKERS", x)
    x <- sub("^WP_EMBRYONIC_STEM_CELL_PLURIPOTENCY_PATHWAYS$", "WP_ESC_PLURIPOTENCY", x)
    x <- sub("^GOBP_SOMATIC_STEM_CELL_POPULATION_MAINTENANCE$", "GOBP_SOMATIC_STEM_MAINT", x)
    x <- sub("^LIM_MAMMARY_STEM_CELL_UP$", "LIM_MAMMARY_STEM_UP", x)
    x <- sub("^PCa_CSC_markers_curated$", "PCa_CSC_curated", x)
    x
}
SHORT <- setNames(short_name(ALL_SETS), ALL_SETS)

# 셀별 mRNAsi: Spearman(weights, log-normalized expr) over shared genes, 청크 처리
compute_mrnasi <- function(obj, w, chunk = 2000L) {
    X <- LayerData(obj, assay = "RNA", layer = "data")
    shared <- intersect(names(w), rownames(X))
    message("  mRNAsi shared genes: ", length(shared), " / ", length(w))
    w  <- w[shared]
    rw <- rank(w)
    n  <- ncol(X); out <- numeric(n); names(out) <- colnames(X)
    idx <- split(seq_len(n), ceiling(seq_len(n) / chunk))
    for (i in idx) {
        m <- as.matrix(X[shared, i, drop = FALSE])
        r <- apply(m, 2, rank)                    # ties = average (Spearman)
        out[i] <- as.numeric(cor(r, rw))
    }
    out
}

score_object <- function(obj, tag) {
    DefaultAssay(obj) <- "RNA"
    if (!"data" %in% Layers(obj[["RNA"]])) obj <- NormalizeData(obj, verbose = FALSE)
    gs_in <- lapply(gene_sets, function(g) intersect(g, rownames(obj)))
    cov <- data.frame(gene_set = names(gene_sets),
                      n_set = lengths(gene_sets), n_detected = lengths(gs_in))
    write.csv(cov, file.path(OUT_DIR, paste0(tag, "_gene_set_coverage.csv")), row.names = FALSE)
    message("[", tag, "] UCell ...")
    obj <- AddModuleScore_UCell(obj, features = gs_in, name = "_UCell",
                                assay = "RNA", slot = "data", ncores = N_CORES)
    message("[", tag, "] AddModuleScore ...")
    for (s in names(gs_in)) {
        obj <- tryCatch(
            AddModuleScore(obj, features = list(gs_in[[s]]), name = "tmpAMS", assay = "RNA"),
            error = function(e) AddModuleScore(obj, features = list(gs_in[[s]]),
                                               name = "tmpAMS", assay = "RNA", nbin = 12))
        obj@meta.data[[paste0(s, "_AMS")]] <- obj$tmpAMS1
        obj$tmpAMS1 <- NULL
    }
    message("[", tag, "] mRNAsi ...")
    mi <- compute_mrnasi(obj, w_vec)
    obj$mRNAsi_raw <- mi
    obj$mRNAsi     <- (mi - min(mi)) / (max(mi) - min(mi))
    obj
}

# ESC-like consensus: 각 ESC set UCell z-score 평균; 증식 보정 잔차
# ESC-like consensus: 각 ESC set UCell z-score 평균.
#   prolifAdj = S/G2M 회귀 잔차, fullAdj = S/G2M + log10(nFeature_RNA) 회귀 잔차
#   (ESC-like 큰 gene set 은 전사체 복잡도(nFeature)와 강하게 공변 → 저복잡도 클러스터
#    (BE 4/BE 6/OE, 메모리: low-count 축) 가 낮게, mRNAsi 는 반대로 높게 잡히는 교란 통제)
add_consensus <- function(md) {
    z <- scale(as.matrix(md[, paste0(ESC_SETS, "_UCell")]))
    md$ESC_consensus <- rowMeans(z)
    md$log_nFeature  <- log10(md$nFeature_RNA)
    md$ESC_consensus_prolifAdj <- resid(lm(ESC_consensus ~ S.Score + G2M.Score, data = md))
    md$ESC_consensus_fullAdj   <- resid(lm(ESC_consensus ~ S.Score + G2M.Score + log_nFeature, data = md))
    md$mRNAsi_nFeatAdj         <- resid(lm(mRNAsi ~ log_nFeature, data = md))
    md$PCa_CSC_fullAdj         <- resid(lm(PCa_CSC_markers_curated_UCell ~ S.Score + G2M.Score + log_nFeature, data = md))
    md
}

# Pseudobulk mRNAsi (Malta 2018 은 bulk 용 설계): annotation × patient 로 count 합산 →
#   log2(CPM+1) → Spearman(weights). 단일세포 sparsity(ties) 문제를 우회한 참조값.
pseudobulk_mrnasi <- function(obj, group_cols, w) {
    cnt <- LayerData(obj, assay = "RNA", layer = "counts")
    grp <- do.call(paste, c(lapply(group_cols, function(g) as.character(obj@meta.data[[g]])), sep = "|"))
    G   <- Matrix::sparse.model.matrix(~ 0 + factor(grp))
    colnames(G) <- levels(factor(grp))
    pb  <- as.matrix(cnt %*% G)
    cpm <- log2(t(t(pb) / colSums(pb)) * 1e6 + 1)
    shared <- intersect(names(w), rownames(cpm))
    data.frame(group = colnames(cpm), n_cells = as.numeric(table(grp)[colnames(cpm)]),
               mRNAsi_pb = as.numeric(cor(cpm[shared, ], w[shared], method = "spearman")))
}

# ============================================================
# 3. Epithelial scoring (load-or-compute) ----
# ============================================================
if (file.exists(OUT_EPI)) {
    message("Loading cached ", OUT_EPI)
    epi_md <- readRDS(OUT_EPI)
    epi <- readRDS(IN_EPI)
    stopifnot(identical(rownames(epi_md), colnames(epi)))
    epi@meta.data <- epi_md
} else {
    epi <- readRDS(IN_EPI)
    epi <- score_object(epi, "epi")
    epi_md <- epi@meta.data

    # Numbat aneuploid 라벨 조인 (있으면)
    epi_md$nb_compartment <- NA_character_; epi_md$nb_pcnv <- NA_real_
    bc <- sub("^P[0-9]+_", "", rownames(epi_md))
    for (S in SAMPLES) {
        fs <- list.files(file.path(NB_DIR, S, "numbat"),
                         pattern = "^clone_post_[0-9]+\\.tsv$", full.names = TRUE)
        if (!length(fs)) next
        it <- as.integer(sub(".*clone_post_([0-9]+)\\.tsv$", "\\1", fs))
        cp <- read.delim(fs[which.max(it)], stringsAsFactors = FALSE)
        i  <- which(epi_md$orig.ident == S)
        m  <- match(bc[i], cp$cell)
        epi_md$nb_compartment[i] <- cp$compartment_opt[m]
        epi_md$nb_pcnv[i]        <- cp$p_cnv[m]
    }
    emb <- Embeddings(epi, "umap")
    epi_md$UMAP_1 <- emb[, 1]; epi_md$UMAP_2 <- emb[, 2]
    epi@meta.data <- epi_md
    saveRDS(epi_md, OUT_EPI)
    write.csv(epi_md[, c("orig.ident", "annotation", grep("_UCell$|_AMS$|mRNAsi|ESC_consensus|nb_", colnames(epi_md), value = TRUE))],
              file.path(OUT_DIR, "epi_stemness_scores_per_cell.csv"))
}
epi_md <- add_consensus(epi_md)
epi_md$annotation <- factor(as.character(epi_md$annotation), levels = LABEL_LEVELS)
epi_md$orig.ident <- factor(as.character(epi_md$orig.ident), levels = SAMPLES)

# pseudobulk mRNAsi (annotation × patient, annotation only)
pb1 <- pseudobulk_mrnasi(epi, c("annotation", "orig.ident"), w_vec) %>%
    separate(group, c("annotation", "orig.ident"), sep = "\\|")
pb2 <- pseudobulk_mrnasi(epi, "annotation", w_vec) %>% rename(annotation = group)
write.csv(pb1, file.path(OUT_DIR, "mRNAsi_pseudobulk_annotation_patient.csv"), row.names = FALSE)
write.csv(pb2, file.path(OUT_DIR, "mRNAsi_pseudobulk_annotation.csv"), row.names = FALSE)
pb1$annotation <- factor(pb1$annotation, levels = LABEL_LEVELS)
p <- ggplot(pb1, aes(annotation, mRNAsi_pb, fill = orig.ident)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_fill_manual(values = setNames(utils_cb_palette(3), SAMPLES), name = NULL) +
    labs(x = NULL, y = "Spearman(weights, log2 CPM)", title = "Pseudobulk mRNAsi (Malta 2018) by annotation × patient") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "Bar_mRNAsi_pseudobulk_annotation_patient.png"), p, width = 11, height = 5, dpi = 200, bg = "white")

UC_COLS   <- paste0(ALL_SETS, "_UCell")
AMS_COLS  <- paste0(ALL_SETS, "_AMS")
EXTRA_COLS <- c("mRNAsi", "mRNAsi_nFeatAdj", "ESC_consensus", "ESC_consensus_prolifAdj",
                "ESC_consensus_fullAdj", "PCa_CSC_fullAdj")
KEY_COLS  <- c(UC_COLS, EXTRA_COLS)
KEY_LAB   <- c(SHORT, mRNAsi = "mRNAsi (Malta2018)",
               mRNAsi_nFeatAdj = "mRNAsi, nFeature-adjusted",
               ESC_consensus = "ESC-like consensus (z-mean)",
               ESC_consensus_prolifAdj = "ESC-like consensus, prolif-adjusted",
               ESC_consensus_fullAdj = "ESC-like consensus, prolif+nFeature-adjusted",
               PCa_CSC_fullAdj = "PCa CSC curated, prolif+nFeature-adjusted")
names(KEY_LAB) <- KEY_COLS

# ============================================================
# 4. Summary tables ----
# ============================================================
long <- epi_md %>%
    mutate(cell = rownames(epi_md)) %>%
    select(cell, orig.ident, annotation, nb_compartment, S.Score, G2M.Score, nFeature_RNA,
           all_of(KEY_COLS)) %>%
    pivot_longer(all_of(KEY_COLS), names_to = "score", values_to = "value") %>%
    mutate(score = factor(score, levels = KEY_COLS))

summ_anno <- long %>% group_by(score, annotation) %>%
    summarise(n = n(), mean = mean(value), median = median(value), sd = sd(value), .groups = "drop") %>%
    group_by(score) %>% mutate(z_of_mean = as.numeric(scale(mean)), rank = rank(-mean)) %>% ungroup()
write.csv(summ_anno, file.path(OUT_DIR, "summary_by_annotation.csv"), row.names = FALSE)

summ_anno_pt <- long %>% group_by(score, orig.ident, annotation) %>%
    summarise(n = n(), mean = mean(value), median = median(value), .groups = "drop")
write.csv(summ_anno_pt, file.path(OUT_DIR, "summary_by_annotation_patient.csv"), row.names = FALSE)

# Kruskal–Wallis (annotation 효과) + epsilon^2
kw <- long %>% group_by(score) %>%
    summarise(kw_p = kruskal.test(value ~ annotation)$p.value,
              eps2 = { h <- kruskal.test(value ~ annotation)$statistic; n <- n(); k <- nlevels(annotation)
                       as.numeric((h - k + 1) / (n - k)) },
              .groups = "drop")
write.csv(kw, file.path(OUT_DIR, "kruskal_by_annotation.csv"), row.names = FALSE)

# Top 5% ESC-consensus cells 구성
top_cut <- quantile(epi_md$ESC_consensus, 0.95)
top_tab <- epi_md %>% mutate(top5 = ESC_consensus >= top_cut) %>%
    group_by(orig.ident, annotation) %>%
    summarise(n = n(), n_top5 = sum(top5), pct_top5 = round(100 * mean(top5), 1), .groups = "drop") %>%
    arrange(desc(pct_top5))
write.csv(top_tab, file.path(OUT_DIR, "top5pct_ESC_consensus_composition.csv"), row.names = FALSE)
top_cut2 <- quantile(epi_md$ESC_consensus_prolifAdj, 0.95)
top_tab2 <- epi_md %>% mutate(top5 = ESC_consensus_prolifAdj >= top_cut2) %>%
    group_by(orig.ident, annotation) %>%
    summarise(n = n(), n_top5 = sum(top5), pct_top5 = round(100 * mean(top5), 1), .groups = "drop") %>%
    arrange(desc(pct_top5))
write.csv(top_tab2, file.path(OUT_DIR, "top5pct_ESC_consensus_prolifAdj_composition.csv"), row.names = FALSE)
top_cut3 <- quantile(epi_md$ESC_consensus_fullAdj, 0.95)
top_tab3 <- epi_md %>% mutate(top5 = ESC_consensus_fullAdj >= top_cut3) %>%
    group_by(orig.ident, annotation) %>%
    summarise(n = n(), n_top5 = sum(top5), pct_top5 = round(100 * mean(top5), 1), .groups = "drop") %>%
    arrange(desc(pct_top5))
write.csv(top_tab3, file.path(OUT_DIR, "top5pct_ESC_consensus_fullAdj_composition.csv"), row.names = FALSE)

# UCell vs AddModuleScore 일치도
ams_cor <- data.frame(gene_set = ALL_SETS,
                      spearman_UCell_vs_AMS = sapply(ALL_SETS, function(s)
                          cor(epi_md[[paste0(s, "_UCell")]], epi_md[[paste0(s, "_AMS")]], method = "spearman")))
write.csv(ams_cor, file.path(OUT_DIR, "UCell_vs_AddModuleScore_concordance.csv"), row.names = FALSE)

# ============================================================
# 5. Figures ----
# ============================================================
theme_set(theme_bw(base_size = 12))

# 5a. Heatmap: annotation × signature (mean UCell z-scored per signature)
hm <- summ_anno %>% mutate(label = KEY_LAB[as.character(score)],
                           label = factor(label, levels = KEY_LAB))
p <- ggplot(hm, aes(annotation, label, fill = z_of_mean)) +
    geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", z_of_mean)), size = 3,
              color = ifelse(abs(hm$z_of_mean) > 1.5, "white", "black")) +
    scale_fill_viridis_c(option = "viridis", name = "z(mean score)\nacross annotations") +
    labs(x = NULL, y = NULL, title = "Stemness signatures across epithelial states",
         subtitle = "Per-signature mean score by annotation, z-scored across annotations (UCell; mRNAsi Spearman)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "Heatmap_meanScore_by_annotation.png"), p, width = 11, height = 7, dpi = 200, bg = "white")

# 5a'. per patient
hm_pt <- summ_anno_pt %>% group_by(score) %>% mutate(z = as.numeric(scale(mean))) %>% ungroup() %>%
    mutate(label = factor(KEY_LAB[as.character(score)], levels = KEY_LAB))
p <- ggplot(hm_pt, aes(annotation, label, fill = z)) + geom_tile(color = "white") +
    facet_wrap(~orig.ident, ncol = 3) +
    scale_fill_viridis_c(option = "viridis", name = "z(mean)") +
    labs(x = NULL, y = NULL, title = "Stemness signatures by annotation × patient") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "Heatmap_meanScore_by_annotation_patient.png"), p, width = 20, height = 7, dpi = 200, bg = "white")

# 5b. Violin by annotation, faceted per signature
p <- ggplot(long, aes(annotation, value, fill = annotation)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", linewidth = 0.2) +
    facet_wrap(~score, scales = "free_y", ncol = 4, labeller = as_labeller(KEY_LAB)) +
    scale_fill_manual(values = ANNO_COLS, guide = "none") +
    labs(x = NULL, y = "score", title = "Stemness scores by epithelial annotation (all patients)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "Violin_by_annotation.png"), p, width = 20, height = 20, dpi = 200, bg = "white")

# 5c. Violin ESC consensus / mRNAsi / prolifAdj by annotation × patient
sub_long <- long %>% filter(score %in% c("ESC_consensus", "ESC_consensus_fullAdj", "mRNAsi", "mRNAsi_nFeatAdj", "PCa_CSC_markers_curated_UCell", "PCa_CSC_fullAdj"))
p <- ggplot(sub_long, aes(annotation, value, fill = annotation)) +
    geom_violin(scale = "width", linewidth = 0.2) +
    geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", linewidth = 0.2) +
    facet_grid(score ~ orig.ident, scales = "free_y", labeller = labeller(score = as_labeller(KEY_LAB))) +
    scale_fill_manual(values = ANNO_COLS, guide = "none") +
    labs(x = NULL, y = "score", title = "Key stemness axes by annotation and patient") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "Violin_keyAxes_by_annotation_patient.png"), p, width = 18, height = 18, dpi = 200, bg = "white")

# 5d. Correlation among signatures + confounders
cor_cols <- c(KEY_COLS, "S.Score", "G2M.Score", "nFeature_RNA", "decontX_contamination")
# within-annotation confound check: annotation 내부에서도 nFeature 와 공변하는지
within <- epi_md %>% group_by(annotation) %>%
    summarise(n = n(),
              rho_ESC_nFeature   = cor(ESC_consensus, nFeature_RNA, method = "spearman"),
              rho_mRNAsi_nFeature = cor(mRNAsi, nFeature_RNA, method = "spearman"),
              rho_ESC_G2M        = cor(ESC_consensus, G2M.Score, method = "spearman"), .groups = "drop")
write.csv(within, file.path(OUT_DIR, "within_annotation_confound_correlations.csv"), row.names = FALSE)
cm <- cor(epi_md[, cor_cols], method = "spearman")
cm_df <- as.data.frame(as.table(cm)); names(cm_df) <- c("a", "b", "rho")
lab_all <- c(KEY_LAB, S.Score = "S.Score", G2M.Score = "G2M.Score",
             nFeature_RNA = "nFeature_RNA", decontX_contamination = "decontX contamination")
cm_df$a <- factor(lab_all[as.character(cm_df$a)], levels = lab_all)
cm_df$b <- factor(lab_all[as.character(cm_df$b)], levels = lab_all)
p <- ggplot(cm_df, aes(a, b, fill = rho)) + geom_tile(color = "white") +
    geom_text(aes(label = sprintf("%.2f", rho)), size = 2.5) +
    scale_fill_viridis_c(option = "viridis", limits = c(-1, 1), name = "Spearman") +
    labs(x = NULL, y = NULL, title = "Signature co-variation and confounders (epithelial cells)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
ggsave(file.path(OUT_DIR, "Correlation_signatures_confounders.png"), p, width = 14, height = 12, dpi = 200, bg = "white")
write.csv(cm, file.path(OUT_DIR, "Correlation_signatures_confounders.csv"))

# 5e. UMAPs (viridis) — all key scores, global + split by patient
for (sc in KEY_COLS) {
    df <- epi_md[order(epi_md[[sc]]), ]
    p <- ggplot(df, aes(UMAP_1, UMAP_2, color = .data[[sc]])) +
        geom_point(size = 0.3) + scale_color_viridis_c(option = "viridis", name = NULL) +
        coord_fixed() + labs(title = KEY_LAB[[sc]]) + theme_classic(base_size = 12)
    ggsave(file.path(OUT_DIR, "UMAP", paste0("UMAP_", sc, ".png")), p, width = 9, height = 8, dpi = 200, bg = "white")
    p <- p + facet_wrap(~orig.ident, ncol = 3) + theme_bw(base_size = 12)
    ggsave(file.path(OUT_DIR, "UMAP", paste0("UMAP_", sc, "_splitBySample.png")), p, width = 22, height = 8, dpi = 200, bg = "white")
}
# annotation reference UMAP
p <- ggplot(epi_md, aes(UMAP_1, UMAP_2, color = annotation)) + geom_point(size = 0.3) +
    scale_color_manual(values = ANNO_COLS) + coord_fixed() + theme_classic(base_size = 12) +
    guides(color = guide_legend(override.aes = list(size = 3))) + labs(title = "Epithelial annotation")
ggsave(file.path(OUT_DIR, "UMAP", "UMAP_annotation_reference.png"), p, width = 9, height = 8, dpi = 200, bg = "white")

# 5f. Numbat aneuploid vs normal (같은 annotation·환자 내 비교)
if (any(!is.na(epi_md$nb_compartment))) {
    nb <- long %>% filter(!is.na(nb_compartment),
                          score %in% c("ESC_consensus", "ESC_consensus_fullAdj", "mRNAsi", "mRNAsi_nFeatAdj",
                                       "PCa_CSC_markers_curated_UCell", "PCa_CSC_fullAdj",
                                       "LIM_MAMMARY_STEM_CELL_UP_UCell"))
    nb_summ <- nb %>% group_by(score, orig.ident, annotation, nb_compartment) %>%
        summarise(n = n(), mean = mean(value), .groups = "drop") %>%
        pivot_wider(names_from = nb_compartment, values_from = c(n, mean)) %>%
        mutate(delta_tumor_minus_normal = mean_tumor - mean_normal)
    nb_test <- nb %>% group_by(score, orig.ident) %>%
        summarise(n_tumor = sum(nb_compartment == "tumor"), n_normal = sum(nb_compartment == "normal"),
                  wilcox_p = if (n_tumor > 5 && n_normal > 5)
                      wilcox.test(value ~ nb_compartment)$p.value else NA_real_,
                  mean_tumor = mean(value[nb_compartment == "tumor"]),
                  mean_normal = mean(value[nb_compartment == "normal"]), .groups = "drop")
    write.csv(nb_summ, file.path(OUT_DIR, "Numbat_tumor_vs_normal_by_annotation.csv"), row.names = FALSE)
    write.csv(nb_test, file.path(OUT_DIR, "Numbat_tumor_vs_normal_wilcox_by_patient.csv"), row.names = FALSE)
    p <- ggplot(nb, aes(annotation, value, fill = nb_compartment)) +
        geom_boxplot(outlier.size = 0.2, linewidth = 0.25, position = position_dodge(width = 0.8)) +
        facet_grid(score ~ orig.ident, scales = "free_y", labeller = labeller(score = as_labeller(KEY_LAB))) +
        scale_fill_manual(values = c(normal = "#0072B2", tumor = "#D55E00"), name = "Numbat") +
        labs(x = NULL, y = "score", title = "Stemness in Numbat aneuploid (tumor) vs diploid (normal) cells") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, "Boxplot_Numbat_tumor_vs_normal.png"), p, width = 18, height = 20, dpi = 200, bg = "white")
}

# 5g. ESC consensus vs proliferation scatter
p <- ggplot(epi_md, aes(S.Score + G2M.Score, ESC_consensus, color = annotation)) +
    geom_point(size = 0.4, alpha = 0.6) + scale_color_manual(values = ANNO_COLS) +
    guides(color = guide_legend(override.aes = list(size = 3))) +
    labs(x = "S.Score + G2M.Score", y = "ESC-like consensus",
         title = sprintf("ESC-like stemness vs proliferation (Spearman rho = %.2f)",
                         cor(epi_md$S.Score + epi_md$G2M.Score, epi_md$ESC_consensus, method = "spearman")))
ggsave(file.path(OUT_DIR, "Scatter_ESCconsensus_vs_proliferation.png"), p, width = 9, height = 7, dpi = 200, bg = "white")

# ============================================================
# 6. All cells (celltype-level overview) ----
# ============================================================
ALL_KEY <- c("BENPORATH_ES_1_UCell", "WONG_EMBRYONIC_STEM_CELL_CORE_UCell",
             "MALTA_CURATED_STEMNESS_MARKERS_UCell", "PCa_CSC_markers_curated_UCell",
             "GOBP_SOMATIC_STEM_CELL_POPULATION_MAINTENANCE_UCell",
             "BENPORATH_PROLIFERATION_UCell", "mRNAsi", "mRNAsi_nFeatAdj",
             "ESC_consensus", "ESC_consensus_fullAdj")
if (file.exists(OUT_ALL)) {
    all_md <- readRDS(OUT_ALL)
} else if (file.exists(IN_ALL)) {
    all_obj <- readRDS(IN_ALL)
    if (!"celltype" %in% colnames(all_obj@meta.data)) stop("combined_CRPC.rds has no celltype")
    all_obj <- score_object(all_obj, "allcells")
    all_md  <- all_obj@meta.data
    if (!all(c("S.Score", "G2M.Score") %in% colnames(all_md))) {
        all_md$S.Score <- 0; all_md$G2M.Score <- 0
    }
    all_md <- add_consensus(all_md)
    emb <- Embeddings(all_obj, "umap"); all_md$UMAP_1 <- emb[, 1]; all_md$UMAP_2 <- emb[, 2]
    saveRDS(all_md, OUT_ALL)
    rm(all_obj); gc()
}
if (exists("all_md")) {
    all_md <- add_consensus(all_md)
    all_md$celltype <- factor(as.character(all_md$celltype))
    # epithelial subtype 라벨 병합
    all_md$group <- as.character(all_md$celltype)
    m <- match(rownames(all_md), rownames(epi_md))
    all_md$group[!is.na(m)] <- paste0("Epi: ", as.character(epi_md$annotation[m[!is.na(m)]]))
    all_md$group <- factor(all_md$group)
    al <- all_md %>% mutate(cell = rownames(all_md)) %>%
        select(cell, orig.ident, celltype, group, all_of(ALL_KEY)) %>%
        pivot_longer(all_of(ALL_KEY), names_to = "score", values_to = "value") %>%
        mutate(score = factor(score, levels = ALL_KEY))
    s_all <- al %>% group_by(score, celltype) %>%
        summarise(n = n(), mean = mean(value), median = median(value), .groups = "drop") %>%
        group_by(score) %>% mutate(rank = rank(-mean)) %>% ungroup()
    write.csv(s_all, file.path(OUT_DIR, "AllCells", "summary_by_celltype.csv"), row.names = FALSE)
    s_grp <- al %>% group_by(score, group) %>%
        summarise(n = n(), mean = mean(value), median = median(value), .groups = "drop") %>%
        group_by(score) %>% mutate(rank = rank(-mean), z = as.numeric(scale(mean))) %>% ungroup()
    write.csv(s_grp, file.path(OUT_DIR, "AllCells", "summary_by_celltype_epiSubtype.csv"), row.names = FALSE)

    ct_cols <- setNames(utils_cb_palette(nlevels(all_md$celltype)), levels(all_md$celltype))
    p <- ggplot(al, aes(celltype, value, fill = celltype)) +
        geom_violin(scale = "width", linewidth = 0.2) +
        geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", linewidth = 0.2) +
        facet_wrap(~score, scales = "free_y", ncol = 4, labeller = as_labeller(KEY_LAB)) +
        scale_fill_manual(values = ct_cols, guide = "none") +
        labs(x = NULL, y = "score", title = "Stemness scores across all cell types") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
    ggsave(file.path(OUT_DIR, "AllCells", "Violin_by_celltype.png"), p, width = 18, height = 12, dpi = 200, bg = "white")

    p <- ggplot(s_grp %>% mutate(label = factor(KEY_LAB[as.character(score)], levels = KEY_LAB)),
                aes(group, label, fill = z)) + geom_tile(color = "white") +
        geom_text(aes(label = sprintf("%.1f", z)), size = 2.8) +
        scale_fill_viridis_c(option = "viridis", name = "z(mean)") +
        labs(x = NULL, y = NULL, title = "Stemness: all cell types + epithelial subtypes (z across groups)") +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
    ggsave(file.path(OUT_DIR, "AllCells", "Heatmap_by_celltype_epiSubtype.png"), p, width = 14, height = 6, dpi = 200, bg = "white")

    for (sc in c("ESC_consensus", "ESC_consensus_fullAdj", "mRNAsi", "mRNAsi_nFeatAdj",
                 "MALTA_CURATED_STEMNESS_MARKERS_UCell", "PCa_CSC_markers_curated_UCell")) {
        df <- all_md[order(all_md[[sc]]), ]
        p <- ggplot(df, aes(UMAP_1, UMAP_2, color = .data[[sc]])) + geom_point(size = 0.3) +
            scale_color_viridis_c(option = "viridis", name = NULL) + coord_fixed() +
            labs(title = paste0(KEY_LAB[[sc]], " — all cells")) + theme_classic(base_size = 12)
        ggsave(file.path(OUT_DIR, "AllCells", paste0("UMAP_", sc, ".png")), p, width = 9, height = 8, dpi = 200, bg = "white")
    }
    p <- ggplot(all_md, aes(UMAP_1, UMAP_2, color = celltype)) + geom_point(size = 0.3) +
        scale_color_manual(values = ct_cols) + coord_fixed() + theme_classic(base_size = 12) +
        guides(color = guide_legend(override.aes = list(size = 3)))
    ggsave(file.path(OUT_DIR, "AllCells", "UMAP_celltype_reference.png"), p, width = 9, height = 8, dpi = 200, bg = "white")
}

# ============================================================
# 7. Console summary ----
# ============================================================
message("\n=== Mean score rank by annotation (1 = highest) ===")
rk <- summ_anno %>% select(score, annotation, rank) %>%
    pivot_wider(names_from = annotation, values_from = rank)
print(as.data.frame(rk), row.names = FALSE)
message("\n=== Kruskal–Wallis (annotation effect) ===")
print(as.data.frame(kw), row.names = FALSE)
message("\n=== UCell vs AddModuleScore Spearman ===")
print(ams_cor, row.names = FALSE)
message("\n=== Top 5% ESC-consensus composition (top 10 rows) ===")
print(head(as.data.frame(top_tab), 10), row.names = FALSE)
message("\n=== Top 5% ESC-consensus (prolif+nFeature-adjusted) composition (top 10 rows) ===")
print(head(as.data.frame(top_tab3), 10), row.names = FALSE)
message("\n=== Pseudobulk mRNAsi by annotation ===")
print(pb2[order(-pb2$mRNAsi_pb), ], row.names = FALSE)
message("\n=== Within-annotation confound correlations ===")
print(as.data.frame(within), row.names = FALSE)
message("\nDone. Outputs: ", OUT_DIR)
