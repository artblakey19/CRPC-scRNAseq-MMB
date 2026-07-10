# Stage 8 — Harmony seed stability of the epithelial reclustering
#
# 목적 (교수 피드백): epithelial reclustering(04번, res 0.4 / 12 cluster / 19411 cell)이
# Harmony 의 RNG seed 선택에 얼마나 의존적인지 검토한다.
#   - 이건 DETERMINISM 검증이 아님. 같은 seed → 같은 결과는 04번의 set.seed(42) 로
#     이미 확보됨(2회 100% 일치).
#   - 여기서 보는 건 ROBUSTNESS: seed 를 바꿔도 클러스터 구조(생물학적 분할)가
#     유지되는가.
#   - 01(integrated) 단계는 lineage 가 명백히 분리돼 검증 불필요; 이 스크립트는
#     오직 04번 epithelial 재클러스터링만 대상으로 한다.
#
# 설계:
#   SCT + PCA 는 결정적이므로 1회만 계산해 고정. Harmony 직전의 set.seed() 만
#   N 개 seed 로 바꿔가며 Harmony → FindNeighbors → FindClusters(res 0.4) 반복.
#   FindClusters 의 random.seed=0 은 고정 → 변동 요인이 오직 Harmony seed 로 격리.
#   각 seed 의 라벨을 canonical(seed 42) 대비 비교:
#     1. ARI / NMI            전역 일치도, seed×seed 히트맵 + reference 대비 bar
#     2. per-cluster Jaccard  Hungarian 매칭, 어떤 클러스터가 견고/취약한가
#     3. per-cell co-assignment frequency  canonical UMAP 오버레이
#     4. cluster 개수 분포 + alluvial(ref → 대표 seed) 흐름
#
# Input:
#   - Results/02_Epithelial_Initial/epi_clustered.rds        (04번과 동일 입력)
#   - Results/03_Epithelial_FilterDecision/filter_decision.csv
#   - Results/04_Epithelial_Filtered/epi_filtered_clustered.rds  (canonical 라벨/UMAP)
# Output:
#   - Results/08_Cluster_Stability/
#       seed_labels.csv                 seed 별 per-cell 라벨 (캐시; 있으면 재계산 skip)
#       ARI_NMI_heatmap.png             seed×seed ARI / NMI
#       ARI_vs_reference_bar.png        각 seed vs canonical(42) ARI / NMI
#       per_cluster_stability.png       클러스터별 평균 Jaccard (생물학 라벨)
#       cell_stability_UMAP.png         per-cell co-assignment frequency (viridis)
#       nclusters_by_seed.png           seed 별 클러스터 개수
#       alluvial_ref_vs_seed.png        canonical → 대표 seed 라벨 흐름
#       summary_ari_nmi.csv             seed 별 ARI/NMI/n_clusters
#       per_cluster_stability.csv       클러스터별 Jaccard 통계

library(Seurat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(glmGamPoi)
library(harmony)
library(aricode)      # ARI, NMI
library(clue)         # solve_LSAP (Hungarian assignment)
library(ggalluvial)
library(pheatmap)
library(viridisLite)

library(future)
plan("multicore", workers = 30)
options(future.globals.maxSize = 128 * 1024^3)

source("scripts/00_utils/scRNA_utils.R")

# ---- Parameters -------------------------------------------------------------
# Seed 목록: 손으로 고른 임의값이 아니라 고정 meta-seed 에서 재현가능하게 무작위
# 추출한다 (cherry-picking 우려 제거 + reproducible). 42 = canonical(04번) 이라
# reference 로 항상 포함. SCTransform/RunPCA 는 각자 내부 seed.use 로 결정적이므로
# 여기서의 set.seed 가 SCT/PCA 재현성을 해치지 않는다(seed42 vs canonical ARI=1 로 검증).
REF_SEED     <- 42
META_SEED    <- 20260708                 # seed 추출 자체의 재현용 (draw 를 고정)
N_RANDOM     <- 20                        # 무작위 seed 개수
set.seed(META_SEED)
RANDOM_SEEDS <- sample.int(1e6, N_RANDOM) # 1..1,000,000 에서 비복원 무작위
SEEDS        <- unique(c(REF_SEED, RANDOM_SEEDS))
RES          <- 0.4
DIMS         <- 1:30
message("Harmony seeds (meta-seed ", META_SEED, "): ",
        paste(SEEDS, collapse = ", "))

IN_RDS       <- "Results/02_Epithelial_Initial/epi_clustered.rds"
DECISION_CSV <- "Results/03_Epithelial_FilterDecision/filter_decision.csv"
CANON_RDS    <- "Results/04_Epithelial_Filtered/epi_filtered_clustered.rds"
OUT_DIR      <- "Results/08_Cluster_Stability"
LABELS_CSV   <- file.path(OUT_DIR, "seed_labels.csv")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# res 0.4 cluster id → 생물학 라벨 (memory: filtered_epi_annotation, 2026-06-30).
# seed 42 재실행은 canonical 과 동일 graph → 동일 번호이므로 이 매핑을 그대로 적용.
# (아래에서 seed 42 vs canonical ARI ≈ 1 로 검증)
CLUSTER_LABELS <- c(
  "0" = "BE1", "7" = "BE2", "10" = "BE3", "1" = "BE4", "2" = "BE5", "8" = "BE6",
  "6" = "OE", "9" = "LE(ARPC)", "3" = "Club", "5" = "Hillock1", "4" = "Hillock2",
  "11" = "Ionocyte"
)
lab_of <- function(cl) ifelse(as.character(cl) %in% names(CLUSTER_LABELS),
                              CLUSTER_LABELS[as.character(cl)],
                              paste0("cl", cl))

# =============================================================================
# 1. Multi-seed clustering (캐시) ----
# =============================================================================
if (!file.exists(LABELS_CSV)) {

    message("== Building fixed SCT + PCA basis (04번과 동일 파이프라인) ==")
    epi <- readRDS(IN_RDS)

    if (!file.exists(DECISION_CSV)) {
        stop("filter_decision.csv 없음: ", DECISION_CSV, call. = FALSE)
    }
    decision <- read.csv(DECISION_CSV, stringsAsFactors = FALSE)
    keep_clusters <- as.character(
        decision$cluster[decision$keep %in% c(TRUE, "TRUE", "true")]
    )
    epi$keep_after_qc <- as.character(epi$seurat_clusters) %in% keep_clusters
    epi <- subset(epi, subset = keep_after_qc)
    message("Epithelial cells after filter: ", ncol(epi))

    DefaultAssay(epi) <- "RNA"
    if ("SCT"     %in% names(epi@assays))     epi[["SCT"]]     <- NULL
    if ("harmony" %in% names(epi@reductions)) epi[["harmony"]] <- NULL
    if ("pca"     %in% names(epi@reductions)) epi[["pca"]]     <- NULL
    if ("umap"    %in% names(epi@reductions)) epi[["umap"]]    <- NULL

    epi[["RNA"]] <- JoinLayers(epi[["RNA"]])
    epi[["RNA"]] <- split(epi[["RNA"]], f = epi$orig.ident)
    epi@meta.data <- dplyr::select(epi@meta.data,
        -matches("_snn_res\\."),
        -any_of(c("seurat_clusters", "keep_after_qc"))
    )

    # 결정적 단계 — seed 무관 (1회만)
    epi <- SCTransform(epi, verbose = FALSE)
    epi <- RunPCA(epi, verbose = FALSE)

    cells  <- colnames(epi)
    labels <- data.frame(cell = cells, stringsAsFactors = FALSE)

    for (s in SEEDS) {
        message("== Harmony seed ", s, " ==")
        if ("harmony" %in% names(epi@reductions)) epi[["harmony"]] <- NULL

        set.seed(s)                     # ← Harmony k-means init 의 유일한 변동 요인
        epi <- IntegrateLayers(
            object = epi, method = HarmonyIntegration,
            orig.reduction = "pca", new.reduction = "harmony",
            normalization.method = "SCT", verbose = FALSE
        )
        epi <- FindNeighbors(epi, reduction = "harmony", dims = DIMS, verbose = FALSE)
        epi <- FindClusters(epi, resolution = RES, verbose = FALSE)   # random.seed=0 고정
        labels[[paste0("seed_", s)]] <- as.character(epi$seurat_clusters)
    }

    write.csv(labels, LABELS_CSV, row.names = FALSE)
    message("Saved per-cell labels → ", LABELS_CSV)
} else {
    message("Loading cached ", LABELS_CSV, " — 지표/그림만 재생성")
    labels <- read.csv(LABELS_CSV, stringsAsFactors = FALSE, check.names = FALSE)
}

seed_cols <- grep("^seed_", names(labels), value = TRUE)
ref_col   <- paste0("seed_", REF_SEED)
L <- lapply(labels[seed_cols], as.character)
names(L) <- seed_cols
ref <- L[[ref_col]]

# ---- Sanity: seed 42 재실행이 canonical(04번 저장본)을 재현하는가 ----
canon_check <- NA_real_
if (file.exists(CANON_RDS)) {
    canon <- readRDS(CANON_RDS)
    canon_lab <- setNames(as.character(canon$SCT_snn_res.0.4), colnames(canon))
    common <- intersect(labels$cell, names(canon_lab))
    if (length(common) > 0) {
        idx <- match(common, labels$cell)
        canon_check <- aricode::ARI(L[[ref_col]][idx], canon_lab[common])
        message(sprintf("Sanity — seed 42 재실행 vs canonical(04번) ARI = %.4f (기대 ≈ 1.0)",
                        canon_check))
    }
} else {
    canon <- NULL
    message("경고: canonical RDS 없음 → UMAP 오버레이/sanity 생략")
}

# =============================================================================
# 2. ARI / NMI matrices ----
# =============================================================================
K <- length(seed_cols)
ari <- matrix(1, K, K, dimnames = list(seed_cols, seed_cols))
nmi <- matrix(1, K, K, dimnames = list(seed_cols, seed_cols))
for (i in seq_len(K - 1)) for (j in (i + 1):K) {
    ari[i, j] <- ari[j, i] <- aricode::ARI(L[[i]], L[[j]])
    nmi[i, j] <- nmi[j, i] <- aricode::NMI(L[[i]], L[[j]])
}
seed_names <- sub("^seed_", "s", seed_cols)
# 표시용 짧은 이름 사본. ari/nmi 자체는 seed_cols 이름을 유지해야
# 뒤쪽의 이름 기반 색인(off_ari[seed_cols], off_ari[alt_cols])이 깨지지 않음.
ari_disp <- ari; nmi_disp <- nmi
dimnames(ari_disp) <- dimnames(nmi_disp) <- list(seed_names, seed_names)

png(file.path(OUT_DIR, "ARI_NMI_heatmap.png"),
    width = 2200, height = 1100, res = 150)
gl <- gridExtra::arrangeGrob(
    pheatmap(ari_disp, display_numbers = TRUE, number_format = "%.2f",
             color = viridis(100), breaks = seq(0, 1, length.out = 101),
             cluster_rows = FALSE, cluster_cols = FALSE,
             main = "Adjusted Rand Index (seed × seed)", silent = TRUE)[[4]],
    pheatmap(nmi_disp, display_numbers = TRUE, number_format = "%.2f",
             color = viridis(100), breaks = seq(0, 1, length.out = 101),
             cluster_rows = FALSE, cluster_cols = FALSE,
             main = "Normalized Mutual Information (seed × seed)", silent = TRUE)[[4]],
    ncol = 2)
grid::grid.draw(gl)
dev.off()

# 각 seed vs canonical(42) — ari/nmi 는 seed_cols 이름 유지 상태
ref_col_i <- which(seed_cols == ref_col)
off_ari <- ari[ref_col, ]
off_nmi <- nmi[ref_col, ]
vs_ref <- data.frame(
    seed = factor(seed_names, levels = seed_names),
    ARI  = as.numeric(off_ari),
    NMI  = as.numeric(off_nmi)
) %>% filter(seed != seed_names[ref_col_i])

p_bar <- vs_ref %>%
    tidyr::pivot_longer(c(ARI, NMI), names_to = "metric", values_to = "value") %>%
    ggplot(aes(reorder(seed, value), value, fill = metric)) +
    geom_col(position = "dodge") +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
    coord_flip(ylim = c(0, 1)) +
    scale_fill_manual(values = utils_cb_palette(2)) +
    labs(x = "Harmony seed", y = "agreement vs canonical (seed 42)",
         title = "Cluster agreement across Harmony seeds",
         subtitle = "dashed = 0.8 (관례적 'stable' 기준)") +
    theme_bw()
ggsave(file.path(OUT_DIR, "ARI_vs_reference_bar.png"),
       p_bar, width = 8, height = 7, bg = "white")

# =============================================================================
# 3. Hungarian 매칭 → per-cell relabel & per-cluster Jaccard ----
# =============================================================================
# alt 라벨을 ref 라벨로 매핑 (총 overlap 최대화)
hungarian_map <- function(alt, ref) {
    ct <- table(alt, ref)
    nr <- nrow(ct); nc <- ncol(ct); n <- max(nr, nc)
    M <- matrix(0, n, n); M[seq_len(nr), seq_len(nc)] <- ct
    asg <- clue::solve_LSAP(max(M) - M)          # 비용 최소화 = overlap 최대화
    map <- setNames(rep(NA_character_, nr), rownames(ct))
    for (i in seq_len(nr)) if (asg[i] <= nc) map[i] <- colnames(ct)[asg[i]]
    map
}

# ref 클러스터 하나가 alt 클러스터와 갖는 best Jaccard
per_cluster_jac <- function(alt, ref) {
    ref_lvls <- sort(unique(ref)); alt_lvls <- sort(unique(alt))
    sapply(ref_lvls, function(r) {
        inR <- ref == r
        max(vapply(alt_lvls, function(c) {
            inC <- alt == c
            sum(inR & inC) / sum(inR | inC)
        }, numeric(1)))
    })
}

alt_cols <- setdiff(seed_cols, ref_col)

# per-cell co-assignment frequency: 각 alt seed 에서 ref 라벨과 동일하게 남는 비율
agree <- vapply(alt_cols, function(sc) {
    m <- hungarian_map(L[[sc]], ref)
    mapped <- m[L[[sc]]]
    as.integer(!is.na(mapped) & mapped == ref)
}, integer(length(ref)))
stability_cell <- rowMeans(agree)

# per-cluster Jaccard 행렬 (행=ref cluster, 열=alt seed)
jac_mat <- vapply(alt_cols, function(sc) per_cluster_jac(L[[sc]], ref),
                  numeric(length(unique(ref))))
rownames(jac_mat) <- sort(unique(ref))

cluster_stab <- data.frame(
    cluster    = rownames(jac_mat),
    label      = lab_of(rownames(jac_mat)),
    mean_jac   = rowMeans(jac_mat),
    min_jac    = apply(jac_mat, 1, min),
    sd_jac     = apply(jac_mat, 1, sd),
    n_cells    = as.integer(table(ref)[rownames(jac_mat)]),
    stringsAsFactors = FALSE
) %>% arrange(desc(mean_jac))

write.csv(cluster_stab, file.path(OUT_DIR, "per_cluster_stability.csv"),
          row.names = FALSE)

p_clus <- cluster_stab %>%
    mutate(label = factor(label, levels = rev(label))) %>%
    ggplot(aes(label, mean_jac, fill = mean_jac)) +
    geom_col() +
    geom_errorbar(aes(ymin = pmax(0, mean_jac - sd_jac),
                      ymax = pmin(1, mean_jac + sd_jac)), width = 0.3) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "grey40") +
    coord_flip(ylim = c(0, 1)) +
    scale_fill_viridis_c(limits = c(0, 1)) +
    labs(x = NULL, y = "mean best-match Jaccard vs canonical (across seeds)",
         title = "어느 클러스터가 seed 에 견고/취약한가",
         subtitle = "높을수록 seed 를 바꿔도 그대로 재현. 오차막대=±SD") +
    theme_bw()
ggsave(file.path(OUT_DIR, "per_cluster_stability.png"),
       p_clus, width = 8, height = 6, bg = "white")

# =============================================================================
# 4. per-cell stability on canonical UMAP ----
# =============================================================================
if (!is.null(canon) && "umap" %in% names(canon@reductions)) {
    emb <- canon@reductions$umap@cell.embeddings
    idx <- match(rownames(emb), labels$cell)
    ok  <- !is.na(idx)
    df <- data.frame(
        UMAP_1 = emb[ok, 1], UMAP_2 = emb[ok, 2],
        stability = stability_cell[idx[ok]],
        cluster   = lab_of(ref[idx[ok]])
    )
    p_um <- ggplot(df, aes(UMAP_1, UMAP_2, color = stability)) +
        geom_point(size = 0.3) +
        scale_color_viridis_c(limits = c(0, 1),
                              name = "co-assignment\nfrequency") +
        labs(title = "Per-cell cluster stability across Harmony seeds",
             subtitle = "각 셀이 canonical 라벨을 유지한 seed 비율 (1 = 항상 동일)") +
        theme_bw() + coord_equal()
    ggsave(file.path(OUT_DIR, "cell_stability_UMAP.png"),
           p_um, width = 9, height = 8, bg = "white")
}

# =============================================================================
# 5. cluster 개수 분포 + alluvial ----
# =============================================================================
nclus <- sapply(seed_cols, function(sc) length(unique(L[[sc]])))
ndf <- data.frame(seed = factor(seed_names, levels = seed_names),
                  n = as.integer(nclus))
p_n <- ggplot(ndf, aes(seed, n, fill = seed == seed_names[ref_col_i])) +
    geom_col() +
    geom_text(aes(label = n), vjust = -0.3, size = 3) +
    scale_fill_manual(values = c("grey60", "#D55E00"), guide = "none") +
    labs(x = "Harmony seed", y = "# clusters (res 0.4)",
         title = "Cluster 개수의 seed 의존성",
         subtitle = "주황 = canonical(seed 42)") +
    theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(OUT_DIR, "nclusters_by_seed.png"),
       p_n, width = 9, height = 5, bg = "white")

# 대표 seed = ref 대비 ARI 가 중앙값인 seed (전형적 케이스)
med_seed <- alt_cols[which.min(abs(as.numeric(off_ari[alt_cols]) -
                                   median(as.numeric(off_ari[alt_cols]))))]
tryCatch({
    adf <- data.frame(
        ref = lab_of(ref),
        alt = paste0("c", L[[med_seed]])
    ) %>% count(ref, alt)
    fill_n <- length(unique(adf$ref))
    p_al <- ggplot(adf, aes(axis1 = ref, axis2 = alt, y = n)) +
        geom_alluvium(aes(fill = ref), alpha = 0.7) +
        geom_stratum(width = 0.35) +
        geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 2.6) +
        scale_x_discrete(limits = c("canonical (seed 42)",
                                    sub("^seed_", "seed ", med_seed)),
                         expand = c(0.1, 0.1)) +
        scale_fill_manual(values = utils_cb_palette(fill_n), guide = "none") +
        labs(y = "cells", x = NULL,
             title = "라벨 흐름: canonical → 대표 seed",
             subtitle = paste0("대표 seed = ", sub("^seed_", "", med_seed),
                               " (ARI 중앙값). split/merge 위치가 보임")) +
        theme_minimal()
    ggsave(file.path(OUT_DIR, "alluvial_ref_vs_seed.png"),
           p_al, width = 8, height = 9, bg = "white")
}, error = function(e) message("alluvial skip: ", conditionMessage(e)))

# =============================================================================
# 6. Summary ----
# =============================================================================
summary_df <- data.frame(
    seed        = sub("^seed_", "", seed_cols),
    n_clusters  = as.integer(nclus),
    ARI_vs_ref  = round(as.numeric(off_ari[seed_cols]), 4),
    NMI_vs_ref  = round(as.numeric(off_nmi[seed_cols]), 4),
    is_canonical = seed_cols == ref_col
)
write.csv(summary_df, file.path(OUT_DIR, "summary_ari_nmi.csv"), row.names = FALSE)

mean_off_ari <- mean(ari[upper.tri(ari)])
mean_off_nmi <- mean(nmi[upper.tri(nmi)])
message("\n================ SUMMARY ================")
message(sprintf("Seeds tested            : %d", K))
message(sprintf("seed42 vs canonical ARI : %.4f", canon_check))
message(sprintf("Mean pairwise ARI       : %.3f", mean_off_ari))
message(sprintf("Mean pairwise NMI       : %.3f", mean_off_nmi))
message(sprintf("Cluster count range     : %d–%d (canonical %d)",
                min(nclus), max(nclus), nclus[ref_col]))
message("Most stable cluster     : ",
        cluster_stab$label[1], sprintf(" (Jaccard %.2f)", cluster_stab$mean_jac[1]))
message("Least stable cluster    : ",
        cluster_stab$label[nrow(cluster_stab)],
        sprintf(" (Jaccard %.2f)", cluster_stab$mean_jac[nrow(cluster_stab)]))
message("Outputs → ", OUT_DIR)
