#!/usr/bin/env Rscript
# =============================================================================
# 13 Henry-ref integration — step 02 : 우리 상피 × Henry 2018 정상 상피 통합·유사도
# -----------------------------------------------------------------------------
# 질문: 우리 CRPC 상피 클러스터가 Henry 정상 BE/Club/Hillock/LE 와 얼마나 비슷한가?
# 세 가지 독립 축으로 본다(서로 가정이 달라 교차검증용).
#   A. 공동 embedding : RNA counts 병합 → batch(CRPC1-3, D17/D27/D35)별 SCT → Harmony.
#        - 우리 세포별 최근접 Henry 세포 k=30 의 Population 구성 (어느 정상 타입 옆에 앉는가)
#        - 전체 kNN 중 Henry 비율 / 기대치 (정상 세포와 실제로 섞이는가; 1=완전 혼합, 0=분리)
#   B. Reference mapping : Henry 를 reference 로 환자별 FindTransferAnchors → MapQuery.
#        - predicted Population + prediction score + MappingScore (reference 가 query 를
#          얼마나 잘 표현하는가; 낮으면 정상 atlas 에 없는 상태)
#   C. Pseudobulk centroid 상관 : dataset 내부 평균으로 centering(10x v2 vs 우리 chemistry
#        offset 제거) 후 HVG Pearson.
# ⚠ 해석: label transfer 는 항상 가장 가까운 정상 라벨을 강제로 붙인다 → 라벨 자체보다
#   score/mixing 을 볼 것. 그리고 "정상 타입과 유사" ≠ benign (memory: CNV-quiet malignant,
#   DNPC plasticity). chemistry batch 는 09/henryref inferCNV 를 망쳤던 그 축이다.
# =============================================================================
suppressPackageStartupMessages({
    library(Seurat); library(Matrix); library(dplyr); library(tidyr)
    library(ggplot2); library(patchwork); library(glmGamPoi); library(harmony)
    library(pheatmap); library(viridisLite); library(RANN)
})
library(future)
N_CORES <- 30
plan("multicore", workers = N_CORES)
options(future.globals.maxSize = 128 * 1024^3)
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS   <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
HENRY_DIR <- "Raw_data/Henry_GSE120716/epi_matrix_all"
# Henry 라벨 모드: "cellxgene" = h5ad Population 원본 그대로 / "corrected" = Club<->Hillock 교정.
#   근거(2026-09-21): 저자 GitLab(StrandLab/sc-TissueMapper_Pr tag 2.0.0) genes.deg.OE1.csv 상위 =
#   SCGB3A1/LCN2/PIGR/PSCA, genes.deg.OE2.csv 1위 = KRT13; 논문 본문 "OE1 SCGB1A1+"(club),
#   "OE2 KRT13+"(hillock). CELLxGENE h5ad 는 'Hillock'(n=1312)=SCGB3A1/PIGR high, 'Club'(n=2530)=KRT13
#   high 로 이름만 뒤바뀌어 있음(count 는 GEO GSE117403 원본과 동일). 두 모드는 출력 폴더가 분리된다.
HENRY_LABELS <- "corrected"
OUT       <- file.path("Results/13_HenryRef_Integration",
                       if (HENRY_LABELS == "corrected") "label_corrected" else "")
VIS       <- file.path(OUT, "Visualization")
dir.create(VIS, showWarnings = FALSE, recursive = TRUE)
K_NN <- 30; DIMS <- 1:30
HENRY_LV <- c("BE", "Club", "Hillock", "LE")

# 1. Load ---------------------------------------------------------------------
epi <- readRDS(EPI_RDS)
ANN_LV <- levels(droplevels(factor(epi$annotation)))
h_counts <- ReadMtx(file.path(HENRY_DIR, "matrix.mtx.gz"),
                    cells = file.path(HENRY_DIR, "barcodes.tsv.gz"),
                    features = file.path(HENRY_DIR, "features.tsv.gz"), feature.column = 2)
h_meta <- read.csv(file.path(HENRY_DIR, "cell_meta.csv"), row.names = 1)
stopifnot(identical(rownames(h_meta), colnames(h_counts)))
h_meta$Population_cellxgene <- h_meta$Population
if (HENRY_LABELS == "corrected") {
    h_meta$Population <- dplyr::recode(h_meta$Population_cellxgene, Club = "Hillock", Hillock = "Club")
}
cat("Henry labels:", HENRY_LABELS, "\n"); print(table(cellxgene = h_meta$Population_cellxgene, used = h_meta$Population))

e_counts <- LayerData(epi, assay = "RNA", layer = "counts")
genes <- intersect(rownames(e_counts), rownames(h_counts))
cat(sprintf("shared genes: %d (ours %d, Henry %d)\n", length(genes), nrow(e_counts), nrow(h_counts)))

meta <- rbind(
    data.frame(row.names = colnames(epi), dataset = "CRPC", batch = as.character(epi$orig.ident),
               label = as.character(epi$annotation)),
    data.frame(row.names = colnames(h_counts), dataset = "Henry", batch = h_meta$donor_id,
               label = paste("Henry", h_meta$Population))
)
LABEL_LV <- c(ANN_LV, paste("Henry", HENRY_LV))
meta$label <- factor(meta$label, levels = LABEL_LV)

# 2. A: joint SCT + Harmony -----------------------------------------------------
JOINT_RDS <- file.path(OUT, "joint_crpc_henry.rds")
if (file.exists(JOINT_RDS)) {
    joint <- readRDS(JOINT_RDS)
} else {
    joint <- CreateSeuratObject(cbind(e_counts[genes, ], h_counts[genes, ]), meta.data = meta)
    joint[["RNA"]] <- split(joint[["RNA"]], f = joint$batch)
    joint <- SCTransform(joint, verbose = FALSE)
    joint <- RunPCA(joint, verbose = FALSE)
    set.seed(42)
    joint <- IntegrateLayers(joint, method = HarmonyIntegration, orig.reduction = "pca",
                             new.reduction = "harmony", normalization.method = "SCT", verbose = FALSE)
    joint[["RNA"]] <- JoinLayers(joint[["RNA"]])
    joint <- RunUMAP(joint, reduction = "harmony", dims = DIMS,
                     n.neighbors = 40, min.dist = 0.3, spread = 1.5, verbose = FALSE)
    saveRDS(joint, JOINT_RDS)
}

emb  <- Embeddings(joint, "harmony")[, DIMS]
is_h <- joint$dataset == "Henry"
# (a) 최근접 Henry 세포 k개의 Population 구성
nn_h  <- nn2(emb[is_h, ], emb[!is_h, ], k = K_NN)$nn.idx
h_pop <- sub("^Henry ", "", as.character(joint$label[is_h]))
nn_frac <- sapply(HENRY_LV, function(p) rowMeans(matrix(h_pop[nn_h] == p, nrow = nrow(nn_h))))
colnames(nn_frac) <- paste0("nnHenry_", HENRY_LV)
# (b) 전체 kNN 중 Henry 비율 / 기대치(전체 Henry 비율); self 제외 위해 k+1
nn_all <- nn2(emb, emb[!is_h, ], k = K_NN + 1)$nn.idx[, -1]
mix <- rowMeans(matrix(is_h[nn_all], nrow = nrow(nn_all))) / mean(is_h)

cell_df <- data.frame(cell = colnames(joint)[!is_h], patient = joint$batch[!is_h],
                      annotation = factor(as.character(joint$label[!is_h]), levels = ANN_LV),
                      nn_frac, henry_mixing_ratio = mix)

# 3. B: reference mapping (환자별) -----------------------------------------------
MAP_CSV <- file.path(OUT, "reference_mapping_per_cell.csv")
if (file.exists(MAP_CSV)) {
    map_df <- read.csv(MAP_CSV)
} else {
    ref <- CreateSeuratObject(h_counts[genes, ], meta.data = h_meta)
    ref <- SCTransform(ref, verbose = FALSE)
    ref <- RunPCA(ref, verbose = FALSE)
    ref <- RunUMAP(ref, dims = DIMS, return.model = TRUE, verbose = FALSE)
    map_df <- bind_rows(lapply(sort(unique(as.character(epi$orig.ident))), function(S) {
        cells <- colnames(epi)[epi$orig.ident == S]
        q <- CreateSeuratObject(e_counts[genes, cells])
        q <- SCTransform(q, verbose = FALSE)
        anc <- FindTransferAnchors(reference = ref, query = q, normalization.method = "SCT",
                                   reference.reduction = "pca", dims = DIMS, verbose = FALSE)
        q <- MapQuery(anchorset = anc, reference = ref, query = q, refdata = list(pop = "Population"),
                      reference.reduction = "pca", reduction.model = "umap", verbose = FALSE)
        ms <- MappingScore(anc, ndim = max(DIMS), verbose = FALSE)
        u <- Embeddings(q, "ref.umap")
        data.frame(cell = colnames(q), predicted_pop = q$predicted.pop,
                   prediction_score = q$predicted.pop.score, mapping_score = ms[colnames(q)],
                   refUMAP_1 = u[, 1], refUMAP_2 = u[, 2])
    }))
    ref_umap <- data.frame(Embeddings(ref, "umap"), Population = ref$Population)
    write.csv(ref_umap, file.path(OUT, "henry_reference_umap.csv"))
    write.csv(map_df, MAP_CSV, row.names = FALSE)
}
ref_umap <- read.csv(file.path(OUT, "henry_reference_umap.csv"), row.names = 1)
colnames(ref_umap)[1:2] <- c("refUMAP_1", "refUMAP_2")
cell_df <- left_join(cell_df, map_df, by = "cell")
write.csv(cell_df, file.path(OUT, "similarity_per_cell.csv"), row.names = FALSE)

# 4. 요약 테이블 ------------------------------------------------------------------
summ <- function(d, ...) d %>% group_by(...) %>% summarise(
    n = n(), across(starts_with("nnHenry_"), mean), henry_mixing_ratio = mean(henry_mixing_ratio),
    prediction_score = median(prediction_score), mapping_score = median(mapping_score), .groups = "drop")
add_pred <- function(s, d, keys) {
    pf <- d %>% count(across(all_of(keys)), predicted_pop) %>% group_by(across(all_of(keys))) %>%
        mutate(f = n / sum(n)) %>% select(-n) %>%
        pivot_wider(names_from = predicted_pop, values_from = f, values_fill = 0, names_prefix = "pred_")
    left_join(s, pf, by = keys)
}
by_ann <- add_pred(summ(cell_df, annotation), cell_df, "annotation")
by_ann_pt <- add_pred(summ(cell_df, annotation, patient), cell_df, c("annotation", "patient"))
write.csv(by_ann, file.path(OUT, "similarity_by_annotation.csv"), row.names = FALSE)
write.csv(by_ann_pt, file.path(OUT, "similarity_by_annotation_patient.csv"), row.names = FALSE)
print(as.data.frame(by_ann), digits = 2)

# 5. C: centered pseudobulk 상관 --------------------------------------------------
pb <- function(counts, grp) {
    m <- t(rowsum(t(as.matrix(counts)), grp))          # genes x groups
    log1p(t(t(m) / colSums(m)) * 1e6)
}
hvg <- intersect(VariableFeatures(joint), genes)
pb_e <- pb(e_counts[hvg, ], as.character(epi$annotation))[, ANN_LV]
pb_h <- pb(h_counts[hvg, ], h_meta$Population)[, HENRY_LV]
cor_mat <- cor(pb_e - rowMeans(pb_e), pb_h - rowMeans(pb_h), method = "pearson")
colnames(cor_mat) <- paste("Henry", HENRY_LV)
write.csv(cor_mat, file.path(OUT, "centroid_correlation_centered.csv"))

# 6. Plots -----------------------------------------------------------------------
joint$label <- factor(joint$label, levels = LABEL_LV)
h_cols <- setNames(c("#0072B2", "#E69F00", "#009E73", "#CC79A7"), paste("Henry", HENRY_LV))
a_cols <- setNames(utils_cb_palette(length(ANN_LV)), ANN_LV)
ud <- data.frame(Embeddings(joint, "umap"), dataset = joint$dataset, label = joint$label, batch = joint$batch)
colnames(ud)[1:2] <- c("U1", "U2")
base <- function(d, col, cols, ttl) ggplot(d, aes(U1, U2, color = .data[[col]])) +
    geom_point(data = ud, color = "grey90", size = 0.1) + geom_point(size = 0.15) +
    scale_color_manual(values = cols, name = NULL) + ggtitle(ttl) + theme_classic() +
    guides(color = guide_legend(override.aes = list(size = 3)))
p1 <- ggplot(ud[sample(nrow(ud)), ], aes(U1, U2, color = dataset)) + geom_point(size = 0.1) +
    scale_color_manual(values = c(CRPC = "#D55E00", Henry = "#0072B2")) + theme_classic() +
    ggtitle("Joint Harmony UMAP — dataset") + guides(color = guide_legend(override.aes = list(size = 3)))
p2 <- base(ud[ud$dataset == "Henry", ], "label", h_cols, "Henry normal (Population)")
p3 <- base(ud[ud$dataset == "CRPC", ], "label", a_cols, "CRPC epithelium (annotation)")
ggsave(file.path(VIS, "joint_umap_overview.png"), p1 + p2 + p3, width = 21, height = 6, dpi = 150, bg = "white")
p4 <- base(ud[ud$dataset == "CRPC", ], "label", a_cols, "CRPC by patient") + facet_wrap(~batch)
ggsave(file.path(VIS, "joint_umap_crpc_by_patient.png"), p4, width = 18, height = 6, dpi = 150, bg = "white")

um <- cbind(ud[!is_h, ], henry_mixing_ratio = pmin(mix, 2))
p5 <- ggplot(um[order(um$henry_mixing_ratio), ], aes(U1, U2, color = henry_mixing_ratio)) +
    geom_point(size = 0.15) + scale_color_viridis_c(option = "viridis", name = "obs/exp") +
    theme_classic() + ggtitle("Henry fraction in kNN (obs/expected, capped 2)")
ggsave(file.path(VIS, "joint_umap_henry_mixing.png"), p5, width = 7.5, height = 6, dpi = 150, bg = "white")

hm <- function(mat, file, main, breaks = seq(0, 1, length.out = 101), w = 5.5) pheatmap(
    mat, cluster_rows = FALSE, cluster_cols = FALSE, display_numbers = TRUE, number_format = "%.2f",
    color = viridis(100), breaks = breaks, main = main, filename = file, width = w, height = 6, number_color = "white")
m_nn <- as.matrix(by_ann[, paste0("nnHenry_", HENRY_LV)]); dimnames(m_nn) <- list(by_ann$annotation, paste("Henry", HENRY_LV))
hm(m_nn, file.path(VIS, "heatmap_nearest_henry_composition.png"), "Nearest-30 Henry cells: Population fraction")
for (cn in paste0("pred_", HENRY_LV)) if (!cn %in% colnames(by_ann)) by_ann[[cn]] <- 0
m_pr <- as.matrix(by_ann[, paste0("pred_", HENRY_LV)]); m_pr[is.na(m_pr)] <- 0
dimnames(m_pr) <- list(by_ann$annotation, paste("Henry", HENRY_LV))
hm(m_pr, file.path(VIS, "heatmap_label_transfer_fraction.png"), "Label transfer: predicted Population fraction")
lim <- max(abs(cor_mat))
hm(cor_mat, file.path(VIS, "heatmap_centroid_correlation.png"), "Centered pseudobulk Pearson r (HVG)",
   breaks = seq(-lim, lim, length.out = 101))

vl <- function(y, ttl) ggplot(cell_df, aes(annotation, .data[[y]], fill = annotation)) +
    geom_violin(scale = "width", linewidth = 0.2) + geom_boxplot(width = 0.12, outlier.shape = NA, fill = "white") +
    scale_fill_manual(values = a_cols, guide = "none") + theme_classic() + ggtitle(ttl) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank())
ggsave(file.path(VIS, "violin_scores_by_annotation.png"),
       vl("prediction_score", "Label-transfer prediction score") / vl("mapping_score", "MappingScore") /
           vl("henry_mixing_ratio", "Henry kNN mixing (obs/exp)"),
       width = 9, height = 12, dpi = 150, bg = "white")

p6 <- ggplot(cell_df, aes(refUMAP_1, refUMAP_2)) +
    geom_point(data = ref_umap, aes(color = paste("Henry", Population)), size = 0.1, alpha = 0.4) +
    scale_color_manual(values = h_cols, name = "Henry reference") +
    geom_point(size = 0.1, color = "black", alpha = 0.5) + facet_wrap(~annotation, ncol = 4) + theme_classic() +
    ggtitle("CRPC cells (black) projected onto Henry reference UMAP") +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(VIS, "refumap_projection_by_annotation.png"), p6, width = 16, height = 11, dpi = 150, bg = "white")
cat("[02] done ->", OUT, "\n")
