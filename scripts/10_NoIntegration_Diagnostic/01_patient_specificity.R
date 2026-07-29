#!/usr/bin/env Rscript
# =============================================================================
# 10 No-Integration Diagnostic — He 2021 식 patient-specificity 검정
# -----------------------------------------------------------------------------
# He et al. 2021 Nat Med (mCRPC scRNA): batch correction 을 쓰지 않고 보면
#   non-malignant 는 환자를 가로질러 co-mingle, cancer 는 환자 특이 클러스터를 형성
#   (patient-specific large somatic CNV 가 원인) → malignancy 신호.
# 우리 상피 임베딩은 Harmony(theta=2) 통합이라 이 신호가 지워져 있다. 여기서는 동일한
#   SCT+PCA 에서 Harmony 만 빼고(=기존 'pca' reduction 재사용) UMAP/클러스터를 다시 만들어,
#   어느 상피 annotation 이 patient-specific(malignant 후보)이고 어느 게 shared(benign
#   후보)인지 진단한다.
# 재사용 근거: epi@reductions$pca 는 Harmony 이전 PCA (SCT 3000 HVG, 3환자 merged,
#   batch 변수/regress 없음) → 그 자체가 "무통합" baseline.
# =============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(patchwork) })
set.seed(42)

OUT  <- "Results/10_NoIntegration_Diagnostic"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
DIMS <- 1:30                                    # 04 파이프라인과 동일

cbP <- c(CRPC1 = "#E69F00", CRPC2 = "#0072B2", CRPC3 = "#009E73")   # 환자 3색(색약친화)
cbA <- c("#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7",
         "#999999","#000000","#332288","#117733","#88CCEE","#AA4499")

epi <- readRDS("Results/05_Epithelial_Downstream/epi_annotated.rds")

# --- Harmony 뺀 임베딩: 기존 pca 재사용 -------------------------------------
epi <- RunUMAP(epi, reduction = "pca", dims = DIMS,
               reduction.name = "umap_noharm", reduction.key = "UMAPnoharm_", verbose = FALSE)
epi <- FindNeighbors(epi, reduction = "pca", dims = DIMS,
                     graph.name = c("nn_noharm", "snn_noharm"), verbose = FALSE)
epi <- FindClusters(epi, graph.name = "snn_noharm", resolution = 0.4, verbose = FALSE)
epi$cluster_noharm <- Idents(epi)

# --- 1) UMAP 비교: Harmony vs no-Harmony (patient & annotation) --------------
gg <- function(r, grp, cols, ttl) DimPlot(epi, reduction = r, group.by = grp, cols = cols,
        shuffle = TRUE, pt.size = 0.2) + ggtitle(ttl)
p <- (gg("umap","orig.ident",cbP,"Harmony · patient") | gg("umap_noharm","orig.ident",cbP,"No-Harmony · patient")) /
     (gg("umap","annotation",cbA,"Harmony · annotation") | gg("umap_noharm","annotation",cbA,"No-Harmony · annotation"))
ggsave(file.path(OUT, "umap_harmony_vs_noharmony.png"), p, width = 16, height = 12, dpi = 150)

# --- 2) patient-specificity 정량 --------------------------------------------
# (a) annotation × patient 구성
comp      <- as.data.frame.matrix(table(epi$annotation, epi$orig.ident))
comp_frac <- round(comp / rowSums(comp), 3)
colnames(comp_frac) <- paste0(colnames(comp_frac), "_frac")
write.csv(cbind(comp, comp_frac), file.path(OUT, "annotation_by_patient.csv"))

# (b) neighbor purity: raw pca 공간에서 각 세포의 kNN 중 같은 환자 비율
emb  <- Embeddings(epi, "pca")[, DIMS]
knn  <- RANN::nn2(emb, k = 31)$nn.idx[, -1]     # self 제외 30-NN
pat  <- as.character(epi$orig.ident)
same <- rowMeans(matrix(pat[knn], ncol = 30) == pat)
pfrac    <- table(pat) / length(pat)
expected <- sum(pfrac^2)                         # 무작위 혼합 기대치
epi$nn_same_patient <- same

pur <- aggregate(nn_same_patient ~ annotation, data = epi@meta.data, mean)
pur$expected_if_mixed <- round(expected, 3)
# segregation index: 0 = 완전 혼합(benign-like), 1 = 완전 환자분리(malignant-like)
pur$segregation_index <- round((pur$nn_same_patient - expected) / (1 - expected), 3)
pur$nn_same_patient   <- round(pur$nn_same_patient, 3)
pur <- pur[order(-pur$segregation_index), ]
write.csv(pur, file.path(OUT, "patient_segregation_by_annotation.csv"), row.names = FALSE)
cat(sprintf("\nexpected same-patient frac if fully mixed = %.3f\n", expected))
cat("=== patient segregation by annotation (He heuristic) ===\n")
print(pur, row.names = FALSE)

# (c) no-Harmony raw cluster × patient / × annotation
rc <- as.data.frame.matrix(table(epi$cluster_noharm, epi$orig.ident))
rc$dominant_frac <- round(apply(rc, 1, max) / rowSums(rc), 3)
write.csv(rc, file.path(OUT, "noharm_cluster_by_patient.csv"))
write.csv(as.data.frame.matrix(table(epi$cluster_noharm, epi$annotation)),
          file.path(OUT, "noharm_cluster_by_annotation.csv"))
cat("\n=== no-Harmony raw clusters: patient composition ===\n"); print(rc)

# --- 3) segregation 을 no-Harmony UMAP 에 투영 ------------------------------
p2 <- FeaturePlot(epi, "nn_same_patient", reduction = "umap_noharm", order = TRUE, pt.size = 0.3) +
      scale_color_viridis_c(option = "viridis", name = "same-patient\nkNN frac") +
      ggtitle("Patient segregation (raw-PCA kNN) on no-Harmony UMAP")
ggsave(file.path(OUT, "umap_noharm_patient_segregation.png"), p2, width = 7.5, height = 6, dpi = 150)

saveRDS(epi@meta.data[, c("orig.ident","annotation","cluster_noharm","nn_same_patient")],
        file.path(OUT, "noharm_meta.rds"))
cat("\n[10] done ->", OUT, "\n")
