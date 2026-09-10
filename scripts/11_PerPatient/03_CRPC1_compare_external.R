#!/usr/bin/env Rscript
# Stage 11-03 — CRPC1 환자 단독 분석(우리 QC) vs 외부 재주석(reannotated_CRPC1.rds) 비교
#
# 외부 객체는 QC 기준(nCount≥500 & nFeature≥400, GEM-X scDblFinder, decontX 없음)이 우리와 다르므로
# (1) 셀 집합 차이, (2) 공통 셀에서 대분류 일치도, (3) 외부 major_celltype 을 우리 단독 UMAP 에 투영,
# (4) 우리 상피 subtype 별 외부 라벨 분포를 낸다.
#
# 입력: Results/11_PerPatient/CRPC1/01_Major/CRPC1_major.rds (celltype 부여됨)
#       Results/11_PerPatient/CRPC1/02_Epithelial/CRPC1_epi.rds (annotation 부여됨; 없으면 상피 비교 생략)
#       reannotated_CRPC1.rds
# 출력: Results/11_PerPatient/CRPC1/03_vs_External/

suppressMessages({ library(Seurat); library(dplyr); library(ggplot2); library(patchwork) })
source("scripts/00_utils/scRNA_utils.R")

OUT_DIR <- "Results/11_PerPatient/CRPC1/03_vs_External"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

major <- readRDS("Results/11_PerPatient/CRPC1/01_Major/CRPC1_major.rds")
ext   <- readRDS("reannotated_CRPC1.rds")

our_bc <- sub("^P1_", "", colnames(major))
ext_bc <- colnames(ext)
common <- intersect(our_bc, ext_bc)
cat(sprintf("our %d / external %d / common %d / our-only %d / external-only %d\n",
            length(our_bc), length(ext_bc), length(common),
            length(setdiff(our_bc, ext_bc)), length(setdiff(ext_bc, our_bc))))
writeLines(c(
    sprintf("our cells (decontX QC): %d", length(our_bc)),
    sprintf("external cells: %d", length(ext_bc)),
    sprintf("common: %d", length(common)),
    sprintf("our-only: %d", length(setdiff(our_bc, ext_bc))),
    sprintf("external-only: %d", length(setdiff(ext_bc, our_bc)))
), file.path(OUT_DIR, "cell_overlap.txt"))

# 외부 라벨을 우리 객체에 붙임
m <- match(our_bc, ext_bc)
major$ext_major_celltype <- as.character(ext$major_celltype)[m]
major$ext_compartment    <- as.character(ext$compartment)[m]

# 셀 집합 차이의 QC 특성 (우리 객체 기준 our-only 셀)
qc_only <- major@meta.data %>%
    mutate(in_external = !is.na(ext_major_celltype)) %>%
    group_by(in_external) %>%
    summarise(n = n(), median_nCount = median(nCount_RNA), median_nFeature = median(nFeature_RNA),
              median_mt = round(median(percent.mt), 2),
              median_decontX = round(median(decontX_contamination), 3), .groups = "drop")
write.csv(qc_only, file.path(OUT_DIR, "our_cells_QC_by_external_membership.csv"), row.names = FALSE)
# 외부-only 셀의 QC (외부 객체 기준) + 외부 라벨
ext_only <- setdiff(ext_bc, our_bc)
ext_only_tab <- ext@meta.data[ext_only, ] %>%
    group_by(major_celltype) %>%
    summarise(n = n(), median_nCount = median(nCount_RNA), median_nFeature = median(nFeature_RNA),
              median_mt = round(median(percent.mt), 2), .groups = "drop") %>% arrange(-n)
write.csv(ext_only_tab, file.path(OUT_DIR, "external_only_cells_by_label.csv"), row.names = FALSE)

# 공통 셀 대분류 교차표
tab <- as.data.frame.matrix(table(ours = major$celltype, external = addNA(factor(major$ext_major_celltype))))
write.csv(tab, file.path(OUT_DIR, "celltype_ours_vs_external.csv"))
print(tab)

# 외부 라벨을 우리 단독 UMAP 에 투영
n1 <- nlevels(major$celltype); n2 <- length(unique(na.omit(major$ext_major_celltype)))
p1 <- DimPlot(major, group.by = "celltype", label = TRUE, repel = TRUE, pt.size = 0.2,
              cols = utils_cb_palette(n1)) + ggtitle("CRPC1 alone — our annotation")
p2 <- DimPlot(major, group.by = "ext_major_celltype", label = TRUE, repel = TRUE, pt.size = 0.2,
              cols = utils_cb_palette(n2), na.value = "grey85") +
    ggtitle("External (reannotated_CRPC1) major_celltype; grey = not in external")
ggsave(file.path(OUT_DIR, "UMAP_ours_vs_external.png"), plot = p1 + p2, width = 22, height = 8, bg = "white")

# 상피 subtype 별 외부 라벨 (외부는 상피 subtype 없음 → Epithelial 인지 / 누락인지만 의미)
EPI_RDS <- "Results/11_PerPatient/CRPC1/02_Epithelial/CRPC1_epi.rds"
if (file.exists(EPI_RDS)) {
    epi <- readRDS(EPI_RDS)
    if ("annotation" %in% colnames(epi@meta.data)) {
        me <- match(sub("^P1_", "", colnames(epi)), ext_bc)
        epi$ext_major_celltype <- as.character(ext$major_celltype)[me]
        tab_e <- epi@meta.data %>%
            group_by(annotation) %>%
            summarise(n = n(), in_external = sum(!is.na(ext_major_celltype)),
                      ext_Epithelial = sum(ext_major_celltype == "Epithelial", na.rm = TRUE),
                      ext_other = sum(!is.na(ext_major_celltype) & ext_major_celltype != "Epithelial"),
                      pct_missing_in_external = round(100 * mean(is.na(ext_major_celltype)), 1),
                      .groups = "drop")
        write.csv(tab_e, file.path(OUT_DIR, "epi_annotation_vs_external.csv"), row.names = FALSE)
        print(as.data.frame(tab_e))
    }
}
message("→ ", OUT_DIR)
