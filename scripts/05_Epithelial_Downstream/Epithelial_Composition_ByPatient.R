# Epithelial_Composition_ByPatient.R
# UMAP_annotation_split.png (split.by = orig.ident) 의 수치판.
# 환자(orig.ident) x annotation 별 cell 수와 해당 샘플 내 비중(%)을 표로 저장.
#
# Reads:  Results/05_Epithelial_Downstream/epi_annotated.rds
# Writes: Results/05_Epithelial_Downstream/Annotation/composition_by_patient_long.csv
#         Results/05_Epithelial_Downstream/Annotation/composition_by_patient_n_wide.csv
#         Results/05_Epithelial_Downstream/Annotation/composition_by_patient_pct_wide.csv

suppressMessages({
    library(Seurat)
    library(dplyr)
    library(tidyr)
})

IN_RDS  <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
OUT_DIR <- "Results/05_Epithelial_Downstream/Annotation"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

epi <- readRDS(IN_RDS)

md <- data.frame(
    sample     = as.character(epi$orig.ident),
    cluster    = as.character(epi$seurat_clusters),
    annotation = epi$annotation,
    stringsAsFactors = FALSE
)

lab_levels <- levels(md$annotation)
sample_levels <- sort(unique(md$sample))

# 0-count 조합도 유지하기 위해 complete grid 로 채움
long <- md %>%
    count(sample, annotation, name = "n", .drop = FALSE) %>%
    complete(sample = sample_levels, annotation = lab_levels, fill = list(n = 0)) %>%
    group_by(sample) %>%
    mutate(sample_total = sum(n),
           pct_of_sample = round(100 * n / sample_total, 2)) %>%
    ungroup() %>%
    group_by(annotation) %>%
    mutate(cluster_total = sum(n),
           pct_of_cluster = round(100 * n / cluster_total, 2)) %>%
    ungroup() %>%
    mutate(annotation = factor(annotation, levels = lab_levels)) %>%
    arrange(sample, annotation)

# annotation -> seurat cluster ID 참고용 컬럼
lab2clu <- md %>%
    distinct(annotation, cluster) %>%
    group_by(annotation) %>%
    summarise(seurat_cluster = paste(sort(as.integer(cluster)), collapse = ","),
              .groups = "drop")

long <- long %>%
    left_join(lab2clu, by = "annotation") %>%
    select(sample, annotation, seurat_cluster,
           n, sample_total, pct_of_sample,
           cluster_total, pct_of_cluster)

write.csv(long, file.path(OUT_DIR, "composition_by_patient_long.csv"), row.names = FALSE)

n_wide <- long %>%
    select(annotation, seurat_cluster, sample, n) %>%
    pivot_wider(names_from = sample, values_from = n) %>%
    arrange(annotation)
n_wide$Total <- rowSums(n_wide[, sample_levels, drop = FALSE])
n_wide <- rbind(
    n_wide,
    data.frame(annotation = "Total", seurat_cluster = NA,
               as.list(colSums(n_wide[, c(sample_levels, "Total"), drop = FALSE])),
               check.names = FALSE)
)
write.csv(n_wide, file.path(OUT_DIR, "composition_by_patient_n_wide.csv"), row.names = FALSE)

pct_wide <- long %>%
    select(annotation, seurat_cluster, sample, pct_of_sample) %>%
    pivot_wider(names_from = sample, values_from = pct_of_sample) %>%
    arrange(annotation)
write.csv(pct_wide, file.path(OUT_DIR, "composition_by_patient_pct_wide.csv"), row.names = FALSE)

message("Wrote composition tables to: ", OUT_DIR)
print(as.data.frame(n_wide))
print(as.data.frame(pct_wide))
