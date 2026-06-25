# _dnpc_evidence_probe.R  (adversarial-verification evidence pack)
# 논쟁적 module의 (a) 유전자별 cluster 평균 발현, (b) He-KRT7 basal-negative 기준 충족,
# (c) POU2F3 tuft 희소 분획 공동발현을 정량화. 기존 캐시 객체 사용.
suppressMessages({ library(Seurat); library(Matrix) })
epi <- readRDS("Results/05_Epithelial_Downstream/epi_publishedDNPC_programs.rds")
DefaultAssay(epi) <- "RNA"
if (!"data" %in% Layers(epi[["RNA"]])) epi <- NormalizeData(epi, verbose = FALSE)
a <- read.csv("Results/05_Epithelial_Downstream/Annotation/annotation_per_cell.csv", stringsAsFactors = FALSE)
lv <- readLines("Results/05_Epithelial_Downstream/Annotation/label_levels.txt")
cl <- factor(a$annotation[match(colnames(epi), a$cell)], levels = lv)
D  <- GetAssayData(epi, assay = "RNA", layer = "data")
OUT <- "Results/05_Epithelial_Downstream/PublishedDNPC_Programs"

gene_mean_by_cluster <- function(genes) {
    genes <- intersect(genes, rownames(D))
    sub <- as.matrix(D[genes, , drop = FALSE])
    t(apply(sub, 1, function(v) tapply(v, cl, mean)))
}
gene_pct_by_cluster <- function(genes) {
    genes <- intersect(genes, rownames(D))
    sub <- as.matrix(D[genes, , drop = FALSE]) > 0
    t(apply(sub, 1, function(v) tapply(v, cl, mean))) * 100
}

contested <- list(
    KRT7_DN   = c("KRT7","S100A9","S100A8","S100A4","S100A6","S100A10","S100A11"),
    Tang_SCL  = c("FOSL1","FOS","FOSB","JUNB","JUN","ATF3","CD44","TACSTD2","CYR61","CTGF","AXL","CCND1","KRT5","KRT14","KRT17","TP63"),
    AP1_IEG   = c("FOS","JUN","JUNB","EGR1","DUSP1","NR4A1","ZFP36"),
    Lab_DNPCp = c("CXCL8","CXCR1","RUNX2","TGFB1","TGFB2"),
    POU2F3    = c("POU2F3","GFI1B","TRPM5","ASCL2","AVIL","LRMP","PLCB2","SOX9"),
    MYC_core  = c("NPM1","NCL","PCNA","ODC1","EIF4E","SRM","CCT2"))
for (nm in names(contested)) {
    write.csv(round(gene_mean_by_cluster(contested[[nm]]), 3),
              file.path(OUT, paste0("probe_genemean_", nm, ".csv")))
}

# (b) He KRT7-subtype 기준: KRT7+ 이면서 basal(KRT5/KRT14/TP63) 음성 세포 분율
get <- function(g) if (g %in% rownames(D)) as.numeric(D[g, ]) else rep(0, ncol(D))
krt7 <- get("KRT7") > 0
basal_neg <- get("KRT5") == 0 & get("KRT14") == 0 & get("TP63") == 0
he_krt7_subtype <- krt7 & basal_neg
df_he <- data.frame(
    cluster   = levels(cl),
    n         = as.integer(table(cl)),
    pct_KRT7pos        = round(100 * tapply(krt7, cl, mean), 1),
    pct_KRT7pos_basalNEG = round(100 * tapply(he_krt7_subtype, cl, mean), 1))
write.csv(df_he, file.path(OUT, "probe_He_KRT7_basalNeg_criterion.csv"), row.names = FALSE)

# (c) POU2F3 tuft 희소 분획: POU2F3+ 및 POU2F3+ ∧ (GFI1B+|TRPM5+) 공동발현 세포 수
pou <- get("POU2F3") > 0
gfi <- get("GFI1B") > 0; trpm5 <- get("TRPM5") > 0
tuft_copos <- pou & (gfi | trpm5)
df_tuft <- data.frame(
    cluster = levels(cl), n = as.integer(table(cl)),
    n_POU2F3pos = as.integer(tapply(pou, cl, sum)),
    n_POU2F3_and_GFI1BorTRPM5 = as.integer(tapply(tuft_copos, cl, sum)))
write.csv(df_tuft, file.path(OUT, "probe_POU2F3_tuft_copositive.csv"), row.names = FALSE)

cat("=== He KRT7+ basal-NEGATIVE criterion (%) ===\n"); print(df_he)
cat("\n=== POU2F3 tuft co-positive cell counts ===\n"); print(df_tuft)
cat("\ntotal POU2F3+ cells:", sum(pou), " | POU2F3+ & (GFI1B+|TRPM5+):", sum(tuft_copos), "\n")
cat("\n=== KRT7_DN per-gene cluster mean ===\n"); print(round(gene_mean_by_cluster(contested$KRT7_DN), 3))
cat("\n=== Tang_SCL per-gene cluster mean ===\n"); print(round(gene_mean_by_cluster(contested$Tang_SCL), 3))
cat("\n=== Lab_DNPCp per-gene cluster mean ===\n"); print(round(gene_mean_by_cluster(contested$Lab_DNPCp), 3))
