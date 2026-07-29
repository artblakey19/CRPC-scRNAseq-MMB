#!/usr/bin/env Rscript
# =============================================================================
# Numbat step 3  —  UMAP 시각화 (HOST / renv, Seurat side)
# -----------------------------------------------------------------------------
# 02b 가 환자별로 낸 numbat 셀 단위 사후확률(clone_post_*.tsv)을 상피 UMAP 에
# 얹어 환자별 split 으로 그린다. numbat 패키지는 필요 없고 TSV 만 읽는다.
#
# 입력 :
#   Results/05_Epithelial_Downstream/epi_annotated.rds        (상피 UMAP + annotation)
#   Results/08_Numbat/<S>/numbat/clone_post_<최종iter>.tsv    (셀 단위 numbat 결과)
# 출력 : Results/08_Numbat/Visualization/
#   UMAP_Numbat_byPatient.png              마스터 그리드(annotation+compartment+p_cnv+clone)
#   UMAP_annotation_byPatient_labeled.png  라벨 붙은 annotation UMAP(환자별)
#   UMAP_Numbat_compartment_byPatient.png  tumor/normal(aneuploid 판정)
#   UMAP_Numbat_pcnv_byPatient.png         p(CNV) 연속형
#   Numbat_compartment_by_annotation.csv   클러스터×compartment 교차표(환자별)
#
# 매핑 : epi 바코드(P1_/P2_/P3_ 접두사)에서 접두사를 떼고 orig.ident 로 환자를 지정,
#   해당 환자 numbat 파일의 cell 과 매칭. numbat query = Epithelial \ Ionocyte 이므로
#   epi 의 Ionocyte 세포는 numbat 결과가 없다 → 회색(n/a)으로 표시(설계상 정상).
#
# 색 규칙(메모리) : discrete = 색약 친화(Okabe-Ito), 연속형 = viridis.
#   compartment normal=파랑 / tumor=주황(vermillion) / n/a(Ionocyte)=회색.
# =============================================================================
suppressMessages({
    library(Seurat); library(ggplot2); library(patchwork)
    library(viridis); library(ggrepel)
})
source("scripts/00_utils/scRNA_utils.R")

EPI_RDS <- "Results/05_Epithelial_Downstream/epi_annotated.rds"
NB_DIR  <- "Results/08_Numbat"
OUT_DIR <- file.path(NB_DIR, "Visualization")
SAMPLES <- c("CRPC1", "CRPC2", "CRPC3")
PT      <- 0.6            # 점을 진하게(크게/불투명) — 매니폴드가 꽉 차 보이도록
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# annotation 라벨 순서(05 와 동일) → 환자 간 색을 고정하기 위한 named 팔레트
LABEL_LEVELS <- c("LE(ARPC)", "Club", "Hillock 1", "Hillock 2",
                  "BE 1", "BE 2", "BE 3", "BE 4", "BE 5", "BE 6",
                  "OE", "Ionocyte")
ANNO_COLS <- setNames(utils_cb_palette(length(LABEL_LEVELS)), LABEL_LEVELS)
COMP_COLS <- c(normal = "grey70", tumor = "#CC3311")    # 회색(diploid) / 빨강(aneuploid)
NA_GREY   <- "grey80"     # p(CNV)·clone 패널의 n/a(Ionocyte) 배경

# --- 1. 상피 UMAP + annotation 을 data.frame 으로 -----------------------------
message("Reading ", EPI_RDS, " ...")
epi <- readRDS(EPI_RDS)
emb <- Embeddings(epi, "umap")
df  <- data.frame(
    UMAP_1     = emb[, 1],
    UMAP_2     = emb[, 2],
    sample     = factor(as.character(epi$orig.ident), levels = SAMPLES),
    annotation = factor(as.character(epi$annotation), levels = LABEL_LEVELS),
    row.names  = colnames(epi),
    stringsAsFactors = FALSE
)
df$bc <- sub("^P[0-9]+_", "", rownames(df))     # numbat cell 과 맞추려 접두사 제거
rm(epi); gc()

# --- 2. numbat 셀 단위 결과(최종 iteration) 조인 -----------------------------
df$nb_compartment <- NA_character_
df$nb_pcnv        <- NA_real_
df$nb_clone       <- NA_integer_

pick_final <- function(S) {              # clone_post_1/2/... 중 최종(최댓값) 선택
    fs  <- list.files(file.path(NB_DIR, S, "numbat"),
                      pattern = "^clone_post_[0-9]+\\.tsv$", full.names = TRUE)
    if (!length(fs)) stop("No clone_post_*.tsv for ", S)
    it  <- as.integer(sub(".*clone_post_([0-9]+)\\.tsv$", "\\1", fs))
    fs[which.max(it)]
}

for (S in SAMPLES) {
    f  <- pick_final(S)
    cp <- read.delim(f, stringsAsFactors = FALSE)
    idx <- which(df$sample == S)
    m   <- match(df$bc[idx], cp$cell)
    df$nb_compartment[idx] <- cp$compartment_opt[m]
    df$nb_pcnv[idx]        <- cp$p_cnv[m]
    df$nb_clone[idx]       <- cp$clone_opt[m]
    message(sprintf("%s: %s | matched %d / %d epi cells (%d n/a = Ionocyte/QC)",
                    S, basename(f), sum(!is.na(m)), length(idx), sum(is.na(m))))
}
df$nb_compartment <- factor(df$nb_compartment, levels = c("normal", "tumor"))

# --- 3. 공통 UMAP 프레임 & 테마 ----------------------------------------------
xr <- range(df$UMAP_1) + 0.03 * c(-1, 1) * diff(range(df$UMAP_1))
yr <- range(df$UMAP_2) + 0.03 * c(-1, 1) * diff(range(df$UMAP_2))

finish <- function(p, ttl) {
    p + coord_equal(xlim = xr, ylim = yr) + ggtitle(ttl) +
        theme_classic(base_size = 11) +
        theme(plot.title   = element_text(size = 11, face = "bold"),
              axis.title   = element_text(size = 8),
              axis.text    = element_text(size = 7),
              legend.title = element_text(size = 9),
              legend.text  = element_text(size = 8),
              legend.key.size = unit(0.4, "cm"))
}

# --- 4. 패널 빌더 -------------------------------------------------------------
# annotation 참조(라벨 = 클러스터 중심). 각 항목 그림 맨 왼쪽에 붙는 "지도" 1장.
#   cells = NULL 이면 전 세포(환자 통합), 특정 환자면 그 환자만.
panel_anno <- function(cells = NULL) {
    d   <- if (is.null(cells)) df else df[df$sample == cells, ]
    d$annotation <- droplevels(d$annotation)
    cen <- aggregate(cbind(UMAP_1, UMAP_2) ~ annotation, data = d, FUN = median)
    p <- ggplot(d, aes(UMAP_1, UMAP_2, colour = annotation)) +
        geom_point(size = PT, stroke = 0) +
        scale_colour_manual(values = ANNO_COLS, drop = TRUE) +
        geom_text_repel(data = cen, aes(label = annotation), colour = "black",
                        size = 3.2, fontface = "bold", seed = 42,
                        bg.colour = "white", bg.r = 0.15,
                        min.segment.length = 0, max.overlaps = Inf,
                        box.padding = 0.3, point.size = NA) +
        guides(colour = "none")
    finish(p, if (is.null(cells)) "annotation (reference)" else paste0(cells, " — annotation"))
}

# n/a(Ionocyte) 세포는 배경 레이어로 먼저 깔아 매니폴드를 온전히 보이게 하고,
# 색 스케일에는 실제 레벨만 남긴다(=> NA 제거 경고 없음). col 로 회색 톤 지정.
bg_na <- function(d, col = NA_GREY) geom_point(data = d[is.na(d$nb_compartment), ],
                                aes(UMAP_1, UMAP_2), colour = col,
                                size = PT, stroke = 0, inherit.aes = FALSE)

# compartment (tumor = aneuploid / normal = diploid / n/a = Ionocyte) ---------
#   n/a 는 normal 과 같은 회색으로 → 그림 전체가 "회색 vs 빨강" 두 톤으로 읽힌다.
panel_comp <- function(s) {
    d  <- df[df$sample == s, ]
    fg <- d[!is.na(d$nb_compartment), ]
    fg <- fg[order(match(as.character(fg$nb_compartment), c("normal", "tumor"))), ]  # tumor 위로
    p <- ggplot() +
        bg_na(d, col = COMP_COLS[["normal"]]) +
        geom_point(data = fg, aes(UMAP_1, UMAP_2, colour = nb_compartment),
                   size = PT, stroke = 0) +
        scale_colour_manual(values = COMP_COLS, name = "Numbat")
    finish(p, paste0(s, " — tumor / normal"))
}

# p_cnv (aneuploid 사후확률, 연속형 viridis) ----------------------------------
panel_pcnv <- function(s) {
    d  <- df[df$sample == s, ]
    fg <- d[!is.na(d$nb_pcnv), ]
    fg <- fg[order(fg$nb_pcnv), ]                        # p 큰 세포를 위로
    p <- ggplot() +
        bg_na(d) +
        geom_point(data = fg, aes(UMAP_1, UMAP_2, colour = nb_pcnv),
                   size = PT, stroke = 0) +
        scale_colour_viridis_c(option = "viridis", limits = c(0, 1), name = "p(CNV)")
    finish(p, paste0(s, " — p(CNV)"))
}

# clone (환자별 subclone; 번호가 환자마다 독립이라 패널별 범례) ----------------
panel_clone <- function(s) {
    d  <- df[df$sample == s, ]
    fg <- d[!is.na(d$nb_clone), ]
    fg$nb_clone <- factor(fg$nb_clone, levels = sort(unique(fg$nb_clone)))
    cols <- setNames(utils_cb_palette(nlevels(fg$nb_clone)), levels(fg$nb_clone))
    p <- ggplot() +
        bg_na(d) +
        geom_point(data = fg[order(fg$nb_clone), ],
                   aes(UMAP_1, UMAP_2, colour = nb_clone), size = PT, stroke = 0) +
        scale_colour_manual(values = cols, name = "clone")
    finish(p, paste0(s, " — clone"))
}

# --- 5. 그림 저장 : 항목별 4장(annotation 참조 + 환자 3장) 한 plot -----------
cap_txt <- paste0("Numbat query = Epithelial \\ Ionocyte; 회색 = normal(diploid)/",
                  "Ionocyte(numbat 미검사). 왼쪽 = annotation 지도(라벨).")

# 오래된(구조가 다른) 산출물 제거 — 사용자 요청대로 대체
invisible(suppressWarnings(file.remove(file.path(OUT_DIR,
    c("UMAP_Numbat_byPatient.png", "UMAP_annotation_byPatient_labeled.png")))))

anno_ref <- panel_anno(NULL)     # 환자 통합 라벨 지도 1장 (모든 그림 맨 왼쪽 공통)

# 항목 정의: 각 항목 = [annotation 지도 | CRPC1 | CRPC2 | CRPC3], guides 는 항목별로 처리
ITEMS <- list(
    list(file = "UMAP_Numbat_compartment_byPatient.png",
         ttl  = "Numbat tumor / normal — patient-split (+ annotation)",
         fn   = panel_comp,  collect = TRUE),
    list(file = "UMAP_Numbat_pcnv_byPatient.png",
         ttl  = "Numbat p(CNV) — patient-split (+ annotation)",
         fn   = panel_pcnv,  collect = TRUE),
    list(file = "UMAP_Numbat_clone_byPatient.png",
         ttl  = "Numbat clone — patient-split (+ annotation)",
         fn   = panel_clone, collect = FALSE)   # clone 번호는 환자마다 독립 → 범례 개별
)

for (it in ITEMS) {
    panels <- c(list(anno_ref), lapply(SAMPLES, it$fn))    # 4장: 지도 + 3환자
    fig <- wrap_plots(panels, nrow = 1)
    if (it$collect) fig <- fig + plot_layout(guides = "collect")
    fig <- fig + plot_annotation(
        title   = it$ttl,
        caption = cap_txt,
        theme   = theme(plot.title   = element_text(size = 16, face = "bold"),
                        plot.caption = element_text(size = 9, colour = "grey35")))
    ggsave(file.path(OUT_DIR, it$file), fig,
           width = 18, height = 5.2, dpi = 200, bg = "white", limitsize = FALSE)
    message("Wrote: ", it$file)
}

# --- 6. annotation × compartment 교차표(환자별) → CSV + 콘솔 ------------------
grid <- expand.grid(sample = SAMPLES, annotation = LABEL_LEVELS,
                    stringsAsFactors = FALSE)
cnt <- function(comp) mapply(function(s, a) {
    sel <- df$sample == s & df$annotation == a
    if (is.na(comp)) sum(sel & is.na(df$nb_compartment))
    else             sum(sel & !is.na(df$nb_compartment) & df$nb_compartment == comp)
}, grid$sample, grid$annotation)
grid$n_normal <- cnt("normal")
grid$n_tumor  <- cnt("tumor")
grid$n_na_iono <- cnt(NA)
grid$n_total  <- grid$n_normal + grid$n_tumor + grid$n_na_iono
grid$pct_tumor <- ifelse(grid$n_total > 0, round(100 * grid$n_tumor / grid$n_total, 1), NA)
grid <- grid[grid$n_total > 0, ]                        # 해당 환자에 없는 클러스터는 제외
grid <- grid[order(grid$sample, match(grid$annotation, LABEL_LEVELS)), ]
write.csv(grid, file.path(OUT_DIR, "Numbat_compartment_by_annotation.csv"),
          row.names = FALSE)
message("Wrote: Numbat_compartment_by_annotation.csv")
message("\n=== tumor fraction by patient (Numbat compartment) ===")
print(round(100 * tapply(df$nb_compartment == "tumor", df$sample, mean, na.rm = TRUE), 1))
message("\n[03] done -> ", OUT_DIR)
