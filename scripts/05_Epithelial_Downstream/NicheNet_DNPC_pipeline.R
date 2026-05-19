#!/usr/bin/env Rscript
# NicheNet: which microenvironment-derived ligands explain the Hillock/Club DNPC terminal state?
# Sender:   non-epi cell types from the full integrated object
# Receiver: Hillock-like 1, Hillock-like 2, Club-like (DNPC terminals)
# Reference (within receiver definition): BE/OE/ARPC = "non-DNPC" epithelial baseline
#
# Strategy: merge non-epi cells (from full) with the epi annotation (from epi_annotated),
# then call nichenet_seuratobj_aggregate per receiver vs the rest of epi.

suppressPackageStartupMessages({
  library(Seurat)
  library(nichenetr)
  library(tidyverse)
})

PROJ <- "/home/MMB/projects/Human CRPC scRNAseq"
NICHENET_DIR <- file.path(PROJ, "Resources/NicheNet")
OUT_DIR <- file.path(PROJ, "Results/05_Epithelial_Downstream/NicheNet")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# 1) Load prior networks --------------------------------------------------
message("[1/6] Loading NicheNet prior networks")
ligand_target_matrix <- readRDS(file.path(NICHENET_DIR, "ligand_target_matrix_nsga2r_final.rds"))
lr_network           <- readRDS(file.path(NICHENET_DIR, "lr_network_human_21122021.rds"))
weighted_networks    <- readRDS(file.path(NICHENET_DIR, "weighted_networks_nsga2r_final.rds"))
lr_network <- lr_network %>% distinct(from, to)

# 2) Build combined object ------------------------------------------------
message("[2/6] Building combined sender(non-epi) + receiver(epi) object")
full <- readRDS(file.path(PROJ, "Results/01_Integrated/combined_CRPC.rds"))
epi  <- readRDS(file.path(PROJ, "Results/05_Epithelial_Downstream/epi_annotated.rds"))

# Non-epi subset from full
nonepi <- subset(full, subset = celltype != "Epithelial")
nonepi$nichenet_label <- as.character(nonepi$celltype)

# Epi subset keeps annotated label
epi$nichenet_label <- as.character(Idents(epi))

# Strip non-RNA assays (SCT merge is brittle across differently-processed objects)
strip_to_rna <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  for (a in setdiff(Assays(obj), "RNA")) obj[[a]] <- NULL
  obj
}
nonepi <- strip_to_rna(nonepi)
epi    <- strip_to_rna(epi)

common_meta <- intersect(colnames(nonepi@meta.data), colnames(epi@meta.data))
nonepi@meta.data <- nonepi@meta.data[, common_meta, drop = FALSE]
epi@meta.data    <- epi@meta.data[,    common_meta, drop = FALSE]

combined <- merge(nonepi, epi, add.cell.ids = c("nonepi","epi"))
DefaultAssay(combined) <- "RNA"
# Ensure data slot exists (NicheNet uses normalized data for receptor expression)
if (length(Layers(combined, assay = "RNA", search = "data")) == 0) {
  combined <- NormalizeData(combined, assay = "RNA", verbose = FALSE)
}
# Join layers if v5 multi-layer
if (length(Layers(combined, assay = "RNA", search = "counts")) > 1) {
  combined <- JoinLayers(combined, assay = "RNA")
}
Idents(combined) <- combined$nichenet_label
message("  combined ncells: ", ncol(combined),
        "   label table:"); print(table(Idents(combined)))

# 3) Define sender / receiver groups --------------------------------------
sender_celltypes <- c("Fibroblast","Endothelial","T/NK cells","Smooth muscle cells",
                      "Stromal","Phagocytes","Immune cells","Mast cells")
receiver_groups  <- list(
  Hillock1 = "Hillock-like 1",
  Hillock2 = "Hillock-like 2",
  Club     = "Club-like"
)
# Reference for DEG = non-DNPC epi
reference_epi <- c("BE 1","BE 2","OE 1","OE 2","OE 3","OE 4","ARPC")

# 4) Run nichenet per receiver -------------------------------------------
all_ligand_activities <- list()
all_target_links      <- list()

save_plot <- function(obj, path, w, h) {
  if (is.null(obj)) return(invisible(NULL))
  tryCatch(ggsave(path, obj, width = w, height = h),
           error = function(e) message("  (skipped ", basename(path), ": ", conditionMessage(e), ")"))
}

for (rname in names(receiver_groups)) {
  receiver <- receiver_groups[[rname]]
  message(sprintf("\n[3/6] NicheNet for receiver_affected=%s  vs  receiver_reference=BE/OE/ARPC", receiver))

  res <- tryCatch(
    nichenet_seuratobj_cluster_de(
      seurat_obj          = combined,
      receiver_affected   = receiver,
      receiver_reference  = reference_epi,
      sender              = sender_celltypes,
      ligand_target_matrix = ligand_target_matrix,
      lr_network          = lr_network,
      weighted_networks   = weighted_networks,
      expression_pct      = 0.10,
      lfc_cutoff          = 0.25,
      geneset             = "up",
      assay_oi            = "RNA"
    ),
    error = function(e) { message("  !! ", conditionMessage(e)); NULL }
  )
  if (is.null(res)) next

  saveRDS(res, file.path(OUT_DIR, sprintf("nichenet_%s.rds", rname)))

  if (!is.null(res$ligand_activities)) {
    la <- res$ligand_activities %>% mutate(receiver = rname)
    all_ligand_activities[[rname]] <- la
  }
  for (slot in c("ligand_target_df","ligand_target_links_df")) {
    if (!is.null(res[[slot]])) {
      all_target_links[[rname]] <- res[[slot]] %>% mutate(receiver = rname)
      break
    }
  }

  message("  saving figures...")
  save_plot(res$ligand_activity_target_heatmap,
            file.path(OUT_DIR, sprintf("%s_ligand_activity_target_heatmap.pdf", rname)), 14, 10)
  save_plot(res$ligand_receptor_heatmap,
            file.path(OUT_DIR, sprintf("%s_ligand_receptor_heatmap.pdf", rname)), 10, 10)
  save_plot(res$ligand_expression_dotplot,
            file.path(OUT_DIR, sprintf("%s_ligand_sender_dotplot.pdf", rname)), 10, 10)
  save_plot(res$ligand_differential_expression_heatmap,
            file.path(OUT_DIR, sprintf("%s_ligand_diff_expr_heatmap.pdf", rname)), 10, 10)
}

# 5) Cross-receiver comparison table -------------------------------------
message("\n[5/6] Writing combined ligand activity table")
la_all <- bind_rows(all_ligand_activities)
write_tsv(la_all, file.path(OUT_DIR, "ligand_activities_all.tsv"))
if (length(all_target_links)) {
  tl_all <- bind_rows(all_target_links)
  write_tsv(tl_all, file.path(OUT_DIR, "ligand_target_links_all.tsv"))
}

# Top-20 ligand by aupr_corrected, per receiver
top20 <- la_all %>%
  group_by(receiver) %>%
  arrange(desc(aupr_corrected), .by_group = TRUE) %>%
  slice_head(n = 20) %>%
  ungroup()
write_tsv(top20, file.path(OUT_DIR, "top20_ligands_per_receiver.tsv"))

# 6) Compact summary heatmap: top-20 ligands × receivers ------------------
message("[6/6] Compact summary heatmap")
top20_long <- la_all %>%
  filter(test_ligand %in% unique(top20$test_ligand)) %>%
  select(receiver, test_ligand, aupr_corrected)

p_sum <- ggplot(top20_long, aes(receiver, reorder(test_ligand, aupr_corrected),
                                fill = aupr_corrected)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "magma") +
  theme_minimal() +
  labs(x = NULL, y = "Ligand", fill = "AUPR\n(corrected)",
       title = "Top NicheNet ligand activity per DNPC receiver")
ggsave(file.path(OUT_DIR, "summary_top20_ligands_heatmap.pdf"), p_sum,
       width = 6, height = 10)

message("\nDone. Outputs in: ", OUT_DIR)
