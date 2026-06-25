#!/usr/bin/env Rscript
# Export epi_annotated.rds raw counts + minimal metadata to a .loom file for pyscenic.

suppressPackageStartupMessages({
  library(Seurat)
  library(SCopeLoomR)
})

PROJ <- "/home/MMB/projects/Human CRPC scRNAseq"
OUT_DIR <- file.path(PROJ, "Results/05_Epithelial_Downstream/SCENIC")
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
LOOM_PATH <- file.path(OUT_DIR, "epi_counts.loom")

message("Loading epi_annotated.rds")
epi <- readRDS(file.path(PROJ, "Results/05_Epithelial_Downstream/epi_annotated.rds"))
DefaultAssay(epi) <- "RNA"

mat <- GetAssayData(epi, assay = "RNA", layer = "counts")
message("counts dim: ", paste(dim(mat), collapse=" x "))

# Filter genes with at least 3 cells expressing
keep_genes <- rowSums(mat > 0) >= 3
mat <- mat[keep_genes, ]
message("after gene filter: ", paste(dim(mat), collapse=" x "))

meta <- data.frame(
  CellID    = colnames(mat),
  CellType  = as.character(Idents(epi)),
  Patient   = as.character(epi$orig.ident),
  row.names = colnames(mat),
  stringsAsFactors = FALSE
)

if (file.exists(LOOM_PATH)) file.remove(LOOM_PATH)
loom <- build_loom(
  file.name        = LOOM_PATH,
  dgem             = mat,
  title            = "CRPC epithelial annotated",
  default.embedding = NULL
)
# SCopeLoomR renamed the annotation API; add_cell_annotation() no longer exists.
# Use add_col_attr() per column (CellType/Patient are categorical -> as.annotation).
add_col_attr(loom = loom, key = "CellType", value = meta$CellType, as.annotation = TRUE)
add_col_attr(loom = loom, key = "Patient",  value = meta$Patient,  as.annotation = TRUE)
close_loom(loom)
message("wrote ", LOOM_PATH)
