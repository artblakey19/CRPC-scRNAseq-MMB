# Rebuild epi2 with EpiCellTypes annotation
# Reproduces Kimetal_NC2024 lines 570-623 from saved epi.rds

library(Seurat)
library(ggplot2)

# Load saved epithelial subset (pre-cell-cycle-regression)
epi <- readRDS("TG Mouse/TriplevDouble.combined1.epi.rds")

# Cell cycle scoring ----
mouse_cc <- readRDS("TG Mouse/mouse_cell_cycle_genes.rds")
s.genes <- mouse_cc$s.genes
g2m.genes <- mouse_cc$g2m.genes

DefaultAssay(epi) <- "RNA"
epi <- ScaleData(epi, features = rownames(epi))
epi <- CellCycleScoring(epi, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

# Cell cycle regression + reclustering ----
epi2 <- epi
DefaultAssay(epi2) <- "integrated"
epi2 <- ScaleData(epi2,
    vars.to.regress = c("S.Score", "G2M.Score"),
    features = rownames(epi2)
)
epi2 <- RunPCA(epi2, features = VariableFeatures(epi2))
epi2 <- FindNeighbors(epi2, reduction = "pca", dims = 1:22)
epi2 <- FindClusters(epi2, resolution = 2.5)
epi2 <- RunUMAP(epi2, reduction = "pca", dims = 1:22)

# Apply EpiCellTypes annotation ----
Idents(epi2) <- "seurat_clusters"
epi2 <- RenameIdents(
    object = epi2,
    "33" = "BE1",
    "6" = "BE2", "16" = "BE2", "29" = "BE2",
    "15" = "BE3", "5" = "BE3",
    "23" = "BE4", "30" = "BE4", "28" = "BE4", "34" = "BE4",
    "35" = "LE1", "21" = "LE1", "12" = "LE1", "13" = "LE1",
    "25" = "LE2", "10" = "LE2", "17" = "LE2",
    "18" = "LE3", "14" = "LE3", "32" = "LE3", "31" = "LE3", "4" = "LE3", "20" = "LE3",
    "11" = "LE4", "1" = "LE4",
    "0" = "LE5", "9" = "LE5",
    "2" = "LE6", "36" = "LE6",
    "3" = "LE7",
    "22" = "LE8", "8" = "LE8", "19" = "LE8",
    "7" = "LE9", "24" = "LE9",
    "26" = "UrLE",
    "27" = "OE"
)
epi2[["EpiCellTypes"]] <- Idents(epi2)

cat("EpiCellTypes distribution:\n")
print(table(epi2$EpiCellTypes))

# Paper color palette (SFig.4f)
epi_cols <- c(
    "salmon", "skyblue1", "olivedrab2", "brown3",
    "deeppink1", "blue", "darkorange", "red", "turquoise3",
    "bisque3", "slategray3", "mediumorchid3",
    "yellow2", "green4", "black"
)

dir.create("TG Mouse/Results", showWarnings = FALSE)

# UMAP colored by EpiCellTypes
p1 <- DimPlot(epi2,
    reduction = "umap", group.by = "EpiCellTypes",
    pt.size = 0.3, label = TRUE, label.size = 4, cols = epi_cols
) +
    ggtitle("EpiCellTypes") +
    theme(plot.title = element_text(hjust = 0.5))

# UMAP split by genotype (DoubleTg vs TripleTg)
p2 <- DimPlot(epi2,
    reduction = "umap", group.by = "EpiCellTypes",
    split.by = "stim", pt.size = 0.3, label = TRUE,
    label.size = 3.5, cols = epi_cols
) +
    ggtitle("EpiCellTypes by Genotype") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("TG Mouse/Results/epi2_EpiCellTypes_UMAP.png",
    p1,
    width = 8, height = 7, dpi = 300
)
ggsave("TG Mouse/Results/epi2_EpiCellTypes_split_UMAP.png",
    p2,
    width = 14, height = 7, dpi = 300
)

# Marker DotPlot (SFig.4h)
DefaultAssay(epi2) <- "RNA"
Idents(epi2) <- "EpiCellTypes"

marker_genes <- c(
    "Tmem171", "Egfl6", "Ncam1", "Clca3a2", "Adm",
    "Krt15", "Palld", "Ctsl", "Tubb6", "Tpm1",
    "Aqp3", "Lgals7", "Col17a1", "Lamb3", "Pvrl1",
    "Sncg", "Ifi202b", "Gpnmb", "Dapl1", "Gpr87",
    "Lars2", "AY036118", "Gm42418", "Gm26917", "Hbb-bs",
    "Defa21", "Gm15293", "Defa5", "Cryba4", "Ltf",
    "Rnf149", "Car2", "Lap3", "Plat", "Tm4sf1",
    "Coch", "Tgfb2", "Dkk2", "Zeb2", "Apoc1",
    "Crip1", "Tspan8", "Ly6c1", "Btc", "B2m",
    "Tgm4", "9530053A07Rik", "Gm5615", "Man1a", "Spink8",
    "Msmb", "Mme", "Apof", "Pcp4", "Agtr1a",
    "Pigr", "Tspan1", "Tnfrsf21", "Dcxr", "Cldn3",
    "Spink1", "Sbpl", "Crabp1", "Col6a3", "Gucy2g",
    "Gsdmc2", "Gsdmc3", "Barx2", "Cxcl15", "Krt4",
    "Serping1", "Igfbp6", "Fbln1", "Serpinf1", "Col1a2"
)
# Keep only genes present in the data
marker_genes <- marker_genes[marker_genes %in% rownames(epi2)]

p3 <- DotPlot(epi2,
    features = marker_genes,
    cols = c("light grey", "red")
) + RotatedAxis() +
    ggtitle("EpiCellTypes Markers") +
    theme(plot.title = element_text(hjust = 0.5))

ggsave("TG Mouse/Results/epi2_EpiCellTypes_markers_DotPlot.png",
    p3,
    width = 20, height = 6, dpi = 300
)

# Save ----
saveRDS(epi2, "TG Mouse/TriplevDouble.combined1.epi2.rds")
cat("\nSaved: TG Mouse/TriplevDouble.combined1.epi2.rds\n")
cat("Figures saved to: TG Mouse/Results/\n")
