#WT ARKO-Gli1 three timepoint SCseq Work Flow

#### Add necessary tools to library ####

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

####E18_CtrlvARKO.combined####
Idents(object = E18_CtrlvARKO.combined) <- "CellTypes"
DimPlot(E18_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3)

Idents(object = E18_CtrlvARKO.combined) <- "stim"
E18_Ctrlonly <- subset(E18_CtrlvARKO.combined, idents = c("E18_Ctrl"))

Idents(object = E18_Ctrlonly) <- "CellTypes"
DimPlot(E18_Ctrlonly, reduction = "umap", pt.size = 0.3)

E18_ARKOonly <- subset(E18_CtrlvARKO.combined, idents = c("E18_ARKO"))
Idents(object = E18_ARKOonly) <- "CellTypes"
DimPlot(E18_ARKOonly, reduction = "umap", pt.size = 0.3)

Idents(object = E18_CtrlvARKO.combined) <- "seurat_clusters"
DimPlot(E18_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3)

Idents(object = E18_Ctrlonly) <- "seurat_clusters"
DimPlot(E18_Ctrlonly, reduction = "umap", pt.size = 0.3)
Idents(object = E18_ARKOonly) <- "seurat_clusters"
DimPlot(E18_ARKOonly, reduction = "umap", pt.size = 0.3)

DefaultAssay(E18_Ctrlonly) <- "RNA"
FeaturePlot(E18_Ctrlonly, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_Ctrlonly, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_Ctrlonly, reduction = "umap", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
DefaultAssay(E18_ARKOonly) <- "RNA"
FeaturePlot(E18_ARKOonly, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_ARKOonly, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(E18_ARKOonly, reduction = "umap", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(E18_Ctrlonly, reduction = "umap", features = c("Ar", "EGFP"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
FeaturePlot(E18_ARKOonly, reduction = "umap", features = c("Ar", "EGFP"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)

####P11_CtrlvARKO.combined####
Idents(object = P11_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P11_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P11_CtrlvARKO.combined) <- "stim"
P11_Ctrlonly <- subset(P11_CtrlvARKO.combined, idents = c("P11_Ctrl"))
P11_ARKOonly <- subset(P11_CtrlvARKO.combined, idents = c("P11_ARKO"))

Idents(object = P11_Ctrlonly) <- "CellTypes"
DimPlot(P11_Ctrlonly, reduction = "umap", pt.size = 0.3)

Idents(object = P11_ARKOonly) <- "CellTypes"
DimPlot(P11_ARKOonly, reduction = "umap", pt.size = 0.3)

table(Idents(object = P11_Ctrlonly))

table(Idents(object = P11_ARKOonly))

Idents(object = P11_Ctrlonly) <- "seurat_clusters"
DimPlot(P11_Ctrlonly, reduction = "umap", pt.size = 0.3)

Idents(object = P11_ARKOonly) <- "seurat_clusters"
DimPlot(P11_ARKOonly, reduction = "umap", pt.size = 0.3)

table(Idents(object = P11_Ctrlonly))

table(Idents(object = P11_ARKOonly))

DefaultAssay(P11_Ctrlonly) <- "RNA"
FeaturePlot(P11_Ctrlonly, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_Ctrlonly, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_Ctrlonly, reduction = "umap", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
DefaultAssay(P11_ARKOonly) <- "RNA"
FeaturePlot(P11_ARKOonly, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_ARKOonly, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P11_ARKOonly, reduction = "umap", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

####P42_CtrlvARKO.combined####
Idents(object = P42_CtrlvARKO.combined) <- "CellTypes"
DimPlot(P42_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = P42_CtrlvARKO.combined) <- "stim"
P42_Ctrlonly <- subset(P42_CtrlvARKO.combined, idents = c("P42_Ctrl"))
P42_ARKOonly <- subset(P42_CtrlvARKO.combined, idents = c("P42_ARKO"))

Idents(object = P42_Ctrlonly) <- "CellTypes"
DimPlot(P42_Ctrlonly, reduction = "umap", pt.size = 0.3)

Idents(object = P42_ARKOonly) <- "CellTypes"
DimPlot(P42_ARKOonly, reduction = "umap", pt.size = 0.3)

table(Idents(object = P42_Ctrlonly))

table(Idents(object = P42_ARKOonly))

Idents(object = P42_Ctrlonly) <- "seurat_clusters"
DimPlot(P42_Ctrlonly, reduction = "umap", pt.size = 0.3)

Idents(object = P42_ARKOonly) <- "seurat_clusters"
DimPlot(P42_ARKOonly, reduction = "umap", pt.size = 0.3)

table(Idents(object = P42_Ctrlonly))

table(Idents(object = P42_ARKOonly))

DefaultAssay(P42_Ctrlonly) <- "RNA"
FeaturePlot(P42_Ctrlonly, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_Ctrlonly, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_Ctrlonly, reduction = "umap", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
DefaultAssay(P42_ARKOonly) <- "RNA"
FeaturePlot(P42_ARKOonly, reduction = "umap", features = c("Ar"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_ARKOonly, reduction = "umap", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
FeaturePlot(P42_ARKOonly, reduction = "umap", features = c("Gli1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

