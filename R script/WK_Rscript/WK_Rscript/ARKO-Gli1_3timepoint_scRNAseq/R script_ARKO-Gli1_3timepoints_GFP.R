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
setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK")

load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/E18.5_only_Analysis/E18only_Rworkspace.RData")
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P42_only_Analysis/P42only_Rworkspace.RData")
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P11_only_Analysis/P11only_Rworkspace.RData")

#### E18.5 Ctrl####
E18_Ctrl <- ScaleData(E18_Ctrl, verbose = FALSE)
E18_Ctrl <- RunPCA(E18_Ctrl, npcs = 50, verbose = FALSE)
ElbowPlot(E18_Ctrl)

E18_Ctrl <- FindNeighbors(E18_Ctrl, reduction = "pca", dims = 15)
E18_Ctrl <- FindClusters(E18_Ctrl, resolution = 0.5)

E18_Ctrl <- RunUMAP(E18_Ctrl, reduction = "pca", dims = 1:15)
DimPlot(E18_Ctrl, reduction = "umap", pt.size = 0.3)

DefaultAssay(E18_Ctrl) <- "RNA"
FeaturePlot(E18_Ctrl, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

#Add ARQ info
E18_Ctrl_EGFPPos <- subset(x=E18_Ctrl, subset = EGFP > 1)
Idents(object = E18_Ctrl_EGFPPos) <- "EGFPPos"
E18_Ctrl_EGFPPos[["EGFPExp"]] <- Idents(object = E18_Ctrl_EGFPPos)
Idents(object = E18_Ctrl_EGFPPos) <- "EGFPExp"
E18_Ctrl$EGFPExp <- Idents(object = E18_Ctrl_EGFPPos)

Idents(object = E18_Ctrl) <- "EGFPExp"
E18_Ctrl_EGFP <- subset(E18_Ctrl, idents = c("EGFPPos"))
DimPlot(E18_Ctrl_EGFP , reduction = "umap", pt.size = 0.3)

#### E18.5 ARKO####
E18_ARKO <- ScaleData(E18_ARKO, verbose = FALSE)
E18_ARKO <- RunPCA(E18_ARKO, npcs = 50, verbose = FALSE)
ElbowPlot(E18_ARKO)

E18_ARKO <- FindNeighbors(E18_ARKO, reduction = "pca", dims = 1:15)
E18_ARKO <- FindClusters(E18_E18_ARKO, resolution = 0.5)
E18_ARKO <- RunUMAP(E18_E18_ARKO, reduction = "pca", dims = 1:15)

#Add ARQ info
DefaultAssay(E18_ARKO) <- "RNA"
E18_ARKO_EGFPPos <- subset(x=E18_ARKO, subset = EGFP > 1)
Idents(object = E18_ARKO_EGFPPos) <- "EGFPPos"
E18_ARKO_EGFPPos[["EGFPExp"]] <- Idents(object = E18_Ctrl_EGFPPos)
Idents(object = E18_ARKO_EGFPPos) <- "EGFPExp"
E18_ARKO$EGFPExp <- Idents(object = E18_ARKO_EGFPPos)

Idents(object = E18_ARKO) <- "EGFPExp"
E18_ARKO_EGFP <- subset(E18_ARKO, idents = c("EGFPPos"))
DimPlot(E18_ARKO_EGFP, reduction = "umap", pt.size = 0.3)

