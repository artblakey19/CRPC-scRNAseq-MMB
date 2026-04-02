#hHGFtg-hMETtg & hHGFtg-hMETtg-Bcat

#Add necessary tools to library
library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)
library(ggpubr)
library(GGally)
library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(scater)
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)
library(irlba)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-Pb/TriplevDouble/TriplevDouble.combined1.epi1")

DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
Triple.combined1.BELEMETPos <- FindVariableFeatures(Triple.combined1.BELEMETPos, selection.method = "vst", nfeatures = 5000)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, verbose = FALSE)
Triple.combined1.BELEMETPos <- RunPCA(Triple.combined1.BELEMETPos, npcs = 50, verbose = FALSE)
ElbowPlot(Triple.combined1.BELEMETPos, ndims = 50)

Triple.combined1.BELEMETPos <- FindNeighbors(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:23)
Triple.combined1.BELEMETPos <- FindClusters(Triple.combined1.BELEMETPos, resolution = 0.5)
Triple.combined1.BELEMETPos <- RunUMAP(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:23)
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap")

Idents(object = Triple.combined1.BELEMETPos) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap", label = TRUE)

#Cell cycle scoring
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
all.genes <- rownames(Triple.combined1.BELEMETPos)
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, features = all.genes)
Triple.combined1.BELEMETPos <- CellCycleScoring(Triple.combined1.BELEMETPos, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = Triple.combined1.BELEMETPos) <- "Phase"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap")

#Cell Cycle regression
DefaultAssay(Triple.combined1.BELEMETPos) <- "RNA"
Triple.combined1.BELEMETPos <- ScaleData(Triple.combined1.BELEMETPos, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Triple.combined1.BELEMETPos))
Triple.combined1.BELEMETPos <- RunPCA(Triple.combined1.BELEMETPos, features = VariableFeatures(Triple.combined1.BELEMETPos))
ElbowPlot(Triple.combined1.BELEMETPos, ndims = 50)

Triple.combined1.BELEMETPos <- FindNeighbors(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:18)
Triple.combined1.BELEMETPos <- FindClusters(Triple.combined1.BELEMETPos, resolution = 0.5)
Triple.combined1.BELEMETPos <- RunUMAP(Triple.combined1.BELEMETPos, reduction = "pca", dims = 1:18)

Idents(object = Triple.combined1.BELEMETPos) <- "EpiCellTypes"
DimPlot(Triple.combined1.BELEMETPos, reduction = "umap", pt.size = 0.3, label = TRUE)

