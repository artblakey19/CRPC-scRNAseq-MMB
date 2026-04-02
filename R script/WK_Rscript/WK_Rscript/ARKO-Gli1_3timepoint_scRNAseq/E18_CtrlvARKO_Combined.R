library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

#### E18.5 ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/Combined_CtrlvARKO/All")


Idents(object = E18_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined) <- "RNA"

tiff(file = "E18_CtrlvARKO.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

tiff(file = "E18_CtrlvARKO.combined Fbln1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Myh11 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Mki67 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Gli1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Ar Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined EGFP Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### E18.5 FBSM ####
#Subcluster FBSM
Idents(object = E18_CtrlvARKO.combined) <- "CellTypes"
E18_CtrlvARKO.combined_FBSM <- subset(E18_CtrlvARKO.combined, idents = c("Fibroblast", "SM", "Prolif Stro"))
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DimPlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)

#Reclustering FibSM
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "integrated"
E18_CtrlvARKO.combined_FBSM <- ScaleData(E18_CtrlvARKO.combined_FBSM, verbose = FALSE)
E18_CtrlvARKO.combined_FBSM <- RunPCA(E18_CtrlvARKO.combined_FBSM, npcs = 50, verbose = FALSE)
ElbowPlot(E18_CtrlvARKO.combined_FBSM, ndims = 50)

E18_CtrlvARKO.combined_FBSM <- FindNeighbors(E18_CtrlvARKO.combined_FBSM, reduction = "pca", dims = 1:22)
E18_CtrlvARKO.combined_FBSM <- FindClusters(E18_CtrlvARKO.combined_FBSM, resolution = 0.5)
E18_CtrlvARKO.combined_FBSM <- RunUMAP(E18_CtrlvARKO.combined_FBSM, reduction = "pca", dims = 1:22)
E18_CtrlvARKO.combined_FBSM <- RunTSNE(E18_CtrlvARKO.combined_FBSM, reduction = "pca", dims = 1:22)
DimPlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E18_CtrlvARKO.combined_FBSM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)
dev.off()

DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_FBSM stim EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Ptch1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ptch1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()


#Cell count
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim"
E18_CtrlvARKO.combined_FBSM$stim.seurat_clusters <- paste(Idents(E18_CtrlvARKO.combined_FBSM), E18_CtrlvARKO.combined_FBSM$seurat_clusters, sep = "_")
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim.seurat_clusters"
table(Idents(E18_CtrlvARKO.combined_FBSM))

#Heatmap
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
E18_CtrlvARKO.combined_FBSM <- ScaleData(E18_CtrlvARKO.combined_FBSM, features = rownames(E18_CtrlvARKO.combined_FBSM))
E18_CtrlvARKO.combined_FBSM.markers <- FindAllMarkers(E18_CtrlvARKO.combined_FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E18_CtrlvARKO.combined_FBSMTop10 <- E18_CtrlvARKO.combined_FBSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "E18_CtrlvARKO.combined_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(E18_CtrlvARKO.combined_FBSM, features = c(E18_CtrlvARKO.combined_FBSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
E18_CtrlvARKO.combined_FBSM.0.1markers <- FindAllMarkers(E18_CtrlvARKO.combined_FBSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E18_CtrlvARKO.combined_FBSM.0.1markers, "E18_CtrlvARKO.combined_FBSM.0.1markers.csv")

#Cluster6
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Aldh1a3 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Aldh1a3"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Aldh1a1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Aldh1a1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Aldh1a2 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Aldh1a2"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Pdgfrb Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Pdgfrb"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Epha3 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Epha3"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Foxl1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Fst Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Fst"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Fgf2 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Fgf2"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cluster10
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Pgr Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Pgr"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Esr1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Esr1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Gata5 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Gata5"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Inhbb Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Inhbb"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Bmp4 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Bmp4"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Wnt5a Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Wnt5a"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()


tiff(file = "E18_CtrlvARKO.combined_FBSM stim Ly6a Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ly6a"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Cd34 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim Thy1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Thy1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

###Add Ar info
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
E18_CtrlvARKO.combined_FBSM_ArPos <- subset(x=E18_CtrlvARKO.combined_FBSM, subset = Ar > 0)
E18_CtrlvARKO.combined_FBSM_ArNeg <- subset(x=E18_CtrlvARKO.combined_FBSM, subset = Ar == 0)
Idents(object = E18_CtrlvARKO.combined_FBSM_ArPos) <- "ArPos"
Idents(object = E18_CtrlvARKO.combined_FBSM_ArNeg) <- "ArNeg"
E18_CtrlvARKO.combined_FBSM_ArPos[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_FBSM_ArPos)
E18_CtrlvARKO.combined_FBSM_ArNeg[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_FBSM_ArNeg)
E18_CtrlvARKO.combined_FBSM_Ar <- merge(x = E18_CtrlvARKO.combined_FBSM_ArPos, y = E18_CtrlvARKO.combined_FBSM_ArNeg)
Idents(object = E18_CtrlvARKO.combined_FBSM_Ar) <- "ArExp"
E18_CtrlvARKO.combined_FBSM$ArExp <- Idents(object = E18_CtrlvARKO.combined_FBSM_Ar)

Idents(object = E18_CtrlvARKO.combined_FBSM) <- "ArExp"
E18_CtrlvARKO.combined_FBSM <- subset(E18_CtrlvARKO.combined_FBSM, idents = c("ArPos", "ArNeg"))
E18_CtrlvARKO.combined_FBSM <- RenameIdents(object = E18_CtrlvARKO.combined_FBSM, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3)

Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim"
E18_CtrlvARKO.combined_FBSM$stim.ArExp <- paste(Idents(E18_CtrlvARKO.combined_FBSM), E18_CtrlvARKO.combined_FBSM$ArExp, sep = "_")
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim.ArExp"

tiff(file = "E18_CtrlvARKO.combined_FBSM stim.ArExp UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", split.by = "stim.ArExp", pt.size = 0.3, label = TRUE)
dev.off()

DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_FBSM stim.ArExp EGFP Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim.ArExp Ar Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim.ArExp Gli1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stim.ArExp Ptch1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ptch1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell count
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim.ArExp"
E18_CtrlvARKO.combined_FBSM$stim.ArExp.seurat_clusters <- paste(Idents(E18_CtrlvARKO.combined_FBSM), E18_CtrlvARKO.combined_FBSM$seurat_clusters, sep = "_")
Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim.ArExp.seurat_clusters"
table(Idents(E18_CtrlvARKO.combined_FBSM))

#Cluster6
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Aldh1a3 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Aldh1a3"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Pdgfrb Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Pdgfrb"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Fst Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Fst"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Fgf2 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Fgf2"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Foxl1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Foxl1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Aldh1a1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Aldh1a1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()


#Cluster10
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Esr1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Esr1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Inhbb Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Inhbb"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Pgr Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Pgr"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Gata5 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Gata5"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Bmp4 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Bmp4"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_FBSM stimAr Wnt5a Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Wnt5a"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Pseudotime of combined_E18_EGFP_FibSM w ProS by E18_ARKO_ArNeg####

Idents(object = E18_CtrlvARKO.combined_FBSM) <- "stim.ArExp"
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
E18EGFPFibSMPseudo <- as.CellDataSet(E18_CtrlvARKO.combined_FBSM)
E18EGFPFibSMPseudo  <- detectGenes(E18EGFPFibSMPseudo, min_expr = 0.1)
print(head(fData(E18EGFPFibSMPseudo)))

expressed_genes <- row.names(subset(fData(E18EGFPFibSMPseudo),
                                    num_cells_expressed >= 10))

pData(E18EGFPFibSMPseudo)$Total_mRNAs <- Matrix::colSums(exprs(E18EGFPFibSMPseudo))
E18EGFPFibSMPseudo <- E18EGFPFibSMPseudo[,pData(E18EGFPFibSMPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(E18EGFPFibSMPseudo), geom =
        "density")

E18EGFPFibSMPseudo <- estimateSizeFactors(E18EGFPFibSMPseudo)
E18EGFPFibSMPseudo <- estimateDispersions(E18EGFPFibSMPseudo)

disp_table <- dispersionTable(E18EGFPFibSMPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, unsup_clustering_genes$gene_id)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2, num_dim = 22,
                                      reduction_method = 'tSNE', verbose = T)
E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, num_clusters = 2)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "stim.ArExp")


diff_test_res <- differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~stim.ArExp")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, ordering_genes)
plot_ordering_genes(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2,
                                      method = 'DDRTree')

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo)

GM_state <- function(E18EGFPFibSMPseudo){
  if (length(unique(pData(E18EGFPFibSMPseudo)$State)) > 1){
    T0_counts <- table(pData(E18EGFPFibSMPseudo)$State, pData(E18EGFPFibSMPseudo)$stim.ArExp)[,"E18_ARKO_ArNeg"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo, root_state = GM_state(E18EGFPFibSMPseudo))

plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Esr1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))


tiff(file = "E18EGFPFibSMPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo seurat_clusters Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim.ArExp Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo CellTypes Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()


tiff(file = "E18EGFPFibSMPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim.ArExp split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)
dev.off()

tiff(file = "E18EGFPFibSMPseudo Ar.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()

#### Pseudotime of combined_E18_EGFP_FibSM by cluster 0####

Idents(object = E18_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_FBSM) <- "RNA"
E18EGFPFibSMPseudo <- as.CellDataSet(E18_CtrlvARKO.combined_FBSM)
E18EGFPFibSMPseudo  <- detectGenes(E18EGFPFibSMPseudo, min_expr = 0.1)
print(head(fData(E18EGFPFibSMPseudo)))

expressed_genes <- row.names(subset(fData(E18EGFPFibSMPseudo),
                                    num_cells_expressed >= 10))

pData(E18EGFPFibSMPseudo)$Total_mRNAs <- Matrix::colSums(exprs(E18EGFPFibSMPseudo))
E18EGFPFibSMPseudo <- E18EGFPFibSMPseudo[,pData(E18EGFPFibSMPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(E18EGFPFibSMPseudo), geom =
        "density")

E18EGFPFibSMPseudo <- estimateSizeFactors(E18EGFPFibSMPseudo)
E18EGFPFibSMPseudo <- estimateDispersions(E18EGFPFibSMPseudo)

disp_table <- dispersionTable(E18EGFPFibSMPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, unsup_clustering_genes$gene_id)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2, num_dim = 22,
                                      reduction_method = 'tSNE', verbose = T)
E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, num_clusters = 2)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "seurat_clusters")


diff_test_res <- differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, ordering_genes)
plot_ordering_genes(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2,
                                      method = 'DDRTree', ncenter = 300)

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo)

GM_state <- function(E18EGFPFibSMPseudo){
  if (length(unique(pData(E18EGFPFibSMPseudo)$State)) > 1){
    T0_counts <- table(pData(E18EGFPFibSMPseudo)$State, pData(E18EGFPFibSMPseudo)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo, root_state = GM_state(E18EGFPFibSMPseudo))

plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)


setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/Combined_CtrlvARKO/All/Pseudo_bycluster0")

tiff(file = "E18EGFPFibSMPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo CellTypes Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()


tiff(file = "E18EGFPFibSMPseudo seurat_clusters Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "E18EGFPFibSMPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()

tiff(file = "E18EGFPFibSMPseudo stim Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim split Pseudotime.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()

tiff(file = "E18EGFPFibSMPseudo stim.ArExp Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim.ArExp split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)
dev.off()


tiff(file = "E18EGFPFibSMPseudo Ar.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo EGFP.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Gli1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Gli1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ptch1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ptch1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Aldh1a3.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Aldh1a3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ptgfrb.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ptgfrb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Fgf2.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Fgf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Esr1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Esr1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Pgr.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Pgr", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Gata5.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Gata5", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Inhbb.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Inhbb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()


#Alternative ordering cells

E18EGFPFibSMPseudo <- detectGenes(E18EGFPFibSMPseudo, min_expr = 0.1)
fData(E18EGFPFibSMPseudo)$use_for_ordering <-
  fData(E18EGFPFibSMPseudo)$num_cells_expressed > 0.05 * ncol(E18EGFPFibSMPseudo)

plot_pc_variance_explained(E18EGFPFibSMPseudo, return_all = F)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo,
                                       max_components = 2,
                                       norm_method = 'log',
                                       num_dim = 22,
                                       reduction_method = 'tSNE',
                                       verbose = T)

E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, verbose = F)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = 'as.factor(seurat_clusters)')


plot_rho_delta(E18EGFPFibSMPseudo, rho_threshold = 2, delta_threshold = 4 )

E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo,
                                    rho_threshold = 2,
                                    delta_threshold = 4,
                                    skip_rho_sigma = T,
                                    verbose = F)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = 'as.factor(seurat_clusters)')


clustering_DEG_genes <-
  differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                       fullModelFormulaStr = '~seurat_clusters',
                       cores = 1)

ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

E18EGFPFibSMPseudo <-
  setOrderingFilter(E18EGFPFibSMPseudo,
                    ordering_genes = ordering_genes)

E18EGFPFibSMPseudo <-
  reduceDimension(E18EGFPFibSMPseudo, method = 'DDRTree')

E18EGFPFibSMPseudo <-
  orderCells(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <-
  orderCells(E18EGFPFibSMPseudo, root_state = GM_state(E18EGFPFibSMPseudo))

plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/Combined_CtrlvARKO/All/Alternative")

tiff(file = "E18EGFPFibSMPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo CellTypes Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "CellTypes", show_branch_points = FALSE)
dev.off()


tiff(file = "E18EGFPFibSMPseudo seurat_clusters Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) 
dev.off()
tiff(file = "E18EGFPFibSMPseudo seurat_clusters split Pseudotime.tiff", width = 12, height = 6, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + facet_wrap(~seurat_clusters, nrow = 2)
dev.off()

tiff(file = "E18EGFPFibSMPseudo stim Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim split Pseudotime.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)
dev.off()

tiff(file = "E18EGFPFibSMPseudo stim.ArExp Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE)
dev.off()
tiff(file = "E18EGFPFibSMPseudo stim.ArExp split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)
dev.off()


tiff(file = "E18EGFPFibSMPseudo Ar.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo EGFP.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "EGFP", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Gli1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Gli1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ptch1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ptch1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Aldh1a3.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Aldh1a3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ptgfrb.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ptgfrb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Fgf2.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Fgf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Esr1.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Esr1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Pgr.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Pgr", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Gata5.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Gata5", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()
tiff(file = "E18EGFPFibSMPseudo Inhbb.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Inhbb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
dev.off()





#### E18.5 Epi ####

tiff(file = "E18_CtrlvARKO.combined Krt15 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Krt15"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Krt14 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Krt14"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Upk3a Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Upk3a"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined Krt18 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined, reduction = "umap", features = c("Krt18"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Subcluster Epi
Idents(object = E18_CtrlvARKO.combined) <- "CellTypes"
E18_CtrlvARKO.combined_Epi <- subset(E18_CtrlvARKO.combined, idents = c("UGE"))
Idents(object = E18_CtrlvARKO.combined_Epi) <- "seurat_clusters"
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

DefaultAssay(E18_CtrlvARKO.combined_Epi) <- "RNA"
E18_CtrlvARKO.combined_Epi <- subset(x=E18_CtrlvARKO.combined_Epi, subset = Epcam > 0)

#Reclustering Epi
Idents(object = E18_CtrlvARKO.combined_Epi) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_Epi) <- "integrated"
E18_CtrlvARKO.combined_Epi <- ScaleData(E18_CtrlvARKO.combined_Epi, verbose = FALSE)
E18_CtrlvARKO.combined_Epi <- RunPCA(E18_CtrlvARKO.combined_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(E18_CtrlvARKO.combined_Epi, ndims = 50)

E18_CtrlvARKO.combined_Epi <- FindNeighbors(E18_CtrlvARKO.combined_Epi, reduction = "pca", dims = 1:15)
E18_CtrlvARKO.combined_Epi <- FindClusters(E18_CtrlvARKO.combined_Epi, resolution = 0.5)
E18_CtrlvARKO.combined_Epi <- RunUMAP(E18_CtrlvARKO.combined_Epi, reduction = "pca", dims = 1:15)
E18_CtrlvARKO.combined_Epi <- RunTSNE(E18_CtrlvARKO.combined_Epi, reduction = "pca", dims = 1:15)
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Epcam"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

tiff(file = "E18_CtrlvARKO.combined_Epi UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)
dev.off()

E18_CtrlvARKO.combined_Epi <- subset(E18_CtrlvARKO.combined_Epi, idents = c("0", "1", "2", "3", "4", "5"))

#Reclustering Epi
Idents(object = E18_CtrlvARKO.combined_Epi) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_Epi) <- "integrated"
E18_CtrlvARKO.combined_Epi <- ScaleData(E18_CtrlvARKO.combined_Epi, verbose = FALSE)
E18_CtrlvARKO.combined_Epi <- RunPCA(E18_CtrlvARKO.combined_Epi, npcs = 50, verbose = FALSE)
ElbowPlot(E18_CtrlvARKO.combined_Epi, ndims = 50)

E18_CtrlvARKO.combined_Epi <- FindNeighbors(E18_CtrlvARKO.combined_Epi, reduction = "pca", dims = 1:15)
E18_CtrlvARKO.combined_Epi <- FindClusters(E18_CtrlvARKO.combined_Epi, resolution = 0.5)
E18_CtrlvARKO.combined_Epi <- RunUMAP(E18_CtrlvARKO.combined_Epi, reduction = "pca", dims = 1:15)
E18_CtrlvARKO.combined_Epi <- RunTSNE(E18_CtrlvARKO.combined_Epi, reduction = "pca", dims = 1:15)
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E18_CtrlvARKO.combined_Epi UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)
dev.off()

#Cell count
Idents(object = E18_CtrlvARKO.combined_Epi) <- "stim"
E18_CtrlvARKO.combined_Epi$stim.seurat_clusters <- paste(Idents(E18_CtrlvARKO.combined_Epi), E18_CtrlvARKO.combined_Epi$seurat_clusters, sep = "_")
Idents(object = E18_CtrlvARKO.combined_Epi) <- "stim.seurat_clusters"
table(Idents(E18_CtrlvARKO.combined_Epi))

#Heatmap
Idents(object = E18_CtrlvARKO.combined_Epi) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_Epi) <- "RNA"
E18_CtrlvARKO.combined_Epi <- ScaleData(E18_CtrlvARKO.combined_Epi, features = rownames(E18_CtrlvARKO.combined_Epi))
E18_CtrlvARKO.combined_Epi.markers <- FindAllMarkers(E18_CtrlvARKO.combined_Epi, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E18_CtrlvARKO.combined_EpiTop10 <- E18_CtrlvARKO.combined_Epi.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "E18_CtrlvARKO.combined_Epi Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(E18_CtrlvARKO.combined_Epi, features = c(E18_CtrlvARKO.combined_EpiTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
E18_CtrlvARKO.combined_Epi.0.1markers <- FindAllMarkers(E18_CtrlvARKO.combined_Epi, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E18_CtrlvARKO.combined_Epi.0.1markers, "E18_CtrlvARKO.combined_Epi.0.1markers.csv")

#FeaturePlot
DefaultAssay(E18_CtrlvARKO.combined_Epi) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_Epi stim EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim Vim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim Zeb1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Zeb1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim Cdh1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim Cdh1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi stim Twist1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Twist1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "E18_CtrlvARKO.combined_Epi Vim Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi Zeb1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Zeb1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Epi Cdh1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Epi, reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()


#### Pseudotime of 6 EGFP_Epi ####
Idents(E18_CtrlvARKO.combined_Epi) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_Epi) <- "RNA"
ARKOvCtrl_EpiPseudo <- as.CellDataSet(E18_CtrlvARKO.combined_Epi)
ARKOvCtrl_EpiPseudo  <- detectGenes(ARKOvCtrl_EpiPseudo, min_expr = 0.1)
print(head(fData(ARKOvCtrl_EpiPseudo)))

expressed_genes <- row.names(subset(fData(ARKOvCtrl_EpiPseudo),
                                    num_cells_expressed >= 10))

pData(ARKOvCtrl_EpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARKOvCtrl_EpiPseudo))
ARKOvCtrl_EpiPseudo <- ARKOvCtrl_EpiPseudo[,pData(ARKOvCtrl_EpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARKOvCtrl_EpiPseudo), geom =
        "density")

ARKOvCtrl_EpiPseudo <- estimateSizeFactors(ARKOvCtrl_EpiPseudo)
ARKOvCtrl_EpiPseudo <- estimateDispersions(ARKOvCtrl_EpiPseudo)

disp_table <- dispersionTable(ARKOvCtrl_EpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARKOvCtrl_EpiPseudo <- setOrderingFilter(ARKOvCtrl_EpiPseudo, unsup_clustering_genes$gene_id)
#plot_ordering_genes(ARKOvCtrl_EpiPseudo)

#ARKOvCtrl_EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
#plot_pc_variance_explained(ARKOvCtrl_EpiPseudo, return_all = F) # norm_method='log'

ARKOvCtrl_EpiPseudo <- reduceDimension(ARKOvCtrl_EpiPseudo, max_components = 2, num_dim = 15,
                                       reduction_method = 'tSNE', verbose = T)
ARKOvCtrl_EpiPseudo <- clusterCells(ARKOvCtrl_EpiPseudo, num_clusters = 2)

plot_cell_clusters(ARKOvCtrl_EpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(ARKOvCtrl_EpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARKOvCtrl_EpiPseudo <- setOrderingFilter(ARKOvCtrl_EpiPseudo, ordering_genes)
plot_ordering_genes(ARKOvCtrl_EpiPseudo)

ARKOvCtrl_EpiPseudo <- reduceDimension(ARKOvCtrl_EpiPseudo, max_components = 2,
                                       method = 'DDRTree')
ARKOvCtrl_EpiPseudo <- orderCells(ARKOvCtrl_EpiPseudo)

GM_state <- function(ARKOvCtrl_EpiPseudo){
  if (length(unique(pData(ARKOvCtrl_EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(ARKOvCtrl_EpiPseudo)$State, pData(ARKOvCtrl_EpiPseudo)$seurat_clusters)[,"2"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARKOvCtrl_EpiPseudo <- orderCells(ARKOvCtrl_EpiPseudo, root_state = GM_state(ARKOvCtrl_EpiPseudo))

plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

plot_cell_trajectory(ARKOvCtrl_EpiPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))


#Alternative ordering cells

ARKOvCtrl_EpiPseudo <- detectGenes(ARKOvCtrl_EpiPseudo, min_expr = 0.1)
fData(ARKOvCtrl_EpiPseudo)$use_for_ordering <-
  fData(ARKOvCtrl_EpiPseudo)$num_cells_expressed > 0.05 * ncol(ARKOvCtrl_EpiPseudo)

plot_pc_variance_explained(ARKOvCtrl_EpiPseudo, return_all = F)

ARKOvCtrl_EpiPseudo <- reduceDimension(ARKOvCtrl_EpiPseudo,
                                       max_components = 2,
                                       norm_method = 'log',
                                       num_dim = 15,
                                       reduction_method = 'tSNE',
                                       verbose = T)

ARKOvCtrl_EpiPseudo <- clusterCells(ARKOvCtrl_EpiPseudo, verbose = F)

plot_cell_clusters(ARKOvCtrl_EpiPseudo, color_by = 'as.factor(seurat_clusters)')


plot_rho_delta(ARKOvCtrl_EpiPseudo, rho_threshold = 2, delta_threshold = 4 )

ARKOvCtrl_EpiPseudo <- clusterCells(ARKOvCtrl_EpiPseudo,
                                    rho_threshold = 2,
                                    delta_threshold = 4,
                                    skip_rho_sigma = T,
                                    verbose = F)

plot_cell_clusters(ARKOvCtrl_EpiPseudo, color_by = 'as.factor(seurat_clusters)')


clustering_DEG_genes <-
  differentialGeneTest(ARKOvCtrl_EpiPseudo[expressed_genes,],
                       fullModelFormulaStr = '~seurat_clusters',
                       cores = 1)

ordering_genes <-
  row.names(clustering_DEG_genes)[order(clustering_DEG_genes$qval)][1:1000]

ARKOvCtrl_EpiPseudo <-
  setOrderingFilter(ARKOvCtrl_EpiPseudo,
                    ordering_genes = ordering_genes)

ARKOvCtrl_EpiPseudo <-
  reduceDimension(ARKOvCtrl_EpiPseudo, method = 'DDRTree')

ARKOvCtrl_EpiPseudo <-
  orderCells(ARKOvCtrl_EpiPseudo)

ARKOvCtrl_EpiPseudo <-
  orderCells(ARKOvCtrl_EpiPseudo, root_state = GM_state(ARKOvCtrl_EpiPseudo))

plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
plot_cell_trajectory(ARKOvCtrl_EpiPseudo, color_by = "stim", show_branch_points = FALSE) + facet_wrap(~stim, nrow = 1)



tiff(file = "ARKOvCtrl_Epi stim EGFP Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi , reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi stim Vim Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi , reduction = "umap", features = c("Vim"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi stim Cdh1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi , reduction = "umap", features = c("Cdh1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "ARKOvCtrl_Epi stim Zeb1 Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARKOvCtrl_Epi , reduction = "umap", features = c("Zeb1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()