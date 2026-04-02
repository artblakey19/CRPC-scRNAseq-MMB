library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

#### P11 ####

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/P11/Combined_CtrlvARKO/All")


Idents(object = P11_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined) <- "RNA"

tiff(file = "P11_CtrlvARKO.combined UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

tiff(file = "P11_CtrlvARKO.combined Fbln1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined Myh11 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined, reduction = "umap", features = c("Myh11"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined Mki67 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined, reduction = "umap", features = c("Mki67"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined Gli1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined Ar Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined EGFP Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### P11 FBSM ####
#Subcluster FBSM
Idents(object = P11_CtrlvARKO.combined) <- "CellTypes"
P11_CtrlvARKO.combined_FBSM <- subset(P11_CtrlvARKO.combined, idents = c("Fibroblast", "SM"))
Idents(object = P11_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DimPlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)

#Reclustering FibSM
Idents(object = P11_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_FBSM) <- "integrated"
P11_CtrlvARKO.combined_FBSM <- ScaleData(P11_CtrlvARKO.combined_FBSM, verbose = FALSE)P11_CtrlvARKO.combined_FBSM <- RunPCA(E18_CtrlvARKO.combined_FBSM, npcs = 50, verbose = FALSE)
ElbowPlot(P11_CtrlvARKO.combined_FBSM, ndims = 50)

P11_CtrlvARKO.combined_FBSM <- FindNeighbors(P11_CtrlvARKO.combined_FBSM, reduction = "pca", dims = 1:22)
P11_CtrlvARKO.combined_FBSM <- FindClusters(P11_CtrlvARKO.combined_FBSM, resolution = 0.5)
P11_CtrlvARKO.combined_FBSM <- RunUMAP(P11_CtrlvARKO.combined_FBSM, reduction = "pca", dims = 1:22)
P11_CtrlvARKO.combined_FBSM <- RunTSNE(P11_CtrlvARKO.combined_FBSM, reduction = "pca", dims = 1:22)
DimPlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "P11_CtrlvARKO.combined_FBSM UMAP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)
dev.off()

DefaultAssay(P11_CtrlvARKO.combined_FBSM) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_FBSM stim EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Ptch1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ptch1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Cell count
Idents(object = P11_CtrlvARKO.combined_FBSM) <- "stim"
P11_CtrlvARKO.combined_FBSM$stim.seurat_clusters <- paste(Idents(P11_CtrlvARKO.combined_FBSM), P11_CtrlvARKO.combined_FBSM$seurat_clusters, sep = "_")
Idents(object = P11_CtrlvARKO.combined_FBSM) <- "stim.seurat_clusters"
table(Idents(P11_CtrlvARKO.combined_FBSM))

#Heatmap
Idents(object = P11_CtrlvARKO.combined_FBSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_FBSM) <- "RNA"
P11_CtrlvARKO.combined_FBSM <- ScaleData(P11_CtrlvARKO.combined_FBSM, features = rownames(P11_CtrlvARKO.combined_FBSM))
P11_CtrlvARKO.combined_FBSM.markers <- FindAllMarkers(P11_CtrlvARKO.combined_FBSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P11_CtrlvARKO.combined_FBSMTop10 <- P11_CtrlvARKO.combined_FBSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P11_CtrlvARKO.combined_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P11_CtrlvARKO.combined_FBSM, features = c(P11_CtrlvARKO.combined_FBSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P11_CtrlvARKO.combined_FBSM.0.1markers <- FindAllMarkers(P11_CtrlvARKO.combined_FBSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P11_CtrlvARKO.combined_FBSM.0.1markers, "P11_CtrlvARKO.combined_FBSM.0.1markers.csv")

#Cluster
DefaultAssay(P11_CtrlvARKO.combined_FBSM) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Aldh1a3 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Aldh1a3"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
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

#Cluster
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Pgr Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Pgr"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Esr1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Esr1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
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


tiff(file = "P11_CtrlvARKO.combined_FBSM stim Ly6a Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Ly6a"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Cd34 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Cd34"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_FBSM stim Thy1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_FBSM, reduction = "umap", features = c("Thy1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
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

###subclustering EGFP+###
Idents(object = P11_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined) <- "RNA"

P11_CtrlvARKO.combined_EGFP <- subset(x=P11_CtrlvARKO.combined, subset = EGFP > 0)

#Add Ar info
DefaultAssay(P11_CtrlvARKO.combined_EGFP) <- "RNA"
P11_CtrlvARKO.combined_EGFP_ArPos <- subset(x=P11_CtrlvARKO.combined_EGFP, subset = Ar > 0)
P11_CtrlvARKO.combined_EGFP_ArNeg <- subset(x=P11_CtrlvARKO.combined_EGFP, subset = Ar == 0)
Idents(object = P11_CtrlvARKO.combined_EGFP_ArPos) <- "ArPos"
Idents(object = P11_CtrlvARKO.combined_EGFP_ArNeg) <- "ArNeg"
P11_CtrlvARKO.combined_EGFP_ArPos[["ArExp"]] <- Idents(object = P11_CtrlvARKO.combined_EGFP_ArPos)
P11_CtrlvARKO.combined_EGFP_ArNeg[["ArExp"]] <- Idents(object = P11_CtrlvARKO.combined_EGFP_ArNeg)
P11_CtrlvARKO.combined_EGFP_Ar <- merge(x = P11_CtrlvARKO.combined_EGFP_ArPos, y = P11_CtrlvARKO.combined_EGFP_ArNeg)
Idents(object = P11_CtrlvARKO.combined_EGFP_Ar) <- "ArExp"
P11_CtrlvARKO.combined_EGFP$ArExp <- Idents(object = P11_CtrlvARKO.combined_EGFP_Ar)

Idents(object = P11_CtrlvARKO.combined_EGFP) <- "ArExp"
P11_CtrlvARKO.combined_EGFP <- subset(P11_CtrlvARKO.combined_EGFP, idents = c("ArPos", "ArNeg"))
P11_CtrlvARKO.combined_EGFP <- RenameIdents(object = P11_CtrlvARKO.combined_EGFP, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(P11_CtrlvARKO.combined_EGFP, reduction = "umap", pt.size = 0.3)

Idents(object = P11_CtrlvARKO.combined_EGFP) <- "stim"
P11_CtrlvARKO.combined_EGFP$stim.ArExp <- paste(Idents(P11_CtrlvARKO.combined_EGFP), P11_CtrlvARKO.combined_EGFP$ArExp, sep = "_")
Idents(object = P11_CtrlvARKO.combined_EGFP) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = P11_CtrlvARKO.combined_EGFP) <- "CellTypes"
P11_CtrlvARKO.combined_EGFP_FibSM <- subset(P11_CtrlvARKO.combined_EGFP, idents = c("Fibroblast", "SM"))
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"
P11_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P11_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(P11_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(P11_CtrlvARKO.combined_EGFP_FibSM, ndims = 50)

P11_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:15)
P11_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(P11_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:15)
P11_CtrlvARKO.combined_EGFP_FibSM <- RunTSNE(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:15)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM stim split UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Arstim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
P11_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(P11_CtrlvARKO.combined_EGFP_FibSM), P11_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(P11_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Ar stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Esr1 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Esr1"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "P11_CtrlvARKO.combined_EGFP_FibSM Aldh1a3 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Aldh1a3"), cols = c("light grey", "blue"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = P11_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(P11_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
P11_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(P11_CtrlvARKO.combined_EGFP_FibSM, features = rownames(P11_CtrlvARKO.combined_EGFP_FibSM))
P11_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(P11_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
P11_CtrlvARKO.combined_EGFP_FibSMTop10 <- P11_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "P11_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(P11_CtrlvARKO.combined_EGFP_FibSM, features = c(P11_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(P11_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "P11_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")


#### P11 Epi ####
Idents(object = P11_CtrlvARKO.combined) <- "stim"
P11_CtrlvARKO.combined_Ctrl <- subset(P11_CtrlvARKO.combined, idents = c("P11_Ctrl"))
P11_CtrlvARKO.combined_ARKO <- subset(P11_CtrlvARKO.combined, idents = c("P11_ARKO"))

DefaultAssay(P11_CtrlvARKO.combined_Ctrl) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_Ctrl EGFP Epcam Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("EGFP", "Epcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

tiff(file = "P11_CtrlvARKO.combined_Ctrl EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("EGFP"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P11_CtrlvARKO.combined_Ctrl Epcam.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("Epcam"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P11_CtrlvARKO.combined_Ctrl Gli1.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("Gli1"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()

DefaultAssay(P11_CtrlvARKO.combined_ARKO) <- "RNA"
tiff(file = "P11_CtrlvARKO.combined_ARKO EGFP Epcam Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("EGFP", "Epcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

tiff(file = "P11_CtrlvARKO.combined_ARKO EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("EGFP"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P11_CtrlvARKO.combined_ARKO Epcam.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("Epcam"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P11_CtrlvARKO.combined_ARKO Gli1.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P11_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("Gli1"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()

#### P42 Epi ####
Idents(object = P42_CtrlvARKO.combined) <- "stim"
P42_CtrlvARKO.combined_Ctrl <- subset(P42_CtrlvARKO.combined, idents = c("P42_Ctrl"))
P42_CtrlvARKO.combined_ARKO <- subset(P42_CtrlvARKO.combined, idents = c("P42_ARKO"))

tiff(file = "P42_CtrlvARKO.combined_Ctrl EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("EGFP"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P42_CtrlvARKO.combined_Ctrl Epcam.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("Epcam"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P42_CtrlvARKO.combined_Ctrl Gli1.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_Ctrl, reduction = "umap", features = c("Gli1"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()

DefaultAssay(P42_CtrlvARKO.combined_ARKO) <- "RNA"
tiff(file = "P42_CtrlvARKO.combined_ARKO EGFP Epcam Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("EGFP", "Epcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

tiff(file = "P42_CtrlvARKO.combined_ARKO EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("EGFP"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P42_CtrlvARKO.combined_ARKO Epcam.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("Epcam"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()
tiff(file = "P42_CtrlvARKO.combined_ARKO Gli1.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(P42_CtrlvARKO.combined_ARKO, reduction = "umap", features = c("Gli1"), max.cutoff = "q90", cols = c("light grey", "red"), pt.size = 0.7, blend.threshold = 0.1)
dev.off()