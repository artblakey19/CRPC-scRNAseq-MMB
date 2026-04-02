#### E18.5_ARKOvCtrl_Gli1+ #### 

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/Combined_CtrlvARKO/Gli1+")

#Add EGFP info
DefaultAssay(E18_CtrlvARKO.combined) <- "RNA"

E18_CtrlvARKO.combined_Gli1 <- subset(x=E18_CtrlvARKO.combined, subset = Gli1 > 0)
Idents(object = E18_CtrlvARKO.combined_Gli1) <- "Gli1Pos"
E18_CtrlvARKO.combined_Gli1[["Gli1Exp"]] <- Idents(object = E18_CtrlvARKO.combined_Gli1)
Idents(object = E18_CtrlvARKO.combined_Gli1) <- "Gli1Exp"
E18_CtrlvARKO.combined$Gli1Exp <- Idents(object = E18_CtrlvARKO.combined_Gli1)

Idents(object = E18_CtrlvARKO.combined) <- "Gli1Exp"
E18_CtrlvARKO.combined_Gli1<- subset(E18_CtrlvARKO.combined, idents = c("Gli1Pos"))
DimPlot(E18_CtrlvARKO.combined_Gli1, reduction = "umap", pt.size = 0.3)

#Add Ar info
DefaultAssay(E18_CtrlvARKO.combined_Gli1) <- "RNA"
E18_CtrlvARKO.combined_Gli1_ArPos <- subset(x=E18_CtrlvARKO.combined_Gli1, subset = Ar > 0)
E18_CtrlvARKO.combined_Gli1_ArNeg <- subset(x=E18_CtrlvARKO.combined_Gli1, subset = Ar == 0)
Idents(object = E18_CtrlvARKO.combined_Gli1_ArPos) <- "ArPos"
Idents(object = E18_CtrlvARKO.combined_Gli1_ArNeg) <- "ArNeg"
E18_CtrlvARKO.combined_Gli1_ArPos[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_Gli1_ArPos)
E18_CtrlvARKO.combined_Gli1_ArNeg[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_Gli1_ArNeg)
E18_CtrlvARKO.combined_Gli1_Ar <- merge(x = E18_CtrlvARKO.combined_Gli1_ArPos, y = E18_CtrlvARKO.combined_Gli1_ArNeg)
Idents(object = E18_CtrlvARKO.combined_Gli1_Ar) <- "ArExp"
E18_CtrlvARKO.combined_Gli1$ArExp <- Idents(object = E18_CtrlvARKO.combined_Gli1_Ar)

Idents(object = E18_CtrlvARKO.combined_Gli1) <- "ArExp"
E18_CtrlvARKO.combined_Gli1 <- subset(E18_CtrlvARKO.combined_Gli1, idents = c("ArPos", "ArNeg"))
E18_CtrlvARKO.combined_Gli1 <- RenameIdents(object = E18_CtrlvARKO.combined_Gli1, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(E18_CtrlvARKO.combined_Gli1, reduction = "umap", pt.size = 0.3)

Idents(object = E18_CtrlvARKO.combined_Gli1) <- "stim"
E18_CtrlvARKO.combined_Gli1$stim.ArExp <- paste(Idents(E18_CtrlvARKO.combined_Gli1), E18_CtrlvARKO.combined_Gli1$ArExp, sep = "_")
Idents(object = E18_CtrlvARKO.combined_Gli1) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = E18_CtrlvARKO.combined_Gli1) <- "CellTypes"
E18_CtrlvARKO.combined_Gli1_FibSM <- subset(E18_CtrlvARKO.combined_Gli1, idents = c("Fibroblast", "SM", "Prolif Stro"))
Idents(object = E18_CtrlvARKO.combined_Gli1_FibSM) <- "seurat_clusters"
DimPlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = E18_CtrlvARKO.combined_Gli1_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_Gli1_FibSM) <- "integrated"
E18_CtrlvARKO.combined_Gli1_FibSM <- ScaleData(E18_CtrlvARKO.combined_Gli1_FibSM, verbose = FALSE)
E18_CtrlvARKO.combined_Gli1_FibSM <- RunPCA(E18_CtrlvARKO.combined_Gli1_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(E18_CtrlvARKO.combined_Gli1_FibSM, ndims = 50)

E18_CtrlvARKO.combined_Gli1_FibSM <- FindNeighbors(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "pca", dims = 1:20)
E18_CtrlvARKO.combined_Gli1_FibSM <- FindClusters(E18_CtrlvARKO.combined_Gli1_FibSM, resolution = 0.5)
E18_CtrlvARKO.combined_Gli1_FibSM <- RunUMAP(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "pca", dims = 1:20)
E18_CtrlvARKO.combined_Gli1_FibSM <- RunTSNE(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "pca", dims = 1:20)
DimPlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E18_CtrlvARKO.combined_Gli1_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", pt.size = 1, label = TRUE)
dev.off()

Idents(object = E18_CtrlvARKO.combined_EGFP2_FibSM) <- "seurat_clusters"
tiff(file = "E18_CtrlvARKO.combined_Gli1_FibSM Arstim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", pt.size = 1, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = E18_CtrlvARKO.combined_Gli1_FibSM) <- "stim.ArExp"
E18_CtrlvARKO.combined_Gli1_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(E18_CtrlvARKO.combined_Gli1_FibSM), E18_CtrlvARKO.combined_Gli1_FibSM$seurat_clusters, sep = "_")
Idents(object = E18_CtrlvARKO.combined_Gli1_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(E18_CtrlvARKO.combined_Gli1_FibSM))

Idents(object = E18_CtrlvARKO.combined_Gli1_FibSM) <- "CellTypes"
DimPlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)

#Featureplots
DefaultAssay(E18_CtrlvARKO.combined_Gli1_FibSM) <- "RNA"

tiff(file = "E18_CtrlvARKO.combined_Gli1_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Gli1_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_Gli1_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_Gli1_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = E18_CtrlvARKO.combined_Gli1_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_Gli1_FibSM) <- "RNA"
E18_CtrlvARKO.combined_Gli1_FibSM <- ScaleData(E18_CtrlvARKO.combined_Gli1_FibSM, features = rownames(E18_CtrlvARKO.combined_Gli1_FibSM))
E18_CtrlvARKO.combined_Gli1_FibSM.markers <- FindAllMarkers(E18_CtrlvARKO.combined_Gli1_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E18_CtrlvARKO.combined_Gli1_FibSMTop10 <- E18_CtrlvARKO.combined_Gli1_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(E18_CtrlvARKO.combined_Gli1_FibSM, features = c(E18_CtrlvARKO.combined_Gli1_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

#DEGs
E18_CtrlvARKO.combined_Gli1_FibSM.0.1markers <- FindAllMarkers(E18_CtrlvARKO.combined_Gli1_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E18_CtrlvARKO.combined_Gli1_FibSM.0.1markers, "E18_CtrlvARKO.combined_Gli1_FibSM.0.1markers.csv")

