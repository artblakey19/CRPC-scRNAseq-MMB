#### E18.5_ARKOvCtrl_EGFP+ #### 

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/GFP+/with ProS")


#subclustering EGFP+
Idents(object = E18_CtrlvARKO.combined) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined) <- "RNA"

E18_CtrlvARKO.combined_EGFP <- subset(x=E18_CtrlvARKO.combined, subset = EGFP > 0)

#Add Ar info
DefaultAssay(E18_CtrlvARKO.combined_EGFP) <- "RNA"
E18_CtrlvARKO.combined_EGFP_ArPos <- subset(x=E18_CtrlvARKO.combined_EGFP, subset = Ar > 0)
E18_CtrlvARKO.combined_EGFP_ArNeg <- subset(x=E18_CtrlvARKO.combined_EGFP, subset = Ar == 0)
Idents(object = E18_CtrlvARKO.combined_EGFP_ArPos) <- "ArPos"
Idents(object = E18_CtrlvARKO.combined_EGFP_ArNeg) <- "ArNeg"
E18_CtrlvARKO.combined_EGFP_ArPos[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_EGFP_ArPos)
E18_CtrlvARKO.combined_EGFP_ArNeg[["ArExp"]] <- Idents(object = E18_CtrlvARKO.combined_EGFP_ArNeg)
E18_CtrlvARKO.combined_EGFP_Ar <- merge(x = E18_CtrlvARKO.combined_EGFP_ArPos, y = E18_CtrlvARKO.combined_EGFP_ArNeg)
Idents(object = E18_CtrlvARKO.combined_EGFP_Ar) <- "ArExp"
E18_CtrlvARKO.combined_EGFP$ArExp <- Idents(object = E18_CtrlvARKO.combined_EGFP_Ar)

Idents(object = E18_CtrlvARKO.combined_EGFP) <- "ArExp"
E18_CtrlvARKO.combined_EGFP <- subset(E18_CtrlvARKO.combined_EGFP, idents = c("ArPos", "ArNeg"))
E18_CtrlvARKO.combined_EGFP <- RenameIdents(object = E18_CtrlvARKO.combined_EGFP, 'ArNeg' = "ArNeg", 'ArPos' = "ArPos")
DimPlot(E18_CtrlvARKO.combined_EGFP, reduction = "umap", pt.size = 0.3)

Idents(object = E18_CtrlvARKO.combined_EGFP) <- "stim"
E18_CtrlvARKO.combined_EGFP$stim.ArExp <- paste(Idents(E18_CtrlvARKO.combined_EGFP), E18_CtrlvARKO.combined_EGFP$ArExp, sep = "_")
Idents(object = E18_CtrlvARKO.combined_EGFP) <- "stim.ArExp"

#Subclustering FibSM
Idents(object = E18_CtrlvARKO.combined_EGFP) <- "CellTypes"
E18_CtrlvARKO.combined_EGFP_FibSM <- subset(E18_CtrlvARKO.combined_EGFP, idents = c("Fibroblast", "SM"))
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3)

#Reclustering FibSM
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "integrated"
E18_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(E18_CtrlvARKO.combined_EGFP_FibSM, verbose = FALSE)
E18_CtrlvARKO.combined_EGFP_FibSM <- RunPCA(E18_CtrlvARKO.combined_EGFP_FibSM, npcs = 50, verbose = FALSE)
ElbowPlot(E18_CtrlvARKO.combined_EGFP_FibSM, ndims = 50)

E18_CtrlvARKO.combined_EGFP_FibSM <- FindNeighbors(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
E18_CtrlvARKO.combined_EGFP_FibSM <- FindClusters(E18_CtrlvARKO.combined_EGFP_FibSM, resolution = 0.5)
E18_CtrlvARKO.combined_EGFP_FibSM <- RunUMAP(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
E18_CtrlvARKO.combined_EGFP_FibSM <- RunTSNE(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "pca", dims = 1:20)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)

tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, label = TRUE)
dev.off()

tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM stim split UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", split.by = "stim", pt.size = 0.3, label = TRUE)
dev.off()

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Arstim UMAP.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3, split.by = "stim.ArExp", label = TRUE)
dev.off()

#Cell count (split)
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
E18_CtrlvARKO.combined_EGFP_FibSM$stim.ArExp.seurat_clusters <- paste(Idents(E18_CtrlvARKO.combined_EGFP_FibSM), E18_CtrlvARKO.combined_EGFP_FibSM$seurat_clusters, sep = "_")
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp.seurat_clusters"
table(Idents(E18_CtrlvARKO.combined_EGFP_FibSM))

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "CellTypes"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Celltypes UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", pt.size = 0.3,label = TRUE)
dev.off()

#Featureplots
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM EGFP Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("EGFP"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ar Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Gli1 Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gli1"), cols = c("light grey", "red"), split.by = "stim.ArExp", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ar stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ar"), cols = c("light grey", "red"), split.by = "stim", min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Heatmap
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18_CtrlvARKO.combined_EGFP_FibSM <- ScaleData(E18_CtrlvARKO.combined_EGFP_FibSM, features = rownames(E18_CtrlvARKO.combined_EGFP_FibSM))
E18_CtrlvARKO.combined_EGFP_FibSM.markers <- FindAllMarkers(E18_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
E18_CtrlvARKO.combined_EGFP_FibSMTop10 <- E18_CtrlvARKO.combined_EGFP_FibSM.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)

tiff(file = "E18_FBSM Heatmap Top10 purple.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(E18_CtrlvARKO.combined_EGFP_FibSM, features = c(E18_CtrlvARKO.combined_EGFP_FibSMTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#DEGs
E18_CtrlvARKO.combined_EGFP_FibSM.0.1markers <- FindAllMarkers(E18_CtrlvARKO.combined_EGFP_FibSM, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E18_CtrlvARKO.combined_EGFP_FibSM.0.1markers, "E18_CtrlvARKO.combined_EGFP_FibSM.0.1markers.csv")


#### Cytotrace ####
Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18_CtrlvARKO.combined_EGFP_FibSMExpdata = data.frame(E18_CtrlvARKO.combined_EGFP_FibSM[["RNA"]]@data)
E18_CtrlvARKO.combined_EGFP_FibSMphenotypedata = FetchData(E18_CtrlvARKO.combined_EGFP_FibSM, c("seurat_clusters"))
write.table(E18_CtrlvARKO.combined_EGFP_FibSMExpdata, file = "FibSMExpdata.txt")
write.table(E18_CtrlvARKO.combined_EGFP_FibSMphenotypedata, file = "FibSMphenotypedata.txt")


#Featureplots
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Clu Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Clu"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Thy1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Thy1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Hand2 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Hand2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ly6e Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ly6e"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Rspo3 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Rspo3"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Sfrp1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Sfrp1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ptn Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ptn"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Aldh1a3 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Aldh1a3"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Pdgfra Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Pdgfra"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Dlk1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Dlk1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Egfr Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Egfr"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Postn Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Postn"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Col6a6 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Col6a6"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Igf2 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Igf2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ccl2 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ccl2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ngf Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ngf"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Csf1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Csf1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Fgf18 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Fgf18"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ctgf Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ctgf"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ptch1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ptch1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Fzd2 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Fzd2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Epha4 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Epha4"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Esr1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Esr1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Bmp2 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Bmp2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Bmp4 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Bmp4"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Bmp7 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Bmp7"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Wnt5a Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Wnt5a"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Dkk2 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Dkk2"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Nkd1 Exp.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Nkd1"), cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Split
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Aldh1a3 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Aldh1a3"), split.by = "stim", cols = c("light grey", "blue"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Pdgfrb stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Pdgfrb"), split.by = "stim", cols = c("light grey", "blue"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Pdgfra stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Pdgfra"), split.by = "stim", cols = c("light grey", "blue"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Fgf2 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Fgf2"), split.by = "stim", cols = c("light grey", "blue"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Esr1 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Esr1"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Pgr stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Pgr"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Gata5 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Gata5"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Inhbb stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Inhbb"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()


tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Ly6a stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Ly6a"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Thy1 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Thy1"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "E18_CtrlvARKO.combined_EGFP_FibSM Cd34 stim Exp.tiff", width = 8, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18_CtrlvARKO.combined_EGFP_FibSM, reduction = "umap", features = c("Cd34"), split.by = "stim", cols = c("light grey", "red"), min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

#### Pseudotime of combined_E18_EGFP_FibSM w ProS by E18_ARKO_ArNeg####

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "stim.ArExp"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18EGFPFibSMPseudo <- as.CellDataSet(E18_CtrlvARKO.combined_EGFP_FibSM)
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

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2, num_dim = 20,
                                      reduction_method = 'tSNE', verbose = T)
E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, num_clusters = 2)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "stim.ArExp")


diff_test_res <- differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~stim.ArExp")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, ordering_genes)
plot_ordering_genes(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2,
                                      method = 'DDRTree', ncenter = 50)

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
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/GFP+/with ProS/Pseudo_by_E18_ARKO")

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
tiff(file = "E18EGFPFibSMPseudo stim.ArExp split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)
dev.off()

tiff(file = "E18EGFPFibSMPseudo Ar.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ar", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Gli1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Gli1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Rspo3.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Rspo3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Mki67.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Mki67", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Sfrp1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Sfrp1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ptn.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ptn", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Igf2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Igf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ngf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ngf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Csf1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Csf1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Fgf2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Fgf2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Fgf18.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Fgf18", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ctgf.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ctgf", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Bmp2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Bmp2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Bmp4.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Bmp4", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Bmp7.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Bmp7", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()

tiff(file = "E18EGFPFibSMPseudo Pdgfrb.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Pdgfrb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Aldh1a3.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Aldh1a3", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Clu.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Clu", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Thy1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Thy1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Hand2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Hand2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ly6e.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ly6e", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple"))
dev.off()

tiff(file = "E18EGFPFibSMPseudo Wnt5a.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Wnt5a", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Ptch1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Ptch1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Bmp4.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Bmp4", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Dkk2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Dkk2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Nkd1.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Nkd1", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()
tiff(file = "E18EGFPFibSMPseudo Axin2.tiff", width = 3, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Axin2", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red"))
dev.off()

#### Pseudotime of combined_E18_EGFP_FibSM w ProS by seurat_clusters####

Idents(object = E18_CtrlvARKO.combined_EGFP_FibSM) <- "seurat_clusters"
DefaultAssay(E18_CtrlvARKO.combined_EGFP_FibSM) <- "RNA"
E18EGFPFibSMPseudo <- as.CellDataSet(E18_CtrlvARKO.combined_EGFP_FibSM)
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

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2, num_dim = 20,
                                      reduction_method = 'tSNE', verbose = T)
E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, num_clusters = 2)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "seurat_clusters")


diff_test_res <- differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, ordering_genes)
plot_ordering_genes(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2,
                                      method = 'DDRTree')

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo)

GM_state <- function(E18EGFPFibSMPseudo){
  if (length(unique(pData(E18EGFPFibSMPseudo)$State)) > 1){
    T0_counts <- table(pData(E18EGFPFibSMPseudo)$State, pData(E18EGFPFibSMPseudo)$seurat_clusters)[,"10"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

E18EGFPFibSMPseudo <- orderCells(E18EGFPFibSMPseudo, root_state = GM_state(E18EGFPFibSMPseudo))

plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)

setwd("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E18.5/GFP+/with ProS/Pseudo_by_seurat_clusters")

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
tiff(file = "E18EGFPFibSMPseudo stim.ArExp split Pseudotime.tiff", width = 12, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, color_by = "stim.ArExp", show_branch_points = FALSE) + facet_wrap(~stim.ArExp, nrow = 1)
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

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2, num_dim = 20,
                                      reduction_method = 'tSNE', verbose = T)
E18EGFPFibSMPseudo <- clusterCells(E18EGFPFibSMPseudo, num_clusters = 2)

plot_cell_clusters(E18EGFPFibSMPseudo, color_by = "seurat_clusters")


diff_test_res <- differentialGeneTest(E18EGFPFibSMPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
E18EGFPFibSMPseudo <- setOrderingFilter(E18EGFPFibSMPseudo, ordering_genes)
plot_ordering_genes(E18EGFPFibSMPseudo)

E18EGFPFibSMPseudo <- reduceDimension(E18EGFPFibSMPseudo, max_components = 2,
                                      method = 'DDRTree')

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
tiff(file = "E18EGFPFibSMPseudo Pdgfrb.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(E18EGFPFibSMPseudo, markers = "Pdgfrb", use_color_gradient = TRUE, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "purple")) + facet_wrap(~stim, nrow = 1)
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
