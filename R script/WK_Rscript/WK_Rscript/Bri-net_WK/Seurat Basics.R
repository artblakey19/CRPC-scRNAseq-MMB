#Set Identities
Idents(object = PINvTumor.combined) 

#TSNE Plot
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 1)

#Expression TSNE
FeaturePlot(PINvTumor.combined, pt.size = 1, features = c("EGFP"), cols = c("light grey", "red"), min.cutoff = 0)

FeaturePlot(PINvTumor.combined, split.by = "ARQExp", pt.size = 1, features = c("ARQ"), cols = c("light grey", "red"), min.cutoff = 0)

#Dot Plot
DotPlot(PINvTumor.combined, features = c("ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

DotPlot(PIN2, features = c("Aqp1", "Cldn5", "Pecam1", "Plvap", "Cdh5", "Spi1", "Tyrobp", "C1qc", "C1qb", "C1qa", "Actg2", "Myh11", "Tagln", "Rgs5", "Acta2", "Lum", "Pdgfra", "Rspo3", "Apod", "Fbln1", "Ccl5", "Cd3e", "Cd28", "Cxcr6", "Cd3d", "Aqp3", "Col17a1", "Krt14", "Krt15", "Krt5", "Stard10", "Alcam", "Cldn3", "Krt19", "Krt18", "Krt8", "ARQ", "EGFP", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

#DEGs
Idents(object = PINvTumor.combined) <- "seurat_clusters.ARQ"
Basal.ARQ.Markers <- FindMarkers(PINvTumor.combined, ident.1 = "2_ARQPos", ident.2 = "2_ARQNeg")
write.table(Basal.ARQ.Markers, file = "BasalARQMarkers.txt")

#Heatmap
BasalARQTop10 <- Basal.ARQ.Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_logFC)
DoHeatmap(PINvTumor.combined, features = c("Ar", BasalARQTop10$gene)) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

#NewLabel
new.cluster.ids <- c("LE", "BE", "BE", "LE", "LE", "LE", "Lym", "Lym", "LE", "Fibro + SM", "LE", "Lym", "Lym", "Leu", "Lym", "Endo", "BE")
names(new.cluster.ids) <- levels(PIN2)
PIN2 <- RenameIdents(PIN2, new.cluster.ids)

# Visualize co-expression of two features simultaneously
FeaturePlot(pbmc, features = c("MS4A1", "CD79A"), blend = TRUE)
FeaturePlot(PIN2, reduction = "tsne", pt.size = 0.5, features = c("ARQ", "Mki67"), blend = TRUE)

#DEG_PIN2
PIN2cluster4v3 <- FindMarkers(PIN2, ident.1 = 4, ident.2 = 3)
> getwd()
[1] "C:/Users/wkyungkim/Documents"
> write.csv(PIN2cluster4v3, "PIN2cluster4v3.csv")

#Clusters the cells
pbmc <- FindNeighbors(pbmc, dims = 1:10)
pbmc <- FindClusters(pbmc, resolution = 0.5)

#Subset
tcells <- subset(cbmc, idents = c("Naive CD4 T", "Memory CD4 T", "CD8 T"))

#Cell count (split)
Idents(object = PINvTumor.combined.Epi) <- "stim"
PINvTumor.combined.Epi$stim.seurat_clusters <- paste(Idents(PINvTumor.combined.Epi), PINvTumor.combined.Epi$seurat_clusters, sep = "_")
Idents(object = PINvTumor.combined.Epi) <- "stim.seurat_clusters"
table(Idents(PINvTumor.combined.Epi))