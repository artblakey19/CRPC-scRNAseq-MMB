

## load in libraries necessary for analysis

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

#PIN PreFiltering Vln
tiff(file = "PIN Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(PIN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

tiff(file = "PIN Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Pre-filteration")
dev.off()

tiff(file = "PIN Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Pre-filteration")
dev.off()

#PIN PostFiltering Vln
tiff(file = "PIN Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(PIN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

tiff(file = "PIN Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN2@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Post-filteration")
dev.off()

tiff(file = "PIN Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(PIN2@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Post-filteration")
dev.off()

tiff(file = "PIN Post-filter Jackstraw.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
JackStrawPlot(PIN2, dims = 1:20)
dev.off()

#Tumor PreFiltering Vln
tiff(file = "Tumor Pre-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

tiff(file = "Tumor Pre-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Pre-filteration")
dev.off()

tiff(file = "Tumor Pre-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Pre-filteration")
dev.off()

#Tumor PostFiltering Vln
tiff(file = "Tumor Post-filter Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(Tumor2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)
dev.off()

tiff(file = "Tumor Post-filter nFeature Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Post-filteration")
dev.off()

tiff(file = "Tumor Post-filter percentMT Hist.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
hist(Tumor2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Post-filteration")
dev.off()

tiff(file = "Tumor Post-filter Jackstraw.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
JackStrawPlot(Tumor2, dims = 1:20)
dev.off()


#Dimplot

#Fig1E
tiff(file = "PIN2 TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PIN2, reduction = "tsne", pt.size = 0.3, group.by = "stim", cols = c("#B71800")) 
dev.off()

tiff(file = "Tumor2 TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Tumor2, reduction = "tsne", pt.size = 0.3, group.by = "stim", cols = c("#009FB7")) 
dev.off()

tiff(file = "PINvTumor.combined TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 0.3, group.by = "stim", cols = c("#B71800", "#009FB7")) 
dev.off()

#NewLabel
new.cluster.ids <- c("LE", "LE", "LE", "BE", "LE", "Lym", "LE", "LE", "LE", "LE", "BE", "LE", "Fib", "LE", "Leu", "Endo", "SM", "BE")
names(new.cluster.ids) <- levels(PINvTumor.combined)
PINvTumor.combined <- RenameIdents(PINvTumor.combined, new.cluster.ids)

#Fig1F
tiff(file = "Combined TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 0.3)
dev.off()

#Fig1G
PINonly <- subset(PINvTumor.combined, idents = c("PIN"))
Idents(object = PINonly) <- "seurat_clusters"
new.cluster.ids <- c("LE", "LE", "LE", "BE", "LE", "Lym", "LE", "LE", "LE", "LE", "BE", "LE", "Fib", "LE", "Leu", "Endo", "SM", "BE")
names(new.cluster.ids) <- levels(PINonly)
PINonly <- RenameIdents(PINonly, new.cluster.ids)

tiff(file = "PINonly TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINonly, reduction = "tsne", pt.size = 0.3) 
dev.off()

Idents(object = PINvTumor.combined) <- "stim"
Tumoronly <- subset(PINvTumor.combined, idents = c("Tumor"))
Idents(object = Tumoronly) <- "seurat_clusters"
new.cluster.ids <- c("LE", "LE", "LE", "BE", "LE", "Lym", "LE", "LE", "LE", "LE", "BE", "LE", "Fib", "LE", "Leu", "Endo", "SM", "BE")
names(new.cluster.ids) <- levels(Tumoronly)
Tumoronly <- RenameIdents(Tumoronly, new.cluster.ids)
DimPlot(Tumoronly, reduction = "tsne", pt.size = 0.3)

tiff(file = "Tumoronly TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Tumoronly, reduction = "tsne", pt.size = 0.3) 
dev.off()

#DotPlot

#Fig1H
Idents(object = PINvTumor.combined) <- "seurat_clusters"
PINvTumor.combined <- RenameIdents(object = PINvTumor.combined, '3' = "BE", '10' = "BE", '17' = "BE", '0' = "LE", '1' = "LE", '2' = "LE", '4' = "LE", '6' = "LE", '7' = "LE", '8' = "LE", '9' = "LE", '11' = "LE", '13' = "LE", '12' = "Fib", '16' = "SM", '14' = "Leu", '15' = "Endo", '5' = "Lym")
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 1, label = TRUE)
PINvTumor.combined[["CellType"]] <- Idents(object = PINvTumor.combined)

DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

tiff(file = "Combined DotPlot.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#FeaturePlot

#Fig1I
Idents(object = PINvTumor.combined) <- "RNA"

tiff(file = "Epcam combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Epcam"), cols = c("light grey", "red"),min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Myh11 combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Myh11"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Fbln1 combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Ar combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "ARQ combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "EGFP combined.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Fig1J
Idents(object = PINvTumor.combined) <- "stim"
PINvTumor.combined$stim.seurat_clusters <- paste(Idents(PINvTumor.combined), PINvTumor.combined$seurat_clusters, sep = "_")
Idents(object = PINvTumor.combined) <- "stim.seurat_clusters"
table(Idents(PINvTumor.combined))
