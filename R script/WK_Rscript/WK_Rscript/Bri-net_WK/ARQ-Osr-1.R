---
  title: "PIN and Tumor Independant Analysis"
author: "Omar Khalid, Ph.D."
date: "9/19/2019"
output: html_document
editor_options: 
  chunk_output_type: console
---
  
  ```{r}

## load in libraries necessary for analysis

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)

PIN.data <- Read10X("./PIN4_outs/raw_feature_bc_matrix/")
PIN <- CreateSeuratObject(counts = PIN.data,  min.cells = 3, min.features = 500, project = "PIN")
PIN <- NormalizeData(PIN)

PIN[["percent.mt"]] <- PercentageFeatureSet(PIN, pattern = "^mt-")

#reanalyze here with new filtering paramaters
PIN2 <- subset(PIN, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)

PIN2 <- FindVariableFeatures(PIN2, selection.method = "vst", nfeatures = 5000)



VlnPlot(PIN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)

hist(PIN@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN pre filteration")

hist(PIN@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN pre filteration")


VlnPlot(PIN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3)


hist(PIN2@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN post filteration")

hist(PIN2@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN post filteration")

PIN2 <- ScaleData(PIN2, verbose = FALSE)
PIN2 <- RunPCA(PIN2, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
PIN2 <- RunUMAP(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- FindNeighbors(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- FindClusters(PIN2, resolution = 0.5)

PIN2 <- RunTSNE(PIN2, reduction = "pca", dims = 1:20)
DimPlot(PIN2, reduction = "umap", label = TRUE)
DimPlot(PIN2, reduction = "tsne", label = TRUE)


PIN2 <- JackStraw(PIN2, num.replicate = 100)
PIN2 <- ScoreJackStraw(PIN2, dims = 1:20)
JackStrawPlot(PIN2, dims = 1:20)


PIN2 <- RunTSNE(PIN2, reduction = "pca", dims = 1:20)
DimPlot(PIN2, reduction = "tsne", label = TRUE)

PIN2 <- FindClusters(PIN2, resolution = 0.5)


DotPlot(PIN2, features = c("Krt8","Krt18","Krt19","Msmb","Acpp","Krt5","Krt14","Krt15","Trp63","Dst","Syp","Grp","Calca","Chga","Tph1", "Scg2", "Fbln1", "Vim", "Apod", "Dcn", "Ptgds","Cfd","Acta2","Tpm2", "Myh11", "Rgs5", "Mt1a","Rgs1","C1qa","C1qb","Tyrobp","C1qc","Ifi27","Ackr1","Sele","Hmox1","Cln5","EGFP","ARQ","Ar"), dot.scale = 8) + RotatedAxis()

PIN2 <- FindClusters(PIN2, resolution = 0.3)


DotPlot(PIN2, features = c("Krt8","Krt18","Krt19","Alcam","Cd24a","Krt5","Krt14","Krt15","Trp63","Col17a1","Syp","Stmn1","Tubb5","Chga","Scg2", "Fbln1",  "Apod", "Fgf2", "Rspo3","Ptgs2","Acta2","Myh11", "Rgs5", "Myl9","Tpm2","C1qa","C1qb","Tyrobp","C1qc","Spi1","Cdh5","Cldn5","Pecam1","Cd200","Sele","Cd28","Cd2","Ccl5","Pdcd1","Ctla4","Ar","ARQ","EGFP"), dot.scale = 8, col.min = 0, assay = "RNA") + RotatedAxis()

```

```{r}

PIN2cluster0 <- FindMarkers(PIN2, ident.1 = 0)
PIN2cluster1 <- FindMarkers(PIN2, ident.1 = 1)
PIN2cluster2 <- FindMarkers(PIN2, ident.1 = 2)
PIN2cluster3 <- FindMarkers(PIN2, ident.1 = 3)
PIN2cluster4 <- FindMarkers(PIN2, ident.1 = 4)
PIN2cluster5 <- FindMarkers(PIN2, ident.1 = 5)
PIN2cluster6 <- FindMarkers(PIN2, ident.1 = 6)
PIN2cluster7 <- FindMarkers(PIN2, ident.1 = 7)
PIN2cluster8 <- FindMarkers(PIN2, ident.1 = 8)
PIN2cluster9 <- FindMarkers(PIN2, ident.1 = 9)
PIN2cluster10 <- FindMarkers(PIN2, ident.1 = 10)
PIN2cluster11 <- FindMarkers(PIN2, ident.1 = 11)
PIN2cluster12 <- FindMarkers(PIN2, ident.1 = 12)
PIN2cluster13 <- FindMarkers(PIN2, ident.1 = 13)
PIN2cluster14 <- FindMarkers(PIN2, ident.1 = 14)
PIN2cluster15 <- FindMarkers(PIN2, ident.1 = 15)
PIN2cluster16 <- FindMarkers(PIN2, ident.1 = 16)



write.csv(PIN2cluster0, "PIN2cluster0.csv")
write.csv(PIN2cluster1, "PIN2cluster1.csv")
write.csv(PIN2cluster2, "PIN2cluster2.csv")
write.csv(PIN2cluster3, "PIN2cluster3.csv")
write.csv(PIN2cluster4, "PIN2cluster4.csv")
write.csv(PIN2cluster5, "PIN2cluster5.csv")
write.csv(PIN2cluster6, "PIN2cluster6.csv")
write.csv(PIN2cluster7, "PIN2cluster7.csv")
write.csv(PIN2cluster8, "PIN2cluster8.csv")
write.csv(PIN2cluster9, "PIN2cluster9.csv")
write.csv(PIN2cluster10, "PIN2cluster10.csv")
write.csv(PIN2cluster11, "PIN2cluster11.csv")
write.csv(PIN2cluster12, "PIN2cluster12.csv")
write.csv(PIN2cluster13, "PIN2cluster13.csv")
write.csv(PIN2cluster14, "PIN2cluster14.csv")
write.csv(PIN2cluster15, "PIN2cluster15.csv")
write.csv(PIN2cluster16, "PIN2cluster16.csv")



```

```{r}

Tumor.data <- Read10X("./Tumor4_outs/raw_feature_bc_matrix/")
Tumor <- CreateSeuratObject(counts = Tumor.data,  min.cells = 3, min.features = 500, project = "Tumor")
Tumor <- NormalizeData(Tumor)


Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^mt-")

Tumor2 <- subset(Tumor, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)

Tumor2 <- FindVariableFeatures(Tumor2, selection.method = "vst", nfeatures = 5000)

VlnPlot(Tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3, col="lightblue")

hist(Tumor@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor pre filteration")

hist(Tumor@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor pre filteration")


VlnPlot(Tumor2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size =0.1, ncol = 3, col="lightblue")

hist(Tumor2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor post filteration")

hist(Tumor2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor post filteration")


Tumor2 <- ScaleData(Tumor2, verbose = FALSE)
Tumor2 <- RunPCA(Tumor2, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
Tumor2 <- RunUMAP(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- FindNeighbors(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- FindClusters(Tumor2, resolution = 0.5)

Tumor2 <- RunTSNE(Tumor2, reduction = "pca", dims = 1:20)
DimPlot(Tumor2, reduction = "umap", label = TRUE)
DimPlot(Tumor2, reduction = "tsne", label = TRUE)

Tumor2 <- JackStraw(Tumor2, num.replicate = 100)
Tumor2 <- ScoreJackStraw(Tumor2, dims = 1:20)
JackStrawPlot(Tumor2, dims = 1:20)



DotPlot(Tumor2, features = c("Krt8","Krt18","Krt19","Msmb","Acpp","Krt5","Krt14","Krt15","Trp63","Dst","Syp","Grp","Calca","Chga","Tph1", "Scg2", "Fbln1", "Vim", "Apod", "Dcn", "Ptgds","Cfd","Acta2","Tpm2", "Myh11", "Rgs5", "Mt1a","Rgs1","C1qa","C1qb","Tyrobp","C1qc","Ifi27","Ackr1","Sele","Hmox1","Cln5","EGFP","ARQ","Ar"), dot.scale = 8) + RotatedAxis()


Tumor2cluster0 <- FindMarkers(Tumor2, ident.1 = 0)
Tumor2cluster1 <- FindMarkers(Tumor2, ident.1 = 1)
Tumor2cluster2 <- FindMarkers(Tumor2, ident.1 = 2)
Tumor2cluster3 <- FindMarkers(Tumor2, ident.1 = 3)
Tumor2cluster4 <- FindMarkers(Tumor2, ident.1 = 4)
Tumor2cluster5 <- FindMarkers(Tumor2, ident.1 = 5)
Tumor2cluster6 <- FindMarkers(Tumor2, ident.1 = 6)
Tumor2cluster7 <- FindMarkers(Tumor2, ident.1 = 7)
Tumor2cluster8 <- FindMarkers(Tumor2, ident.1 = 8)
Tumor2cluster9 <- FindMarkers(Tumor2, ident.1 = 9)
Tumor2cluster10 <- FindMarkers(Tumor2, ident.1 = 10)
Tumor2cluster11 <- FindMarkers(Tumor2, ident.1 = 11)
Tumor2cluster12 <- FindMarkers(Tumor2, ident.1 = 12)
Tumor2cluster13 <- FindMarkers(Tumor2, ident.1 = 13)

write.csv(Tumor2cluster0, "Tumor2cluster0.csv")
write.csv(Tumor2cluster1, "Tumor2cluster1.csv")
write.csv(Tumor2cluster2, "Tumor2cluster2.csv")
write.csv(Tumor2cluster3, "Tumor2cluster3.csv")
write.csv(Tumor2cluster4, "Tumor2cluster4.csv")
write.csv(Tumor2cluster5, "Tumor2cluster5.csv")
write.csv(Tumor2cluster6, "Tumor2cluster6.csv")
write.csv(Tumor2cluster7, "Tumor2cluster7.csv")
write.csv(Tumor2cluster8, "Tumor2cluster8.csv")
write.csv(Tumor2cluster9, "Tumor2cluster9.csv")
write.csv(Tumor2cluster10, "Tumor2cluster10.csv")
write.csv(Tumor2cluster11, "Tumor2cluster11.csv")
write.csv(Tumor2cluster12, "Tumor2cluster12.csv")
write.csv(Tumor2cluster13, "Tumor2cluster13.csv")


```

```{r}

#Add ARQ Exp Data
PINARQPos <- subset(x=PIN2, subset = ARQ > 0, slot = 'counts')
PINARQNeg <- subset(x=PIN2, subset = ARQ == 0, slot = 'counts')
Idents(object = PINARQPos) <- "ARQPos"
Idents(object = PINARQNeg) <- "ARQNeg"
PINARQPos[["ARQExp"]] <- Idents(object = PINARQPos)
PINARQNeg[["ARQExp"]] <- Idents(object = PINARQNeg)
PINARQ <- merge(x = PINARQPos, y = PINARQNeg)
Idents(object = PINARQ) <- "ARQExp"
PIN2$ARQExp <- Idents(object = PINARQ)

TumorARQPos <- subset(x=Tumor2, subset = ARQ > 0, slot = 'counts')
TumorARQNeg <- subset(x=Tumor2, subset = ARQ == 0, slot = 'counts')
Idents(object = TumorARQPos) <- "ARQPos"
Idents(object = TumorARQNeg) <- "ARQNeg"
TumorARQPos[["ARQExp"]] <- Idents(object = TumorARQPos)
TumorARQNeg[["ARQExp"]] <- Idents(object = TumorARQNeg)
TumorARQ <- merge(x = TumorARQPos, y = TumorARQNeg)
Idents(object = TumorARQ) <- "ARQExp"
Tumor2$ARQExp <- Idents(object = TumorARQ)

#Merging Datasets

#Stash old idents
Tumor2[["orig.clusters"]] <- Idents(object = Tumor2)
PIN2[["orig.clusters"]] <- Idents(object = PIN2)

#Set Current idents
Idents(object = PIN2) <- "seurat_clusters"

PIN2$stim <- "PIN"
Tumor2$stim <- "Tumor"

PINvTumor.anchors <- FindIntegrationAnchors(object.list = list(PIN2, Tumor2), dims = 1:20)

PINvTumor.combined <- IntegrateData(anchorset = PINvTumor.anchors, dims = 1:20)

DefaultAssay(PINvTumor.combined) <- "integrated"

# Run the standard workflow for visualization and clustering
PINvTumor.combined <- ScaleData(PINvTumor.combined, verbose = FALSE)
PINvTumor.combined <- RunPCA(PINvTumor.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
PINvTumor.combined <- FindNeighbors(PINvTumor.combined, reduction = "pca", dims = 1:20)
PINvTumor.combined <- FindClusters(PINvTumor.combined, resolution = 0.5)
PINvTumor.combined <- RunTSNE(PINvTumor.combined, reduction = "pca", dims = 1:20)
p1 <- DimPlot(PINvTumor.combined, reduction = "tsne", label = TRUE)
p2 <- DimPlot(PINvTumor.combined, reduction = "tsne", group.by = "stim")
plot_grid(p1, p2)

DefaultAssay(PINvTumor.combined) <- "RNA"

PINvTumor.combined.markers <- FindAllMarkers(PINvTumor.combined, min.pct = 0.25, logfc.threshold = 0.25)
PINvTumor.combined.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(PINvTumor.combined.markers, file = "PINvTumorcombinedmarkers.txt")
PINvTumorTop5 <- PINvTumor.combined.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_logFC)


FeaturePlot(PINvTumor.combined, features = c("EGFP", "ARQ"), cols = c("light grey", "red"), min.cutoff = 0, label = TRUE)
FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c( "ARQ", "EGFP"), cols = c("light grey", "red"), min.cutoff = 0, label = TRUE, split.by = "stim")


```


```{r}

#Cell Type Groups

#Grouping Cells (SplitGroups)

Idents(object = PINvTumor.combined) <- "seurat_clusters"
PINvTumor.combined <- RenameIdents(object = PINvTumor.combined, '2' = "BE1", '15' = "BE2", '0' = "LE1", '1' = "LE2", '3' = "LE3", '4' = "LE4", '5' = "LE5", '6' = "LE6", '9' = "LE7", '10' = "LE8", '8' = "Fibro", '14' = "SM", '11' = "Leu", '13' = "Endo", '7' = "I1", '12' = "I2", '16' = "I3", '17' = "I4")
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 1)
PINvTumor.combined[["SplitCellClass"]] <- Idents(object = PINvTumor.combined)

DotPlot(PINvTumor.combined, features = c( "Cd2", "Cd28",	"Sele",	"Plvap",	"Cdh5", "Pecam1", "Cldn5", "Spi1",	"Tyrobp",	"C1qc",	 "C1qb",	"C1qa",	"Myl9", "Tagln",	"Rgs5",	 "Myh11",	"Acta2",	"Ptgs2", "Rspo3", "Fgf2",	 "Apod", "Fbln1", "Cd24a", "Alcam", "Krt19",	 "Krt8",	"Trp63",	"Col17a1", "Krt15",	 "Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()


#Combined Groups

Idents(object = PINvTumor.combined) <- "seurat_clusters"
PINvTumor.combined <- RenameIdents(object = PINvTumor.combined, '2' = "BE", '15' = "BE", '0' = "LE", '1' = "LE", '3' = "LE", '4' = "LE", '5' = "LE", '6' = "LE", '9' = "LE", '10' = "LE", '8' = "Fibro", '14' = "SM", '11' = "Leu", '13' = "Endo", '7' = "I", '12' = "I", '16' = "I", '17' = "I")
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 1)
PINvTumor.combined[["CellClass"]] <- Idents(object = PINvTumor.combined)
DotPlot(PINvTumor.combined, features = c( "Cd2", "Cd28",	"Sele",	"Plvap",	"Cdh5", "Pecam1", "Cldn5", "Spi1",	"Tyrobp",	"C1qc",	 "C1qb",	"C1qa",	"Myl9", "Tagln",	"Rgs5",	 "Myh11",	"Acta2",	"Ptgs2", "Rspo3", "Fgf2",	 "Apod", "Fbln1", "Cd24a", "Alcam", "Krt19",	 "Krt8",	"Trp63",	"Col17a1", "Krt15",	 "Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

Idents(object = PINvTumor.combined) <- "CellClass"
LE.markers <- FindMarkers(PINvTumor.combined, ident.1 = "LE", min.pct = 0.25, logfc.threshold = 0.25)
LE.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.table(LE.markers, file = "LEMarkers.txt")


```

```{r}

#Separate By ARQExp

PINvTumor.combined$SplitCellClass.ARQ <- paste(Idents(PINvTumor.combined), PINvTumor.combined$ARQExp, sep = "_")
PINvTumor.combined$SplitCellClass <- Idents(PINvTumor.combined)
Idents(PINvTumor.combined) <- "SplitCellClass.ARQ"