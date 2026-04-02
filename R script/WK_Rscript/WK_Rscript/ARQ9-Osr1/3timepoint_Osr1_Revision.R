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

### E18.5 ###
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/E18.5_only_Analysis.RData")
load("W:/ARKO-Gli1_3timepoint_scRNAseq_WK/E17.5/E17_Ctrl_ARKO.RData")
load("//isi-dcnl/user_data/zjsun/group/Adam Olson/AO Files/ARKO Gli1/Single Cell/WT/WTE17 Single Cell Workspace.RData")

### P11 ###
load("//isi-dcnl/user_data/zjsun/group/Sun_Lab_Sequencing_Data/Gli1-ARKO_3timepoints/SunLab_Analysis/P11_only_Analysis.RData")

### P35 ###
load("E:/DHL - SC analysis/P3P6to P35/P35_32414")

#### Initial Filtering and Clustering ####

#Setup workspace to make file calling & saving easy

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Osr1")

#E18.5

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

E181.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/counts/36982_6-E18-WT-Gli1/outs/filtered_feature_bc_matrix")
E181 <- CreateSeuratObject(counts = E181.data,  min.cells = 3, min.features = 200, project = "E181")
E181 <- NormalizeData(E181)

E182.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/20190417_Tgen/28126_scRNA/counts/28126_count/outs/filtered_feature_bc_matrix")
E182 <- CreateSeuratObject(counts = E182.data,  min.cells = 3, min.features = 200, project = "E182")
E182 <- NormalizeData(E182)

E183.data <- Read10X("//isi-dcnl/user_data/zjsun/group/Adam Olson/AO Files/ARKO Gli1/Single Cell/WT/filtered_feature_bc_matrix")
E183 <- CreateSeuratObject(counts = E183.data,  min.cells = 3, min.features = 200, project = "E183")
E183 <- NormalizeData(E183)

P111.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191214_32848_33580_test/counts/32848_ARCreER_P9/outs/filtered_feature_bc_matrix")
P111 <- CreateSeuratObject(counts = P111.data,  min.cells = 3, min.features = 200, project = "P111")
P111 <- NormalizeData(P111)

P112.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/counts/36977_1-P11-WT-Gli1/outs/filtered_feature_bc_matrix")
P112 <- CreateSeuratObject(counts = P112.data,  min.cells = 3, min.features = 200, project = "P112")
P112 <- NormalizeData(P112)

P113.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/190927/31923_ARCrER/counts/31923_ARCrER/outs/filtered_feature_bc_matrix")
P113 <- CreateSeuratObject(counts = P113.data,  min.cells = 3, min.features = 200, project = "P113")
P113 <- NormalizeData(P113)

P351.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/190219_scRNA/27347_WT_count_EGFPmm10/outs/filtered_feature_bc_matrix")
P351 <- CreateSeuratObject(counts = P351.data,  min.cells = 3, min.features = 200, project = "P351")
P351 <- NormalizeData(P351)

P352.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/191011_32414_test/32414_ARCreER_P35/32414_ARCreER_P35/outs/filtered_feature_bc_matrix")
P352 <- CreateSeuratObject(counts = P352.data,  min.cells = 3, min.features = 200, project = "P352")
P352 <- NormalizeData(P352)

P353.data <- Read10X("//isi-dcnl/user_data/zjsun/seq/200513/36980_4-P42-WT-Gli1/36980_count_ref5/outs/filtered_feature_bc_matrix")
P353 <- CreateSeuratObject(counts = P353.data,  min.cells = 3, min.features = 200, project = "P353")
P353 <- NormalizeData(P353)

#Initial processing & filtering

E181[["percent.mt"]] <- PercentageFeatureSet(E181, pattern = "^mt-")
VlnPlot(E181, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(E181@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E181 Pre-filteration")

hist(E181@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E181 Pre-filteration")

E1811 <- subset(E181, subset = nFeature_RNA > 500 & nFeature_RNA < 4000 & percent.mt < 10)
VlnPlot(E1811, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(E1811@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E1811 Post-filteration")

hist(E1811@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E1811 Post-filteration")

E1811 <- FindVariableFeatures(E1811, selection.method = "vst", nfeatures = 5000)

all.genes <- rownames(E1811)
E1811 <- ScaleData(E1811, features = all.genes)
E1811 <- RunPCA(E1811, features = VariableFeatures(object = E1811))
print(E1811[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E1811, dims = 1:2, reduction = "pca")
DimPlot(E1811, reduction = "pca")
ElbowPlot(E1811, ndims = 50)

E1811 <- FindNeighbors(E1811, dims = 1:20)
E1811 <- FindClusters(E1811, resolution = 0.5)
head(Idents(E1811), 5)
E1811 <- RunTSNE(E1811, dims = 1:20)
E1811 <- RunUMAP(E1811, dims = 1:20)
DimPlot(E1811, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(E1811) <- "RNA"
E1811Osr1Pos <- subset(x=E1811, subset = Osr1 > 0)
E1811Osr1Neg <- subset(x=E1811, subset = Osr1 == 0)
Idents(object = E1811Osr1Pos) <- "Osr1Pos"
Idents(object = E1811Osr1Neg) <- "Osr1Neg"
E1811Osr1Pos[["Osr1Exp"]] <- Idents(object = E1811Osr1Pos)
E1811Osr1Neg[["Osr1Exp"]] <- Idents(object = E1811Osr1Neg)
E1811Osr1 <- merge(x = E1811Osr1Pos, y = E1811Osr1Neg)
Idents(object = E1811Osr1) <- "Osr1Exp"
E1811$Osr1Exp <- Idents(object = E1811Osr1)

#Add Prom1 info
DefaultAssay(E1811) <- "RNA"
E1811Prom1Pos <- subset(x=E1811, subset = Prom1 > 0)
E1811Prom1Neg <- subset(x=E1811, subset = Prom1 == 0)
Idents(object = E1811Prom1Pos) <- "Prom1Pos"
Idents(object = E1811Prom1Neg) <- "Prom1Neg"
E1811Prom1Pos[["Prom1Exp"]] <- Idents(object = E1811Prom1Pos)
E1811Prom1Neg[["Prom1Exp"]] <- Idents(object = E1811Prom1Neg)
E1811Prom1 <- merge(x = E1811Prom1Pos, y = E1811Prom1Neg)
Idents(object = E1811Prom1) <- "Prom1Exp"
E1811$Prom1Exp <- Idents(object = E1811Prom1)

Idents(object = E1811) <- "Osr1Exp"
E1811$Osr1Exp.Prom1Exp <- paste(Idents(E1811), E1811$Prom1Exp, sep = "_")
Idents(object = E1811) <- "Osr1Exp.Prom1Exp"
table(Idents(E1811))

#Add Itga4 info
DefaultAssay(E1811) <- "RNA"
E1811Itga4Pos <- subset(x=E1811, subset = Itga4 > 0)
E1811Itga4Neg <- subset(x=E1811, subset = Itga4 == 0)
Idents(object = E1811Itga4Pos) <- "Itga4Pos"
Idents(object = E1811Itga4Neg) <- "Itga4Neg"
E1811Itga4Pos[["Itga4Exp"]] <- Idents(object = E1811Itga4Pos)
E1811Itga4Neg[["Itga4Exp"]] <- Idents(object = E1811Itga4Neg)
E1811Itga4 <- merge(x = E1811Itga4Pos, y = E1811Itga4Neg)
Idents(object = E1811Itga4) <- "Itga4Exp"
E1811$Itga4Exp <- Idents(object = E1811Itga4)

Idents(object = E1811) <- "Osr1Exp"
E1811$Osr1Exp.Itga4Exp <- paste(Idents(E1811), E1811$Itga4Exp, sep = "_")
Idents(object = E1811) <- "Osr1Exp.Itga4Exp"
table(Idents(E1811))

#Add Ly6a info
DefaultAssay(E1811) <- "RNA"
E1811Ly6aPos <- subset(x=E1811, subset = Ly6a > 0)
E1811Ly6aNeg <- subset(x=E1811, subset = Ly6a == 0)
Idents(object = E1811Ly6aPos) <- "Ly6aPos"
Idents(object = E1811Ly6aNeg) <- "Ly6aNeg"
E1811Ly6aPos[["Ly6aExp"]] <- Idents(object = E1811Ly6aPos)
E1811Ly6aNeg[["Ly6aExp"]] <- Idents(object = E1811Ly6aNeg)
E1811Ly6a <- merge(x = E1811Ly6aPos, y = E1811Ly6aNeg)
Idents(object = E1811Ly6a) <- "Ly6aExp"
E1811$Ly6aExp <- Idents(object = E1811Ly6a)

Idents(object = E1811) <- "Osr1Exp"
E1811$Osr1Exp.Ly6aExp <- paste(Idents(E1811), E1811$Ly6aExp, sep = "_")
Idents(object = E1811) <- "Osr1Exp.Ly6aExp"
table(Idents(E1811))

#Add Tacstd2 info
DefaultAssay(E1811) <- "RNA"
E1811Tacstd2Pos <- subset(x=E1811, subset = Tacstd2 > 0)
E1811Tacstd2Neg <- subset(x=E1811, subset = Tacstd2 == 0)
Idents(object = E1811Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = E1811Tacstd2Neg) <- "Tacstd2Neg"
E1811Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = E1811Tacstd2Pos)
E1811Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = E1811Tacstd2Neg)
E1811Tacstd2 <- merge(x = E1811Tacstd2Pos, y = E1811Tacstd2Neg)
Idents(object = E1811Tacstd2) <- "Tacstd2Exp"
E1811$Tacstd2Exp <- Idents(object = E1811Tacstd2)

Idents(object = E1811) <- "Osr1Exp"
E1811$Osr1Exp.Tacstd2Exp <- paste(Idents(E1811), E1811$Tacstd2Exp, sep = "_")
Idents(object = E1811) <- "Osr1Exp.Tacstd2Exp"
table(Idents(E1811))

#Add Trp63 info
DefaultAssay(E1811) <- "RNA"
E1811Trp63Pos <- subset(x=E1811, subset = Trp63 > 0)
E1811Trp63Neg <- subset(x=E1811, subset = Trp63 == 0)
Idents(object = E1811Trp63Pos) <- "Trp63Pos"
Idents(object = E1811Trp63Neg) <- "Trp63Neg"
E1811Trp63Pos[["Trp63Exp"]] <- Idents(object = E1811Trp63Pos)
E1811Trp63Neg[["Trp63Exp"]] <- Idents(object = E1811Trp63Neg)
E1811Trp63 <- merge(x = E1811Trp63Pos, y = E1811Trp63Neg)
Idents(object = E1811Trp63) <- "Trp63Exp"
E1811$Trp63Exp <- Idents(object = E1811Trp63)

Idents(object = E1811) <- "Osr1Exp"
E1811$Osr1Exp.Trp63Exp <- paste(Idents(E1811), E1811$Trp63Exp, sep = "_")
Idents(object = E1811) <- "Osr1Exp.Trp63Exp"
table(Idents(E1811))

#Add Psca info
DefaultAssay(E1811) <- "RNA"
E1811PscaPos <- subset(x=E1811, subset = Psca > 0)
E1811PscaNeg <- subset(x=E1811, subset = Psca == 0)
Idents(object = E1811PscaPos) <- "PscaPos"
Idents(object = E1811PscaNeg) <- "PscaNeg"
E1811PscaPos[["PscaExp"]] <- Idents(object = E1811PscaPos)
E1811PscaNeg[["PscaExp"]] <- Idents(object = E1811PscaNeg)
E1811Psca <- merge(x = E1811PscaPos, y = E1811PscaNeg)
Idents(object = E1811Psca) <- "PscaExp"
E1811$PscaExp <- Idents(object = E1811Psca)

Idents(object = E1811) <- "Osr1Exp"
E1811$Osr1Exp.PscaExp <- paste(Idents(E1811), E1811$PscaExp, sep = "_")
Idents(object = E1811) <- "Osr1Exp.PscaExp"
table(Idents(E1811))

#Initial processing & filtering

E182[["percent.mt"]] <- PercentageFeatureSet(E182, pattern = "^mt-")
VlnPlot(E182, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(E182@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E182 Pre-filteration")

hist(E182@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E182 Pre-filteration")

plot1 <- FeatureScatter(E182, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E182, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E1821 <- subset(E182, subset = nFeature_RNA> 450 & nFeature_RNA < 8000 & percent.mt < 10)
VlnPlot(E1821, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(E1821@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E1821 Post-filteration")

hist(E1821@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E1821 Post-filteration")

plot1 <- FeatureScatter(E1821, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E1821, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E1821 <- FindVariableFeatures(E1821, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(E1821)
E1821 <- ScaleData(E1821, features = all.genes)
E1821 <- RunPCA(E1821, features = VariableFeatures(object = E1811))
print(E1821[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E1821, dims = 1:2, reduction = "pca")
DimPlot(E1821, reduction = "pca")
ElbowPlot(E1821, ndims = 50)

E1821 <- FindNeighbors(E1821, dims = 1:20)
E1821 <- FindClusters(E1821, resolution = 0.5)
head(Idents(E1821), 5)
E1821 <- RunTSNE(E1821, dims = 1:20)
E1821 <- RunUMAP(E1821, dims = 1:20)
DimPlot(E1821, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(E1821) <- "RNA"
E1821Osr1Pos <- subset(x=E1821, subset = Osr1 > 0)
E1821Osr1Neg <- subset(x=E1821, subset = Osr1 == 0)
Idents(object = E1821Osr1Pos) <- "Osr1Pos"
Idents(object = E1821Osr1Neg) <- "Osr1Neg"
E1821Osr1Pos[["Osr1Exp"]] <- Idents(object = E1821Osr1Pos)
E1821Osr1Neg[["Osr1Exp"]] <- Idents(object = E1821Osr1Neg)
E1821Osr1 <- merge(x = E1821Osr1Pos, y = E1821Osr1Neg)
Idents(object = E1821Osr1) <- "Osr1Exp"
E1821$Osr1Exp <- Idents(object = E1821Osr1)

#Add Prom1 info
DefaultAssay(E1821) <- "RNA"
E1821Prom1Pos <- subset(x=E1821, subset = Prom1 > 0)
E1821Prom1Neg <- subset(x=E1821, subset = Prom1 == 0)
Idents(object = E1821Prom1Pos) <- "Prom1Pos"
Idents(object = E1821Prom1Neg) <- "Prom1Neg"
E1821Prom1Pos[["Prom1Exp"]] <- Idents(object = E1821Prom1Pos)
E1821Prom1Neg[["Prom1Exp"]] <- Idents(object = E1821Prom1Neg)
E1821Prom1 <- merge(x = E1821Prom1Pos, y = E1821Prom1Neg)
Idents(object = E1821Prom1) <- "Prom1Exp"
E1821$Prom1Exp <- Idents(object = E1821Prom1)

Idents(object = E1821) <- "Osr1Exp"
E1821$Osr1Exp.Prom1Exp <- paste(Idents(E1821), E1821$Prom1Exp, sep = "_")
Idents(object = E1821) <- "Osr1Exp.Prom1Exp"
table(Idents(E1821))

#Add Itga4 info
DefaultAssay(E1821) <- "RNA"
E1821Itga4Pos <- subset(x=E1821, subset = Itga4 > 0)
E1821Itga4Neg <- subset(x=E1821, subset = Itga4 == 0)
Idents(object = E1821Itga4Pos) <- "Itga4Pos"
Idents(object = E1821Itga4Neg) <- "Itga4Neg"
E1821Itga4Pos[["Itga4Exp"]] <- Idents(object = E1821Itga4Pos)
E1821Itga4Neg[["Itga4Exp"]] <- Idents(object = E1821Itga4Neg)
E1821Itga4 <- merge(x = E1821Itga4Pos, y = E1821Itga4Neg)
Idents(object = E1821Itga4) <- "Itga4Exp"
E1821$Itga4Exp <- Idents(object = E1821Itga4)

Idents(object = E1821) <- "Osr1Exp"
E1821$Osr1Exp.Itga4Exp <- paste(Idents(E1821), E1821$Itga4Exp, sep = "_")
Idents(object = E1821) <- "Osr1Exp.Itga4Exp"
table(Idents(E1821))

#Add Ly6a info
DefaultAssay(E1821) <- "RNA"
E1821Ly6aPos <- subset(x=E1821, subset = Ly6a > 0)
E1821Ly6aNeg <- subset(x=E1821, subset = Ly6a == 0)
Idents(object = E1821Ly6aPos) <- "Ly6aPos"
Idents(object = E1821Ly6aNeg) <- "Ly6aNeg"
E1821Ly6aPos[["Ly6aExp"]] <- Idents(object = E1821Ly6aPos)
E1821Ly6aNeg[["Ly6aExp"]] <- Idents(object = E1821Ly6aNeg)
E1821Ly6a <- merge(x = E1821Ly6aPos, y = E1821Ly6aNeg)
Idents(object = E1821Ly6a) <- "Ly6aExp"
E1821$Ly6aExp <- Idents(object = E1821Ly6a)

Idents(object = E1821) <- "Osr1Exp"
E1821$Osr1Exp.Ly6aExp <- paste(Idents(E1821), E1821$Ly6aExp, sep = "_")
Idents(object = E1821) <- "Osr1Exp.Ly6aExp"
table(Idents(E1821))

#Add Tacstd2 info
DefaultAssay(E1821) <- "RNA"
E1821Tacstd2Pos <- subset(x=E1821, subset = Tacstd2 > 0)
E1821Tacstd2Neg <- subset(x=E1821, subset = Tacstd2 == 0)
Idents(object = E1821Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = E1821Tacstd2Neg) <- "Tacstd2Neg"
E1821Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = E1821Tacstd2Pos)
E1821Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = E1821Tacstd2Neg)
E1821Tacstd2 <- merge(x = E1821Tacstd2Pos, y = E1821Tacstd2Neg)
Idents(object = E1821Tacstd2) <- "Tacstd2Exp"
E1821$Tacstd2Exp <- Idents(object = E1821Tacstd2)

Idents(object = E1821) <- "Osr1Exp"
E1821$Osr1Exp.Tacstd2Exp <- paste(Idents(E1821), E1821$Tacstd2Exp, sep = "_")
Idents(object = E1821) <- "Osr1Exp.Tacstd2Exp"
table(Idents(E1821))

#Add Trp63 info
DefaultAssay(E1821) <- "RNA"
E1821Trp63Pos <- subset(x=E1821, subset = Trp63 > 0)
E1821Trp63Neg <- subset(x=E1821, subset = Trp63 == 0)
Idents(object = E1821Trp63Pos) <- "Trp63Pos"
Idents(object = E1821Trp63Neg) <- "Trp63Neg"
E1821Trp63Pos[["Trp63Exp"]] <- Idents(object = E1821Trp63Pos)
E1821Trp63Neg[["Trp63Exp"]] <- Idents(object = E1821Trp63Neg)
E1821Trp63 <- merge(x = E1821Trp63Pos, y = E1821Trp63Neg)
Idents(object = E1821Trp63) <- "Trp63Exp"
E1821$Trp63Exp <- Idents(object = E1821Trp63)

Idents(object = E1821) <- "Osr1Exp"
E1821$Osr1Exp.Trp63Exp <- paste(Idents(E1821), E1821$Trp63Exp, sep = "_")
Idents(object = E1821) <- "Osr1Exp.Trp63Exp"
table(Idents(E1821))

#Add Psca info
DefaultAssay(E1821) <- "RNA"
E1821PscaPos <- subset(x=E1821, subset = Psca > 0)
E1821PscaNeg <- subset(x=E1821, subset = Psca == 0)
Idents(object = E1821PscaPos) <- "PscaPos"
Idents(object = E1821PscaNeg) <- "PscaNeg"
E1821PscaPos[["PscaExp"]] <- Idents(object = E1821PscaPos)
E1821PscaNeg[["PscaExp"]] <- Idents(object = E1821PscaNeg)
E1821Psca <- merge(x = E1821PscaPos, y = E1821PscaNeg)
Idents(object = E1821Psca) <- "PscaExp"
E1821$PscaExp <- Idents(object = E1821Psca)

Idents(object = E1821) <- "Osr1Exp"
E1821$Osr1Exp.PscaExp <- paste(Idents(E1821), E1821$PscaExp, sep = "_")
Idents(object = E1821) <- "Osr1Exp.PscaExp"
table(Idents(E1821))

#E183
#Initial processing & filtering

E183[["percent.mt"]] <- PercentageFeatureSet(E183, pattern = "^mt-")
VlnPlot(E183, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(E183@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E183 Pre-filteration")

hist(E183@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E183 Pre-filteration")

plot1 <- FeatureScatter(E183, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E183, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

E1831 <- subset(E183, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 10)
VlnPlot(E1831, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(E1831@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "E1831 Post-filteration")

hist(E1831@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "E1831 Post-filteration")

plot1 <- FeatureScatter(E1831, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(E1831, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

DefaultAssay(E1831) <- "RNA"
E1831 <- NormalizeData(E1831, verbose = FALSE)
E1831 <- FindVariableFeatures(E1831, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(E1831)
E1831 <- ScaleData(E1831, features = all.genes)
E1831 <- RunPCA(E1831, features = VariableFeatures(E1831))
print(E1831[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(E1831, dims = 1:2, reduction = "pca")
DimPlot(E1831, reduction = "pca")
ElbowPlot(E1831, ndims = 50)

E1831 <- FindNeighbors(E1831, dims = 1:20)
E1831 <- FindClusters(E1831, resolution = 0.5)
head(Idents(E1831), 5)
E1831 <- RunTSNE(E1831, dims = 1:20)
E1831 <- RunUMAP(E1831, dims = 1:20)
DimPlot(E1831, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(E1831) <- "RNA"
E1831Osr1Pos <- subset(x=E1831, subset = Osr1 > 0)
E1831Osr1Neg <- subset(x=E1831, subset = Osr1 == 0)
Idents(object = E1831Osr1Pos) <- "Osr1Pos"
Idents(object = E1831Osr1Neg) <- "Osr1Neg"
E1831Osr1Pos[["Osr1Exp"]] <- Idents(object = E1831Osr1Pos)
E1831Osr1Neg[["Osr1Exp"]] <- Idents(object = E1831Osr1Neg)
E1831Osr1 <- merge(x = E1831Osr1Pos, y = E1831Osr1Neg)
Idents(object = E1831Osr1) <- "Osr1Exp"
E1831$Osr1Exp <- Idents(object = E1831Osr1)

#Add Prom1 info
DefaultAssay(E1831) <- "RNA"
E1831Prom1Pos <- subset(x=E1831, subset = Prom1 > 0)
E1831Prom1Neg <- subset(x=E1831, subset = Prom1 == 0)
Idents(object = E1831Prom1Pos) <- "Prom1Pos"
Idents(object = E1831Prom1Neg) <- "Prom1Neg"
E1831Prom1Pos[["Prom1Exp"]] <- Idents(object = E1831Prom1Pos)
E1831Prom1Neg[["Prom1Exp"]] <- Idents(object = E1831Prom1Neg)
E1831Prom1 <- merge(x = E1831Prom1Pos, y = E1831Prom1Neg)
Idents(object = E1831Prom1) <- "Prom1Exp"
E1831$Prom1Exp <- Idents(object = E1831Prom1)

Idents(object = E1831) <- "Osr1Exp"
E1831$Osr1Exp.Prom1Exp <- paste(Idents(E1831), E1831$Prom1Exp, sep = "_")
Idents(object = E1831) <- "Osr1Exp.Prom1Exp"
table(Idents(E1831))

#Add Itga4 info
DefaultAssay(E1831) <- "RNA"
E1831Itga4Pos <- subset(x=E1831, subset = Itga4 > 0)
E1831Itga4Neg <- subset(x=E1831, subset = Itga4 == 0)
Idents(object = E1831Itga4Pos) <- "Itga4Pos"
Idents(object = E1831Itga4Neg) <- "Itga4Neg"
E1831Itga4Pos[["Itga4Exp"]] <- Idents(object = E1831Itga4Pos)
E1831Itga4Neg[["Itga4Exp"]] <- Idents(object = E1831Itga4Neg)
E1831Itga4 <- merge(x = E1831Itga4Pos, y = E1831Itga4Neg)
Idents(object = E1831Itga4) <- "Itga4Exp"
E1831$Itga4Exp <- Idents(object = E1831Itga4)

Idents(object = E1831) <- "Osr1Exp"
E1831$Osr1Exp.Itga4Exp <- paste(Idents(E1831), E1831$Itga4Exp, sep = "_")
Idents(object = E1831) <- "Osr1Exp.Itga4Exp"
table(Idents(E1831))

#Add Ly6a info
DefaultAssay(E1831) <- "RNA"
E1831Ly6aPos <- subset(x=E1831, subset = Ly6a > 0)
E1831Ly6aNeg <- subset(x=E1831, subset = Ly6a == 0)
Idents(object = E1831Ly6aPos) <- "Ly6aPos"
Idents(object = E1831Ly6aNeg) <- "Ly6aNeg"
E1831Ly6aPos[["Ly6aExp"]] <- Idents(object = E1831Ly6aPos)
E1831Ly6aNeg[["Ly6aExp"]] <- Idents(object = E1831Ly6aNeg)
E1831Ly6a <- merge(x = E1831Ly6aPos, y = E1831Ly6aNeg)
Idents(object = E1831Ly6a) <- "Ly6aExp"
E1831$Ly6aExp <- Idents(object = E1831Ly6a)

Idents(object = E1831) <- "Osr1Exp"
E1831$Osr1Exp.Ly6aExp <- paste(Idents(E1831), E1831$Ly6aExp, sep = "_")
Idents(object = E1831) <- "Osr1Exp.Ly6aExp"
table(Idents(E1831))

#Add Tacstd2 info
DefaultAssay(E1831) <- "RNA"
E1831Tacstd2Pos <- subset(x=E1831, subset = Tacstd2 > 0)
E1831Tacstd2Neg <- subset(x=E1831, subset = Tacstd2 == 0)
Idents(object = E1831Tacstd2Pos) <- "Tacstd2Pos"
Idents(object = E1831Tacstd2Neg) <- "Tacstd2Neg"
E1831Tacstd2Pos[["Tacstd2Exp"]] <- Idents(object = E1831Tacstd2Pos)
E1831Tacstd2Neg[["Tacstd2Exp"]] <- Idents(object = E1831Tacstd2Neg)
E1831Tacstd2 <- merge(x = E1831Tacstd2Pos, y = E1831Tacstd2Neg)
Idents(object = E1831Tacstd2) <- "Tacstd2Exp"
E1831$Tacstd2Exp <- Idents(object = E1831Tacstd2)

Idents(object = E1831) <- "Osr1Exp"
E1831$Osr1Exp.Tacstd2Exp <- paste(Idents(E1831), E1831$Tacstd2Exp, sep = "_")
Idents(object = E1831) <- "Osr1Exp.Tacstd2Exp"
table(Idents(E1831))

#Add Trp63 info
DefaultAssay(E1831) <- "RNA"
E1831Trp63Pos <- subset(x=E1831, subset = Trp63 > 0)
E1831Trp63Neg <- subset(x=E1831, subset = Trp63 == 0)
Idents(object = E1831Trp63Pos) <- "Trp63Pos"
Idents(object = E1831Trp63Neg) <- "Trp63Neg"
E1831Trp63Pos[["Trp63Exp"]] <- Idents(object = E1831Trp63Pos)
E1831Trp63Neg[["Trp63Exp"]] <- Idents(object = E1831Trp63Neg)
E1831Trp63 <- merge(x = E1831Trp63Pos, y = E1831Trp63Neg)
Idents(object = E1831Trp63) <- "Trp63Exp"
E1831$Trp63Exp <- Idents(object = E1831Trp63)

Idents(object = E1831) <- "Osr1Exp"
E1831$Osr1Exp.Trp63Exp <- paste(Idents(E1831), E1831$Trp63Exp, sep = "_")
Idents(object = E1831) <- "Osr1Exp.Trp63Exp"
table(Idents(E1831))

#Add Psca info
DefaultAssay(E1831) <- "RNA"
E1831PscaPos <- subset(x=E1831, subset = Psca > 0)
E1831PscaNeg <- subset(x=E1831, subset = Psca == 0)
Idents(object = E1831PscaPos) <- "PscaPos"
Idents(object = E1831PscaNeg) <- "PscaNeg"
E1831PscaPos[["PscaExp"]] <- Idents(object = E1831PscaPos)
E1831PscaNeg[["PscaExp"]] <- Idents(object = E1831PscaNeg)
E1831Psca <- merge(x = E1831PscaPos, y = E1831PscaNeg)
Idents(object = E1831Psca) <- "PscaExp"
E1831$PscaExp <- Idents(object = E1831Psca)

Idents(object = E1831) <- "Osr1Exp"
E1831$Osr1Exp.PscaExp <- paste(Idents(E1831), E1831$PscaExp, sep = "_")
Idents(object = E1831) <- "Osr1Exp.PscaExp"
table(Idents(E1831))


#P111
#Initial processing & filtering

P111[["percent.mt"]] <- PercentageFeatureSet(P111, pattern = "^mt-")
VlnPlot(P111, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P111@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P111 Pre-filteration")

hist(P111@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P111 Pre-filteration")

P1111 <- subset(P111, subset =  nFeature_RNA < 5000 & percent.mt < 10)
VlnPlot(P1111, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P1111@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P1111 Post-filteration")

hist(P1111@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P1111 Post-filteration")

plot1 <- FeatureScatter(P1111, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P1111, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P1111 <- FindVariableFeatures(P1111, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P1111)
P1111 <- ScaleData(P1111, features = all.genes)
P1111 <- RunPCA(P1111, features = VariableFeatures(object = E1811))
print(P1111[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P1111, dims = 1:2, reduction = "pca")
DimPlot(P1111, reduction = "pca")
ElbowPlot(P1111, ndims = 50)

P1111 <- FindNeighbors(P1111, dims = 1:20)
P1111 <- FindClusters(P1111, resolution = 0.5)
head(Idents(P1111), 5)
P1111 <- RunTSNE(P1111, dims = 1:20)
P1111 <- RunUMAP(P1111, dims = 1:20)
DimPlot(P1111, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(P1111) <- "RNA"
P1111Osr1Pos <- subset(x=P1111, subset = Osr1 > 0)
P1111Osr1Neg <- subset(x=P1111, subset = Osr1 == 0)
Idents(object = P1111Osr1Pos) <- "Osr1Pos"
Idents(object = P1111Osr1Neg) <- "Osr1Neg"
P1111Osr1Pos[["Osr1Exp"]] <- Idents(object = P1111Osr1Pos)
P1111Osr1Neg[["Osr1Exp"]] <- Idents(object = P1111Osr1Neg)
P1111Osr1 <- merge(x = P1111Osr1Pos, y = P1111Osr1Neg)
Idents(object = P1111Osr1) <- "Osr1Exp"
P1111$Osr1Exp <- Idents(object = P1111Osr1)

Idents(object = P1111) <- "Osr1Exp"
table(Idents(P1111))

#P112
#Initial processing & filtering

P112[["percent.mt"]] <- PercentageFeatureSet(P112, pattern = "^mt-")
VlnPlot(P112, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P112@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P112 Pre-filteration")

hist(P112@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P112 Pre-filteration")

P1121 <- subset(P112, subset =  nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(P1121, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P1121@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P1121 Post-filteration")

hist(P1121@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P1121 Post-filteration")

plot1 <- FeatureScatter(P1121, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P1121, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

DefaultAssay(P1121) <- "RNA"
P1121 <- NormalizeData(P1121, verbose = FALSE)
P1121 <- FindVariableFeatures(P1121, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P1121)
P1121 <- ScaleData(P1121, features = all.genes)
P1121 <- RunPCA(P1121, features = VariableFeatures(P1121))
print(P1121[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P1121, dims = 1:2, reduction = "pca")
DimPlot(P1121, reduction = "pca")

ElbowPlot(P1121, ndims = 50)

P1121 <- FindNeighbors(P1121, dims = 1:20)
P1121 <- FindClusters(P1121, resolution = 0.5)
head(Idents(P1121), 5)
P1121 <- RunTSNE(P1121, dims = 1:20)
P1121 <- RunUMAP(P1121, dims = 1:20)
DimPlot(P1121, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(P1121) <- "RNA"
P1121Osr1Pos <- subset(x=P1121, subset = Osr1 > 0)
P1121Osr1Neg <- subset(x=P1121, subset = Osr1 == 0)
Idents(object = P1121Osr1Pos) <- "Osr1Pos"
Idents(object = P1121Osr1Neg) <- "Osr1Neg"
P1121Osr1Pos[["Osr1Exp"]] <- Idents(object = P1121Osr1Pos)
P1121Osr1Neg[["Osr1Exp"]] <- Idents(object = P1121Osr1Neg)
P1121Osr1 <- merge(x = P1121Osr1Pos, y = P1121Osr1Neg)
Idents(object = P1121Osr1) <- "Osr1Exp"
P1121$Osr1Exp <- Idents(object = P1121Osr1)

Idents(object = P1121) <- "Osr1Exp"
table(Idents(P1121))


#P113
#Initial processing & filtering

P113[["percent.mt"]] <- PercentageFeatureSet(P113, pattern = "^mt-")
VlnPlot(P113, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P113@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P113 Pre-filteration")

hist(P113@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P113 Pre-filteration")

P1131 <- subset(P113, subset =  nFeature_RNA < 2000 & percent.mt < 10)
VlnPlot(P1131, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P1131@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P1131 Post-filteration")

hist(P1131@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P1131 Post-filteration")

plot1 <- FeatureScatter(P1131, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P1131, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P1131 <- FindVariableFeatures(P1131, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P1131)
P1131 <- ScaleData(P1131, features = all.genes)
P1131 <- RunPCA(P1131, features = VariableFeatures(object = E1811))
print(P1131[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P1131, dims = 1:2, reduction = "pca")
DimPlot(P1131, reduction = "pca")
ElbowPlot(P1131, ndims = 50)

P1131 <- FindNeighbors(P1131, dims = 1:20)
P1131 <- FindClusters(P1131, resolution = 0.5)
head(Idents(P1131), 5)
P1131 <- RunTSNE(P1131, dims = 1:20)
P1131 <- RunUMAP(P1131, dims = 1:20)
DimPlot(P1131, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(P1131) <- "RNA"
P1131Osr1Pos <- subset(x=P1131, subset = Osr1 > 0)
P1131Osr1Neg <- subset(x=P1131, subset = Osr1 == 0)
Idents(object = P1131Osr1Pos) <- "Osr1Pos"
Idents(object = P1131Osr1Neg) <- "Osr1Neg"
P1131Osr1Pos[["Osr1Exp"]] <- Idents(object = P1131Osr1Pos)
P1131Osr1Neg[["Osr1Exp"]] <- Idents(object = P1131Osr1Neg)
P1131Osr1 <- merge(x = P1131Osr1Pos, y = P1131Osr1Neg)
Idents(object = P1131Osr1) <- "Osr1Exp"
P1131$Osr1Exp <- Idents(object = P1131Osr1)

Idents(object = P1131) <- "Osr1Exp"
table(Idents(P1131))

#P351
#Initial processing & filtering

P351[["percent.mt"]] <- PercentageFeatureSet(P351, pattern = "^mt-")
VlnPlot(P351, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P351@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P351 Pre-filteration")

hist(P351@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P351 Pre-filteration")

plot1 <- FeatureScatter(P351, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P351, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P3511 <- subset(P351, subset = nFeature_RNA > 500 & nFeature_RNA < 8500 & percent.mt < 15)
VlnPlot(P3511, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P3511@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P3511 Post-filteration")

hist(P3511@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P3511 Post-filteration")

plot1 <- FeatureScatter(P3511, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3511, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P3511 <- FindVariableFeatures(P3511, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P3511)
P3511 <- ScaleData(P3511, features = all.genes)
P3511 <- RunPCA(P3511, features = VariableFeatures(object = E1811))
print(P3511[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P3511, dims = 1:2, reduction = "pca")
DimPlot(P3511, reduction = "pca")
ElbowPlot(P3511, ndims = 50)

P3511 <- FindNeighbors(P3511, dims = 1:20)
P3511 <- FindClusters(P3511, resolution = 0.5)
head(Idents(P3511), 5)
P3511 <- RunTSNE(P3511, dims = 1:20)
P3511 <- RunUMAP(P3511, dims = 1:20)
DimPlot(P3511, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(P3511) <- "RNA"
P3511Osr1Pos <- subset(x=P3511, subset = Osr1 > 0)
P3511Osr1Neg <- subset(x=P3511, subset = Osr1 == 0)
Idents(object = P3511Osr1Pos) <- "Osr1Pos"
Idents(object = P3511Osr1Neg) <- "Osr1Neg"
P3511Osr1Pos[["Osr1Exp"]] <- Idents(object = P3511Osr1Pos)
P3511Osr1Neg[["Osr1Exp"]] <- Idents(object = P3511Osr1Neg)
P3511Osr1 <- merge(x = P3511Osr1Pos, y = P3511Osr1Neg)
Idents(object = P3511Osr1) <- "Osr1Exp"
P3511$Osr1Exp <- Idents(object = P3511Osr1)

Idents(object = P3511) <- "Osr1Exp"
table(Idents(P3511))

#P352
#Initial processing & filtering

P352[["percent.mt"]] <- PercentageFeatureSet(P352, pattern = "^mt-")
VlnPlot(P352, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P352@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P352 Pre-filteration")

hist(P352@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P352 Pre-filteration")

plot1 <- FeatureScatter(P352, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P352, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P3521 <- subset(P352, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 & percent.mt < 15)
VlnPlot(P3521, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P3521@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P3521 Post-filteration")

hist(P3521@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P3521 Post-filteration")

plot1 <- FeatureScatter(P3521, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3521, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

DefaultAssay(P3521) <- "RNA"
P3521 <- NormalizeData(P3521, verbose = FALSE)
P3521 <- FindVariableFeatures(P3521, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P3521)
P3521 <- ScaleData(P3521, features = all.genes)
P3521 <- RunPCA(P3521, features = VariableFeatures(P3521))
print(P3521[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P3521, dims = 1:2, reduction = "pca")
DimPlot(P3521, reduction = "pca")

ElbowPlot(P3521, ndims = 50)

P3521 <- FindNeighbors(P3521, dims = 1:20)
P3521 <- FindClusters(P3521, resolution = 0.5)
head(Idents(P3521), 5)
P3521 <- RunTSNE(P3521, dims = 1:20)
P3521 <- RunUMAP(P3521, dims = 1:20)
DimPlot(P3521, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(P3521) <- "RNA"
P3521Osr1Pos <- subset(x=P3521, subset = Osr1 > 0)
P3521Osr1Neg <- subset(x=P3521, subset = Osr1 == 0)
Idents(object = P3521Osr1Pos) <- "Osr1Pos"
Idents(object = P3521Osr1Neg) <- "Osr1Neg"
P3521Osr1Pos[["Osr1Exp"]] <- Idents(object = P3521Osr1Pos)
P3521Osr1Neg[["Osr1Exp"]] <- Idents(object = P3521Osr1Neg)
P3521Osr1 <- merge(x = P3521Osr1Pos, y = P3521Osr1Neg)
Idents(object = P3521Osr1) <- "Osr1Exp"
P3521$Osr1Exp <- Idents(object = P3521Osr1)

Idents(object = P3521) <- "Osr1Exp"
table(Idents(P3521))

#P353
#Initial processing & filtering

P353[["percent.mt"]] <- PercentageFeatureSet(P353, pattern = "^mt-")
VlnPlot(P353, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P353@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P353 Pre-filteration")

hist(P353@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P353 Pre-filteration")

plot1 <- FeatureScatter(P353, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P353, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P3531 <- subset(P353, subset = nFeature_RNA < 2500 & percent.mt < 15)
VlnPlot(P3531, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

hist(P3531@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "P3531 Post-filteration")

hist(P3531@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "P3531 Post-filteration")

plot1 <- FeatureScatter(P3531, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(P3531, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2

P3531 <- FindVariableFeatures(P3531, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(P3531)
P3531 <- ScaleData(P3531, features = all.genes)
P3531 <- RunPCA(P3531, features = VariableFeatures(object = E1811))
print(P3531[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(P3531, dims = 1:2, reduction = "pca")
DimPlot(P3531, reduction = "pca")
ElbowPlot(P3531, ndims = 50)

P3531 <- FindNeighbors(P3531, dims = 1:20)
P3531 <- FindClusters(P3531, resolution = 0.5)
head(Idents(P3531), 5)
P3531 <- RunTSNE(P3531, dims = 1:20)
P3531 <- RunUMAP(P3531, dims = 1:20)
DimPlot(P3531, reduction = "umap", pt.size = 1, label = TRUE)

#Add Osr1 info
DefaultAssay(P3531) <- "RNA"
P3531Osr1Pos <- subset(x=P3531, subset = Osr1 > 0)
P3531Osr1Neg <- subset(x=P3531, subset = Osr1 == 0)
Idents(object = P3531Osr1Pos) <- "Osr1Pos"
Idents(object = P3531Osr1Neg) <- "Osr1Neg"
P3531Osr1Pos[["Osr1Exp"]] <- Idents(object = P3531Osr1Pos)
P3531Osr1Neg[["Osr1Exp"]] <- Idents(object = P3531Osr1Neg)
P3531Osr1 <- merge(x = P3531Osr1Pos, y = P3531Osr1Neg)
Idents(object = P3531Osr1) <- "Osr1Exp"
P3531$Osr1Exp <- Idents(object = P3531Osr1)

Idents(object = P3531) <- "Osr1Exp"
table(Idents(P3531))

#### Osr1 & Stem cell markers ####

#E1821 (E17.5-1)
Idents(object = E1821) <- "seurat_clusters"

#Add Osr1 info
DefaultAssay(E1821) <- "RNA"
E1821Osr1Pos <- subset(x=E1821, subset = Osr1 > 0)
E1821Osr1Neg <- subset(x=E1821, subset = Osr1 == 0)
Idents(object = E1821Osr1Pos) <- "Osr1Pos"
Idents(object = E1821Osr1Neg) <- "Osr1Neg"
E1821Osr1Pos[["Osr1Exp"]] <- Idents(object = E1821Osr1Pos)
E1821Osr1Neg[["Osr1Exp"]] <- Idents(object = E1821Osr1Neg)
E1821Osr1 <- merge(x = E1821Osr1Pos, y = E1821Osr1Neg)
Idents(object = E1821Osr1) <- "Osr1Exp"
E1821$Osr1Exp <- Idents(object = E1821Osr1)

Idents(object = E1821) <- "Osr1Exp"
E1822 <- subset(E1821, idents = c("Osr1Pos", "Osr1Neg"))
E1822 <- RenameIdents(object = E1822, 'Osr1Neg' = "Osr1Neg", 'Osr1Pos' = "Osr1Pos")
DimPlot(E1822, reduction = "umap", pt.size = 0.3)

#DEGs
DefaultAssay(E1822) <- "RNA"
E1822.0.1.Markers <- FindMarkers(E1822, ident.1 = "Osr1Pos", ident.2 = "Osr1Neg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E1822.0.1.Markers, "E1822.0.1.Markers.csv")

DefaultAssay(E1821) <- "RNA"
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Oct4"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
tiff(file = "PINonly EGFP Myh11 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonly, reduction = "umap", features = c("EGFP", "Myh11"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()

#E1831 (E17.5-2)
Idents(object = E1831) <- "seurat_clusters"

#Add Osr1 info
DefaultAssay(E1831) <- "RNA"
E1831Osr1Pos <- subset(x=E1831, subset = Osr1 > 0)
E1831Osr1Neg <- subset(x=E1831, subset = Osr1 == 0)
Idents(object = E1831Osr1Pos) <- "Osr1Pos"
Idents(object = E1831Osr1Neg) <- "Osr1Neg"
E1831Osr1Pos[["Osr1Exp"]] <- Idents(object = E1831Osr1Pos)
E1831Osr1Neg[["Osr1Exp"]] <- Idents(object = E1831Osr1Neg)
E1831Osr1 <- merge(x = E1831Osr1Pos, y = E1831Osr1Neg)
Idents(object = E1831Osr1) <- "Osr1Exp"
E1831$Osr1Exp <- Idents(object = E1831Osr1)

Idents(object = E1831) <- "Osr1Exp"
E1832 <- subset(E1831, idents = c("Osr1Pos", "Osr1Neg"))
E1832 <- RenameIdents(object = E1832, 'Osr1Neg' = "Osr1Neg", 'Osr1Pos' = "Osr1Pos")
DimPlot(E1832, reduction = "umap", pt.size = 0.3)

#DEGs
DefaultAssay(E1832) <- "RNA"
E1832.0.1.Markers <- FindMarkers(E1832, ident.1 = "Osr1Pos", ident.2 = "Osr1Neg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(E1832.0.1.Markers, "E1832.0.1.Markers.csv")

DefaultAssay(E1831) <- "RNA"

tiff(file = "E1831 Osr1 Ly6a Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Ly6a"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Psca Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Psca"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Prom1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Prom1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Bmi1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Bmi1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Trop2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Tacstd2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 p63 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Trp63"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

tiff(file = "E1831 Osr1 Sox2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Sox2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Oct4 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Pou5f1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Tbx3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Tbx3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()


tiff(file = "E1831 Osr1 Zeb1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Zeb1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Twist1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Twist1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Twist2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Twist2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()


tiff(file = "E1831 Osr1 Ar Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Ar"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Gli1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Gli1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

tiff(file = "E1831 Osr1 Lef1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Lef1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Foxf1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Foxf1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Gata2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Gata2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Hoxd13 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Hoxd13"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1831 Osr1 Hoxa10 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Hoxa10"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

tiff(file = "E1831 Osr1 Irf6 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Irf6"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

tiff(file = "E1831 Osr1 Irf6 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Irf6"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

tiff(file = "E1831 Osr1 Msx2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1831, reduction = "umap", features = c("Osr1", "Msx2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

DefaultAssay(E1821) <- "RNA"

tiff(file = "E1821 Osr1 Ly6a Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Ly6a"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Psca Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Psca"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Prom1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Prom1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Bmi1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Bmi1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Trop2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Tacstd2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 p63 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Trp63"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

tiff(file = "E1821 Osr1 Sox2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Sox2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Oct4 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Oct4"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Tbx3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Tbx3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()


tiff(file = "E1821 Osr1 Zeb1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Zeb1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Twist1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Twist1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E1821 Osr1 Twist2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Twist2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()


tiff(file = "E1821 Osr1 Aldh1a1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E1821, reduction = "umap", features = c("Osr1", "Aldh1a1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

#### Merging Datasets ####

#Stash old idents
E1831[["orig.clusters"]] <- Idents(object = E1831)
P1121[["orig.clusters"]] <- Idents(object = P1121)
P3521[["orig.clusters"]] <- Idents(object = P3521)

#Set Current idents
Idents(object = E1831) <- "seurat_clusters"
Idents(object = P1121) <- "seurat_clusters"
Idents(object = P3521) <- "seurat_clusters"
E1831$stim <- "E18.5"
P1121$stim <- "P14"
P3521$stim <- "P35"

Osr1.combined.anchors <- FindIntegrationAnchors(object.list = list(E1831, P1121, P3521), dims = 1:20)
Osr1.combined <- IntegrateData(anchorset = Osr1.combined.anchors, dims = 1:20)
DefaultAssay(Osr1.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
Osr1.combined <- ScaleData(Osr1.combined, verbose = FALSE)
Osr1.combined <- RunPCA(Osr1.combined, npcs = 30, verbose = FALSE)
ElbowPlot(Osr1.combined, ndims = 50)
# t-SNE and Clustering
Osr1.combined <- FindNeighbors(Osr1.combined, reduction = "pca", dims = 1:16)
Osr1.combined <- FindClusters(Osr1.combined, resolution = 0.5)
Osr1.combined <- RunTSNE(Osr1.combined, reduction = "pca", dims = 1:16)
Osr1.combined <- RunUMAP(Osr1.combined, reduction = "pca", dims = 1:16)
DimPlot(Osr1.combined, reduction = "umap", pt.size = 0.3, group.by = "stim", cols = c("darkred", "darkblue", "grey")) 

#Cell identification
DefaultAssay(Osr1.combined) <- "RNA"

tiff(file = "Osr1.combined Epcam Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Vim Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Fbln1 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Myh11 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Myh11"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Krt5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Krt8 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Cdh5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Rgs1 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Rgs1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Mpz Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Mpz"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Tyrobp Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Tyrobp"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Rgs5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Pate4 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Pate4"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined Chga Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Chga"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#DEGs
DefaultAssay(Osr1.combined) <- "RNA"
Osr1.combined <- ScaleData(Osr1.combined, features = rownames(Osr1.combined))
Osr1.combined.markers <- FindAllMarkers(Osr1.combined, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(Osr1.combined.markers, "Osr1.combined.markers.csv")

#New labeling
Idents(object = Osr1.combined) <- "seurat_clusters"
Osr1.combined <- RenameIdents(object = Osr1.combined, '2' = "BE", '6' = "BE", '10' = "BE", '5' = "LE", '14' = "NE", '15' = "SV", '0' = "SM", '4' = "SM", '1' = "FB", '3' = "FB", '7' = "FB", '9' = "ProS", '8' = "Glia", '12' = "Leu", '13' = "Peri", '11' = "VE")
Osr1.combined[["CellType"]] <- Idents(object = Osr1.combined)
DimPlot(Osr1.combined, reduction = "umap", pt.size = 0.3)

tiff(file = "Osr1.combined CellType UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined, reduction = "umap", pt.size = 0.3) 
dev.off()
tiff(file = "Osr1.combined CellType split UMAP.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined, reduction = "umap", split.by = "stim", pt.size = 0.3) 
dev.off()
tiff(file = "Osr1.combined Osr1 split Exp.tiff", width = 18, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined, reduction = "umap", features = c("Osr1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 1.5, max.cutoff = "q90")
dev.off()

#Subset E18.5
Idents(object = Osr1.combined1) <- "stim"
E18only <- subset(Osr1.combined, idents = c("E18.5"))

DefaultAssay(E18only) <- "RNA"

tiff(file = "E18.5only Osr1 Prom1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Prom1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Irf6 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Irf6"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Msx2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Msx2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Foxf1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Foxf1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Tbx3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Tbx3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Hoxa10 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Hoxa10"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Hoxd13 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Hoxd13"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Ly6a Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Ly6a"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Tacstd2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Tacstd2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Psca Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Psca"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = .5, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Prom1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Prom1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Alcam Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Alcam"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Ly6a Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Ly6a"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Tacstd2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Tacstd2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Itga6 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Itga6"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Tbx3 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Tbx3"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Psca Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Psca"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Bmi1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Bmi1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Trp63 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Trp63"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 1, blend.threshold = 0.1)
dev.off()

#### Cell Cycle Regression ####
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")

s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

DefaultAssay(Osr1.combined) <- "RNA"
all.genes <- rownames(Osr1.combined)
Osr1.combined <- ScaleData(Osr1.combined, features = all.genes)

Osr1.combined <- CellCycleScoring(Osr1.combined, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

Idents(object = Osr1.combined) <- "Phase"
DimPlot(Osr1.combined, reduction = "umap")
tiff(file = "Osr1.combined Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined, reduction = "umap", pt.size = 0.3, cols = c("purple", "#FFD966","#1D762E"))
dev.off()

Idents(object = Osr1.combined) <- "Phase"
DimPlot(Osr1.combined, reduction = "umap")
tiff(file = "Osr1.combined Cell Cyle stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined, reduction = "umap", split.by = "stim", pt.size = 0.3, cols = c("purple", "#FFD966","#1D762E"))
dev.off()

#Take Cell cycle out-1 
mouse_cell_cycle_genes <- readRDS(file = "//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARKOvCtrl_3timepoint/E18.5/mouse_cell_cycle_genes.rds")

s.genes <- mouse_cell_cycle_genes$s.genes
g2m.genes <- mouse_cell_cycle_genes$g2m.genes

Osr1.combined1 <- Osr1.combined

DefaultAssay(Osr1.combined1) <- "RNA"
Osr1.combined1 <- FindVariableFeatures(Osr1.combined1, selection.method = "vst")
Osr1.combined1 <- ScaleData(Osr1.combined1, features = rownames(Osr1.combined1))
Osr1.combined1 <- RunPCA(Osr1.combined1, features = VariableFeatures(Osr1.combined1))

Osr1.combined1 <- CellCycleScoring(Osr1.combined1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
Osr1.combined1 <- ScaleData(Osr1.combined1, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Osr1.combined1))
Osr1.combined1 <- RunPCA(Osr1.combined1, features = VariableFeatures(Osr1.combined1))
print(Osr1.combined1[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Osr1.combined1, dims = 1:2, reduction = "pca")
DimPlot(Osr1.combined1, reduction = "pca")
ElbowPlot(Osr1.combined1, ndims = 50)
Osr1.combined1 <- FindNeighbors(Osr1.combined1, dims = 1:16)
Osr1.combined1 <- FindClusters(Osr1.combined1, resolution = 0.5)
Osr1.combined1 <- RunTSNE(Osr1.combined1, dims = 1:16)
Osr1.combined1 <- RunUMAP(Osr1.combined1, dims = 1:16)

Idents(object = Osr1.combined1) <- "Phase"
tiff(file = "Osr1.combined1 Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined1, reduction = "umap", pt.size = 0.3, cols = c("purple", "#FFD966","#1D762E"))
dev.off()
tiff(file = "Osr1.combined1 Cell Cyle stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined1, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("purple", "#FFD966","#1D762E"))
dev.off()

#Take Cell cycle out-2 
Osr1.combined2 <- Osr1.combined1
DefaultAssay(Osr1.combined2) <- "integrated"
Osr1.combined2 <- ScaleData(Osr1.combined2, vars.to.regress = c("S.Score", "G2M.Score"), features = rownames(Osr1.combined1))
Osr1.combined2 <- RunPCA(Osr1.combined2, features = VariableFeatures(Osr1.combined2))
ElbowPlot(Osr1.combined2, ndims = 30)
Osr1.combined2 <- FindNeighbors(Osr1.combined2, reduction = "pca", dims = 1:16)
Osr1.combined2 <- FindClusters(Osr1.combined2, resolution = 0.5)
Osr1.combined2 <- RunUMAP(Osr1.combined2, reduction = "pca", dims = 1:16)
Osr1.combined2 <- RunTSNE(Osr1.combined2, reduction = "pca", dims = 1:16)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3, label = TRUE)

Idents(object = Osr1.combined2) <- "Phase"
tiff(file = "Osr1.combined2 Cell Cyle UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3, cols = c("purple", "#FFD966","#1D762E"))
dev.off()
tiff(file = "Osr1.combined2 Cell Cyle stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3, split.by = "stim", cols = c("purple", "#FFD966","#1D762E"))
dev.off()
tiff(file = "Osr1.combined2 ProS highlight after cellcyleregression UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3, cols = c("grey", "grey", "grey", "grey", "grey", "grey", "purple", "grey", "grey", "grey", "grey"))
dev.off()

#Cell type identification
DefaultAssay(Osr1.combined2) <- "RNA"
tiff(file = "Osr1.combined2 Epcam Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Vim Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Vim"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Fbln1 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Fbln1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Myh11 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Myh11"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Krt5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Krt8 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Cdh5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Cdh5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Mpz Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Mpz"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Tyrobp Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Tyrobp"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Rgs5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Rgs5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Pate4 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Pate4"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "Osr1.combined2 Chga Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Chga"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#New labeling
Idents(object = Osr1.combined2) <- "seurat_clusters"
Osr1.combined2 <- RenameIdents(object = Osr1.combined2, '2' = "BE", '5' = "BE", '4' = "LE", '1' = "FB", '3' = "FB", '8' = "FB", '0' = "SM", '6' = "SM", '11' = "Peri", '10' = "Leu", '9' = "VE", '7' = "Gli1", '12' = "SV", '13' = "NE")
Osr1.combined2[["CellType"]] <- Idents(object = Osr1.combined2)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.3)

#UMAP
Idents(object = Osr1.combined2) <- "CellType"
tiff(file = "Osr1.combined2 CellType UMAP.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.5, cols = c("#1D762E", "#E06666",  "red", "#FF9933", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
dev.off()
tiff(file = "Osr1.combined2 CellType stim UMAP.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
DimPlot(Osr1.combined2, reduction = "umap", pt.size = 0.5, split.by = "stim", cols = c("#1D762E", "#E06666",  "red", "#FF9933", "purple", "#FFD966", "#28CC44", "#3399FF", "#3333FF", "#51CCC9"))
dev.off()

#Featureplots
DefaultAssay(Osr1.combined2) <- "RNA"
tiff(file = "Osr1.combined2 Osr1 stim Exp.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(Osr1.combined2, reduction = "umap", features = c("Osr1"), split.by = "stim", cols = c("light grey", "red"), pt.size = 0.5, max.cutoff = "q90")
dev.off()

#Dotplots
Idents(object = Osr1.combined2) <- "CellType"
tiff(file = "Osr1.combined2 CellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(Osr1.combined2, features = c("Krt14", "Krt15", "Krt17", "Krt5", "Krt6a", "Krt19", "Tspan1", "Krt8", "Clu", "Cldn3", "Igfbp3", "Fbln1", "Lum", "Bgn", "Col15a1", "Myl9", "Acta2", "Tagln", "Actg2", "Myh11", "Rgs5", "Ndufa4l2", "Apold1", "Mustn1",  "Lrrc32", "Cd74", "Ccl4", "C1qa", "C1qb", "Tyrobp", "Cldn5", "Cdh5", "Pecam1", "Plvap", "Aqp1", "Mpz", "Plp1", "Fabp7", "Cryab", "Gpm6b", "Pate4", "Svs2", "Svs4", "Wfdc15b", "Plac8", "Scg2", "Chgb", "Chga", "Syp", "Scg3"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Subset E18.5
Idents(object = Osr1.combined2) <- "stim"
E18only <- subset(Osr1.combined2, idents = c("E18.5"))

DefaultAssay(E18only) <- "RNA"

tiff(file = "E18.5only Osr1 Prom1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Prom1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Ly6a Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Ly6a"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Tacstd2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Tacstd2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Psca Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Psca"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Itga6 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Itga6"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
dev.off()
tiff(file = "E18.5only Osr1 Trp63 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(E18only, reduction = "umap", features = c("Osr1", "Trp63"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "olivedrab4", "tomato3"), pt.size = 2, blend.threshold = 0.1)
dev.off()
