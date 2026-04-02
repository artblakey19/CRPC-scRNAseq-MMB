####SU2C/PCF polyA####
library(reshape2)
library(vegan)
library(rgl)
library(gplots)
library(grid)
library(gridExtra)
library(GenomicFeatures)
library(ggplot2)
library(statmod)
library(edgeR)
library(Seurat)
library(SeuratWrappers)
library(dplyr)
library(tidyverse)

####load data####
setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF_revision")

#RNAseq data
countdata <- read.table(file = "SU2C_PCF_CRPC_mRNAseq_polyA.csv", sep = ",", header = T, row.names=1, as.is=T)

a <- as.matrix(countdata)
#Log Transform
b=log2(a/10+1)
c=log2(a+1)
mCRPC_log =  CreateSeuratObject(counts = c, project = "2019SU2C")

#Clustering
mCRPC_log <- FindVariableFeatures(mCRPC_log, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(mCRPC_log)
mCRPC_log <- ScaleData(mCRPC_log, features = all.genes)
mCRPC_log <- RunPCA(mCRPC_log, features = VariableFeatures(object = mCRPC_log))

ElbowPlot(mCRPC_log, ndims = 50)

mCRPC_log <- FindNeighbors(mCRPC_log, reduction = "pca", dims = 1:30)
mCRPC_log <- FindClusters(mCRPC_log, resolution = 0.5)
mCRPC_log <- RunUMAP(mCRPC_log, reduction = "pca", dims = 1:30)
DimPlot(mCRPC_log, reduction = "umap", pt.size = 0.3) 

#Add classification
ARNE_meta_log1 <- mCRPC_log@meta.data
ARNE_meta_log1<- rownames_to_column(ARNE_meta_log1)
write_csv(ARNE_meta_log1, "ARNE_meta_log1.csv")
ARNE_meta_log <- read_csv("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF_revision/ARNE_meta_log1_1.csv")
ARNE_meta_log <- column_to_rownames(ARNE_meta_log, var = "rowname")
mCRPC_log <- AddMetaData(mCRPC_log,ARNE_meta_log)
#Rename CellTypes
Idents(object = mCRPC_log) <- "ARNE"
DimPlot(mCRPC_log, reduction = "umap", pt.size = 0.3, label=TRUE)
mCRPC_log <- RenameIdents(object = mCRPC_log, 
                          'ARpos_NEneg'="ARPC",'ARneg_NEpos'="NEPC",'ARneg_NEneg'="DNPC",
                          'ARpos_NEpos'="APNP",'ARneg_NElow'="ANNL", 'ARlow_NEneg'="ALNN")  
mCRPC_log[["ARNE"]] <- Idents(object = mCRPC_log)
DimPlot(mCRPC_log, reduction = "umap", pt.size = 1, label=TRUE)

table(Idents(mCRPC_log))

####Subset ARPC NEPC DNPC####
Idents(object = mCRPC_log) <- "ARNE"
ARNEDNPC <- subset(mCRPC_log, idents = c("ARPC", "NEPC", "DNPC"))

#Heatmap for AR signature
ARlist <- list(c("NAP1L2", "TMPRSS2", "CHRNA2", "TARP", "STEAP4", "ALDH1A3", "PART1", "PLPP1", "KLK3", "KLK2", 'NKX3-1', "AR", "PMEPA1", "SLC45A3", "FKBP5"))
ARNEDNPC <- AddModuleScore(object = ARNEDNPC, features = ARlist, name = "ARDownstream") 
ARNEDNPC[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = ARNEDNPC, vars = 'ARDownstream1')))
all.genes <- rownames(ARNEDNPC)
ARNEDNPC <- ScaleData(ARNEDNPC, features = all.genes)
tiff(file = "ARNEDNPC ARscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = ARNEDNPC, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()
tiff(file = "ARNEDNPC ARscore Heatmap-1.tiff", width = 5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = ARNEDNPC, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#Heatmap for NE signature
NElist <- list(c('NKX2-1', "INSM1", "SOX2", "LMO3", "POU3F2", "ASCL1", "CHGA", "CHGB", "SCN3A", "CELF3", "ELAVL4", "PCSK1", "SRRM4", "SNAP25", "SYP", "CHRNB2", "ACTL6B", "SCG3", "ENO2"))
ARNEDNPC <- AddModuleScore(object = ARNEDNPC, features = NElist, name = "NEmodule") 
ARNEDNPC[['NEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = ARNEDNPC, vars = 'NEmodule1')))
all.genes <- rownames(ARNEDNPC)
ARNEDNPC <- ScaleData(ARNEDNPC, features = all.genes)
tiff(file = "ARNEDNPC NEscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = ARNEDNPC, features = 'NEmodule1', assay = 'NEmodule', slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#MET with AR downstream Heatmap
DefaultAssay(ARNEDNPC) <- "RNA"
ARNEDNPC <- ScaleData(ARNEDNPC, features = rownames(ARNEDNPC))
tiff(file = "ARNEDNPC MET with AR all downstream Heatmap.tiff", width = 5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "AR", "NAP1L2", "TMPRSS2", "CHRNA2", "TARP", "STEAP4", "ALDH1A3", "PART1", "PLPP1", "KLK3", "KLK2", 'NKX3-1', "PMEPA1", "SLC45A3", "FKBP5"), draw.lines = TRUE, size = 3, disp.max = 2, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with AR downstream Heatmap-1
DefaultAssay(ARNEDNPC) <- "RNA"
ARNEDNPC <- ScaleData(ARNEDNPC, features = rownames(ARNEDNPC))
tiff(file = "ARNEDNPC MET with AR all downstream Heatmap-1.tiff", width = 5, height = 4, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", 'NKX3-1', "ALDH1A3", "PMEPA1", "STEAP4"), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected HGF/MET signature
tiff(file = "ARNEDNPC HGF selected downstream Heatmap without HGF-1.tiff", width = 5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "AQP9","RUNX1", "SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3"), draw.lines = TRUE, size = 3, disp.max = 3.4, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected HGF/MET signature
tiff(file = "ARNEDNPC HGF selected downstream Heatmap without HGF-1.tiff", width = 5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "AQP9","RUNX1", "HGF",
                                 "ABCC3",  "SLC43A3",
                                 "GPNMB",
                                 "CHI3L1", "VNN1", "CXCL1",
                                 "IL1B", "CCL13", "IL1R2",
                                 "CLEC5A", "CCL8", "KYNU", "LY75"), draw.lines = TRUE, size = 3, disp.max = 3.4, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected HGF/MET signature final
tiff(file = "ARNEDNPC HGF selected downstream Heatmap final.tiff", width = 5, height = 2.7, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "CHI3L1", 
                                 "CXCL1", "IL1B","CCL13", "CLEC5A", "KYNU", "LY75", "SLC43A3"), draw.lines = TRUE, size = 2.8, disp.max = 3, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()



#MET with selected WNT signature
tiff(file = "ARNEDNPC WNT selected downstream Heatmap without HGF-1.tiff", width = 5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "LEF1", "TCF7L2", 
                                 "PLAUR", "MMP7", "LGR5", "FST", "BMP2"), draw.lines = TRUE, size = 2, disp.max = 2.5, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected WNT signature-1
tiff(file = "ARNEDNPC WNT selected downstream Heatmap without HGF-2.tiff", width = 5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", 
                                 "PLAUR", "MMP7","PLAU", "FST",
                                 "ID2", "PTTG1", "CD44",
                                 "S100A8", "JAK2"), draw.lines = TRUE, size = 2, disp.max = 3, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected WNT signature-2
tiff(file = "ARNEDNPC stanford Wnt Heatmap without HGF-3.tiff", width = 5, height = 10.5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET", "MYC",	"CCND1",	"TCF7",	"LEF1",	"PPARD",
                                 "JUN",	"FOSL1",	"PLAUR",	"MMP7",	"AXIN2",	"NRCAM",	
                                 "TCF4",	"GAST",	"CD44",	"EFNB1",	"CLDN1",	"BIRC5",
                                 "VEGFA",	"FGF18",	"ATOH1",	"MET",	"EDN1",	"MYCBP",	
                                 "L1CAM",	"ID2",	"JAG1",	"TIAM1",	"DKK1",	"FGF9",
                                 "FGF20",	"LGR5",	"SOX17",	"PTTG1",	"FZD7",	"FST",	
                                 "EN2",	"GJA1",	"STRA6",	"RHOU",	"TWIST1",	"TCF20",
                                 "CCN4",	"CDX4",	"SFRP2",	"EGFR",	"CTLA4",	"VCAN",	"TNFRSF19"
), draw.lines = TRUE, size = 2, disp.max = 2.8, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected WNT signature-3
tiff(file = "ARNEDNPC selected stanford Wnt Heatmap.tiff", width = 5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET",	"PLAUR",	"MMP7",	"CD44",	"CLDN1",
                                 "EDN1",	"ID2",	"PTTG1",	"FST",	"CTLA4"
), draw.lines = TRUE, size = 2, disp.max = 2.8, disp.min = -2.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected WNT signature-4
tiff(file = "ARNEDNPC selected FEVR Wnt Heatmap.tiff", width = 5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET",	"PLAUR",	"MMP7",	"FST",
                                 "PKLR", "HAL",
                                 "S100A8", "PCSK9", "FABP6", "SERPINB1",
                                 "TNFSF13B", "JAK2", "IRF1","TRIM34"
), draw.lines = TRUE, size = 2, disp.max = 2.8, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected WNT signature-final
tiff(file = "ARNEDNPC selected FEVR Wnt Heatmap final.tiff", width = 5, height = 2.6, units = "in", compression = "lzw", res =200)
DoHeatmap(ARNEDNPC, features = c("MET",	"PLAUR",	"MMP7",	
                                 "PKLR", "HAL", "PCSK9", "SERPINB1",
                                 "TNFSF13B", "JAK2"
), draw.lines = TRUE, size = 2, disp.max = 3, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

####DEGs####
#DEGs_DNPCvARPC
DefaultAssay(ARNEDNPC) <- "RNA"
all.genes <- rownames(ARNEDNPC)
ARNEDNPC <- ScaleData(ARNEDNPC, features = all.genes)
DNPCvARPC.Markers <- FindMarkers(ARNEDNPC, ident.1 = c("DNPC"), 
                                 ident.2 = c("ARPC"), min.pct = 0, logfc.threshold = 0)
write.csv(DNPCvARPC.Markers, "DNPCvARPC.Markers.csv")

#p.adjust
DEG_DNPCvARPC <- read.csv("DNPCvARPC.Markers.csv") 
DEG_DNPCvARPC_pvalue <- DEG_DNPCvARPC$p_val
DEG_DNPCvARPC_pvalue=as.numeric(DEG_DNPCvARPC_pvalue)
DEG_DNPCvARPC_BH = p.adjust(DEG_DNPCvARPC_pvalue, "BH")
write.csv(DEG_DNPCvARPC_BH, "DEG_DNPCvARPC_BH.csv")

#DEGs_NEPCvARPC
DefaultAssay(ARNEDNPC) <- "RNA"
all.genes <- rownames(ARNEDNPC)
ARNEDNPC <- ScaleData(ARNEDNPC, features = all.genes)
NEPCvARPC.Markers <- FindMarkers(ARNEDNPC, ident.1 = c("NEPC"), 
                                 ident.2 = c("ARPC"), min.pct = 0, logfc.threshold = 0)
write.csv(NEPCvARPC.Markers, "NEPCvARPC.Markers.csv")

#p.adjust
DEG_NEPCvARPC <- read.csv("NEPCvARPC.Markers.csv") 
DEG_NEPCvARPC_pvalue <- DEG_NEPCvARPC$p_val
DEG_NEPCvARPC_pvalue=as.numeric(DEG_NEPCvARPC_pvalue)
DEG_NEPCvARPC_BH = p.adjust(DEG_NEPCvARPC_pvalue, "BH")
write.csv(DEG_NEPCvARPC_BH, "DEG_NEPCvARPC_BH.csv")

####DEGs_limma####
#DEGs_DNPCvARPC
DefaultAssay(ARNEDNPC) <- "RNA"
all.genes <- rownames(ARNEDNPC)
ARNEDNPC <- ScaleData(ARNEDNPC, features = all.genes)
DNPCvARPC.Markers1 <- FindMarkers(ARNEDNPC, ident.1 = c("DNPC"), 
                                 ident.2 = c("ARPC"),  min.pct = 0, logfc.threshold = 0)
write.csv(DNPCvARPC.Markers1, "DNPCvARPC.Markers1.csv")

#p.adjust
DEG_DNPCvARPC <- read.csv("DNPCvARPC.Markers1.csv") 
DEG_DNPCvARPC_pvalue <- DEG_DNPCvARPC$p_val
DEG_DNPCvARPC_pvalue=as.numeric(DEG_DNPCvARPC_pvalue)
DEG_DNPCvARPC_BH = p.adjust(DEG_DNPCvARPC_pvalue, "BH")
write.csv(DEG_DNPCvARPC_BH, "DEG_DNPCvARPC_BH1.csv")

#DEGs_NEPCvARPC
DefaultAssay(ARNEDNPC) <- "RNA"
all.genes <- rownames(ARNEDNPC)
ARNEDNPC <- ScaleData(ARNEDNPC, features = all.genes)
NEPCvARPC.Markers <- FindMarkers(ARNEDNPC, ident.1 = c("NEPC"), 
                                 ident.2 = c("ARPC"), min.pct = 0, logfc.threshold = 0)
write.csv(NEPCvARPC.Markers, "NEPCvARPC.Markers.csv")

#p.adjust
DEG_NEPCvARPC <- read.csv("NEPCvARPC.Markers.csv") 
DEG_NEPCvARPC_pvalue <- DEG_NEPCvARPC$p_val
DEG_NEPCvARPC_pvalue=as.numeric(DEG_NEPCvARPC_pvalue)
DEG_NEPCvARPC_BH = p.adjust(DEG_NEPCvARPC_pvalue, "BH")
write.csv(DEG_NEPCvARPC_BH, "DEG_NEPCvARPC_BH.csv")

####without log####
setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF_revision")
countdata <- read.table(file = "SU2C_PCF_CRPC_mRNAseq_polyA.csv", sep = ",", header = T, row.names=1, as.is=T)

mCRPC <- CreateSeuratObject(counts = countdata,  min.cells = 3, min.features = 200, project = "SU2C")

#Clustering
mCRPC <- FindVariableFeatures(mCRPC, selection.method = "vst", nfeatures = 5000)

mCRPC <- ScaleData(mCRPC, verbose = FALSE)
mCRPC <- RunPCA(mCRPC, npcs = 50, verbose = FALSE)
ElbowPlot(mCRPC, ndims = 50)

mCRPC <- FindNeighbors(mCRPC, reduction = "pca", dims = 1:20)
mCRPC <- FindClusters(mCRPC, resolution = 0.5)
mCRPC <- RunTSNE(mCRPC, reduction = "pca", dims = 1:20)
mCRPC <- RunUMAP(mCRPC, reduction = "pca", dims = 1:20)
DimPlot(mCRPC, reduction = "umap", pt.size = 0.3) 

#Add classification
ARNE_meta1 <- mCRPC@meta.data
ARNE_meta1<- rownames_to_column(ARNE_meta1)
write_csv(ARNE_meta1, "ARNE_meta1.csv")
ARNE_meta <- read_csv("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF_revision/ARNE_meta1_1.csv")
ARNE_meta <- column_to_rownames(ARNE_meta, var = "rowname")
mCRPC <- AddMetaData(mCRPC,ARNE_meta)
#Rename CellTypes
Idents(object = mCRPC) <- "ARNE"
DimPlot(mCRPC, reduction = "umap", pt.size = 0.3, label=TRUE)
mCRPC <- RenameIdents(object = mCRPC, 
                          'ARpos_NEneg'="ARPC",'ARneg_NEpos'="NEPC",'ARneg_NEneg'="DNPC",
                          'ARpos_NEpos'="APNP",'ARneg_NElow'="ANNL", 'ARlow_NEneg'="ALNN")  
mCRPC[["ARNE"]] <- Idents(object = mCRPC)
DimPlot(mCRPC, reduction = "umap", pt.size = 1, label=TRUE)

table(Idents(mCRPC))

####Subset ARPC NEPC DNPC####
setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF_revision/wo log")

Idents(object = mCRPC) <- "ARNE"
mCRPC_ARNEDNPC <- subset(mCRPC, idents = c("ARPC", "NEPC", "DNPC"))

#Heatmap for AR signature
ARlist <- list(c("NAP1L2", "TMPRSS2", "CHRNA2", "TARP", "STEAP4", "ALDH1A3", "PART1", "PLPP1", "KLK3", "KLK2", 'NKX3-1', "AR", "PMEPA1", "SLC45A3", "FKBP5"))
mCRPC_ARNEDNPC <- AddModuleScore(object = mCRPC_ARNEDNPC, features = ARlist, name = "ARDownstream") 
mCRPC_ARNEDNPC[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC_ARNEDNPC, vars = 'ARDownstream1')))
all.genes <- rownames(mCRPC_ARNEDNPC)
mCRPC_ARNEDNPC <- ScaleData(mCRPC_ARNEDNPC, features = all.genes)
tiff(file = "mCRPC_ARNEDNPC ARscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_ARNEDNPC, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE, disp.max = 1, disp.min = -1) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()
tiff(file = "mCRPC_ARNEDNPC ARscore Heatmap-1.tiff", width = 5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_ARNEDNPC, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE, disp.max = 1, disp.min = -1) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#Heatmap for NE signature
NElist <- list(c('NKX2-1', "INSM1", "SOX2", "LMO3", "POU3F2", "ASCL1", "CHGA", "CHGB", "SCN3A", "CELF3", "ELAVL4", "PCSK1", "SRRM4", "SNAP25", "SYP", "CHRNB2", "ACTL6B", "SCG3", "ENO2"))
mCRPC_ARNEDNPC <- AddModuleScore(object = mCRPC_ARNEDNPC, features = NElist, name = "NEmodule") 
mCRPC_ARNEDNPC[['NEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC_ARNEDNPC, vars = 'NEmodule1')))
all.genes <- rownames(mCRPC_ARNEDNPC)
mCRPC_ARNEDNPC <- ScaleData(mCRPC_ARNEDNPC, features = all.genes)
tiff(file = "mCRPC_ARNEDNPC NEscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_ARNEDNPC, features = 'NEmodule1', assay = 'NEmodule', slot = 'data', draw.lines = TRUE, disp.max = 1, disp.min = -1) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#MET with AR downstream Heatmap
DefaultAssay(mCRPC_ARNEDNPC) <- "RNA"
mCRPC_ARNEDNPC <- ScaleData(mCRPC_ARNEDNPC, features = rownames(mCRPC_ARNEDNPC))
tiff(file = "mCRPC_ARNEDNPC MET with AR all downstream Heatmap-1.tiff", width = 4, height = 4, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_ARNEDNPC, features = c("MET", "AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", 'NKX3-1', "ALDH1A3", "PMEPA1", "STEAP4"), draw.lines = TRUE, size = 3, disp.max = 1.3, disp.min = -1.4) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected HGF/MET signature
tiff(file = "mCRPC_ARNEDNPC HGF selected downstream Heatmap.tiff", width = 5, height = 4, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_ARNEDNPC, features = c("MET", "AQP9","RUNX1", "SERPINE1",
                                 "ABCC3",  "SLC43A3",
                                 "GPNMB","RASGRP3"), draw.lines = TRUE, size = 3, disp.max = 1.9, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#MET with selected WNT signature
tiff(file = "mCRPC_ARNEDNPC WNT selected downstream Heatmap.tiff", width = 5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_ARNEDNPC, features = c("MET", "LEF1", "TCF7L2", 
                                 "PLAUR", "MMP7", "LGR5", "FST", "BMP2"), draw.lines = TRUE, size = 2, disp.max = 1.9, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

####DEGs####
#DEGs_DNPCvARPC
DefaultAssay(mCRPC_ARNEDNPC) <- "RNA"
all.genes <- rownames(mCRPC_ARNEDNPC)
mCRPC_ARNEDNPC <- ScaleData(mCRPC_ARNEDNPC, features = all.genes)
DNPCvARPC.Markers <- FindMarkers(mCRPC_ARNEDNPC, ident.1 = c("DNPC"), 
                                 ident.2 = c("ARPC"), min.pct = 0, logfc.threshold = 0)
write.csv(DNPCvARPC.Markers, "DNPCvARPC.Markers.csv")

#p.adjust
DEG_DNPCvARPC <- read.csv("DNPCvARPC.Markers.csv") 
DEG_DNPCvARPC_pvalue <- DEG_DNPCvARPC$p_val
DEG_DNPCvARPC_pvalue=as.numeric(DEG_DNPCvARPC_pvalue)
DEG_DNPCvARPC_BH = p.adjust(DEG_DNPCvARPC_pvalue, "BH")
write.csv(DEG_DNPCvARPC_BH, "DEG_DNPCvARPC_BH.csv")

#DEGs_NEPCvARPC
DefaultAssay(mCRPC_ARNEDNPC) <- "RNA"
all.genes <- rownames(mCRPC_ARNEDNPC)
mCRPC_ARNEDNPC <- ScaleData(mCRPC_ARNEDNPC, features = all.genes)
NEPCvARPC.Markers <- FindMarkers(mCRPC_ARNEDNPC, ident.1 = c("NEPC"), 
                                 ident.2 = c("ARPC"), min.pct = 0, logfc.threshold = 0)
write.csv(NEPCvARPC.Markers, "NEPCvARPC.Markers.csv")

#p.adjust
DEG_NEPCvARPC <- read.csv("NEPCvARPC.Markers.csv") 
DEG_NEPCvARPC_pvalue <- DEG_NEPCvARPC$p_val
DEG_NEPCvARPC_pvalue=as.numeric(DEG_NEPCvARPC_pvalue)
DEG_NEPCvARPC_BH = p.adjust(DEG_NEPCvARPC_pvalue, "BH")
write.csv(DEG_NEPCvARPC_BH, "DEG_NEPCvARPC_BH.csv")

####Correlation Analysis####

