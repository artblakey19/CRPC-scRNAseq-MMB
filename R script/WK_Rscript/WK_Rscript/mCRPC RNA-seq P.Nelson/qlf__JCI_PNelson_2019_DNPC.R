####Load required packages####
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

###
setwd("//Users/won/Desktop/Sun Lab/Data/HGF-MET-Bcat/Human/2019_JCI_PNelson")
countdata <- read.table(file = "COUNTS.csv", sep = ",", header = T, row.names=1, as.is=T)

a <- as.matrix(countdata)
#Log Transform
b=log2(a/10+1)
c=log2(a+1)
CRPC2019 =  CreateSeuratObject(counts = c, project = "2019PNelson")

ARNE_meta1 <- CRPC2019@meta.data
ARNE_meta1<- rownames_to_column(ARNE_meta1)
write_csv(ARNE_meta1, "ARNE_meta1.csv")
ARNE_meta <- column_to_rownames(ARNE_meta, var = "rowname")
CRPC2019 <- AddMetaData(CRPC2019,ARNE_meta)

#Clustering
CRPC2019 <- FindVariableFeatures(CRPC2019, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(CRPC2019)
CRPC2019 <- ScaleData(CRPC2019, features = all.genes)
CRPC2019 <- RunPCA(CRPC2019, features = VariableFeatures(object = CRPC2019))

ElbowPlot(CRPC2019, ndims = 50)

CRPC2019 <- FindNeighbors(CRPC2019, reduction = "pca", dims = 1:20)
CRPC2019 <- FindClusters(CRPC2019, resolution = 0.5)
CRPC2019 <- RunUMAP(CRPC2019, reduction = "pca", dims = 1:20)
DimPlot(CRPC2019, reduction = "umap", pt.size = 0.3) 

Idents(object = CRPC2019) <- "ARNE"
DimPlot(CRPC2019, reduction = "umap", pt.size = 0.3) 

#DEGs_
DefaultAssay(CRPC2019) <- "RNA"
all.genes <- rownames(CRPC2019)
CRPC2019 <- ScaleData(CRPC2019, features = all.genes)
DNNEPCvARPC.Markers <- FindMarkers(CRPC2019, ident.1 = c("DNPC1", "DNPC2", "NEPC"), 
                                   ident.2 = c("ARPC"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(DNNEPCvARPC.Markers, "DNNEPCvARPC.Markers.csv")











####SU2C_PCF_mCRPC####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/mCRPC RNA-seq P.Nelson/RNAseq")
countdata <- read.table(file = "ARPCNEPCDNPC_FPKM_Final.csv", sep = ",", header = T, row.names=1, as.is=T) #duplicate 'row.names' are not allowed

mCRPC <- CreateSeuratObject(counts = countdata,  min.cells = 3, min.features = 200, project = "2019_PNelson")

#Clustering
mCRPC <- FindVariableFeatures(mCRPC, selection.method = "vst", nfeatures = 2500)

mCRPC <- ScaleData(mCRPC, verbose = FALSE)
mCRPC <- RunPCA(mCRPC, npcs = 50, verbose = FALSE, approx=FALSE)
ElbowPlot(mCRPC, ndims = 50)

mCRPC <- FindNeighbors(mCRPC, reduction = "pca", dims = 1:30)
mCRPC <- FindClusters(mCRPC, resolution = 0.5)
mCRPC <- RunUMAP(mCRPC, reduction = "pca", dims = 1:30)
DimPlot(mCRPC, reduction = "umap", pt.size = 0.3) 




#AR score NE score assignment
ARScore <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
AR.genes <- ARScore[[1]]
NEScore <- list(c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
NE.genes <- NEScore[[1]]

DefaultAssay(mCRPC) <- "RNA"
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
mCRPC <- CellCycleScoring(mCRPC, s.features = AR.genes, g2m.features = NE.genes, set.ident = TRUE)

Idents(object = mCRPC) <- "Phase"
DimPlot(mCRPC, reduction = "umap")

mCRPC <- RenameIdents(object = mCRPC, 'S' = "ARPC", 'G2M' = "NEPC", 'G1' = "DNPC")
mCRPC[["ARNE"]] <- Idents(object = mCRPC)

Idents(object = mCRPC) <- "ARNE"

tiff(file = "mCRPC AR NE Score UMAP.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(mCRPC, reduction = "umap", pt.size = 1, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

tiff(file = "mCRPC MET expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(mCRPC, reduction = "umap", features = c("MET"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

tiff(file = "mCRPC MET expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(mCRPC, reduction = "umap", features = c("MET"), cols = c("light grey", "red"), pt.size = 1, max.cutoff = "q90")
dev.off()

tiff(file = "mCRPC MET Vln.tiff", width = 5, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC, features = "MET", pt.size = 0.3, cols = c("lightseagreen", "orange", "magenta2"))
dev.off()

tiff(file = "mCRPC AR WNT MET downstream Heatmap.tiff", width = 5, height = 6, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("AR", "SYP", "HGF", "MET", "CTNNB1", "MYC", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR", "TCF7L2", "XPO1"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

table(Idents(mCRPC))

#Heatmap for AR signature
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
mCRPC <- AddModuleScore(object = mCRPC, features = ARlist, name = "ARDownstream") 
mCRPC[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC, vars = 'ARDownstream1')))
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
tiff(file = "mCRPC ARscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()
tiff(file = "mCRPC ARscore Heatmap-1.tiff", width = 5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#Heatmap for NE signature
NElist <- list(c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
mCRPC <- AddModuleScore(object = mCRPC, features = NElist, name = "NEmodule") 
mCRPC[['NEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC, vars = 'NEmodule1')))
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
tiff(file = "mCRPC NEscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC, features = 'NEmodule1', assay = 'NEmodule', slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#MET with AR downstream Heatmap
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
png(file = "mCRPC MET with AR all downstream Heatmap.png", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("MET", "AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1",  "ALDH1A3", "STEAP4"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

###MET signaling
#Heatmap for all HGF/MET signature
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
tiff(file = "mCRPC HGF all downstream Heatmap-1.tiff", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("HGF",  "MET", "AQP9", "CCL16", "CXCL1",
                               "ABCC3", "IL1R2", "IRF4", "LY75", "TNFAIP6",
                               "SNTB1", "APOC1", "SPHK1", "MMP9", "RUNX1",
                               "SLC43A3", "SNX10", "RAB38", "SERPINF1", "ADAM28",
                               "SERPINE1", "AKR1B1", "SPP1", "GPNMB", "DAB2", 
                               "RASGRP3", "MMP2", "LPXN", "TNIP3", "RFTN1",
                               "ST3GAL1", "FSCN1", "SOD2"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Heatmap for selected HGF/MET signature
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
tiff(file = "mCRPC HGF selected downstream Heatmap.tiff", width = 5.5, height = 4, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("HGF", "MET", "AQP9","RUNX1", "CCL16","SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3","MMP2","RFTN1"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

tiff(file = "mCRPC HGF selected downstream Heatmap without HGF.tiff", width = 5.5, height = 4, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("MET", "AQP9","RUNX1", "CCL16","SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3","MMP2","RFTN1"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

tiff(file = "mCRPC HGF selected downstream Heatmap-1.tiff", width = 5.5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("HGF", "MET", "AQP9","RUNX1", "SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

tiff(file = "mCRPC HGF selected downstream Heatmap without HGF-1.tiff", width = 5.5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("MET", "AQP9","RUNX1", "SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Correlation_DNPC_MET with AR downstream
GOI <- c( "MET", "KLK3", "KLK2", "NKX3-1", "PMEPA1", "AR")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('DNPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR downstream spearman correlation DNPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_ARPC_MET with AR downstream
GOI <- c("MET", "KLK3", "KLK2", "NKX3-1", "PMEPA1", "AR")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR downstream spearman correlation ARPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_NEPC_MET with AR downstream
GOI <- c("MET", "KLK3", "KLK2", "NKX3-1", "PMEPA1", "AR")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('NEPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR downstream spearman correlation NEPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

###WNT signaling
#Correlation for all WNT downstream
GOI <- c( "MET", "MYC", "CCND1", "AXIN2",
          "LEF1","TCF7", "TCF7L2", "TCF4",
          "TNFAIP3", "BMP2", "LGR5", "LGR6",
          "RSPO1", "RSPO2", "RSPO3", "RSPO4", 
          "FGFR2", "TLR2", "SERPINB1", "NKD1",
          "CLDN1", "CTLA4", "DKK1", "EDN1", "FOSL1",
          "FST", "GJA1",  "MMP7", "NRCAM", "PLAU", 
          "PLAUR", "PTTG1", "VCAN")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)

#Correlation for selected WNT downstream
GOI <- c( "AR", "MET", "MYC", "CCND1", "AXIN2", "LEF1", "TCF4", "TCF7L2",
          "LGR6", "FST",  "BMP2", 
          "MMP7", "SERPINB1", "TNFAIP3")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)


#Heatmap for all WNT signature
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
tiff(file = "mCRPC MET and all WNT downstream Heatmap.tiff", width = 5.5, height = 4, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c( "MET", "HGF", "LEF1", "TCF7L2", "AXIN2", "PLAUR", "MMP7", "SNAI1", "TNFAIP3",
                               "LGR6", "FST", "BMP2"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

tiff(file = "mCRPC MET and selected WNT downstream Heatmap.tiff", width = 5.5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c( "MET", "HGF", "LEF1", "TCF7L2", "PLAUR", "MMP7",
                               "LGR6", "FST", "BMP2" ), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Correlation_DNPC_MET with Bcat downstream
GOI <- c( "HGF", "MET", "LGR6", "FST", "BMP2", "MMP7")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('DNPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with WNT downstream spearman correlation DNPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_ARPC_MET with Bcat downstream
GOI <- c( "HGF", "MET", "LGR6", "FST", "BMP2", "MMP7")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with WNT downstream spearman correlation ARPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_NEPC_MET with Bcat downstream
GOI <- c( "HGF", "MET", "LGR6", "FST", "BMP2", "MMP7")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('NEPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with WNT downstream spearman correlation NEPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_DNPC_MET with AR and Bcat downstream
GOI <- c( "HGF", "MET", "LGR6", "FST", "BMP2", "MMP7", "KLK3", "KLK2", "NKX3-1", "PMEPA1", "AR")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('NEPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR WNT downstream spearman correlation NEPC.tiff", width = 9, height = 9, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_ARPC_MET with AR and Bcat downstream
GOI <- c( "HGF", "MET", "LGR6", "FST", "BMP2", "MMP7", "KLK3", "KLK2", "NKX3-1", "PMEPA1", "AR")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR WNT downstream spearman correlation ARPC.tiff", width = 9, height = 9, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_NEPC_MET with AR and Bcat downstream
GOI <- c( "HGF", "MET", "LGR6", "FST", "BMP2", "MMP7", "KLK3", "KLK2", "NKX3-1", "PMEPA1", "AR")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('NEPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR WNT downstream spearman correlation NEPC.tiff", width = 9, height = 9, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

##MET and MET downstream
#Correlation for all MET downstream
GOI <- c( "AR", "MET", "AQP9", "RUNX1",
          "SERPINE1", "ABCC3",
          "SLC43A3", "GPNMB", "RASGRP3")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)

#Correlation_DNPC_MET with AR and Bcat downstream
GOI <- c( "AR", "MET", "RUNX1", "AQP9", "SLC43A3", "ABCC3")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('DNPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR MET downstream spearman correlation DNPC.tiff", width = 9, height = 9, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()

#Correlation_DNPC_MET with AR and Bcat downstream
GOI <- c( "AR", "MET", "RUNX1", "AQP9", "SLC43A3", "ABCC3")  
GOI_index <- is.element(rownames(mCRPC),GOI)
Cell_index <- is.element(Idents(mCRPC), c('ARPC'))
expr_GOI <- mCRPC@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- mCRPC@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR MET downstream spearman correlation ARPC.tiff", width = 9, height = 9, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()


#Metadata after ARNE scoring
write.csv(mCRPC@meta.data, file = "mCRPC@meta.data.csv")

####AR+ARPC ARlow DNPC####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/Human PCa/SU2C_PCF/AR+ARPCvAR-DNPC")

Idents(object = mCRPC) <- "ARNE"
DNPCvARPC <- subset(mCRPC, idents = c("ARPC", "DNPC"))

DefaultAssay(DNPCvARPC) <- "RNA"
tiff(file = "DNPCvARPC AR expression plots.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(DNPCvARPC, reduction = "umap", features = c("AR"), cols = c("light grey", "red"), split.by = "ARNE", pt.size = 1, max.cutoff = "q90", keep.scale = "all")
dev.off()

#Add AR expression info.
#Add hMETtg information
DefaultAssay(DNPCvARPC) <- "RNA"
DNPCvARPCARPos <- subset(x=DNPCvARPC,  subset = `AR` > 20)
DNPCvARPCARNeg <- subset(x=DNPCvARPC,  subset = `AR` < 20)
Idents(object = DNPCvARPCARPos) <- "ARPos"
Idents(object = DNPCvARPCARNeg) <- "ARNeg"
DNPCvARPCARPos[["ARExp"]] <- Idents(object = DNPCvARPCARPos)
DNPCvARPCARNeg[["ARExp"]] <- Idents(object = DNPCvARPCARNeg)
DNPCvARPCAR <- merge(x = DNPCvARPCARPos, y = DNPCvARPCARNeg)
Idents(object = DNPCvARPCAR) <- "ARExp"
DNPCvARPC$ARExp <- Idents(object = DNPCvARPCAR)
Idents(object = DNPCvARPC) <- "ARExp"
DimPlot(DNPCvARPC, reduction = "umap", pt.size = 0.3, label = TRUE)

#Counts
Idents(object = DNPCvARPC) <- "ARNE"
DNPCvARPC$ARExp.ARNE <- paste(Idents(DNPCvARPC), DNPCvARPC$ARExp, sep = "_")
Idents(object = DNPCvARPC) <- "ARExp.ARNE"
table(Idents(DNPCvARPC))

#Subset
Idents(object = DNPCvARPC) <- "ARExp.ARNE"
DNPCvARPC_ARExp <- subset(DNPCvARPC, idents = c("ARPC_ARPos", "DNPC_ARNeg"))

#Heatmap for AR signature
Idents(object = DNPCvARPC_ARExp) <- "ARNE"
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
DNPCvARPC_ARExp <- AddModuleScore(object = DNPCvARPC_ARExp, features = ARlist, name = "ARDownstream") 
DNPCvARPC_ARExp[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = DNPCvARPC_ARExp, vars = 'ARDownstream1')))
all.genes <- rownames(DNPCvARPC_ARExp)
DNPCvARPC_ARExp <- ScaleData(DNPCvARPC_ARExp, features = all.genes)
tiff(file = "DNPCvARPC_ARExp ARscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = DNPCvARPC_ARExp, features = 'ARDownstream1', assay = 'ARmodule', disp.max = 6, disp.min = -6, slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()
tiff(file = "DNPCvARPC_ARExp ARscore Heatmap-1.tiff", width = 5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = DNPCvARPC_ARExp, features = 'ARDownstream1', assay = 'ARmodule', disp.max = 6, disp.min = -6, slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#Heatmap for NE signature
NElist <- list(c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
DNPCvARPC_ARExp <- AddModuleScore(object = DNPCvARPC_ARExp, features = NElist, name = "NEmodule") 
DNPCvARPC_ARExp[['NEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = DNPCvARPC_ARExp, vars = 'NEmodule1')))
all.genes <- rownames(DNPCvARPC_ARExp)
DNPCvARPC_ARExp <- ScaleData(DNPCvARPC_ARExp, features = all.genes)
tiff(file = "DNPCvARPC_ARExp NEscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = DNPCvARPC_ARExp, features = 'NEmodule1', assay = 'NEmodule', disp.max = 10, disp.min = -1, slot = 'data', draw.lines = TRUE) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off

#MET with AR downstream Heatmap
DefaultAssay(DNPCvARPC_ARExp) <- "RNA"
DNPCvARPC_ARExp <- ScaleData(DNPCvARPC_ARExp, features = rownames(DNPCvARPC_ARExp))
tiff(file = "DNPCvARPC_ARExp MET with AR all downstream Heatmap.png", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(DNPCvARPC_ARExp, features = c("MET", "AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1",  "ALDH1A3", "STEAP4"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Correlation_ARPC_MET with AR downstream
Idents(object = DNPCvARPC_ARExp) <- "ARNE"
GOI <- c("MET", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4")  
GOI_index <- is.element(rownames(DNPCvARPC_ARExp),GOI)
Cell_index <- is.element(Idents(DNPCvARPC_ARExp), c('ARPC'))
expr_GOI <- DNPCvARPC_ARExp@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- DNPCvARPC_ARExp@assays$RNA@counts[GOI_index,Cell_index]

tiff(file = "MET with AR downstream spearman correlation ARPC.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE, label_size = 6)
dev.off()



####SU2C_PCF_mCRPC_ABITreatment####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF")
countdata1 <- read.table(file = "SU2C_PCF_CRPC_mRNAseq_polyA_1.csv", sep = ",", header = T, row.names=1, as.is=T)

mCRPC1 <- CreateSeuratObject(counts = countdata1,  min.cells = 3, min.features = 200, project = "SU2C")

#Clustering
mCRPC1 <- FindVariableFeatures(mCRPC1, selection.method = "vst", nfeatures = 5000)

mCRPC1 <- ScaleData(mCRPC1, verbose = FALSE)
mCRPC1 <- RunPCA(mCRPC1, npcs = 50, verbose = FALSE)
ElbowPlot(mCRPC1, ndims = 50)

mCRPC1 <- FindNeighbors(mCRPC1, reduction = "pca", dims = 1:20)
mCRPC1 <- FindClusters(mCRPC1, resolution = 0.5)
mCRPC1 <- RunTSNE(mCRPC1, reduction = "pca", dims = 1:20)
mCRPC1 <- RunUMAP(mCRPC1, reduction = "pca", dims = 1:20)
DimPlot(mCRPC1, reduction = "umap", pt.size = 0.3) 

#AR score NE score assignment
ARScore <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
AR.genes <- ARScore[[1]]
NEScore <- list(c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
NE.genes <- NEScore[[1]]

DefaultAssay(mCRPC1) <- "RNA"
all.genes <- rownames(mCRPC1)
mCRPC1 <- ScaleData(mCRPC1, features = all.genes)
mCRPC1 <- CellCycleScoring(mCRPC1, s.features = AR.genes, g2m.features = NE.genes, set.ident = TRUE)

Idents(object = mCRPC1) <- "Phase"
DimPlot(mCRPC1, reduction = "umap")

mCRPC1 <- RenameIdents(object = mCRPC1, 'S' = "ARPC", 'G2M' = "NEPC", 'G1' = "DNPC")
mCRPC1[["ARNE"]] <- Idents(object = mCRPC1)

Idents(object = mCRPC1) <- "ARNE"

#Add Ar information
DefaultAssay(mCRPC1) <- "RNA"
mCRPC1_ARPC <- subset(mCRPC1, idents = c("ARPC"))
DefaultAssay(mCRPC1_ARPC) <- "RNA"
mCRPC1_ARPC <- ScaleData(mCRPC1_ARPC, features = all.genes)
FeaturePlot(mCRPC1_ARPC, reduction = "umap", features = c("AR"), cols = c("light grey", "red"), pt.size = 0.3)

mCRPC1ARPos <- subset(x=mCRPC1,  subset = `AR` > 250)
mCRPC1ARNeg <- subset(x=mCRPC1,  subset = `AR` < 250)
Idents(object = mCRPC1ARPos) <- "ARPos"
Idents(object = mCRPC1ARNeg) <- "ARNeg"
mCRPC1ARPos[["ARExp"]] <- Idents(object = mCRPC1ARPos)
mCRPC1ARNeg[["ARExp"]] <- Idents(object = mCRPC1ARNeg)
mCRPC1AR <- merge(x = mCRPC1ARPos, y = mCRPC1ARNeg)
Idents(object = mCRPC1AR) <- "ARExp"
mCRPC1$ARExp <- Idents(object = mCRPC1AR)
Idents(object = mCRPC1) <- "ARExp"
DimPlot(mCRPC1, reduction = "umap", pt.size = 0.3, label = TRUE)

#Cell Counts
Idents(object = mCRPC1) <- "ARNE"
mCRPC1$ARExp.ARNE <- paste(Idents(mCRPC1), mCRPC1$ARExp, sep = "_")
Idents(object = mCRPC1) <- "ARExp.ARNE"
table(Idents(mCRPC1))

#Subset
Idents(object = mCRPC1) <- "ARExp.ARNE"
mCRPC1_1 <- subset(mCRPC1, idents = c("ARPC_ARPos", "DNPC_ARNeg", "DNPC_ARPos", "NEPC_ARNeg"))

Idents(object = mCRPC1_1) <- "ARNE"
mCRPC1_1 <- ScaleData(mCRPC1_1, features = all.genes)

DoHeatmap(mCRPC1_1, features = c("AR", "CHGA", "SYP", "HGF", "MET", "CTNNB1", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR", "TCF7L2", "XPO1"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

tiff(file = "TriplevDouble.combined1.epi2 AR WNT MET downstream Heatmap.tiff", width = 10, height = 6, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC1_1, features = c("AR", "CHGA", "SYP", "HGF", "MET", "CTNNB1", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR", "TCF7L2", "XPO1"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#subset
Idents(object = mCRPC1) <- "ARExp.ARNE"
mCRPC1_2 <- subset(mCRPC1, idents = c("ARPC_ARPos", "DNPC_ARNeg", "NEPC_ARNeg"))

Idents(object = mCRPC1_2) <- "ARNE"
mCRPC1_2 <- ScaleData(mCRPC1_2, features = all.genes)

DoHeatmap(mCRPC1_2, features = c("AR", "CHGA", "SYP", "HGF", "MET", "CTNNB1", "MYC", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR", "TCF7L2", "XPO1"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))

tiff(file = "TriplevDouble.combined1.epi2 AR WNT MET downstream Heatmap.tiff", width = 5, height = 6, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC1_2, features = c("AR", "CHGA", "SYP", "HGF", "MET", "CTNNB1", "MYC", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR", "TCF7L2", "XPO1"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

table(Idents(mCRPC1_2))

####SU2C_PCF_mCRPC_original ABI metadata####
meta.data = read.csv("Treatment.csv")
mCRPC[["Treatment"]] <- meta.data$Treatment[match(rownames(mCRPC@meta.data), meta.data$X)]

#Cell Counts
Idents(object = mCRPC) <- "ARNE"
mCRPC$Treatment.ARNE <- paste(Idents(mCRPC), mCRPC$Treatment, sep = "_")
Idents(object = mCRPC) <- "Treatment.ARNE"
table(Idents(mCRPC))

mCRPC_ABI <- subset(mCRPC, idents = c("ARPC_Naive", "ARPC_Exposed",
                                     "NEPC_Naive", "NEPC_Exposed",
                                     "DNPC_Naive", "DNPC_Exposed"))
DimPlot(mCRPC_ABI, reduction = "umap")

mCRPC_ABI <- RenameIdents(object = mCRPC_ABI, 'ARPC_Naive' = "ARPC", 'ARPC_Exposed' = "ARPC", 
                          'DNPC_Naive' = "DNPC", 'DNPC_Exposed' = "DNPC", 
                          'NEPC_Naive' = "NEPC", 'NEPC_Exposed' = "NEPC")


mCRPC_ABI <- RenameIdents(object = mCRPC_ABI, 'ARPC_Naive' = "ARPC_Naive", 'ARPC_Exposed' = "ARPC_Exposed", 
                      'DNPC_Naive' = "DNPC_Naive", 'DNPC_Exposed' = "DNPC_Exposed", 
                      'NEPC_Naive' = "NEPC_Naive", 'NEPC_Exposed' = "NEPC_Exposed")
mCRPC_ABI[["Treatment.ARNE"]] <- Idents(object = mCRPC_ABI)

table(Idents(mCRPC_ABI))

tiff(file = "mCRPC ABI AR WNT MET Ribosome G2M downstream Heatmap.tiff", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_ABI, features = c("AR", "CHGA", "SYP", "MKI67", "PCNA", "CENPF", "CENPE", "TOP1", "NASP", 
                                  "HGF", "MET", "RAP1B", "ITGB1", "UBA52",
                                  "CTNNB1", "MYC", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR",  "TCF7L2", 
                                  "XPO1", "RPL12", "RPS16", "RPL3", "RPS5", "RPS20"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Exposed
mCRPC_ABI_Exposed <- subset(mCRPC, idents = c("ARPC_Exposed",
                                     "NEPC_Exposed",
                                     "DNPC_Exposed"))
DefaultAssay(mCRPC_ABI_Exposed) <- "RNA"
all.genes <- rownames(mCRPC_ABI_Exposed)
mCRPC_ABI_Exposed <- ScaleData(mCRPC_ABI_Exposed, features = all.genes)
tiff(file = "mCRPC_ABI_Exposed AR WNT MET Ribosome G2M downstream Heatmap.tiff", width = 4.5, height = 6, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_ABI_Exposed, features = c("AR", "CHGA", "SYP", "MKI67", "PCNA", "CENPF", "CENPE", "TOP1", "NASP", 
                                  "HGF", "MET", "RAP1B", "ITGB1", "UBA52",
                                  "CTNNB1", "MYC", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR",  "TCF7L2", 
                                  "XPO1", "RPL12", "RPS16", "RPL3", "RPS5", "RPS20"), draw.lines = TRUE) + scale_fill_gradientn(colors = c("purple", "black", "yellow"))
dev.off()

#Correlation
GOI <- c('AR','HGF','MET', "CTNNB1", "MYC", "AXIN2", "TCF4", "LEF1", "ZEB2", "FST", "PLAUR",  "TCF7L2", 'XPO1', "RPL12", 'RPL16')  
GOI_index <- is.element(rownames(E1831),GOI)
Cell_index <- is.element(Idents(E1831),c('1'))
expr_GOI <- E1831@assays$RNA@data[GOI_index,Cell_index] 
expr_GOI <- E1831@assays$RNA@counts[GOI_index,Cell_index]
ggcorr(t(expr_GOI), method = c("pairwise", "spearman"), label = TRUE)


##DEG Analysis
#Limma
log10_fpkm = log10(fpkm_table+1)
log10_fpkmData = data.frame(log10_fpkm)
lmfit <- lmFit(logCPM, design)
lmfit <- eBayes(lmfit, trend=TRUE)
topTable(lmfit, coef = 2)
library(MAST)
mCRPC_ABI[["compare.ARNE"]] <- Idents(object = mCRPC_ABI)
ARPC_vs_DNPC_markers <- FindMarkers(mCRPC_ABI, ident.1 = "ARPC", ident.2 = "DNPC",
                                    test.use = "MAST",min.pct = 0.1, only.pos = F)

###RV-Script_FPKM_to_DEG####
#Create log2+1 normalized fpkm for ARPC/DNPC

library(scran)
library(tidyverse)
#_________________________________________________
#Subset
Idents(object = mCRPC) <- "Treatment.ARNE"
mARDNPC_ABI_od <- subset(mCRPC_ABI, idents = c("ARPC", "DNPC"))
table(Idents(mCRPC_ABI))
##extract count matrix
mARDNPC_count <- mARDNPC_ABI@assays$RNA@counts
a <- as.matrix(mARDNPC_count)
#Log Transform
b=log2(a/10+1)
c=log2(a+1)
#Check Gen Expression profile
b[1:5,1:5]
b['MET', ]
c['MET', ]
mARDNPC_count_log2 <- as.data.frame(c)
#Metadata
head(mARDNPC_ABI@meta.data)
mARDNPC_ABI_meta<- data.frame(mARDNPC_ABI_od@meta.data)
mARDNPC_ABI_meta<- rownames_to_column(mARDNPC_ABI_meta)
mARDNPC_ABI_meta <- select(mARDNPC_ABI_meta,rowname,orig.ident,Treatment,ARNE,
                           Treatment.ARNE,compare.ARNE,S.Score,G2M.Score, Phase)
mARDNPC_ABI_meta <- column_to_rownames(mARDNPC_ABI_meta, var = "rowname")
#create seurat object
mARDNPC_ABI =  CreateSeuratObject(counts = mARDNPC_count_log2, project = "mARDNPC")
#add metadata
mARDNPC_ABI <- AddMetaData(mARDNPC_ABI,mARDNPC_ABI_meta)
Idents(object = mARDNPC_ABI) <- "Treatment.ARNE"
table(Idents(mARDNPC_ABI))
mARDNPC_ABI <- FindVariableFeatures(mARDNPC_ABI, selection.method = "vst", nfeatures = 6000)
#Intigrate#Fast MNN..........
library(SeuratWrappers)
#Cluster
mARDNPC_ABI <- FindVariableFeatures(mARDNPC_ABI, selection.method = "vst", nfeatures = 6000)
all.genes <- rownames(mARDNPC_ABI)
mARDNPC_ABI <- ScaleData(mARDNPC_ABI, features = all.genes)
mARDNPC_ABI<- RunPCA(mARDNPC_ABI, features = VariableFeatures(object = mARDNPC_ABI))
ElbowPlot(mARDNPC_ABI)
mARDNPC_ABI <- RunUMAP(mARDNPC_ABI, dims = 1:20)
mARDNPC_ABI <- FindNeighbors(mARDNPC_ABI, dims = 1:30)
mARDNPC_ABI <- FindClusters(mARDNPC_ABI)
Idents(object = mARDNPC_ABI) <- "compare.ARNE"
DimPlot(mARDNPC_ABI, reduction = "umap", pt.size = 1, label = TRUE)
#DEGs
library(MAST)
ARPC_vs_DNPC_markers <- FindAllMarkers(mARDNPC_ABI,test.use = "wilcox", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.table(ARPC_vs_DNPC_markers, "ARPC_vs_DNPC_markers.txt")

##DEGs DNPC_Exposed vs ARPC_Naive
Idents(object = mARDNPC_ABI) <- "Treatment.ARNE"
DimPlot(mARDNPC_ABI, reduction = "umap", pt.size = 1, label = TRUE)
DNPC_Exposed_vs_ARPC_Naive_markers <- FindMarkers(mARDNPC_ABI, test.use = "wilcox", ident.1 = "DNPC_Exposed", ident.2 = 
                                                    "ARPC_Naive", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.table(DNPC_Exposed_vs_ARPC_Naive_markers, "DNPC_Exposed_vs_ARPC_Naive_markers.txt")

##DEGs DNPC_Exposed vs ARPC_Exposed
DNPC_Exposed_vs_ARPC_Exposed_markers <- FindMarkers(mARDNPC_ABI, test.use = "wilcox", ident.1 = "DNPC_Exposed", ident.2 = 
                                                    "ARPC_Exposed", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.table(DNPC_Exposed_vs_ARPC_Exposed_markers, "DNPC_Exposed_vs_ARPC_Exposed_markers.txt")

##DEGs DNPC vs ARPC
Idents(object = mARDNPC_ABI) <- "ARNE"
DimPlot(mARDNPC_ABI, reduction = "umap", pt.size = 1, label = TRUE)
DNPC_vs_ARPC_markers <- FindMarkers(mARDNPC_ABI, test.use = "wilcox", ident.1 = "DNPC", ident.2 = 
                                                    "ARPC", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.table(DNPC_vs_ARPC_markers, "DNPC_vs_ARPC_markers.txt")

####mCRPC_Treatment####

#ARNE
#Heatmap
Idents(object = mCRPC) <- "ARNE"
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
tiff(file = "mCRPC ARNE XPO RP family.tiff", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                                        "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3","GPNMB","RASGRP3",
                                        "LEF1", "TCF7L2", "PLAUR", "MMP7","LGR6", "FST", "BMP2", 
                                        "XPO1", "XPO4", "XPO5", "XPO6", "XPO7",
                                        "RPL12", "RPL39L", 
                                        "RPS16", "RPS20", "RPS29", 
                              "EIF4A1", "EIF4E2", "EIF4EBP1", "EIF5A2", "EIF1B", "EIF5B"
                              ), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Heatmap
Idents(object = mCRPC) <- "ARNE"
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
tiff(file = "mCRPC ARNE XPO RP family selected.tiff", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                              "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3","GPNMB","RASGRP3",
                              "LEF1", "TCF7L2", "PLAUR", "MMP7","LGR6", "FST", "BMP2", 
                              "XPO1", "XPO7",
                             "RPL39L", 
                              "RPS16", "RPS20", "RPS29", 
                               "EIF4EBP1", "EIF5A2", "EIF1B", "EIF5B"
), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Treatment.ARNE
Idents(object = mCRPC) <- "ARNE"
mCRPC$Treatment.ARNE <- paste(Idents(mCRPC), mCRPC$Treatment, sep = "_")
Idents(object = mCRPC) <- "Treatment.ARNE"
table(Idents(mCRPC))

Idents(object = mCRPC) <- "Treatment.ARNE"
mCRPC <- RenameIdents(object = mCRPC, 
                             'ARPC_Naive'="ARPC_Naive", 'ARPC_Exposed'="ARPC_Exposed",'ARPC_Unknown'="ARPC_Unknown",
                      'NEPC_Naive'="NEPC_Naive", 'NEPC_Exposed'="NEPC_Exposed",'NEPC_Unknown'="NEPC_Unknown",
                      'DNPC_Naive'="DNPC_Naive", 'DNPC_Exposed'="DNPC_Exposed",'DNPC_Unknown'="DNPC_Unknown")  
mCRPC[["Treatment.ARNE"]] <- Idents(object = mCRPC)

#Heatmap
Idents(object = mCRPC) <- "Treatment.ARNE"
DefaultAssay(mCRPC) <- "RNA"
mCRPC <- ScaleData(mCRPC, features = rownames(mCRPC))
tiff(file = "mCRPC Treatment.ARNE XPO RP family selected.tiff", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                              "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3","GPNMB","RASGRP3",
                              "LEF1", "TCF7L2", "PLAUR", "MMP7","LGR6", "FST", "BMP2", 
                              "XPO1", "XPO7",
                              "RPL39L", 
                              "RPS16", "RPS20", "RPS29", 
                              "EIF4EBP1", "EIF5A2", "EIF1B", "EIF5B"
), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#DNPC_ExposedvARPC_all
DefaultAssay(mCRPC) <- "RNA"
Idents(object = mCRPC) <- "Treatment.ARNE"
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
DNPC_ExposedvARPC.Markers <- FindMarkers(mCRPC, ident.1 = c("DNPC_Exposed"), 
                                         ident.2 = c("ARPC_Naive", "ARPC_Unknown", "ARPC_Exposed"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(DNPC_ExposedvARPC.Markers, "DNPC_ExposedvARPC.Markers.csv")

#Heatmap
Idents(object = mCRPC) <- "Treatment.ARNE"
mCRPC_Treatment <- subset(mCRPC, idents = c("ARPC_Naive", "ARPC_Unknown", "ARPC_Exposed", "DNPC_Exposed"))

Idents(object = mCRPC_Treatment) <- "ARNE"
DefaultAssay(mCRPC_Treatment) <- "RNA"
mCRPC_Treatment <- ScaleData(mCRPC_Treatment, features = rownames(mCRPC_Treatment))
tiff(file = "mCRPC_Treatment XPO RP family.tiff", width = 6, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_Treatment, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                                        "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3","GPNMB","RASGRP3",
                                        "LEF1", "TCF7L2", "PLAUR", "MMP7","LGR6", "FST", "BMP2", 
                                        "XPO1", "XPO4", "XPO5", "XPO6", "XPO7",
                                        "RPL12", "RPL39L", "RPL11", 
                                        "RPS16", "RPS20", "RPS29", "RPS24", "RPS27L"
), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#DNPC_ExposedvARPC_Naive
DefaultAssay(mCRPC) <- "RNA"
Idents(object = mCRPC) <- "Treatment.ARNE"
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
DNPC_ExposedvARPC_Naive.Markers <- FindMarkers(mCRPC, ident.1 = c("DNPC_Exposed"), 
                                               ident.2 = c("ARPC_Naive"), min.pct = 0.01, logfc.threshold = 0.01)
write.csv(DNPC_ExposedvARPC_Naive.Markers, "DNPC_ExposedvARPC_Naive.Markers.csv")

Idents(object = mCRPC) <- "Treatment.ARNE"
mCRPC_Treatment1 <- subset(mCRPC, idents = c("ARPC_Naive", "DNPC_Exposed"))

Idents(object = mCRPC_Treatment1) <- "Treatment.ARNE"
DefaultAssay(mCRPC_Treatment1) <- "RNA"
mCRPC_Treatment1 <- ScaleData(mCRPC_Treatment1, features = rownames(mCRPC_Treatment1))
tiff(file = "mCRPC_Treatment1 Treatment.ARNE XPO RP family.tiff", width = 4.5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_Treatment1, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                                         "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3","GPNMB","RASGRP3",
                                         "LEF1", "TCF7L2", "PLAUR", "MMP7","LGR6", "FST", "BMP2", 
                                         "XPO1", "XPO4", "XPO5", "XPO6", "XPO7",
                                         "RPL12", "RPL39L",  
                                         "RPS16", "RPS27L", "RPS29", "RPS20", "RPS24",
                                         "EIF4E3", "EIF4EBP2", "EIF1B", "EIF4EBP1"
), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

Idents(object = mCRPC_Treatment1) <- "Treatment.ARNE"
DefaultAssay(mCRPC_Treatment1) <- "RNA"
mCRPC_Treatment1 <- ScaleData(mCRPC_Treatment1, features = rownames(mCRPC_Treatment1))
tiff(file = "mCRPC_Treatment1 Treatment.ARNE XPO RP family selected.tiff", width = 4.5, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_Treatment1, features = c("AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "ALDH1A3", "STEAP4",
                                         "MET", "AQP9","RUNX1", "SERPINE1","ABCC3","SLC43A3","GPNMB","RASGRP3",
                                         "LEF1", "TCF7L2", "PLAUR", "MMP7","LGR6", "FST", "BMP2", 
                                         "XPO1", "XPO7",
                                         "RPL39L",  
                                         "RPS27L", "RPS29", "RPS20", "RPS24",
                                         "EIF4E3", "EIF4EBP2", "EIF1B"
), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -1.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()