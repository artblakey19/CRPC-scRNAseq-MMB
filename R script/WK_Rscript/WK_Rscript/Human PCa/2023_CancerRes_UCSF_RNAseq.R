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
library(dplyr)
library(tidyverse)

####load data####
setwd("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/2023_CancerRes_UCSF")

#AR score NE score assignment
ARScore <- list(c("AR", "NKX3-1", "KLK3", "CHRNA2", "SLC45A3", "CD3G", "NAP1L2"))
AR.gene.list <- ARScore[[1]]
NEScore <- list(c("CHGA", "SYP", "ACTL6B", "SNAP25", "INSM1", "ASCL1", "CHRNB2", "SRRM4", "CELF3", "PCSK1", "SOX2", "POU3F2", "LMO3", "NKX2-1"))
NEPC.gene.list <- NEScore[[1]]

#RNAseq data
RNA.TPM.df <- read.table(file = "Lundberg_CR_2023_TPM_210_samples.csv", sep = ",", header = T, row.names=1, as.is=T)

a <- as.matrix(RNA.TPM.df)
#Log Transform
b=log2(a/10+1)
c=log2(a+1)
mCRPC_log =  CreateSeuratObject(counts = c, project = "2023UCSF")

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

DefaultAssay(mCRPC_log) <- "RNA"
all.genes <- rownames(mCRPC_log)
mCRPC_log <- ScaleData(mCRPC_log, features = all.genes)
mCRPC_log <- CellCycleScoring(mCRPC_log, s.features = AR.gene.list, g2m.features = NEPC.gene.list, set.ident = TRUE)

ARNE_meta_log1 <- mCRPC_log@meta.data
ARNE_meta_log1<- rownames_to_column(ARNE_meta_log1)
write_csv(ARNE_meta_log1, "ARNE_meta_log1.csv")
ARNE_meta_log <- read_csv("//isi-dcnl/user_data/zjsun/group/DO NOT USE-Lab Members/Won Kyung Kim/Human PCa/2023_CancerRes_UCSF/ARNE_meta_log.csv")
ARNE_meta_log <- column_to_rownames(ARNE_meta_log, var = "rowname")
mCRPC_log <- AddMetaData(mCRPC_log,ARNE_meta_log)
#Rename CellTypes
Idents(object = mCRPC_log) <- "ARNE"
DimPlot(mCRPC_log, reduction = "umap", pt.size = 0.3, label=TRUE)
mCRPC_log <- RenameIdents(object = mCRPC_log, 
                                        'ARPC'="ARPC",'ALNN'="ALNN",'APNP'="APNP",
                                        'NEPC'="NEPC", 'DNPC'="DNPC")  
mCRPC_log[["ARNE"]] <- Idents(object = mCRPC_log)
DimPlot(mCRPC_log, reduction = "umap", pt.size = 1, label=TRUE)

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Human PCa/2023_CancerRes_UCSF/NEW")

#Heatmap for AR signature
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
mCRPC_log <- AddModuleScore(object = mCRPC_log, features = ARlist, name = "ARDownstream") 
mCRPC_log[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC_log, vars = 'ARDownstream1')))
all.genes <- rownames(mCRPC_log)
mCRPC_log <- ScaleData(mCRPC_log, features = all.genes)
tiff(file = "mCRPC_log ARscore Heatmap.tiff", width = 7, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_log, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE, disp.max = 5, disp.min = -4) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()
tiff(file = "mCRPC_log ARscore Heatmap-1.tiff", width = 7, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_log, features = 'ARDownstream1', assay = 'ARmodule', slot = 'data', draw.lines = TRUE, disp.max = 5, disp.min = -4) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#Heatmap for NE signature
NElist <- list(c("SYP", "CHGA", "CHGB", "ENO2", "CHRNB2", "SCG3", "SCN3A", "PCSK1", "ELAVL4", "NKX2-1"))
mCRPC_log <- AddModuleScore(object = mCRPC_log, features = NElist, name = "NEmodule") 
mCRPC_log[['NEmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC_log, vars = 'NEmodule1')))
all.genes <- rownames(mCRPC_log)
mCRPC_log <- ScaleData(mCRPC_log, features = all.genes)
tiff(file = "mCRPC_log NEscore Heatmap.tiff", width = 7, height = 1.5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_log, features = 'NEmodule1', assay = 'NEmodule', slot = 'data', draw.lines = TRUE, disp.max = 4, disp.min = -4) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()
tiff(file = "mCRPC_log NEscore Heatmap-1.tiff", width = 7, height = 5, units = "in", compression = "lzw", res =200)
DoHeatmap(object = mCRPC_log, features = 'NEmodule1', assay = 'NEmodule', slot = 'data', draw.lines = TRUE, disp.max = 4, disp.min = -4) + scale_fill_gradientn(colors = c("darkblue", "white", "darkred"))
dev.off()

#Heatmap for selected HGF/MET signature
DefaultAssay(mCRPC_log) <- "RNA"
mCRPC_log <- ScaleData(mCRPC_log, features = rownames(mCRPC_log))
tiff(file = "mCRPC_log HGF selected downstream Heatmap.tiff", width = 5.5, height = 3.2, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_log, features = c("MET", "AQP9","RUNX1", "SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3"), draw.lines = TRUE, size = 3, disp.max = 3.2, disp.min = -2.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

tiff(file = "mCRPC_log MET and WNT downstream Heatmap.tiff", width = 5.5, height = 3.2, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_log, features = c( "MET", "LEF1", "TCF4", "PLAUR", "MMP7",
                               "LGR6", "FST", "BMP2" ), draw.lines = TRUE, size = 3, disp.max = 3.2, disp.min = -2.5) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()


tiff(file = "mCRPC_log MET and ARscore Heatmap.tiff", width = 5.5, height = 3.2, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC_log, features = c("MET", "AR", "KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "ALDH1A3", "PMEPA1", "STEAP4"), draw.lines = TRUE, size = 3, disp.max = 1.5, disp.min = -2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Vlnplots for AR & AR score



#Vlnplots

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

DefaultAssay(mCRPC_log) <- "RNA"
Idents(object = mCRPC_log) <- "ARNE"
tiff(file = "mCRPC_log Xpo1 Vln.tiff", width = 5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "XPO1", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log AR Vln.tiff", width = 5, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "AR", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log MET Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "AQP9", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

tiff(file = "mCRPC_log RPL3L Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "RPL3L", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log RPL39L Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "RPL39L", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log RPS9 Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "RPS9", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log MRPL1 Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "MRPL1", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log MRPL52 Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "MRPL52", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()
tiff(file = "mCRPC_log MRPS5 Vln.tiff", width = 5, height = 4.5, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "MRPS5", pt.size = 0) + NoLegend() +
  stat_summary(fun.data=data_summary) + stat_summary(fun = mean, geom='point', size = 25, colour = "green2", shape = 95)
dev.off()

#AR signature vlnplot
ARlist <- list(c("AR", "NKX3-1", "KLK3", "CHRNA2", "SLC45A3", "CD3G", "NAP1L2"))
mCRPC_log <- AddModuleScore(object = mCRPC_log, features = ARlist, name = "ARscore")
Idents(object = mCRPC_log) <- "ARNE"
VlnPlot(mCRPC_log, features = "ARscore1", pt.size = 0.01) 

tiff(file = "mCRPC_log AR signature Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "ARscore1", pt.size = 0.01) + NoLegend()
dev.off()

#NE signature vlnplot
NElist <- list(c("CHGA", "SYP", "ACTL6B", "SNAP25", "INSM1", "ASCL1", "CHRNB2", "SRRM4", "CELF3", "PCSK1", "SOX2", "POU3F2", "LMO3", "NKX2-1"))
mCRPC_log <- AddModuleScore(object = mCRPC_log, features = NElist, name = "NEscore")
Idents(object = mCRPC_log) <- "ARNE"
VlnPlot(mCRPC_log, features = "NEscore1", pt.size = 0.01)

tiff(file = "mCRPC_log NE signature Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "NEscore1", pt.size = 0.01) + NoLegend()
dev.off()

#HGF MET signature vlnplot
HGFlist <- list(c("MET", "AQP9","RUNX1", "SERPINE1",
                  "ABCC3",  "SLC43A3",
                  "GPNMB","RASGRP3"))
mCRPC_log <- AddModuleScore(object = mCRPC_log, features = HGFlist, name = "HGFscore")
Idents(object = mCRPC_log) <- "ARNE"
VlnPlot(mCRPC_log, features = "HGFscore1", pt.size = 0.01)

tiff(file = "mCRPC_log HGF signature Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "HGFscore1", pt.size = 0.01) + NoLegend()
dev.off()

#WNT signature vlnplot
WNTlist <- list(c("LEF1", "TCF4", "PLAUR", "MMP7",
                  "LGR6", "FST", "BMP2"))
mCRPC_log <- AddModuleScore(object = mCRPC_log, features = WNTlist, name = "WNTscore")
Idents(object = mCRPC_log) <- "ARNE"
VlnPlot(mCRPC_log, features = "WNTscore1", pt.size = 0.01)

tiff(file = "mCRPC_log WNT signature Vln.tiff", width = 6, height = 3, units = "in", compression = "lzw", res = 800)
VlnPlot(mCRPC_log, features = "WNTscore1", pt.size = 0.01) + NoLegend()
dev.off()

mCRPC_log.celltype <- mCRPC_log@meta.data$ARNE
write.csv(mCRPC_log.celltype, file = "mCRPC_log.celltype.csv")

mCRPC_log.ARscore_meta <- mCRPC_log@meta.data$ARscore1
write.csv(mCRPC_log.ARscore_meta, file = "mCRPC_log.ARscore_meta.csv")
mCRPC_log.NEscore_meta <- mCRPC_log@meta.data$NEscore1
write.csv(mCRPC_log.NEscore_meta, file = "mCRPC_log.NEscore_meta.csv")
mCRPC_log.HGFscore_meta <- mCRPC_log@meta.data$HGFscore1
write.csv(mCRPC_log.HGFscore_meta, file = "mCRPC_log.HGFscore_meta.csv")
mCRPC_log.WNTscore_meta <- mCRPC_log@meta.data$WNTscore1
write.csv(mCRPC_log.WNTscore_meta, file = "mCRPC_log.WNTscore_meta.csv")

####Wilcoxon test####

library(dplyr)

my_data <- read.csv(file.choose())
print(my_data)

group_by(my_data, Celltype) %>%
  summarise(
    count = n(),
    median = median(WNT, na.rm = TRUE),
    IQR = IQR(WNT, na.rm = TRUE)
  )

res <- wilcox.test(WNT ~ Celltype, data = my_data,
                   exact = FALSE)
res


#KStest
x <- c(1.200428215,
       0.189479643,
       0.458719555,
       1.94483904,
       0.589226421,
       1.564496762,
       -0.13651079,
       1.384169534,
       -0.03064624,
       -0.061832951,
       0.826284361,
       1.004069852,
       0.148243117
       
)

y <- c(0.578094301,
       1.05147356,
       -0.292621735,
       0.187384679,
       0.241355269,
       -0.715938361,
       0.156688268
)

  ks.test(x, y, alternative = c("two.sided", "less", "greater"),
          exact = NULL)
  
  
####DEGs####
##DEGs DNPC vs ARPC
Idents(object = mCRPC_log) <- "ARNE"
DimPlot(mCRPC_log, reduction = "umap", pt.size = 1, label = TRUE)
DNPC_vs_ARPC_markers <- FindMarkers(mCRPC_log, test.use = "wilcox", ident.1 = "DNPC", ident.2 = 
                                      "ARPC", min.pct = 0.0, logfc.threshold = 0.0, only.pos = F)
write.csv(DNPC_vs_ARPC_markers, "DNPC_vs_ARPC_Markers_NEW.csv")

#p.adjust
DEG_DNPCvARPC <- read.csv("DNPC_vs_ARPC_Markers_NEW.csv") 
DEG_DNPCvARPC_pvalue <- DEG_DNPCvARPC$p_val
DEG_DNPCvARPC_pvalue=as.numeric(DEG_DNPCvARPC_pvalue)
DEG_DNPCvARPC_BH = p.adjust(DEG_DNPCvARPC_pvalue, "BH")
write.csv(DEG_DNPCvARPC_BH, "DEG_DNPCvARPC_BH.csv")

##DEGs NEPC vs ARPC
Idents(object = mCRPC_log) <- "ARNE"
DimPlot(mCRPC_log, reduction = "umap", pt.size = 1, label = TRUE)
NEPC_vs_ARPC_markers <- FindMarkers(mCRPC_log, test.use = "wilcox", ident.1 = "NEPC", ident.2 = 
                                      "ARPC", min.pct = 0.0, logfc.threshold = 0.0, only.pos = F)
write.csv(NEPC_vs_ARPC_markers, "NEPC_vs_ARPC_Markers_NEW.csv")

#p.adjust
DEG_NEPCvARPC <- read.csv("NEPC_vs_ARPC_Markers_NEW.csv") 
DEG_NEPCvARPC_pvalue <- DEG_NEPCvARPC$p_val
DEG_NEPCvARPC_pvalue=as.numeric(DEG_NEPCvARPC_pvalue)
DEG_NEPCvARPC_BH = p.adjust(DEG_NEPCvARPC_pvalue, "BH")
write.csv(DEG_NEPCvARPC_BH, "DEG_NEPCvARPC_BH.csv")

##DEGs DNPC vs NEPC
Idents(object = mCRPC_log) <- "ARNE"
DimPlot(mCRPC_log, reduction = "umap", pt.size = 1, label = TRUE)
DNPC_vs_NEPC_markers <- FindMarkers(mCRPC_log, test.use = "wilcox", ident.1 = "DNPC", ident.2 = 
                                      "NEPC", min.pct = 0.0, logfc.threshold = 0.0, only.pos = F)
write.csv(DNPC_vs_NEPC_markers, "DNPC_vs_NEPC_Markers_NEW.csv")

#p.adjust
DEG_DNPCvNEPC <- read.csv("DNPC_vs_NEPC_Markers_NEW.csv") 
DEG_DNPCvNEPC_pvalue <- DEG_DNPCvNEPC$p_val
DEG_DNPCvNEPC_pvalue=as.numeric(DEG_DNPCvNEPC_pvalue)
DEG_DNPCvNEPC_BH = p.adjust(DEG_DNPCvNEPC_pvalue, "BH")
write.csv(DEG_DNPCvNEPC_BH, "DEG_DNPCvNEPC_BH.csv")

##DEGs ALNN vs ARPC
Idents(object = mCRPC_log) <- "ARNE"
DimPlot(mCRPC_log, reduction = "umap", pt.size = 1, label = TRUE)
ALNN_vs_ARPC_markers <- FindMarkers(mCRPC_log, test.use = "wilcox", ident.1 = "ALNN", ident.2 = 
                                      "ARPC", min.pct = 0.0, logfc.threshold = 0.0, only.pos = F)
write.csv(ALNN_vs_ARPC_markers, "ALNN_vs_ARPC_markers_NEW.csv")

#p.adjust
DEG_ALNNvARPC <- read.csv("ALNN_vs_ARPC_markers_NEW.csv") 
DEG_ALNNvARPC_pvalue <- DEG_ALNNvARPC$p_val
DEG_ALNNvARPC_pvalue=as.numeric(DEG_ALNNvARPC_pvalue)
DEG_ALNNvARPC_BH = p.adjust(DEG_ALNNvARPC_pvalue, "BH")
write.csv(DEG_ALNNvARPC_BH, "DEG_ALNNvARPC_BH.csv")

##DEGs APNP vs ARPC
Idents(object = mCRPC_log) <- "ARNE"
APNP_vs_ARPC_markers <- FindMarkers(mCRPC_log, test.use = "wilcox", ident.1 = "APNP", ident.2 = 
                                      "ARPC", min.pct = 0.0, logfc.threshold = 0.0, only.pos = F)
write.csv(APNP_vs_ARPC_markers, "APNP_vs_ARPC_markers_NEW.csv")

#p.adjust
DEG_APNPvARPC <- read.csv("APNP_vs_ARPC_markers_NEW.csv") 
DEG_APNPvARPC_pvalue <- DEG_APNPvARPC$p_val
DEG_APNPvARPC_pvalue=as.numeric(DEG_APNPvARPC_pvalue)
DEG_APNPvARPC_BH = p.adjust(DEG_APNPvARPC_pvalue, "BH")
write.csv(DEG_APNPvARPC_BH, "DEG_APNPvARPC_BH.csv")

####DEGs-wolog####

#RNAseq data
RNA.TPM.df <- read.table(file = "Lundberg_CR_2023_TPM_210_samples.csv", sep = ",", header = T, row.names=1, as.is=T)

mCRPC =  CreateSeuratObject(counts = RNA.TPM.df, project = "2023UCSF")

#Clustering
mCRPC <- FindVariableFeatures(mCRPC, selection.method = "vst", nfeatures = 5000)
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
mCRPC <- RunPCA(mCRPC, features = VariableFeatures(object = mCRPC))

ElbowPlot(mCRPC, ndims = 50)

mCRPC <- FindNeighbors(mCRPC, reduction = "pca", dims = 1:30)
mCRPC <- FindClusters(mCRPC, resolution = 0.5)
mCRPC <- RunUMAP(mCRPC, reduction = "pca", dims = 1:30)
DimPlot(mCRPC, reduction = "umap", pt.size = 0.3) 

DefaultAssay(mCRPC) <- "RNA"
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
mCRPC <- CellCycleScoring(mCRPC, s.features = AR.gene.list, g2m.features = NEPC.gene.list, set.ident = TRUE)

ARNE_meta1 <- mCRPC@meta.data
ARNE_meta1<- rownames_to_column(ARNE_meta1)
write_csv(ARNE_meta1, "ARNE_meta1.csv")
ARNE_meta <- read_csv("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Human PCa/2023_CancerRes_UCSF/ARNE_meta.csv")
ARNE_meta <- column_to_rownames(ARNE_meta, var = "rowname")
mCRPC <- AddMetaData(mCRPC,ARNE_meta)
#Rename CellTypes
Idents(object = mCRPC) <- "ARNE"
DimPlot(mCRPC, reduction = "umap", pt.size = 0.3, label=TRUE)
mCRPC <- RenameIdents(object = mCRPC, 
                          'ARPC'="ARPC",'ALNN'="ALNN",'APNP'="APNP",
                          'NEPC'="NEPC", 'DNPC'="DNPC")  
mCRPC[["ARNE"]] <- Idents(object = mCRPC)
DimPlot(mCRPC, reduction = "umap", pt.size = 1, label=TRUE)

#AR signature vlnplot
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
mCRPC <- AddModuleScore(object = mCRPC, features = ARlist, name = "ARscore")
Idents(object = mCRPC) <- "ARNE"
VlnPlot(mCRPC, features = "ARscore1", pt.size = 0.01) 

tiff(file = "GSM4089151_P1_2.epi2 AR signature Vln.tiff", width = 6, height = 4, units = "in", compression = "lzw", res = 800)
VlnPlot(GSM4089151_P1_2.epi2, features = "ARscore1", pt.size = 0.01, cols = c("steelblue1", "blue", "bisque3", "blueviolet", "aquamarine3","deeppink2", "plum4", "darkorange1", "black", "darkkhaki", "darkgoldenrod1", "darkolivegreen3", "red")) + NoLegend()
dev.off()

####DEGs DNPC vs ARPC####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Human PCa/2023_CancerRes_UCSF/mCRPC")

Idents(object = mCRPC) <- "ARNE"
DimPlot(mCRPC, reduction = "umap", pt.size = 1, label = TRUE)
DNPC_vs_ARPC_markers <- FindMarkers(mCRPC, test.use = "wilcox", ident.1 = "DNPC", ident.2 = 
                                      "ARPC", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.csv(DNPC_vs_ARPC_markers, "DNPC_vs_ARPC_Markers_NEW.csv")

#p.adjust
DEG_DNPCvARPC <- read.csv("DNPC_vs_ARPC_Markers_NEW.csv") 
DEG_DNPCvARPC_pvalue <- DEG_DNPCvARPC$p_val
DEG_DNPCvARPC_pvalue=as.numeric(DEG_DNPCvARPC_pvalue)
DEG_DNPCvARPC_BH = p.adjust(DEG_DNPCvARPC_pvalue, "BH")
write.csv(DEG_DNPCvARPC_BH, "DEG_DNPCvARPC_BH.csv")

##DEGs NEPC vs ARPC
Idents(object = mCRPC) <- "ARNE"
DimPlot(mCRPC, reduction = "umap", pt.size = 1, label = TRUE)
NEPC_vs_ARPC_markers <- FindMarkers(mCRPC, test.use = "wilcox", ident.1 = "NEPC", ident.2 = 
                                      "ARPC", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.csv(NEPC_vs_ARPC_markers, "NEPC_vs_ARPC_Markers_NEW.csv")

#p.adjust
DEG_NEPCvARPC <- read.csv("NEPC_vs_ARPC_Markers_NEW.csv") 
DEG_NEPCvARPC_pvalue <- DEG_NEPCvARPC$p_val
DEG_NEPCvARPC_pvalue=as.numeric(DEG_NEPCvARPC_pvalue)
DEG_NEPCvARPC_BH = p.adjust(DEG_NEPCvARPC_pvalue, "BH")
write.csv(DEG_NEPCvARPC_BH, "DEG_NEPCvARPC_BH.csv")

##DEGs DNPC vs NEPC
Idents(object = mCRPC) <- "ARNE"
DimPlot(mCRPC, reduction = "umap", pt.size = 1, label = TRUE)
DNPC_vs_NEPC_markers <- FindMarkers(mCRPC, test.use = "wilcox", ident.1 = "DNPC", ident.2 = 
                                      "NEPC", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.csv(DNPC_vs_NEPC_markers, "DNPC_vs_NEPC_Markers_NEW.csv")

#p.adjust
DEG_DNPCvNEPC <- read.csv("DNPC_vs_NEPC_Markers_NEW.csv") 
DEG_DNPCvNEPC_pvalue <- DEG_DNPCvNEPC$p_val
DEG_DNPCvNEPC_pvalue=as.numeric(DEG_DNPCvNEPC_pvalue)
DEG_DNPCvNEPC_BH = p.adjust(DEG_DNPCvNEPC_pvalue, "BH")
write.csv(DEG_DNPCvNEPC_BH, "DEG_DNPCvNEPC_BH.csv")

##DEGs ALNN vs ARPC
Idents(object = mCRPC) <- "ARNE"
ALNN_vs_ARPC_markers <- FindMarkers(mCRPC, test.use = "wilcox", ident.1 = "ALNN", ident.2 = 
                                      "ARPC", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.csv(ALNN_vs_ARPC_markers, "ALNN_vs_ARPC_markers_NEW.csv")

#p.adjust
DEG_ALNNvARPC <- read.csv("ALNN_vs_ARPC_markers_NEW.csv") 
DEG_ALNNvARPC_pvalue <- DEG_ALNNvARPC$p_val
DEG_ALNNvARPC_pvalue=as.numeric(DEG_ALNNvARPC_pvalue)
DEG_ALNNvARPC_BH = p.adjust(DEG_ALNNvARPC_pvalue, "BH")
write.csv(DEG_ALNNvARPC_BH, "DEG_ALNNvARPC_BH.csv")

##DEGs APNP vs ARPC
Idents(object = mCRPC) <- "ARNE"
APNP_vs_ARPC_markers <- FindMarkers(mCRPC, test.use = "wilcox", ident.1 = "APNP", ident.2 = 
                                      "ARPC", min.pct = 0.01, logfc.threshold = 0.01, only.pos = F)
write.csv(APNP_vs_ARPC_markers, "APNP_vs_ARPC_markers_NEW.csv")

#p.adjust
DEG_APNPvARPC <- read.csv("APNP_vs_ARPC_markers_NEW.csv") 
DEG_APNPvARPC_pvalue <- DEG_APNPvARPC$p_val
DEG_APNPvARPC_pvalue=as.numeric(DEG_APNPvARPC_pvalue)
DEG_APNPvARPC_BH = p.adjust(DEG_APNPvARPC_pvalue, "BH")
write.csv(DEG_APNPvARPC_BH, "DEG_APNPvARPC_BH.csv")
