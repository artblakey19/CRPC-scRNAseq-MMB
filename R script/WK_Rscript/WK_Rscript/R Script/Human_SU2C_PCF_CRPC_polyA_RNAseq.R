#Load data
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Human PCa/SU2C_PCF")

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

#Heatmap for AR signature
ARlist <- list(c("KLK3", "KLK2", "TMPRSS2", "FKBP5", "NKX3-1", "PMEPA1", "AR", "ALDH1A3", "STEAP4"))
mCRPC <- AddModuleScore(object = mCRPC, features = ARlist, name = "ARDownstream") 
mCRPC[['ARmodule']] <- CreateAssayObject(data = t(x = FetchData(object = mCRPC, vars = 'ARDownstream1')))
all.genes <- rownames(mCRPC)
mCRPC <- ScaleData(mCRPC, features = all.genes)
tiff(file = "mCRPC ARscore Heatmap.tiff", width = 5, height = 1.5, units = "in", compression = "lzw", res =200)
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

#Heatmap for MET downstream
tiff(file = "mCRPC HGF selected downstream Heatmap.tiff", width = 5.5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("MET", "AQP9","RUNX1", "SERPINE1",
                              "ABCC3",  "SLC43A3",
                              "GPNMB","RASGRP3"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

#Heatmap for WNT downstream
tiff(file = "mCRPC WNT selected downstream Heatmap.tiff", width = 5.5, height = 3.5, units = "in", compression = "lzw", res =200)
DoHeatmap(mCRPC, features = c("LEF1", "TCF7L2", "PLAUR", "MMP7",
                              "LGR5", "FST", "BMP2"), draw.lines = TRUE, size = 3, disp.max = 1.2, disp.min = -1.2) + scale_fill_gradientn(colors = c("purple", "black", "yellow"), na.value = "white")
dev.off()

