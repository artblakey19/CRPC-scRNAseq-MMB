####mGFP+AR+vmGFP+AR-vmGFP-AR+vmGFP-AR- in ArCre at E17.5####

#Setup workspace to make file calling & saving easy

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ArCreER/ArCreE17")

Idents(object = ARCre2) <- "OverAllCellType"
tiff(file = "ARCre2 tsne OverAllCellType.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARCre2, reduction = "tsne", pt.size = 0.3) 
dev.off()

Idents(object = ARCre2) <- "ArEGFPExp"
tiff(file = "ARCre2 tsne ArEGFPExp.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARCre2, reduction = "tsne", pt.size = 0.3) 
dev.off()

Idents(object = ARCre2) <- "OverAllCellType"
UGE <- subset(ARCre2, idents = c("UGE"))

Idents(object = UGE) <- "ArEGFPExp"
tiff(file = "UGE tsne ArEGFPExp.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(UGE, reduction = "tsne", pt.size = 0.3) 
dev.off()

table(Idents(UGE))

Idents(object = ARCre2) <- "OverAllCellType"
FBSM <- subset(ARCre2, idents = c("Fb", "SM"))

Idents(object = FBSM) <- "ArEGFPExp"
tiff(file = "FBSM tsne ArEGFPExp.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(FBSM, reduction = "tsne", pt.size = 0.3) 
dev.off()

table(Idents(FBSM))

#### AR EGFP Exp Analysis ####
DefaultAssay(ARCre2) <- "RNA"
DoublePos1 <- subset(x=ARCre2, subset = Ar > 0 & EGFP > 1)
Aronly1 <- subset(x=ARCre2, subset = Ar > 0 & EGFP < 1)
EGFPonly1 <- subset(x=ARCre2, subset = Ar == 0 & EGFP > 1)
DoubleNeg1 <- subset(x=ARCre2, subset = Ar == 0 & EGFP < 1)
Idents(object = DoublePos1) <- "DoublePos1"
Idents(object = Aronly1) <- "Aronly1"
Idents(object = EGFPonly1) <- "EGFPonly1"
Idents(object = DoubleNeg1) <- "DoubleNeg1"
DoublePos1[["EGFPArExp"]] <- Idents(object = DoublePos1)
Aronly1[["EGFPArExp"]] <- Idents(object = Aronly1)
EGFPonly1[["EGFPArExp"]] <- Idents(object = EGFPonly1)
DoubleNeg1[["EGFPArExp"]] <- Idents(object = DoubleNeg1)
Temp4 <- list(DoublePos1, DoubleNeg1, Aronly1, EGFPonly1)
EGFPAr1 <- merge(x = DoublePos1, y = Aronly1)
EGFPAr2 <- merge(x = EGFPonly1, y = DoubleNeg1)
EGFPAr <- merge(x = EGFPAr1, y = EGFPAr2)
ARCre2[["EGFPArExp"]] <- Idents(object = EGFPAr)
DimPlot(ARCre2, reduction = "tsne", group.by = "EGFPArExp")

Idents(object = ARCre2) <- "OverAllCellType"
UGE <- subset(ARCre2, idents = c("UGE"))

Idents(object = UGE) <- "EGFPArExp"
tiff(file = "UGE tsne ArEGFPExp.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(UGE, reduction = "tsne", pt.size = 0.3) 
dev.off()

table(Idents(UGE))

Idents(object = ARCre2) <- "OverAllCellType"
FBSM <- subset(ARCre2, idents = c("Fb", "SM"))

Idents(object = FBSM) <- "EGFPArExp"
tiff(file = "FBSM tsne EGFPArExp.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(FBSM, reduction = "tsne", pt.size = 0.3) 
dev.off()

table(Idents(FBSM))

####Subset/Rename FBSM for SCENIC####
Idents(object = FBSM) <- "ArEGFPExp"
FBSM2 <- subset(FBSM, idents = c("DoubleNeg", "DoublePos", "EGFPonly"))
DimPlot(FBSM2, reduction = "tsne", pt.size = 0.3) 

Idents(object = FBSM2) <- "ArEGFPExp"
FBSM2 <- RenameIdents(object = FBSM2, 'DoublePos' = "mGFPPosArPos", 'EGFPonly' = "mGFPPosArNeg",
                                     'DoubleNeg' = "mGFPNegArNeg")
FBSM2[["CellType"]] <- Idents(object = FBSM2)

Idents(object = FBSM2) <- "CellType"
DimPlot(FBSM2, reduction = "tsne", pt.size = 0.3) 

####DEGs####

#mGFPPosArPosvsmGFPPosArNeg
Idents(object = FBSM2) <- "CellType"
DefaultAssay(FBSM2) <- "RNA"
mGFPPosArPosvArNeg.0.1.Markers <- FindMarkers(FBSM2, ident.1 = "mGFPPosArPos", ident.2 = "mGFPPosArNeg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(mGFPPosArPosvArNeg.0.1.Markers, "mGFPPosArPosvArNeg.0.1.Markers.csv")
mGFPPosArPosvArNeg.0.Markers <- FindMarkers(FBSM2, ident.1 = "mGFPPosArPos", ident.2 = "mGFPPosArNeg", min.pct = 0, logfc.threshold = 0)
write.csv(mGFPPosArPosvArNeg.0.Markers, "mGFPPosArPosvArNeg.0.Markers.csv")

#p.adjust
DEG_mGFPPosArPosvArNeg <- read.csv("mGFPPosArPosvArNeg.0.Markers.csv") 
DEG_mGFPPosArPosvArNeg_pvalue <- DEG_mGFPPosArPosvArNeg$p_val
DEG_mGFPPosArPosvArNeg_pvalue=as.numeric(DEG_mGFPPosArPosvArNeg_pvalue)
DEG_mGFPPosArPosvArNeg_pvalue_BH = p.adjust(DEG_mGFPPosArPosvArNeg_pvalue, "BH")
write.csv(DEG_mGFPPosArPosvArNeg_pvalue_BH, "DEG_mGFPPosArPosvArNeg_pvalue_BH.csv")

#mGFPPosArPosvsmGFPNegArNeg
Idents(object = FBSM2) <- "CellType"
DefaultAssay(FBSM2) <- "RNA"
mGFPPosArPosvmGFPNegArNeg.0.1.Markers <- FindMarkers(FBSM2, ident.1 = "mGFPPosArPos", ident.2 = "mGFPNegArNeg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(mGFPPosArPosvmGFPNegArNeg.0.1.Markers, "mGFPPosArPosvmGFPNegArNeg.0.1.Markers.csv")
mGFPPosArPosvmGFPNegArNeg.0.Markers <- FindMarkers(FBSM2, ident.1 = "mGFPPosArPos", ident.2 = "mGFPNegArNeg", min.pct = 0, logfc.threshold = 0)
write.csv(mGFPPosArPosvmGFPNegArNeg.0.Markers, "mGFPPosArPosvmGFPNegArNeg.0.Markers.csv")

#p.adjust
DEG_mGFPPosArPosvmGFPNegArNeg <- read.csv("mGFPPosArPosvmGFPNegArNeg.0.Markers.csv") 
DEG_mGFPPosArPosvmGFPNegArNeg_pvalue <- DEG_mGFPPosArPosvmGFPNegArNeg$p_val
DEG_mGFPPosArPosvmGFPNegArNeg_pvalue=as.numeric(DEG_mGFPPosArPosvmGFPNegArNeg_pvalue)
DEG_mGFPPosArPosvmGFPNegArNeg_pvalue_BH = p.adjust(DEG_mGFPPosArPosvmGFPNegArNeg_pvalue, "BH")
write.csv(DEG_mGFPPosArPosvmGFPNegArNeg_pvalue_BH, "DEG_mGFPPosArPosvmGFPNegArNeg_BH.csv")

#mGFPPosArNegvsmGFPNegArNeg
Idents(object = FBSM2) <- "CellType"
DefaultAssay(FBSM2) <- "RNA"
mGFPPosArNegvmGFPNegArNeg.0.1.Markers <- FindMarkers(FBSM2, ident.1 = "mGFPPosArNeg", ident.2 = "mGFPNegArNeg", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(mGFPPosArNegvmGFPNegArNeg.0.1.Markers, "mGFPPosArNegvmGFPNegArNeg.0.1.Markers.csv")
mGFPPosArNegvmGFPNegArNeg.0.Markers <- FindMarkers(FBSM2, ident.1 = "mGFPPosArNeg", ident.2 = "mGFPNegArNeg", min.pct = 0, logfc.threshold = 0)
write.csv(mGFPPosArNegvmGFPNegArNeg.0.Markers, "mGFPPosArNegvmGFPNegArNeg.0.Markers.csv")

#p.adjust
DEG_mGFPPosArNegvmGFPNegArNeg <- read.csv("mGFPPosArNegvmGFPNegArNeg.0.Markers.csv") 
DEG_mGFPPosArNegvmGFPNegArNeg_pvalue <- DEG_mGFPPosArNegvmGFPNegArNeg$p_val
DEG_mGFPPosArNegvmGFPNegArNeg_pvalue=as.numeric(DEG_mGFPPosArNegvmGFPNegArNeg_pvalue)
DEG_mGFPPosArNegvmGFPNegArNeg_pvalue_BH = p.adjust(DEG_mGFPPosArNegvmGFPNegArNeg_pvalue, "BH")
write.csv(DEG_mGFPPosArNegvmGFPNegArNeg_pvalue_BH, "DEG_mGFPPosArNegvmGFPNegArNeg_BH.csv")
