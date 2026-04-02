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
