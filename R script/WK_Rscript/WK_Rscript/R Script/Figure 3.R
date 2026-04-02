#Fig3A

#Tsne Plots
tiff(file = "Combined Epi Highlited TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 0.3, cols = c("blue", "blue", "grey", "grey", "grey", "grey", "grey"))
dev.off()

#Fig3B
tiff(file = "Epirecluster TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINvTumor.combined.Epi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

#Fig3C
tiff(file = "PINonlyEpi TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PINonlyEpi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

tiff(file = "TumoronlyEpi TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(TumoronlyEpi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

#Fig3E

#DotPlot
tiff(file = "CellType DotPlot.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
DotPlot(PINvTumor.combined.Epi, features = c("Cd53", "Sh2d2a", "Cytip", "Srgn", "Rgs1", "Gjb2", "Timp4", "Cxcl17", "Oit1", "Cyp2f2", "Cdca3", "Birc5", "Ube2c", "Stmn1", "Mki67", "Myof", "Mgat4a", "Bace2", "Cyba", "Wfdc2", "Hoxb13", "Hmgcs2", "Ceacam2", "Mme", "Apof", "Pnliprp1", "C1s2", "C1rb", "Pbsn", "Tgm4", "Pmaip1", "Crip1", "Cd55", "Kctd14", "Sbspon", "Apoc4", "Nkain4", "Gulo", "Npl", "Mt3", "Ptgds", "Svs4", "Gpx3", "Svs3a", "Wfdc15b", "Nxf7", "Syngr1", "Msx2", "Defa20", "Wif1", "Lamb3", "Htra1", "Lbp", "Ltbp4", "Sult5a1", "Col17a1", "Cldn1", "Tpm2", "Krt17", "Krt14", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()
dev.off()

#Fig3F

#FeaturePlot
DefaultAssay(PINonlyEpi) <- "RNA"

tiff(file = "PIN hAR Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "PIN Cdh1 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Cdh1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "PIN CK5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "PIN CK8 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "PIN Pbsn Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "PIN Mki67 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

DefaultAssay(TumoronlyEpi) <- "RNA"

tiff(file = "Tumor hAR Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Tumor Cdh1 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Cdh1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Tumor CK5 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Tumor CK8 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Tumor Pbsn Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "Tumor Mki67 Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()




