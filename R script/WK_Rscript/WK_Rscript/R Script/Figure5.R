####Figure 5#### 

#Fig5A
tiff(file = "PIN2 Epi Highlited TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PIN2, reduction = "tsne", pt.size = 0.3, cols = c("blue", "blue", "blue", "blue", "blue", "blue", "grey", "grey", "blue", "grey", "blue", "grey", "grey", "grey", "grey", "grey", "blue"))
dev.off()

#Fig5B
tiff(file = "EpionlyPIN TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(EpionlyPIN, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

#Fig5C
tiff(file = "EpionlyPIN hAR Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(EpionlyPIN, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "EpionlyPIN EGFP Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(EpionlyPIN, reduction = "tsne", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Fig5D
tiff(file = "EpionlyPIN Trajectory Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

#Fig5E
tiff(file = "EpionlyPIN Trajectory seurat_clusters.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

#Fig5F
tiff(file = "EpionlyPIN Trajectory ARQ.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, markers = "ARQ", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "EpionlyPIN Trajectory EGFP.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, markers = "EGFP", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "EpionlyPIN Trajectory Krt5.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "EpionlyPIN Trajectory Krt8.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "EpionlyPIN Trajectory Krt19.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "EpionlyPIN Trajectory Pbsn.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(PINEpiPseudo, markers = "Pbsn", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#Fig5G
tiff(file = "PIN2 BE Highlited TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(PIN2, reduction = "tsne", pt.size = 0.3, cols = c("grey", "blue", "blue", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))
dev.off()

#Fig5H
tiff(file = "BEonlyPIN TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))
dev.off()

# Cytotrace
ARQPosBEARQNegBEPIN2Expdata = data.frame(ARQPosBEARQNegBEPIN2[["RNA"]]@data)
ARQPosBEARQNegBEPIN2phenotypedata = FetchData(ARQPosBEARQNegBEPIN2, c("seurat_clusters"))
write.table(ARQPosBEARQNegBEPIN2Expdata, file = "ARQPosBEARQNegBEPIN2Expdata.txt")
write.table(ARQPosBEARQNegBEPIN2phenotypedata, file = "ARQPosBEARQNegBEPIN2phenotypedata.txt")

#Fig5I

#FeaturePlot
tiff(file = "BEonlyPIN ARQ.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

tiff(file = "BEonlyPIN EGFP.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Fig5J

#Heatmap
tiff(file = "Cluster4v2 Heatmap50.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(ARQPosBE5ARQNegBE3PIN2, features = c(ARQPosBE5ARQNegBE3PIN2Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

#Fig5N

#vlnPlot
tiff(file = "Cluster4v2BE Igf1r Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Igf1r", pt.size = 0)
dev.off()

tiff(file = "Cluster4v2BE Fos Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Fos", pt.size = 0)
dev.off()

tiff(file = "Cluster4v2BE Mapk13 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Mapk13", pt.size = 0)
dev.off()

tiff(file = "Cluster4v2BE Jak2 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Jak2", pt.size = 0)
dev.off()




#Etcs

tiff(file = "Cluster4v2BE Wnt4 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Wnt4", pt.size = 0)
dev.off()

tiff(file = "Cluster4v2BE Nek7 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Nek7", pt.size = 0)
dev.off()

tiff(file = "Cluster4v2BE Scd1 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Scd1", pt.size = 0)
dev.off()

tiff(file = "Cluster4v2BE Hpgd Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Hpgd", pt.size = 0)
dev.off()

to_be_tested <- row.names(subset(fData(PINEpiPseudo), gene_short_name %in% c("ARQ")))
cds_subset <- PINEpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "EpionlyPIN plot genes Pseudotime ARQ.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

to_be_tested <- row.names(subset(fData(PINEpiPseudo), gene_short_name %in% c("EGFP")))
cds_subset <- PINEpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "EpionlyPIN plot genes Pseudotime EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

to_be_tested <- row.names(subset(fData(PINEpiPseudo), gene_short_name %in% c("Krt5")))
cds_subset <- PINEpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "EpionlyPIN plot genes Pseudotime Krt5.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

to_be_tested <- row.names(subset(fData(PINEpiPseudo), gene_short_name %in% c("Krt8")))
cds_subset <- PINEpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "EpionlyPIN plot genes Pseudotime Krt8.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

to_be_tested <- row.names(subset(fData(PINEpiPseudo), gene_short_name %in% c("Krt19")))
cds_subset <- PINEpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "EpionlyPIN plot genes Pseudotime Krt19.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

to_be_tested <- row.names(subset(fData(PINEpiPseudo), gene_short_name %in% c("Pbsn")))
cds_subset <- PINEpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "EpionlyPIN plot genes Pseudotime Pbsn.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))
dev.off()

