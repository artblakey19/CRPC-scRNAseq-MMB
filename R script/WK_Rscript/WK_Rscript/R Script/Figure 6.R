###Subset ARQ Positive cells in PIN2###

#Fig6A
tiff(file = "PIN2 Epcam Exp.tiff", width = 6, height = 6, units = "in", compression = "lzw", res = 800)
FeaturePlot(PIN2, reduction = "tsne", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()

#Fig6B
tiff(file = "EpionlyPINARQPos TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(EpionlyPINARQPos, reduction = "tsne", pt.size = 0.3) 
dev.off()

#Fig6C
tiff(file = "ARQPosPIN2 TSNE.tiff", width = 6.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQonlyPIN2, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF")) 
dev.off()

#Fig6D
tiff(file = "ARQPosPINPseudo Trajectory Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

#Fig6E
tiff(file = "ARQPosPINPseudo Trajectory Seurat.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE)
dev.off()

#Fig6F
tiff(file = "ARQPosPINPseudo Trajectory ARQ.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, markers = "ARQ", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQPosPINPseudo Trajectory EGFP.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, markers = "EGFP", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQPosPINPseudo Trajectory Krt5.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQPosPINPseudo Trajectory Krt8.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQPosPINPseudo Trajectory Krt19.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQPosPINPseudo Trajectory Pbsn.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQPINEpiPseudo, markers = "Pbsn", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#Fig6G
tiff(file = "EpiARQPos TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQPosEpi, reduction = "tsne", pt.size = 0.3) + NoLegend()
dev.off()

#Fig6H
tiff(file = "REcluster EpiARQPos TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQPosEpi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF"))
dev.off()

#Fig6I
tiff(file = "EpiARQPos stim TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(ARQPosEpi, reduction = "tsne", group.by = "stim", cols = c("#B71800", "#009FB7"))
dev.off()

#Fig6J
tiff(file = "EpiARQPosPIN PINonly TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(EpiARQPosPIN, reduction = "tsne", cols = c("#B71800")) + NoLegend()
dev.off()
tiff(file = "EpiARQPosTumor Tumoronly TSNE.tiff", width = 7.5, height = 6, units = "in", compression = "lzw", res = 800)
DimPlot(EpiARQPosTumor, reduction = "tsne", cols = c("#009FB7")) + NoLegend()
dev.off()

#Fig6K
tiff(file = "ARQEpiPseudo Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

#Fig6L
tiff(file = "ARQEpiPseudo Seurat.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF"))
dev.off()

#Fig6M
tiff(file = "ARQEpiPseudo Trajectory Krt5.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Krt8.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Krt18.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Krt18", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Krt19.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#Fig6N
tiff(file = "ARQEpiPseudo Trajectory Igf1r.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Igf1r", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory fos.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Fos", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Mapk13.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Mapk13", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Jak2.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Jak2", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "ARQEpiPseudo Trajectory Axin2.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Axin2", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Ccnd1.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Ccdn1", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Tcf4.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Tcf4", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()
tiff(file = "ARQEpiPseudo Trajectory Lgr5.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(ARQEpiPseudo, markers = "Lgr5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()


####Etcs####

#plot_genes_in_pseudotime

to_be_tested <- row.names(subset(fData(PINPseudo), gene_short_name %in% c("Igf1r")))
cds_subset <- PINPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))

tiff(file = "ARQPosPIN2 plot genes Pseudotime Igf1r.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))
dev.off()

to_be_tested <- row.names(subset(fData(PINPseudo), gene_short_name %in% c("Fos")))
cds_subset <- PINPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "ARQPosPIN2 plot genes Pseudotime Fos.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))
dev.off()

to_be_tested <- row.names(subset(fData(PINPseudo), gene_short_name %in% c("Mapk13")))
cds_subset <- PINPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "ARQPosPIN2 plot genes Pseudotime Mapk13.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))
dev.off()

to_be_tested <- row.names(subset(fData(PINPseudo), gene_short_name %in% c("Jak2")))
cds_subset <- PINPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "ARQPosPIN2 plot genes Pseudotime Jak2.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))
dev.off()

#FeaturePlot
FeaturePlot(ARQPosPIN2, reduction = "tsne", features = c("Igf1r"), cols = c("light grey", "purple"),  max.cutoff = "q90", pt.size = 0.5)



#Plot_genes_branched_pseudotime
PIN_genes <- row.names(subset(fData(PINPseudo),
                              gene_short_name %in% c("Igf1r", "Fos")))

tiff(file = "ARQPosPIN2 plot genes branched Pseudotime Igf1r Fos.tiff", width = 4, height = 8, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(PINPseudo[PIN_genes,],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)
dev.off()

PIN_genes <- row.names(subset(fData(PINPseudo),
                              gene_short_name %in% c("Mapk13", "Jak2")))
tiff(file = "ARQPosPIN2 plot genes branched Pseudotime Mapk13 Jak2.tiff", width = 4, height = 8, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(PINPseudo[PIN_genes,],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)
dev.off()




#Heatmap
Idents(object = ARQPosPIN2) <- "seurat_clusters"
Cluster5v4 <- subset(ARQPosPIN2, idents = c("5", "4"))
DefaultAssay(Cluster5v4) <- "RNA"
Cluster5v4 <- ScaleData(Cluster5v4, features = rownames(Cluster5v4))
Cluster5v4.markers <- FindAllMarkers(Cluster5v4, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Cluster5v4Top50 <- Cluster5v4.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

tiff(file = "Cluster4v2 Heatmap50.tiff", width = 8, height = 10, units = "in", compression = "lzw", res = 200)
DoHeatmap(Cluster5v4, features = c(Cluster5v4Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))
dev.off()

#vlnPlot
tiff(file = "Cluster5v2 Nucb2 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Nucb2", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Atf3 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Atf3", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Hsp90b1 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Hsp90b1", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Fos Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Fos", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Mapk13 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Mapk13", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Scd1 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Scd1", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Pgap2 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Pgap2", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Hpgd Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Hpgd", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Serinc3 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Serinc3", pt.size = 0)
dev.off()
tiff(file = "Cluster5v2 Hsp90b1 Vln.tiff", width = 3, height = 2, units = "in", compression = "lzw", res = 800)
VlnPlot(Cluster5v4, features = "Hsp90b1", pt.size = 0)
dev.off

#
PIN_genes <- row.names(subset(fData(PINPseudo),
                              gene_short_name %in% c("Nucb2", "Atf3", "Hsp90b1")))
plot_genes_branched_pseudotime(PINPseudo[PIN_genes,],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)

PIN_genes <- row.names(subset(fData(PINPseudo),
                              gene_short_name %in% c("Fos", "Mapk13", "Scd1")))
plot_genes_branched_pseudotime(PINPseudo[PIN_genes,],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)

PIN_genes <- row.names(subset(fData(PINPseudo),
                              gene_short_name %in% c("Pgap2", "Hpgd", "Serinc3")))
plot_genes_branched_pseudotime(PINPseudo[PIN_genes,],
                               branch_point = 1,
                               color_by = "Pseudotime",
                               ncol = 1)


