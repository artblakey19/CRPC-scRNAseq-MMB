#Fig4A

tiff(file = "PINvTumor Trajectory Pseudotime.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)
dev.off()

#Fig4B

tiff(file = "PINvTumor Trajectory seurat_clusters.tiff", width = 8, height = 8, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

#Fig4C
tiff(file = "PINvTumor Trajectory ARQ.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "ARQ", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory EGFP.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "EGFP", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Krt5.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Krt8.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Krt19.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Pbsn.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Pbsn", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#Fig4D
tiff(file = "branched heatmap EpiPseudo.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_genes_branched_heatmap(EpiPseudo[row.names(BEAM_res)[1:50]], branch_point = 1, num_clusters = 6, cores=4, use_gene_short_name=TRUE, show_rownames=TRUE)
dev.off()

#Fig4F
tiff(file = "PINvTumor Trajectory Axin2.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Axin2", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Tcf4.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Tcf4", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Ccnd1.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Ccnd1", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

tiff(file = "PINvTumor Trajectory Lgr5.tiff", width = 4, height = 5, units = "in", compression = "lzw", res = 800)
plot_cell_trajectory(EpiPseudo, markers = "Lgr5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
dev.off()

#Fig4G
EpiPseudo_genes <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Tcf4", "Axin2")))

tiff(file = "PINvTumor plot_genes_branched_pseudotime Tcf4 Axin2.tiff", width = 4, height = 6, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(EpiPseudo[EpiPseudo_genes,],
                               branch_point = 2,
                               color_by = "seurat_clusters",
                               ncol = 1)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

EpiPseudo_genes <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Ccnd1", "Lgr5")))

tiff(file = "PINvTumor plot_genes_branched_pseudotime Ccnd1 Lgr5.tiff", width = 4, height = 6, units = "in", compression = "lzw", res = 800)
plot_genes_branched_pseudotime(EpiPseudo[EpiPseudo_genes,],
                               branch_point = 2,
                               color_by = "seurat_clusters",
                               ncol = 1)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

#Fig4H

#Co-expression FeaturePlot
DefaultAssay(PINvTumor.combined.Epi) <- "RNA"
tiff(file = "PINvTumor ARQ Tcf4 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi, reduction = "tsne", features = c("ARQ", "Tcf4"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "tomato3", "olivedrab4"), pt.size = .5, blend.threshold = 0.1)
dev.off()

tiff(file = "PINvTumor ARQ Axin2 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi, reduction = "tsne", features = c("ARQ", "Axin2"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "tomato3", "olivedrab4"), pt.size = .5, blend.threshold = 0.1)
dev.off()

tiff(file = "PINvTumor ARQ Ccnd1 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi, reduction = "tsne", features = c("ARQ", "Ccnd1"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "tomato3", "olivedrab4"), pt.size = .5, blend.threshold = 0.1)
dev.off()

tiff(file = "PINvTumor ARQ Lgr5 Blend.tiff", width = 12, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(PINvTumor.combined.Epi, reduction = "tsne", features = c("ARQ", "Lgr5"), blend = TRUE, max.cutoff = "q90", cols = c("light grey", "tomato3", "olivedrab4"), pt.size = .5, blend.threshold = 0.1)
dev.off()


#etc

#Cytotrace
PINvTumor.combined.EpiExpdata = data.frame(PINvTumor.combined.Epi[["RNA"]]@data)
PINvTumor.combined.Epiphenotypedata = FetchData(PINvTumor.combined.Epi, c("seurat_clusters"))
write.table(PINvTumor.combined.EpiExpdata, file = "PINvTumor.combined.EpiExpdata.txt")
write.table(PINvTumor.combined.Epiphenotypedata, file = "PINvTumor.combined.Epiphenotypedata.txt")

#plot_genes_in_pseudotime

to_be_tested <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("ARQ")))
cds_subset <- EpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "PINvTumor plot genes Pseudotime ARQ.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

to_be_tested <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("EGFP")))
cds_subset <- EpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "PINvTumor plot genes Pseudotime EGFP.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

to_be_tested <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Krt5")))
cds_subset <- EpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "PINvTumor plot genes Pseudotime Krt5.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

to_be_tested <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Krt8")))
cds_subset <- EpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "PINvTumor plot genes Pseudotime Krt8.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

to_be_tested <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Krt19")))
cds_subset <- EpiPseudo[to_be_tested,]
diff_test_res <- differentialGeneTest(cds_subset,fullModelFormulaStr = "~sm.ns(Pseudotime)")

tiff(file = "PINvTumor plot genes Pseudotime Krt19.tiff", width = 4, height = 4, units = "in", compression = "lzw", res = 800)
plot_genes_in_pseudotime(cds_subset, color_by = "seurat_clusters") + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))
dev.off()

