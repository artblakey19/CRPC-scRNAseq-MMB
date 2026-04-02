#WT ARQ-Osr1 SCseq Work Flow


#### Add necessary tools to library ####

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(RColorBrewer)

#### Initial Filtering and Clustering ####

#Setup workspace to make file calling & saving easy

setwd("W:/ARQ9_scRnAseq_WK")

#PIN 

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

PIN.data <- Read10X("//isi-dcnl/user_data/zjsun/BIC/AR Transgene Tumorigenesis/ARQ9-Osr1 SingleCell Seq/PCa_AR_transgene_analysis/PIN4_outs/filtered_feature_bc_matrix")
PIN <- CreateSeuratObject(counts = PIN.data,  min.cells = 3, min.features = 500, project = "PIN")
PIN <- NormalizeData(PIN)

#Start here if loading an R object

#Initial processing & filtering

PIN[["percent.mt"]] <- PercentageFeatureSet(PIN, pattern = "^mt-")
VlnPlot(PIN, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

hist(PIN@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Pre-filteration")

hist(PIN@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Pre-filteration")

PIN2 <- subset(PIN, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
VlnPlot(PIN2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

hist(PIN2@meta.data$nFeature_RNA, breaks = 100, col = "pink", xlab = "nFeature_RNA", main = "PIN Post-filteration")

hist(PIN2@meta.data$percent.mt, breaks = 100, col = "pink", xlab = "percent.mt", main = "PIN Post-filteration")

PIN2 <- FindVariableFeatures(PIN2, selection.method = "vst", nfeatures = 5000)

#Clustering
PIN2 <- ScaleData(PIN2, verbose = FALSE)
PIN2 <- RunPCA(PIN2, npcs = 30, verbose = FALSE)
PIN2 <- RunUMAP(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- FindNeighbors(PIN2, reduction = "pca", dims = 1:20)
PIN2 <- FindClusters(PIN2, resolution = 0.5)
PIN2 <- JackStraw(PIN2, num.replicate = 100)
PIN2 <- ScoreJackStraw(PIN2, dims = 1:20)
JackStrawPlot(PIN2, dims = 1:20)

PIN2 <- RunTSNE(PIN2, reduction = "pca", dims = 1:20)
DimPlot(PIN2, reduction = "tsne", label = TRUE)

#Tumor

#Following steps are for starting with filtered_feature_bc_matrix
# Files in matrix 1) barcodes.tsv.gz 2) features.tsv.gz 3) matrix.mtx.gz

Tumor.data <- Read10X("//isi-dcnl/user_data/zjsun/BIC/AR Transgene Tumorigenesis/ARQ9-Osr1 SingleCell Seq/PCa_AR_transgene_analysis/Tumor4_outs/filtered_feature_bc_matrix")
Tumor <- CreateSeuratObject(counts = Tumor.data,  min.cells = 3, min.features = 500, project = "Tumor")
Tumor <- NormalizeData(Tumor)

#Start here if loading an R object

#Initial processing & filtering

Tumor[["percent.mt"]] <- PercentageFeatureSet(Tumor, pattern = "^mt-")
VlnPlot(Tumor, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

hist(Tumor@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Pre-filteration")

hist(Tumor@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Pre-filteration")

Tumor2 <- subset(Tumor, subset = nFeature_RNA > 1000 & nFeature_RNA < 7000 & percent.mt < 10)
VlnPlot(Tumor2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

hist(Tumor2@meta.data$nFeature_RNA, breaks = 100, col = "lightblue", xlab = "nFeature_RNA", main = "Tumor Post-filteration")

hist(Tumor2@meta.data$percent.mt, breaks = 100, col = "lightblue", xlab = "percent.mt", main = "Tumor Post-filteration")

Tumor2 <- FindVariableFeatures(Tumor2, selection.method = "vst", nfeatures = 5000)
Tumor2 <- ScaleData(Tumor2, verbose = FALSE)
Tumor2 <- RunPCA(Tumor2, npcs = 30, verbose = FALSE)
Tumor2 <- RunUMAP(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- FindNeighbors(Tumor2, reduction = "pca", dims = 1:20)
Tumor2 <- FindClusters(Tumor2, resolution = 0.5)
Tumor2 <- JackStraw(Tumor2, num.replicate = 100)
Tumor2 <- ScoreJackStraw(Tumor2, dims = 1:20)
JackStrawPlot(Tumor2, dims = 1:20)

Tumor2 <- RunTSNE(Tumor2, reduction = "pca", dims = 1:20)
DimPlot(Tumor2, reduction = "tsne", label = TRUE)

#### Merging Datasets ####

#Stash old idents
Tumor2[["orig.clusters"]] <- Idents(object = Tumor2)
PIN2[["orig.clusters"]] <- Idents(object = PIN2)

#Set Current idents
Idents(object = PIN2) <- "seurat_clusters"
Idents(object = Tumor2) <- "seurat_clusters"
PIN2$stim <- "PIN"
Tumor2$stim <- "Tumor"
PINvTumor.anchors <- FindIntegrationAnchors(object.list = list(PIN2, Tumor2), dims = 1:20)
PINvTumor.combined <- IntegrateData(anchorset = PINvTumor.anchors, dims = 1:20)
DefaultAssay(PINvTumor.combined) <- "integrated"

#Run the standard workflow for visualization and clustering
PINvTumor.combined <- ScaleData(PINvTumor.combined, verbose = FALSE)
PINvTumor.combined <- RunPCA(PINvTumor.combined, npcs = 30, verbose = FALSE)
# t-SNE and Clustering
PINvTumor.combined <- FindNeighbors(PINvTumor.combined, reduction = "pca", dims = 1:20)
PINvTumor.combined <- FindClusters(PINvTumor.combined, resolution = 0.5)
PINvTumor.combined <- RunTSNE(PINvTumor.combined, reduction = "pca", dims = 1:20)
DimPlot(PINvTumor.combined, reduction = "tsne", label = TRUE)

DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 0.3, group.by = "stim", cols = c("#B71800", "#009FB7")) 

#New labeling
new.cluster.ids <- c("LE", "LE", "LE", "BE", "LE", "Lym", "LE", "LE", "LE", "LE", "BE", "LE", "Fib", "LE", "Leu", "Endo", "SM", "BE")
names(new.cluster.ids) <- levels(PINvTumor.combined)
PINvTumor.combined <- RenameIdents(PINvTumor.combined, new.cluster.ids)
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 0.3)

Idents(object = PINvTumor.combined) <- "stim"
PINonly <- subset(PINvTumor.combined, idents = c("PIN"))
Idents(object = PINonly) <- "seurat_clusters"
new.cluster.ids <- c("LE", "LE", "LE", "BE", "LE", "Lym", "LE", "LE", "LE", "LE", "BE", "LE", "Fib", "LE", "Leu", "Endo", "SM", "BE")
names(new.cluster.ids) <- levels(PINonly)
PINonly <- RenameIdents(PINonly, new.cluster.ids)
DimPlot(PINonly, reduction = "tsne", pt.size = 0.3)

Idents(object = PINvTumor.combined) <- "stim"
Tumoronly <- subset(PINvTumor.combined, idents = c("Tumor"))
Idents(object = Tumoronly) <- "seurat_clusters"
new.cluster.ids <- c("LE", "LE", "LE", "BE", "LE", "Lym", "LE", "LE", "LE", "LE", "BE", "LE", "Fib", "LE", "Leu", "Endo", "SM", "BE")
names(new.cluster.ids) <- levels(Tumoronly)
Tumoronly <- RenameIdents(Tumoronly, new.cluster.ids)
DimPlot(Tumoronly, reduction = "tsne", pt.size = 0.3)

#Dotplot
Idents(object = PINvTumor.combined) <- "seurat_clusters"
PINvTumor.combined <- RenameIdents(object = PINvTumor.combined, '3' = "BE", '10' = "BE", '17' = "BE", '0' = "LE", '1' = "LE", '2' = "LE", '4' = "LE", '6' = "LE", '7' = "LE", '8' = "LE", '9' = "LE", '11' = "LE", '13' = "LE", '12' = "Fib", '16' = "SM", '14' = "Leu", '15' = "Endo", '5' = "Lym")
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 1, label = TRUE)
PINvTumor.combined[["CellType"]] <- Idents(object = PINvTumor.combined)

DotPlot(PINvTumor.combined, features = c("Cxcr6", "Cd28", "Cd3e", "Cd3d", "Ccl5", "Pecam1", "Cldn5", "Cdh5", "Plvap", "Aqp1", "C1qc", "C1qb", "C1qa", "Spi1", "Tyrobp", "Tagln", "Actg2", "Rgs5", "Myh11", "Acta2", "Fbln1", "Rspo3", "Pdgfra", "Apod", "Lum", "Cldn3", "Stard10", "Alcam", "Krt19", "Krt8", "Aqp3", "Col17a1", "Krt15" ,"Krt14", "Krt5", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

#FeaturePlot

Idents(object = PINvTumor.combined) <- "RNA"

FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Ar"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Epcam"), cols = c("light grey", "red"),min.cutoff = 0, pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Myh11"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINvTumor.combined, reduction = "tsne", features = c("Fbln1"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Cell counts
Idents(object = PINvTumor.combined) <- "stim"
PINvTumor.combined$stim.seurat_clusters <- paste(Idents(PINvTumor.combined), PINvTumor.combined$seurat_clusters, sep = "_")
Idents(object = PINvTumor.combined) <- "stim.seurat_clusters"
table(Idents(PINvTumor.combined))

#Reclustering Epithelial cells
Idents(object = PINvTumor.combined) <- "CellType"
DimPlot(PINvTumor.combined, reduction = "tsne", pt.size = 0.3, cols = c("blue", "blue", "grey", "grey", "grey", "grey", "grey"))

PINvTumor.combined.Epi <- subset(PINvTumor.combined, idents = c("0", "1", "2", "3", "4", "6", "7", "8", "9", "10", "11", "13", "17"))

Idents(object = PINvTumor.combined.Epi) <- "seurat_clusters"
DefaultAssay(PINvTumor.combined.Epi) <- "integrated"

# Run the standard workflow for visualization and clustering
PINvTumor.combined.Epi <- ScaleData(PINvTumor.combined.Epi, verbose = FALSE)
PINvTumor.combined.Epi <- RunPCA(PINvTumor.combined.Epi, npcs = 30, verbose = FALSE)

# t-SNE and Clustering
PINvTumor.combined.Epi <- FindNeighbors(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
PINvTumor.combined.Epi <- FindClusters(PINvTumor.combined.Epi, resolution = 0.5)
PINvTumor.combined.Epi <- RunTSNE(PINvTumor.combined.Epi, reduction = "pca", dims = 1:20)
DimPlot(PINvTumor.combined.Epi, reduction = "tsne", label = TRUE)

Idents(object = PINvTumor.combined.Epi) <- "stim"

PINonlyEpi <- subset(PINvTumor.combined.Epi, idents = c("PIN"))
Idents(object = PINonlyEpi) <- "seurat_clusters"

DimPlot(PINonlyEpi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))

TumoronlyEpi <- subset(PINvTumor.combined.Epi, idents = c("Tumor"))
Idents(object = TumoronlyEpi) <- "seurat_clusters"

DimPlot(TumoronlyEpi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))

#Cell counts
Idents(object = PINvTumor.combined.Epi) <- "stim"
PINvTumor.combined.Epi$stim.seurat_clusters <- paste(Idents(PINvTumor.combined.Epi), PINvTumor.combined.Epi$seurat_clusters, sep = "_")
Idents(object = PINvTumor.combined.Epi) <- "stim.seurat_clusters"
table(Idents(PINvTumor.combined.Epi))

#DotPlot
Idents(object = PINvTumor.combined.Epi) <- "seurat_clusters"
PINvTumor.combined.Epi <- RenameIdents(object = PINvTumor.combined.Epi, '4' = "BE1", '7' = "BE2", '0' = "LE1", '1' = "LE2", '2' = "LE3", '3' = "LE4", '5' = "LE5", '8' = "LE6", '9' = "LE7", '6' = "ProE", '10' = "Club", '11' = "OE")
DimPlot(PINvTumor.combined.Epi, reduction = "tsne", pt.size = 1, label = TRUE)
PINvTumor.combined.Epi[["CellType"]] <- Idents(object = PINvTumor.combined.Epi)

DotPlot(PINvTumor.combined.Epi, features = c("Cd53", "Sh2d2a", "Cytip", "Srgn", "Rgs1", "Gjb2", "Timp4", "Cxcl17", "Oit1", "Cyp2f2", "Cdca3", "Birc5", "Ube2c", "Stmn1", "Mki67", "Myof", "Mgat4a", "Bace2", "Cyba", "Wfdc2", "Hoxb13", "Hmgcs2", "Ceacam2", "Mme", "Apof", "Pnliprp1", "C1s2", "C1rb", "Pbsn", "Tgm4", "Pmaip1", "Crip1", "Cd55", "Kctd14", "Sbspon", "Apoc4", "Nkain4", "Gulo", "Npl", "Mt3", "Ptgds", "Svs4", "Gpx3", "Svs3a", "Wfdc15b", "Nxf7", "Syngr1", "Msx2", "Defa20", "Wif1", "Lamb3", "Htra1", "Lbp", "Ltbp4", "Sult5a1", "Col17a1", "Cldn1", "Tpm2", "Krt17", "Krt14", "EGFP", "ARQ", "Ar"), cols = c("light grey", "red")) + RotatedAxis()

#FeaturePlots
DefaultAssay(PINonlyEpi) <- "RNA"

FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Cdh1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(PINonlyEpi, reduction = "tsne", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

DefaultAssay(TumoronlyEpi) <- "RNA"

FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Cdh1"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Krt5"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Krt8"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Pbsn"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(TumoronlyEpi, reduction = "tsne", features = c("Mki67"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

#### PINvTumor Epithelial Pseudotime ####

DefaultAssay(PINvTumor.combined.Epi) <- "RNA"
EpiPseudo <- as.CellDataSet(PINvTumor.combined.Epi)
EpiPseudo <- detectGenes(EpiPseudo, min_expr = 0.1)
print(head(fData(EpiPseudo)))

expressed_genes <- row.names(subset(fData(EpiPseudo),
                                    num_cells_expressed >= 10))

pData(EpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(EpiPseudo))
EpiPseudo <- EpiPseudo[,pData(EpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(EpiPseudo), geom =
        "density")

EpiPseudo <- estimateSizeFactors(EpiPseudo)
EpiPseudo <- estimateDispersions(EpiPseudo)

disp_table <- dispersionTable(EpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
EpiPseudo <- setOrderingFilter(EpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(EpiPseudo)

# EpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(EpiPseudo, return_all = F) # norm_method='log'

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2, num_dim = 20,
                             reduction_method = 'tSNE', verbose = T)
EpiPseudo <- clusterCells(EpiPseudo, num_clusters = 2)

plot_cell_clusters(EpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(EpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
EpiPseudo <- setOrderingFilter(EpiPseudo, ordering_genes)
plot_ordering_genes(EpiPseudo)

EpiPseudo <- reduceDimension(EpiPseudo, max_components = 2,
                             method = 'DDRTree')

EpiPseudo <- orderCells(EpiPseudo)

GM_state <- function(EpiPseudo){
  if (length(unique(pData(EpiPseudo)$State)) > 1){
    T0_counts <- table(pData(EpiPseudo)$State, pData(EpiPseudo)$seurat_clusters)[,"7"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

EpiPseudo <- orderCells(EpiPseudo, root_state = GM_state(EpiPseudo))

#Visualization

plot_cell_trajectory(EpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

plot_cell_trajectory(EpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))

plot_cell_trajectory(EpiPseudo, markers = "ARQ", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "EGFP", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Pbsn", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

#Analyzing branches in single cell trajectory
BEAM_res <- BEAM(EpiPseudo, branch_point = 2, cores = 4)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
my_branched2_heatmap <- plot_genes_branched_heatmap(EpiPseudo[row.names(subset(BEAM_res, qval < 1e-4)),],
                                                    branch_point = 2,
                                                    num_clusters = 6,
                                                    cores = 4,
                                                    use_gene_short_name = TRUE,
                                                    show_rownames = TRUE,
                                                    return_heatmap = TRUE)
plot_genes_branched_heatmap(EpiPseudo[row.names(BEAM_res)[1:50]], branch_point = 2, num_clusters = 6, cores=4, use_gene_short_name=TRUE, show_rownames=TRUE)

table(my_branched2_heatmap$annotation_row$Cluster)
my_row2_6 <- my_branched2_heatmap$annotation_row
my_row2_6 <- data.frame(cluster = my_row2_6$Cluster,
                        gene = row.names(my_row2_6),
                        stringsAsFactors = FALSE)

Branch2.6cluster1.markers <- table(my_row2_6[my_row2_6$cluster == 1,'gene'])
write.table(Branch2.6cluster1.markers, file = "Branch2.6cluster1.markers")
Branch2.6cluster2.markers <- table(my_row2_6[my_row2_6$cluster == 2,'gene'])
write.table(Branch2.6cluster2.markers, file = "Branch2.6cluster2.markers")
Branch2.6cluster3.markers <- table(my_row2_6[my_row2_6$cluster == 3,'gene'])
write.table(Branch2.6cluster3.markers, file = "Branch2.6cluster3.markers")
Branch2.6cluster4.markers <- table(my_row2_6[my_row2_6$cluster == 4,'gene'])
write.table(Branch2.6cluster4.markers, file = "Branch2.6cluster4.markers")
Branch2.6cluster5.markers <- table(my_row2_6[my_row2_6$cluster == 5,'gene'])
write.table(Branch2.6cluster5.markers, file = "Branch2.6cluster5.markers")
Branch2.6cluster6.markers <- table(my_row2_6[my_row2_6$cluster == 6,'gene'])
write.table(Branch2.6cluster6.markers, file = "Branch2.6cluster6.markers")

#Expression level of Wnt downstreams in plot_cell_trajectory

plot_cell_trajectory(EpiPseudo, markers = "Tcf4", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Axin2", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Ccnd1", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(EpiPseudo, markers = "Lgr5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

#Expression level of Wnt downstreams in plot_genes_branched_pseudotime

EpiPseudo_genes <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Tcf4", "Axin2")))
plot_genes_branched_pseudotime(EpiPseudo[EpiPseudo_genes,],
                               branch_point = 2,
                               color_by = "seurat_clusters",
                               ncol = 1)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))

EpiPseudo_genes <- row.names(subset(fData(EpiPseudo), gene_short_name %in% c("Ccnd1", "Lgr5")))
plot_genes_branched_pseudotime(EpiPseudo[EpiPseudo_genes,],
                               branch_point = 2,
                               color_by = "seurat_clusters",
                               ncol = 1)+ scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF", "#FF5CFF"))

#### Epithelial cells from PIN ####
DimPlot(PIN2, reduction = "tsne", pt.size = 0.3, cols = c("blue", "blue", "blue", "blue", "blue", "blue", "grey", "grey", "blue", "grey", "blue", "grey", "grey", "grey", "grey", "grey", "blue"))

EpionlyPIN <- subset(PIN2, idents = c("0", "1", "2", "3", "4", "5", "8", "10", "16"))

# RE-clustering EpionlyPIN
EpionlyPIN <- ScaleData(EpionlyPIN, verbose = FALSE)
EpionlyPIN <- RunPCA(EpionlyPIN, npcs = 30, verbose = FALSE)
ElbowPlot(EpionlyPIN)

# t-SNE and Clustering
EpionlyPIN <- FindNeighbors(EpionlyPIN, reduction = "pca", dims = 1:10)
EpionlyPIN <- FindClusters(EpionlyPIN, resolution = 0.5)
EpionlyPIN <- RunTSNE(EpionlyPIN, reduction = "pca", dims = 1:10)
DimPlot(EpionlyPIN, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))

#FeaturePlot
FeaturePlot(EpionlyPIN, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(EpionlyPIN, reduction = "tsne", features = c("EGFP"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

#### EpionlyPIN Pseudotime ####

DefaultAssay(EpionlyPIN) <- "RNA"

PINEpiPseudo <- as.CellDataSet(EpionlyPIN)

PINEpiPseudo <- detectGenes(PINEpiPseudo, min_expr = 0.1)
print(head(fData(PINEpiPseudo)))
expressed_genes <- row.names(subset(fData(PINEpiPseudo),
                                    num_cells_expressed >= 10))

pData(PINEpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(PINEpiPseudo))
PINEpiPseudo <- PINEpiPseudo[,pData(PINPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(PINEpiPseudo), geom =
        "density")

PINEpiPseudo <- estimateSizeFactors(PINEpiPseudo)
PINEpiPseudo <- estimateDispersions(PINEpiPseudo)

disp_table <- dispersionTable(PINEpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
PINEpiPseudo <- setOrderingFilter(PINEpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(PINEpiPseudo)

# PINPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(PINEpiPseudo, return_all = F) # norm_method='log'
plot_pc_variance_explained(PINEpiPseudo, return_all = F)
PINEpiPseudo <- reduceDimension(PINEpiPseudo, max_components = 2, num_dim = 10,
                                reduction_method = 'tSNE', verbose = T)
PINEpiPseudo <- clusterCells(PINEpiPseudo, num_clusters = 2)

plot_cell_clusters(PINEpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(PINEpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
PINEpiPseudo <- setOrderingFilter(PINEpiPseudo, ordering_genes)
plot_ordering_genes(PINEpiPseudo)

PINEpiPseudo <- reduceDimension(PINEpiPseudo, max_components = 2,
                                method = 'DDRTree')

PINEpiPseudo <- orderCells(PINEpiPseudo)

GM_state <- function(PINEpiPseudo){
  if (length(unique(pData(PINEpiPseudo)$State)) > 1){
    T0_counts <- table(pData(PINEpiPseudo)$State, pData(PINEpiPseudo)$seurat_clusters)[,"1"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

PINEpiPseudo <- orderCells(PINEpiPseudo, root_state = GM_state(PINEpiPseudo))

#Cytotrace
EpionlyPINExpdata = data.frame(EpionlyPIN[["RNA"]]@data)
EpionlyPINphenotypedata = FetchData(EpionlyPIN, c("seurat_clusters"))
write.table(EpionlyPINExpdata, file = "EpionlyPINExpdata.txt")
write.table(EpionlyPINphenotypedata, file = "EpionlyPINphenotypedata.txt")

#Visualization

plot_cell_trajectory(PINEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

plot_cell_trajectory(PINEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "grey50"))

plot_cell_trajectory(PINEpiPseudo, markers = "ARQ", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(PINEpiPseudo, markers = "EGFP", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(PINEpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(PINEpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(PINEpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(PINEpiPseudo, markers = "Pbsn", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))


#### Basel cells in PIN ####
Idents(object = PIN2) <- "seurat_clusters"
DimPlot(PIN2, reduction = "tsne", pt.size = 0.3, cols = c("grey", "blue", "blue", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey", "grey"))

BEonlyPIN2 <- subset(PIN2, idents = c("1", "2"))

#Adding Expression Info
DefaultAssay(BEonlyPIN2) <- "RNA"
BEonlyPIN2ARQPos <- subset(x=BEonlyPIN2, subset = ARQ > 1)
BEonlyPIN2ARQNeg <- subset(x=BEonlyPIN2, subset = ARQ == 0)
Idents(object = BEonlyPIN2ARQPos) <- "ARQPos"
Idents(object = BEonlyPIN2ARQNeg) <- "ARQNeg"
BEonlyPIN2ARQPos[["ARQExp"]] <- Idents(object = BEonlyPIN2ARQPos)
BEonlyPIN2ARQNeg[["ARQExp"]] <- Idents(object = BEonlyPIN2ARQNeg)
BEonlyPIN2ARQ <- merge(x = BEonlyPIN2ARQPos, y = BEonlyPIN2ARQNeg)
Idents(object = BEonlyPIN2ARQ) <- "ARQExp"
BEonlyPIN2$ARQExp <- Idents(object = BEonlyPIN2ARQ)

Idents(object = BEonlyPIN2) <- "ARQExp"
ARQPosBEARQNegBEPIN2 <- subset(BEonlyPIN, idents = c("ARQPos", "ARQNeg"))
DimPlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", pt.size = 0.3, label = TRUE)

#Reclustering ARQ+ BE & ARQ- BE
ARQPosBEARQNegBEPIN2 <- ScaleData(ARQPosBEARQNegBEPIN2, verbose = FALSE)
ARQPosBEARQNegBEPIN2 <- RunPCA(ARQPosBEARQNegBEPIN2, npcs = 30, verbose = FALSE)
ElbowPlot(ARQPosBEARQNegBEPIN2)

#t-SNE and Clustering
ARQPosBEARQNegBEPIN2 <- FindNeighbors(ARQPosBEARQNegBEPIN2, reduction = "pca", dims = 1:10)
ARQPosBEARQNegBEPIN2 <- FindClusters(ARQPosBEARQNegBEPIN2, resolution = 0.5)
ARQPosBEARQNegBEPIN2 <- RunTSNE(ARQPosBEARQNegBEPIN2, reduction = "pca", dims = 1:10)
DimPlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF"))

#FeaturePlot
Idents(object = ARQPosBEARQNegBEPIN2) <- "RNA"
FeaturePlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", features = c("ARQ"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

FeaturePlot(ARQPosBEARQNegBEPIN2, reduction = "tsne", features = c("EGFP"), cols = c("light grey", "red"), pt.size = 0.3, max.cutoff = "q90")

#Cytotrace
ARQPosBEARQNegBEPIN2Expdata = data.frame(ARQPosBEARQNegBEPIN2[["RNA"]]@data)
ARQPosBEARQNegBEPIN2phenotypedata = FetchData(ARQPosBEARQNegBEPIN2, c("seurat_clusters"))
write.table(ARQPosBEARQNegBEPIN2Expdata, file = "ARQPosBEARQNegBEPIN2Expdata.txt")
write.table(ARQPosBEARQNegBEPIN2phenotypedata, file = "ARQPosBEARQNegBEPIN2phenotypedata.txt")

#Heatmap
Idents(object = ARQPosBEARQNegBEPIN2) <- "seurat_clusters"
ARQPosBE5ARQNegBE3PIN2 <- subset(ARQPosBEARQNegBEPIN2, idents = c("2", "4"))
DefaultAssay(ARQPosBE5ARQNegBE3PIN2) <- "RNA"
ARQPosBE5ARQNegBE3PIN2 <- ScaleData(ARQPosBE5ARQNegBE3PIN2, features = rownames(ARQPosBE5ARQNegBE3PIN2))
ARQPosBE5ARQNegBE3PIN2.markers <- FindAllMarkers(ARQPosBE5ARQNegBE3PIN2, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
ARQPosBE5ARQNegBE3PIN2PIN2Top50 <- ARQPosBE5ARQNegBE3PIN2.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)
DoHeatmap(ARQPosBE5ARQNegBE3PIN2, features = c(ARQPosBE5ARQNegBE3PIN2Top50$gene)) + NoLegend() + scale_fill_gradientn(colors = c("blue", "white", "red"))

#DEGs
ARQPosBE5ARQNegBE3PIN2.all.Markers <- FindMarkers(ARQPosBEARQNegBEPIN2, ident.1 = "4", ident.2 = "2")
write.csv(ARQPosBE5ARQNegBE3PIN2.all.Markers, "ARQPosBE5ARQNegBE3PIN2.all.Markers.csv")

#ViolinPlot
VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Igf1r", pt.size = 0)

VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Fos", pt.size = 0)

VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Mapk13", pt.size = 0)

VlnPlot(ARQPosBE5ARQNegBE3PIN2, features = "Jak2", pt.size = 0)

#### ARQ+ EpionlyPIN Pseudotime ####
Idents(object = PIN2) <- "seurat_clusters"
DefaultAssay(PIN2) <- "RNA"
FeaturePlot(PIN2, reduction = "tsne", features = c("Epcam"), cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")

#Adding Expression Info
DefaultAssay(PIN2) <- "RNA"
PINARQPos <- subset(x=PIN2, subset = ARQ > 1)
PINARQNeg <- subset(x=PIN2, subset = ARQ == 0)
Idents(object = PINARQPos) <- "ARQPos"
Idents(object = PINARQNeg) <- "ARQNeg"
PINARQPos[["ARQExp"]] <- Idents(object = PINARQPos)
PINARQNeg[["ARQExp"]] <- Idents(object = PINARQNeg)
PINARQ <- merge(x = PINARQPos, y = PINARQNeg)
Idents(object = PINARQ) <- "ARQExp"
PIN2$ARQExp <- Idents(object = PINARQ)

Idents(object = PIN2) <- "ARQExp"
ARQonlyPIN <- subset(PIN2, idents = c("ARQPos"))
DimPlot(ARQonlyPIN, reduction = "tsne", pt.size = 0.3) 

#RE-clustering ARQEpionlyPIN
DefaultAssay(ARQonlyPIN) <- "RNA"
ARQonlyPIN <- ScaleData(ARQonlyPIN, verbose = FALSE)
ARQonlyPIN <- RunPCA(ARQonlyPIN, npcs = 30, verbose = FALSE)
ElbowPlot(ARQonlyPIN)

# t-SNE and Clustering
ARQonlyPIN<- FindNeighbors(ARQonlyPIN, reduction = "pca", dims = 1:10)
ARQonlyPIN <- FindClusters(ARQonlyPIN, resolution = 0.5)
ARQonlyPIN <- RunTSNE(ARQonlyPIN, reduction = "pca", dims = 1:10)

ARQonlyPIN2 <- subset(ARQonlyPIN, idents = c("0", "1", "2", "3", "4", "5"))
DimPlot(ARQonlyPIN2, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF")) 

####ARQonlyPIN Pseudotime####
DefaultAssay(ARQonlyPIN2) <- "RNA"

ARQPINEpiPseudo <- as.CellDataSet(ARQonlyPIN2)

ARQPINEpiPseudo <- detectGenes(ARQPINEpiPseudo, min_expr = 0.1)
print(head(fData(ARQPINEpiPseudo)))
expressed_genes <- row.names(subset(fData(ARQPINEpiPseudo),
                                    num_cells_expressed >= 10))

pData(ARQPINEpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARQPINEpiPseudo))
ARQPINEpiPseudo <- ARQPINEpiPseudo[,pData(ARQPINEpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARQPINEpiPseudo), geom =
        "density")

ARQPINEpiPseudo <- estimateSizeFactors(ARQPINEpiPseudo)
ARQPINEpiPseudo <- estimateDispersions(ARQPINEpiPseudo)

disp_table <- dispersionTable(ARQPINEpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARQPINEpiPseudo <- setOrderingFilter(ARQPINEpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(ARQPINEpiPseudo)

# PINPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(ARQPINEpiPseudo, return_all = F) # norm_method='log'

ARQPINEpiPseudo <- reduceDimension(ARQPINEpiPseudo, max_components = 2, num_dim = 10,
                                   reduction_method = 'tSNE', verbose = T)
ARQPINEpiPseudo <- clusterCells(ARQPINEpiPseudo, num_clusters = 4)

plot_cell_clusters(ARQPINEpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(ARQPINEpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARQPINEpiPseudo <- setOrderingFilter(ARQPINEpiPseudo, ordering_genes)
plot_ordering_genes(ARQPINEpiPseudo)

ARQPINEpiPseudo <- reduceDimension(ARQPINEpiPseudo, max_components = 2,
                                   method = 'DDRTree')

ARQPINEpiPseudo <- orderCells(ARQPINEpiPseudo)

GM_state <- function(ARQPINEpiPseudo){
  if (length(unique(pData(ARQPINEpiPseudo)$State)) > 1){
    T0_counts <- table(pData(ARQPINEpiPseudo)$State, pData(ARQPINEpiPseudo)$seurat_clusters)[,"0"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARQPINEpiPseudo <- orderCells(ARQPINEpiPseudo, root_state = GM_state(ARQPINEpiPseudo))

#Visualization
plot_cell_trajectory(ARQPINEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

plot_cell_trajectory(ARQPINEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE)

####ARQPosPINvTumor Pseudotime####

###Subset ARQ Positive cells in integrated Epi###
Idents(object = PINvTumor.combined.Epi) <- "seurat_clusters"
DefaultAssay(PINvTumor.combined.Epi) <- "RNA"

#Adding Expression Info
EpiARQPos <- subset(x=PINvTumor.combined.Epi, subset = ARQ > 1)
EpiARQNeg <- subset(x=PINvTumor.combined.Epi, subset = ARQ == 0)
Idents(object = EpiARQPos) <- "ARQPos"
Idents(object = EpiARQNeg) <- "ARQNeg"
EpiARQPos[["ARQExp"]] <- Idents(object = EpiARQPos)
EpiARQNeg[["ARQExp"]] <- Idents(object = EpiARQNeg)
EpiARQ <- merge(x = EpiARQPos, y = EpiARQNeg)
Idents(object = EpiARQ) <- "ARQExp"
PINvTumor.combined.Epi$ARQExp <- Idents(object = EpiARQ)
Idents(object = PINvTumor.combined.Epi) <- "ARQExp"
ARQPosEpi <- subset(PINvTumor.combined.Epi, idents = c("ARQPos"))
DimPlot(ARQPosEpi, reduction = "tsne", pt.size = 0.3) + NoLegend()

#Recluster EpiARQPos
DefaultAssay(ARQPosEpi) <- "integrated"
ARQPosEpi <- ScaleData(ARQPosEpi, verbose = FALSE)
ARQPosEpi <- RunPCA(ARQPosEpi, npcs = 30, verbose = FALSE)
ElbowPlot(ARQPosEpi)

ARQPosEpi <- FindNeighbors(ARQPosEpi, reduction = "pca", dims = 1:10)
ARQPosEpi <- FindClusters(ARQPosEpi, resolution = 0.5)
ARQPosEpi <- RunTSNE(ARQPosEpi, reduction = "pca", dims = 1:10)
DimPlot(ARQPosEpi, reduction = "tsne", pt.size = 0.3, cols = c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF"))

DimPlot(ARQPosEpi, reduction = "tsne", group.by = "stim", cols = c("#B71800", "#009FB7"))

Idents(object = ARQPosEpi) <- "stim"
EpiARQPosPIN <- subset(ARQPosEpi, idents = c("PIN"))
DimPlot(EpiARQPosPIN, reduction = "tsne", cols = c("#B71800")) + NoLegend()

EpiARQPosTumor <- subset(EpiARQPos, idents = c("Tumor"))
DimPlot(EpiARQPosTumor, reduction = "tsne", cols = c("#009FB7")) + NoLegend()

#### ARQ Positive Epi Pseudotime ####
DefaultAssay(ARQPosEpi) <- "RNA"
ARQEpiPseudo <- as.CellDataSet(ARQPosEpi)

ARQEpiPseudo <- detectGenes(ARQEpiPseudo, min_expr = 0.1)
print(head(fData(ARQEpiPseudo)))
expressed_genes <- row.names(subset(fData(ARQEpiPseudo),
                                    num_cells_expressed >= 10))

pData(ARQEpiPseudo)$Total_mRNAs <- Matrix::colSums(exprs(ARQEpiPseudo))
ARQEpiPseudo <- ARQEpiPseudo[,pData(ARQEpiPseudo)$Total_mRNAs < 1e6]
qplot(Total_mRNAs, data = pData(ARQEpiPseudo), geom =
        "density")

ARQEpiPseudo <- estimateSizeFactors(ARQEpiPseudo)
ARQEpiPseudo <- estimateDispersions(ARQEpiPseudo)

disp_table <- dispersionTable(ARQEpiPseudo)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
ARQEpiPseudo <- setOrderingFilter(ARQEpiPseudo, unsup_clustering_genes$gene_id)
plot_ordering_genes(ARQEpiPseudo)

# ARQEpiPseudo@auxClusteringData[["tSNE"]]$variance_explained <- NULL
plot_pc_variance_explained(ARQEpiPseudo, return_all = F) 

ARQEpiPseudo <- reduceDimension(ARQEpiPseudo, max_components = 2, num_dim = 10,
                                reduction_method = 'tSNE', verbose = T)
ARQEpiPseudo <- clusterCells(ARQEpiPseudo, num_clusters = 2)
plot_cell_clusters(ARQEpiPseudo, color_by = "seurat_clusters")

diff_test_res <- differentialGeneTest(ARQEpiPseudo[expressed_genes,],
                                      fullModelFormulaStr = "~seurat_clusters")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
ARQEpiPseudo <- setOrderingFilter(ARQEpiPseudo, ordering_genes)
plot_ordering_genes(ARQEpiPseudo)

ARQEpiPseudo <- reduceDimension(ARQEpiPseudo, max_components = 2,
                                method = 'DDRTree')

ARQEpiPseudo <- orderCells(ARQEpiPseudo)
GM_state <- function(ARQEpiPseudo){
  if (length(unique(pData(ARQEpiPseudo)$State)) > 1){
    T0_counts <- table(pData(ARQEpiPseudo)$State, pData(ARQEpiPseudo)$seurat_clusters)[,"4"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}

ARQEpiPseudo <- orderCells(ARQEpiPseudo, root_state = GM_state(ARQEpiPseudo))

#Visualization
plot_cell_trajectory(ARQEpiPseudo, color_by = "Pseudotime", show_branch_points = FALSE)

plot_cell_trajectory(ARQEpiPseudo, color_by = "seurat_clusters", show_branch_points = FALSE) + scale_color_manual(breaks = c("X", "Y", "Z"), values=c("#E06666", "#FF9933", "#FFD966", "#A7EAA7", "#1D762E", "#3399FF", "#C7BEE1", "#761D65", "#C67E3D", "grey50", "#3333FF"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Krt5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Krt8", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Krt18", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Krt19", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Igf1r", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Fos", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Mapk13", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Jak2", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Axin2", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Tcf4", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Ccnd1", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))

plot_cell_trajectory(ARQEpiPseudo, markers = "Lgr5", use_color_gradient = TRUE, cell_size = 0.5, show_branch_points = FALSE) + scale_color_gradientn(colours = c("gray80", "red3"))
