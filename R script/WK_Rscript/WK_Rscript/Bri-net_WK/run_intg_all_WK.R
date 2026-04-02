library(Seurat)
#library(devtools)
library(dplyr)
library(data.table)
library(Matrix)
library(ggplot2)
library(patchwork)


options(future.globals.maxSize = 8000 * 1024^2)
base_dir="/net/isi-dcnl/ifs/user_data/zjsun/group/Sun_Lab_Sequencing_Data/"
out_dir = paste0(base_dir,"anal_hjcho/intg_All/")
setwd(out_dir)
print(out_dir)

#Cas-WT
data_info1=paste0(base_dir,"Gli1-ARKO_Castration/30044/filtered_feature_bc_matrix")
sc.data1 <- Read10X(data.dir =data_info1)
sc.obj1 <- CreateSeuratObject(counts = sc.data1, project = "INTG", min.cells = 3, min.features = 200)
sc.obj1[["percent.mt"]] <- PercentageFeatureSet(sc.obj1, pattern = "^mt-")
group1 <- sample("CAS", size = length(colnames(sc.obj1)), replace = TRUE)
names(group1) <- colnames(sc.obj1)
sc.obj1 <- AddMetaData(object = sc.obj1, metadata = group1, col.name = "groups")
type1 <- sample("WT", size = length(colnames(sc.obj1)), replace = TRUE)
names(type1) <- colnames(sc.obj1)
sc.obj1 <- AddMetaData(object = sc.obj1, metadata = type1, col.name = "types")


#Cas-Knockout
data_info2=paste0(base_dir,"Gli1-ARKO_Castration/30045/filtered_feature_bc_matrix")

sc.data2 <- Read10X(data.dir =data_info2)
sc.obj2 <- CreateSeuratObject(counts = sc.data2, project = "INTG", min.cells = 3, min.features = 200)
sc.obj2[["percent.mt"]] <- PercentageFeatureSet(sc.obj2, pattern = "^mt-")
group2 <- sample("CAS", size = length(colnames(sc.obj2)), replace = TRUE)
names(group2) <- colnames(sc.obj2)
sc.obj2 <- AddMetaData(object = sc.obj2, metadata = group2, col.name = "groups")
type2 <- sample("KO", size = length(colnames(sc.obj2)), replace = TRUE)
names(type2) <- colnames(sc.obj2)
sc.obj2 <- AddMetaData(object = sc.obj2, metadata = type2, col.name = "types")

#Regen-WT
data_info3=paste0(base_dir,"Gli1-ARKO_Regen/29979/filtered_feature_bc_matrix")
sc.data3 <- Read10X(data.dir =data_info3)
sc.obj3 <- CreateSeuratObject(counts = sc.data3, project = "INTG", min.cells = 3, min.features = 200)
sc.obj3[["percent.mt"]] <- PercentageFeatureSet(sc.obj3, pattern = "^mt-")
group3 <- sample("Regen", size = length(colnames(sc.obj3)), replace = TRUE)
names(group3) <- colnames(sc.obj3)
sc.obj3 <- AddMetaData(object = sc.obj3, metadata = group3, col.name = "groups")
type3 <- sample("WT", size = length(colnames(sc.obj3)), replace = TRUE)
names(type3) <- colnames(sc.obj3)
sc.obj3 <- AddMetaData(object = sc.obj3, metadata = type3, col.name = "types")

#Cas-Knockout
data_info4=paste0(base_dir,"Gli1-ARKO_Regen/29980/filtered_feature_bc_matrix")

sc.data4 <- Read10X(data.dir =data_info4)
sc.obj4 <- CreateSeuratObject(counts = sc.data4, project = "INTG", min.cells = 3, min.features = 200)
sc.obj4[["percent.mt"]] <- PercentageFeatureSet(sc.obj4, pattern = "^mt-")
group4 <- sample("Regen", size = length(colnames(sc.obj4)), replace = TRUE)
names(group4) <- colnames(sc.obj4)
sc.obj4 <- AddMetaData(object = sc.obj4, metadata = group4, col.name = "groups")
type4 <- sample("KO", size = length(colnames(sc.obj4)), replace = TRUE)
names(type4) <- colnames(sc.obj4)
sc.obj4 <- AddMetaData(object = sc.obj4, metadata = type4, col.name = "types")


sc.list = c(sc.obj1, sc.obj2, sc.obj3, sc.obj4)
names(sc.list) = c('Cas-WT', 'Cas-KO', 'Regen-WT','Regen-KO')
for (i in 1:length(sc.list)) {
   sc.list[[i]] <- SCTransform(sc.list[[i]], verbose = FALSE)
}
sc.features <- SelectIntegrationFeatures(object.list = sc.list, nfeatures = 3000)
sc.list <- PrepSCTIntegration(object.list = sc.list, anchor.features = sc.features,
    verbose = FALSE)

#identify anchors
sc.anchors <- FindIntegrationAnchors(object.list = sc.list, normalization.method = "SCT", anchor.features = sc.features, verbose = FALSE)
sc.raw <- IntegrateData(anchorset = sc.anchors, normalization.method = "SCT", verbose = FALSE)
#get plot befor filtration
print("before filter")
print(dim(sc.raw))
png("vlnplot_before.png", width=800, height=700)
VlnPlot(sc.raw, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()

plot1 <- FeatureScatter(sc.raw, feature1 = "nCount_RNA", feature2 = "percent.mt") + geom_hline(yintercept=15, linetype="dashed", color = "red") + annotate("text", x = 340000, y = 15, label = "15%", vjust = -0.5, color="red", size=5)
plot2 <- FeatureScatter(sc.raw, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_hline(yintercept=750, linetype="dashed", color = "red") + annotate("text", x = 340000, y = 750, label = "750", vjust = -0.5, color="red", size=5)+geom_hline(yintercept=8500, linetype="dashed", color = "red") + annotate("text", x = 340000, y = 8500, label = "8500", vjust = -0.5, color="red", size=5)
png("featureScatter.png", width=1000, height=1000)
plot1 + plot2
dev.off()

#filter later after integrate

sc.obj <-subset(x=sc.raw, subset = nFeature_RNA > 750 & nFeature_RNA < 8500 & percent.mt < 15)

png("vlnplot_after.png", width=800, height=700)
VlnPlot(sc.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
dev.off()
print("after filter")
print(dim(sc.obj))


DefaultAssay(sc.obj) <- "RNA"
all.genes <- rownames(sc.obj)
DefaultAssay(sc.obj) <- "integrated"
sc.obj <- ScaleData(sc.obj, features = all.genes)
sc.obj <- RunPCA(sc.obj, features = VariableFeatures(object = sc.obj))
png("elbowPlot.png", width=800, height=800)
ElbowPlot(sc.obj, ndims = 40)
dev.off()

#need to find dimention from elbow plot
# Determine percent of variation associated with each PC
pct <- sc.obj[["pca"]]@stdev / sum(sc.obj[["pca"]]@stdev) * 100

# Calculate cumulative percents for each PC
cumu <- cumsum(pct)

# Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
co1 <- which(cumu > 90 & pct < 5)[1]

co1
# Determine the difference between variation of PC and subsequent PC
co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1

# Minimum of the two calculation
pcs <- min(co1, co2)


#cluster
sc.obj <- FindNeighbors(sc.obj, dims = 1:pcs)
sc.obj <- FindClusters(sc.obj, resolution = 0.5) #19 -> #18 clusters changed

DefaultAssay(sc.obj) <- "SCT"

sc.keep <- sc.obj

sc.obj <- RunUMAP(sc.keep, reduction = "pca", dims = 1:30)

save.image()

png("dimPlot.png", width=800, height=800)
DimPlot(sc.obj, reduction = "pca")
dev.off()

png("dimplot_split.png", width=2000, height=700)
DimPlot(sc.obj, reduction = "umap", label = TRUE, split.by="groups")
dev.off()


#Find markers
deg <- FindAllMarkers(sc.obj, min.pct = 0.25, logfc.threshold = 0.25)
#deg.sorted <- deg[order(abs(deg$avg_logFC)), ]
write.table(deg, file=paste0("all_markers_",out_nm,".txt"), sep="\t",row.names=FALSE,quote=FALSE)

#get volcano to compare
library(EnhancedVolcano)
DefaultAssay(sc.obj) <- "SCT"
m <- FindMarkers(sc.obj, ident.1 = grp_nm1, ident.2 = grp_nm2, group.by = 'groups', subset.ident = 0, logfc.threshold = 0)  #all genes
m <- m[!is.infinite(m$avg_logFC),]
m.sorted <- m[order(abs(m$avg_logFC)), ]

png(paste0("volcano_cluster0_scale.png"), width=800, height=800)
EnhancedVolcano(m.sorted, lab = rownames(m.sorted), x = 'avg_logFC', y = 'p_val')
dev.off()
   
DefaultAssay(sc.obj) <- "RNA"
mb <- FindMarkers(sc.obj, ident.1 = grp_nm1, ident.2 = grp_nm2, group.by = 'groups', subset.ident = 0, logfc.threshold = 0)   #all genes
mb <- mb[!is.infinite(mb$avg_logFC),]
mb.sorted <- mb[order(abs(mb$avg_logFC)), ]
png(paste0("volcano_cluster0_no_scale.png"), width=800, height=800)
#png("volcano_before_scale.png", width=800, height=800)
EnhancedVolcano(mb.sorted, lab = rownames(mb.sorted), x = 'avg_logFC', y = 'p_val')
dev.off()

#Cell counts (split)
Idents(object = sc.obj) <- "groups"
sc.obj$groups.types <- paste(Idents(sc.obj), sc.obj$types, sep = "_")
Idents(object = sc.obj) <- "groups.types"
sc.obj$groups.types.seurat_clusters <- paste(Idents(sc.obj), sc.obj$seurat_clusters, sep = "_")
Idents(object = sc.obj) <- "groups.types.seurat_clusters"
table(Idents(sc.obj))

#FeaturePlots (split)

setwd("W:/ARKO-Gli1_Cas_scRNAseq/intg_All")

DefaultAssay(sc.obj) <- "RNA"
Idents(object = sc.obj) <- "seurat_clusters"

tiff(file = "intg_All Ar Exp.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(sc.obj, reduction = "umap", features = c("Ar"), split.by = 'groups.types', cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "intg_All EGFP Exp.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(sc.obj, reduction = "umap", features = c("EGFP"), split.by = 'groups.types', cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()
tiff(file = "intg_All Gli1 Exp.tiff", width = 16, height = 4, units = "in", compression = "lzw", res = 800)
FeaturePlot(sc.obj, reduction = "umap", features = c("Gli1"), split.by = 'groups.types', cols = c("light grey", "purple"), pt.size = 0.3, max.cutoff = "q90")
dev.off()




