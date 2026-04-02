#Running SCENIC on HPC

library(Seurat)
library(devtools)
library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(RColorBrewer)
library(SeuratData)
library(BiocManager)
library(SeuratWrappers)
library(monocle3)
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/Cell Type Identification/E18.5")

load("E18_Ctrl1.RData")

Idents(object = E18_Ctrl1) <- "seurat_clusters"
DimPlot(E18_Ctrl1, reduction = "umap", pt.size = 1)
cellInfo <- data.frame(Class=Idents(E18_Ctrl1))
exprMat <- as.matrix(GetAssayData(E18_Ctrl1, slot = "counts"))
cellInfo$nGene <- colSums(exprMat>0)
head(cellInfo)

cellInfo <- data.frame(cellInfo)
cellTypeColumn <- "Class"
colnames(cellInfo)[which(colnames(cellInfo)==cellTypeColumn)] <- "seurat_clusters"
cbind(table(cellInfo$seurat_clusters))

dir.create("int")
saveRDS(cellInfo, file="int/cellInfo.Rds")

colVars <- list(seurat_clusters=c("0"="#CCCCCC", 
                           "1"="#666666",
                           "2"="#000000","3"="#FFCC00","4"="#FF9900","5"="#FF6600","6"="#FF3300",
                           "7"="#33FF00","8"="#339900","9"="#336600","10"="#CCFF00","11"="#00FF66",
                           "12"="#996600","13"="#660000","14"="#FF0000","15"="#CC0033","16"="#FF3399",
                           "17"="#FF9999","18"="#FFCCFF","19"="#99FFFF","20"="#33CCCC","21"="#006666",
                           "22"="#0000FF","23"="#00CCFF","24"="#9900CC","25"="#6666FF","26"="#3300CC"
                           
))

colVars$seurat_clusters <- colVars$seurat_clusters[intersect(names(colVars$seurat_clusters), cellInfo$seurat_clusters)]
saveRDS(colVars, file="int/colVars.Rds")
plot.new(); legend(0,1, fill=colVars$seurat_clusters, legend=names(colVars$seurat_clusters))

org="mgi"
dbDir="cisTarget_databases" # RcisTarget databases location
myDatasetTitle="E18_Ctrl1 SCENIC RUN" # choose a name for your analysis
data(defaultDbNames)
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, datasetTitle=myDatasetTitle, nCores=10) 
scenicOptions@inputDatasetInfo$cellInfo <- "int/cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "int/colVars.Rds"
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions, minCountsPerGene=3*.01*ncol(exprMat), minSamples=ncol(exprMat)*.01)

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)

logMat <- exprMat # Better if it is logged/normalized

save.image("E18_Ctrl1.RData")

runCorrelation(exprMat_filtered, scenicOptions)

runGenie3(exprMat_filtered, scenicOptions)

scenicOptions <- readRDS("int/scenicOptions.Rds")
scenicOptions@settings$verbose <- TRUE
scenicOptions@settings$nCores <- 10
scenicOptions@settings$seed <- 123

scenicOptions@settings$dbs <- scenicOptions@settings$dbs["10kb"]

runSCENIC_1_coexNetwork2modules(scenicOptions)
runSCENIC_2_createRegulons(scenicOptions)
runSCENIC_3_scoreCells(scenicOptions, logMat)

nPcs <- c(5,15,50)

scenicOptions@settings$seed <- 123 # same seed for all of them
# Run t-SNE with different settings:
tSNE1 <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50))
tSNE2 <- tsneAUC(scenicOptions, aucType="AUC", nPcs=nPcs, perpl=c(5,15,50), onlyHighConf=TRUE, filePrefix="int/tSNE_oHC")
# Plot as pdf (individual files in int/):
tSNE3 <- paste0("int/",grep(".Rds", grep("tSNE_", list.files("int"), value=T), value=T))

par(mfrow=c(length(nPcs), 3))
tSNE4 <- paste0("int/",grep(".Rds", grep("tSNE_AUC", list.files("int"), value=T, perl = T), value=T))
plotTsne_compareSettings(tSNE4, scenicOptions, showLegend=FALSE, varName="seurat_clusters", cex=.5)

scenicOptions@settings$defaultTsne$aucType <- "AUC"
scenicOptions@settings$defaultTsne$dims <- 5
scenicOptions@settings$defaultTsne$perpl <- 15
saveRDS(scenicOptions, file="int/scenicOptions.Rds")

save.image("E18_Ctrl1.RData")

load("E18_Ctrl1.RData")
aucellApp <- plotTsne_AUCellApp(scenicOptions, logMat) # default t-SNE
tSNE_scenic <- readRDS(tsneFileName(scenicOptions))
aucell_regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
par(mfrow=c(2,3))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, exprMat, aucell_regulonAUC[onlyNonDuplicatedExtended(rownames(aucell_regulonAUC))[c("Ar", "Foxl1", "Foxf1","Gli1", "Gli2")],], plots="Expression")
Cairo::CairoPDF("output/Step4_BinaryRegulonActivity_tSNE_colByAUC.pdf", width=20, height=15)
par(mfrow=c(4,6))
AUCell::AUCell_plotTSNE(tSNE_scenic$Y, cellsAUC=aucell_regulonAUC, plots="AUC")
dev.off()
regulons <- loadInt(scenicOptions, "regulons")
regulons <- loadInt(scenicOptions, "aucell_regulons")
head(cbind(onlyNonDuplicatedExtended(names(regulons))))
regulonAUC <- loadInt(scenicOptions, "aucell_regulonAUC")
regulonAUC <- regulonAUC[onlyNonDuplicatedExtended(rownames(regulonAUC)),]
regulonActivity_byseurat_clusters <- sapply(split(rownames(cellInfo), cellInfo$seurat_clusters),
                                     function(cells) rowMeans(getAUC(regulonAUC)[,cells]))
regulonActivity_byseurat_clusters_Scaled <- t(scale(t(regulonActivity_byseurat_clusters), center = T, scale=T))
write.csv(regulonActivity_byseurat_clusters_Scaled, file = "E18_Ctrl1_TFs_SCENIC.csv")
pheatmap::pheatmap(regulonActivity_byCseurat_clusters_Scaled, #fontsize_row=3, 
                   color=colorRampPalette(c("blue","white","red"))(100), breaks=seq(0, 1, length.out = 100),
                   treeheight_row=10, treeheight_col=10, border_color=NA)

topRegulators <- reshape2::melt(regulonActivity_byseurat_clusters_Scaled)
colnames(topRegulators) <- c("Regulon", "seurat_clusters", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators, file = "E18_Ctrl1_topRegulators.csv")
save.image("E18_Ctrl1.RData")

#not scaled...
topRegulators <- reshape2::melt(regulonActivity_byseurat_clusters)
colnames(topRegulators) <- c("Regulon", "seurat_clusters", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>0),]
viewTable(topRegulators)
write.csv(topRegulators, file = "E18_Ctrl1_topRegulators.csv")

#Binarized version
minPerc <- .7
binaryRegulonActivity <- loadInt(scenicOptions, "aucell_binary_nonDupl")
cellInfo_binarizedCells <- cellInfo[which(rownames(cellInfo)%in% colnames(binaryRegulonActivity)),, drop=FALSE]
regulonActivity_byseurat_clusters_Binarized <- sapply(split(rownames(cellInfo_binarizedCells), cellInfo_binarizedCells$seurat_clusters), 
                                               function(cells) rowMeans(binaryRegulonActivity[,cells, drop=FALSE]))
binaryActPerc_subset <- regulonActivity_byseurat_clusters_Binarized[which(rowSums(regulonActivity_byseurat_clusters_Binarized>minPerc)>0),]
ComplexHeatmap::Heatmap(binaryActPerc_subset, name="Regulon activity (%)", col = c("white","pink","red"))
topRegulators <- reshape2::melt(regulonActivity_byseurat_clusters_Binarized)
colnames(topRegulators) <- c("Regulon", "seurat_clusters", "RelativeActivity")
topRegulators <- topRegulators[which(topRegulators$RelativeActivity>minPerc),]
viewTable(topRegulators)