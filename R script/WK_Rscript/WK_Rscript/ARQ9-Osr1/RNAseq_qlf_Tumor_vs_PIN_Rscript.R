install.packages('reshape2')
install.packages('vegan')
install.packages('rgl')
install.packages('gplots')
install.packages('grid')
install.packages('gridExtra')
BiocManager::install("GenomicFeatures")

library(reshape2)
library(vegan)
library(rgl)
library(gplots)
library(grid)
library(gridExtra)
library(GenomicFeatures)
library(ggplot2)
library(statmod)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/RNAseq")

####CPM####
wk_RAW <- cbind(genes_df$Sample_31885_s9364, genes_df$Sample_31886_s6562, genes_df$Sample_31892_s9369, genes_df$Sample_31888_s7112, genes_df$Sample_31889_s7181)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("31885_s9364", "31886_s6562", "31892_s9369", "31888_s7112", "31889_s7181")
rownames(wk_RAW) <- row.names(genes_df)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("PIN", "PIN", "PIN", "Tumor", "Tumor"), genes=data.frame(Chr=as.character(genes_df$seqnames), Start=genes_df$start, End=genes_df$end, Strand=as.character(genes_df$strand), ID=genes_df$gene_id, Symbol=genes_df$symbol, Length=genes_df$total_exon_length))
colnames(wk_RAW_DGEList) <- c("31885_s9364", "31886_s6562", "31892_s9369", "31888_s7112", "31889_s7181")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "upperquartile")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

####TPM####
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/genes_df$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/genes_df$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/genes_df$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/genes_df$total_exon_length
normalized_RPKM[,5] <- normalized_RPKM[,5]*1000/genes_df$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

####Further analysis####
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("31885", "31886", "31892", "31888", "31889")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan","pink",colors())
class(fixed.color.palatte)
length(fixed.color.palatte)
color.palette <- fixed.color.palatte[1:num.sample.types]
color.palette
labelColors = rep("black",times=ncol(normalized_RPKM.keep.log2))
labelColors
for (i in 1:num.sample.types){
  labelColors[qc.grp == as.character(groups[i])] = color.palette[i]
}
labelColors
pca.file = paste("pca_by_","Group",".png",sep="")
png(file=pca.file)

plot(pc.values$PC1, pc.values$PC2, col = labelColors, xlab = paste("PC1 (",round(100* variance.explained[1] , digits = 2),"%)", sep = ""),ylab = paste("PC2 (",round(100* variance.explained[2] , digits = 2),"%)", sep = ""), pch=19)
text(pc.values$PC1, pc.values$PC2, c("31885", "31886", "31892", "31888", "31889"), pos = 1)
legend("topright", legend=groups, col=labelColors, pch=19)
dev.off()
groups

class(groups)
treatmentR3 <- factor(as.character(qc.grp))
treatmentR3
comp_10 <- relevel(treatmentR3, ref="PIN")
comp_10
comp_10_design <- model.matrix(~comp_10)
comp_10_design
rownames(comp_10_design) <- colnames(w1)
comp_10_design
w1_comp_10 <- estimateDisp(w1, comp_10_design, robust = TRUE)
w1_comp_10$common.dispersion
plotBCV(w1_comp_10)
fit_comp_10 <- glmQLFit(w1_comp_10, comp_10_design, robust=TRUE)
plotQLDisp(fit_comp_10)
names(fit_comp_10)
head(fit_comp_10$coefficients)

fit_comp_10$design
qlf_Tumor_vs_PIN <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Tumor_vs_PIN)

qlf_Tumor_vs_PIN$table$normRPKM_31885_s9364 <- normalized_RPKM[row.names(qlf_Tumor_vs_PIN$table),1]
qlf_Tumor_vs_PIN$table$normRPKM_31886_s6562 <- normalized_RPKM[row.names(qlf_Tumor_vs_PIN$table),2]
qlf_Tumor_vs_PIN$table$normRPKM_31892_s9369 <- normalized_RPKM[row.names(qlf_Tumor_vs_PIN$table),3]
qlf_Tumor_vs_PIN$table$normRPKM_31888_s7112 <- normalized_RPKM[row.names(qlf_Tumor_vs_PIN$table),4]
qlf_Tumor_vs_PIN$table$normRPKM_31889_s7181 <- normalized_RPKM[row.names(qlf_Tumor_vs_PIN$table),5]

qlf_Tumor_vs_PIN$table$RAW_31885_s9364 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_PIN$table),1]
qlf_Tumor_vs_PIN$table$RAW_31886_s6562 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_PIN$table),2]
qlf_Tumor_vs_PIN$table$RAW_31892_s9369 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_PIN$table),3]
qlf_Tumor_vs_PIN$table$RAW_31888_s7112 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_PIN$table),4]
qlf_Tumor_vs_PIN$table$RAW_31889_s7181 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_PIN$table),5]

qlf_Tumor_vs_PIN_data_frame <- as.data.frame(topTags(qlf_Tumor_vs_PIN, n=30000))
sum(ifelse((qlf_Tumor_vs_PIN_data_frame$logFC >= 1)&(qlf_Tumor_vs_PIN_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Tumor_vs_PIN_data_frame$logFC <= -1)&(qlf_Tumor_vs_PIN_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Tumor_vs_PIN_data_frame$logFC >= log2(1.5))&(qlf_Tumor_vs_PIN_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Tumor_vs_PIN_data_frame$logFC <= -log2(1.5))&(qlf_Tumor_vs_PIN_data_frame$PValue < 0.05),-1,0))

qlf_Tumor_vs_PIN_data_frame$UP <- ifelse((qlf_Tumor_vs_PIN_data_frame$logFC >= log2(1.5))&(qlf_Tumor_vs_PIN_data_frame$PValue < 0.05),1,0)
qlf_Tumor_vs_PIN_data_frame$DOWN <- ifelse((qlf_Tumor_vs_PIN_data_frame$logFC <= -log2(1.5))&(qlf_Tumor_vs_PIN_data_frame$PValue < 0.05),-1,0)

sum(qlf_Tumor_vs_PIN_data_frame$UP)
sum(qlf_Tumor_vs_PIN_data_frame$DOWN)
qlf_Tumor_vs_PIN_data_frame$UP_OR_DOWN <- qlf_Tumor_vs_PIN_data_frame$UP + qlf_Tumor_vs_PIN_data_frame$DOWN
write.table(qlf_Tumor_vs_PIN_data_frame, file="qlf_Tumor_vs_PIN_data_frame_final.txt", quote=F, sep="\t", col.names=NA)

####Vendiagram####
venn.plot <- draw.triple.venn(area1 = 267, area2 = 1571, area3 = 948,
                              n12 = 169, n23= 948, n13 = 53, n123 = 53)
grid.draw(venn.plot);
grid.newpage();

library(Vennerable)
Vdemo <- Venn(SetNames = c("A", "B", "C"), Weight = c('000' = 0, '100' = 171, '010' = 122, '110' = 0, '001' = 1080, '101' = 135, '011' = 1702, '111' = 63))

tiff(file = "RNAseq Venn Diagram.tiff", width = 5, height = 5, units = "in", compression = "lzw", res = 800)
plot(Vdemo, doWeights = TRUE, type = "circles")
dev.off()
###RV-Script_FPKM_to_DGE####
#Create log2+1 normalized fpkm for ARPC/DNPC
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.16")
library(scran)
#_________________________________________________
#Subset
mARDNPC_ABI <- subset(mCRPC_ABI, idents = c("ARPC", "DNPC"))
table(Idents(mCRPC_ABI))
##extract count matrix
mARDNPC_count <- mARDNPC_ABI@assays$RNA@counts
a <- as.matrix(mARDNPC_count)
#Log Transform
b=log2(a/10+1)
c=log2(a+1)

#Check Gen Expression profile
b[1:5,1:5]
b['MET', ]
c['MET', ]
mARDNPC_count_log2 <- as.data.frame(c)
#Metadata
head(mARDNPC_ABI@meta.data)
mARDNPC_ABI_meta<- data.frame(mARDNPC_ABI@meta.data)
mARDNPC_ABI_meta<- rownames_to_column(mARDNPC_ABI_meta)
mARDNPC_ABI_meta <- select(mARDNPC_ABI_meta,rowname,orig.ident,Treatment,ARNE,
                           Treatment.ARNE,compare.ARNE,S.Score,G2M.Score, Phase)
mARDNPC_ABI_meta <- column_to_rownames(mARDNPC_ABI_meta, var = "rowname")
#create seurat object
mARDNPC_ABI =  CreateSeuratObject(counts = c, project = "mARDNPC")
#add metadata
mARDNPC_ABI <- AddMetaData(mARDNPC_ABI,mARDNPC_ABI_meta)
Idents(object = mARDNPC_ABI) <- "Treatment.ARNE"
table(Idents(mARDNPC_ABI))
mARDNPC_ABI <- FindVariableFeatures(mARDNPC_ABI, selection.method = "vst", nfeatures = 6000)
saveRDS(mARDNPC_ABI,"mARDNPC_ABI.rds")
#Intigrate#Fast MNN..........
library(SeuratWrappers)
mARDNPC_ABI <- RunFastMNN(object.list = SplitObject( mARDNPC_ABI,  split.by = "Treatment.ARNE"), nfeatures = 6000)
mARDNPC_ABI
#Cluster
mARDNPC_ABI <- FindVariableFeatures(mARDNPC_ABI, selection.method = "vst", nfeatures = 6000)
all.genes <- rownames(mARDNPC_ABI)
mARDNPC_ABI <- ScaleData(mARDNPC_ABI, features = all.genes)
mARDNPC_ABI<- RunPCA(int_obj, features = VariableFeatures(object = int_obj))
ElbowPlot(int_obj)
mARDNPC_ABI <- RunUMAP(int_obj, reduction = "mnn", dims = 1:30)
mARDNPC_ABI <- FindNeighbors(mARDNPC_ABI, reduction = "mnn", dims = 1:30)
imARDNPC_ABI <- FindClusters(mARDNPC_ABI)
#DEGs
ARPC_vs_DNPC_markers <- FindMarkers(mCRPC_ABI, ident.1 = "ARPC", ident.2 = "DNPC",
                                    +             test.use = "wilcox",min.pct = 0.1, only.pos = F)