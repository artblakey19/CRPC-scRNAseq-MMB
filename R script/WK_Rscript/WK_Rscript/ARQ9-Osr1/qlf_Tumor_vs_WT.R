setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/RNAseq/TumorvsWT")

####CPM####
wk_RAW <- cbind(genes_df$Sample_18716_WT_prostate, genes_df$Sample_28466_13_s8126, genes_df$Sample_28467_14_s7722, genes_df$Sample_31888_s7112, genes_df$Sample_31889_s7181)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("18716_WT_prostate", "28466_13_s8126", "28467_14_s7722", "31888_s7112", "31889_s7181")
rownames(wk_RAW) <- row.names(genes_df)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("WT", "WT", "WT", "Tumor", "Tumor"), genes=data.frame(Chr=as.character(genes_df$seqnames), Start=genes_df$start, End=genes_df$end, Strand=as.character(genes_df$strand), ID=genes_df$gene_id, Symbol=genes_df$symbol, Length=genes_df$total_exon_length))
colnames(wk_RAW_DGEList) <- c("18716_WT_prostate", "28466_13_s8126", "28467_14_s7722", "31888_s7112", "31889_s7181")
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

qc.grp <- c("18716", "28466", "28467", "31888", "31889")
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
text(pc.values$PC1, pc.values$PC2, c("18716", "28466", "28467", "31888", "31889"), pos = 1)
legend("topright", legend=groups, col=labelColors, pch=19)
dev.off()
groups

class(groups)
treatmentR3 <- factor(as.character(qc.grp))
treatmentR3
comp_10 <- relevel(treatmentR3, ref="WT")
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
qlf_Tumor_vs_WT <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Tumor_vs_WT)

qlf_Tumor_vs_WT$table$normRPKM_18716_WT_prostate <- normalized_RPKM[row.names(qlf_Tumor_vs_WT$table),1]
qlf_Tumor_vs_WT$table$normRPKM_28466_13_s8126 <- normalized_RPKM[row.names(qlf_Tumor_vs_WT$table),2]
qlf_Tumor_vs_WT$table$normRPKM_28467_14_s7722 <- normalized_RPKM[row.names(qlf_Tumor_vs_WT$table),3]
qlf_Tumor_vs_WT$table$normRPKM_31888_s7112 <- normalized_RPKM[row.names(qlf_Tumor_vs_WT$table),4]
qlf_Tumor_vs_WT$table$normRPKM_31889_s7181 <- normalized_RPKM[row.names(qlf_Tumor_vs_WT$table),5]

qlf_Tumor_vs_WT$table$RAW_18716_WT_prostate <- w1_comp_10$counts[row.names(qlf_Tumor_vs_WT$table),1]
qlf_Tumor_vs_WT$table$RAW_28466_13_s8126 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_WT$table),2]
qlf_Tumor_vs_WT$table$RAW_28467_14_s7722 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_WT$table),3]
qlf_Tumor_vs_WT$table$RAW_31888_s7112 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_WT$table),4]
qlf_Tumor_vs_WT$table$RAW_31889_s7181 <- w1_comp_10$counts[row.names(qlf_Tumor_vs_WT$table),5]

qlf_Tumor_vs_WT_data_frame <- as.data.frame(topTags(qlf_Tumor_vs_WT, n=30000))
sum(ifelse((qlf_Tumor_vs_WT_data_frame$logFC >= 1)&(qlf_Tumor_vs_WT_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Tumor_vs_WT_data_frame$logFC <= -1)&(qlf_Tumor_vs_WT_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Tumor_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Tumor_vs_WT_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Tumor_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Tumor_vs_WT_data_frame$PValue < 0.05),-1,0))

qlf_Tumor_vs_WT_data_frame$UP <- ifelse((qlf_Tumor_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Tumor_vs_WT_data_frame$PValue < 0.05),1,0)
qlf_Tumor_vs_WT_data_frame$DOWN <- ifelse((qlf_Tumor_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Tumor_vs_WT_data_frame$PValue < 0.05),-1,0)

sum(qlf_Tumor_vs_WT_data_frame$UP)
sum(qlf_Tumor_vs_WT_data_frame$DOWN)
qlf_Tumor_vs_WT_data_frame$UP_OR_DOWN <- qlf_Tumor_vs_WT_data_frame$UP + qlf_Tumor_vs_WT_data_frame$DOWN
write.table(qlf_Tumor_vs_WT_data_frame, file="qlf_Tumor_vs_WT_data_frame_final.txt", quote=F, sep="\t", col.names=NA)