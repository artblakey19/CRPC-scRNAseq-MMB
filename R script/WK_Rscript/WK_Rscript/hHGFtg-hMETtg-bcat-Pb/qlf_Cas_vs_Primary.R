
library(reshape2)
library(vegan)
library(rgl)
library(gplots)
library(grid)
library(gridExtra)
library(GenomicFeatures)
library(ggplot2)
library(statmod)
library(edgeR)
library(Clipper)

####Primary Vs WT (2 vs 3)####
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/1. PrimaryVsWT")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48627, countdata$COHP_48628, countdata$COHP_18716, countdata$COHP_28466, countdata$COHP_28467)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48627_Primary", "48628_Primary", "18716_WT", "28466_WT", "28467_WT")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("Primary", "Primary", "WT", "WT", "WT"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("Primary", "Primary", "WT", "WT", "WT"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48627_Primary", "48628_Primary", "18716_WT", "28466_WT", "28467_WT")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
normalized_RPKM[,5] <- normalized_RPKM[,5]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48627_Primary", "48628_Primary", "18716_WT", "28466_WT", "28467_WT")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan","pink", colors())
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
text(pc.values$PC1, pc.values$PC2, c("48627_Primary", "48628_Primary", "18716_WT", "28466_WT", "28467_WT"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
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
qlf_Primary_vs_WT <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Primary_vs_WT)

qlf_Primary_vs_WT$table$normRPKM_48627_Primary <- normalized_RPKM[row.names(qlf_Primary_vs_WT$table),1]
qlf_Primary_vs_WT$table$normRPKM_48628_Primary <- normalized_RPKM[row.names(qlf_Primary_vs_WT$table),2]
qlf_Primary_vs_WT$table$normRPKM_18716_WT <- normalized_RPKM[row.names(qlf_Primary_vs_WT$table),3]
qlf_Primary_vs_WT$table$normRPKM_28466_WT <- normalized_RPKM[row.names(qlf_Primary_vs_WT$table),4]
qlf_Primary_vs_WT$table$normRPKM_28467_WT <- normalized_RPKM[row.names(qlf_Primary_vs_WT$table),5]

qlf_Primary_vs_WT$table$RAW_48627_Primary <- w1_comp_10$counts[row.names(qlf_Primary_vs_WT$table),1]
qlf_Primary_vs_WT$table$RAW_48628_Primary <- w1_comp_10$counts[row.names(qlf_Primary_vs_WT$table),2]
qlf_Primary_vs_WT$table$RAW_18716_WT <- w1_comp_10$counts[row.names(qlf_Primary_vs_WT$table),3]
qlf_Primary_vs_WT$table$RAW_28466_WT <- w1_comp_10$counts[row.names(qlf_Primary_vs_WT$table),4]
qlf_Primary_vs_WT$table$RAW_28467_WT <- w1_comp_10$counts[row.names(qlf_Primary_vs_WT$table),5]

qlf_Primary_vs_WT_data_frame <- as.data.frame(topTags(qlf_Primary_vs_WT, n=30000))
sum(ifelse((qlf_Primary_vs_WT_data_frame$logFC >= 1)&(qlf_Primary_vs_WT_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Primary_vs_WT_data_frame$logFC <= -1)&(qlf_Primary_vs_WT_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Primary_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Primary_vs_WT_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Primary_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Primary_vs_WT_data_frame$PValue < 0.05),-1,0))

qlf_Primary_vs_WT_data_frame$UP <- ifelse((qlf_Primary_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Primary_vs_WT_data_frame$PValue < 0.05),1,0)
qlf_Primary_vs_WT_data_frame$DOWN <- ifelse((qlf_Primary_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Primary_vs_WT_data_frame$PValue < 0.05),-1,0)

sum(qlf_Primary_vs_WT_data_frame$UP)
sum(qlf_Primary_vs_WT_data_frame$DOWN)
qlf_Primary_vs_WT_data_frame$UP_OR_DOWN <- qlf_Primary_vs_WT_data_frame$UP + qlf_Primary_vs_WT_data_frame$DOWN
write.table(qlf_Primary_vs_WT_data_frame, file="qlf_Primary_vs_WT_data_frame_final.txt", quote=F, sep="\t", col.names=NA)

####Cas Vs WT (3 vs 3)####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/CasVsWT")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48624, countdata$COHP_48625, countdata$COHP_48626, countdata$COHP_18716, countdata$COHP_28466, countdata$COHP_28467)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48624_Castration", "48625_Castration", "48626_Castration", "18716_WT", "28466_WT", "28467_WT")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:6], group=c("Castration", "Castration", "Castration", "WT", "WT", "WT"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:6], group=c("Castration", "Castration", "Castration", "WT", "WT", "WT"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48624_Castration", "48625_Castration", "48626_Castration", "18716_WT", "28466_WT", "28467_WT")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
normalized_RPKM[,5] <- normalized_RPKM[,5]*1000/countdata$total_exon_length
normalized_RPKM[,6] <- normalized_RPKM[,6]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48624_Castration", "48625_Castration", "48626_Castration", "18716_WT", "28466_WT", "28467_WT")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan","pink", "red", colors())
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
text(pc.values$PC1, pc.values$PC2, c("48624_Castration", "48625_Castration", "48626_Castration", "18716_WT", "28466_WT", "28467_WT"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
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
qlf_Cas_vs_WT <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Cas_vs_WT)

qlf_Cas_vs_WT$table$normRPKM_48624_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_WT$table),1]
qlf_Cas_vs_WT$table$normRPKM_48625_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_WT$table),2]
qlf_Cas_vs_WT$table$normRPKM_48626_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_WT$table),3]
qlf_Cas_vs_WT$table$normRPKM_28466_WT <- normalized_RPKM[row.names(qlf_Cas_vs_WT$table),4]
qlf_Cas_vs_WT$table$normRPKM_28467_WT <- normalized_RPKM[row.names(qlf_Cas_vs_WT$table),5]
qlf_Cas_vs_WT$table$normRPKM_28467_WT <- normalized_RPKM[row.names(qlf_Cas_vs_WT$table),6]

qlf_Cas_vs_WT$table$RAW_48624_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_WT$table),1]
qlf_Cas_vs_WT$table$RAW_48625_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_WT$table),2]
qlf_Cas_vs_WT$table$RAW_48626_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_WT$table),3]
qlf_Cas_vs_WT$table$RAW_28466_WT <- w1_comp_10$counts[row.names(qlf_Cas_vs_WT$table),4]
qlf_Cas_vs_WT$table$RAW_28467_WT <- w1_comp_10$counts[row.names(qlf_Cas_vs_WT$table),5]
qlf_Cas_vs_WT$table$RAW_28467_WT <- w1_comp_10$counts[row.names(qlf_Cas_vs_WT$table),6]

qlf_Cas_vs_WT_data_frame <- as.data.frame(topTags(qlf_Cas_vs_WT, n=30000))
sum(ifelse((qlf_Cas_vs_WT_data_frame$logFC >= 1)&(qlf_Cas_vs_WT_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_WT_data_frame$logFC <= -1)&(qlf_Cas_vs_WT_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Cas_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_WT_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_WT_data_frame$PValue < 0.05),-1,0))

qlf_Cas_vs_WT_data_frame$UP <- ifelse((qlf_Cas_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_WT_data_frame$PValue < 0.05),1,0)
qlf_Cas_vs_WT_data_frame$DOWN <- ifelse((qlf_Cas_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_WT_data_frame$PValue < 0.05),-1,0)

sum(qlf_Cas_vs_WT_data_frame$UP)
sum(qlf_Cas_vs_WT_data_frame$DOWN)
qlf_Cas_vs_WT_data_frame$UP_OR_DOWN <- qlf_Cas_vs_WT_data_frame$UP + qlf_Cas_vs_WT_data_frame$DOWN
write.table(qlf_Cas_vs_WT_data_frame, file="qlf_Cas_vs_WT_data_frame_final.txt", quote=F, sep="\t", col.names=NA)

####Lung vs WT (2 vs 3)####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/LungVsWT")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48629, countdata$COHP_48630, countdata$COHP_18716, countdata$COHP_28466, countdata$COHP_28467)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48629_Lung", "48630_Lung", "18716_WT", "28466_WT", "28467_WT")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("Lung", "Lung", "WT", "WT", "WT"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("Lung", "Lung", "WT", "WT", "WT"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48629_Lung", "48630_Lung", "18716_WT", "28466_WT", "28467_WT")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
normalized_RPKM[,5] <- normalized_RPKM[,5]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48629_Lung", "48630_Lung", "18716_WT", "28466_WT", "28467_WT")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan","pink", colors())
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
text(pc.values$PC1, pc.values$PC2, c("48629_Lung", "48630_Lung", "18716_WT", "28466_WT", "28467_WT"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
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
qlf_Lung_vs_WT <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Lung_vs_WT)

qlf_Lung_vs_WT$table$normRPKM_48629_Lung <- normalized_RPKM[row.names(qlf_Lung_vs_WT$table),1]
qlf_Lung_vs_WT$table$normRPKM_48630_Lung <- normalized_RPKM[row.names(qlf_Lung_vs_WT$table),2]
qlf_Lung_vs_WT$table$normRPKM_18716_WT <- normalized_RPKM[row.names(qlf_Lung_vs_WT$table),3]
qlf_Lung_vs_WT$table$normRPKM_28466_WT <- normalized_RPKM[row.names(qlf_Lung_vs_WT$table),4]
qlf_Lung_vs_WT$table$normRPKM_28467_WT <- normalized_RPKM[row.names(qlf_Lung_vs_WT$table),5]

qlf_Lung_vs_WT$table$RAW_48629_Lung <- w1_comp_10$counts[row.names(qlf_Lung_vs_WT$table),1]
qlf_Lung_vs_WT$table$RAW_48630_Lung <- w1_comp_10$counts[row.names(qlf_Lung_vs_WT$table),2]
qlf_Lung_vs_WT$table$RAW_18716_WT <- w1_comp_10$counts[row.names(qlf_Lung_vs_WT$table),3]
qlf_Lung_vs_WT$table$RAW_28466_WT <- w1_comp_10$counts[row.names(qlf_Lung_vs_WT$table),4]
qlf_Lung_vs_WT$table$RAW_28467_WT <- w1_comp_10$counts[row.names(qlf_Lung_vs_WT$table),5]

qlf_Lung_vs_WT_data_frame <- as.data.frame(topTags(qlf_Lung_vs_WT, n=30000))
sum(ifelse((qlf_Lung_vs_WT_data_frame$logFC >= 1)&(qlf_Lung_vs_WT_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Lung_vs_WT_data_frame$logFC <= -1)&(qlf_Lung_vs_WT_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Lung_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Lung_vs_WT_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Lung_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Lung_vs_WT_data_frame$PValue < 0.05),-1,0))

qlf_Lung_vs_WT_data_frame$UP <- ifelse((qlf_Lung_vs_WT_data_frame$logFC >= log2(1.5))&(qlf_Lung_vs_WT_data_frame$PValue < 0.05),1,0)
qlf_Lung_vs_WT_data_frame$DOWN <- ifelse((qlf_Lung_vs_WT_data_frame$logFC <= -log2(1.5))&(qlf_Lung_vs_WT_data_frame$PValue < 0.05),-1,0)

sum(qlf_Lung_vs_WT_data_frame$UP)
sum(qlf_Lung_vs_WT_data_frame$DOWN)
qlf_Lung_vs_WT_data_frame$UP_OR_DOWN <- qlf_Lung_vs_WT_data_frame$UP + qlf_Lung_vs_WT_data_frame$DOWN
write.table(qlf_Lung_vs_WT_data_frame, file="qlf_Lung_vs_WT_data_frame_final.txt", quote=F, sep="\t", col.names=NA)


####Cas Vs Intact (3 vs 2)####

setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/CasVsPrimary")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48624, countdata$COHP_48625, countdata$COHP_48626, countdata$COHP_48627, countdata$COHP_48628)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48624_Castration", "48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("Castration", "Castration", "Castration", "Primary", "Primary"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48624_Castration", "48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
normalized_RPKM[,5] <- normalized_RPKM[,5]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48624_Castration", "48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan","pink", colors())
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
text(pc.values$PC1, pc.values$PC2, c("48624_Castration", "48625_Castration", "48626_Castration", "48627_Primary", "48628_Primary"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
dev.off()
groups

class(groups)
treatmentR3 <- factor(as.character(qc.grp))
treatmentR3
comp_10 <- relevel(treatmentR3, ref="Primary")
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
qlf_Cas_vs_Primary <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Cas_vs_Primary)

qlf_Cas_vs_Primary$table$normRPKM_48624_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),1]
qlf_Cas_vs_Primary$table$normRPKM_48625_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),2]
qlf_Cas_vs_Primary$table$normRPKM_48626_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),3]
qlf_Cas_vs_Primary$table$normRPKM_48627_Primary <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),4]
qlf_Cas_vs_Primary$table$normRPKM_48628_Primary <- normalized_RPKM[row.names(qlf_Cas_vs_Primary$table),5]


qlf_Cas_vs_Primary$table$RAW_48624_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),1]
qlf_Cas_vs_Primary$table$RAW_48625_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),2]
qlf_Cas_vs_Primary$table$RAW_48626_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),3]
qlf_Cas_vs_Primary$table$RAW_48627_Primary <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),4]
qlf_Cas_vs_Primary$table$RAW_48628_Primary <- w1_comp_10$counts[row.names(qlf_Cas_vs_Primary$table),5]

qlf_Cas_vs_Primary_data_frame <- as.data.frame(topTags(qlf_Cas_vs_Primary, n=30000))
sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC >= 1)&(qlf_Cas_vs_Primary_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC <= -1)&(qlf_Cas_vs_Primary_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_Primary_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),-1,0))

qlf_Cas_vs_Primary_data_frame$UP <- ifelse((qlf_Cas_vs_Primary_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),1,0)
qlf_Cas_vs_Primary_data_frame$DOWN <- ifelse((qlf_Cas_vs_Primary_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_Primary_data_frame$PValue < 0.05),-1,0)

sum(qlf_Cas_vs_Primary_data_frame$UP)
sum(qlf_Cas_vs_Primary_data_frame$DOWN)
qlf_Cas_vs_Primary_data_frame$UP_OR_DOWN <- qlf_Cas_vs_Primary_data_frame$UP + qlf_Cas_vs_Primary_data_frame$DOWN
write.table(qlf_Cas_vs_Primary_data_frame, file="qlf_Cas_vs_Primary_data_frame_final.txt", quote=F, sep="\t", col.names=NA)

#heatmap
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

#load dataframe
setwd("//isi-dcnl/user_data/zjsun/group/Lab Members/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/CasVsPrimary")
Cas_vs_Primary_Heatmap <- read.csv("Cas_vs_Primary_Heatmap.csv")

#create matrix file
mat <- as.matrix(Heatmap_CasVsIntact_4) 
my_colors = c("darkblue","white","darkred")
my_colors = colorRampPalette(my_colors)(50)
fpkm_log2 = log2(mat + 1)
pheatmap(fpkm_log2[,], scale = "row", color = my_colors,show_column_dend = FALSE, show_column_names = F,show_row_names = F)


#create matrix file
mat <- as.matrix(Heatmap_CasVIntact_combined) 
my_colors = c("darkblue","white","darkred")
my_colors = colorRampPalette(my_colors)(50)
fpkm_log2 = log2(mat + 1)
pheatmap(mat, scale = "row", color = my_colors,show_column_dend = FALSE, show_column_names = F,show_row_names = F)


Heatmap(mat)
htmap <- Heatmap(mat,show_column_dend = FALSE,show_column_names = F,show_row_names = F,
                 col = colorRamp2(c(-3,0,3), c("darkblue","white","darkred")),name = "Fold change",
                 width = unit (5,"cm"), height = unit(10,"cm"),resolution = 300,compression = "lzw",
                 heatmap_legend_param = list(legend_direction = "horizontal",
                                             legend_positin = "bottom",legend_width = unit(2, "cm"))) 



####Cas vs Lung (3 vs 2)####
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/CasVsLung")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48624, countdata$COHP_48625, countdata$COHP_48626, countdata$COHP_48629, countdata$COHP_48630)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48624_Castration", "48625_Castration", "48626_Castration", "48629_Lung", "48630_Lung")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:5], group=c("Castration", "Castration", "Castration", "Lung", "Lung"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48624_Castration", "48625_Castration", "48626_Castration", "48629_Lung", "48630_Lung")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
normalized_RPKM[,5] <- normalized_RPKM[,5]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48624_Castration", "48625_Castration", "48626_Castration", "48629_Lung", "48630_Lung")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan","pink", colors())
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
text(pc.values$PC1, pc.values$PC2, c("48624_Castration", "48625_Castration", "48626_Castration", "48629_Lung", "48630_Lung"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
dev.off()
groups

class(groups)
treatmentR3 <- factor(as.character(qc.grp))
treatmentR3
comp_10 <- relevel(treatmentR3, ref="Lung")
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
qlf_Cas_vs_Lung <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Cas_vs_Lung)

qlf_Cas_vs_Lung$table$normRPKM_48624_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Lung$table),1]
qlf_Cas_vs_Lung$table$normRPKM_48625_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Lung$table),2]
qlf_Cas_vs_Lung$table$normRPKM_48626_Castration <- normalized_RPKM[row.names(qlf_Cas_vs_Lung$table),3]
qlf_Cas_vs_Lung$table$normRPKM_48629_Lung <- normalized_RPKM[row.names(qlf_Cas_vs_Lung$table),4]
qlf_Cas_vs_Lung$table$normRPKM_48630_Lung <- normalized_RPKM[row.names(qlf_Cas_vs_Lung$table),5]


qlf_Cas_vs_Lung$table$RAW_48624_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Lung$table),1]
qlf_Cas_vs_Lung$table$RAW_48625_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Lung$table),2]
qlf_Cas_vs_Lung$table$RAW_48626_Castration <- w1_comp_10$counts[row.names(qlf_Cas_vs_Lung$table),3]
qlf_Cas_vs_Lung$table$RAW_48629_Lung <- w1_comp_10$counts[row.names(qlf_Cas_vs_Lung$table),4]
qlf_Cas_vs_Lung$table$RAW_48626_Lung <- w1_comp_10$counts[row.names(qlf_Cas_vs_Lung$table),5]

qlf_Cas_vs_Lung_data_frame <- as.data.frame(topTags(qlf_Cas_vs_Lung, n=30000))
sum(ifelse((qlf_Cas_vs_Lung_data_frame$logFC >= 1)&(qlf_Cas_vs_Lung_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_Lung_data_frame$logFC <= -1)&(qlf_Cas_vs_Lung_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Cas_vs_Lung_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_Lung_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Cas_vs_Lung_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_Lung_data_frame$PValue < 0.05),-1,0))

qlf_Cas_vs_Lung_data_frame$UP <- ifelse((qlf_Cas_vs_Lung_data_frame$logFC >= log2(1.5))&(qlf_Cas_vs_Lung_data_frame$PValue < 0.05),1,0)
qlf_Cas_vs_Lung_data_frame$DOWN <- ifelse((qlf_Cas_vs_Lung_data_frame$logFC <= -log2(1.5))&(qlf_Cas_vs_Lung_data_frame$PValue < 0.05),-1,0)

sum(qlf_Cas_vs_Lung_data_frame$UP)
sum(qlf_Cas_vs_Lung_data_frame$DOWN)
qlf_Cas_vs_Lung_data_frame$UP_OR_DOWN <- qlf_Cas_vs_Lung_data_frame$UP + qlf_Cas_vs_Lung_data_frame$DOWN
write.table(qlf_Cas_vs_Lung_data_frame, file="qlf_Cas_vs_Lung_data_frame_final.txt", quote=F, sep="\t", col.names=NA)

####Lung Vs Intact (2 vs 2)####

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/hHGFtg-hMETtg-b-cat-Pb/RNAseq/LungVsPrimary")
countdata <- read.csv("countdata.csv")

##CPM
wk_RAW <- cbind(countdata$COHP_48629, countdata$COHP_48630, countdata$COHP_48627, countdata$COHP_48628)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("48629_Lung", "48630_Lung", "48627_Primary", "48628_Primary")
rownames(wk_RAW) <- row.names(countdata)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:4], group=c("Lung", "Lung", "Primary", "Primary"), genes=data.frame(Symbol=countdata$symbol, Length=countdata$total_exon_length))
colnames(wk_RAW_DGEList) <- c("48629_Lung", "48630_Lung", "48627_Primary", "48628_Primary")
names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList, method = "TMM")
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w)

##TPM
normalized_RPKM <- CPM_w
normalized_RPKM[,1] <- normalized_RPKM[,1]*1000/countdata$total_exon_length
normalized_RPKM[,2] <- normalized_RPKM[,2]*1000/countdata$total_exon_length
normalized_RPKM[,3] <- normalized_RPKM[,3]*1000/countdata$total_exon_length
normalized_RPKM[,4] <- normalized_RPKM[,4]*1000/countdata$total_exon_length
keep_wk <- rowSums(normalized_RPKM >= 1) >= 1.4
sum(keep_wk)
w1 <- w[keep_wk, , keep.lib.sizes=FALSE]
w1$samples
colSums(w1$counts)
plotMDS(w1)
normalized_RPKM.keep <- normalized_RPKM[keep_wk,]
normalized_RPKM.keep.log2 <- log2(normalized_RPKM.keep+1)
write.table(normalized_RPKM.keep.log2, "normalized_RPKM.keep.log2.txt", quote=FALSE, sep="\t", col.names=NA)

##Further analysis
class(normalized_RPKM.keep.log2)
pca.values <- prcomp(na.omit(data.matrix(normalized_RPKM.keep.log2)))
pc.values <- data.frame(pca.values$rotation)
pc.values
variance.explained <- (pca.values $sdev)^2 / sum(pca.values $sdev^2)
pca.table <- data.frame(PC = 1:length(variance.explained), percent.variation = variance.explained, t(pc.values))
pca.text.file = paste("Group","_pca_values.txt",sep="")
write.table(pca.table, pca.text.file, quote=F, row.names=F, sep="\t")

qc.grp <- c("48629_Lung", "48630_Lung", "48627_Primary", "48628_Primary")
qc.grp
qc.grp <- as.character(w1$samples$group)
qc.grp
groups = levels(as.factor(as.character(qc.grp)))
groups
num.sample.types = length(groups)
num.sample.types
fixed.color.palatte = c("green","orange","purple","cyan", colors())
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
text(pc.values$PC1, pc.values$PC2, c("48629_Lung", "48630_Lung", "48627_Primary", "48628_Primary"), pos = 1)
legend("bottomright", legend=groups, col=color.palette, pch=19)
dev.off()
groups

class(groups)
treatmentR3 <- factor(as.character(qc.grp))
treatmentR3
comp_10 <- relevel(treatmentR3, ref="Primary")
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
qlf_Lung_vs_Primary <- glmQLFTest(fit_comp_10, coef=2)
topTags(qlf_Lung_vs_Primary)

qlf_Lung_vs_Primary$table$normRPKM_48629_Lung <- normalized_RPKM[row.names(qlf_Lung_vs_Primary$table),1]
qlf_Lung_vs_Primary$table$normRPKM_48630_Lung <- normalized_RPKM[row.names(qlf_Lung_vs_Primary$table),2]
qlf_Lung_vs_Primary$table$normRPKM_48627_Primary <- normalized_RPKM[row.names(qlf_Lung_vs_Primary$table),3]
qlf_Lung_vs_Primary$table$normRPKM_48628_Primary <- normalized_RPKM[row.names(qlf_Lung_vs_Primary$table),4]


qlf_Lung_vs_Primary$table$RAW_48629_Lung <- w1_comp_10$counts[row.names(qlf_Lung_vs_Primary$table),1]
qlf_Lung_vs_Primary$table$RAW_48630_Lung <- w1_comp_10$counts[row.names(qlf_Lung_vs_Primary$table),2]
qlf_Lung_vs_Primary$table$RAW_48627_Primary <- w1_comp_10$counts[row.names(qlf_Lung_vs_Primary$table),3]
qlf_Lung_vs_Primary$table$RAW_48628_Primary <- w1_comp_10$counts[row.names(qlf_Lung_vs_Primary$table),4]

qlf_Lung_vs_Primary_data_frame <- as.data.frame(topTags(qlf_Lung_vs_Primary, n=30000))
sum(ifelse((qlf_Lung_vs_Primary_data_frame$logFC >= 1)&(qlf_Lung_vs_Primary_data_frame$FDR < 0.05),1,0))
sum(ifelse((qlf_Lung_vs_Primary_data_frame$logFC <= -1)&(qlf_Lung_vs_Primary_data_frame$FDR < 0.05),-1,0))

sum(ifelse((qlf_Lung_vs_Primary_data_frame$logFC >= log2(1.5))&(qlf_Lung_vs_Primary_data_frame$PValue < 0.05),1,0))
sum(ifelse((qlf_Lung_vs_Primary_data_frame$logFC <= -log2(1.5))&(qlf_Lung_vs_Primary_data_frame$PValue < 0.05),-1,0))

qlf_Lung_vs_Primary_data_frame$UP <- ifelse((qlf_Lung_vs_Primary_data_frame$logFC >= log2(1.5))&(qlf_Lung_vs_Primary_data_frame$PValue < 0.05),1,0)
qlf_Lung_vs_Primary_data_frame$DOWN <- ifelse((qlf_Lung_vs_Primary_data_frame$logFC <= -log2(1.5))&(qlf_Lung_vs_Primary_data_frame$PValue < 0.05),-1,0)

sum(qlf_Lung_vs_Primary_data_frame$UP)
sum(qlf_Lung_vs_Primary_data_frame$DOWN)
qlf_Lung_vs_Primary_data_frame$UP_OR_DOWN <- qlf_Lung_vs_Primary_data_frame$UP + qlf_Lung_vs_Primary_data_frame$DOWN
write.table(qlf_Lung_vs_Primary_data_frame, file="qlf_Lung_vs_Primary_data_frame_final.txt", quote=F, sep="\t", col.names=NA)


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
