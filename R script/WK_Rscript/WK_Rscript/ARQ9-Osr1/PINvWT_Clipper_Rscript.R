setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/RNAseq/PINvsWT/Clipper-1")

library(devtools)

install_github("JSB-UCLA/Clipper")
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
library(edgeR)
library(Clipper)

#Practice
dataurl = 'http://www.stat.ucla.edu/~jingyi.li/data/Clipper/'
count1 = readRDS(url(paste0(dataurl, 'simcount1.rds'))) 
count2 = readRDS(url(paste0(dataurl, 'simcount2.rds')))

r1 = ncol(count1)
r2 = ncol(count2)
cond_idx = rep(2, r1 + r2)
cond_idx[1:r1] = 1
dat = cbind(count1, count2)
cond_idx = factor(cond_idx)
y <- DGEList(counts=dat,group=cond_idx)
y <- calcNormFactors(y)
count_norm = cpm(y)

view(count_norm)
view(count_norm[,1:r1])
view(count_norm[,-(1:r1)])

re_clipper = Clipper(score.exp = log(base = 2, count_norm[,1:r1] + 1), 
                     score.back = log(base = 2, count_norm[,-(1:r1)] + 1), 
                     FDR = 0.5,
                     analysis = "differential")

####Method1####
wk_RAW1 <- cbind(genes_df$Sample_31885_s9364, genes_df$Sample_31886_s6562, genes_df$Sample_31892_s9369)
wk_RAW1 <- as.data.frame(wk_RAW1)
colnames(wk_RAW1) <- c("31885_s9364", "31886_s6562", "31892_s9369")
rownames(wk_RAW1) <- row.names(genes_df)

wk_RAW2 <- cbind(genes_df$Sample_18716_WT_prostate, genes_df$Sample_28466_13_s8126, genes_df$Sample_28467_14_s7722)
wk_RAW2 <- as.data.frame(wk_RAW2)
colnames(wk_RAW2) <- c("18716_WT_prostate", "28466_13_s8126", "28467_14_s7722")
rownames(wk_RAW2) <- row.names(genes_df)


r1 = ncol(wk_RAW1)
r2 = ncol(wk_RAW2)
cond_idx = rep(2, r1 + r2)
cond_idx[1:r1] = 1
dat = cbind(wk_RAW1, wk_RAW2)
cond_idx = factor(cond_idx)
y <- DGEList(counts=dat,group=cond_idx)
y <- calcNormFactors(y)
count_norm = cpm(y)

re_clipper = Clipper(score.exp = log(base = 2, count_norm[,1:r1] + 0.25), 
                     score.back = log(base = 2, count_norm[,-(1:r1)] + 0.25), 
                     FDR = 0.25,
                     analysis = "differential")

head(re_clipper$contrast.score.value)
write.csv(re_clipper$contrast.score.value, "PINvWT.clipper.contrast.score.value.csv")

head(re_clipper$q)
write.csv(re_clipper$q, "PINvWT.clipper.q-1.csv")

head(re_clipper$discoveries)
write.csv(re_clipper$discoveries, "PINvWT.clipper.discoveries.csv")

###Method2###
wk_RAW <- cbind(genes_df$Sample_28466_13_s8126, genes_df$Sample_28467_14_s7722, genes_df$Sample_31885_s9364, genes_df$Sample_31892_s9369)
wk_RAW <- as.data.frame(wk_RAW)
colnames(wk_RAW) <- c("28466_13_s8126", "28467_14_s7722", "31885_s9364", "31892_s9369")
rownames(wk_RAW) <- row.names(genes_df)
View(wk_RAW)

wk_RAW_DGEList <- DGEList(counts=wk_RAW[,1:4], group=c("WT", "WT", "PIN", "PIN"), genes=data.frame(Chr=as.character(genes_df$seqnames), Start=genes_df$start, End=genes_df$end, Strand=as.character(genes_df$strand), ID=genes_df$gene_id, Symbol=genes_df$symbol, Length=genes_df$total_exon_length))
colnames(wk_RAW_DGEList) <- c("28466_13_s8126", "28467_14_s7722", "31885_s9364", "31892_s9369")
write.csv(wk_RAW_DGEList, "wk_RAW_DGEList.csv")

names(wk_RAW_DGEList)
names(wk_RAW)
head(wk_RAW_DGEList$counts)

wk_RAW_DGEList$samples
head(wk_RAW_DGEList$genes)
w <- calcNormFactors(wk_RAW_DGEList)
w$samples
CPM_w <- cpm(w)
dim(CPM_w)
View(CPM_w[,3:4])
View(CPM_w[,1:2])

re_clipper = Clipper(score.exp = log(base = 2, CPM_w[,3:4] + 1), 
                     score.back = log(base = 2, CPM_w[,1:2] + 1), 
                     FDR = 0.1,
                     analysis = "differential")

head(re_clipper$contrast.score.value)
write.csv(re_clipper$contrast.score.value, "PINvWT.clipper.contrast.score.value.csv")

head(re_clipper$q)
write.csv(re_clipper$q, "PINvWT.clipper.q.csv")

head(re_clipper$discoveries)
write.csv(re_clipper$discoveries, "PINvWT.clipper.discoveries.csv")

