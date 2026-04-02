setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/RNAseq/PINvsWT/Clipper")

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

re_clipper = Clipper(score.exp = log(base = 2, count_norm[,1:r1] + 1), 
                     score.back = log(base = 2, count_norm[,-(1:r1)] + 1), 
                     FDR = 0.5,
                     analysis = "differential")

#My RNAseq data
wk_count2 <- cbind(genes_df$Sample_31888_s7112, genes_df$Sample_31889_s7181)
wk_count2 <- as.data.frame(wk_count2)
colnames(wk_count2) <- c("31888_s7112", "31889_s7181")
rownames(wk_count2) <- row.names(genes_df)
View(wk_count2)

wk_count1 <- cbind(genes_df$Sample_31885_s9364, genes_df$Sample_31886_s6562, genes_df$Sample_31892_s9369)
wk_count1 <- as.data.frame(wk_count1)
colnames(wk_count1) <- c("31885_s9364", "31886_s6562", "31892_s9369")
rownames(wk_count1) <- row.names(genes_df)
View(wk_count1)

r1 = ncol(wk_count1)
r2 = ncol(wk_count2)
cond_idx = rep(2, r1 + r2)
cond_idx[1:r1] = 1
dat = cbind(wk_count1, wk_count2)
cond_idx = factor(cond_idx)
y <- DGEList(counts=dat,group=cond_idx)
y <- calcNormFactors(y)
count_norm = cpm(y)

re_clipper = Clipper(score.exp = log(base = 2, count_norm[,1:r1] + 1), 
                     score.back = log(base = 2, count_norm[,-(1:r1)] + 1), 
                     FDR = 0.5,
                     analysis = "differential")

saveRDS(re_clipper, file = "TumorvWT_clipper.rds")
saveRDS(genes_df, file = "genes_df.rds")

