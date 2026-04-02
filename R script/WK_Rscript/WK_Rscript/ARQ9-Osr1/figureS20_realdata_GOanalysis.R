######## huaman monocyte #######
setwd('../data_DEG_analysis/')
library(DESeq2)
library(edgeR)
library(ChIPpeakAnno)
library(Clipper)
library(clusterProfiler)
library(org.Hs.eg.db)
library("AnnotationDbi")

mydeseq2 = function(count1, count2){
  
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  dds <- DESeqDataSetFromMatrix(dat, DataFrame(cond_idx), ~ cond_idx)
  dds <- DESeq(dds)
  count_norm = counts(dds,normalized = T, replaced = F)
  res <- results(dds, alpha = 0.05)
  discovery = rownames(res)[which(res$padj <= 0.05)]
  return(discovery)
}
myedger = function(count1, count2){
  n1 = ncol(count1)
  n2 = ncol(count2)
  cond_idx = rep(2, n1 + n2)
  cond_idx[1:n1] = 1
  dat = cbind(count1, count2)
  if(!is.factor(cond_idx)){
    cond_idx = factor(cond_idx)
  }
  y <- DGEList(counts=dat,group=cond_idx)
  keep <- filterByExpr(y)
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  count_norm = cpm(y)
  design <- model.matrix(~cond_idx)
  y <- estimateDisp(y,design)
  fit <- glmQLFit(y,design)
  qlf <- glmQLFTest(fit,coef=2)
  
  qlf <- topTags(qlf, n = nrow(dat), p.value = 0.05)
  discovery = rownames(qlf)
  return(discovery)
}
count1 = readRDS('count1.rds')
count2 = readRDS('count2.rds')

count1 = count1[,c(2,3,4)]
count2 = count2[,c(2,3,4)]
re_deseq2 = suppressMessages(mydeseq2(count1 = count1, count2 = count2))
re_edger = myedger(count1 = count1, count2 = count2)


dat = cbind(count1, count2)
cond_idx = rep(c(1,2), each = 3)
y <- DGEList(counts=dat,group=cond_idx)
y <- calcNormFactors(y)
count_norm_edger = cpm(y) #### normalized matrix
count1 = count_norm_edger[, (1:3)]
count2 = count_norm_edger[,-(1:3)]
re_clipper <- Clipper(log(count1+1), log(count2+1), 'd')
re_clipper <- rownames(count1)[re_clipper$discoveries[[1]]]

ensemblid <- matrix(unlist(strsplit(rownames(count1), '[.]')), nrow = 2)[1,]
allIDs = convert2EntrezID(IDs=ensemblid, orgAnn="org.Hs.eg.db",
                          ID_type="ensembl_gene_id")
geneID <- function(genelist, ont){
  ensemblid <- matrix(unlist(strsplit(genelist, '[.]')), nrow = 2)[1,]
  entrezIDs = convert2EntrezID(IDs=ensemblid, orgAnn="org.Hs.eg.db",
                               ID_type="ensembl_gene_id")
  return(entrezIDs)
}
re_clipper <- geneID(re_clipper)
re_edger <- geneID(re_edger)
re_deseq2 <- geneID(re_deseq2)


ego_edger <- enrichGO(re_edger, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
ego_deseq2 <- enrichGO(re_deseq2, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
ego_clipper <- enrichGO(re_clipper, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
sum(ego_edger[,7] <= 0.01)
sum(ego_clipper[,7] <= 0.01)
sum(ego_deseq2[,7] <= 0.01)
ego_edger_mf <- gene2GO(re_edger, "MF")
ego_deseq2_mf <- gene2GO(re_deseq2, "MF")
ego_clipper_mf <- gene2GO(re_clipper, "MF")
sum(ego_edger_mf[,7] <= 0.01)
nrow(ego_clipper_mf)
nrow(ego_deseq2_mf)

re_diff <- setdiff(re_clipper, re_edger)
ego_diff <- enrichGO(re_diff, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
sum(ego_diff[,7] <= 0.01)
re_diff_e <- setdiff(re_edger, re_clipper)
ego_diff_e <- enrichGO(re_diff_e, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
sum(ego_diff_e[,7] <= 0.01)
re_diff2 <- setdiff(re_clipper, re_deseq2)
ego_diff2 <- enrichGO(re_diff2, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
re_diff2_d <- setdiff(re_deseq2, re_clipper)
ego_diff2_d <- enrichGO(re_diff2_d, OrgDb  = org.Hs.eg.db, ont = "BP", universe = allIDs)
sum(ego_diff2[,7] <= 0.01)
sum(ego_diff2_d[,7] <= 0.01)



######## GO term ###########
pdf(file = paste0('GOterm_clipper_edger.pdf'), height = 4, width = 7)
tab_GO1 <- as.data.frame(ego_diff[1:8,c(1,2,7)])
tab_GO1$`GO term (ID)` <- paste0(tab_GO1$Description, ' (', tab_GO1$ID, ')')
tab_GO1$`qvalue (Clipper)` <- ego_clipper[match(tab_GO1$ID, ego_clipper[,1]),7]
p_GO <- ggtexttable(tab_GO1[,4:5], rows = NULL)
p_GO <- p_GO %>% tab_add_title(text = "GO terms enriched in Clipper-specific DEGs in Clipper vs. edgeR comparison", face = "plain", size = 11)
p_GO
dev.off()

pdf(file = paste0('GOterm_clipper_deseq2.pdf'), height = 4, width = 7)
tab_GO2 <- as.data.frame(ego_diff2[1:5,c(1,2,7)])
tab_GO2$`GO term (ID)` <- paste0(tab_GO2$Description, ' (', tab_GO2$ID, ')')
tab_GO2$`qvalue (Clipper)` <- ego_clipper[match(tab_GO2$ID, ego_clipper[,1]),7]
p_GO <- ggtexttable(tab_GO2[,4:5], rows = NULL)
p_GO <- p_GO %>% tab_add_title(text = "GO terms enriched in Clipper-specific DEGs in Clipper vs. DESeq2 comparison", face = "plain", size = 11)
p_GO
dev.off()


