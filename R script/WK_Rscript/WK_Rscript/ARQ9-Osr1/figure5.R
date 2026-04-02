######## huaman monocyte #######
setwd('~/Dropbox/clipper/DE analysis')
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

######## veen diagram #########
library(ggVennDiagram)
venn_list <- list(Clipper = re_clipper,
                  edgeR = re_edger,
                  DESeq2 = re_deseq2)
saveRDS(venn_list, "GOanalysis/genelist.rds")
venn_list <- readRDS("GOanalysis/genelist.rds")
print(sapply(venn_list, length))

pdf(file = paste0('Veen.pdf'), height = 4, width = 4)
venn_slingshot <- ggVennDiagram(venn_list, label = "count", label_alpha = 0) + scale_fill_gradient(low="white", high = "#2ca25f")
venn_slingshot
dev.off()
######## bar plot ############
theme_set(theme_bw() +
            theme(
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank()))

GOnum <- data.frame(pair = c("Total", "Total","Total", "Clipper vs edgeR", "Clipper vs edgeR", "Clipper vs DESeq2", "Clipper vs DESeq2"), 
                    method = c("Clipper", "edgeR", "DESeq2", "Clipper", "edgeR", "Clipper", "DESeq2"), 
                    num = c(110, 87, 22, 8,0,5,0))
color_vec <- c("#619BC3", "#A75E4C", "#D8A64E", "#898B11", "#BBBD46")
pdf(file = paste0('Barplot_GO.pdf'), height = 4, width = 4)
p_GOnum <- ggplot(data = GOnum, aes(x=pair, y=num, fill=method)) +
  geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values=color_vec[1:3]) + ylab("Number of significant GO terms") +
  theme(axis.title.x=element_blank(), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(angle = 0, size = 9, colour = "black")
        #      axis.text.x=element_blank(),        axis.ticks.x=element_blank()
  ) + theme(aspect.ratio = 1, legend.position = c(0.2, 0.5), legend.title = element_blank()) + geom_text(aes(label= num), vjust=1.4, color="black", size=3.5, position=position_dodge(0.9))

p_GOnum
dev.off()

######## GO term ###########
tab_GO <- as.data.frame(ego_deseq2[1:4,c(1,2,7)])
tab_GO$`GO term (ID)` <- paste0(tab_GO$Description, ' (', tab_GO$ID, ')')
tab_GO$`qvalue (DESeq2)` <- tab_GO$qvalue
tab_GO$`qvalue (edgeR)` <- ego_edger[match(tab_GO$ID, ego_edger[,1]),7]
tab_GO$`qvalue (Clipper)` <- ego_clipper[match(tab_GO$ID, ego_clipper[,1]),7]


library(ggpubr)
pdf(file = paste0('GOterm.pdf'), height = 4, width = 8)
p_GO <- ggtexttable(tab_GO[4:7], rows = NULL)
p_GO <- p_GO %>% tab_add_title(text = "Top GO terms enriched in all three sets of identified DEGs", face = "plain", size = 12)
p_GO
dev.off()
