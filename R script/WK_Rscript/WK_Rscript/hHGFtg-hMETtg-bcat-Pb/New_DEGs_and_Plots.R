## Bulk RNA-seq analysis
#DEGS
BiocManager::install("DESeq2")
#Packages
library(tidyverse)
library(edgeR)

#Import Data
cts <- countdata_no_lung
cts <- dplyr::select(cts,!(total_exon_length))
cts <-column_to_rownames(cts, var = "symbol")
cts <- as.matrix(cts)
head(cts [1:5,])
#Groups
Groups <- c("Cas", "Cas", "Cas", "Pri", "Pri", "WT","WT","WT")
#dataset for 2vs2 comparison-----------------------------
cts2 <- countdata_no_lung
cts2 <- dplyr::select(cts2,!c(COHP_48624,COHP_28466,total_exon_length))
cts2 <-column_to_rownames(cts2, var = "symbol")
cts2 <- as.matrix(cts2)
head(cts2 [1:5,])
Groups2 <- c("Cas", "Cas","Pri", "Pri", "WT","WT")
##-------------------------------------------------
d <- DGEList(counts=cts,group=factor(Groups))
dim(d)
d.full <- d
head(d$counts)
head(cpm(d))
apply(d$counts, 2, sum)
keep <- rowSums(cpm(d)>100) >= 2
d <- d[keep,]
dim(d)
d <- calcNormFactors(d)
d
d$samples$lib.size <- colSums(d$counts)
d$samples
plotMDS(d, method="bcv", col=as.numeric(d$samples$group))
legend("bottomleft", as.character(unique(d$samples$group)), col=1:3, pch=20)
#Inter-library variation
d1 <- estimateCommonDisp(d, verbose=T)
names(d1)
d1 <- estimateTagwiseDisp(d1)
names(d1)
#GLM Estimates
design.mat <- model.matrix(~ 0 + d$samples$group)
colnames(design.mat) <- levels(d$samples$group)
d2 <- estimateGLMCommonDisp(d,design.mat)
d2 <- estimateGLMTrendedDisp(d2,design.mat, method="power")
d2 <- estimateGLMTagwiseDisp(d2,design.mat)
plotBCV(d2)
# Compare with DESeq2-------------------------------
library(DESeq2)
vignette("DESeq2")
cds <- DESeqDataSetFromMatrix(data.frame(d$counts),colData = Groups,design = d$samples$group)
cds <- estimateSizeFactors( cds )
sizeFactors( cds )
cds <- estimateDispersions( cds , method="blind")
##-----------------------------------------------------
plotDispEsts(cds)
#Differential Expression#Exact-test
et12 <- exactTest(d1, pair=c(1,2)) # compare groups Cas and Pri
et13 <- exactTest(d1, pair=c(1,3)) # compare groups Cas and WT
et14 <- exactTest(d1, pair=c(2,1)) # compare groups Pri and WT
et12
#Differential Expression#GLM
design.mat
fit <- glmFit(d2, design.mat)
lrt12 <- glmLRT(fit, contrast=c(1,-1,0)) # compare groups Cas and Pri
lrt13 <- glmLRT(fit, contrast=c(1,0,-1)) # compare groups Cas and wt
lrt23 <- glmLRT(fit, contrast=c(0,1,-1)) # compare groups pri and wt
topTags(lrt12,n=10)
##Set FDR 0.05/p-val
de_casvspri_3vs2 <- as.data.frame(topTags(lrt12, n=30000))
tst <- sum(ifelse((de_casvspri_3vs2$logFC >= 1)&(de_casvspri_3vs2$FDR < 0.05),1,0))
sum(ifelse((de_casvspri_3vs2$logFC <= -1)&(de_casvspri_3vs2$FDR < 0.05),-1,0))
#UP/Down
de_casvspri_3vs2$UP <- ifelse((de_casvspri_3vs2$logFC >= log2(1.5))&(de_casvspri_3vs2$PValue < 0.05),1,0)
de_casvspri_3vs2$DOWN <- ifelse((de_casvspri_3vs2$logFC <= -log2(1.5))&(de_casvspri_3vs2$PValue < 0.05),-1,0)
de_casvspri_3vs2$UP_OR_DOWN <- de_casvspri_3vs2$UP + de_casvspri_3vs2$DOWN
sum(de_casvspri_3vs2$UP)
sum(de_casvspri_3vs2$DOWN)
#Write
write.table(de_casvspri_3vs2,"cas_vs_pri_3vs2.txt")
write.table(de_casvspri_2vs2,"cas_vs_pri_2vs2.txt")
write.table(de_casvswt_3vs3,"cas_vs_wt_3vs3.txt")
write.table(de_casvswt_2vs3,"Primary_vs_wt_2vs3.txt")
write.table(de_casvswt_2vs2,"Primary_vs_wt_2vs2.txt")
#Mouse to human ortholog conversion
library(babelgene)
cas_vs_pri_3vs2_gsea_rank <- rename(cas_vs_pri_3vs2_gsea_rank, "symbol" = "gene")
genelist = cas_vs_pri_3vs2_gsea_rank$symbol
ortho_cas_pri_3vs2 <- orthologs(genes = genelist, species = "mouse", human = F, min_support = 5, top = TRUE)

casvspri_3vs2_hu_rank <- inner_join(cas_vs_pri_3vs2_gsea_rank,ortho_cas_pri_3vs2,by = "symbol") %>% select(contains (c("symbol","logFC")))

write.table(casvspri_3vs2_hu_rank,"cas_vs_pri_3vs2_hu_gsea_rank.txt")

df1 <- as.data.frame(lrt12)
print(et12)
et23 <-  as.data.frame(et23)
write.table(topTags(et12), "cas_vs_pri.txt")
topTags(et13)
topTags(et23)
et12
# combine
#remove
rm(et23)
####Plot####
#heatmap
library(ComplexHeatmap)
library(circlize)
library(pheatmap)

#create a dataframe of common genes
common_genes <- inner_join(Cas_vs_WT_3vs3,cas_pri_3vs2,by = "gene", suffix = c("_casvswt", "_casvspri")) %>%
  inner_join(., Primary_vs_WT_2vs3, by='gene',suffix = c(".,", "_privswt"))
common_genes <- unique(common_genes)
common_genes <- common_genes %>% select(starts_with (c("gene","logFC")))
common_genes <- rename(common_genes, logFC_privswt = logFC) %>% column_to_rownames( var = "gene")
#create matrix file
mat <- as.matrix(common_genes) 
pheatmap::pheatmap(common_genes, cluster_cols = F, cluster_rows = F)
rownames(mat) = common_genes$gene
colnames(mat) = colnames(common_genes)
mat <- log10(mat+10)

Heatmap(mat)
htmap <- Heatmap(mat,show_column_dend = FALSE,show_column_names = F,show_row_names = F,
        col = colorRamp2(c(-3,0,3), c("darkblue","white","darkred")),name = "Fold change",
        width = unit (5,"cm"), height = unit(10,"cm"),resolution = 300,compression = "lzw",
        heatmap_legend_param = list(legend_direction = "horizontal",
        legend_positin = "bottom",legend_width = unit(2, "cm"))) 


# library
#venn diagram
library(hrbrthemes)
library(tm)
library(proustr)
library(VennDiagram)
library(eulerr)

#Make the plot

x = list(Cas = cas_pri_3vs2$gene ,Pri = Cas_vs_WT_3vs3$gene,WT = Primary_vs_WT_2vs3$gene)

eulvenn <- plot(euler(x), shape = "ellipse",
     quantities = T,fintsize = 20,labels = F,imagetype="pdf",
     fill = "transparent",edges = c("red","green","purple"),
     height = 480 , width = 480 ,resolution = 300,compression = "lzw",
     lwd = 4,rotation = 10,cex = 0.5,cat.cex = 0.3,
cat.pos = c(-27, 27, 135),cat.dist = c(0.055, 0.055, 0.085))
eulvenn

#spiderplot
library(fmsb)
# To use the fmsb package, add 2 lines to the data frame: the max and min of each variable to show on the plot!
data <- rbind(rep(20,5) , rep(0,5) , data)

# Color vector
colors_border=c( rgb(0.2,0.5,0.5,0.9), rgb(0.8,0.2,0.5,0.9) , rgb(0.7,0.5,0.1,0.9) )
colors_in=c( rgb(0.2,0.5,0.5,0.4), rgb(0.8,0.2,0.5,0.4) , rgb(0.7,0.5,0.1,0.4) )

# plot with default options:
radarchart( data  , axistype=1 , 
            #custom polygon
            pcol=colors_border , pfcol=colors_in , plwd=4 , plty=1,
            #custom the grid
            cglcol="grey", cglty=1, axislabcol="grey", caxislabels=seq(0,20,5), cglwd=0.8,
            #custom labels
            vlcex=0.8 
)

# Add a legend
legend(x=0.7, y=1, legend = rownames(data[-c(1,2),]), bty = "n", pch=20 , col=colors_in , text.col = "grey", cex=1.2, pt.cex=3)