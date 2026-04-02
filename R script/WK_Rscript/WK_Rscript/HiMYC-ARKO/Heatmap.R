#Heatmap_Comaprison_human_data

library(tidyr)
library(dplyr)
library(ggplot2)
library(ggrepel)
library(gplots)
library(ggdendro)
library(RColorBrewer)

setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/Human data")

#DEGs
DefaultAssay(CtrlvARKO.combined.MYCBE) <- "RNA"
Idents(object = CtrlvARKO.combined.MYCBE) <- "stim"
CtrlvARKO.combined.MYCBE <- ScaleData(CtrlvARKO.combined.MYCBE, features = rownames(CtrlvARKO.combined.MYCBE))

MYCBE_ARKOvCtrl.marker <- FindMarkers(CtrlvARKO.combined.MYCBE, ident.1 = "ARKO", ident.2 = "Ctrl", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(MYCBE_ARKOvCtrl.marker, "MYCBE_ARKOvCtrl.marker.csv")

MYCBE_CtrlvARKO.marker <- FindMarkers(CtrlvARKO.combined.MYCBE, ident.1 = "Ctrl", ident.2 = "ARKO", min.pct = 0.1, logfc.threshold = 0.1)
write.csv(MYCBE_CtrlvARKO.marker, "MYCBE_CtrlvARKO.marker.csv")

#Import data from files
setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/HiMYC-ARKO/Final/Human data")

BE_vs_LE <- as.data.frame(read.csv("Tang_et_al_BEvLE_DEGs.csv"))
ARKO_vs_Ctrl <- as.data.frame(read.csv("MYCBE_ARKOvCtrl.csv"))
Ctrl_vs_ARKO <- as.data.frame(read.csv("MYCBE_CtrlvARKO.csv"))

#Make gene table with adjusted P-values and Log2-FoldChanges
make_volcano_tab <- function(df, ixsymbol, ixpv, ixfc, cutpv, cutfc) {
  ng <- dim(df)[1]
  colnames(df)[ixsymbol] <- "Symbol"
  colnames(df)[ixpv] <- "pval"
  colnames(df)[ixfc] <- "log2FC"
  df <- df %>% select(Symbol, pval, log2FC)
  df %>% mutate(degflag = as.numeric(pval<0.05 & abs(log2FC)>=0.1) * sign(log2FC),
                significance = -log10(pval)) -> df
  df$significance[is.infinite(df$significance)] <- 300
  df$significance[df$significance>300] <- 300
  top10symbol <- df %>% filter(degflag==1) %>% slice_min(pval, n=10) %>% select(Symbol)
  bottom10symbol <- df %>% filter(degflag==-1) %>% slice_min(pval, n=10) %>% select(Symbol)
  top10ix <- which(df$Symbol %in% unlist(top10symbol))
  bottom10ix <- which(df$Symbol %in% unlist(bottom10symbol))
  df$top10 <- ""
  df$top10[top10ix] <- df$Symbol[top10ix]
  df$top10[bottom10ix] <- df$Symbol[bottom10ix]
  return(df)
}

voltab_ARKO_vs_Ctrl <- make_volcano_tab(ARKO_vs_Ctrl, 1, 2, 3, 0.05, 0.1)
voltab_Ctrl_vs_ARKO <- make_volcano_tab(Ctrl_vs_ARKO, 1, 2, 3, 0.05, 0.1)
voltab_BE_vs_LE <- make_volcano_tab(BE_vs_LE, 1, 7, 6, 0.05, 0.1)

#Make foldchange matrix for heatmap
fcmat <- merge(voltab_ARKO_vs_Ctrl[,c(1,3)], voltab_Ctrl_vs_ARKO[,c(1,3)], 
               by.x = "Symbol", by.y = "Symbol", all.x = T) 
fcmat <- merge(fcmat, voltab_BE_vs_LE[,c(1,3)], by.x = "Symbol", by.y = "Symbol", all.x = T)
fcmat <- na.omit(fcmat)
colnames(fcmat) <- c("Symbol","ARKO vs Ctrl", "Ctrl vs ARKO", "BE vs LE")

# subset of foldchange matrix with DEGs
degs <- unique(c(voltab_ARKO_vs_Ctrl$Symbol[voltab_ARKO_vs_Ctrl$degflag != 0],
                 voltab_Ctrl_vs_ARKO$Symbol[voltab_Ctrl_vs_ARKO$degflag != 0],
                 voltab_BE_vs_LE$Symbol[voltab_BE_vs_LE$degflag != 0]))

degmat <- as.matrix(fcmat[fcmat$Symbol %in% degs,2:4])
rownames(degmat) <- fcmat$Symbol[fcmat$Symbol %in% degs]
write.csv(degmat, "degmat.csv")

#Draw Heatmap using qplot package
mypal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(11)
mypal[c(4:8)] <- "#F7F7F7"
heatmap.2(degmat, col=mypal, hclustfun = function(x) hclust(x,method = 'ward.D2'),
          breaks = seq(-1,1,length.out=12), Colv = NA, dendrogram = "row",
          colsep = c(1,2), margins = c(3,10), srtCol = 0, adjCol = c(0.5,0.5),
          offsetRow = 0, offsetCol = 0, key.xlab = "Log2 Fold Change", key.title = "",
          density.info = "none", trace = "none", cexCol = 1, labRow = NA)

#Draw Heatmap using certain genes using qplot package
df_heatmap <- as.data.frame(read.csv("17genes.csv"))
mypal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(11)
mypal[c(5:7)] <- "#F7F7F7"
df_mat <- as.matrix(df_heatmap[,2:4])
colnames(df_mat) <- c("ARKO vs Ctrl","Ctrl vs ARKO","BE vs LE")
heatmap.2(df_mat, col=mypal, distfun = function(x) dist(x, method ='euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          breaks = seq(-1.5,1.5,length.out=12), Colv = NA, dendrogram = "none",
          colsep = c(1:8), margins = c(3,6), srtCol = 45, adjCol = c(0.5,0.5),
          offsetRow = 0, offsetCol = 1, key.xlab = "Log2 Fold Change", key.title = "",
          density.info = "none", trace = "none", cexCol = 1, labRow = df_heatmap[,1])

#Draw Heatmap using certain genes using qplot package
df_heatmap1 <- as.data.frame(read.csv("common_genes.csv"))
mypal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(11)
mypal[c(5:7)] <- "#F7F7F7"
df_mat <- as.matrix(df_heatmap1[,2:4])
colnames(df_mat) <- c("ARKO vs Ctrl","Ctrl vs ARKO","BE vs LE")
heatmap.2(df_mat, col=mypal, distfun = function(x) dist(x, method ='euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          breaks = seq(-2,2,length.out=12), Colv = NA, dendrogram = "none",
          colsep = c(1:8), margins = c(3,6), srtCol = 45, adjCol = c(0.5,0.5),
          offsetRow = 0, offsetCol = 1, key.xlab = "Log2 Fold Change", key.title = "",
          density.info = "none", trace = "none", cexCol = 1, labRow = df_heatmap1[,1])

#Draw Heatmap using certain genes using qplot package
df_heatmap1 <- as.data.frame(read.csv("common_genes.csv"))
mypal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(11)
mypal[c(5:7)] <- "#F7F7F7"
df_mat <- as.matrix(cbind(-df_heatmap1[,4:4], df_heatmap1[,4:4], df_heatmap1[2:3]))
colnames(df_mat) <- c("LE vs BE","BE vs LE", "ARKO vs Ctrl", "Ctrl vs ARKO")

tiff(file = "human Heatmap common -1.5_1.5.tiff", width = 5, height = 15, units = "in", compression = "lzw", res = 200)
heatmap.2(df_mat, col=mypal, distfun = function(x) dist(x, method ='euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          breaks = seq(-1.5,1.5,length.out=12), Colv = NA, dendrogram = "none",
          colsep = c(1:8), margins = c(3,6), srtCol = 45, adjCol = c(0.5,0.5),
          offsetRow = 0, offsetCol = 1, key.xlab = "Log2 Fold Change", key.title = "",
          density.info = "none", trace = "none", cexCol = 1, labRow = df_heatmap1[,1])
dev.off()

#Draw Heatmap using certain genes using qplot package
df_heatmap <- as.data.frame(read.csv("16genes.csv"))
mypal <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(11)
mypal[c(5:7)] <- "#F7F7F7"
df_mat <- as.matrix(cbind(-df_heatmap[,4:4], df_heatmap[,4:4], df_heatmap[2:3]))
colnames(df_mat) <- c("LE vs BE","BE vs LE", "ARKO vs Ctrl", "Ctrl vs ARKO")

tiff(file = "human Heatmap intersting common -1.5_1.5.tiff", width = 5, height = 15, units = "in", compression = "lzw", res = 200)
heatmap.2(df_mat, col=mypal, distfun = function(x) dist(x, method ='euclidean'),
          hclustfun = function(x) hclust(x,method = 'ward.D2'),
          breaks = seq(-1.5,1.5,length.out=12), Colv = NA, dendrogram = "none",
          colsep = c(1:8), margins = c(3,6), srtCol = 45, adjCol = c(0.5,0.5),
          offsetRow = 0, offsetCol = 1, key.xlab = "Log2 Fold Change", key.title = "",
          density.info = "none", trace = "none", cexCol = 1, labRow = df_heatmap[,1])
dev.off()
