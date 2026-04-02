#Env
BiocManager::install("ComplexHeatmap")
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ggplotify)

df <- Heatmap_Triple_and_Double_LE_EpiCellTypes_4
df <- as.data.frame(df)
#Subset
df2 <- select(df, !(Fkbp5:Bmpr1b))
#
df <- column_to_rownames(df, var = "Clusters")
df <- as.matrix(df)
d <- pheatmap::pheatmap(df, cluster_cols = F, cluster_rows = F)

RColorBrewer::brewer.pal(n=3,name = "PuBlye")
df2 <- column_to_rownames(df2, var = "Clusters")
df2 <- as.matrix(df2)
plt2 <- pheatmap(df, cluster_cols = F,cluster_rows = F,
                   color = colorRampPalette(c("purple","grey", "yellow"))(101),
                 cellwidth = 20,fontsize = 10,angle_col = 45)
plt2
# Adding a pseudo-count of 10 to avoid strong color jumps with just 1 cell.
# Change to log10
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("purple","grey0", "yellow"))(101),
                   cluster_cols = F, cluster_rows = F)
pheatmap::pheatmap(log10(df+10),scale="column",
                   fontsize = 10,color = colorRampPalette(c("darkblue","white", "darkred"))(101),
                   cluster_cols = F, cluster_rows = F)


pheatmap(log10(df+10), color=colorRampPalette(c("darkblue","white", "darkred"))(101),
         cluster_cols = F, cluster_rows = F)


