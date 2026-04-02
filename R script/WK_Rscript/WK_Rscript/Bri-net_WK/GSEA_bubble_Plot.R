
# import package
library("ggplot2")
library(tidyverse)
library(reshape2)

##Try#Plot
ggplot(data=df_melted, aes(x=Location, y=value, color=Treatment)) +
  geom_point(position=position_dodge(width=0.3))

# plot: dot plot
data <-  as.data.frame(GSEA_bubble_Plot)
ggplot(data = data, aes(x = data$Group, y =`Enriched Pathways`, 
                        color =`p-value`, size = "NES")) + 
  geom_point() +
  scale_color_gradient(low = "green", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA enrichment analysis")

# Merged gsea bubble plot
GSEA_bubble_Plot_byGroup1 <- read_csv("//bri-net/Cancer Biology/ZJsun grp/WK/HGF-hMET_Bcat/scRNAseq/GSEA_bubble_Plot_byGroup1.csv")
data <-  as.data.frame(GSEA_bubble_Plot_byGroup1)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`FDR`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot Glandular Solid.tiff", width = 2.2, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = Group, y =`GSEA pathways`, 
                        color =`FDR`, size = `NES`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  ggtitle("GSEA")
dev.off()

# Glandular gsea bubble plot
GSEA_bubble_Plot_Glandular1 <- read_csv("//bri-net/Cancer Biology/ZJsun grp/WK/HGF-hMET_Bcat/scRNAseq/GSEA_bubble_Plot_Glandular1.csv")
data <-  as.data.frame(GSEA_bubble_Plot_Glandular1)
ggplot(data = data, aes(x = NES, y =`GSEA pathways`, 
                        color =`FDR`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(0.5,2.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot Glandular.tiff", width = 2.5, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = NES, y =`GSEA pathways`, 
                        color =`FDR`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(0.5,2.5)) +
  ggtitle("GSEA")
dev.off()

# Solid gsea bubble plot
GSEA_bubble_Plot_Solid1 <- read_csv("//bri-net/Cancer Biology/ZJsun grp/WK/HGF-hMET_Bcat/scRNAseq/GSEA_bubble_Plot_Solid1.csv")
data <-  as.data.frame(GSEA_bubble_Plot_Solid1)
ggplot(data = data, aes(x = NES, y =`GSEA pathways`, 
                        color =`FDR`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(0.5,2.5)) +
  ggtitle("GSEA")

tiff(file = "GSEA bubble plot Solid.tiff", width = 2.5, height = 4, units = "in", compression = "lzw", res = 800)
ggplot(data = data, aes(x = NES, y =`GSEA pathways`, 
                        color =`FDR`, size = `size`)) + 
  geom_point() +
  scale_color_gradient(low = "blue", high = "red") +
  theme_bw() + 
  ylab("") + 
  xlab("") + 
  scale_x_continuous(limits=c(0.5,2.5)) +
  ggtitle("GSEA")
dev.off()