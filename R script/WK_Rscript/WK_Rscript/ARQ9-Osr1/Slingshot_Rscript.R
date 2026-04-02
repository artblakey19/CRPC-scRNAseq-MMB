#Slingshot

library(slingshot)
library(BUSpaRse)
library(tidyverse)
library(tidymodels)
library(Seurat)
library(scales)
library(viridis)
library(Matrix)

DefaultAssay(BEcombined) <- "RNA"
Idents(object = BEcombined) <- "seurat_clusters"

BE_sds <- slingshot(Embeddings(BEcombined, "umap"), clusterLabels = BEcombined$seurat_clusters, 
                    start.clus = 6, stretch = 2)

cell_pal <- function(cell_vars, pal_fun,...) {
  if (is.numeric(cell_vars)) {
    pal <- pal_fun(100, ...)
    return(pal[cut(cell_vars, breaks = 100)])
  } else {
    categories <- sort(unique(cell_vars))
    pal <- setNames(pal_fun(length(categories), ...), categories)
    return(pal[cell_vars])
  }
}

cell_colors <- cell_pal(BEcombined$cell_type, brewer_pal("qual", "Set2"))
cell_colors_clust <- cell_pal(BEcombined$seurat_clusters, hue_pal())

plot(reducedDim(BE_sds), col = cell_colors_clust, pch = 16, cex = 0.5)
lines(BE_sds, lwd = 2, type = 'lineages', col = 'black')

nc <- 3
pt <- slingPseudotime(BE_sds)
nms <- colnames(pt)
nr <- ceiling(length(nms)/nc)
pal <- viridis(100, end = 0.95)
par(mfrow = c(nr, nc))
for (i in nms) {
  colors <- pal[cut(pt[,i], breaks = 100)]
  plot(reducedDim(BE_sds), col = colors, pch = 16, cex = 0.5, main = i)
  lines(BE_sds, lwd = 2, col = 'black', type = 'lineages')
}

# Get top highly variable genes
top_hvg <- HVFInfo(BEcombined) %>% 
  mutate(., bc = rownames(.)) %>% 
  arrange(desc(residual_variance)) %>% 
  top_n(300, residual_variance) %>% 
  pull(bc)

top_hvg <- HVFInfo(object = BEcombined, assay = "RNA")[1:300, ]
top_hvg <- HVFInfo(object = BEcombined[["RNA"]], selection.method = 'vst')[1:300, ]

# Prepare data for random forest
dat_use <- t(GetAssayData(BEcombined, slot = "data")[top_hvg,])
dat_use_df <- cbind(slingPseudotime(BE_sds)[,2], dat_use) # Do curve 2, so 2nd columnn
colnames(dat_use_df)[1] <- "pseudotime"
dat_use_df <- as.data.frame(dat_use_df[!is.na(dat_use_df[,1]),])

count.data <- GetAssayData(object = BEcombined[["RNA"]], slot = "counts")
count.data <- as.matrix(x = count.data + 1)
new.seurat.object <- SetAssayData(
  object = BEcombined,
  slot = "counts",
  new.data = count.data,
  assay = "RNA"
)

