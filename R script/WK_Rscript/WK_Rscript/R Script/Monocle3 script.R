
library(Seurat)
library(devtools)
library(R.utils)
library(SeuratWrappers)
library(monocle3)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(remotes) 

setwd("//isi-dcnl/user_data/zjsun/group/Aly Buckley/HGF hMET Bcat Pb")

#Convert Seurat to Monocle3 cell data set class
cds <- as.cell_data_set(TriplevDouble.combined1.epi1)

###Necessary?
cds.partition <- c(rep(1,length(cds@colData@rownames)))
names(cds.partition) <- cds@colData@rownames
cds.partition <- as.factor(cds.partition)

cds@clusters$UMAP$partitions <- cds.partition
list_cluster <- TriplevDouble.combined1.epi1@active.ident
cds@int_colData@listData$reducedDims$UMAP <- TriplevDouble.combined1.epi1@reductions$umap@cell.embeddings


#Construction single cell trajectories
cds <- cluster_cells(cds = cds, reduction_method = "UMAP") #necessary?
cds <- learn_graph(cds, use_partition = FALSE)

#plot single cell trajectory
plot_cells(cds = cds,
           color_cells_by = "cluster",
           show_trajectory_graph = TRUE,
           label_branch_points = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           group_label_size = 4)

plot_cells(cds = cds,
           color_cells_by = "EpiCellTypes",
           show_trajectory_graph = TRUE,
           label_branch_points = FALSE,
           label_leaves=FALSE,
           label_roots = FALSE,
           group_label_size = 4)

#plot single cell trajectory by pseudotime
cds <- order_cells(cds, reduction_method = "UMAP", root_cells = colnames(cds[clusters(cds)==2]))##looks incorrect, manual may be better ? 
cds <- order_cells(cds, reduction_method = "UMAP") #works
plot_cells(cds=cds,
           color_cells_by = "pseudotime",
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           label_roots = FALSE,
           graph_label_size=4)

pseudotime(cds)
cds$monocle3_pseudotime <- pseudotime(cds)
data.pseudo <- as.data.frame(colData(cds))
ggplot(data.pseudo, aes(monocle3_pseudotime, reorder(EpiCellTypes, monocle3_pseudotime, median), fill = EpiCellTypes)) + geom_boxplot()

#DEGs Monocle3
deg.cds <- graph_test(cds, neighbor_graph="principal_graph", cores=4) #correct parameter 4?
write.csv(deg.cds, file = "deg.cds.csv")

#Plot cells gene expression
target.genes <- c("hMETtg","hHGFtg","Ar","Xpo1","Ctnnb1","Myc")
plot_cells(cds=cds,
           genes = target.genes,
           show_trajectory_graph = FALSE,
           label_cell_groups=FALSE)


