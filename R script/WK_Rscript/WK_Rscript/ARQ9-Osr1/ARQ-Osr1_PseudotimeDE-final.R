setwd("//isi-dcnl/user_data/zjsun/group/Won Kyung Kim/ARQ9-Osr1/Nat Comm Revision/Fig.6/pseudotime/PseudotimeDE")

library(PseudotimeDE)
library(SingleCellExperiment)
library(slingshot)
library(tibble)
library(dplyr)
library(scales)
library(irlba)
library(Seurat)

list.files()
ARQPosEpi <- readRDS("ARQPosEpi.rds")
ARQNegEpi <- readRDS("ARQNegEpi.rds")

ARQPosEpi1 <- as.SingleCellExperiment(ARQPosEpi)
ARQNegEpi1 <- as.SingleCellExperiment(ARQNegEpi)

##ARQPosEpi

####Perform pseudotime inference on the original dataset####
rd <- irlba::prcomp_irlba(t(logcounts(ARQPosEpi1)), scale. = FALSE)$x[, 1:2]

reducedDims(ARQPosEpi1) <- SimpleList(UMAP = rd)
colData(ARQPosEpi1)$cl <- 1

ARQPosEpi_ori <- slingshot(ARQPosEpi1, reducedDim = 'UMAP', clusterLabels = "cl")
ARQPosEpi_ori_tbl <- tibble(cell = colnames(ARQPosEpi1), pseudotime = rescale(colData(ARQPosEpi_ori)$slingPseudotime_1))

head(ARQPosEpi_ori_tbl)

####Perform pseudotime inference on subsamples####
set.seed(123)
options(mc.cores = 2)
n = 1000

## Ganerate random subsamples
ARQPosEpi_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(ARQPosEpi1)[2]), size = 0.8*dim(ARQPosEpi1)[2], replace = FALSE)
})

ARQPosEpi_sub_tbl <- mclapply(ARQPosEpi_index, function(x, sce) {
  sce <- sce[, x]
  rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(UMAP = rd)
  
  fit <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = "cl")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, ARQPosEpi_ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = ARQPosEpi1)

PseudotimeDE::plotUncertainty(ARQPosEpi_ori_tbl, ARQPosEpi_sub_tbl)

####Perform DE test####
system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5", "Krt5", "Krt19"),
                                                 ori.tbl = ARQPosEpi_ori_tbl,
                                                 sub.tbl = ARQPosEpi_sub_tbl[1:100], ## To save time, use 100 subsamples
                                                 mat = ARQPosEpi1, ## You can also use a matrix or SeuratObj as the input
                                                 model = "nb",
                                                 mc.cores = 2))

print(res)

####Visualization####
ARQPosEpi2 <- GetAssayData(object = ARQPosEpi, slot = "counts")
ARQPosEpi2 <- as.matrix(ARQPosEpi2)
#PseudotimeDE::plotCurve showed error like below. So I newly generated matrix file like upper...
#Error in t.default(count_mat) : argument is not a matrix

tiff(file = "ARQPosEpi PseudotimeDE ARQ, Axin2, Ccnd1, Krt19, Krt5, Lgr5, Tcf4.tiff", width = 6, height = 12, units = "in", compression = "lzw", res = 800)
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi2,
                        model.fit = res$gam.fit)
dev.off()

##ARQNegEpi

####Perform pseudotime inference on the original dataset####
rd <- irlba::prcomp_irlba(t(logcounts(ARQNegEpi1)), scale. = FALSE)$x[, 1:2]

reducedDims(ARQNegEpi1) <- SimpleList(UMAP = rd)
colData(ARQNegEpi1)$cl <- 1

ARQNegEpi_ori <- slingshot(ARQNegEpi1, reducedDim = 'UMAP', clusterLabels = "cl")
ARQNegEpi_ori_tbl <- tibble(cell = colnames(ARQNegEpi1), pseudotime = rescale(colData(ARQNegEpi_ori)$slingPseudotime_1))

head(ARQNegEpi_ori_tbl)

####Perform pseudotime inference on subsamples####
set.seed(123)
options(mc.cores = 2)
n = 100

## Ganerate random subsamples
ARQNegEpi_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(ARQNegEpi1)[2]), size = 0.8*dim(ARQNegEpi1)[2], replace = FALSE)
})

ARQNegEpi_sub_tbl <- mclapply(ARQNegEpi_index, function(x, sce) {
  sce <- sce[, x]
  rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(UMAP = rd)
  
  fit <- slingshot(sce, reducedDim = 'UMAP', clusterLabels = "cl")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, ARQNegEpi_ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = ARQNegEpi1)

PseudotimeDE::plotUncertainty(ARQNegEpi_ori_tbl, ARQNegEpi_sub_tbl)

####Perform DE test####
system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5", "Krt5", "Krt19"),
                                                 ori.tbl = ARQNegEpi_ori_tbl,
                                                 sub.tbl = ARQNegEpi_sub_tbl[1:100], ## To save time, use 100 subsamples
                                                 mat = ARQNegEpi1, ## You can also use a matrix or SeuratObj as the input
                                                 model = "nb",
                                                 mc.cores = 2))

print(res)

####Visualization####
ARQNegEpi2 <- GetAssayData(object = ARQNegEpi, slot = "counts")
ARQNegEpi2 <- as.matrix(ARQNegEpi2)
#PseudotimeDE::plotCurve showed error like below. So I newly generated matrix file like upper...
#Error in t.default(count_mat) : argument is not a matrix

tiff(file = "ARQNegEpi PseudotimeDE ARQ, Axin2, Ccnd1, Krt19, Krt5, Lgr5, Tcf4.tiff", width = 6, height = 12, units = "in", compression = "lzw", res = 800)
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi2,
                        model.fit = res$gam.fit)
dev.off()