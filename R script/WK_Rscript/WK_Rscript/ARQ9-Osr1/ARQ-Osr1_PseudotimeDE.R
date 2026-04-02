#Running PseudotimeDE on HPC

install.packages(c("SingleCellExperiment", "tibble", "irlba"))

library(Seurat)
library(devtools)

devtools::install_github("SONGDONGYUAN1994/PseudotimeDE")

library(dplyr)
library(Matrix)
library(cowplot)
library(ggplot2)
library(reticulate)
library(monocle)
library(slingshot)
library(RColorBrewer)
library(SCENIC)
library(GENIE3)
library(AUCell)
library(RcisTarget)

library(tidyverse)
library(tidymodels)
library(scales)
library(viridis)
library(dplyr)
devtools::install_github("SONGDONGYUAN1994/PseudotimeDE")
library(PseudotimeDE)
library(SingleCellExperiment)
library(tibble)
library(irlba)

list.files()
ARQPosEpi <- readRDS("ARQPosEpi.rds")
ARQNegEpi <- readRDS("ARQNegEpi.rds")

ARQPosEpi1 <- as.SingleCellExperiment(ARQPosEpi)
ARQNegEpi1 <- as.SingleCellExperiment(ARQNegEpi)

##ARQPosEpi

####Dataset####
data(ARQPosEpi1, package = "PseudotimeDE")

####Perform pseudotime inference on the original dataset####
rd <- irlba::prcomp_irlba(t(logcounts(ARQPosEpi1)), scale. = FALSE)$x[, 1:2]

reducedDims(ARQPosEpi1) <- SimpleList(PCA = rd)
colData(ARQPosEpi1)$cl <- 3

fit_ori <- slingshot(ARQPosEpi1, reducedDim = 'PCA', clusterLabels = "cl")
ARQPosEpi_ori_tbl <- tibble(cell = colnames(ARQPosEpi1), pseudotime = rescale(colData(fit_ori)$slingPseudotime_1))

head(ARQPosEpi_ori_tbl)

####Perform pseudotime inference on subsamples####
set.seed(123)
## Set the cores for parallelization. Note that mclapply doesnot work on Windows.
options(mc.cores = 2)

## Number of subsmaples
n = 100

## Ganerate random subsamples
ARQPosEpi_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(ARQPosEpi1)[2]), size = 0.8*dim(ARQPosEpi1)[2], replace = FALSE)
})

ARQPosEpi_sub_tbl <- mclapply(ARQPosEpi_index, function(x, sce) {
  sce <- sce[, x]
  rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(PCA = rd)
  
  fit <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "cl")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, ARQPosEpi_ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = ARQPosEpi1)

data(ARQPosEpi_ori_tbl, package = "PseudotimeDE")
data(ARQPosEpi_sub_tbl, package = "PseudotimeDE")

PseudotimeDE::plotUncertainty(ARQPosEpi_ori_tbl, ARQPosEpi_sub_tbl)

Cairo::CairoPDF("output/ARQPosEpi_plotUncertainty.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotUncertainty(ARQPosEpi_ori_tbl, ARQPosEpi_sub_tbl)
dev.off()

####Perform DE test####
system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5"),
                                                 ori.tbl = ARQPosEpi_ori_tbl,
                                                 sub.tbl = ARQPosEpi_sub_tbl[1:100], ## To save time, use 100 subsamples
                                                 mat = ARQPosEpi1, ## You can also use a matrix or SeuratObj as the input
                                                 model = "nb",
                                                 mc.cores = 2))

print(res)

PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi1,
                        model.fit = res$gam.fit)

Cairo::CairoPDF("output/ARQPosEpi_plotCurve_res.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi1,
                        model.fit = res$gam.fit)
dev.off()

system.time(res_zinb <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5"),
                                                      ori.tbl = ARQPosEpi_ori_tbl,
                                                      sub.tbl = ARQPosEpi_sub_tbl[1:100],
                                                      mat = ARQPosEpi1,
                                                      model = "zinb"))

print(res_zinb)

PseudotimeDE::plotCurve(gene.vec = res_zinb$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi1,
                        model.fit = res_zinb$gam.fit)

Cairo::CairoPDF("output/ARQPosEpi_plotCurve_res_zinb.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi1,
                        model.fit = res$gam.fit)
dev.off()

system.time(res_gaussian <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5"),
                                                          ori.tbl = ARQPosEpi_ori_tbl,
                                                          sub.tbl = ARQPosEpi_sub_tbl[1:100],
                                                          mat = ARQPosEpi1,
                                                          model = "gaussian", 
                                                          assay.use = "logcounts"))

print(res_gaussian)

PseudotimeDE::plotCurve(gene.vec = res_gaussian$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi1,
                        model.fit = res_gaussian$gam.fit, assay.use = "logcounts")

Cairo::CairoPDF("output/ARQPosEpi_plotCurve_res_gaussian.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQPosEpi_ori_tbl,
                        mat = ARQPosEpi1,
                        model.fit = res$gam.fit)
dev.off()

save.image("ARQ9-Osr1_PseudotimeDE.RData")

##ARQNegEpi

####Dataset####
data(ARQNegEpi1, package = "PseudotimeDE")

####Perform pseudotime inference on the original dataset####
rd <- irlba::prcomp_irlba(t(logcounts(ARQNegEpi1)), scale. = FALSE)$x[, 1:2]

reducedDims(ARQNegEpi1) <- SimpleList(PCA = rd)
colData(ARQNegEpi1)$cl <- 1

fit_ori1 <- slingshot(ARQNegEpi1, reducedDim = 'PCA', clusterLabels = "cl")
ARQNegEpi_ori_tbl <- tibble(cell = colnames(ARQNegEpi1), pseudotime = rescale(colData(fit_ori1)$slingPseudotime_1))

head(ARQNegEpi_ori_tbl)

####Perform pseudotime inference on subsamples####
set.seed(123)
## Set the cores for parallelization. Note that mclapply doesnot work on Windows.
options(mc.cores = 2)

## Number of subsmaples
n = 100

## Ganerate random subsamples
ARQNegEpi_index <- mclapply(seq_len(n), function(x) {
  sample(x = c(1:dim(ARQNegEpi1)[2]), size = 0.8*dim(ARQNegEpi1)[2], replace = FALSE)
})

ARQNegEpi_sub_tbl <- mclapply(ARQNegEpi_index, function(x, sce) {
  sce <- sce[, x]
  rd <- irlba::prcomp_irlba(t(logcounts(sce)), scale. = FALSE)$x[, 1:2]
  reducedDims(sce) <- SimpleList(PCA = rd)
  
  fit <- slingshot(sce, reducedDim = 'PCA', clusterLabels = "cl")
  tbl <- tibble(cell = colnames(sce), pseudotime = rescale(colData(fit)$slingPseudotime_1))
  
  ## Make sure the direction of pseudotime is the same as the original pseudotime
  merge.tbl <- left_join(tbl, ARQNegEpi_ori_tbl, by = "cell")
  
  if(cor(merge.tbl$pseudotime.x, merge.tbl$pseudotime.y) < 0) {
    tbl <- dplyr::mutate(tbl, pseudotime = 1-pseudotime)
  }
  tbl
}, sce = ARQNegEpi1)

data(ARQNegEpi_ori_tbl, package = "PseudotimeDE")
data(ARQNegEpi_sub_tbl, package = "PseudotimeDE")

PseudotimeDE::plotUncertainty(ARQNegEpi_ori_tbl, ARQNegEpi_sub_tbl)

Cairo::CairoPDF("output/ARQNegEpi_plotUncertainty.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotUncertainty(ARQNegEpi_ori_tbl, ARQNegEpi_sub_tbl)
dev.off()

####Perform DE test####
system.time(res <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5"),
                                                 ori.tbl = ARQNegEpi_ori_tbl,
                                                 sub.tbl = ARQNegEpi_sub_tbl[1:100], ## To save time, use 100 subsamples
                                                 mat = ARQNegEpi1, ## You can also use a matrix or SeuratObj as the input
                                                 model = "nb",
                                                 mc.cores = 2))

print(res)

PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi1,
                        model.fit = res$gam.fit)

Cairo::CairoPDF("output/ARQNegEpi_plotCurve_res.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi1,
                        model.fit = res$gam.fit)
dev.off()

system.time(res_zinb <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5"),
                                                      ori.tbl = ARQNegEpi_ori_tbl,
                                                      sub.tbl = ARQNegEpi_sub_tbl[1:100],
                                                      mat = ARQNegEpi1,
                                                      model = "zinb"))

print(res_zinb)

PseudotimeDE::plotCurve(gene.vec = res_zinb$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi1,
                        model.fit = res_zinb$gam.fit)

Cairo::CairoPDF("output/ARQNegEpi_plotCurve_res_zinb.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi1,
                        model.fit = res$gam.fit)
dev.off()

system.time(res_gaussian <- PseudotimeDE::runPseudotimeDE(gene.vec = c("ARQ", "Tcf4", "Axin2", "Ccnd1", "Lgr5"),
                                                          ori.tbl = ARQNegEpi_ori_tbl,
                                                          sub.tbl = ARQNegEpi_sub_tbl[1:100],
                                                          mat = ARQNegEpi1,
                                                          model = "gaussian", 
                                                          assay.use = "logcounts"))

print(res_gaussian)

PseudotimeDE::plotCurve(gene.vec = res_gaussian$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi1,
                        model.fit = res_gaussian$gam.fit, assay.use = "logcounts")

Cairo::CairoPDF("output/ARQNegEpi_plotCurve_res_gaussian.pdf", width=15, height=10)
par(mfrow=c(4,6))
PseudotimeDE::plotCurve(gene.vec = res$gene,
                        ori.tbl = ARQNegEpi_ori_tbl,
                        mat = ARQNegEpi1,
                        model.fit = res$gam.fit)
dev.off()
