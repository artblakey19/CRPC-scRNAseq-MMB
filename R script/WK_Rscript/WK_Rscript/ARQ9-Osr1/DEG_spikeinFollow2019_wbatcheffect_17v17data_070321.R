remove(list = ls())
library(parallel)
library(DESeq2)
library(edgeR)
library(IHW)
library(scales)
library(preprocessCore)
source('clipper042520_w2sided.R')
source('auxiliaryfunctions_070321.R', echo=TRUE)
filename = "DEG_spikeinFollow2019_wbatcheffect_17vs17data_070321"
FDR = (1:10)/100
# dat_tot = readRDS(file = paste0('trueDE_FDR', 100*FDR,'.rds'))
r1 = 3
r2 = 3
r = 17
count1 = readRDS('count1.rds')
count2 = readRDS('count2.rds')
# cor(count1[, c(-3, -10)])
# sum(cor(count1[,-c(1,2,3,6,10, 11, 13)]) <= 0.8)
# which(cor(count1[, -c(1,2,3,6,10, 11, 13)]) <= 0.8, arr.ind = T)
# 
# cor(count2)
# count2 = count2[, c(-2, -6)]
# count1 = count1[, -c(1,2,3,6,10, 11, 13)]

# idx2rm = (rowSums(count2)==0)  | (rowSums(count1)==0)
# count1 = count1[!idx2rm, ]
# count2 = count2[!idx2rm, ]
# sum(rowSums(count2) == 0)
# sum(rowSums(count1) ==  0)
# geneid = rownames(count1)
# sample_idx1 = colnames(count1)
# sample_idx2 = colnames(count2)
# 
# prop0_1 = apply(count1, 1, function(x){
#   mean(x==0)
# })
# prop0_2 = apply(count2, 1, function(x){
#   mean(x==0)
# })
# hist(prop0_1)
# hist(prop0_2)
# 
# idx2rm = (prop0_2 > 0.3 )
# # idx2rm = rep(F, nrow(count1))
# sum(idx2rm)
# count1 = count1[!idx2rm, ]
# count2 = count2[!idx2rm, ]
# mean(count1 == 0)
# mean(count2 == 0)
# sum(rowSums(count2) == 0)
# sum(rowSums(count1) ==  0)
# 
# rownames(count1) = geneid[!idx2rm]
# rownames(count2) = geneid[!idx2rm]
########

# ############ quantile normalize ############
# count1 = normalize.quantiles(as.matrix(count1),copy = TRUE)
# count2 = normalize.quantiles(as.matrix(count2),copy = T)

# prop0_1 = apply(count1, 1, function(x){
#   mean(x==0)
# })
# prop0_2 = apply(count2, 1, function(x){
#   mean(x==0)
# })
# 
# idx2rm = (prop0_2 > 0.3)
# sum(idx2rm)
# count1 = count1[!idx2rm, ]
# count2 = count2[!idx2rm, ]
# mean(count1 == 0)
# mean(count2 == 0)
# sum(rowSums(count2) == 0)
# sum(rowSums(count1) ==  0)
# 
# rownames(count1) = geneid[!idx2rm]
# rownames(count2) = geneid[!idx2rm]



m1 = rowMeans(count1)
m2 = rowMeans(count2)
idx0_1 = apply(count1, 1, function(x){
  any(x == 0)
})
idx0_2 = apply(count2, 1, function(x){
  any(x == 0)
})

idx2fd = which( (!idx0_1) & (!idx0_2))
# idx2rm = which(m1<=1 | m2<=1)
# sum(m2<=1)
# fd = readRDS(file= 'foldchange_fromGierlinski2015.rds')
fd =  ((m1 + 1)/(m2 + 1))
sd_ls = apply(count2, 1, sd)
hist(fd/sd_ls, breaks = 100)
median(fd/sd_ls)
# plot(y = fd, x = sd_ls)
fd = fd[ fd >= 16 ]
length(fd) ## 205
hist(fd, breaks = 1000)
# mu = mean(logfd)
# sd = sd(logfd)
set.seed(12345)
trueDE_idx = sample(idx2fd, min(0.3*nrow(count1), length(idx2fd)))
# trueDE_idx = which(0 > 1)
trueDE = rownames(count1)[trueDE_idx]
fd_trueDE = sample_fd(fd = fd, n = length(trueDE_idx))
hist(fd_trueDE)
# *(as.numeric(runif(length(trueDE_idx)) < 0.5)*2 - 1)
# sd_ls = apply(count2, 1, sd)
# saveRDS(sd_ls[-trueDE_idx], file = 'sd_ls_old.rds')
fd_tot = rep(1, nrow(count1))
fd_tot[trueDE_idx] = fd_trueDE
saveRDS(list(fd_tot = fd_tot, trueDE_idx = trueDE_idx), file = paste0(filename, '_trueFD&DEidx.rds'))
count2_wfd = apply(count2, 2, function(x){x*fd_tot})
re = list(count_wofd = count2, count_wfd = count2_wfd, trueDE_idx = trueDE_idx)
saveRDS(re, file = 'realdata_human_2019.rds')




ncores = 34
iter = 100
########
set.seed(1)

re = Compare(spikeInTarget = 'DEG')


saveRDS(re, file = paste0(filename,'_fd16.rds'))

################################
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
re = readRDS(paste0(filename,'_fd16.rds'))

re_test_uniformpval = test_uniformpval(re, trueDE_idx)
pval_test_uniformpval = lapply(re_test_uniformpval, function(x){x[2,]})
naprop_test_uniformpval = lapply(re_test_uniformpval, function(x){x[1,]})
saveRDS(list(pval = pval_test_uniformpval, 
             naprop = naprop_test_uniformpval), 'pval_NDEGifuniform_human_2019.rds')

# pdf(file = paste0('boxplot_',filename,'_nDEG_uniformtest.pdf'), height = 5, width = 7)
# re_test_uniformpval = lapply(re_test_uniformpval, function(x){x[!is.na(x)]})
# vioplot( re_test_uniformpval,
#          # as.matrix(as.data.frame(re_test_uniformpval)),
#          main = filename)
# dev.off()
# saveRDS(re_test_uniformpval, file = 'pval_NDEGifuniform_human_2019.rds')
# 

p_ls = coarseplot(re)
p1 = p_ls$p1
p2 = p_ls$p2
pdf(paste0(filename,'_fd16.pdf'),  width = 10, height = 5)
p = cbind(p1, p2)
grid.arrange(p)
dev.off()


# dat = re
# print(length(re[[1]]$fdp))
# fdp = lapply(dat, '[[','fdp')
# fdp = lapply(1:length(re[[1]]$fdp), function(i){
#   rowMeans(sapply(fdp, '[[',i))
# })
# 
# pow = lapply(dat, '[[','pow')
# pow = lapply(1:length(re[[1]]$fdp), function(i){
#   rowMeans(sapply(pow, '[[',i))
# })
# 
# fdp
# pow
# ################### lineplots ###################
# library(ggplot2)
# library(grid)
# library(gridExtra)
# library(egg)
# library(gtable)
# 
# method_names = names(re[[1]]$fdp)
# dat = cbind.data.frame(fdr = unlist(fdp),
#                        pow = unlist(pow),
#                        methods = rep(method_names, each = 10),
#                        q = rep(FDR, length(method_names)))
# p1 = ggplot(dat,aes(x = q, y = fdr, group = methods)) +
#   geom_line(aes(q, y = fdr, col = methods), lwd = 1) + 
#   geom_abline(slope = 1, lty = 'dashed', color ="#525252", size = 0.5)  
# 
# p2 = ggplot(dat,aes(x = q, y = pow, group = methods)) +
#   geom_line(aes(q, y = pow, col = methods), lwd = 1) 
# 
# pdf(file = 'fdr_DEG_spikeinFollow2019_17vs17data_061821.pdf',  width = 5, height = 5)
# p1
# dev.off()
# pdf(file = 'pow_DEG_spikeinFollow2019_17vs17data_061821.pdf',  width = 5, height = 5)
# p2
# dev.off()
# 
# 
# 
# 
