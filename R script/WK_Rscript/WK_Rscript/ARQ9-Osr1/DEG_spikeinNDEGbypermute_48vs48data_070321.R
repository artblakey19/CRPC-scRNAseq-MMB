remove(list = ls())
library(parallel)
library(DESeq2)
library(edgeR)
library(IHW)
source('clipper042520_w2sided.R')
source('auxiliaryfunctions_070321.R', echo=TRUE)
filename = "DEG_spikeinNDEGbypermute_48vs48data_070321"
FDR = (1:10)/100
r1 = 3
r2 = 3
count = readRDS('count_matrices_gierlinski.rds')
count1 = count$count1
count2 = count$count2
r1_tot = ncol(count1)
r2_tot = ncol(count2)
####################
generate_subcount = function( count1, count2,r1, r2){
  
  
  # ## generate subcounts from condition 1 by random sampling
  # subcount1 = t(apply(count1, 1, sample, size = r1, replace = T))
  # subcount2 = t(apply(count1, 1, sample, size = r2, replace = T))
  # 
  # ## add spike-in true DEG
  # subcount1[trueDE_idx,] = t(apply(count1_trueDE,1, sample, size = r1, replace = T))
  # subcount2[trueDE_idx,] = t(apply(count2_trueDE,1, sample, size = r2, replace = T))
  # re = list(subcount1 = subcount1, subcount2 = subcount2)
  
  subcount1 = count1[,sample(r1_tot, r1, replace = F)]
  subcount2 = count2[,sample(r2_tot, r2, replace = F)]
  
  subcount_tot_nDE = cbind(subcount1[-trueDE_idx,], subcount2[-trueDE_idx,])
  # subcount_tot_nDE = count1[-trueDE_idx, ]
  subcount_tot_nDE = t(apply(subcount_tot_nDE, 1, sample, size = r1 + r2, replace = F))
  subcount1[-trueDE_idx,] = subcount_tot_nDE[,1:r1]
  subcount2[-trueDE_idx,] = subcount_tot_nDE[,-(1:r1)]
  
  # lgm1 = log(base =2, rowMeans(subcount1) + 1 )
  # lgm2 = log(base =2, rowMeans(subcount2) + 1 )
  # plot(lgm1, lgm2, cex = 0.3, col = 1)
  # points(lgm1[trueDE_idx], lgm2[trueDE_idx], col = 2, cex = 0.3)
  # points(lgm1[-trueDE_idx], lgm2[-trueDE_idx], col = 3, cex = 0.3)
  re = list(subcount1 = as.matrix(subcount1), subcount2 = as.matrix(subcount2))
  return(re)
}


############ use edgeR to normalize ############
dat = cbind(count1, count2)
cond_idx = rep(c(1,2), c(r1_tot, r2_tot))
y <- DGEList(counts = dat, group = cond_idx)
y <- calcNormFactors(y)
count_norm_edger = cpm(y) #### normalized matrix
count1 = count_norm_edger[, (1:r1_tot)]
count2 = count_norm_edger[,-(1:r1_tot)]



m1 = rowMeans(count1)
m2 = rowMeans(count2)
logfd = log(base = 2, (m1 + 1)/(m2+1))

# idx_m1_na = m1 == 0.0
# idx_m2_na = m2 == 0.0
# hist(logfd)
# sum(logfd > 2, na.rm = T)
trueDE_idx = which(abs(logfd) > 1.5)
trueDE = rownames(count1)[trueDE_idx]
count1_trueDE = count1[trueDE_idx,]
count2_trueDE = count2[trueDE_idx,]
fd_tot = rep(1, nrow(count1))
fd_tot[trueDE_idx] = (2^logfd)[trueDE_idx]
saveRDS(list(fd_tot = fd_tot, trueDE_idx = trueDE_idx), file = paste0(filename, '_trueFD&DEidx.rds'))
re = list(count1 = count2, count2 = count2, trueDE_idx = trueDE_idx)
saveRDS(re, file = 'realdata_yeast_perm.rds')
# hist(log(base = 2, (m1[trueDE_idx] + 1)/(m2[trueDE_idx]+ 1)), breaks = 100)
# plot(log(base = 2, m1+1), log(base = 2, m2 + 1), cex =0.3)
# points(log(base = 2, m1[trueDE_idx]+1), log(base = 2, m2[trueDE_idx] + 1), col = 2, cex = 0.3)
# points(log(base = 2, m1[-trueDE_idx]+1), log(base = 2, m2[-trueDE_idx] + 1), col = 3, cex = 0.3)

ncores = 34
iter = 100
#########
set.seed(1)

re = Compare(spikeInTarget = 'NDEG')

saveRDS(re, file = 'DEG_spikeinNDEGbypermute_48vs48data_070321.rds')


############ coarse plotting
library(ggplot2)
library(grid)
library(gridExtra)
library(egg)
library(gtable)
re =  readRDS(file = 'DEG_spikeinNDEGbypermute_48vs48data_070321.rds')


re_test_uniformpval = test_uniformpval(re, trueDE_idx)
pval_test_uniformpval = lapply(re_test_uniformpval, function(x){x[2,]})
naprop_test_uniformpval = lapply(re_test_uniformpval, function(x){x[1,]})
saveRDS(list(pval = pval_test_uniformpval, 
             naprop = naprop_test_uniformpval), 'pval_NDEGifuniform_yeast_perm.rds')


p_ls = coarseplot(re)
p1 = p_ls$p1
p2 = p_ls$p2
pdf(file = 'DEG_spikeinNDEGbypermute_48vs48data_070321.pdf',  width = 10, height = 5)
p = cbind(p1, p2)
grid.arrange(p)
dev.off()



