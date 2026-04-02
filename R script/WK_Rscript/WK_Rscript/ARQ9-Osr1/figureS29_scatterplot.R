remove(list = ls())
library(multiHiCcompare)
library(diffHic)
library(lattice)
library(edgeR)
library(FIND)
library(csaw)
library(truncnorm)
library(gridExtra)
source('functions.R')
source('~/Dropbox/research projects/clipper/source042520/clipper042520_w2sided.R', echo=F)
rep1 = readRDS(file = 'primary_rep_processed.rds') 
rep2 = readRDS(file = 'replicate_rep_processed.rds')
ncores = 20-1
iter = 100
resolution = 10^6

FDR = (1:10)/100

rep1_binid = convert2binned(rep1)
rep2_binid = convert2binned(rep2)


############### setting parameters ###########
set.seed(1)
trueup_binidx = rbind(expand.grid(list(69:78, 143:157)), 
                      expand.grid(list(109:122, 187:197)),
                      expand.grid(list(143:152, 189: 198))
)

truedn_binidx = rbind(expand.grid(list(47:60, 112:121)), 
                      expand.grid(list(166: 177, 229: 238)),
                      expand.grid(list(51:80, 113:118)),
                      expand.grid(list(188:198, 240:249)))

true_binidx = rbind(trueup_binidx, truedn_binidx)
true_m = matrix(0, nrow = 250, ncol = 250)
for(i in 1:nrow(trueup_binidx)){
  true_idx_i = as.matrix(trueup_binidx[i,])
  true_m[true_idx_i[1], true_idx_i[2]] = 1
  true_m[true_idx_i[2], true_idx_i[1]] = 1
}
for(i in 1:nrow(truedn_binidx)){
  true_idx_i = as.matrix(truedn_binidx[i,])
  true_m[true_idx_i[1], true_idx_i[2]] = -1
  true_m[true_idx_i[2], true_idx_i[1]] = -1
}
levelplot(true_m)

true_binidx_char = convert2char(true_binidx[,1], true_binidx[,2])
totrow_idx = convert2char(rep1_binid[,1],rep1_binid[,2])
truerow_idx = match(true_binidx_char, totrow_idx)
mean(convert2char(trueup_binidx[,1],trueup_binidx[,2]) %in% totrow_idx)
mean(convert2char(truedn_binidx[,1],truedn_binidx[,2]) %in% totrow_idx)

#### sample log-fold change 
lgfc_up = sapply(1:nrow(trueup_binidx), function(i){
  true_idx_i = as.matrix(trueup_binidx[i,])
  mu_norm = 1/(true_idx_i[2] - true_idx_i[1])*100
  # mu_norm
  rtruncnorm(1, mean = mu_norm, sd = 0.5, a = 0.5)
})
lgfc_dn = sapply(1:nrow(truedn_binidx), function(i){
  true_idx_i = as.matrix(truedn_binidx[i,])
  mu_norm = 1/(true_idx_i[2] - true_idx_i[1])*100
  # mu_norm
  rtruncnorm(1, mean = - mu_norm, sd = 0.5, b = -0.5)
})
# lgfc_up = rnorm(nrow(trueup_binidx), mean = 1.5, sd = 0.5)
# lgfc_dn = -rnorm(nrow(truedn_binidx), mean = 1.5, sd = 0.5)
#### assign log-fc according to how far the indices are from the diagonal
true_binidx_char = convert2char(start1 = true_binidx[,1], start2 = true_binidx[,2])
lgfc = c(lgfc_up, lgfc_dn)

n_exp = 2
n_back = 2
mean_exp = compute_spikein_mean(rep1_binid, true_binidx, lgfc)
mean_back = rep1_binid
# levelplot(log10(convert2matrix(mean_back)))
# # which(m_exp >3.5*10^6, arr.ind = T)
# levelplot(log10(convert2matrix(mean_exp)))
# levelplot((convert2matrix(mean_back)))
# levelplot((convert2matrix(mean_exp)))
saveRDS(convert2matrix(mean_back), file = 'mean_back.rds')
saveRDS(convert2matrix(mean_exp), file = 'mean_exp.rds')
dat_exp = reshape2::melt(mean_exp, c("x", "y"), value.name = "z")
dat_back = reshape2::melt(mean_back, c("x", "y"), value.name = "z")


set.seed(1)
dat_ls_binned = generatedata(mean_exp = mean_exp, 
                             mean_back = mean_back,
                             n_exp = n_exp, n_back = n_back,
                             totrow_idx = totrow_idx)

exp1 = dat_ls_binned$exp[[1]]
exp2 = dat_ls_binned$exp[[2]]
exp1 = convert2matrix(exp1)
exp2 = convert2matrix(exp2)
n_maxd = 250
indices = as.matrix(expand.grid(list(1:n_maxd, 1:n_maxd)))
indices = indices[ indices[,1]>= indices[,2],]
n_pair = nrow(indices)

dist = t(sapply(1:n_pair, function(i){
  dist_ij = abs(indices[i,1] - indices[i,2])
  inter_ij = log10(exp1[indices[i,1],indices[i,2]]+1)
  return(c(dist_ij, inter_ij))
}))
cor(dist)
cor(dist, method = 'spearman')
pdf('scatterplot_hic_bindistVSintensity.png', height = 5, width = 5)
plot(dist, 
     xlab = 'Distance between two genomic regions', 
     ylab = 'Contact intensity between two genomic regions',
     pch = 16, col = alpha(1, 0.2), cex = 0.5)
text(x = 150, y = 6*10^5, label = 'Pearson cor: -0.16\n Spearman Cor: -0.44')
dev.off()

