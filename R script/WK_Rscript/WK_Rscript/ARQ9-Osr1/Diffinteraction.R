remove(list = ls())
library(multiHiCcompare)
library(diffHic)
library(lattice)
library(edgeR)
library(FIND)
library(csaw)
library(truncnorm)
library(gridExtra)
library(Clipper)
source('functions.R')
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


#### assign log-fc according to how far the indices are from the diagonal
true_binidx_char = convert2char(start1 = true_binidx[,1], start2 = true_binidx[,2])
lgfc = c(lgfc_up, lgfc_dn)

n_exp = 2
n_back = 2
mean_exp = compute_spikein_mean(rep1_binid, true_binidx, lgfc)
mean_back = rep1_binid
levelplot(log10(convert2matrix(mean_back)))
# which(m_exp >3.5*10^6, arr.ind = T)
levelplot(log10(convert2matrix(mean_exp)))
levelplot((convert2matrix(mean_back)))
levelplot((convert2matrix(mean_exp)))
saveRDS(convert2matrix(mean_back), file = 'mean_back.rds')


hic = mclapply(1:iter, function(it){
  set.seed(it)
  dat_ls_binned = generatedata(mean_exp = mean_exp, 
                               mean_back = mean_back,
                               n_exp = n_exp, n_back = n_back,
                               totrow_idx = totrow_idx)
  ####### clipper without log #######
  clipper_input = format2clipperinput(dat_ls_binned, n_exp, n_back)
  
  ####### clipper withl log #######
  re_clipper = clipper2sided(score_exp = log(clipper_input$exp),
                             score_back = log(clipper_input$back),
                             FDR = FDR, ifpowerful = F, nknockoff = 1)
  discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
  discovery_ls = lapply(discovery_ls, function(x){
    clipper_input$feature[x]
  })
  fdppow_clipperlog = compute_fdppow(discovery_ls = discovery_ls, trueidx = true_binidx_char )
  fdppow_clipperlog
  
  ################# hiccompare #################
  

  re_hiccompare = my_multihiccompare(score_exp_ls = dat_ls_binned$exp, 
                                     score_back_ls = dat_ls_binned$back,
                                     FDR = FDR,true_binidx_char = true_binidx_char, 
                                     resolution = resolution,
                                     ifnormalize = F)
  # fdppow_hiccompare = compute_fdppow(discovery_ls = re_hiccompare$discovery_ls, trueidx = true_binidx_char)
  re_hiccompare$fdppow
  ################# diffhic #################
  
  re_diffhic = my_diffhic(score_exp_ls = dat_ls_binned$exp,
                          score_back_ls = dat_ls_binned$back,
                          design = c(rep(1,n_exp), rep(2, n_back)),
                          FDR = FDR,
                          truerow_idx = truerow_idx,
                          ifnormalize = F)
  
  re_diffhic$fdppow
  ################# FIND #################
  ### normalize before input into Find
  re_find = my_find(score_exp_ls = dat_ls_binned$exp,
                    score_back_ls = dat_ls_binned$back,
                    FDR = FDR,
                    true_binidx_char = true_binidx_char)
  re_find
  
  fdp = list(clipper_log = fdppow_clipperlog[1,],
             hiccompare = re_hiccompare$fdppow[1,],
             diffhic = re_diffhic$fdppow[1,],
             find = re_find[1,])
  pow = list(clipper_log = fdppow_clipperlog[2,],
             hiccompare = re_hiccompare$fdppow[2,],
             diffhic = re_diffhic$fdppow[2,],
             find = re_find[2,])
  
  return(list(fdp = fdp, pow = pow))
  
  
},mc.cores = ncores)


saveRDS(hic, file = 'Diffinteraction_spikein_nbatcheffect.rds')
