remove(list = ls())
#source('~/Dropbox/clipper/source042520/clipper042520_w2sided.R', echo=F)
source('clipper042520_w2sided.R', echo=F)

library(parallel)
library(locfdr)
library(truncdist)
library(ks)
library(qvalue)
ncores = detectCores() - 2
iter = 200

n <- 10000
q_tot = (1:10)/100
out_prop <- c(0, 0.1, 0.2)


####### functions ######

compute_pooled_pval = function(back, exp){
  re = sapply(exp, function(x){
    mean(back >= x, na.rm = T)
  })
  return(re)
  
}

compute_pooled_fdr = function(back, exp, pit = 1){
  re = sapply(exp, function(x){
    sum(back >= x, na.rm = T)*pit/max(1, sum(exp >= x, na.rm = T))
  })
  return(re)
}

compute_pair_fdr = function(back, exp){
  diff1 <- exp - back
  diff2 <- back - exp
  re = sapply(diff1, function(x){
    min(1, sum(diff2 >= x, na.rm = T)/max(1, sum(diff1 >= x, na.rm = T)))
  })
  return(re)
}

find_discovery_wpval = function(pval, q){
  pval.adj = p.adjust(pval, method = 'BH')
  discovery_ls = lapply(q, function(q_i){
    discovery = which(pval.adj <= q_i)
  })
  return(discovery_ls)
}

find_discovery_wqval = function(pval, q){
  qval = try(qvalue(pval), silent = T)
  if(class(qval) == "try-error"){
    qval = qvalue(pval, pi0 = 1)
  }else{
    discovery_ls = lapply(q, function(q_i){
      discovery = which(qval$qvalues <= q_i)
    })
  }
  return(discovery_ls)
}

find_discovery_wfdr = function(fdr, q){
  discovery_ls = lapply(q, function(q_i){
    discovery = which(fdr<= q_i)
  })
  return(discovery_ls)
}

compute_fdppow = function(discovery_ls, trueidx){
  sapply(discovery_ls, function(discovery){
    fdp = sum(!discovery %in% trueidx )/max(length(discovery),1)
    pow = sum(discovery %in% trueidx)/length(trueidx)
    return(c(fdp, pow))
  })
}

compute_paired_pval_norm = function(exp, back, sd =1){
  pnorm(exp-back, mean = 0, sd = sd, lower.tail = F)
}


find_discovery_wlocfdr = function(loc.fdr, q ){
  if(class(loc.fdr)== "try-error"){
    discovery_ls = rep(list(which(1<0)), length(q))
  }else{
    discovery_ls = lapply(q, function(q_i){
      discovery = which(loc.fdr <= q_i)
    })
  }
  
  return(discovery_ls)
}


compute_paired_pval_nb_c = function(x,y, nb.p = 1/2){
  k = x-y
  size_hat = mean(c(x,y))
  y_ls = 0:100
  tailprob =  sapply(y_ls, function(y_i){
    # print(y)
    prod1 = pnbinom(y_i + k - 1, size = size_hat, prob = nb.p, lower.tail = F)
    prod2 = dnbinom(y_i, size = size_hat, prob = nb.p)
    if(identical(prod1, 0) | identical(prod2, 0)){
      return(0)
    }else{
      return(prod1*prod2)
    }
  })
  sum(tailprob)
}



compute_paired_pval_pois_mis = function(exp, back, sd = 1){
  logexp = log(exp + 0.01)
  logback =  log(back + 0.01)
  pnorm( logexp - logback, mean = 0, sd = sd, lower.tail = F)
}

localfdr_multrep_p = function(exp, back){
  
  Z = 0.5* (exp[,1] + exp[,2] - back[,1] - back[,2])
  z0 = 0.5* (exp[,1] + back[,1]) - 0.5 * (exp[,2] + back[,2])
  model1 <- kde(Z)
  model2 <- kde(z0)
  Z_p <- pmin(1, (1-pkde(Z, model2))/(1-pkde(Z, model1)))
  return(Z_p)
}

localfdr_multrep_p_miss = function(exp, back){
  zandz0 <- sapply(1:nrow(exp), function(j) {
    exp0 <- na.omit(exp[j,])
    back0 <- na.omit(back[j,])
    if(length(exp0 > 0 & length(back0 >0))) {
      Z <- 0.5* (exp0[1] + exp0[2] - back0[1] - back0[2])
    } else {Z <- NA}
    if(length(exp0 > 1 & length(back0 >1))) {
      z0 <- 0.5* (exp0[1] + back0[1] - exp0[2] - back0[2])
    } else {z0 <- NA}
    return(c(Z, z0))
  })
  Z <- zandz0[1,]
  z0 <- zandz0[2,]
  model1 <- kde(Z[!is.na(Z)])
  model2 <- kde(z0[!is.na(z0)])
  Z[is.na(Z)] <- 0
  z0[is.na(z0)] <- 0
  Z_p <- pmin(1, dkde(Z, model2)/dkde(Z, model1))
  return(Z_p)
}

########## distribution of difference between two nb ###########
twosample_nb = function(x,s1,s2,p1,p2)  {
  p = 0
  for (y in 0:(s2*3)) {
    p = p + pnbinom(y + x -1, s1, p1, lower.tail = F) * dnbinom(y, s2, p1 )
  }
  return(p)
}


########################################################################################
######################   Result 1  #########################
######################  Fig1 multiple replicates, generated from norm   #############
########################################################################################

##############outlier################
outlier_norm = function(n, mean.norm){
  a = qnorm(0.99, mean.norm )
  rtrunc(n, spec = 'norm', mean = mean.norm, a = a)
}

set.seed(1)

nonnull_prop <- c(0.2, 0.3, 0.4)


sim_multirep_norm_homo <- lapply(1:length(nonnull_prop), function(i_prop){
  set.seed(1)
  prop =  nonnull_prop[i_prop]
  trueidx = 1:(prop*n)
  
  #simulate mean of each point from Gaussian
  mean_back = rep(0, n*3)
  mean_exp = rep(0, n*3)
  mean_exp[1:(length(trueidx)*3)] = rep(rnorm(length(trueidx), mean = 5, sd = 1), each = 3)
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnorm(n*3, mean = 0, sd = 1) + mean_back
    score_exp <- rnorm(n*3, mean = 0, sd = 1) + mean_exp
 
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    pval_pooled[is.na(pval_pooled)] <- 1
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "greater", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], mu = mean(score_back[j,!(is.na(score_back[j,]))]), alternative = "greater")$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "greater", var.equal = F)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method  = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method  = 'max', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method = 'diff', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})
    
    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pair_c = fdppow_pair_c[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pair_c = fdppow_pair_c[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')
  
  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr, 
                             q = rep(q_tot, length(methods)), 
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_multirep_norm_homo) = nonnull_prop
saveRDS(sim_multirep_norm_homo, 'sim_multirep_norm_homo_prop_interesting.rds')





sim_multirep_norm_hete <- lapply(1:length(nonnull_prop), function(i_prop){
  
  prop =  nonnull_prop[i_prop]
  set.seed(1)
  trueidx = 1:(prop*n)
  #simulate mean of each point from Gaussian
  mean_back = rep(rnorm(n, mean = 0, sd = 2), each = 3)
  mean_exp = mean_back
  mean_exp[1:(length(trueidx)*3)] = rep(rnorm(length(trueidx), mean = 5, sd = 1), each = 3)
  
  # print(paste0('divergence: ', v, '; N = ',n))
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnorm(n*3, mean = 0, sd = 1) + mean_back
    score_exp <- rnorm(n*3, mean = 0, sd = 1) + mean_exp
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "greater", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], mu = mean(score_back[j,!(is.na(score_back[j,]))]), alternative = "greater")$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (sum(is.na(score_back[j,]))>1|sum(is.na(score_exp[j,]))>1) {
        return(1)} else {
          return(t.test(score_exp[j,!(is.na(score_exp[j,]))], score_back[j,!(is.na(score_back[j,]))], alternative = "greater", var.equal = F)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method  = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method  = 'max', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method = 'diff', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
 
    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})
    
    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pair_c = fdppow_pair_c[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pair_c = fdppow_pair_c[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')
  
  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr, 
                             q = rep(q_tot, length(methods)), 
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})
names(sim_multirep_norm_hete) = out_prop
saveRDS(sim_multirep_norm_hete, 'sim_multirep_norm_hete_prop_interesting.rds')


######################  Fig1 multiple replicates, generated from Pois   #############
########################################################################################

outlier_pois = function(n, mean.pois){
  a = qpois(0.99, mean.pois )
  rtrunc(n, spec = 'pois', lambda = mean.pois, a = a)
}


#simulate mean of each point from Poisson
sim_multirep_pois_homo <- lapply(1:length(nonnull_prop), function(i_prop){
  
  prop =  nonnull_prop[i_prop]
  set.seed(1)
  trueidx = 1:(prop*n)
  #simulate mean of each point from Poisson
  mean_back = rep(20, n*3)
  mean_exp = rep(20, n*3)
  mean_exp[1:(length(trueidx)*3)] = rep(rpois(length(trueidx), lambda = 40), each = 3)
  
  
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rpois(n*3, lambda = mean_back)
    score_exp <- rpois(n*3, lambda = mean_exp)

    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
   
    ##### use BH paired correct specification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "greater")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(1 - ppois(sum(score_exp[j,], na.rm = T)-1, lambda = length(na.omit(score_exp))/length(na.omit(score_back))*sum(score_back[j,], na.rm = T)))
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (length(unique(na.omit(score_exp[j,])))<2|length(unique(na.omit(score_back[j,])))<2) {
        return(1)} else {
          return(t.test(-log(score_exp[j,!(is.na(score_exp[j,]))]+0.01), -log(score_back[j,!(is.na(score_back[j,]))]+0.01), alternative = "greater", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method  = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method  = 'max', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method = 'diff', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})
    
    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pair_c = fdppow_pair_c[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pair_c = fdppow_pair_c[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')
  
  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr, 
                             q = rep(q_tot, length(methods)), 
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_multirep_pois_homo) = out_prop
saveRDS(sim_multirep_pois_homo, 'sim_multirep_pois_homo_prop_interesting.rds')




sim_multirep_pois_hete <- lapply(1:length(nonnull_prop), function(i_prop){
  
  prop =  nonnull_prop[i_prop]
  set.seed(1)
  trueidx = 1:(prop*n)
  mean_back = rep(rnbinom(n, size = 20, prob = 0.5), each = 3)
  mean_exp = mean_back
  mean_exp[1:(length(trueidx)*3)] = rep(rpois(length(trueidx), lambda = 40), each = 3)
  
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rpois(n*3, lambda = mean_back)
    score_exp <- rpois(n*3, lambda = mean_exp)
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired correct specification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "greater")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(1 - ppois(sum(score_exp[j,], na.rm = T)-1, lambda = length(na.omit(score_exp))/length(na.omit(score_back))*sum(score_back[j,], na.rm = T)))
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m1 misspecification use different variance
    p <- sapply(1:n, function(j) {
      if (length(unique(na.omit(score_exp[j,])))<2|length(unique(na.omit(score_back[j,])))<2) {
        return(1)} else {
          return(t.test(-log(score_exp[j,!(is.na(score_exp[j,]))]+0.01), -log(score_back[j,!(is.na(score_back[j,]))]+0.01), alternative = "greater", var.equal = T)$p.value)}})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method  = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method  = 'max', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method = 'diff', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})
    
    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pair_c = fdppow_pair_c[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pair_c = fdppow_pair_c[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')
  
  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr, 
                             q = rep(q_tot, length(methods)), 
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_multirep_pois_hete) = out_prop
saveRDS(sim_multirep_pois_hete, 'sim_multirep_pois_hete_prop_interesting.rds')

######################  Fig1 single replicates, generated from NB   #############
########################################################################################
nb.p = 1/2
##############outlier################
outlier_nb = function(n, size.nb){
  a = qnbinom(0.99, size.nb, nb.p )
  rtrunc(n, spec = 'nbinom', size = size.nb,prob = nb.p, a = a)
}

sim_multirep_nb_homo <- lapply(1:length(nonnull_prop), function(i_prop){
  
  prop =  nonnull_prop[i_prop]
  set.seed(1)
  trueidx = 1:(prop*n)
  mean_back = rep(20, n*3)
  mean_exp = rep(20, n*3)
  mean_exp[1:(length(trueidx)*3)] = rep(rnbinom(length(trueidx), size = 45, prob = nb.p), each = 3)
  
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnbinom(n*3, size = mean_back, prob = nb.p)
    score_exp <- rnbinom(n*3, size = mean_exp, prob = nb.p)

    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      score_exp0 <- na.omit(score_exp[j,])
      score_back0 <- na.omit(score_back[j,])
      if (length(score_exp0) == 0|length(score_back0) == 0) {return(1)} else {
        size_hat <- mean(c(score_back0, score_exp0))
        return(twosample_nb(sum(score_exp0)-sum(score_back0), size_hat*length(score_exp0), size_hat*length(score_back0), 0.5, 0.5))
      }
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,])) > 2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(1 - pnbinom(sum(score_exp[j,], na.rm = T)-1, size = length(na.omit(score_exp))/length(na.omit(score_back))*sum(score_back[j,], na.rm = T), prob = nb.p))
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use model misspecification
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "greater")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method  = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method  = 'max', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method = 'diff', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})
    
    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pair_c = fdppow_pair_c[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pair_c = fdppow_pair_c[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')
  
  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr, 
                             q = rep(q_tot, length(methods)), 
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_multirep_nb_homo) = out_prop
saveRDS(sim_multirep_nb_homo, 'sim_multirep_nb_homo_prop_interesting.rds')




sim_multirep_nb_hete <- lapply(1:length(nonnull_prop), function(i_prop){
  
  prop =  nonnull_prop[i_prop]
  set.seed(1)
  trueidx = 1:(prop*n)
  mean_back = rep(rnbinom(n, size = 20, prob = 0.5), each = 3)
  mean_exp = mean_back
  mean_exp[1:(length(trueidx)*3)] = rep(rnbinom(length(trueidx), size = 45, prob = nb.p), each = 3)
  
  fdppow_ls <- mclapply(1:iter, function (it){
    
    set.seed(it)
    score_back <- rnbinom(n*3, size = mean_back, prob = nb.p)
    score_exp <- rnbinom(n*3, size = mean_exp, prob = nb.p)
    score_back = matrix(score_back, ncol = 3, byrow = T)
    score_exp = matrix(score_exp, ncol = 3, byrow = T)
    
    ##### use BH pooled
    pval_pooled = compute_pooled_pval(back = rowMeans(score_back, na.rm = T), exp = rowMeans(score_exp, na.rm = T))
    discovery_ls = find_discovery_wpval(pval = pval_pooled, q = q_tot)
    fdppow_pooled = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = pval_pooled, q = q_tot)
    fdppow_pooled_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)

    ##### use BH paired correct specification
    p <- sapply(1:n, function(j) {
      score_exp0 <- na.omit(score_exp[j,])
      score_back0 <- na.omit(score_back[j,])
      if (length(score_exp0) == 0|length(score_back0) == 0) {return(1)} else {
        size_hat <- mean(c(score_back0, score_exp0))
        return(twosample_nb(sum(score_exp0)-sum(score_back0), size_hat*length(score_exp0), size_hat*length(score_back0), 0.5, 0.5))
      }
    })
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_c = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_c_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use BH paired m2 misspecification as variance 1
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,])) > 2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(1 - pnbinom(sum(score_exp[j,], na.rm = T)-1, size = length(na.omit(score_exp))/length(na.omit(score_back))*sum(score_back[j,], na.rm = T), prob = nb.p))
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_2as1 = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_2as1_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    
    p <- sapply(1:n, function(j){
      if (sum(is.na(score_back[j,]))>2|sum(is.na(score_exp[j,]))>2) {
        return(1)} else { return(poisson.test(c(sum(score_exp[j,], na.rm = T),sum(score_back[j,], na.rm = T)), c(sum(!is.na(score_exp[j,])), sum(!is.na(score_back[j,]))), alternative = "greater")$p.value)
        }})
    discovery_ls = find_discovery_wpval(pval = p, q = q_tot)
    fdppow_pair_mis = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    discovery_ls = find_discovery_wqval(pval = p, q = q_tot)
    fdppow_pair_mis_q = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### clipper 
    re_clipper <- clipper1sided(score_exp = score_exp, 
                                score_back = score_back, 
                                FDR = q_tot, 
                                importanceScore_method  = 'diff', 
                                FDR_control_method = 'BC',
                                ifpowerful = T)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper max BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method  = 'max', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clippermax_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    ##### use clipper diff BC
    re_clipper <- clipper1sided( score_exp = score_exp, 
                                 score_back = score_back, 
                                 FDR = q_tot, 
                                 importanceScore_method = 'diff', 
                                 FDR_control_method = 'BC',
                                 ifpowerful = F)
    discovery_ls = lapply(re_clipper$results,'[[', 'discovery')
    fdppow_clipperdiff_bc = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    

    ##### emp null
    t1 <- sapply(1:n, function(i){
      if(length(unique(na.omit(score_exp[i,])))<=1|length(na.omit(unique(score_back[i,])))<=1) {return(NA)}else {
        test.out <- t.test(score_exp[i,],score_back[i,], var.equal = TRUE)
        return(qnorm(pt(test.out$statistic, test.out$parameter)))}})
    
    loc.fdr <- try(locfdr(t1[!is.na(t1)],  plot = 0)$fdr, silent = T)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdremp = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    loc.fdr = localfdr_multrep_p_miss(exp = score_exp, back = score_back)
    discovery_ls = find_discovery_wlocfdr(loc.fdr = loc.fdr, q = q_tot)
    fdppow_locfdrper = compute_fdppow(discovery_ls = discovery_ls, trueidx = trueidx)
    
    fdp = list(pooled = fdppow_pooled[1,], 
               pair_2as1 = fdppow_pair_2as1[1,],
               pair_mis = fdppow_pair_mis[1,],
               pair_c = fdppow_pair_c[1,],
               pooled_q = fdppow_pooled_q[1,], 
               pair_2as1_q = fdppow_pair_2as1_q[1,],
               pair_mis_q = fdppow_pair_mis_q[1,],
               pair_c_q = fdppow_pair_c_q[1,],
               clipper = fdppow_clipper[1,],
               clippermax_bc = fdppow_clippermax_bc[1,],
               clipperdiff_bc = fdppow_clipperdiff_bc[1,],
               locfdremp = fdppow_locfdremp[1,],
               locfdrper = fdppow_locfdrper[1,])
    
    pow = list(pooled = fdppow_pooled[2,], 
               pair_2as1 = fdppow_pair_2as1[2,],
               pair_mis = fdppow_pair_mis[2,],
               pair_c = fdppow_pair_c[2,],
               pooled_q = fdppow_pooled_q[2,], 
               pair_2as1_q = fdppow_pair_2as1_q[2,],
               pair_mis_q = fdppow_pair_mis_q[2,],
               pair_c_q = fdppow_pair_c_q[2,],
               clipper = fdppow_clipper[2,],
               clippermax_bc = fdppow_clippermax_bc[2,],
               clipperdiff_bc = fdppow_clipperdiff_bc[2,],
               locfdremp = fdppow_locfdremp[2,],
               locfdrper = fdppow_locfdrper[2,])
    
    return(list(fdp = fdp, pow = pow))
  },mc.cores = ncores)
  
  fdr = lapply( fdppow_ls, '[[', 'fdp')
  pow = lapply(fdppow_ls, '[[', 'pow')
  
  methods = names(fdr[[1]])
  re_fdr = lapply(methods, function(m){
    temp = sapply(fdr, '[[',m)
    fdr = rowMeans(temp, na.rm = T)
    fdr_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(fdr = fdr, fdr_sd = fdr_sd))
  })
  names(re_fdr) =methods
  re_fdr = Reduce('rbind', re_fdr)
  re_fdr = cbind.data.frame( re_fdr, 
                             q = rep(q_tot, length(methods)), 
                             methods = rep(methods, each = length(q_tot))
  )
  re_pow = lapply(methods, function(m){
    temp = sapply(pow, '[[',m)
    pow = rowMeans(temp, na.rm = T)
    pow_sd = apply(temp, 1, sd, na.rm = T)
    return(cbind(pow = pow, pow_sd = pow_sd))
  })
  names(re_pow) =methods
  re_pow = Reduce('rbind', re_pow)
  re = cbind.data.frame(re_fdr, re_pow)
  return(re)
})

names(sim_multirep_nb_hete) = out_prop
saveRDS(sim_multirep_nb_hete, 'sim_multirep_nb_hete_prop_interesting.rds')

########################################################################################





